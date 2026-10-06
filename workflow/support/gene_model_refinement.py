#!/usr/bin/env python3
"""Restartable, opt-in synteny-guided coding model refinement.

Original annotations are immutable. Accepted predictions and representative
adoption are separate decisions; every effective sequence has an exact GFF path.
"""
import argparse
import bisect
import contextlib
import copy
import functools
import gzip
import hashlib
import importlib.metadata
import json
import math
import re
import resource
import shutil
import sqlite3
import sys
import time
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote

try:
    import rescue_gene_models as rescue
    from fasta_sequence_store import exclusive_lock, open_text
    from format_species_annotation.common import parse_gff_attributes
    from gene_model_catalog import build_catalog, indexed_genome, validate_candidate, write_catalog
    from gene_model_selection import pair_score, select_representatives
    from gene_model_species_profiles import parameters_for, read_profiles
    from gene_model_store import _connection as store_connection
    from gene_model_store import build_store, iter_loci, iter_locus_keys, load_locus, select_from_store
    from input_generation_array_state import atomic_json, digest, digest_paths
except ImportError:
    from . import rescue_gene_models as rescue
    from .fasta_sequence_store import exclusive_lock, open_text
    from .format_species_annotation.common import parse_gff_attributes
    from .gene_model_catalog import build_catalog, indexed_genome, validate_candidate, write_catalog
    from .gene_model_selection import pair_score, select_representatives
    from .gene_model_species_profiles import parameters_for, read_profiles
    from .gene_model_store import _connection as store_connection
    from .gene_model_store import build_store, iter_loci, iter_locus_keys, load_locus, select_from_store
    from .input_generation_array_state import atomic_json, digest, digest_paths

SCHEMA = 1
MAP_FIELDS = ('species', 'gene_id', 'candidate_id', 'source_transcript_id', 'status', 'score', 'margin', 'reason')
DEFAULTS = dict(policy='conserved', mode='conservative', min_margin=0.10, min_support=2,
                minimum_coverage=0.95, minimum_identity=0.50, max_intron=20000,
                padding=2000, max_interval=200000, candidate_limit=32)


# Invocation-local reuse avoids re-reading every catalog for every target worker.
# Fresh CLI boundaries and each publication still verify frozen source bytes.
_INVOCATION_CACHE = None


def invocation_cached(function):
    @functools.wraps(function)
    def wrapped(root, value, *args, **kwargs):
        if _INVOCATION_CACHE is None:
            return function(root, value, *args, **kwargs)
        semantic = args[:1] if function.__name__ in {'catalog_species', 'predict_species'} else ()
        if function.__name__ in {'select', 'catalog_index'}:
            semantic = (kwargs.get('predictions', args[0] if args else False),)
        key = (function.__name__, str(root), semantic)
        if key not in _INVOCATION_CACHE:
            _INVOCATION_CACHE[key] = function(root, value, *args, **kwargs)
        return _INVOCATION_CACHE[key]
    return wrapped


def implementation():
    modules = (sys.modules[build_catalog.__module__], sys.modules[select_representatives.__module__],
               sys.modules[rescue.validate_model.__module__], sys.modules[atomic_json.__module__],
               sys.modules[build_store.__module__])
    support = Path(__file__).parent
    dependencies = ('cds_model_normalisation.py', 'gff_feature_structure.py', 'gff_attribute_syntax.py',
                    'fasta_sequence_store.py', 'species_labeling.py', 'pairwise_synteny.py',
                    'representative_selection.py', 'rescue_anchor_admission.py', 'gene_model_species_profiles.py',
                    'format_species_writers.py', 'format_species_common.py', 'format_species_constants.py',
                    'format_species_provider_config.py', 'format_species_taxonomy.py')
    files = [support / name for name in dependencies] + list((support / 'format_species_annotation').rglob('*.py'))
    return {str(Path(m.__file__).name): digest(m.__file__) for m in modules} | {Path(__file__).name: digest(__file__)} | {str(p.relative_to(support)): digest(p) for p in files}


def dependency_identities():
    from Bio.Align import _pairwisealigner
    return {'versions': {name: importlib.metadata.version(name) for name in ('biopython', 'pysam')},
            'python': sys.version, 'sqlite': sqlite3.sqlite_version,
            'pairwise_aligner_sha256': digest(_pairwisealigner.__file__)}


def read_table(path):
    return rescue.table(path)


def plan(output, inputs=None, rescue_output=None, edges=None, rna=None, species_profiles=None, **parameters):
    root = Path(output).resolve()
    root.mkdir(parents=True, exist_ok=True)
    with exclusive_lock(root / '.plan.lock'):
        if rescue_output:
            if inputs:
                raise ValueError('Choose --inputs or --rescue-output; explicit inputs must not be silently ignored')
            anchor_root = Path(rescue_output).resolve()
            anchors = rescue.load(anchor_root)
            sources = copy.deepcopy(anchors['request']['sources'])
            files = [str(anchor_root / 'plan.json')]
            augmented = anchor_root / 'augmented'
            if augmented.exists():
                receipts = {}
                for name in anchors['species']:
                    directory = anchor_root / 'rescued' / name
                    if not rescue.verified(directory, rescue.rescue_key(anchor_root, anchors, name)):
                        raise ValueError('Missing-gene rescue incomplete or corrupt: ' + name)
                    receipts[name] = digest(directory / 'receipt.json')
                    files.append(str(directory / 'receipt.json'))
                augmented_key = {'plan': digest(anchor_root / 'plan.json'), 'rescue_receipts': receipts}
                if not rescue.verified(augmented, augmented_key):
                    raise ValueError('Augmented rescue inputs incomplete or corrupt')
                rows = read_table(augmented / 'inputs.tsv')
                if len(rows) != len(sources) or {row['species'] for row in rows} != set(sources):
                    raise ValueError('Augmented rescue input species disagree with the anchor plan')
                for row in rows:
                    source = sources[row['species']]
                    for field, role in (('cds', 'fasta'), ('gff', 'gff'), ('genome', 'genome')):
                        path = Path(row[field])
                        source[role] = str((path if path.is_absolute() else augmented / path).resolve())
                files += [str(augmented / 'receipt.json'), str(augmented / 'inputs.tsv')]
        else:
            anchor_root = None
            if not inputs:
                raise ValueError('Provide --inputs or --rescue-output')
            sources = {}
            for row in read_table(inputs):
                name = row['species']
                rescue.safe_token(name, 'species')
                if name in sources:
                    raise ValueError('Duplicate species: ' + name)
                base = Path(inputs).resolve().parent
                def resolve_path(value, base=base):
                    p = Path(value)
                    return str((p if p.is_absolute() else base / p).resolve())
                sources[name] = {'species': name, 'fasta': resolve_path(row.get('cds', row.get('fasta', ''))),
                                 'gff': resolve_path(row['gff']), 'genome': resolve_path(row['genome']),
                                 'genetic_code': int(row.get('genetic_code') or 1)}
            files = [str(Path(inputs).resolve())]
        if not sources:
            raise ValueError('No species')
        params = {**DEFAULTS, **parameters}
        if params['policy'] not in {'longest', 'conserved'} or params['mode'] not in {'off', 'audit', 'conservative'}:
            raise ValueError('Unknown selection policy/refinement mode')
        if any(not math.isfinite(float(params[k])) or not 0 <= float(params[k]) <= 1
               for k in ('min_margin', 'minimum_coverage', 'minimum_identity')):
            raise ValueError('Invalid probability/margin')
        if any(int(params[k]) < 1 for k in ('min_support', 'max_intron', 'max_interval', 'candidate_limit')) or params['padding'] < 0:
            raise ValueError('Invalid resource/evidence bound')
        files += [s[k] for s in sources.values() for k in ('fasta', 'gff', 'genome')]
        profiles = (read_profiles(species_profiles, sources) if species_profiles else
                    copy.deepcopy(anchors['request'].get('species_profiles', {})) if anchor_root else {})
        for optional in (edges, rna, species_profiles):
            if optional:
                files.append(str(Path(optional).resolve()))
        tool = shutil.which('miniprot') if params['mode'] != 'off' else None
        request = {'schema': SCHEMA, 'sources': sources, 'files': digest_paths(files), 'parameters': params,
                   'implementation': implementation(), 'dependencies': dependency_identities(), 'rescue_output': str(anchor_root) if anchor_root else None,
                   'edges': str(Path(edges).resolve()) if edges else None,
                   'rna': str(Path(rna).resolve()) if rna else None, 'species_profiles': profiles,
                   'miniprot': {'path': tool, 'sha256': digest(tool)} if tool else None}
        if not edges and not anchor_root:
            raise ValueError('A frozen synteny rescue plan or explicit trusted correspondence table is required')
        if params['mode'] != 'off' and not tool:
            raise ValueError('miniprot is required for prediction')
        path = root / 'plan.json'
        if path.exists():
            previous = json.loads(path.read_text())
            if previous['request'] != request:
                raise ValueError('Frozen refinement plan differs; use a new output directory')
            return previous
        value = {'request': request, 'species': sorted(sources)}
        atomic_json(path, value, immutable=True)
        return value


def load(root, names=None):
    root = Path(root)
    value = json.loads((root / 'plan.json').read_text())
    request = value['request']
    if request['schema'] != SCHEMA or request['implementation'] != implementation() or request['dependencies'] != dependency_identities():
        raise ValueError('Refinement implementation changed; use a new output directory')
    species_paths = {source[k] for source in request['sources'].values() for k in ('fasta', 'gff', 'genome')}
    files = request['files'] if names is None else {p: value for p, value in request['files'].items() if p not in species_paths}
    if names is not None:
        files.update({request['sources'][n][k]: request['files'][request['sources'][n][k]] for n in names for k in ('fasta', 'gff', 'genome')})
    if digest_paths(files) != files:
        raise ValueError('Frozen refinement input changed')
    tool = request['miniprot']
    if tool and digest(tool['path']) != tool['sha256']:
        raise ValueError('Frozen predictor changed')
    return value


def stage(root, relative, key, builder, names=None):
    dependencies = key.get('dependencies', {})
    snapshots = {root / 'catalog' / name / 'receipt.json': expected
                 for name, expected in dependencies.get('catalog', {}).items()}
    snapshots.update({root / directory / 'receipt.json': expected
                      for directory, expected in dependencies.get('index', {}).items()})
    for field, directory in (('correspondence', 'correspondence'), ('initial', 'selection_initial'), ('selection', 'selection_final')):
        if field in dependencies:
            snapshots[root / directory / 'receipt.json'] = dependencies[field]
    if 'anchors' in dependencies or 'prepared' in dependencies:
        anchor_root = Path(json.loads((root / 'plan.json').read_text())['request']['rescue_output'])
        for field, directory in (('anchors', 'synteny'), ('prepared', 'prepared')):
            snapshots.update({anchor_root / directory / name / 'receipt.json': expected
                              for name, expected in dependencies.get(field, {}).items()})
    snapshots.update({root / 'predictions' / name / 'receipt.json': expected
                      for name, expected in dependencies.get('predictions', {}).items()})
    snapshots = {str(path): expected for path, expected in snapshots.items()}

    def guard():
        load(root, names)
        rescue.require_same_key(snapshots, digest_paths(snapshots))
        for path in snapshots:
            directory = Path(path).parent
            receipt = json.loads(Path(path).read_text())
            if not rescue.verified(directory, receipt['key']):
                raise ValueError('Stage dependency content changed: ' + str(directory))

    def measured(tmp):
        start, cpu = time.perf_counter(), time.process_time()
        children_before = resource.getrusage(resource.RUSAGE_CHILDREN)
        builder(tmp)
        children = resource.getrusage(resource.RUSAGE_CHILDREN)
        self_usage = resource.getrusage(resource.RUSAGE_SELF)
        rss_factor = 1024 ** 2 if sys.platform == 'darwin' else 1024
        atomic_json(tmp / 'performance.json', {'builder_wall_seconds': time.perf_counter() - start,
                    'builder_process_cpu_seconds': time.process_time() - cpu,
                    'builder_child_cpu_seconds': children.ru_utime + children.ru_stime - children_before.ru_utime - children_before.ru_stime,
                    'process_lifetime_peak_rss_mib': self_usage.ru_maxrss / rss_factor,
                    'child_lifetime_peak_rss_mib': children.ru_maxrss / rss_factor,
                    'note': 'Builder excludes source/receipt verification; RSS is a lifetime high-water mark, not an isolated stage peak.'})
    return rescue.stage(root, Path(relative), {'plan': digest(root / 'plan.json'), **key}, measured, guard)


@invocation_cached
def catalog_species(root, value, name):
    source = value['request']['sources'][name]
    def build(tmp):
        catalog = build_catalog(name, source['fasta'], source['gff'], source['genome'], source['genetic_code'])
        if value['request']['rna']:
            rna = rna_path_index(read_table(value['request']['rna']))
            for gene in catalog['loci']:
                for candidate in gene['candidates']:
                    paths = rna_support(candidate, name, rna)
                    candidate['quality']['rna_supported'] = bool(paths)
                    candidate.setdefault('support', {})['rna_paths'] = paths
        write_catalog(catalog, tmp)
        # The integration contract is independent of catalog writer filenames.
        atomic_json(tmp / 'catalog.json', catalog)
    return stage(root, Path('catalog') / name, {}, build, [name])


@invocation_cached
def catalog_index(root, value, predictions=False):
    directories = [catalog_species(root, value, name) for name in value['species']]
    dependencies = {'catalog': {n: digest(root / 'catalog' / n / 'receipt.json') for n in value['species']}}
    if predictions:
        for name in value['species']:
            predict_species(root, value, name)
        dependencies['predictions'] = {n: digest(root / 'predictions' / n / 'receipt.json') for n in value['species']}
    def build(tmp):
        db = tmp / 'loci.sqlite3'
        summary = build_store(directories, db)
        summary['database'] = str((root / ('catalog_index_final' if predictions else 'catalog_index') / 'loci.sqlite3').resolve())
        if predictions:
            accepted_count = 0
            with sqlite3.connect(db) as connection:
                for name in value['species']:
                    additions = json.loads((root / 'predictions' / name / 'predictions.json').read_text())
                    by_gene = defaultdict(list)
                    for model in additions:
                        if model['status'] == 'accepted':
                            by_gene[model['gene_id']].append(model['candidate'])
                    for gene, candidates in by_gene.items():
                        locus = json.loads(connection.execute('SELECT json FROM loci WHERE species=? AND gene_id=?', (name, gene)).fetchone()[0])
                        locus['candidates'].extend(candidates)
                        for candidate in candidates:
                            connection.execute('INSERT INTO candidate_owners(species,candidate_id,gene_id) VALUES(?,?,?)', (name, candidate['candidate_id'], gene))
                        serialized = json.dumps(locus, sort_keys=True)
                        connection.execute('UPDATE loci SET json=? WHERE species=? AND gene_id=?', (serialized, name, gene))
                        accepted_count += len(candidates)
                        summary['max_locus_json_bytes'] = max(summary['max_locus_json_bytes'], len(serialized.encode()))
            summary['accepted_prediction_count'] = accepted_count
            summary['candidate_count'] += accepted_count
        atomic_json(tmp / 'summary.json', summary)
    directory = stage(root, 'catalog_index_final' if predictions else 'catalog_index', {'dependencies': dependencies}, build)
    return directory / 'loci.sqlite3'


def infer_flanked_loci(db, species_a, species_b, anchor_pairs, position_index, params):
    """Nominate existing nonanchor loci only between two conserved anchors.

    A rank correspondence requires equal bounded counts and independent coding
    similarity. Copy conflicts remain ambiguous when all block evidence merges.
    """
    result = []
    def interior(species, left, right):
        positions, tracks = position_index[species]
        a, b = positions[left], positions[right]
        if a['seqid'] != b['seqid']:
            return []
        low, high = sorted((a, b), key=lambda g: g['start'])
        start, end = low['end'], high['start']
        if end <= start or end - start > params['max_interval']:
            return []
        starts, genes = tracks[a['seqid']]
        found = [g for g in genes[bisect.bisect_left(starts, start):bisect.bisect_left(starts, end)] if g['end'] <= end]
        return found if a['start'] < b['start'] else list(reversed(found))
    for left, right in zip(anchor_pairs, anchor_pairs[1:], strict=False):
        if None in (*left, *right):
            continue
        a = interior(species_a, left[0], right[0])
        b = interior(species_b, left[1], right[1])
        if not a or len(a) != len(b) or len(a) > 8:
            continue
        for ga, gb in zip(a, b, strict=True):
            la, lb = load_locus(db, species_a, ga['gene_id']), load_locus(db, species_b, gb['gene_id'])
            ca = max(la['candidates'], key=lambda c: c.get('corrected_cds_length', len(c['cds'])))
            cb = max(lb['candidates'], key=lambda c: c.get('corrected_cds_length', len(c['cds'])))
            score = pair_score(ca, cb)
            if score.bounded or score.identity < params['minimum_identity'] or min(score.coverage_a, score.coverage_b) < 0.75:
                continue
            result.append(dict(species_a=species_a, gene_a=ga['gene_id'], species_b=species_b, gene_b=gb['gene_id'],
                               weight=0.8, ambiguous=False,
                               evidence={'kind': 'two_flanking_anchors', 'left': list(left), 'right': list(right),
                                         'identity': score.identity, 'coverage_a': score.coverage_a, 'coverage_b': score.coverage_b}))
    return result


@invocation_cached
def correspondence(root, value, cpus=1, comparison_cache=None):
    db = catalog_index(root, value)
    request = value['request']
    dependencies = {'catalog': {n: digest(root / 'catalog' / n / 'receipt.json') for n in value['species']},
                    'index': {db.parent.name: digest(db.parent / 'receipt.json')}}
    anchor_root = Path(request['rescue_output']) if request['rescue_output'] else None
    jobs = []
    anchor_plan = None
    if anchor_root and not request['edges']:
        anchor_plan = rescue.load(anchor_root)
        jobs = [j for j in anchor_plan['synteny_jobs'] if j['a'] != j['b']]
        for job in jobs:
            rescue.synteny(anchor_root, anchor_plan, job['index'], cpus, comparison_cache)
        dependencies['anchors'] = {j['id']: digest(anchor_root / 'synteny' / j['id'] / 'receipt.json') for j in jobs}
        dependencies['prepared'] = {}
        anchor_hash = digest(anchor_root / 'plan.json')
        for name in value['species']:
            directory = rescue.prepared(anchor_root, anchor_plan, name)
            if not rescue.verified(directory, {'plan': anchor_hash, 'species': name}):
                raise ValueError('Prepared annotation incomplete or corrupted: ' + name)
            dependencies['prepared'][name] = digest(directory / 'receipt.json')
    def build(tmp):
        edges = []
        valid = defaultdict(set)
        for species, gene in iter_locus_keys(db):
            valid[species].add(gene)
        if request['edges']:
            for row in read_table(request['edges']):
                a, b, ga, gb = (row[k] for k in ('species_a', 'species_b', 'gene_a', 'gene_b'))
                if a == b or ga not in valid.get(a, set()) or gb not in valid.get(b, set()):
                    raise ValueError('Correspondence references absent locus or same species: ' + str(row))
                weight = float(row.get('weight') or 1)
                if not math.isfinite(weight) or weight <= 0 or weight > 1:
                    raise ValueError('Invalid correspondence weight')
                edges.append(dict(species_a=a, gene_a=ga, species_b=b, gene_b=gb, weight=weight,
                                  evidence=row.get('evidence') or 'user_frozen_synteny',
                                  ambiguous=str(row.get('ambiguous', '')).lower() in {'1', 'true', 'yes'}))
        else:
            aliases = {}
            for n in value['species']:
                lookup = {}
                for g in iter_loci(db, n):
                    for token in {g['gene_id'], g['gene_id'].removeprefix(n + '_'),
                                  *(c['source_gene_id'] for c in g['candidates']), *(c.get('gene_token', '') for c in g['candidates'])}:
                        if token in lookup and lookup[token] != g['gene_id']:
                            raise ValueError('Ambiguous gene alias: ' + token)
                        lookup[token] = g['gene_id']
                aliases[n] = {r['jcvi_id']: lookup.get(r['locus_id']) or lookup.get(r['original_id'])
                              for r in read_table(anchor_root / 'prepared' / n / 'genes.id_map.tsv') if r['status'] == 'selected'}
            position_index = {}
            for n in value['species']:
                positions, tracks = {}, defaultdict(list)
                for g in iter_loci(db, n):
                    blocks = [b for c in g['candidates'] for b in c['blocks']]
                    if not blocks:
                        continue
                    position = dict(gene_id=g['gene_id'], seqid=g['seqid'], start=min(b[0] for b in blocks), end=max(b[1] for b in blocks))
                    positions[g['gene_id']] = position
                    tracks[g['seqid']].append(position)
                ordered = {seqid: sorted(genes, key=lambda g: g['start']) for seqid, genes in tracks.items()}
                position_index[n] = positions, {seqid: ([g['start'] for g in genes], genes) for seqid, genes in ordered.items()}
            with store_connection(db) as connection:
                for job in jobs:
                    if not rescue.verified(anchor_root / 'synteny' / job['id'], rescue.comparison_key(anchor_root, job)):
                        raise ValueError('Unverified synteny comparison')
                    blocks = json.loads((anchor_root / 'synteny' / job['id'] / 'blocks.json').read_text())
                    for block in blocks:
                        anchor_pairs = [(aliases[job['a']].get(a), aliases[job['b']].get(b)) for a, b in block]
                        edges.extend(infer_flanked_loci(connection, job['a'], job['b'], anchor_pairs, position_index, request['parameters']))
                        for i, (a, b) in enumerate(block):
                            ga, gb = aliases[job['a']].get(a), aliases[job['b']].get(b)
                            if not ga or not gb:
                                continue
                            # Boundary anchors lack two independent flanks and are audit-only.
                            edges.append(dict(species_a=job['a'], gene_a=ga, species_b=job['b'], gene_b=gb,
                                              weight=1.0, evidence={'comparison': job['id'], 'anchors': len(block),
                                                                  'left_flank': i > 0, 'right_flank': i + 1 < len(block)},
                                              ambiguous=i == 0 or i + 1 == len(block)))
        unique = {}
        for e in edges:
            key = tuple(sorted(((e['species_a'], e['gene_a']), (e['species_b'], e['gene_b']))))
            old = unique.get(key)
            if old is None or old['ambiguous'] and not e['ambiguous']:
                unique[key] = e
        edges = [unique[k] for k in sorted(unique)]
        neighbors = defaultdict(set)
        for e in edges:
            if not e['ambiguous']:
                neighbors[(e['species_a'], e['gene_a'], e['species_b'])].add(e['gene_b'])
                neighbors[(e['species_b'], e['gene_b'], e['species_a'])].add(e['gene_a'])
        for e in edges:
            if len(neighbors[e['species_a'], e['gene_a'], e['species_b']]) > 1 or len(neighbors[e['species_b'], e['gene_b'], e['species_a']]) > 1:
                e['ambiguous'] = True
        atomic_json(tmp / 'edges.json', edges)
        atomic_json(tmp / 'summary.json', {'edges': len(edges), 'ambiguous': sum(e['ambiguous'] for e in edges)})
    return stage(root, 'correspondence', {'dependencies': dependencies}, build)


@invocation_cached
def select(root, value, predictions=False):
    db = catalog_index(root, value, predictions)
    directory = correspondence(root, value)
    edges = json.loads((directory / 'edges.json').read_text())
    params = value['request']['parameters']
    dependencies = {'correspondence': digest(directory / 'receipt.json'),
                    'catalog': {n: digest(root / 'catalog' / n / 'receipt.json') for n in value['species']},
                    'index': {db.parent.name: digest(db.parent / 'receipt.json')}}
    if predictions:
        dependencies['predictions'] = {n: digest(root / 'predictions' / n / 'receipt.json') for n in value['species']}
    def build(tmp):
        selection = select_from_store(db, edges, policy=params['policy'], min_margin=params['min_margin'],
                                           min_support=params['min_support'], candidate_limit=params['candidate_limit'])
        atomic_json(tmp / 'selection.json', selection)
        rescue.write_tsv(tmp / 'representative_map.tsv', MAP_FIELDS,
                         [[r.get(k, '') for k in MAP_FIELDS] for r in selection['selections']])
    return stage(root, 'selection_final' if predictions else 'selection_initial', {'dependencies': dependencies}, build)


def rna_path_index(rows):
    """Index complete coding paths once; junctions alone never supply support."""
    if isinstance(rows, dict):
        return rows
    index = defaultdict(set)
    for row in rows:
        blocks = json.loads(row['cds_blocks'])
        if not isinstance(blocks, list) or not blocks or row['strand'] not in {'+', '-'}:
            raise ValueError('Invalid whole coding RNA path')
        if any(not isinstance(b, list) or len(b) != 2 or any(type(v) is not int for v in b) or b[0] < 0 or b[1] <= b[0] for b in blocks):
            raise ValueError('RNA coding intervals require integer zero-based half-open coordinates')
        if any((a[1] > b[0] if row['strand'] == '+' else b[1] > a[0]) for a, b in zip(blocks, blocks[1:], strict=False)):
            raise ValueError('RNA coding intervals must be disjoint and in transcript order')
        if int(row.get('count') or 1) > 0:
            key = row['species'], row['seqid'], row['strand'], tuple(map(tuple, blocks))
            index[key].add(row.get('transcript_id') or 'rna_path')
    return {key: sorted(values) for key, values in index.items()}


def rna_support(candidate, name, rows):
    key = name, candidate['seqid'], candidate['strand'], tuple(tuple(b[:2]) for b in candidate['blocks'])
    return list(rna_path_index(rows).get(key, ()))


def annotation_ownership_spans(gff_path, catalog):
    owners = {c[key]: g['gene_id'] for g in catalog['loci'] for c in g['candidates']
              for key in ('source_gene_id', 'source_transcript_id', 'gene_token') if c.get(key)}
    spans = {}
    with open_text(Path(gff_path)) as handle:
        for line in handle:
            if line.rstrip('\r\n') == '##FASTA':
                break
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip().split('\t')
            if len(f) != 9:
                raise ValueError('Invalid annotation row')
            parsed = parse_gff_attributes(f[8])
            attr = {key: values[0] for key, values in parsed.items() if values}
            identifier = attr.get('ID', attr.get('transcript_id', attr.get('gene_id', '')))
            declared = parsed.get('Parent', ()) + parsed.get('gene_id', ()) + parsed.get('transcript_id', ())
            if (f[2].lower().endswith(('gene', 'rna', 'transcript')) or f[2] in {'CDS', 'exon'}
                    or identifier in owners or parsed.get('gene_id') or parsed.get('transcript_id')
                    or any(token in owners for token in declared)):
                fallback = attr.get('gene_id') or next(iter(parsed.get('Parent', ())), identifier)
                owner = owners.get(identifier) or next((owners[p] for p in declared if p in owners), 'annotation:' + fallback)
                key = f[0], owner
                start, end = int(f[3]) - 1, int(f[4])
                previous = spans.get(key, (start, end))
                spans[key] = min(previous[0], start), max(previous[1], end)
    return [{'seqid': seqid, 'start': start, 'end': end, 'gene_id': owner}
            for (seqid, owner), (start, end) in sorted(spans.items())]


def classify_predictions(models, catalog, edges, params, rna_rows, genome_hash):
    """Require target ownership, intact genomic ORF and independent support."""
    rna_rows = rna_path_index(rna_rows)
    loci = {g['gene_id']: g for g in catalog['loci']}
    # Freeze the same species/locus predicates once. Scanning a cohort-wide
    # correspondence graph for every predicted coding path is quadratic.
    ambiguous_loci, trusted_donors = set(), defaultdict(set)
    for edge in edges:
        for side, other in (('a', 'b'), ('b', 'a')):
            if edge['species_' + side] != catalog['species']:
                continue
            gene = edge['gene_' + side]
            if edge.get('ambiguous'):
                ambiguous_loci.add(gene)
            else:
                trusted_donors[gene].add(edge['species_' + other])
    spans = []
    for g in loci.values():
        blocks = [b for c in g['candidates'] for b in c['blocks']]
        if blocks:
            spans.append((g['seqid'], min(b[0] for b in blocks), max(b[1] for b in blocks), g['gene_id']))
    spans.extend((r['seqid'], r['start'], r['end'], r['gene_id']) for r in catalog.get('annotation_spans', []))
    tracks = defaultdict(list)
    for seqid, start, end, owner in spans:
        tracks[seqid].append((start, end, owner))
    indexes = {}
    for seqid, rows in tracks.items():
        ordered = sorted(rows)
        maxima, maximum = [], 0
        for _, end, _ in ordered:
            maximum = max(maximum, end)
            maxima.append(maximum)
        indexes[seqid] = [r[0] for r in ordered], ordered, maxima
    def foreign_overlap(seqid, start, end, owner):
        starts, rows, maxima = indexes.get(seqid, ([], [], []))
        i = bisect.bisect_left(starts, end) - 1
        while i >= 0 and maxima[i] > start:
            a, b, gene = rows[i]
            if gene != owner and a < end and b > start:
                return True
            i -= 1
        return False
    grouped = {}
    for m in models:
        g = loci[m['gene_id']]
        candidate = {'seqid': m['seqid'], 'strand': m['strand'], 'blocks': m['cds'], 'cds': m['sequence'],
                     'origin': 'predicted', 'source_gene_id': g['candidates'][0]['source_gene_id'],
                     'gene_token': g.get('gene_token', g['candidates'][0].get('gene_token', g['gene_id'])), 'junctions': []}
        shape = (candidate['seqid'], candidate['strand'], tuple(tuple(b) for b in candidate['blocks']))
        original_shapes = {(c.get('seqid', g['seqid']), c.get('strand', g['strand']), tuple(tuple(b) for b in c['blocks'])) for c in g['candidates']}
        # Alignment quality belongs to this donor, not the shared coding path.
        # Failed alignments cannot veto a supported path or count as support.
        alignment_problems = sorted(set(m['problems']) & {'low_coverage', 'low_identity'})
        problems = sorted(set(m['problems']) - set(alignment_problems))
        if g.get('ambiguous_coordinates'):
            problems.append('ambiguous_locus_coordinates')
        if m['gene_id'] in ambiguous_loci:
            problems.append('ambiguous_locus_correspondence')
        if shape in original_shapes:
            problems.append('existing_coding_path')
        if candidate['strand'] != g['strand'] or candidate['seqid'] != g['seqid']:
            problems.append('wrong_locus_strand')
        start, end = min((b[0] for b in m['cds']), default=0), max((b[1] for b in m['cds']), default=0)
        if foreign_overlap(m['seqid'], start, end, m['gene_id']):
            problems.append('overlap_other_locus')
        owners = [c for c in g['candidates'] if any(a[0] < b[1] and a[1] > b[0] for a in c['blocks'] for b in m['cds'])]
        if not owners:
            problems.append('no_overlap_owned_locus')
        # Protected exceptions cannot become normal genes through donor transfer.
        if any(c.get('quality', {}).get('annotated_pseudogene') or c.get('quality', {}).get('translation_exception') or c.get('quality', {}).get('annotated_exception') or c.get('quality', {}).get('sequence_exception') for c in g['candidates']):
            problems.append('protected_annotation_exception')
        identity = {'assembly': genome_hash, 'gene': m['gene_id'], 'shape': shape, 'cds': candidate['cds']}
        identifier = catalog['species'] + '_ggrefine_' + hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:20]
        candidate.update(candidate_id=identifier, source_transcript_id=identifier,
                         coding_key=hashlib.sha256(json.dumps(shape).encode()).hexdigest())
        candidate['quality'] = validate_candidate(candidate, catalog.get('genetic_code', 1))
        if not candidate['quality']['valid_orf']:
            problems.append('invalid_predicted_coding_path')
        if m['donor_species'] not in trusted_donors[m['gene_id']]:
            problems.append('untrusted_donor_correspondence')
        candidate['quality']['rna_supported'] = bool(rna_support(candidate, catalog['species'], rna_rows))
        if not candidate['quality']['rna_supported'] and any(c.get('quality', {}).get('sequence_mismatch') for c in owners):
            problems.append('unresolved_source_sequence_mismatch')
        key = m['gene_id'], shape
        if key not in grouped:
            grouped[key] = {'gene_id': m['gene_id'], 'candidate': candidate, 'problems': problems,
                            'donors': [], 'alignments': [], 'rna_paths': rna_support(candidate, catalog['species'], rna_rows)}
        row = grouped[key]
        supported = not alignment_problems and 'untrusted_donor_correspondence' not in problems
        if supported:
            row['donors'].append(m['donor_species'])
        row['alignments'].append({'donor_species': m['donor_species'], 'donor_candidate': m['donor_candidate'],
                                  'identity': m['identity'], 'coverage': m['coverage'],
                                  'supports_path': supported, 'problems': alignment_problems})
        row['problems'] = sorted(set(row['problems'] + problems))
    result = []
    for row in grouped.values():
        row['donors'] = sorted(set(row['donors']))
        if not row['donors']:
            row['problems'].append('no_qualifying_donor_alignment')
        if len(row['donors']) < params['min_support'] and not row['rna_paths']:
            row['problems'].append('insufficient_independent_support')
        row['status'] = 'accepted' if not row['problems'] and params['mode'] == 'conservative' else 'proposal'
        valid_original = any(c['quality'].get('valid_orf') for c in loci[row['gene_id']]['candidates'])
        row['change_type'] = 'isoform_addition' if valid_original else 'model_revision'
        row['evidence_class'] = 'rna_path_supported' if row['rna_paths'] else 'homology_only_predicted'
        row['candidate']['quality']['representative_eligible'] = row['status'] == 'accepted' and (bool(row['rna_paths']) or not valid_original)
        row['candidate']['support'] = {'donors': row['donors'], 'rna_paths': row['rna_paths'], 'class': row['evidence_class']}
        result.append(row)
    active = []
    ordered = sorted(result, key=lambda r: (r['candidate']['seqid'], min(b[0] for b in r['candidate']['blocks'])))
    for row in ordered:
        c = row['candidate']
        if not c['quality'].get('valid_orf') or set(row['problems']) - {'insufficient_independent_support'}:
            continue
        start, end = min(b[0] for b in c['blocks']), max(b[1] for b in c['blocks'])
        active = [(s, e, other) for s, e, other in active if other['candidate']['seqid'] == c['seqid'] and e > start]
        for _, _, other in active:
            if other['gene_id'] != row['gene_id']:
                for conflicting in (row, other):
                    conflicting['status'] = 'proposal'
                    if 'competing_predictions_overlap' not in conflicting['problems']:
                        conflicting['problems'].append('competing_predictions_overlap')
                    conflicting['candidate']['quality']['representative_eligible'] = False
        active.append((start, end, row))
    return sorted(result, key=lambda r: (r['gene_id'], r['candidate']['candidate_id']))


@invocation_cached
def predict_species(root, value, name, cpus=1):
    initial = select(root, value)
    cat_dir = catalog_species(root, value, name)
    db = catalog_index(root, value)
    catalog = json.loads((cat_dir / 'catalog_metadata.json').read_text())
    catalog['loci'] = list(iter_loci(db, name))
    corr = correspondence(root, value)
    edges = json.loads((corr / 'edges.json').read_text())
    params, source = parameters_for(value['request'], name), value['request']['sources'][name]
    dependencies = {'initial': digest(initial / 'receipt.json'), 'catalog': {n: digest(root / 'catalog' / n / 'receipt.json') for n in value['species']},
                    'correspondence': digest(corr / 'receipt.json'),
                    'index': {db.parent.name: digest(db.parent / 'receipt.json')}}
    def build(tmp):
        if params['mode'] == 'off':
            atomic_json(tmp / 'predictions.json', [])
            return
        catalog['annotation_spans'] = annotation_ownership_spans(source['gff'], catalog)
        loci = {(name, g['gene_id']): g for g in catalog['loci']}
        decisions = {(r['species'], r['gene_id']): r for r in json.loads((initial / 'selection.json').read_text())['selections']}
        proteins, windows, queries = {}, defaultdict(list), {}
        proposals = []
        with store_connection(db) as connection:
            for edge in edges:
                if edge['ambiguous'] or name not in {edge['species_a'], edge['species_b']}:
                    continue
                target_side = 'a' if edge['species_a'] == name else 'b'
                donor_side = 'b' if target_side == 'a' else 'a'
                target_gene = edge['gene_' + target_side]
                donor = edge['species_' + donor_side]
                target = loci[name, target_gene]
                donor_locus = load_locus(connection, donor, edge['gene_' + donor_side])
                if target.get('ambiguous_coordinates') or donor_locus.get('ambiguous_coordinates') or decisions[name, target_gene]['status'] == 'ambiguous_correspondence' or decisions[donor, donor_locus['gene_id']]['status'] == 'ambiguous_correspondence':
                    proposals.append({'gene_id': target_gene, 'status': 'ambiguous_locus_ownership', 'donor': donor})
                    continue
                if not any(c['cds'] and c['blocks'] for c in target['candidates']):
                    proposals.append({'gene_id': target_gene, 'status': 'unreconstructible_original', 'donor': donor})
                    continue
                if len(target['candidates']) > params['candidate_limit'] or len(donor_locus['candidates']) > params['candidate_limit']:
                    proposals.append({'gene_id': target_gene, 'status': 'candidate_limit_exceeded', 'donor': donor})
                    continue
                baseline = next(c for c in target['candidates'] if c['candidate_id'] == decisions[name, target_gene]['candidate_id'])
                # Intact concordant single coding paths need no predictor process.
                for candidate in donor_locus['candidates'][:params['candidate_limit']]:
                    if not candidate['quality'].get('usable') or not candidate['protein']:
                        continue
                    discrepancy = not baseline['quality'].get('valid_orf') or len(target['candidates']) < len(donor_locus['candidates']) or abs(len(candidate['protein']) - len(baseline['protein'])) > max(3, len(candidate['protein']) * .05)
                    if not discrepancy:
                        score = pair_score(baseline, candidate)
                        discrepancy = not score.bounded and score.identity >= params['minimum_identity'] and (min(score.coverage_a, score.coverage_b) < params['minimum_coverage'] or (len(baseline['blocks']) > 1 and len(candidate['blocks']) > 1 and score.junction_similarity < 0.8))
                    if not discrepancy:
                        continue
                    represented = False
                    for known in target['candidates']:
                        if not known['quality'].get('valid_orf'):
                            continue
                        match = pair_score(known, candidate)
                        if not match.bounded and min(match.coverage_a, match.coverage_b) >= params['minimum_coverage'] and match.identity >= params['minimum_identity'] and match.junction_similarity >= 0.95:
                            represented = True
                            break
                    if represented:
                        continue
                    blocks = [b for c in target['candidates'] for b in c['blocks']]
                    start, end = max(0, min(b[0] for b in blocks) - params['padding']), max(b[1] for b in blocks) + params['padding']
                    if end - start > params['max_interval']:
                        proposals.append({'gene_id': target_gene, 'status': 'interval_too_large', 'donor': donor})
                        continue
                    qid = 'q' + hashlib.sha256((target_gene + donor + candidate['candidate_id']).encode()).hexdigest()[:24]
                    region = {'id': qid, 'seqid': target['seqid'], 'start': start, 'end': end, 'query': candidate['candidate_id'],
                              'donor': donor, 'gene_id': target_gene, 'donor_candidate': candidate['candidate_id']}
                    queries[qid] = region
                    proteins.setdefault(donor, {})[candidate['candidate_id']] = candidate['protein']
                    windows[target['seqid'], start, end].append(region)
        validated = []
        # Nominate first: a species with no search windows needs no genome I/O.
        with indexed_genome(source['genome']) if windows else contextlib.nullcontext() as genome:
            bounded = {}
            for (seqid, start, end), regions in windows.items():
                end = min(end, genome.get_reference_length(seqid))
                for region in regions:
                    region['end'] = end
                bounded[seqid, start, end] = regions
            models = rescue.search_intervals(tmp, bounded, proteins, genome, source['genetic_code'], params['max_intron'], cpus) if windows else []
            for model in models:
                region = queries[model['query']]
                model['seqid'] = region['seqid']
                model['cds'] = [[a + region['start'], b + region['start'], p] for a, b, p in model['cds']]
                model.update(gene_id=region['gene_id'], donor_species=region['donor'], donor_candidate=region['donor_candidate'])
                checked = rescue.validate_model(model, genome, source['genetic_code'], params)
                if any(a < region['start'] or b > region['end'] for a, b, _ in checked['cds']):
                    checked['problems'].append('outside_search_window')
                validated.append(checked)
        rna = read_table(value['request']['rna']) if value['request']['rna'] else []
        predictions = classify_predictions(validated, catalog, edges, params, rna, value['request']['files'][source['genome']])
        atomic_json(tmp / 'predictions.json', predictions)
        atomic_json(tmp / 'nominations.json', {'queries': list(queries.values()), 'proposals': proposals})
        atomic_json(tmp / 'summary.json', {'queries': len(queries), 'windows': len(windows), 'alignments': len(validated),
                                         'accepted': sum(r['status'] == 'accepted' for r in predictions),
                                         'proposals': sum(r['status'] != 'accepted' for r in predictions) + len(proposals)})
    donor_names = {name}
    for edge in edges:
        if not edge['ambiguous'] and name in {edge['species_a'], edge['species_b']}:
            donor_names.update((edge['species_a'], edge['species_b']))
    return stage(root, Path('predictions') / name, {'dependencies': dependencies}, build, sorted(donor_names))


def candidate_gff(gene, candidate, name, *, gene_id=None, gene_token=None):
    blocks = candidate['blocks']
    def attr(text):
        return quote(str(text), safe='._-:')
    gid, tid = attr(gene['gene_id'] if gene_id is None else gene_id), attr(candidate['source_transcript_id'])
    token = attr(gene['gene_id'] if gene_token is None else gene_token)
    start, end = min(b[0] for b in blocks), max(b[1] for b in blocks)
    seqid, strand = candidate.get('seqid', gene['seqid']), candidate.get('strand', gene['strand'])
    def row(kind, a, b, phase, attributes):
        return f'{seqid}\tGeneGalleon\t{kind}\t{a + 1}\t{b}\t.\t{strand}\t{phase}\t{attributes}\n'
    lines = [row('gene', start, end, '.', 'ID=' + gid), row('mRNA', start, end, '.', f'ID={tid};Parent={gid};gene_id={token}')]
    for i, (a, b, phase) in enumerate(blocks, 1):
        lines += [row('exon', a, b, '.', f'ID={tid}.exon{i};Parent={tid}'),
                  row('CDS', a, b, phase, f'ID={tid}.cds{i};Parent={tid};gene_id={token}')]
    return ''.join(lines)


def selected_gff_rows(original, transcript_ids, *, retain_all=False):
    """Canonicalize source relationships, retaining selected or all features."""
    records, parents, implicit = [], {}, {}
    for line in original.splitlines(keepends=True):
        if line.strip() == '##FASTA':
            break
        if line.startswith('#') or not line.strip():
            records.append((line, None, {}, []))
            continue
        fields = line.rstrip().split('\t')
        if len(fields) != 9:
            raise ValueError('Invalid source GFF row')
        is_gtf = bool(re.match(r'^\s*[^\s=;]+\s+', fields[8]))
        parsed = parse_gff_attributes(fields[8])
        attributes = {key: values[0] for key, values in parsed.items() if values}
        parent_ids = list(parsed.get('Parent', ()))
        tid, gid = attributes.get('transcript_id'), attributes.get('gene_id')
        # GTF encodes its graph on every child rather than requiring ID/Parent rows.
        if tid and not attributes.get('ID') and fields[2].lower().endswith(('rna', 'transcript')):
            attributes['ID'] = tid
        elif gid and not tid and fields[2].lower().endswith('gene') and not attributes.get('ID'):
            attributes['ID'] = gid
        if not parent_ids and tid:
            parent_ids = [gid] if attributes.get('ID') == tid else [tid]
            parent_ids = list(filter(None, parent_ids))
        if tid and gid:
            parents.setdefault(tid, set()).add(gid)
            implicit.setdefault(tid, (gid, fields[0], fields[6], []) )[3].append((int(fields[3]), int(fields[4])))
        identifier = attributes.get('ID')
        if identifier:
            parents.setdefault(identifier, set()).update(parent_ids)
        # Produce GFF3 for every effective view; exact source syntax remains archived.
        if parsed and (is_gtf or bool(re.fullmatch(r'[^\s=;]+', fields[8].strip()))):
            values = dict(parsed)
            for key in ('ID',):
                if key in attributes:
                    values[key] = (attributes[key],)
            if parent_ids:
                values['Parent'] = tuple(parent_ids)
            fields[8] = ';'.join(key + '=' + ','.join(quote(value, safe='._-:') for value in vals) for key, vals in values.items())
            line = '\t'.join(fields) + '\n'
        else:
            # Preserve structured GFF3 fields such as Target and Gap. Add only
            # graph attributes missing from transcript/gene identity evidence.
            additions = []
            if attributes.get('ID') and not parsed.get('ID'):
                additions.append('ID=' + quote(attributes['ID'], safe='._-:'))
            if parent_ids and not parsed.get('Parent'):
                additions.append('Parent=' + ','.join(quote(parent, safe='._-:') for parent in parent_ids))
            if additions:
                fields[8] = fields[8].rstrip(';') + ';' + ';'.join(additions)
                line = '\t'.join(fields) + '\n'
        records.append((line, fields, attributes, parent_ids))
    if retain_all:
        transcript_ids = set(transcript_ids) | set(implicit)
    ancestors, pending = set(), list(transcript_ids)
    while pending:
        identifier = pending.pop()
        for parent in parents.get(identifier, set()) - ancestors - {''}:
            ancestors.add(parent)
            pending.append(parent)
    coding_owners = {parent for _, fields, _, pids in records if fields is not None and fields[2] == 'CDS' for parent in pids}
    excluded_coding_owners = coding_owners - transcript_ids
    selected_ancestors = transcript_ids | ancestors
    descendants = set(transcript_ids)
    changed = True
    while changed:
        changed = False
        for _, fields, attributes, parent_ids in records:
            if fields is None or attributes.get('ID') in ancestors or attributes.get('ID') in excluded_coding_owners:
                continue
            if descendants.intersection(parent_ids):
                identifier = attributes.get('ID')
                if identifier and identifier not in descendants:
                    descendants.add(identifier)
                    changed = True
    structural_parents = descendants | ancestors
    lines, known_ids = [], {attrs['ID'] for _, f, attrs, _ in records if f is not None and 'ID' in attrs}
    gene_spans = defaultdict(list)
    for tid in sorted(transcript_ids & implicit.keys()):
        gid, seqid, strand, spans = implicit[tid]
        start, end = min(a for a, _ in spans), max(b for _, b in spans)
        gene_spans[(gid, seqid, strand)].append((start, end))
        if tid not in known_ids:
            kind = 'mRNA' if tid in coding_owners else 'transcript'
            lines.append(f'{seqid}\tGeneGalleon\t{kind}\t{start}\t{end}\t.\t{strand}\t.\tID={quote(tid, safe="._-:")};Parent={quote(gid, safe="._-:")};gene_id={quote(gid, safe="._-:")}\n')
    for _, fields, attributes, parent_ids in records:
        if (fields is not None and fields[2].lower().endswith(('rna', 'transcript'))
                and (retain_all or attributes.get('ID') in selected_ancestors)):
            for gid in parent_ids:
                gene_spans[(gid, fields[0], fields[6])].append((int(fields[3]), int(fields[4])))
    synthesized_genes = []
    for (gid, seqid, strand), spans in sorted(gene_spans.items()):
        if gid not in known_ids:
            start, end = min(a for a, _ in spans), max(b for _, b in spans)
            synthesized_genes.append(f'{seqid}\tGeneGalleon\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={quote(gid, safe="._-:")}\n')
    lines = ['##gff-version 3\n'] + synthesized_genes + lines
    for line, fields, attributes, parent_ids in records:
        if fields is None:
            if not line.startswith('##gff-version'):
                lines.append(line)
            continue
        identifier = attributes.get('ID', '')
        structural = identifier in transcript_ids or identifier in ancestors
        keep = retain_all or structural or identifier in descendants or (identifier not in coding_owners and bool(descendants.intersection(parent_ids)))
        if keep:
            if parent_ids and not retain_all:
                retained = [parent for parent in parent_ids if parent in (structural_parents if structural else descendants)]
                if not retained:
                    continue
                parts = fields[8].split(';')
                fields[8] = ';'.join('Parent=' + ','.join(quote(parent, safe='._-:') for parent in retained) if part.startswith('Parent=') else part for part in parts)
                line = '\t'.join(fields) + '\n'
            lines.append(line)
    return ''.join(lines)


def extend_gene_bounds(original, bounds):
    lines, changes = [], []
    for line in original.splitlines(keepends=True):
        if line.startswith('#') or not line.strip():
            lines.append(line)
            continue
        fields = line.rstrip().split('\t')
        if len(fields) == 9 and fields[2].lower().endswith('gene'):
            parsed = parse_gff_attributes(fields[8])
            identifier = next(iter(parsed.get('ID', parsed.get('gene_id', ()))), None)
            if identifier in bounds:
                start, end = bounds[identifier]
                a, b = min(int(fields[3]), start + 1), max(int(fields[4]), end)
                if (a, b) != (int(fields[3]), int(fields[4])):
                    changes.append({'source_gene_id': identifier, 'before': [int(fields[3]), int(fields[4])], 'after': [a, b]})
                    fields[3], fields[4] = str(a), str(b)
                    line = '\t'.join(fields) + '\n'
        lines.append(line)
    return ''.join(lines), changes


def analysis_coding_candidate(candidate, genetic_code=None):
    """Clip only already-known partial codon bases for the CDS analysis view."""
    from Bio.Seq import Seq
    quality = candidate.get('quality', {})
    expected = str(candidate.get('protein', '')).upper()
    if not quality.get('usable') or not expected or '*' in expected:
        raise ValueError('Analysis CDS requires an admitted source translation')
    offset = int(quality.get('translation_offset', 0))
    if offset not in {0, 1, 2}:
        raise ValueError('Invalid known translation offset')
    strand = candidate.get('strand')
    if strand not in {'+', '-'}:
        raise ValueError('Analysis CDS requires a known strand')
    sequence = str(candidate['cds']).upper()
    blocks = sorted((list(b) for b in candidate['blocks']), key=lambda b: (b[0], b[1]), reverse=strand == '-')
    if not blocks or any(b[0] < 0 or b[1] <= b[0] for b in blocks) or sum(b[1] - b[0] for b in blocks) != len(sequence):
        raise ValueError('Analysis CDS geometry differs from source sequence')
    coding = sequence[offset:]
    trailing = len(coding) % 3
    coding = coding[:len(coding) - trailing]
    code = int(genetic_code if genetic_code is not None else quality.get('genetic_code', 1))
    if not coding or str(Seq(coding).translate(table=code)).removesuffix('*') != expected:
        raise ValueError('Analysis CDS does not reproduce the admitted protein')
    def clip(amount, at_start):
        while amount:
            index = 0 if at_start else -1
            width = blocks[index][1] - blocks[index][0]
            removed = min(width, amount)
            if removed == width:
                blocks.pop(index)
            elif (strand == '+') == at_start:
                blocks[index][0] += removed
            else:
                blocks[index][1] -= removed
            amount -= removed
    clip(offset, True)
    clip(trailing, False)
    cumulative = 0
    junctions = []
    for index, block in enumerate(blocks):
        block[2] = (3 - cumulative % 3) % 3
        cumulative += block[1] - block[0]
        if index + 1 < len(blocks):
            following = blocks[index + 1]
            donor, acceptor = (block[1], following[0]) if strand == '+' else (block[0], following[1])
            junctions.append([donor, acceptor, cumulative % 3])
    result = copy.deepcopy(candidate)
    result.update(cds=coding, protein=expected, blocks=blocks, junctions=junctions,
                  analysis={'removed_first_bases': offset, 'removed_final_bases': trailing,
                            'source_cds_length': len(sequence), 'analysis_cds_length': len(coding)})
    result['quality'].update(translation_offset=0, incomplete_codon=False)
    return result


def analysis_gff_rows(selected_gff_text, analysis_candidates):
    """Keep the selected source graph and clip CDS rows per admitted owner."""
    candidates = list(analysis_candidates.values()) if isinstance(analysis_candidates, dict) else list(analysis_candidates)
    by_owner = {}
    indexes = {}
    for candidate in candidates:
        owner = candidate['source_transcript_id']
        if owner in by_owner:
            raise ValueError('Duplicate admitted transcript identity')
        if not candidate.get('quality', {}).get('usable') or not candidate.get('protein'):
            raise ValueError('Analysis GFF requires admitted coding paths')
        by_owner[owner] = candidate
        ordered = sorted(tuple(b) for b in candidate['blocks'])
        indexes[owner] = [b[0] for b in ordered], ordered
    source = selected_gff_rows(selected_gff_text, set(by_owner))
    records, used_ids, owner_shapes, exported = [], set(), defaultdict(lambda: defaultdict(list)), defaultdict(list)
    for line in source.splitlines(keepends=True):
        if line.startswith('#') or not line.strip():
            records.append((line, None, {}, {}))
            continue
        fields = line.rstrip('\r\n').split('\t')
        attributes = parse_gff_attributes(fields[8])
        identifier = next(iter(attributes.get('ID', ())), None)
        if identifier:
            used_ids.add(identifier)
        shapes = {}
        if fields[2] == 'CDS':
            owners = attributes.get('Parent', ()) or attributes.get('transcript_id', ()) or attributes.get('ID', ())
            start, end = int(fields[3]) - 1, int(fields[4])
            for owner in sorted(set(owners) & by_owner.keys()):
                candidate = by_owner[owner]
                if (fields[0], fields[6]) != (candidate['seqid'], candidate['strand']):
                    raise ValueError('Analysis GFF source axis differs from admitted path')
                starts, blocks = indexes[owner]
                position = bisect.bisect_left(starts, start)
                matching = []
                while position < len(blocks) and blocks[position][0] < end:
                    block = blocks[position]
                    if block[1] > end:
                        raise ValueError('Analysis CDS extends beyond its source block')
                    matching.append(block)
                    position += 1
                shapes[owner] = matching
                exported[owner].extend(matching)
                if identifier:
                    owner_shapes[identifier][owner].extend(matching)
        records.append((line, fields, attributes, shapes))
    for owner, candidate in by_owner.items():
        if sorted(exported[owner]) != sorted(tuple(b) for b in candidate['blocks']):
            raise ValueError('Analysis GFF does not contain the exact admitted CDS path')
    feature_ids = {}
    replacements = defaultdict(set)
    for identifier, paths in sorted(owner_shapes.items()):
        groups = defaultdict(list)
        for owner, shape in sorted(paths.items()):
            groups[tuple(sorted(shape))].append(owner)
        for shape, owners in sorted(groups.items()):
            assigned = identifier
            if len(groups) > 1 and shape:
                key = json.dumps([identifier, shape, owners], sort_keys=True)
                assigned = identifier + '_gganalysis_' + hashlib.sha256(key.encode()).hexdigest()[:16]
                suffix = 0
                while assigned in used_ids:
                    suffix += 1
                    assigned = identifier + '_gganalysis_' + hashlib.sha256(key.encode()).hexdigest()[:16] + '_' + str(suffix)
                used_ids.add(assigned)
            for owner in owners:
                feature_ids[identifier, owner] = assigned
            if shape:
                replacements[identifier].add(assigned)
        replacements.setdefault(identifier, set())
    def replace_attribute(text, key, values):
        replacement = key + '=' + ','.join(quote(str(value), safe='._-:') for value in values)
        parts = text.split(';')
        if any(part.startswith(key + '=') for part in parts):
            return ';'.join(replacement if part.startswith(key + '=') else part for part in parts)
        return text.rstrip(';') + ';' + replacement
    output, noncoding = [], []
    for line, fields, attributes, shapes in records:
        if fields is None:
            output.append(line)
        elif fields[2] != 'CDS':
            noncoding.append((fields, attributes))
            output.append((fields, attributes))
        else:
            identifier = next(iter(attributes.get('ID', ())), None)
            groups = defaultdict(list)
            for owner, blocks in shapes.items():
                assigned = feature_ids.get((identifier, owner), identifier)
                for block in blocks:
                    groups[block, assigned].append(owner)
            for (block, assigned), owners in sorted(groups.items()):
                revised = list(fields)
                revised[3], revised[4], revised[7] = str(block[0] + 1), str(block[1]), str(block[2])
                if attributes.get('Parent'):
                    revised[8] = replace_attribute(revised[8], 'Parent', sorted(owners))
                if identifier and assigned != identifier:
                    revised[8] = replace_attribute(revised[8], 'ID', [assigned])
                output.append('\t'.join(revised) + '\n')
    changed = True
    while changed:
        changed = False
        retained = []
        for fields, attributes in noncoding:
            parents = attributes.get('Parent', ())
            mapped = sorted({replacement for parent in parents for replacement in replacements.get(parent, {parent})})
            if parents and not mapped:
                identifier = next(iter(attributes.get('ID', ())), None)
                if identifier and replacements.get(identifier) != set():
                    replacements[identifier] = set()
                    changed = True
                continue
            retained.append((fields, attributes))
        noncoding = retained
    retained = {id(fields) for fields, _ in noncoding}
    lines = []
    for item in output:
        if isinstance(item, str):
            lines.append(item)
        else:
            fields, attributes = item
            if id(fields) not in retained:
                continue
            parents = attributes.get('Parent', ())
            mapped = sorted({replacement for parent in parents for replacement in replacements.get(parent, {parent})})
            if parents and list(parents) != mapped:
                fields[8] = replace_attribute(fields[8], 'Parent', mapped)
            lines.append('\t'.join(fields) + '\n')
    return ''.join(lines)


@invocation_cached
def finalize(root, value):
    final = select(root, value, predictions=True)
    db = catalog_index(root, value, predictions=True)
    selections = json.loads((final / 'selection.json').read_text())['selections']
    dependencies = {'selection': digest(final / 'receipt.json'),
                    'catalog': {n: digest(root / 'catalog' / n / 'receipt.json') for n in value['species']},
                    'index': {db.parent.name: digest(db.parent / 'receipt.json')}}
    def build(tmp):
        for role in ('species_cds', 'species_protein', 'species_gff', 'species_genome', 'analysis_cds', 'analysis_gff', 'full_annotation', 'source_annotation', 'source_cds', 'all_candidates'):
            (tmp / role).mkdir()
        shutil.copyfile(final / 'representative_map.tsv', tmp / 'representative_map.tsv')
        rows, changes, translation_audit, coding_audit, exclusions, gene_bounds_audit = [], [], [], [], [], []
        for name in value['species']:
            catalog = json.loads((root / 'catalog' / name / 'catalog_metadata.json').read_text())
            catalog['loci'] = iter_loci(db, name)
            source = value['request']['sources'][name]
            selected = {r['gene_id']: r for r in selections if r['species'] == name}
            paths = {'cds': tmp / 'species_cds' / (name + '.fa'), 'protein': tmp / 'species_protein' / (name + '.fa'),
                     'gff': tmp / 'species_gff' / (name + '.gff3'), 'genome': tmp / 'species_genome' / (name + '.fa'),
                     'analysis_cds': tmp / 'analysis_cds' / (name + '.fa'), 'analysis_gff': tmp / 'analysis_gff' / (name + '.gff3')}
            # Preserve the exact supplied CDS bytes independently of genomic reconstruction.
            cds_suffix = '.fa.gz' if Path(source['fasta']).name.lower().endswith('.gz') else '.fa'
            shutil.copyfile(source['fasta'], tmp / 'source_cds' / (name + cds_suffix))
            # Copy to make a self-contained verifiable effective input view.
            with open_text(Path(source['genome'])) as handle, paths['genome'].open('w') as out:
                shutil.copyfileobj(handle, out)
            full_path = tmp / 'full_annotation' / (name + '.gff3')
            source_annotation = tmp / 'source_annotation' / (name + '.gff3')
            if Path(source['gff']).name.lower().endswith('.gz'):
                with gzip.open(source['gff'], 'rb') as handle, source_annotation.open('wb') as out:
                    shutil.copyfileobj(handle, out)
            else:
                shutil.copyfile(source['gff'], source_annotation)
            original = source_annotation.read_bytes().decode('utf-8')
            directive = re.search(r'^##FASTA(?:\r?\n|$)', original, flags=re.MULTILINE)
            before_fasta, sep, embedded = (original[:directive.start()], '##FASTA', original[directive.start() + len('##FASTA'):]) if directive else (original, '', '')
            additions, predicted_rows, original_transcripts, expanded_bounds, analysis_candidates = [], [], set(), {}, []
            with paths['cds'].open('w') as cds, paths['protein'].open('w') as pep, paths['gff'].open('w') as gff, \
                    paths['analysis_cds'].open('w') as analysis_cds, \
                    (tmp / 'all_candidates' / (name + '.fa')).open('w') as all_cds:
                for gene in catalog['loci']:
                    chosen = selected[gene['gene_id']]
                    candidates = {c['candidate_id']: c for c in gene['candidates']}
                    candidate = candidates[chosen['candidate_id']]
                    exportable = bool(candidate['cds']) and bool(candidate['blocks']) and not gene.get('ambiguous_coordinates') and not candidate['quality'].get('sequence_mismatch') and not candidate['quality'].get('structure_problem')
                    if exportable:
                        cds.write(f">{gene['gene_id']}\n{candidate['cds']}\n")
                    else:
                        problem = candidate['quality'].get('structure_problem')
                        reason = 'invalid_annotation_structure:' + str(problem) if problem else 'source_cds_sequence_mismatch' if candidate['quality'].get('sequence_mismatch') else 'unreconstructible_or_ambiguous_coordinates'
                        exclusions.append((name, gene['gene_id'], candidate['candidate_id'], reason))
                    admitted = exportable and candidate['quality'].get('usable', False) and bool(candidate['protein']) and '*' not in candidate['protein']
                    if admitted:
                        pep.write(f">{gene['gene_id']}\n{candidate['protein']}\n")
                        coding_candidate = analysis_coding_candidate(candidate)
                        analysis_cds.write(f">{gene['gene_id']}\n{coding_candidate['cds']}\n")
                        analysis_candidates.append(coding_candidate)
                        head = candidate['quality'].get('translation_offset', 0)
                        tail = len(candidate['cds']) - head - len(coding_candidate['cds'])
                        coding_audit.append((name, gene['gene_id'], candidate['candidate_id'], 'included', head, tail,
                                             len(candidate['cds']), len(coding_candidate['cds'])))
                    else:
                        coding_audit.append((name, gene['gene_id'], candidate['candidate_id'], 'excluded', '', '', len(candidate['cds']), 0))
                    translation_audit.append((name, gene['gene_id'], candidate['candidate_id'], 'included' if admitted else 'excluded', json.dumps(candidate['quality'], sort_keys=True)))
                    if not exportable:
                        pass
                    elif candidate['origin'] == 'original':
                        original_transcripts.add(candidate['source_transcript_id'])
                    else:
                        predicted_rows.append(candidate_gff(gene, candidate, name))
                    for c in candidates.values():
                        all_cds.write(f">{c['candidate_id']}\n{c['cds']}\n")
                        if c['origin'] == 'predicted':
                            # Preserve the original gene and attach only new transcript/CDS rows.
                            lines = candidate_gff(gene, c, name, gene_id=c['source_gene_id'],
                                                  gene_token=c.get('gene_token', gene.get('gene_token', gene['gene_id']))).splitlines(keepends=True)[1:]
                            start, end = min(b[0] for b in c['blocks']), max(b[1] for b in c['blocks'])
                            old = expanded_bounds.get(c['source_gene_id'], (start, end))
                            expanded_bounds[c['source_gene_id']] = min(old[0], start), max(old[1], end)
                            additions.extend(lines)
                    changes.append({**chosen, 'selected_origin': candidate['origin']})
                gff.write(selected_gff_rows(original, original_transcripts) + ''.join(predicted_rows))
            paths['analysis_gff'].write_text(analysis_gff_rows(paths['gff'].read_text(), analysis_candidates))
            derived_original, bound_changes = extend_gene_bounds(selected_gff_rows(before_fasta, set(), retain_all=True), expanded_bounds)
            gene_bounds_audit.extend(dict(species=name, **row) for row in bound_changes)
            full_path.write_text(derived_original.rstrip() + '\n' + ''.join(additions) + (sep + embedded if sep else ''))
            row = {'species': name, **{k: str((root / 'effective' / p.relative_to(tmp)).resolve()) for k, p in paths.items()},
                   'representative_map': str((root / 'effective' / 'representative_map.tsv').resolve()),
                   'genetic_code': source['genetic_code'], 'representative_map_sha256': digest(tmp / 'representative_map.tsv')}
            for k, path in paths.items():
                row[k + '_sha256'] = digest(path)
            rows.append(row)
        rescue.write_tsv(tmp / 'effective_exclusions.tsv', ('species', 'gene_id', 'candidate_id', 'reason'), exclusions)
        atomic_json(tmp / 'gene_bounds_changes.json', gene_bounds_audit)
        rescue.write_tsv(tmp / 'translation_admission.tsv', ('species', 'gene_id', 'candidate_id', 'status', 'quality'), translation_audit)
        rescue.write_tsv(tmp / 'coding_admission.tsv', ('species', 'gene_id', 'candidate_id', 'status', 'head_bases_removed', 'tail_bases_removed', 'source_cds_length', 'analysis_cds_length'), coding_audit)
        rescue.write_tsv(tmp / 'species_genetic_code.tsv', ('species', 'genetic_code'), [(r['species'], r['genetic_code']) for r in rows])
        fields = list(rows[0])
        rescue.write_tsv(tmp / 'inputs.tsv', fields, [[r[k] for k in fields] for r in rows])
        atomic_json(tmp / 'changes.json', changes)
        atomic_json(tmp / 'summary.json', {'species': len(rows), 'loci': len(selections),
                                         'changed_representatives': sum(r['status'] == 'conserved' for r in selections),
                                         'predicted_selected': sum(r['selected_origin'] == 'predicted' for r in changes)})
    return stage(root, 'effective', {'dependencies': dependencies}, build)


def verify_inputs(inputs, field=None):
    def parsed_identity(path):
        info = Path(path).stat()
        return (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns,
                str(Path(path).resolve(strict=True)))

    input_identity = parsed_identity(inputs)
    rows = read_table(inputs)
    if not rows:
        raise ValueError('Empty effective input manifest')
    roots = set()
    for row in rows:
        roots.add(str(Path(row['cds']).parent.parent))
    if len(roots) != 1:
        raise ValueError('Effective inputs must share one published view')
    root = Path(next(iter(roots)))
    receipt_path = root / 'receipt.json'
    receipt_identity = parsed_identity(receipt_path)
    receipt = json.loads(receipt_path.read_text())
    if (not isinstance(receipt, dict) or 'key' not in receipt or not isinstance(receipt.get('files'), dict)
            or not receipt['files'] or any(not isinstance(path, str) or not isinstance(value, str)
                                          or not (root / path).is_file() for path, value in receipt['files'].items())):
        raise ValueError('Effective view is incomplete or corrupt')
    # Fresh full hashes are shared only inside this boundary. Include every
    # receipt member, including unselected annotations and admitted views, as
    # well as external role paths supported by existing manifests. digest_paths
    # deduplicates resolved aliases and fences mutations across the whole batch.
    roles = ('cds', 'protein', 'gff', 'genome', 'representative_map')
    paths = [inputs, receipt_path, *(row[role] for row in rows for role in roles),
             *(root / path for path in receipt['files'])]
    hashes = digest_paths(paths)
    # Parsing happened before the batch began; bind the parsed rows/receipt to
    # those same bytes even if they changed before their fresh hash was taken.
    if parsed_identity(inputs) != input_identity or parsed_identity(receipt_path) != receipt_identity:
        raise OSError('Effective manifest or receipt changed while verifying')
    for row in rows:
        for role in roles:
            if hashes[str(row[role])] != row[role + '_sha256']:
                message = ('Effective representative map changed' if role == 'representative_map'
                           else 'Effective input changed: ' + row[role])
                raise ValueError(message)
    if any(hashes[str(root / path)] != value for path, value in receipt['files'].items()):
        raise ValueError('Effective view is incomplete or corrupt')
    # A copied TSV must still be the published manifest. Role hashes alone do
    # not bind its species labels, genetic codes, row membership, or paths.
    if hashes[str(inputs)] != receipt['files'].get('inputs.tsv'):
        raise ValueError('Effective input manifest differs from the published view')
    if field in {'layout', 'analysis_layout', 'coding_layout'}:
        coding = field != 'layout'
        names = ('analysis_cds' if coding else 'species_cds', 'species_protein', 'analysis_gff' if coding else 'species_gff',
                 'species_genome', 'representative_map.tsv', 'species_genetic_code.tsv')
        if coding and any(not (root / name).exists() for name in names[:3]):
            raise ValueError('Published view lacks admitted coding inputs; generate a new refinement output')
        layout = '\t'.join(str(root / name) for name in names)
        if field == 'coding_layout':
            codes = {int(row['genetic_code']) for row in rows if Path(row['analysis_cds']).stat().st_size}
            if len(codes) != 1:
                raise ValueError('CDS analysis requires one common genetic code; use input_sequence_mode=protein for mixed codes')
            layout += '\t' + str(codes.pop())
        return layout
    if field:
        return str(root / ('representative_map.tsv' if field == 'representative_map' else field if field.startswith('analysis_') else 'species_' + field))
    return rows


def inspect_stages(root, value):
    """Check the required dependency graph, including absent publications.

    A receipt's own key is not evidence that it belongs to this run or its
    current upstream stages. Reconstruct those keys from the frozen plan and
    the current dependency receipts without building or repairing anything.
    """
    plan_hash = digest(root / 'plan.json')

    def receipt_hash(relative):
        try:
            return digest(root / relative / 'receipt.json')
        except OSError:
            return None

    catalogs = {name: receipt_hash(Path('catalog') / name) for name in value['species']}
    predictions = {name: receipt_hash(Path('predictions') / name) for name in value['species']}
    initial_index = {'catalog_index': receipt_hash(Path('catalog_index'))}
    final_index = {'catalog_index_final': receipt_hash(Path('catalog_index_final'))}
    expected = {str(Path('catalog') / name): {} for name in value['species']}
    expected['catalog_index'] = {'dependencies': {'catalog': catalogs}}
    correspondence_dependencies = {'catalog': catalogs, 'index': initial_index}
    anchors = {}
    request = value['request']
    if request['rescue_output'] and not request['edges']:
        anchor_root = Path(request['rescue_output'])
        anchor_plan = rescue.load(anchor_root)
        anchor_hashes = {}
        for job in anchor_plan['synteny_jobs']:
            if job['a'] == job['b']:
                continue
            directory = anchor_root / 'synteny' / job['id']
            try:
                anchor_hashes[job['id']] = digest(directory / 'receipt.json')
                anchors['anchors/' + job['id']] = rescue.verified(directory, rescue.comparison_key(anchor_root, job))
            except (OSError, ValueError, KeyError):
                anchor_hashes[job['id']] = None
                anchors['anchors/' + job['id']] = False
        correspondence_dependencies['anchors'] = anchor_hashes
        prepared_hashes = {}
        anchor_hash = digest(anchor_root / 'plan.json')
        for name in value['species']:
            directory = anchor_root / 'prepared' / name
            try:
                prepared_hashes[name] = digest(directory / 'receipt.json')
                anchors['prepared/' + name] = rescue.verified(directory, {'plan': anchor_hash, 'species': name})
            except (OSError, ValueError):
                prepared_hashes[name] = None
                anchors['prepared/' + name] = False
        correspondence_dependencies['prepared'] = prepared_hashes
    expected['correspondence'] = {'dependencies': correspondence_dependencies}
    selection_dependencies = {'correspondence': receipt_hash(Path('correspondence')), 'catalog': catalogs,
                              'index': initial_index}
    expected['selection_initial'] = {'dependencies': selection_dependencies}
    prediction_dependencies = {'initial': receipt_hash(Path('selection_initial')), **selection_dependencies}
    expected.update({str(Path('predictions') / name): {'dependencies': prediction_dependencies} for name in value['species']})
    expected['catalog_index_final'] = {'dependencies': {'catalog': catalogs, 'predictions': predictions}}
    expected['selection_final'] = {'dependencies': {**selection_dependencies, 'predictions': predictions, 'index': final_index}}
    expected['effective'] = {'dependencies': {'selection': receipt_hash(Path('selection_final')), 'catalog': catalogs,
                                            'index': final_index}}
    return {relative: rescue.verified(root / relative, {'plan': plan_hash, **key})
            for relative, key in expected.items()} | anchors


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    prep = sub.add_parser('plan')
    prep.add_argument('--output', required=True, type=Path)
    prep.add_argument('--inputs', type=Path)
    prep.add_argument('--rescue-output', type=Path)
    prep.add_argument('--edges', type=Path)
    prep.add_argument('--rna', type=Path, help='Whole coding RNA paths: species,seqid,strand,cds_blocks(JSON),transcript_id,count TSV')
    prep.add_argument('--species-profiles', type=Path, help='Explicit target-species prediction parameter overrides (TSV)')
    for key, default in DEFAULTS.items():
        prep.add_argument('--' + key.replace('_', '-'), type=type(default), default=default)
    for command in ('catalog', 'correspondence', 'select', 'predict', 'finalize', 'run', 'status', 'qc'):
        s = sub.add_parser(command)
        s.add_argument('--output', type=Path, required=True)
        s.add_argument('--task-index', type=int)
        s.add_argument('--cpus', type=int, default=1)
        s.add_argument('--comparison-cache', type=Path)
    check = sub.add_parser('verify-inputs')
    check.add_argument('--inputs', required=True, type=Path)
    check.add_argument('--field', choices=['cds', 'protein', 'gff', 'genome', 'analysis_cds', 'analysis_gff', 'representative_map', 'layout', 'analysis_layout', 'coding_layout'])
    return p


def main():
    global _INVOCATION_CACHE
    _INVOCATION_CACHE = {}
    args = parser().parse_args()
    if args.command == 'verify-inputs':
        result = verify_inputs(args.inputs, args.field)
        print(result if isinstance(result, str) else json.dumps(result, sort_keys=True))
        return
    if args.command == 'plan':
        plan(args.output, args.inputs, args.rescue_output, args.edges, args.rna,
             species_profiles=args.species_profiles, **{k: getattr(args, k) for k in DEFAULTS})
        return
    if args.cpus < 1:
        raise ValueError('--cpus must be positive')
    root = args.output.resolve()
    value = load(root)
    names = value['species']
    if args.task_index is not None:
        if not 1 <= args.task_index <= len(names):
            raise ValueError('Invalid task index')
        names = [names[args.task_index - 1]]
    if args.command in {'catalog', 'run'}:
        for name in names:
            catalog_species(root, value, name)
    if args.command in {'correspondence', 'run'}:
        correspondence(root, value, args.cpus, args.comparison_cache)
    if args.command in {'select', 'run'}:
        select(root, value)
    if args.command in {'predict', 'run'}:
        for name in names:
            predict_species(root, value, name, args.cpus)
    if args.command in {'finalize', 'run'}:
        finalize(root, value)
    if args.command in {'status', 'qc'}:
        if args.command == 'qc':
            verify_inputs(root / 'effective' / 'inputs.tsv')
        result = {'species': len(value['species']), 'stages': inspect_stages(root, value)}
        if args.command == 'qc' and not all(result['stages'].values()):
            invalid = [name for name, valid in result['stages'].items() if not valid]
            raise ValueError('Corrupt refinement stage or missing publication: ' + ', '.join(invalid))
        print(json.dumps(result, sort_keys=True))
    load(root)


if __name__ == '__main__':
    main()

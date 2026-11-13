#!/usr/bin/env python3
"""Restartable coding-model refinement and conserved isoform selection.

Original annotations are immutable. Genomic representatives have exact GFF
paths, and assembly representatives retain their original transcript IDs.
"""
import argparse
import bisect
import contextlib
import copy
import csv
import functools
import gzip
import hashlib
import importlib.metadata
import io
import json
import math
import os
import re
import resource
import shutil
import sqlite3
import sys
import tempfile
import time
from collections import defaultdict, deque
from concurrent.futures import ProcessPoolExecutor
from itertools import chain, islice
from multiprocessing import get_context
from pathlib import Path
from urllib.parse import quote

from Bio.Seq import Seq
from cdskit.translate import translate_sequence_string

try:
    import gene_model_transcriptome as transcriptome
    import refinement_rna
    import refinement_transcriptome_edges
    import rescue_gene_models as rescue
    from fasta_sequence_store import exclusive_lock, open_text
    from filesystem_io import scratch_build_directory
    from format_species_annotation.common import parse_gff_attributes
    from gene_model_catalog import (
        SOURCE_RESCUE_LOCUS_ATTRIBUTES,
        _formatted_id,
        build_catalog,
        indexed_genome,
        validate_candidate,
        write_catalog,
    )
    from gene_model_selection import _baseline, pair_score, select_representatives
    from gene_model_species_profiles import parameters_for, read_profiles
    from gene_model_store import _connection as store_connection
    from gene_model_store import build_store, iter_loci, iter_locus_keys, load_locus, select_from_store
    from input_generation_array_state import FreshDigestBatch, atomic_json, digest, digest_paths
    from refinement_input_snapshot import RefinementInputSnapshot
    from refinement_receipt_snapshot import ReceiptSnapshot
    from rescue_coding_paths import coding_overlap
    from rescue_homology_policy import homology_supported
    from rescue_model_store import frozen_model_store_key, iter_revision_models, verify_model_store_key
    from rescue_prediction_cache import stream_json_array
    from rescue_prepared_snapshot import PreparedSnapshot
except ImportError:
    from . import gene_model_transcriptome as transcriptome
    from . import refinement_rna, refinement_transcriptome_edges
    from . import rescue_gene_models as rescue
    from .fasta_sequence_store import exclusive_lock, open_text
    from .filesystem_io import scratch_build_directory
    from .format_species_annotation.common import parse_gff_attributes
    from .gene_model_catalog import (
        SOURCE_RESCUE_LOCUS_ATTRIBUTES,
        _formatted_id,
        build_catalog,
        indexed_genome,
        validate_candidate,
        write_catalog,
    )
    from .gene_model_selection import _baseline, pair_score, select_representatives
    from .gene_model_species_profiles import parameters_for, read_profiles
    from .gene_model_store import _connection as store_connection
    from .gene_model_store import build_store, iter_loci, iter_locus_keys, load_locus, select_from_store
    from .input_generation_array_state import FreshDigestBatch, atomic_json, digest, digest_paths
    from .refinement_input_snapshot import RefinementInputSnapshot
    from .refinement_receipt_snapshot import ReceiptSnapshot
    from .rescue_coding_paths import coding_overlap
    from .rescue_homology_policy import homology_supported
    from .rescue_model_store import frozen_model_store_key, iter_revision_models, verify_model_store_key
    from .rescue_prediction_cache import stream_json_array
    from .rescue_prepared_snapshot import PreparedSnapshot

SCHEMA = 1
MAP_FIELDS = ('species', 'gene_id', 'candidate_id', 'source_transcript_id', 'status', 'score', 'margin', 'reason')
CONCRETE_CDS = re.compile('[ACGT]+')
OWNERSHIP_ATTRIBUTES = frozenset(('ID', 'Parent', 'gene_id', 'transcript_id'))
DEFAULTS = dict(policy='conserved', mode='conservative', isoform_adoption='rna_required', min_margin=0.10, min_support=2,
                minimum_coverage=0.70, minimum_identity=0.50, minimum_broad_coverage=0.90,
                minimum_lower_coverage_identity=0.70, max_intron=20000,
                padding=2000, max_interval=200000, candidate_limit=32)


# Invocation-local reuse avoids re-reading every catalog for every target worker.
# Fresh CLI boundaries and each publication still verify frozen source bytes.
_INVOCATION_CACHE = None
_RECEIPT_SNAPSHOT = None
_INPUT_SNAPSHOT = None


def _receipt_verifier(cpus):
    """Bound parallel full-byte member hashing by the allocated CPUs."""
    return functools.partial(rescue.verified, hash_workers=min(cpus, 4))


@contextlib.contextmanager
def invocation_context(cpus=1):
    """Retain immutable dependency proofs only within an explicit command."""
    global _INVOCATION_CACHE, _RECEIPT_SNAPSHOT, _INPUT_SNAPSHOT
    previous = _INVOCATION_CACHE, _RECEIPT_SNAPSHOT, _INPUT_SNAPSHOT
    snapshot = ReceiptSnapshot(_receipt_verifier(cpus))
    inputs = RefinementInputSnapshot()
    _INVOCATION_CACHE, _RECEIPT_SNAPSHOT = {}, snapshot
    _INPUT_SNAPSHOT = inputs
    try:
        with inputs, snapshot:
            yield
    finally:
        _INVOCATION_CACHE, _RECEIPT_SNAPSHOT, _INPUT_SNAPSHOT = previous


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
                    'format_species_provider_config.py', 'format_species_taxonomy.py',
                    'rescue_terminal_completion.py', 'rescue_coding_paths.py',
                    'rescue_homology_policy.py', 'rescue_path_compatibility.py', 'rescue_boundary_revalidation.py',
                    'rescue_model_store.py', 'rescue_prediction_cache.py', 'rescue_raw_validation.py',
                    'rescue_prepared_snapshot.py', 'refinement_receipt_snapshot.py',
                    'refinement_input_snapshot.py', 'annevo_gene_models.py', 'oriongeno_runtime.py', 'oriongeno_model_profiles.py', 'transcript_utrs.py', 'genome_sequence_identity.py',
                    'gene_model_source_retention.py', 'gene_model_transcriptome.py', 'rename_rnabloom_transcripts.py',
                    'model_protein_correspondence.py', 'refinement_rna.py', 'refinement_transcriptome_edges.py')
    files = [support / name for name in dependencies] + list((support / 'format_species_annotation').rglob('*.py'))
    return {str(Path(m.__file__).name): digest(m.__file__) for m in modules} | {Path(__file__).name: digest(__file__)} | {str(p.relative_to(support)): digest(p) for p in files}


def dependency_identities():
    from Bio.Align import _pairwisealigner
    return {'versions': {name: importlib.metadata.version(name) for name in ('biopython', 'pysam')},
            'python': sys.version, 'sqlite': sqlite3.sqlite_version,
            'pairwise_aligner_sha256': digest(_pairwisealigner.__file__)}


def read_table(path):
    return rescue.table(path)


def revision_worker(request, name):
    """Resolve either frozen revision publication without inventing a JSON path."""
    stored = request.get('revision_model_stores', {}).get(name)
    if stored:
        return Path(stored['directory'])
    legacy = request.get('revision_candidates', {}).get(name)
    return Path(legacy).parent if legacy else None


def revision_models(request, name):
    directory = revision_worker(request, name)
    if directory is None:
        return iter(())
    stored = request.get('revision_model_stores', {}).get(name)
    if stored:
        return iter_revision_models(directory, frozen_key=stored['key'])
    return iter_revision_models(directory)


def read_model_publication(path, receipt_sha, *, required=True):
    receipt = (path.parent / 'receipt.json').read_bytes()
    if hashlib.sha256(receipt).hexdigest() != receipt_sha:
        raise ValueError('Rescue worker receipt changed during planning: ' + str(path.parent))
    expected = json.loads(receipt)['files'].get(path.name)
    if expected is None:
        if required:
            raise ValueError('Rescue has no verified candidate publication: ' + str(path))
        return [], None
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError('Candidates differ from the verified rescue publication: ' + str(path))
    return json.loads(data), expected


def freeze_model_reviews(approvals, anchor_root, receipts, publication, label):
    approvals = Path(approvals).resolve()
    data = approvals.read_bytes()
    rows = list(csv.DictReader(io.StringIO(data.decode(), newline=''), delimiter='\t'))
    if not rows or any(set(row) != {'species', 'candidate_id'} or not all(row.values()) for row in rows):
        raise ValueError(label + ' approvals require species and candidate_id columns and at least one row')
    requested = defaultdict(set)
    for row in rows:
        name, identifier = row['species'], row['candidate_id']
        if name not in receipts:
            raise ValueError(label + ' approval requires completed rescue for species: ' + name)
        if identifier in requested[name]:
            raise ValueError('Duplicate ' + label.lower() + ' approval: ' + identifier)
        requested[name].add(identifier)
    publications = {str(approvals): hashlib.sha256(data).hexdigest()}
    selected = {}
    for name, identifiers in sorted(requested.items()):
        path = anchor_root / 'rescued' / name / publication
        rows, expected = read_model_publication(path, receipts[name])
        candidates = {row['candidate_id']: row for row in rows}
        claimed = set()
        selected[name] = []
        for identifier in sorted(identifiers):
            if identifier not in candidates:
                raise ValueError('Unknown ' + label.lower() + ' candidate for ' + name + ': ' + identifier)
            candidate = candidates[identifier]
            owners = {row['gene_id'] for row in candidate['owner_overlap']}
            if owners & claimed:
                raise ValueError('Approved ' + label.lower() + ' candidates share a source gene: ' + identifier)
            claimed.update(owners)
            selected[name].append(candidate)
        publications[str(path)] = expected
    return selected, publications


def plan(output, inputs=None, rescue_output=None, edges=None, rna=None, species_profiles=None, annevo_inputs=None, oriongeno_inputs=None, *, split_model_merges=None, model_repairs=None, transcriptome_inputs=None, rna_bam_dir=None, transcriptome_assembly_dir=None,
         assembly_genetic_codes=None, assembly_cpus=1, donor_plan=None, target=None, **parameters):
    root = Path(output).resolve()
    root.mkdir(parents=True, exist_ok=True)
    with exclusive_lock(root / '.plan.lock'):
        transcriptome_sources = {}
        transcriptome_proof = FreshDigestBatch()
        if transcriptome_inputs:
            transcriptome_inputs = Path(transcriptome_inputs).resolve()
            transcriptome_proof.read([transcriptome_inputs])
            transcriptome_sources = transcriptome.read_sources(transcriptome_inputs)
        receipts = {}
        pool = None
        genomic_targets = []
        if rescue_output:
            if inputs:
                raise ValueError('Choose --inputs or --rescue-output; explicit inputs must not be silently ignored')
            anchor_root = Path(rescue_output).resolve()
            # Consume the frozen producer; its schema, receipts and inputs are
            # verified below. A downstream catalogue update does not rerun it.
            anchors = rescue.load(anchor_root, inspect_only=True)
            sources = copy.deepcopy(anchors['request']['sources'])
            files = [str(anchor_root / 'plan.json')]
            augmented = anchor_root / 'augmented'
            pool = anchors['request'].get('fixed_donor_pool')
            if donor_plan and (not pool or Path(donor_plan).resolve() != Path(pool['root'])):
                raise ValueError('Refinement donor plan must match the rescue donor plan')
            if pool:
                genomic_targets = anchors['species']
                if target and target not in genomic_targets:
                    raise ValueError('Requested refinement target does not match the rescue target')
                if len(genomic_targets) != 1 or not augmented.is_dir():
                    raise ValueError('Fixed-panel refinement requires completed rescue for exactly one genomic target')
                target, = genomic_targets
                used = set(anchors['donors'][target])
                if used != set(pool['species']) - {target}:
                    raise ValueError('Rescue donors differ from the fixed panel')
                sources = {name: sources[name] for name in [target, *sorted(used)]}
            revision_candidates, revision_model_stores, revision_donor_maps = {}, {}, {}
            if augmented.exists():
                for name in anchors['species']:
                    directory = anchor_root / 'rescued' / name
                    if not rescue.verified(directory, rescue.rescue_key(anchor_root, anchors, name)):
                        raise ValueError('Missing-gene rescue incomplete or corrupt: ' + name)
                    receipts[name] = digest(directory / 'receipt.json')
                    files.append(str(directory / 'receipt.json'))
                    revision_path = directory / 'revision_candidates.json'
                    if (directory / 'model_store/manifest.json').is_file():
                        stored_key = frozen_model_store_key(directory, kind='revision')
                        verify_model_store_key(directory, stored_key)
                        revision_model_stores[name] = {'directory': str(directory.resolve()), 'key': stored_key}
                        files.extend(str((directory / member).resolve()) for member in stored_key['files'])
                    elif revision_path.is_file():
                        worker_receipt = json.loads((directory / 'receipt.json').read_text())
                        if worker_receipt['files'].get('revision_candidates.json') != digest(revision_path):
                            raise ValueError('Rescue revisions are not a verified worker publication: ' + name)
                        revision_candidates[name] = str(revision_path.resolve())
                        files.append(str(revision_path.resolve()))
                augmented_key = {'plan': digest(anchor_root / 'plan.json'), 'rescue_receipts': receipts}
                if not rescue.verified(augmented, augmented_key):
                    raise ValueError('Augmented rescue inputs incomplete or corrupt')
                rows = read_table(augmented / 'inputs.tsv')
                if len(rows) != len(anchors['species']) or {row['species'] for row in rows} != set(anchors['species']):
                    raise ValueError('Augmented rescue input species disagree with the anchor plan')
                for row in rows:
                    source = sources[row['species']]
                    for field, role in (('cds', 'fasta'), ('gff', 'gff'), ('genome', 'genome')):
                        path = Path(row[field])
                        source[role] = str((path if path.is_absolute() else augmented / path).resolve())
                files += [str(augmented / 'receipt.json'), str(augmented / 'inputs.tsv')]
                donors = set()
                for name in sorted(set(revision_candidates) | set(revision_model_stores)):
                    request_revisions = {'revision_candidates': revision_candidates,
                                         'revision_model_stores': revision_model_stores}
                    for model in revision_models(request_revisions, name):
                        donors.update(e['donor'] for e in [model.get('evidence', {}), *model.get('support', [])]
                                      if e.get('donor'))
                for donor in sorted(donors):
                    if donor not in sources:
                        raise ValueError('Rescue revision donor is absent from the frozen cohort: ' + donor)
                    directory = anchor_root / 'prepared' / donor
                    if not rescue.verified(directory, {'plan': digest(anchor_root / 'plan.json'), 'species': donor}):
                        raise ValueError('Rescue revision donor preparation is incomplete: ' + donor)
                    metadata = {role: str((directory / filename).resolve()) for role, filename in
                                (('mapping', 'genes.id_map.tsv'), ('protein', 'genes.pep'), ('receipt', 'receipt.json'))}
                    revision_donor_maps[donor] = metadata
                    files.extend(metadata.values())
        else:
            anchor_root = None
            revision_candidates, revision_model_stores, revision_donor_maps = {}, {}, {}
            sources = {}
            for row in read_table(inputs) if inputs else ():
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
            files = [str(Path(inputs).resolve())] if inputs else []
            if donor_plan:
                if inputs:
                    raise ValueError('Fixed-panel genomic refinement requires a completed fixed-panel rescue')
                donor_root = Path(donor_plan).resolve(strict=True)
                donors = rescue.load(donor_root, inspect_only=True)
                if donors['request'].get('fixed_donor_pool'):
                    raise ValueError('A target rescue plan cannot define a donor panel')
                pool = {'root': str(donor_root), 'species': donors['species'],
                        'plan_sha256': rescue.plan_digest(donor_root, donors)}
                sources = copy.deepcopy(donors['request']['sources'])
        if transcriptome_assembly_dir and not transcriptome_inputs and not genomic_targets:
            if pool and not target:
                raise ValueError('Use --target to select one assembly from a fixed-panel assembly directory')
            transcriptome_sources, prepared_files = transcriptome.prepare_assemblies(
                transcriptome_assembly_dir, {} if pool else sources, root / 'assembly_inputs', assembly_cpus, assembly_genetic_codes,
                transcriptome_proof, **({'only': target} if target else {}))
            files.extend(prepared_files)
        if pool:
            if not genomic_targets:
                for name in transcriptome_sources:
                    sources.pop(name, None)
            fixed_donors = sorted(sources.keys() - set(genomic_targets))
            targets = sorted(set(genomic_targets) | transcriptome_sources.keys())
            if len(targets) != 1 or (genomic_targets and transcriptome_sources):
                raise ValueError('Use a separate fixed-panel refinement output for each target')
            if target and targets != [target]:
                raise ValueError('Requested refinement target is absent from the inputs')
            if edges:
                raise ValueError('Fixed-panel refinement uses the rescue or protein correspondence')
            donor_root = Path(pool['root'])
            if root.is_relative_to(donor_root) or donor_root.is_relative_to(root):
                raise ValueError('Use separate refinement and donor plan directories')
            donor_plan_path = str(donor_root / 'plan.json')
            donor_record = rescue.load(donor_root, inspect_only=True)
            if rescue.plan_digest(donor_root, donor_record) != pool['plan_sha256']:
                raise ValueError('Frozen donor plan changed')
            expected_donors = {sources[name][role]: donor_record['request']['files'][sources[name][role]]
                               for name in fixed_donors for role in ('fasta', 'gff', 'genome')}
            transcriptome_proof.read(expected_donors, expected=expected_donors)
            transcriptome_proof.read([donor_plan_path], expected={donor_plan_path: pool['plan_sha256']})
            files.append(donor_plan_path)
        elif target:
            raise ValueError('--target requires a fixed donor plan')
        if sources.keys() & transcriptome_sources.keys():
            raise ValueError('A species may occur in only one of the genomic and transcriptome inputs')
        if not sources and not transcriptome_sources:
            raise ValueError('Provide species through --inputs, --rescue-output or --transcriptome-inputs')
        if transcriptome_sources:
            if transcriptome_inputs:
                files.append(str(transcriptome_inputs))
            files.extend(s[k] for s in transcriptome_sources.values() for k in ('fasta', 'transcript_gene_map'))
        inherited = {k: v for k, v in anchors['request']['parameters'].items() if k.startswith('terminal_')} if anchor_root else {}
        params = {**DEFAULTS, **inherited, **parameters}
        if params['policy'] not in {'longest', 'conserved'} or params['mode'] not in {'off', 'audit', 'conservative'}:
            raise ValueError('Unknown selection policy/refinement mode')
        if params['isoform_adoption'] not in {'rna_required', 'conservation_supported'}:
            raise ValueError('Unknown isoform adoption policy')
        if params['isoform_adoption'] == 'conservation_supported' and params['policy'] != 'conserved':
            raise ValueError('Conservation-supported isoform adoption requires the conserved selection policy')
        if any(not math.isfinite(float(params[k])) or not 0 <= float(params[k]) <= 1
               for k in ('min_margin', 'minimum_coverage', 'minimum_identity', 'minimum_broad_coverage',
                         'minimum_lower_coverage_identity')):
            raise ValueError('Invalid probability/margin')
        if any(int(params[k]) < 1 for k in ('min_support', 'max_intron', 'max_interval', 'candidate_limit')) or params['padding'] < 0:
            raise ValueError('Invalid resource/evidence bound')
        files += [s[k] for s in sources.values() for k in ('fasta', 'gff', 'genome')]
        prediction_sources = {n: sources[n] for n in genomic_targets} if pool else sources
        profiles = (read_profiles(species_profiles, prediction_sources) if species_profiles else
                    copy.deepcopy(anchors['request'].get('species_profiles', {})) if anchor_root else {})
        for optional in (edges, rna, species_profiles):
            if optional:
                files.append(str(Path(optional).resolve()))
        rna_bams = refinement_rna.discover(rna_bam_dir, prediction_sources) if rna_bam_dir else {}
        files.extend(path for entries in rna_bams.values() for entry in entries for path in entry.values() if path)
        import annevo_gene_models
        inherited_annevo = anchors['request'].get('annevo') if anchor_root else None
        if inherited_annevo:
            annevo_gene_models.verify_inputs(inherited_annevo)
        annevo_table = annevo_inputs or (inherited_annevo['table'] if inherited_annevo else None)
        annevo = annevo_gene_models.freeze_inputs(annevo_table, prediction_sources) if annevo_table else None
        if annevo:
            annevo_gene_models.verify_inputs(annevo)
            files.extend(annevo['files'])
        inherited_oriongeno = anchors['request'].get('oriongeno') if anchor_root else None
        if inherited_oriongeno:
            annevo_gene_models.verify_inputs(inherited_oriongeno)
        oriongeno = (annevo_gene_models.freeze_inputs(oriongeno_inputs, prediction_sources, engine='oriongeno')
                     if oriongeno_inputs else copy.deepcopy(inherited_oriongeno))
        if oriongeno:
            files.extend(oriongeno['files'])
        genomic_prediction = bool(genomic_targets if pool else sources) and params['mode'] != 'off'
        tool = shutil.which('miniprot') if genomic_prediction else None
        reviews, publications, claimed = {}, {}, set()
        for key, approvals, publication, label in (
            ('split_model_merges', split_model_merges, 'split_model_candidates.json', 'Split model'),
            ('model_repairs', model_repairs, 'model_repair_candidates.json', 'Model repair'),
        ):
            if not approvals:
                continue
            if not receipts or params['mode'] != 'conservative':
                raise ValueError(label + ' approvals require completed --rescue-output and --mode conservative')
            selected, proof = freeze_model_reviews(approvals, anchor_root, receipts, publication, label)
            owners = {(name, owner['gene_id']) for name, candidates in selected.items()
                      for candidate in candidates for owner in candidate['owner_overlap']}
            if owners & claimed:
                raise ValueError('Approved repair and merge share a source gene: ' + str(sorted(owners & claimed)))
            claimed.update(owners)
            reviews[key] = selected
            publications.update(proof)
            files.extend(proof)
        if params['mode'] == 'conservative':
            for key, publication in (('automatic_model_repairs', 'model_repair_candidates.json'),
                                     ('automatic_split_merges', 'split_model_candidates.json')):
                automatic = {}
                for name, receipt_sha in sorted(receipts.items()):
                    path = anchor_root / 'rescued' / name / publication
                    rows, expected = read_model_publication(path, receipt_sha, required=False)
                    candidates = []
                    for row in rows:
                        if any((name, owner['gene_id']) in claimed for owner in row['owner_overlap']):
                            continue
                        support = row['support']
                        if key == 'automatic_model_repairs':
                            if row['kind'] != 'coding_path':
                                continue
                            support = [hit for hit in support if hit.get('syntenic_locus') or hit.get('reciprocal_locus')]
                        if len({hit['donor'] for hit in support} - {name}) >= 2:
                            candidates.append(row)
                    if candidates:
                        automatic[name] = candidates
                        publications[str(path)] = expected
                        files.append(str(path))
                if automatic:
                    reviews[key] = automatic
        request = {'schema': SCHEMA, 'sources': sources, 'files': transcriptome_proof.read(files), 'parameters': params,
                   'implementation': implementation(), 'dependencies': dependency_identities(), 'rescue_output': str(anchor_root) if anchor_root else None,
                   'edges': str(Path(edges).resolve()) if edges else None,
                   'rna': str(Path(rna).resolve()) if rna else None, 'species_profiles': profiles,
                   'revision_candidates': revision_candidates,
                   'revision_model_stores': revision_model_stores,
                   'revision_donor_maps': revision_donor_maps,
                   'miniprot': {'path': tool, 'sha256': digest(tool)} if tool else None}
        if any(request['files'][path] != expected for path, expected in publications.items()):
            raise ValueError('Reviewed candidate publication changed during planning')
        request.update(reviews)
        if rna_bams:
            request['rna_bams'] = rna_bams
        if annevo:
            request['annevo'] = annevo
        if oriongeno:
            request['oriongeno'] = oriongeno
        if transcriptome_sources:
            request['transcriptome_sources'] = transcriptome_sources
        if pool:
            request['target_species'] = targets
            request['fixed_donors'] = fixed_donors
            request['donor_plan_path'] = donor_plan_path
        if not edges and transcriptome_sources and len(sources) + len(transcriptome_sources) > 1:
            diamond = shutil.which('diamond')
            if not diamond:
                raise ValueError('DIAMOND is required to infer transcriptome correspondence')
            request['transcriptome_diamond'] = {'path': diamond, 'sha256': digest(diamond)}
        if not edges and not anchor_root and not transcriptome_sources:
            raise ValueError('A frozen synteny rescue plan or explicit trusted correspondence table is required')
        if genomic_prediction and not tool:
            raise ValueError('miniprot is required for prediction')
        transcriptome_proof.check()
        path = root / 'plan.json'
        if path.exists():
            previous = json.loads(path.read_text())
            if previous['request'] != request:
                raise ValueError('Frozen refinement plan differs; use a new output directory')
            return previous
        value = {'request': request, 'species': sorted(sources.keys() | transcriptome_sources.keys())}
        atomic_json(path, value, immutable=True)
        return value


def target_species(value):
    return value['request'].get('target_species', value['species'])


def donor_cache(request):
    return Path(request['donor_plan_path']).parent.parent / 'gene_model_refinement_donor_cache'


def donor_catalog_key(request, name):
    source = request['sources'][name]
    return {'schema': 1, 'species': name, 'genetic_code': source['genetic_code'],
            'files': {role: request['files'][source[role]] for role in ('fasta', 'gff', 'genome')},
            'implementation': request['implementation'], 'dependencies': request['dependencies']}


def donor_catalog_relative(request, name):
    key = json.dumps(donor_catalog_key(request, name), sort_keys=True).encode()
    return Path('catalog') / hashlib.sha256(key).hexdigest()


def catalog_directory(root, value, name):
    request = value['request']
    if name in request.get('fixed_donors', ()):
        return donor_cache(request) / donor_catalog_relative(request, name)
    return root / 'catalog' / name


def catalog_hashes(root, value):
    return {name: digest(catalog_directory(root, value, name) / 'receipt.json') for name in value['species']}


def load(root, names=None):
    root = Path(root)
    try:
        value = json.loads((root / 'plan.json').read_text())
    except (OSError, ValueError) as exc:
        if _INPUT_SNAPSHOT is not None:
            _INPUT_SNAPSHOT.reject_plan_read(exc)
        raise
    if _INPUT_SNAPSHOT is not None:
        _INPUT_SNAPSHOT.verify_plan(root, value)
    request = value['request']
    if request['schema'] != SCHEMA or request['implementation'] != implementation() or request['dependencies'] != dependency_identities():
        raise ValueError('Refinement implementation changed; use a new output directory')
    source_paths = {name: [source[k] for k in ('fasta', 'gff', 'genome')]
                    for name, source in request['sources'].items()}
    source_paths.update({name: [source[k] for k in ('fasta', 'transcript_gene_map')]
                         for name, source in request.get('transcriptome_sources', {}).items()})
    species_paths = set(chain.from_iterable(source_paths.values()))
    native_groups = [request.get(engine, {}).get('predictions', {}) for engine in ('annevo', 'oriongeno')]
    species_paths.update(prediction[k] for native in native_groups for prediction in native.values() for k in ('gff', 'receipt'))
    bam_paths = {n: {path for entry in entries for path in entry.values() if path}
                 for n, entries in request.get('rna_bams', {}).items()}
    species_paths.update(path for paths in bam_paths.values() for path in paths)
    files = request['files'] if names is None else {p: value for p, value in request['files'].items() if p not in species_paths}
    if names is not None:
        files.update({path: request['files'][path] for name in names for path in source_paths[name]})
        files.update({native[n][k]: request['files'][native[n][k]] for native in native_groups for n in names if n in native for k in ('gff', 'receipt')})
        files.update({path: request['files'][path] for n in names for path in bam_paths.get(n, ())})
    if _INPUT_SNAPSHOT is None:
        if digest_paths(files) != files:
            raise ValueError('Frozen refinement input changed')
    else:
        _INPUT_SNAPSHOT.verify_files(files)
    for key in ('miniprot', 'transcriptome_diamond'):
        tool = request.get(key)
        if not tool:
            continue
        if _INPUT_SNAPSHOT is not None:
            _INPUT_SNAPSHOT.verify_predictor(tool['path'], tool['sha256'])
        elif digest(tool['path']) != tool['sha256']:
            raise ValueError('Frozen predictor changed: ' + key)
    return value


def stage(root, relative, key, builder, names=None, prepared_snapshot=None):
    dependencies = key.get('dependencies', {})
    catalog_deps = dependencies.get('catalog', {})
    value = json.loads((root / 'plan.json').read_text()) if catalog_deps else None
    snapshots = {catalog_directory(root, value, name) / 'receipt.json': expected
                 for name, expected in catalog_deps.items()}
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
    if dependencies.get('rescue_workers'):
        request = json.loads((root / 'plan.json').read_text())['request']
        snapshots.update({revision_worker(request, name) / 'receipt.json': expected
                          for name, expected in dependencies['rescue_workers'].items()})
    snapshots = {str(path): expected for path, expected in snapshots.items()}
    prepared_receipts = ({str(prepared_snapshot.root / 'prepared' / name / 'receipt.json'): name
                          for name in dependencies.get('prepared', {})} if prepared_snapshot else {})
    ordinary = {path: expected for path, expected in snapshots.items() if path not in prepared_receipts}
    receipt_snapshot = _RECEIPT_SNAPSHOT

    def guard():
        load(root, names)
        if receipt_snapshot is None:
            rescue.require_same_key(ordinary, digest_paths(ordinary))
        else:
            receipt_snapshot.verify(ordinary)
        if prepared_snapshot is not None:
            try:
                prepared_snapshot.check_all()
                for path, name in prepared_receipts.items():
                    rescue.require_same_key(snapshots[path], prepared_snapshot.receipt_digest(name))
            except (OSError, ValueError) as exc:
                raise ValueError('Stage dependency content changed: prepared annotation') from exc
        for path in snapshots:
            if path in prepared_receipts or receipt_snapshot is not None:
                continue
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
    def publication_guard():
        if _INPUT_SNAPSHOT is not None:
            load(root, names)
        if receipt_snapshot is not None:
            receipt_snapshot.check(ordinary)
        if prepared_snapshot is not None:
            try:
                prepared_snapshot.check_all()
            except (OSError, ValueError) as exc:
                raise ValueError('Stage dependency content changed: prepared annotation') from exc

    options = {'publication_guard': publication_guard} if prepared_snapshot or receipt_snapshot or _INPUT_SNAPSHOT else {}
    if receipt_snapshot is not None:
        options['verify'] = receipt_snapshot.verify_reusable
    return rescue.stage(root, Path(relative), {'plan': digest(root / 'plan.json'), **key}, measured, guard, **options)


@invocation_cached
def catalog_species(root, value, name):
    request = value['request']
    if name in request.get('fixed_donors', ()):
        source = request['sources'][name]

        def build_donor(tmp):
            catalog = build_catalog(name, source['fasta'], source['gff'], source['genome'], source['genetic_code'])
            for locus in catalog['loci']:
                locus['fixed_representative'] = True
            write_catalog(catalog, tmp)

        def guard():
            load(root, [name])

        options = {'verify': _RECEIPT_SNAPSHOT.verify_reusable} if _RECEIPT_SNAPSHOT is not None else {}
        return rescue.stage(donor_cache(request), donor_catalog_relative(request, name),
                            donor_catalog_key(request, name), build_donor, guard, publication_guard=guard, **options)

    def build(tmp):
        transcriptome_source = value['request'].get('transcriptome_sources', {}).get(name)
        if transcriptome_source:
            write_catalog(transcriptome.build_catalog(name, transcriptome_source), tmp)
            return
        source = value['request']['sources'][name]
        catalog = build_catalog(name, source['fasta'], source['gff'], source['genome'], source['genetic_code'])
        entries = value['request'].get('rna_bams', {}).get(name, [])
        if value['request']['rna'] or entries:
            evidence = read_table(value['request']['rna']) if value['request']['rna'] else []
            evidence.extend(refinement_rna.evidence_rows(
                name, source['genome'], entries,
                (c for gene in catalog['loci'] if len(gene['candidates']) > 1 or any(
                    c.get('rescue_alternative_coding_path') for c in gene['candidates']) for c in gene['candidates']), tmp))
            rna = rna_path_index(evidence)
            for gene in catalog['loci']:
                for candidate in gene['candidates']:
                    paths = rna_support(candidate, name, rna)
                    candidate['quality']['rna_supported'] = bool(paths)
                    candidate.setdefault('support', {})['rna_paths'] = paths
        params = parameters_for(value['request'], name)
        for gene in catalog['loci']:
            for candidate in gene['candidates']:
                if candidate.get('rescue_alternative_coding_path'):
                    paths = candidate.get('support', {}).get('rna_paths', [])
                    status = 'accepted' if candidate['quality'].get('valid_orf') else 'proposal'
                    set_representative_admission({'candidate': candidate, 'rna_paths': paths, 'status': status}, True, params)
        write_catalog(catalog, tmp)
    return stage(root, Path('catalog') / name, {}, build, [name])


@invocation_cached
def catalog_index(root, value, predictions=False):
    dependencies = {}
    if predictions:
        initial_db = catalog_index(root, value)
        initial_receipt = (initial_db.parent / 'receipt.json').read_bytes()
        initial_record = json.loads(initial_receipt)
        initial_files = initial_record['files']
        dependencies['index'] = {'catalog_index': hashlib.sha256(initial_receipt).hexdigest()}
        for name in target_species(value):
            predict_species(root, value, name)
        dependencies['predictions'] = {n: digest(root / 'predictions' / n / 'receipt.json') for n in target_species(value)}
    else:
        directories = [catalog_species(root, value, name) for name in value['species']]
    dependencies['catalog'] = catalog_hashes(root, value)
    if predictions and initial_record['key']['dependencies']['catalog'] != dependencies['catalog']:
        raise ValueError('Initial refinement index catalog receipts changed before final indexing')

    def build(tmp):
        if predictions:
            raw_summary = (initial_db.parent / 'summary.json').read_bytes()
            if hashlib.sha256(raw_summary).hexdigest() != initial_files['summary.json']:
                raise ValueError('Initial refinement index summary changed during copying')
            summary = json.loads(raw_summary)
            overlays = catalog_prediction_overlays(root, target_species(value))
            first_overlay = next(overlays, None)
            summary['accepted_prediction_count'] = 0
        # Without accepted predictions a verified sequential copy is all the
        # final index needs. Keep SQLite updates on local scratch when present.
        scratch = (contextlib.nullcontext(tmp) if predictions and first_overlay is None
                   else scratch_build_directory(tmp))
        with scratch as work:
            db = work / 'loci.sqlite3'
            if predictions:
                copy_refinement_index(initial_db, db, initial_files[db.name])
                if first_overlay is not None:
                    apply_catalog_predictions(db, summary, chain((first_overlay,), overlays))
            else:
                summary = build_store(directories, db)
            if work != tmp:
                copy_refinement_index(db, tmp / db.name)
            if work != tmp or predictions:
                with (tmp / db.name).open('rb') as handle:
                    os.fsync(handle.fileno())
        summary['database'] = str((root / ('catalog_index_final' if predictions else 'catalog_index') / 'loci.sqlite3').resolve())
        atomic_json(tmp / 'summary.json', summary)
    directory = stage(root, 'catalog_index_final' if predictions else 'catalog_index', {'dependencies': dependencies}, build)
    return directory / 'loci.sqlite3'


def catalog_prediction_overlays(root, species):
    for name in species:
        by_gene = defaultdict(list)
        for model in stream_json_array(root / 'predictions' / name / 'predictions.json'):
            if model['status'] == 'accepted':
                by_gene[model['gene_id']].append(model['candidate'])
        for gene, candidates in by_gene.items():
            yield name, gene, candidates


def apply_catalog_predictions(db, summary, overlays):
    accepted_count = 0
    with contextlib.closing(sqlite3.connect(db)) as connection, connection:
        for name, gene, candidates in overlays:
            locus = load_locus(connection, name, gene)
            locus['candidates'].extend(candidates)
            connection.executemany('INSERT INTO candidate_owners(species,candidate_id,gene_id) VALUES(?,?,?)',
                                   ((name, candidate['candidate_id'], gene) for candidate in candidates))
            serialized = json.dumps(locus, sort_keys=True)
            connection.execute('UPDATE loci SET json=? WHERE species=? AND gene_id=?', (serialized, name, gene))
            accepted_count += len(candidates)
            summary['max_locus_json_bytes'] = max(summary['max_locus_json_bytes'], len(serialized.encode()))
    summary['accepted_prediction_count'] = accepted_count
    summary['candidate_count'] += accepted_count


def copy_refinement_index(source, destination, expected_sha256=None):
    """Copy with a checksum in one pass and verify the expected source bytes."""
    checksum = hashlib.sha256()
    with source.open('rb') as original, destination.open('wb') as copy_file:
        for chunk in iter(lambda: original.read(1024 * 1024), b''):
            copy_file.write(chunk)
            checksum.update(chunk)
    if expected_sha256 is not None and checksum.hexdigest() != expected_sha256:
        raise OSError('Copied refinement SQLite index differs from its source database')


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


def uses_rescue_correspondence(request):
    return bool(request['rescue_output'] and (not request['edges'] or request.get('transcriptome_sources')))


@invocation_cached
def correspondence(root, value, cpus=1, comparison_cache=None):
    db = catalog_index(root, value)
    request = value['request']
    anchor_root = Path(request['rescue_output']) if request['rescue_output'] else None
    if uses_rescue_correspondence(request):
        anchor_plan = rescue.load(anchor_root, inspect_only=True)
        with PreparedSnapshot(anchor_root, anchor_plan, rescue) as snapshot:
            return _correspondence(root, value, db, cpus, comparison_cache, anchor_plan, snapshot)
    return _correspondence(root, value, db, cpus, comparison_cache)


def _correspondence(root, value, db, cpus, comparison_cache, anchor_plan=None, prepared_snapshot=None):
    request = value['request']
    dependencies = {'catalog': catalog_hashes(root, value),
                    'index': {db.parent.name: digest(db.parent / 'receipt.json')}}
    anchor_root = Path(request['rescue_output']) if request['rescue_output'] else None
    jobs = []
    if anchor_plan is not None:
        jobs = [j for j in anchor_plan['synteny_jobs'] if j['a'] != j['b']]
        dependencies['anchors'] = {}
        for job in jobs:
            directory = anchor_root / 'synteny' / job['id']
            if 'fixed_donors' not in request:
                rescue.synteny(anchor_root, anchor_plan, job['index'], cpus, comparison_cache,
                               prepared_snapshot=prepared_snapshot)
            receipt = (directory / 'receipt.json').read_bytes()
            if json.loads(receipt).get('key') != rescue.comparison_key(anchor_root, job, prepared_snapshot):
                raise ValueError('Unverified synteny comparison: ' + job['id'])
            dependencies['anchors'][job['id']] = hashlib.sha256(receipt).hexdigest()
        dependencies['prepared'] = {}
        anchor_species = sorted({job[side] for job in jobs for side in ('a', 'b')})
        if 'fixed_donors' not in request:
            anchor_species = sorted(request['sources'])
        for name in anchor_species:
            prepared_snapshot.prepared(name)
            dependencies['prepared'][name] = prepared_snapshot.receipt_digest(name)
    def build(tmp):
        edges = []
        sequence_edges = set()
        if request['edges']:
            valid = defaultdict(set)
            for species, gene in iter_locus_keys(db):
                valid[species].add(gene)
            for row in read_table(request['edges']):
                a, b, ga, gb = (row[k] for k in ('species_a', 'species_b', 'gene_a', 'gene_b'))
                if a == b or ga not in valid.get(a, set()) or gb not in valid.get(b, set()):
                    raise ValueError('Correspondence references absent locus or same species: ' + str(row))
                weight = float(row.get('weight') or 1)
                if not math.isfinite(weight) or weight <= 0 or weight > 1:
                    raise ValueError('Invalid correspondence weight')
                sequence_only = any(n in request.get('transcriptome_sources', {}) for n in (a, b))
                ambiguity = (row.get('ambiguous') or '').strip().lower()
                if sequence_only:
                    pair = tuple(sorted(((a, ga), (b, gb))))
                    if pair in sequence_edges or ambiguity not in {'', '0', '1', 'false', 'true', 'no', 'yes'}:
                        raise ValueError('Duplicate or malformed transcriptome correspondence: ' + str(row))
                    sequence_edges.add(pair)
                edges.append(dict(species_a=a, gene_a=ga, species_b=b, gene_b=gb, weight=weight,
                                  evidence=row.get('evidence') or ('user_frozen_orthology' if sequence_only else 'user_frozen_synteny'),
                                  ambiguous=ambiguity in {'1', 'true', 'yes'}))
                if sequence_only:
                    edges[-1]['comparison_basis'] = 'protein_alignment'
        if anchor_plan is not None:
            aliases, position_index = {}, {}
            for n in anchor_species:
                prepared_snapshot.check([n])
                lookup, positions, tracks = {}, {}, defaultdict(list)
                for g in iter_loci(db, n):
                    for token in {g['gene_id'], g['gene_id'].removeprefix(n + '_'),
                                  *(c['source_gene_id'] for c in g['candidates']), *(c.get('gene_token', '') for c in g['candidates'])}:
                        if token in lookup and lookup[token] != g['gene_id']:
                            raise ValueError('Ambiguous gene alias: ' + token)
                        lookup[token] = g['gene_id']
                    blocks = [b for c in g['candidates'] for b in c['blocks']]
                    if not blocks:
                        continue
                    position = dict(gene_id=g['gene_id'], seqid=g['seqid'], start=min(b[0] for b in blocks),
                                    end=max(b[1] for b in blocks))
                    positions[g['gene_id']] = position
                    tracks[g['seqid']].append(position)
                aliases[n] = {r['jcvi_id']: lookup.get(r['locus_id']) or lookup.get(r['original_id'])
                              for r in read_table(anchor_root / 'prepared' / n / 'genes.id_map.tsv') if r['status'] == 'selected'}
                ordered = {seqid: sorted(genes, key=lambda g: g['start']) for seqid, genes in tracks.items()}
                position_index[n] = positions, {seqid: ([g['start'] for g in genes], genes) for seqid, genes in ordered.items()}
            with store_connection(db) as connection:
                for job in jobs:
                    directory = anchor_root / 'synteny' / job['id']
                    if _RECEIPT_SNAPSHOT is None:
                        blocks = json.loads((directory / 'blocks.json').read_text())
                    else:
                        present, blocks = _RECEIPT_SNAPSHOT.read_json_member(
                            directory / 'receipt.json', dependencies['anchors'][job['id']], 'blocks.json')
                        if not present:
                            raise ValueError('Synteny comparison has no receipted blocks')
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
        if request.get('transcriptome_diamond'):
            with contextlib.ExitStack() as stack:
                scratch = Path(stack.enter_context(tempfile.TemporaryDirectory(prefix='transcriptome-proteins-', dir=tmp)))
                executable = request['transcriptome_diamond']['path']
                if 'target_species' in request:
                    target, = request['target_species']
                    donors = request['fixed_donors']
                    edges.extend(refinement_transcriptome_edges.fixed_panel_edges(
                        db, target, donors, executable, request['parameters'], scratch, cpus, donor_cache(request),
                        {n: donor_catalog_key(request, n) for n in donors}, lambda: load(root)))
                else:
                    edges.extend(refinement_transcriptome_edges.infer_edges(
                        db, value['species'], request['transcriptome_sources'], executable,
                        request['parameters'], scratch, cpus))
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
        atomic_json(tmp / 'donor_species.json', correspondence_donor_species(edges))
        atomic_json(tmp / 'summary.json', {'edges': len(edges), 'ambiguous': sum(e['ambiguous'] for e in edges)})
    return stage(root, 'correspondence', {'dependencies': dependencies}, build,
                 prepared_snapshot=prepared_snapshot)


def correspondence_donor_species(edges):
    """Project the final nonambiguous edges without changing source guards."""
    donors = defaultdict(set)
    for edge in edges:
        if not edge['ambiguous']:
            a, b = edge['species_a'], edge['species_b']
            donors[a].add(b)
            donors[b].add(a)
    return {name: sorted(values) for name, values in sorted(donors.items())}


def prediction_donor_names(directory, name, species):
    """Use a verified additive projection; old publications retain edge scan."""
    path = directory / 'receipt.json'
    expected = digest(path)
    if _RECEIPT_SNAPSHOT is not None:
        present, donors = _RECEIPT_SNAPSHOT.read_json_member(path, expected, 'donor_species.json')
    else:
        # Outside an explicit command, keep ordinary full verification and
        # discard this single-use read fence immediately.
        with ReceiptSnapshot(rescue.verified) as proof:
            present, donors = proof.read_json_member(path, expected, 'donor_species.json')
    if not present:
        donors = correspondence_donor_species(json.loads((directory / 'edges.json').read_text()))
    if (not isinstance(donors, dict) or not set(donors) <= set(species)
            or any(not isinstance(values, list) or any(not isinstance(value, str) for value in values)
                   or values != sorted(set(values)) or not set(values) <= set(species)
                   for values in donors.values())):
        raise ValueError('Invalid correspondence donor species projection')
    return sorted({name, *donors.get(name, [])})


@invocation_cached
def select(root, value, predictions=False):
    db = catalog_index(root, value, predictions)
    directory = correspondence(root, value)
    params = value['request']['parameters']
    dependencies = {'correspondence': digest(directory / 'receipt.json'),
                    'catalog': catalog_hashes(root, value),
                    'index': {db.parent.name: digest(db.parent / 'receipt.json')}}
    if predictions:
        dependencies['predictions'] = {n: digest(root / 'predictions' / n / 'receipt.json') for n in target_species(value)}
        # Assembly-only targets have no genomic prediction or repair.
        if set(target_species(value)) <= value['request'].get('transcriptome_sources', {}).keys():
            initial = select(root, value)
            dependencies['initial'] = digest(initial / 'receipt.json')

    def build(tmp):
        if 'initial' in dependencies:
            if json.loads((db.parent / 'summary.json').read_text())['accepted_prediction_count'] != 0:
                raise ValueError('Assembly-only target has an unexpected genomic prediction')
            for filename in ('selection.json', 'representative_map.tsv'):
                shutil.copyfile(initial / filename, tmp / filename)
            return
        edges = json.loads((directory / 'edges.json').read_text())
        # Graph components read loci in an order that is expensive on shared
        # storage. Verify one sequential copy before opening it for selection.
        scratch_context = (scratch_build_directory(db.parent) if 'fixed_donors' in value['request']
                           else contextlib.nullcontext(db.parent))
        with scratch_context as scratch:
            selection_db = db
            if scratch != db.parent:
                receipt = (db.parent / 'receipt.json').read_bytes()
                if hashlib.sha256(receipt).hexdigest() != dependencies['index'][db.parent.name]:
                    raise ValueError('Refinement index receipt changed before selection')
                selection_db = scratch / db.name
                copy_refinement_index(db, selection_db, json.loads(receipt)['files'][db.name])
            selection = select_from_store(selection_db, edges, policy=params['policy'], min_margin=params['min_margin'],
                                          min_support=params['min_support'], candidate_limit=params['candidate_limit'])
        if predictions:
            repairs = (row for name in target_species(value)
                       for row in stream_json_array(root / 'predictions' / name / 'predictions.json')
                       if row.get('same_coordinate_source_repair'))
            selection = apply_same_coordinate_repairs(selection, repairs)
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
        blocks = row['cds_blocks']
        if isinstance(blocks, str):
            blocks = json.loads(blocks)
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
            parsed = parse_gff_attributes(f[8], keys=OWNERSHIP_ATTRIBUTES)
            attr = {key: values[0] for key, values in parsed.items() if values}
            identifier = attr.get('ID', attr.get('transcript_id', attr.get('gene_id', '')))
            declared = parsed.get('Parent', ()) + parsed.get('gene_id', ()) + parsed.get('transcript_id', ())
            if (f[2].lower().endswith(('gene', 'rna', 'transcript')) or f[2] in {'CDS', 'exon'}
                    or identifier in owners or parsed.get('gene_id') or parsed.get('transcript_id')
                    or any(token in owners for token in declared)):
                fallback = attr.get('gene_id') or next(iter(parsed.get('Parent', ())), identifier)
                owner = owners.get(identifier) or next((owners[p] for p in declared if p in owners), 'annotation:' + fallback)
                # Opposite-strand annotation is not ownership of this coding
                # path. Unknown strand retains the conservative span gate.
                strand = f[6] if f[6] in {'+', '-'} else '.'
                key = f[0], owner, strand
                start, end = int(f[3]) - 1, int(f[4])
                previous = spans.get(key, (start, end))
                spans[key] = min(previous[0], start), max(previous[1], end)
    coding = defaultdict(list)
    for gene in catalog['loci']:
        if gene.get('ambiguous_coordinates') or any(c.get('quality', {}).get('structure_problem') for c in gene['candidates']):
            continue
        for candidate in gene['candidates']:
            coding[gene['gene_id'], candidate.get('seqid', gene['seqid']),
                   candidate.get('strand', gene['strand'])].extend(candidate['blocks'])
    result = []
    for (seqid, owner, strand), (start, end) in sorted(spans.items()):
        row = {'seqid': seqid, 'start': start, 'end': end, 'gene_id': owner, 'strand': strand}
        # For a reconstructed coding locus, an intron/UTR overlap alone does
        # not fuse coding models. Noncoding or unbound owners keep full spans.
        blocks = coding.get((owner, seqid, strand), []) if strand != '.' else []
        if blocks:
            row['cds_blocks'] = sorted({tuple(b[:2]) for b in blocks})
        result.append(row)
    return result


def set_representative_admission(row, valid_original, params):
    """Set the independent adoption gate without changing prediction evidence."""
    conservation_adoption = (valid_original and not row['rna_paths']
                             and params.get('isoform_adoption', 'rna_required') == 'conservation_supported'
                             and params['policy'] == 'conserved')
    quality = row['candidate']['quality']
    quality['representative_eligible'] = row['status'] == 'accepted' and (bool(row['rna_paths']) or not valid_original or conservation_adoption)
    quality.pop('representative_admission', None)
    if conservation_adoption:
        quality['representative_admission'] = 'conservation_supported'


def same_coordinate_source_owner(model, gene, genetic_code):
    """Identify one source-invalid, otherwise complete exact genomic path."""
    if gene.get('ambiguous_coordinates') or not isinstance(genetic_code, int) or isinstance(genetic_code, bool):
        return None
    shape = (model['seqid'], model['strand'], tuple(tuple(b) for b in model['cds']))
    matches = [c for c in gene['candidates'] if
               (c.get('seqid', gene['seqid']), c.get('strand', gene['strand']),
                tuple(tuple(b) for b in c['blocks'])) == shape]
    if len(matches) != 1:
        return None
    original = matches[0]
    quality = original.get('quality', {})
    sources = original.get('source_cds', [])
    if (original.get('origin') != 'original' or not quality.get('sequence_mismatch')
            or original.get('cds') != model.get('sequence') or len(sources) != 1
            or not sources[0].get('cds') or sources[0].get('sequence_agreement') is not False
            or sources[0]['cds'] == original['cds']):
        return None
    prohibited = ('partial', 'annotated_partial', 'phase_conflict', 'phase_unknown',
                  'phase_unresolved', 'internal_stop', 'ambiguous', 'invalid_base',
                  'incomplete_codon', 'translation_uncertain', 'structure_problem',
                  'annotated_exception', 'annotated_pseudogene', 'translation_exception',
                  'sequence_exception', 'frameshift', 'ownership_conflict')
    if any(quality.get(flag) for flag in prohibited):
        return None
    probe = copy.deepcopy(original)
    probe['quality']['sequence_mismatch'] = False
    try:
        fresh = validate_candidate(probe, genetic_code)
    except (TypeError, ValueError):
        return None
    if not fresh.get('valid_orf') or any(fresh.get(flag) for flag in prohibited):
        return None
    return original


def apply_same_coordinate_repairs(selection, prediction_rows):
    """Update a newly built selection in place, retaining bounded graph audit."""
    repairs = defaultdict(list)
    for row in prediction_rows:
        proof = row.get('same_coordinate_source_repair')
        candidate = row.get('candidate', {})
        if (row.get('status') != 'accepted' or not isinstance(proof, dict)
                or proof.get('schema') != 1 or not proof.get('baseline_candidate_id')
                or proof.get('minimum_external_species', 0) < 2
                or len(set(proof.get('independent_external_species', []))) < proof['minimum_external_species']
                or not candidate.get('quality', {}).get('valid_orf')
                or not candidate.get('quality', {}).get('representative_eligible')
                or hashlib.sha256(candidate.get('cds', '').encode()).hexdigest() != proof.get('genomic_cds_sha256')
                or candidate.get('source_transcript_id') != proof.get('source_transcript_id')):
            continue
        shape = (candidate.get('seqid'), candidate.get('strand'),
                 tuple(tuple(b) for b in candidate.get('blocks', [])))
        if hashlib.sha256(json.dumps(shape).encode()).hexdigest() != proof.get('coding_key'):
            continue
        probe = copy.deepcopy(candidate)
        try:
            fresh = validate_candidate(probe, proof.get('genetic_code'))
        except (TypeError, ValueError):
            continue
        if not fresh.get('valid_orf') or not fresh.get('usable'):
            continue
        repairs[proof['species'], row['gene_id']].append(row)
    applied = {}
    for decision in selection['selections']:
        matches = repairs.get((decision['species'], decision['gene_id']), [])
        if len(matches) != 1:
            continue
        row = matches[0]
        proof = row['same_coordinate_source_repair']
        if decision['candidate_id'] not in {proof['original_candidate_id'], row['candidate']['candidate_id']}:
            continue
        previous = copy.deepcopy(decision)
        baseline = proof['baseline_candidate_id']
        applied[decision['species'], decision['gene_id']] = {
            'species': decision['species'], 'gene_id': decision['gene_id'],
            'original_candidate_id': decision['candidate_id'],
            'final_candidate_id': row['candidate']['candidate_id'],
            'baseline_candidate_id': baseline, 'source_selected_flags': [],
            'change_count_delta': (int(row['candidate']['candidate_id'] != baseline)
                                   - int(decision['candidate_id'] != baseline))}
        decision.update(candidate_id=row['candidate']['candidate_id'], status='sequence_repaired',
                        reason='trusted_same_coordinate_genome_cds',
                        source_correspondence_decision=previous,
                        same_coordinate_source_repair=copy.deepcopy(proof),
                        orthology='unassigned', expected_copy='unassigned')
    if not applied:
        return selection
    # Retain only original flags for repaired loci, not a second scores table.
    for score in selection.get('scores', []):
        audit = applied.get((score['species'], score['gene_id']))
        if audit is not None:
            audit['source_selected_flags'].append(
                {'candidate_id': score['candidate_id'], 'selected': score['selected']})
            score['selected'] = score['candidate_id'] == audit['final_candidate_id']
    if isinstance(selection.get('metrics'), dict):
        selection['source_correspondence_metrics'] = dict(selection['metrics'])
        selection['metrics'] = dict(selection['metrics'])
        if 'changed_representatives' in selection['metrics']:
            selection['metrics']['changed_representatives'] += sum(
                audit['change_count_delta'] for audit in applied.values())
        selection['metrics']['same_coordinate_sequence_repairs'] = len(applied)
    selection['same_coordinate_repair_audit'] = list(applied.values())
    return selection


def classify_predictions(models, catalog, edges, params, rna_rows, genome_hash):
    """Require target ownership, intact genomic ORF and independent support."""
    rna_rows = rna_path_index(rna_rows)
    models = iter(models)
    try:
        first_model = next(models)
    except StopIteration:
        return []
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
            spans.append((g['seqid'], min(b[0] for b in blocks), max(b[1] for b in blocks), g['gene_id'],
                          '.' if g.get('ambiguous_coordinates') else g.get('strand', '.'), (() if g.get('ambiguous_coordinates') or
                            any(c.get('quality', {}).get('structure_problem') for c in g['candidates'])
                            else tuple(tuple(b[:2]) for b in blocks))))
    spans.extend((r['seqid'], r['start'], r['end'], r['gene_id'], r.get('strand', '.'),
                  tuple(tuple(b[:2]) for b in r.get('cds_blocks', ()))) for r in catalog.get('annotation_spans', []))
    tracks = defaultdict(list)
    for seqid, start, end, owner, strand, blocks in spans:
        tracks[seqid].append((start, end, owner, strand, blocks))
    indexes = {}
    for seqid, rows in tracks.items():
        ordered = sorted(rows)
        maxima, maximum = [], 0
        for _, end, *_ in ordered:
            maximum = max(maximum, end)
            maxima.append(maximum)
        indexes[seqid] = [r[0] for r in ordered], ordered, maxima
    def foreign_overlap(seqid, strand, blocks, start, end, owner):
        starts, rows, maxima = indexes.get(seqid, ([], [], []))
        i = bisect.bisect_left(starts, end) - 1
        while i >= 0 and maxima[i] > start:
            a, b, gene, other_strand, other_blocks = rows[i]
            if (gene != owner and a < end and b > start and (other_strand not in {'+', '-'} or other_strand == strand)
                    and (other_strand not in {'+', '-'} or not other_blocks
                         or any(x[0] < y[1] and x[1] > y[0] for x in blocks for y in other_blocks))):
                return True
            i -= 1
        return False
    grouped = {}
    for m in chain((first_model,), models):
        g = loci[m['gene_id']]
        candidate = {'seqid': m['seqid'], 'strand': m['strand'], 'blocks': m['cds'], 'cds': m['sequence'],
                     'origin': 'predicted', 'source_gene_id': g['candidates'][0]['source_gene_id'],
                     'gene_token': g.get('gene_token', g['candidates'][0].get('gene_token', g['gene_id'])), 'junctions': []}
        shape = (candidate['seqid'], candidate['strand'], tuple(tuple(b) for b in candidate['blocks']))
        original_shapes = {(c.get('seqid', g['seqid']), c.get('strand', g['strand']), tuple(tuple(b) for b in c['blocks'])) for c in g['candidates']}
        # Alignment quality belongs to this donor, not the shared coding path.
        # Failed alignments cannot veto a supported path or count as support.
        alignment_problems = sorted(set(m['problems']) & {'low_coverage', 'low_identity',
                                                        'low_reciprocal_coverage', 'annevo_alignment_bounded',
                                                        'annevo_unconfirmed_coding_structure', 'oriongeno_alignment_bounded',
                                                        'oriongeno_unconfirmed_coding_structure'})
        problems = sorted(set(m['problems']) - set(alignment_problems))
        if g.get('ambiguous_coordinates'):
            problems.append('ambiguous_locus_coordinates')
        if m['gene_id'] in ambiguous_loci:
            problems.append('ambiguous_locus_correspondence')
        source_repair_owner = same_coordinate_source_owner(m, g, catalog.get('genetic_code'))
        if m.get('same_coordinate_repair_only') and source_repair_owner is None:
            problems.append('same_coordinate_repair_required')
        if shape in original_shapes:
            problems.append('existing_coding_path')
        if candidate['strand'] != g['strand'] or candidate['seqid'] != g['seqid']:
            problems.append('wrong_locus_strand')
        start, end = min((b[0] for b in m['cds']), default=0), max((b[1] for b in m['cds']), default=0)
        if foreign_overlap(m['seqid'], m['strand'], m['cds'], start, end, m['gene_id']):
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
        rna_confirmed_structure = (candidate['quality']['rna_supported']
                                   and any(p.endswith('_unconfirmed_coding_structure') for p in alignment_problems))
        if rna_confirmed_structure:
            alignment_problems = [p for p in alignment_problems if not p.endswith('_unconfirmed_coding_structure')]
        from gene_model_source_retention import source_structure_retention
        retention = source_structure_retention(candidate, g['candidates'], catalog.get('genetic_code', 1))
        explicit_intact_adoption = (params.get('isoform_adoption') == 'conservation_supported'
                                   and any(c['quality'].get('valid_orf') for c in g['candidates']))
        if not candidate['quality']['rna_supported'] and not retention['preserved'] and not explicit_intact_adoption:
            problems.append('unconfirmed_source_structure_change')
        if not candidate['quality']['rna_supported'] and any(c.get('quality', {}).get('sequence_mismatch') for c in owners):
            problems.append('unresolved_source_sequence_mismatch')
        key = m['gene_id'], shape
        if key not in grouped:
            grouped[key] = {'gene_id': m['gene_id'], 'candidate': candidate, 'problems': problems,
                            'donors': [], 'alignments': [], 'rna_paths': rna_support(candidate, catalog['species'], rna_rows),
                            'source_structure_retention': retention}
        row = grouped[key]
        supported = not alignment_problems and 'untrusted_donor_correspondence' not in problems
        if supported:
            row['donors'].append(m['donor_species'])
        if (supported and source_repair_owner is not None
                and m['donor_species'] != catalog['species']
                and 'minimum_identity' in params and 'minimum_coverage' in params
                and homology_supported(m['coverage'], m['identity'], params)):
            row.setdefault('_same_coordinate_repair_support', []).append(
                (m['donor_species'], source_repair_owner))
        row['alignments'].append({'donor_species': m['donor_species'], 'donor_candidate': m['donor_candidate'],
                                  'identity': m['identity'], 'coverage': m['coverage'],
                                  'supports_path': supported, 'problems': alignment_problems,
                                  **({'terminal_completion': m['terminal_completion']} if m.get('terminal_completion') else {}),
                                  **({'rescue_revision_model_id': m['rescue_revision_model_id']} if m.get('rescue_revision_model_id') else {})})
        for engine in ('annevo', 'oriongeno'):
            if m.get('evidence', {}).get(engine):
                row['alignments'][-1][engine] = m['evidence'][engine]
                if rna_confirmed_structure:
                    row['alignments'][-1]['structure_confirmation'] = 'whole_coding_RNA'
        row['problems'] = sorted(set(row['problems'] + problems))
    result = []
    for row in grouped.values():
        row['donors'] = sorted(set(row['donors']))
        repair_support = row.pop('_same_coordinate_repair_support', [])
        external = sorted({donor for donor, _ in repair_support})
        minimum = max(2, params['min_support'])
        exceptional = {'existing_coding_path', 'unresolved_source_sequence_mismatch',
                       'ambiguous_locus_correspondence'}
        if (len(external) >= minimum and params['mode'] == 'conservative'
                and params['policy'] == 'conserved'
                and not set(row['problems']) - exceptional
                and row['candidate']['quality'].get('valid_orf')
                and row['candidate']['cds'] == repair_support[0][1]['cds']):
            owner = repair_support[0][1]
            row['same_coordinate_source_repair'] = {
                'schema': 1, 'species': catalog['species'],
                'original_candidate_id': owner['candidate_id'],
                'baseline_candidate_id': _baseline(loci[row['gene_id']], loci[row['gene_id']]['candidates'])['candidate_id'],
                'source_transcript_id': owner['source_transcript_id'],
                'supplied_cds_sha256': hashlib.sha256(owner['source_cds'][0]['cds'].encode()).hexdigest(),
                'genomic_cds_sha256': hashlib.sha256(row['candidate']['cds'].encode()).hexdigest(),
                'assembly_sha256': genome_hash, 'minimum_external_species': minimum,
                'genetic_code': catalog['genetic_code'], 'coding_key': row['candidate']['coding_key'],
                'independent_external_species': external,
                'correspondence_ambiguity': 'ambiguous_locus_correspondence' in row['problems'],
                'orthology': 'unassigned', 'expected_copy': 'unassigned'}
            row['candidate']['source_transcript_id'] = owner['source_transcript_id']
            row['candidate']['same_coordinate_source_repair'] = copy.deepcopy(row['same_coordinate_source_repair'])
            row['candidate']['gene_id'] = row['gene_id']
            row['problems'] = sorted(set(row['problems']) - exceptional)
        if not row['donors']:
            row['problems'].append('no_qualifying_donor_alignment')
        if len(row['donors']) < params['min_support'] and not row['rna_paths']:
            row['problems'].append('insufficient_independent_support')
        row['status'] = 'accepted' if not row['problems'] and params['mode'] == 'conservative' else 'proposal'
        valid_original = any(c['quality'].get('valid_orf') for c in loci[row['gene_id']]['candidates'])
        row['change_type'] = ('model_revision' if row.get('same_coordinate_source_repair')
                              else 'isoform_addition' if valid_original else 'model_revision')
        row['evidence_class'] = ('homology_supported_same_coordinate_repair'
                                 if row.get('same_coordinate_source_repair') else
                                 'rna_path_supported' if row['rna_paths'] else 'homology_only_predicted')
        set_representative_admission(row, valid_original, params)
        row['candidate']['support'] = {'donors': row['donors'], 'rna_paths': row['rna_paths'], 'class': row['evidence_class']}
        for engine in ('annevo', 'oriongeno'):
            if any(a.get(engine) for a in row['alignments']):
                row['candidate']['support'][engine] = [a[engine] for a in row['alignments'] if a.get(engine)]
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
            if (other['gene_id'] != row['gene_id'] and other['candidate']['strand'] == c['strand']
                    and any(a[0] < b[1] and a[1] > b[0] for a in c['blocks'] for b in other['candidate']['blocks'])):
                for conflicting in (row, other):
                    conflicting['status'] = 'proposal'
                    if 'competing_predictions_overlap' not in conflicting['problems']:
                        conflicting['problems'].append('competing_predictions_overlap')
                    conflicting['candidate']['quality']['representative_eligible'] = False
        active.append((start, end, row))
    return sorted(result, key=lambda r: (r['gene_id'], r['candidate']['candidate_id']))


def rescue_revision_donor_resolver(request, connection):
    """Resolve a frozen rescue protein query to its exact catalog locus/path."""
    lookups, proteins, resolved = {}, {}, {}
    def resolve(donor, query):
        key = donor, query
        if key in resolved:
            return resolved[key]
        metadata = request.get('revision_donor_maps', {}).get(donor)
        if not metadata:
            return None
        if donor not in lookups:
            lookups[donor] = {row['jcvi_id']: row for row in read_table(metadata['mapping']) if row['status'] == 'selected'}
            proteins[donor] = {identifier: sequence for identifier, _, sequence in rescue.fasta_records(Path(metadata['protein']))}
        row = lookups[donor].get(query)
        reference = proteins[donor].get(query)
        found = {}
        if row and reference:
            for token in {row.get('locus_id', ''), row.get('original_id', '')} - {''}:
                try:
                    locus = load_locus(connection, donor, _formatted_id(donor, token))
                except KeyError:
                    continue
                found[locus['gene_id']] = locus
        result = None
        if len(found) == 1:
            locus = next(iter(found.values()))
            agreeing = [c for c in locus['candidates'] if c['quality'].get('usable') and c.get('protein') == reference]
            bound = [c for c in agreeing if row['original_id'] in c.get('source_fasta_ids', [])]
            agreeing = bound or agreeing
            baseline = locus.get('source_baseline_candidate_id') or locus.get('source_baseline_coding_candidate_id')
            selected = next((c for c in agreeing if c['candidate_id'] == baseline), None)
            if selected is None and agreeing and len({c['coding_key'] for c in agreeing}) == 1:
                selected = min(agreeing, key=lambda c: c['candidate_id'])
            if selected is not None and not locus.get('ambiguous_coordinates'):
                result = locus['gene_id'], selected['candidate_id']
        resolved[key] = result
        return result
    return resolve


def import_rescue_revisions(models, catalog, edges, params, genome, resolve_donor):
    """Nominate existing loci without converting missing-gene support to trust.

    Only the original ownership-routing decision is removed. The current
    genome, protected annotations, correspondence, donor-specific alignment
    metrics and ordinary refinement acceptance/adoption gates still apply.
    """
    owners = defaultdict(set)
    for gene in catalog['loci']:
        for token in {gene['gene_id'], gene.get('source_gene_id', ''), gene.get('gene_token', ''),
                      *(c.get('source_gene_id', '') for c in gene['candidates']),
                      *(c.get('gene_token', '') for c in gene['candidates'])} - {''}:
            owners[token].add(gene['gene_id'])
    trusted = defaultdict(set)
    for edge in edges:
        if edge.get('ambiguous'):
            continue
        for side, other in (('a', 'b'), ('b', 'a')):
            if edge['species_' + side] == catalog['species']:
                trusted[edge['gene_' + side], edge['species_' + other]].add(edge['gene_' + other])
    validated, proposals = [], []
    for raw in models:
        nomination = {'model_id': raw.get('model_id', ''), 'status': 'rescue_revision_nomination',
                      'revision_owner_ids': raw.get('revision_owner_ids', [])}
        declared = raw.get('revision_owner_ids', [])
        mapped = set().union(*(owners.get(token, set()) for token in declared)) if declared else set()
        if len(declared) != 1 or len(mapped) != 1:
            proposals.append({**nomination, 'reason': 'ambiguous_or_absent_rescue_revision_owner'})
            continue
        gene_id = next(iter(mapped))
        original_problems = set(raw.get('problems', [])) - {'overlap_existing_annotation', 'low_coverage', 'low_identity'}
        evidence_rows = [(raw.get('evidence', {}), raw)]
        primary = raw.get('evidence', {})
        for evidence in raw.get('support', []):
            if (evidence.get('donor'), evidence.get('query')) == (primary.get('donor'), primary.get('query')):
                continue
            metrics = evidence.get('alignment', evidence)
            if not all(key in metrics for key in ('coverage', 'identity', 'problems')):
                proposals.append({**nomination, 'gene_id': gene_id, 'donor': evidence.get('donor'),
                                  'reason': 'rescue_support_without_alignment_metrics'})
                continue
            evidence_rows.append((evidence, metrics))
        for evidence, metrics in evidence_rows:
            donor, query = evidence.get('donor'), evidence.get('query')
            correspondence = resolve_donor(donor, query) if donor and query else None
            if correspondence is None or correspondence[0] not in trusted[gene_id, donor]:
                proposals.append({**nomination, 'gene_id': gene_id, 'donor': donor, 'query': query,
                                  'reason': 'untrusted_rescue_donor_query_correspondence'})
                continue
            model = copy.deepcopy(raw)
            model.update(gene_id=gene_id, donor_species=donor, donor_candidate=correspondence[1],
                         coverage=metrics['coverage'], identity=metrics['identity'])
            # Reconstruct from this genome even when the rescue publication
            # already recorded successful QC. Do not accept supplied DNA.
            checked = rescue.validate_model(model, genome, catalog.get('genetic_code', 1), params)
            preserved = original_problems | (set(metrics.get('problems', [])) - {'overlap_existing_annotation'})
            if raw.get('sequence') != checked['sequence']:
                preserved.add('imported_rescue_sequence_mismatch')
            checked['problems'] = sorted(set(checked['problems']) | preserved)
            checked['rescue_revision_model_id'] = raw.get('model_id', '')
            validated.append(checked)
    return validated, proposals


def search_refinement_intervals(root, request, tmp, windows, proteins, genome, source, params, cpus):
    if not windows:
        return []
    code = source['genetic_code']
    intron = params['max_intron']
    executable = request['miniprot']['path']
    if 'fixed_donors' not in request:
        return rescue.search_intervals(tmp, windows, proteins, genome, code, intron, cpus,
                                       compact=True, executable=executable)
    key = {'schema': 1, 'genome': request['files'][source['genome']], 'code': code, 'max_intron': intron,
           'miniprot': request['miniprot']['sha256'], 'implementation': request['implementation'],
           'windows': [[list(window), [[region['id'], region['donor'], region['query'],
                         hashlib.sha256(proteins[region['donor']][region['query']].encode()).hexdigest()]
                        for region in regions]] for window, regions in windows.items()]}
    relative = Path('miniprot') / hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()

    def build_raw(work):
        models = rescue.search_intervals(work, windows, proteins, genome, code, intron, cpus,
                                         compact=True, executable=executable)
        atomic_json(work / 'raw_models.json', models)

    def guard():
        load(root)

    cache = rescue.stage(donor_cache(request), relative, key, build_raw, guard, publication_guard=guard)
    receipt_path = cache / 'receipt.json'
    receipt_sha = digest(receipt_path)
    with ReceiptSnapshot(rescue.verified) as proof:
        if not rescue.verified_dependency(proof, cache, key):
            raise ValueError('Unverified refinement predictor cache')
        present, models = proof.read_json_member(receipt_path, receipt_sha, 'raw_models.json')
        if not present:
            raise ValueError('Refinement predictor cache lacks its models')
        record = json.loads(receipt_path.read_bytes())
        proof.check({receipt_path: receipt_sha})
        for member in (rescue.INTERVAL_ARCHIVE, rescue.LOCAL_RECORDS):
            shutil.copyfile(cache / member, tmp / member)
            if digest(tmp / member) != record['files'][member]:
                raise ValueError('Refinement predictor diagnostics changed during copy')
    return models


def _nomination_indexes(targets, baseline, donors, params):
    """Find donor coding paths that warrant a new predictor search."""
    searches = []
    target_baseline = targets[baseline]
    for index, candidate in enumerate(donors):
        if not candidate['quality'].get('usable') or not candidate['protein']:
            continue
        baseline_score = None
        discrepancy = (not target_baseline['quality'].get('valid_orf')
                       or len(targets) < len(donors)
                       or abs(len(candidate['protein']) - len(target_baseline['protein']))
                       > max(3, len(candidate['protein']) * .05))
        if not discrepancy:
            baseline_score = pair_score(target_baseline, candidate)
            discrepancy = (not baseline_score.bounded and baseline_score.identity >= params['minimum_identity']
                           and (not homology_supported(min(baseline_score.coverage_a, baseline_score.coverage_b),
                                                      baseline_score.identity, params)
                                or (len(target_baseline['blocks']) > 1 and len(candidate['blocks']) > 1
                                    and baseline_score.junction_similarity < 0.8)))
        if not discrepancy:
            continue
        for known in targets:
            if not known['quality'].get('valid_orf'):
                continue
            match = (baseline_score if known is target_baseline and baseline_score is not None
                     else pair_score(known, candidate))
            if not match.bounded and homology_supported(min(match.coverage_a, match.coverage_b), match.identity, params) and match.junction_similarity >= 0.95:
                break
        else:
            searches.append(index)
    return searches


def _nomination_batch(rows, params):
    return [_nomination_indexes(*row, params) if row is not None else [] for row in rows]


def _ordered_nominations(inputs, params, cpus):
    batch_size = 64
    minimum_parallel_edges = 512
    initial = list(islice(inputs, minimum_parallel_edges)) if cpus > 1 else []
    remaining = chain(initial, inputs)
    if len(initial) < minimum_parallel_edges:
        for context, row in remaining:
            yield context, _nomination_indexes(*row, params) if row is not None else []
        return

    def consume(batch):
        contexts, future = batch
        yield from zip(contexts, future.result(), strict=True)

    # Genome/SQLite handles stay in the caller. Results are applied in edge
    # order so query insertion and proposal order remain stable.
    workers = min(cpus, 8)
    pending = deque()
    with ProcessPoolExecutor(max_workers=workers, mp_context=get_context('spawn')) as executor:
        while batch := list(islice(remaining, batch_size)):
            contexts, rows = zip(*batch, strict=True)
            pending.append((contexts, executor.submit(_nomination_batch, rows, params)))
            if len(pending) >= workers * 2:
                yield from consume(pending.popleft())
        while pending:
            yield from consume(pending.popleft())


@invocation_cached
def predict_species(root, value, name, cpus=1):
    if name not in target_species(value):
        raise ValueError('Donor-only species has no refinement prediction: ' + name)
    initial = select(root, value)
    cat_dir = catalog_species(root, value, name)
    db = catalog_index(root, value)
    corr = correspondence(root, value)
    donor_names = prediction_donor_names(corr, name, value['species'])
    params, source = parameters_for(value['request'], name), value['request']['sources'].get(name)
    dependencies = {'initial': digest(initial / 'receipt.json'), 'catalog': catalog_hashes(root, value),
                    'correspondence': digest(corr / 'receipt.json'),
                    'index': {db.parent.name: digest(db.parent / 'receipt.json')}}
    revision_directory = revision_worker(value['request'], name)
    if revision_directory is not None:
        dependencies['rescue_workers'] = {name: digest(revision_directory / 'receipt.json')}
    def build(tmp):
        if params['mode'] == 'off' or source is None:
            atomic_json(tmp / 'predictions.json', [])
            return
        catalog = json.loads((cat_dir / 'catalog_metadata.json').read_text())
        catalog['loci'] = list(iter_loci(db, name))
        edges = json.loads((corr / 'edges.json').read_text())
        edges = [edge for edge in edges if edge.get('comparison_basis') != 'protein_alignment'
                 and name in (edge['species_a'], edge['species_b'])]
        loci = {g['gene_id']: g for g in catalog['loci']}
        needed = {(edge['species_' + side], edge['gene_' + side])
                  for edge in edges if not edge['ambiguous'] for side in ('a', 'b')}
        # The receipted map contains every selection's ID and status without
        # the cohort-wide scoring and audit payloads from selection.json.
        with (initial / 'representative_map.tsv').open(newline='', encoding='utf-8') as handle:
            decisions = {(row['species'], row['gene_id']): (row['candidate_id'], row['status'])
                         for row in csv.DictReader(handle, delimiter='\t')
                         if (row['species'], row['gene_id']) in needed}
        proteins, windows, queries = {}, defaultdict(list), {}
        proposals = []
        def scoring_inputs(connection):
            for edge in edges:
                if edge['ambiguous']:
                    continue
                target_side = 'a' if edge['species_a'] == name else 'b'
                donor_side = 'b' if target_side == 'a' else 'a'
                target_gene = edge['gene_' + target_side]
                donor = edge['species_' + donor_side]
                target = loci[target_gene]
                donor_locus = load_locus(connection, donor, edge['gene_' + donor_side])
                target_candidate_id, target_status = decisions[name, target_gene]
                baseline = next(i for i, c in enumerate(target['candidates'])
                                if c['candidate_id'] == target_candidate_id)
                target_baseline = target['candidates'][baseline]
                repair_only = (target_status == 'ambiguous_correspondence'
                               and same_coordinate_source_owner(
                                   {'seqid': target['seqid'], 'strand': target['strand'],
                                    'cds': target_baseline['blocks'], 'sequence': target_baseline['cds']},
                                   target, catalog.get('genetic_code')) is not None)
                if (target.get('ambiguous_coordinates') or donor_locus.get('ambiguous_coordinates')
                        or (target_status == 'ambiguous_correspondence' and not repair_only)
                        or decisions[donor, donor_locus['gene_id']][1] == 'ambiguous_correspondence'):
                    yield ({'gene_id': target_gene, 'status': 'ambiguous_locus_ownership', 'donor': donor}, None)
                    continue
                if not any(c['cds'] and c['blocks'] for c in target['candidates']):
                    yield ({'gene_id': target_gene, 'status': 'unreconstructible_original', 'donor': donor}, None)
                    continue
                if len(target['candidates']) > params['candidate_limit'] or len(donor_locus['candidates']) > params['candidate_limit']:
                    yield ({'gene_id': target_gene, 'status': 'candidate_limit_exceeded', 'donor': donor}, None)
                    continue
                yield ((target_gene, donor, target, donor_locus['candidates'], repair_only),
                       (target['candidates'], baseline, donor_locus['candidates']))

        with store_connection(db) as connection:
            for context, searches in _ordered_nominations(scoring_inputs(connection), params, cpus):
                if isinstance(context, dict):
                    proposals.append(context)
                    continue
                if not searches:
                    continue
                target_gene, donor, target, candidates, repair_only = context
                blocks = [b for c in target['candidates'] for b in c['blocks']]
                start, end = max(0, min(b[0] for b in blocks) - params['padding']), max(b[1] for b in blocks) + params['padding']
                for index in searches:
                    candidate = candidates[index]
                    if end - start > params['max_interval']:
                        proposals.append({'gene_id': target_gene, 'status': 'interval_too_large', 'donor': donor})
                        continue
                    qid = 'q' + hashlib.sha256((target_gene + donor + candidate['candidate_id']).encode()).hexdigest()[:24]
                    region = {'id': qid, 'seqid': target['seqid'], 'start': start, 'end': end, 'query': candidate['candidate_id'],
                              'donor': donor, 'gene_id': target_gene, 'donor_candidate': candidate['candidate_id'],
                              **({'same_coordinate_repair_only': True} if repair_only else {})}
                    queries[qid] = region
                    proteins.setdefault(donor, {})[candidate['candidate_id']] = candidate['protein']
                    windows[target['seqid'], start, end].append(region)
        validated = []
        # Nominate first: a species with no search windows needs no genome I/O.
        revisions = iter(revision_models(value['request'], name))
        first_revision = next(revisions, None)
        revision_count = 0
        def counted_revisions():
            nonlocal revision_count
            for model in chain((first_revision,), revisions):
                revision_count += 1
                yield model
        native_engines = [e for e in ('annevo', 'oriongeno') if name in value['request'].get(e, {}).get('predictions', {})]
        native_enabled = bool(native_engines)
        with indexed_genome(source['genome']) if windows or first_revision is not None or native_enabled else contextlib.nullcontext() as genome:
            raw_validation = rescue.RawValidationMemo(rescue.validate_model, genome, source['genetic_code'], params,
                                                      candidate_homology=rescue.homology_problems)
            bounded = {}
            for (seqid, start, end), regions in windows.items():
                end = min(end, genome.get_reference_length(seqid))
                for region in regions:
                    region['end'] = end
                bounded[seqid, start, end] = regions
            models = search_refinement_intervals(root, value['request'], tmp, bounded, proteins,
                                                 genome, source, params, cpus)
            for model in models:
                region = queries[model['query']]
                model['seqid'] = region['seqid']
                model['cds'] = [[a + region['start'], b + region['start'], p] for a, b, p in model['cds']]
                model.update(gene_id=region['gene_id'], donor_species=region['donor'], donor_candidate=region['donor_candidate'])
                if region.get('same_coordinate_repair_only'):
                    model['same_coordinate_repair_only'] = True
                checked = raw_validation.validate(model)
                checked = rescue.complete_terminals(checked, genome, source['genetic_code'], params,
                    proteins[region['donor']][region['donor_candidate']], rescue.validate_model)
                if any(a < region['start'] or b > region['end'] for a, b, _ in checked['cds']):
                    checked['problems'].append('outside_search_window')
                validated.append(checked)
            for engine in native_engines:
                import annevo_gene_models
                native_regions = [{**region, 'query': region['donor_candidate']} for region in queries.values()]
                native_models, native_audit = annevo_gene_models.nominate(
                    value['request'][engine], name, native_regions, proteins, genome, params, cpus=cpus, engine=engine)
                native_checked = []
                for model in native_models:
                    region = queries[model['query']]
                    model.update(gene_id=region['gene_id'], donor_species=region['donor'], donor_candidate=region['donor_candidate'])
                    if region.get('same_coordinate_repair_only'):
                        model['same_coordinate_repair_only'] = True
                    checked = raw_validation.validate(model)
                    checked['problems'] = sorted(set(checked['problems']) | set(model[engine + '_problems']))
                    native_checked.append(checked)
                annevo_gene_models.confirm_transferred_structure(native_checked, validated)
                validated.extend(native_checked)
                atomic_json(tmp / (engine + '_candidates.json'), native_audit)
            if first_revision is not None:
                with store_connection(db) as connection:
                    imported, import_proposals = import_rescue_revisions(
                        counted_revisions(), catalog, edges, params, genome,
                        rescue_revision_donor_resolver(value['request'], connection))
                validated.extend(imported)
                proposals.extend(import_proposals)
        rna = read_table(value['request']['rna']) if value['request']['rna'] else []
        rna.extend(refinement_rna.evidence_rows(name, source['genome'],
                                               value['request'].get('rna_bams', {}).get(name, []), validated, tmp))
        if validated:
            catalog['annotation_spans'] = annotation_ownership_spans(source['gff'], catalog)
        predictions = classify_predictions(validated, catalog, edges, params, rna, value['request']['files'][source['genome']])
        atomic_json(tmp / 'predictions.json', predictions)
        atomic_json(tmp / 'nominations.json', {'queries': list(queries.values()), 'proposals': proposals})
        atomic_json(tmp / 'summary.json', {'queries': len(queries), 'windows': len(windows), 'alignments': len(validated),
                                         'rescue_revision_models': revision_count,
                                         'accepted': sum(r['status'] == 'accepted' for r in predictions),
                                         'proposals': sum(r['status'] != 'accepted' for r in predictions) + len(proposals)})
    return stage(root, Path('predictions') / name, {'dependencies': dependencies}, build, donor_names)


def source_rescue_locus_provenance(gene):
    """Preserve historical locus evidence without borrowing old path support."""
    fields = {key for _, key in SOURCE_RESCUE_LOCUS_ATTRIBUTES}
    result = {}
    records = [gene.get('source_rescue_locus_provenance', {})]
    for candidate in gene.get('candidates', ()):
        records.extend((candidate.get('source_rescue_locus_provenance', {}), candidate.get('rescue_path_selection', {})))
    for record in records:
        for key in fields & record.keys():
            value = str(record[key])
            if key in result and result[key] != value:
                raise ValueError('Conflicting historical rescue locus provenance: ' + gene['gene_id'])
            result[key] = value
    return result


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
    history = source_rescue_locus_provenance(gene)
    historical_attributes = ''.join(f';{attribute}={attr(history[key])}' for attribute, key in SOURCE_RESCUE_LOCUS_ATTRIBUTES
                                    if key in history)
    path_attributes = ''
    donors = candidate.get('support', {}).get('donors')
    if donors is not None:
        if not isinstance(donors, list) or any(not isinstance(donor, str) or not donor for donor in donors):
            raise ValueError('Invalid current coding-path donor support: ' + candidate['candidate_id'])
        path_attributes = ';path_independent_donor_species_count=' + str(len(set(donors) - {name}))
    import annevo_gene_models
    path_attributes += annevo_gene_models.gff_attributes(
        [{engine: origin} for engine in ('annevo', 'oriongeno') for origin in candidate.get('support', {}).get(engine, [])])
    lines = [row('gene', start, end, '.', 'ID=' + gid + historical_attributes),
             row('mRNA', start, end, '.', f'ID={tid};Parent={gid};gene_id={token}' + path_attributes)]
    for i, (a, b, phase) in enumerate(blocks, 1):
        lines += [row('exon', a, b, '.', f'ID={tid}.exon{i};Parent={tid}'),
                  row('CDS', a, b, phase, f'ID={tid}.cds{i};Parent={tid};gene_id={token}')]
    return ''.join(lines)


def _gff_source_graph(original, *, retain_attributes=True):
    """Parse source relationships, retaining decoded attributes only for callers that need them."""
    records, parents, implicit = [], {}, {}
    lines = original.splitlines(keepends=True) if isinstance(original, str) else original
    for line in lines:
        if line.strip() == '##FASTA':
            break
        if line.startswith('#') or not line.strip():
            attributes, parent_ids = ({}, []) if retain_attributes else (None, ())
            records.append((line, None, attributes, parent_ids))
            continue
        fields = line.rstrip().split('\t')
        if len(fields) != 9:
            raise ValueError('Invalid source GFF row')
        is_gtf = bool(re.match(r'^\s*[^\s=;]+\s+', fields[8]))
        parsed = parse_gff_attributes(fields[8])
        attributes = {key: values[0] for key, values in parsed.items() if values}
        parent_ids = parsed.get('Parent', ())
        tid, gid = attributes.get('transcript_id'), attributes.get('gene_id')
        # GTF encodes its graph on every child rather than requiring ID/Parent rows.
        if tid and not attributes.get('ID') and fields[2].lower().endswith(('rna', 'transcript')):
            attributes['ID'] = tid
        elif gid and not tid and fields[2].lower().endswith('gene') and not attributes.get('ID'):
            attributes['ID'] = gid
        if not parent_ids and tid:
            parent_ids = [gid] if attributes.get('ID') == tid else [tid]
            parent_ids = tuple(filter(None, parent_ids))
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
        records.append((line, fields, attributes if retain_attributes else identifier,
                        list(parent_ids) if retain_attributes else parent_ids))
    return records, parents, implicit


def _selected_gff_source_rows(source_graph, transcript_ids, *, retain_all=False, synthetic_parents_only=False):
    """Render a selected or complete view without mutating the parsed graph."""
    if synthetic_parents_only and not retain_all:
        raise ValueError('Synthetic parent export requires the complete source graph')
    records, parents, implicit = source_graph

    def record_id(value):
        return value.get('ID') if isinstance(value, dict) else value

    if retain_all:
        transcript_ids = set(transcript_ids) | set(implicit)
    ancestors, pending = set(), list(transcript_ids)
    while pending and not retain_all:
        identifier = pending.pop()
        for parent in parents.get(identifier, set()) - ancestors - {''}:
            ancestors.add(parent)
            pending.append(parent)
    coding_owners = {parent for _, fields, _, pids in records if fields is not None and fields[2] == 'CDS' for parent in pids}
    excluded_coding_owners = coding_owners - transcript_ids
    selected_ancestors = transcript_ids | ancestors
    descendants = set(transcript_ids)
    changed = True
    while changed and not retain_all:
        changed = False
        for _, fields, identifier, parent_ids in records:
            identifier = record_id(identifier)
            if fields is None or identifier in ancestors or identifier in excluded_coding_owners:
                continue
            if descendants.intersection(parent_ids):
                if identifier and identifier not in descendants:
                    descendants.add(identifier)
                    changed = True
    structural_parents = descendants | ancestors
    lines, known_ids = [], {identifier for _, fields, value, _ in records if fields is not None
                           and (identifier := record_id(value)) is not None}
    gene_spans = defaultdict(list)
    for tid in sorted(transcript_ids & implicit.keys()):
        gid, seqid, strand, spans = implicit[tid]
        if gid in known_ids and tid in known_ids:
            continue
        start, end = min(a for a, _ in spans), max(b for _, b in spans)
        if gid not in known_ids:
            gene_spans[(gid, seqid, strand)].append((start, end))
        if tid not in known_ids:
            kind = 'mRNA' if tid in coding_owners else 'transcript'
            lines.append(f'{seqid}\tGeneGalleon\t{kind}\t{start}\t{end}\t.\t{strand}\t.\tID={quote(tid, safe="._-:")};Parent={quote(gid, safe="._-:")};gene_id={quote(gid, safe="._-:")}\n')
    for _, fields, identifier, parent_ids in records:
        if (fields is not None and parent_ids and fields[2].lower().endswith(('rna', 'transcript'))
                and (retain_all or record_id(identifier) in selected_ancestors)):
            span = int(fields[3]), int(fields[4])
            for gid in parent_ids:
                if gid not in known_ids:
                    gene_spans[(gid, fields[0], fields[6])].append(span)
    synthesized_genes = []
    for (gid, seqid, strand), spans in sorted(gene_spans.items()):
        start, end = min(a for a, _ in spans), max(b for _, b in spans)
        synthesized_genes.append(f'{seqid}\tGeneGalleon\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={quote(gid, safe="._-:")}\n')
    lines = ['##gff-version 3\n'] + synthesized_genes + lines
    if synthetic_parents_only:
        return ''.join(lines)
    for line, fields, identifier, parent_ids in records:
        if fields is None:
            if not line.startswith('##gff-version'):
                lines.append(line)
            continue
        identifier = record_id(identifier) or ''
        structural = identifier in transcript_ids or identifier in ancestors
        keep = retain_all or structural or identifier in descendants or (identifier not in coding_owners and bool(descendants.intersection(parent_ids)))
        if keep:
            if parent_ids and not retain_all:
                retained = [parent for parent in parent_ids if parent in (structural_parents if structural else descendants)]
                if not retained:
                    continue
                fields = fields[:]
                parts = fields[8].split(';')
                fields[8] = ';'.join('Parent=' + ','.join(quote(parent, safe='._-:') for parent in retained) if part.startswith('Parent=') else part for part in parts)
                line = '\t'.join(fields) + '\n'
            lines.append(line)
    return ''.join(lines)


def selected_gff_rows(original, transcript_ids, *, retain_all=False):
    """Canonicalize source relationships, retaining selected or all features."""
    return _selected_gff_source_rows(_gff_source_graph(original), transcript_ids, retain_all=retain_all)


def extend_gene_bounds(original, bounds):
    if not bounds:
        return original, []
    lines, changes = [], []
    for line in original.splitlines(keepends=True):
        if line.startswith('#') or not line.strip():
            lines.append(line)
            continue
        fields = line.rstrip().split('\t')
        if len(fields) == 9 and fields[2].lower().endswith('gene'):
            parsed = parse_gff_attributes(fields[8], keys={'ID', 'gene_id'})
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
    # Catalog-admitted ordinary DNA uses the same translation kernel as
    # catalog construction. Preserve Biopython's handling of other spellings
    # and code overrides for direct callers.
    if code == 1 and CONCRETE_CDS.fullmatch(coding):
        translated = translate_sequence_string(coding, code, False)
    else:
        translated = str(Seq(coding).translate(table=code)) if coding else ''
    if translated.removesuffix('*') != expected:
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


def source_phase_repair(candidate, species):
    """A complete, source-matching CDS fixes conflicting numeric phases uniquely."""
    quality = candidate['quality']
    blocks = candidate['blocks']
    sequence = candidate['cds']
    if (not quality.get('phase_conflict') or not quality.get('has_start') or not quality.get('has_stop')
            or quality.get('partial') or not blocks or blocks[0][2] != 0
            or any(block[2] not in (0, 1, 2) for block in blocks)
            or not any(record['cds'].upper() == sequence for record in candidate.get('source_cds', []))
            or rescue.phase_conventions([{'start': a, 'end': b, 'phase': str(p)} for a, b, p in blocks])):
        return None
    fixed = rescue.coding_phases(blocks)
    shape = candidate['seqid'], candidate['strand'], fixed
    identifier = species + '_ggphase_' + hashlib.sha256(repr((candidate['source_gene_id'], shape)).encode()).hexdigest()[:16]
    return {'candidate_id': identifier, 'kind': 'phase_only', 'adoption': 'automatic',
            'basis': 'exact_source_cds', 'seqid': candidate['seqid'], 'strand': candidate['strand'],
            'cds': fixed, 'cds_sha256': hashlib.sha256(sequence.encode()).hexdigest(), 'support': []}


def automatic_coding_gain(original, checked, *, syntenic):
    """Keep the intact source CDS or extend a partial model in the same frame."""
    old, new = original['cds'], checked['sequence']
    if len(new) - len(old) < 3 * max(15, .10 * len(old) / 3):
        return False
    _, in_frame = coding_overlap({**original, 'cds': original['blocks']}, checked)
    quality = original['quality']
    if quality['valid_orf'] or not syntenic:
        return in_frame == len(old) and (new.endswith(old) or new.startswith(old))
    blocks, repaired = original['blocks'], checked['cds']
    start_side, stop_side = (0, 1) if original['strand'] == '+' else (1, 0)
    same_start = blocks[0][start_side] == repaired[0][start_side]
    same_stop = blocks[-1][stop_side] == repaired[-1][stop_side]
    return (quality['partial'] and in_frame * 3 >= len(old) * 2
            and ((quality['has_start'] and same_start) or (quality['has_stop'] and same_stop)))


def repair_coordinates(review):
    blocks = review['cds']
    return {'seqid': review['seqid'], 'strand': review['strand'], 'cds': blocks,
            'start': min(b[0] for b in blocks), 'end': max(b[1] for b in blocks)}


def automatic_model_choices(options, originals, ownership, genome, code, params):
    """Choose a unique qualifying path and keep competing qualifying paths unresolved."""
    choices = []
    for owner, rows in options.items():
        original, chosen = originals[owner]
        supported = []
        for review in rows:
            coordinates = repair_coordinates(review)
            geometry, exact = rescue.model_repair_geometry(coordinates, ownership)
            if [row['gene_id'] for row in geometry] != [owner]:
                continue
            basis = review.get('basis', 'exact_source_cds')
            if review['kind'] == 'phase_only':
                if not exact:
                    continue
                # This path has exact source evidence rather than alignment metrics.
                metrics = {'coverage': 0., 'identity': 0.}
                checks = {**params, 'minimum_coverage': 0., 'minimum_identity': 0., 'minimum_broad_coverage': 0.}
                support = []
            else:
                if exact or geometry != review['owner_overlap']:
                    continue
                support = [row for row in review['support']
                           if owner in (row.get('syntenic_locus'), row.get('reciprocal_locus'))
                           and row['donor'] != original['species']
                           and homology_supported(row['coverage'], row['identity'], params)]
                if len({row['donor'] for row in support}) < max(2, params['min_support']):
                    continue
                syntenic = len({row['donor'] for row in support if row.get('syntenic_locus') == owner}) >= max(2, params['min_support'])
                basis = 'syntenic_cds_gain' if syntenic else 'reciprocal_cds_extension'
                metrics, checks = support[0], params
            checked = rescue.validate_model({**coordinates, 'coverage': metrics['coverage'],
                                             'identity': metrics['identity'], 'frameshift': False},
                                            genome, code, checks)
            if (checked['problems'] or checked['cds'] != review['cds']
                    or hashlib.sha256(checked['sequence'].encode()).hexdigest() != review['cds_sha256']):
                continue
            if review['kind'] == 'coding_path' and not automatic_coding_gain(original, checked, syntenic=syntenic):
                continue
            supported.append(({**review, 'owner_overlap': geometry, 'support': support,
                               'adoption': 'automatic', 'basis': basis}, checked))
        if len(supported) != 1:
            continue
        review, checked = supported[0]
        if chosen['candidate_id'] == original['candidate_id'] or (
                chosen['cds'] == checked['sequence'] and chosen['blocks'] == review['cds']
                and chosen['seqid'] == review['seqid'] and chosen['strand'] == review['strand']):
            choices.append((review, checked))
    return choices


def automatic_split_choices(options, originals, ownership, genome, code, params, species):
    """Merge incomplete fragments only when independent donors agree on one CDS."""
    choices = []
    for review in options:
        members = {row['gene_id'] for row in review['owner_overlap']}
        if not members <= originals.keys():
            continue
        coordinates = repair_coordinates(review)
        if rescue.split_model_geometry(coordinates, ownership) != review['owner_overlap']:
            continue
        support = [row for row in review['support'] if row['donor'] != species
                   and homology_supported(row['coverage'], row['identity'], params)]
        if len({row['donor'] for row in support}) < max(2, params['min_support']):
            continue
        checked = rescue.validate_model({**coordinates, 'coverage': support[0]['coverage'],
                                         'identity': support[0]['identity'], 'frameshift': False},
                                        genome, code, params)
        if (checked['problems'] or checked['cds'] != review['cds']
                or hashlib.sha256(checked['sequence'].encode()).hexdigest() != review['cds_sha256']):
            continue
        for member in sorted(members):
            original = originals[member]
            _, in_frame = coding_overlap({**original, 'cds': original['blocks']}, checked)
            if in_frame * 3 < 2 * len(original['cds']):
                break
        else:
            choices.append(({**review, 'support': support, 'adoption': 'automatic',
                             'basis': 'independent_donors_complementary_fragments'}, checked))
    return choices


def independent_source_phases(source_gff, options, originals):
    """Return phase repairs whose CDS rows can be edited without affecting other transcripts."""
    expected = {}
    for owner, rows in options.items():
        if rows[0]['kind'] == 'phase_only':
            model = originals[owner][0]
            expected[owner] = {(model['seqid'], model['strand'], a, b, model['source_transcript_id'])
                               for a, b, _ in model['blocks']}
    if not expected:
        return set()
    pending = set().union(*expected.values())
    positions = {key[:4] for key in pending}
    seen, ambiguous = set(), set()
    with open_text(Path(source_gff)) as handle:
        for line in handle:
            if line.strip() == '##FASTA':
                break
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip().split('\t')
            if len(fields) != 9 or fields[2] != 'CDS':
                continue
            locus = fields[0], fields[6], int(fields[3]) - 1, int(fields[4])
            if locus not in positions:
                continue
            attrs = parse_gff_attributes(fields[8], keys={'Parent', 'transcript_id'})
            parents = attrs.get('Parent') or attrs.get('transcript_id', ())
            matches = {(*locus, parent) for parent in parents} & pending
            if len(parents) != 1:
                ambiguous.update(matches)
            ambiguous.update(matches & seen)
            seen.update(matches)
    return {owner for owner, keys in expected.items() if keys <= seen and not keys & ambiguous}


def prepare_reviewed_models(request, name, db, selected):
    """Recheck approved merge and repair paths against the frozen target inputs."""
    merges = request.get('split_model_merges', {}).get(name, [])
    repairs = request.get('model_repairs', {}).get(name, [])
    source = request['sources'][name]
    parameters = parameters_for(request, name)
    targets = {row['gene_id'] for review in chain(merges, repairs) for row in review['owner_overlap']}
    nominations = defaultdict(list)
    for row in request.get('automatic_model_repairs', {}).get(name, []):
        nominations[row['owner_overlap'][0]['gene_id']].append(row)
    split_options = request.get('automatic_split_merges', {}).get(name, [])
    split_targets = {row['gene_id'] for review in split_options for row in review['owner_overlap']}
    owners, auto_originals, auto_options, split_originals = {}, {}, {}, {}
    seen_automatic, duplicate_automatic = set(), set()
    for gene in iter_loci(db, name):
        originals = [c for c in gene['candidates'] if c['origin'] == 'original']
        ids = {row['source_gene_id'] for row in originals}
        duplicate_automatic.update(ids & seen_automatic)
        seen_automatic.update(ids)
        if ids & targets:
            if len(ids) != 1 or gene.get('ambiguous_coordinates') or ids & owners.keys():
                raise ValueError('Ambiguous catalog ownership for reviewed model: ' + gene['gene_id'])
            owners[next(iter(ids))] = gene
            continue
        if parameters['mode'] != 'conservative' or gene.get('ambiguous_coordinates') or len(originals) != 1:
            continue
        original = originals[0]
        owner = original['source_gene_id']
        if any(original['quality'].get(flag) for flag in (
                'annotated_exception', 'structure_problem', 'internal_stop', 'ambiguous',
                'invalid_base', 'translation_uncertain', 'phase_unresolved', 'sequence_mismatch')):
            continue
        chosen_id = selected[gene['gene_id']]['candidate_id']
        if (owner in split_targets and chosen_id == original['candidate_id']
                and original['quality']['partial'] and original['quality']['usable']):
            split_originals[owner] = original
            owners[owner] = gene
        phase = source_phase_repair(original, name)
        options = [phase] if phase else nominations.get(owner, []) if original['quality']['usable'] else []
        if options:
            if owner in auto_originals:
                auto_options.pop(owner, None)
                continue
            chosen = next(c for c in gene['candidates'] if c['candidate_id'] == chosen_id)
            auto_originals[owner], auto_options[owner], owners[owner] = (original, chosen), options, gene
    manual_owners = {owner: owners[owner] for owner in targets if owner in owners}
    if targets != manual_owners.keys() or len({gene['gene_id'] for gene in manual_owners.values()}) != len(targets):
        raise ValueError('Approved model owners do not match distinct catalog loci: ' + name)
    for owner in duplicate_automatic:
        auto_options.pop(owner, None)
        split_originals.pop(owner, None)
    independent = independent_source_phases(source['gff'], auto_options, auto_originals)
    auto_options = {owner: rows for owner, rows in auto_options.items()
                    if rows[0]['kind'] != 'phase_only' or owner in independent}
    if not (targets or auto_options or split_originals):
        return [], {}, [], []
    with open_text(Path(source['gff'])) as handle:
        ownership = rescue.OwnershipIndex(rescue.original_ownership(handle, rescue.attributes))
    code = source['genetic_code']
    new_genes, repaired, merge_audit, repair_audit, footprints = [], {}, [], [], []
    with indexed_genome(source['genome']) as genome:
        automatic = automatic_model_choices(auto_options, auto_originals, ownership, genome, code, parameters)
        automatic.extend(automatic_split_choices(split_options, split_originals, ownership, genome, code, parameters, name))
        automatic_paths = rescue.OwnershipIndex([{**repair_coordinates(row), 'gene_id': row['candidate_id']}
                                                 for row, _ in automatic])
        manual_paths = rescue.OwnershipIndex([{**repair_coordinates(row), 'gene_id': row['candidate_id']}
                                              for row in chain(merges, repairs)])
        for review, checked in chain(((row, None) for row in chain(merges, repairs)), automatic):
            identifier = review['candidate_id']
            kind = review.get('kind')
            blocks = review['cds']
            coordinates = repair_coordinates(review)
            if checked is not None:
                if (automatic_paths.overlapping(coordinates, owner=identifier)
                        or manual_paths.overlapping(coordinates)):
                    continue
                geometry, support = review['owner_overlap'], review['support']
            else:
                if kind:
                    geometry, exact = rescue.model_repair_geometry(coordinates, ownership)
                    if kind != ('phase_only' if exact else 'coding_path'):
                        raise ValueError('Approved model repair type changed: ' + identifier)
                else:
                    geometry = rescue.split_model_geometry(coordinates, ownership)
                if not geometry or geometry != review['owner_overlap']:
                    raise ValueError('Approved model ownership changed: ' + identifier)
                support = [row for row in review['support']
                           if homology_supported(row['coverage'], row['identity'], parameters)]
                if not support or (kind and len({r['donor'] for r in support} - {name}) < (
                        1 if kind == 'phase_only' else max(2, parameters['min_support']))):
                    raise ValueError('Approved model lacks required independent homology: ' + identifier)
                checked = rescue.validate_model({**coordinates, 'coverage': support[0]['coverage'],
                                                 'identity': support[0]['identity'], 'frameshift': False},
                                                genome, code, parameters)
            sequence = checked['sequence']
            donors = sorted({row['donor'] for row in support})
            if (checked['problems'] or checked['cds'] != blocks
                    or hashlib.sha256(sequence.encode()).hexdigest() != review['cds_sha256']):
                raise ValueError('Approved model no longer has a complete validated CDS: ' + identifier)

            source_id, gene_id, gene_token, transcript = identifier, identifier, identifier, identifier + '.t1'
            audit = {**review, 'species': name}
            if kind:
                source_id = geometry[0]['gene_id']
                gene = owners[source_id]
                gene_id = gene['gene_id']
                gene_token = gene.get('gene_token', gene_id)
                originals = [c for c in gene['candidates'] if c['origin'] == 'original']
                if any(c['quality'].get('annotated_exception') for c in originals):
                    raise ValueError('Approved model repair targets a protected annotation: ' + identifier)
                if kind == 'phase_only':
                    matches = [c for c in originals if
                               [b[:2] for b in c['blocks']] == [b[:2] for b in blocks]
                               and c['quality'].get('phase_conflict') and c['cds'] == sequence
                               and any(row['cds'].upper() == sequence for row in c.get('source_cds', []))]
                    if len(matches) != 1:
                        raise ValueError('Phase-only repair requires one conflicting source model with identical CDS: ' + identifier)
                    transcript = matches[0]['source_transcript_id']
                audit['previous_candidate_id'] = selected[gene_id]['candidate_id']
            elif identifier in selected:
                raise ValueError('Approved split model ID collides with an existing catalog locus: ' + identifier)
            automatic_adoption = review.get('adoption') == 'automatic'
            origin = ('automatic_model_repair' if automatic_adoption else 'reviewed_model_repair') if kind else 'split_model_merge'
            candidate = {'candidate_id': identifier, 'gene_id': gene_id, 'species': name,
                         'source_gene_id': source_id, 'source_transcript_id': transcript, 'gene_token': gene_token,
                         'seqid': review['seqid'], 'strand': review['strand'], 'blocks': blocks,
                         'cds': sequence, 'origin': origin, 'support': {'donors': donors}}
            candidate['quality'] = validate_candidate(candidate, code)
            if not candidate['quality']['valid_orf'] or not candidate['protein']:
                if review.get('adoption') == 'automatic':
                    continue
                raise ValueError('Approved model cannot be admitted for translation: ' + identifier)
            audit.update(gene_id=gene_id, source_transcript_id=transcript)
            if kind:
                candidate['repair_kind'] = kind
                repaired[gene_id] = candidate
                repair_audit.append(audit)
            else:
                new_genes.append({'gene_id': identifier, 'seqid': review['seqid'], 'strand': review['strand'],
                                  'candidates': [candidate]})
                audit['replaced_gene_ids'] = [owners[row['gene_id']]['gene_id'] for row in geometry]
                merge_audit.append(audit)
            selected[gene_id] = {'species': name, 'gene_id': gene_id, 'candidate_id': identifier,
                                 'source_transcript_id': transcript,
                                 'status': 'repaired_model' if kind else 'merged_split_model',
                                 'reason': ('automatic_model_repair' if kind else 'automatic_split_model_merge')
                                 if automatic_adoption else
                                 ('reviewed_rescue_model_repair' if kind else 'reviewed_rescue_split_model')}
            footprints.append({**coordinates, 'gene_id': identifier})
    approved_paths = rescue.OwnershipIndex(footprints)
    for path in footprints:
        if approved_paths.overlapping(path, owner=path['gene_id']):
            raise ValueError('Approved model coding paths overlap: ' + path['gene_id'])
    return new_genes, repaired, merge_audit, repair_audit


def apply_reviewed_phase_rows(graph, repairs):
    """Change only the selected transcript's CDS phase fields in derived views."""
    expected = {(c['source_transcript_id'], c['seqid'], c['strand'], a, b): phase
                for c in repairs if c['repair_kind'] == 'phase_only' for a, b, phase in c['blocks']}
    if not expected:
        return
    seen = set()
    for i, (_, fields, identifier, parents) in enumerate(graph[0]):
        if fields is None or fields[2] != 'CDS':
            continue
        keys = [(p, fields[0], fields[6], int(fields[3]) - 1, int(fields[4])) for p in parents]
        matches = expected.keys() & set(keys)
        if not matches:
            continue
        if len(matches) != 1 or len(parents) != 1 or matches & seen:
            raise ValueError('Cannot isolate reviewed phase correction to one transcript: ' + str(matches))
        seen.update(matches)
        fields = fields[:]
        fields[7] = str(expected[next(iter(matches))])
        graph[0][i] = ('\t'.join(fields) + '\n', fields, identifier, parents)
    if seen != expected.keys():
        raise ValueError('Reviewed phase correction did not match source CDS rows: '
                         + str(sorted(expected.keys() - seen)))


def write_sequence_view(tmp, db, request, selected_by_species, retired_by_species, coding_audit):
    """Publish paired coding sequences for the complete mixed cohort."""
    output = tmp / 'sequence'
    for directory in ('species_cds', 'species_protein'):
        (output / directory).mkdir(parents=True)
    genome_admission = {(name, gene): status for name, gene, _, status, *_ in coding_audit}
    transcriptome_sources = request['transcriptome_sources']
    inputs, decisions = [], []
    for name in sorted(selected_by_species):
        rows = selected_by_species[name]
        cds_path, pep_path = (output / role / (name + '.fa') for role in ('species_cds', 'species_protein'))
        if name in request['sources']:
            code = request['sources'][name]['genetic_code']
            os.link(tmp / 'analysis_cds' / cds_path.name, cds_path)
            os.link(tmp / 'species_protein' / pep_path.name, pep_path)
            for row in rows:
                gene = row['gene_id']
                if gene not in retired_by_species.get(name, ()):
                    decisions.append((row, 'genome_annotation', '', genome_admission[name, gene]))
            kind = 'genome_annotation'
        else:
            code = transcriptome_sources[name]['genetic_code']
            selected = {row['gene_id']: row for row in rows}
            with cds_path.open('w') as cds, pep_path.open('w') as protein:
                for locus in iter_loci(db, name):
                    gene = locus['gene_id']
                    row = selected[gene]
                    records = locus['candidates']
                    candidate = next(c for c in records if c['candidate_id'] == row['candidate_id'])
                    baseline = _baseline(locus, records)['source_transcript_id']
                    usable = candidate['quality']['usable']
                    if usable:
                        sequence = candidate['cds']
                        cds.write(f'>{gene}\n{sequence[:len(sequence) // 3 * 3]}\n')
                        protein.write(f'>{gene}\n{candidate["protein"]}\n')
                    decisions.append((row, 'transcriptome_assembly', baseline, 'included' if usable else 'excluded'))
            kind = 'transcriptome_assembly'
        inputs.append((name, kind, cds_path.relative_to(output), pep_path.relative_to(output), code))
    rescue.write_tsv(output / 'inputs.tsv', ('species', 'source_kind', 'cds', 'protein', 'genetic_code'), inputs)
    rescue.write_tsv(output / 'species_genetic_code.tsv', ('species', 'genetic_code'),
                     [(name, code) for name, _, _, _, code in inputs])
    rescue.write_tsv(output / 'representative_map.tsv', (*MAP_FIELDS, 'source_kind', 'baseline_transcript_id', 'sequence_admission'),
                     [[*(row.get(k, '') for k in MAP_FIELDS), kind, baseline, admission]
                      for row, kind, baseline, admission in decisions])
    assembly_rows = [row for row, kind, _, _ in decisions if kind == 'transcriptome_assembly']
    return {'transcriptome_species': len(transcriptome_sources), 'transcriptome_loci': len(assembly_rows),
            'transcriptome_changed_representatives': sum(row['status'] == 'conserved' for row in assembly_rows)}


def finalize(root, value):
    final = select(root, value, predictions=True)
    db = catalog_index(root, value, predictions=True)
    dependencies = {'selection': digest(final / 'receipt.json'),
                    'catalog': catalog_hashes(root, value),
                    'index': {db.parent.name: digest(db.parent / 'receipt.json')}}
    def build(tmp):
        targets = set(target_species(value))
        selections = [row for row in json.loads((final / 'selection.json').read_text())['selections']
                      if row['species'] in targets]
        selected_by_species = defaultdict(list)
        for selected_row in selections:
            selected_by_species[selected_row['species']].append(selected_row)
        merged_by_species, repaired_by_species, retired_by_species = {}, {}, {}
        merges, repairs = [], []
        request = value['request']
        for name in sorted(targets & request['sources'].keys()):
            if parameters_for(request, name)['mode'] != 'conservative':
                continue
            selected = {row['gene_id']: row for row in selected_by_species[name]}
            new_genes, repaired, merge_audit, repair_audit = prepare_reviewed_models(request, name, db, selected)
            selected_by_species[name] = list(selected.values())
            merged_by_species[name], repaired_by_species[name] = new_genes, repaired
            retired_by_species[name] = {gene_id for row in merge_audit for gene_id in row['replaced_gene_ids']}
            merges.extend(merge_audit)
            repairs.extend(repair_audit)
        for role in ('species_cds', 'species_protein', 'species_gff', 'species_genome', 'analysis_cds', 'analysis_gff', 'full_annotation', 'source_annotation', 'source_cds', 'all_candidates', 'species_transcript', 'utr_audit'):
            (tmp / role).mkdir()
        map_path = tmp / 'representative_map.tsv'
        if merges or repairs or request.get('transcriptome_sources') or 'target_species' in request:
            rescue.write_tsv(map_path, MAP_FIELDS, [[row.get(k, '') for k in MAP_FIELDS]
                             for name in sorted(targets & request['sources'].keys()) for row in selected_by_species[name]
                             if row['gene_id'] not in retired_by_species.get(name, ())])
        else:
            shutil.copyfile(final / 'representative_map.tsv', map_path)
        if merges:
            atomic_json(tmp / 'split_model_merges.json', merges)
        if repairs:
            atomic_json(tmp / 'model_repairs.json', repairs)
        representative_map_proof = FreshDigestBatch()
        rows, changes, translation_audit, coding_audit, exclusions, gene_bounds_audit = [], [], [], [], [], []
        unpublished_candidates, utr_summaries = [], []
        for name in sorted(targets & request['sources'].keys()):
            catalog = json.loads((root / 'catalog' / name / 'catalog_metadata.json').read_text())
            catalog['loci'] = iter_loci(db, name)
            source = value['request']['sources'][name]
            selected = {r['gene_id']: r for r in selected_by_species[name]}
            paths = {'cds': tmp / 'species_cds' / (name + '.fa'), 'protein': tmp / 'species_protein' / (name + '.fa'),
                     'gff': tmp / 'species_gff' / (name + '.gff3'), 'genome': tmp / 'species_genome' / (name + '.fa'),
                     'analysis_cds': tmp / 'analysis_cds' / (name + '.fa'), 'analysis_gff': tmp / 'analysis_gff' / (name + '.gff3'),
                     'transcript': tmp / 'species_transcript' / (name + '.fa')}
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
            source_graph = _gff_source_graph(original, retain_attributes=False)
            retired = retired_by_species.get(name, set())
            new_genes = merged_by_species.get(name, [])
            repaired = repaired_by_species.get(name, {})
            apply_reviewed_phase_rows(source_graph, repaired.values())
            new_transcripts = [gene['candidates'][0] for gene in new_genes]
            new_transcripts.extend(c for c in repaired.values() if c['repair_kind'] == 'coding_path')
            if new_transcripts:
                occupied = source_graph[1].keys() | source_graph[2].keys()
                occupied.update(row['source_transcript_id'] for row in selected_by_species[name]
                                if row['status'] not in {'merged_split_model', 'repaired_model'})
                for candidate in new_transcripts:
                    tid = candidate['source_transcript_id']
                    identifiers = {tid, *(f'{tid}.{kind}{i}' for kind in ('exon', 'cds')
                                   for i in range(1, len(candidate['blocks']) + 1))}
                    if candidate['origin'] == 'split_model_merge':
                        identifiers.add(candidate['gene_id'])
                    if identifiers & occupied:
                        raise ValueError('Approved model GFF ID collides: ' + candidate['candidate_id'])
                    occupied.update(identifiers)
            directive = re.search(r'^##FASTA(?:\r?\n|$)', original, flags=re.MULTILINE)
            sep, embedded = ('##FASTA', original[directive.start() + len('##FASTA'):]) if directive else ('', '')
            additions, predicted_rows, original_transcripts, expanded_bounds, analysis_candidates = [], [], set(), {}, []
            with paths['cds'].open('w') as cds, paths['protein'].open('w') as pep, paths['gff'].open('w') as gff, \
                    paths['analysis_cds'].open('w') as analysis_cds, \
                    (tmp / 'all_candidates' / (name + '.fa')).open('w') as all_cds:
                for gene in chain(catalog['loci'], new_genes):
                    chosen = selected[gene['gene_id']]
                    candidates = {c['candidate_id']: c for c in gene['candidates']}
                    if gene['gene_id'] in repaired:
                        repair = repaired[gene['gene_id']]
                        if repair['candidate_id'] in candidates:
                            raise ValueError('Approved model repair ID collides: ' + repair['candidate_id'])
                        candidates[repair['candidate_id']] = repair
                    candidate = candidates[chosen['candidate_id']]
                    for c in candidates.values():
                        if (c['origin'] == 'predicted' and not c.get('rescue_alternative_coding_path')
                                and not c.get('same_coordinate_source_repair') and not c['quality'].get('rna_supported')
                                and chosen['status'] == 'ambiguous_correspondence'):
                            # Keep copy-unresolved homology in the prediction audit.
                            unpublished_candidates.append({'species': name, 'gene_id': gene['gene_id'],
                                'candidate_id': c['candidate_id'], 'reason': chosen['reason']})
                            continue
                        all_cds.write(f">{c['candidate_id']}\n{c['cds']}\n")
                        if c['origin'] == 'split_model_merge':
                            additions.append(candidate_gff(gene, c, name))
                        elif (c.get('repair_kind') == 'coding_path' or
                              (c['origin'] == 'predicted' and not c.get('rescue_alternative_coding_path')
                               and not c.get('same_coordinate_source_repair'))):
                            # Preserve the original gene and attach only new transcript/CDS rows.
                            lines = candidate_gff(gene, c, name, gene_id=c['source_gene_id'],
                                                  gene_token=c.get('gene_token', gene.get('gene_token', gene['gene_id']))).splitlines(keepends=True)[1:]
                            start, end = min(b[0] for b in c['blocks']), max(b[1] for b in c['blocks'])
                            old = expanded_bounds.get(c['source_gene_id'], (start, end))
                            expanded_bounds[c['source_gene_id']] = min(old[0], start), max(old[1], end)
                            additions.extend(lines)
                    if gene['gene_id'] in retired:
                        continue
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
                    elif (candidate['origin'] == 'original' or candidate.get('rescue_alternative_coding_path')
                          or candidate.get('same_coordinate_source_repair') or candidate.get('repair_kind') == 'phase_only'):
                        original_transcripts.add(candidate['source_transcript_id'])
                    else:
                        predicted_rows.append(candidate_gff(gene, candidate, name))
                    changes.append({**chosen, 'selected_origin': candidate['origin']})
                gff.write(_selected_gff_source_rows(source_graph, original_transcripts) + ''.join(predicted_rows))
            for identifier in retired:
                del selected[identifier]
            import transcript_utrs
            native_prediction = value['request'].get('oriongeno', {}).get('predictions', {}).get(name)
            native_models = transcript_utrs.load_models(native_prediction['gff']) if native_prediction else None
            with indexed_genome(source['genome']) as genome:
                augmented, utr_audit, transcripts = transcript_utrs.complete(
                    paths['gff'].read_text(), original, genome, name, native_models)
            paths['gff'].write_text(augmented)
            with paths['transcript'].open('w') as stream:
                for tid, sequence in sorted(transcripts.items()):
                    stream.write('>' + quote(tid, safe='._-:') + ' utr_policy=projected_longest_per_side\n' + sequence + '\n')
            audited_transcripts = {r['transcript_id'] for r in utr_audit}
            for gene_id, chosen in selected.items():
                if chosen['source_transcript_id'] not in audited_transcripts:
                    utr_audit.append(dict(species=name, gene_id=gene_id,
                                          transcript_id=chosen['source_transcript_id'], status='unexportable'))
            atomic_json(tmp / 'utr_audit' / (name + '.json'), utr_audit)
            utr_summaries.append((name, len(selected), len(transcripts),
                                  sum(bool(r.get('utr5', {}).get('length')) for r in utr_audit),
                                  sum(bool(r.get('utr3', {}).get('length')) for r in utr_audit),
                                  sum(bool(r.get('utr5', {}).get('length')) and bool(r.get('utr3', {}).get('length')) for r in utr_audit)))
            paths['analysis_gff'].write_text(analysis_gff_rows(augmented, analysis_candidates))
            derived_original, bound_changes = extend_gene_bounds(_selected_gff_source_rows(source_graph, set(), retain_all=True), expanded_bounds)
            gene_bounds_audit.extend(dict(species=name, **row) for row in bound_changes)
            full_path.write_text(derived_original.rstrip() + '\n' + ''.join(additions) + (sep + embedded if sep else ''))
            updates = {}
            selected_models = transcript_utrs.read_models(augmented, strict=False)
            for record in utr_audit:
                if 'utr5' not in record:
                    continue
                tid = record['transcript_id']
                updates[tid] = dict(model=selected_models[tid], exons=selected_models[tid]['exons'],
                                    choices={'5': record['utr5'], '3': record['utr3']})
            full_path.write_text(transcript_utrs.render(full_path.read_text(), updates))
            del source_graph
            map_hashes = representative_map_proof.read([tmp / 'representative_map.tsv'])
            row = {'species': name, **{k: str((root / 'effective' / p.relative_to(tmp)).resolve()) for k, p in paths.items()},
                   'representative_map': str((root / 'effective' / 'representative_map.tsv').resolve()),
                   'genetic_code': source['genetic_code'], 'representative_map_sha256': map_hashes[str(tmp / 'representative_map.tsv')]}
            for k, path in paths.items():
                row[k + '_sha256'] = digest(path)
            rows.append(row)
        rescue.write_tsv(tmp / 'effective_exclusions.tsv', ('species', 'gene_id', 'candidate_id', 'reason'), exclusions)
        atomic_json(tmp / 'gene_bounds_changes.json', gene_bounds_audit)
        atomic_json(tmp / 'unpublished_candidates.json', unpublished_candidates)
        rescue.write_tsv(tmp / 'utr_summary.tsv', ('species', 'representatives', 'exported_transcripts', 'with_5prime_utr', 'with_3prime_utr', 'with_both_utrs'), utr_summaries)
        rescue.write_tsv(tmp / 'translation_admission.tsv', ('species', 'gene_id', 'candidate_id', 'status', 'quality'), translation_audit)
        rescue.write_tsv(tmp / 'coding_admission.tsv', ('species', 'gene_id', 'candidate_id', 'status', 'head_bases_removed', 'tail_bases_removed', 'source_cds_length', 'analysis_cds_length'), coding_audit)
        rescue.write_tsv(tmp / 'species_genetic_code.tsv', ('species', 'genetic_code'), [(r['species'], r['genetic_code']) for r in rows])
        if rows:
            fields = list(rows[0])
            rescue.write_tsv(tmp / 'inputs.tsv', fields, [[r[k] for k in fields] for r in rows])
        sequence_summary = (write_sequence_view(tmp, db, request, selected_by_species, retired_by_species, coding_audit)
                            if request.get('transcriptome_sources') else {})
        atomic_json(tmp / 'changes.json', changes)
        atomic_json(tmp / 'summary.json', {'species': len(rows), 'loci': len(changes),
                                         'changed_representatives': sum(r['status'] in {'conserved', 'sequence_repaired', 'merged_split_model'} for r in changes),
                                         'same_coordinate_sequence_repairs': sum(r['status'] == 'sequence_repaired' for r in changes),
                                         'predicted_selected': sum(r['selected_origin'] == 'predicted' for r in changes),
                                         **sequence_summary})
        representative_map_proof.check()
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


def inspect_stages(root, value, cpus=1):
    """Check the required dependency graph, including absent publications.

    A receipt's own key is not evidence that it belongs to this run or its
    current upstream stages. Reconstruct those keys from the frozen plan and
    the current dependency receipts without building or repairing anything.
    """
    plan_hash = digest(root / 'plan.json')
    verify = _receipt_verifier(cpus)

    def receipt_hash(relative):
        try:
            return digest(root / relative / 'receipt.json')
        except OSError:
            return None

    catalog_dirs = {name: catalog_directory(root, value, name) for name in value['species']}
    catalogs = {name: receipt_hash(path) for name, path in catalog_dirs.items()}
    predictions = {name: receipt_hash(Path('predictions') / name) for name in target_species(value)}
    initial_index = {'catalog_index': receipt_hash(Path('catalog_index'))}
    final_index = {'catalog_index_final': receipt_hash(Path('catalog_index_final'))}
    fixed = set(value['request'].get('fixed_donors', ()))
    expected = {str(Path('catalog') / name): {} for name in value['species'] if name not in fixed}
    expected['catalog_index'] = {'dependencies': {'catalog': catalogs}}
    correspondence_dependencies = {'catalog': catalogs, 'index': initial_index}
    request = value['request']
    anchors = {'donor_catalog/' + name: verify(catalog_dirs[name], donor_catalog_key(request, name))
               for name in sorted(fixed)}
    if uses_rescue_correspondence(request):
        anchor_root = Path(request['rescue_output'])
        anchor_plan = rescue.load(anchor_root, inspect_only=True)
        jobs = [job for job in anchor_plan['synteny_jobs'] if job['a'] != job['b']]
        prepared_hashes = {}
        anchor_hash = rescue.plan_digest(anchor_root, anchor_plan)
        prepared_names = ({job[side] for job in jobs for side in ('a', 'b')}
                          if 'fixed_donors' in request else request['sources'])
        with ReceiptSnapshot(verify) as prepared_proof:
            for name in sorted(prepared_names):
                directory = anchor_root / 'prepared' / name
                member = directory / 'receipt.json'
                try:
                    valid = rescue.verified_dependency(prepared_proof, directory, {'plan': anchor_hash, 'species': name})
                    prepared_hashes[name] = digest(member) if valid else None
                    if valid:
                        prepared_proof.check({member: prepared_hashes[name]})
                except (OSError, ValueError, KeyError):
                    prepared_hashes[name] = None
                anchors['prepared/' + name] = prepared_hashes[name] is not None
            anchor_hashes = {}
            for job in jobs:
                directory = anchor_root / 'synteny' / job['id']
                try:
                    prepared = {name: prepared_hashes[name] for name in {job['a'], job['b']}}
                    if not all(prepared.values()):
                        raise ValueError('Prepared comparison inputs are unavailable')
                    anchor_hashes[job['id']] = digest(directory / 'receipt.json')
                    key = {'plan': anchor_hash, 'job': job, 'prepared': prepared}
                    anchors['anchors/' + job['id']] = rescue.verified_dependency(prepared_proof, directory, key)
                    if anchors['anchors/' + job['id']]:
                        prepared_proof.check({directory / 'receipt.json': anchor_hashes[job['id']]})
                except (OSError, ValueError, KeyError):
                    anchor_hashes[job['id']] = None
                    anchors['anchors/' + job['id']] = False
            # Ordinary missing/corrupt results never acquire a proof. A change
            # to a generation already verified during QC must still fail closed.
            prepared_proof.check_all()
        correspondence_dependencies['anchors'] = anchor_hashes
        correspondence_dependencies['prepared'] = prepared_hashes
    expected['correspondence'] = {'dependencies': correspondence_dependencies}
    selection_dependencies = {'correspondence': receipt_hash(Path('correspondence')), 'catalog': catalogs,
                              'index': initial_index}
    expected['selection_initial'] = {'dependencies': selection_dependencies}
    prediction_dependencies = {'initial': receipt_hash(Path('selection_initial')), **selection_dependencies}
    for name in target_species(value):
        dependencies = dict(prediction_dependencies)
        revision_directory = revision_worker(request, name)
        if revision_directory is not None:
            try:
                worker_hash = digest(revision_directory / 'receipt.json')
            except OSError:
                worker_hash = None
            dependencies['rescue_workers'] = {name: worker_hash}
        expected[str(Path('predictions') / name)] = {'dependencies': dependencies}
    expected['catalog_index_final'] = {'dependencies': {'catalog': catalogs, 'index': initial_index,
                                                        'predictions': predictions}}
    final_selection_dependencies = {**selection_dependencies, 'predictions': predictions, 'index': final_index}
    if set(target_species(value)) <= request.get('transcriptome_sources', {}).keys():
        final_selection_dependencies['initial'] = prediction_dependencies['initial']
    expected['selection_final'] = {'dependencies': final_selection_dependencies}
    expected['effective'] = {'dependencies': {'selection': receipt_hash(Path('selection_final')), 'catalog': catalogs,
                                            'index': final_index}}
    return {relative: verify(root / relative, {'plan': plan_hash, **key})
            for relative, key in expected.items()} | anchors


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    prep = sub.add_parser('plan')
    prep.add_argument('--output', required=True, type=Path)
    prep.add_argument('--inputs', type=Path)
    prep.add_argument('--transcriptome-assembly-dir', type=Path, help='Discover assemblies in standard transcriptome generation output')
    prep.add_argument('--assembly-genetic-codes', type=Path, help='Species genetic code TSV for discovered assemblies')
    prep.add_argument('--assembly-cpus', type=int, default=1)
    prep.add_argument('--transcriptome-inputs', type=Path, help='Assembly species TSV: species,cds,transcript_gene_map[,genetic_code]')
    prep.add_argument('--rescue-output', type=Path)
    prep.add_argument('--donor-plan', type=Path, help='Frozen rescue donor plan for one transcriptome target; genomic targets inherit it from rescue')
    prep.add_argument('--target', help='Species to discover in a shared transcriptome assembly directory')
    prep.add_argument('--edges', type=Path)
    prep.add_argument('--rna-bam-dir', type=Path, help='Directory of species-named, genome-aligned RNA BAMs')
    prep.add_argument('--rna', type=Path, help='Whole coding RNA paths: species,seqid,strand,cds_blocks(JSON),transcript_id,count TSV')
    prep.add_argument('--model-repairs', type=Path, help='Reviewed single-locus rescue candidates TSV: species,candidate_id')
    prep.add_argument('--split-model-merges', type=Path, help='Reviewed rescue merge candidates TSV: species,candidate_id')
    prep.add_argument('--oriongeno-inputs', type=Path, help='Native OrionGeno TSV: species,gff,receipt; also supplies UTRs')
    prep.add_argument('--annevo-inputs', type=Path, help='Source-bound native ANNEVO candidate TSV: species,gff,receipt')
    prep.add_argument('--species-profiles', type=Path, help='Explicit target-species prediction parameter overrides (TSV)')
    for key, default in DEFAULTS.items():
        prep.add_argument('--' + key.replace('_', '-'), type=type(default), default=default)
    assembly = sub.add_parser('prepare-transcriptome', help='Extract one CDS per assembled transcript without collapsing isoforms')
    assembly.add_argument('--assembly', type=Path, required=True)
    assembly.add_argument('--species', required=True)
    assembly.add_argument('--output', type=Path, required=True)
    grouping = assembly.add_mutually_exclusive_group(required=True)
    grouping.add_argument('--assembler', choices=['trinity', 'rnaspades', 'rna-bloom2'])
    grouping.add_argument('--gene-map', type=Path, help='Explicit transcript_id, gene_id[,baseline] TSV')
    assembly.add_argument('--corset-clusters', type=Path)
    assembly.add_argument('--genetic-code', type=int, default=1)
    assembly.add_argument('--cpus', type=int, default=1)
    for command in ('catalog', 'correspondence', 'select', 'predict', 'finalize', 'run', 'status', 'qc'):
        s = sub.add_parser(command)
        s.add_argument('--output', type=Path, required=True)
        s.add_argument('--task-index', type=int)
        s.add_argument('--cpus', type=int, default=1)
        s.add_argument('--comparison-cache', type=Path)
    check = sub.add_parser('verify-inputs')
    check.add_argument('--inputs', required=True, type=Path)
    check.add_argument('--field', choices=['cds', 'protein', 'gff', 'genome', 'analysis_cds', 'analysis_gff', 'representative_map', 'transcript', 'layout', 'analysis_layout', 'coding_layout'])
    utr = sub.add_parser('reproject-utrs', help='Re-evaluate UTRs from a completed, immutable coding-model view')
    utr.add_argument('--source-output', required=True, type=Path)
    utr.add_argument('--output', required=True, type=Path)
    utr.add_argument('--source-receipt-sha256', required=True)
    utr.add_argument('--species', action='append', help='Repeat for a subset; defaults to the complete source cohort')
    return p


def _main(args):
    if args.command == 'prepare-transcriptome':
        transcriptome.prepare(args.assembly, args.species, args.output, args.genetic_code, args.cpus,
                              args.assembler, args.corset_clusters, args.gene_map)
        return
    if args.command == 'reproject-utrs':
        import transcript_utr_reassessment
        transcript_utr_reassessment.reproject(args.source_output, args.output, args.source_receipt_sha256, args.species)
        return
    if args.command == 'verify-inputs':
        result = verify_inputs(args.inputs, args.field)
        print(result if isinstance(result, str) else json.dumps(result, sort_keys=True))
        return
    if args.command == 'plan':
        plan(args.output, args.inputs, args.rescue_output, args.edges, args.rna,
             species_profiles=args.species_profiles, annevo_inputs=args.annevo_inputs, oriongeno_inputs=args.oriongeno_inputs,
             split_model_merges=args.split_model_merges, model_repairs=args.model_repairs,
             transcriptome_inputs=args.transcriptome_inputs, rna_bam_dir=args.rna_bam_dir,
             transcriptome_assembly_dir=args.transcriptome_assembly_dir,
             assembly_genetic_codes=args.assembly_genetic_codes, assembly_cpus=args.assembly_cpus, donor_plan=args.donor_plan, target=args.target,
             **{k: getattr(args, k) for k in DEFAULTS})
        return
    if args.cpus < 1:
        raise ValueError('--cpus must be positive')
    root = args.output.resolve()
    value = load(root)
    names = target_species(value) if args.command in {'predict', 'run'} else value['species']
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
        if args.command == 'qc' and set(target_species(value)) & value['request']['sources'].keys():
            verify_inputs(root / 'effective' / 'inputs.tsv')
        result = {'species': len(target_species(value)), 'stages': inspect_stages(root, value, args.cpus)}
        if args.command == 'qc' and not all(result['stages'].values()):
            invalid = [name for name, valid in result['stages'].items() if not valid]
            raise ValueError('Corrupt refinement stage or missing publication: ' + ', '.join(invalid))
        print(json.dumps(result, sort_keys=True))
    load(root)


def main():
    args = parser().parse_args()
    with invocation_context(cpus=max(1, getattr(args, 'cpus', 1))):
        _main(args)


if __name__ == '__main__':
    main()

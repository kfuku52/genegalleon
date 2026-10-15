#!/usr/bin/env python3
"""Source-bound ANNEVO candidates for the existing conservative model gates.

ANNEVO supplies coding paths, never independent biological support. Protein
support is measured against each nominated donor; callers retain ownership,
genomic QC, independent-species, admission and representative-selection rules.
"""
import argparse
import bisect
import csv
import json
import os
import re
import tempfile
from collections import defaultdict, deque
from concurrent.futures import ProcessPoolExecutor
from itertools import chain, islice
from multiprocessing import get_context
from pathlib import Path

from Bio.Seq import Seq

try:
    from fasta_sequence_store import exclusive_lock
except ImportError:
    # Establish the shared support import path for package callers before
    # importing the model selector's transitive support dependencies.
    from .fasta_sequence_store import exclusive_lock

try:
    from annevo_model_selection import PROFILES
    from format_species_annotation.common import parse_gff_attributes
    from gene_model_selection import pair_score
    from genome_sequence_identity import genome_sequence_digest
    from input_generation_array_state import atomic_json, digest, digest_paths
except ImportError:
    from .annevo_model_selection import PROFILES
    from .format_species_annotation.common import parse_gff_attributes
    from .gene_model_selection import pair_score
    from .genome_sequence_identity import genome_sequence_digest
    from .input_generation_array_state import atomic_json, digest, digest_paths


def prediction_identity(gff, receipt, genome_sha256, genetic_code=1, genome=None):
    """Verify actual GG native prediction provenance before importing a path."""
    receipt = Path(receipt)
    raw_sha = digest(receipt)
    value = json.loads(receipt.read_text())
    if digest(receipt) != raw_sha:
        raise ValueError('ANNEVO receipt changed during reading')
    if not isinstance(value, dict):
        raise ValueError('ANNEVO native receipt must be an object')
    gpu = value.get('gpu_probe')
    if (type(value.get('schema_version')) is not int or value['schema_version'] != 1
            or value.get('kind') != 'genegalleon-annevo-native-prediction'
            or type(value.get('native_returncode')) is not int or value['native_returncode'] != 0
            or not isinstance(gpu, dict) or gpu.get('gpu_verified') is not True
            or value.get('model_profile') not in PROFILES or genetic_code != 1):
        raise ValueError('ANNEVO candidates require a successful supported native GPU receipt and genetic code 1')
    settings = value.get('inference_settings')
    # Old Magnoliopsida receipts predate inference_settings. Newly supported
    # profiles must retain the producer's explicit model/lineage pairing.
    if ((settings is not None or value['model_profile'] != 'Magnoliopsida')
            and (not isinstance(settings, dict) or settings.get('lineage') != value['model_profile'])):
        raise ValueError('ANNEVO candidate model profile differs from its native inference lineage')
    genome_match = 'raw_identity'
    if value.get('genome_sha256') != genome_sha256:
        expected_sequence = value.get('genome_sequence_sha256')
        if (genome is None or not re.fullmatch('[0-9a-f]{64}', str(expected_sequence))
                or genome_sequence_digest(genome) != expected_sequence or digest(genome) != genome_sha256):
            raise ValueError('ANNEVO candidate genome differs from the target assembly bytes or contig sequences')
        genome_match = 'contig_sequence_identity'
    gff = Path(gff)
    if (type(value.get('GFF_bytes')) is not int or value.get('GFF_sha256') != digest(gff)
            or value['GFF_bytes'] != gff.stat().st_size):
        raise ValueError('ANNEVO candidate GFF differs from its native receipt')
    for field, size in (('runtime_digest', 64), ('model_sha256', 64), ('source_commit', 40), ('cache_key', 64),
                        ('genome_sha256', 64)):
        if not re.fullmatch('[0-9a-f]{' + str(size) + '}', str(value.get(field, ''))):
            raise ValueError('ANNEVO receipt lacks its exact ' + field)
    return {key: value[key] for key in ('genome_sha256', 'GFF_sha256', 'GFF_bytes', 'runtime_digest',
                                      'model_sha256', 'source_commit', 'model_profile', 'cache_key')} | {
                                          'receipt_sha256': raw_sha, 'target_genome_sha256': genome_sha256,
                                          'genome_match': genome_match,
                                          **({'genome_sequence_sha256': value['genome_sequence_sha256']}
                                             if value.get('genome_sequence_sha256') else {})}


def freeze_inputs(path, sources):
    """Freeze a subset of target species; missing rows leave their route unchanged."""
    table = Path(path).resolve()
    table_sha = digest(table)
    result, files = {}, {str(table): table_sha}
    with table.open(newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if reader.fieldnames != ['species', 'gff', 'receipt']:
            raise ValueError('ANNEVO input TSV requires species, gff, receipt columns in that order')
        for row in reader:
            if None in row or any(not value for value in row.values()):
                raise ValueError('Incomplete ANNEVO input row')
            name = row['species']
            if name not in sources or name in result:
                raise ValueError('Unknown or duplicate ANNEVO target species: ' + name)
            locations = {}
            for role in ('gff', 'receipt'):
                location = Path(row[role])
                locations[role] = str((location if location.is_absolute() else table.parent / location).resolve())
            source = sources[name]
            genome_sha = digest(source['genome'])
            identity = prediction_identity(**locations, genome_sha256=genome_sha,
                                           genetic_code=source.get('genetic_code', 1), genome=source['genome'])
            result[name] = {**locations, 'identity': identity}
            files.update({locations['gff']: identity['GFF_sha256'], locations['receipt']: identity['receipt_sha256']})
    if not result:
        raise ValueError('ANNEVO input TSV contains no predictions')
    if digest_paths(files) != files:
        raise ValueError('ANNEVO inputs changed during planning')
    return {'table': str(table), 'predictions': result, 'files': files}


def verify_inputs(frozen, names=None):
    if not frozen:
        return
    selected = frozen['predictions'] if names is None else {n: frozen['predictions'][n] for n in names if n in frozen['predictions']}
    paths = [frozen['table'], *(p[role] for p in selected.values() for role in ('gff', 'receipt'))]
    expected = {p: frozen['files'][p] for p in paths}
    if digest_paths(expected) != expected:
        raise ValueError('Frozen ANNEVO candidate input changed')


def read_paths(path, genome):
    """Read native GFF relationships and zero-based CDS blocks without repair."""
    features, transcripts, cds = {}, {}, defaultdict(list)
    references = set(genome.references)
    with Path(path).open() as stream:
        for number, line in enumerate(stream, 1):
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\r\n').split('\t')
            if len(f) != 9 or f[2] not in {'gene', 'mRNA', 'exon', 'CDS'}:
                raise ValueError('Unsupported native ANNEVO GFF row: ' + str(number))
            parsed = parse_gff_attributes(f[8])
            if f[2] == 'gene' and parsed.get('Parent'):
                raise ValueError('ANNEVO gene has an unexpected parent')
            if any(len(parsed.get(key, ())) != 1 for key in ('ID',) + (() if f[2] == 'gene' else ('Parent',))):
                raise ValueError('ANNEVO GFF requires unique feature and parent identifiers')
            attrs = {key: values[0] for key, values in parsed.items() if len(values) == 1}
            identifier = attrs.get('ID')
            if not identifier or identifier in features:
                raise ValueError('Missing or duplicate ANNEVO feature ID')
            start, end = int(f[3]) - 1, int(f[4])
            if (f[0] not in references or f[6] not in {'+', '-'}
                    or not 0 <= start < end <= genome.get_reference_length(f[0])):
                raise ValueError('ANNEVO GFF feature outside target genome or unknown strand')
            features[identifier] = (f[2], f[0], start, end, f[6], attrs.get('Parent'))
            if f[2] == 'mRNA':
                transcripts[identifier] = {'native_transcript_id': identifier, 'native_gene_id': attrs.get('Parent'),
                                           'seqid': f[0], 'strand': f[6]}
            elif f[2] == 'CDS':
                if f[7] not in {'0', '1', '2'} or not attrs.get('Parent'):
                    raise ValueError('ANNEVO CDS lacks an explicit phase or parent')
                cds[attrs['Parent']].append([start, end, int(f[7])])
    for identifier, (kind, seqid, start, end, strand, parent) in features.items():
        if kind == 'gene':
            if parent:
                raise ValueError('ANNEVO gene has an unexpected parent')
            continue
        owner = features.get(parent)
        required = 'gene' if kind == 'mRNA' else 'mRNA'
        if not owner or owner[0] != required or owner[1] != seqid or owner[4] != strand or not owner[2] <= start < end <= owner[3]:
            raise ValueError('ANNEVO GFF has an invalid feature ownership relationship: ' + identifier)
    paths = []
    for identifier, transcript in transcripts.items():
        blocks = sorted(cds.get(identifier, []), reverse=transcript['strand'] == '-')
        if not blocks:
            raise ValueError('ANNEVO transcript has no coding path')
        ordered = sorted(blocks)
        if any(a[1] > b[0] for a, b in zip(ordered, ordered[1:], strict=False)):
            raise ValueError('ANNEVO coding blocks overlap')
        paths.append({**transcript, 'cds': blocks,
                      'start': min(b[0] for b in blocks), 'end': max(b[1] for b in blocks)})
    return paths


def _score_protein_batch(pairs):
    """Only strings cross the process boundary; genome handles stay in caller."""
    return [pair_score({'protein': a}, {'protein': b}) for a, b in pairs]


def _scored_pairs(pairs, cpus):
    """Preserve pair order with at most two small pending batches per CPU."""
    if cpus == 1:
        for path, region, native, donor in pairs:
            yield path, region, pair_score({'protein': native}, {'protein': donor})
        return
    batch_size = 64
    initial = list(islice(pairs, max(512, cpus * batch_size)))
    if len(initial) < 512:
        for path, region, native, donor in initial:
            yield path, region, pair_score({'protein': native}, {'protein': donor})
        return
    workers = min(cpus, (len(initial) + batch_size - 1) // batch_size)
    pending = deque()
    remaining = chain(initial, pairs)
    # Spawn prevents inheriting the caller's indexed genome or SQLite handles.
    with ProcessPoolExecutor(max_workers=workers, mp_context=get_context('spawn')) as executor:
        while batch := list(islice(remaining, batch_size)):
            future = executor.submit(_score_protein_batch, [(a, b) for _, _, a, b in batch])
            pending.append((batch, future))
            if len(pending) >= workers * 2:
                ready, result = pending.popleft()
                for (path, region, _, _), score in zip(ready, result.result(), strict=True):
                    yield path, region, score
        for ready, result in pending:
            for (path, region, _, _), score in zip(ready, result.result(), strict=True):
                yield path, region, score


def nominate(frozen, name, regions, proteins, genome, params, *, cpus=1):
    """Align only paths fully inside nominated windows to their actual donors."""
    if type(cpus) is not int or cpus < 1:
        raise ValueError('ANNEVO comparison CPUs must be a positive integer')
    if not frozen or name not in frozen['predictions']:
        return [], {'native_paths': 0, 'candidate_pairs': 0, 'proposals': []}
    verify_inputs(frozen, [name])
    prediction = frozen['predictions'][name]
    paths = read_paths(prediction['gff'], genome)
    tracks = defaultdict(list)
    for path in paths:
        tracks[path['seqid']].append(path)
    indexes = {}
    for seqid, rows in tracks.items():
        rows.sort(key=lambda p: (p['start'], p['end'], p['native_transcript_id']))
        indexes[seqid] = [p['start'] for p in rows], rows
    results, proposals, nominated, protein_cache = [], [], set(), {}
    cap = int(params.get('candidate_limit', 30))
    def candidate_pairs():
        for region in regions:
            # A whole-genome nomination cannot assign an ANNEVO hit to a locus.
            # Existing miniprot genomic placement supplies a bounded target span.
            if region.get('genome_only') and not region.get('annevo_bounded_placement'):
                continue
            starts, rows = indexes.get(region['seqid'], ([], []))
            hits = [p for p in rows[bisect.bisect_left(starts, region['start']):bisect.bisect_left(starts, region['end'])]
                    if p['end'] <= region['end']]
            if len(hits) > cap:
                proposals.append({'query': region['id'], 'reason': 'annevo_candidate_limit_exceeded', 'paths': len(hits)})
                continue
            for path in hits:
                nominated.add(path['native_transcript_id'])
                token = path['native_transcript_id']
                if token not in protein_cache:
                    sequence = ''.join(str(Seq(genome.fetch(path['seqid'], a, b)).reverse_complement())
                                       if path['strand'] == '-' else genome.fetch(path['seqid'], a, b)
                                       for a, b, _ in path['cds']).upper()
                    protein_cache[token] = str(Seq(sequence).translate()).removesuffix('*') if len(sequence) % 3 == 0 else ''
                donor = proteins[region['donor']][region['query']]
                yield path, region, protein_cache[token], donor
    for path, region, score in _scored_pairs(candidate_pairs(), cpus):
        alignment = score.as_dict()
        # Neither input supplies donor exon coordinates. Equal empty
        # junction lists are not evidence of conserved target splicing.
        # Retain the existing raw PairScore fields for saved-data readers,
        # while explicitly marking their missing structural input.
        alignment['junction_information_available'] = False
        alignment['score_basis'] = 'protein_only_without_donor_junctions'
        problems = []
        if score.bounded:
            problems.append('annevo_alignment_bounded')
        if min(score.coverage_a, score.coverage_b) < params['minimum_coverage']:
            problems.append('low_reciprocal_coverage')
        evidence = {**region, 'annevo': {**prediction['identity'],
                    'native_transcript_id': path['native_transcript_id'], 'native_gene_id': path['native_gene_id'],
                    'alignment': alignment}}
        results.append({'seqid': path['seqid'], 'strand': path['strand'], 'cds': [list(b) for b in path['cds']],
                        'query': region['id'], 'coverage': score.coverage_b, 'identity': score.identity,
                        'frameshift': False, 'evidence': evidence, 'search': 'annevo_nominated_path',
                        'annevo_problems': problems})
    proposals.extend({'native_transcript_id': p['native_transcript_id'], 'reason': 'no_bounded_donor_nomination'}
                     for p in paths if p['native_transcript_id'] not in nominated)
    verify_inputs(frozen, [name])
    return results, {'native_paths': len(paths), 'candidate_pairs': len(results), 'proposals': proposals}


def gff_attributes(evidence):
    """Keep engine provenance distinct from independent donor-species support."""
    from urllib.parse import quote
    origins = [r['annevo'] for r in evidence if r.get('annevo')]
    if not origins:
        return ''
    ids = sorted({r['native_transcript_id'] for r in origins})
    return ';prediction_engine=ANNEVO;annevo_native_transcripts=' + ','.join(quote(i, safe='._-:') for i in ids)


def confirm_transferred_structure(native_models, transferred):
    """Protein similarity alone does not corroborate a target exon chain.

    Require the same donor query's actual genomic mapping to retain every
    internal junction and compatible coding coordinates. Partial start/stop
    transfer can still corroborate internal structure. Another ANNEVO row
    cannot supply this check, and this never adds a donor species.
    """
    from rescue_coding_paths import coding_shape, compatible_paths
    def key(model):
        evidence = model.get('evidence', {})
        return (model.get('gene_id', evidence.get('target')), model.get('donor_species', evidence.get('donor')),
                model.get('donor_candidate', evidence.get('query')), model['seqid'], model['strand'])
    def junctions(model):
        ordered = sorted(model['cds'])
        return tuple((a[1], b[0]) for a, b in zip(ordered, ordered[1:], strict=False))
    def checked_support(model):
        donor, query = key(model)[1:3]
        return {**model, 'problems': [], 'support': [],
                'evidence': {**model.get('evidence', {}), 'donor': donor, 'query': query}}
    allowed = {'missing_start', 'missing_stop', 'unanchored_genome_search'}
    by_query = defaultdict(list)
    for model in transferred:
        if (not model.get('evidence', {}).get('annevo') and model.get('sequence') and model.get('cds')
                and not set(model.get('problems', [])) - allowed):
            by_query[key(model)].append(model)
    for model in native_models:
        model['problems'] = [p for p in model['problems'] if p != 'annevo_unconfirmed_coding_structure']
        matches = [other for other in by_query[key(model)] if junctions(other) == junctions(model)
                   and compatible_paths(checked_support(model), checked_support(other))]
        if model.get('evidence', {}).get('genome_only'):
            # No existing target model corroborates either terminus here.
            # Do not turn a partial/paralog placement into a new complete gene.
            matches = [other for other in matches if coding_shape(other) == coding_shape(model)
                       and not {'missing_start', 'missing_stop'} & set(other.get('problems', []))]
        model['evidence']['annevo']['genomic_structure_transfer'] = {
            'confirmed': bool(matches), 'counting_unit': 'same_donor_query_genomic_mapping',
            'scope': 'whole_coding_path' if model.get('evidence', {}).get('genome_only') else 'internal_junctions',
            'matching_paths': [{'query': row['query'], 'cds': row['cds']} for row in matches]}
        if not matches:
            model['problems'] = sorted(set(model['problems']) | {'annevo_unconfirmed_coding_structure'})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--species', required=True)
    parser.add_argument('--genome', required=True, type=Path)
    parser.add_argument('--gff', required=True, type=Path)
    parser.add_argument('--receipt', required=True, type=Path)
    parser.add_argument('--genetic-code', type=int, default=1)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if not re.fullmatch('[A-Za-z0-9_.-]+', args.species):
        raise ValueError('Invalid species identifier')
    identity = prediction_identity(args.gff, args.receipt, digest(args.genome),
                                   genetic_code=args.genetic_code, genome=args.genome)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with exclusive_lock(args.output.with_suffix(args.output.suffix + '.lock')):
        with tempfile.NamedTemporaryFile(mode='w', dir=args.output.parent, newline='', delete=False) as stream:
            temporary = Path(stream.name)
            writer = csv.writer(stream, delimiter='\t')
            writer.writerow(['species', 'gff', 'receipt'])
            writer.writerow([args.species, str(args.gff.resolve()), str(args.receipt.resolve())])
            stream.flush()
            os.fsync(stream.fileno())
        try:
            if args.output.is_symlink() or (args.output.exists() and args.output.read_bytes() != temporary.read_bytes()):
                raise ValueError('ANNEVO candidate table already exists with different inputs')
            provenance = args.output.with_suffix(args.output.suffix + '.provenance.json')
            if provenance.is_symlink():
                raise ValueError('ANNEVO candidate provenance is redirected')
            if not provenance.exists() or json.loads(provenance.read_text()) != identity:
                atomic_json(provenance, identity)
            if not args.output.exists():
                temporary.replace(args.output)
        finally:
            temporary.unlink(missing_ok=True)


if __name__ == '__main__':
    main()

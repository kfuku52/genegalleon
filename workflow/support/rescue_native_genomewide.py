"""Candidate-first native rescue outside the original coding annotation.

Protein search nominates real external donors, never a genomic placement or
orthology. All returned paths still require the ordinary genomic, terminal,
unanchored placement, ownership and coding-path admission checks.
"""
import copy
import hashlib
import json
import time
from array import array
from collections import Counter, defaultdict
from pathlib import Path

from Bio.Seq import Seq

SEARCH = 'oriongeno_genomewide_path'
CONTRACT = {'schema': 1, 'scope': 'all_native_paths_without_original_cds_overlap',
            'donors': 'all_prepared_proteins_of_frozen_selected_external_species',
            'search': 'diamond_ultra_sensitive_unmasked_evalue_1000_no_target_cap',
            'admission': 'existing_genomic_and_unique_external_donor_gates'}


def _r():
    try:
        import rescue_gene_models
    except ImportError:
        from . import rescue_gene_models
    return rescue_gene_models


def _native():
    try:
        import annevo_gene_models
    except ImportError:
        from . import annevo_gene_models
    return annevo_gene_models


def _candidate_protein_pairs(path, queries, subjects, minimum):
    """Retain hit IDs compactly and yield unique pairs in their original sort order."""
    subject_ids = sorted(subjects)
    subject_indexes = {token: index for index, token in enumerate(subject_ids)}
    subject_lengths = [len(subjects[token].upper().removesuffix('*')) for token in subject_ids]
    query_lengths = {token: len(protein) for token, protein in queries.items()}
    candidates = defaultdict(lambda: array('I'))
    rejected = Counter()
    with path.open() as stream:
        for line in stream:
            fields = line.split()
            if len(fields) != 2:
                raise ValueError('Unknown genomewide protein search identifier')
            q, s = fields
            a, index = query_lengths.get(q), subject_indexes.get(s)
            if a is None or index is None:
                raise ValueError('Unknown genomewide protein search identifier')
            b = subject_lengths[index]
            # Even perfect pairing cannot align more residues than the shorter sequence.
            if (a / b if a < b else b / a) < minimum:
                rejected[q] += 1
            else:
                candidates[q].append(index)
    # Exact short peptides need not seed a DIAMOND alignment.
    by_sequence = {protein: subject_indexes[token] for token, protein in subjects.items()}
    for q, protein in queries.items():
        if protein in by_sequence:
            candidates[q].append(by_sequence[protein])

    def comparisons():
        for q in sorted(candidates):
            # Deduplicate only the current query, including hits also added by
            # the exact-match pass. Never keep a Python set for every query.
            for index in sorted(set(candidates.pop(q))):
                s = subject_ids[index]
                yield q, s, queries[q], subjects[s]

    return rejected, comparisons()


def run_protein_search(queries, subjects, search_dir, cpus):
    """Produce the complete DIAMOND pair table before alignment and validation."""
    r = _r()
    for file, sequences in (('native.pep', queries), ('donors.pep', subjects)):
        with (search_dir / file).open('w') as out:
            for token, protein in sequences.items():
                out.write(f'>{token}\n{protein}\n')
    r.run(['diamond', 'makedb', '--in', search_dir / 'donors.pep',
           '--db', search_dir / 'donors', '--threads', cpus], search_dir, 'makedb')
    r.run(['diamond', 'blastp', '--db', search_dir / 'donors', '--query', search_dir / 'native.pep',
           '--out', search_dir / 'pairs.tsv', '--outfmt', '6', 'qseqid', 'sseqid',
           '--ultra-sensitive', '--masking', '0', '--comp-based-stats', '0', '--evalue', '1000',
           '--max-target-seqs', '0', '--max-hsps', '1', '--threads', cpus], search_dir, 'blastp')
    (search_dir / 'donors.dmnd').unlink(missing_ok=True)


def nominate(frozen, name, proteins, genome, ownership, params, directory, *, cpus=1, validate=None,
             search=run_protein_search):
    """Search each eligible native isoform, retaining identical-protein aliases.

    There is no donor representation filter, candidate window or query cap. Search
    negatives are explicitly search negatives, not proof that no homolog exists.
    Only sequence-identical aliases share an alignment; loci/isoforms stay distinct.
    Validate each result in the caller while the bounded scoring pool runs,
    without retaining separate raw and validated lists.
    """
    if type(cpus) is not int or cpus < 1:
        raise ValueError('Native genomewide CPUs must be a positive integer')
    if not frozen or name not in frozen['predictions']:
        return [], {'contract': CONTRACT, 'native_paths': 0, 'paths': [], 'candidate_pairs': 0,
                    'native_prediction_available': False}
    native = _native()
    native.verify_inputs(frozen, [name])
    prediction = frozen['predictions'][name]
    paths = native.read_paths(prediction['gff'], genome, engine='oriongeno')
    audit = {'contract': CONTRACT, 'native_paths': len(paths), 'paths': [], 'candidate_pairs': 0}
    groups, rows = defaultdict(list), {}
    started = time.perf_counter()
    for path in paths:
        token = path['native_transcript_id']
        if token in rows:
            raise ValueError('Duplicate native transcript in genomewide rescue')
        owners = ownership.overlapping(path)
        row = {'native_transcript_id': token, 'native_gene_id': path['native_gene_id'],
               'seqid': path['seqid'], 'strand': path['strand'], 'cds': path['cds']}
        rows[token] = row
        audit['paths'].append(row)
        if owners:
            row.update(status='original_cds_overlap', original_owner_ids=sorted({p['gene_id'] for p in owners}))
            continue
        dna = ''.join(str(Seq(genome.fetch(path['seqid'], a, b)).reverse_complement())
                      if path['strand'] == '-' else genome.fetch(path['seqid'], a, b)
                      for a, b, _ in path['cds']).upper()
        if len(dna) % 3 or set(dna) - set('ACGT'):
            row.update(status='invalid_native_cds', reasons=['incomplete_frame_or_ambiguous_cds'])
            continue
        protein = str(Seq(dna).translate()).removesuffix('*')
        if not protein or '*' in protein:
            row.update(status='invalid_native_cds', reasons=['empty_translation_or_internal_stop'])
            continue
        row['status'] = 'no_protein_search_hit'
        groups[protein].append(path)
    donor_aliases = defaultdict(list)
    for donor in sorted(set(proteins) - {name}):
        for query, protein in sorted(proteins[donor].items()):
            token = 'orion_genome_query_' + hashlib.sha256(
                json.dumps([name, donor, query], separators=(',', ':')).encode()).hexdigest()[:20]
            donor_aliases[protein].append((donor, query, token))
    if not donor_aliases:
        raise ValueError('Genomewide native rescue has no external donor proteins')
    queries = {f'n{i}': protein for i, protein in enumerate(sorted(groups))}
    subjects = {f'd{i}': protein for i, protein in enumerate(sorted(donor_aliases))}
    search_dir = Path(directory) / 'native_genomewide_search'
    search_dir.mkdir()
    rejected_hits = Counter()
    comparisons = ()
    if queries:
        search(queries, subjects, search_dir, cpus)
        rejected_hits, comparisons = _candidate_protein_pairs(search_dir / 'pairs.tsv', queries, subjects,
                                                               params['minimum_coverage'])
        for token, count in rejected_hits.items():
            for path in groups[queries[token]]:
                rows[path['native_transcript_id']].update(status='evaluated', length_rejected_protein_hits=count)
    aligned_pairs = 0
    results = []
    validation_seconds = 0.0
    for q, s, score in native._scored_pairs(iter(comparisons), cpus):
        aligned_pairs += 1
        alignment = {**score.as_dict(), 'junction_information_available': False,
                     'score_basis': 'protein_only_without_donor_junctions'}
        problems = (['oriongeno_alignment_bounded'] if score.bounded else [])
        if min(score.coverage_a, score.coverage_b) < params['minimum_coverage']:
            problems.append('low_reciprocal_coverage')
        for path in groups[queries[q]]:
            rows[path['native_transcript_id']]['status'] = 'evaluated'
            for donor, query, token in donor_aliases[subjects[s]]:
                evidence = {'target': name, 'donor': donor, 'query': query, 'genome_only': True,
                            'placement': 'unanchored', 'orthology': 'unassigned', 'expected_copy': 'unassigned',
                            'native_search_scope': CONTRACT['scope'], 'oriongeno': {
                                **prediction['identity'], 'native_transcript_id': path['native_transcript_id'],
                                'native_gene_id': path['native_gene_id'], 'alignment': alignment}}
                evidence['id'] = token
                model = {'seqid': path['seqid'], 'strand': path['strand'], 'cds': copy.deepcopy(path['cds']),
                         'query': evidence['id'], 'coverage': score.coverage_b, 'identity': score.identity,
                         'frameshift': False, 'evidence': evidence, 'search': SEARCH,
                         'oriongeno_problems': list(problems)}
                if validate is not None:
                    validation_started = time.perf_counter()
                    model = validate(model)
                    validation_seconds += time.perf_counter() - validation_started
                results.append(model)
    native.verify_inputs(frozen, [name])
    audit.update(unique_query_proteins=len(queries), unique_donor_proteins=len(subjects),
                 external_donor_records=sum(map(len, donor_aliases.values())),
                 aligned_unique_protein_pairs=aligned_pairs, candidate_pairs=len(results),
                 length_rejected_protein_hits=sum(rejected_hits.values()),
                 search_and_alignment_seconds=time.perf_counter() - started - validation_seconds)
    return results, audit


def finish_audit(audit, models):
    """Account for every native path after the actual shared admission decision."""
    outcomes = {}
    for model in models:
        if model.get('search') != SEARCH:
            continue
        evidence = model['evidence']
        token = evidence['oriongeno']['native_transcript_id']
        if token not in outcomes:
            outcomes[token] = Counter(), set(), set(), set()
        statuses, reasons, donors, identifiers = outcomes[token]
        statuses[model['status']] += 1
        reasons.update(model['problems'])
        if model['status'] in {'accepted', 'duplicate_support', 'accepted_alternative_path'}:
            donors.add(evidence['donor'])
            if model.get('model_id'):
                identifiers.add(model['model_id'])
    for row in audit['paths']:
        if row['status'] != 'evaluated':
            continue
        statuses, reasons, donors, identifiers = outcomes.get(row['native_transcript_id'], ({}, set(), (), ()))
        if row.get('length_rejected_protein_hits'):
            reasons.add('low_reciprocal_coverage')
        elif not statuses:
            raise ValueError('Evaluated native path disappeared before admission')
        row.update(status='admitted' if donors else 'not_admitted',
                   candidate_records=sum(statuses.values()), admission_statuses=dict(statuses),
                   reasons=sorted(reasons), external_donor_species=sorted(donors), model_ids=sorted(identifiers))
    audit['path_status_counts'] = dict(Counter(row['status'] for row in audit['paths']))
    audit['native_gene_counts'] = dict(Counter(
        status for status in audit['path_status_counts']
        for _ in {row['native_gene_id'] for row in audit['paths'] if row['status'] == status}))
    audit['native_gene_count_note'] = 'A gene with multiple isoforms can occur in multiple status categories.'
    if sum(audit['path_status_counts'].values()) != audit['native_paths']:
        raise ValueError('Native path inventory does not reconcile')
    return audit


def transfer_structure(models, proteins, genome, code, params, directory, *, cpus=1):
    """Map nominated donors to real candidate loci without inventing anchors.

    The ordinary native admission requires an independently transferred coding
    structure. Previously unsearched/represented donors need this actual Miniprot
    mapping as well. Candidate coordinates only define a search window; every
    returned path stays unanchored and must pass normal genomic/locus checks.
    """
    r = _r()
    placement_flags = {'unanchored_genome_search', 'outside_expected_synteny_interval'}
    windows, evidence = defaultdict(dict), {}
    for model in models:
        if set(model['problems']) - placement_flags:
            continue
        start, end = min(b[0] for b in model['cds']), max(b[1] for b in model['cds'])
        if end - start > params['max_interval']:
            model['problems'].append('native_span_exceeds_max_interval')
            continue
        padding = params['terminal_max_extension']
        window = (model['seqid'], max(0, start - padding),
                  min(genome.get_reference_length(model['seqid']), end + padding))
        base = model['evidence']
        token = 'orion_genome_placement_' + hashlib.sha256(
            json.dumps([base['id'], *window], separators=(',', ':')).encode()).hexdigest()[:20]
        if token in evidence:
            ids = evidence[token]['native_genomewide_placement']['native_transcript_ids']
            if base['oriongeno']['native_transcript_id'] not in ids:
                ids.append(base['oriongeno']['native_transcript_id'])
            continue
        region = {key: value for key, value in base.items() if key != 'oriongeno'}
        region.update(id=token, seqid=window[0], start=window[1], end=window[2],
                      native_genomewide_placement={'basis': 'native_coding_coordinates_only',
                                                  'synteny_assignment': False,
                                                  'native_identity': {k: v for k, v in base['oriongeno'].items()
                                                                      if k not in {'alignment', 'genomic_structure_transfer'}},
                                                  'native_transcript_ids': [base['oriongeno']['native_transcript_id']]})
        evidence[token] = region
        windows[window][token] = region
    if not windows:
        return []
    target = Path(directory) / 'native_genomewide_transfer'
    target.mkdir()
    predictions = r.search_intervals(target, {key: list(rows.values()) for key, rows in windows.items()},
                                     proteins, genome, code, params['max_intron'], cpus, compact=True)
    for model in predictions:
        region = evidence[model['query']]
        model['seqid'] = region['seqid']
        model['cds'] = [[a + region['start'], b + region['start'], p] for a, b, p in model['cds']]
        model['evidence'] = region
        model['search'] = 'native_genomewide_miniprot'
    return predictions


def build_plan(previous_output, output, *, producer_source=None, completed_output=None,
               completed_producer_source=None, compact_diagnostics=False):
    """Prefer verified coding evidence, retaining ordinary raw cache replay."""
    r = _r()
    import rescue_admission_replay
    import rescue_native_genomewide_reuse

    if completed_producer_source and not completed_output:
        raise ValueError('A completed producer source requires its output directory')
    if producer_source:
        return rescue_admission_replay.build_plan(
            previous_output, producer_source, output, native_extension=True,
            completed_output=completed_output, completed_producer_source=completed_producer_source,
            compact_diagnostics=compact_diagnostics)
    previous, root = Path(previous_output).resolve(strict=True), Path(output).resolve()
    reason = 'Archived producer source unavailable'
    for source in rescue_native_genomewide_reuse.archived_sources(previous):
        try:
            return rescue_admission_replay.build_plan(
                previous, source, root, native_extension=True, automatic=True,
                completed_output=completed_output, completed_producer_source=completed_producer_source,
                compact_diagnostics=compact_diagnostics)
        except rescue_native_genomewide_reuse.IncompatibleEvidence as exc:
            reason = str(exc)
    if root.exists() or root == previous or previous in root.parents or root in previous.parents:
        raise ValueError('Native genomewide extension requires a new disjoint output directory')
    prior = r.load(previous, inspect_only=True)
    if (not prior['request'].get('oriongeno') or prior['request'].get('admission_cache')
            or prior['request'].get('oriongeno_genomewide')):
        raise ValueError('Native extension requires an ordinary completed native rescue')
    frozen = r.frozen_prediction_cache_key(previous)
    if set(frozen['species']) != set(prior['species']):
        raise ValueError('Native extension requires every completed previous worker')
    plan = copy.deepcopy(prior)
    plan['output'] = str(root)
    plan['request']['tools'] = r.identities()
    if any(prior['request']['tools'][key] != plan['request']['tools'][key]
           for key in ('miniprot', 'miniprot_sha256')):
        raise ValueError('Native extension requires the original Miniprot binary for raw cache reuse')
    plan['request']['prediction_cache'] = frozen
    plan['request']['oriongeno_genomewide'] = CONTRACT
    plan['request']['native_genomewide_parent'] = {
        'root': str(previous), 'plan_sha256': frozen['plan_sha256'],
        'validated_reuse_unavailable': reason,
        'prepared': {n: r.digest(previous / 'prepared' / n / 'receipt.json') for n in prior['species']},
        'synteny': {j['id']: r.digest(previous / 'synteny' / j['id'] / 'receipt.json') for j in prior['synteny_jobs']}}
    if compact_diagnostics:
        plan['request']['native_genomewide_parent']['compact_diagnostics'] = True
        plan['request']['output_storage'] = {**r.output_storage(prior), 'diagnostics': 'summary'}
    if completed_output:
        plan['request']['native_genomewide_worker_imports'] = rescue_native_genomewide_reuse.freeze_completed(
            plan, completed_output, completed_producer_source)
    metadata = {name: (previous / name).read_bytes()
                for name in ('selection.nwk', 'quality.tsv', 'references.nwk', 'references.tsv')}
    plan['request']['native_genomewide_parent']['metadata_sha256'] = {
        name: hashlib.sha256(raw).hexdigest() for name, raw in metadata.items()}
    r.require_same_key(frozen, r.frozen_prediction_cache_key(previous))
    r.require_same_key(frozen['plan_sha256'], r.plan_digest(previous, prior))
    r.refresh_source_versions(plan)
    root.mkdir(parents=True)
    for name, raw in metadata.items():
        (root / name).write_bytes(raw)
    r.atomic_json(root / 'plan.json', plan, immutable=True)
    return plan


def prepare(root, plan, cpus):
    """Import complete receipted stages, with full byte validation during copy."""
    if plan['request'].get('admission_cache'):
        import rescue_admission_replay
        rescue_admission_replay.prepare(root, plan, cpus)
    else:
        _prepare_raw(root, plan, cpus)
    completed = plan['request'].get('native_genomewide_worker_imports', {}).get('workers', {})
    if completed:
        import rescue_native_genomewide_reuse
        r = _r()
        with r.ReceiptSnapshot(lambda directory, key: r.verified(directory, key, hash_workers=cpus)) as deps:
            for name in completed:
                rescue_native_genomewide_reuse.import_worker(root, plan, name, cpus, dependencies=deps)


def _prepare_raw(root, plan, cpus):
    r = _r()
    parent = plan['request']['native_genomewide_parent']
    previous = Path(parent['root'])
    prior = r.load(previous, inspect_only=True)
    old_sha, new_sha = r.plan_digest(previous, prior), r.plan_digest(root, plan)
    r.require_same_key(parent['plan_sha256'], old_sha)
    expected = copy.deepcopy(prior)
    expected['output'] = str(root)
    expected['request']['tools'] = plan['request']['tools']
    if parent.get('compact_diagnostics'):
        expected['request']['output_storage'] = {**r.output_storage(prior), 'diagnostics': 'summary'}
    for key in ('prediction_cache', 'oriongeno_genomewide', 'native_genomewide_parent',
                'native_genomewide_worker_imports'):
        if key not in plan['request']:
            continue
        expected['request'][key] = plan['request'][key]
    # Generation hints are refreshed by the planner. The source SHA-256 values
    # in the inherited request remain part of this exact comparison.
    expected.pop('source_versions', None)
    r.require_same_key(expected, {k: v for k, v in plan.items() if k != 'source_versions'})
    r.require_same_key(CONTRACT, plan['request']['oriongeno_genomewide'])
    batch = r.RescueSourceBatch(plan)
    def guard():
        batch.check()
        r.require_same_key(old_sha, r.plan_digest(previous, prior))
        r.require_same_key(new_sha, r.plan_digest(root, plan))
    for name, sha in parent['metadata_sha256'].items():
        if batch.read([root / name, previous / name]) != {str(root / name): sha, str(previous / name): sha}:
            raise ValueError('Native extension metadata changed')
    for name in plan['species']:
        r.verify_sources(plan, [name], ['fasta', 'gff', 'genome'], verification=batch)
        relative = Path('prepared') / name
        if batch.read([previous / relative / 'receipt.json']) != {str(previous / relative / 'receipt.json'): parent['prepared'][name]}:
            raise ValueError('Native extension prepared generation changed')
        def build(tmp, relative=relative, name=name):
            r.copy_verified_comparison(previous / relative, tmp, {'plan': old_sha, 'species': name}, cpus)
        r.stage(root, relative, {'plan': new_sha, 'species': name}, build, guard, hash_workers=cpus)
    for job in plan['synteny_jobs']:
        relative = Path('synteny') / job['id']
        if batch.read([previous / relative / 'receipt.json']) != {str(previous / relative / 'receipt.json'): parent['synteny'][job['id']]}:
            raise ValueError('Native extension comparison generation changed')
        def build(tmp, relative=relative, job=job):
            r.copy_verified_comparison(previous / relative, tmp, r.comparison_key(previous, job), cpus)
        r.stage(root, relative, r.comparison_key(root, job), build, guard, hash_workers=cpus)

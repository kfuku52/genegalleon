"""Candidate-first native rescue outside the original coding annotation.

Protein search nominates real external donors, never a genomic placement or
orthology. All returned paths still require the ordinary genomic, terminal,
unanchored placement, ownership and coding-path admission checks.
"""
import copy
import hashlib
import json
import time
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


def nominate(frozen, name, proteins, genome, ownership, params, directory, *, cpus=1):
    """Search each eligible native isoform, retaining identical-protein aliases.

    There is no donor representation filter, candidate window or query cap. Search
    negatives are explicitly search negatives, not proof that no homolog exists.
    Only sequence-identical aliases share an alignment; loci/isoforms stay distinct.
    """
    if type(cpus) is not int or cpus < 1:
        raise ValueError('Native genomewide CPUs must be a positive integer')
    if not frozen or name not in frozen['predictions']:
        return [], {'contract': CONTRACT, 'native_paths': 0, 'paths': [], 'candidate_pairs': 0,
                    'native_prediction_available': False}
    native, r = _native(), _r()
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
            donor_aliases[protein].append((donor, query))
    if not donor_aliases:
        raise ValueError('Genomewide native rescue has no external donor proteins')
    queries = {f'n{i}': protein for i, protein in enumerate(sorted(groups))}
    subjects = {f'd{i}': protein for i, protein in enumerate(sorted(donor_aliases))}
    candidates = set()
    search_dir = Path(directory) / 'native_genomewide_search'
    search_dir.mkdir()
    if queries:
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
        with (search_dir / 'pairs.tsv').open() as stream:
            for line in stream:
                fields = line.split()
                if len(fields) != 2 or fields[0] not in queries or fields[1] not in subjects:
                    raise ValueError('Unknown genomewide protein search identifier')
                candidates.add(tuple(fields))
        # Exact short peptides need not seed a DIAMOND alignment. Preserve all
        # exact matches independently of the search program's seeding heuristic.
        by_sequence = {protein: token for token, protein in subjects.items()}
        candidates.update((token, by_sequence[protein]) for token, protein in queries.items() if protein in by_sequence)
    def comparisons():
        for q, s in sorted(candidates):
            yield q, s, queries[q], subjects[s]
    results = []
    for q, s, score in native._scored_pairs(comparisons(), cpus):
        alignment = {**score.as_dict(), 'junction_information_available': False,
                     'score_basis': 'protein_only_without_donor_junctions'}
        problems = (['oriongeno_alignment_bounded'] if score.bounded else [])
        if min(score.coverage_a, score.coverage_b) < params['minimum_coverage']:
            problems.append('low_reciprocal_coverage')
        for path in groups[queries[q]]:
            rows[path['native_transcript_id']]['status'] = 'evaluated'
            for donor, query in donor_aliases[subjects[s]]:
                evidence = {'target': name, 'donor': donor, 'query': query, 'genome_only': True,
                            'placement': 'unanchored', 'orthology': 'unassigned', 'expected_copy': 'unassigned',
                            'native_search_scope': CONTRACT['scope'], 'oriongeno': {
                                **prediction['identity'], 'native_transcript_id': path['native_transcript_id'],
                                'native_gene_id': path['native_gene_id'], 'alignment': alignment}}
                evidence['id'] = 'orion_genome_query_' + hashlib.sha256(
                    json.dumps([name, donor, query], separators=(',', ':')).encode()).hexdigest()[:20]
                results.append({'seqid': path['seqid'], 'strand': path['strand'], 'cds': copy.deepcopy(path['cds']),
                                'query': evidence['id'], 'coverage': score.coverage_b, 'identity': score.identity,
                                'frameshift': False, 'evidence': evidence, 'search': SEARCH,
                                'oriongeno_problems': list(problems)})
    native.verify_inputs(frozen, [name])
    audit.update(unique_query_proteins=len(queries), unique_donor_proteins=len(subjects),
                 external_donor_records=sum(map(len, donor_aliases.values())),
                 aligned_unique_protein_pairs=len(candidates), candidate_pairs=len(results),
                 search_and_alignment_seconds=time.perf_counter() - started)
    (search_dir / 'donors.dmnd').unlink(missing_ok=True)
    return results, audit


def finish_audit(audit, models):
    """Account for every native path after the actual shared admission decision."""
    outcomes = defaultdict(list)
    for model in models:
        if model.get('search') == SEARCH:
            outcomes[model['evidence']['oriongeno']['native_transcript_id']].append(model)
    for row in audit['paths']:
        if row['status'] != 'evaluated':
            continue
        records = outcomes[row['native_transcript_id']]
        if not records:
            raise ValueError('Evaluated native path disappeared before admission')
        admitted = [m for m in records if m['status'] in {'accepted', 'duplicate_support', 'accepted_alternative_path'}]
        row.update(status='admitted' if admitted else 'not_admitted',
                   candidate_records=len(records), admission_statuses=dict(Counter(m['status'] for m in records)),
                   reasons=sorted({p for m in records for p in m['problems']}),
                   external_donor_species=sorted({m['evidence']['donor'] for m in admitted}),
                   model_ids=sorted({m['model_id'] for m in admitted if m.get('model_id')}))
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
                                     proteins, genome, code, params['max_intron'], cpus)
    for model in predictions:
        region = evidence[model['query']]
        model['seqid'] = region['seqid']
        model['cds'] = [[a + region['start'], b + region['start'], p] for a, b, p in model['cds']]
        model['evidence'] = region
        model['search'] = 'native_genomewide_miniprot'
    return predictions


def build_plan(previous_output, output, *, producer_source=None, completed_output=None, completed_producer_source=None):
    """Freeze a new ordinary raw replay plus the explicit native search scope.

    Prepared donor proteins and synteny are imported as frozen scientific inputs,
    not asserted to be newly produced by the current formatter. Prediction-cache
    verification still checks their bytes, all source inputs and Miniprot identity.
    """
    r = _r()
    if producer_source:
        import rescue_admission_replay
        return rescue_admission_replay.build_plan(previous_output,producer_source,output,native_extension=True,
            completed_output=completed_output,completed_producer_source=completed_producer_source)
    if completed_output or completed_producer_source:
        raise ValueError('Completed native workers require the explicit validated-evidence producer source')
    previous, root = Path(previous_output).resolve(strict=True), Path(output).resolve()
    if root.exists() or root == previous or previous in root.parents or root in previous.parents:
        raise ValueError('Native genomewide extension requires a new disjoint output directory')
    prior = r.load(previous, inspect_only=True)
    if not prior['request'].get('oriongeno') or prior['request'].get('admission_cache'):
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
        'prepared': {n: r.digest(previous / 'prepared' / n / 'receipt.json') for n in prior['species']},
        'synteny': {j['id']: r.digest(previous / 'synteny' / j['id'] / 'receipt.json') for j in prior['synteny_jobs']}}
    metadata = {name: (previous / name).read_bytes()
                for name in ('selection.nwk', 'quality.tsv', 'references.nwk', 'references.tsv')}
    plan['request']['native_genomewide_parent']['metadata_sha256'] = {
        name: hashlib.sha256(raw).hexdigest() for name, raw in metadata.items()}
    r.require_same_key(frozen, r.frozen_prediction_cache_key(previous))
    r.require_same_key(frozen['plan_sha256'], r.plan_digest(previous, prior))
    root.mkdir(parents=True)
    for name, raw in metadata.items():
        (root / name).write_bytes(raw)
    r.atomic_json(root / 'plan.json', plan, immutable=True)
    return plan


def prepare(root, plan, cpus):
    """Import complete receipted stages, with full byte validation during copy."""
    r = _r()
    if plan['request'].get('admission_cache'):
        import rescue_admission_replay
        import rescue_native_genomewide_reuse
        rescue_admission_replay.prepare(root,plan,cpus)
        for name in plan['request'].get('native_genomewide_worker_imports',{}).get('workers',{}):
            rescue_native_genomewide_reuse.import_worker(root,plan,name,cpus)
        return
    parent = plan['request']['native_genomewide_parent']
    previous = Path(parent['root'])
    prior = r.load(previous, inspect_only=True)
    old_sha, new_sha = r.plan_digest(previous, prior), r.plan_digest(root, plan)
    r.require_same_key(parent['plan_sha256'], old_sha)
    expected = copy.deepcopy(prior)
    expected['output'] = str(root)
    expected['request']['tools'] = plan['request']['tools']
    for key in ('prediction_cache', 'oriongeno_genomewide', 'native_genomewide_parent'):
        expected['request'][key] = plan['request'][key]
    r.require_same_key(expected, plan)
    r.require_same_key(CONTRACT, plan['request']['oriongeno_genomewide'])
    batch = r.FreshDigestBatch()
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

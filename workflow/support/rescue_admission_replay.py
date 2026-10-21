"""Standard, source-proven admission replay of completed rescue evidence.

This is deliberately narrower than prediction-cache replay: searches, raw QC,
boundary alternatives and native donor/structure evidence must be unchanged.
Only placement, ownership and coding-path admission run again. Unknown producer
implementations or changed scientific settings require ordinary raw replay.
"""
import ast
import copy
import hashlib
import json
import shutil
import time
from collections import Counter
from pathlib import Path

try:
    import rescue_model_store as store
    from rescue_model_store_patch import patch_model_store
except ImportError:
    from . import rescue_model_store as store
    from .rescue_model_store_patch import patch_model_store

POST = frozenset({'overlap_existing_annotation', 'competing_new_models',
                  'ambiguous_coding_path_representative', 'ambiguous_annevo_coding_path_representative',
                  'ambiguous_oriongeno_coding_path_representative', 'outside_expected_synteny_interval',
                  'unanchored_genome_search'})
MUTABLE = frozenset({'status', 'model_id', 'parent_model_id', 'alternative_coding_paths', 'path_selection',
                     'locus_support', 'support', 'revision_owner_ids', 'revision_disposition',
                     'placement_evidence', 'quality_evidence', 'presence_evidence', 'integration_conflict',
                     'retained_coding_path'})
MODULES = {
    'rescue_raw_validation.py': 'raw_validation_implementation',
    'rescue_homology_policy.py': 'homology_problems_implementation',
    'rescue_path_compatibility.py': 'protein_path_support_implementation',
    'rescue_terminal_completion.py': 'complete_terminals_implementation',
    'rescue_boundary_revalidation.py': 'revalidate_boundaries_implementation',
    'rescue_additional_candidates.py': 'nominate_genome_only_candidates_implementation',
    'gene_model_selection.py': 'annevo_alignment_implementation',
    'gene_model_species_profiles.py': 'species_profiles_implementation',
    'rescue_model_quality.py': 'quality_implementation',
    'rescue_anchor_admission.py': 'anchor_admission_implementation',
    'cds_model_normalisation.py': 'cds_normalisation_implementation',
    'pairwise_synteny.py': 'mapping_implementation',
    'fasta_sequence_store.py': 'reader_implementation',
    'gene_model_catalog.py': 'genome_index_implementation',
    'gff_attribute_syntax.py': 'attribute_syntax_implementation',
    'rescue_prediction_cache.py': 'verify_prediction_cache_implementation',
    'rescue_search_inputs.py': 'search_inputs_implementation',
    'transcript_utrs.py': 'transcript_utr_implementation',
    'genome_sequence_identity.py': 'annevo_genome_identity_implementation',
    'format_species_annotation/common.py': 'ownership_attribute_parser_implementation',
}
FUNCTIONS = {
    'rescue_gene_models.py': ('validate_model', 'prepared', 'build_comparison', 'align_self', 'parse_anchors',
                              'run', 'candidate_intervals', 'local_flanks', 'candidates', 'iter_miniprot'),
    'annevo_gene_models.py': ('read_paths', '_score_protein_batch', '_scored_pairs', 'nominate',
                            'confirm_transferred_structure'),
    'rescue_coding_paths.py': ('coding_shape', 'coding_overlap', 'alignment_support', 'compatible_paths', 'original_ownership'),
    'input_generation_array_state.py': ('FreshDigestBatch', 'digest', 'digest_paths', 'atomic_json'),
}
ALLOWED_TOOL_CHANGES = frozenset({'implementation', 'model_store_implementation', 'resolve_coding_paths_implementation',
                                'annevo_candidate_implementation', 'state_implementation',
                                'comparison_copy_proof_implementation', 'output_hashes_implementation',
                                'admission_replay_implementation', 'admission_patch_implementation'})


def _r():
    try:
        import rescue_gene_models
    except ImportError:
        from . import rescue_gene_models
    return rescue_gene_models


def _definitions(path, names):
    definitions = {node.name: ast.dump(node, include_attributes=False)
                   for node in ast.parse(Path(path).read_text()).body if hasattr(node, 'name')}
    if set(names) - definitions.keys():
        raise ValueError('Missing rescue replay function: ' + str(path))
    return {name: definitions[name] for name in names}


def prove_contract(previous_plan, producer_source):
    """Bind unchanged strict science to actual archived producer source bytes."""
    r = _r()
    old = Path(producer_source).resolve(strict=True) / 'workflow/support'
    current = Path(r.__file__).resolve().parent
    before, after = previous_plan['request']['tools'], r.identities()
    bindings = {**MODULES, 'rescue_gene_models.py': 'implementation',
                'annevo_gene_models.py': 'annevo_candidate_implementation',
                'rescue_coding_paths.py': 'resolve_coding_paths_implementation',
                'input_generation_array_state.py': 'state_implementation'}
    hashes = {}
    batch = r.FreshDigestBatch()
    observed = batch.read([*[old / filename for filename in bindings],
                           *[current / filename for filename in MODULES.keys() | FUNCTIONS.keys()]])
    for filename, field in bindings.items():
        value = observed[str(old / filename)]
        if before.get(field) != value:
            raise ValueError('Producer source is not bound to its plan: ' + filename)
        hashes[str(old / filename)] = value
    for filename in MODULES:
        if hashes[str(old / filename)] != observed[str(current / filename)]:
            raise ValueError('Strict rescue evidence implementation changed: ' + filename)
    for filename, names in FUNCTIONS.items():
        if _definitions(old / filename, names) != _definitions(current / filename, names):
            raise ValueError('Rescue replay function changed: ' + filename)
    permitted = set(ALLOWED_TOOL_CHANGES)
    # This route is absent from legacy plans. Once used, its implementation is
    # strict evidence and cannot change during an admission-only replay.
    if not previous_plan['request'].get('oriongeno_genomewide'):
        permitted.add('native_genomewide_implementation')
    # The worker dispatcher was absent from older raw producers and is not
    # executed by admission replay. Once declared, retain its exact identity.
    if 'cached_validation_implementation' not in before:
        permitted.add('cached_validation_implementation')
    if 'genomewide_reuse_implementation' not in before:
        permitted.add('genomewide_reuse_implementation')
    if not previous_plan['request'].get('oriongeno'):
        permitted.update(('oriongeno_candidate_implementation', 'oriongeno_profiles_implementation'))
    if any(before.get(k) != after.get(k) for k in before.keys() | after.keys() if k not in permitted):
        raise ValueError('Rescue evidence runtime/tool identities changed')
    if previous_plan['request'].get('gemoma_jar'):
        raise ValueError('Admission replay does not yet qualify GeMoMa evidence')
    batch.check()
    return {'schema': 1, 'producer_source': str(old.parent.parent), 'source_files': hashes,
            'current_implementation': after['implementation'],
            'strict_qc_unchanged': True, 'native_alignment_and_structure_unchanged': True,
            'searches_reexecuted': False, 'native_prediction_reexecuted': False}


def _legacy_annotation_history(plan):
    """Keep a producer's separately declared annotation source, without QC reuse."""
    value = plan['request'].get('integration_reclassification')
    if value is None:
        inherited = plan['request'].get('admission_cache', {}).get('annotation_history')
        if inherited is None:
            return None
        value = {'source_root': inherited['root'], 'source_plan_sha256': inherited['plan_sha256'],
                 'source_worker_receipts': inherited['worker_receipts']}
    if not isinstance(value, dict):
        raise ValueError('Malformed prior annotation-source declaration')
    root, sha, workers = (value.get(k) for k in ('source_root', 'source_plan_sha256', 'source_worker_receipts'))
    def valid_sha(s):
        return isinstance(s, str) and len(s) == 64 and not set(s) - set('0123456789abcdef')
    if (not isinstance(root, str) or not Path(root).is_absolute() or '..' in Path(root).parts
            or not valid_sha(sha) or not isinstance(workers, dict) or set(workers) != set(plan['species'])
            or not all(valid_sha(s) for s in workers.values())):
        raise ValueError('Prior annotation source lacks complete frozen identities')
    return {'root': str(Path(root).resolve(strict=True)), 'plan_sha256': sha, 'worker_receipts': workers}


def build_plan(previous_output, producer_source, output, *, native_extension=False,
               completed_output=None, completed_producer_source=None):
    r = _r()
    previous, root = Path(previous_output).resolve(strict=True), Path(output).resolve()
    if root == previous or previous in root.parents or root in previous.parents or root.exists():
        raise ValueError('Admission replay requires a new, disjoint output directory')
    prior = r.load(previous, inspect_only=True)
    if native_extension:
        if prior['request'].get('oriongeno_genomewide'):
            raise ValueError('Native evidence extension requires the pre-genomewide producer')
        import rescue_native_genomewide_reuse
        proof=rescue_native_genomewide_reuse.prove_contract(prior,producer_source)
    else:
        if completed_output or completed_producer_source:
            raise ValueError('Completed imports require the explicit native extension route')
        proof = prove_contract(prior, producer_source)
    frozen = r.frozen_prediction_cache_key(previous)
    if set(frozen['species']) != set(prior['species']):
        raise ValueError('Admission replay requires all completed producer workers')
    if any(row.get('model_store', {}).get('format') != store.FORMAT for row in frozen['species'].values()):
        raise ValueError('Admission replay requires compact raw publications')
    plan = copy.deepcopy(prior)
    plan['output'] = str(root)
    plan['request']['tools'] = r.identities()
    plan['request']['prediction_cache'] = frozen
    plan['request']['admission_cache'] = {'root': str(previous), 'plan_sha256': frozen['plan_sha256'],
                                         'contract': proof, 'species': frozen['species'],
                                         'annotation_history': _legacy_annotation_history(prior)}
    plan['request'].pop('integration_reclassification', None)
    if native_extension:
        import rescue_native_genomewide
        plan['request']['oriongeno_genomewide']=rescue_native_genomewide.CONTRACT
        if bool(completed_output)!=bool(completed_producer_source):
            raise ValueError('Completed native workers require their archived producer source')
        if completed_output:
            plan['request']['native_genomewide_worker_imports']=rescue_native_genomewide_reuse.freeze_completed(
                plan,completed_output,completed_producer_source)
    metadata = {}
    for name in ('selection.nwk', 'quality.tsv', 'references.nwk', 'references.tsv'):
        metadata[name] = (previous / name).read_bytes()
    plan['request']['admission_cache']['metadata_sha256'] = {
        name: hashlib.sha256(raw).hexdigest() for name, raw in metadata.items()}
    r.require_same_key(frozen, r.frozen_prediction_cache_key(previous))
    if any(r.digest(previous / name) != sha for name, sha in plan['request']['admission_cache']['metadata_sha256'].items()):
        raise ValueError('Admission metadata changed during plan creation')
    root.mkdir(parents=True)
    for name, raw in metadata.items():
        (root / name).write_bytes(raw)
    r.atomic_json(root / 'plan.json', plan, immutable=True)
    return plan


def previous_context(root, plan):
    r = _r()
    cache = plan['request']['admission_cache']
    previous = Path(cache['root'])
    r.require_same_key(cache['species'], plan['request']['prediction_cache']['species'])
    prior = r.load(previous, inspect_only=True)
    r.require_same_key(cache['annotation_history'], _legacy_annotation_history(prior))
    if r.plan_digest(previous, prior) != cache['plan_sha256']:
        raise ValueError('Admission producer plan changed')
    native=cache['contract'].get('scope')=='native_genomewide_extension'
    if native:
        import rescue_native_genomewide
        import rescue_native_genomewide_reuse
        proof=rescue_native_genomewide_reuse.prove_contract(prior,cache['contract']['producer_source'])
        rescue_native_genomewide_reuse.verify_completed(plan)
    else:
        proof=prove_contract(prior,cache['contract']['producer_source'])
    r.require_same_key(cache['contract'],proof)
    expected = copy.deepcopy(prior)
    expected['output'] = str(root)
    expected['request']['tools'] = plan['request']['tools']
    expected['request']['prediction_cache'] = plan['request']['prediction_cache']
    expected['request']['admission_cache'] = cache
    expected['request'].pop('integration_reclassification', None)
    if native:
        expected['request']['oriongeno_genomewide']=rescue_native_genomewide.CONTRACT
        if 'native_genomewide_worker_imports' in plan['request']:
            expected['request']['native_genomewide_worker_imports']=plan['request']['native_genomewide_worker_imports']
    r.require_same_key(expected, plan)
    r.require_same_key(plan['request']['prediction_cache'], r.frozen_prediction_cache_key(previous))
    for name, sha in cache['metadata_sha256'].items():
        if r.digest(root / name) != sha or r.digest(previous / name) != sha:
            raise ValueError('Admission reference metadata changed')
    return previous, prior


def prepare(root, plan, cpus):
    """Import fully verified, semantically unchanged prepared/comparison stages."""
    r = _r()
    previous, prior = previous_context(root, plan)
    old_hash, new_hash = r.plan_digest(previous, prior), r.plan_digest(root, plan)
    inputs = r.FreshDigestBatch()
    def guard():
        inputs.check()
        r.require_same_key(old_hash, r.plan_digest(previous, prior))
        r.require_same_key(new_hash, r.plan_digest(root, plan))
    for name in plan['species']:
        r.verify_sources(plan, [name], ['fasta', 'gff', 'genome'], verification=inputs)
        old_key, new_key = {'plan': old_hash, 'species': name}, {'plan': new_hash, 'species': name}
        def build(tmp, name=name, old_key=old_key):
            r.copy_verified_comparison(previous / 'prepared' / name, tmp, old_key, cpus)
        r.stage(root, Path('prepared') / name, new_key, build, guard, hash_workers=cpus)
    for job in plan['synteny_jobs']:
        old_key, new_key = r.comparison_key(previous, job), r.comparison_key(root, job)
        def build(tmp, job=job, old_key=old_key):
            r.copy_verified_comparison(previous / 'synteny' / job['id'], tmp, old_key, cpus)
        r.stage(root, Path('synteny') / job['id'], new_key, build, guard, hash_workers=cpus)


def eligible(context):
    return not set(context.get('problems', ())) - POST


def reset_admission(model):
    result = copy.deepcopy({key: value for key, value in model.items() if key not in MUTABLE})
    result['problems'] = sorted(set(result.get('problems', ())) - POST)
    return result


def read_retained(root, plan, name, previous_worker, receipt_sha, batch, history):
    r = _r()
    retained = r.read_retained_coding_paths(previous_worker, receipt_sha, verification=batch, history_audit=history)
    extra = plan['request']['admission_cache']['annotation_history']
    if extra is None:
        return retained
    parent = Path(extra['root']) / 'rescued' / name
    expected = {str(Path(extra['root']) / 'plan.json'): extra['plan_sha256'],
                str(parent / 'receipt.json'): extra['worker_receipts'][name]}
    if batch.read(expected) != expected:
        raise ValueError('Separately bound annotation history changed')
    producer = json.loads((Path(extra['root']) / 'plan.json').read_text())
    before, current = producer['request'], plan['request']
    a, b = before['sources'][name], current['sources'][name]
    if before['files'][a['genome']] != current['files'][b['genome']] or a['genetic_code'] != b['genetic_code']:
        raise ValueError('Separately bound annotation genome/code differs')
    already = [row for row in history if row.get('status') == 'verified'
               and Path(row['directory']).resolve() == parent.resolve()]
    if already:
        if any(row['receipt_sha256'] != extra['worker_receipts'][name] for row in already):
            raise ValueError('Annotation history branch has conflicting identities')
        return retained
    inherited = r.read_retained_coding_paths(parent, extra['worker_receipts'][name], verification=batch,
                                             history_audit=history, _ancestors=(previous_worker.resolve(),))
    for shape, source in inherited.items():
        if shape in retained and retained[shape]['sequence_sha256'] != source['sequence_sha256']:
            raise ValueError('Previous annotations disagree on genomic CDS sequence')
        retained.setdefault(shape, source)
    batch.check()
    return retained


def _rewrite_table(source, destination, updates, appended, make_row, count, *, allow_tail=False):
    import csv
    with source.open(newline='') as inp, destination.open('x', newline='') as out:
        reader, writer = csv.reader(inp, delimiter='\t'), csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(next(reader))
        total = 0
        inserted = False
        for ordinal, row in enumerate(reader):
            if ordinal == count:
                writer.writerows(make_row(m) for m in appended)
                inserted = True
            if ordinal in updates:
                if row[0] != updates[ordinal]['query']:
                    raise ValueError('Admission audit order changed')
                row = make_row(updates[ordinal])
            writer.writerow(row)
            total += 1
        if total < count or (not allow_tail and total != count):
            raise ValueError('Admission audit/model count differs')
        if not inserted:
            writer.writerows(make_row(m) for m in appended)


def rescue(root, plan, name, cpus):
    r = _r()
    from gene_model_catalog import indexed_genome
    from input_generation_array_state import FreshDigestBatch
    previous, prior = previous_context(root, plan)
    source, params = plan['request']['sources'][name], r.parameters_for(plan['request'], name)
    old_worker = previous / 'rescued' / name
    frozen = plan['request']['admission_cache']['species'][name]
    if r.digest(old_worker / 'receipt.json') != frozen['receipt_sha256']:
        raise ValueError('Admission worker generation changed')
    old_receipt, observed_sha = store._receipt(old_worker)
    r.require_same_key(frozen['receipt_sha256'], observed_sha)
    r.require_same_key(old_receipt['key'], r.rescue_key(previous, prior, name))
    raw_key = store.frozen_model_store_key(old_worker)
    batch = FreshDigestBatch()
    expected = {str(old_worker / member): sha for member, sha in old_receipt['files'].items()}
    expected[str(old_worker / 'receipt.json')] = frozen['receipt_sha256']
    expected.update(plan['request']['admission_cache']['contract']['source_files'])
    key = r.rescue_key(root, plan, name)
    def guard():
        batch.check()
        r.require_same_key(key, r.rescue_key(root, plan, name))
        store._check_key(old_worker, raw_key, 'models')
        r.verify_sources(plan, [name], ['genome', 'fasta', 'gff'])
    def build(tmp):
        started = time.perf_counter()
        if batch.read(expected) != expected:
            raise ValueError('Admission worker is incomplete or corrupted')
        history = []
        retained = read_retained(root, plan, name, old_worker, frozen['receipt_sha256'], batch, history)
        counts = Counter()
        def select(context):
            counts['raw_records'] += 1
            return eligible(context)
        selected = dict(store.iter_model_items(old_worker, context_filter=select, frozen_key=raw_key, verification=batch))
        for ordinal, primary in store._iter(old_worker, 'accepted', verification=batch, with_ordinals=True):
            if ordinal not in selected or store._json(selected[ordinal]) != store._json(primary):
                raise ValueError('Prior primary is not an eligible, identical raw candidate')
        fingerprints = {i: hashlib.sha256(store._json(m)).digest() for i, m in selected.items()}
        partials = {i: copy.deepcopy(m.get('partial_evidence')) for i, m in selected.items()}
        current = {i: reset_admission(m) for i, m in selected.items()}
        selected.clear()
        if any(m.get('evidence', {}).get(engine) and not plan['request'].get(engine)
               for m in current.values() for engine in ('annevo', 'oriongeno')):
            raise ValueError('Native evidence is not enabled by its frozen producer plan')
        counts['admission_rechecked'] = len(current)
        counts['unchanged_strict_rejections'] = counts['raw_records'] - len(current)
        with r.open_text(Path(source['gff'])) as handle:
            existing = r.original_ownership(handle, r.attributes)
        proteins = {donor: {i: s for i, _, s in r.fasta_records(root / 'prepared' / donor / 'genes.pep')}
                    for donor in [name, *plan['donors'][name]]}
        native_audit = None
        with indexed_genome(source['genome']) as genome:
            extra, historical = r.revalidate_retained_coding_paths(retained, list(current.values()), genome,
                                                                  source['genetic_code'], params, proteins, name)
            models = [*current.values(), *extra]
            for model in models:
                r.check_interval(model)
            if plan['request']['admission_cache']['contract'].get('scope') == 'native_genomewide_extension':
                import rescue_native_genomewide_reuse
                additions, native_audit = rescue_native_genomewide_reuse.fresh_candidates(
                    root, plan, name, proteins, genome, r.OwnershipIndex(existing), params, tmp, cpus, models)
                extra.extend(additions)
                models.extend(additions)
            placement = r.reassess_unanchored_models(models, params['unanchored_min_species'],
                                                     ownership=r.OwnershipIndex(existing))
            r.consolidate(models, existing, name, plan.get('nearest_references', {}).get(name, ()),
                          retained_paths=retained)
            for model in models:
                if model['status'] in {'accepted', 'accepted_alternative_path'}:
                    fresh = r.validate_model(model, genome, source['genetic_code'], params)
                    if fresh['problems'] or fresh['cds'] != model['cds'] or fresh['sequence'] != model['sequence']:
                        raise ValueError('Admitted historical CDS failed fresh genomic QC')
                model['quality_evidence'] = r.model_quality(model, source['genetic_code'])
        if native_audit is not None:
            import rescue_native_genomewide
            rescue_native_genomewide.finish_audit(native_audit, models)
            r.atomic_json(tmp / 'oriongeno_genomewide.json', native_audit)
        if any(m.get('partial_evidence') != partials[i] for i, m in current.items()):
            raise ValueError('Admission replay changed partial evidence')
        revisions = {}
        for m in models:
            if len(m.get('revision_owner_ids', [])) == 1 and not set(m['problems']) - {'overlap_existing_annotation'}:
                shape = (*r.coding_shape(m), m['revision_owner_ids'][0])
                evidence = {**m['evidence'], 'alignment': {'coverage': m['coverage'], 'identity': m['identity'],
                                                          'problems': list(m['problems'])}}
                revisions.setdefault(shape, {**m, 'support': []})['support'].append(evidence)
        changed = {i: m for i, m in current.items() if hashlib.sha256(store._json(m)).digest() != fingerprints[i]}
        manifest = patch_model_store(old_worker, tmp, changed, appended=extra,
                                     revisions=revisions.values(), frozen_key=raw_key,
                                     append_partial_evidence=native_audit is not None)
        # Copy every other receipt member without truncating a shared inode.
        mutable = {'audit.tsv', 'quality_flags.tsv', 'search_inputs.json', 'placement_audit.json',
                   'retained_path_revalidation.json', 'retained_annotation_sources.json', 'admission_replay.json'}
        if native_audit is not None:
            mutable.add('oriongeno_genomewide.json')
        for member in old_receipt['files']:
            if member.startswith('model_store/') or member in mutable:
                continue
            source_path, target = store._member(old_worker, member), tmp / member
            target.parent.mkdir(parents=True, exist_ok=True)
            with source_path.open('rb') as inp, target.open('xb') as out:
                shutil.copyfileobj(inp, out)
            if r.digest(target) != old_receipt['files'][member]:
                raise ValueError('Admission sidecar changed during copying')
        def audit(m):
            return (m['query'], m['evidence']['query'], m['status'], ','.join(m['problems']),
                    m['coverage'], m['identity'], m.get('model_id', ''))
        def quality(m):
            q = m['quality_evidence']
            return (m['query'], m.get('model_id', ''), m['status'], q['start_codon'],
                    q['terminal_alignment']['n_aligned'], q['terminal_alignment']['c_aligned'],
                    ','.join(q['donor_species']), ','.join(q['flags']), q['query_alignment']['aligned_query_fraction'],
                    q['query_alignment']['internal_unaligned_query_fraction'])
        _rewrite_table(old_worker / 'audit.tsv', tmp / 'audit.tsv', current, extra, audit, counts['raw_records'], allow_tail=True)
        _rewrite_table(old_worker / 'quality_flags.tsv', tmp / 'quality_flags.tsv', current, extra, quality, counts['raw_records'])
        if 'search_inputs.json' in old_receipt['files']:
            search = json.loads((old_worker / 'search_inputs.json').read_text())
        elif counts['raw_records'] == 0:
            # Ordinary no-window workers publish an empty store without search
            # sidecars. Their native extension starts a real search below.
            search = {'schema': 1, 'empty_producer_search': True,
                      'prepared_receipts': {n: None for n in [name, *plan['donors'][name]]}}
        else:
            raise ValueError('Nonempty admission producer lacks search input identities')
        search.update(plan_sha256=r.plan_digest(root, plan), prediction_cache=plan['request']['prediction_cache'],
                      prepared_receipts={n: r.digest(root / 'prepared' / n / 'receipt.json') for n in search['prepared_receipts']})
        r.atomic_json(tmp / 'search_inputs.json', search)
        r.atomic_json(tmp / 'placement_audit.json', placement)
        r.atomic_json(tmp / 'retained_path_revalidation.json', historical)
        r.atomic_json(tmp / 'retained_annotation_sources.json', history)
        r.atomic_json(tmp / 'admission_replay.json', {'schema': 1, 'contract': plan['request']['admission_cache']['contract'],
                      'source_worker': str(old_worker), 'source_receipt_sha256': frozen['receipt_sha256'],
                      'counts': dict(counts), 'changed_records': len(changed), 'appended_records': len(extra),
                      'published_counts': manifest['counts'], 'full_stream_verified': True,
                      'admitted_cds_freshly_validated': True, 'wall_seconds': time.perf_counter() - started})
    return r.stage(root, Path('rescued') / name, key, build, guard, hash_workers=cpus)

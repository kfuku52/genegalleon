"""Actual protein search, genomic QC and shared admission on synthetic loci."""
import copy
import importlib
import json
import random
import shutil
import sys
from collections import Counter
from pathlib import Path

import pytest
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
native = importlib.import_module('rescue_native_genomewide')
annevo = importlib.import_module('annevo_gene_models')
rescue = importlib.import_module('rescue_gene_models')
refinement = importlib.import_module('gene_model_refinement')
catalog = importlib.import_module('gene_model_catalog')
native_fixture = importlib.import_module('test_annevo_gene_models').native_fixture
runtime_fixtures = importlib.import_module('test_rescue_gene_models_runtime')
hidden_models = runtime_fixtures.hidden_models


def fixture(tmp_path, *, strand='+', aliases=False, second_locus=False, coding_fault=None):
    genome, gff, receipt, table, sources, coding, blocks = native_fixture(tmp_path, strand=strand)
    codons = ['GCT', 'CGT', 'AAT', 'GAT', 'TGT', 'CAA', 'GAA', 'GGT', 'CAT', 'ATT',
              'CTG', 'AAA', 'ATG', 'TTT', 'CCT', 'TCT', 'ACT', 'TGG', 'TAT', 'GTT']
    rng = random.Random(19)
    coding = 'ATG' + ''.join(rng.choices(codons, k=100)) + 'TAA'
    dna = 'C' * 30 + coding[:120] + 'GT' + 'AC' * 12 + 'AG' + coding[120:] + 'C' * 30
    genome.write_text('>chr\n' + (str(Seq(dna).reverse_complement()) if strand == '-' else dna) + '\n')
    if coding_fault:
        text = genome.read_text()
        text = text.replace('GAA', 'TAA', 1) if coding_fault == 'stop' else text.replace('GTAC', 'CTAC', 1)
        genome.write_text(text)
    if aliases:
        rows = gff.read_text().splitlines()
        gff.write_text('\n'.join(rows + [line.replace('native.t', 'native.alt').replace('ID=c', 'ID=altc').replace('ID=e', 'ID=alte')
                                         for line in rows if '\tgene\t' not in line and not line.startswith('#')]) + '\n')
    if second_locus:
        genome.write_text(genome.read_text() + '>other\n' + genome.read_text().splitlines()[1] + '\n')
        gff.write_text(gff.read_text() + '\n'.join(
            line.replace('chr\t', 'other\t').replace('native', 'paralog').replace('ID=c', 'ID=pc').replace('ID=e', 'ID=pe')
            for line in gff.read_text().splitlines() if not line.startswith('#')) + '\n')
    value = json.loads(receipt.read_text())
    value.update(kind='genegalleon-oriongeno-native-prediction', model_profile='plants',
                 inference_settings={'lineage': 'plants'}, genome_sha256=annevo.digest(genome),
                 GFF_sha256=annevo.digest(gff), GFF_bytes=gff.stat().st_size)
    receipt.write_text(json.dumps(value))
    frozen = annevo.freeze_inputs(table, sources, engine='oriongeno')
    protein = str(Seq(coding).translate()).removesuffix('*')
    return genome, frozen, protein, blocks


def evaluate(tmp_path, *, owners=(), minimum_species=1, **kwargs):
    genome, frozen, protein, blocks = fixture(tmp_path, **kwargs)
    ownership = rescue.OwnershipIndex(list(owners))
    params = {**refinement.DEFAULTS, 'unanchored_min_species': minimum_species, 'terminal_max_extension': 300}
    proteins = {'Target_species': {'self': protein}, 'Donor_one': {'represented': protein, 'identical_alias': protein}}
    with catalog.indexed_genome(genome) as sequence:
        models, audit = native.nominate(frozen, 'Target_species', proteins, sequence, ownership, params, tmp_path, cpus=1)
        models = [rescue.check_interval(rescue.validate_model(m, sequence, 1, params)) for m in models]
        transfers = native.transfer_structure(models, proteins, sequence, 1, params, tmp_path, cpus=1)
        transfers = [rescue.check_interval(rescue.validate_model(m, sequence, 1, params)) for m in transfers]
        annevo.confirm_transferred_structure(models, transfers)
        models += transfers
        rescue.reassess_unanchored_models(models, minimum_species, ownership=ownership)
        models = rescue.consolidate(models, list(owners), 'Target_species', ['Donor_one'])
    return models, native.finish_audit(audit, models), blocks


@pytest.mark.parametrize('strand', ['+', '-'])
def test_outside_any_window_and_represented_donors_are_admitted_with_real_search(tmp_path, strand):
    models, audit, blocks = evaluate(tmp_path, strand=strand)
    accepted = [m for m in models if m['status'] == 'accepted' and m['search'] == native.SEARCH]
    assert len(accepted) == 1
    assert accepted[0]['cds'] == blocks and accepted[0]['problems'] == []
    assert accepted[0]['placement_evidence']['independent_donor_species'] == ['Donor_one']
    assert audit['path_status_counts'] == {'admitted': 1}
    assert audit['candidate_pairs'] == 2 and audit['aligned_unique_protein_pairs'] == 1
    assert not any(key.startswith('expected_') and key != 'expected_copy' for key in accepted[0]['evidence'])
    assert accepted[0]['evidence']['oriongeno']['alignment']['junction_information_available'] is False


@pytest.mark.parametrize('boundary', [False, True])
def test_coverage_prefilter_preserves_threshold_and_accounts_for_rejected_hits(tmp_path, monkeypatch, boundary):
    genome, frozen, protein, _ = fixture(tmp_path)
    donor = protein[:50] + '*'
    minimum = 50 / len(protein)
    if not boundary:
        minimum += 0.000001
    params = {**refinement.DEFAULTS, 'minimum_coverage': minimum}
    aligned = []
    run = rescue.run
    score = annevo._scored_pairs

    def one_hit(command, directory, label, *args, **kwargs):
        if label == 'blastp':
            (directory / 'pairs.tsv').write_text('n0\td0\n')
        else:
            return run(command, directory, label, *args, **kwargs)

    def observed(pairs, cpus):
        pairs = list(pairs)
        aligned.extend(pairs)
        return score(pairs, cpus)

    monkeypatch.setattr(rescue, 'run', one_hit)
    monkeypatch.setattr(annevo, '_scored_pairs', observed)
    with catalog.indexed_genome(genome) as sequence:
        models, audit = native.nominate(frozen, 'Target_species', {'Other_species': {'donor': donor}},
                                       sequence, rescue.OwnershipIndex([]), params, tmp_path)
    assert len(aligned) == len(models) == int(boundary)
    assert audit['length_rejected_protein_hits'] == int(not boundary)
    if not boundary:
        path, = native.finish_audit(audit, models)['paths']
        assert path['status'] == 'not_admitted'
        assert path['reasons'] == ['low_reciprocal_coverage']
        assert path['candidate_records'] == 0


def test_compact_pairs_keep_sorted_unique_ids_and_rejection_multiplicity(tmp_path):
    hits = tmp_path / 'pairs.tsv'
    query = {'n2': 'M', 'n10': 'MMM', 'n3': 'MAG'}
    subject = {'d2': 'MMM', 'd10': 'MAM', 'd3': 'M'}
    hits.write_text('n2 d2\nn10 d2\nn2 d2\nn10 d10\nn10 d2\n')
    rejected, pairs = native._candidate_protein_pairs(hits, query, subject, 0.7)
    assert rejected == {'n2': 2}
    assert list(pairs) == [('n10', 'd10', 'MMM', 'MAM'), ('n10', 'd2', 'MMM', 'MMM'),
                           ('n2', 'd3', 'M', 'M')]


@pytest.mark.parametrize('bad_hit', ['n00 d0', 'n0 d00', 'n0', 'n0 d0 d1'])
def test_compact_pairs_reject_malformed_or_unknown_hit_tokens(tmp_path, bad_hit):
    hits = tmp_path / 'pairs.tsv'
    hits.write_text(bad_hit + '\n')
    with pytest.raises(ValueError, match='Unknown genomewide protein search identifier'):
        native._candidate_protein_pairs(hits, {'n0': 'M'}, {'d0': 'M'}, 0.7)


def test_identical_isoforms_keep_both_native_ids_but_one_gene_model(tmp_path):
    models, audit, _ = evaluate(tmp_path, aliases=True)
    assert sum(m['status'] == 'accepted' for m in models) == 1
    assert audit['native_paths'] == 2 and audit['path_status_counts'] == {'admitted': 2}
    assert audit['native_gene_counts'] == {'admitted': 1}
    assert {m['evidence']['oriongeno']['native_transcript_id'] for m in models if m['search'] == native.SEARCH} == {'native.t', 'native.alt'}


def test_identical_protein_at_two_distinct_new_loci_does_not_invent_copy_assignment(tmp_path):
    models, audit, _ = evaluate(tmp_path, second_locus=True)
    assert not any(m['status'] == 'accepted' for m in models)
    assert audit['path_status_counts'] == {'not_admitted': 2}
    assert {m['placement_evidence']['reason'] for m in models} == {'insufficient_unique_donor_support'}


@pytest.mark.parametrize('fault,reason', [('stop', 'invalid_native_cds'), ('splice', 'not_admitted')])
def test_disrupted_native_genome_does_not_pass_homology_only(tmp_path, fault, reason):
    models, audit, _ = evaluate(tmp_path, coding_fault=fault)
    assert not any(m['status'] == 'accepted' for m in models)
    assert audit['path_status_counts'] == {reason: 1}
    if fault == 'splice':
        assert all('noncanonical_splice' in m['problems'] for m in models)


def test_original_cds_overlap_is_accounted_without_duplicating_gene(tmp_path):
    owner = {'gene_id': 'original', 'seqid': 'chr', 'strand': '+', 'start': 30, 'end': 364,
             'cds': [[30, 150], [178, 364]]}
    models, audit, _ = evaluate(tmp_path, owners=[owner])
    assert models == [] and audit['path_status_counts'] == {'original_cds_overlap': 1}
    assert audit['paths'][0]['original_owner_ids'] == ['original']


def test_donor_aliases_and_native_predictions_are_not_independent_species(tmp_path):
    models, audit, _ = evaluate(tmp_path, minimum_species=2)
    assert not any(m['status'] == 'accepted' for m in models)
    assert audit['path_status_counts'] == {'not_admitted': 1}


def test_validation_during_nomination_preserves_order_and_audit(tmp_path):
    genome, frozen, protein, _ = fixture(tmp_path, aliases=True, second_locus=True)
    params = {**refinement.DEFAULTS, 'terminal_max_extension': 300}
    proteins = {'Donor_one': {'p': protein, 'alias': protein}}
    baseline, streamed = (tmp_path / name for name in ('baseline', 'streamed'))
    baseline.mkdir()
    streamed.mkdir()
    with catalog.indexed_genome(genome) as sequence:
        def validate(model):
            return rescue.check_interval(rescue.validate_model(model, sequence, 1, params))

        ownership = rescue.OwnershipIndex([])
        raw, baseline_audit = native.nominate(frozen, 'Target_species', proteins, sequence,
                                              ownership, params, baseline)
        checked, streamed_audit = native.nominate(frozen, 'Target_species', proteins, sequence,
                                                   ownership, params, streamed, validate=validate)
        assert checked == [validate(model) for model in raw]
    for audit in (baseline_audit, streamed_audit):
        audit.pop('search_and_alignment_seconds')
    assert streamed_audit == baseline_audit
    assert streamed_audit['candidate_pairs'] == len(checked) == 8


def test_frozen_gff_change_fails_before_search(tmp_path):
    genome, frozen, protein, _ = fixture(tmp_path)
    Path(frozen['predictions']['Target_species']['gff']).write_text('changed')
    with catalog.indexed_genome(genome) as sequence, pytest.raises(ValueError, match='Frozen ANNEVO'):
        native.nominate(frozen, 'Target_species', {'Donor': {'p': protein}}, sequence,
                        rescue.OwnershipIndex([]), refinement.DEFAULTS, tmp_path)


def test_genomewide_native_cache_never_reuses_native_alignment_or_old_admission():
    cache = importlib.import_module('rescue_prediction_cache')
    model = {'query': 'q', 'search': native.SEARCH, 'evidence': {'oriongeno': {'native_transcript_id': 't'}}}
    assert cache._raw_prediction(copy.deepcopy(model), {'q': model['evidence']}) is None


@pytest.mark.parametrize('search,engine', [('annevo_nominated_path','annevo'),
    ('oriongeno_nominated_path','oriongeno'), ('oriongeno_genomewide_path','oriongeno')])
def test_native_raw_snapshot_precedes_derived_structure_confirmation(search, engine):
    cache = importlib.import_module('rescue_prediction_cache')
    model = dict(query='q',search=search,seqid='chr',strand='+',
        evidence={'donor':'Donor','query':'protein',engine:{'native_transcript_id':'t','alignment':{'identity':.9}}})
    model['raw_prediction'] = copy.deepcopy(model)
    model['evidence'][engine]['genomic_structure_transfer'] = {'confirmed':True,'scope':'whole_coding_path'}
    original = copy.deepcopy(model)
    assert cache._raw_prediction(model, {'q':{'id':'q'}}) is None
    assert model == original
    for key,value in [('donor','Other'),('query','other protein')]:
        changed = copy.deepcopy(model)
        changed['evidence'][key] = value
        with pytest.raises(ValueError,match='model provenance'):
            cache._raw_prediction(changed, {'q':{'id':'q'}})
    changed = copy.deepcopy(model)
    changed['evidence'][engine]['alignment']['identity'] = .1
    with pytest.raises(ValueError,match='model provenance'):
        cache._raw_prediction(changed, {'q':{'id':'q'}})


def test_miniprot_evidence_cannot_use_native_confirmation_exception():
    cache = importlib.import_module('rescue_prediction_cache')
    model = dict(query='q',search='genome_fallback',seqid='chr',strand='+',
        evidence={'donor':'Donor','query':'protein','oriongeno':{'native_transcript_id':'t'}})
    model['raw_prediction'] = copy.deepcopy(model)
    model['evidence']['oriongeno']['genomic_structure_transfer'] = {'confirmed':True}
    with pytest.raises(ValueError,match='model provenance'):
        cache._raw_prediction(model, {'q':{'id':'q'}})


def test_inventory_cannot_silently_lose_evaluated_path():
    with pytest.raises(ValueError, match='disappeared'):
        native.finish_audit({'native_paths': 1, 'paths': [{'status': 'evaluated', 'native_transcript_id': 't'}]}, [])


def test_missing_native_species_leaves_other_rescue_routes_unchanged(tmp_path):
    models, audit = native.nominate({'predictions': {}}, 'Not_in_native_subset', {}, None, None, {}, tmp_path)
    assert models == [] and audit['native_prediction_available'] is False
    assert native.finish_audit(audit, [])['path_status_counts'] == {}
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("strand,coding_fault", [("+", None), ("-", None), ("+", "splice")])
def test_donor_free_native_uses_genomic_qc_and_one_model_per_gene(tmp_path, strand, coding_fault):
    admission = importlib.import_module("rescue_oriongeno_admission")
    genome, frozen, _, blocks = fixture(tmp_path, strand=strand, aliases=True, coding_fault=coding_fault)
    models = []
    with catalog.indexed_genome(genome) as sequence:
        audit = admission.admit(frozen, "Target_species", sequence, 1, refinement.DEFAULTS,
                                rescue.OwnershipIndex([]), models)
        if coding_fault:
            assert models == []
            assert audit["genes"]["native"] == {"status": "invalid_cds", "reasons": ["noncanonical_splice"]}
            return
        assert len(models) == 1
        model = models[0]
        assert model["cds"] == blocks
        assert model["coverage"] is None and model["identity"] is None
        assert rescue.alignment_support(model) == []
        assert audit["genes"]["native"]["native_transcript_id"] == "native.alt"
        assert audit["genes"]["native"]["status"] == "admitted"
        # A previous rescue is a locus owner even when it has another ID.
        second = admission.admit(frozen, "Target_species", sequence, 1, refinement.DEFAULTS,
                                 rescue.OwnershipIndex([]), models)
        assert second["genes"]["native"]["status"] == "already_rescued"
        assert models == [model]


def test_completed_extension_reuses_old_searches_and_publishes_new_native_rescue(hidden_models, monkeypatch):
    directory, names, coding = hidden_models
    target = names[0]
    native_gff = directory / 'native.gff3'
    records = ['##gff-version 3']
    for index in (0, 8):
        a, b = index * (len(coding[index]) + 60), index * (len(coding[index]) + 60) + len(coding[index])
        records += [f'chr1\tOrionGeno\tgene\t{a+1}\t{b}\t.\t+\t.\tID=n{index}',
                    f'chr1\tOrionGeno\tmRNA\t{a+1}\t{b}\t.\t+\t.\tID=n{index}.t;Parent=n{index}',
                    f'chr1\tOrionGeno\tCDS\t{a+1}\t{b}\t.\t+\t0\tParent=n{index}.t']
    native_gff.write_text('\n'.join(records) + '\n')
    genome = directory / 'genome' / (target + '.genome.fa')
    receipt = directory / 'native.json'
    receipt.write_text(json.dumps(dict(schema_version=1, kind='genegalleon-oriongeno-native-prediction',
        native_returncode=0, gpu_probe={'gpu_verified': True}, model_profile='plants',
        inference_settings={'lineage': 'plants'}, runtime_digest='a'*64, model_sha256='b'*64,
        source_commit='c'*40, cache_key='d'*64, genome_sha256=rescue.digest(genome),
        GFF_sha256=rescue.digest(native_gff), GFF_bytes=native_gff.stat().st_size)))
    table = directory / 'native.tsv'
    table.write_text(f'species\tgff\treceipt\n{target}\t{native_gff}\t{receipt}\n')
    old = directory / 'previous'
    args = rescue.parser().parse_args(['plan', '--cds-dir', str(directory/'cds'), '--gff-dir', str(directory/'gff'),
        '--genome-dir', str(directory/'genome'), '--busco-dir', str(directory/'busco'), '--genome-fallback', '0',
        '--tree', str(directory/'tree.nwk'), '--output', str(old), '--oriongeno-inputs', str(table),
        '--diagnostics', 'legacy'])
    def build_empty_comparison(tmp, kind):
        rescue.atomic_json(tmp / 'blocks.json', [])
        if kind == 'pair':
            # A completed pair comparison retains its raw alignment even with no hits.
            (tmp / 'target.query.last').touch()

    def publish_workers(root):
        plan = rescue.build_plan(args)
        for name in names:
            rescue.prepared(root, plan, name)
        for job in plan['synteny_jobs']:
            rescue.stage(root, Path('synteny')/job['id'], rescue.comparison_key(root, job),
                         lambda tmp, kind=job['kind']: build_empty_comparison(tmp, kind))
        for name in names:
            rescue.rescue(root, plan, name, 1)
        return plan

    plan = publish_workers(old)
    store = importlib.import_module('rescue_model_store')
    assert not list(store.iter_accepted_models(old/'rescued'/target))
    configured = directory / 'configured'
    args.output, args.oriongeno_without_donor = configured, True
    configured_plan = publish_workers(configured)
    rescue.finalize(configured, configured_plan)
    accepted = list(store.iter_accepted_models(configured/'rescued'/target))
    assert len(accepted) == 1 and accepted[0]['sequence'] == coding[8]
    assert accepted[0]['quality_evidence']['donor_species'] == []
    audit = json.loads((configured/'rescued'/target/'oriongeno_without_donor.json').read_text())
    assert audit['genes']['n0']['status'] == 'original_cds_overlap'
    assert audit['genes']['n8']['status'] == 'admitted'
    gff = (configured/'augmented/species_gff'/f'{target}.rescue.gff3').read_text()
    assert 'rescue_admission=oriongeno_without_donor;locus_independent_donor_species_count=0' in gff
    report = importlib.import_module('gene_model_refinement_busco')
    methods = Counter()
    counts, per_locus = report.classify_rescue_support(
        accepted, target, {accepted[0]['model_id']}, configured_plan, method_intersections=methods)
    assert methods == {'oriongeno': 1}
    assert sum(counts.values()) == 0
    grouped = report.group_support_evidence(per_locus, target, configured_plan)
    assert grouped['oriongeno_without_donor'] == 1
    admission = importlib.import_module('rescue_oriongeno_admission')
    overlay = directory / 'overlay'
    with monkeypatch.context() as patch:
        def forbid_search(*args, **kwargs):
            raise AssertionError('OrionGeno additions must reuse the completed donor search')
        patch.setattr(rescue, 'run', forbid_search)
        admission.extend(old, overlay, 2)
        admission.extend(old, overlay, 2)
    assert list(store.iter_accepted_models(overlay/'rescued'/target)) == accepted
    assert (overlay/'augmented/species_gff'/f'{target}.rescue.gff3').read_text() == gff
    original = (old/'rescued'/target/'receipt.json').read_bytes()
    new = directory / 'extension'
    extended = native.build_plan(old, new)
    native.prepare(new, extended, 2)
    assert extended['request']['parameters'] == plan['request']['parameters']
    assert extended['request']['species_profiles'] == plan['request']['species_profiles']
    assert extended['donors'] == plan['donors']
    run = rescue.run
    def guarded(command, root, label, stdout=None):
        assert label not in {'miniprot_index', 'miniprot_genome'}
        return run(command, root, label, stdout)
    monkeypatch.setattr(rescue, 'run', guarded)
    with monkeypatch.context() as interrupted:
        def stop_after_search(*args):
            raise ValueError('Simulated interruption after DIAMOND')

        interrupted.setattr(native, '_candidate_protein_pairs', stop_after_search)
        with pytest.raises(ValueError, match='Simulated interruption after DIAMOND'):
            rescue.rescue(new, extended, target, 2)
    miniprot_checkpoint = new / 'miniprot_predictions' / target / 'receipt.json'
    diamond_checkpoint = new / 'oriongeno_protein_search' / target / 'receipt.json'
    frozen_searches = {path: path.read_bytes() for path in (miniprot_checkpoint, diamond_checkpoint)}
    with monkeypatch.context() as retry:
        def reuse_search(command, root, label, stdout=None):
            assert label not in {'makedb', 'blastp'}, 'Retry repeated the completed native protein search'
            return guarded(command, root, label, stdout)

        retry.setattr(rescue, 'run', reuse_search)
        rescue.rescue(new, extended, target, 2)
    assert {path: path.read_bytes() for path in frozen_searches} == frozen_searches
    accepted = list(store.iter_accepted_models(new/'rescued'/target))
    assert len(accepted) == 1 and accepted[0]['sequence'] == coding[8]
    assert (old/'rescued'/target/'receipt.json').read_bytes() == original
    audit = json.loads((new/'rescued'/target/'oriongeno_genomewide.json').read_text())
    assert audit['path_status_counts'] == {'admitted': 1, 'original_cds_overlap': 1}
    assert rescue.verified(new/'rescued'/target, rescue.rescue_key(new, extended, target))
    # A fresh invocation reads the ordinary verified publication idempotently.
    rescue.rescue(new, rescue.load(new), target, 2)
    # Validated evidence import must produce the same complete native decisions,
    # while running fresh search/mapping and never replaying the old raw stream.
    source = directory / 'source'
    shutil.copytree(Path(rescue.__file__).parent, source / 'workflow/support')
    fast = directory / 'validated_extension'
    fast_plan = native.build_plan(old, fast, compact_diagnostics=True)
    assert fast_plan['request']['admission_cache']['contract']['producer_source'] == str(source.resolve())
    native.prepare(fast, fast_plan, 2)
    def forbidden(*args, **kwargs):
        raise AssertionError('validated evidence must not repeat old raw QC')
    with monkeypatch.context() as patch:
        patch.setattr(rescue, 'verify_prediction_cache', forbidden)
        rescue.rescue(fast, fast_plan, target, 2)
    assert list(store.iter_models(fast/'rescued'/target)) == list(store.iter_models(new/'rescued'/target))
    fast_audit = json.loads((fast/'rescued'/target/'oriongeno_genomewide.json').read_text())
    assert fast_audit['path_status_counts'] == audit['path_status_counts']
    fast_worker = fast / 'rescued' / target
    assert not (fast_worker / 'quality_flags.tsv').exists()
    detailed = rescue.export_worker_diagnostics(fast, fast_plan, target, directory / 'detailed', 2)
    for filename in ('audit.tsv', 'quality_flags.tsv'):
        assert (detailed / filename).read_text() == (new / 'rescued' / target / filename).read_text()

    reuse = importlib.import_module('rescue_native_genomewide_reuse')
    prove = reuse.prove_contract
    with monkeypatch.context() as patch:
        def changed_parser(plan, source):
            return {**prove(plan, source), 'native_parser_equivalence_required': True}

        def compare_parser(plan, name, genome):
            if name == target:
                raise reuse.IncompatibleEvidence('Frozen native coding paths changed under the current GFF parser')
            return {'native_prediction_available': False}

        patch.setattr(reuse, 'prove_contract', changed_parser)
        patch.setattr(reuse, 'parser_equivalence', compare_parser)
        mixed = directory / 'mixed_extension'
        mixed_plan = native.build_plan(old, mixed)
        native.prepare(mixed, mixed_plan, 2)
        rescue.rescue(mixed, mixed_plan, target, 2)
        rescue.rescue(mixed, mixed_plan, names[1], 2)
        validation = json.loads((mixed / 'rescued' / target / 'cached_prediction_validation.json').read_text())
        assert validation['reuse'] == 'raw_prediction_cache'
        assert 'coding paths changed' in validation['reason']
        assert (mixed / 'rescued' / names[1] / 'admission_replay.json').is_file()
        assert list(store.iter_accepted_models(mixed / 'rescued' / target)) == accepted

    archived_file = source / 'workflow/support/rescue_raw_validation.py'
    archived_file.write_text(archived_file.read_text() + '\n# incompatible archived identity\n')
    fallback = native.build_plan(old, directory / 'automatic_raw_extension')
    assert 'producer source differs' in fallback['request']['native_genomewide_parent']['validated_reuse_unavailable']
    archived_file.write_bytes((Path(rescue.__file__).parent / archived_file.name).read_bytes())

    # A new qualified generation may independently import a completed worker,
    # retaining its real producer receipt and performing no fresh science.
    imported = directory / 'imported_extension'
    imported_plan = native.build_plan(old, imported, completed_output=new)
    frozen = imported_plan['request']['native_genomewide_worker_imports']['workers']
    assert set(frozen) == {target}
    # An unrelated array task can finish after the snapshot was taken.
    rescue.rescue(new, extended, names[1], 2)
    native.prepare(imported, imported_plan, 2)
    assert not (imported / 'rescued' / names[1]).exists()
    monkeypatch.setattr(native, 'nominate', forbidden)
    rescue.rescue(imported, imported_plan, target, 2)
    assert list(store.iter_models(imported/'rescued'/target)) == list(store.iter_models(new/'rescued'/target))
    provenance = json.loads((imported/'rescued'/target/'completed_native_worker_import.json').read_text())
    assert provenance['source_receipt_sha256'] == rescue.digest(new/'rescued'/target/'receipt.json')
    assert not provenance['new_scientific_computation']

    worker_receipt = new / 'rescued' / target / 'receipt.json'
    worker_receipt.write_bytes(worker_receipt.read_bytes() + b'\n')
    with pytest.raises(ValueError, match='worker generations changed'):
        reuse.verify_completed(imported_plan)

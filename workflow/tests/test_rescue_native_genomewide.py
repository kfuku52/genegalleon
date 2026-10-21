"""Actual protein search, genomic QC and shared admission on synthetic loci."""
import copy
import importlib
import json
import random
import sys
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


def test_inventory_cannot_silently_lose_evaluated_path():
    with pytest.raises(ValueError, match='disappeared'):
        native.finish_audit({'native_paths': 1, 'paths': [{'status': 'evaluated', 'native_transcript_id': 't'}]}, [])


def test_missing_native_species_leaves_other_rescue_routes_unchanged(tmp_path):
    models, audit = native.nominate({'predictions': {}}, 'Not_in_native_subset', {}, None, None, {}, tmp_path)
    assert models == [] and audit['native_prediction_available'] is False
    assert native.finish_audit(audit, [])['path_status_counts'] == {}
    assert list(tmp_path.iterdir()) == []


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
        '--genome-dir', str(directory/'genome'), '--busco-dir', str(directory/'busco'),
        '--tree', str(directory/'tree.nwk'), '--output', str(old), '--oriongeno-inputs', str(table)])
    plan = rescue.build_plan(args)
    for name in names:
        rescue.prepared(old, plan, name)
    for job in plan['synteny_jobs']:
        rescue.stage(old, Path('synteny')/job['id'], rescue.comparison_key(old, job),
                     lambda tmp: rescue.atomic_json(tmp/'blocks.json', []))
    # This fixture has no original nominated windows. Native rescue therefore
    # supplies a real missing locus without relying on another search's hit.
    monkeypatch.setattr(rescue, 'candidates', lambda *_: [])
    monkeypatch.setattr(rescue, 'nominate_genome_only_candidates', lambda *_: ([], {}))
    for name in names:
        rescue.rescue(old, plan, name, 1)
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
    rescue.rescue(new, extended, target, 2)
    store = importlib.import_module('rescue_model_store')
    accepted = list(store.iter_accepted_models(new/'rescued'/target))
    assert len(accepted) == 1 and accepted[0]['sequence'] == coding[8]
    assert (old/'rescued'/target/'receipt.json').read_bytes() == original
    audit = json.loads((new/'rescued'/target/'oriongeno_genomewide.json').read_text())
    assert audit['path_status_counts'] == {'admitted': 1, 'original_cds_overlap': 1}
    assert rescue.verified(new/'rescued'/target, rescue.rescue_key(new, extended, target))
    # A fresh invocation reads the ordinary verified publication idempotently.
    rescue.rescue(new, rescue.load(new), target, 2)

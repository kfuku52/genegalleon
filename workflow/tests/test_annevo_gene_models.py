"""Conservative admission/immutability tests with explicitly synthetic native paths."""
import copy
import importlib
import json
import sys
from pathlib import Path

import pytest
from Bio.Seq import Seq

SUPPORT = Path(__file__).resolve().parents[1] / 'support'
sys.path.insert(0, str(SUPPORT))
annevo = importlib.import_module('annevo_gene_models')
catalog_module = importlib.import_module('gene_model_catalog')
rescue = importlib.import_module('rescue_gene_models')
refinement = importlib.import_module('gene_model_refinement')


def native_fixture(tmp_path, *, strand='+'):
    coding = 'ATG' + 'GAA' * 100 + 'TAA'
    dna = 'C' * 30 + coding[:120] + 'GT' + 'AC' * 12 + 'AG' + coding[120:] + 'C' * 30
    blocks = [[30, 150, 0], [178, 364, 0]]
    if strand == '-':
        size = len(dna)
        dna = str(Seq(dna).reverse_complement())
        blocks = [[size - b, size - a, phase] for a, b, phase in blocks]
    genome = tmp_path / 'genome.fa'
    genome.write_text('>chr\n' + dna + '\n')
    start, end = min(b[0] for b in blocks), max(b[1] for b in blocks)
    lines = [f'chr\tANNEVO\tgene\t{start + 1}\t{end}\t.\t{strand}\t.\tID=native',
             f'chr\tANNEVO\tmRNA\t{start + 1}\t{end}\t.\t{strand}\t.\tID=native.t;Parent=native']
    for i, (a, b, phase) in enumerate(blocks):
        lines.extend([f'chr\tANNEVO\tCDS\t{a + 1}\t{b}\t.\t{strand}\t{phase}\tID=c{i};Parent=native.t',
                      f'chr\tANNEVO\texon\t{a + 1}\t{b}\t.\t{strand}\t.\tID=e{i};Parent=native.t'])
    gff = tmp_path / 'native.gff3'
    gff.write_text('##gff-version 3\n' + '\n'.join(lines) + '\n')
    receipt = tmp_path / 'native.json'
    value = dict(schema_version=1, kind='genegalleon-annevo-native-prediction', native_returncode=0,
                 gpu_probe={'gpu_verified': True}, model_profile='Magnoliopsida', runtime_digest='a' * 64,
                 inference_settings={'lineage': 'Magnoliopsida'},
                 model_sha256='b' * 64, source_commit='c' * 40, cache_key='d' * 64,
                 genome_sha256=annevo.digest(genome), GFF_sha256=annevo.digest(gff), GFF_bytes=gff.stat().st_size)
    receipt.write_text(json.dumps(value))
    table = tmp_path / 'inputs.tsv'
    table.write_text('species\tgff\treceipt\nTarget_species\tnative.gff3\tnative.json\n')
    sources = {'Target_species': {'genome': str(genome), 'genetic_code': 1}}
    return genome, gff, receipt, table, sources, coding, blocks


@pytest.mark.parametrize('strand', ['+', '-'])
def test_actual_native_path_coordinates_and_real_donor_alignment(tmp_path, strand):
    genome, _gff, _receipt, table, sources, coding, blocks = native_fixture(tmp_path, strand=strand)
    frozen = annevo.freeze_inputs(table, sources)
    regions = [dict(id='q' + donor, donor=donor, query='protein', seqid='chr', start=0, end=394)
               for donor in ('Donor_one', 'Donor_two')]
    protein = str(Seq(coding).translate()).removesuffix('*')
    proteins = {r['donor']: {'protein': protein} for r in regions}
    with catalog_module.indexed_genome(genome) as sequence:
        models, audit = annevo.nominate(frozen, 'Target_species', regions, proteins, sequence, refinement.DEFAULTS)
        checked = [rescue.validate_model(m, sequence, 1, refinement.DEFAULTS) for m in models]
    assert len(models) == 2 and audit['native_paths'] == 1 and audit['proposals'] == []
    assert all(m['sequence'] == coding and m['cds'] == blocks and m['problems'] == [] for m in checked)
    assert all(m['coverage'] == m['identity'] == 1 for m in models)
    assert all(m['evidence']['annevo']['native_transcript_id'] == 'native.t' for m in models)


@pytest.mark.parametrize('profile', annevo.PROFILES)
def test_all_qualified_model_lineages_freeze_without_changing_candidate_gates(tmp_path, profile):
    _genome, _gff, receipt, table, sources, *_ = native_fixture(tmp_path)
    document = json.loads(receipt.read_text())
    document.update(model_profile=profile, inference_settings={'lineage': profile})
    receipt.write_text(json.dumps(document))
    frozen = annevo.freeze_inputs(table, sources)
    assert frozen['predictions']['Target_species']['identity']['model_profile'] == profile
    assert frozen['predictions']['Target_species']['identity']['model_sha256'] == document['model_sha256']


def test_legacy_magnoliopsida_receipt_remains_compatible(tmp_path):
    _genome, _gff, receipt, table, sources, *_ = native_fixture(tmp_path)
    document = json.loads(receipt.read_text())
    document.pop('inference_settings')
    receipt.write_text(json.dumps(document))
    assert annevo.freeze_inputs(table, sources)['predictions']['Target_species']['identity']['model_profile'] == 'Magnoliopsida'


@pytest.mark.parametrize('profile,settings', [('Mammalia', None), ('Actinopteri', {'lineage': 'Magnoliopsida'}),
                                             ('Magnoliopsida', {'lineage': 'Fungi'})])
def test_mismatched_or_missing_new_inference_lineage_is_rejected(tmp_path, profile, settings):
    _genome, _gff, receipt, table, sources, *_ = native_fixture(tmp_path)
    document = json.loads(receipt.read_text())
    document.update(model_profile=profile, inference_settings=settings)
    receipt.write_text(json.dumps(document))
    with pytest.raises(ValueError, match='inference lineage'):
        annevo.freeze_inputs(table, sources)


@pytest.mark.parametrize('field,value', [('kind', 'other'), ('native_returncode', 1), ('model_profile', 'unsupported'),
                                       ('source_commit', 'main'), ('runtime_digest', 'invalid'),
                                       ('genome_sha256', 'f' * 64), ('GFF_sha256', 'f' * 64), ('GFF_bytes', 0),
                                       ('gpu_probe', {'gpu_verified': False})])
def test_unqualified_or_mismatched_native_receipt_rejected(tmp_path, field, value):
    _genome, _gff, receipt, table, sources, *_ = native_fixture(tmp_path)
    document = json.loads(receipt.read_text())
    document[field] = value
    receipt.write_text(json.dumps(document))
    with pytest.raises(ValueError, match='ANNEVO'):
        annevo.freeze_inputs(table, sources)


@pytest.mark.parametrize('role', ['gff', 'receipt', 'table'])
def test_frozen_native_input_changes_fail_at_fresh_boundaries(tmp_path, role):
    _genome, gff, receipt, table, sources, *_ = native_fixture(tmp_path)
    frozen = annevo.freeze_inputs(table, sources)
    path = {'gff': gff, 'receipt': receipt, 'table': table}[role]
    path.write_text(path.read_text() + '\n')
    with pytest.raises(ValueError, match='Frozen ANNEVO'):
        annevo.verify_inputs(frozen, ['Target_species'])


@pytest.mark.parametrize('row', ['Unknown_species\tnative.gff3\tnative.json\n',
                               'Target_species\tnative.gff3\tnative.json\nTarget_species\tnative.gff3\tnative.json\n'])
def test_unknown_or_duplicate_species_fail(tmp_path, row):
    _genome, _gff, _receipt, table, sources, *_ = native_fixture(tmp_path)
    table.write_text('species\tgff\treceipt\n' + row)
    with pytest.raises(ValueError, match='Unknown or duplicate'):
        annevo.freeze_inputs(table, sources)


def test_no_native_input_leaves_prediction_boundary_empty(tmp_path):
    assert annevo.nominate(None, 'Target_species', [], {}, None, {}) == (
        [], {'native_paths': 0, 'candidate_pairs': 0, 'proposals': []})
    assert list(tmp_path.iterdir()) == []


def test_reciprocal_coverage_guards_gene_fusions_and_bounded_alignment(tmp_path):
    genome, _gff, _receipt, table, sources, coding, _blocks = native_fixture(tmp_path)
    frozen = annevo.freeze_inputs(table, sources)
    protein = str(Seq(coding).translate()).removesuffix('*')
    region = dict(id='query', donor='Donor_species', query='short', seqid='chr', start=0, end=394)
    with catalog_module.indexed_genome(genome) as sequence:
        models, _audit = annevo.nominate(frozen, 'Target_species', [region],
                                         {'Donor_species': {'short': protein[:40]}}, sequence, refinement.DEFAULTS)
    assert models[0]['coverage'] == 1  # Donor-only coverage would falsely pass.
    assert 'low_reciprocal_coverage' in models[0]['annevo_problems']


def test_whole_genome_placement_does_not_erase_independent_support_requirement(tmp_path):
    genome, _gff, _receipt, table, sources, coding, _blocks = native_fixture(tmp_path)
    frozen = annevo.freeze_inputs(table, sources)
    protein = str(Seq(coding).translate()).removesuffix('*')
    region = dict(id='query', donor='Donor_species', query='protein', seqid='chr', start=0, end=394,
                  genome_only=True, annevo_bounded_placement=True)
    with catalog_module.indexed_genome(genome) as sequence:
        models, _audit = annevo.nominate(frozen, 'Target_species', [region],
                                         {'Donor_species': {'protein': protein}}, sequence, refinement.DEFAULTS)
        model = rescue.check_interval(rescue.validate_model(models[0], sequence, 1, refinement.DEFAULTS))
    assert model['evidence']['genome_only'] is True
    rescue.reassess_unanchored_models([model], 2)
    assert model['problems']  # ANNEVO prediction itself is not another donor.


def test_native_model_outside_a_nominated_window_remains_audited(tmp_path):
    genome, _gff, _receipt, table, sources, *_ = native_fixture(tmp_path)
    with catalog_module.indexed_genome(genome) as sequence:
        models, audit = annevo.nominate(annevo.freeze_inputs(table, sources), 'Target_species', [], {}, sequence, refinement.DEFAULTS)
    assert not models and audit['proposals'] == [{'native_transcript_id': 'native.t', 'reason': 'no_bounded_donor_nomination'}]


@pytest.mark.parametrize('replacement', [('Parent=native.t', 'Parent=missing'), ('ID=e0', 'ID=c0'),
                                       ('\t+\t0\tID=c0', '\t+\t.\tID=c0')])
def test_invalid_native_gff_graph_is_not_silently_repaired(tmp_path, replacement):
    genome, gff, _receipt, _table, _sources, *_ = native_fixture(tmp_path)
    gff.write_text(gff.read_text().replace(*replacement))
    with catalog_module.indexed_genome(genome) as sequence, pytest.raises(ValueError, match='ANNEVO'):
        annevo.read_paths(gff, sequence)


def test_duplicate_predictors_do_not_count_as_independent_species():
    fixture = importlib.import_module('test_gene_model_refinement')
    catalog, models, edges = fixture.classification_fixture()
    native = copy.deepcopy(models[0])
    native['evidence'] = {'annevo': {'native_transcript_id': 'native.t'}}
    result = fixture.classify(catalog, [models[0], native], edges)[0]
    assert result['status'] == 'proposal' and result['donors'] == ['Species_donor1']
    assert 'insufficient_independent_support' in result['problems']


def test_weak_annevo_donor_cannot_veto_two_independently_supported_donors():
    fixture = importlib.import_module('test_gene_model_refinement')
    catalog, models, edges = fixture.classification_fixture()
    native = copy.deepcopy(models[0])
    native['problems'] = ['low_reciprocal_coverage']
    native['evidence'] = {'annevo': {'native_transcript_id': 'native.t'}}
    result = fixture.classify(catalog, models + [native], edges)[0]
    assert result['status'] == 'accepted' and len(result['donors']) == 2
    assert result['alignments'][-1]['supports_path'] is False


def test_scoped_override_is_forwarded_by_input_generation():
    import subprocess
    result = subprocess.run(['bash', '-c', 'source "$1"; gg_print_entrypoint_config_vars gg_input_generation_entrypoint.sh',
                             'vars', str(SUPPORT / 'gg_entrypoint_config_vars.sh')], capture_output=True, text=True, check=True)
    assert 'gene_model_annevo_inputs' in result.stdout.splitlines()
    core = (SUPPORT.parent / 'core/gg_input_generation_core.sh').read_text()
    assert 'rescue_args+=(--annevo-inputs' in core and 'refinement_args+=(--annevo-inputs' in core


def test_fasta_serialisation_changes_require_sequence_bound_receipt(tmp_path):
    import gzip
    genome, _gff, receipt, table, sources, *_ = native_fixture(tmp_path)
    original = genome.read_text()
    compressed = tmp_path / 'normalised.fa.gz'
    compressed.write_bytes(gzip.compress(('>chr different description\n' + original.splitlines()[1].lower() + '\n').encode()))
    sources['Target_species']['genome'] = str(compressed)
    with pytest.raises(ValueError, match='candidate genome differs'):
        annevo.freeze_inputs(table, sources)  # Legacy receipts still require exact bytes.
    value = json.loads(receipt.read_text())
    value['genome_sequence_sha256'] = annevo.genome_sequence_digest(genome)
    receipt.write_text(json.dumps(value))
    frozen = annevo.freeze_inputs(table, sources)
    assert frozen['predictions']['Target_species']['identity']['genome_match'] == 'contig_sequence_identity'
    assert annevo.genome_sequence_digest(compressed) == value['genome_sequence_sha256']
    compressed.write_bytes(gzip.compress(original.replace('>chr', '>another_contig').encode()))
    with pytest.raises(ValueError, match='candidate genome differs'):
        annevo.freeze_inputs(table, sources)


@pytest.mark.parametrize('text', ['>a\nATG\n>a\nATG\n', '>a\n', 'ATG\n', '>a\nAXG\n'])
def test_invalid_semantic_genome_identity_fails(tmp_path, text):
    path = tmp_path / 'invalid.fa'
    path.write_text(text)
    with pytest.raises(ValueError):
        annevo.genome_sequence_digest(path)


def test_native_cached_paths_are_skipped_without_weakening_raw_provenance_checks():
    cache = importlib.import_module('rescue_prediction_cache')
    region = {'id': 'q', 'donor': 'Donor_species', 'query': 'protein'}
    native = {'query': 'q', 'search': 'annevo_nominated_path', 'seqid': 'chr', 'strand': '+',
              'evidence': {**region, 'annevo': {'native_transcript_id': 'native.t'}}}
    native['raw_prediction'] = copy.deepcopy(native)
    assert cache._raw_prediction(native, {'q': region}) is None
    native['raw_prediction']['seqid'] = 'other'
    with pytest.raises(ValueError, match='does not match its model provenance'):
        cache._raw_prediction(native, {'q': region})

"""Conservative admission/immutability tests with explicitly synthetic native paths."""
import copy
import importlib
import json
import subprocess
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


def test_package_import_works_without_support_on_python_path():
    root = SUPPORT.parents[1]
    code = 'import sys;sys.path.insert(0,' + repr(str(root)) + ');import workflow.support.annevo_gene_models'
    result = subprocess.run([sys.executable, '-I', '-c', code], capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr


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


@pytest.mark.parametrize('engine', ['annevo', 'oriongeno'])
@pytest.mark.parametrize('strand', ['+', '-'])
def test_actual_native_path_coordinates_and_real_donor_alignment(tmp_path, strand, engine):
    genome, _gff, _receipt, table, sources, coding, blocks = native_fixture(tmp_path, strand=strand)
    if engine == 'oriongeno':
        value = json.loads(_receipt.read_text())
        value.update(kind='genegalleon-oriongeno-native-prediction', model_profile='plants', inference_settings={'lineage': 'plants'})
        _receipt.write_text(json.dumps(value))
    frozen = annevo.freeze_inputs(table, sources, engine=engine)
    regions = [dict(id='q' + donor, donor=donor, query='protein', seqid='chr', start=0, end=394)
               for donor in ('Donor_one', 'Donor_two')]
    protein = str(Seq(coding).translate()).removesuffix('*')
    proteins = {r['donor']: {'protein': protein} for r in regions}
    with catalog_module.indexed_genome(genome) as sequence:
        models, audit = annevo.nominate(frozen, 'Target_species', regions, proteins, sequence, refinement.DEFAULTS, engine=engine)
        checked = [rescue.validate_model(m, sequence, 1, refinement.DEFAULTS) for m in models]
    assert len(models) == 2 and audit['native_paths'] == 1 and audit['proposals'] == []
    assert all(m['sequence'] == coding and m['cds'] == blocks and m['problems'] == [] for m in checked)
    assert all(m['coverage'] == m['identity'] == 1 for m in models)
    assert all(m['evidence'][engine]['native_transcript_id'] == 'native.t' for m in models)
    assert all(m['evidence'][engine]['alignment']['junction_information_available'] is False for m in models)
    assert all('score' in m['evidence'][engine]['alignment'] and 'junction_similarity' in m['evidence'][engine]['alignment'] for m in models)


@pytest.mark.parametrize('strand', ['+', '-'])
def test_parallel_nomination_retains_exact_order_scores_audit_and_inputs(tmp_path, strand):
    genome, _gff, _receipt, table, sources, coding, _blocks = native_fixture(tmp_path, strand=strand)
    frozen = annevo.freeze_inputs(table, sources)
    regions = [dict(id=f'q{i}', donor='Donor_one' if i % 2 else 'Donor_two', query='protein',
                    seqid='chr', start=0, end=394) for i in range(600)]
    regions += [dict(id='unbounded', donor='Donor_one', query='protein', seqid='chr',
                     start=0, end=394, genome_only=True),
                dict(id='no_hit', donor='Donor_two', query='protein', seqid='chr', start=0, end=20)]
    protein = str(Seq(coding).translate()).removesuffix('*')
    proteins = {'Donor_one': {'protein': protein}, 'Donor_two': {'protein': 'M' + 'K' * 100}}
    before = copy.deepcopy((frozen, regions, proteins))
    with catalog_module.indexed_genome(genome) as sequence:
        serial = annevo.nominate(frozen, 'Target_species', regions, proteins, sequence, refinement.DEFAULTS)
        parallel = annevo.nominate(frozen, 'Target_species', regions, proteins, sequence, refinement.DEFAULTS, cpus=2)
    assert parallel == serial
    assert [row['query'] for row in parallel[0]] == [f'q{i}' for i in range(600)]
    assert (frozen, regions, proteins) == before
    parallel[0][0]['cds'][0][0] += 1
    assert parallel[0][0]['cds'] != parallel[0][1]['cds']


@pytest.mark.parametrize('cpus', [0, -1, True, 1.5, '2'])
def test_nomination_rejects_invalid_cpu_budget(cpus):
    with pytest.raises(ValueError, match='CPUs must be a positive integer'):
        annevo.nominate(None, 'Target_species', [], {}, None, {}, cpus=cpus)


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


@pytest.mark.parametrize('code', [1, 6, 12])
def test_candidate_export_cli_checks_the_actual_genetic_code(tmp_path, code):
    genome, gff, receipt, _table, sources, *_ = native_fixture(tmp_path)
    output = tmp_path / 'export.tsv'
    if code != 1:
        output.write_text('existing candidate export')
    result = subprocess.run([sys.executable, str(SUPPORT / 'annevo_gene_models.py'),
        '--species', 'Target_species', '--genome', str(genome), '--gff', str(gff),
        '--receipt', str(receipt), '--genetic-code', str(code), '--output', str(output)],
        capture_output=True, text=True, timeout=30)
    if code == 1:
        assert result.returncode == 0, result.stderr
        assert annevo.freeze_inputs(output, sources)['predictions']['Target_species']['identity']['model_profile'] == 'Magnoliopsida'
    else:
        assert result.returncode != 0 and 'genetic code 1' in result.stderr
        assert output.read_text() == 'existing candidate export'
        assert not output.with_suffix('.tsv.provenance.json').exists()


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


@pytest.mark.parametrize('strand', ['+', '-'])
def test_transfer_confirmation_needs_actual_same_query_and_junctions(tmp_path, strand):
    genome, _gff, _receipt, table, sources, coding, _blocks = native_fixture(tmp_path, strand=strand)
    region = dict(id='q', donor='Donor_species', query='protein', seqid='chr', start=0, end=394)
    proteins = {'Donor_species': {'protein': str(Seq(coding).translate()).removesuffix('*')}}
    with catalog_module.indexed_genome(genome) as sequence:
        models, _ = annevo.nominate(annevo.freeze_inputs(table, sources), 'Target_species', [region], proteins, sequence, refinement.DEFAULTS)
        checked = rescue.validate_model(models[0], sequence, 1, refinement.DEFAULTS)
    transferred = {**copy.deepcopy(checked), 'evidence': region, 'search': 'synteny_interval'}
    positive = copy.deepcopy(checked)
    annevo.confirm_transferred_structure([positive], [transferred])
    assert positive['evidence']['annevo']['genomic_structure_transfer']['confirmed']
    for bad in [copy.deepcopy(checked), {**copy.deepcopy(transferred), 'evidence': {**region, 'query': 'other'}},
                {**copy.deepcopy(transferred), 'problems': ['frameshift']},
                {**copy.deepcopy(transferred), 'cds': [[30, 364, 0]]}]:
        negative = copy.deepcopy(checked)
        annevo.confirm_transferred_structure([negative], [bad])
        assert 'annevo_unconfirmed_coding_structure' in negative['problems']


@pytest.mark.parametrize('strand', ['+', '-'])
def test_source_structure_retains_exons_and_complete_termini(strand):
    from gene_model_source_retention import source_structure_retention
    blocks = [[100, 130, 0], [200, 233, 0], [300, 330, 0]]
    if strand == '-':
        blocks.reverse()
    original = dict(candidate_id='source', seqid='chr', strand=strand, blocks=blocks,
                    cds='ATG' + 'AAA' * 29 + 'TAA', quality={})
    candidate = dict(seqid='chr', strand=strand, blocks=[*blocks])
    assert source_structure_retention(candidate, [original], 1)['preserved']
    # A frame-changing missing exon does not invalidate the retained endpoints.
    missing = {**original, 'blocks': [blocks[0], blocks[-1]], 'cds': 'ATG' + 'AAA' * 17 + 'A' + 'TAA'}
    assert source_structure_retention(candidate, [missing], 1)['preserved']
    for altered in ([[100, 130, 0], [300, 330, 0]],
                    [[103, 130, 0], [200, 233, 0], [300, 330, 0]],
                    [[97, 130, 0], [200, 233, 0], [300, 330, 0]],
                    [[100, 130, 0], [200, 236, 0], [300, 330, 0]]):
        assert not source_structure_retention({**candidate, 'blocks': altered}, [original], 1)['preserved']


def test_incomplete_terminal_extension_remains_possible():
    from gene_model_source_retention import source_structure_retention
    original = dict(candidate_id='source', seqid='chr', strand='+', blocks=[[100, 130, 0]],
                    cds='ATG' + 'AAA' * 9, quality={})
    assert source_structure_retention(dict(seqid='chr', strand='+', blocks=[[100, 160, 0]]), [original], 1)['preserved']


def test_new_locus_needs_whole_path_transfer_and_cannot_borrow_partial_termini(tmp_path):
    genome, _gff, _receipt, table, sources, coding, _blocks = native_fixture(tmp_path)
    region = dict(id='q', target='Target_species', donor='Donor_species', query='protein', seqid='chr', start=0, end=394,
                  genome_only=True, annevo_bounded_placement=True)
    proteins = {'Donor_species': {'protein': str(Seq(coding).translate()).removesuffix('*')}}
    with catalog_module.indexed_genome(genome) as sequence:
        models, _ = annevo.nominate(annevo.freeze_inputs(table, sources), 'Target_species', [region], proteins, sequence, refinement.DEFAULTS)
        native = rescue.check_interval(rescue.validate_model(models[0], sequence, 1, refinement.DEFAULTS))
    transferred = {**copy.deepcopy(native), 'evidence': region, 'search': 'genome_fallback'}
    good = copy.deepcopy(native)
    annevo.confirm_transferred_structure([good], [transferred])
    assert good['evidence']['annevo']['genomic_structure_transfer']['confirmed']
    partial = copy.deepcopy(transferred)
    partial['cds'][0][0] += 3
    partial['problems'].append('missing_start')
    proposed = copy.deepcopy(native)
    annevo.confirm_transferred_structure([proposed], [partial])
    assert 'annevo_unconfirmed_coding_structure' in proposed['problems']


def test_homology_cannot_discard_retained_source_exon_without_whole_path_rna():
    fixture = importlib.import_module('test_gene_model_refinement')
    catalog, models, edges = fixture.classification_fixture()
    original = catalog['loci'][0]['candidates'][0]
    original['blocks'].append([200, 240, 0])
    original['cds'] += 'AAA' * 13 + 'A'
    row = fixture.classify(catalog, models, edges)[0]
    assert row['status'] == 'proposal' and 'unconfirmed_source_structure_change' in row['problems']
    assert not row['candidate']['quality']['representative_eligible']
    rna = [dict(species=catalog['species'], seqid='chr1', strand='+', cds_blocks=json.dumps([b[:2] for b in models[0]['cds']]),
                transcript_id='whole_rna', count=1)]
    admitted = fixture.classify(catalog, models, edges, rna=rna)[0]
    assert admitted['status'] == 'accepted' and admitted['candidate']['quality']['representative_eligible']


def test_whole_path_rna_can_confirm_structure_without_miniprot_junction_agreement():
    fixture = importlib.import_module('test_gene_model_refinement')
    catalog, models, edges = fixture.classification_fixture()
    for model in models:
        model['problems'] = ['annevo_unconfirmed_coding_structure']
        model['evidence'] = {'annevo': {'native_transcript_id': 'native.t'}}
    assert fixture.classify(catalog, models, edges)[0]['status'] == 'proposal'
    rna = [dict(species=catalog['species'], seqid='chr1', strand='+', cds_blocks=json.dumps([b[:2] for b in models[0]['cds']]),
                transcript_id='actual_whole_RNA_fixture', count=1)]
    row = fixture.classify(catalog, models, edges, rna=rna)[0]
    assert row['status'] == 'accepted'
    assert all(a['structure_confirmation'] == 'whole_coding_RNA' for a in row['alignments'])


def test_mask_sampling_is_reproducible_across_process_hash_seeds():
    import os
    import subprocess
    benchmark_dir = SUPPORT.parent / 'benchmarks'
    code = '''import json,sys
sys.path.insert(0,sys.argv[1])
import benchmark_annevo_masked_models as b
ids=['gene_'+str(i) for i in range(80)]
choices={g:[(100,131),(200,232),(300,334)] for g in reversed(ids)}
genes,exons=b.sampled_masks(dict.fromkeys(ids),choices,.25,20261014)
print(json.dumps([sorted(genes),{g:sorted(v) for g,v in exons.items()}],sort_keys=True))
'''
    results = [subprocess.check_output([sys.executable, '-c', code, str(benchmark_dir)], text=True,
                                      env={**os.environ, 'PYTHONHASHSEED': value}) for value in ('1', '79')]
    assert results[0] == results[1]
    genes, exons = json.loads(results[0])
    assert len(genes) == len(exons) == 20


def test_benchmark_separates_nonexact_hidden_locus_from_off_reference_addition(tmp_path):
    sys.path.insert(0, str(SUPPORT.parent / 'benchmarks'))
    benchmark = importlib.import_module('benchmark_annevo_masked_models')
    fixture = importlib.import_module('test_rescue_coding_paths')
    models = []
    for blocks in ([[0, 63, 0]], [[100, 160, 0]]):
        row = fixture.model('toy', blocks, ('Donor_one', 'Donor_two'))
        row['evidence']['target'] = 'Target_species'
        models.append(row)
    gff = tmp_path / 'masked.gff3'
    gff.write_text('##gff-version 3\n')
    truth = {'reference_gene': {'seqid': 'chr1', 'strand': '+', 'blocks': [[0, 60, 0]]}}
    result = benchmark.evaluate_rescue(models, gff, truth, {'reference_gene'}, 'Target_species')
    assert result['accepted_additions'] == 2
    assert result['nonexact_additions_at_hidden_loci'] == result['off_reference_additions_unknown'] == 1


def genome_benchmark():
    sys.path.insert(0, str(SUPPORT.parent / 'benchmarks'))
    return importlib.import_module('benchmark_annevo_genome_models')


def test_genome_benchmark_matches_all_isoforms_and_separates_unknown_loci():
    b = genome_benchmark()
    first = dict(seqid='chr', strand='+', blocks=[[10, 20, 0], [40, 51, 2]])
    second = {**first, 'blocks': [[10, 23, 0], [40, 51, 2]]}
    other = dict(seqid='chr', strand='+', blocks=[[80, 110, 0]])
    catalog = {'loci': [dict(gene_id='gene', candidates=[first, second]), dict(gene_id='other', candidates=[other])]}
    reference = b.reference_index(catalog)
    assert b.match_category(second, *reference)[:2] == ('exact_any_reference_isoform', ['gene'])
    unknown = dict(seqid='chr', strand='+', cds=[[130, 160, 0]])
    assert b.match_category(unknown, *reference) == ('off_reference_unknown', [], [])
    fusion = dict(seqid='chr', strand='+', cds=[[10, 20, 0], [80, 110, 0]])
    assert b.match_category(fusion, *reference)[0] == 'multiple_reference_loci'
    # Exonic overlap is required: an intron-nested gene has no coding overlap.
    nested = dict(seqid='chr', strand='+', cds=[[25, 35, 0]])
    assert b.match_category(nested, *reference)[0] == 'off_reference_unknown'


def external_experiment_fixture(tmp_path):
    b = genome_benchmark()
    genome, native = tmp_path / 'genome.fa', tmp_path / 'native.gff3'
    genome.write_text('>chr\nATGAAATAA\n')
    native.write_text('##gff-version 3\n')
    inventory = dict(root='/sealed', source_commit='a' * 40,
                     files=[dict(path='model.pt', sha256='b' * 64)])
    ipath = tmp_path / 'runtime_inventory.json'
    ipath.write_text(json.dumps(inventory))
    manifest = dict(schema_version=1, kind='predictor-native-package-v1', tool='ANNEVO', species='Target_species',
                    inputs={'genome': '/sealed/genome.fa'}, options={'lineage': 'Magnoliopsida'}, root='/output',
                    runtime={'inventory': '/sealed/inventory.json', 'source_commit': 'a' * 40, 'model': '/sealed/model.pt'},
                    files={'/sealed/genome.fa': b.rescue.digest(genome), '/sealed/inventory.json': b.rescue.digest(ipath)})
    mpath = tmp_path / 'manifest.json'
    mpath.write_text(json.dumps(manifest))
    receipt = dict(manifest_sha256=b.rescue.digest(mpath), upstream_returncode=0, integration_returncode=0,
                   gpu_probe={'gpu_verified': True, 'tool': 'ANNEVO'}, upstream_completion_reported=True,
                   artifact_sha256=b.rescue.digest(native), artifact_bytes=native.stat().st_size,
                   command=['annotation.py', '-g', '/sealed/genome.fa', '-m', '/sealed/model.pt',
                            '-l', 'Magnoliopsida', '-o', '/output/prediction.gff3'])
    rpath = tmp_path / 'upstream_exit.json'
    rpath.write_text(json.dumps(receipt))
    return b, genome, native, rpath, receipt


def test_external_benchmark_keeps_original_provenance_separate_from_native_gg_receipts(tmp_path):
    b, genome, native, rpath, receipt = external_experiment_fixture(tmp_path)
    frozen = b.external_prediction(tmp_path, 'Target_species', genome, native)
    identity = frozen['predictions']['Target_species']['identity']
    assert identity['provenance_kind'] == 'verified_external_native_experiment'
    assert 'cache_key' not in identity
    with pytest.raises(ValueError, match='successful.*native'):
        annevo.prediction_identity(native, rpath, b.rescue.digest(genome))
    assert json.loads(rpath.read_text()) == receipt


@pytest.mark.parametrize('change', ['gpu', 'exit', 'manifest', 'model', 'genome', 'gff'])
def test_external_benchmark_rejects_failed_or_changed_original_evidence(tmp_path, change):
    b, genome, native, rpath, receipt = external_experiment_fixture(tmp_path)
    if change == 'gpu':
        receipt['gpu_probe']['gpu_verified'] = False
    elif change == 'exit':
        receipt['upstream_returncode'] = False  # bool cannot masquerade as integer zero
    elif change == 'manifest':
        receipt['manifest_sha256'] = '0' * 64
    elif change == 'model':
        receipt['command'][receipt['command'].index('-m') + 1] = '/other/model.pt'
    elif change == 'genome':
        genome.write_text('>chr\nATGCCCTAA\n')
    else:
        native.write_text('##gff-version 3\n# changed\n')
    rpath.write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match='binding|sealed inputs'):
        b.external_prediction(tmp_path, 'Target_species', genome, native)


@pytest.mark.parametrize('strand', ['+', '-'])
def test_completed_actual_donor_transfer_confirms_native_path_without_waiving_placement(tmp_path, strand):
    genome, _gff, _receipt, table, sources, coding, blocks = native_fixture(tmp_path, strand=strand)
    region = dict(id='q', target='Target_species', donor='Donor_species', query='protein',
                  seqid='chr', start=0, end=394, genome_only=True, annevo_bounded_placement=True)
    protein = str(Seq(coding).translate()).removesuffix('*')
    proteins = {'Donor_species': {'protein': protein}}
    with catalog_module.indexed_genome(genome) as sequence:
        raw_native, _ = annevo.nominate(annevo.freeze_inputs(table, sources), 'Target_species', [region], proteins,
                                        sequence, refinement.DEFAULTS)
        native = rescue.check_interval(rescue.validate_model(raw_native[0], sequence, 1, refinement.DEFAULTS))
        partial_blocks = copy.deepcopy(blocks)
        if strand == '+':
            partial_blocks[0][0] += 3
        else:
            partial_blocks[0][1] -= 3
        raw = dict(seqid='chr', strand=strand, cds=partial_blocks, query='q', evidence=region, frameshift=False,
                   coverage=(len(protein) - 1) / len(protein), identity=1., query_start=1,
                   query_end=len(protein), query_length=len(protein), search='genome_fallback')
        partial = rescue.validate_model(raw, sequence, 1, refinement.DEFAULTS)
        before = copy.deepcopy(native)
        annevo.confirm_transferred_structure([before], [rescue.check_interval(copy.deepcopy(partial))])
        assert 'annevo_unconfirmed_coding_structure' in before['problems']
        completed = rescue.complete_terminals(partial, sequence, 1, refinement.DEFAULTS, protein, rescue.validate_model)
        completed = rescue.check_interval(completed)
        annevo.confirm_transferred_structure([native], [completed])
    assert completed['terminal_completion']['status'] == 'completed'
    assert completed['cds'] == blocks
    assert completed['problems'] == native['problems'] == ['unanchored_genome_search']
    assert native['evidence']['annevo']['genomic_structure_transfer']['scope'] == 'whole_coding_path'
    # Genomic confirmation cannot resolve copy/placement ownership on its own.
    assert native['evidence']['genome_only'] is True


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('phase', [-1, None])
def test_unspecified_phase_does_not_erase_observed_source_initiation(strand, phase):
    from gene_model_source_retention import source_structure_retention
    blocks = [[100, 130, phase], [200, 232, -1], [300, 330, -1]]
    if strand == '-':
        blocks.reverse()
        blocks[0][2] = phase
    # Missing internal coding sequence prevents whole-ORF phase inference,
    # but the retained physical start and stop codons still constrain repair.
    source = dict(candidate_id='source', seqid='chr', strand=strand, blocks=[blocks[0], blocks[-1]],
                  cds='ATG' + 'AAA' * 18 + 'A' + 'TAA', quality={'translation_offset': 0})
    exact = dict(seqid='chr', strand=strand, blocks=copy.deepcopy(blocks))
    assert source_structure_retention(exact, [source], 1)['preserved']
    extended = copy.deepcopy(exact)
    if strand == '+':
        extended['blocks'][0][0] -= 12
    else:
        extended['blocks'][0][1] += 12
    audit = source_structure_retention(extended, [source], 1)
    assert not audit['preserved'] and not audit['sources'][0]['complete_start_preserved']
    # Explicit partial starts or a known nonzero first phase keep their
    # existing incomplete-terminus behavior; no complete ORF is inferred here.
    for partial in [{**source, 'quality': {'annotated_partial': True}},
                    {**source, 'blocks': [[*source['blocks'][0][:2], 1], source['blocks'][-1]]}]:
        assert source_structure_retention(extended, [partial], 1)['preserved']

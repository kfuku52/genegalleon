"""Selected CDSs stay fixed while transcript ends supply projected UTRs."""
import importlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
utrs = importlib.import_module('transcript_utrs')


class Genome:
    references = ('chr',)

    def __init__(self, sequence='ACGT'*300):
        self.sequence = sequence

    def fetch(self, seqid, start, end):
        assert seqid == 'chr'
        return self.sequence[start:end]

    def get_reference_length(self, _seqid):
        return len(self.sequence)


def model(tid='t', gene='g', cds=((100, 400),), exons=((50, 450),), strand='+'):
    def row(kind, a, b, attrs, phase='.'):
        return f'chr\tsource\t{kind}\t{a+1}\t{b}\t.\t{strand}\t{phase}\t{attrs}\n'
    text = row('gene', exons[0][0], exons[-1][1], 'ID='+gene)
    text += row('mRNA', exons[0][0], exons[-1][1], f'ID={tid};Parent={gene}')
    for i, (a, b) in enumerate(exons):
        text += row('exon', a, b, f'ID={tid}.e{i};Parent={tid}')
    for i, (a, b) in enumerate(cds):
        text += row('CDS', a, b, f'ID={tid}.c{i};Parent={tid}', '0')
    return text


@pytest.mark.parametrize('strand', ['+', '-'])
def test_boundary_shifts_choose_longest_independently_and_preserve_cds(strand):
    selected = model(cds=((150, 350),), exons=((150, 350),), strand=strand)
    original = model(cds=((100, 400),), exons=((30, 420),), strand=strand)
    native = utrs.read_models(model('native', 'ng', cds=((90, 390),), exons=((50, 480),), strand=strand))
    result, audit, sequences = utrs.complete(selected, original, Genome(), 'Species_test', native)
    assert [r for r in result.splitlines() if '\tCDS\t' in r] == [r for r in selected.splitlines() if '\tCDS\t' in r]
    low, high = ('utr5', 'utr3') if strand == '+' else ('utr3', 'utr5')
    assert audit[0][low]['source'] == 'existing' and audit[0][low]['length'] == 120
    assert audit[0][high]['source'] == 'oriongeno' and audit[0][high]['length'] == 130
    expected = Genome().sequence[30:480]
    assert sequences['t'] == (str(Seq(expected).reverse_complement()) if strand == '-' else expected)
    reconstructed = utrs.read_models(result)['t']
    assert reconstructed['exons'] == [[30, 480]]
    assert reconstructed['cds'] == [(150, 350, 0)]


def test_spliced_length_ranking_does_not_count_introns():
    selected = model(cds=((300, 400),), exons=((300, 400),))
    original = model(cds=((300, 400),), exons=((0, 20), (280, 400)))
    native = utrs.read_models(model('n', 'ng', cds=((300, 400),), exons=((240, 400),)))
    _, audit, _ = utrs.complete(selected, original, Genome(), 'S', native)
    assert audit[0]['utr5']['source'] == 'oriongeno'
    assert audit[0]['utr5']['length'] == 60


@pytest.mark.parametrize('strand', ['+', '-'])
def test_cds_extending_past_one_outer_end_only_removes_that_side(strand):
    selected = model(cds=((20, 350),), exons=((20, 350),), strand=strand)
    native = utrs.read_models(model('n', 'ng', strand=strand))
    _, audit, _ = utrs.complete(selected, selected, Genome(), 'S', native)
    low, high = ('utr5', 'utr3') if strand == '+' else ('utr3', 'utr5')
    assert audit[0][low]['length'] == 0
    assert audit[0][high]['length'] == 100


def test_internal_coding_exon_differences_do_not_prevent_utrs():
    selected = model(cds=((100, 160), (220, 400)), exons=((100, 160), (220, 400)))
    native = utrs.read_models(model('n', 'ng', cds=((100, 180), (250, 400)), exons=((40, 180), (250, 470))))
    result, audit, _ = utrs.complete(selected, selected, Genome(), 'S', native)
    assert audit[0]['utr5']['length'] == 60 and audit[0]['utr3']['length'] == 70
    assert utrs.read_models(result)['t']['cds'] == [(100, 160, 0), (220, 400, 0)]


def test_utrs_outside_declared_exons_are_not_silently_ignored():
    text = model() + 'chr\tsource\tfive_prime_UTR\t11\t20\t.\t+\t.\tParent=t\n'
    with pytest.raises(ValueError, match='UTR is not contained'):
        utrs.read_models(text)
    assert utrs.read_models(text, strict=False)['t']['problem'] == 'utr_outside_exons'


@pytest.mark.parametrize('strand', ['+', '-'])
def test_mislabeled_source_utr_is_replaced_with_correct_side(strand):
    text = model(strand=strand) + f'chr\tsource\tfive_prime_UTR\t401\t450\t.\t{strand}\t.\tID=wrong;Parent=t\n'
    result, _, _ = utrs.complete(text, text, Genome(), 'S')
    rows = [line.split('\t') for line in result.splitlines() if '\tUTR\t' in line or '_UTR\t' in line]
    assert len(rows) == 2
    high = next(row for row in rows if row[3:5] == ['401', '450'])
    assert high[2] == ('three_prime_UTR' if strand == '+' else 'five_prime_UTR')


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('motif,canonical', [('GTAG', True), ('AAAA', False)])
def test_boundary_in_donor_intron_keeps_projected_utr_and_audits_connection(motif, canonical, strand):
    selected = model(cds=((100, 400),), exons=((100, 400),), strand=strand)
    native = utrs.read_models(model('n', 'ng', cds=((120, 400),), exons=((20, 50), (120, 450)), strand=strand))
    sequence = list('A'*1200)
    genomic_motif = str(Seq(motif).reverse_complement()) if strand == '-' else motif
    sequence[50:52], sequence[98:100] = genomic_motif[:2], genomic_motif[2:]
    _, audit, _ = utrs.complete(selected, selected, Genome(''.join(sequence)), 'S', native)
    chosen = audit[0]['utr5' if strand == '+' else 'utr3']
    assert chosen['length'] == 30 and chosen['new_junction']
    assert chosen['junction_canonical'] == canonical
    assert chosen['junction_motif'] == motif[:2]+'-'+motif[2:]


def test_native_model_overlapping_two_loci_is_not_assigned_by_longest_utr():
    selected = model('a', 'ga', cds=((100, 200),), exons=((100, 200),)) + model('b', 'gb', cds=((300, 400),), exons=((300, 400),))
    native = utrs.read_models(model('n', 'ng', cds=((100, 400),), exons=((10, 490),)))
    _, audit, _ = utrs.complete(selected, selected, Genome(), 'S', native)
    assert all(r['utr5']['length'] == 0 and r['utr3']['length'] == 0 for r in audit)
    assert all(r['native_status'] == 'no_unique_coding_overlap' for r in audit)


@pytest.mark.parametrize('strand', ['+', '-'])
def test_unique_exact_native_cds_keeps_utrs_when_a_neighbor_also_overlaps(strand):
    selected = model('a', 'ga', cds=((100, 250),), exons=((100, 250),), strand=strand)
    selected += model('b', 'gb', cds=((200, 400),), exons=((200, 400),), strand=strand)
    native = utrs.read_models(model('n', 'ng', cds=((100, 250),), exons=((50, 300),), strand=strand))
    _, audit, _ = utrs.complete(selected, selected, Genome(), 'S', native)
    by_id = {row['transcript_id']: row for row in audit}
    assert by_id['a']['utr5']['length'] == by_id['a']['utr3']['length'] == 50
    assert by_id['a']['native_status'] == 'mapped'
    assert by_id['b']['native_status'] == 'no_unique_coding_overlap'


def test_same_length_prefers_existing_and_full_annotation_keeps_other_isoforms():
    selected = model()
    other = model('other', cds=((100, 400),), exons=((50, 450),)).splitlines(keepends=True)[1:]
    full = selected + ''.join(other)
    augmented, audit, _ = utrs.complete(selected, full, Genome(), 'S', utrs.read_models(model('n', 'ng')))
    assert audit[0]['utr5']['source'] == 'existing'
    updates = {'t': dict(model=utrs.read_models(augmented)['t'], exons=[[50, 450]], choices={'5': audit[0]['utr5'], '3': audit[0]['utr3']})}
    rendered = utrs.render(full+'##FASTA\n>chr\nACGT\n', updates)
    assert 'ID=other.e0;Parent=other' in rendered
    assert rendered.endswith('##FASTA\n>chr\nACGT\n')


def test_oriongeno_receipt_and_gff_are_bound_to_target_assembly(tmp_path):
    native = importlib.import_module('annevo_gene_models')
    genome, gff, receipt, table = [tmp_path/p for p in ('genome.fa', 'native.gff3', 'receipt.json', 'inputs.tsv')]
    genome.write_text('>chr\n'+'A'*600+'\n')
    gff.write_text(model())
    value = dict(schema_version=1, kind='genegalleon-oriongeno-native-prediction', native_returncode=0,
                 gpu_probe={'gpu_verified': True}, model_profile='plants', inference_settings={'lineage': 'plants'},
                 runtime_digest='a'*64, cache_key='b'*64, source_commit='c'*40, model_sha256='d'*64,
                 genome_sha256=native.digest(genome), GFF_sha256=native.digest(gff), GFF_bytes=gff.stat().st_size)
    receipt.write_text(json.dumps(value))
    table.write_text('species\tgff\treceipt\nS\tnative.gff3\treceipt.json\n')
    frozen = native.freeze_inputs(table, {'S': {'genome': str(genome)}}, engine='oriongeno')
    assert frozen['predictions']['S']['identity']['model_profile'] == 'plants'
    gff.write_text(gff.read_text()+'# changed\n')
    with pytest.raises(ValueError, match='changed'):
        native.verify_inputs(frozen)


def test_whole_cohort_finalization_preserves_cds_and_binds_utrs(tmp_path):
    from test_gene_model_refinement import tiny_inputs

    refinement = importlib.import_module('gene_model_refinement')
    native_helper = importlib.import_module('annevo_gene_models')
    inputs, edges, rows = tiny_inputs(tmp_path)
    table = tmp_path/'oriongeno.tsv'
    lines = ['species\tgff\treceipt\n']
    originals = {}
    for row in rows:
        name = row['species']
        Path(row['genome']).write_text('>chr\n'+'A'*10+'ATGAAACCCTAA'+'A'*28+'\n')
        source_text = model(cds=((10, 22),), exons=((0, 24),))
        Path(row['gff']).write_text(source_text)
        originals[name] = source_text
        gff, receipt = tmp_path/(name+'.native.gff3'), tmp_path/(name+'.native.json')
        gff.write_text(model('n', 'ng', cds=((10, 22),), exons=((2, 35),)))
        value = dict(schema_version=1, kind='genegalleon-oriongeno-native-prediction', native_returncode=0,
                     gpu_probe={'gpu_verified': True}, model_profile='plants', inference_settings={'lineage': 'plants'},
                     runtime_digest='a'*64, cache_key='b'*64, source_commit='c'*40, model_sha256='d'*64,
                     genome_sha256=native_helper.digest(row['genome']), GFF_sha256=native_helper.digest(gff), GFF_bytes=gff.stat().st_size)
        receipt.write_text(json.dumps(value))
        lines.append(f'{name}\t{gff}\t{receipt}\n')
    table.write_text(''.join(lines))
    before = {row['species']: Path(row['cds']).read_bytes() for row in rows}
    baseline_root = tmp_path/'baseline'
    baseline = refinement.finalize(baseline_root, refinement.plan(baseline_root, inputs=inputs, edges=edges, mode='off'))
    root = tmp_path/'with-utrs'
    result = refinement.finalize(root, refinement.plan(root, inputs=inputs, edges=edges, mode='off', oriongeno_inputs=table))
    assert len(refinement.verify_inputs(result/'inputs.tsv')) == len(rows)
    assert (result/'representative_map.tsv').read_bytes() == (baseline/'representative_map.tsv').read_bytes()
    for row in rows:
        name = row['species']
        assert Path(row['gff']).read_text() == originals[name]
        assert Path(row['cds']).read_bytes() == before[name]
        for role in ('species_cds', 'species_protein', 'analysis_cds'):
            assert (result/role/(name+'.fa')).read_bytes() == (baseline/role/(name+'.fa')).read_bytes()
        audit = json.loads((result/'utr_audit'/(name+'.json')).read_text())[0]
        assert audit['utr5']['source'] == 'existing' and audit['utr5']['length'] == 10
        assert audit['utr3']['source'] == 'oriongeno' and audit['utr3']['length'] == 13
        assert utrs.load_models(result/'species_gff'/(name+'.gff3'))['t']['exons'] == [[0, 35]]
        assert utrs.load_models(result/'full_annotation'/(name+'.gff3'))['t']['exons'] == [[0, 35]]
    assert len((result/'utr_summary.tsv').read_text().splitlines()) == len(rows)+1
    fasta = result/'species_transcript'/(rows[0]['species']+'.fa')
    fasta.write_text(fasta.read_text()+'A')
    with pytest.raises(ValueError, match='corrupt'):
        refinement.verify_inputs(result/'inputs.tsv')


def test_directory_discovery_freezes_only_the_selected_cohort(tmp_path):
    from test_annevo_gene_models import native_fixture

    helper = importlib.import_module('annevo_gene_models')
    _genome, gff, receipt, _table, sources, *_ = native_fixture(tmp_path)
    value = json.loads(receipt.read_text())
    value.update(kind='genegalleon-oriongeno-native-prediction', model_profile='plants', inference_settings={'lineage': 'plants'})
    receipt.write_text(json.dumps(value))
    root = tmp_path/'predictions'
    for species in ('Target_species', 'Other_species'):
        directory = root/species
        directory.mkdir(parents=True)
        (directory/'candidate_inputs.tsv').write_text(f'species\tgff\treceipt\n{species}\t{gff}\t{receipt}\n')
    frozen = helper.freeze_inputs(root, sources, engine='oriongeno')
    assert set(frozen['predictions']) == {'Target_species'}
    helper.verify_inputs(frozen)
    (root/'Target_species/candidate_inputs.tsv').write_text('changed')
    with pytest.raises(ValueError, match='changed'):
        helper.verify_inputs(frozen)


def test_shared_prediction_directory_without_current_cohort_is_a_frozen_empty_subset(tmp_path):
    helper = importlib.import_module('annevo_gene_models')
    folder = tmp_path/'Other_species'
    folder.mkdir()
    (folder/'candidate_inputs.tsv').write_text('species\tgff\treceipt\nOther_species\tmissing\tmissing\n')
    frozen = helper.freeze_inputs(tmp_path, {'Current_species': {}}, engine='oriongeno')
    assert frozen['predictions'] == {} and frozen['tables'] == [] and frozen['files'] == {}
    helper.verify_inputs(frozen)
    current = tmp_path/'Current_species'
    current.mkdir()
    (current/'candidate_inputs.tsv').write_text('new prediction after planning')
    helper.verify_inputs(frozen)
    assert frozen['predictions'] == {}


@pytest.mark.parametrize('setting,present,expected', [('auto', True, 'root'), ('auto', False, ''), ('', True, ''), ('explicit.tsv', True, 'explicit.tsv')])
def test_input_stage_automatic_connection_and_explicit_disable(tmp_path, setting, present, expected):
    core = Path(__file__).resolve().parents[1]/'core/gg_input_generation_core.sh'
    text = core.read_text()
    start = text.index('resolve_oriongeno_candidate_inputs() {')
    function = text[start:text.index('\n}\n', start)+3]
    root = tmp_path/'genome_prediction/oriongeno'
    if present:
        (root/'Species_a').mkdir(parents=True)
        (root/'Species_a/candidate_inputs.tsv').write_text('fixture')
    result = subprocess.run(['bash', '-c', function+'\nresolve_oriongeno_candidate_inputs'],
                            env={**os.environ, 'gene_model_oriongeno_inputs': setting,
                                 'gg_workspace_output_dir': str(tmp_path)}, capture_output=True, text=True, check=True)
    assert result.stdout.strip() == (str(root) if expected == 'root' else expected)

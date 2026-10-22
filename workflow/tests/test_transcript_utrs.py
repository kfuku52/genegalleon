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
def test_shared_terminal_geometry_excludes_introns_and_uses_explicit_cds_endpoints(strand):
    donor = utrs.read_models(model(cds=((100, 180), (220, 300)),
                                  exons=((20, 40), (80, 180), (220, 320), (340, 370)), strand=strand))['t']
    low, high = ('5', '3') if strand == '+' else ('3', '5')
    assert utrs.coding_boundary(donor, low) == 100 and utrs.coding_boundary(donor, high) == 300
    assert utrs.terminal_utr_blocks(donor, low) == [[20, 40], [80, 100]]
    assert utrs.terminal_utr_blocks(donor, high) == [[300, 320], [340, 370]]
    assert utrs.terminal_utr_blocks(donor, low, boundary=30) == [[20, 30]]
    assert utrs.terminal_utr_blocks(donor, high, boundary=380) == []
    assert utrs.exons_within_reference(donor, 370)
    assert not utrs.exons_within_reference(donor, 369)
    assert not utrs.exons_within_reference(donor, None)


@pytest.mark.parametrize('strand,side', [('.', '5'), ('+', 'other')])
def test_shared_geometry_rejects_an_unresolved_side_or_strand(strand, side):
    donor = {'strand': strand, 'cds': [(100, 300, 0)], 'exons': [[50, 350]]}
    with pytest.raises(ValueError, match='side 5 or 3'):
        utrs.terminal_utr_blocks(donor, side)


@pytest.mark.parametrize('changed', ['seqid', 'strand'])
def test_projection_rejects_cross_contig_or_cross_strand_donors(changed):
    selected = utrs.read_models(model())['t']
    donor = {**selected, changed: 'other' if changed == 'seqid' else '-'}
    with pytest.raises(ValueError, match='same contig and strand'):
        utrs.project(donor, selected, '5', Genome())


@pytest.mark.parametrize('strand', ['+', '-'])
def test_shorter_selected_cds_does_not_manufacture_utrs_from_a_cds_only_donor(strand):
    donor = utrs.read_models(model(cds=((100, 400),), exons=((100, 400),), strand=strand))['t']
    selected = utrs.read_models(model(cds=((150, 350),), exons=((150, 350),), strand=strand))['t']
    for side in ('5', '3'):
        result = utrs.project(donor, selected, side, Genome())
        assert result['blocks'] == [] and result['length'] == 0 and result['reason'] == 'no_source_utr'


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
    assert audit[0][low]['source_transcript'] == 't' and audit[0][high]['source_transcript'] == 'native'
    assert audit[0]['annotation_basis'] == 'projected_longest_per_side'
    attributes = [utrs.parse_gff_attributes(r.split('\t')[8]) for r in result.splitlines() if '_prime_UTR\t' in r]
    assert {a['utr_source_transcript'] for a in attributes} == {('t',), ('native',)}
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


def test_fragmented_genome_reference_metadata_is_read_once_per_completion():
    class FragmentedGenome(Genome):
        accesses = 0

        @property
        def references(self):
            self.accesses += 1
            return (*('unused'+str(i) for i in range(16000)), 'chr')

    genome = FragmentedGenome()
    selected = model() + model('t2', 'g2', cds=((600, 900),), exons=((550, 950),))
    native = utrs.read_models(model('n', 'ng'))
    _, audit, sequences = utrs.complete(selected, selected, genome, 'S', native)
    assert genome.accesses == 1
    assert all(row['utr5']['length'] == row['utr3']['length'] == 50 for row in audit)
    assert len(sequences) == 2


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


def coge_segments(strand='+'):
    text = model(cds=((100, 160), (220, 400)), exons=((400, 450), (500, 520)), strand=strand)
    rows = []
    for line in text.splitlines():
        fields = line.replace('\tsource\t', '\tCoGe\t').split('\t')
        if fields[2] in {'gene', 'mRNA'}:
            fields[3] = '101'
        if fields[2] != 'gene':
            fields[8] += ';coge_fid=' + ('9' if fields[2] == 'CDS' else '7')
        rows.append('\t'.join(fields)+'\n')
    return ''.join(rows)


@pytest.mark.parametrize('strand', ['+', '-'])
def test_coge_coding_and_terminal_exon_segments_are_interpreted_automatically(strand):
    text = coge_segments(strand)
    parsed = utrs.read_models(text)['t']
    assert parsed['exons'] == [[100, 160], [220, 450], [500, 520]]
    assert parsed['exon_path_basis'] == 'coge_coding_plus_terminal_noncoding_segments'
    result, audit, seq = utrs.complete(text, text, Genome(), 'S')
    side = 'utr3' if strand == '+' else 'utr5'
    assert audit[0][side]['length'] == 70
    assert audit[0]['source_exon_interpretations'][0]['declared_exons'] == [[400, 450], [500, 520]]
    assert len(seq['t']) == 310
    assert [r for r in result.splitlines() if '\tCDS\t' in r] == [r for r in text.splitlines() if '\tCDS\t' in r]
    assert utrs.read_models(result)['t']['exons'] == parsed['exons']


@pytest.mark.parametrize('change', ['source', 'identity', 'bounds', 'internal_segment', 'shared_parent', 'rna_identity'])
def test_coge_interpretation_requires_source_identity_and_terminal_structure(change):
    text = coge_segments()
    if change == 'source':
        text = text.replace('\tCoGe\t', '\tother\t')
    elif change == 'identity':
        text = text.replace('Parent=t;coge_fid=7', 'Parent=t;coge_fid=8')
    elif change == 'bounds':
        text = text.replace('mRNA\t101\t520', 'mRNA\t91\t520')
    elif change == 'internal_segment':
        text = text.replace('exon\t401\t450', 'exon\t181\t190')
    elif change == 'shared_parent':
        text = text.replace('Parent=t;coge_fid=7', 'Parent=t,other;coge_fid=7')
    else:
        text = text.replace('Parent=g;coge_fid=7', 'Parent=g')
    with pytest.raises(ValueError):
        utrs.read_models(text)
    assert utrs.read_models(text, strict=False)['t'].get('problem')


def test_unique_full_coding_intron_chain_resolves_terminal_neighbour_overlap():
    selected = model('a', 'ga', cds=((100, 200), (300, 400)), exons=((100, 200), (300, 400)))
    selected += model('b', 'gb', cds=((420, 450),), exons=((420, 450),))
    native = utrs.read_models(model('n', 'ng', cds=((90, 200), (300, 470)), exons=((30, 200), (300, 500))))
    _, audit, _ = utrs.complete(selected, selected, Genome(), 'S', native)
    by_id = {r['transcript_id']: r for r in audit}
    assert by_id['a']['utr5']['length'] == 70 and by_id['a']['utr3']['length'] == 100
    match = by_id['a']['native_correspondence']['candidates'][0]
    assert match['basis'] == 'identical_coding_intron_chain' and match['assigned_gene'] == 'ga'
    assert by_id['b']['native_status'] == 'no_unique_coding_overlap'


def test_shared_full_intron_chain_remains_ambiguous():
    selected = model('a', 'ga', cds=((100, 200), (300, 400)), exons=((100, 200), (300, 400)))
    selected += model('b', 'gb', cds=((110, 200), (300, 390)), exons=((110, 200), (300, 390)))
    native = utrs.read_models(model('n', 'ng', cds=((90, 200), (300, 470)), exons=((30, 200), (300, 500))))
    _, audit, _ = utrs.complete(selected, selected, Genome(), 'S', native)
    assert all(r['native_status'] == 'no_unique_coding_overlap' for r in audit)
    assert all(r['native_correspondence']['status'] == 'ambiguous_coding_overlap' for r in audit)


def test_unmapped_native_diagnostics_distinguish_absence_from_opposite_strand():
    selected = utrs.read_models(model('a', 'ga', cds=((100, 200),), exons=((100, 200),))
                                + model('b', 'gb', cds=((600, 700),), exons=((600, 700),)))
    native = utrs.read_models(model('n', 'ng', cds=((100, 200),), exons=((100, 200),), strand='-'))
    _, reasons = utrs.native_correspondence(selected, native)
    assert reasons['a']['status'] == 'opposite_strand_coding_overlap'
    assert reasons['b']['status'] == 'no_coding_overlap'
    _, reasons = utrs.native_correspondence(selected, {})
    assert reasons['a']['status'] == 'no_native_prediction_on_contig'


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


@pytest.mark.parametrize('explicit_native', [False, True])
def test_whole_cohort_finalization_preserves_cds_and_binds_utrs(tmp_path, explicit_native):
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
        if explicit_native:
            import oriongeno_outputs
            staged = tmp_path / (name + '.explicit.gff3')
            oriongeno_outputs.write_explicit_utrs(gff, staged, lengths={'chr': 50})
            gff.write_bytes(staged.read_bytes())
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


@pytest.mark.parametrize('kind', ['tRNA', 'ncRNA', 'transcript'])
def test_standalone_noncoding_features_do_not_block_coding_completion(kind):
    selected = model()
    original = selected + f'chr\tsource\t{kind}\t501\t550\t.\t+\t.\tID=noncoding\n'
    original += 'chr\tsource\texon\t501\t550\t.\t+\t.\tParent=noncoding\n'
    result, audit, sequences = utrs.complete(selected, original, Genome(), 'S')
    assert utrs.read_models(result)['t']['cds'] == [(100, 400, 0)]
    assert audit[0]['utr5']['length'] == audit[0]['utr3']['length'] == 50
    assert sequences['t'] == Genome().sequence[50:450]


@pytest.mark.parametrize('with_exons', [False, True])
def test_completion_without_utrs_keeps_every_gff_byte(with_exons):
    text = model(cds=((100, 400),), exons=((100, 400),))
    if not with_exons:
        text = ''.join(line for line in text.splitlines(keepends=True) if '\texon\t' not in line)
    result, audit, sequences = utrs.complete(text, text, Genome(), 'S')
    assert result == text
    assert sequences['t'] == Genome().sequence[100:400]
    assert audit[0]['utr5']['length'] == audit[0]['utr3']['length'] == 0


@pytest.mark.parametrize('ownership', ['missing', 'ambiguous', 'duplicate'])
def test_coding_transcript_ownership_still_fails_strict_import(ownership):
    text = model()
    if ownership == 'missing':
        text = text.replace('ID=t;Parent=g', 'ID=t')
    elif ownership == 'ambiguous':
        text = text.replace('ID=t;Parent=g', 'ID=t;Parent=g,other')
    else:
        text += 'chr\tsource\tmRNA\t51\t450\t.\t+\t.\tID=t;Parent=g\n'
    with pytest.raises(ValueError, match='unique identity and one gene owner'):
        utrs.read_models(text)
    _, audit, sequences = utrs.complete(text, text, Genome(), 'S')
    assert audit[0]['status'] == 'unsupported_structure'
    assert sequences == {}

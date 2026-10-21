"""UTR-only publications bind consumed inputs and retain all coding rows."""
import hashlib
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
import transcript_utr_reassessment as reassess
from test_transcript_utrs import coge_segments, model


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def completed(tmp_path):
    root = tmp_path / 'source'
    view = root / 'effective'
    view.mkdir(parents=True)
    (root / 'plan.json').write_text(json.dumps(dict(species=['S'], request={})))
    text = coge_segments()
    for role, suffix, value in (
            ('species_gff', '.gff3', text), ('analysis_gff', '.gff3', text.replace('CDS\t101\t160', 'CDS\t104\t160')),
            ('full_annotation', '.gff3', text), ('source_annotation', '.gff3', text),
            ('species_cds', '.fa', '>g\n'+'A'*240+'\n'), ('species_genome', '.fa', '>chr\n'+'A'*1200+'\n'),
            ('utr_audit', '.json', json.dumps([dict(species='S', transcript_id='t', gene_id='g', status='unsupported_structure')]))):
        (view / role).mkdir()
        (view / role / ('S'+suffix)).write_text(value)
    (view / 'representative_map.tsv').write_text('species\tgene_id\tsource_transcript_id\nS\tg\tt\n')
    roles = dict(gff='species_gff/S.gff3', analysis_gff='analysis_gff/S.gff3', cds='species_cds/S.fa', genome='species_genome/S.fa')
    (view / 'inputs.tsv').write_text('\t'.join(['species', *roles, 'genetic_code'])+'\n'
        +'\t'.join(['S', *(str(view/p) for p in roles.values()), '1'])+'\n')
    files = {str(p.relative_to(view)): sha(p) for p in view.rglob('*') if p.is_file()}
    (view / 'receipt.json').write_text(json.dumps(dict(key=dict(plan=sha(root/'plan.json')), files=files)))
    return root, tmp_path / 'reassessment', sha(view/'receipt.json')


def test_reassessment_corrects_coge_without_modifying_source_or_coding_rows(completed):
    source, output, expected = completed
    before = {p: sha(p) for p in source.rglob('*') if p.is_file()}
    result = reassess.reproject(source, output, expected)
    assert all(sha(p) == digest for p, digest in before.items())
    summary = json.loads((result / 'summary.json').read_text())[0]
    assert summary['final_gene_models'] == summary['with_any_utr'] == summary['interpreted_coge_models'] == 1
    assert summary['cds_unchanged']
    for role, old_role in [('species_gff', 'species_gff'), ('full_annotation', 'full_annotation'), ('analysis_gff', 'analysis_gff')]:
        assert reassess._coding_rows((result/role/'S.gff3').read_text()) == reassess._coding_rows((source/'effective'/old_role/'S.gff3').read_text())
    assert (result / 'receipt.json').is_file()
    audit = json.loads((result / 'utr_audit/S.json').read_text())
    assert audit[0]['final_gene_id'] == 'g' and audit[0]['final_exported']


def test_reassessment_reuses_same_generation(completed, monkeypatch):
    source, output, expected = completed
    first = reassess.reproject(source, output, expected)
    monkeypatch.setattr(reassess.utrs, 'complete', lambda *_a, **_kw: pytest.fail('completed publication was recomputed'))
    assert reassess.reproject(source, output, expected) == first


@pytest.mark.parametrize('invalid', ['receipt', 'plan', 'member', 'species', 'inplace'])
def test_reassessment_rejects_changed_or_unbound_inputs(completed, invalid):
    source, output, expected = completed
    species = None
    if invalid == 'receipt':
        expected = '0'*64
    elif invalid == 'plan':
        (source/'plan.json').write_text('{"species":["S"],"request":{"changed":true}}')
    elif invalid == 'member':
        with (source/'effective/species_gff/S.gff3').open('a') as stream:
            stream.write('# changed\n')
    elif invalid == 'species':
        species = ['S', 'S']
    else:
        output = source / 'new'
    with pytest.raises(ValueError):
        reassess.reproject(source, output, expected, species)
    assert not (output/'effective/receipt.json').exists()


def test_mid_command_input_change_cannot_publish(completed, monkeypatch):
    source, output, expected = completed
    complete = reassess.utrs.complete
    def changing(*args, **kwargs):
        result = complete(*args, **kwargs)
        with (source/'effective/species_gff/S.gff3').open('a') as stream:
            stream.write('# changed during UTR computation\n')
        return result
    monkeypatch.setattr(reassess.utrs, 'complete', changing)
    with pytest.raises(OSError, match='changed'):
        reassess.reproject(source, output, expected)
    assert not (output/'effective/receipt.json').exists()


@pytest.mark.parametrize('invalid', [None, 'gff', 'assembly'])
def test_reassessment_binds_native_prediction_and_target_assembly(completed, invalid):
    source, output, expected = completed
    gff, receipt = source/'native.gff3', source/'native.json'
    genome = source/'effective/species_genome/S.fa'
    gff.write_text(model('n', 'ng', cds=((100, 160), (220, 400)), exons=((80, 160), (220, 550))))
    value = dict(schema_version=1, kind='genegalleon-oriongeno-native-prediction', native_returncode=0,
        gpu_probe={'gpu_verified': True}, model_profile='plants', inference_settings={'lineage': 'plants'},
        runtime_digest='a'*64, cache_key='b'*64, source_commit='c'*40, model_sha256='d'*64,
        genome_sha256='0'*64 if invalid == 'assembly' else sha(genome), GFF_sha256=sha(gff), GFF_bytes=gff.stat().st_size)
    receipt.write_text(json.dumps(value))
    identity = {**value, 'receipt_sha256': sha(receipt)}
    plan_path = source/'plan.json'
    plan = json.loads(plan_path.read_text())
    plan['request'] = dict(files={str(gff): sha(gff), str(receipt): sha(receipt)},
        oriongeno=dict(predictions=dict(S=dict(gff=str(gff), receipt=str(receipt), identity=identity))))
    plan_path.write_text(json.dumps(plan))
    effective = source/'effective/receipt.json'
    frozen = json.loads(effective.read_text())
    frozen['key']['plan'] = sha(plan_path)
    effective.write_text(json.dumps(frozen))
    expected = sha(effective)
    if invalid == 'gff':
        gff.write_text(gff.read_text()+'# changed native bytes\n')
    if invalid:
        with pytest.raises(ValueError):
            reassess.reproject(source, output, expected)
        assert not (output/'effective/receipt.json').exists()
    else:
        result = reassess.reproject(source, output, expected)
        audit = json.loads((result/'utr_audit/S.json').read_text())[0]
        assert audit['utr5']['source'] == audit['utr3']['source'] == 'oriongeno'
        assert audit['utr5']['length'] == 20 and audit['utr3']['length'] == 150


def test_standard_cli_exposes_utr_only_reassessment(completed, monkeypatch):
    import gene_model_refinement as refinement
    source, output, expected = completed
    monkeypatch.setattr(sys, 'argv', ['gene_model_refinement.py', 'reproject-utrs',
        '--source-output', str(source), '--source-receipt-sha256', expected, '--output', str(output), '--species', 'S'])
    refinement.main()
    assert (output/'effective/summary.json').is_file()


def test_reassessment_does_not_reuse_previous_projected_utr_as_independent_donor(completed):
    from test_transcript_utrs import Genome
    source, output, expected = completed
    original = (source/'effective/source_annotation/S.gff3').read_text()
    old, _, _ = reassess.utrs.complete(original, original, Genome(), 'S',
        reassess.utrs.read_models(model('n', 'ng', cds=((100, 160), (220, 400)), exons=((80, 160), (220, 700)))))
    view = source/'effective'
    for role in ('species_gff', 'analysis_gff', 'full_annotation'):
        (view/role/'S.gff3').write_text(old)
    frozen = json.loads((view/'receipt.json').read_text())
    frozen['files'] = {name: sha(view/name) for name in frozen['files']}
    (view/'receipt.json').write_text(json.dumps(frozen))
    result = reassess.reproject(source, output, sha(view/'receipt.json'))
    audit = json.loads((result/'utr_audit/S.json').read_text())[0]
    assert audit['utr3']['source'] == 'existing' and audit['utr3']['length'] == 70
    assert audit['utr5']['length'] == 0
    full = reassess.utrs.load_models(result/'full_annotation/S.gff3')['t']
    assert full['exons'] == [[100, 160], [220, 450], [500, 520]]

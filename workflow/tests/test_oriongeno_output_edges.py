"""Boundary regressions for retained native annotations and their verification."""
import importlib
import random
from pathlib import Path

import pytest
from test_oriongeno_outputs import utr_regions
from test_transcript_utrs import Genome, model

outputs = importlib.import_module('oriongeno_outputs')
utrs = importlib.import_module('transcript_utrs')
native = importlib.import_module('annevo_runtime')


def test_escaped_parent_ids_survive_native_export_and_final_completion(tmp_path):
    text = model(tid='t%3b1')
    source, target = tmp_path/'source.gff3', tmp_path/'target.gff3'
    source.write_text(text)
    outputs.write_explicit_utrs(source, target)
    final, _, _ = utrs.complete(target.read_text(), text, Genome(), 'S')
    target.write_text(final)
    assert native.gff_counts(target)['CDS'] == 1


@pytest.mark.parametrize('newline', [b'\n', b'\r\n', b'\r'])
def test_native_export_preserves_original_rows_and_embedded_fasta_bytes(tmp_path, newline):
    annotation = model().encode().replace(b'\n', newline)
    fasta = b'##FASTA'+newline+b'>chr'+newline+b'ACGT'+newline
    source, target = tmp_path/'source.gff3', tmp_path/'target.gff3'
    source.write_bytes(annotation + fasta)
    outputs.write_explicit_utrs(source, target)
    assert source.read_bytes() == annotation + fasta
    assert target.read_bytes().startswith(annotation)
    assert target.read_bytes().endswith(fasta)
    assert b'\tfive_prime_UTR\t' in target.read_bytes().split(b'##FASTA')[0]


def test_generated_utrs_precede_their_gff_group_terminator(tmp_path):
    text = model() + '###\n' + model('t2', 'g2', cds=((600, 900),), exons=((550, 950),)) + '###\n'
    source, target = tmp_path/'source.gff3', tmp_path/'target.gff3'
    source.write_text(text)
    outputs.write_explicit_utrs(source, target)
    groups = target.read_text().split('###\n')
    assert 'Parent=t;' in groups[0] and '\tfive_prime_UTR\t' in groups[0]
    assert 'Parent=t2;' in groups[1] and '\tfive_prime_UTR\t' in groups[1]
    assert not groups[2].strip()
    again = tmp_path/'again.gff3'
    outputs.write_explicit_utrs(target, again)
    assert again.read_bytes() == target.read_bytes()
    final, _, _ = utrs.complete(text, text, Genome(), 'S')
    groups = final.split('###\n')
    assert '\tfive_prime_UTR\t' in groups[0] and '\tfive_prime_UTR\t' in groups[1]
    assert not groups[2].strip()


@pytest.mark.parametrize('rna_row', [False, True])
def test_projected_native_ends_stay_in_the_selected_group(rna_row):
    selected = model(cds=((100, 400),), exons=((100, 400),))
    if not rna_row:
        selected = ''.join(r for r in selected.splitlines(keepends=True) if '\tmRNA\t' not in r and '\texon\t' not in r)
        selected = selected.replace('Parent=t', 'Parent=g')
    selected += '###\n'
    donor = utrs.read_models(model('n', 'ng'))
    final, _, _ = utrs.complete(selected, selected, Genome(), 'S', donor)
    assert final.endswith('###\n') and '\tfive_prime_UTR\t' in final
    assert utrs.complete(final, selected, Genome(), 'S', donor)[0] == final


def test_final_added_features_are_not_swallowed_by_an_unterminated_comment():
    text = model() + '# trailing comment'
    final, _, _ = utrs.complete(text, text, Genome(), 'S')
    assert '# trailing comment\n' in final
    assert sum('_prime_UTR\t' in r for r in final.splitlines() if not r.startswith('#')) == 2


def test_canonical_duplicate_id_cannot_hide_behind_different_percent_case(tmp_path):
    text = model(tid='t%3b1') + 'chr\tsource\texon\t51\t450\t.\t+\t.\tID=t%3B1;Parent=t%3b1\n'
    source = tmp_path/'source.gff3'
    source.write_text(text)
    with pytest.raises(RuntimeError, match='duplicate'):
        native.gff_counts(source)


def test_native_gff_validation_stops_before_embedded_fasta(tmp_path):
    source = tmp_path/'source.gff3'
    source.write_text(model() + '##FASTA\n>chr\nACGT\n')
    assert native.gff_counts(source)['CDS'] == 1


@pytest.mark.parametrize('text', [
    'ID=t;Parent=g;Name=unrelated;Note=ignored',
    ' ID = t%3b1 ; Parent = g%2C1,g2 ;',
    'ID=t;ID=second;Parent=g;Parent=g2',
    'ID="t";Parent=\'g\'',
    'ID=t%252C1;Parent=g%3B1',
    'ID=t;ID=;Parent=,g,,',
    'gene_id "g"; transcript_id "t"; ID "native";',
    'gene_id "a=b"; ID "native";',
    '.',
    'bare_identifier',
])
def test_native_identity_scan_agrees_with_the_full_attribute_reader(text):
    try:
        expected = utrs.parse_gff_attributes(text)
    except ValueError:
        with pytest.raises(ValueError):
            native._gff_identities(text)
        return
    ids, parents = native._gff_identities(text)
    assert tuple(ids) == expected.get('ID', ()) and tuple(parents) == expected.get('Parent', ())


@pytest.mark.parametrize('alias', ['self', 'symlink', 'hardlink'])
def test_streamed_export_cannot_truncate_its_input(tmp_path, alias):
    import os
    source, output = tmp_path/'source.gff3', tmp_path/'out.gff3'
    source.write_text(model())
    before = source.read_bytes()
    if alias == 'self':
        output = source
    elif alias == 'symlink':
        output.symlink_to(source)
    else:
        os.link(source, output)
    with pytest.raises(ValueError, match='overwrite'):
        outputs.write_explicit_utrs(source, output)
    assert source.read_bytes() == before


def test_transcript_cannot_cross_a_closed_gff_group(tmp_path):
    source, output = tmp_path/'source.gff3', tmp_path/'out.gff3'
    source.write_text(model() + '###\nchr\tsource\tUTR\t51\t100\t.\t+\t.\tParent=t\n')
    output.write_bytes(b'previous output')
    with pytest.raises(ValueError, match='closed GFF groups'):
        outputs.write_explicit_utrs(source, output)
    assert output.read_bytes() == b'previous output'


def test_streamed_export_fences_source_changes_between_passes(tmp_path, monkeypatch):
    source, output = tmp_path/'source.gff3', tmp_path/'out.gff3'
    source.write_text(model())
    output.write_bytes(b'previous output')
    read = outputs.read_model_rows
    def changed_after_read(stream):
        result = read(stream)
        source.write_text(source.read_text() + '# changed\n')
        return result
    monkeypatch.setattr(outputs, 'read_model_rows', changed_after_read)
    with pytest.raises(ValueError, match='changed'):
        outputs.write_explicit_utrs(source, output)
    assert output.read_bytes() == b'previous output'


@pytest.mark.parametrize('seed', range(50))
def test_generated_valid_models_keep_consistent_utrs_and_cds_across_both_stages(tmp_path, seed):
    rng = random.Random(seed)
    position, exons = rng.randrange(10), []
    for _ in range(rng.randint(1, 8)):
        length = rng.randint(12, 70)
        exons.append((position, position + length))
        position += length + rng.randint(10, 50)
    start, end = exons[0][0] + 4, exons[-1][1] - 4
    cds = tuple((max(a, start), min(b, end)) for a, b in exons if a < end and b > start)
    strand = '+' if seed % 2 else '-'
    tid = 't%3b' + str(seed)
    text = model(tid=tid, cds=cds, exons=exons, strand=strand)
    if seed % 3 == 0:
        text += '###\n'
    raw = text.replace('\n', '\r\n').encode() if seed % 4 == 0 else text.encode()
    source, explicit, again = [tmp_path / p for p in ['source.gff3', 'explicit.gff3', 'again.gff3']]
    source.write_bytes(raw)
    outputs.write_explicit_utrs(source, explicit, lengths={'chr': 1200})
    outputs.write_explicit_utrs(explicit, again, lengths={'chr': 1200})
    assert source.read_bytes() == raw and explicit.read_bytes() == again.read_bytes()
    expected = {}
    for kind, blocks in [('five_prime_UTR' if strand == '+' else 'three_prime_UTR',
                          [[a, min(b, start)] for a, b in exons if a < start]),
                         ('three_prime_UTR' if strand == '+' else 'five_prime_UTR',
                          [[max(a, end), b] for a, b in exons if b > end])]:
        expected[('t;' + str(seed), kind, strand)] = blocks
    assert utr_regions(explicit.read_text()) == expected
    result, _, _ = utrs.complete(explicit.read_text(), text, Genome(), 'S', utrs.load_models(explicit))
    assert utr_regions(result) == expected
    assert [r for r in result.splitlines() if '\tCDS\t' in r] == [r for r in text.splitlines() if '\tCDS\t' in r]
    again.write_text(result)
    assert native.gff_counts(again)['CDS'] == len(cds)


@pytest.mark.parametrize('changed_after', ['self', 'later_companion'])
def test_companion_hash_verification_fences_concurrent_changes(tmp_path, monkeypatch, changed_after):
    gff = tmp_path/'genes.gff3'
    paths = {k: v for k, v in outputs.companion_paths(gff).items() if k in {'native_gff', 'repeat_gff'}}
    for path in paths.values():
        path.write_bytes(b'sealed companion')
    value = {'inference_settings': {'output_utrs': True, 'output_repeat': True},
             'annotation_outputs': {'schema_version': 1, 'output_utrs': True, 'output_repeat': True,
                'utr_policy': outputs.UTR_POLICY, 'files': {k: outputs.artifact_record(p, p) for k, p in paths.items()}}}
    original = outputs.file_sha
    def changed_after_read(path):
        result = original(path)
        if Path(path) == paths['native_gff'] and changed_after == 'self' or Path(path) == paths['repeat_gff'] and changed_after == 'later_companion':
            paths['native_gff'].write_bytes(b'edited companion')
        return result
    monkeypatch.setattr(outputs, 'file_sha', changed_after_read)
    with pytest.raises((OSError, ValueError), match='companion|changed'):
        outputs.verify_companions(value, gff)

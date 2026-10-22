"""Explicit UTRs preserve coding paths; native repeat companions remain verifiable."""
import gzip
import hashlib
import importlib
import sys
from pathlib import Path

import pytest
from test_transcript_utrs import model

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
outputs = importlib.import_module('oriongeno_outputs')
runtime = importlib.import_module('annevo_runtime')


def exported(tmp_path, text, *, lengths=None):
    source, target = tmp_path / 'source.gff3', tmp_path / 'explicit.gff3'
    source.write_text(text)
    proof = outputs.write_explicit_utrs(source, target, lengths=lengths)
    return source, target, proof


@pytest.mark.parametrize('strand', ['+', '-'])
def test_spliced_utrs_preserve_all_original_rows_and_coding_phase(tmp_path, strand):
    text = model(cds=((100, 180), (220, 300)), exons=((20, 40), (80, 180), (220, 320), (340, 370)), strand=strand)
    source, target, proof = exported(tmp_path, text, lengths={'chr': 400})
    assert source.read_text() == text
    result = target.read_text()
    assert result.startswith(text)
    utrs = [r.split('\t') for r in result.splitlines() if '\tUTR\t' in r or '_prime_UTR\t' in r]
    low, high = ('five_prime_UTR', 'three_prime_UTR') if strand == '+' else ('three_prime_UTR', 'five_prime_UTR')
    assert {(r[2], int(r[3]), int(r[4])) for r in utrs} == {
        (low, 21, 40), (low, 81, 100), (high, 301, 320), (high, 341, 370)}
    assert proof['added_features'] == {'five_prime_UTR': 2, 'three_prime_UTR': 2}
    assert all(r[6] == strand and r[7] == '.' for r in utrs)
    outputs.write_explicit_utrs(target, tmp_path / 'again.gff3')
    assert (tmp_path / 'again.gff3').read_bytes() == target.read_bytes()


def test_cds_only_has_no_invented_utrs_and_embedded_fasta_is_retained(tmp_path):
    text = model(cds=((100, 300),), exons=((100, 300),)) + '##FASTA\n>chr\nACGT\n'
    _, target, proof = exported(tmp_path, text)
    assert target.read_text() == text
    assert sum(proof['added_features'].values()) == 0


def test_added_features_precede_embedded_fasta_and_do_not_lose_comments(tmp_path):
    text = '# a note containing ##FASTA\n' + model() + '##FASTA\n>chr\nACGT\n'
    _, target, _ = exported(tmp_path, text)
    result = target.read_text()
    assert result.startswith('# a note containing ##FASTA\n')
    assert result.endswith('##FASTA\n>chr\nACGT\n')
    assert result.index('\tfive_prime_UTR\t') < result.index('\n##FASTA\n')


def test_partial_existing_utrs_are_completed_without_duplicates_or_id_collisions(tmp_path):
    text = model() + ('chr\tsource\tfive_prime_UTR\t51\t70\t.\t+\t.\tID=t.ggutr.5.1;Parent=t\n')
    _, target, proof = exported(tmp_path, text)
    assert target.read_text().startswith(text)
    assert 'ID=t.ggutr.5.1_;Parent=t' in target.read_text()
    assert '\tfive_prime_UTR\t71\t100\t' in target.read_text()
    assert proof['added_features']['five_prime_UTR'] == 1
    outputs.write_explicit_utrs(target, tmp_path / 'again.gff3')
    assert (tmp_path / 'again.gff3').read_bytes() == target.read_bytes()


def test_parent_identifier_retains_the_original_percent_encoding(tmp_path):
    text = model(tid='t%3b1')
    _, target, _ = exported(tmp_path, text)
    assert ';Parent=t%3b1;annotation_basis=' in target.read_text()
    assert runtime.gff_counts(target)['five_prime_UTR'] == 1


@pytest.mark.parametrize('feature,start,end', [('five_prime_UTR', 101, 120), ('three_prime_UTR', 51, 70), ('UTR', 95, 105)])
def test_conflicting_explicit_utrs_fail_instead_of_changing_cds(tmp_path, feature, start, end):
    text = model() + f'chr\tsource\t{feature}\t{start}\t{end}\t.\t+\t.\tID=bad;Parent=t\n'
    with pytest.raises(ValueError, match='conflicts'):
        exported(tmp_path, text)


@pytest.mark.parametrize('lengths', [{'other': 1000}, {'chr': 420}])
def test_export_rejects_out_of_genome_models(tmp_path, lengths):
    with pytest.raises(ValueError, match='outside'):
        exported(tmp_path, model(), lengths=lengths)


@pytest.mark.parametrize('strand', ['+', '-', '.'])
def test_native_repeat_validation_and_zero_predictions(tmp_path, strand):
    path = tmp_path / 'repeat.gff3'
    path.write_text('##gff-version 3\n')
    assert outputs.validate_repeat_gff(path, {'chr': 30}) == {'repeat_region': 0}
    path.write_text(f'##gff-version 3\nchr\tOrionGeno\trepeat_region\t1\t30\t.\t{strand}\t.\tID=r1\n')
    assert outputs.validate_repeat_gff(path, {'chr': 30}) == {'repeat_region': 1}


@pytest.mark.parametrize('change', ['length', 'seqid', 'kind', 'duplicate', 'header', 'parent'])
def test_invalid_repeat_outputs_fail(tmp_path, change):
    row = 'chr\tOrionGeno\trepeat_region\t1\t30\t.\t+\t.\tID=r1\n'
    text = '##gff-version 3\n' + row
    text = {'length': text.replace('\t30\t', '\t31\t'), 'seqid': text.replace('chr\t', 'absent\t'),
            'kind': text.replace('repeat_region', 'CDS'), 'duplicate': text + row,
            'header': row, 'parent': text.replace('ID=r1', 'ID=r1;Parent=g')}[change]
    path = tmp_path / 'repeat.gff3'
    path.write_text(text)
    with pytest.raises(ValueError):
        outputs.validate_repeat_gff(path, {'chr': 30})


@pytest.mark.parametrize('previous', [False, True])
@pytest.mark.parametrize('failed', ['gene.hmm_boundaries.jsonl.gz', 'gene.utr_repairs.jsonl.gz', 'receipt.json'])
def test_bundle_publication_failure_restores_every_existing_file(tmp_path, monkeypatch, failed, previous):
    destinations = [tmp_path / name for name in ['gene.gff3', 'gene.native.gff3', 'gene.repeat.gff3',
                    'gene.hmm_boundaries.jsonl.gz', 'gene.utr_repairs.jsonl.gz', 'receipt.json']]
    sources = [tmp_path / ('new' + str(i)) for i in range(len(destinations) - 1)]
    if previous:
        for i, path in enumerate(destinations):
            path.write_text('old' + str(i))
    for path in sources:
        path.write_text('new')
    original = runtime.os.replace
    def fail(source, destination):
        if '.staging-' in str(source) and destination.name == failed:
            raise OSError('bundle publication failed')
        return original(source, destination)
    monkeypatch.setattr(runtime.os, 'replace', fail)
    with pytest.raises(OSError, match='bundle publication'):
        runtime.publish_pair(sources[0], destinations[0], destinations[-1], {},
                             companions=list(zip(sources[1:], destinations[1:-1], strict=True)))
    if previous:
        assert [p.read_text() for p in destinations] == ['old' + str(i) for i in range(len(destinations))]
    else:
        assert not any(p.exists() for p in destinations)
    assert not list(tmp_path.glob('*.backup-*')) and not list(tmp_path.glob('*.staging-*'))


@pytest.mark.parametrize('damage', ['missing', 'changed', 'symlink', 'contract'])
def test_companions_are_required_for_new_receipts_but_not_invented_for_old_ones(tmp_path, damage):
    gff = tmp_path / 'gene.gff3'
    paths = outputs.companion_paths(gff)
    paths = {key: paths[key] for key in ('native_gff', 'repeat_gff')}
    for path in paths.values():
        path.write_text('sealed companion')
    value = {'inference_settings': {'output_repeat': True, 'output_utrs': True},
             'annotation_outputs': {'schema_version': 1, 'output_repeat': True, 'output_utrs': True,
                                   'utr_policy': outputs.UTR_POLICY,
                                   'files': {k: outputs.artifact_record(p, p) for k, p in paths.items()}}}
    outputs.verify_companions(value, gff)
    if damage == 'missing':
        paths['repeat_gff'].unlink()
    elif damage == 'changed':
        paths['repeat_gff'].write_text('changed')
    elif damage == 'symlink':
        paths['repeat_gff'].unlink()
        paths['repeat_gff'].symlink_to(paths['native_gff'])
    else:
        value['annotation_outputs']['files'].pop('repeat_gff')
    with pytest.raises(ValueError, match='companion'):
        outputs.verify_companions(value, gff)
    outputs.verify_companions({}, gff)


def test_new_output_flags_cannot_be_passed_off_as_a_historical_receipt(tmp_path):
    with pytest.raises(ValueError, match='missing'):
        outputs.verify_companions({'inference_settings': {'output_repeat': True, 'output_utrs': True}}, tmp_path/'x.gff3')


@pytest.mark.parametrize('text', [b'', b'{"status":"joined"}\n{"status":"repaired","note":"\xe6\xa4\x8d\xe7\x89\xa9"}'])
def test_audit_gzip_retains_exact_bytes_and_is_independent_of_temporary_names(tmp_path, text):
    source, compressed, again = [tmp_path / name for name in ('audit.jsonl', 'audit.jsonl.gz', 'other.jsonl.gz')]
    source.write_bytes(text)
    metadata = outputs.compress_audit_jsonl(source, compressed)
    source.rename(tmp_path / 'renamed.jsonl')
    assert outputs.compress_audit_jsonl(tmp_path / 'renamed.jsonl', again) == metadata
    assert compressed.read_bytes() == again.read_bytes()
    assert gzip.decompress(compressed.read_bytes()) == text
    assert metadata == dict(compression='gzip', uncompressed_sha256=hashlib.sha256(text).hexdigest(),
                            uncompressed_bytes=len(text), jsonl_records=len(text.splitlines()))


@pytest.mark.parametrize('damage', ['missing', 'symlink', 'invalid-json', 'nonobject', 'blank'])
def test_audit_compression_rejects_missing_and_invalid_sources(tmp_path, damage):
    source = tmp_path / 'audit.jsonl'
    if damage == 'symlink':
        target = tmp_path / 'target.jsonl'
        target.write_text('{}\n')
        source.symlink_to(target)
    elif damage != 'missing':
        source.write_bytes({'invalid-json': b'{broken}\n', 'nonobject': b'[]\n', 'blank': b'\n'}[damage])
    with pytest.raises(ValueError, match='audit'):
        outputs.compress_audit_jsonl(source, tmp_path / 'audit.jsonl.gz')


@pytest.mark.parametrize('name', ['hmm_boundaries', 'utr_repairs'])
@pytest.mark.parametrize('damage', ['missing', 'changed', 'contract', 'downgrade', 'metadata'])
def test_audit_companions_are_bound_to_versioned_receipts(tmp_path, name, damage):
    gff = tmp_path / 'genes.gff3'
    paths = outputs.companion_paths(gff)
    records = {}
    for key, path in paths.items():
        if key in outputs.AUDIT_SUFFIXES:
            source = tmp_path / (key + '.jsonl')
            source.write_text('{}\n')
            metadata = outputs.compress_audit_jsonl(source, path)
        else:
            path.write_text('sealed companion')
            metadata = {}
        records[key] = outputs.artifact_record(path, path) | metadata
    settings = {'output_repeat': True, 'output_utrs': True, 'audit_log_compression': 'gzip'}
    value = {'inference_settings': settings.copy(), 'annotation_outputs': {
        'schema_version': 2, **settings, 'utr_policy': outputs.UTR_POLICY, 'files': records}}
    outputs.verify_companions(value, gff)
    if damage == 'missing':
        paths[name].unlink()
    elif damage == 'changed':
        paths[name].write_bytes(b'changed gzip')
    elif damage == 'contract':
        records.pop(name)
    elif damage == 'metadata':
        records[name]['jsonl_records'] = -1
    else:
        value['annotation_outputs']['schema_version'] = 1
        for key in outputs.AUDIT_SUFFIXES:
            records.pop(key)
    with pytest.raises(ValueError, match='companion|audit'):
        outputs.verify_companions(value, gff)


def test_cli_defaults_and_stage_flags_are_explicit(monkeypatch, capsys):
    import oriongeno_runtime as native
    captured = []
    monkeypatch.setattr(native, 'run', lambda args: captured.append(args) or {})
    monkeypatch.setattr(sys, 'argv', ['oriongeno_runtime.py', 'run', '--receipt', 'r', '--genome', 'g',
                                   '--output', 'o', '--output-receipt', 'p'])
    native.main()
    assert captured[0].output_repeat == captured[0].output_utrs == 1
    capsys.readouterr()
    from test_annevo_genome_annotation import source_function
    result = source_function('gg_validate_genome_annotation_backend; printf "%s %s" "$oriongeno_output_utrs" "$oriongeno_output_repeat"', env={'genome_annotation_backend': 'oriongeno'})
    assert result.returncode == 0 and result.stdout == '1 1'
    for flag in ['oriongeno_output_utrs', 'oriongeno_output_repeat']:
        result = source_function('gg_validate_genome_annotation_backend', env={'genome_annotation_backend': 'oriongeno', flag: '2'})
        assert result.returncode != 0

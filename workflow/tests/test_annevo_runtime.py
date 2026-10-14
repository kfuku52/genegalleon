"""Optional runtime integrity and failure-safe publication, with small offline assets."""
import argparse
import hashlib
import importlib
import io
import json
import os
import platform
import subprocess
import sys
import tarfile
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / 'support'
sys.path.insert(0, str(SUPPORT))
cache = importlib.import_module('annevo_runtime_cache')
download = importlib.import_module('annevo_runtime_download')
runtime = importlib.import_module('annevo_runtime')


def test_actual_cli_cold_offline_prepare_reports_original_error_without_network(tmp_path, monkeypatch, capsys):
    command = [str(SUPPORT / 'annevo_runtime.py'), 'prepare', '--cache-root', str(tmp_path / 'cli-cache'), '--no-download']
    result = subprocess.run([sys.executable, '-E', '-s', '-B', *command], capture_output=True, text=True, timeout=15)
    expected = 'ANNEVO: ANNEVO runtime manifest is not cached and downloads are disabled\n'
    assert result.returncode == 1
    assert result.stdout == '' and result.stderr == expected
    assert not list((tmp_path / 'cli-cache').rglob('READY.json'))

    # The same unmodified CLI path must reject before its resolver can access
    # the registry, even though the default locator names a public OCI artifact.
    monkeypatch.setattr(cache, 'resolve_manifest', lambda *a: pytest.fail('offline CLI touched the network resolver'))
    monkeypatch.setattr(sys, 'argv', [*command[:3], str(tmp_path / 'in-process-cache'), '--no-download'])
    assert runtime.main() == 1
    assert capsys.readouterr().err == expected


@pytest.mark.parametrize('failure', [OSError('native download failed'), tarfile.ReadError('invalid runtime archive')])
def test_cli_reports_prepare_io_and_tar_failures_without_replacing_them(tmp_path, monkeypatch, capsys, failure):
    def fail(*a, **kwargs):
        raise failure
    monkeypatch.setattr(runtime, 'prepare', fail)
    monkeypatch.setattr(sys, 'argv', [str(SUPPORT / 'annevo_runtime.py'), 'prepare', '--cache-root', str(tmp_path / 'cache')])
    assert runtime.main() == 1
    captured = capsys.readouterr()
    assert captured.out == ''
    assert captured.err == 'ANNEVO: ' + str(failure) + '\n'


def bundle(tmp_path, monkeypatch):
    builder_spec = importlib.util.spec_from_file_location('bundle_builder_fixture', SUPPORT.parents[1] / 'container/scripts/build_annevo_runtime_bundle.py')
    builder = importlib.util.module_from_spec(builder_spec)
    builder_spec.loader.exec_module(builder)
    fixture_spec = importlib.util.spec_from_file_location('bundle_fixture', Path(__file__).with_name('test_annevo_runtime_bundle.py'))
    fixtures = importlib.util.module_from_spec(fixture_spec)
    fixture_spec.loader.exec_module(fixtures)
    source = fixtures.official_fixture(tmp_path)
    monkeypatch.setattr(builder, 'inspect_environment', lambda *a: {'machine': platform.machine(), 'torch_version': 'fixture', 'compiled_cuda': 'fixture'})
    def pack(prefix, output, executable):
        fixtures.environment_archive(output)
        return {'strict_missing_files': True}
    monkeypatch.setattr(builder, 'pack_environment', pack)
    abi = {'system': 'Linux', 'machine': platform.machine(), 'libc': dict(zip(('name', 'version'), platform.libc_ver(), strict=True)), 'image_digest': 'sha256:'+'a'*64}
    path, manifest = builder.build_bundle(tmp_path/'env', source, tmp_path/'distribution', abi, 'fixture')
    monkeypatch.setattr(download, 'platform_name', lambda: manifest['platform'])
    original_fetch = cache.fetch_verified
    def fetch(url, target, sha, size, **kwargs):
        if 'raw.githubusercontent.com' in url:
            url = (source/'saved_model'/Path(url).name).as_uri()
        return original_fetch(url, target, sha, size, **kwargs)
    monkeypatch.setattr(cache, 'fetch_verified', fetch)
    monkeypatch.setattr(cache.subprocess, 'run', lambda *a, **k: subprocess.CompletedProcess(a[0], 0, '', ''))
    monkeypatch.setattr(cache, '_cpu_check', lambda *a: manifest['cuda_profile'])
    return path, manifest


def test_cold_install_warm_reuse_and_selected_receipt_stability(tmp_path, monkeypatch):
    path, manifest = bundle(tmp_path, monkeypatch)
    root = tmp_path/'cache'
    first = cache.prepare(root, 'Magnoliopsida', str(path))
    selected = Path(first['receipt_path'])
    before = selected.read_bytes(), selected.stat().st_mtime_ns
    monkeypatch.setattr(cache, 'resolve_manifest', lambda *a: pytest.fail('warm access touched network'))
    monkeypatch.setattr(cache.subprocess, 'run', lambda *a, **k: pytest.fail('warm access reinstalled environment'))
    second = cache.prepare(root, 'Magnoliopsida', str(path), no_download=True)
    assert first['reused'] is False and second['reused'] is True
    assert first['generation'] == second['generation']
    assert (selected.read_bytes(), selected.stat().st_mtime_ns) == before
    assert len(list((root/'runtimes').rglob('g-*'))) == 2  # generation + inventory
    Path(first['model_path']).write_bytes(b'corrupt')
    with pytest.raises(RuntimeError, match='checksum'):
        cache.prepare(root, 'Magnoliopsida', str(path), no_download=True)


def test_failed_install_never_publishes_ready_and_can_retry(tmp_path, monkeypatch):
    path, _ = bundle(tmp_path, monkeypatch)
    def fail(*a):
        raise RuntimeError('dependency import failed')
    monkeypatch.setattr(cache, '_cpu_check', fail)
    root = tmp_path/'cache'
    with pytest.raises(RuntimeError, match='dependency'):
        cache.prepare(root, 'Magnoliopsida', str(path))
    assert not list(root.rglob('READY.json'))
    with pytest.raises(RuntimeError, match='not installed'):
        cache.prepare(root, 'Magnoliopsida', str(path), no_download=True)
    monkeypatch.setattr(cache, '_cpu_check', lambda *a: {})
    assert cache.prepare(root, 'Magnoliopsida', str(path))['reused'] is False


@pytest.mark.parametrize('failed_entrypoint', [None, 'prediction.py', 'decoding.py'])
def test_relocated_cpu_readiness_checks_every_native_cli_and_inference_dependencies(tmp_path, monkeypatch, failed_entrypoint):
    source = tmp_path / 'source'
    source.mkdir()
    receipt = {'python': str(tmp_path / 'environment/bin/python'), 'native_entrypoint': str(source / 'annotation.py')}
    profile = {'torch_version': 'OFFLINE', 'compiled_cuda': 'OFFLINE'}
    calls = []

    def run(command, **kwargs):
        calls.append((command, kwargs))
        if '-c' in command:
            code = command[-1]
            assert 'torchmetrics' not in code and 'Bio,tqdm' in code
            return subprocess.CompletedProcess(command, 0, json.dumps(profile), '')
        assert kwargs['cwd'] == source.parent
        assert command[1:4] == ['-E', '-s', '-B']
        if Path(command[-2]).name == failed_entrypoint:
            return subprocess.CompletedProcess(command, 1, '', 'ModuleNotFoundError: missing native model')
        return subprocess.CompletedProcess(command, 0, 'OFFLINE native help', '')

    monkeypatch.setattr(cache.subprocess, 'run', run)
    if failed_entrypoint:
        with pytest.raises(RuntimeError, match=failed_entrypoint):
            cache._cpu_check(receipt, {'cuda_profile': profile})
    else:
        assert cache._cpu_check(receipt, {'cuda_profile': profile}) == profile
        assert [Path(command[-2]).name for command, _ in calls[1:]] == ['annotation.py', 'prediction.py', 'decoding.py']


def test_runtime_modification_and_redirected_manifest_are_rejected(tmp_path, monkeypatch):
    path, _ = bundle(tmp_path, monkeypatch)
    root = tmp_path/'cache'
    receipt = cache.prepare(root, 'Magnoliopsida', str(path))
    changed = Path(receipt['native_entrypoint'])
    changed.chmod(0o600)
    changed.write_bytes(b'changed')
    with pytest.raises(RuntimeError, match='member changed'):
        runtime.read_selected(receipt['receipt_path'])
    alias = next((root/'locators').glob('*.json'))
    raw = alias.read_bytes()
    alias.unlink()
    external = tmp_path/'foreign.json'
    external.write_bytes(raw)
    alias.symlink_to(external)
    with pytest.raises(RuntimeError, match='redirected'):
        cache.prepare(root, 'Magnoliopsida', str(path))


@pytest.mark.parametrize('failure', ['stale_pax_path', 'normalized_duplicate'])
def test_archive_rejection_identifies_actual_member_without_accepting_malformed_payload(tmp_path, failure):
    archive = tmp_path / 'malformed.tar.gz'
    with tarfile.open(archive, 'w:gz', format=tarfile.PAX_FORMAT) as stream:
        if failure == 'stale_pax_path':
            member = tarfile.TarInfo('environment/lib/original.so')
            member.pax_headers = {'path': 'lib/original.so'}
            stream.addfile(member)
        else:
            stream.addfile(tarfile.TarInfo('source/a'))
            stream.addfile(tarfile.TarInfo('source/./a'))
    output = tmp_path / 'unpacked'
    output.mkdir()
    with pytest.raises(RuntimeError) as error:
        cache._safe_extract(archive, output)
    if failure == 'stale_pax_path':
        assert str(error.value) == "Unsafe ANNEVO archive member: name='lib/original.so', normalized='lib/original.so'"
        assert not list(output.iterdir())
    else:
        assert str(error.value) == "Duplicate ANNEVO archive member: name='source/./a', normalized='source/a'"


@pytest.mark.parametrize('name,link', [('../escape', None), ('source/escape', '../../outside'), ('source/a/../b', None)])
def test_archive_traversal_and_external_links_are_rejected(tmp_path, name, link):
    archive = tmp_path/'bad.tar.gz'
    with tarfile.open(archive, 'w:gz') as stream:
        member = tarfile.TarInfo(name)
        if link:
            member.type, member.linkname = tarfile.SYMTYPE, link
            stream.addfile(member)
        else:
            member.size = 1
            stream.addfile(member, io.BytesIO(b'x'))
    with pytest.raises((RuntimeError, tarfile.TarError)):
        cache._safe_extract(archive, tmp_path/'out')
    assert not (tmp_path/'escape').exists()


def test_pair_publication_failure_restores_both_existing_outputs(tmp_path, monkeypatch):
    prediction, output, receipt = [tmp_path/n for n in ('new', 'output.gff3', 'receipt.json')]
    prediction.write_bytes(b'new prediction')
    output.write_bytes(b'old prediction')
    receipt.write_bytes(b'old receipt')
    original = os.replace
    def fail(source, destination):
        if '.staging-' in str(source) and destination == receipt:
            raise OSError('publication failure')
        return original(source, destination)
    monkeypatch.setattr(runtime.os, 'replace', fail)
    with pytest.raises(OSError, match='publication failure'):
        runtime.publish_pair(prediction, output, receipt, {'fresh': True})
    assert output.read_bytes() == b'old prediction'
    assert receipt.read_bytes() == b'old receipt'
    assert sorted(p.name for p in tmp_path.iterdir()) == ['new', 'output.gff3', 'receipt.json']


def test_gpu_unavailable_does_not_invoke_native_predictor(monkeypatch):
    monkeypatch.setattr(runtime.subprocess, 'run', lambda *a, **k: subprocess.CompletedProcess(a[0], 1, '', 'CUDA is unavailable'))
    with pytest.raises(RuntimeError, match='allocated GPU probe failed'):
        runtime.probe_gpu({'python': sys.executable}, {})


def test_oci_blob_identity_and_part_mapping(monkeypatch):
    part = {'name': 'part.000', 'sha256': 'b'*64, 'size': 100}
    raw = json.dumps({'parts': [part]}).encode()
    sha = hashlib.sha256(raw).hexdigest()
    index = {'layers': [{'mediaType': download.MANIFEST_MEDIA, 'digest': 'sha256:'+sha, 'size': len(raw)}, {'digest': 'sha256:'+'b'*64, 'size': 100}]}
    monkeypatch.setattr(download, 'registry_headers', lambda *a: ('https://ghcr.io/v2/fixture/runtime', 'latest', {'Authorization': 'ephemeral'}))
    monkeypatch.setattr(download, '_read', lambda url, *a: json.dumps(index).encode() if '/manifests/' in url else raw)
    document, _ = download.resolve_manifest('oci://ghcr.io/fixture/runtime:latest')
    assert document['parts'][0]['url'] == 'https://ghcr.io/v2/fixture/runtime/blobs/sha256:'+'b'*64
    index['layers'][1]['size'] = 101
    with pytest.raises(RuntimeError, match='missing'):
        download.resolve_manifest('oci://ghcr.io/fixture/runtime:latest')


def test_native_failure_retains_old_pair_and_logs(tmp_path, monkeypatch):
    genome, selected, output, receipt = [tmp_path/n for n in ('genome.fa', 'selected.json', 'result.gff3', 'result.json')]
    for path in (genome, selected, output, receipt):
        path.write_bytes(b'original')
    genome.write_bytes(b'>seq\nATGAAATAA\n')
    value = {'model_profile': 'Magnoliopsida', 'python': sys.executable, 'model_path': str(tmp_path/'model'), 'generation': str(tmp_path/'installed'), 'native_entrypoint': str(tmp_path/'annotation.py')}
    monkeypatch.setattr(runtime, 'read_selected', lambda *a: value)
    monkeypatch.setattr(runtime, 'probe_gpu', lambda *a: {'gpu_verified': True})
    monkeypatch.setattr(runtime, 'inference_settings', lambda *a: {'overlap_pred': False})
    monkeypatch.setattr(runtime.subprocess, 'run', lambda *a, **k: subprocess.CompletedProcess(a[0], 1))
    args = argparse.Namespace(receipt=str(selected), genome=str(genome), output=str(output), output_receipt=str(receipt), model_profile='Magnoliopsida', threads=4, batch_size=32, sequence_threshold=25)
    with pytest.raises(RuntimeError, match='native prediction failed'):
        runtime.run(args)
    assert output.read_bytes() == receipt.read_bytes() == b'original'
    assert output.with_name(output.name+'.logs').is_dir()


def test_two_parallel_callers_install_one_generation(tmp_path, monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    path, manifest = bundle(tmp_path, monkeypatch)
    checks = []
    monkeypatch.setattr(cache, '_cpu_check', lambda *a: checks.append(True) or manifest['cuda_profile'])
    def prepare():
        return cache.prepare(tmp_path/'shared-cache', 'Magnoliopsida', str(path))
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(lambda _: prepare(), range(2)))
    assert checks == [True]
    assert len({r['generation'] for r in results}) == 1
    assert sorted(r['reused'] for r in results) == [False, True]


@pytest.mark.parametrize("profile", ["Mammalia", "Aves", "Actinopteri", "Insecta", "Fungi", "Magnoliopsida"])
def test_gzip_genome_is_unpacked_in_short_native_work_and_original_is_bound(tmp_path, monkeypatch, profile):
    import gzip
    genome, selected, output, receipt = [tmp_path/n for n in ('genome.fa.gz', 'selected.json', 'result.gff3', 'result.json')]
    fasta = b'>seq\nATGAAATAA\n'
    genome.write_bytes(gzip.compress(fasta))
    selected.write_bytes(b'original selected')
    value = {'model_profile': profile, 'python': sys.executable, 'model_path': str(tmp_path/'model'), 'generation': str(tmp_path/'installed'), 'native_entrypoint': str(tmp_path/'annotation.py'), 'runtime_digest': 'a'*64, 'cache_key': 'b'*64, 'source_commit': 'c'*40, 'model_sha256': 'd'*64}
    monkeypatch.setattr(runtime, 'read_selected', lambda *a: value)
    monkeypatch.setattr(runtime, 'probe_gpu', lambda *a: {'gpu_verified': True})
    (tmp_path / 'prediction.py').write_text('def main():\n length_config=' + repr({profile: [30720, 5120, 32]}) + '\n')
    def predict(argv, **kwargs):
        assert argv[argv.index("-l") + 1] == profile
        assert ("--overlap_pred" in argv) == (profile in {"Mammalia", "Actinopteri"})
        assert argv[argv.index("-m") + 1] == value["model_path"]
        native_genome = Path(argv[argv.index('-g')+1])
        assert native_genome.read_bytes() == fasta
        assert str(native_genome).startswith('/tmp/gg-annevo-')
        assert len(str(native_genome)) < 100
        Path(argv[argv.index('-o')+1]).write_text('seq\tANNEVO\tgene\t1\t9\t.\t+\t.\tID=g1\nseq\tANNEVO\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\nseq\tANNEVO\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=t1\n')
        return subprocess.CompletedProcess(argv, 0)
    monkeypatch.setattr(runtime.subprocess, 'run', predict)
    args = argparse.Namespace(receipt=str(selected), genome=str(genome), output=str(output), output_receipt=str(receipt), model_profile=profile, threads=4, batch_size=32, sequence_threshold=25)
    result = runtime.run(args)
    assert result['genome_sha256'] == download.file_sha(genome)
    assert result['native_genome_sha256'] == hashlib.sha256(fasta).hexdigest()
    assert result['genome_gzip'] is True
    assert result['GFF_sha256'] == download.file_sha(output)


def test_mixed_profiles_share_environment_and_fetch_only_selected_weights(tmp_path, monkeypatch):
    path, manifest = bundle(tmp_path, monkeypatch)
    # The tiny offline bundle remains unchanged; model objects are separate.
    source = tmp_path / 'source/saved_model'
    for profile in ('Mammalia', 'Fungi'):
        weight = source / ('ANNEVO_' + profile + '.pt')
        weight.write_bytes(('fixture ' + profile).encode())
        manifest['models'][profile] = {
            'path': 'source/saved_model/' + weight.name,
            'official_url': 'https://raw.githubusercontent.com/xjtu-omics/ANNEVO/' + manifest['source']['commit'] + '/saved_model/' + weight.name,
            'sha256': download.file_sha(weight), 'size': weight.stat().st_size,
        }
    path.write_text(json.dumps(manifest))
    from concurrent.futures import ThreadPoolExecutor
    checks = []
    monkeypatch.setattr(cache, '_cpu_check', lambda *a: checks.append(True) or manifest['cuda_profile'])
    def selected(profile):
        return cache.prepare(tmp_path / 'shared', profile, str(path), receipt_path=tmp_path / profile / 'runtime.json')
    with ThreadPoolExecutor(2) as pool:
        results = list(pool.map(selected, ('Mammalia', 'Fungi')))
    assert len(checks) == 1
    assert len({r['generation'] for r in results}) == 1
    assert {p.name for p in (tmp_path / 'shared/models').glob('*/*.pt')} == {'ANNEVO_Mammalia.pt', 'ANNEVO_Fungi.pt'}
    for result in results:
        profile = result['model_profile']
        assert json.loads(Path(result['receipt_path']).read_text())['model_profile'] == profile
        warm = cache.prepare(tmp_path / 'shared', profile, str(path), receipt_path=result['receipt_path'], no_download=True)
        assert warm['generation'] == result['generation']

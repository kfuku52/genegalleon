"""Bounded provider-boundary tests; these do not qualify native CUDA inference."""
import argparse
import gzip
import importlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
native = importlib.import_module('oriongeno_runtime')


@pytest.fixture
def local_runtime(tmp_path, monkeypatch):
    package, model = tmp_path/'package', tmp_path/'checkpoint'
    package.mkdir()
    (package/'__init__.py').write_text('# native package fixture\n')
    model.mkdir()
    (model/'config.json').write_text('{}')
    (model/'model.safetensors').write_bytes(b'fixture weights')
    def execute(argv, _env, **_kwargs):
        if '-c' in argv:
            return subprocess.CompletedProcess(argv, 0, json.dumps({'package': str(package)}))
        return subprocess.CompletedProcess(argv, 0, 'oriongeno @ fixture\ntorch==fixture\n')
    monkeypatch.setattr(native, 'execute', execute)
    args = argparse.Namespace(cache_root=str(tmp_path/'cache'), receipt=str(tmp_path/'selected.json'),
                              python=sys.executable, checkpoint=str(model), source_commit='c'*40,
                              model_profile='plants', no_download=True)
    return args, package, model


def test_external_runtime_prepare_and_warm_offline_verification(local_runtime):
    args, _package, _model = local_runtime
    first = native.prepare(args)
    before = Path(args.receipt).read_bytes(), Path(args.receipt).stat().st_mtime_ns
    assert native.prepare(args) == first
    assert (Path(args.receipt).read_bytes(), Path(args.receipt).stat().st_mtime_ns) == before


def test_upstream_pytorch_bin_checkpoint_is_supported(local_runtime):
    args, _package, model = local_runtime
    (model/'model.safetensors').rename(model/'pytorch_model.bin')
    assert native.prepare(args)['model_path'] == str(model)


@pytest.mark.parametrize('option,value', [('source_commit', 'd'*40), ('source_ref', 'another-branch'),
                                         ('cache_root', 'another-cache')])
def test_warm_runtime_rejects_changed_provenance_and_cache(local_runtime, tmp_path, option, value):
    args, _package, _model = local_runtime
    native.prepare(args)
    original = Path(args.receipt).read_bytes()
    setattr(args, option, str(tmp_path/value) if option == 'cache_root' else value)
    with pytest.raises(ValueError, match='differs'):
        native.prepare(args)
    assert Path(args.receipt).read_bytes() == original


def test_runtime_cache_key_cannot_diverge_from_recorded_digest(local_runtime):
    args, _package, _model = local_runtime
    native.prepare(args)
    value = json.loads(Path(args.receipt).read_text())
    value['cache_key'] = 'e'*64
    Path(args.receipt).write_text(json.dumps(value))
    with pytest.raises(ValueError, match='changed'):
        native.read_selected(args.receipt)


def test_runtime_environment_keeps_native_caches_local_and_drops_python_injection(tmp_path, monkeypatch):
    for key in ('PYTHONHOME', 'PYTHONPATH', 'PYTHONPYCACHEPREFIX', 'HF_HUB_CACHE', 'TRITON_CACHE_DIR'):
        monkeypatch.setenv(key, '/unrelated')
    env = native.environment(tmp_path, 3)
    assert all(key not in env for key in ('PYTHONHOME', 'PYTHONPATH', 'PYTHONPYCACHEPREFIX'))
    for key in ('HF_HUB_CACHE', 'TRITON_CACHE_DIR', 'PIP_CACHE_DIR', 'CONDA_PKGS_DIRS', 'MAMBA_ROOT_PREFIX'):
        assert Path(env[key]).is_relative_to(tmp_path)
    assert env['NUMBA_NUM_THREADS'] == '3'


@pytest.mark.parametrize('changed', ['package', 'model', 'receipt'])
def test_runtime_snapshot_rejects_changes(local_runtime, changed):
    args, package, model = local_runtime
    native.prepare(args)
    if changed == 'package':
        (package/'__init__.py').write_text('# changed\n')
    elif changed == 'model':
        (model/'model.safetensors').write_bytes(b'changed weights')
    else:
        value = json.loads(Path(args.receipt).read_text())
        value['model_profile'] = 'mammals'
        Path(args.receipt).write_text(json.dumps(value))
    with pytest.raises((ValueError, RuntimeError), match='changed'):
        native.read_selected(args.receipt)


def test_cold_offline_runtime_does_not_install_or_download(tmp_path):
    args = argparse.Namespace(cache_root=str(tmp_path/'cache'), receipt=str(tmp_path/'selected.json'),
                              python='', checkpoint='', source_commit='', model_profile='plants', no_download=True)
    with pytest.raises(ValueError, match='cold'):
        native.prepare(args)
    assert not Path(args.receipt).exists()


def test_collect_tables_resolves_relative_paths_and_reuses_identical_output(tmp_path):
    root = tmp_path/'predictions'
    for name in ('Species_a', 'Species_b'):
        folder = root/name
        folder.mkdir(parents=True)
        (folder/'candidate_inputs.tsv').write_text('species\tgff\treceipt\n'+name+'\tnative.gff3\tnative.json\n')
    output = tmp_path/'combined.tsv'
    assert native.collect(root, output)['species'] == 2
    assert str(root/'Species_a/native.gff3') in output.read_text()
    before = output.stat().st_mtime_ns
    native.collect(root, output)
    assert output.stat().st_mtime_ns == before


@pytest.mark.parametrize('target', ['candidate_inputs.tsv', 'native.gff3', 'native.json', 'gff-alias'])
def test_collect_does_not_overwrite_native_inputs(tmp_path, target):
    root = tmp_path/'predictions'
    folder = root/'Species_a'
    folder.mkdir(parents=True)
    table = folder/'candidate_inputs.tsv'
    table.write_text('species\tgff\treceipt\nSpecies_a\tnative.gff3\tnative.json\n')
    for name in ('native.gff3', 'native.json'):
        (folder/name).write_text('preserved input')
    if target == 'gff-alias':
        os.link(folder/'native.gff3', folder/target)
    before = {p: p.read_bytes() for p in folder.iterdir()}
    with pytest.raises(ValueError, match='overwrite'):
        native.collect(root, folder/target)
    assert {p: p.read_bytes() for p in folder.iterdir()} == before


def test_native_log_redirect_cannot_modify_runtime(local_runtime, tmp_path):
    setup, package, _model = local_runtime
    native.prepare(setup)
    genome = tmp_path/'genome.fa'
    genome.write_text('>chr\nATGAAATAA\n')
    gff = tmp_path/'native.gff3'
    gff.with_name(gff.name+'.logs').symlink_to(package, target_is_directory=True)
    args = argparse.Namespace(receipt=setup.receipt, genome=str(genome), output=str(gff),
                              output_receipt=str(tmp_path/'native.json'), threads=1, batch_size=1)
    with pytest.raises(ValueError, match='runtime'):
        native.run(args)
    assert sorted(p.name for p in package.iterdir()) == ['__init__.py']


def test_native_failure_preserves_previous_complete_pair(local_runtime, tmp_path, monkeypatch):
    setup, _package, _model = local_runtime
    runtime = native.prepare(setup)
    genome, gff, receipt = [tmp_path/name for name in ('genome.fa', 'prediction.gff3', 'native.json')]
    genome.write_text('>chr\nATGAAATAA\n')
    gff.write_text('previous GFF')
    receipt.write_text('previous receipt')
    monkeypatch.setattr(native, 'read_selected', lambda _path: runtime)
    def execute(argv, _env, **_kwargs):
        if '-c' in argv:
            return subprocess.CompletedProcess(argv, 0, json.dumps({'gpu_verified': True}))
        raise subprocess.CalledProcessError(1, argv)
    monkeypatch.setattr(native, 'execute', execute)
    args = argparse.Namespace(receipt=setup.receipt, genome=str(genome), output=str(gff),
                              output_receipt=str(receipt), threads=2, batch_size=1)
    with pytest.raises(subprocess.CalledProcessError):
        native.run(args)
    assert gff.read_text() == 'previous GFF' and receipt.read_text() == 'previous receipt'


@pytest.mark.parametrize('compressed', [False, True])
def test_successful_native_publication_binds_actual_gff_and_genome(local_runtime, tmp_path, monkeypatch, compressed):
    from test_transcript_utrs import model

    setup, _package, _model = local_runtime
    runtime = native.prepare(setup)
    original_execute = native.execute
    genome = tmp_path/('genome.fa.gz' if compressed else 'genome.fa')
    sequence = b'>chr\nCCCATGAAATAACCC\n'
    genome.write_bytes(gzip.compress(sequence) if compressed else sequence)
    gff, receipt = tmp_path/'native.gff3', tmp_path/'separate-receipts/native.json'
    def execute(argv, env, **kwargs):
        if '--output' in argv:
            assert Path(argv[argv.index('--genome')+1]).read_bytes() == sequence
            Path(argv[argv.index('--output')+1]).write_text(model(cds=((3, 12),), exons=((0, 15),)))
            return subprocess.CompletedProcess(argv, 0, '')
        if '-c' in argv and 'gpu_verified' in argv[argv.index('-c')+1]:
            return subprocess.CompletedProcess(argv, 0, json.dumps({'gpu_verified': True}))
        return original_execute(argv, env, **kwargs)
    monkeypatch.setattr(native, 'execute', execute)
    args = argparse.Namespace(receipt=setup.receipt, genome=str(genome), output=str(gff), output_receipt=str(receipt), threads=2, batch_size=1)
    value = native.run(args)
    assert json.loads(receipt.read_text()) == value
    assert value['kind'] == 'genegalleon-oriongeno-native-prediction'
    assert value['genome_sha256'] == native.file_sha(genome)
    assert value['GFF_sha256'] == native.file_sha(gff)
    assert value['runtime_digest'] == runtime['runtime_digest']
    assert value['feature_counts']['CDS'] == 1

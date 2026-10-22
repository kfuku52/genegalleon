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


def write_prediction(argv, text):
    path = Path(argv[argv.index('--output')+1])
    path.write_text(text)
    if argv[argv.index('--output-repeat')+1] == 'True':
        path.with_name(path.stem+'.repeat.gff3').write_text('##gff-version 3\nchr\tOrionGeno\trepeat_region\t1\t3\t.\t+\t.\tID=r1\n')


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


def test_production_conda_c_compiler_is_exposed_to_triton(tmp_path, monkeypatch):
    monkeypatch.delenv('CC', raising=False)
    monkeypatch.setattr(native.shutil, 'which', lambda name: '/opt/conda/bin/'+name
                        if name == 'x86_64-conda-linux-gnu-cc' else None)
    assert native.environment(tmp_path)['CC'] == '/opt/conda/bin/x86_64-conda-linux-gnu-cc'
    monkeypatch.setenv('CC', '/explicit/compiler')
    assert native.environment(tmp_path)['CC'] == '/explicit/compiler'


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


@pytest.mark.parametrize('shared', [False, True])
@pytest.mark.parametrize('changed', ['package', 'model', 'receipt'])
def test_read_selected_fences_changes_during_package_inventory(local_runtime, monkeypatch, shared, changed):
    args, package, model = local_runtime
    native.prepare(args)
    targets = {'package': package / '__init__.py', 'model': model / 'model.safetensors',
               'receipt': Path(args.receipt)}
    original = native.execute
    def changed_during_freeze(argv, env, **kwargs):
        result = original(argv, env, **kwargs)
        target = targets[changed]
        target.write_bytes(target.read_bytes() + b'changed during pip freeze\n')
        return result
    monkeypatch.setattr(native, 'execute', changed_during_freeze)
    verification = native.RuntimeVerification() if shared else None
    with pytest.raises((OSError, RuntimeError), match='changed'):
        native.read_selected(args.receipt, verification=verification)


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
    monkeypatch.setattr(native, 'read_selected', lambda _path, **kw: runtime)
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


@pytest.mark.parametrize('compression', ['plain', 'gzip-suffix', 'gzip-magic', 'gzip-upper'])
@pytest.mark.parametrize('auto_selected', [False, True])
def test_successful_native_publication_binds_actual_gff_and_genome(local_runtime, tmp_path, monkeypatch, compression, auto_selected):
    from test_transcript_utrs import model

    setup, _package, _model = local_runtime
    runtime = native.prepare(setup)
    original_execute = native.execute
    compressed = compression != 'plain'
    filename = {'gzip-suffix': 'genome.fa.gz', 'gzip-upper': 'genome.fa.GZ'}.get(compression, 'genome.fa')
    genome = tmp_path/filename
    sequence = b'>chr\nCCCATGAAATAACCC\n'
    genome.write_bytes(gzip.compress(sequence) if compressed else sequence)
    gff, receipt = tmp_path/'native.gff3', tmp_path/'separate-receipts/native.json'
    def execute(argv, env, **kwargs):
        if '--output' in argv:
            assert Path(argv[argv.index('--genome')+1]).read_bytes() == sequence
            write_prediction(argv, model(cds=((3, 12),), exons=((0, 15),)))
            return subprocess.CompletedProcess(argv, 0, '')
        if '-c' in argv and 'gpu_verified' in argv[argv.index('-c')+1]:
            return subprocess.CompletedProcess(argv, 0, json.dumps({'gpu_verified': True}))
        return original_execute(argv, env, **kwargs)
    monkeypatch.setattr(native, 'execute', execute)
    monkeypatch.setenv('ORIONGENO_SPECIES_NAME', 'Unrelated_species')
    args = argparse.Namespace(receipt=setup.receipt, genome=str(genome), output=str(gff), output_receipt=str(receipt), threads=2, batch_size=1)
    if auto_selected:
        from test_oriongeno_model_selection import fixture, selection

        selected = fixture(tmp_path, 3702)
        # The fixture creates a genome too; retain the sequence under test.
        genome.write_bytes(gzip.compress(sequence) if compressed else sequence)
        selected.genome = str(genome)
        record = selection.select(selected)
        args.selection_receipt, args.species = selected.output, selected.species
    value = native.run(args)
    assert json.loads(receipt.read_text()) == value
    assert value['kind'] == 'genegalleon-oriongeno-native-prediction'
    assert value['genome_sha256'] == native.file_sha(genome)
    assert value['GFF_sha256'] == native.file_sha(gff)
    assert value['runtime_digest'] == runtime['runtime_digest']
    assert value['feature_counts']['CDS'] == 1
    assert value['feature_counts']['five_prime_UTR'] == value['feature_counts']['three_prime_UTR'] == 1
    assert value['annotation_outputs']['output_repeat'] is True
    assert value['annotation_outputs']['output_utrs'] is True
    for name, artifact in value['annotation_outputs']['files'].items():
        path = native.companion_paths(gff)[name]
        assert artifact['sha256'] == native.file_sha(path)
        assert artifact['bytes'] == path.stat().st_size
    assert '\tUTR\t' not in native.companion_paths(gff)['native_gff'].read_text()
    from oriongeno_outputs import verify_companions
    verify_companions(value, gff)
    from annevo_gene_models import prediction_identity
    assert prediction_identity(gff, receipt, value['genome_sha256'], engine='oriongeno')['GFF_sha256'] == value['GFF_sha256']
    assert value['argv'][value['argv'].index('--species-name')+1] == ''
    if auto_selected:
        assert value['model_selection'] == record
        assert value['model_selection_sha256'] == native.file_sha(args.selection_receipt)
        assert value['species'] == 'Test_species'


@pytest.mark.parametrize('changed', ['mapping', 'selection', 'taxonomy'])
def test_auto_selection_change_during_native_prediction_preserves_previous_pair(local_runtime, tmp_path, monkeypatch, changed):
    import sqlite3

    from test_oriongeno_model_selection import fixture, selection
    from test_transcript_utrs import model

    setup, _package, _model = local_runtime
    native.prepare(setup)
    selected = fixture(tmp_path, 3702)
    selection.select(selected)
    gff, receipt = tmp_path/'native.gff3', tmp_path/'native.json'
    gff.write_text('previous GFF')
    receipt.write_text('previous receipt')
    original_execute = native.execute
    def execute(argv, env, **kwargs):
        if '--output' in argv:
            write_prediction(argv, model(cds=((0, 9),), exons=((0, 9),)))
            if changed == 'taxonomy':
                with sqlite3.connect(selected.taxonomy_db) as conn:
                    conn.execute("INSERT INTO synonym VALUES(3702,'late synonym')")
            else:
                target = Path(selected.taxid_map if changed == 'mapping' else selected.output)
                target.write_text(target.read_text()+'\n')
            return subprocess.CompletedProcess(argv, 0, '')
        if '-c' in argv and 'gpu_verified' in argv[argv.index('-c')+1]:
            return subprocess.CompletedProcess(argv, 0, json.dumps({'gpu_verified': True}))
        return original_execute(argv, env, **kwargs)
    monkeypatch.setattr(native, 'execute', execute)
    args = argparse.Namespace(receipt=setup.receipt, genome=selected.genome, output=str(gff), output_receipt=str(receipt),
                              threads=2, batch_size=1, selection_receipt=selected.output, species=selected.species)
    with pytest.raises((OSError, ValueError, RuntimeError), match='changed'):
        native.run(args)
    assert gff.read_text() == 'previous GFF' and receipt.read_text() == 'previous receipt'


@pytest.mark.parametrize('utrs,repeats', [(0, 0), (0, 1), (1, 0)])
def test_native_output_controls_are_recorded_and_preserve_original_models(local_runtime, tmp_path, monkeypatch, utrs, repeats):
    from test_transcript_utrs import model
    setup, _, _ = local_runtime
    native.prepare(setup)
    original = native.execute
    text = model(cds=((3, 12),), exons=((0, 15),))
    def execute(argv, env, **kwargs):
        if '--output' in argv:
            write_prediction(argv, text)
            return subprocess.CompletedProcess(argv, 0, '')
        if '-c' in argv and 'gpu_verified' in argv[argv.index('-c')+1]:
            return subprocess.CompletedProcess(argv, 0, json.dumps({'gpu_verified': True}))
        return original(argv, env, **kwargs)
    monkeypatch.setattr(native, 'execute', execute)
    genome = tmp_path/'genome.fa'
    genome.write_text('>chr\nCCCATGAAATAACCC\n')
    output = tmp_path/'genes.gff3'
    args = argparse.Namespace(receipt=setup.receipt, genome=str(genome), output=str(output),
                              output_receipt=str(tmp_path/'prediction.json'), threads=1, batch_size=1,
                              output_utrs=utrs, output_repeat=repeats)
    value = native.run(args)
    paths = native.companion_paths(output)
    assert paths['native_gff'].read_text() == text
    assert paths['repeat_gff'].is_file() == bool(repeats)
    assert ('five_prime_UTR' in output.read_text()) == bool(utrs)
    assert [r for r in output.read_text().splitlines() if '\tCDS\t' in r] == [r for r in text.splitlines() if '\tCDS\t' in r]
    from oriongeno_outputs import verify_companions
    verify_companions(value, output)


@pytest.mark.parametrize('changed', ['mapping', 'selection', 'taxonomy', 'taxonomy-wal', 'genome'])
def test_final_runtime_inventory_cannot_hide_changed_selection_inputs(local_runtime, tmp_path, monkeypatch, changed):
    from test_oriongeno_model_selection import fixture, selection
    from test_transcript_utrs import model

    setup, _package, _model = local_runtime
    native.prepare(setup)
    selected = fixture(tmp_path, 3702)
    selection.select(selected)
    gff, receipt = tmp_path/'native.gff3', tmp_path/'native.json'
    gff.write_text('previous GFF')
    receipt.write_text('previous receipt')
    original_execute = native.execute
    inventories = 0
    def execute(argv, env, **kwargs):
        nonlocal inventories
        if '--output' in argv:
            write_prediction(argv, model(cds=((0, 9),), exons=((0, 9),)))
            return subprocess.CompletedProcess(argv, 0, '')
        if '-c' in argv and 'gpu_verified' in argv[argv.index('-c')+1]:
            return subprocess.CompletedProcess(argv, 0, json.dumps({'gpu_verified': True}))
        result = original_execute(argv, env, **kwargs)
        if 'freeze' in argv:
            inventories += 1
            if inventories == 2:
                if changed == 'taxonomy-wal':
                    Path(selected.taxonomy_db+'-wal').write_bytes(b'active writer')
                else:
                    target = Path({'mapping': selected.taxid_map, 'selection': selected.output,
                                   'taxonomy': selected.taxonomy_db, 'genome': selected.genome}[changed])
                    target.write_bytes(target.read_bytes()+b'\n')
        return result
    monkeypatch.setattr(native, 'execute', execute)
    args = argparse.Namespace(receipt=setup.receipt, genome=selected.genome, output=str(gff), output_receipt=str(receipt),
                              threads=2, batch_size=1, selection_receipt=selected.output, species=selected.species)
    with pytest.raises((OSError, ValueError, RuntimeError), match='changed|active transaction'):
        native.run(args)
    assert gff.read_text() == 'previous GFF' and receipt.read_text() == 'previous receipt'


def test_external_virtualenv_keeps_its_own_package_environment(tmp_path):
    import venv

    environment = tmp_path/'venv'
    venv.EnvBuilder(system_site_packages=True, with_pip=False, symlinks=True).create(environment)
    python = environment/'bin/python'
    code = "import sysconfig; print(sysconfig.get_path('purelib'))"
    site = Path(subprocess.check_output([python, '-c', code], text=True).strip())
    package = site/'oriongeno'
    package.mkdir()
    (package/'__init__.py').write_text('# isolated native package fixture\n')
    model = tmp_path/'checkpoint'
    model.mkdir()
    (model/'config.json').write_text('{}')
    (model/'model.safetensors').write_bytes(b'fixture weights')
    args = argparse.Namespace(cache_root=str(tmp_path/'cache'), receipt=str(tmp_path/'selected.json'),
                              python=str(python), checkpoint=str(model), source_commit='c'*40,
                              model_profile='plants', no_download=True)
    value = native.prepare(args)
    assert Path(value['python']) == python
    assert native.read_selected(args.receipt) == value


def test_native_log_hardlink_cannot_truncate_checkpoint(local_runtime, tmp_path, monkeypatch):
    setup, _package, model = local_runtime
    runtime = native.prepare(setup)
    monkeypatch.setattr(native, 'read_selected', lambda _path, **kw: runtime)
    original_execute = native.execute
    def execute(argv, env, **kwargs):
        if '-c' in argv and 'gpu_verified' in argv[argv.index('-c')+1]:
            return subprocess.CompletedProcess(argv, 0, json.dumps({'gpu_verified': True}))
        if '--output' in argv:
            from test_transcript_utrs import model
            Path(argv[argv.index('--output')+1]).write_text(model(cds=((0, 9),), exons=((0, 9),)))
            return subprocess.CompletedProcess(argv, 0, '')
        return original_execute(argv, env, **kwargs)
    monkeypatch.setattr(native, 'execute', execute)
    genome, output = tmp_path/'genome.fa', tmp_path/'prediction.gff3'
    genome.write_text('>chr\nATGAAATAA\n')
    logs = output.with_name(output.name+'.logs')
    logs.mkdir()
    weights = model/'model.safetensors'
    os.link(weights, logs/'native.stdout.txt')
    before = weights.read_bytes()
    args = argparse.Namespace(receipt=setup.receipt, genome=str(genome), output=str(output),
                              output_receipt=str(tmp_path/'native.json'), threads=1, batch_size=1)
    with pytest.raises(ValueError, match='overwrite|runtime'):
        native.run(args)
    assert weights.read_bytes() == before


@pytest.mark.parametrize('link_kind', ['symlink', 'hardlink'])
@pytest.mark.parametrize('target_role', ['package', 'weights', 'previous_gff', 'previous_receipt'])
def test_log_aliases_preserve_runtime_and_previous_publication(local_runtime, tmp_path, link_kind, target_role):
    setup, package, checkpoint = local_runtime
    genome, output, receipt = [tmp_path/name for name in ('genome.fa', 'prediction.gff3', 'native.json')]
    genome.write_text('>chr\nATGAAATAA\n')
    output.write_text('previous GFF')
    receipt.write_text('previous receipt')
    targets = {'package': package/'__init__.py', 'weights': checkpoint/'model.safetensors',
               'previous_gff': output, 'previous_receipt': receipt}
    target = targets[target_role]
    logs = output.with_name(output.name+'.logs')
    logs.mkdir()
    log = logs/'native.stderr.txt'
    if link_kind == 'symlink':
        log.symlink_to(target)
    else:
        os.link(target, log)
    native.prepare(setup)
    before = {p: p.read_bytes() for p in [*targets.values(), genome, Path(setup.receipt)]}
    args = argparse.Namespace(receipt=setup.receipt, genome=str(genome), output=str(output),
                              output_receipt=str(receipt), threads=1, batch_size=1)
    with pytest.raises(ValueError, match='overwrite'):
        native.run(args)
    assert {p: p.read_bytes() for p in before} == before


def test_warm_override_distinguishes_venvs_with_the_same_base_interpreter(local_runtime, tmp_path):
    setup, _package, _model = local_runtime
    pythons = []
    for name in ('venv-a', 'venv-b'):
        folder = tmp_path/name/'bin'
        folder.mkdir(parents=True)
        python = folder/'python'
        python.symlink_to(sys.executable)
        pythons.append(python)
    setup.python = str(pythons[0])
    native.prepare(setup)
    before = Path(setup.receipt).read_bytes()
    setup.python = str(pythons[1])
    with pytest.raises(ValueError, match='differs'):
        native.prepare(setup)
    assert Path(setup.receipt).read_bytes() == before


def test_legacy_python_receipt_accepts_an_alias_in_the_same_environment(local_runtime, tmp_path):
    import shutil

    setup, _package, _model = local_runtime
    folder = tmp_path/'external/bin'
    folder.mkdir(parents=True)
    binary, alias = folder/'python3', folder/'python'
    shutil.copy2(sys.executable, binary)
    alias.symlink_to(binary.name)
    # Older receipts recorded the dereferenced executable name.
    setup.python = str(binary)
    expected = native.prepare(setup)
    before = Path(setup.receipt).read_bytes()
    setup.python = str(alias)
    assert native.prepare(setup) == expected
    assert Path(setup.receipt).read_bytes() == before


def test_runtime_preparation_respects_shared_filesystem_publication_lock(local_runtime):
    from concurrent.futures import ThreadPoolExecutor, TimeoutError
    from threading import Event

    from shared_namespace_lock import namespace_lock

    setup, _package, _model = local_runtime
    expected = native.prepare(setup)
    started = Event()
    def prepare():
        started.set()
        return native.prepare(setup)
    with ThreadPoolExecutor(max_workers=1) as executor:
        with namespace_lock(Path(setup.cache_root)/'.prepare.lock', exclusive=True):
            pending = executor.submit(prepare)
            assert started.wait(timeout=5)
            with pytest.raises(TimeoutError):
                pending.result(timeout=.2)
        assert pending.result(timeout=5) == expected

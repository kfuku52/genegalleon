"""The read-only preparation path cannot start assembly or claim workflow completion."""
import importlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / 'support'
sys.path.insert(0, str(SUPPORT))
prep = importlib.import_module('transcriptome_read_preparation')


def file_snapshot(path):
    info = path.stat()
    # Reading may advance atime; retain file identity, write metadata and bytes.
    identity = (info.st_dev, info.st_ino, info.st_mode, info.st_nlink, info.st_uid,
                info.st_gid, info.st_size, info.st_mtime_ns, info.st_ctime_ns)
    return identity, path.read_bytes()


def fixture(tmp_path):
    root = tmp_path.resolve()
    for directory in ('workspace/input/amalgkit_metadata', 'downloads', 'runtime', 'reads', 'receipts', 'task-locks'):
        (root / directory).mkdir(parents=True, exist_ok=True)
    metadata = root / 'workspace/input/amalgkit_metadata/Plant_one_metadata.tsv'
    metadata.write_text('run\tscientific_name\tsample_group\nSRR100\tPlant one\tflower\nSRR101\tPlant one\tleaf\nSRR102\tPlant one\troot\n')
    image = root / 'image.sif'
    image.write_bytes(b'qualified fixture image')
    core = 'runtime/workflow/core/gg_transcriptome_generation_core.sh'
    (root / core).parent.mkdir(parents=True)
    (root / core).write_text('native core fixture')
    row = dict(task_id=1, species='Plant_one', runs=['SRR100', 'SRR101', 'SRR102'],
               files=[dict(path=str(metadata), basename=metadata.name, sha256=prep.file_sha(metadata), size=metadata.stat().st_size)])
    manifest = dict(root=str(root), project_root=str(root), workspace_root=str(root / 'workspace'),
                    downloads_root=str(root / 'downloads'), fastq_source_root=str(root / 'reads'),
                    busco_lineage='embryophyta_odb12', image_path=str(image), image_sha256=prep.file_sha(image),
                    container_binary='/fixture/apptainer', input_mode=prep.INPUT_MODE, species=[row],
                    managed_files={core: dict(sha256=prep.file_sha(root / core))})
    from test_genome_annotation_cohort import owning_release
    release = owning_release(manifest, 'gg_transcriptome_generation')
    manifest['runtime_release_root'] = str(release)
    manifest['runtime_release_manifest_sha256'] = prep.file_sha(release / '.kfauto-release.json')
    return manifest, row


def test_inspection_is_exact_and_does_not_write_metadata(tmp_path):
    manifest, row = fixture(tmp_path)
    path = Path(row['files'][0]['path'])
    before = file_snapshot(path)
    assert prep.inspect_metadata(path, row['files'][0]['sha256'], row['runs'])['run_ids'] == row['runs']
    assert before == file_snapshot(path)
    with pytest.raises(ValueError, match='membership'):
        prep.inspect_metadata(path, row['files'][0]['sha256'], row['runs'][:2])
    with pytest.raises(ValueError, match='changed'):
        prep.inspect_metadata(path, '0' * 64, row['runs'])


def test_actual_inspection_cli_reads_the_selected_metadata(tmp_path):
    _, row = fixture(tmp_path)
    record = row['files'][0]
    result = subprocess.run([sys.executable, '-S', '-B', str(SUPPORT / 'transcriptome_read_preparation.py'), '--inspect-json',
                             json.dumps(dict(path=record['path'], sha256=record['sha256'], runs=row['runs']))],
                            capture_output=True, text=True, check=True)
    assert json.loads(result.stdout) == dict(metadata_sha256=record['sha256'], run_ids=row['runs'], files=[record['basename']])


def test_container_preflight_binds_actual_metadata_without_requiring_future_fastqs(tmp_path, monkeypatch):
    manifest, row = fixture(tmp_path)
    commands = []
    expected = prep.inspect_metadata(row['files'][0]['path'], row['files'][0]['sha256'], row['runs'])
    def run(command, **kwargs):
        commands.append(command)
        return subprocess.CompletedProcess(command, 0, json.dumps(expected), '')
    monkeypatch.setattr(prep.subprocess, 'run', run)
    result = prep.container_preflight(manifest, row, 'f' * 64)
    assert result['schema_version'] == 2 and result['input_mode'] == 'metadata-read-preparation'
    assert result['metadata_sha256'] == row['files'][0]['sha256']
    assert commands[0][0] == manifest['container_binary']
    assert '--inspect-json' in commands[0]
    assert not list((tmp_path / 'reads').iterdir())


def test_only_getfastq_can_run_and_offline_second_pass_is_required(tmp_path, monkeypatch):
    manifest, row = fixture(tmp_path)
    calls = []
    def entrypoint(native, name, env):
        calls.append(env)
        assert name == 'gg_transcriptome_generation'
        assert env['GG_TRANSCRIPTOME_MODE_TRANSCRIPTOME_ASSEMBLY'] == 'metadata'
        assert env['GG_COMMON_BUSCO_LINEAGE'] == 'embryophyta_odb12'
        assert env['GG_TRANSCRIPTOME_REMOVE_AMALGKIT_FASTQ_AFTER_COMPLETION'] == '0'
        assert env['GG_TRANSCRIPTOME_TRANSCRIPTOME_GETFASTQ_CACHE_DIR'] == manifest['fastq_source_root']
        assert all(env['GG_TRANSCRIPTOME_RUN_' + flag] == '0' for flag in prep.DISABLED)
        assert env['GG_ARRAY_TASK_ID'] == '1'
    monkeypatch.setattr(prep, 'entrypoint', entrypoint)
    completion = Path(manifest['fastq_source_root']) / row['species'] / 'getfastq_completion.json'
    completion.parent.mkdir()
    completion.write_text('native completion fixture')
    provenance = Path(manifest['workspace_root']) / 'output/transcriptome_assembly/artifact_provenance/Plant_one.getfastq.json'
    provenance.parent.mkdir(parents=True)
    provenance.write_text('native provenance fixture')
    monkeypatch.setattr(prep, 'verify_native', lambda *args: dict(native_provenance_sha256=prep.file_sha(provenance),
                                                              completion_manifest_sha256=prep.file_sha(completion)))
    result = prep.prepare(manifest, row)
    assert [env['GG_TRANSCRIPTOME_RUN_AMALGKIT_GETFASTQ'] for env in calls] == ['1', '0']
    assert result['status'] == 'read-preparation-completed'
    assert result['assembly_complete'] is False
    assert result['completion_manifest_sha256'] == prep.file_sha(completion)


def test_zero_exit_without_native_completion_does_not_count_as_prepared(tmp_path, monkeypatch):
    manifest, row = fixture(tmp_path)
    monkeypatch.setattr(prep, 'entrypoint', lambda *args: None)
    with pytest.raises((FileNotFoundError, ValueError, prep.ProvenanceError)):
        prep.prepare(manifest, row)


def recovered_fixture(tmp_path, monkeypatch):
    manifest, row = fixture(tmp_path)
    cache = Path(manifest['fastq_source_root']) / row['species']
    cache.mkdir(parents=True)
    (cache / 'getfastq_completion.json').write_text('completed native manifest')
    parameters = dict(rrna_filter='no', rrna_filter_jobs='1', rrna_filter_chunk_spots='5000000',
        rrna_filter_memory_limit='32G', contam_filter='no', contam_filter_rank='domain',
        read_name='trinity', remove_sra='yes', remove_tmp='yes', dump_print='yes',
        aws='yes', ncbi='yes', ena='yes', redo='no')
    provenance = Path(manifest['workspace_root']) / 'output/transcriptome_assembly/artifact_provenance/Plant_one.getfastq.json'
    provenance.parent.mkdir(parents=True)
    provenance.write_text(json.dumps(dict(parameters=parameters)))
    contract = cache / '.getfastq_cache_contract.json'
    contract.write_text(json.dumps(dict(schema_version=1, species=row['species'],
        metadata_sha256=row['files'][0]['sha256'], metadata_size_bytes=row['files'][0]['size'],
        completion_manifest='getfastq_completion.json', content_stat_fingerprint='[]', parameters=parameters)))
    def verify(*args):
        return dict(native_provenance_sha256=prep.file_sha(provenance),
                    getfastq_provenance_sha256=prep.file_sha(provenance),
                    completion_manifest_sha256=prep.file_sha(cache / 'getfastq_completion.json'))
    monkeypatch.setattr(prep, 'verify_native', verify)
    monkeypatch.setattr(prep, 'entrypoint', lambda *args: pytest.fail('recovery must not run the workflow'))
    def run(command, **kwargs):
        assert command[-2] == '--inspect-cache-json'
        value = json.loads(command[-1])
        return subprocess.CompletedProcess(command, 0, json.dumps({
            **prep.inspect_metadata(value['path'], value['sha256'], value['runs']),
            'prepared_cache': value['expected']}), '')
    monkeypatch.setattr(prep.subprocess, 'run', run)
    return manifest, row, contract


def test_recovery_only_verifies_native_cache_and_preserves_every_existing_byte(tmp_path, monkeypatch):
    manifest, row, _ = recovered_fixture(tmp_path, monkeypatch)
    before = {str(path): file_snapshot(path) for path in tmp_path.rglob('*') if path.is_file()}
    result = prep.recover_prepared_cache(manifest, row)
    assert result['status'] == 'read-preparation-completed' and result['assembly_complete'] is False
    assert result['runs'] == row['runs']
    assert before == {str(path): file_snapshot(path) for path in tmp_path.rglob('*') if path.is_file()}


@pytest.mark.parametrize('key,value', [('metadata_sha256', '0' * 64), ('species', 'Other_plant'),
    ('metadata_size_bytes', 0), ('parameters', {}), ('completion_manifest', 'other.json')])
def test_recovery_rejects_cache_parameter_or_metadata_drift(tmp_path, monkeypatch, key, value):
    manifest, row, contract = recovered_fixture(tmp_path, monkeypatch)
    original = json.loads(contract.read_text())
    contract.write_text(json.dumps({**original, key: value}))
    with pytest.raises(ValueError, match='contract'):
        prep.recover_prepared_cache(manifest, row)


def test_recovery_rejects_failed_native_audit_and_changed_cache_view(tmp_path, monkeypatch):
    manifest, row, _ = recovered_fixture(tmp_path, monkeypatch)
    monkeypatch.setattr(prep.subprocess, 'run', lambda command, **kwargs: subprocess.CompletedProcess(command, 0, '{}', ''))
    with pytest.raises(ValueError, match='changed'):
        prep.recover_prepared_cache(manifest, row)
    monkeypatch.setattr(prep, 'verify_native', lambda *args: (_ for _ in ()).throw(ValueError('native audit failed')))
    with pytest.raises(ValueError, match='native audit'):
        prep.recover_prepared_cache(manifest, row)


@pytest.mark.parametrize('stale_atime', [False, True])
def test_recovery_receipt_retains_source_failure_identity_and_separate_verification_job(tmp_path, monkeypatch, stale_atime):
    manifest, row, _ = recovered_fixture(tmp_path, monkeypatch)
    for name, record in manifest['managed_files'].items():
        record['sha256'] = prep.file_sha(tmp_path / name)
    path = tmp_path / 'manifest.json'
    path.write_text(json.dumps(manifest))
    monkeypatch.setenv('SLURM_JOB_ID', '778')
    monkeypatch.setenv('SLURM_ARRAY_JOB_ID', '777')
    monkeypatch.setenv('SLURM_ARRAY_TASK_ID', '1')
    monkeypatch.setattr(prep, 'container_preflight', lambda *args: dict(native_input_verified=True))
    monkeypatch.setattr(sys, 'argv', ['native-prep', '--package-root', str(tmp_path), '--manifest-sha256',
        prep.file_sha(path), '--task-id', '1', '--recover-source-job-id', '123'])
    prep.main()
    receipt = tmp_path / 'receipts/Plant_one.json'
    value = json.loads(receipt.read_text())
    assert value['job_id'] == value['recovery_job_id'] == '777' and value['source_job_id'] == '123'
    assert value['recovery_kind'] == 'native-read-cache-verification-v1'
    if stale_atime:
        info = receipt.stat()
        os.utime(receipt, ns=(info.st_mtime_ns - 86_400_000_000_000, info.st_mtime_ns))
    before = file_snapshot(receipt)
    prep.main()
    assert before == file_snapshot(receipt)
    original = {key: item for key, item in value.items()
                if key not in {'source_job_id', 'recovery_job_id', 'recovery_kind', 'recovery_helper_sha256'}}
    original['job_id'] = '123'
    receipt.write_text(json.dumps(original))
    if stale_atime:
        info = receipt.stat()
        os.utime(receipt, ns=(info.st_mtime_ns - 86_400_000_000_000, info.st_mtime_ns))
    before = file_snapshot(receipt)
    prep.main()
    assert before == file_snapshot(receipt)


def test_dispatch_preserves_owning_transcriptome_release_through_real_bootstrap(tmp_path, monkeypatch):
    manifest, _ = fixture(tmp_path)
    native, env = prep.native_environment(manifest, generate=True)
    root = Path(manifest['runtime_release_root'])
    original = subprocess.run
    calls = []
    def boundary(argv, *, cwd, env, check):
        assert cwd == root
        assert env['KFAUTO_RUNTIME_RELEASE_ROOT'] == str(root)
        command = ('source "$1"; gg_validate_explicit_runtime_release "$KFAUTO_RUNTIME_RELEASE_ROOT" '
                   'gg_transcriptome_generation "$1" "$KFAUTO_RUNTIME_RELEASE_SNAPSHOT_SHA256"')
        result = original(['bash', '-c', command, 'bootstrap', str(root / 'workflow/support/gg_entrypoint_bootstrap.sh')],
                          env=env, capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        calls.append(argv)
    monkeypatch.setattr(prep.subprocess, 'run', boundary)
    prep.entrypoint(native, 'gg_transcriptome_generation', env)
    assert len(calls) == 1


def test_cached_assembly_never_downloads_metadata_reads_or_runs_quantification(tmp_path, monkeypatch):
    manifest, row = fixture(tmp_path)
    manifest.update(input_mode=prep.ASSEMBLY_MODE, cache_parameters={'amalgkit_rrna_filter': 'no', 'amalgkit_contam_filter': 'no'})
    calls = []
    def entrypoint(native, name, env):
        assert env['GG_COMMON_TMP_ROOT'] == 'auto'
        assert env['GG_TRANSCRIPTOME_RUN_AMALGKIT_GETFASTQ'] == '0'
        assert env['GG_TRANSCRIPTOME_RUN_AMALGKIT_METADATA_OR_INTEGRATE'] == '0'
        assert env['GG_TRANSCRIPTOME_RUN_AMALGKIT_QUANT'] == env['GG_TRANSCRIPTOME_RUN_AMALGKIT_MERGE'] == '0'
        assert env['GG_TRANSCRIPTOME_AMALGKIT_RRNA_FILTER'] == 'no'
        assert all(env['GG_TRANSCRIPTOME_RUN_' + flag] == '1' for flag in prep.ASSEMBLY_FLAGS)
        calls.append(env)
    monkeypatch.setattr(prep, 'entrypoint', entrypoint)
    monkeypatch.setattr(prep, 'verify_native', lambda *args: dict(verified=True))
    completion = Path(manifest['fastq_source_root']) / row['species'] / 'getfastq_completion.json'
    completion.parent.mkdir()
    completion.write_text('native completion fixture')
    assert prep.prepare(manifest, row)['assembly_complete'] is True
    assert len(calls) == 1


def test_cached_preflight_rechecks_the_native_index_and_verified_handoff_identity(tmp_path):
    import gzip
    manifest, row = fixture(tmp_path)
    cache = Path(manifest['fastq_source_root']) / row['species']
    cache.mkdir()
    entries = []
    for run in row['runs']:
        directory = cache / run
        directory.mkdir()
        fastq = directory / (run + '.amalgkit.fastq.gz')
        with gzip.open(fastq, 'wb') as stream:
            stream.write(b'@read\nACGT\n+\n!!!!\n')
        stats = directory / 'getfastq_stats.tsv'
        stats.write_text('run\tnum_written\tbp_fastp_in\n' + run + '\t1\t4\n')
        def contract(path):
            return dict(path=str(path.relative_to(cache)), size=path.stat().st_size, sha256=prep.file_sha(path))
        entries.append(dict(run=run, status='complete', files=[contract(fastq)], stats=contract(stats)))
    completion = cache / 'getfastq_completion.json'
    completion.write_text(json.dumps(dict(schema_version=4, read_source='amalgkit', status='complete', run_count=3, runs=entries)))
    contract = cache / '.getfastq_cache_contract.json'
    contract.write_text('native cache parameter proof')
    value = dict(path=row['files'][0]['path'], sha256=row['files'][0]['sha256'], runs=row['runs'], cache=str(cache),
                 expected=dict(completion_sha256=prep.file_sha(completion), contract_sha256=prep.file_sha(contract)))
    core = str(SUPPORT.parent / 'core/gg_transcriptome_generation_core.sh')
    command = ['bash', '-c', 'GG_CORE_SOURCE_ONLY=1 source "$3"; validate_amalgkit_getfastq_completion_manifest_index "$1" "$2" yes',
               'fixture', str(completion), value['path'], core]
    value['expected']['content_stat_fingerprint'] = json.loads(subprocess.check_output(command, text=True))
    assert prep.inspect_cache(value)['prepared_cache'] == value['expected']
    content = fastq.read_bytes()
    fastq.write_bytes(content[:-1] + bytes([content[-1] ^ 1]))
    with pytest.raises((subprocess.CalledProcessError, ValueError)):
        prep.inspect_cache(value)


def test_cached_assembly_cannot_claim_completion_without_all_native_proofs(tmp_path, monkeypatch):
    manifest, row = fixture(tmp_path)
    manifest.update(input_mode=prep.ASSEMBLY_MODE, cache_parameters={})
    monkeypatch.setattr(prep, 'entrypoint', lambda *args: None)
    with pytest.raises((FileNotFoundError, ValueError)):
        prep.prepare(manifest, row)

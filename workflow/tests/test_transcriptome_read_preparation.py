"""The read-only preparation path cannot start assembly or claim workflow completion."""
import importlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / 'support'
sys.path.insert(0, str(SUPPORT))
prep = importlib.import_module('transcriptome_read_preparation')


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
    return manifest, row


def test_inspection_is_exact_and_does_not_write_metadata(tmp_path):
    manifest, row = fixture(tmp_path)
    path = Path(row['files'][0]['path'])
    before = path.stat(), path.read_bytes()
    assert prep.inspect_metadata(path, row['files'][0]['sha256'], row['runs'])['run_ids'] == row['runs']
    assert before == (path.stat(), path.read_bytes())
    with pytest.raises(ValueError, match='membership'):
        prep.inspect_metadata(path, row['files'][0]['sha256'], row['runs'][:2])
    with pytest.raises(ValueError, match='changed'):
        prep.inspect_metadata(path, '0' * 64, row['runs'])


def test_actual_inspection_cli_reads_the_selected_metadata(tmp_path):
    _, row = fixture(tmp_path)
    record = row['files'][0]
    result = subprocess.run([sys.executable, str(SUPPORT / 'transcriptome_read_preparation.py'), '--inspect-json',
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

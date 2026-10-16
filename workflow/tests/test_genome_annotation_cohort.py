"""Frozen per-species handoff, native provenance and CPU-only completeness gates."""
import csv
import hashlib
import http.server
import importlib
import json
import subprocess
import sys
import threading
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / 'support'
sys.path.insert(0, str(SUPPORT))
cohort = importlib.import_module('genome_annotation_cohort')
provenance = importlib.import_module('artifact_provenance')


def fixture(tmp_path):
    project = tmp_path.resolve()
    for directory in ('work', 'assets', 'runtime'):
        (project / directory).mkdir()
    (project / 'image.sif').write_bytes(b'fixture')
    genome = project / 'source.fa'
    genome.write_text('>contig\nATGAAATAA\n')
    member = dict(task=1, species='Plant_one', filename='Plant_one.genome.fa',
                  source=dict(local_path=str(genome), sha256=cohort.file_sha(genome)))
    manifest = dict(schema_version=1, busco_lineage='embryophyta_odb12',
                    project_root=str(project), work_root=str(project / 'work'),
                    downloads_root=str(project / 'assets'), runtime_root=str(project / 'runtime'),
                    image_path=str(project / 'image.sif'), members=[member])
    workspace = project / 'work/Plant_one'
    workspace.mkdir()
    return manifest, member, workspace


def record(path, workspace, step, inputs, outputs, parameters=None):
    logical = workspace / 'output/.gg_global_artifacts'
    payload = dict(schema_version=1, step=step, family_id='Plant_one',
                   inputs=provenance.describe_paths(inputs, logical, workspace),
                   outputs=provenance.describe_paths(outputs, logical, workspace),
                   optional_outputs=[], parameters=parameters or {})
    cohort.atomic_json(path, payload)


def prediction(workspace):
    root = workspace / 'output/genome_prediction/oriongeno/Plant_one'
    formatted = root / 'formatted/generation'
    for role in ('cds', 'gff', 'genome'):
        directory = formatted / ('species_' + role)
        directory.mkdir(parents=True)
        filename = 'Plant_one.gff3' if role == 'gff' else 'Plant_one.' + role + '.fa'
        (directory / filename).write_text('formatted fixture ' + role)
        (directory / (filename + '.audit.json')).write_text('formatter audit fixture')
    gff = root / 'Plant_one.gff3'
    gff.write_text('native fixture')
    genome = workspace / 'input/species_genome/Plant_one.genome.fa'
    record(root / 'prediction.provenance.json', workspace, 'genome_annotation_oriongeno_prediction',
           [('genome', str(genome))], [('gff', str(gff))])
    record(root / 'format.provenance.json', workspace, 'genome_annotation_oriongeno_format',
           [('genome', str(genome)), ('gff', str(gff))], [('formatted', str(formatted))])
    return root


def test_unique_frozen_membership_and_exact_lineage(tmp_path):
    manifest, member, _ = fixture(tmp_path)
    assert cohort.select(manifest, 1) == member
    for bad in ({'busco_lineage': 'eukaryota_odb12'}, {'members': [member, member]},
                {'members': [{**member, 'task': 2}]}, {'members': [{**member, 'species': '../escape'}]},
                {'members': [{**member, 'filename': 'Other.genome.fa'}]}):
        with pytest.raises(ValueError):
            cohort.select({**manifest, **bad}, 1)
    with pytest.raises(ValueError):
        cohort.select(manifest, 2)


def test_local_acquisition_is_physical_reusable_and_hash_bound(tmp_path):
    manifest, member, workspace = fixture(tmp_path)
    cohort.acquire(manifest, member, workspace)
    target = workspace / 'input/species_genome/Plant_one.genome.fa'
    assert target.samefile(member['source']['local_path'])
    assert not target.is_symlink()
    first = (workspace / 'genome.acquired.json').read_bytes()
    cohort.acquire(manifest, member, workspace)
    assert first == (workspace / 'genome.acquired.json').read_bytes()
    Path(member['source']['local_path']).write_text('changed')
    with pytest.raises(ValueError, match='changed'):
        cohort.acquire(manifest, member, workspace)


def test_public_acquisition_calls_native_download_only_and_rejects_bad_md5(tmp_path, monkeypatch):
    manifest, _, workspace = fixture(tmp_path)
    member = dict(task=1, species='Plant_one', filename='Plant_one.genome.fa',
                  source=dict(url='https://example.invalid/genome.fa', md5='0' * 32))
    calls = []
    def native(manifest, name, env):
        calls.append((name, env))
        rows = list(csv.DictReader(open(env['GG_INPUT_DOWNLOAD_MANIFEST']), delimiter='\t'))
        assert rows == [dict(provider='direct', id='Plant_one', species_key='Plant_one',
                             genome_url=member['source']['url'], genome_filename=member['filename'], genome_archive_member='')]
        root = cohort.provider_raw_dir('direct', workspace / 'raw', member['species'])
        root.mkdir(parents=True)
        (root / member['filename']).write_text('download fixture')
    monkeypatch.setattr(cohort, 'entrypoint', native)
    with pytest.raises(ValueError, match='published MD5'):
        cohort.acquire(manifest, member, workspace)
    assert not (workspace / 'genome.acquired.json').exists()
    assert not (workspace / 'input/species_genome/Plant_one.genome.fa').exists()
    assert calls[0][0] == 'gg_input_generation'
    assert calls[0][1]['GG_INPUT_DOWNLOAD_ONLY'] == '1'
    assert calls[0][1]['GG_INPUT_RUN_SPECIES_BUSCO'] == '0'


def test_native_public_download_and_warm_reuse_preserve_genome_bytes(tmp_path, monkeypatch):
    manifest, member, workspace = fixture(tmp_path)
    data = b'>contig\nATGAAATAA\n'
    requests = []
    class Handler(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            requests.append(self.path)
            self.send_response(200)
            self.send_header('Content-Length', str(len(data)))
            self.end_headers()
            self.wfile.write(data)

        def log_message(self, *args):
            pass
    server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    member['source'] = dict(url=f'http://127.0.0.1:{server.server_port}/genome.fa',
                            md5=hashlib.md5(data, usedforsecurity=False).hexdigest())
    def native(manifest, name, env):
        subprocess.run([sys.executable, str(SUPPORT / 'format_species_inputs.py'), '--provider', 'direct',
                        '--download-only', '--download-manifest', env['GG_INPUT_DOWNLOAD_MANIFEST'],
                        '--resolved-manifest-output', str(workspace / 'resolved.tsv'),
                        '--download-dir', env['GG_INPUT_DOWNLOAD_DIR'], '--jobs', '1'], env=env, check=True)
    monkeypatch.setattr(cohort, 'entrypoint', native)
    try:
        cohort.acquire(manifest, member, workspace)
        first = (workspace / 'genome.acquired.json').read_bytes()
        count = len(requests)
        cohort.acquire(manifest, member, workspace)
        assert first == (workspace / 'genome.acquired.json').read_bytes()
        assert len(requests) == count == 1
        assert (workspace / 'input/species_genome/Plant_one.genome.fa').read_bytes() == data
    finally:
        server.shutdown()
        thread.join()


@pytest.mark.parametrize('damage', ['missing', 'changed', 'wrong_lineage', 'none'])
def test_busco_requires_current_native_prediction_and_successful_native_proof(tmp_path, monkeypatch, damage):
    manifest, member, workspace = fixture(tmp_path)
    cohort.acquire(manifest, member, workspace)
    root = prediction(workspace)
    calls = []
    def native(manifest, name, env):
        calls.append(env)
        assert env['GG_GENOME_ANNOTATION_GENOME_ANNOTATION_BACKEND'] == 'none'
        assert env['GG_GENOME_ANNOTATION_RUN_BUSCO_CDS'] == '1'
        assert env['GG_COMMON_BUSCO_LINEAGE'] == 'embryophyta_odb12'
        assert Path(env['gg_workspace_dir']) == workspace / 'busco'
        assert env['GG_ARRAY_TASK_ID'] == '1'
        if damage in {'wrong_lineage', 'none'}:
            phase_workspace = workspace / 'busco'
            output = phase_workspace / 'output/short.txt'
            output.parent.mkdir(parents=True)
            output.write_text('BUSCO fixture')
            record(phase_workspace / 'output/artifact_provenance/genome_annotation/Plant_one.busco_cds.json',
                   phase_workspace, 'genome_annotation_busco_cds', [], [('short_summary', str(output))],
                   {'lineage': 'eukaryota_odb12' if damage == 'wrong_lineage' else 'embryophyta_odb12'})
    monkeypatch.setattr(cohort, 'entrypoint', native)
    if damage == 'changed':
        (root / 'Plant_one.gff3').write_text('tampered')
    if damage == 'none':
        cohort.annotate(manifest, member, workspace, 'busco')
    else:
        with pytest.raises((ValueError, FileNotFoundError, provenance.ProvenanceError)):
            cohort.annotate(manifest, member, workspace, 'busco')
    assert bool(calls) == (damage != 'changed')


def test_wrong_manifest_hash_rejected_without_creating_workspaces(tmp_path):
    manifest, _, workspace = fixture(tmp_path)
    path = tmp_path / 'cohort.json'
    path.write_text(json.dumps(manifest))
    workspace.rmdir()
    result = subprocess.run([sys.executable, str(SUPPORT / 'genome_annotation_cohort.py'), '--manifest', str(path),
                             '--manifest-sha256', '0' * 64, '--phase', 'acquire', '--task', '1'], capture_output=True, text=True)
    assert result.returncode != 0
    assert 'manifest changed' in result.stderr
    assert not workspace.exists()

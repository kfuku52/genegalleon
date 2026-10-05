import csv
import gzip
import json
import subprocess
import sys
import time
from pathlib import Path

import pytest

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
PLAN_SCRIPT = SUPPORT_DIR / "plan_input_generation_tasks.py"
RUN_TASK_SCRIPT = SUPPORT_DIR / "run_input_generation_task.py"
MERGE_SCRIPT = SUPPORT_DIR / "merge_input_generation_shards.py"
STAGE_SCRIPT = SUPPORT_DIR / "stage_input_generation_downloads.py"
REQUIRE_GENOMES_SCRIPT = SUPPORT_DIR / "validate_required_genomes.py"
REQUIRE_OUTPUTS_SCRIPT = SUPPORT_DIR / "validate_required_species_outputs.py"


@pytest.fixture
def bound_staging_plan(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(SUPPORT_DIR))
    species = 'Arabidopsis_thaliana'
    raw = tmp_path / 'raw'
    write_direct_species_fixture(raw, species)
    roles = {role: raw / species / (species + suffix) for role, suffix in
             (('cds', '.cds.fa'), ('gff', '.gff'), ('genome', '.genome.fa'))}
    manifest = tmp_path / 'manifest.tsv'
    manifest.write_text('provider\tid\tspecies_key\tbind_local_sources\tcds_url\tgff_url\tgenome_url\n'
                        + f'direct\tfixture\t{species}\t1\t'
                        + '\t'.join(path.as_uri() for path in roles.values()) + '\n')
    plan = tmp_path / 'plan.json'
    result = run_python(PLAN_SCRIPT, '--provider', 'all', '--download-manifest', str(manifest),
                        '--download-dir', str(tmp_path / 'downloads'), '--stage-downloads', '--outfile', str(plan))
    assert result.returncode == 0, result.stderr
    return plan, roles


def test_staging_reads_unique_sources_per_boundary_and_resume_keeps_receipts(bound_staging_plan, monkeypatch):
    import input_generation_array_state as state
    import stage_input_generation_downloads as staging

    plan, roles = bound_staging_plan
    original_digest = state.digest
    calls = []

    def counted(path):
        calls.append(str(path))
        return original_digest(path)

    monkeypatch.setattr(state, 'digest', counted)
    monkeypatch.setattr(staging, 'digest', counted)
    staging.stage_downloads(plan)
    for path in roles.values():
        assert calls.count(str(path)) == 2  # Fresh preflight plus independent publication.
    receipt = Path(str(plan) + '.tasks/1.json')
    frozen = receipt.read_bytes()
    calls.clear()
    staging.stage_downloads(plan)
    assert receipt.read_bytes() == frozen
    for path in roles.values():
        assert calls.count(str(path)) == 1
    payload = json.loads(frozen)
    payload['task']['input_sha256'][str(roles['cds'])] = '0' * 64
    receipt.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match='Staged raw input changed'):
        staging.stage_downloads(plan)


@pytest.fixture
def reusable_staging_plan(bound_staging_plan):
    import input_generation_array_state as state
    import stage_input_generation_downloads as staging

    old_plan, roles = bound_staging_plan
    staging.stage_downloads(old_plan)
    task = json.loads(old_plan.read_text())['tasks'][0]
    row = dict(task['manifest_row'])
    row.update(reuse_staged_plan=str(old_plan), reuse_staged_plan_sha256=state.digest(old_plan),
               reuse_staged_workspace=str(old_plan.parent), reuse_staged_task_index='1',
               reuse_staged_receipt_sha256=state.digest(Path(str(old_plan) + '.tasks/1.json')))
    manifest = old_plan.parent / 'reuse.tsv'
    with manifest.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter='\t')
        writer.writeheader()
        writer.writerow(row)
    new_plan = old_plan.with_name('new.json')
    return old_plan, new_plan, manifest, roles


def plan_reuse(new_plan, manifest):
    return run_python(PLAN_SCRIPT, '--provider', 'all', '--download-manifest', str(manifest),
                      '--download-dir', str(new_plan.parent / 'new-downloads'),
                      '--stage-downloads', '--outfile', str(new_plan))


def test_cross_plan_reuse_reads_once_keeps_outputs_and_same_plan_resume(reusable_staging_plan, monkeypatch):
    import input_generation_array_state as state
    import plan_input_generation_tasks as planner
    import stage_input_generation_downloads as staging

    old_plan, new_plan, manifest, roles = reusable_staging_plan
    original_digest = state.digest
    calls = []

    def counted(path):
        calls.append(str(path))
        return original_digest(path)

    monkeypatch.setattr(state, 'digest', counted)
    monkeypatch.setattr(sys, 'argv', ['plan', '--provider', 'all', '--download-manifest', str(manifest),
                                    '--download-dir', str(new_plan.parent / 'new-downloads'),
                                    '--stage-downloads', '--outfile', str(new_plan)])
    assert planner.main() == 0
    assert not any(str(path) in calls for path in roles.values())

    def unexpected_gzip(*args, **kwargs):
        pytest.fail('Identical successfully staged inputs should reuse gzip validation')

    monkeypatch.setattr(staging, 'validate_gzip_with_cache', unexpected_gzip)
    staging.stage_downloads(new_plan)
    for path in roles.values():
        assert calls.count(str(path)) == 1
    original = json.loads(Path(str(old_plan) + '.tasks/1.json').read_text())['task']
    migrated = json.loads(Path(str(new_plan) + '.tasks/1.json').read_text())['task']
    assert {k: v for k, v in migrated.items() if k != 'staged_input_reuse'} == original
    frozen = Path(str(new_plan) + '.tasks/1.json').read_bytes()
    calls.clear()
    staging.stage_downloads(new_plan)
    assert Path(str(new_plan) + '.tasks/1.json').read_bytes() == frozen
    for path in roles.values():
        assert calls.count(str(path)) == 1


@pytest.mark.parametrize('fault', ['plan_sha', 'receipt_sha', 'species', 'index', 'role', 'partial', 'unbound'])
def test_cross_plan_reuse_rejects_wrong_proof(reusable_staging_plan, fault):
    old_plan, new_plan, manifest, roles = reusable_staging_plan
    with manifest.open() as handle:
        row = next(csv.DictReader(handle, delimiter='\t'))
    if fault == 'plan_sha':
        row['reuse_staged_plan_sha256'] = '0' * 64
    elif fault == 'receipt_sha':
        row['reuse_staged_receipt_sha256'] = '0' * 64
    elif fault == 'species':
        row['species_key'] = 'Arabidopsis_lyrata'
    elif fault == 'index':
        row['reuse_staged_task_index'] = '2'
    elif fault == 'role':
        row['cds_url'] = roles['genome'].as_uri()
    elif fault == 'partial':
        row['reuse_staged_workspace'] = ''
    else:
        row['bind_local_sources'] = '0'
    with manifest.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter='\t')
        writer.writeheader()
        writer.writerow(row)
    result = plan_reuse(new_plan, manifest)
    assert result.returncode != 0
    assert 'Staging reuse' in result.stderr or 'staging reuse' in result.stderr
    assert not new_plan.exists()


@pytest.mark.parametrize('moment', ['before_staging', 'binding', 'publication', 'after_publication'])
@pytest.mark.parametrize('mutation', ['restore_mtime', 'replace'])
def test_cross_plan_reuse_rejects_changed_bytes_and_late_changes(reusable_staging_plan, monkeypatch, moment, mutation):
    import os

    import stage_input_generation_downloads as staging

    _old_plan, new_plan, manifest, roles = reusable_staging_plan
    assert plan_reuse(new_plan, manifest).returncode == 0
    source = roles['cds']

    def mutate():
        before = source.stat()
        time.sleep(0.02)
        data = source.read_bytes().replace(b'ATG', b'ACG')
        if mutation == 'replace':
            replacement = source.with_suffix('.new')
            replacement.write_bytes(data)
            os.replace(replacement, source)
        else:
            source.write_bytes(data)
        os.utime(source, ns=(before.st_atime_ns, before.st_mtime_ns))

    if moment == 'before_staging':
        mutate()
    elif moment in ('binding', 'publication'):
        original = staging.bound_local_manifest_task

        def binding(task, **kwargs):
            if moment == 'binding':
                mutate()
            result = original(task, **kwargs)
            if moment == 'publication':
                mutate()
            return result

        monkeypatch.setattr(staging, 'bound_local_manifest_task', binding)
    else:
        original = staging.atomic_json

        def publish(path, *args, **kwargs):
            result = original(path, *args, **kwargs)
            mutate()
            return result

        monkeypatch.setattr(staging, 'atomic_json', publish)
    with pytest.raises(ValueError, match='(Input changed during staging reuse|Staging reuse input SHA256 mismatch)'):
        staging.stage_downloads(new_plan)
    assert not Path(str(new_plan) + '.prepared.json').exists()


def test_cross_plan_reuse_rechecks_sealed_evidence_on_staging(reusable_staging_plan):
    import stage_input_generation_downloads as staging

    old_plan, new_plan, manifest, _roles = reusable_staging_plan
    assert plan_reuse(new_plan, manifest).returncode == 0
    receipt = Path(str(old_plan) + '.tasks/1.json')
    receipt.write_text(receipt.read_text() + ' ')
    with pytest.raises(ValueError, match='evidence SHA256 mismatch'):
        staging.stage_downloads(new_plan)


def test_cross_plan_reuse_maps_old_container_paths(reusable_staging_plan):
    import input_generation_array_state as state
    import stage_input_generation_downloads as staging

    old_plan, new_plan, manifest, _roles = reusable_staging_plan
    receipt_path = Path(str(old_plan) + '.tasks/1.json')
    receipt = json.loads(receipt_path.read_text())
    task = receipt['task']
    mapping = {path: '/workspace/' + str(Path(path).relative_to(old_plan.parent)) for path in task['input_sha256']}
    task['input_sha256'] = {mapping[path]: sha for path, sha in task['input_sha256'].items()}
    for role in ('cds', 'gff', 'genome'):
        task[role + '_path'] = mapping[task[role + '_path']]
    receipt_path.write_text(json.dumps(receipt))
    text = manifest.read_text().replace(next(csv.DictReader(manifest.open(), delimiter='\t'))['reuse_staged_receipt_sha256'], state.digest(receipt_path))
    manifest.write_text(text)
    result = plan_reuse(new_plan, manifest)
    assert result.returncode == 0, result.stderr
    staging.stage_downloads(new_plan)


def test_cross_plan_reuse_keeps_worker_fresh_content_verification(reusable_staging_plan):
    from argparse import Namespace

    import run_input_generation_task as worker
    import stage_input_generation_downloads as staging

    _old_plan, new_plan, manifest, roles = reusable_staging_plan
    assert plan_reuse(new_plan, manifest).returncode == 0
    staging.stage_downloads(new_plan)
    task = json.loads(new_plan.read_text())['tasks'][0]
    args = Namespace(task_plan=str(new_plan), task_index=1)
    assert worker.resolve_manifest_task(task, args)['cds_path'] == roles['cds']
    roles['cds'].write_bytes(roles['cds'].read_bytes().replace(b'ATG', b'ACG'))
    with pytest.raises(ValueError, match='Local manifest input changed'):
        worker.resolve_manifest_task(task, args)


def test_cross_plan_reuse_keeps_current_coge_validation(reusable_staging_plan, monkeypatch):
    import input_generation_array_state as state
    import stage_input_generation_downloads as staging

    old_plan, new_plan, manifest, _roles = reusable_staging_plan
    old = json.loads(old_plan.read_text())
    old['tasks'][0]['provider'] = 'coge'
    old_plan.write_text(json.dumps(old))
    receipt_path = Path(str(old_plan) + '.tasks/1.json')
    receipt = json.loads(receipt_path.read_text())
    receipt['task']['provider'] = 'coge'
    receipt['plan_sha256'] = state.digest(old_plan)
    receipt_path.write_text(json.dumps(receipt))
    with manifest.open() as handle:
        row = next(csv.DictReader(handle, delimiter='\t'))
    row.update(provider='coge', reuse_staged_plan_sha256=state.digest(old_plan),
               reuse_staged_receipt_sha256=state.digest(receipt_path))
    with manifest.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter='\t')
        writer.writeheader()
        writer.writerow(row)
    result = plan_reuse(new_plan, manifest)
    assert result.returncode == 0, result.stderr

    def current_check(*args, **kwargs):
        raise ValueError('Current CoGe payload check')

    monkeypatch.setattr(staging, 'validate_coge_export_gff_file', current_check)
    with pytest.raises(ValueError, match='Current CoGe payload check'):
        staging.stage_downloads(new_plan)
    assert not Path(str(new_plan) + '.tasks/1.json').exists()


@pytest.mark.parametrize('boundary', ['binding', 'publication', 'resume'])
@pytest.mark.parametrize('mutation', ['write', 'restore_mtime', 'replace'])
def test_staging_rejects_changes_between_boundaries(bound_staging_plan, monkeypatch, boundary, mutation):
    import os

    import stage_input_generation_downloads as staging

    plan, roles = bound_staging_plan
    source = roles['cds']

    def mutate():
        before = source.stat()
        changed = source.read_bytes().replace(b'ATG', b'ACG')
        if mutation == 'replace':
            replacement = source.with_suffix('.replacement')
            replacement.write_bytes(changed)
            os.replace(replacement, source)
        else:
            source.write_bytes(changed)
        if mutation in ('restore_mtime', 'replace'):
            os.utime(source, ns=(before.st_atime_ns, before.st_mtime_ns))

    if boundary == 'resume':
        staging.stage_downloads(plan)
        mutate()
        with pytest.raises(ValueError, match='Local manifest input changed'):
            staging.stage_downloads(plan)
    else:
        original_binding = staging.bound_local_manifest_task

        def changed_binding(task, **kwargs):
            if boundary == 'binding':
                mutate()
            actual = original_binding(task, **kwargs)
            if boundary == 'publication':
                mutate()
            return actual

        monkeypatch.setattr(staging, 'bound_local_manifest_task', changed_binding)
        with pytest.raises(ValueError, match='(Bound local source|Local input changed|Input changed during staging reuse)'):
            staging.stage_downloads(plan)
        assert not Path(str(plan) + '.tasks/1.json').exists()


def test_planning_and_binding_hash_aliased_roles_once(bound_staging_plan, monkeypatch):
    import input_generation_array_state as state
    import plan_input_generation_tasks as planner
    import stage_input_generation_downloads as staging

    plan, roles = bound_staging_plan
    manifest = plan.parent / 'manifest.tsv'
    manifest.write_text(manifest.read_text().replace(roles['cds'].as_uri(), roles['genome'].as_uri()))
    aliased_plan = plan.with_name('aliased.json')
    original_digest = state.digest
    calls = []

    def counted(path):
        calls.append(str(path))
        return original_digest(path)

    monkeypatch.setattr(state, 'digest', counted)
    monkeypatch.setattr(sys, 'argv', ['planner', '--provider', 'all', '--download-manifest', str(manifest),
                        '--download-dir', str(plan.parent / 'downloads'), '--stage-downloads', '--outfile', str(aliased_plan)])
    assert planner.main() == 0
    assert calls.count(str(roles['genome'])) == 1
    calls.clear()
    task = json.loads(aliased_plan.read_text())['tasks'][0]
    actual = staging.bound_local_manifest_task(task)
    assert actual['cds_path'] == actual['genome_path'] == roles['genome']
    assert calls.count(str(roles['genome'])) == 1
    calls.clear()
    staging.stage_downloads(aliased_plan)
    assert calls.count(str(roles['genome'])) == 2


def test_digest_rejects_path_replacement_while_reading(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(SUPPORT_DIR))
    import input_generation_array_state as state

    source = tmp_path / 'source.fa'
    source.write_bytes(b'>chr1\nATG\n')
    replacement = tmp_path / 'replacement.fa'
    replacement.write_bytes(b'>chr1\nACG\n')
    original_fstat = state.os.fstat
    calls = 0

    def replace_after_fstat(fd):
        nonlocal calls
        info = original_fstat(fd)
        calls += 1
        if calls == 2:
            state.os.replace(replacement, source)
        return info

    monkeypatch.setattr(state.os, 'fstat', replace_after_fstat)
    with pytest.raises(OSError, match='File changed while hashing'):
        state.digest(source)


def test_unique_hash_batch_rejects_source_changed_after_its_read(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(SUPPORT_DIR))
    import input_generation_array_state as state

    first, second = tmp_path / 'first.fa', tmp_path / 'second.fa'
    first.write_bytes(b'>chr1\nATG\n')
    second.write_bytes(b'>chr2\nATG\n')
    original_digest = state.digest

    def change_earlier_file(path):
        if path == str(second):
            before = first.stat()
            # Lustre may expose second-resolution ctime. Cross that boundary
            # so this test specifically exercises the ctime fence.
            time.sleep(1.1)
            first.write_bytes(b'>chr1\nACG\n')
            state.os.utime(first, ns=(before.st_atime_ns, before.st_mtime_ns))
        return original_digest(path)

    monkeypatch.setattr(state, 'digest', change_earlier_file)
    with pytest.raises(OSError, match='File changed while hashing'):
        state.digest_paths([first, second, first])


def test_worker_verifies_bound_inputs_once_and_rejects_conflicting_receipt(bound_staging_plan, monkeypatch):
    import input_generation_array_state as state
    import run_input_generation_task as runner
    import stage_input_generation_downloads as staging

    plan, roles = bound_staging_plan
    staging.stage_downloads(plan)
    original_digest = state.digest
    calls = []

    def counted(path):
        calls.append(str(path))
        return original_digest(path)

    monkeypatch.setattr(state, 'digest', counted)
    monkeypatch.setattr(runner, 'digest', counted)
    metadata = plan.parent / 'metadata.json'
    argv = ['runner', '--task-plan', str(plan), '--task-index', '1', '--describe-only',
            '--species-cds-dir', str(plan.parent / 'cds'), '--species-gff-dir', str(plan.parent / 'gff'),
            '--species-genome-dir', str(plan.parent / 'genomes'), '--task-meta-output', str(metadata)]
    monkeypatch.setattr(sys, 'argv', argv)
    assert runner.main() == 0
    for path in roles.values():
        assert calls.count(str(path)) == 1
    frozen = metadata.read_bytes()
    receipt = Path(str(plan) + '.tasks/1.json')
    payload = json.loads(receipt.read_text())
    payload['task']['input_sha256'][str(roles['cds'])] = '0' * 64
    receipt.write_text(json.dumps(payload))
    with pytest.raises(SystemExit) as error:
        runner.main()
    assert error.value.code == 2
    assert metadata.read_bytes() == frozen
    with pytest.raises(ValueError, match='conflicts with the frozen plan'):
        state.frozen_input_hashes(plan, state.load_plan(plan), 1)


def test_legacy_download_rejects_original_changed_before_receipt(bound_staging_plan, monkeypatch):
    from types import SimpleNamespace

    import run_input_generation_task as runner

    plan, roles = bound_staging_plan
    payload = json.loads(plan.read_text())
    payload.pop('download_mode')
    payload['tasks'][0]['manifest_row']['bind_local_sources'] = '0'
    plan.write_text(json.dumps(payload))
    original_download = runner.fsi.download_from_manifest

    def changed_download(*args, **kwargs):
        report = original_download(*args, **kwargs)
        roles['cds'].write_text('>gene1\nACGAAATTT\n')
        return report

    monkeypatch.setattr(runner.fsi, 'download_from_manifest', changed_download)
    args = SimpleNamespace(task_plan=plan, task_index=1, download_timeout=5,
                           http_header=[], auth_bearer_token_env='')
    with pytest.raises(ValueError, match='Local manifest input changed'):
        runner.resolve_manifest_task(runner.deserialize_task(payload['tasks'][0]), args)
    assert not Path(str(plan) + '.tasks/1.json').exists()


def test_bound_sources_preserve_formatted_outputs_without_raw_copy_and_reject_changes(tmp_path):
    raw = tmp_path / "raw"
    species = "Arabidopsis_thaliana"
    write_direct_species_fixture(raw, species)
    outputs = []
    for mode in ("0", "1"):
        root = tmp_path / mode
        root.mkdir()
        manifest = root / "manifest.tsv"
        roles = {role: raw / species / (species + suffix) for role, suffix in
                 (("cds", ".cds.fa"), ("gff", ".gff"), ("genome", ".genome.fa"))}
        fields = ["provider", "id", "species_key", "bind_local_sources", *[r + "_url" for r in roles]]
        with manifest.open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerow({"provider": "direct", "id": "fixture", "species_key": species,
                             "bind_local_sources": mode, **{r + "_url": p.as_uri() for r, p in roles.items()}})
        plan = root / "plan.json"
        planned = run_python(PLAN_SCRIPT, "--provider", "all", "--download-manifest", str(manifest),
                             "--download-dir", str(root / "downloads"), "--stage-downloads", "--outfile", str(plan))
        assert planned.returncode == 0, planned.stderr
        staged = run_python(STAGE_SCRIPT, "--task-plan", str(plan), "--require-gff", "--require-genome")
        assert staged.returncode == 0, staged.stdout + staged.stderr
        task = json.loads(Path(str(plan) + ".tasks/1.json").read_text())["task"]
        if mode == "1":
            assert all(task[r + "_path"] == str(p) for r, p in roles.items())
            assert not (root / "downloads").exists()
        else:
            assert all(task[r + "_path"] != str(p) for r, p in roles.items())
        args = ("--task-plan", str(plan), "--task-index", "1", "--species-cds-dir", str(root / "cds"),
                "--species-gff-dir", str(root / "gff"), "--species-genome-dir", str(root / "genome"))
        result = run_python(RUN_TASK_SCRIPT, *args)
        assert result.returncode == 0, result.stdout + result.stderr
        outputs.append({role: [(p.name, gzip.open(p, "rt").read()) for p in sorted((root / role).glob("*.gz"))]
                        for role in roles})
        assert all(outputs[-1].values())
    assert outputs[0] == outputs[1]
    roles["cds"].write_text(">changed\nATG\n")
    rejected = run_python(RUN_TASK_SCRIPT, *args)
    assert rejected.returncode != 0 and "Local manifest input changed" in rejected.stderr
    restaged = run_python(STAGE_SCRIPT, "--task-plan", str(plan))
    assert restaged.returncode != 0 and "Local manifest input changed" in restaged.stderr


@pytest.mark.parametrize('invalid', ['html', 'header_only', 'corrupt_gzip', 'remote_role', 'archive_member'])
def test_bound_coge_sources_do_not_bypass_payload_or_source_guards(tmp_path, invalid):
    sys.path.insert(0, str(SUPPORT_DIR))
    from input_generation_array_state import digest
    from stage_input_generation_downloads import bound_local_manifest_task

    cds, gff, genome = (tmp_path / name for name in ('cds.fa', 'source.gff', 'genome.fa'))
    cds.write_text('>gene1\nATG\n')
    genome.write_text('>chr1\nATG\n')
    gff.write_text('##gff-version 3\nchr1\tCoGe\tCDS\t1\t3\t.\t+\t0\tID=c;Parent=t\n')
    if invalid == 'html':
        gff.write_text('<html>failure</html>')
    elif invalid == 'header_only':
        gff.write_text('##gff-version 3\n')
    elif invalid == 'corrupt_gzip':
        gff = tmp_path / 'source.gff.gz'
        gff.write_bytes(gzip.compress(b'##gff-version 3\n')[:-5])
    row = {'provider': 'coge', 'id': '123', 'species_key': 'Example_species', 'bind_local_sources': '1',
           'cds_url': cds.as_uri(), 'gff_url': gff.as_uri(), 'genome_url': genome.as_uri()}
    if invalid == 'remote_role':
        row['genome_url'] = 'https://example.invalid/genome.fa'
    elif invalid == 'archive_member':
        row['gff_archive_member'] = 'source.gff'
    task = {'provider': 'coge', 'species_key': 'Example_species', 'species_prefix': 'Example_species',
            'manifest_row': row, 'input_sha256': {str(p): digest(p) for p in (cds, gff, genome)}}
    with pytest.raises(ValueError):
        bound_local_manifest_task(task)


def test_staged_http_inputs_run_without_server_and_reject_missing_or_changed_cache(tmp_path):
    import functools
    import threading
    from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

    raw = tmp_path / "raw"
    species = "Arabidopsis_thaliana"
    write_direct_species_fixture(raw, species)
    server = ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(SimpleHTTPRequestHandler, directory=str(raw)))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    manifest = tmp_path / "source.tsv"
    base = f"http://127.0.0.1:{server.server_port}/{species}/{species}"
    manifest.write_text("provider\tid\tspecies_key\tcds_url\tgff_url\tgenome_url\n"
                        f"direct\tfixture\t{species}\t{base}.cds.fa\t{base}.gff\t{base}.genome.fa\n")
    plan = tmp_path / "plan.json"
    planned = run_python(PLAN_SCRIPT, "--provider", "all", "--download-manifest", str(manifest),
                         "--download-dir", str(tmp_path / "downloads"), "--stage-downloads", "--outfile", str(plan))
    assert planned.returncode == 0, planned.stderr
    args = ("--task-plan", str(plan), "--task-index", "1", "--species-cds-dir", str(tmp_path / "cds"),
            "--species-gff-dir", str(tmp_path / "gff"), "--species-genome-dir", str(tmp_path / "genome"))
    try:
        missing = run_python(RUN_TASK_SCRIPT, *args)
        assert missing.returncode != 0 and "Staged download receipt is missing" in missing.stderr
        staged = run_python(STAGE_SCRIPT, "--task-plan", str(plan), "--jobs", "3")
        assert staged.returncode == 0, staged.stdout + staged.stderr
    finally:
        server.shutdown()
        server.server_close()
        thread.join(3)
    manifest.write_text("no longer valid\n")
    completed = run_python(RUN_TASK_SCRIPT, *args)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert list((tmp_path / "cds").glob("*.fa.gz"))
    restaged = run_python(STAGE_SCRIPT, "--task-plan", str(plan))
    assert restaged.returncode == 0 and "no downloads needed" in restaged.stdout, restaged.stderr
    cached = json.loads(Path(str(plan) + ".tasks/1.json").read_text())["task"]
    Path(cached["cds_path"]).write_text(">changed\nATG\n")
    rejected = run_python(RUN_TASK_SCRIPT, *args)
    assert rejected.returncode != 0 and "Raw input changed" in rejected.stderr
    rejected_stage = run_python(STAGE_SCRIPT, "--task-plan", str(plan))
    assert rejected_stage.returncode != 0 and "Staged raw input changed" in rejected_stage.stderr


def test_staged_manifest_roles_override_nonstandard_direct_and_ncbi_filenames(tmp_path):
    cases = (
        ("direct", "fixture", "Chrysanthemum_morifolium", "Cmo.fasta"),
        ("ncbi", "GCA_000000001.1", "Nicotiana_benthamiana", "Nbe_scf.fa"),
    )
    for provider, source_id, species, genome_name in cases:
        root = tmp_path / species
        root.mkdir()
        cds = root / "genes.cds.fa"
        gff = root / "genes.gff"
        genome = root / genome_name
        cds.write_text(">gene1\nATGAAATTT\n", encoding="utf-8")
        gff.write_text("##gff-version 3\nchr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1\n", encoding="utf-8")
        genome.write_text(">chr1\nATGAAATTT\n", encoding="utf-8")
        manifest = root / "manifest.tsv"
        manifest.write_text(
            "provider\tid\tspecies_key\tcds_url\tgff_url\tgenome_url\tcds_filename\tgff_filename\tgenome_filename\n"
            f"{provider}\t{source_id}\t{species}\t{cds.as_uri()}\t{gff.as_uri()}\t{genome.as_uri()}\t{cds.name}\t{gff.name}\t{genome.name}\n",
            encoding="utf-8",
        )
        plan = root / "plan.json"
        planned = run_python(PLAN_SCRIPT, "--provider", "all", "--download-manifest", str(manifest),
                             "--download-dir", str(root / "downloads"), "--stage-downloads", "--outfile", str(plan),
                             "--require-gff", "--require-genome")
        assert planned.returncode == 0, planned.stderr
        staged = run_python(STAGE_SCRIPT, "--task-plan", str(plan), "--require-gff", "--require-genome")
        assert staged.returncode == 0, staged.stdout + staged.stderr
        task = json.loads(Path(str(plan) + ".tasks/1.json").read_text())["task"]
        assert Path(task["genome_path"]).name == genome_name
        assert Path(task["cds_path"]).name == cds.name


@pytest.mark.parametrize("provider,species,source_prefix", [
    ("ensemblplants", "Oryza_sativa_tropical_japonica_subgroup", "Oryza_sativa_azucena.AzucenaRS1"),
    ("ensemblplants", "Oryza_sativa_aus_subgroup", "Oryza_sativa_n22.OsN22RS2"),
    ("citrusgenomedb", "Citrus_clementina_x_Citrus_x_tangelo", "Citrus_reticulata_fairchild"),
])
@pytest.mark.parametrize("staged", [True, False])
def test_manifest_roles_preserve_provider_alias_and_hybrid_identity(tmp_path, provider, species, source_prefix, staged):
    roles = {role: tmp_path / (source_prefix + suffix) for role, suffix in
             (("cds", ".cds.fa"), ("gff", ".gff3"), ("genome", ".dna.toplevel.fa"))}
    roles["cds"].write_text(">gene1\nATGAAATTT\n")
    roles["gff"].write_text("##gff-version 3\nchr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1\n")
    roles["genome"].write_text(">chr1\nATGAAATTT\n")
    fields = ["provider", "id", "species_key", *[r + suffix for r in roles for suffix in ("_url", "_filename")]]
    manifest = tmp_path / "manifest.tsv"
    with manifest.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow({"provider": provider, "id": "fixture", "species_key": species,
                         **{r + "_url": p.as_uri() for r, p in roles.items()},
                         **{r + "_filename": p.name if provider == "ensemblplants" else "" for r, p in roles.items()}})
    plan = tmp_path / "plan.json"
    planned = run_python(PLAN_SCRIPT, "--provider", "all", "--download-manifest", str(manifest),
                         "--download-dir", str(tmp_path / "downloads"), "--outfile", str(plan),
                         *(["--stage-downloads"] if staged else []))
    assert planned.returncode == 0, planned.stderr
    args = ("--task-plan", str(plan), "--task-index", "1", "--describe-only",
            "--task-meta-output", str(tmp_path / "meta.json"), "--species-cds-dir", str(tmp_path / "cds"),
            "--species-gff-dir", str(tmp_path / "gff"), "--species-genome-dir", str(tmp_path / "genome"))
    if staged:
        result = run_python(STAGE_SCRIPT, "--task-plan", str(plan), "--require-gff", "--require-genome")
        assert result.returncode == 0, result.stdout + result.stderr
    result = run_python(RUN_TASK_SCRIPT, *args)
    assert result.returncode == 0, result.stdout + result.stderr
    receipt = Path(str(plan) + ".tasks/1.json")
    task = json.loads(receipt.read_text())["task"]
    assert task["species_key"] == task["species_prefix"] == species
    assert task["provider"] == provider
    for role, source in roles.items():
        assert Path(task[role + "_path"]).read_bytes() == source.read_bytes()
    if provider == "ensemblplants":
        assert Path(task["cds_path"]).name == roles["cds"].name
        assert Path(task["cds_path"]).parent.name == "original_files"
    else:
        assert Path(task["cds_path"]).parent.name == species
    before = receipt.read_bytes()
    if staged:
        resumed = run_python(STAGE_SCRIPT, "--task-plan", str(plan), "--require-gff", "--require-genome")
        assert resumed.returncode == 0 and "no downloads needed" in resumed.stdout, resumed.stderr
    else:
        assert run_python(RUN_TASK_SCRIPT, *args).returncode == 0
    assert receipt.read_bytes() == before


@pytest.mark.parametrize("invalid", ["traversal", "symlink", "extension", "missing", "header_only_coge"])
def test_provider_manifest_roles_reject_invalid_inputs(tmp_path, monkeypatch, invalid):
    monkeypatch.syspath_prepend(str(SUPPORT_DIR))
    from format_species_provider_resolvers import provider_raw_dir
    from stage_input_generation_downloads import explicit_manifest_task

    provider = "coge" if invalid == "header_only_coge" else "ensemblplants"
    task = {"provider": provider, "species_key": "Example_species", "species_prefix": "Example_species"}
    raw = provider_raw_dir(provider, tmp_path, task["species_key"])
    raw.mkdir(parents=True)
    cds = raw / "alias.cds.fa"
    cds.write_text(">gene1\nATG\n")
    row = {"id": "123", "cds_url": "https://example.invalid/alias.cds.fa", "cds_filename": cds.name}
    if invalid == "traversal":
        row["cds_filename"] = "../alias.cds.fa"
    elif invalid == "symlink":
        link = raw / "link.fa"
        link.symlink_to(cds)
        row["cds_filename"] = link.name
    elif invalid == "extension":
        row["cds_filename"] = "alias.html"
    elif invalid == "missing":
        cds.unlink()
    else:
        gff = raw / "alias.gff3"
        gff.write_text("##gff-version 3\n")
        row.update(gff_url="https://example.invalid/alias.gff3", gff_filename=gff.name)
    with pytest.raises(ValueError):
        explicit_manifest_task(task, row, tmp_path)


def test_prepare_resources_are_separate_from_compute_array(tmp_path):
    helper = SUPPORT_DIR.parent / "gg_input_generation_array.py"
    result = run_python(helper, "--task-plan", str(tmp_path / "plan.json"), "--cpus", "4", "--memory", "32G",
                        "--prepare-cpus", "8", "--prepare-memory", "8G", "--partition", "compute",
                        "--prepare-partition", "network")
    assert result.returncode == 0, result.stderr
    prepare = next(line for line in result.stdout.splitlines() if "MODE=array_prepare " in line)
    worker = next(line for line in result.stdout.splitlines() if "MODE=array_worker " in line)
    assert "--cpus-per-task=8" in prepare and "--mem=8G" in prepare and "--partition=network" in prepare
    assert "--cpus-per-task=4" in worker and "--mem=32G" in worker and "--partition=compute" in worker


def test_required_genome_is_opt_in_for_local_array_planning(tmp_path):
    source = tmp_path / "Direct" / "species_wise_original"
    write_direct_species_fixture(source, "Arabidopsis_thaliana")
    (source / "Arabidopsis_thaliana" / "Arabidopsis_thaliana.genome.fa").unlink()
    args = ("--provider", "direct", "--input-dir", str(source), "--outfile", str(tmp_path / "plan.json"))
    default = run_python(PLAN_SCRIPT, *args)
    assert default.returncode == 0, default.stderr
    required = run_python(PLAN_SCRIPT, *args, "--require-genome")
    assert required.returncode != 0
    assert "Required genome input is missing for Arabidopsis_thaliana" in required.stderr


def test_required_gff_is_opt_in_for_local_array_planning(tmp_path):
    source = tmp_path / "Direct" / "species_wise_original"
    write_direct_species_fixture(source, "Arabidopsis_thaliana")
    (source / "Arabidopsis_thaliana" / "Arabidopsis_thaliana.gff").unlink()
    args = ("--provider", "direct", "--input-dir", str(source), "--outfile", str(tmp_path / "plan.json"))
    default = run_python(PLAN_SCRIPT, *args)
    assert default.returncode == 0, default.stderr
    required = run_python(PLAN_SCRIPT, *args, "--require-gff")
    assert required.returncode != 0
    assert "Required GFF input is missing for Arabidopsis_thaliana" in required.stderr


def test_required_genome_rejects_partial_staged_download_without_receipt(tmp_path):
    import functools
    import threading
    from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

    raw = tmp_path / "raw"
    species = "Arabidopsis_thaliana"
    write_direct_species_fixture(raw, species)
    (raw / species / f"{species}.genome.fa").unlink()
    server = ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(SimpleHTTPRequestHandler, directory=str(raw)))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        manifest = tmp_path / "manifest.tsv"
        base = f"http://127.0.0.1:{server.server_port}/{species}/{species}"
        manifest.write_text("provider\tid\tspecies_key\tcds_url\tgff_url\tgenome_url\n"
                            f"direct\tfixture\t{species}\t{base}.cds.fa\t{base}.gff\t{base}.genome.fa\n")
        plan = tmp_path / "plan.json"
        planned = run_python(PLAN_SCRIPT, "--provider", "all", "--download-manifest", str(manifest),
                             "--download-dir", str(tmp_path / "downloads"), "--stage-downloads", "--outfile", str(plan))
        assert planned.returncode == 0, planned.stderr
        required = run_python(STAGE_SCRIPT, "--task-plan", str(plan), "--require-genome")
        assert required.returncode != 0
        assert "Required genome input is missing for Arabidopsis_thaliana" in required.stderr
        assert not Path(str(plan) + ".tasks/1.json").exists()
        optional = run_python(STAGE_SCRIPT, "--task-plan", str(plan))
        assert optional.returncode == 0, optional.stderr
        assert json.loads(Path(str(plan) + ".tasks/1.json").read_text())["task"]["genome_path"] is None
        required_cached = run_python(STAGE_SCRIPT, "--task-plan", str(plan), "--require-genome")
        assert required_cached.returncode != 0
        assert "Required genome input is missing for Arabidopsis_thaliana" in required_cached.stderr
    finally:
        server.shutdown()
        server.server_close()
        thread.join(3)


def test_required_gff_rejects_partial_staged_download_without_receipt(tmp_path):
    import functools
    import threading
    from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

    raw = tmp_path / "raw"
    species = "Arabidopsis_thaliana"
    write_direct_species_fixture(raw, species)
    (raw / species / f"{species}.gff").unlink()
    server = ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(SimpleHTTPRequestHandler, directory=str(raw)))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        manifest = tmp_path / "manifest.tsv"
        base = f"http://127.0.0.1:{server.server_port}/{species}/{species}"
        manifest.write_text("provider\tid\tspecies_key\tcds_url\tgff_url\tgenome_url\n"
                            f"direct\tfixture\t{species}\t{base}.cds.fa\t{base}.gff\t{base}.genome.fa\n")
        plan = tmp_path / "plan.json"
        planned = run_python(PLAN_SCRIPT, "--provider", "all", "--download-manifest", str(manifest),
                             "--download-dir", str(tmp_path / "downloads"), "--stage-downloads", "--outfile", str(plan))
        assert planned.returncode == 0, planned.stderr
        required = run_python(STAGE_SCRIPT, "--task-plan", str(plan), "--require-gff")
        assert required.returncode != 0
        assert "Required GFF input is missing for Arabidopsis_thaliana" in required.stderr
        assert not Path(str(plan) + ".tasks/1.json").exists()
        optional = run_python(STAGE_SCRIPT, "--task-plan", str(plan))
        assert optional.returncode == 0, optional.stderr
        cached = json.loads(Path(str(plan) + ".tasks/1.json").read_text())["task"]
        assert cached["gff_path"] is None
        required_cached = run_python(STAGE_SCRIPT, "--task-plan", str(plan), "--require-gff")
        assert required_cached.returncode != 0
    finally:
        server.shutdown()
        server.server_close()
        thread.join(3)


def test_required_genome_summary_rejects_missing_or_empty_output(tmp_path):
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nATG\n")
    summary = tmp_path / "species.tsv"
    summary.write_text("species_prefix\tgenome_output_path\nArabidopsis_thaliana\t" + str(genome) + "\n")
    ok = run_python(REQUIRE_GENOMES_SCRIPT, "--species-summary", str(summary), "--expected-task-count", "1")
    assert ok.returncode == 0, ok.stderr
    genome.write_text("")
    missing = run_python(REQUIRE_GENOMES_SCRIPT, "--species-summary", str(summary), "--expected-task-count", "1")
    assert missing.returncode != 0
    assert "Required formatted genome is missing" in missing.stderr
    genome.write_text("not a FASTA\nATG\n")
    invalid = run_python(REQUIRE_GENOMES_SCRIPT, "--species-summary", str(summary))
    assert invalid.returncode != 0
    genome_gz = tmp_path / "genome.fa.gz"
    with gzip.open(genome_gz, "wt") as handle:
        handle.write(">chr1\nATG\n")
    summary.write_text("species_prefix\tgenome_output_path\nArabidopsis_thaliana\t" + str(genome_gz) + "\n")
    compressed = run_python(REQUIRE_GENOMES_SCRIPT, "--species-summary", str(summary))
    assert compressed.returncode == 0, compressed.stderr


def test_required_cds_and_gff_outputs_are_independent(tmp_path):
    cds = tmp_path / "cds.fa.gz"
    gff = tmp_path / "genes.gff.gz"
    with gzip.open(cds, "wt") as handle:
        handle.write(">gene1\nATG\n")
    with gzip.open(gff, "wt") as handle:
        handle.write("##gff-version 3\nchr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=gene1\n")
    summary = tmp_path / "species.tsv"
    summary.write_text("species_prefix\tcds_output_path\tgff_output_path\n"
                       f"Arabidopsis_thaliana\t{cds}\t{gff}\n")
    flags = ("--species-summary", str(summary), "--require-cds", "--require-gff")
    valid = run_python(REQUIRE_OUTPUTS_SCRIPT, *flags)
    assert valid.returncode == 0, valid.stderr
    with gzip.open(gff, "wt") as handle:
        handle.write("##gff-version 3\n")
    invalid_gff = run_python(REQUIRE_OUTPUTS_SCRIPT, *flags)
    assert invalid_gff.returncode != 0
    assert "Required formatted GFF is missing or invalid" in invalid_gff.stderr
    cds_only = run_python(REQUIRE_OUTPUTS_SCRIPT, "--species-summary", str(summary), "--require-cds")
    assert cds_only.returncode == 0, cds_only.stderr
    with gzip.open(cds, "wt") as handle:
        handle.write(">gene1\n")
    invalid_cds = run_python(REQUIRE_OUTPUTS_SCRIPT, "--species-summary", str(summary), "--require-cds")
    assert invalid_cds.returncode != 0
    assert "Required formatted CDS is missing or invalid" in invalid_cds.stderr


def test_required_outputs_reject_late_corruption_and_truncated_gzip(tmp_path):
    genome = tmp_path / "genome.fa.gz"
    cds = tmp_path / "cds.fa.gz"
    gff = tmp_path / "genes.gff.gz"
    summary = tmp_path / "species.tsv"
    summary.write_text(
        "species_prefix\tgenome_output_path\tcds_output_path\tgff_output_path\n"
        f"Arabidopsis_thaliana\t{genome}\t{cds}\t{gff}\n"
    )
    valid_fasta = ">chr1\nATG\n>chr2\nGCC\n"
    valid_gff = "##gff-version 3\nchr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=gene1\n"
    genome.write_bytes(gzip.compress(valid_fasta.encode()))
    cds.write_bytes(gzip.compress(valid_fasta.encode()))
    gff.write_bytes(gzip.compress(valid_gff.encode()))
    args = ("--species-summary", str(summary), "--require-genome", "--require-cds", "--require-gff")
    assert run_python(REQUIRE_OUTPUTS_SCRIPT, *args).returncode == 0

    genome.write_bytes(gzip.compress(valid_fasta.encode())[:-4])
    invalid_genome = run_python(REQUIRE_OUTPUTS_SCRIPT, *args)
    assert "Required formatted genome is missing or invalid" in invalid_genome.stderr
    genome.write_bytes(gzip.compress(valid_fasta.encode()))

    cds.write_bytes(gzip.compress(">gene1\nATG\n>gene2\n".encode()))
    invalid_cds = run_python(REQUIRE_OUTPUTS_SCRIPT, *args)
    assert "Required formatted CDS is missing or invalid" in invalid_cds.stderr
    cds.write_bytes(gzip.compress(valid_fasta.encode()))

    gff.write_bytes(gzip.compress(valid_gff.encode())[:-4])
    invalid_gff = run_python(REQUIRE_OUTPUTS_SCRIPT, *args)
    assert "Required formatted GFF is missing or invalid" in invalid_gff.stderr
    gff.write_bytes(gzip.compress((valid_gff + "chr1\tbroken\n").encode()))
    invalid_late_gff = run_python(REQUIRE_OUTPUTS_SCRIPT, *args)
    assert "Required formatted GFF is missing or invalid" in invalid_late_gff.stderr


def test_required_outputs_reject_invalid_sequence_symbols_and_gff_fields(tmp_path):
    cds = tmp_path / "cds.fa.gz"
    gff = tmp_path / "genes.gff.gz"
    summary = tmp_path / "species.tsv"
    summary.write_text(
        "species_prefix\tcds_output_path\tgff_output_path\n"
        f"Species_a\t{cds}\t{gff}\n",
        encoding="utf-8",
    )
    args = ("--species-summary", str(summary), "--require-cds", "--require-gff")
    valid_gff = "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=g1\n"
    gff.write_bytes(gzip.compress(valid_gff.encode()))
    for sequence in ("???", "ATG@", "---"):
        cds.write_bytes(gzip.compress(f">g1\n{sequence}\n".encode()))
        rejected = run_python(REQUIRE_OUTPUTS_SCRIPT, *args)
        assert rejected.returncode != 0
        assert "Required formatted CDS is missing or invalid" in rejected.stderr

    cds.write_bytes(gzip.compress(b">g1\nATGNRYS\n"))
    assert run_python(REQUIRE_OUTPUTS_SCRIPT, *args).returncode == 0
    invalid_features = (
        "chr1\tsrc\tgene\tabc\t-2\t.\t?\t9\tID=g1\n",
        "chr1\tsrc\tgene\t9\t1\t.\t+\t.\tID=g1\n",
        "chr1\tsrc\tgene\t1\t9\tNaN\t+\t.\tID=g1\n",
        "chr1\tsrc\tgene\t1\t9\t.\tinvalid\t.\tID=g1\n",
    )
    for feature in invalid_features:
        gff.write_bytes(gzip.compress(feature.encode()))
        rejected = run_python(REQUIRE_OUTPUTS_SCRIPT, *args)
        assert rejected.returncode != 0
        assert "Required formatted GFF is missing or invalid" in rejected.stderr


def run_python(script: Path, *args):
    return subprocess.run(
        [sys.executable, str(script), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def write_direct_species_fixture(root: Path, species_name: str) -> None:
    species_dir = root / species_name
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / f"{species_name}.cds.fa").write_text(
        f">{species_name}.gene1.t1\nATGAAATTT\n",
        encoding="utf-8",
    )
    (species_dir / f"{species_name}.gff").write_text(
        "\n".join(
            [
                "##gff-version 3",
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=gene1.t1;Parent=gene1",
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=gene1.t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    (species_dir / f"{species_name}.genome.fa").write_text(
        ">chr1\nATGAAATTT\n",
        encoding="utf-8",
    )


def test_plan_input_generation_tasks_discovers_direct_species(tmp_path: Path):
    input_root = tmp_path / "Direct" / "species_wise_original"
    write_direct_species_fixture(input_root, "Arabidopsis_thaliana")
    write_direct_species_fixture(input_root, "Oryza_sativa")
    task_plan = tmp_path / "task_plan.json"

    completed = run_python(
        PLAN_SCRIPT,
        "--provider",
        "direct",
        "--input-dir",
        str(input_root),
        "--outfile",
        str(task_plan),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    payload = json.loads(task_plan.read_text(encoding="utf-8"))
    assert payload["task_count"] == 2
    assert payload["species"] == ["Arabidopsis_thaliana", "Oryza_sativa"]
    assert payload["tasks"][0]["provider"] == "direct"
    assert payload["tasks"][0]["gene_grouping_mode"] == "rescue_overlap"
    assert payload["tasks"][0]["gff_repair_mode"] == "safe"
    assert payload["tasks"][0]["format_strict"] is False


def test_run_input_generation_task_and_merge_shards(tmp_path: Path):
    input_root = tmp_path / "Direct" / "species_wise_original"
    write_direct_species_fixture(input_root, "Arabidopsis_thaliana")
    write_direct_species_fixture(input_root, "Oryza_sativa")
    task_plan = tmp_path / "task_plan.json"
    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    shard_dir = tmp_path / "species_summary_shards"
    stats_dir = tmp_path / "task_stats_shards"
    meta_dir = tmp_path / "task_meta_shards"
    shard_dir.mkdir(parents=True, exist_ok=True)
    stats_dir.mkdir(parents=True, exist_ok=True)
    meta_dir.mkdir(parents=True, exist_ok=True)

    completed = run_python(
        PLAN_SCRIPT,
        "--provider",
        "direct",
        "--input-dir",
        str(input_root),
        "--outfile",
        str(task_plan),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    for task_index in (1, 2):
        completed = run_python(
            RUN_TASK_SCRIPT,
            "--task-plan",
            str(task_plan),
            "--task-index",
            str(task_index),
            "--species-cds-dir",
            str(out_cds),
            "--species-gff-dir",
            str(out_gff),
            "--species-genome-dir",
            str(out_genome),
            "--species-summary-output",
            str(shard_dir / f"{task_index}.tsv"),
            "--stats-output",
            str(stats_dir / f"{task_index}.json"),
            "--task-meta-output",
            str(meta_dir / f"{task_index}.json"),
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    (stats_dir / "1.mapping.json").write_text(json.dumps({
        "phase_conflicts_total": 3,
        "utr_conflicts_total": 1,
    }), encoding="utf-8")

    aggregate_stats = tmp_path / "aggregate_stats.json"
    mapping_qc = tmp_path / "species_mapping_qc.tsv"
    merged_species_summary = tmp_path / "gg_input_generation_species.tsv"
    completed = run_python(
        MERGE_SCRIPT,
        "--species-summary-shard-dir",
        str(shard_dir),
        "--species-summary-output",
        str(merged_species_summary),
        "--task-stats-dir",
        str(stats_dir),
        "--aggregate-stats-output",
        str(aggregate_stats),
        "--mapping-qc-output", str(mapping_qc),
        "--expected-task-count",
        "2",
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    cds_outputs = sorted(out_cds.glob("*.fa.gz"))
    assert len(cds_outputs) == 2
    with gzip.open(cds_outputs[0], "rt", encoding="utf-8") as handle:
        assert handle.read().startswith(">")

    with open(merged_species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 2
    assert {row["species_prefix"] for row in rows} == {"Arabidopsis_thaliana", "Oryza_sativa"}

    payload = json.loads(aggregate_stats.read_text(encoding="utf-8"))
    assert payload["task_stats_files"] == 2
    assert payload["mapping_qc_files"] == 1
    assert payload["phase_conflicts_total"] == 3
    assert payload["utr_conflicts_total"] == 1
    assert payload["num_species_cds_files"] == 2
    assert payload["num_species_gff_files"] == 2
    assert payload["num_species_genome_files"] == 2
    with mapping_qc.open(newline="") as handle:
        qc_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [(row["species_prefix"], row["qc_status"]) for row in qc_rows] == [
        ("Arabidopsis_thaliana", "available"), ("Oryza_sativa", "not_recorded")
    ]
    assert qc_rows[0]["phase_conflicts_total"] == "3"
    assert payload["cds_gff_records_mapped"] == 2
    assert payload["cds_gff_records_unmapped"] == 0


@pytest.mark.parametrize("qc_kind", ["mapping", "longest"])
def test_merge_requires_validation_qc_to_be_bound_to_worker_receipt(tmp_path: Path, qc_kind):
    raw = tmp_path / "Direct" / "species_wise_original"
    write_direct_species_fixture(raw, "Arabidopsis_thaliana")
    plan = tmp_path / "plan.json"
    state = SUPPORT_DIR / "input_generation_array_state.py"
    assert run_python(PLAN_SCRIPT, "--provider", "direct", "--input-dir", str(raw),
                      "--outfile", str(plan)).returncode == 0
    shard_dir = tmp_path / "species_summary_shards"
    stats_dir = tmp_path / "task_stats_shards"
    shard_dir.mkdir()
    stats_dir.mkdir()
    summary = shard_dir / "1.tsv"
    stats = stats_dir / "1.json"
    generated = run_python(
        RUN_TASK_SCRIPT, "--task-plan", str(plan), "--task-index", "1",
        "--species-cds-dir", str(tmp_path / "cds"),
        "--species-gff-dir", str(tmp_path / "gff"),
        "--species-genome-dir", str(tmp_path / "genome"),
        "--species-summary-output", str(summary), "--stats-output", str(stats),
    )
    assert generated.returncode == 0, generated.stderr
    receipt_args = ("complete", "--task-plan", str(plan), "--task-index", "1",
                    "--file", str(summary), "--file", str(stats))
    assert run_python(state, *receipt_args).returncode == 0
    qc = stats_dir / ("1." + qc_kind + ".json")
    payload = ({"phase_conflicts_total": 2, "utr_conflicts_total": 1} if qc_kind == "mapping" else
               {"source_gene_validation": {"Arabidopsis_thaliana": {
                   "source_gene_complete": False, "source_gene_check": "unresolved",
                   "source_gene_unresolved_records": 2}}})
    qc.write_text(json.dumps(payload), encoding="utf-8")
    merge_args = (
        "--species-summary-shard-dir", str(shard_dir),
        "--species-summary-output", str(tmp_path / "merged.tsv"),
        "--task-stats-dir", str(stats_dir),
        "--aggregate-stats-output", str(tmp_path / "aggregate.json"),
        "--expected-task-count", "1", "--task-plan", str(plan),
    )
    unbound = run_python(MERGE_SCRIPT, *merge_args)
    assert unbound.returncode != 0
    assert ("Mapping QC" if qc_kind == "mapping" else "Source ownership QC") + " is not bound" in unbound.stderr
    assert run_python(state, *receipt_args, "--file", str(qc)).returncode == 0
    merged = run_python(MERGE_SCRIPT, *merge_args)
    assert merged.returncode == 0, merged.stderr
    if qc_kind == "longest":
        aggregate = json.loads((tmp_path / "aggregate.json").read_text())
        assert aggregate["ownership_qc_files"] == 1
        assert aggregate["source_gene_validation"] == payload["source_gene_validation"]
    qc.write_text(json.dumps({"phase_conflicts_total": 0}), encoding="utf-8")
    assert run_python(MERGE_SCRIPT, *merge_args).returncode != 0


def test_manifest_planning_defers_downloads_and_freezes_inputs(tmp_path):
    source = tmp_path / "sources"
    write_direct_species_fixture(source, "Arabidopsis_thaliana")
    manifest = tmp_path / "manifest.tsv"
    species = "Arabidopsis_thaliana"
    raw = source / species
    fields = ["provider", "id", "species_key", "cds_url", "gff_url", "genome_url"]
    row = ["direct", "fixture", species, (raw / (species + ".cds.fa")).as_uri(),
           (raw / (species + ".gff")).as_uri(), (raw / (species + ".genome.fa")).as_uri()]
    manifest.write_text("\t".join(fields) + "\n" + "\t".join(row) + "\n")
    plan = tmp_path / "plan.json"
    downloads = tmp_path / "downloads"
    args = ["--provider", "all", "--download-manifest", str(manifest), "--download-dir", str(downloads), "--outfile", str(plan)]
    prepared = run_python(PLAN_SCRIPT, *args)
    assert prepared.returncode == 0, prepared.stderr
    assert not downloads.exists()
    frozen = plan.read_bytes()
    manifest.unlink()  # The worker must consume the frozen row, not mutable input.
    worker = run_python(RUN_TASK_SCRIPT, "--task-plan", str(plan), "--task-index", "1",
                        "--species-cds-dir", str(tmp_path / "cds"), "--species-gff-dir", str(tmp_path / "gff"),
                        "--species-genome-dir", str(tmp_path / "genome"), "--describe-only",
                        "--task-meta-output", str(tmp_path / "meta.json"))
    assert worker.returncode == 0, worker.stderr + worker.stdout
    assert downloads.exists()
    meta = json.loads((tmp_path / "meta.json").read_text())
    assert meta["species_prefix"] == species
    assert Path(meta["cds_path"]).is_file()
    assert plan.read_bytes() == frozen


def test_duplicate_species_and_changed_plan_are_rejected(tmp_path):
    raw = tmp_path / "raw"
    write_direct_species_fixture(raw, "Arabidopsis_thaliana")
    plan = tmp_path / "plan.json"
    args = ["--provider", "direct", "--input-dir", str(raw), "--outfile", str(plan)]
    assert run_python(PLAN_SCRIPT, *args).returncode == 0
    frozen = plan.read_bytes()
    assert run_python(PLAN_SCRIPT, *args).returncode == 0
    write_direct_species_fixture(raw, "Oryza_sativa")
    assert run_python(PLAN_SCRIPT, *args).returncode != 0
    assert plan.read_bytes() == frozen
    manifest = tmp_path / "duplicate.tsv"
    manifest.write_text("provider\tid\tspecies_key\n" + "direct\tx\tArabidopsis_thaliana\n" * 2)
    result = run_python(PLAN_SCRIPT, "--provider", "all", "--download-manifest", str(manifest),
                        "--download-dir", str(tmp_path / "downloads"), "--outfile", str(tmp_path / "bad.json"))
    assert result.returncode != 0
    assert "Duplicate species" in result.stderr


def test_receipts_detect_changed_outputs_and_retry_selects_only_pending(tmp_path):
    state = SUPPORT_DIR / "input_generation_array_state.py"
    raw = tmp_path / "raw"
    write_direct_species_fixture(raw, "Arabidopsis_thaliana")
    write_direct_species_fixture(raw, "Oryza_sativa")
    plan = tmp_path / "plan.json"
    assert run_python(PLAN_SCRIPT, "--provider", "direct", "--input-dir", str(raw), "--outfile", str(plan)).returncode == 0
    assert run_python(state, "configure", "--task-plan", str(plan), "--prepare").returncode == 0
    assert run_python(state, "prepared", "--task-plan", str(plan)).returncode == 0
    output = tmp_path / "result"
    output.write_text("valid")
    assert run_python(state, "complete", "--task-plan", str(plan), "--task-index", "1", "--file", str(output)).returncode == 0
    pending = run_python(state, "pending", "--task-plan", str(plan))
    assert pending.stdout.strip() == "2"
    helper = SUPPORT_DIR.parent / "gg_input_generation_array.py"
    preview = run_python(helper, "--task-plan", str(plan), "--retry", "--cpus", "3", "--max-running", "2")
    assert preview.returncode == 0, preview.stderr
    assert "--array=2%2" in preview.stdout
    assert "--dependency=afterok:WORKER_JOB_ID" in preview.stdout
    assert "--cpus-per-task=3" in preview.stdout
    output.write_text("corrupt")
    assert run_python(state, "pending", "--task-plan", str(plan)).stdout.strip() == "1,2"


def test_retry_submission_refuses_any_active_legacy_worker_array(tmp_path, monkeypatch):
    import os
    state = SUPPORT_DIR / "input_generation_array_state.py"
    raw = tmp_path / "raw"
    write_direct_species_fixture(raw, "Arabidopsis_thaliana")
    plan = tmp_path / "plan.json"
    assert run_python(PLAN_SCRIPT, "--provider", "direct", "--input-dir", str(raw), "--outfile", str(plan)).returncode == 0
    assert run_python(state, "configure", "--task-plan", str(plan), "--prepare").returncode == 0
    assert run_python(state, "prepared", "--task-plan", str(plan)).returncode == 0
    fake = tmp_path / "squeue"
    fake.write_text("#!/bin/sh\necho 40457_29\n")
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path) + os.pathsep + os.environ["PATH"])
    helper = SUPPORT_DIR.parent / "gg_input_generation_array.py"
    result = run_python(helper, "--task-plan", str(plan), "--retry", "--submit")
    assert result.returncode != 0
    assert "40457_29" in result.stderr


def test_slurm_helper_waits_for_prepare_and_submits_afterok(tmp_path, monkeypatch):
    import os
    plan = tmp_path / "plan.json"
    calls = tmp_path / "calls.jsonl"
    fake = tmp_path / "sbatch"
    fake.write_text("#!" + sys.executable + "\n" + '''import hashlib, json, os, sys
from pathlib import Path
log = Path(os.environ["TEST_CALLS"])
mode = os.environ["GG_INPUT_INPUT_GENERATION_MODE"]
with log.open("a") as handle:
    handle.write(json.dumps({"mode": mode, "argv": sys.argv[1:]}) + "\\n")
if mode == "array_prepare":
    if os.environ.get("TEST_PREPARE_FAIL"):
        print("simulated Slurm rejection", file=sys.stderr)
        sys.exit(1)
    Path(os.environ["GG_INPUT_TASK_PLAN_OUTPUT"]).write_text(json.dumps({"task_count": 2, "tasks": [{"species_prefix": "A_b"}, {"species_prefix": "C_d"}]}))
    plan = Path(os.environ["GG_INPUT_TASK_PLAN_OUTPUT"])
    Path(str(plan) + ".settings.json").write_text("{}")
    Path(str(plan) + ".prepared.json").write_text(json.dumps({"plan_sha256": hashlib.sha256(plan.read_bytes()).hexdigest(), "settings_sha256": hashlib.sha256(b"{}").hexdigest()}))
print({"array_prepare": "101", "array_worker": "102", "array_finalize": "103"}[mode])
''')
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path) + os.pathsep + os.environ["PATH"])
    monkeypatch.setenv("TEST_CALLS", str(calls))
    helper = SUPPORT_DIR.parent / "gg_input_generation_array.py"
    result = run_python(helper, "--task-plan", str(plan), "--submit", "--max-running", "5")
    assert result.returncode == 0, result.stderr
    submissions = [json.loads(line) for line in calls.read_text().splitlines()]
    assert [call["mode"] for call in submissions] == ["array_prepare", "array_worker", "array_finalize"]
    assert "--wait" in submissions[0]["argv"]
    assert "--array=1-2%5" in submissions[1]["argv"]
    assert "--dependency=afterok:102" in submissions[2]["argv"]
    assert not any("--mem-per-cpu" in arg for call in submissions for arg in call["argv"])
    assert all(any(arg.startswith("--wrap=") for arg in call["argv"]) for call in submissions)
    calls.unlink()
    monkeypatch.setenv("TEST_PREPARE_FAIL", "1")
    result = run_python(helper, "--task-plan", str(plan), "--submit")
    assert result.returncode != 0
    assert len(calls.read_text().splitlines()) == 1
    assert "simulated Slurm rejection" in result.stderr


@pytest.mark.parametrize("rescue_resources", [False, True])
def test_slurm_rescue_chain_uses_frozen_comparison_and_species_arrays(tmp_path, monkeypatch, rescue_resources):
    import os
    plan = tmp_path / "tmp" / "plan.json"
    plan.parent.mkdir()
    calls = tmp_path / "calls.jsonl"
    fake = tmp_path / "sbatch"
    fake.write_text("#!" + sys.executable + "\n" + '''import hashlib, json, os
from pathlib import Path
mode = os.environ["GG_INPUT_INPUT_GENERATION_MODE"]
with Path(os.environ["TEST_CALLS"]).open("a") as handle:
    handle.write(json.dumps({"mode": mode, "rescue": os.environ.get("GG_INPUT_RUN_GENE_MODEL_RESCUE"), "argv": __import__("sys").argv[1:]}) + "\\n")
plan = Path(os.environ["GG_INPUT_TASK_PLAN_OUTPUT"])
if mode == "array_prepare":
    plan.write_text(json.dumps({"task_count": 2, "tasks": [{"species_prefix": "A_b"}, {"species_prefix": "C_d"}]}))
    Path(str(plan) + ".settings.json").write_text("{}")
    Path(str(plan) + ".prepared.json").write_text(json.dumps({"plan_sha256": hashlib.sha256(plan.read_bytes()).hexdigest(), "settings_sha256": hashlib.sha256(b"{}").hexdigest()}))
if mode == "array_finalize":
    root = Path(os.environ["GG_INPUT_GENE_MODEL_RESCUE_DIR"])
    root.mkdir(exist_ok=True)
    (root / "plan.json").write_text(json.dumps({"species": ["A_b", "C_d"], "synteny_jobs": [{"index": i, "id": str(i)} for i in range(1, 4)]}))
print({"array_prepare": 101, "array_worker": 102, "array_finalize": 103, "rescue_synteny": 104, "rescue_models": 105, "rescue_finalize": 106}[mode])
''')
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path) + os.pathsep + os.environ["PATH"])
    monkeypatch.setenv("TEST_CALLS", str(calls))
    helper = SUPPORT_DIR.parent / "gg_input_generation_array.py"
    extra = ["--rescue-max-running", "6", "--rescue-cpus", "8", "--rescue-memory", "64G"] if rescue_resources else []
    result = run_python(helper, "--task-plan", str(plan), "--rescue", "--submit", "--max-running", "5", *extra)
    assert result.returncode == 0, result.stderr
    submissions = [json.loads(line) for line in calls.read_text().splitlines()]
    assert [c["mode"] for c in submissions] == ["array_prepare", "array_worker", "array_finalize", "rescue_synteny", "rescue_models", "rescue_finalize"]
    assert all(c["rescue"] == "1" for c in submissions)
    assert "--wait" in submissions[2]["argv"]
    assert "--array=1-3%5" in submissions[3]["argv"]
    assert "--dependency=afterok:104" in submissions[4]["argv"]
    assert "--array=1-2%" + ("6" if rescue_resources else "5") in submissions[4]["argv"]
    for index in [4, 5]:
        assert "--cpus-per-task=" + ("8" if rescue_resources else "4") in submissions[index]["argv"]
        assert "--mem=" + ("64G" if rescue_resources else "32G") in submissions[index]["argv"]
    assert "--mem=32G" in submissions[3]["argv"]
    assert "--dependency=afterok:105" in submissions[5]["argv"]


def test_slurm_rescue_retry_requires_qc_worker_receipt(tmp_path):
    import hashlib
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({"task_count": 1, "tasks": [{"species_prefix": "A_b"}]}))
    Path(str(plan) + ".settings.json").write_text("{}")
    Path(str(plan) + ".prepared.json").write_text(json.dumps({
        "plan_sha256": hashlib.sha256(plan.read_bytes()).hexdigest(),
        "settings_sha256": hashlib.sha256(b"{}").hexdigest()}))
    rescue = tmp_path / "rescue"
    rescue.mkdir()
    (rescue / "plan.json").write_text(json.dumps({"species": ["A_b"], "synteny_jobs": []}))
    models = rescue / "rescued" / "A_b"
    models.mkdir(parents=True)
    (models / "models.json").write_text("[]")
    (models / "receipt.json").write_text(json.dumps({"key": {"plan": hashlib.sha256((rescue / "plan.json").read_bytes()).hexdigest()},
                                                    "files": {"models.json": hashlib.sha256(b"[]").hexdigest()}}))
    helper = SUPPORT_DIR.parent / "gg_input_generation_array.py"
    result = run_python(helper, "--task-plan", str(plan), "--rescue", "--rescue-output", str(rescue), "--retry")
    assert result.returncode == 0, result.stderr
    assert "GG_INPUT_INPUT_GENERATION_MODE=rescue_models" in result.stdout and "--array=1%8" in result.stdout


def test_slurm_rescue_active_retry_rejects_before_submitting_any_job(tmp_path, monkeypatch):
    import os
    for name, contents in (("squeue", "#!/bin/sh\nprintf '123_1\\n'\n"),
                           ("sbatch", "#!/bin/sh\ntouch \"$TEST_SUBMISSION\"\n")):
        path = tmp_path / name
        path.write_text(contents)
        path.chmod(0o755)
    submitted = tmp_path / "submitted"
    monkeypatch.setenv("PATH", str(tmp_path) + os.pathsep + os.environ["PATH"])
    monkeypatch.setenv("TEST_SUBMISSION", str(submitted))
    helper = SUPPORT_DIR.parent / "gg_input_generation_array.py"
    result = run_python(helper, "--task-plan", str(tmp_path / "plan.json"), "--rescue", "--retry", "--submit")
    assert result.returncode != 0 and "rescue jobs are active" in result.stderr
    assert not submitted.exists()


@pytest.mark.parametrize("corruption", ["foreign_job", "foreign_species", "malformed_files", "prepared",
                                       "repaired_comparison", "repaired_prepared", "effective", "none"])
def test_slurm_rescue_retry_checks_full_job_identity_and_artifacts(tmp_path, corruption):
    import hashlib
    plan = tmp_path / "task_plan.json"
    plan.write_text(json.dumps({"task_count": 2, "tasks": [{"species_prefix": "A_b"}, {"species_prefix": "C_d"}]}))
    Path(str(plan) + ".settings.json").write_text("{}")
    Path(str(plan) + ".prepared.json").write_text(json.dumps({
        "plan_sha256": hashlib.sha256(plan.read_bytes()).hexdigest(),
        "settings_sha256": hashlib.sha256(b"{}").hexdigest()}))
    root = tmp_path / "rescue"
    root.mkdir()
    job = {"id": "comparison_000001", "index": 1, "a": "A_b", "b": "C_d", "kind": "pair"}
    (root / "plan.json").write_text(json.dumps({"species": ["A_b", "C_d"], "synteny_jobs": [job],
                                               "donors": {"A_b": ["C_d"], "C_d": ["A_b"]}}))
    fingerprint = hashlib.sha256((root / "plan.json").read_bytes()).hexdigest()
    def complete(directory, key):
        directory.mkdir(parents=True)
        (directory / "data").write_text("ok")
        receipt = {"key": key, "files": {"data": hashlib.sha256(b"ok").hexdigest()}}
        (directory / "receipt.json").write_text(json.dumps(receipt))
        return directory / "receipt.json"
    for n in ["A_b", "C_d"]:
        complete(root / "prepared" / n, {"plan": fingerprint, "species": n})
    comparison = complete(root / "synteny" / job["id"], {"plan": fingerprint, "job": job,
                           "prepared": {n: hashlib.sha256((root / "prepared" / n / "receipt.json").read_bytes()).hexdigest()
                                        for n in ["A_b", "C_d"]}})
    for n in ["A_b", "C_d"]:
        models = complete(root / "rescued" / n, {"plan": fingerprint, "species": n,
                          "comparisons": {job["id"]: hashlib.sha256(comparison.read_bytes()).hexdigest()},
                          "prepared": {p: hashlib.sha256((root / "prepared" / p / "receipt.json").read_bytes()).hexdigest()
                                       for p in ["A_b", "C_d"]}})
        effective = complete(root / "effective" / n, {"plan": fingerprint,
                             "rescue_receipts": {n: hashlib.sha256(models.read_bytes()).hexdigest()}})
        complete(root / "workers" / n, {"plan": fingerprint, "species": n})
        worker = root / "workers" / n / "receipt.json"
        payload = json.loads(worker.read_text())
        payload["files"].update({f"../../rescued/{n}/receipt.json": hashlib.sha256(models.read_bytes()).hexdigest(),
                                 f"../../effective/{n}/receipt.json": hashlib.sha256(effective.read_bytes()).hexdigest()})
        worker.write_text(json.dumps(payload))
    if corruption == "foreign_job":
        payload = json.loads(comparison.read_text())
        payload["key"]["job"]["id"] = "comparison_000002"
        comparison.write_text(json.dumps(payload))
    elif corruption == "foreign_species":
        (root / "workers" / "C_d" / "receipt.json").write_bytes((root / "workers" / "A_b" / "receipt.json").read_bytes())
    elif corruption == "malformed_files":
        comparison.write_text(json.dumps({"key": {"plan": fingerprint, "job": job}, "files": ["invalid"]}))
    elif corruption == "prepared":
        (root / "prepared" / "A_b" / "data").write_text("corrupted")
    elif corruption != "none":
        directory = {"repaired_comparison": comparison.parent,
                     "repaired_prepared": root / "prepared" / "A_b",
                     "effective": root / "effective" / "A_b"}[corruption]
        (directory / "data").write_text("valid new result")
        payload = json.loads((directory / "receipt.json").read_text())
        payload["files"]["data"] = hashlib.sha256(b"valid new result").hexdigest()
        (directory / "receipt.json").write_text(json.dumps(payload))
    helper = SUPPORT_DIR.parent / "gg_input_generation_array.py"
    result = run_python(helper, "--task-plan", str(plan), "--rescue", "--rescue-output", str(root), "--retry")
    assert result.returncode == 0, result.stderr
    lines = result.stdout.splitlines()
    pair_lines = [line for line in lines if "GG_INPUT_INPUT_GENERATION_MODE=rescue_synteny" in line]
    worker_lines = [line for line in lines if "GG_INPUT_INPUT_GENERATION_MODE=rescue_models" in line]
    if corruption == "foreign_species":
        assert pair_lines == [] and len(worker_lines) == 1 and "--array=2%8" in worker_lines[0]
    elif corruption == "effective":
        assert pair_lines == [] and len(worker_lines) == 1 and "--array=1%8" in worker_lines[0]
    elif corruption == "repaired_comparison":
        assert pair_lines == [] and len(worker_lines) == 1 and "--array=1-2%8" in worker_lines[0]
    elif corruption == "none":
        assert pair_lines == [] and worker_lines == []
    else:
        assert len(pair_lines) == 1 and "--array=1%8" in pair_lines[0]
        assert len(worker_lines) == 1 and "--array=1-2%8" in worker_lines[0]


def test_workspace_rejects_second_plan_and_malformed_receipt_is_pending(tmp_path):
    state = SUPPORT_DIR / "input_generation_array_state.py"
    raw = tmp_path / "raw"
    write_direct_species_fixture(raw, "Arabidopsis_thaliana")
    first = tmp_path / "first.json"
    second = tmp_path / "second.json"
    for plan in (first, second):
        assert run_python(PLAN_SCRIPT, "--provider", "direct", "--input-dir", str(raw), "--outfile", str(plan)).returncode == 0
    workspace = tmp_path / "workspace"
    claim = ["claim-workspace", "--workspace", str(workspace), "--prepare", "--task-plan"]
    assert run_python(state, *claim, str(first)).returncode == 0
    assert run_python(state, *claim, str(second)).returncode != 0
    assert run_python(state, "index", "--task-plan", str(first), "--task-index", "01").stdout.strip() == "1"
    output = tmp_path / "output"
    output.write_text("ok")
    assert run_python(state, "complete", "--task-plan", str(first), "--task-index", "1", "--file", str(output)).returncode == 0
    receipt = Path(str(first) + ".completed") / "1.json"
    content = json.loads(receipt.read_text())
    content["files"] = ["bad structure"]
    receipt.write_text(json.dumps(content))
    pending = run_python(state, "pending", "--task-plan", str(first))
    assert pending.returncode == 0, pending.stderr
    assert pending.stdout.strip() == "1"


def test_manifest_source_change_and_foreign_resolved_cache_are_rejected(tmp_path):
    species = "Arabidopsis_thaliana"
    raw_root = tmp_path / "raw"
    write_direct_species_fixture(raw_root, species)
    raw = raw_root / species
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("provider\tid\tspecies_key\nlocal\t" + str(raw) + "\t" + species + "\n")
    plan = tmp_path / "plan.json"
    downloads = tmp_path / "downloads"
    assert run_python(PLAN_SCRIPT, "--provider", "local", "--download-manifest", str(manifest),
                      "--download-dir", str(downloads), "--outfile", str(plan)).returncode == 0
    worker_args = ["--task-plan", str(plan), "--task-index", "1", "--species-cds-dir", str(tmp_path / "cds"),
                   "--species-gff-dir", str(tmp_path / "gff"), "--species-genome-dir", str(tmp_path / "genome"), "--describe-only"]
    source = raw / (species + ".cds.fa")
    original = source.read_bytes()
    source.write_bytes(original + b"AAA\n")
    changed = run_python(RUN_TASK_SCRIPT, *worker_args)
    assert changed.returncode != 0 and "changed after planning" in changed.stderr
    assert not downloads.exists()
    source.write_bytes(original)
    assert run_python(RUN_TASK_SCRIPT, *worker_args).returncode == 0
    cache = Path(str(plan) + ".tasks") / "1.json"
    content = json.loads(cache.read_text())
    content["plan_sha256"] = "belongs to another plan"
    cache.write_text(json.dumps(content))
    foreign = run_python(RUN_TASK_SCRIPT, *worker_args)
    assert foreign.returncode != 0 and "another plan/task" in foreign.stderr
    assert not (tmp_path / "cds").exists()


def test_completion_rejects_raw_changes_during_worker(tmp_path):
    state = SUPPORT_DIR / "input_generation_array_state.py"
    raw = tmp_path / "raw"
    write_direct_species_fixture(raw, "Arabidopsis_thaliana")
    plan = tmp_path / "plan.json"
    assert run_python(PLAN_SCRIPT, "--provider", "direct", "--input-dir", str(raw), "--outfile", str(plan)).returncode == 0
    source = next(raw.glob("*/*.cds.fa"))
    source.write_text(source.read_text() + "AAA\n")
    output = tmp_path / "output"
    output.write_text("would be stale")
    completed = run_python(state, "complete", "--task-plan", str(plan), "--task-index", "1", "--file", str(output))
    assert completed.returncode != 0
    assert not (Path(str(plan) + ".completed") / "1.json").exists()


def test_array_worker_reports_empty_cds_before_formatting_genome(tmp_path):
    raw = tmp_path / "raw"
    write_direct_species_fixture(raw, "Arabidopsis_thaliana")
    next(raw.glob("*/*.cds.fa")).unlink()
    next(raw.glob("*/*.gff*")).write_text(
        "##gff-version 3\nchr1\t.\tregion\t1\t9\t.\t+\t.\tID=chr1\n")
    plan = tmp_path / "plan.json"
    planned = run_python(PLAN_SCRIPT, "--provider", "direct", "--input-dir", str(raw),
                         "--outfile", str(plan))
    assert planned.returncode == 0, planned.stderr
    genome_dir = tmp_path / "formatted_genome"
    result = run_python(RUN_TASK_SCRIPT, "--task-plan", str(plan), "--task-index", "1",
                        "--species-cds-dir", str(tmp_path / "cds"),
                        "--species-gff-dir", str(tmp_path / "gff"),
                        "--species-genome-dir", str(genome_dir))
    assert result.returncode != 0
    assert "No nuclear CDS records for Arabidopsis_thaliana" in result.stderr
    assert not list(genome_dir.glob("*.fa.gz"))


def test_completion_hashes_raw_receipt_alias_once(tmp_path, monkeypatch):
    import input_generation_array_state as state

    raw = tmp_path / "raw.fa"
    raw.write_text(">Species_a_gene1\nATG\n")
    alias = tmp_path / "alias.fa"
    alias.symlink_to(raw)
    output = tmp_path / "output"
    output.write_text("valid output\n")
    plan = tmp_path / "plan.json"
    expected = state.digest(raw)
    plan.write_text(json.dumps({"task_count": 1, "tasks": [{
        "species_prefix": "Species_a", "input_sha256": {str(alias): expected}}]}))
    calls = []
    original = state.digest

    def counted(path):
        calls.append(str(path))
        return original(path)

    monkeypatch.setattr(state, "digest", counted)
    monkeypatch.setattr(sys, "argv", ["state", "complete", "--task-plan", str(plan),
                        "--task-index", "1", "--file", str(raw), "--file", str(output)])
    state.main()
    assert calls.count(str(raw.resolve())) == 1
    receipt = json.loads(state.receipt_path(plan, 1).read_text())
    assert receipt["files"] == {str(raw.resolve()): expected, str(alias): expected,
                               str(output.resolve()): original(output)}
    assert state.verify_receipt(plan, 1)


def test_completion_rejects_source_mutation_while_hashing_later_output(tmp_path, monkeypatch):
    import input_generation_array_state as state

    raw = tmp_path / "raw"
    raw.write_bytes(b"original")
    output = tmp_path / "output"
    output.write_bytes(b"result")
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({"task_count": 1, "tasks": [{
        "species_prefix": "Species_a", "input_sha256": {str(raw): state.digest(raw)}}]}))
    original = state.digest

    def mutate(path):
        value = original(path)
        if Path(path) == output:
            # Cross Lustre's ctime/mtime boundary before a same-size rewrite.
            time.sleep(1.1)
            raw.write_bytes(b"modified")
        return value

    monkeypatch.setattr(state, "digest", mutate)
    monkeypatch.setattr(sys, "argv", ["state", "complete", "--task-plan", str(plan),
                        "--task-index", "1", "--file", str(raw), "--file", str(output)])
    with pytest.raises(OSError, match="changed while hashing"):
        state.main()
    assert not state.receipt_path(plan, 1).exists()


def test_digest_paths_rejects_alias_target_replacement(tmp_path, monkeypatch):
    import input_generation_array_state as state

    first, second, alias = (tmp_path / name for name in ("first", "second", "alias"))
    first.write_bytes(b"identical")
    second.write_bytes(b"identical")
    alias.symlink_to(first)
    original = state.digest

    def replace(path):
        result = original(path)
        alias.unlink()
        alias.symlink_to(second)
        return result

    monkeypatch.setattr(state, "digest", replace)
    with pytest.raises(OSError, match="changed while hashing"):
        state.digest_paths([alias])


def test_custom_output_directory_cannot_be_claimed_by_two_workspaces(tmp_path):
    state = SUPPORT_DIR / "input_generation_array_state.py"
    raw = tmp_path / "raw"
    write_direct_species_fixture(raw, "Arabidopsis_thaliana")
    first = tmp_path / "first.json"
    second = tmp_path / "second.json"
    shared = tmp_path / "custom_species_cds"
    for plan in (first, second):
        assert run_python(PLAN_SCRIPT, "--provider", "direct", "--input-dir", str(raw), "--outfile", str(plan)).returncode == 0
    assert run_python(state, "claim-workspace", "--prepare", "--task-plan", str(first), "--workspace", str(tmp_path / "ws1"), "--file", str(shared)).returncode == 0
    second_claim = run_python(state, "claim-workspace", "--prepare", "--task-plan", str(second), "--workspace", str(tmp_path / "ws2"), "--file", str(shared))
    assert second_claim.returncode != 0 and "another array plan" in second_claim.stderr
    assert not (tmp_path / "ws2" / ".array-plan.json").exists()


def test_prepared_marker_rejects_changed_shared_lineage(tmp_path):
    state = SUPPORT_DIR / "input_generation_array_state.py"
    raw = tmp_path / "raw"
    write_direct_species_fixture(raw, "Arabidopsis_thaliana")
    plan = tmp_path / "plan.json"
    assert run_python(PLAN_SCRIPT, "--provider", "direct", "--input-dir", str(raw), "--outfile", str(plan)).returncode == 0
    assert run_python(state, "configure", "--prepare", "--task-plan", str(plan)).returncode == 0
    lineage = tmp_path / "lineage.txt"
    lineage.write_text("eukaryota_odb12")
    assert run_python(state, "prepared", "--task-plan", str(plan), "--file", str(lineage)).returncode == 0
    assert run_python(state, "check-prepared", "--task-plan", str(plan)).returncode == 0
    lineage.write_text("different_lineage")
    assert run_python(state, "check-prepared", "--task-plan", str(plan)).returncode != 0


def test_nonregular_receipt_input_fails_without_blocking(tmp_path):
    import os
    raw = tmp_path / "raw"
    write_direct_species_fixture(raw, "Arabidopsis_thaliana")
    plan = tmp_path / "plan.json"
    assert run_python(PLAN_SCRIPT, "--provider", "direct", "--input-dir", str(raw), "--outfile", str(plan)).returncode == 0
    fifo = tmp_path / "pipe"
    os.mkfifo(fifo)
    result = subprocess.run([sys.executable, str(SUPPORT_DIR / "input_generation_array_state.py"), "complete",
                             "--task-plan", str(plan), "--task-index", "1", "--file", str(fifo)],
                            capture_output=True, text=True, timeout=3)
    assert result.returncode != 0


def test_array_plan_rejects_escaping_species_and_download_filenames(tmp_path):
    for index, (species, filename) in enumerate((("../Arabidopsis_thaliana", "cds.fa"),
                                               ("Arabidopsis_thaliana", "../../outside.fa"),
                                               (".hidden_species", "cds.fa"))):
        manifest = tmp_path / f"bad{index}.tsv"
        manifest.write_text(f"provider\tid\tspecies_key\tcds_filename\ndirect\tx\t{species}\t{filename}\n")
        output = tmp_path / f"bad{index}.json"
        result = run_python(PLAN_SCRIPT, "--provider", "all", "--download-manifest", str(manifest),
                            "--download-dir", str(tmp_path / "downloads"), "--outfile", str(output))
        assert result.returncode != 0 and "filename components" in result.stderr
        assert not output.exists()


@pytest.mark.parametrize('combined', [False, True])
def test_slurm_refinement_chain_routes_sparse_pairs_and_species_with_dependencies(tmp_path, monkeypatch, combined):
    import os
    plan = tmp_path / 'tmp' / 'plan.json'
    plan.parent.mkdir()
    calls = tmp_path / 'calls.jsonl'
    fake = tmp_path / 'sbatch'
    fake.write_text('#!' + sys.executable + '\n' + '''import hashlib, json, os, sys
from pathlib import Path
mode = os.environ['GG_INPUT_INPUT_GENERATION_MODE']
with Path(os.environ['TEST_CALLS']).open('a') as handle:
    handle.write(json.dumps({'mode': mode, 'argv': sys.argv[1:], 'refinement': os.environ.get('GG_INPUT_RUN_GENE_MODEL_REFINEMENT'), 'anchors': os.environ.get('GG_INPUT_GENE_MODEL_RESCUE_DIR')}) + '\\n')
plan = Path(os.environ['GG_INPUT_TASK_PLAN_OUTPUT'])
if mode == 'array_prepare':
    plan.write_text(json.dumps({'task_count': 2, 'tasks': [{'species_prefix': 'A_b'}, {'species_prefix': 'C_d'}]}))
    Path(str(plan) + '.settings.json').write_text('{}')
    Path(str(plan) + '.prepared.json').write_text(json.dumps({'plan_sha256': hashlib.sha256(plan.read_bytes()).hexdigest(), 'settings_sha256': hashlib.sha256(b'{}').hexdigest()}))
if mode == 'array_finalize':
    if os.environ.get('GG_INPUT_RUN_GENE_MODEL_RESCUE') == '1':
        anchors = Path(os.environ['GG_INPUT_GENE_MODEL_RESCUE_DIR'])
        anchors.mkdir(exist_ok=True)
        (anchors / 'plan.json').write_text(json.dumps({'species': ['A_b', 'C_d'], 'donors': {'A_b': ['C_d'], 'C_d': ['A_b']}, 'synteny_jobs': [{'a': 'A_b', 'b': 'C_d', 'index': 1, 'id': 'comparison_1'}]}))
if mode == 'refinement_prepare' or mode == 'array_finalize' and os.environ.get('GG_INPUT_RUN_GENE_MODEL_RESCUE') != '1':
    root = Path(os.environ['GG_INPUT_GENE_MODEL_REFINEMENT_DIR'])
    root.mkdir(exist_ok=True)
    anchors = Path(os.environ['GG_INPUT_GENE_MODEL_RESCUE_DIR']) if os.environ.get('GG_INPUT_RUN_GENE_MODEL_RESCUE') == '1' else root.parent / 'anchors'
    anchors.mkdir(exist_ok=True)
    if not (anchors / 'plan.json').exists():
        (anchors / 'plan.json').write_text(json.dumps({'synteny_jobs': [{'a': 'A_b', 'b': 'C_d', 'index': 1}, {'a': 'A_b', 'b': 'A_b', 'index': 2}]}))
    (root / 'plan.json').write_text(json.dumps({'species': ['A_b', 'C_d'], 'request': {'rescue_output': str(anchors), 'edges': None}}))
print({'array_prepare': 101, 'array_worker': 102, 'array_finalize': 103, 'rescue_synteny': 104, 'refinement_catalog': 105, 'refinement_correspondence': 106, 'refinement_predict': 107, 'refinement_finalize': 108, 'rescue_models': 109, 'rescue_finalize': 110, 'refinement_prepare': 111}[mode])
''')
    fake.chmod(0o755)
    monkeypatch.setenv('PATH', str(tmp_path) + os.pathsep + os.environ['PATH'])
    monkeypatch.setenv('TEST_CALLS', str(calls))
    helper = SUPPORT_DIR.parent / 'gg_input_generation_array.py'
    result = run_python(helper, '--task-plan', str(plan), '--refinement', '--submit', '--max-running', '5',
                        '--rescue-max-running', '3', '--rescue-cpus', '8', '--rescue-memory', '64G', *(['--rescue'] if combined else []))
    assert result.returncode == 0, result.stderr
    submissions = [json.loads(line) for line in calls.read_text().splitlines()]
    if combined:
        assert [row['mode'] for row in submissions[3:7]] == ['rescue_synteny', 'rescue_models', 'rescue_finalize', 'refinement_prepare']
        assert '--wait' in submissions[5]['argv']
        assert '--dependency=afterok:110' in submissions[6]['argv'] and '--wait' in submissions[6]['argv']
        submissions = submissions[:3] + submissions[7:]
    assert [row['mode'] for row in submissions] == ['array_prepare', 'array_worker', 'array_finalize',
           'rescue_synteny', 'refinement_catalog', 'refinement_correspondence', 'refinement_predict', 'refinement_finalize']
    assert all(row['refinement'] == '1' for row in submissions)
    assert '--wait' in submissions[2]['argv']
    assert '--array=1%5' in submissions[3]['argv']  # Self comparisons are not isoform correspondence.
    assert submissions[3]['anchors'].endswith('/gene_model_rescue' if combined else '/anchors')
    assert '--dependency=afterok:104:105' in submissions[5]['argv']
    assert '--dependency=afterok:106' in submissions[6]['argv']
    assert '--array=1-2%3' in submissions[6]['argv']
    assert '--cpus-per-task=8' in submissions[6]['argv'] and '--mem=64G' in submissions[6]['argv']
    assert '--dependency=afterok:107' in submissions[7]['argv']


def test_slurm_refinement_retry_rejects_active_workers_before_submission(tmp_path, monkeypatch):
    import os
    for name, text in [('squeue', '#!/bin/sh\necho 900_1\n'), ('sbatch', '#!/bin/sh\ntouch "$TEST_SUBMISSION"\n')]:
        path = tmp_path / name
        path.write_text(text)
        path.chmod(0o755)
    marker = tmp_path / 'submitted'
    monkeypatch.setenv('PATH', str(tmp_path) + os.pathsep + os.environ['PATH'])
    monkeypatch.setenv('TEST_SUBMISSION', str(marker))
    helper = SUPPORT_DIR.parent / 'gg_input_generation_array.py'
    result = run_python(helper, '--task-plan', str(tmp_path / 'plan.json'), '--refinement', '--retry', '--submit')
    assert result.returncode != 0 and 'refinement workers are active' in result.stderr
    assert not marker.exists()

"""Behavior tests for external, disposable workflow scratch storage."""
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / 'support'
RUNNER = SUPPORT / 'task_tmp.py'


def run_env(tmp_path, **overrides):
    root = tmp_path / 'scratch'
    root.mkdir(exist_ok=True)
    return dict(os.environ, GG_TMP_MOUNT=str(root), GG_TMP_WORKSPACE_ID='/project/a',
                GG_ARRAY_TASK_ID='1', **overrides)


def run(env, script, workflow='gg_gene_evolution'):
    return subprocess.run([sys.executable, str(RUNNER), '--workflow', workflow,
                           '--', 'bash', '-c', script], env=env, text=True,
                          capture_output=True)


def records(tmp_path):
    return list((tmp_path / 'scratch').glob('genegalleon-*/*/*/run-*/record.json'))


def test_success_cleanup_and_failure_retention(tmp_path):
    env = run_env(tmp_path)
    result = run(env, 'echo data > "$GG_TMP_TASK_ROOT/result"; exit 2')
    assert result.returncode == 2, result.stderr
    record = records(tmp_path)[0]
    assert (record.parent / 'work/result').read_text().strip() == 'data'
    assert json.loads(record.read_text())['exit_code'] == 2
    result = run(env, 'test -d "$TMPDIR"; test ! -e "$GG_TMP_TASK_ROOT/result"')
    assert result.returncode == 0, result.stderr
    assert not records(tmp_path)


def test_reuse_and_keep(tmp_path):
    env = run_env(tmp_path, delete_preexisting_tmp_dir='0', delete_tmp_dir='0')
    assert run(env, 'echo resume > "$GG_TMP_TASK_ROOT/result"; exit 1').returncode == 1
    result = run(env, 'test "$(cat "$GG_TMP_TASK_ROOT/result")" = resume')
    assert result.returncode == 0, result.stderr
    assert len(records(tmp_path)) == 1


def test_retention_is_scoped_and_preserves_unmanaged_files(tmp_path):
    env = run_env(tmp_path, gene_family_tmp_max_dirs='1')
    unrelated = tmp_path / 'scratch/unrelated'
    unrelated.write_text('keep')
    for task in ('1', '2', '3'):
        env['GG_ARRAY_TASK_ID'] = task
        assert run(env, 'exit 1').returncode == 1
    assert len(records(tmp_path)) == 1
    assert unrelated.read_text() == 'keep'
    env['GG_TMP_WORKSPACE_ID'] = '/project/b'
    assert run(env, 'exit 1').returncode == 1
    assert len(records(tmp_path)) == 2


def test_active_task_not_reused_or_deleted(tmp_path):
    env = run_env(tmp_path, delete_preexisting_tmp_dir='0')
    ready = tmp_path / 'ready'
    child = subprocess.Popen([sys.executable, str(RUNNER), '--workflow', 'gg_gene_evolution',
                              '--', 'bash', '-c', f'touch "{ready}"; sleep 30'], env=env,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        for _ in range(200):
            if ready.exists():
                break
            time.sleep(.02)
        assert ready.exists()
        assert run(env, 'exit 1').returncode == 1
        assert len(records(tmp_path)) == 2
        child.send_signal(signal.SIGTERM)
        child.communicate(timeout=10)
        assert child.returncode == 143
        assert len(records(tmp_path)) == 2
    finally:
        if child.poll() is None:
            child.kill()
            child.communicate()


def test_symlink_user_root_rejected(tmp_path):
    env = run_env(tmp_path)
    target = tmp_path / 'unrelated'
    target.mkdir()
    (tmp_path / 'scratch' / f'genegalleon-{os.getuid()}').symlink_to(target)
    result = run(env, 'true')
    assert result.returncode == 1
    assert 'symlink' in result.stderr
    assert not list(target.iterdir())


def shell(script, env):
    return subprocess.run(['bash', '-c', f'source "{SUPPORT}/gg_util.sh"; {script}'],
                          env=env, text=True, capture_output=True)


def test_mount_and_mapping(tmp_path):
    env = run_env(tmp_path, GG_COMMON_TMP_ROOT=str(tmp_path / 'scratch'),
                  gg_workspace_dir='/workspace')
    result = shell('gg_configure_task_tmp_mount && '
                   'test "$SINGULARITYENV_GG_TMP_MOUNT" = /gg_tmp && '
                   'test "$APPTAINERENV_GG_TMP_HOST_ROOT" = "$GG_COMMON_TMP_ROOT"', env)
    assert result.returncode == 0, result.stderr
    env.update(GG_TMP_TASK_ROOT='/gg_tmp/run/work')
    result = shell('gg_task_tmp_path /workspace/output/tmp/task', env)
    assert result.stdout.strip() == '/gg_tmp/run/work/output/tmp/task'
    assert shell('gg_task_tmp_path /other/tmp', env).returncode != 0
    env['GG_COMMON_TMP_ROOT'] = 'workspace'
    result = shell('gg_task_tmp_path /workspace/output/tmp/task', env)
    assert result.stdout.strip() == '/workspace/output/tmp/task'


def test_env_resolved_at_launch_and_invalid_root_fails(tmp_path):
    env = run_env(tmp_path, GG_COMMON_TMP_ROOT='env', gg_workspace_dir='/workspace')
    env.pop('TMPDIR', None)
    assert shell('gg_configure_task_tmp_mount', env).returncode != 0
    env['TMPDIR'] = str(tmp_path / 'scratch')
    assert shell('gg_configure_task_tmp_mount', env).returncode == 0
    env['GG_COMMON_TMP_ROOT'] = str(tmp_path / 'missing')
    assert shell('gg_configure_task_tmp_mount', env).returncode != 0


def test_publication_failure_preserves_scratch(tmp_path):
    env = run_env(tmp_path, gg_workspace_dir=str(tmp_path))
    # Existing cp_out stages at the destination before renaming. A failed copy
    # must propagate and keep scratch, without replacing an existing output.
    destination = tmp_path / 'result'
    destination.write_text('original')
    result = run(env, f'set -e; source "{SUPPORT}/gg_util.sh"; '
                      f'cp_out "$GG_TMP_TASK_ROOT/missing" "{destination}"')
    assert result.returncode != 0
    assert destination.read_text() == 'original'
    assert len(records(tmp_path)) == 1


@pytest.mark.parametrize("metrics", [False, True])
def test_signal_is_not_success_even_when_child_trap_exits_zero(tmp_path, metrics):
    env = run_env(tmp_path)
    ready = tmp_path / 'ready'
    script = f'trap "exit 0" TERM; touch "{ready}"; while :; do sleep .1; done'
    command = ['bash', '-c', script]
    if metrics:
        command = [sys.executable, str(SUPPORT / 'resource_metrics.py'),
                   '--directory', str(tmp_path / 'metrics'), '--workflow', 'gg_gene_evolution',
                   '--runtime-id', 'test', '--server-id', 'test', '--', *command]
    child = subprocess.Popen([sys.executable, str(RUNNER), '--workflow', 'gg_gene_evolution',
                              '--', *command], env=env,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        for _ in range(200):
            if ready.exists():
                break
            time.sleep(.02)
        assert ready.exists()
        child.send_signal(signal.SIGTERM)
        stdout, stderr = child.communicate(timeout=10)
        assert child.returncode == 143, stdout + stderr
        assert len(records(tmp_path)) == 1
    finally:
        if child.poll() is None:
            child.kill()
            child.communicate()


@pytest.mark.parametrize("suffix", [":other", ",other", "\n"])
def test_canonical_mount_path_is_validated(tmp_path, suffix):
    unsafe = tmp_path / ("scratch" + suffix)
    unsafe.mkdir()
    alias = tmp_path / 'safe-alias'
    alias.symlink_to(unsafe)
    env = run_env(tmp_path, GG_COMMON_TMP_ROOT=str(alias), gg_workspace_dir='/workspace')
    result = shell('gg_configure_task_tmp_mount', env)
    assert result.returncode != 0


def test_other_gene_mode_scratch_is_not_deleted(tmp_path):
    env = run_env(tmp_path, mode_gene_evolution='orthogroup')
    assert run(env, 'exit 1').returncode == 1
    record = records(tmp_path)[0]
    env['mode_gene_evolution'] = 'query2family'
    assert run(env, 'true').returncode == 0
    assert record.exists()


def test_corrupt_record_is_preserved_without_blocking_new_job(tmp_path):
    env = run_env(tmp_path)
    assert run(env, 'exit 1').returncode == 1
    record = records(tmp_path)[0]
    record.write_text('[]')
    result = run(env, 'true')
    assert result.returncode == 0, result.stderr
    assert record.read_text() == '[]'


def test_non_exec_adapter_rejects_external_scratch(tmp_path):
    script = tmp_path / 'gg_probe_core.sh'
    script.write_text('exit 0\n')
    runtime = tmp_path / 'runtime'
    called = tmp_path / 'called'
    runtime.write_text(f'#!/bin/sh\ntouch "{called}"\nexit 0\n')
    runtime.chmod(0o755)
    env = run_env(tmp_path, GG_COMMON_TMP_ROOT=str(tmp_path / 'scratch'))
    result = shell(f'singularity_command=("{runtime}" shell); '
                   f'gg_run_container_shell_script image "{script}"', env)
    assert result.returncode != 0
    assert not called.exists()


@pytest.mark.parametrize('key,value', [
    ('delete_tmp_dir', 'invalid'), ('delete_preexisting_tmp_dir', 'invalid'),
    ('GG_ARRAY_TASK_ID', '0'), ('mode_gene_evolution', 'invalid'),
    ('gene_family_tmp_max_bytes', '-1'),
])
def test_invalid_config_never_prunes_previous_scratch(tmp_path, key, value):
    env = run_env(tmp_path)
    assert run(env, 'exit 1').returncode == 1
    record = records(tmp_path)[0]
    before = record.read_bytes()
    env[key] = value
    assert run(env, 'true').returncode != 0
    assert record.read_bytes() == before


@pytest.mark.parametrize('limit,value', [
    ('gene_family_tmp_max_bytes', '1'), ('gene_family_tmp_max_files', '1'),
])
def test_retention_budget_removes_only_owned_idle_runs(tmp_path, limit, value):
    env = run_env(tmp_path)
    assert run(env, 'echo payload > "$GG_TMP_TASK_ROOT/data"; exit 1').returncode == 1
    previous = records(tmp_path)[0]
    env.update({limit: value, 'GG_ARRAY_TASK_ID': '2'})
    assert run(env, 'true').returncode == 0
    assert not previous.exists()


def test_retention_age_expires_old_inactive_run(tmp_path):
    env = run_env(tmp_path)
    assert run(env, 'exit 1').returncode == 1
    record = records(tmp_path)[0]
    data = json.loads(record.read_text())
    data['updated'] = time.time() - 10 * 86400
    record.write_text(json.dumps(data))
    env['GG_ARRAY_TASK_ID'] = '2'
    assert run(env, 'true').returncode == 0
    assert not record.exists()


def test_killed_metrics_parent_cannot_unlock_live_core_scratch(tmp_path):
    env = run_env(tmp_path)
    ready = tmp_path / 'pids'
    script = f'echo "$PPID $$" > "{ready}"; while :; do sleep .1; done'
    command = [sys.executable, str(RUNNER), '--workflow', 'gg_gene_evolution', '--',
               sys.executable, str(SUPPORT / 'resource_metrics.py'),
               '--directory', str(tmp_path / 'metrics'), '--workflow', 'gg_gene_evolution',
               '--runtime-id', 'test', '--server-id', 'test', '--', 'bash', '-c', script]
    core_pid = None
    with (tmp_path / 'log').open('w') as log:
        child = subprocess.Popen(command, env=env, stdout=log, stderr=log)
        try:
            for _ in range(200):
                if ready.exists() and ready.read_text().strip():
                    break
                time.sleep(.02)
            metrics_pid, core_pid = map(int, ready.read_text().split())
            record = records(tmp_path)[0]
            os.kill(metrics_pid, signal.SIGKILL)
            assert child.wait(timeout=10) == 137
            # The computation survived its monitor and still owns its scratch.
            os.kill(core_pid, 0)
            assert run(env, 'true').returncode == 0
            assert record.exists()
        finally:
            if core_pid is not None:
                try:
                    os.killpg(core_pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
            if child.poll() is None:
                child.kill()
                child.wait(timeout=10)


@pytest.mark.parametrize('mode', ['query2family', 'orthogroup'])
def test_changed_family_assignment_preserves_other_family_scratch(tmp_path, mode):
    workspace = tmp_path / 'workspace'
    env = run_env(tmp_path, gg_workspace_dir=str(workspace), mode_gene_evolution=mode)
    if mode == 'query2family':
        inputs = workspace / 'input/query_gene'
        inputs.mkdir(parents=True)
        (inputs / 'family_b').write_text('query')
    else:
        table = workspace / 'output/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv'
        table.parent.mkdir(parents=True)
        table.write_text('Orthogroup\tSpecies\nOG0001\t1\n')
    assert run(env, 'exit 1').returncode == 1
    record = records(tmp_path)[0]
    if mode == 'query2family':
        (inputs / 'family_a').write_text('new query')
    else:
        table.write_text('Orthogroup\tSpecies\nOG0002\t1\n')
    assert run(env, 'true').returncode == 0
    assert record.exists()


@pytest.mark.parametrize("interrupt", [False, True])
def test_success_cleanup_waits_for_retention_scan(tmp_path, interrupt):
    import fcntl

    env = run_env(tmp_path)
    ready, finish = tmp_path / 'ready', tmp_path / 'finish'
    command = [sys.executable, str(RUNNER), '--workflow', 'gg_gene_evolution', '--',
               'bash', '-c', f'touch "{ready}"; while test ! -e "{finish}"; do sleep .02; done']
    child = subprocess.Popen(command, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        for _ in range(200):
            if ready.exists():
                break
            time.sleep(.02)
        assert ready.exists()
        record = records(tmp_path)[0]
        # Another invocation is scanning idle runs under the shared scope lock.
        with (record.parent.parent / 'scope.lock').open('rb') as scope:
            fcntl.flock(scope, fcntl.LOCK_EX)
            finish.touch()
            time.sleep(.3)
            assert record.exists(), 'Completion removed a run during a retention scan'
            if interrupt:
                child.send_signal(signal.SIGTERM)
                time.sleep(.1)
        stdout, stderr = child.communicate(timeout=10)
        assert child.returncode == (143 if interrupt else 0), stdout + stderr
        assert bool(records(tmp_path)) == interrupt
    finally:
        finish.touch()
        if child.poll() is None:
            child.kill()
        child.communicate(timeout=10)

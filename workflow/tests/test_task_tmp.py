"""Behavior tests for external, disposable workflow scratch storage."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

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


def test_input_array_workers_with_separate_scratch_roots(tmp_path, monkeypatch):
    # Exercise the full prepare/parallel workers/finalize core flow with the
    # existing deterministic toolchain fixture. Each job has isolated scratch,
    # as on separate compute nodes; shared plans and shards must still work.
    import shlex
    import test_gg_input_generation_end_to_end as integration

    original_core = integration.CORE_PATH
    original_env = integration._core_env
    wrapper = tmp_path / 'external-core.sh'
    wrapper.write_text(
        'exec ' + shlex.join([sys.executable, str(RUNNER), '--workflow',
                             'gg_input_generation', '--', 'bash', str(original_core)]) + '\n'
    )

    def external_env(*args, **kwargs):
        env = original_env(*args, **kwargs)
        scratch = tmp_path / ('scratch-' + env['input_generation_mode'] + '-' + env.get('GG_ARRAY_TASK_ID', '1'))
        scratch.mkdir(exist_ok=True)
        env.update(GG_COMMON_TMP_ROOT=str(scratch), GG_TMP_MOUNT=str(scratch),
                   GG_TMP_WORKSPACE_ID=env['gg_workspace_dir'])
        return env

    monkeypatch.setattr(integration, 'CORE_PATH', wrapper)
    monkeypatch.setattr(integration, '_core_env', external_env)
    integration.test_gg_input_generation_array_mode_end_to_end_with_parallel_workers(tmp_path)
    assert not list(tmp_path.glob('scratch-*/genegalleon-*/*/*/run-*'))

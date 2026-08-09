import os
import shutil
import subprocess
from pathlib import Path

import pytest


WORKFLOW_DIR = Path(__file__).resolve().parents[1]
LOCK_HELPER = WORKFLOW_DIR / "support" / "gg_shared_lock.sh"


def run_bash(script: str, env_updates=None) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    for key in (
        "GG_ARRAY_TASK_COUNT",
        "SLURM_ARRAY_TASK_COUNT",
        "PBS_ARRAY_TASK_COUNT",
        "SGE_TASK_FIRST",
        "SGE_TASK_LAST",
        "SGE_TASK_STEPSIZE",
    ):
        env.pop(key, None)
    env.update(env_updates or {})
    return subprocess.run(
        ["bash", "-c", script],
        text=True,
        capture_output=True,
        check=False,
        env=env,
    )


def run_finalizer(state_root: Path, task_id: int) -> subprocess.CompletedProcess[str]:
    script = f'''
source "{LOCK_HELPER}"
if gg_array_finalizer_claim "{state_root}" summary 3; then
  echo CLAIMED
  gg_array_finalizer_complete
else
  status=$?
  echo "WAIT:${{status}}"
fi
'''
    env = os.environ.copy()
    env.update(
        {
            "GG_JOB_ID": "job-42",
            "GG_ARRAY_TASK_ID": str(task_id),
            "GG_ARRAY_TASK_COUNT": "3",
        }
    )
    return subprocess.run(
        ["bash", "-c", script],
        text=True,
        capture_output=True,
        check=False,
        env=env,
    )


@pytest.mark.parametrize(
    ("environment", "expected"),
    [
        ({"SLURM_ARRAY_TASK_COUNT": "7"}, "7"),
        ({"PBS_ARRAY_TASK_COUNT": "5"}, "5"),
        ({"SGE_TASK_FIRST": "2", "SGE_TASK_LAST": "10", "SGE_TASK_STEPSIZE": "2"}, "5"),
    ],
)
def test_array_expected_task_count_uses_scheduler_metadata(environment, expected):
    completed = run_bash(
        f'source "{LOCK_HELPER}"\ngg_array_expected_task_count 99\n',
        environment,
    )
    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == expected


@pytest.mark.skipif(shutil.which("flock") is None, reason="flock is required for array finalization")
def test_only_last_ready_array_task_claims_finalizer(tmp_path: Path):
    state_root = tmp_path / "finalizers"
    first = run_finalizer(state_root, 1)
    second = run_finalizer(state_root, 2)
    last = run_finalizer(state_root, 3)
    repeated = run_finalizer(state_root, 3)

    assert first.returncode == 0, first.stderr
    assert second.returncode == 0, second.stderr
    assert last.returncode == 0, last.stderr
    assert repeated.returncode == 0, repeated.stderr
    assert "WAIT:1" in first.stdout
    assert "WAIT:1" in second.stdout
    assert "CLAIMED" in last.stdout
    assert "WAIT:1" in repeated.stdout
    assert (state_root / "summary" / "job-42" / "done").is_file()


@pytest.mark.skipif(shutil.which("flock") is None, reason="flock is required for array finalization")
def test_local_single_task_finalizer_does_not_reuse_completed_state(tmp_path: Path):
    script = f'''
source "{LOCK_HELPER}"
if gg_array_finalizer_claim "{tmp_path / 'finalizers'}" summary 1; then
  echo CLAIMED
  gg_array_finalizer_complete
fi
'''
    first = run_bash(script)
    second = run_bash(script)

    assert first.returncode == 0, first.stderr
    assert second.returncode == 0, second.stderr
    assert "CLAIMED" in first.stdout
    assert "CLAIMED" in second.stdout

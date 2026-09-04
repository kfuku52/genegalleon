import os
import subprocess
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"


def run_shell(script, *args, env=None):
    return subprocess.run(
        ["bash", "-c", 'set -euo pipefail\nsource "$1/gg_util.sh"\n' + script,
         "audit", str(SUPPORT), *map(str, args)],
        env={"PATH": os.environ["PATH"], **(env or {})},
        capture_output=True, text=True, timeout=15,
    )


def test_cleanup_preserves_output_published_before_removal(tmp_path):
    target = tmp_path / "stat_branch"
    target.mkdir()
    script = r'''
audit_target="$2/stat_branch"
publish() { printf 'result\n' > "$audit_target/new_family.tsv"; }
find() {
  if [[ "$1" == "$audit_target" && "$*" == *"-print -quit"* ]]; then
    local snapshot
    snapshot=$(command find "$@")
    publish
    printf '%s' "$snapshot"
  else
    command find "$@"
  fi
}
rmdir() { publish; command rmdir "$@"; }
remove_empty_subdirs "$2"
'''
    result = run_shell(script, tmp_path)
    assert result.returncode == 0, result.stderr
    assert (target / "new_family.tsv").read_text() == "result\n"


def test_cleanup_removes_only_empty_directories(tmp_path):
    (tmp_path / "empty").mkdir()
    (tmp_path / "kept").mkdir()
    (tmp_path / "kept" / ".hidden").write_text("result")
    result = run_shell('remove_empty_subdirs "$2"', tmp_path)
    assert result.returncode == 0, result.stderr
    assert not (tmp_path / "empty").exists()
    assert (tmp_path / "kept" / ".hidden").read_text() == "result"


def test_slurm_array_finalizes_once_using_forwarded_common_job_id(tmp_path):
    script = r'''
gg_normalize_scheduler_env >/dev/null
gg_workspace_dir="$2"
gg_workflow_dir="$1/.."
gg_container_image_path="$2/runtime.sif"
set_singularityenv >/dev/null
[[ "$GG_JOB_ID" == "$SLURM_JOB_ID" ]]
[[ "$GG_ARRAY_JOB_ID" == 101 ]]
[[ "$APPTAINERENV_GG_ARRAY_JOB_ID" == 101 ]]
GG_ARRAY_JOB_ID=$SINGULARITYENV_GG_ARRAY_JOB_ID
unset SLURM_ARRAY_JOB_ID
if gg_array_finalizer_claim "$2/finalizers" summary 3; then
  printf 'summary\n' >> "$2/summary-runs"
  gg_array_finalizer_complete
fi
'''
    (tmp_path / "runtime.sif").touch()
    for task, job in ((1, 102), (2, 103), (3, 101), (3, 101)):
        result = run_shell(script, tmp_path, env={
            "SLURM_JOB_ID": str(job), "SLURM_ARRAY_JOB_ID": "101",
            "SLURM_ARRAY_TASK_ID": str(task), "SLURM_ARRAY_TASK_COUNT": "3",
            "SLURM_CPUS_PER_TASK": "1",
        })
        assert result.returncode == 0, result.stderr
    assert (tmp_path / "summary-runs").read_text() == "summary\n"
    assert len(list((tmp_path / "finalizers").rglob("done"))) == 1


@pytest.mark.parametrize("memory, overrides, expected", [
    ({"SLURM_MEM_PER_NODE": "8192"}, {}, 8),
    ({"SLURM_MEM_PER_NODE": "12345"}, {}, 12),
    ({"SLURM_MEM_PER_CPU": "1536"}, {}, 12),
    ({"SLURM_MEM_PER_NODE": "8192"}, {"GG_MEM_TOTAL_GB": "6"}, 6),
    ({"SLURM_MEM_PER_NODE": "8192"}, {"GG_MEM_PER_CPU_GB": "2"}, 16),
    ({"SLURM_MEM_PER_NODE": "8192"}, {"MEM_PER_HOST": "7"}, 7),
])
def test_slurm_memory_uses_allocation_and_respects_explicit_overrides(memory, overrides, expected):
    result = run_shell(
        'gg_normalize_scheduler_env >/dev/null\nprintf "%s %s\\n" "$GG_MEM_TOTAL_GB" "$GG_MEM_TOOL_GB"',
        env={"SLURM_JOB_ID": "123", "SLURM_CPUS_PER_TASK": "8", **memory, **overrides},
    )
    assert result.returncode == 0, result.stderr
    total, tool = map(int, result.stdout.split())
    assert total == expected
    assert 0 < tool <= total


@pytest.mark.parametrize("status", [0, 23])
def test_disabled_semaphore_preserves_command_status(tmp_path, status):
    result = run_shell('gg_run_with_shared_semaphore "$2" 0 audit bash -c "exit $3"', tmp_path, status)
    assert result.returncode == status, result.stderr

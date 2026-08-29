import os
import shutil
import subprocess
from pathlib import Path

import pytest

WORKFLOW_DIR = Path(__file__).resolve().parents[1]
LOCK_HELPER = WORKFLOW_DIR / "support" / "gg_shared_lock.sh"
TRANSCRIPTOME_CORE = WORKFLOW_DIR / "core" / "gg_transcriptome_generation_core.sh"


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


@pytest.mark.skipif(shutil.which("flock") is None, reason="flock is required for array finalization")
def test_different_singleton_jobs_serialize_one_shared_stage_transaction(tmp_path: Path):
    state_root = tmp_path / "finalizers"
    artifact = tmp_path / "annotation_summary" / "summary.pdf"
    trace = tmp_path / "trace.log"
    script = f'''
source "{LOCK_HELPER}"
claimed=0
transaction_locked=0
cleanup() {{
  status=$?
  if [[ $transaction_locked -eq 1 ]]; then
    gg_stage_transaction_lock_release || true
  fi
  if [[ $claimed -eq 1 ]]; then
    gg_array_finalizer_release || true
  fi
  exit "$status"
}}
if gg_array_finalizer_claim "{state_root}" summary 1; then
  claimed=1
  trap cleanup EXIT
else
  exit $?
fi
gg_stage_transaction_lock_acquire "{state_root}" summary
transaction_locked=1
if [[ -s "{artifact}" ]]; then
  printf 'skip:%s\n' "$GG_JOB_ID" >> "{trace}"
else
  stage=$(mktemp -d "{tmp_path}/summary-stage.XXXXXX")
  mkdir -p "$stage/output"
  printf 'start:%s\n' "$GG_JOB_ID" >> "{trace}"
  sleep 0.2
  printf 'complete\n' > "$stage/output/summary.pdf"
  mv -- "$stage/output" "{artifact.parent}"
  rmdir -- "$stage"
  printf 'finish:%s\n' "$GG_JOB_ID" >> "{trace}"
fi
gg_stage_transaction_lock_release
transaction_locked=0
gg_array_finalizer_complete
claimed=0
trap - EXIT
'''

    processes = []
    for job_id, task_id in (("replacement-a", "102"), ("replacement-b", "131")):
        env = os.environ.copy()
        env.update(
            {
                "GG_SCHEDULER_KIND": "slurm",
                "GG_JOB_ID": job_id,
                "GG_ARRAY_TASK_ID": task_id,
                "GG_ARRAY_TASK_COUNT": "1",
                "GG_LOCK_POLL_SECONDS": "1",
            }
        )
        processes.append(
            subprocess.Popen(
                ["bash", "-c", script],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=env,
            )
        )

    completed = [process.communicate(timeout=10) for process in processes]

    assert all(process.returncode == 0 for process in processes), completed
    assert artifact.read_text(encoding="utf-8") == "complete\n"
    trace_lines = trace.read_text(encoding="utf-8").splitlines()
    assert len([line for line in trace_lines if line.startswith("start:")]) == 1
    assert len([line for line in trace_lines if line.startswith("finish:")]) == 1
    assert len([line for line in trace_lines if line.startswith("skip:")]) == 1
    assert trace_lines[0].startswith("start:")
    assert trace_lines[1].startswith("finish:")
    assert trace_lines[2].startswith("skip:")


def test_transcriptome_summary_uses_private_staging_and_transactional_publication():
    core = TRANSCRIPTOME_CORE.read_text(encoding="utf-8")
    summary_start = core.index("task='Multispecies summary'")
    summary_end = core.index("if [[ ${remove_amalgkit_fastq_after_completion}", summary_start)
    summary = core[summary_start:summary_end]

    lock_position = summary.index("gg_stage_transaction_lock_acquire")
    provenance_position = summary.index("gg_artifact_prepare_stage")
    publish_position = summary.index("mv_out_replace_dir")
    record_position = summary.index("gg_artifact_record")

    assert lock_position < provenance_position < publish_position < record_position
    assert ".annotation_summary.gg-work.XXXXXX" in summary
    assert 'cd "${transcriptome_summary_stage_dir}"' in summary
    assert 'rm -rf -- "${transcriptome_summary_output_dir}"' not in summary

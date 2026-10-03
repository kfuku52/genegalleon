import os
import shlex
import subprocess
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"


def run_shell(script, *args, env=None, timeout=15):
    return subprocess.run(
        ["bash", "-c", 'set -euo pipefail\nsource "$1/gg_util.sh"\n' + script,
         "audit", str(SUPPORT), *map(str, args)],
        env={"PATH": os.environ["PATH"], **(env or {})},
        capture_output=True, text=True, timeout=timeout,
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



def test_parameter_snapshot_updates_equal_size_edits_and_keeps_identical_files(tmp_path):
    core = (SUPPORT.parent / "core/gg_gene_evolution_core.sh").read_text()
    start = core.index("# Copy parameter files and codes to ")
    end = core.index('\ncd "${gg_workspace_dir}"', start)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    traits, tree, pruned = [inputs / name for name in ("trait.tsv", "species.nwk", "pruned.nwk")]
    traits.write_text("name\ttrait\nA\t0\n")
    tree.write_text("(A:1,B:1);\n")
    pruned.write_text("(A:1,B:1);\n")
    saved = tmp_path / "parameters"
    setup = "\n".join(f"{key}={shlex.quote(str(value))}" for key, value in {
        "file_og_parameters_dir": saved, "file_sp_trait": traits,
        "species_tree": tree, "species_tree_pruned": pruned,
    }.items()) + "\n"
    script = setup + core[start:end]
    for _ in range(2):
        result = run_shell(script, timeout=30)
        assert result.returncode == 0, result.stdout + result.stderr
        traits.write_text("name\ttrait\nA\t1\n")
        tree.write_text("(A:2,B:1);\n")
    assert (saved / traits.name).read_bytes() == traits.read_bytes()
    assert (saved / tree.name).read_bytes() == tree.read_bytes()
    before = {path: path.stat().st_mtime_ns for path in saved.iterdir() if path.is_file()}
    result = run_shell(script, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    assert before == {path: path.stat().st_mtime_ns for path in before}


def test_failed_site_summary_rebuild_preserves_completed_result_and_manifest(tmp_path):
    core_text = (SUPPORT.parent / "core/gg_gene_summary_core.sh").read_text()
    start = core_text.index("run_csubst_site_convergence_summary_for_source() {")
    end = core_text.index("\n}\n", start) + 3
    workspace = tmp_path / "workspace"
    families = workspace / "output/query2family"
    families.mkdir(parents=True)
    trait = workspace / "input/trait.tsv"
    trait.parent.mkdir()
    trait.write_text("name\ttrait\nA\t0\n")
    orthofinder = workspace / "output/orthofinder"
    orthofinder.mkdir()
    (orthofinder / "inputs.tsv").write_text("data\n")
    child_dir = tmp_path / "child"
    child_dir.mkdir()
    child = child_dir / "gg_convergent_sites_core.sh"
    child.write_text('set -euo pipefail\nprintf "valid\\n" > "$dir_out/prior.tsv"\n')
    out = workspace / "output/summary/sites"
    variables = dict(gg_workspace_dir=workspace, gg_workspace_input_dir=workspace / "input",
                     gg_workspace_output_dir=workspace / "output", gg_support_dir=SUPPORT,
                     gg_core_dir=child_dir, summary_output_dir=workspace / "output/summary",
                     gene_family_source="query2family", dir_gene_family=families,
                     csubst_site_output_dir=out, csubst_site_trait_file=trait,
                     csubst_site_orthofinder_dir=orthofinder, artifact_stale_policy="rebuild",
                     run_csubst_site_convergence_summary=1, csubst_site_arity_range=2,
                     csubst_site_trait="all", csubst_site_skip_lower_order="yes",
                     csubst_site_min_fg_stem_ratio=0, csubst_site_min_ocn_any2spe=1,
                     csubst_site_min_omega_c_any2spe=1, csubst_site_min_ocn_cod=1,
                     csubst_site_max_candidates_per_arity=1, csubst_site_nonsyn_recode="no")
    setup = "\n".join(f"{key}={shlex.quote(str(value))}" for key, value in variables.items()) + "\n"
    script = setup + core_text[start:end] + "\nrun_csubst_site_convergence_summary_for_source\n"
    result = run_shell(script, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    manifest = workspace / "output/summary/artifact_provenance/query2family.csubst_site.json"
    before = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in [out / "prior.tsv", manifest]}
    trait.write_text("name\ttrait\nA\t1\n")
    child.write_text('printf "partial\\n" > "$dir_out/partial.tsv"\nexit 23\n')
    result = run_shell(script, timeout=60)
    assert result.returncode == 23, result.stdout + result.stderr
    assert before == {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in before}
    assert not (out / "partial.tsv").exists()
    assert not list(out.parent.glob("sites.rebuild.*"))
    child.write_text('printf "new\\n" > "$dir_out/current.tsv"\n')
    result = run_shell(script, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert not (out / "prior.tsv").exists()
    assert (out / "current.tsv").read_text() == "new\n"
    assert manifest.read_bytes() != before[manifest][0]
    assert not list(out.parent.glob("sites.rebuild.*"))

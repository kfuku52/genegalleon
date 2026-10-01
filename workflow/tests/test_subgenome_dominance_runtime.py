import json
import os
import subprocess
from pathlib import Path

from workflow.tests.test_subgenome_dominance import fixture_inputs

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_standalone_stage_reuse_and_stale_input_preserve_results(tmp_path):
    workspace = tmp_path / "workspace"
    fixture_inputs(workspace / "input")
    env = {**os.environ, "gg_workspace_dir": str(workspace), "GG_COMMON_TMP_ROOT": "workspace",
           "GG_ARRAY_TASK_ID": "1", "GG_JOB_ID": "subgenome_test", "GG_TASK_CPUS": "1", "GG_MEM_PER_CPU_GB": "2",
           "genome_evolution_mode": "subgenome", "artifact_stale_policy": "stop",
           "subgenome_manifest": str(workspace / "input/manifest.tsv"), "subgenome_bootstrap_replicates": "200"}
    tree = workspace / "output/species_tree/precious.txt"
    tree.parent.mkdir(parents=True)
    tree.write_text("preserve unrelated stages")
    def run():
        return subprocess.run(["bash", str(REPO_ROOT / "workflow/core/gg_genome_evolution_core.sh")],
                              cwd=REPO_ROOT, env=env, capture_output=True, text=True, timeout=90)
    first = run()
    assert first.returncode == 0, first.stdout + first.stderr
    output = workspace / "output/genome_evolution/subgenome_dominance/fixture/summary.json"
    assert json.loads(output.read_text())["statistics"][1]["effect"] == 1
    saved = output.read_bytes()
    for name in ("contrasts.png", "contrasts.svg", "contrasts_absolute.png", "contrasts_absolute.svg"):
        assert (output.parent / name).stat().st_size > 0
    mtime = output.stat().st_mtime_ns
    second = run()
    assert second.returncode == 0, second.stdout + second.stderr
    assert output.stat().st_mtime_ns == mtime
    assert tree.read_text() == "preserve unrelated stages"
    expression = workspace / "input/expr.tsv"
    expression.write_text(expression.read_text().replace("g0A\t2\t2\t2", "g0A\t4\t4\t4"))
    stale = run()
    assert stale.returncode != 0
    assert output.read_bytes() == saved
    env["artifact_stale_policy"] = "rebuild"
    rebuilt = run()
    assert rebuilt.returncode == 0, rebuilt.stdout + rebuilt.stderr
    assert output.read_bytes() != saved

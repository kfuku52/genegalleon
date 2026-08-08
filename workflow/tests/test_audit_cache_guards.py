import subprocess
import sys
from pathlib import Path

import pandas as pd

WORKFLOW_DIR = Path(__file__).resolve().parents[1]
SCRIPT = WORKFLOW_DIR / "support" / "audit_cache_guards.py"
BASELINE = Path(__file__).resolve().parent / "data" / "artifact_cache_guard_baseline.txt"


def run_cli(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *map(str, args)],
        text=True,
        capture_output=True,
        check=False,
    )


def test_cache_guard_audit_detects_existence_only_stage(tmp_path):
    core = tmp_path / "workflow" / "core"
    core.mkdir(parents=True)
    script = core / "example_core.sh"
    script.write_text(
        'task="Example"\n'
        'if [[ ! -s "${file_example}" && ${run_example} -eq 1 ]]; then\n'
        '  gg_step_start "${task}"\n'
        "else\n"
        '  gg_step_skip "${task}"\n'
        "fi\n",
        encoding="utf-8",
    )
    report = tmp_path / "report.tsv"

    completed = run_cli(
        core,
        "--source-root",
        tmp_path,
        "--output-tsv",
        report,
        "--fail-on-findings",
    )
    assert completed.returncode == 1
    frame = pd.read_csv(report, sep="\t", keep_default_na=False)
    assert frame.loc[0, "task"] == "Example"
    assert frame.loc[0, "guard_type"] == "existence_only"
    assert frame.loc[0, "output_variables"] == "file_example"


def test_cache_guard_audit_does_not_flag_provenance_controlled_stage(tmp_path):
    core = tmp_path / "workflow" / "core"
    core.mkdir(parents=True)
    script = core / "example_core.sh"
    script.write_text(
        'task="Example"\n'
        'if gg_artifact_needs_run "${provenance_args[@]}"; then\n'
        "  example_needs_update=1\n"
        "fi\n"
        "if [[ ${example_needs_update} -eq 1 && ${run_example} -eq 1 ]]; then\n"
        '  gg_step_start "${task}"\n'
        "fi\n",
        encoding="utf-8",
    )
    report = tmp_path / "report.tsv"

    completed = run_cli(
        core,
        "--source-root",
        tmp_path,
        "--output-tsv",
        report,
        "--fail-on-findings",
    )
    assert completed.returncode == 0
    frame = pd.read_csv(report, sep="\t", keep_default_na=False)
    assert frame.empty


def test_repository_has_no_new_unprovenanced_cache_guards(tmp_path):
    report = tmp_path / "report.tsv"
    completed = run_cli(
        WORKFLOW_DIR / "core",
        "--source-root",
        WORKFLOW_DIR.parent,
        "--output-tsv",
        report,
        "--baseline",
        BASELINE,
    )
    assert completed.returncode == 0, completed.stderr
    assert "unapproved_findings=0" in completed.stderr

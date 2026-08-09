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


def test_cache_guard_audit_detects_mtime_and_count_only_stages(tmp_path):
    core = tmp_path / "workflow" / "core"
    core.mkdir(parents=True)
    script = core / "example_core.sh"
    script.write_text(
        'task="Mtime cache"\n'
        'if [[ "${file_input}" -nt "${file_output}" ]]; then\n'
        '  gg_step_start "${task}"\n'
        "fi\n"
        'task="Count cache"\n'
        'if [[ ${num_inputs} -ne ${num_outputs} && ${run_count} -eq 1 ]]; then\n'
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
    assert completed.returncode == 1
    frame = pd.read_csv(report, sep="\t", keep_default_na=False)
    assert frame[["task", "guard_type"]].to_dict(orient="records") == [
        {"task": "Mtime cache", "guard_type": "mtime_only"},
        {"task": "Count cache", "guard_type": "count_only"},
    ]


def test_cache_guard_audit_detects_early_exit_on_existing_output(tmp_path):
    core = tmp_path / "workflow" / "core"
    core.mkdir(parents=True)
    script = core / "example_core.sh"
    script.write_text(
        'task="Internal cache"\n'
        'if [[ -s "${file_output}" ]]; then\n'
        '  echo "Skipping cached output"\n'
        "  return 0\n"
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
    assert frame.loc[0, "guard_type"] == "early_exit_existence_only"
    assert frame.loc[0, "output_variables"] == "file_output"


def test_cache_guard_audit_allows_explicitly_audited_early_exit(tmp_path):
    core = tmp_path / "workflow" / "core"
    core.mkdir(parents=True)
    script = core / "example_core.sh"
    script.write_text(
        '# gg-cache-guard: audited - source cache is content-addressed\n'
        'if [[ -s "${file_output}" ]]; then\n'
        "  return 0\n"
        "fi\n",
        encoding="utf-8",
    )

    completed = run_cli(
        core,
        "--source-root",
        tmp_path,
        "--fail-on-findings",
    )
    assert completed.returncode == 0


def test_cache_guard_audit_detects_python_early_exit_on_existing_output(tmp_path):
    support = tmp_path / "workflow" / "support"
    support.mkdir(parents=True)
    script = support / "example.py"
    script.write_text(
        "def build(output_path, overwrite):\n"
        "    if output_path.exists() and not overwrite:\n"
        "        return {'status': 'skip'}\n",
        encoding="utf-8",
    )
    report = tmp_path / "report.tsv"

    completed = run_cli(
        support,
        "--source-root",
        tmp_path,
        "--output-tsv",
        report,
        "--fail-on-findings",
    )
    assert completed.returncode == 1
    frame = pd.read_csv(report, sep="\t", keep_default_na=False)
    assert frame.loc[0, "guard_type"] == "python_early_exit_existence_only"
    assert frame.loc[0, "task"] == "build"
    assert frame.loc[0, "output_variables"] == "output_path"


def test_repository_has_no_new_unprovenanced_cache_guards(tmp_path):
    report = tmp_path / "report.tsv"
    completed = run_cli(
        WORKFLOW_DIR / "core",
        WORKFLOW_DIR / "support",
        "--source-root",
        WORKFLOW_DIR.parent,
        "--output-tsv",
        report,
        "--baseline",
        BASELINE,
    )
    assert completed.returncode == 0, completed.stderr
    assert "unapproved_findings=0" in completed.stderr

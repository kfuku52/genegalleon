"""Exercise the refinement finisher with stubbed tools and strict Bash."""

import os
import subprocess
from pathlib import Path

import pytest


@pytest.fixture
def run_finisher(tmp_path):
    core = Path(__file__).resolve().parents[1] / "core/gg_input_generation_core.sh"
    signature = "finish_gene_model_refinement() {"
    body = core.read_text().split(signature, 1)[1].split("\n}\n", 1)[0]
    log = tmp_path / "commands.txt"
    script = """set -euo pipefail
python() { printf '%s\\n' "$*" >> "$COMMAND_LOG"; }
run_refinement_omark() { printf 'omark\\n' >> "$COMMAND_LOG"; return "$OMARK_STATUS"; }
annotate_gene_model_rescue_swissprot() { printf 'audit %s\\n' "$1" >> "$COMMAND_LOG"; }
ensure_shared_busco_lineage_ready() { busco_lineage_resolved=embryophyta_odb12; }
ensure_busco_download_path() { printf '%s\\n' /db; }
gg_memory_parallel_job_cap() { printf '%s\\n' 2; }
""" + signature + body + "\n}\nfinish_gene_model_refinement\n"

    def run(**settings):
        log.write_text("")
        env = dict(
            PATH=os.environ["PATH"], COMMAND_LOG=str(log), OMARK_STATUS="0",
            run_species_busco="1", run_species_omark="0", gg_support_dir="/support",
            run_gene_model_rescue_swissprot="0", gene_model_rescue_dir=str(tmp_path / "rescue"),
            gene_model_refinement_rescue_dir="", gene_model_refinement_inputs="", gene_model_rescue_swissprot_dir="",
            gene_model_refinement_dir=str(tmp_path / "refinement"), species_cds_dir="/all-cds",
            GG_TASK_CPUS="8", species_busco_parallel_jobs="auto", task_plan_output="/plan",
            gg_workspace_dir="/workspace", GG_MEM_TOOL_GB="64", species_busco_memory_gb_per_job="16",
            busco_lineage_resolved="",
        )
        result = subprocess.run(["bash", "-c", script], env={**env, **settings}, capture_output=True, text=True)
        return result, log.read_text().splitlines()

    return run


@pytest.mark.parametrize("enabled", ["0", "1"], ids=["busco-disabled", "busco-enabled"])
@pytest.mark.parametrize("trailing_slash", ["", "/"], ids=["no-slash", "trailing-slash"])
def test_finisher_wires_review_and_bounds_busco_resources(tmp_path, run_finisher, enabled, trailing_slash):
    refinement = tmp_path / "refinement"
    result, commands = run_finisher(
        run_species_busco=enabled, gene_model_refinement_dir=str(refinement) + trailing_slash,
    )
    assert result.returncode == 0, result.stderr
    assert " finalize " in commands[0] and " qc " in commands[1]
    assert "plot_gene_model_refinement.py" in commands[2] and "--cds-dir /all-cds" in commands[2]
    assert f"--report {refinement}.review " in commands[2]
    if enabled == "1":
        assert "gene_model_refinement_busco.py" in commands[3]
        assert "--jobs 2 --cpus 4" in commands[3]
        assert f"--report {refinement}.review/busco" in commands[3]
    else:
        assert len(commands) == 3


@pytest.mark.parametrize("enabled", [False, True])
def test_finisher_reuses_swissprot_annotation_and_forwards_verified_audit(tmp_path, run_finisher, enabled):
    anchor = tmp_path / "rescue"
    (anchor / "augmented").mkdir(parents=True)
    (anchor / "augmented/receipt.json").write_text("{}")
    result, commands = run_finisher(run_gene_model_rescue_swissprot=str(int(enabled)))
    assert result.returncode == 0, result.stderr
    assert any(c.startswith("audit ") for c in commands) == enabled
    comparison = next(c for c in commands if "gene_model_refinement_busco.py" in c)
    assert (f"--rescue-swissprot-dir {anchor}.swissprot" in comparison) == enabled
    assert ("--no-rescue-swissprot" in comparison) == (not enabled)


@pytest.mark.parametrize("status", [0, 17])
def test_finisher_runs_omark_after_qc_and_propagates_failure(run_finisher, status):
    result, commands = run_finisher(run_species_omark="1", OMARK_STATUS=str(status))
    assert result.returncode == status, result.stderr
    assert " finalize " in commands[0] and " qc " in commands[1]
    assert commands[2] == "omark"
    if status == 0:
        assert "plot_gene_model_refinement.py" in commands[3]
    else:
        assert len(commands) == 3

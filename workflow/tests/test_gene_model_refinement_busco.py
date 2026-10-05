"""Paired BUSCO scores bind counts, provenance, species coverage and cache bytes."""

import json
import os
import subprocess
from importlib import import_module
from pathlib import Path
from xml.etree import ElementTree

import pytest

from workflow.tests.test_plot_gene_model_refinement import completed

busco = import_module("gene_model_refinement_busco")
review = import_module("plot_gene_model_refinement")


def summary(complete=8, duplicated=1, version="6.1.0"):
    return f"""# BUSCO version is: {version}
# The lineage dataset is: embryophyta_odb12 (Creation date: 2026-05-22, number of genomes: 78, number of BUSCOs: 10)
# BUSCO was run in mode: euk_tran
C:{complete * 10:.1f}%[S:{(complete - duplicated) * 10:.1f}%,D:{duplicated * 10:.1f}%],F:10.0%,M:{(9-complete)*10:.1f}%,n:10
{complete} Complete BUSCOs (C)
{complete - duplicated} Complete and single-copy BUSCOs (S)
{duplicated} Complete and duplicated BUSCOs (D)
1 Fragmented BUSCOs (F)
{9 - complete} Missing BUSCOs (M)
10 Total BUSCO groups searched
Dependencies and versions:
 hmmsearch: 3.4
 metaeuk: 7.bba0d80
"""


def test_all_species_coverage_marks_cds_only_species_unassessed(tmp_path):
    root = completed(tmp_path)
    cds_dir = tmp_path / "all_species"
    cds_dir.mkdir()
    extra = cds_dir / "Drosophyllum_lusitanicum.fa"
    extra.write_text(">Dros_g1\nATGAAATAA\n")
    pairs = busco.input_pairs(root, cds_dir)
    row = next(r for r in pairs if r["species"] == "Drosophyllum_lusitanicum")
    assert len(pairs) == 4
    assert row["refinement_status"] == "not_analysed"
    assert row["before"] == row["after"] == str(extra)
    data = review.collect(root, max_loci=1, cds_dir=cds_dir)
    assert data["species"][row["species"]]["refinement_status"] == "not_analysed"
    assert "accepted_repair_paths" not in data["species"][row["species"]]
    output = tmp_path / "review"
    output.mkdir()
    review.plot_summary(data, output)
    review.write_review(data, output)
    assert "Drosophyllum lusitanicum [not analysed]" in (output / "summary.svg").read_text()
    assert "No matching genome/GFF" in (output / "review.html").read_text()
    assert "Not analysed</td>" in (output / "review.html").read_text()


def test_exact_counts_drive_deltas_and_plot_includes_excluded_species(tmp_path):
    before, after = tmp_path / "before.txt", tmp_path / "after.txt"
    before.write_text(summary())
    after.write_text(summary(9, 2))
    pair = {"species": "Drosophyllum_lusitanicum", "refinement_status": "not_analysed", "reason": "No genome/GFF"}
    row = busco.paired_result(pair, busco.read_result(before), busco.read_result(after))
    assert row["delta_complete"] == 1
    assert row["delta_complete_pp"] == 10
    assert row["delta_duplicated"] == 1
    busco.plot_comparison([row], tmp_path)
    ElementTree.parse(tmp_path / "busco_comparison.svg")
    assert "Drosophyllum lusitanicum  [not analysed]" in (tmp_path / "busco_comparison.svg").read_text()


@pytest.mark.parametrize("field,value", [("busco_version", "6.0"), ("mode", "proteins"),
                                         ("lineage_creation_date", "2025-01-01"), ("total", 100),
                                         ("dependencies", {"metaeuk": "different"})])
def test_noncomparable_pairs_are_rejected(tmp_path, field, value):
    path = tmp_path / "summary.txt"
    path.write_text(summary())
    before = busco.read_result(path)
    after = dict(before, **{field: value})
    with pytest.raises(ValueError, match="Noncomparable"):
        busco.paired_result({"species": "Species_a"}, before, after)


def test_bad_counts_or_absent_metadata_are_not_zero_scores(tmp_path):
    path = tmp_path / "summary.txt"
    path.write_text(summary().replace("10 Total", "11 Total"))
    with pytest.raises(ValueError, match="Inconsistent"):
        busco.read_result(path)
    path.write_text(summary().replace("Creation date: 2026-05-22,", ""))
    with pytest.raises(ValueError, match="comparability metadata"):
        busco.read_result(path)


def test_same_input_reuses_bound_before_and_detects_tampering(tmp_path, monkeypatch):
    source = tmp_path / "Species_a.fa"
    source.write_text(">a\nATGAAATAA\n")
    pair = {"species": "Species_a", "before": str(source), "after": str(source)}
    contract = {"lineage_path": "lineage", "download_path": "db"}
    calls = []
    def fake_run(command, **kwargs):
        calls.append(command)
        target = Path(kwargs["cwd"]) / "busco"
        target.mkdir()
        (target / "short_summary.specific.embryophyta_odb12.busco.txt").write_text(summary())
    monkeypatch.setattr(busco.subprocess, "run", fake_run)
    first = busco.run_one(pair, "before", tmp_path, contract, 2)
    second = busco.run_one(pair, "after", tmp_path, contract, 2)
    assert first["complete"] == second["complete"] == 8
    assert len(calls) == 1
    receipt = tmp_path / "runs/Species_a/after/receipt.json"
    assert json.loads(receipt.read_text())["reused_identical_before"]
    source.write_text(">a\nATGCCCTAA\n")
    with pytest.raises(ValueError, match="cache changed"):
        busco.run_one(pair, "after", tmp_path, contract, 2)


def test_busco_failure_does_not_publish_complete_receipt(tmp_path, monkeypatch):
    source = tmp_path / "Species_a.fa"
    source.write_text(">a\nATGAAATAA\n")
    def fail(*args, **kwargs):
        raise busco.subprocess.CalledProcessError(1, "busco")
    monkeypatch.setattr(busco.subprocess, "run", fail)
    with pytest.raises(busco.subprocess.CalledProcessError):
        busco.run_one({"species": "Species_a", "before": str(source)}, "before", tmp_path,
                      {"lineage_path": "lineage", "download_path": "db"}, 2)
    assert not (tmp_path / "runs/Species_a/before/receipt.json").exists()


@pytest.mark.parametrize("enabled", [0, 1])
def test_input_generation_finisher_wires_review_and_bounds_busco_resources(tmp_path, enabled):
    core = Path(__file__).resolve().parents[1] / "core/gg_input_generation_core.sh"
    function = core.read_text().split("finish_gene_model_refinement() {", 1)[1].split("\n}\n", 1)[0]
    log = tmp_path / "commands.txt"
    script = """set -euo pipefail
python() { printf '%s\\n' "$*" >> "$COMMAND_LOG"; }
ensure_shared_busco_lineage_ready() { busco_lineage_resolved=embryophyta_odb12; }
ensure_busco_download_path() { printf '%s\\n' /db; }
gg_memory_parallel_job_cap() { printf '%s\\n' 2; }
finish_gene_model_refinement() {""" + function + "\n}\nfinish_gene_model_refinement\n"
    env = dict(os.environ, COMMAND_LOG=str(log), run_species_busco=str(enabled), gg_support_dir="/support",
               gene_model_refinement_dir="/refinement", species_cds_dir="/all-cds", GG_TASK_CPUS="8",
               species_busco_parallel_jobs="auto", task_plan_output="/plan", gg_workspace_dir="/workspace",
               GG_MEM_TOOL_GB="64", species_busco_memory_gb_per_job="16", busco_lineage_resolved="")
    result = subprocess.run(["bash", "-c", script], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    commands = log.read_text().splitlines()
    assert " finalize " in commands[0] and " qc " in commands[1]
    assert "plot_gene_model_refinement.py" in commands[2] and "--cds-dir /all-cds" in commands[2]
    if enabled:
        assert "gene_model_refinement_busco.py" in commands[3]
        assert "--jobs 2 --cpus 4" in commands[3]
        assert "--report /refinement.review/busco" in commands[3]
    else:
        assert len(commands) == 3

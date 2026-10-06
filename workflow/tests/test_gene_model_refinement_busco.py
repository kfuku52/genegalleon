"""Paired BUSCO scores bind counts, provenance, species coverage and cache bytes."""

import gzip
import json
import os
import subprocess
from importlib import import_module
from pathlib import Path
from xml.etree import ElementTree

import pytest

from workflow.tests.test_gene_model_refinement import refinement, tiny_inputs
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


def test_rescue_counts_use_gene_features_not_transcripts_and_verify_sources(tmp_path):
    inputs, edges, sources = tiny_inputs(tmp_path)
    gff = Path(sources[0]["gff"])
    # One rescued gene with two transcripts and three CDS features is one locus.
    gff.write_text(gff.read_text().replace("\ts\t", "\tgenegalleon_rescue\t"))
    root = tmp_path / "refinement"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    refinement.finalize(root, value)
    pairs = busco.input_pairs(root)
    changes = busco.collect_model_changes(root, pairs)
    assert changes["species"]["Species_target"]["prior_rescued_loci"] == 1
    assert changes["species"]["Species_donor1"]["prior_rescued_loci"] == 0
    assert sum(s["accepted_repair_paths"] for s in changes["species"].values()) == 0
    assert changes["plan_sha256"] == busco.digest(root / "plan.json")
    mismatched = [dict(p) for p in pairs]
    mismatched[0]["after"] = str(tmp_path / "another.fa")
    with pytest.raises(ValueError, match="different BUSCO inputs"):
        busco.collect_model_changes(root, mismatched)
    gff.write_text(gff.read_text() + "# changed\n")
    with pytest.raises(ValueError, match="source annotation changed"):
        busco.collect_model_changes(root, pairs)


def test_combined_figure_preserves_palette_count_units_and_unavailable_species(tmp_path):
    path = tmp_path / "summary.txt"
    path.write_text(summary())
    result = busco.read_result(path)
    rows = [busco.paired_result({"species": name, "refinement_status": status}, result, result)
            for name, status in [("Species_a", "analysed"), ("Drosophyllum_lusitanicum", "not_analysed")]]
    changes = {"species": {
        "Species_a": {"refinement_status": "analysed", "prior_rescued_loci": 924,
                      "accepted_repair_paths": 9, "accepted_isoform_paths": 28},
        "Drosophyllum_lusitanicum": {"refinement_status": "not_analysed", "prior_rescued_loci": None,
                                   "accepted_repair_paths": None, "accepted_isoform_paths": None},
    }}
    busco.plot_comparison(rows, tmp_path, changes)
    svg = (tmp_path / "busco_comparison.svg").read_text()
    ElementTree.fromstring(svg)
    for color in ("#000000", "#b22222", "#666666", "#cccccc"):
        assert color in svg.lower()
    assert "Missing-gene rescue" in svg and "Accepted coding paths" in svg
    assert ">924<" in svg and ">9 / 28<" in svg
    assert svg.count(">Not analysed<") == 2
    assert "may share a locus" in svg
    changes["species"]["Drosophyllum_lusitanicum"]["prior_rescued_loci"] = 0
    with pytest.raises(ValueError, match="unavailable, not zero"):
        busco.plot_comparison(rows, tmp_path, changes)
    changes["species"].pop("Drosophyllum_lusitanicum")
    with pytest.raises(ValueError, match="membership differ"):
        busco.plot_comparison(rows, tmp_path, changes)


def test_rescue_support_classification_uses_all_support_and_frozen_overlap():
    name = "Species_target"
    plan = {"nearest_references": {name: ["Near", "Overlap"]},
            "common_references": ["Balanced", "Overlap", name],
            "donors": {name: ["Near", "Balanced", "Overlap"]},
            "synteny_jobs": [{"id": "self", "a": name, "b": name, "kind": "self"}]}
    def model(identifier, donors, status="accepted"):
        return {"model_id": identifier, "status": status,
                "evidence": {"donor": donors[0]},
                "support": [{"donor": d, "target": name} for d in donors]}
    models = [model("n", ["Near", "Near"]), model("c", ["Balanced"]),
              model("b", ["Near", "Balanced"]), model("o", ["Overlap"]),
              model("duplicate", ["Near"], "duplicate_support")]
    models += [{"model_id": "s", "status": "accepted", "support": [{"donor": name, "target": name, "comparison": "self"}]},
               {"model_id": "mixed", "status": "accepted", "support": [{"donor": name, "target": name, "comparison": "self"},
                                                                         {"donor": "Near", "target": name}]}]
    counts, evidence = busco.classify_rescue_support(iter(models), name, {"n", "c", "b", "o", "s", "mixed"}, plan)
    assert counts == {"nearest_only": 2, "balanced_only": 1, "both": 2}
    assert evidence["s"]["category"] == "self_only"
    assert evidence["mixed"]["category"] == "nearest_only"
    assert evidence["n"]["supporting_donors"] == ["Near"]
    assert evidence["o"]["category"] == "both"
    mapped, evidence = busco.classify_rescue_support([model("tx", ["Near"])], name, {"gene": {"gene", "tx"}}, plan)
    assert mapped["nearest_only"] == 1 and evidence["gene"]["source_model_id"] == "tx"
    with pytest.raises(ValueError, match="gene IDs differ"):
        busco.classify_rescue_support(models, name, {"n"}, plan)
    with pytest.raises(ValueError, match="Duplicate accepted"):
        busco.classify_rescue_support([models[0], models[0]], name, {"n"}, plan)
    with pytest.raises(ValueError, match="lacks supporting"):
        busco.classify_rescue_support([dict(models[0], support=[])], name, {"n"}, plan)
    with pytest.raises(ValueError, match="donor/target differs"):
        busco.classify_rescue_support([model("x", ["Unknown"])], name, {"x"}, plan)
    with pytest.raises(ValueError, match="donor/target differs"):
        busco.classify_rescue_support([model("x", [name])], name, {"x"}, plan)


def test_imported_rescue_support_verifies_receipts_models_and_gff(tmp_path):
    inputs, edges, sources = tiny_inputs(tmp_path)
    gff = Path(sources[0]["gff"])
    gff.write_text(gff.read_text().replace("\ts\t", "\tgenegalleon_rescue\t"))
    root = tmp_path / "refinement"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    refinement.finalize(root, value)
    rescue = tmp_path / "original_rescue"
    rescue.mkdir()
    names = [r["species"] for r in sources]
    busco.atomic_json(rescue / "plan.json", {
        "nearest_references": {n: ["Species_donor1"] for n in names},
        "common_references": ["Species_donor2"],
        "donors": {n: sorted({"Species_donor1", "Species_donor2"} - {n}) for n in names},
    })
    plan_hash = busco.digest(rescue / "plan.json")
    receipts, files = {}, {}
    for n, source in zip(names, sources, strict=True):
        models = [{"status": "accepted", "model_id": "t1", "support": [
            {"donor": "Species_donor1", "target": n}, {"donor": "Species_donor2", "target": n}]}] if n == names[0] else []
        worker = rescue / "rescued" / n
        busco.atomic_json(worker / "models.json", models)
        busco.atomic_json(worker / "receipt.json", {"key": {"plan": plan_hash, "species": n},
                                                   "files": {"models.json": busco.digest(worker / "models.json")}})
        receipts[n] = busco.digest(worker / "receipt.json")
        files["species_gff/" + n + ".rescue.gff3"] = busco.digest(source["gff"])
    augmented = {"key": {"plan": plan_hash, "rescue_receipts": receipts}, "files": files}
    busco.atomic_json(rescue / "augmented/receipt.json", augmented)
    pairs = busco.input_pairs(root)
    changes = busco.collect_model_changes(root, pairs, rescue)
    assert changes["species"][names[0]]["rescue_support_counts"] == {"nearest_only": 0, "balanced_only": 0, "both": 1}
    assert changes["evidence"][names[0]]["rescued_loci_support"]["g"]["supporting_donors"] == names[1:]
    assert changes["evidence"][names[0]]["rescued_loci_support"]["g"]["source_model_id"] == "t1"
    assert changes["rescue_reference_selection"]["plan_sha256"] == plan_hash
    altered = rescue / "rescued" / names[0] / "models.json"
    original = altered.read_text()
    altered.write_text(original + "\n")
    with pytest.raises(ValueError, match="models changed"):
        busco.collect_model_changes(root, pairs, rescue)
    altered.write_text(original)
    augmented["files"]["species_gff/" + names[0] + ".rescue.gff3"] = "wrong"
    busco.atomic_json(rescue / "augmented/receipt.json", augmented)
    with pytest.raises(ValueError, match="differs from the source annotation"):
        busco.collect_model_changes(root, pairs, rescue)


def test_rescue_stacks_render_all_three_groups_and_reject_missing_or_wrong_totals(tmp_path):
    path = tmp_path / "summary.txt"
    path.write_text(summary())
    result = busco.read_result(path)
    rows = [busco.paired_result({"species": n, "refinement_status": status}, result, result)
            for n, status in [("Species_a", "analysed"), ("Drosophyllum_lusitanicum", "not_analysed")]]
    changes = {"species": {
        "Species_a": {"refinement_status": "analysed", "prior_rescued_loci": 11, "accepted_repair_paths": 2,
                      "accepted_isoform_paths": 3, "rescue_self_only_loci": 1,
                      "rescue_support_counts": {"nearest_only": 3, "balanced_only": 2, "both": 5}},
        "Drosophyllum_lusitanicum": {"refinement_status": "not_analysed", "prior_rescued_loci": None,
                                   "accepted_repair_paths": None, "accepted_isoform_paths": None},
    }}
    busco.plot_comparison(rows, tmp_path, changes)
    svg = (tmp_path / "busco_comparison.svg").read_text()
    for color in busco.RESCUE_SUPPORT_COLOURS:
        assert color in svg
    for label in busco.RESCUE_SUPPORT_LABELS:
        assert label in svg
    assert "donor belonging to both lists" in svg
    assert ">10 / 1<" in svg and "self-only loci are outside the three groups" in svg
    assert svg.count(">Not analysed<") == 2
    changes["species"]["Species_a"]["rescue_support_counts"]["both"] = 6
    with pytest.raises(ValueError, match="sum to the rescued"):
        busco.plot_comparison(rows, tmp_path, changes)


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


@pytest.mark.parametrize("compressed", [False, True])
def test_same_input_reuses_bound_before_and_detects_tampering(tmp_path, monkeypatch, compressed):
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
        full = target / "run_lineage/full_table.tsv"
        full.parent.mkdir()
        payload = b"# BUSCO id\tStatus\tSequence\nBUSCO_1\tComplete\ta\n"
        if compressed:
            with gzip.open(str(full) + ".gz", "wb") as handle:
                handle.write(payload)
        else:
            full.write_bytes(payload)
    monkeypatch.setattr(busco.subprocess, "run", fake_run)
    first = busco.run_one(pair, "before", tmp_path, contract, 2)
    second = busco.run_one(pair, "after", tmp_path, contract, 2)
    assert first["complete"] == second["complete"] == 8
    assert len(calls) == 1
    receipt = tmp_path / "runs/Species_a/after/receipt.json"
    assert json.loads(receipt.read_text())["reused_identical_before"]
    assert (receipt.parent / "full_table.tsv").read_bytes() == (receipt.parent.parent / "before/full_table.tsv").read_bytes()
    full_table = receipt.parent / "full_table.tsv"
    original = full_table.read_bytes()
    full_table.write_text("tampered\n")
    with pytest.raises(ValueError, match="cache changed"):
        busco.run_one(pair, "after", tmp_path, contract, 2)
    full_table.write_bytes(original)
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


def test_plot_only_verifies_historical_scores_and_rejects_tampering(tmp_path, monkeypatch):
    contract = {"historic_evaluation": "fixture"}
    pair = {"species": "Species_a", "refinement_status": "analysed", "reason": ""}
    for phase, count in (("before", 8), ("after", 9)):
        source = tmp_path / (phase + ".fa")
        source.write_text(">a\nATG" + ("AAA" if phase == "before" else "CCC") + "TAA\n")
        pair[phase] = str(source)
        directory = tmp_path / "runs/Species_a" / phase
        directory.mkdir(parents=True)
        (directory / "summary.txt").write_text(summary(count))
        busco.atomic_json(directory / "receipt.json", {
            "key": {"contract": contract, "source_sha256": busco.digest(source)},
            "summary_sha256": busco.digest(directory / "summary.txt"),
        })
    row = busco.paired_result(pair, busco.read_result(tmp_path / "runs/Species_a/before/summary.txt"),
                              busco.read_result(tmp_path / "runs/Species_a/after/summary.txt"))
    busco.atomic_json(tmp_path / "contract.json", {"contract": contract, "pairs": [pair]})
    value = {"contract": contract, "species": [row]}
    busco.atomic_json(tmp_path / "busco_comparison.json", value)
    def fail(*args, **kwargs):
        raise AssertionError("Predictor must not execute while redrawing")
    monkeypatch.setattr(busco, "run_one", fail)
    assert busco.render_existing(tmp_path)[0]["delta_complete"] == 1
    provenance = json.loads((tmp_path / "rendering_provenance.json").read_text())
    assert not provenance["predictor_executed"] and provenance["verified_full_tables"] == 0
    row["delta_complete"] = 99
    busco.atomic_json(tmp_path / "busco_comparison.json", value)
    with pytest.raises(ValueError, match="delta changed"):
        busco.render_existing(tmp_path)
    row["delta_complete"] = 1
    busco.atomic_json(tmp_path / "busco_comparison.json", value)
    Path(pair["after"]).write_text(">a\nATGTTTTAA\n")
    with pytest.raises(ValueError, match="input or score changed"):
        busco.render_existing(tmp_path)


@pytest.mark.parametrize("enabled", [0, 1])
@pytest.mark.parametrize("trailing_slash", [False, True])
def test_input_generation_finisher_wires_review_and_bounds_busco_resources(tmp_path, enabled, trailing_slash):
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
               gene_model_refinement_dir="/refinement/" if trailing_slash else "/refinement", species_cds_dir="/all-cds", GG_TASK_CPUS="8",
               species_busco_parallel_jobs="auto", task_plan_output="/plan", gg_workspace_dir="/workspace",
               GG_MEM_TOOL_GB="64", species_busco_memory_gb_per_job="16", busco_lineage_resolved="")
    result = subprocess.run(["bash", "-c", script], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    commands = log.read_text().splitlines()
    assert " finalize " in commands[0] and " qc " in commands[1]
    assert "plot_gene_model_refinement.py" in commands[2] and "--cds-dir /all-cds" in commands[2]
    assert "--report /refinement.review " in commands[2]
    if enabled:
        assert "gene_model_refinement_busco.py" in commands[3]
        assert "--jobs 2 --cpus 4" in commands[3]
        assert "--report /refinement.review/busco" in commands[3]
    else:
        assert len(commands) == 3

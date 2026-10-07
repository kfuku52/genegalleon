"""Paired BUSCO scores bind counts, provenance, species coverage and cache bytes."""

import gzip
import json
import os
import shutil
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


@pytest.mark.parametrize("species_count", [1, 24])
def test_three_stage_busco_separates_rescue_and_refinement_and_preserves_palette(tmp_path, monkeypatch, species_count):
    from matplotlib.figure import Figure
    original_save = Figure.savefig
    def save(fig, *args, **kwargs):
        fig.canvas.draw()
        assert len(fig.axes) == 2
        # Each species has exactly three four-category stacks. Check phase
        # order and percentages, rather than only the presence of SVG colours.
        patches = fig.axes[0].patches[:12 * species_count]
        assert len(patches) == 12 * species_count
        for phase_index, (phase, offset) in enumerate((("pre_rescue", -.26), ("before", 0), ("after", .26))):
            for i, row in enumerate(rows):
                stack = patches[(phase_index * species_count + i) * 4:(phase_index * species_count + i + 1) * 4]
                assert sum(p.get_width() for p in stack) == pytest.approx(100)
                for patch, category in zip(stack, busco.STATUS, strict=True):
                    assert patch.get_height() == pytest.approx(.18)
                    assert patch.get_y() + patch.get_height() / 2 == pytest.approx(i + offset)
                    assert patch.get_width() == pytest.approx(100 * row[phase + "_result"][category] / row[phase + "_result"]["total"])
        renderer = fig.canvas.get_renderer()
        assert all(ax.get_legend() is None for ax in fig.axes)
        assert len(fig.legends) == 2
        for legend in fig.legends:
            bounds = legend.get_window_extent(renderer)
            assert all(not bounds.overlaps(ax.bbox) for ax in fig.axes)
            assert all(not bounds.overlaps(ax.xaxis.label.get_window_extent(renderer)) for ax in fig.axes)
            assert not bounds.overlaps(fig.texts[-1].get_window_extent(renderer))
            start, end = (0, 0) if "Single-copy" in [t.get_text() for t in legend.get_texts()] else (1, 1)
            assert fig.axes[start].bbox.x0 - 1 <= bounds.x0 < bounds.x1 <= fig.axes[end].bbox.x1 + 1
        return original_save(fig, *args, **kwargs)
    monkeypatch.setattr(Figure, "savefig", save)
    summaries = []
    for i, complete in enumerate((7, 9, 8)):
        path = tmp_path / f"s{i}.txt"
        path.write_text(summary(complete))
        summaries.append(busco.read_result(path))
    pair = {"species": "Drosophyllum_lusitanicum", "pre_rescue": "pre.fa",
            "before": "rescued.fa", "after": "refined.fa",
            "refinement_status": "not_analysed", "reason": "No genome"}
    row = busco.staged_result(pair, summaries[1], summaries[2], summaries[0])
    assert (row["delta_rescue_complete"], row["delta_complete"], row["delta_total_complete"]) == (2, -1, 1)
    rows = [dict(row, species=f"Species_{i}", refinement_status="analysed") for i in range(species_count - 1)] + [row]
    busco.plot_three_stage(rows, tmp_path)
    svg = (tmp_path / "busco_three_stage.svg").read_text()
    assert "Drosophyllum lusitanicum (not analysed)" in svg
    assert all(c.lower() in svg.lower() for c in busco.STATUS_COLOURS)
    assert "Top to bottom within each species" in svg
    summaries[0]["lineage"] = "insecta_odb12"
    with pytest.raises(ValueError, match="Noncomparable"):
        busco.staged_result(pair, summaries[1], summaries[2], summaries[0])


def test_swissprot_diagnostic_plot_records_thresholds_and_excluded_species(tmp_path, monkeypatch):
    import copy

    from matplotlib.figure import Figure
    original_save = Figure.savefig
    def save(fig, *args, **kwargs):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        assert len(fig.legends) == 2
        for ax, legend in zip(fig.axes, fig.legends, strict=True):
            bounds = legend.get_window_extent(renderer)
            assert ax.bbox.x0 - 1 <= bounds.x0 < bounds.x1 <= ax.bbox.x1 + 1
            assert all(not bounds.overlaps(a.bbox) for a in fig.axes)
            assert not bounds.overlaps(ax.xaxis.label.get_window_extent(renderer))
            assert not bounds.overlaps(fig.texts[-1].get_window_extent(renderer))
        return original_save(fig, *args, **kwargs)
    monkeypatch.setattr(Figure, "savefig", save)
    path = tmp_path / "summary.txt"
    path.write_text(summary())
    result = busco.read_result(path)
    rows = [busco.paired_result({"species": "Species_a", "refinement_status": "analysed", "reason": ""}, result, result),
            busco.paired_result({"species": "Drosophyllum_lusitanicum", "refinement_status": "not_analysed", "reason": ""}, result, result)]
    from rescue_swissprot_evidence import DEFAULTS, NO_SUPPORT_REASONS
    changes = {"swissprot_evidence": {"parameters": DEFAULTS}, "species": {
        "Species_a": {"refinement_status": "analysed", "prior_rescued_loci": 6,
                      "rescue_swissprot_groups": dict(zip(busco.SWISSPROT_GROUPS, (1, 0, 0, 5, 0), strict=True)),
                      "rescue_partial_te_groups": dict(zip(("primary_te_support", "partial_te_only", "no_te_support", "not_assessed"),
                                                         (1, 2, 3, 0), strict=True)),
                      "rescue_no_support_reasons": dict.fromkeys(NO_SUPPORT_REASONS, 1)},
        "Drosophyllum_lusitanicum": {"refinement_status": "not_analysed", "prior_rescued_loci": None}}}
    busco.plot_swissprot_diagnostics(rows, tmp_path, changes)
    svg = (tmp_path / "rescue_swissprot_diagnostics.svg").read_text()
    assert "Partial TE flag only" in svg
    assert "without a competing-score filter" in svg
    assert "short-protein thresholds are unchanged" in svg
    assert "Not analysed" in svg
    for field, replacement in (
        ("rescue_partial_te_groups", None),
        ("rescue_no_support_reasons", {}),
        ("rescue_swissprot_groups", dict(zip(busco.SWISSPROT_GROUPS, (0, 1, 0, 5, 0), strict=True))),
        ("prior_rescued_loci", 7),
    ):
        broken = copy.deepcopy(changes)
        broken["species"]["Species_a"][field] = replacement
        with pytest.raises(ValueError, match="Swiss-Prot diagnostic"):
            busco.plot_swissprot_diagnostics(rows, tmp_path, broken)
    broken = copy.deepcopy(changes)
    broken["species"]["Species_a"]["rescue_no_support_reasons"]["no_returned_hits"] = -1
    with pytest.raises(ValueError, match="invalid Swiss-Prot diagnostic"):
        busco.plot_swissprot_diagnostics(rows, tmp_path, broken)
    broken = copy.deepcopy(changes)
    broken["species"]["Drosophyllum_lusitanicum"]["rescue_no_support_reasons"] = dict.fromkeys(NO_SUPPORT_REASONS, 0)
    with pytest.raises(ValueError, match="Unanalysed Swiss-Prot diagnostics"):
        busco.plot_swissprot_diagnostics(rows, tmp_path, broken)


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


@pytest.mark.parametrize(('record', 'category'), [
    ({'status': 'available', 'masked_fraction': .8, 'te_annotated_fraction': .1, 'classes': ['LTR/Gypsy', 'Simple_repeat']}, 'te_hit'),
    ({'status': 'available', 'masked_fraction': .2, 'te_annotated_fraction': 0, 'classes': ['Unknown']}, 'other_repeat_hit'),
    ({'status': 'available', 'masked_fraction': 0, 'te_annotated_fraction': 0, 'classes': []}, 'no_repeat_hit'),
    ({'status': 'not_provided', 'masked_fraction': None, 'te_annotated_fraction': None, 'classes': None}, 'not_assessed'),
])
def test_repeat_categories_keep_unassessed_separate_from_negative(record, category):
    assert busco.repeat_group(record) == category
    with pytest.raises(ValueError, match='repeat'):
        busco.repeat_group(dict(record, masked_fraction=-.1))


def test_repeat_audit_join_binds_models_counts_and_snapshot(tmp_path, monkeypatch):
    name = 'Species_a'
    directory = tmp_path / name
    source = {'rescue_models_sha256': 'models_hash', 'rescue_receipt_sha256': 'receipt_hash',
              'rescued_loci_support': {'gene': {'source_model_id': 'tx'}}}
    changes = {'rescue_reference_selection': {'plan_sha256': 'plan_hash'}, 'evidence': {name: source},
               'species': {name: {'refinement_status': 'analysed', 'prior_rescued_loci': 1},
                           'Drosophyllum_lusitanicum': {'refinement_status': 'not_analysed', 'prior_rescued_loci': None}}}
    busco.collect_rescue_repeat_evidence(changes)
    assert changes['species'][name]['rescue_repeat_groups']['not_assessed'] == 1
    assert changes['species']['Drosophyllum_lusitanicum']['rescue_repeat_groups'] is None
    records = [{'model_id': 'tx', 'rescue_status': 'accepted', 'repeat': {
        'status': 'available', 'masked_fraction': .7, 'te_annotated_fraction': .6, 'classes': ['LTR/Gypsy']}}]
    busco.atomic_json(directory / 'evidence.json', records)
    receipt = {'key': {'schema': 1, 'species': name, 'inputs': {
        '/rescue/plan.json': 'plan_hash', '/rescue/rescued/Species_a/models.json': 'models_hash',
        '/rescue/rescued/Species_a/receipt.json': 'receipt_hash'}},
        'files': {'evidence.json': busco.digest(directory / 'evidence.json')}}
    busco.atomic_json(directory / 'receipt.json', receipt)
    busco.collect_rescue_repeat_evidence(changes, tmp_path)
    assert changes['species'][name]['rescue_repeat_groups']['te_hit'] == 1
    assert changes['evidence'][name]['rescued_loci_repeat']['loci']['gene']['te_annotated_fraction'] == .6
    receipt['key']['inputs']['/rescue/rescued/Species_a/models.json'] = 'other_models'
    busco.atomic_json(directory / 'receipt.json', receipt)
    with pytest.raises(ValueError, match='different rescue models'):
        busco.collect_rescue_repeat_evidence(changes, tmp_path)
    receipt['key']['inputs']['/rescue/rescued/Species_a/models.json'] = 'models_hash'
    records[0]['model_id'] = 'other_tx'
    busco.atomic_json(directory / 'evidence.json', records)
    receipt['files']['evidence.json'] = busco.digest(directory / 'evidence.json')
    busco.atomic_json(directory / 'receipt.json', receipt)
    with pytest.raises(ValueError, match='membership differs'):
        busco.collect_rescue_repeat_evidence(changes, tmp_path)
    records[0]['model_id'] = 'tx'
    busco.atomic_json(directory / 'evidence.json', records)
    receipt['files']['evidence.json'] = busco.digest(directory / 'evidence.json')
    busco.atomic_json(directory / 'receipt.json', receipt)
    import rescue_model_evidence
    read = rescue_model_evidence.read_json_snapshot
    def changed(path):
        result = read(path)
        if path.name == 'evidence.json':
            path.write_text('[]')
        return result
    monkeypatch.setattr(rescue_model_evidence, 'read_json_snapshot', changed)
    with pytest.raises(ValueError, match='changed while loading'):
        busco.collect_rescue_repeat_evidence(changes, tmp_path)


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


def test_accepted_path_support_uses_verified_donors_not_target_rna():
    selection = {"nearest_references": {"Target": ["Near", "Overlap"]}, "common_references": ["Balanced", "Overlap"]}
    allowed = {"Target", "Near", "Balanced", "Overlap", "Other"}
    def model(identifier, donors):
        return {"status": "accepted", "change_type": "isoform_addition", "gene_id": "same_gene", "donors": donors,
                "candidate": {"candidate_id": identifier, "source_transcript_id": identifier, "support": {"donors": donors, "rna_paths": ["target_RNA"]}},
                "alignments": [{"donor_species": d} for d in donors]}
    models = [model("n", ["Near", "Near", "Other"]), model("b", ["Balanced"]),
              model("both", ["Near", "Balanced"]), model("overlap", ["Overlap"]),
              model("s", ["Target"]), model("mixed", ["Target", "Near"]), model("other", ["Target", "Other"])]
    counts, evidence = busco.classify_accepted_path_support(models, "Target", selection, allowed)
    assert counts == {"nearest_only": 2, "balanced_only": 1, "both": 2, "self_only": 1, "other_interspecies": 1}
    assert len(evidence) == 7 and evidence["n"]["other_supporting_donors"] == ["Other"]
    assert evidence["b"]["category"] == "balanced_only"  # Target RNA does not count as self homology.
    assert busco.group_support_evidence(evidence, "Target", selection, allow_other=True) == {
        "self_only": 2, "relative_only": 1, "phylogenetic_only": 1, "multiple": 3, "other_interspecies": 0,
    }
    for broken in [model("empty", []), model("unknown", ["Unknown"])]:
        with pytest.raises(ValueError, match="valid supporting donors"):
            busco.classify_accepted_path_support([broken], "Target", selection, allowed)
    with pytest.raises(ValueError, match="records disagree"):
        busco.classify_accepted_path_support([models[0], models[0]], "Target", selection, allowed)
    broken = model("bad", ["Near"])
    broken["alignments"] = [{"donor_species": "Balanced"}]
    with pytest.raises(ValueError, match="records disagree"):
        busco.classify_accepted_path_support([broken], "Target", selection, allowed)


def test_saved_model_summary_adds_path_support_and_rejects_changed_receipt(tmp_path):
    inputs, edges, sources = tiny_inputs(tmp_path)
    root = tmp_path / "refinement"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    refinement.finalize(root, value)
    changes = busco.collect_model_changes(root, busco.input_pairs(root))
    names = list(changes["species"])
    changes["rescue_reference_selection"] = {"nearest_references": {n: [] for n in names}, "common_references": names}
    busco.collect_accepted_path_support(root, changes)
    assert all(s["accepted_path_support_counts"] == dict.fromkeys(busco.PATH_SUPPORT, 0) for s in changes["species"].values())
    assert all(e["accepted_paths_support"] == {} for e in changes["evidence"].values())
    assert all(s["accepted_path_support_groups"] == dict.fromkeys(busco.PATH_SUPPORT_GROUPS, 0) for s in changes["species"].values())
    receipt = root / "predictions" / names[0] / "receipt.json"
    receipt.write_text(receipt.read_text() + "\n")
    with pytest.raises(ValueError, match="Saved prediction receipt changed"):
        busco.collect_accepted_path_support(root, changes)


def test_srp_regrouping_preserves_legacy_counts_and_requires_complete_evidence():
    import copy
    selection = {"nearest_references": {"Target": ["Target", "Near", "Near2", "Overlap"]},
                 "common_references": ["Target", "Balanced", "Overlap"]}
    evidence = {key: {"supporting_donors": donors} for key, donors in {
        "s": ["Target"], "r": ["Near", "Near2"], "p": ["Balanced"],
        "sr": ["Target", "Near"], "sp": ["Target", "Balanced"], "rp": ["Overlap"],
        "srp": ["Target", "Near", "Balanced"],
    }.items()}
    paths = {**evidence, "other": {"supporting_donors": ["Other"]}}
    changes = {"rescue_reference_selection": selection, "species": {
        "Target": {"refinement_status": "analysed", "prior_rescued_loci": 7,
                   "accepted_repair_paths": 2, "accepted_isoform_paths": 6,
                   "rescue_support_counts": {"nearest_only": 2, "balanced_only": 2, "both": 2},
                   "rescue_self_only_loci": 1,
                   "accepted_path_support_counts": {"nearest_only": 2, "balanced_only": 2, "both": 2,
                                                    "self_only": 1, "other_interspecies": 1}},
        "Excluded": {"refinement_status": "not_analysed", "prior_rescued_loci": None,
                     "accepted_repair_paths": None, "accepted_isoform_paths": None},
    }, "evidence": {"Target": {"rescued_loci_support": evidence, "accepted_paths_support": paths}}}
    original = copy.deepcopy(changes)
    busco.regroup_model_support(changes)
    for key, value in original["species"]["Target"].items():
        assert changes["species"]["Target"][key] == value
    assert changes["evidence"] == original["evidence"]
    assert changes["species"]["Target"]["rescue_support_groups"] == {
        "self_only": 1, "relative_only": 1, "phylogenetic_only": 1, "multiple": 4,
    }
    assert changes["species"]["Target"]["accepted_path_support_groups"] == {
        "self_only": 1, "relative_only": 1, "phylogenetic_only": 1, "multiple": 4, "other_interspecies": 1,
    }
    assert changes["species"]["Excluded"] == original["species"]["Excluded"]
    broken = copy.deepcopy(original)
    del broken["evidence"]["Target"]["accepted_paths_support"]["other"]
    with pytest.raises(ValueError, match="complete per-model evidence"):
        busco.regroup_model_support(broken)
    assert "rescue_support_groups" not in broken["species"]["Target"]  # No partial updates.
    changes["species"]["Target"]["rescue_support_groups"]["multiple"] -= 1
    with pytest.raises(ValueError, match="disagree with per-model evidence"):
        busco.regroup_model_support(changes)
    for donors in [[], [None], [""]]:
        with pytest.raises(ValueError, match="Missing supporting donor evidence"):
            busco.group_support_evidence({"bad": {"supporting_donors": donors}}, "Target", selection)
    with pytest.raises(ValueError, match="differs from the frozen"):
        busco.group_support_evidence({"bad": {"supporting_donors": ["Other"]}}, "Target", selection)


@pytest.mark.parametrize("grouped", [False, True])
@pytest.mark.parametrize("swissprot", [False, True])
@pytest.mark.parametrize("staged", [False, True])
def test_rescue_and_two_path_stacks_include_all_support_and_reject_wrong_totals(tmp_path, monkeypatch, grouped, swissprot, staged):
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    original_barh = Axes.barh
    original_save = Figure.savefig
    bars = []
    def barh(ax, y, width, *args, **kwargs):
        bars.append((ax, width, kwargs.get("left", 0), kwargs.get("color"), y))
        return original_barh(ax, y, width, *args, **kwargs)
    monkeypatch.setattr(Axes, "barh", barh)
    def save(fig, *args, **kwargs):
        fig.canvas.draw()
        if len(fig.axes) == 2:  # Separate three-stage figure has its own geometry check.
            return original_save(fig, *args, **kwargs)
        assert len(fig.axes) == (4 if staged else 5)
        renderer = fig.canvas.get_renderer()
        legends = [legend.get_window_extent(renderer) for legend in fig.legends]
        if staged:
            assert all(not ax.title.get_window_extent(renderer).overlaps(text.get_window_extent(renderer))
                       for ax in fig.axes for text in fig.texts[:2])
        labels = [ax.xaxis.label.get_window_extent(renderer) for ax in fig.axes]
        note = fig.texts[-1].get_window_extent(renderer)
        assert all(not a.overlaps(b) for a in legends for b in labels)
        assert all(not a.overlaps(note) for a in legends)
        assert all(not a.overlaps(b) for i, a in enumerate(legends) for b in legends[i + 1:])
        for legend, bounds in zip(fig.legends, legends, strict=True):
            title = legend.get_title().get_text().replace("\n", " ")
            if title == "BUSCO status":
                left, right = fig.axes[0].bbox.x0, fig.axes[0 if staged else 1].bbox.x1
            else:
                ax = fig.axes[-1] if "coding-path" in title or "Coding-path" in title else fig.axes[-2]
                left, right = ax.bbox.x0, ax.bbox.x1
            assert left - 1 <= bounds.x0 < bounds.x1 <= right + 1
        return original_save(fig, *args, **kwargs)
    monkeypatch.setattr(Figure, "savefig", save)
    path = tmp_path / "summary.txt"
    path.write_text(summary())
    result = busco.read_result(path)
    rows = [busco.paired_result({"species": n, "refinement_status": status}, result, result)
            for n, status in [("Species_a", "analysed"), ("Drosophyllum_lusitanicum", "not_analysed")]]
    if staged:
        rows = [busco.staged_result({k: v for k, v in row.items() if not k.endswith("_result")}, result, result, result)
                for row in rows]
    changes = {"species": {
        "Species_a": {"refinement_status": "analysed", "prior_rescued_loci": 11, "accepted_repair_paths": 2,
                      "accepted_isoform_paths": 3, "rescue_self_only_loci": 1,
                      "rescue_repeat_groups": dict(zip(busco.REPEAT_GROUPS, [2, 3, 4, 2], strict=True)),
                      "rescue_support_counts": {"nearest_only": 3, "balanced_only": 2, "both": 5},
                      "accepted_path_support_counts": dict.fromkeys(busco.PATH_SUPPORT, 1)},
        "Drosophyllum_lusitanicum": {"refinement_status": "not_analysed", "prior_rescued_loci": None,
                                   "accepted_repair_paths": None, "accepted_isoform_paths": None},
    }}
    if grouped:
        changes["species"]["Species_a"].update(
            rescue_support_groups={"self_only": 1, "relative_only": 2, "phylogenetic_only": 1, "multiple": 7},
            accepted_path_support_groups=dict.fromkeys(busco.PATH_SUPPORT_GROUPS, 1),
        )
    if swissprot:
        changes["species"]["Species_a"]["rescue_swissprot_groups"] = dict(zip(busco.SWISSPROT_GROUPS, [2, 3, 1, 4, 1], strict=True))
        changes["swissprot_evidence"] = {"parameters": {
            "evalue": 1e-7, "query_coverage": .6, "target_coverage": .65, "minimum_alignment": 55,
            "score_fraction": .95, "max_hits": 40, "sensitivity": 7,
        }}
    busco.plot_comparison(rows, tmp_path, changes)
    svg = (tmp_path / "busco_comparison.svg").read_text()
    from xml.etree import ElementTree
    displayed = " ".join(" ".join(ElementTree.fromstring(svg).itertext()).split())
    for color in (*busco.RESCUE_SUPPORT_COLOURS, busco.RESCUE_SELF_COLOUR):
        assert color in svg
    for label in busco.SUPPORT_GROUP_LABELS if grouped else (*busco.RESCUE_SUPPORT_LABELS, busco.RESCUE_SELF_LABEL):
        assert label in displayed
    assert ("at least two support types" if grouped else "donor belonging to both lists") in svg
    if grouped:
        assert "S/R/P" not in svg
    assert ">11<" in svg
    assert ("Target RNA is separate from self-species homology" if grouped else "mixed self/interspecies support uses the interspecies group") in svg
    self_bar = next(b for b in bars if b[3] == busco.RESCUE_SELF_COLOUR)
    assert self_bar[1:3] == (1, 0 if grouped else 10)
    rescue_bars = [b for b in bars if b[0] is self_bar[0]]
    rescue_upper = [b for b in rescue_bars if b[4] == -.20]
    rescue_lower = [b for b in rescue_bars if b[4] == .20]
    assert len(rescue_upper) == 4 and sum(b[1] for b in rescue_upper) == 11
    assert len(rescue_lower) == (5 if swissprot else 4) and sum(b[1] for b in rescue_lower) == 11
    assert [b[1] for b in rescue_lower] == ([2, 3, 1, 4, 1] if swissprot else [2, 3, 4, 2])
    assert [b[2] for b in rescue_lower] == ([0, 2, 5, 6, 10] if swissprot else [0, 2, 5, 9])
    for label in busco.SWISSPROT_LABELS if swissprot else busco.REPEAT_GROUP_LABELS:
        assert label in displayed
    assert ("Upper: donors; lower: Swiss-Prot" if swissprot else "Upper: support; lower: repeats") in svg
    assert ("no support does not exclude TE origin" if swissprot else "no hit does not establish a true gene") in svg
    if swissprot:
        for description in ("E-value &lt;= 1e-07", "paired residues &gt;= 55 aa", "query coverage &gt;= 60%", "target coverage &gt;= 65%",
                            "Bit score &gt;= 95% of the best qualifying hit", "sensitivity 7; max hits 40", "no sequence-identity cutoff",
                            "Transposable element keyword", "TE silencing/regulation", "each locus counts once"):
            assert description in svg
    coding_ax = next(b[0] for b in bars if b[3] == "#187d97")
    upper = [b for b in bars if b[0] is coding_ax and b[4] == -.20]
    lower = [b for b in bars if b[0] is coding_ax and b[4] == .20]
    assert len(upper) == 2 and sum(b[1] for b in upper) == 5
    assert len(lower) == 5 and sum(b[1] for b in lower) == 5
    assert [b[2] for b in lower] == [0, 1, 2, 3, 4]
    assert [b[3] for b in lower] == list(busco.PATH_SUPPORT_GROUP_COLOURS if grouped else busco.PATH_SUPPORT_COLOURS)
    assert "Upper: repair / isoform; lower: support" in svg and "Other interspecies only" in svg
    assert svg.count(">Not analysed<") == 2
    changes["species"]["Species_a"]["rescue_repeat_groups"]["not_assessed"] = 3
    with pytest.raises(ValueError, match="Repeat groups must sum"):
        busco.plot_comparison(rows, tmp_path, changes)
    changes["species"]["Species_a"]["rescue_repeat_groups"]["not_assessed"] = 2
    changes["species"]["Species_a"]["rescue_support_counts"]["both"] = 6
    with pytest.raises(ValueError, match="sum to the rescued"):
        busco.plot_comparison(rows, tmp_path, changes)
    changes["species"]["Species_a"]["rescue_support_counts"]["both"] = 5
    changes["species"]["Species_a"]["accepted_path_support_counts"]["self_only"] = 2
    with pytest.raises(ValueError, match="sum to the accepted"):
        busco.plot_comparison(rows, tmp_path, changes)
    changes["species"]["Drosophyllum_lusitanicum"]["accepted_path_support_counts"] = dict.fromkeys(busco.PATH_SUPPORT, 0)
    changes["species"]["Species_a"]["accepted_path_support_counts"]["self_only"] = 1
    with pytest.raises(ValueError, match="coding-path support counts must be unavailable"):
        busco.plot_comparison(rows, tmp_path, changes)
    if grouped:
        del changes["species"]["Drosophyllum_lusitanicum"]["accepted_path_support_counts"]
        changes["species"]["Species_a"]["rescue_support_groups"]["multiple"] += 1
        with pytest.raises(ValueError, match="S/R/P support groups must sum"):
            busco.plot_comparison(rows, tmp_path, changes)
        changes["species"]["Species_a"]["rescue_support_groups"]["multiple"] -= 1
        changes["species"]["Drosophyllum_lusitanicum"]["accepted_path_support_groups"] = dict.fromkeys(busco.PATH_SUPPORT_GROUPS, 0)
        with pytest.raises(ValueError, match="groups must be unavailable"):
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
    relocated = tmp_path / "relocated_report"
    relocated.mkdir()
    shutil.copytree(tmp_path / "runs", relocated / "runs")
    for name in ("busco_comparison.json", "contract.json"):
        shutil.copyfile(tmp_path / name, relocated / name)
    original_comparison_hash = busco.digest(relocated / "busco_comparison.json")
    assert busco.render_existing(relocated) == [row]
    assert busco.digest(relocated / "busco_comparison.json") == original_comparison_hash
    assert row["before_result"]["source"] == str(tmp_path / "runs/Species_a/before/summary.txt")
    original_read = busco.read_result
    before_source = Path(pair["before"])
    before_bytes = before_source.read_bytes()
    def mutate_earlier_input(path):
        result = original_read(path)
        if path == tmp_path / "runs/Species_a/after/summary.txt":
            before_source.write_text(">a\nATGTTTTAA\n")
        return result
    with monkeypatch.context() as patch:
        patch.setattr(busco, "read_result", mutate_earlier_input)
        with pytest.raises(OSError, match="File changed while hashing"):
            busco.render_existing(tmp_path)
    before_source.write_bytes(before_bytes)
    row["delta_complete"] = 99
    busco.atomic_json(tmp_path / "busco_comparison.json", value)
    with pytest.raises(ValueError, match="delta changed"):
        busco.render_existing(tmp_path)
    row["delta_complete"] = 1
    busco.atomic_json(tmp_path / "busco_comparison.json", value)
    Path(pair["after"]).write_text(">a\nATGTTTTAA\n")
    with pytest.raises(ValueError, match="input or score changed"):
        busco.render_existing(tmp_path)


@pytest.mark.parametrize("relation", ["same", "report_inside_output", "output_inside_report"])
@pytest.mark.parametrize("plot_only", [False, True])
def test_busco_cli_keeps_reports_outside_immutable_publications(tmp_path, monkeypatch, capsys, relation, plot_only):
    import sys
    root, report = tmp_path / "refinement", tmp_path / "refinement"
    if relation == "report_inside_output":
        report /= "review"
    elif relation == "output_inside_report":
        root /= "effective"
    args = ["gene_model_refinement_busco.py", "--output", str(root), "--report", str(report)]
    args += ["--plot-only"] if plot_only else ["--lineage", "/lineage", "--download-path", "/db"]
    monkeypatch.setattr(sys, "argv", args)
    def fail(*args, **kwargs):
        raise AssertionError("Unsafe report path must fail before reading or writing a publication")
    monkeypatch.setattr(busco, "render_existing", fail)
    monkeypatch.setattr(busco, "input_pairs", fail)
    with pytest.raises(SystemExit) as exc:
        busco.main()
    assert exc.value.code == 2
    assert "Report must be separate from the immutable refinement tree" in capsys.readouterr().err


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
               run_gene_model_rescue_swissprot="0", gene_model_rescue_dir="/rescue", gene_model_refinement_rescue_dir="",
               gene_model_refinement_inputs="",
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

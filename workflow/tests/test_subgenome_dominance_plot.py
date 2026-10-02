import csv
import json
import shutil
import xml.etree.ElementTree as ET

import pytest

from workflow.support import subgenome_dominance_plot as plot


def analyses():
    rows = [{"metric": metric, "group_id": "G01", "subgenome_a": "A", "subgenome_b": b,
             "tissue": "leaf" if metric.startswith("expression") else "", "n_loci": 60, "n_blocks": 6,
             "effect": -.2, "ci_low": -.35, "ci_high": .1, "p_value": .005,
             "q_value": .02, "multiple_testing_n": 10, "status": "estimated"}
            for metric in plot.METRICS[:2] for b in ("B", "C")]
    return [{"analysis_id": f"a{i}", "species": species, "reference": "Beta", "pair_set": "all_pairs",
             "assignment_scope": "local", "expression_unit": "TPM", "statistics": [dict(row) for row in rows]}
            for i, species in enumerate(("Species_one", "Species_two", "Species_three"))]


def test_three_species_uniform_circles_filters_and_absolute_bounds(tmp_path):
    config = {"metrics": list(plot.METRICS[:2]), "reference": "Beta", "pair_set": "all_pairs", "tissue": "leaf",
              "species_order": ["Species_three", "Species_one", "Species_two"],
              "font_size": 8, "width_pt": 650, "height_pt": 320}
    metadata = plot.comparison_plot(analyses(), tmp_path, config=config)
    assert metadata["point_count"] == 12
    assert metadata["figure_size_pt"] == [650, 320]
    assert metadata["point_unit"] == "group-specific subgenome contrast"
    rows = list(csv.DictReader((tmp_path / "comparison_absolute_points.tsv").open(), delimiter="\t"))
    assert rows[0]["species"] == "Species_three"
    assert {r["q_value"] for r in rows} == {"0.02"}
    assert {r["multiple_testing_n"] for r in rows} == {"10"}
    assert all(float(r["plotted_ci_low"]) == 0 for r in rows)
    assert rows[0]["plotted_effect"] == "20.0" and rows[0]["plotted_unit"] == "percentage_points"
    assert rows[-1]["plotted_effect"] == "0.2" and rows[-1]["plotted_unit"] == "log2_ratio"
    assert (tmp_path / "comparison_absolute.pdf").read_bytes().startswith(b"%PDF")
    tree = ET.parse(tmp_path / "comparison_absolute.svg")
    texts = [el for el in tree.iter() if el.tag.endswith("text")]
    assert any("Benjamini–Hochberg" in "".join(el.itertext()) for el in texts)
    # All point markers reference one circle path; caps and ticks use separate paths.
    curved_paths = {"#" + el.get("id") for el in tree.iter() if el.tag.endswith("path")
                    and el.get("id") and "C" in el.get("d", "")}
    uses = [el for el in tree.iter() if el.tag.endswith("use")
            and el.get("{http://www.w3.org/1999/xlink}href") in curved_paths]
    assert len(uses) == 12
    assert len({el.get("{http://www.w3.org/1999/xlink}href") for el in uses}) == 1


def test_display_filter_does_not_readjust_q_and_limits_are_shared(tmp_path):
    result = plot.comparison_plot(analyses(), tmp_path, config={"formats": ["svg"],
                                  "contrast_anchors": {"Species_one": "B"}, "tissue": "leaf"})
    assert result["point_count"] == 10
    assert result["shared_x_limits"]["retention_difference"] == pytest.approx([-1.4, 39.55])
    assert "unchanged" in result["tests"]


@pytest.mark.parametrize("config", [{"font_size": 0}, {"formats": ["eps"]}, {"typo": True},
                                  {"significance_threshold": float("nan")}, {"species_order": ["x", "x"]}])
def test_invalid_config_fails(config):
    with pytest.raises(ValueError):
        plot.validate_config(config)


def test_missing_font_fails_instead_of_silent_fallback(tmp_path):
    with pytest.raises(ValueError, match="Failed to find font"):
        plot.comparison_plot(analyses(), tmp_path, config={"font_family": "NoSuchFont_GeneGalleon"})


def test_supplied_font_files_take_precedence_over_installed_faces(tmp_path):
    from matplotlib import font_manager
    fonts = []
    for style, weight in (("normal", "normal"), ("italic", "normal"), ("normal", "bold")):
        original = font_manager.findfont(font_manager.FontProperties(family="DejaVu Sans", style=style, weight=weight))
        supplied = tmp_path / f"supplied-{style}-{weight}.ttf"
        shutil.copyfile(original, supplied)
        fonts.append(str(supplied))
    result = plot.comparison_plot(analyses(), tmp_path, config={"font_files": fonts, "formats": ["svg"]})
    assert set(result["font_hashes"]) == set(fonts)
    config = tmp_path / "plot.json"
    config.write_text(json.dumps({"font_files": ["supplied-normal-normal.ttf"]}))
    resolved, hashes = plot.load_config(config)
    assert resolved["font_files"] == [fonts[0]] and hashes[fonts[0]] == plot.sha256(fonts[0])


def test_clipped_intervals_and_empty_filters_fail(tmp_path):
    with pytest.raises(ValueError, match="clips"):
        plot.comparison_plot(analyses(), tmp_path, config={"absolute_x_limits": {"retention_difference": [0, 10]}})
    with pytest.raises(ValueError, match="No analyses"):
        plot.comparison_plot(analyses(), tmp_path, config={"reference": "Missing"})


def test_report_cli_input_contract_and_source_hashes(tmp_path):
    entry = analyses()[0]
    statistics = tmp_path / "statistics.tsv"
    with statistics.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=entry["statistics"][0], delimiter="\t")
        writer.writeheader()
        writer.writerows(entry["statistics"])
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("analysis_id\tspecies\tstatistics_file\treference\tpair_set\n"
                        "a0\tSpecies_one\tstatistics.tsv\tBeta\tall_pairs\n")
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"reference": "Beta", "formats": ["svg", "pdf"]}))
    signed, folded = plot.report(manifest, tmp_path / "plots", config)
    assert not signed["absolute"] and folded["absolute"]
    assert folded["input_hashes"][str(statistics)] == plot.sha256(statistics)
    assert str(config) in folded["input_hashes"]


@pytest.mark.parametrize("key", ["font_size", "point_size", "error_bar_width", "dpi", "width_pt", "height_pt",
                                "individual_width_pt", "individual_height_pt", "capsize", "significance_threshold"])
def test_boolean_numeric_config_is_rejected(key):
    with pytest.raises(ValueError):
        plot.validate_config({key: True})


def test_config_results_cannot_mutate_shared_defaults():
    first = plot.validate_config({})
    first["metrics"].clear()
    first["contrast_anchors"]["species"] = "A"
    second = plot.validate_config({})
    assert second["metrics"] == list(plot.METRICS) and second["contrast_anchors"] == {}


@pytest.mark.parametrize("absolute", [False, True])
@pytest.mark.parametrize("effect,low,high", [(.2, -.1, .1), (-.2, -.1, .1), (.2, .3, .4), (-.2, -.4, -.3)])
def test_confidence_set_need_not_contain_point_estimate(tmp_path, absolute, effect, low, high):
    entry = analyses()[0]
    entry["statistics"] = [{**entry["statistics"][2], "effect": effect, "ci_low": low, "ci_high": high}]
    result = plot.comparison_plot([entry], tmp_path, absolute=absolute,
                                  config={"metrics": ["expression_log2_ratio"], "formats": ["svg"]})
    assert result["point_count"] == 1
    points = list(csv.DictReader((tmp_path / ("comparison_absolute_points.tsv" if absolute else "comparison_points.tsv")).open(), delimiter="\t"))
    expected = plot.display_interval(entry["statistics"][0], absolute)
    assert tuple(float(points[0][k]) for k in ("plotted_effect", "plotted_ci_low", "plotted_ci_high")) == expected
    swapped = {**entry["statistics"][0], "effect": -effect, "ci_low": -high, "ci_high": -low}
    assert plot.display_interval(swapped, True) == plot.display_interval(entry["statistics"][0], True)


def test_unestimable_species_and_metrics_preserve_panel_topology(tmp_path):
    entries = analyses()
    entries[1]["statistics"] = []
    for row in entries[2]["statistics"]:
        row.update(effect=None, ci_low=None, ci_high=None, p_value=None, q_value=None, n_loci=0, n_blocks=0)
    result = plot.comparison_plot(entries, tmp_path, config={"metrics": list(plot.METRICS[:2]), "formats": ["svg"],
                                                           "species_order": [a["species"] for a in entries]})
    assert result["point_count"] == 4 and len(result["panels"]) == 6
    assert sum(p["status"] == "not_estimable" for p in result["panels"]) == 4
    empty = plot.comparison_plot([entries[1]], tmp_path / "empty", config={"formats": ["svg"]})
    assert empty["point_count"] == 0 and len(empty["panels"]) == 3
    assert (tmp_path / "empty/comparison_absolute_points.tsv").read_text().startswith("analysis_id\t")
    with (tmp_path / "empty/comparison_absolute_points.tsv").open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        assert {"n_loci", "n_blocks", "p_value", "plotted_y", "plotted_unit"} <= set(reader.fieldnames)
        assert not list(reader)


@pytest.mark.parametrize("changes", [{"q_value": float("nan")}, {"effect": float("inf")}, {"n_blocks": -1},
                                     {"n_blocks": 61}, {"n_loci": None}, {"q_value": -.1}, {"q_value": .001},
                                     {"ci_low": .3, "ci_high": .1}, {"ci_low": None}, {"metric": "unknown"},
                                     {"subgenome_b": "A"}, {"n_opportunities": 2}])
def test_invalid_statistics_fail_before_render(tmp_path, changes):
    entry = analyses()[0]
    entry["statistics"][0].update(changes)
    with pytest.raises(ValueError):
        plot.comparison_plot([entry], tmp_path)
    assert not list(tmp_path.iterdir())


def report_fixture(tmp_path, statistics_name="statistics.tsv"):
    entry = analyses()[0]
    statistics = tmp_path / statistics_name
    with statistics.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=entry["statistics"][0], delimiter="\t")
        writer.writeheader()
        writer.writerows(entry["statistics"])
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("analysis_id\tspecies\tstatistics_file\n" + f"a0\tSpecies_one\t{statistics_name}\n")
    config = tmp_path / "config.json"
    config.write_text('{"formats": ["svg"]}')
    return manifest, statistics, config


def test_report_render_failure_preserves_previous_report(tmp_path, monkeypatch):
    manifest, _, config = report_fixture(tmp_path)
    output = tmp_path / "plots"
    plot.report(manifest, output, config)
    saved = {p.name: p.read_bytes() for p in output.iterdir()}
    original = plot.comparison_plot
    def fail_absolute(*args, **kwargs):
        if kwargs["absolute"]:
            raise ValueError("forced absolute render failure")
        return original(*args, **kwargs)
    monkeypatch.setattr(plot, "comparison_plot", fail_absolute)
    with pytest.raises(ValueError, match="forced"):
        plot.report(manifest, output, config)
    assert {p.name: p.read_bytes() for p in output.iterdir()} == saved


def test_report_changed_input_preserves_previous_report(tmp_path, monkeypatch):
    manifest, statistics, config = report_fixture(tmp_path)
    output = tmp_path / "plots"
    plot.report(manifest, output, config)
    saved = {p.name: p.read_bytes() for p in output.iterdir()}
    original = plot.comparison_plot
    def change_input(*args, **kwargs):
        result = original(*args, **kwargs)
        if kwargs["absolute"]:
            statistics.write_text(statistics.read_text() + "\n")
        return result
    monkeypatch.setattr(plot, "comparison_plot", change_input)
    with pytest.raises(ValueError, match="input changed"):
        plot.report(manifest, output, config)
    assert {p.name: p.read_bytes() for p in output.iterdir()} == saved


def test_report_publication_failure_rolls_back_all_files(tmp_path, monkeypatch):
    manifest, _, config = report_fixture(tmp_path)
    output = tmp_path / "plots"
    plot.report(manifest, output, config)
    saved = {p.name: p.read_bytes() for p in output.iterdir()}
    original = plot.os.replace
    failed = False
    def fail_one(source, target):
        nonlocal failed
        if not failed and source.name == "comparison_absolute.svg" and source.parent.name != "plots":
            failed = True
            raise OSError("forced publication failure")
        return original(source, target)
    monkeypatch.setattr(plot.os, "replace", fail_one)
    with pytest.raises(OSError, match="forced"):
        plot.report(manifest, output, config)
    assert failed and {p.name: p.read_bytes() for p in output.iterdir()} == saved


def test_failed_rollback_retains_previous_files_for_recovery(tmp_path, monkeypatch):
    manifest, _, config = report_fixture(tmp_path)
    output = tmp_path / "plots"
    plot.report(manifest, output, config)
    saved = (output / "comparison.svg").read_bytes()
    original = plot.os.replace
    def persistent_failure(source, target):
        if source.name == "comparison_absolute.svg" and source.parent.name != "plots" and not source.parent.name.startswith(".subgenome-report-recovery-"):
            raise OSError("forced publish failure")
        if source.name == "comparison.svg" and source.parent.name.startswith(".subgenome-report-recovery-"):
            raise OSError("forced rollback failure")
        return original(source, target)
    monkeypatch.setattr(plot.os, "replace", persistent_failure)
    with pytest.raises(RuntimeError, match="rollback incomplete.*recovery"):
        plot.report(manifest, output, config)
    recoveries = list(tmp_path.glob(".subgenome-report-recovery-*"))
    assert len(recoveries) == 1 and (recoveries[0] / "comparison.svg").read_bytes() == saved


def test_report_removes_obsolete_formats_and_preserves_unrelated_files(tmp_path):
    manifest, _, config = report_fixture(tmp_path)
    output = tmp_path / "plots"
    config.write_text('{"formats": ["svg", "pdf"]}')
    plot.report(manifest, output, config)
    (output / "notes.txt").write_text("preserve")
    config.write_text('{"formats": ["svg"]}')
    plot.report(manifest, output, config)
    assert not list(output.glob("*.pdf")) and (output / "notes.txt").read_text() == "preserve"
    completion = json.loads((output / "report_manifest.json").read_text())
    assert all(plot.sha256(output / name) == digest for name, digest in completion["output_hashes"].items())


def test_report_input_output_collision_is_rejected_before_writing(tmp_path):
    manifest, statistics, config = report_fixture(tmp_path, "comparison_points.tsv")
    saved = statistics.read_bytes()
    with pytest.raises(ValueError, match="overwrite an input"):
        plot.report(manifest, tmp_path, config)
    assert statistics.read_bytes() == saved


def test_report_respects_output_lock(tmp_path):
    import fcntl
    import hashlib
    manifest, _, config = report_fixture(tmp_path)
    output = tmp_path / "plots"
    key = hashlib.sha256(str(output.resolve()).encode()).hexdigest()[:20]
    with (tmp_path / f".subgenome-report-{key}.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(ValueError, match="owns the output lock"):
            plot.report(manifest, output, config)
    assert not output.exists()


@pytest.mark.parametrize("corruption", ["duplicate_header", "short_row", "extra_cell", "invalid_number"])
def test_report_malformed_tables_fail_cleanly(tmp_path, corruption):
    manifest, statistics, config = report_fixture(tmp_path)
    lines = statistics.read_text().splitlines()
    if corruption == "duplicate_header":
        lines[0] = lines[0].replace("group_id", "metric")
    elif corruption == "short_row":
        lines[1] = "\t".join(lines[1].split("\t")[:-1])
    elif corruption == "extra_cell":
        lines[1] += "\textra"
    else:
        lines[1] = lines[1].replace("0.02", "nan")
    statistics.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError):
        plot.report(manifest, tmp_path / "plots", config)
    assert not (tmp_path / "plots").exists()


def test_reversed_duplicate_contrast_is_rejected(tmp_path):
    entry = analyses()[0]
    row = entry["statistics"][0]
    entry["statistics"].append({**row, "subgenome_a": row["subgenome_b"], "subgenome_b": row["subgenome_a"],
                               "effect": -row["effect"], "ci_low": -row["ci_high"], "ci_high": -row["ci_low"]})
    with pytest.raises(ValueError, match="Repeated comparison"):
        plot.comparison_plot([entry], tmp_path)


def test_missing_tissue_is_checked_within_selected_analysis_metadata(tmp_path):
    entries = analyses()
    entries[0]["reference"] = "Other"
    for row in entries[0]["statistics"]:
        if row["tissue"]:
            row["tissue"] = "flower"
    with pytest.raises(ValueError, match="tissue is absent"):
        plot.comparison_plot(entries, tmp_path, config={"reference": "Beta", "tissue": "flower"})


@pytest.mark.parametrize("anchors", [{"Typo_species": "A"}, {"Species_one": "D"}])
def test_unknown_anchor_species_or_label_fails(tmp_path, anchors):
    with pytest.raises(ValueError, match="unknown species|anchor is absent"):
        plot.comparison_plot(analyses(), tmp_path, config={"contrast_anchors": anchors})


def test_legacy_statistics_without_new_inference_columns_remain_readable(tmp_path):
    manifest, statistics, config = report_fixture(tmp_path)
    with statistics.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    rows = [{k: v for k, v in r.items() if k not in {"multiple_testing_n", "p_value"}} for r in rows]
    with statistics.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    _, result = plot.report(manifest, tmp_path / "plots", config)
    assert result["point_count"] == 4

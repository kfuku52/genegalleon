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
             "assignment_scope": "local", "expression_unit": "TPM", "statistics": rows}
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
    with pytest.raises(ValueError, match="No estimable"):
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

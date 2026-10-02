"""Real NWKIT rendering and reconciliation-guided root selection in the GG runtime."""

import json
import shlex
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / "workflow" / "support"


def run(*args):
    return subprocess.run([str(value) for value in args], capture_output=True, text=True)


def test_dated_tree_rejects_duplicate_tip_names(tmp_path):
    from workflow.support.dated_tree_presentation import render_dated_tree

    tree, plot = tmp_path / "tree.nwk", tmp_path / "plot.svg"
    tree.write_text("(A:10,A:10);")
    plot.write_text("previous plot")
    with pytest.raises(ValueError, match="unique"):
        render_dated_tree(tree, plot, geological_background="none")
    assert plot.read_text() == "previous plot"


@pytest.mark.parametrize("table, content", [
    ("tip_order", "species_id\tspecies_id\nwrong\tA\nwrong\tB\n"),
    ("tip_annotations", "species_id\tcolour\tcolour\nA\tinvalid\tblack\n"),
    ("age_clades", "descendant_species\tdescendant_species\nwrong\tA,B\n"),
    ("busco_summary", "Species\tbusco_cds_single\tbusco_cds_single\tbusco_cds_duplicated\tbusco_cds_fragmented\tbusco_cds_missing\tbusco_cds_total\nA\t99\t1\t0\t0\t0\t1\nB\t99\t1\t0\t0\t0\t1\n"),
    ("tip_order", "species_id\nA\tunlabelled extra field\nB\tunlabelled extra field\n"),
])
def test_dated_tree_rejects_ambiguous_metadata_without_replacing_plot(tmp_path, table, content):
    from workflow.support.dated_tree_presentation import render_dated_tree

    tree, plot, metadata = tmp_path / "tree.nwk", tmp_path / "plot.svg", tmp_path / "metadata.tsv"
    tree.write_text("(A:10,B:10);")
    metadata.write_text(content)
    plot.write_text("previous plot")
    with pytest.raises(ValueError, match="header|row width"):
        render_dated_tree(tree, plot, geological_background="none", **{table: metadata})
    assert plot.read_text() == "previous plot"


def test_dated_tree_restores_matplotlib_style_after_render(tmp_path):
    import matplotlib

    from workflow.support.dated_tree_presentation import render_dated_tree

    tree = tmp_path / "tree.nwk"
    tree.write_text("(A:10,B:10);")
    with matplotlib.rc_context({"font.family": "serif", "font.size": 13, "pdf.fonttype": 3, "svg.fonttype": "path"}):
        before = {key: matplotlib.rcParams[key] for key in ("font.family", "font.size", "pdf.fonttype", "svg.fonttype")}
        render_dated_tree(tree, tmp_path / "plot.svg", geological_background="none")
        assert {key: matplotlib.rcParams[key] for key in before} == before


def test_dated_tree_mixed_interval_legends_fit_narrow_page(tmp_path):
    from workflow.support.dated_tree_presentation import render_dated_tree

    tree = tmp_path / "tree.nhx"
    tree.write_text("((A:10,B:10):20[&&NHX:age=10:age_ci_low=8:age_ci_high=12:age_ci_kind=equal-tail:age_ci_level=0.95],C:30)[&95%HPD={28,32}];")
    report = render_dated_tree(tree, tmp_path / "plot.svg", figure_width=3.6, geological_background="none")
    assert report["credible_interval_count"] == 2
    assert report["legend_rows"] == 2
    assert "equal-tailed credible intervals" in report["credible_interval_label"]


def test_dated_tree_multiline_branch_legends_do_not_overlap(tmp_path):
    from workflow.support.dated_tree_presentation import render_dated_tree

    tree, annotations = tmp_path / "tree.nwk", tmp_path / "branches.tsv"
    tree.write_text("(A:10,B:10);")
    annotations.write_text('descendant_species\tlabel\tsymbol\nA\t"First event\nwith explanation\nand third line"\t^\nB\t"Second event with\na long explanation\nand another line"\tx\n')
    report = render_dated_tree(tree, tmp_path / "plot.svg", branch_annotations=annotations, figure_width=3.6,
                               geological_background="none")
    boxes = report["legend_bbox_points"]
    assert all(box[3] < other[1] or box[1] > other[3]
               for index, box in enumerate(boxes) for other in boxes[index + 1:])


def test_dated_tree_cannot_replace_any_consumed_busco_header(tmp_path):
    from workflow.support.dated_tree_presentation import render_dated_tree

    tree, summary = tmp_path / "tree.nwk", tmp_path / "summary.tsv"
    tree.write_text("(A:10,B:10);")
    summary.write_text("Species\tbusco_cds_single\tbusco_cds_duplicated\tbusco_cds_fragmented\tbusco_cds_missing\tbusco_cds_total\nA\t1\t0\t0\t0\t1\nB\t1\t0\t0\t0\t1\n")
    results = tmp_path / "busco"
    results.mkdir()
    original = "# The lineage dataset is: embryophyta_odb12\n"
    for species in "AB":
        for suffix in ("short.txt", "full.tsv"):
            (results / f"{species}.busco.{suffix}").write_text(original)
    header = results / "A.busco.short.txt"
    with pytest.raises(ValueError, match="replace|input"):
        render_dated_tree(tree, tmp_path / "plot.svg", busco_summary=summary, busco_results=results,
                          layout_report=header, geological_background="none")
    assert header.read_text() == original


def test_dated_tree_cannot_replace_geological_boundary_table(tmp_path, monkeypatch):
    from workflow.support import dated_tree_presentation as presentation

    tree, table = tmp_path / "tree.nwk", tmp_path / "periods.tsv"
    tree.write_text("(A:10,B:10);")
    original = presentation.GEOLOGICAL_DATA.read_bytes()
    table.write_bytes(original)
    monkeypatch.setattr(presentation, "GEOLOGICAL_DATA", table)
    with pytest.raises(ValueError, match="replace|input"):
        presentation.render_dated_tree(tree, tmp_path / "plot.svg", layout_report=table)
    assert table.read_bytes() == original


def test_dated_tree_branch_symbols_near_present_are_complete_and_separate_from_age_labels(tmp_path, monkeypatch):
    from matplotlib.figure import Figure

    from workflow.support.dated_tree_presentation import render_dated_tree

    tree, annotations = tmp_path / "tree.nwk", tmp_path / "branches.tsv"
    tree.write_text("((A:10,B:10):20,C:30);")
    annotations.write_text("descendant_species\tlabel\tsymbol\tbranch_fraction\nA\tRecent event\tx\t0.001\nA,B\tEarlier event\t^\t0.02\n")
    save = Figure.savefig
    observed = []

    def inspect(figure, *args, **kwargs):
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        markers = [line for line in figure.axes[0].lines if (line.get_gid() or "").startswith("branch-event-")]
        assert len(markers) == 2 and all(not line.get_clip_on() for line in markers)
        for text in figure.axes[0].texts:
            if text.get_text() in {"10.0", "30.0"}:
                assert all(not text.get_window_extent(renderer).overlaps(line.get_window_extent(renderer))
                           for line in markers)
        observed.append(True)
        return save(figure, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", inspect)
    report = render_dated_tree(tree, tmp_path / "plot.svg", branch_annotations=annotations, node_ages="all",
                               geological_background="none")
    assert observed == [True] and report["mean_age_label_count"] == 2


def test_dated_tree_rejects_overlapping_branch_symbols_without_replacing_outputs(tmp_path):
    from workflow.support.dated_tree_presentation import render_dated_tree

    tree, annotations, plot = tmp_path / "tree.nwk", tmp_path / "branches.tsv", tmp_path / "plot.svg"
    tree.write_text("(A:10,B:10);")
    annotations.write_text("descendant_species\tlabel\tsymbol\tbranch_fraction\nA\tGain\t^\t0.5\nA\tLoss\tx\t0.501\n")
    plot.write_text("previous plot")
    with pytest.raises(ValueError, match="Branch symbols overlap"):
        render_dated_tree(tree, plot, branch_annotations=annotations, geological_background="none")
    assert plot.read_text() == "previous plot"


@pytest.mark.parametrize("width", [3.6, 4.8])
def test_dated_tree_deep_time_source_footer_wraps_within_page(tmp_path, width):
    from workflow.support.dated_tree_presentation import render_dated_tree

    tree = tmp_path / "tree.nwk"
    tree.write_text("(A:4500,B:4500);")
    report = render_dated_tree(tree, tmp_path / "plot.svg", figure_width=width, show_geological_source=True)
    assert report["geological_source_credit_visible"] is True


@pytest.mark.parametrize("width", [3.6, 4.8, 6.0, 7.2])
def test_dated_tree_publication_preserves_intervals_dataset_and_species_rows(tmp_path, width):
    import xml.etree.ElementTree as ET

    tree = tmp_path / "tree.nwk"
    tree.write_text("((A_a:10,B_b:10):140,C_c:150)[&95%HPD={145,155}];")
    summary = tmp_path / "annotation_summary.tsv"
    summary.write_text(
        "Species\tbusco_cds_single\tbusco_cds_duplicated\tbusco_cds_fragmented\tbusco_cds_missing\tbusco_cds_total\tbusco_cds_lineage\n"
        "A a\t8\t1\t0\t1\t10\tembryophyta_odb12\n"
        "C c\t6\t2\t1\t1\t10\tembryophyta_odb12\n"
        "B b\t7\t0\t2\t1\t10\tembryophyta_odb12\n"
    )
    plot, report = tmp_path / "plot.svg", tmp_path / "report.json"
    result = run(
        sys.executable,
        SUPPORT / "plot_dated_tree.py",
        "--infile",
        tree,
        "--outfile",
        plot,
        "--busco-summary",
        summary,
        "--layout-report",
        report,
        "--geological-background",
        "period",
        "--figure-width",
        width,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    texts = [item.text for item in ET.parse(plot).iter("{http://www.w3.org/2000/svg}text")]
    assert "95% highest posterior density intervals" in texts
    assert "Number of BUSCO genes" in texts and "(embryophyta_odb12)" in texts
    data = json.loads(report.read_text())
    assert data["geological_label_placement"] == "above_tree"
    assert data["geological_source_credit_visible"] is False
    assert not any(text and text.startswith("Geological periods:") for text in texts)
    assert len(data["geological_labels"]) == len(data["geological_intervals"])
    for period, label in zip(data["geological_intervals"], data["geological_labels"], strict=True):
        assert label["name"] == period["name"]
        assert texts.count(period["name"]) == 1
        assert label["anchor_Ma"] == (period["young_Ma"] + period["old_Ma"]) / 2
        assert label["bbox_points"][1] > data["tree_plot_bbox_points"][3]
        assert label["bbox_points"][1] - data["tree_plot_bbox_points"][3] == pytest.approx(4)
    colours = [p["colour"] for p in data["geological_intervals"]]
    assert colours == (["#F7F7F7", "#E7E7E7"] * len(colours))[:len(colours)]
    assert data["credible_interval_count"] == 1
    assert data["busco_counts"]["C_c"] == [6, 2, 1, 1]
    assert data["tip_y_coordinates"] == {"A_a": 2, "B_b": 1, "C_c": 0}
    assert data["geological_intervals"][-1]["young_Ma"] == 143.1
    assert data["all_ages_Ma"][0]["mean"] == 150
    assert data["all_ages_Ma"][0]["low"] == 145
    assert data["all_ages_Ma"][0]["high"] == 155
    assert data["figure_size_inches"][0] == width
    assert data["font_size_points"] == 8
    assert data["tree_x_axis_position"] == "bottom"
    assert data["tree_x_axis_y_points"] == pytest.approx(data["busco_plot_bbox_points"][1])
    assert data["tree_x_axis_y_points"] < data["busco_percentage_axis_y_points"]
    assert data["mean_age_label_count"] == 0
    assert data["credible_interval_style"]["alpha"] is None
    assert data["credible_interval_style"]["zorder"] < data["credible_interval_style"]["tree_zorder"]
    assert data["credible_interval_style"]["bar_capstyle"] == "butt"
    assert data["credible_interval_style"]["cap_linewidth_points"] == data["credible_interval_style"]["bar_linewidth_points"]
    groups = [item.attrib.get("id", "") for item in ET.parse(plot).iter("{http://www.w3.org/2000/svg}g")]
    assert groups.index("age-interval-0") < groups.index("tree-branch-0")
    svg_groups = {item.attrib.get("id"): item for item in ET.parse(plot).iter("{http://www.w3.org/2000/svg}g")}
    bar = next(svg_groups["age-interval-0"].iter("{http://www.w3.org/2000/svg}path"))
    bar_style = dict(entry.strip().split(": ", 1) for entry in bar.attrib["style"].split(";") if entry.strip())
    assert bar_style.get("stroke-linecap", "butt") == "butt"
    for cap in svg_groups["age-interval-cap-0"].iter("{http://www.w3.org/2000/svg}use"):
        cap_style = dict(entry.strip().split(": ", 1) for entry in cap.attrib["style"].split(";") if entry.strip())
        assert cap_style["stroke-width"] == bar_style["stroke-width"]
    for item in ET.parse(plot).iter():
        if item.attrib.get("id", "").startswith("age-interval-"):
            assert all("opacity" not in value for child in item.iter() for value in child.attrib.values())
    legend_boxes = data["legend_bbox_points"]
    assert all(box[3] < other[1] or box[1] > other[3]
               for index, box in enumerate(legend_boxes) for other in legend_boxes[index + 1:])
    assert tree.read_text().endswith("[&95%HPD={145,155}];")


def test_dated_tree_source_footer_is_opt_in_and_reclaims_physical_space(tmp_path):
    import xml.etree.ElementTree as ET

    tree = tmp_path / "tree.nwk"
    tree.write_text("(" + ",".join(f"Species_{i}:200" for i in range(20)) + ");")
    reports = []
    for visible in (False, True):
        plot, report = tmp_path / f"credit-{visible}.svg", tmp_path / f"credit-{visible}.json"
        command = [sys.executable, SUPPORT / "plot_dated_tree.py", "--infile", tree,
                   "--outfile", plot, "--geological-background", "period", "--layout-report", report]
        if visible:
            command.append("--show-geological-source")
        result = run(*command)
        assert result.returncode == 0, result.stdout + result.stderr
        texts = [item.text for item in ET.parse(plot).iter("{http://www.w3.org/2000/svg}text")]
        assert any(text and text.startswith("Geological periods:") for text in texts) is visible
        data = json.loads(report.read_text())
        assert data["geological_source_credit_visible"] is visible
        assert data["geological_source"].endswith("ChronostratChart2026-06.pdf")
        reports.append(data)
    assert (reports[1]["figure_size_inches"][1] - reports[0]["figure_size_inches"][1]) * 72 == pytest.approx(20)
    assert [report["row_spacing_points"] for report in reports] == pytest.approx([9, 9])
    assert reports[0]["all_ages_Ma"] == reports[1]["all_ages_Ma"]


def test_dated_tree_branch_events_select_exact_stems_and_share_legends(tmp_path):
    import xml.etree.ElementTree as ET

    tree, annotations = tmp_path / "tree.nwk", tmp_path / "branches.tsv"
    original = "((A_a:10,B_b:10):140,C_c:150)[&95%HPD={145,155}];"
    tree.write_text(original)
    annotations.write_text(
        "event_id\tdescendant_species\tlabel\tsymbol\tbranch_fraction\n"
        "gain-ab\tB_b,A_a\tCarnivory gain\t^\t\n"
        "loss-a\tA_a\tCarnivory loss\tx\t0.75\n"
        "loss-b\tB_b\tCarnivory loss\tx\t\n")
    plot, report = tmp_path / "plot.svg", tmp_path / "report.json"
    result = run(sys.executable, SUPPORT / "plot_dated_tree.py", "--infile", tree,
                 "--outfile", plot, "--branch-annotations", annotations, "--layout-report", report)
    assert result.returncode == 0, result.stdout + result.stderr
    root = ET.parse(plot)
    texts = [item.text for item in root.iter("{http://www.w3.org/2000/svg}text")]
    assert texts.count("Carnivory gain") == texts.count("Carnivory loss") == 1
    groups = {item.attrib.get("id"): item for item in root.iter("{http://www.w3.org/2000/svg}g")}
    assert all("branch-event-" + name in groups for name in ("gain-ab", "loss-a", "loss-b"))
    data = json.loads(report.read_text())
    assert [event["display_position_Ma"] for event in data["branch_annotations"]] == [80, 7.5, 5]
    assert [event["y"] for event in data["branch_annotations"]] == [1.5, 2, 1]
    assert data["credible_interval_count"] == 1
    assert data["all_ages_Ma"][0]["mean"] == 150
    boxes = data["legend_bbox_points"]
    assert len(boxes) == 2 and boxes[0][3] < boxes[1][1]
    assert tree.read_text() == original


@pytest.mark.parametrize("rows, error", [
    ("A_a,C_c\tGain\t^\t0.5", "exact non-root"),
    ("A_a,B_b,C_c\tGain\t^\t0.5", "exact non-root"),
    ("A_a\tGain\t^\t0", "Branch fraction"),
    ("A_a\tGain\t^\tnan", "Branch fraction"),
    ("A_a\tGain\tbad\t0.5", "Invalid branch"),
    ("A_a\tGain\t^\t0.5\nB_b\tGain\tx\t0.5", "conflicting symbols"),
    ("A_a\tGain\t^\t0.5\nA_a\tLoss\tx\t0.5", "same display position"),
])
def test_dated_tree_rejects_bad_branch_annotations_preserving_outputs(tmp_path, rows, error):
    tree, annotations = tmp_path / "tree.nwk", tmp_path / "branches.tsv"
    tree.write_text("((A_a:10,B_b:10):140,C_c:150);")
    annotations.write_text("descendant_species\tlabel\tsymbol\tbranch_fraction\n" + rows + "\n")
    plot, report = tmp_path / "plot.svg", tmp_path / "report.json"
    plot.write_text("old plot")
    report.write_text("old report")
    result = run(sys.executable, SUPPORT / "plot_dated_tree.py", "--infile", tree,
                 "--outfile", plot, "--branch-annotations", annotations, "--layout-report", report)
    assert result.returncode != 0 and error in result.stderr
    assert plot.read_text() == "old plot" and report.read_text() == "old report"


def test_dated_tree_busco_header_discovery_and_mixed_dataset_rejection(tmp_path):
    sys.path.insert(0, str(SUPPORT))
    from dated_tree_presentation import read_busco

    summary_dir = tmp_path / "annotation_summary"
    summary_dir.mkdir()
    summary = summary_dir / "annotation_summary.tsv"
    summary.write_text(
        "Species\tbusco_cds_single\tbusco_cds_duplicated\tbusco_cds_fragmented\tbusco_cds_missing\tbusco_cds_total\n"
        "A\t1\t0\t0\t0\t1\nB\t1\t0\t0\t0\t1\n"
    )
    results = tmp_path / "species_cds_busco_full"
    results.mkdir()
    for species in "AB":
        (results / (species + ".busco.full.tsv")).write_text(
            "# The lineage dataset is: embryophyta_odb12 (number of BUSCOs: 1)\n"
        )
    counts, dataset, sources = read_busco(summary, ["A", "B"])
    assert dataset == "embryophyta_odb12" and len(sources) == 2 and counts["A"] == [1, 0, 0, 0]
    (results / "B.busco.full.tsv").write_text("# The lineage dataset is: embryophyta_odb10\n")
    tree, plot, report = tmp_path / "tree.nwk", tmp_path / "plot.pdf", tmp_path / "report.json"
    tree.write_text("(A:10,B:10);")
    plot.write_bytes(b"old plot")
    report.write_text("old report")
    result = run(
        sys.executable,
        SUPPORT / "plot_dated_tree.py",
        "--infile",
        tree,
        "--outfile",
        plot,
        "--busco-summary",
        summary,
        "--layout-report",
        report,
    )
    assert result.returncode != 0 and "Mixed BUSCO" in result.stderr
    assert plot.read_bytes() == b"old plot" and report.read_text() == "old report"


@pytest.mark.parametrize("missing", ["NA", "", "n/a"])
def test_dated_tree_busco_missing_lineage_never_becomes_a_dataset_name(tmp_path, missing):
    from workflow.support.dated_tree_presentation import read_busco

    summary = tmp_path / "summary.tsv"
    summary.write_text("Species\tbusco_cds_single\tbusco_cds_duplicated\tbusco_cds_fragmented\tbusco_cds_missing\tbusco_cds_total\tbusco_cds_lineage\n"
                       f"A\t1\t0\t0\t0\t1\t{missing}\nB\t1\t0\t0\t0\t1\t{missing}\n")
    assert read_busco(summary, ["A", "B"])[1] is None
    for species in "AB":
        (tmp_path / f"{species}.busco.short.txt").write_text("# The lineage dataset is: embryophyta_odb12\n")
    assert read_busco(summary, ["A", "B"], tmp_path)[1] == "embryophyta_odb12"


def test_dated_tree_busco_empty_header_cannot_consume_the_following_line(tmp_path):
    from workflow.support.dated_tree_presentation import read_busco

    summary = tmp_path / "summary.tsv"
    summary.write_text("Species\tbusco_cds_single\tbusco_cds_duplicated\tbusco_cds_fragmented\tbusco_cds_missing\tbusco_cds_total\nA\t1\t0\t0\t0\t1\nB\t1\t0\t0\t0\t1\n")
    for species in "AB":
        (tmp_path / f"{species}.busco.short.txt").write_text("# The lineage dataset is:\n# Busco id\tStatus\n")
    assert read_busco(summary, ["A", "B"], tmp_path)[1] is None


def test_dated_tree_geological_header_handles_narrow_periods_on_deep_time_tree(tmp_path):
    tree, plot, report = tmp_path / "tree.nwk", tmp_path / "plot.svg", tmp_path / "report.json"
    tree.write_text("(A:4500,B:4500);")
    result = run(sys.executable, SUPPORT / "plot_dated_tree.py", "--infile", tree, "--outfile", plot,
                 "--geological-background", "period", "--layout-report", report)
    assert result.returncode == 0, result.stdout + result.stderr
    data = json.loads(report.read_text())
    assert len(data["geological_labels"]) == 13
    assert any(abs(label["horizontal_offset_points"]) > 1 for label in data["geological_labels"])
    boxes = [label["bbox_points"] for label in data["geological_labels"]]
    assert all(box[1] > data["tree_plot_bbox_points"][3] for box in boxes)
    for index, box in enumerate(boxes):
        assert all(box[2] <= other[0] or box[0] >= other[2] for other in boxes[index + 1:])
    assert tree.read_text() == "(A:4500,B:4500);"


@pytest.mark.parametrize("background, has_periods", [("none", False), (None, True)])
def test_dated_tree_layout_report_selects_presentation_with_optional_background(tmp_path, background, has_periods):
    tree, plot, report = tmp_path / "tree.nwk", tmp_path / "plot.pdf", tmp_path / "report.json"
    tree.write_text("(A:10,B:10);")
    command = [
        sys.executable,
        SUPPORT / "plot_dated_tree.py",
        "--infile",
        tree,
        "--outfile",
        plot,
        "--layout-report",
        report,
    ]
    if background is not None:
        command.extend(["--geological-background", background])
    result = run(*command)
    assert result.returncode == 0, result.stdout + result.stderr
    assert bool(json.loads(report.read_text())["geological_intervals"]) is has_periods
    assert json.loads(report.read_text())["figure_size_inches"][0] == 4.8


@pytest.mark.parametrize("text", [
    "(A:1,B:1)", "(A:1,B:1);garbage", "(A:1,B:1);(C:1,D:1);", "();",
    "(A:1,B:1));", "#NEXUS\nBEGIN TREES;\nEND;", "(A:-1,B:1);",
])
def test_dated_tree_plot_rejects_invalid_inputs_without_replacing_plot(tmp_path, text):
    source, plot = tmp_path / "tree.nwk", tmp_path / "tree.pdf"
    source.write_text(text)
    plot.write_bytes(b"previous plot")
    result = run(sys.executable, SUPPORT / "plot_dated_tree.py", "--infile", source, "--outfile", plot)
    assert result.returncode != 0
    assert plot.read_bytes() == b"previous plot"


@pytest.mark.parametrize("text", [
    "(A:10,B:10);",
    "((A:10,B:10)named:20[&&NHX:age=10:age_ci_low=8:age_ci_high=12:age_ci_kind=HPD:age_ci_level=0.95],C:30);",
    "#NEXUS\nBEGIN TREES;\nTREE dated = [&R] (A:10,B:10)[&95%HPD={8,12}];\nEND;",
])
def test_dated_tree_plot_accepts_public_unit_nhx_and_figtree(tmp_path, text):
    source, plot = tmp_path / "tree.nwk", tmp_path / "tree.pdf"
    source.write_text(text)
    result = run(sys.executable, SUPPORT / "plot_dated_tree.py", "--infile", source, "--outfile", plot)
    assert result.returncode == 0, result.stdout + result.stderr
    assert plot.read_bytes().startswith(b"%PDF-")
    assert b"NWKIT" in plot.read_bytes()
    assert source.read_text() == text


def test_dated_tree_plot_drops_ci_when_rounded_tree_is_slightly_non_ultrametric(tmp_path):
    source, plot = tmp_path / "rounded.nhx", tmp_path / "rounded.pdf"
    source.write_text(
        "((A:10.000000,B:9.999000)named:20.000000[&&NHX:age=10:age_ci_low=8:"
        "age_ci_high=12:age_ci_kind=HPD:age_ci_level=0.95],C:30.000000);"
    )
    result = run(sys.executable, SUPPORT / "plot_dated_tree.py", "--infile", source, "--outfile", plot)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "slightly non-ultrametric" in result.stderr
    assert plot.read_bytes().startswith(b"%PDF-")
    assert source.read_text().startswith("((A:10.000000")


@pytest.mark.parametrize("candidate_method,reason", [
    ("mad", "mad_compatible_with_reconciliation"),
    ("midpoint", "midpoint_compatible_with_reconciliation"),
    ("outgroup", "first_reconciliation_candidate"),
])
def test_nwkit_root_selection_keeps_reconciliation_priority_and_branch_lengths(tmp_path, candidate_method, reason):
    from nwkit.clade_mapping import projected_root_split
    from nwkit.util import read_tree

    source = tmp_path / "tree.nwk"
    source.write_text("((A:7.8,B:0.949):0.244,(C:7.8,D:1.439):2.87);")
    candidates = tmp_path / "candidates.nwk"
    candidate = candidates
    if candidate_method:
        args = ["nwkit", "root", "--method", candidate_method, "--infile", source, "--outfile", candidate]
        if candidate_method == "outgroup":
            args += ["--outgroup", "D"]
        result = run(*args)
        assert result.returncode == 0, result.stderr
    output, table, plot = tmp_path / "selected.nwk", tmp_path / "roots.tsv", tmp_path / "roots.pdf"
    result = run(sys.executable, SUPPORT / "species_tree_guided_gene_tree_rooting.py",
                 "--in-tree", source, "--candidate-trees", candidates, "--out-tree", output,
                 "--comparison-table", table, "--comparison-plot", plot)
    assert result.returncode == 0, result.stdout + result.stderr
    assert f"Selected root: {reason}" in result.stdout
    before, after = read_tree(str(source), "auto", True), read_tree(str(output), "auto", True)
    assert set(before.leaf_names()) == set(after.leaf_names()) == set("ABCD")
    for a in "ABCD":
        for b in "ABCD":
            assert before.get_distance(a, b) == pytest.approx(after.get_distance(a, b))
    if candidate_method:
        candidate_tree = read_tree(str(candidate), "auto", True)
        assert projected_root_split(after, frozenset("ABCD")) == projected_root_split(candidate_tree, frozenset("ABCD"))
    assert table.read_text().startswith("method\tstatus\t")
    assert plot.read_bytes().startswith(b"%PDF-")


@pytest.mark.parametrize("mode", ["dna", "pep"])
def test_root_core_stage_invalidates_engine_and_preserves_previous_bundle_on_failure(tmp_path, mode):
    core = (ROOT / "workflow/core/gg_genome_evolution_core.sh").read_text()
    begin = core.index("busco_species_tree_assisted_gene_tree_rooting() {")
    function = core[begin:core.index("\nbusco_grampa()", begin)]
    title = "DNA" if mode == "dna" else "protein"
    begin = core.index(f'task="Species-tree-guided gene tree rooting of duplicate-containing BUSCO {title} trees"')
    stage = core[begin:core.index('\ntask=', begin + 1)]
    original = "((A:7.8,B:0.949):0.244,(C:7.8,D:1.439):2.87);"
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "OG1.busco.nwk").write_text(original)
    archives = tmp_path / "archives"
    archives.mkdir()
    archive = archives / "OG1.busco.roots.nwk"

    def write_candidate(text):
        archive.write_text(text)

    write_candidate(original)
    output = tmp_path / "output"
    script = tmp_path / "stage.sh"
    script.write_text(f'''set -euo pipefail
gg_support_dir={shlex.quote(str(SUPPORT))}
source "${{gg_support_dir}}/gg_util.sh"
gg_workspace_dir={shlex.quote(str(tmp_path))}
gg_workspace_output_dir={shlex.quote(str(output))}
genome_evolution_provenance_dir="${{gg_workspace_output_dir}}/provenance"
genome_nwkit_identity="${{1:-test-engine-one}}"
file_dated_species_tree={shlex.quote(str(inputs / "OG1.busco.nwk"))}
dir_busco_iqtree_{mode}={shlex.quote(str(inputs))}
dir_busco_reconciliation_{mode}={shlex.quote(str(archives))}
dir_busco_rooted_txt_{mode}="${{gg_workspace_output_dir}}/reports"
dir_busco_rooted_nwk_{mode}="${{gg_workspace_output_dir}}/trees"
run_busco_dupaware_root_{mode}=1
species_label_parser=taxonomic
GG_GENOME_PARALLEL_JOBS=1
artifact_stale_policy=rebuild
gg_step_start() {{ echo run >> stage-runs.txt; }}
gg_step_skip() {{ echo skip >> stage-runs.txt; }}
{function}
{stage}
''')

    def execute(identity="test-engine-one"):
        return subprocess.run(["bash", str(script), identity], cwd=tmp_path, capture_output=True, text=True)

    for _ in range(2):
        result = execute()
        assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "stage-runs.txt").read_text().splitlines() == ["run", "skip"]
    result = execute("test-engine-two")
    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "stage-runs.txt").read_text().splitlines() == ["run", "skip", "run"]
    saved = {p: p.read_bytes() for folder in (output / "reports", output / "trees") for p in folder.glob("OG1*")}
    assert len(saved) == 4
    manifest = output / "provenance" / f"busco.root_{mode}.json"
    assert json.loads(manifest.read_text())["step"] == f"genome_evolution_busco_root_{mode}"
    saved_manifest = manifest.read_bytes()
    write_candidate("((A:1,B:1):1,(C:1,WrongTip:1):1);")
    failed = execute("test-engine-two")
    assert failed.returncode != 0
    assert {p: p.read_bytes() for p in saved} == saved
    assert manifest.read_bytes() == saved_manifest


@pytest.mark.parametrize("invalid", ["missing_file", "empty_file", "duplicate_tips", "different_tips", "unresolved_root"])
def test_root_selection_rejects_bad_candidates_and_preserves_outputs(tmp_path, invalid):
    source = tmp_path / "tree.nwk"
    source.write_text("((A:7.8,B:0.949):0.244,(C:7.8,D:1.439):2.87);")
    candidates = tmp_path / "candidates.nwk"
    if invalid != "missing_file":
        texts = {
            "empty_file": "",
            "duplicate_tips": "((A:1,A:1,B:1):1,(C:1,D:1):1);",
            "different_tips": "((A:1,B:1):1,(C:1,WrongTip:1):1);",
            "unresolved_root": "(A:1,B:1,C:1,D:1);",
        }
        candidates.write_text(texts[invalid])
    outputs = [tmp_path / name for name in ("selected.nwk", "roots.tsv", "roots.pdf")]
    for path in outputs:
        path.write_bytes(b"previous output")
    result = run(sys.executable, SUPPORT / "species_tree_guided_gene_tree_rooting.py",
                 "--in-tree", source, "--candidate-trees", candidates, "--out-tree", outputs[0],
                 "--comparison-table", outputs[1], "--comparison-plot", outputs[2])
    assert result.returncode != 0, result.stdout + result.stderr
    assert all(path.read_bytes() == b"previous output" for path in outputs)


def test_dated_tree_stage_keeps_pdf_pair_on_failed_install_and_tracks_engine(tmp_path):
    core = (ROOT / "workflow/core/gg_genome_evolution_core.sh").read_text()
    begin = core.index('task="Dated species tree plotting"')
    stage = core[begin:core.index('\n# Species taxonomy', begin)]
    source = tmp_path / "mcmctree_95CI.nhx"
    source.write_text("(A:10,B:10)Root[&&NHX:age=10:age_ci_low=8:age_ci_high=12:age_ci_kind=HPD:age_ci_level=0.95];")
    output = tmp_path / "output"
    script = tmp_path / "stage.sh"
    script.write_text(f'''set -euo pipefail
gg_support_dir={shlex.quote(str(SUPPORT))}
source "${{gg_support_dir}}/gg_util.sh"
gg_workspace_dir={shlex.quote(str(tmp_path))}
gg_workspace_output_dir={shlex.quote(str(output))}
genome_evolution_provenance_dir="${{gg_workspace_output_dir}}/provenance"
genome_nwkit_identity="${{1:-engine-one}}"
dir_mcmctree2={shlex.quote(str(tmp_path))}
file_mcmctree_dated_nwk={shlex.quote(str(tmp_path / "absent.nwk"))}
file_dated_species_tree={shlex.quote(str(source))}
file_plot_mcmctree_pdf="${{gg_workspace_output_dir}}/dated.pdf"
file_dated_species_tree_pdf="${{gg_workspace_output_dir}}/summary.pdf"
run_plot_mcmctreer=1
artifact_stale_policy=rebuild
gg_step_start() {{ echo run >> runs.txt; }}
gg_step_skip() {{ echo skip >> runs.txt; }}
mv() {{
  if [[ ${{2:-}} == *.gg-stage.* && ${{3:-}} == "${{file_dated_species_tree_pdf}}" && ${{FAIL_INSTALL:-0}} == 1 ]]; then
    return 1
  fi
  command mv "$@"
}}
{stage}
''')
    import os

    def execute(identity="engine-one", fail=False):
        return subprocess.run(["bash", str(script), identity], cwd=tmp_path,
                              env=dict(os.environ, FAIL_INSTALL="1" if fail else "0"), capture_output=True, text=True)

    for _ in range(2):
        result = execute()
        assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "runs.txt").read_text().splitlines() == ["run", "skip"]
    paths = [output / "dated.pdf", output / "summary.pdf"]
    old = [p.read_bytes() for p in paths]
    assert old[0] == old[1]
    failed = execute("engine-two", fail=True)
    assert failed.returncode != 0
    assert [p.read_bytes() for p in paths] == old
    result = execute("engine-two")
    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "runs.txt").read_text().splitlines() == ["run", "skip", "run", "run"]
    assert paths[0].read_bytes() == paths[1].read_bytes()


def test_root_adapter_consumes_real_nwkit_candidates(tmp_path):
    source = tmp_path / "gene.nwk"
    source.write_text("((A_a_1:1,B_b_1:1):1,(A_a_2:1,B_b_2:1):1);")
    species = tmp_path / "species.nwk"
    species.write_text("(A_a:1,B_b:1);")
    candidates = tmp_path / "candidates.nwk"
    result = run("nwkit", "root", "--method", "reconciliation", "--infile", source,
                 "--species-tree", species, "--outfile", tmp_path / "first.nwk",
                 "--candidates-out", candidates)
    assert result.returncode == 0, result.stderr
    output, table, plot = tmp_path / "selected.nwk", tmp_path / "roots.tsv", tmp_path / "roots.pdf"
    result = run(sys.executable, SUPPORT / "species_tree_guided_gene_tree_rooting.py",
                 "--in-tree", source, "--candidate-trees", candidates, "--out-tree", output,
                 "--comparison-table", table, "--comparison-plot", plot)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Selected root: mad_compatible_with_reconciliation" in result.stdout
    assert plot.read_bytes().startswith(b"%PDF-")

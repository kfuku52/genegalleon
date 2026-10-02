"""Container-only rendering and gene-summary integration for saved-family plots."""
import os
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas
import pytest

from workflow.tests.test_presence_absence_inputs import saved_manifest

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / "workflow/support"


def run_r(args, directory):
    assert shutil.which("Rscript"), "Run this test in a GeneGalleon container"
    return subprocess.run(["Rscript", str(SUPPORT / "plot_query2family_presence_absence.R"), *args],
                          cwd=directory, capture_output=True, text=True, check=True)


def summary_fixture(tmp_path):
    tree = tmp_path / "species.nwk"
    tree.write_text("(Species_one:1,Species_two:1);\n")
    families = [f"family_{i}" for i in range(70)]
    long_table = tmp_path / "long.tsv"
    pandas.DataFrame([
        dict(species=species, species_display=species.replace("_", " "), query=family,
             query_order=i, presence=1, copy_number=1, status="complete")
        for i, family in enumerate(families, 1) for species in ("Species_one", "Species_two")
    ]).to_csv(long_table, sep="\t", index=False)
    return [f"--species_tree={tree}", f"--long_table={long_table}"]


def test_empty_outputs_auto_width_label_map_and_focus_rows(tmp_path):
    args = summary_fixture(tmp_path)
    labels = tmp_path / "labels.tsv"
    display_name = "Readable_family_name_with_wide_WWWWWWWWWWWWWWWWWWWWW"
    labels.write_text(f"kind\tid\tlabel\nfamily\tfamily_0\t{display_name}\n")
    svg = tmp_path / "summary.svg"
    run_r([*args, "--width=auto", f"--label_map={labels}", "--focus_species=Species_one",
           "--out_pdf=", f"--out_svg={svg}"], tmp_path)
    assert not (tmp_path / "1").exists()
    root = ET.parse(svg).getroot()
    width, height = [float(value) for value in root.attrib["viewBox"].split()[2:]]
    assert width == pytest.approx((4.8 + 70 * 0.14) * 72, abs=0.02)
    text = svg.read_text()
    assert display_name in text
    assert ">family_0</text>" not in text
    label_boxes = [(float(y), float(length)) for y, length in re.findall(
        r"translate\([0-9.]+,([0-9.]+)\) rotate\(-90\).*?textLength='([0-9.]+)px'", text)]
    assert len(label_boxes) == 70
    assert max(y + length for y, length in label_boxes) < height - 3
    # One row outline, distinct from the matrix tiles and their white borders.
    assert any(element.tag.endswith("rect") and "stroke-width: 0.64" in element.attrib.get("style", "")
               and "stroke: none" not in element.attrib.get("style", "")
               and float(element.attrib.get("width", "0")) > width * 0.8 for element in root.iter())


def test_empty_svg_and_both_empty_outputs_do_not_write_sentinel_file(tmp_path):
    args = summary_fixture(tmp_path)
    pdf = tmp_path / "summary.pdf"
    run_r([*args, f"--out_pdf={pdf}", "--out_svg="], tmp_path)
    assert pdf.read_bytes().startswith(b"%PDF")
    result = run_r([*args, "--out_pdf=", "--out_svg="], tmp_path)
    assert "Wrote PDF" not in result.stdout and "Wrote SVG" not in result.stdout
    assert not (tmp_path / "1").exists()


@pytest.mark.parametrize("source", ["query2family", "orthogroup"])
def test_gene_summary_manifest_only_reuses_saved_trees_and_evidence(tmp_path, source):
    assert shutil.which("Rscript"), "Run this test in a GeneGalleon container"
    _store, manifest = saved_manifest(tmp_path)
    tree = tmp_path / "species.nwk"
    tree.write_text("((Anchor_one:1,Anchor_two:1):1,Missing_species:2);\n")
    workspace = tmp_path / "workspace"
    output = tmp_path / "summary"
    env = os.environ.copy()
    env.update(gg_support_dir=str(SUPPORT), gg_workspace_dir=str(workspace),
               gene_family_source=source, summary_output_dir=str(output),
               run_species_taxonomy="0", run_family_completion_summary="0",
               run_presence_absence_summary="1", presence_absence_species_tree=str(tree),
               presence_absence_species_tree_ci="", presence_absence_species_tree_support="",
               presence_absence_busco_table="", presence_absence_ortholog_basis="query_gene",
               presence_absence_family_manifest=str(manifest), presence_absence_query_label="label",
               presence_absence_focus_species="Anchor_one", presence_absence_plot_width="auto",
               presence_absence_legend_columns="3" if source == "orthogroup" else "auto")
    result = subprocess.run(["bash", str(ROOT / "workflow/core/gg_gene_summary_core.sh")],
                            env=env, cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    prefix = output / "query2family_query_gene_orthologs"
    assert Path(str(prefix) + ".pdf").read_bytes().startswith(b"%PDF")
    assert ">HOG1 Anchor_one_a</text>" in Path(str(prefix) + ".svg").read_text()
    columns = pandas.read_csv(str(prefix) + ".columns.tsv", sep="\t")
    assert len(columns) == 2 and set(columns.family_id) == {"Combined"}
    mapping = pandas.read_csv(str(prefix) + ".query_map.tsv", sep="\t")
    assert set(mapping.hog_ids) == {"HOG1", "HOG2"}
    assert Path(str(prefix) + ".selection.tsv").is_file()
    assert Path(str(prefix) + ".overlap.tsv").is_file()
    assert "Skipping gene-family presence/absence summary because no" not in result.stdout
    svg_path = Path(str(prefix) + ".svg")
    text = svg_path.read_text()
    label_boxes = [(float(y), float(length)) for y, length in re.findall(
        r"translate\([0-9.]+,([0-9.]+)\) rotate\(-90\).*?textLength='([0-9.]+)px'", text)]
    svg = ET.parse(svg_path).getroot()
    legend = next(element for element in svg.iter()
                  if element.tag.endswith("text") and element.text == "Query-gene orthologs")
    assert max(y + length for y, length in label_boxes) + 6 < float(legend.attrib["y"])
    if source == "orthogroup":
        legend_text = [element for element in svg.iter()
                       if element.tag.endswith("text") and element.text in (
                           "undetected", "query-anchor-specific", "shared across query anchors", "D#: mapped duplication")]
        assert len({element.attrib["x"] for element in legend_text}) > 1
        for left in legend_text:
            for right in legend_text:
                if left.attrib["y"] == right.attrib["y"] and float(left.attrib["x"]) < float(right.attrib["x"]):
                    assert float(left.attrib["x"]) + float(left.attrib["textLength"].removesuffix("px")) + 6 < float(right.attrib["x"])

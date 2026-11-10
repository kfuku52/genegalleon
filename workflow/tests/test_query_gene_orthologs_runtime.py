"""R rendering contracts run in the GeneGalleon runtime lane."""
import subprocess
import xml.etree.ElementTree as ET

import pandas
import pytest
from test_query_gene_orthologs import (
    SUPPORT_DIR,
    collect_weak_fixture,
)
from test_query_gene_orthologs import (
    weak_duplication_fixture as fixture_provider,
)

weak_duplication_fixture = fixture_provider


def test_strict_ortholog_copy_numbers_fit_between_evidence_bands(weak_duplication_fixture, tmp_path):
    fixture = weak_duplication_fixture
    out_dir = collect_weak_fixture(fixture, tmp_path / "strict_band", 0)
    command = ["Rscript", str(SUPPORT_DIR / "plot_query2family_presence_absence.R"),
               f"--species_tree={fixture['species_tree']}", f"--species_mapping_tree={fixture['species_tree']}",
               f"--long_table={fixture['long_table']}", "--reference_species=Reference_species",
               "--width=7.2", "--evidence_layout=band", f"--out_svg={out_dir / 'strict.svg'}"]
    command += [f"--ortholog_{kind}_table={out_dir / (name + '.tsv')}" for kind, name in
                (("column", "columns"), ("glyph", "glyphs"), ("tree", "tree"), ("synteny", "synteny"), ("ufboot", "ufboot"))]
    subprocess.run(command, check=True, capture_output=True, text=True)
    elements = list(ET.parse(out_dir / "strict.svg").getroot().iter())
    labels = [element for element in elements if element.tag.endswith("text") and element.text == "1"]
    cells = [element for element in elements if element.tag.endswith("rect") and
             "fill: #2166AC" in element.attrib.get("style", "") and
             any(float(element.attrib["x"]) <= float(label.attrib.get("x", -1)) <=
                 float(element.attrib["x"]) + float(element.attrib["width"]) and
                 float(element.attrib["y"]) <= float(label.attrib.get("y", -1)) <=
                 float(element.attrib["y"]) + float(element.attrib["height"])
                 for label in labels)]
    assert cells
    assert min(float(cell.attrib["height"]) * 0.64 for cell in cells) >= 8
    assert "Additional ortholog candidate" not in (out_dir / "strict.svg").read_text()
    too_short = subprocess.run([*command, "--height=1"], check=False, capture_output=True, text=True)
    assert too_short.returncode != 0
    assert "height is too small for copy-number labels" in too_short.stderr


def test_fill_rows_and_loss_marks_with_a_noncontiguous_glyph(tmp_path):
    tree = tmp_path / "species.nwk"
    tree.write_text("(Species_one:1,Species_two:1);\n")
    long_table = tmp_path / "long.tsv"
    columns = tmp_path / "columns.tsv"
    glyphs = tmp_path / "glyphs.tsv"
    losses = tmp_path / "losses.tsv"
    basis = "query_gene"
    families = {"F1": ("a", "b", "c"), "F2": ("d", "e")}
    pandas.DataFrame([
        dict(species=species, species_display=species.replace("_", " "), query=family,
             query_order=order, presence=1, copy_number=1, status="complete")
        for order, family in enumerate(families, 1) for species in ("Species_one", "Species_two")
    ]).to_csv(long_table, sep="\t", index=False)
    pandas.DataFrame([
        dict(column_order=column, family_id=family, family_order=1 if family == "F1" else 2,
             basis=basis, query_ids=anchor, query_count=1, anchor_species="Species_one",
             anchor_cds_fasta_id=f"Species_one_{anchor}", anchor_gene_id=anchor,
             plot_label=anchor, anchor_tip_branch_id=column)
        for column, (family, anchor) in enumerate(
            ((family, anchor) for family, anchors in families.items() for anchor in anchors), 1
        )
    ]).to_csv(columns, sep="\t", index=False)
    pandas.DataFrame([
        dict(species="Species_two", family_id=family, family_order=1 if family == "F1" else 2,
             basis=basis, relation=relation, anchor_cds_fasta_ids=";".join(f"Species_one_{a}" for a in anchors),
             anchor_gene_ids=";".join(anchors), anchor_count=len(anchors), copy_number=copies,
             start_order=start, end_order=end, lane_index=index, lane_count=count)
        for family, relation, anchors, copies, start, end, index, count in (
            ("F1", "ambiguous", ("a", "c"), 3, 1, 3, 1, 1),
            ("F2", "specific", ("d",), 2, 4, 4, 1, 2),
            ("F2", "shared_ancestral", ("d", "e"), 1, 4, 5, 2, 2),
        )
    ]).to_csv(glyphs, sep="\t", index=False)
    losses.write_text(
        "family_id\tspecies\tanchor_cds_fasta_id\n"
        "F1\tSpecies_two\tSpecies_one_a\n"
        "F1\tSpecies_two\tSpecies_one_c\n"
        "F1\tSpecies_two\tSpecies_one_b\n"
    )
    (tmp_path / "species_label_utils.r").symlink_to(SUPPORT_DIR / "species_label_utils.r")
    wrapper = tmp_path / "inspect.R"
    wrapper.write_text(
        f'source("{SUPPORT_DIR / "plot_query2family_presence_absence.R"}", local=TRUE)\n'
        'write.table(glyph_rect_df, "geometry.tsv", sep="\\t", quote=FALSE, row.names=FALSE)\n'
        'write.table(loss_label_df, "loss-geometry.tsv", sep="\\t", quote=FALSE, row.names=FALSE)\n'
    )
    svg = tmp_path / "view.svg"
    cmd = ["Rscript", str(wrapper), f"--species_tree={tree}", f"--long_table={long_table}",
           "--ortholog_basis=query_gene", f"--ortholog_column_table={columns}", f"--ortholog_glyph_table={glyphs}",
           "--evidence_layout=off", f"--out_svg={svg}"]
    subprocess.run(cmd, cwd=tmp_path, check=True, capture_output=True, text=True)
    compact = pandas.read_csv(tmp_path / "geometry.tsv", sep="\t")
    subprocess.run([*cmd, "--glyph_height=fill", f"--loss_table={losses}"],
                   cwd=tmp_path, check=True, capture_output=True, text=True)
    filled = pandas.read_csv(tmp_path / "geometry.tsv", sep="\t")
    loss = pandas.read_csv(tmp_path / "loss-geometry.tsv", sep="\t")
    single = filled[filled.family_id == "F1"].iloc[0]
    stacked = filled[filled.family_id == "F2"]
    assert single.ymin == pytest.approx(stacked.ymin.min())
    assert single.ymax == pytest.approx(stacked.ymax.max())
    assert single.ymax - single.ymin > (compact.ymax - compact.ymin).max()
    assert single.text_y == pytest.approx((single.ymin + single.ymax) / 2)
    assert single.copy_label == "3†"
    assert loss.y.item() < single.ymin
    assert ">0†</text>" in svg.read_text()


@pytest.mark.parametrize("case", ["strict", "candidates", "subset_species", "legacy", "partial", "wrong_genes"])
def test_duplication_bars_follow_displayed_gene_membership(weak_duplication_fixture, tmp_path, case):
    fixture = weak_duplication_fixture
    threshold = 0 if case == "strict" else 0.05
    out_dir = collect_weak_fixture(fixture, tmp_path / "bars", threshold, "query_gene")
    tree_path = out_dir / "tree.tsv"
    tree = pandas.read_csv(tree_path, sep="\t", keep_default_na=False)
    metadata = ["displayed_gene_ids", "displayed_child1_gene_ids", "displayed_child2_gene_ids"]
    if case == "legacy":
        tree.drop(columns=metadata).to_csv(tree_path, sep="\t", index=False)
    elif case == "partial":
        tree.drop(columns=metadata[:1]).to_csv(tree_path, sep="\t", index=False)
    elif case == "wrong_genes":
        tree.loc[tree.displayed_gene_ids != "", "displayed_gene_ids"] += ";Missing_gene"
        tree.to_csv(tree_path, sep="\t", index=False)
    long_table = fixture["long_table"]
    if case == "subset_species":
        long_table = out_dir / "subset.long.tsv"
        original = pandas.read_csv(fixture["long_table"], sep="\t")
        original.loc[original.species != "Shared_species"].to_csv(long_table, sep="\t", index=False)
    wrapper = out_dir / "export.R"
    wrapper.write_text(
        f'source("{SUPPORT_DIR / "plot_query2family_presence_absence.R"}", local=TRUE)\n'
        f'write.table(duplication_bar_nodes, "{out_dir / "bar_nodes.tsv"}", sep="\\t", quote=FALSE, row.names=FALSE)\n'
    )
    (out_dir / "species_label_utils.r").symlink_to(SUPPORT_DIR / "species_label_utils.r")
    result = subprocess.run([
        "Rscript", str(wrapper), f"--species_tree={fixture['species_tree']}",
        f"--species_mapping_tree={fixture['species_tree']}", f"--long_table={long_table}",
        "--ortholog_basis=query_gene", f"--dup_conf_score_threshold={threshold}",
        f"--ortholog_column_table={out_dir / 'columns.tsv'}", f"--ortholog_glyph_table={out_dir / 'glyphs.tsv'}",
        f"--ortholog_tree_table={tree_path}", f"--ortholog_dup_conf_table={out_dir / 'dup_conf.tsv'}",
        "--evidence_layout=off", f"--out_svg={out_dir / 'bars.svg'}",
    ], check=False, capture_output=True, text=True)
    if case in ("partial", "wrong_genes"):
        assert result.returncode != 0
        assert ("incomplete" if case == "partial" else "disagree") in result.stderr
        return
    assert result.returncode == 0, result.stdout + result.stderr
    bars = pandas.read_csv(out_dir / "bar_nodes.tsv", sep="\t")
    assert (fixture["weak_branch"] in set(bars.node_id)) == (case in ("candidates", "legacy"))
    svg = (out_dir / "bars.svg").read_text()
    assert ("Bar height = full-family duplication count" if case == "legacy" else
            "Bar height = displayed-gene duplication count") in svg
    if case == "legacy":
        assert "Legacy ortholog tree table" in result.stderr
    assert len(set(bars.node_id)) == len(bars)

import gzip
import re
import shutil
import subprocess
import xml.etree.ElementTree as ElementTree
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace

import pandas
import pytest

from workflow.support.gene_family_output_store import (
    archive_completed_outputs,
    family_context,
)

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"


def svg_number(element: ElementTree.Element, attribute: str, default: float) -> float:
    try:
        return float(element.attrib[attribute])
    except (KeyError, ValueError):
        return default


def count_svg_matrix_edge_bands(svg_path: Path, species_label: str, column_count: int) -> int:
    root = ElementTree.parse(svg_path).getroot()
    species_baselines = [
        float(element.attrib["y"])
        for element in root.iter()
        if element.tag.endswith("text")
        and (element.text or "") == species_label
        and "y" in element.attrib
    ]
    assert len(species_baselines) == 1
    species_baseline = species_baselines[0]

    rectangles = []
    for element in root.iter():
        if not element.tag.endswith("rect"):
            continue
        try:
            rectangle = {
                "x": float(element.attrib["x"]),
                "y": float(element.attrib["y"]),
                "width": float(element.attrib["width"]),
                "height": float(element.attrib["height"]),
                "style": element.attrib.get("style", ""),
            }
        except (KeyError, ValueError):
            continue
        rectangles.append(rectangle)

    row_tiles = [
        rectangle
        for rectangle in rectangles
        if "stroke: #D9D9D9" in rectangle["style"]
        and rectangle["y"] <= species_baseline <= rectangle["y"] + rectangle["height"]
    ]
    assert len(row_tiles) == column_count
    row_ymin = min(rectangle["y"] for rectangle in row_tiles)
    row_ymax = max(rectangle["y"] + rectangle["height"] for rectangle in row_tiles)
    matrix_xmin = min(rectangle["x"] for rectangle in row_tiles)
    matrix_xmax = max(rectangle["x"] + rectangle["width"] for rectangle in row_tiles)
    tile_height = min(rectangle["height"] for rectangle in row_tiles)
    tolerance = 0.05
    return sum(
        rectangle["height"] < tile_height * 0.3
        and rectangle["x"] >= matrix_xmin - tolerance
        and rectangle["x"] + rectangle["width"] <= matrix_xmax + tolerance
        and rectangle["y"] >= row_ymin - tolerance
        and rectangle["y"] + rectangle["height"] <= row_ymax + tolerance
        for rectangle in rectangles
    )


def assert_svg_copy_numbers_above_edge_bands(svg_path: Path, species_label: str) -> None:
    root = ElementTree.parse(svg_path).getroot()
    elements = list(root.iter())
    element_order = {id(element): index for index, element in enumerate(elements)}
    species_baselines = [
        float(element.attrib["y"])
        for element in elements
        if element.tag.endswith("text")
        and (element.text or "") == species_label
        and "y" in element.attrib
    ]
    assert len(species_baselines) == 1
    species_baseline = species_baselines[0]

    row_tiles = [
        element
        for element in elements
        if element.tag.endswith("rect")
        and "stroke: #D9D9D9" in element.attrib.get("style", "")
        and svg_number(element, "y", float("inf"))
        <= species_baseline
        <= svg_number(element, "y", float("-inf"))
        + svg_number(element, "height", 0.0)
    ]
    assert row_tiles
    row_ymin = min(svg_number(element, "y", float("inf")) for element in row_tiles)
    row_ymax = max(
        svg_number(element, "y", float("-inf"))
        + svg_number(element, "height", 0.0)
        for element in row_tiles
    )
    matrix_xmin = min(svg_number(element, "x", float("inf")) for element in row_tiles)
    matrix_xmax = max(
        svg_number(element, "x", float("-inf"))
        + svg_number(element, "width", 0.0)
        for element in row_tiles
    )
    tile_height = min(svg_number(element, "height", float("inf")) for element in row_tiles)
    edge_bands = [
        element
        for element in elements
        if element.tag.endswith("rect")
        and svg_number(element, "height", float("inf")) < tile_height * 0.3
        and matrix_xmin <= svg_number(element, "x", float("-inf")) <= matrix_xmax
        and row_ymin
        <= svg_number(element, "y", float("-inf"))
        <= row_ymax
    ]
    assert edge_bands

    copy_numbers = [
        element
        for element in elements
        if element.tag.endswith("text")
        and re.fullmatch(r"[0-9]+", element.text or "")
        and matrix_xmin <= svg_number(element, "x", float("-inf")) <= matrix_xmax
        and row_ymin - 1.0
        <= svg_number(element, "y", float("-inf"))
        <= row_ymax + 1.0
    ]
    assert copy_numbers
    assert min(element_order[id(element)] for element in copy_numbers) > max(
        element_order[id(element)] for element in edge_bands
    )


def load_module(name: str):
    path = SUPPORT_DIR / name
    spec = spec_from_file_location(path.stem, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def stat_row(
    branch_id,
    parent,
    child1,
    child2,
    event,
    name,
    species="",
    marker="",
    support="",
    generax_support="",
):
    return {
        "branch_id": branch_id,
        "parent": parent,
        "child1": child1,
        "child2": child2,
        "so_event": event,
        "node_name": name,
        "spnode_coverage": species,
        "query_marker_source": marker,
        "support_unrooted": support,
        "support_generax_ufboot": generax_support,
    }


def write_family(query_dir, output_root, family_id, query_records, rows, cds_ids, gzip_cds=False):
    query_text = "".join(
        f">{query_id} description | {query_name} | focal species\nAAAA\n"
        for query_id, query_name in query_records
    )
    (query_dir / family_id).write_text(query_text, encoding="utf-8")
    stat_dir = output_root / "stat_branch"
    cds_dir = output_root / "cds_fasta"
    stat_dir.mkdir(parents=True, exist_ok=True)
    cds_dir.mkdir(parents=True, exist_ok=True)
    pandas.DataFrame(rows).to_csv(stat_dir / f"{family_id}_stat.branch.tsv", sep="\t", index=False)
    cds_text = "".join(f">{cds_id}\nATG\n" for cds_id in cds_ids)
    if gzip_cds:
        with gzip.open(cds_dir / f"{family_id}_cds.fa.gz", "wt", encoding="utf-8") as handle:
            handle.write(cds_text)
    else:
        (cds_dir / f"{family_id}_cds.fasta").write_text(cds_text, encoding="utf-8")


def write_query_only_family(query_dir, output_root, family_id, species, query_count):
    leaf_ids = list(range(query_count))
    rows_by_id = {
        leaf_id: stat_row(
            leaf_id,
            -1,
            -1,
            -1,
            "L",
            f"{species}_GENE_{family_id}_{leaf_id + 1}",
            species,
            f"direct:{family_id}_q{leaf_id + 1}",
        )
        for leaf_id in leaf_ids
    }
    current_root = leaf_ids[0]
    next_node_id = query_count
    for leaf_id in leaf_ids[1:]:
        rows_by_id[current_root]["parent"] = next_node_id
        rows_by_id[leaf_id]["parent"] = next_node_id
        rows_by_id[next_node_id] = stat_row(
            next_node_id,
            -1,
            current_root,
            leaf_id,
            "D",
            f"{family_id}_duplication_{next_node_id}",
            species,
        )
        current_root = next_node_id
        next_node_id += 1
    query_records = [
        (f"{family_id}_q{index + 1}", f"Friendly_{family_id}_{index + 1}")
        for index in leaf_ids
    ]
    cds_ids = [f"{species}_GENE_{family_id}_{index + 1}" for index in leaf_ids]
    write_family(
        query_dir,
        output_root,
        family_id,
        query_records,
        [rows_by_id[node_id] for node_id in sorted(rows_by_id)],
        cds_ids,
        gzip_cds=query_count % 2 == 0,
    )


def create_multifamily_fixture(tmp_path: Path):
    query_dir = tmp_path / "input" / "query_gene"
    output_root = tmp_path / "output" / "query2family"
    query_dir.mkdir(parents=True)

    write_family(
        query_dir,
        output_root,
        "TWO",
        [("qC", "FriendlyC"), ("qB", "FriendlyB"), ("qA", "FriendlyA")],
        [
            stat_row(0, -1, 1, 2, "D", "root", "n0"),
            stat_row(1, 0, 3, 4, "D", "AB_duplication", "Reference_species"),
            stat_row(2, 0, -1, -1, "L", "Reference_species_GENE_C", "Reference_species", "direct:qC"),
            stat_row(3, 1, -1, -1, "L", "Reference_species_GENE_A", "Reference_species", "direct:qA"),
            stat_row(4, 1, -1, -1, "L", "Reference_species_GENE_B", "Reference_species", "direct:qB"),
        ],
        ["Reference_species_GENE_A", "Reference_species_GENE_B", "Reference_species_GENE_C"],
        gzip_cds=True,
    )
    write_family(
        query_dir,
        output_root,
        "ONE",
        [("qY", "FriendlyY"), ("qX", "FriendlyX")],
        [
            stat_row(0, -1, 1, 2, "S", "root"),
            stat_row(1, 0, 5, 6, "D", "ancestral_duplication_1", "Ancestor_species"),
            stat_row(2, 0, 3, 4, "D", "XY_duplication", "Reference_species"),
            stat_row(3, 2, 10, 11, "S", "X_speciation_outer", support=87),
            stat_row(4, 2, -1, -1, "L", "Reference_species_GENE_Y", "Reference_species", "direct:qY"),
            stat_row(5, 1, 7, 8, "D", "ancestral_duplication_2", "Ancestor_species"),
            stat_row(6, 1, -1, -1, "L", "Ancestor_species_copy_1", "Ancestor_species"),
            stat_row(7, 5, -1, -1, "L", "Ancestor_species_copy_2", "Ancestor_species"),
            stat_row(8, 5, -1, -1, "L", "Ancestor_species_copy_3", "Ancestor_species"),
            stat_row(9, 11, -1, -1, "L", "Reference_species_GENE_X", "Reference_species", "direct:qX"),
            stat_row(10, 3, -1, -1, "L", "Ancestor_species_X_sister", "Ancestor_species"),
            stat_row(11, 3, 9, 12, "S", "X_speciation_inner", support=96),
            stat_row(12, 11, -1, -1, "L", "Ancestor_species_X_sister_2", "Ancestor_species"),
        ],
        ["Reference_species_GENE_X", "Reference_species_GENE_Y"],
    )
    write_family(
        query_dir,
        output_root,
        "SOLO",
        [("qSolo", "FriendlySolo")],
        [
            stat_row(
                0,
                -1,
                -1,
                -1,
                "L",
                "Reference_species_LONG_GENE_ID_001",
                "Reference_species",
                "direct:qSolo",
            )
        ],
        ["Reference_species_LONG_GENE_ID_001"],
    )
    selection = tmp_path / "selection.tsv"
    selection.write_text("family_id\nTWO\nONE\nSOLO\n", encoding="utf-8")
    return query_dir, output_root, selection


def run_multifamily_collection(tmp_path: Path):
    mod = load_module("query_gene_orthologs.py")
    query_dir, output_root, selection = create_multifamily_fixture(tmp_path)
    synteny_dir = output_root / "synteny"
    synteny_dir.mkdir(parents=True, exist_ok=True)
    pandas.DataFrame(
        [
            ["Reference_species_GENE_X", "Reference_species", "upstream", -2, "rx1", "GX1", 3],
            ["Reference_species_GENE_X", "Reference_species", "downstream", 2, "rx2", "GX2", 2],
            ["Reference_species_GENE_Y", "Reference_species", "upstream", -2, "ry1", "GY1", 2],
            ["Reference_species_GENE_Y", "Reference_species", "downstream", 2, "ry2", "GY2", 2],
            ["Ancestor_species_copy_1", "Ancestor_species", "upstream", -1, "a1x", "GX1", 3],
            ["Ancestor_species_copy_1", "Ancestor_species", "downstream", 1, "a1x2", "GX2", 2],
            ["Ancestor_species_copy_1", "Ancestor_species", "downstream", 2, "a1y", "GY1", 2],
            ["Ancestor_species_copy_1", "Ancestor_species", "downstream", 3, "a1y2", "GY2", 2],
            ["Ancestor_species_copy_2", "Ancestor_species", "upstream", -1, "a2", "GX1", 3],
            ["Ancestor_species_copy_3", "Ancestor_species", "upstream", -1, "a3", "UNIQUE", 1],
        ],
        columns=[
            "node_name", "species", "direction", "offset", "neighbor_gene", "group_id", "group_size"
        ],
    ).to_csv(synteny_dir / "ONE_synteny.tsv", sep="\t", index=False)
    out_dir = tmp_path / "summary"
    args = SimpleNamespace(
        dir_gene_family=str(output_root),
        dir_query_gene=str(query_dir),
        family_file=str(selection),
        reference_species="Reference_species",
        out_columns=str(out_dir / "columns.tsv"),
        out_glyphs=str(out_dir / "glyphs.tsv"),
        out_tree=str(out_dir / "tree.tsv"),
        out_synteny=str(out_dir / "synteny.tsv"),
        out_ufboot=str(out_dir / "ufboot.tsv"),
    )
    mod.run(args)
    return out_dir


def test_multifamily_columns_glyphs_and_trees_remain_family_scoped(tmp_path: Path):
    out_dir = run_multifamily_collection(tmp_path)
    columns = pandas.read_csv(out_dir / "columns.tsv", sep="\t")
    glyphs = pandas.read_csv(out_dir / "glyphs.tsv", sep="\t")
    tree = pandas.read_csv(out_dir / "tree.tsv", sep="\t")
    synteny = pandas.read_csv(out_dir / "synteny.tsv", sep="\t")
    ufboot = pandas.read_csv(out_dir / "ufboot.tsv", sep="\t")

    assert columns["family_id"].tolist() == ["TWO", "TWO", "TWO", "ONE", "ONE", "SOLO"]
    assert columns["gene_id"].tolist() == [
        "GENE_A",
        "GENE_B",
        "GENE_C",
        "GENE_X",
        "GENE_Y",
        "LONG_GENE_ID_001",
    ]
    assert set(columns["reference_species"]) == {"Reference_species"}
    assert columns["plot_label"].tolist() == columns["gene_id"].tolist()
    assert columns["column_order"].tolist() == list(range(1, 7))

    ranges = {
        family_id: (int(group["column_order"].min()), int(group["column_order"].max()))
        for family_id, group in columns.groupby("family_id", sort=False)
    }
    for row in glyphs.itertuples(index=False):
        family_min, family_max = ranges[row.family_id]
        assert family_min <= row.start_order <= row.end_order <= family_max
    ancestral = glyphs.loc[
        (glyphs["family_id"] == "ONE") &
        (glyphs["species"] == "Ancestor_species") &
        (glyphs["relation"] == "shared_ancestral")
    ].iloc[0]
    assert ancestral["relation"] == "shared_ancestral"
    assert (ancestral["start_order"], ancestral["end_order"]) == (4, 5)

    assert tree.groupby("family_id")["is_tip"].sum().to_dict() == {"ONE": 2, "SOLO": 1, "TWO": 3}
    assert len(set(zip(tree["family_id"], tree["node_id"], strict=True))) == len(tree)
    assert tree["node_id"].duplicated().any(), "Fixture must exercise node-ID collisions across families"
    solo = tree.loc[tree["family_id"] == "SOLO"].iloc[0]
    assert solo["gene_id"] == "LONG_GENE_ID_001"
    assert solo["node_height"] == 0
    two_duplications = tree.loc[(tree["family_id"] == "TWO") & (tree["event"] == "D")]
    assert two_duplications["mapped_species_node"].tolist() == ["n0", "Reference_species"]
    assert sorted(two_duplications["duplication_index"].astype(int).tolist()) == [1, 2]
    one_duplications = tree.loc[(tree["family_id"] == "ONE") & (tree["event"] == "D")]
    assert one_duplications["mapped_species_node"].tolist() == ["Ancestor_species", "Reference_species", "Ancestor_species"]
    assert one_duplications["duplication_index"].astype(int).tolist() == [1, 2, 3]
    assert one_duplications["in_reference_tree"].astype(int).tolist() == [0, 1, 0]
    one_ancestor_synteny = synteny.loc[
        (synteny["family_id"] == "ONE") &
        (synteny["species"] == "Ancestor_species")
    ]
    assert set(one_ancestor_synteny["synteny_status"]) >= {
        "supported", "single_anchor", "no_support"
    }
    evaluated_ufboot = ufboot.loc[ufboot["orthology_ufboot_status"] == "evaluated"]
    assert evaluated_ufboot.shape[0] == 2
    assert set(evaluated_ufboot["candidate_cds_fasta_id"]) == {
        "Ancestor_species_X_sister",
        "Ancestor_species_X_sister_2",
    }
    assert set(evaluated_ufboot["reference_cds_fasta_id"]) == {
        "Reference_species_GENE_X"
    }
    assert set(evaluated_ufboot["orthology_mrca_event"]) == {"S"}
    assert set(evaluated_ufboot["ufboot_support_source"]) == {"support_unrooted"}
    assert sorted(evaluated_ufboot["decisive_branch_ufboot"].tolist()) == [87, 96]
    assert "mrca_is_root" in set(
        ufboot.loc[
            ufboot["orthology_ufboot_status"] == "not_evaluable",
            "orthology_ufboot_unavailable_reason",
        ]
    )


def test_multifamily_tables_render_to_pdf_and_svg(tmp_path: Path):
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript is unavailable")
    out_dir = run_multifamily_collection(tmp_path)
    tree_path = tmp_path / "species.nwk"
    tree_path.write_text(
        "(Ancestor_species:1,Reference_species:0)dated_root;\n",
        encoding="utf-8",
    )
    mapping_tree_path = tmp_path / "species_mapping.nwk"
    mapping_tree_path.write_text(
        "(Ancestor_species:1,Reference_species:1)n0;\n",
        encoding="utf-8",
    )
    glyph_table = pandas.read_csv(out_dir / "glyphs.tsv", sep="\t")
    glyph_table.loc[glyph_table.index[0], "relation"] = "ambiguous"
    glyph_table.to_csv(out_dir / "glyphs.tsv", sep="\t", index=False)
    species = ["Ancestor_species", "Reference_species"]
    long_table = tmp_path / "long.tsv"
    pandas.DataFrame(
        {
            "species": species,
            "species_display": [value.replace("_", " ") for value in species],
            "query": ["ONE", "TWO"],
            "query_order": [1, 2],
            "presence": [1, 1],
            "copy_number": [1, 1],
            "status": ["complete"] * 2,
        }
    ).to_csv(long_table, sep="\t", index=False)
    pdf_path = out_dir / "multifamily.pdf"
    svg_path = out_dir / "multifamily.svg"
    plot_command = [
            "Rscript",
                str(SUPPORT_DIR / "plot_query2family_presence_absence.R"),
                f"--species_tree={tree_path}",
                f"--species_mapping_tree={mapping_tree_path}",
                f"--long_table={long_table}",
            f"--ortholog_column_table={out_dir / 'columns.tsv'}",
            f"--ortholog_glyph_table={out_dir / 'glyphs.tsv'}",
            f"--ortholog_tree_table={out_dir / 'tree.tsv'}",
            f"--ortholog_synteny_table={out_dir / 'synteny.tsv'}",
            f"--ortholog_ufboot_table={out_dir / 'ufboot.tsv'}",
            "--reference_species=Reference_species",
            f"--out_pdf={pdf_path}",
            f"--out_svg={svg_path}",
        ]
    subprocess.run(
        plot_command,
        check=True,
        capture_output=True,
        text=True,
    )

    assert pdf_path.stat().st_size > 1000
    svg_text = svg_path.read_text(encoding="utf-8")
    for expected in ("TWO", "ONE", "SOLO", "GENE_A", "GENE_X", "LONG_GENE_ID_001"):
        assert expected in svg_text
    assert "FriendlyA" not in svg_text
    assert "Reference species orthologs" in svg_text
    assert "Reference species genes" not in svg_text
    assert "Duplication bar colors" in svg_text
    assert "Bar height = duplication count" in svg_text
    assert "non-contiguous orthology" in svg_text
    assert "Evidence band (top): local synteny anchors" in svg_text
    assert "Local synteny anchors (A)" not in svg_text
    assert "Fill: highest per-copy anchor count (0 = none)" in svg_text
    assert "Color: highest per-copy A in cell" not in svg_text
    assert "Window: 3 upstream + 3 downstream gene models" in svg_text
    assert "Evidence band (bottom): Gene tree UFBoot (%)" in svg_text
    assert "Fill: lowest per-copy MRCA-branch support" in svg_text
    assert "Numeric fill requires support for every copy" in svg_text
    assert "Evidence band states" in svg_text
    assert "unavailable" in svg_text
    assert "reference self" in svg_text
    assert "No band: reference self" not in svg_text
    assert "No circle: unavailable or reference self" not in svg_text
    assert count_svg_matrix_edge_bands(
        svg_path,
        "Reference species",
        column_count=6,
    ) == 12
    assert count_svg_matrix_edge_bands(
        svg_path,
        "Ancestor species",
        column_count=6,
    ) > 0
    assert_svg_copy_numbers_above_edge_bands(svg_path, "Reference species")
    assert_svg_copy_numbers_above_edge_bands(svg_path, "Ancestor species")
    glyph_svg_path = out_dir / "multifamily.glyph.svg"
    glyph_command = [
        argument for argument in plot_command if not argument.startswith("--out_")
    ] + ["--evidence_layout=glyph", f"--out_svg={glyph_svg_path}"]
    subprocess.run(
        glyph_command,
        check=True,
        capture_output=True,
        text=True,
    )
    glyph_svg_text = glyph_svg_path.read_text(encoding="utf-8")
    assert "Local synteny anchors" in glyph_svg_text
    assert "Color: highest per-copy anchor count" in glyph_svg_text
    assert "No circle: unavailable or reference self" in glyph_svg_text
    assert "Gene tree UFBoot (%)" in glyph_svg_text
    assert "Evidence band (top):" not in glyph_svg_text

    rail_svg_path = out_dir / "multifamily.rail.svg"
    rail_command = [
        argument for argument in plot_command if not argument.startswith("--out_")
    ] + ["--evidence_layout=rail", f"--out_svg={rail_svg_path}"]
    subprocess.run(
        rail_command,
        check=True,
        capture_output=True,
        text=True,
    )
    rail_svg_text = rail_svg_path.read_text(encoding="utf-8")
    assert "Evidence rail (top): local synteny anchors" in rail_svg_text
    assert "Evidence rail (bottom): Gene tree UFBoot (%)" in rail_svg_text
    assert "Evidence rail states" in rail_svg_text
    assert "reference self" in rail_svg_text
    assert "No band: reference self" not in rail_svg_text

    off_svg_path = out_dir / "multifamily.off.svg"
    off_command = [
        argument for argument in plot_command if not argument.startswith("--out_")
    ] + ["--evidence_layout=off", f"--out_svg={off_svg_path}"]
    subprocess.run(
        off_command,
        check=True,
        capture_output=True,
        text=True,
    )
    off_svg_text = off_svg_path.read_text(encoding="utf-8")
    assert "Local synteny anchors" not in off_svg_text
    assert "Gene tree UFBoot" not in off_svg_text
    assert "Evidence band states" not in off_svg_text
    expected_min_ufboot_color = subprocess.run(
        [
            "Rscript",
            "-e",
            'cat(grDevices::hcl.colors(256, palette="viridis")[[1 + round(87 / 100 * 255)]])',
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    assert svg_text.count(f"fill: {expected_min_ufboot_color};") == 1
    expected_zero_support_color = subprocess.run(
        [
            "Rscript",
            "-e",
            'cat(grDevices::hcl.colors(256, palette="viridis")[[1]])',
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    assert svg_text.count(f"fill: {expected_zero_support_color};") >= 2
    expected_glyph_ufboot_color = subprocess.run(
        [
            "Rscript",
            "-e",
            'cat(grDevices::hcl.colors(256, palette="Inferno")[[1 + round(87 / 100 * 255)]])',
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    assert glyph_svg_text.count(f"fill: {expected_glyph_ufboot_color};") == 1
    assert "Local synteny (n = copy count)" not in svg_text
    assert "single anchor only" not in svg_text
    assert "Circle area = duplication count" not in svg_text
    bar_title_y = float(
        re.search(
            r"<text x='[0-9.]+' y='([0-9.]+)'.*?>Bar height = duplication count</text>",
            svg_text,
        ).group(1)
    )
    family_legend_y = [
        float(value)
        for family_id in ("TWO", "ONE", "SOLO")
        for value in re.findall(
            rf"<text x='[0-9.]+' y='([0-9.]+)'.*?>{family_id}</text>",
            svg_text,
        )
    ]
    assert family_legend_y
    assert bar_title_y > max(family_legend_y) + 8
    assert "D#: mapped duplication" not in svg_text
    assert not re.search(r">D[0-9]+</text>", svg_text)
    assert "Branch length" in svg_text
    assert "Million years ago" not in svg_text
    assert re.search(r">1</text>", svg_text)

    synteny_table = pandas.read_csv(out_dir / "synteny.tsv", sep="\t")
    is_reference_self = (
        synteny_table["candidate_cds_fasta_id"]
        == synteny_table["reference_cds_fasta_id"]
    )
    invalid_status_cases = [
        ("nonself_as_self", synteny_table.index[~is_reference_self][0], "reference_self"),
        ("self_as_nonself", synteny_table.index[is_reference_self][0], "no_support"),
    ]
    for case_name, row_index, invalid_status in invalid_status_cases:
        invalid_synteny = synteny_table.copy()
        invalid_synteny.loc[row_index, "synteny_status"] = invalid_status
        invalid_synteny_path = out_dir / f"synteny.{case_name}.tsv"
        invalid_synteny.to_csv(invalid_synteny_path, sep="\t", index=False)
        invalid_result = subprocess.run(
            [
                argument
                for argument in plot_command
                if not argument.startswith("--out_")
                and not argument.startswith("--ortholog_synteny_table=")
            ]
            + [
                f"--ortholog_synteny_table={invalid_synteny_path}",
                f"--out_svg={out_dir / f'invalid_synteny.{case_name}.svg'}",
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        assert invalid_result.returncode != 0
        assert (
            "Ortholog synteny reference_self status disagrees with "
            "candidate/reference IDs"
        ) in (invalid_result.stdout + invalid_result.stderr)

    tree_table = pandas.read_csv(out_dir / "tree.tsv", sep="\t")
    mapped_duplications = tree_table.loc[
        (tree_table["event"] == "D") & tree_table["mapped_species_node"].notna()
    ]
    expected_bar_count = mapped_duplications.groupby(
        ["family_id", "mapped_species_node"], sort=False
    ).ngroups
    colored_rectangles = [
        tuple(float(value) for value in match[:4])
        for match in re.findall(
            r"<rect x='([0-9.]+)' y='([0-9.]+)' width='([0-9.]+)' height='([0-9.]+)' "
            r"style='stroke-width: 0.34;[^']*fill: #[0-9A-F]{6};' />",
            svg_text,
        )[:expected_bar_count]
    ]
    assert len(colored_rectangles) == expected_bar_count
    rectangles_by_baseline = {}
    for x, y, width, height in colored_rectangles:
        rectangles_by_baseline.setdefault(round(y + height, 1), []).append((x, width))
    same_branch_groups = [group for group in rectangles_by_baseline.values() if len(group) > 1]
    assert same_branch_groups, "Fixture must place multiple family bars on one zero-length branch"
    for group in same_branch_groups:
        ordered = sorted(group)
        assert all(
            left_x + left_width <= right_x + 0.02
            for (left_x, left_width), (right_x, _) in zip(
                ordered, ordered[1:], strict=False
            )
        )

    ambiguous_tree = tmp_path / "ambiguous_species_mapping.nwk"
    ambiguous_tree.write_text(
        "(Ancestor_species:1,Reference_species:0)Reference_species;\n",
        encoding="utf-8",
    )
    ambiguous_result = subprocess.run(
        [
            "Rscript",
            str(SUPPORT_DIR / "plot_query2family_presence_absence.R"),
            f"--species_tree={ambiguous_tree}",
            f"--long_table={long_table}",
            f"--ortholog_column_table={out_dir / 'columns.tsv'}",
            f"--ortholog_glyph_table={out_dir / 'glyphs.tsv'}",
            f"--ortholog_tree_table={out_dir / 'tree.tsv'}",
            "--reference_species=Reference_species",
            f"--out_svg={out_dir / 'ambiguous.svg'}",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert ambiguous_result.returncode != 0
    assert "Mapped duplication label occurs more than once" in (
        ambiguous_result.stdout + ambiguous_result.stderr
    )


def test_eight_families_and_twenty_three_query_columns_render_without_cross_family_failure(
    tmp_path: Path,
):
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript is unavailable")
    mod = load_module("query_gene_orthologs.py")
    query_dir = tmp_path / "input" / "query_gene"
    output_root = tmp_path / "output" / "query2family"
    query_dir.mkdir(parents=True)
    query_counts = [1, 2, 3, 4, 5, 1, 3, 4]
    family_ids = [f"FAM{index + 1:02d}" for index in range(len(query_counts))]
    reference_species = "Reference_species"
    for family_id, query_count in zip(family_ids, query_counts, strict=True):
        write_query_only_family(
            query_dir,
            output_root,
            family_id,
            reference_species,
            query_count,
        )
    selection = tmp_path / "selection.tsv"
    selection.write_text(
        "family_id\n" + "\n".join(reversed(family_ids)) + "\n",
        encoding="utf-8",
    )
    out_dir = tmp_path / "summary"
    mod.run(
        SimpleNamespace(
            dir_gene_family=str(output_root),
            dir_query_gene=str(query_dir),
            family_file=str(selection),
            reference_species=reference_species,
            out_columns=str(out_dir / "columns.tsv"),
            out_glyphs=str(out_dir / "glyphs.tsv"),
            out_tree=str(out_dir / "tree.tsv"),
            out_synteny=str(out_dir / "synteny.tsv"),
            out_ufboot=str(out_dir / "ufboot.tsv"),
        )
    )
    columns = pandas.read_csv(out_dir / "columns.tsv", sep="\t")
    assert columns.shape[0] == sum(query_counts) == 23
    assert columns["family_id"].drop_duplicates().tolist() == list(reversed(family_ids))

    rendered_family_ids = [
        f"FAMILY_WITH_A_VERY_LONG_INPUT_FILENAME_{index + 1:02d}"
        for index in range(len(family_ids))
    ]
    rendered_family_by_source = dict(zip(family_ids, rendered_family_ids, strict=True))
    for table_name in ("columns.tsv", "glyphs.tsv", "tree.tsv", "synteny.tsv", "ufboot.tsv"):
        table_path = out_dir / table_name
        table = pandas.read_csv(table_path, sep="\t")
        table["family_id"] = table["family_id"].map(rendered_family_by_source)
        table.to_csv(table_path, sep="\t", index=False)

    species_tree = tmp_path / "species.nwk"
    species = [reference_species, "Other_species"]
    species_tree.write_text("(Reference_species:1,Other_species:1);\n", encoding="utf-8")
    long_table = tmp_path / "long.tsv"
    pandas.DataFrame(
        {
            "species": species,
            "species_display": [value.replace("_", " ") for value in species],
            "query": ["FAM01", "FAM01"],
            "query_order": [1, 1],
            "presence": [1, 0],
            "copy_number": [1, 0],
            "status": ["complete"] * len(species),
        }
    ).to_csv(long_table, sep="\t", index=False)
    pdf_path = out_dir / "stress.pdf"
    svg_path = out_dir / "stress.svg"
    subprocess.run(
        [
            "Rscript",
            str(SUPPORT_DIR / "plot_query2family_presence_absence.R"),
            f"--species_tree={species_tree}",
            f"--long_table={long_table}",
            f"--ortholog_column_table={out_dir / 'columns.tsv'}",
            f"--ortholog_glyph_table={out_dir / 'glyphs.tsv'}",
            f"--ortholog_tree_table={out_dir / 'tree.tsv'}",
            f"--ortholog_synteny_table={out_dir / 'synteny.tsv'}",
            f"--ortholog_ufboot_table={out_dir / 'ufboot.tsv'}",
            f"--reference_species={reference_species}",
            f"--out_pdf={pdf_path}",
            f"--out_svg={svg_path}",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert pdf_path.stat().st_size > 1000
    svg_text = svg_path.read_text(encoding="utf-8")
    for family_id in rendered_family_ids:
        assert family_id in svg_text
    assert "Friendly_FAM" not in svg_text
    assert "Branch length" in svg_text
    assert "Million years ago" not in svg_text
    assert "Duplication bar colors" in svg_text
    assert "Bar height = duplication count" in svg_text
    assert "Circle area = duplication count" not in svg_text
    assert "D#: query-gene duplication" not in svg_text
    assert "D#: mapped duplication" not in svg_text
    assert not re.search(r">D[0-9]+</text>", svg_text)
    for family_id in rendered_family_ids:
        assert re.search(rf">{re.escape(family_id)}</text>", svg_text)
    assert not re.search(r">F[0-9]+:", svg_text)
    query_label_boxes = [
        (float(y), float(length))
        for y, length in re.findall(
            r"translate\([0-9.]+,([0-9.]+)\) rotate\(-90\).*?textLength='([0-9.]+)px'.*?>GENE_",
            svg_text,
        )
    ]
    legend_y_match = re.search(
        r"<text x='[0-9.]+' y='([0-9.]+)'.*?>Reference species orthologs</text>",
        svg_text,
    )
    assert query_label_boxes and legend_y_match
    assert float(legend_y_match.group(1)) > max(y + length for y, length in query_label_boxes) + 8
    assert "Evidence band (top): local synteny anchors" not in svg_text
    assert "Evidence band (bottom): Gene tree UFBoot (%)" not in svg_text
    assert "unavailable" not in svg_text
    assert "Evidence band states" in svg_text
    assert "reference self" in svg_text
    assert "No band: reference self" not in svg_text
    assert count_svg_matrix_edge_bands(
        svg_path,
        "Reference species",
        column_count=sum(query_counts),
    ) == 2 * sum(query_counts)
    assert_svg_copy_numbers_above_edge_bands(svg_path, "Reference species")

    vertical_family_title_boxes = []
    family_legend_x = []
    for family_id in rendered_family_ids:
        title_match = re.search(
            rf"translate\([0-9.]+,([0-9.]+)\) rotate\(-90\).*?"
            rf"textLength='([0-9.]+)px'.*?>{re.escape(family_id)}</text>",
            svg_text,
        )
        assert title_match
        vertical_family_title_boxes.append(
            (float(title_match.group(1)), float(title_match.group(2)))
        )
        legend_match = re.search(
            rf"<text x='([0-9.]+)' y='[0-9.]+' .*?>{re.escape(family_id)}</text>",
            svg_text,
        )
        assert legend_match
        family_legend_x.append(float(legend_match.group(1)))
    assert vertical_family_title_boxes
    assert min(y - length for y, length in vertical_family_title_boxes) > 5
    assert len({round(value, 2) for value in family_legend_x}) == 1


@pytest.mark.parametrize(
    ("cds_fasta_id", "species", "expected"),
    [
        ("Arabidopsis_thaliana_AT3G26790", "Arabidopsis_thaliana", "AT3G26790"),
        ("Species.one.GENE.1", "Species.one", "GENE.1"),
        ("Species-one-GENE-1", "Species-one", "GENE-1"),
        ("unprefixed_gene", "Species_one", "unprefixed_gene"),
    ],
)
def test_gene_id_is_derived_from_exact_cds_fasta_id(cds_fasta_id, species, expected):
    mod = load_module("query_gene_orthologs.py")
    assert mod.gene_id_from_cds_fasta_id(cds_fasta_id, species) == expected


def test_missing_query_tip_in_cds_fasta_is_a_hard_error():
    mod = load_module("query_gene_orthologs.py")
    by_id = {1: {"node_name": "Species_GENE1", "spnode_coverage": "Species"}}
    with pytest.raises(ValueError, match="not found in the family CDS FASTA"):
        mod.resolve_query_cds_definitions(
            by_id=by_id,
            query_tip_by_id={"query1": 1},
            definitions=[{"query_id": "query1", "query_label": "Friendly"}],
            cds_fasta_ids=["Species_OTHER"],
            family_id="FAM",
        )


def test_missing_query_marker_is_a_hard_error():
    mod = load_module("query_gene_orthologs.py")
    rows = [stat_row(0, -1, -1, -1, "L", "Species_GENE1", "Species")]
    with pytest.raises(ValueError, match="no gene-tree tip marker"):
        mod.select_query_tip_nodes(rows, [{"query_id": "query1", "query_label": "Friendly"}])


def test_reference_species_selection_includes_unmarked_tips_and_excludes_other_species():
    mod = load_module("query_gene_orthologs.py")
    rows = [
        stat_row(0, 3, -1, -1, "L", "Reference_species_GENE1", "Reference_species", "direct:q1"),
        stat_row(1, 3, -1, -1, "L", "Reference_species_GENE2", "Reference_species"),
        stat_row(2, 3, -1, -1, "L", "Other_species_GENE3", "Other_species"),
        stat_row(3, -1, 0, 1, "S", "root", "n0"),
    ]

    selected, definitions = mod.select_reference_tip_nodes(rows, "Reference species")

    assert selected == {
        "Reference_species_GENE1": 0,
        "Reference_species_GENE2": 1,
    }
    assert [definition["query_id"] for definition in definitions] == list(selected)


def test_stat_branch_and_cds_fasta_can_both_be_read_from_zip_storage(tmp_path: Path):
    mod = load_module("query_gene_orthologs.py")
    query_dir = tmp_path / "query_gene"
    output_root = tmp_path / "query2family"
    query_dir.mkdir()
    write_family(
        query_dir,
        output_root,
        "ZIP",
        [("q1", "Friendly")],
        [stat_row(0, -1, -1, -1, "L", "Zip_species_GENE_Z", "Zip_species", "direct:q1")],
        ["Zip_species_GENE_Z"],
    )
    synteny_dir = output_root / "synteny"
    synteny_dir.mkdir(parents=True)
    pandas.DataFrame(
        [
            [
                "Zip_species_GENE_Z",
                "Zip_species",
                "upstream",
                -1,
                "Zip_species_NEIGHBOR",
                "G1",
                1,
            ]
        ],
        columns=[
            "node_name", "species", "direction", "offset", "neighbor_gene", "group_id", "group_size"
        ],
    ).to_csv(synteny_dir / "ZIP_synteny.tsv", sep="\t", index=False)
    for subdir, suffix in (("stat_tree", "_stat.tree.tsv"), ("tree_plot", "_tree_plot.pdf")):
        path = output_root / subdir / f"ZIP{suffix}"
        path.parent.mkdir(parents=True)
        path.write_text(f"{subdir}\n", encoding="utf-8")
    family_ids, family_from_name = family_context("query2family", query_dir=query_dir)
    archive_completed_outputs(output_root, "query2family", family_ids, family_from_name)

    columns, glyphs, _ = mod.collect_query_gene_orthologs(
        output_root,
        query_dir,
        reference_species="Zip_species",
    )
    evidence = mod.collect_reference_synteny_evidence(
        mod.GeneFamilyOutputStore(output_root), columns, glyphs
    )
    ufboot_evidence = mod.collect_reference_ufboot_evidence(
        mod.GeneFamilyOutputStore(output_root), columns, glyphs
    )

    assert columns[0]["gene_id"] == "GENE_Z"
    assert evidence[0]["synteny_status"] == "reference_self"
    assert ufboot_evidence[0]["orthology_ufboot_status"] == "reference_self"
    assert not (output_root / "stat_branch").exists()
    assert not (output_root / "cds_fasta").exists()
    assert not (output_root / "synteny").exists()

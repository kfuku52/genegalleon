import shutil
import subprocess
import xml.etree.ElementTree as ET
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace

import pandas
import pytest

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
CORE_SCRIPT = Path(__file__).resolve().parents[1] / "core" / "gg_gene_summary_core.sh"
ENTRYPOINT_SCRIPT = Path(__file__).resolve().parents[1] / "gg_gene_summary_entrypoint.sh"
CONFIG_REGISTRY = SUPPORT_DIR / "gg_entrypoint_config_vars.sh"


def load_module(name: str):
    path = SUPPORT_DIR / name
    spec = spec_from_file_location(f"workflow.support.{path.stem}", path)
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


@pytest.fixture
def weak_duplication_fixture(tmp_path):
    """One extra Shared copy crosses a D with one shared species out of twenty."""
    query_dir = tmp_path / "queries"
    output_root = tmp_path / "families"
    (output_root / "stat_branch").mkdir(parents=True)
    (output_root / "cds_fasta").mkdir()
    query_dir.mkdir()
    (query_dir / "FAM").write_text(
        ">q1\nAAAA\n>q2\nAAAA\n", encoding="utf-8"
    )
    reference_clade = (
        "D", "reference_duplication",
        ("Reference_species", "Reference_species_REF1", "direct:q1"),
        ("Reference_species", "Reference_species_REF2", "direct:q2"),
    )
    left_clade = (
        "S", "shared_speciation", reference_clade,
        ("Shared_species", "Shared_species_strict", ""),
    )
    for number in range(18):
        species = f"Other_species_{number:02d}"
        left_clade = ("S", f"speciation_{number}", left_clade, (species, f"{species}_gene", ""))
    tree = (
        "S", "root", ("Outside_species", "Outside_species_gene", ""),
        ("D", "weak_duplication", left_clade, ("Shared_species", "Shared_species_candidate", "")),
    )
    rows = []

    def visit(node, parent):
        branch_id = len(rows)
        rows.append(None)
        if len(node) == 3:
            species, name, marker = node
            rows[branch_id] = stat_row(branch_id, parent, -1, -1, "L", name, species, marker)
        else:
            event, name, left, right = node
            child1 = visit(left, branch_id)
            child2 = visit(right, branch_id)
            mapping = "Reference_species" if name == "reference_duplication" else "n0"
            rows[branch_id] = stat_row(branch_id, parent, child1, child2, event, name, mapping, generax_support=100)
        return branch_id

    visit(tree, -1)
    stat_path = output_root / "stat_branch" / "FAM_stat.branch.tsv"
    pandas.DataFrame(rows).to_csv(stat_path, sep="\t", index=False)
    tips = [row for row in rows if row["so_event"] == "L"]
    (output_root / "cds_fasta" / "FAM_cds.fasta").write_text(
        "".join(f">{row['node_name']}\nATG\n" for row in tips), encoding="utf-8"
    )
    species = sorted({row["spnode_coverage"] for row in tips})
    species_tree = tmp_path / "species.nwk"
    species_tree.write_text("(" + ",".join(f"{sp}:1" for sp in species) + ")n0;\n", encoding="utf-8")
    long_table = tmp_path / "long.tsv"
    pandas.DataFrame([
        dict(species=sp, species_display=sp.replace("_", " "), query="FAM", query_order=1,
             presence=1, copy_number=sum(row["spnode_coverage"] == sp for row in tips), status="complete")
        for sp in species
    ]).to_csv(long_table, sep="\t", index=False)
    return dict(query_dir=query_dir, output_root=output_root, species_tree=species_tree,
                long_table=long_table, weak_branch=next(row["branch_id"] for row in rows if row["node_name"] == "weak_duplication"))


def collect_weak_fixture(fixture, out_dir, threshold, basis="reference_species"):
    mod = load_module("query_gene_orthologs.py")
    mod.run(SimpleNamespace(
        basis=basis, dir_gene_family=str(fixture["output_root"]), dir_query_gene=str(fixture["query_dir"]),
        family_file="", reference_species="Reference_species", dup_conf_score_threshold=threshold,
        out_columns=str(out_dir / "columns.tsv"), out_glyphs=str(out_dir / "glyphs.tsv"),
        out_tree=str(out_dir / "tree.tsv"), out_synteny=str(out_dir / "synteny.tsv"),
        out_ufboot=str(out_dir / "ufboot.tsv"), out_dup_conf=str(out_dir / "dup_conf.tsv"),
        out_query_map=str(out_dir / "query_map.tsv"),
    ))
    return out_dir


@pytest.mark.parametrize("basis", ["reference_species", "query_gene"])
def test_weak_duplication_cutoff_retains_strict_calls_events_and_anchor_identity(
    weak_duplication_fixture, tmp_path, basis,
):
    baseline = collect_weak_fixture(weak_duplication_fixture, tmp_path / "strict", 0, basis)
    below = collect_weak_fixture(weak_duplication_fixture, tmp_path / "below", 0.049, basis)
    flagged = collect_weak_fixture(weak_duplication_fixture, tmp_path / "flagged", 0.05, basis)
    assert (flagged / "columns.tsv").read_bytes() == (baseline / "columns.tsv").read_bytes()
    # Candidate membership changes bar provenance, while every original tree
    # node, S/D call, map and index remains identical.
    provenance = ["displayed_gene_ids", "displayed_child1_gene_ids", "displayed_child2_gene_ids"]
    original_tree = pandas.read_csv(baseline / "tree.tsv", sep="\t", keep_default_na=False)
    candidate_tree = pandas.read_csv(flagged / "tree.tsv", sep="\t", keep_default_na=False)
    pandas.testing.assert_frame_equal(original_tree.drop(columns=provenance), candidate_tree.drop(columns=provenance))
    original_d = original_tree.loc[original_tree.node_id == weak_duplication_fixture["weak_branch"]].iloc[0]
    candidate_d = candidate_tree.loc[candidate_tree.node_id == weak_duplication_fixture["weak_branch"]].iloc[0]
    assert not original_d.displayed_child1_gene_ids or not original_d.displayed_child2_gene_ids
    assert candidate_d.displayed_child1_gene_ids and candidate_d.displayed_child2_gene_ids
    assert (below / "glyphs.tsv").read_bytes() == (baseline / "glyphs.tsv").read_bytes()
    assert pandas.read_csv(baseline / "dup_conf.tsv", sep="\t").empty
    strict = pandas.read_csv(baseline / "glyphs.tsv", sep="\t", keep_default_na=False)
    glyphs = pandas.read_csv(flagged / "glyphs.tsv", sep="\t", keep_default_na=False)
    old_columns = [col for col in strict.columns if col not in ("lane_index", "lane_count")]
    pandas.testing.assert_frame_equal(
        strict[old_columns].sort_values(["species", "start_order"]).reset_index(drop=True),
        glyphs.loc[glyphs.relation != "weak_duplication", old_columns]
        .sort_values(["species", "start_order"]).reset_index(drop=True),
    )
    weak = glyphs[glyphs.relation == "weak_duplication"]
    assert len(weak) == 1
    assert weak.iloc[0].gene_ids == "Shared_species_candidate"
    assert weak.iloc[0].copy_number == 1
    assert (weak.iloc[0].start_order, weak.iloc[0].end_order) == (1, 2)
    assert set(glyphs.loc[glyphs.species == "Shared_species", "gene_ids"]) == {
        "Shared_species_strict", "Shared_species_candidate",
    }
    evidence = pandas.read_csv(flagged / "dup_conf.tsv", sep="\t")
    assert len(evidence) == 2
    assert set(evidence.mrca_branch_id) == {weak_duplication_fixture["weak_branch"]}
    assert set(evidence.shared_species_count) == {1}
    assert set(evidence.union_species_count) == {20}
    assert set(evidence.dup_conf_score) == {0.05}
    assert set(evidence.branch_ufboot) == {100}
    anchor_field = "anchor_cds_fasta_id" if basis == "query_gene" else "reference_cds_fasta_id"
    assert set(evidence[anchor_field]) == {"Reference_species_REF1", "Reference_species_REF2"}
    if basis == "query_gene":
        assert weak.iloc[0].anchor_query_ids == "q1;q2"
    ufboot = pandas.read_csv(flagged / "ufboot.tsv", sep="\t", keep_default_na=False)
    weak_support = ufboot[ufboot.relation == "weak_duplication"]
    assert len(weak_support) == 2
    assert set(weak_support.orthology_mrca_event) == {"D"}
    assert set(weak_support.orthology_ufboot_status) == {"not_evaluable"}
    assert set(weak_support.orthology_ufboot_unavailable_reason) == {"weak_duplication"}
    assert set(weak_support.decisive_branch_ufboot) == {""}
    permissive = collect_weak_fixture(weak_duplication_fixture, tmp_path / "permissive", 1, basis)
    all_evidence = pandas.read_csv(permissive / "dup_conf.tsv", sep="\t")
    assert "Reference_species" not in set(all_evidence.species)


@pytest.mark.parametrize("threshold", ["nan", "inf", "-inf", "-0.01", "1.01", "bad"])
def test_invalid_duplication_cutoff_fails_before_reading_inputs(threshold):
    mod = load_module("query_gene_orthologs.py")
    with pytest.raises(ValueError, match="finite number between 0 and 1"):
        mod.run(SimpleNamespace(dup_conf_score_threshold=threshold))


def test_positive_duplication_cutoff_requires_its_audit_output_before_reading_inputs():
    mod = load_module("query_gene_orthologs.py")
    with pytest.raises(ValueError, match="out_dup_conf is required"):
        mod.run(SimpleNamespace(dup_conf_score_threshold=0.05))


@pytest.mark.parametrize("basis", ["reference_species", "query_gene"])
def test_additional_candidate_requires_a_saved_cds_sequence(weak_duplication_fixture, tmp_path, basis):
    fixture = weak_duplication_fixture
    cds = fixture["output_root"] / "cds_fasta" / "FAM_cds.fasta"
    cds.write_text(cds.read_text().replace(">Shared_species_candidate\nATG\n", ""))
    with pytest.raises(ValueError, match="candidate.*CDS FASTA"):
        collect_weak_fixture(fixture, tmp_path / "missing_cds", 0.05, basis)


def test_additional_candidates_reject_inconsistent_saved_species_overlap(weak_duplication_fixture, tmp_path):
    fixture = weak_duplication_fixture
    path = fixture["output_root"] / "stat_branch/FAM_stat.branch.tsv"
    table = pandas.read_csv(path, sep="\t", keep_default_na=False)
    # The root splits disjoint species sets; a saved D label cannot have score zero.
    table.loc[table.node_name == "root", "so_event"] = "D"
    table.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="Saved species-overlap event disagrees"):
        collect_weak_fixture(fixture, tmp_path / "inconsistent", 0.05)


@pytest.mark.parametrize("basis", ["reference_species", "query_gene"])
@pytest.mark.parametrize("corruption", [
    "ufboot_relation", "ufboot_gene", "synteny_gene", "cutoff", "score", "counts", "mrca",
    "missing_pair", "missing_support_pair", "missing_table", "glyph_gene_count", "zero_overlap",
    "mrca_without_support", "wrong_D_without_support", "mixed_mrcas_without_support", "wrong_descendants",
])
def test_plot_rejects_mixed_additional_candidate_evidence(
    weak_duplication_fixture, tmp_path, basis, corruption,
):
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript is unavailable")
    fixture = weak_duplication_fixture
    out_dir = collect_weak_fixture(fixture, tmp_path / "mixed", 0.05, basis)
    if corruption == "wrong_descendants":
        path = out_dir / "tree.tsv"
        table = pandas.read_csv(path, sep="\t", keep_default_na=False)
        weak = table.node_id == fixture["weak_branch"]
        table.loc[weak, "displayed_child1_gene_ids"] += ";Shared_species_candidate"
        table.loc[weak, "displayed_child2_gene_ids"] = ""
        table.to_csv(path, sep="\t", index=False)
    elif corruption in ("score", "counts", "mrca", "missing_pair", "zero_overlap",
                        "mrca_without_support", "wrong_D_without_support", "mixed_mrcas_without_support"):
        path = out_dir / "dup_conf.tsv"
        table = pandas.read_csv(path, sep="\t", keep_default_na=False)
        if corruption == "score":
            table.loc[0, "dup_conf_score"] = 0.01
        elif corruption == "counts":
            table.loc[0, "union_species_count"] = 19
        elif corruption == "mrca":
            table.loc[0, "mrca_branch_id"] = 0
        elif corruption == "zero_overlap":
            table.loc[:, "shared_species_count"] = 0
            table.loc[:, "dup_conf_score"] = 0
        elif corruption == "mrca_without_support":
            table.loc[:, "mrca_branch_id"] = 0
        elif corruption == "mixed_mrcas_without_support":
            table.loc[0, "mrca_branch_id"] = 0
        elif corruption == "wrong_D_without_support":
            tree = pandas.read_csv(out_dir / "tree.tsv", sep="\t")
            other_d = tree.loc[(tree.event == "D") & (tree.node_id != fixture["weak_branch"]), "node_id"].iloc[0]
            table.loc[:, "mrca_branch_id"] = other_d
        else:
            table = table.iloc[1:]
        table.to_csv(path, sep="\t", index=False)
    elif corruption == "glyph_gene_count":
        path = out_dir / "glyphs.tsv"
        table = pandas.read_csv(path, sep="\t", keep_default_na=False)
        table.loc[table.relation == "weak_duplication", "gene_ids"] = "Shared_species_candidate;Shared_species_wrong_gene"
        table.to_csv(path, sep="\t", index=False)
    elif corruption not in ("cutoff", "missing_table"):
        name = "synteny" if corruption == "synteny_gene" else "ufboot"
        path = out_dir / f"{name}.tsv"
        table = pandas.read_csv(path, sep="\t", keep_default_na=False)
        weak = table.relation == "weak_duplication"
        if corruption == "missing_support_pair":
            table = table.drop(table[weak].index[0])
        elif corruption == "ufboot_relation":
            table.loc[weak, "relation"] = "shared_ancestral"
            table.loc[weak, "orthology_mrca_event"] = "S"
            table.loc[weak, "orthology_ufboot_status"] = "evaluated"
            table.loc[weak, "decisive_branch_ufboot"] = "100"
            table.loc[weak, "orthology_ufboot_unavailable_reason"] = ""
        else:
            table.loc[weak, "candidate_cds_fasta_id"] = "Shared_species_wrong_gene"
        table.to_csv(path, sep="\t", index=False)
    command = [
        "Rscript", str(SUPPORT_DIR / "plot_query2family_presence_absence.R"),
        f"--species_tree={fixture['species_tree']}", f"--species_mapping_tree={fixture['species_tree']}",
        f"--long_table={fixture['long_table']}", f"--ortholog_basis={basis}",
        "--reference_species=Reference_species",
        f"--dup_conf_score_threshold={'0.04' if corruption == 'cutoff' else '0.05'}",
        f"--ortholog_column_table={out_dir / 'columns.tsv'}", f"--ortholog_glyph_table={out_dir / 'glyphs.tsv'}",
        f"--ortholog_tree_table={out_dir / 'tree.tsv'}", f"--ortholog_synteny_table={out_dir / 'synteny.tsv'}",
        f"--ortholog_ufboot_table={out_dir / 'ufboot.tsv'}",
        f"--ortholog_dup_conf_table={out_dir / 'dup_conf.tsv'}", f"--out_svg={out_dir / 'mixed.svg'}",
    ]
    if corruption == "missing_table":
        command = [arg for arg in command if not arg.startswith("--ortholog_dup_conf_table=")]
    if corruption.endswith("without_support") or corruption == "wrong_descendants":
        command = [arg for arg in command if not arg.startswith("--ortholog_ufboot_table=")]
    result = subprocess.run(command, check=False, capture_output=True, text=True)
    assert result.returncode != 0, f"Mixed {corruption} evidence was accepted"
    assert ("require --ortholog_dup_conf_table" if corruption == "missing_table" else "disagree") in result.stderr


@pytest.mark.parametrize("basis", ["reference_species", "query_gene"])
@pytest.mark.parametrize("legend_columns,evidence_layout", [(1, "band"), (2, "rail"), (3, "glyph"), (3, "off")])
def test_weak_duplication_candidates_render_in_orange_without_false_ufboot(
    weak_duplication_fixture, tmp_path, basis, legend_columns, evidence_layout,
):
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript is unavailable")
    fixture = weak_duplication_fixture
    out_dir = collect_weak_fixture(fixture, tmp_path / "plot", 0.05, basis)
    plot_command = [
        "Rscript", str(SUPPORT_DIR / "plot_query2family_presence_absence.R"),
        f"--species_tree={fixture['species_tree']}", f"--species_mapping_tree={fixture['species_tree']}",
        f"--long_table={fixture['long_table']}", f"--ortholog_basis={basis}",
        "--reference_species=Reference_species", "--dup_conf_score_threshold=0.05",
        f"--ortholog_column_table={out_dir / 'columns.tsv'}", f"--ortholog_glyph_table={out_dir / 'glyphs.tsv'}",
        f"--ortholog_tree_table={out_dir / 'tree.tsv'}", f"--ortholog_synteny_table={out_dir / 'synteny.tsv'}",
        f"--ortholog_ufboot_table={out_dir / 'ufboot.tsv'}",
        f"--ortholog_dup_conf_table={out_dir / 'dup_conf.tsv'}",
        f"--legend_columns={legend_columns}", f"--evidence_layout={evidence_layout}",
        f"--out_svg={out_dir / 'presence.svg'}", f"--out_pdf={out_dir / 'presence.pdf'}",
    ]
    subprocess.run(plot_command, check=True, capture_output=True, text=True)
    svg = (out_dir / "presence.svg").read_text(encoding="utf-8")
    assert "Additional ortholog candidate" in svg
    assert "duplication confidence score &lt;=0.05" in svg
    assert "weak-D candidate" not in svg
    assert "#FDBA74" in svg.upper()
    assert "#2166AC" in svg.upper()
    assert (out_dir / "presence.pdf").stat().st_size > 1000
    # Copy numbers need physical room inside both stacked lanes, including the
    # central 64% left between top/bottom evidence bands.
    elements = list(ET.parse(out_dir / "presence.svg").getroot().iter())
    numeric_labels = [element for element in elements if element.tag.endswith("text") and element.text == "1"]
    candidate_cells = [element for element in elements if element.tag.endswith("rect") and
                       "fill: #FDBA74" in element.attrib.get("style", "") and
                       any(float(element.attrib["x"]) <= float(label.attrib.get("x", -1)) <=
                           float(element.attrib["x"]) + float(element.attrib["width"]) and
                           float(element.attrib["y"]) <= float(label.attrib.get("y", -1)) <=
                           float(element.attrib["y"]) + float(element.attrib["height"])
                           for label in numeric_labels)]
    assert candidate_cells
    available_fraction = 0.64 if evidence_layout == "band" else 1
    assert min(float(cell.attrib["height"]) * available_fraction for cell in candidate_cells) >= 8
    if legend_columns == 1:
        too_short = subprocess.run([*plot_command, "--height=1"], check=False, capture_output=True, text=True)
        assert too_short.returncode != 0
        assert "height is too small for copy-number labels" in too_short.stderr
    support = pandas.read_csv(out_dir / "ufboot.tsv", sep="\t", keep_default_na=False)
    weak = support.relation == "weak_duplication"
    support.loc[weak, "orthology_ufboot_status"] = "evaluated"
    support.loc[weak, "decisive_branch_ufboot"] = "100"
    support.loc[weak, "orthology_ufboot_unavailable_reason"] = ""
    support.to_csv(out_dir / "ufboot.tsv", sep="\t", index=False)
    result = subprocess.run(plot_command, check=False, capture_output=True, text=True)
    assert result.returncode != 0
    assert "cannot report speciation-based orthology UFBoot" in result.stderr


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


def test_weak_duplication_candidates_do_not_paint_intervening_unassigned_anchors(
    weak_duplication_fixture,
):
    mod = load_module("query_gene_orthologs.py")
    fixture = weak_duplication_fixture
    columns, _, _ = mod.collect_query_gene_orthologs(
        fixture["output_root"], fixture["query_dir"], "Reference_species"
    )
    rows = mod.read_stat_branch(mod.GeneFamilyOutputStore(fixture["output_root"]), "FAM")
    outside = next(row for row in rows if row["node_name"] == "Outside_species_gene")
    # A query basis can mix anchor species; the middle anchor has a strict S MRCA.
    middle = dict(columns[0], column_order=2, cds_fasta_id=outside["node_name"],
                  gene_id="gene", reference_tip_branch_id=int(outside["branch_id"]))
    columns[1]["column_order"] = 3
    columns.insert(1, middle)
    glyphs, evidence = mod.add_weak_duplication_candidates(
        mod.GeneFamilyOutputStore(fixture["output_root"]), columns, [], 0.05
    )
    candidate_glyphs = [row for row in glyphs if row["gene_ids"] == "Shared_species_candidate"]
    assert {(row["start_order"], row["end_order"]) for row in candidate_glyphs} == {(1, 1), (3, 3)}
    assert len(evidence) == 2


def test_query_gene_orthologs_connect_preduplication_copy_and_count_paralogs(tmp_path: Path):
    mod = load_module("query_gene_orthologs.py")
    query_dir = tmp_path / "input" / "query_gene"
    stat_dir = tmp_path / "output" / "query2family" / "stat_branch"
    cds_dir = stat_dir.parent / "cds_fasta"
    query_dir.mkdir(parents=True)
    stat_dir.mkdir(parents=True)
    cds_dir.mkdir(parents=True)
    (query_dir / "AFL").write_text(
        ">qLEC description | LEC2 | Arabidopsis thaliana\nAAAA\n"
        ">qFUS description | FUS3 | Arabidopsis thaliana\nAAAA\n",
        encoding="utf-8",
    )
    rows = [
        stat_row(0, -1, 1, 2, "S", "root"),
        stat_row(1, 0, -1, -1, "L", "Ancestor_gene", "Ancestor_species"),
        stat_row(2, 0, 3, 4, "D", "FUS_LEC_duplication", "n1"),
        stat_row(3, 2, 5, 6, "S", "FUS_clade"),
        stat_row(4, 2, 9, 10, "S", "LEC_clade"),
        stat_row(5, 3, -1, -1, "L", "Arabidopsis_thaliana_AT3G26790", "Arabidopsis_thaliana", "direct:qFUS"),
        stat_row(6, 3, 7, 8, "D", "Beta_duplication"),
        stat_row(7, 6, -1, -1, "L", "Beta_FUS_a", "Beta_species"),
        stat_row(8, 6, -1, -1, "L", "Beta_FUS_b", "Beta_species"),
        stat_row(9, 4, -1, -1, "L", "Arabidopsis_thaliana_AT1G28300", "Arabidopsis_thaliana", "best:qLEC"),
        stat_row(10, 4, -1, -1, "L", "Gamma_LEC", "Gamma_species"),
    ]
    rows[2]["spnode_generax"] = "n1_generax"
    rows[6]["spnode_generax"] = "Beta_species"
    rows[5]["spnode_generax"] = "Wrong_tip_species"
    pandas.DataFrame(rows).to_csv(stat_dir / "AFL_stat.branch.tsv", sep="\t", index=False)
    (cds_dir / "AFL_cds.fasta").write_text(
        ">Arabidopsis_thaliana_AT3G26790\nATG\n"
        ">Arabidopsis_thaliana_AT1G28300\nATG\n",
        encoding="utf-8",
    )

    columns, glyphs, tree_nodes = mod.collect_query_gene_orthologs(
        dir_gene_family=stat_dir.parent,
        dir_query_gene=query_dir,
        reference_species="Arabidopsis_thaliana",
    )

    assert {row["reference_species"] for row in columns} == {"Arabidopsis_thaliana"}
    assert [row["gene_id"] for row in columns] == ["AT3G26790", "AT1G28300"]
    assert [row["cds_fasta_id"] for row in columns] == [
        "Arabidopsis_thaliana_AT3G26790",
        "Arabidopsis_thaliana_AT1G28300",
    ]
    shared = [row for row in glyphs if row["species"] == "Ancestor_species"]
    assert len(shared) == 1
    assert shared[0]["relation"] == "shared_ancestral"
    assert shared[0]["reference_gene_ids"] == "AT3G26790;AT1G28300"
    assert (shared[0]["start_order"], shared[0]["end_order"]) == (1, 2)
    assert shared[0]["copy_number"] == 1

    beta = [row for row in glyphs if row["species"] == "Beta_species"]
    assert len(beta) == 1
    assert beta[0]["relation"] == "specific"
    assert beta[0]["reference_gene_ids"] == "AT3G26790"
    assert beta[0]["copy_number"] == 2

    tree_by_node = {row["node_id"]: row for row in tree_nodes}
    assert len(tree_nodes) == 4
    fus_tip = next(row for row in tree_nodes if row["gene_id"] == "AT3G26790")
    lec_tip = next(row for row in tree_nodes if row["gene_id"] == "AT1G28300")
    assert fus_tip["parent_node_id"] == lec_tip["parent_node_id"] == 2
    assert tree_by_node[2]["event"] == "D"
    assert tree_by_node[2]["node_height"] == 1
    assert tree_by_node[2]["plot_order"] == 1.5
    assert tree_by_node[2]["mapped_species_node"] == "n1_generax"
    assert tree_by_node[2]["duplication_index"] == 1
    assert tree_by_node[2]["in_reference_tree"] == 1
    assert tree_by_node[6]["mapped_species_node"] == "Beta_species"
    assert tree_by_node[6]["in_reference_tree"] == 0
    assert fus_tip["mapped_species_node"] == "Arabidopsis_thaliana"


def test_query_gene_ortholog_glyph_lanes_only_split_overlapping_spans():
    mod = load_module("query_gene_orthologs.py")
    glyphs = [
        {"species": "Sp", "start_order": 1, "end_order": 2, "family_order": 1},
        {"species": "Sp", "start_order": 2, "end_order": 3, "family_order": 2},
        {"species": "Sp", "start_order": 4, "end_order": 4, "family_order": 3},
    ]

    mod.assign_lanes(glyphs)

    assert [glyph["lane_index"] for glyph in glyphs] == [1, 2, 1]
    assert {glyph["lane_count"] for glyph in glyphs} == {2}


def test_query_gene_ortholog_glyph_lanes_are_independent_between_families():
    mod = load_module("query_gene_orthologs.py")
    glyphs = [
        {"species": "Sp", "family_id": "A", "start_order": 1, "end_order": 2, "family_order": 1},
        {"species": "Sp", "family_id": "A", "start_order": 2, "end_order": 2, "family_order": 1},
        {"species": "Sp", "family_id": "B", "start_order": 3, "end_order": 3, "family_order": 2},
    ]

    mod.assign_lanes(glyphs)

    assert [glyph["lane_count"] for glyph in glyphs] == [2, 2, 1]


def test_local_synteny_uses_distinct_shared_groups_and_retains_reverse_order():
    mod = load_module("query_gene_orthologs.py")
    reference_rows = [
        {"offset": -4, "group_id": "G1", "neighbor_gene": "r1"},
        {"offset": -2, "group_id": "G2", "neighbor_gene": "r2"},
        {"offset": 3, "group_id": "G3", "neighbor_gene": "r3"},
    ]
    candidate_rows = [
        {"offset": -3, "group_id": "G3", "neighbor_gene": "c3"},
        {"offset": 1, "group_id": "G2", "neighbor_gene": "c2a"},
        {"offset": 2, "group_id": "G2", "neighbor_gene": "c2b"},
        {"offset": 5, "group_id": "G1", "neighbor_gene": "c1"},
    ]

    metrics = mod.local_synteny_metrics(reference_rows, candidate_rows, window_radius=5)

    assert metrics["shared_anchor_count"] == 3
    assert metrics["local_synteny_score"] == pytest.approx(0.3)
    assert metrics["collinear_anchor_count"] == 3
    assert metrics["collinearity_ratio"] == pytest.approx(1.0)
    assert metrics["collinear_orientation"] == "reverse"
    assert metrics["shared_group_ids"] == "G1;G2;G3"


def test_reference_synteny_evidence_calls_two_anchors_supported_and_one_anchor_single(
    tmp_path: Path,
):
    mod = load_module("query_gene_orthologs.py")
    output_root = tmp_path / "query2family"
    synteny_dir = output_root / "synteny"
    synteny_dir.mkdir(parents=True)
    pandas.DataFrame(
        [
            ["Reference_species_REF1", "Reference_species", "upstream", -2, "r1", "G1", 3],
            ["Reference_species_REF1", "Reference_species", "downstream", 2, "r2", "G2", 2],
            ["Other_species_COPY_A", "Other_species", "upstream", -1, "a1", "G1", 3],
            ["Other_species_COPY_A", "Other_species", "downstream", 1, "a2", "G2", 2],
            ["Other_species_COPY_B", "Other_species", "upstream", -1, "b1", "G1", 3],
        ],
        columns=[
            "node_name", "species", "direction", "offset", "neighbor_gene", "group_id", "group_size"
        ],
    ).to_csv(synteny_dir / "FAM_synteny.tsv", sep="\t", index=False)
    columns = [
        {
            "column_order": 1,
            "family_id": "FAM",
            "family_order": 1,
            "reference_species": "Reference_species",
            "cds_fasta_id": "Reference_species_REF1",
            "gene_id": "REF1",
            "plot_label": "REF1",
            "reference_tip_branch_id": 1,
        }
    ]
    glyphs = [
        {
            "species": "Reference_species",
            "family_id": "FAM",
            "family_order": 1,
            "reference_species": "Reference_species",
            "relation": "specific",
            "reference_cds_fasta_ids": "Reference_species_REF1",
            "reference_gene_ids": "REF1",
            "reference_gene_count": 1,
            "copy_number": 1,
            "gene_ids": "Reference_species_REF1",
            "start_order": 1,
            "end_order": 1,
            "is_contiguous": 1,
            "lane_index": 1,
            "lane_count": 1,
        },
        {
            "species": "Other_species",
            "family_id": "FAM",
            "family_order": 1,
            "reference_species": "Reference_species",
            "relation": "specific",
            "reference_cds_fasta_ids": "Reference_species_REF1",
            "reference_gene_ids": "REF1",
            "reference_gene_count": 1,
            "copy_number": 2,
            "gene_ids": "Other_species_COPY_A;Other_species_COPY_B",
            "start_order": 1,
            "end_order": 1,
            "is_contiguous": 1,
            "lane_index": 1,
            "lane_count": 1,
        },
    ]

    evidence = mod.collect_reference_synteny_evidence(
        mod.GeneFamilyOutputStore(output_root), columns, glyphs
    )
    by_candidate = {row["candidate_cds_fasta_id"]: row for row in evidence}

    assert by_candidate["Reference_species_REF1"]["synteny_status"] == "reference_self"
    assert by_candidate["Other_species_COPY_A"]["synteny_status"] == "supported"
    assert by_candidate["Other_species_COPY_A"]["shared_anchor_count"] == 2
    assert by_candidate["Other_species_COPY_A"]["synteny_window_radius"] == 2
    assert by_candidate["Other_species_COPY_B"]["synteny_status"] == "single_anchor"
    assert by_candidate["Other_species_COPY_B"]["shared_anchor_count"] == 1


def test_reference_ufboot_evidence_uses_nonroot_speciation_mrca_branch(
    tmp_path: Path,
):
    mod = load_module("query_gene_orthologs.py")
    output_root = tmp_path / "query2family"
    stat_dir = output_root / "stat_branch"
    stat_dir.mkdir(parents=True)
    rows = [
        stat_row(0, -1, 1, 4, "S", "root"),
        stat_row(1, 0, 2, 3, "S", "orthology_mrca", support=0.93),
        stat_row(2, 1, -1, -1, "L", "Reference_species_REF1", "Reference_species"),
        stat_row(3, 1, -1, -1, "L", "Other_species_COPY", "Other_species"),
        stat_row(4, 0, -1, -1, "L", "Outgroup_species_COPY", "Outgroup_species"),
    ]
    pandas.DataFrame(rows).to_csv(
        stat_dir / "FAM_stat.branch.tsv", sep="\t", index=False
    )
    columns = [
        {
            "column_order": 1,
            "family_id": "FAM",
            "family_order": 1,
            "reference_species": "Reference_species",
            "cds_fasta_id": "Reference_species_REF1",
            "gene_id": "REF1",
            "plot_label": "REF1",
            "reference_tip_branch_id": 2,
        }
    ]
    glyphs = [
        {
            "species": "Reference_species",
            "family_id": "FAM",
            "family_order": 1,
            "reference_species": "Reference_species",
            "relation": "specific",
            "reference_cds_fasta_ids": "Reference_species_REF1",
            "copy_number": 1,
            "gene_ids": "Reference_species_REF1",
            "start_order": 1,
            "end_order": 1,
            "lane_index": 1,
            "lane_count": 1,
        },
        {
            "species": "Other_species",
            "family_id": "FAM",
            "family_order": 1,
            "reference_species": "Reference_species",
            "relation": "specific",
            "reference_cds_fasta_ids": "Reference_species_REF1",
            "copy_number": 1,
            "gene_ids": "Other_species_COPY",
            "start_order": 1,
            "end_order": 1,
            "lane_index": 1,
            "lane_count": 1,
        },
    ]

    evidence = mod.collect_reference_ufboot_evidence(
        mod.GeneFamilyOutputStore(output_root), columns, glyphs
    )
    by_candidate = {row["candidate_cds_fasta_id"]: row for row in evidence}

    assert by_candidate["Reference_species_REF1"]["orthology_ufboot_status"] == "reference_self"
    candidate = by_candidate["Other_species_COPY"]
    assert candidate["orthology_mrca_branch_id"] == 1
    assert candidate["orthology_mrca_event"] == "S"
    assert candidate["decisive_branch_ufboot"] == pytest.approx(93)
    assert candidate["ufboot_support_source"] == "support_unrooted"
    assert candidate["orthology_ufboot_status"] == "evaluated"
    assert candidate["orthology_ufboot_unavailable_reason"] == ""


def test_reference_ufboot_prefers_explicit_generax_support_and_keeps_one_percent(
    tmp_path: Path,
):
    mod = load_module("query_gene_orthologs.py")
    output_root = tmp_path / "query2family"
    stat_dir = output_root / "stat_branch"
    stat_dir.mkdir(parents=True)
    rows = [
        stat_row(0, -1, 1, 4, "S", "root"),
        stat_row(
            1,
            0,
            2,
            3,
            "S",
            "orthology_mrca",
            support=0.93,
            generax_support=1,
        ),
        stat_row(2, 1, -1, -1, "L", "Reference_species_REF1", "Reference_species"),
        stat_row(3, 1, -1, -1, "L", "Other_species_COPY", "Other_species"),
        stat_row(4, 0, -1, -1, "L", "Outgroup_species_COPY", "Outgroup_species"),
    ]
    pandas.DataFrame(rows).to_csv(
        stat_dir / "FAM_stat.branch.tsv", sep="\t", index=False
    )
    columns = [
        {
            "column_order": 1,
            "family_id": "FAM",
            "family_order": 1,
            "reference_species": "Reference_species",
            "cds_fasta_id": "Reference_species_REF1",
            "gene_id": "REF1",
            "plot_label": "REF1",
            "reference_tip_branch_id": 2,
        }
    ]
    glyphs = [
        {
            "species": "Other_species",
            "family_id": "FAM",
            "family_order": 1,
            "reference_species": "Reference_species",
            "relation": "specific",
            "reference_cds_fasta_ids": "Reference_species_REF1",
            "copy_number": 1,
            "gene_ids": "Other_species_COPY",
            "start_order": 1,
            "end_order": 1,
            "lane_index": 1,
            "lane_count": 1,
        }
    ]

    evidence = mod.collect_reference_ufboot_evidence(
        mod.GeneFamilyOutputStore(output_root), columns, glyphs
    )

    assert evidence[0]["ufboot_support_source"] == "support_generax_ufboot"
    assert evidence[0]["decisive_branch_ufboot"] == pytest.approx(1)


def test_reference_ufboot_rejects_multiple_mrca_branches_within_one_glyph(
    tmp_path: Path,
):
    mod = load_module("query_gene_orthologs.py")
    output_root = tmp_path / "query2family"
    stat_dir = output_root / "stat_branch"
    stat_dir.mkdir(parents=True)
    rows = [
        stat_row(0, -1, 1, 5, "S", "root"),
        stat_row(1, 0, 3, 2, "S", "outer_mrca", support=87),
        stat_row(2, 1, 4, 6, "S", "inner_mrca", support=96),
        stat_row(3, 1, -1, -1, "L", "Other_species_COPY_A", "Other_species"),
        stat_row(4, 2, -1, -1, "L", "Reference_species_REF1", "Reference_species"),
        stat_row(5, 0, -1, -1, "L", "Outgroup_species_COPY", "Outgroup_species"),
        stat_row(6, 2, -1, -1, "L", "Other_species_COPY_B", "Other_species"),
    ]
    pandas.DataFrame(rows).to_csv(
        stat_dir / "FAM_stat.branch.tsv", sep="\t", index=False
    )
    columns = [
        {
            "column_order": 1,
            "family_id": "FAM",
            "family_order": 1,
            "reference_species": "Reference_species",
            "cds_fasta_id": "Reference_species_REF1",
            "gene_id": "REF1",
            "plot_label": "REF1",
            "reference_tip_branch_id": 4,
        }
    ]
    glyphs = [
        {
            "species": "Other_species",
            "family_id": "FAM",
            "family_order": 1,
            "reference_species": "Reference_species",
            "relation": "specific",
            "reference_cds_fasta_ids": "Reference_species_REF1",
            "copy_number": 2,
            "gene_ids": "Other_species_COPY_A;Other_species_COPY_B",
            "start_order": 1,
            "end_order": 1,
            "lane_index": 1,
            "lane_count": 1,
        }
    ]

    with pytest.raises(
        ValueError,
        match="do not share one orthology-defining speciation branch",
    ):
        mod.collect_reference_ufboot_evidence(
            mod.GeneFamilyOutputStore(output_root), columns, glyphs
        )


def test_stat_branch_duplicate_branch_id_is_a_hard_error():
    mod = load_module("query_gene_orthologs.py")
    rows = [
        stat_row(0, -1, -1, -1, "L", "root"),
        stat_row(0, -1, -1, -1, "L", "duplicate"),
    ]

    with pytest.raises(ValueError, match="duplicate branch_id: 0"):
        mod.build_tree_index(rows)


def test_stat_branch_missing_parent_is_a_hard_error():
    mod = load_module("query_gene_orthologs.py")
    rows = [
        stat_row(0, -1, -1, -1, "L", "root"),
        stat_row(1, 99, -1, -1, "L", "orphan"),
    ]

    with pytest.raises(ValueError, match="parent does not exist"):
        mod.build_tree_index(rows)


def test_stat_branch_explicit_children_must_match_parent_links():
    mod = load_module("query_gene_orthologs.py")
    rows = [
        stat_row(0, -1, 1, -1, "S", "root"),
        stat_row(1, 0, -1, -1, "L", "left"),
        stat_row(2, 0, -1, -1, "L", "right"),
    ]

    with pytest.raises(ValueError, match="child columns disagree with parent links"):
        mod.build_tree_index(rows)


def test_stat_branch_without_explicit_child_columns_uses_parent_links():
    mod = load_module("query_gene_orthologs.py")
    rows = [
        {"branch_id": 2, "parent": 0, "so_event": "L", "node_name": "right"},
        {"branch_id": 0, "parent": -1, "so_event": "S", "node_name": "root"},
        {"branch_id": 1, "parent": 0, "so_event": "L", "node_name": "left"},
    ]

    by_id, children, root = mod.build_tree_index(rows)

    assert root == 0
    assert children[0] == [1, 2]
    assert mod.depth_first_tip_order(by_id, children, root) == [1, 2]


def test_stat_branch_disconnected_cycle_is_a_hard_error():
    mod = load_module("query_gene_orthologs.py")
    rows = [
        stat_row(0, -1, -1, -1, "L", "root"),
        stat_row(1, 2, 2, -1, "D", "cycle_a"),
        stat_row(2, 1, 1, -1, "D", "cycle_b"),
    ]

    with pytest.raises(ValueError, match="disconnected from its root"):
        mod.build_tree_index(rows)


def test_family_file_selection_is_deduplicated_without_reordering(tmp_path: Path):
    mod = load_module("query_gene_orthologs.py")
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    for family_id in ("AHA", "AFL", "YABBY"):
        (query_dir / family_id).write_text("query\n", encoding="utf-8")
    family_file = tmp_path / "families.tsv"
    family_file.write_text(
        "family_id\nYABBY\nAFL\nYABBY\nAHA\nAFL\n",
        encoding="utf-8",
    )

    assert mod.read_family_ids(query_dir, family_file) == ["YABBY", "AFL", "AHA"]


def test_reference_species_row_must_remain_one_copy_per_reference_gene():
    mod = load_module("query_gene_orthologs.py")
    rows = [
        stat_row(0, -1, 1, 2, "S", "incorrect_speciation"),
        stat_row(1, 0, -1, -1, "L", "Reference_species_GENE1", "Reference_species"),
        stat_row(2, 0, -1, -1, "L", "Reference_species_GENE2", "Reference_species"),
    ]

    with pytest.raises(ValueError, match="one-to-one identity row"):
        mod.collect_family_orthologs(
            rows=rows,
            cds_fasta_ids=["Reference_species_GENE1", "Reference_species_GENE2"],
            reference_species="Reference_species",
            family_id="BAD",
            family_order=1,
            first_column_order=1,
        )


def test_duplicate_cds_fasta_ids_are_a_hard_error(tmp_path: Path):
    mod = load_module("query_gene_orthologs.py")
    output_root = tmp_path / "query2family"
    cds_dir = output_root / "cds_fasta"
    cds_dir.mkdir(parents=True)
    (cds_dir / "DUP_cds.fasta").write_text(
        ">Reference_species_GENE1\nATG\n"
        ">Reference_species_GENE1\nATG\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="duplicate sequence IDs"):
        mod.read_family_cds_fasta_ids(
            mod.GeneFamilyOutputStore(output_root),
            "DUP",
        )


def test_query_gene_basis_coalesces_shared_tips_and_retains_query_mapping():
    mod = load_module("query_gene_orthologs.py")
    rows = [
        stat_row(0, -1, 1, 2, "S", "root"),
        stat_row(1, 0, -1, -1, "L", "Other_species_COPY", "Other_species"),
        stat_row(2, 0, 3, 4, "D", "query_duplication", "n1"),
        stat_row(
            3,
            2,
            -1,
            -1,
            "L",
            "Species_A_GENE_A",
            "Species_A",
            "direct:q1|best:q2",
        ),
        stat_row(4, 2, -1, -1, "L", "Species_B_GENE_B", "Species_B", "best:q3"),
    ]
    definitions = [
        {"query_id": "q1", "query_label": "Query one"},
        {"query_id": "q2", "query_label": "Query two"},
        {"query_id": "q3", "query_label": "Query three"},
    ]

    columns, glyphs, tree_nodes, query_map = mod.collect_family_query_anchor_orthologs(
        rows=rows,
        cds_fasta_ids=["Species_A_GENE_A", "Species_B_GENE_B", "Other_species_COPY"],
        query_definitions=definitions,
        family_id="FAM",
        family_order=1,
        first_column_order=1,
        first_query_order=1,
    )

    assert len(columns) == 2
    assert columns[0]["query_ids"] == "q1;q2"
    assert columns[0]["query_count"] == 2
    assert columns[0]["anchor_source"] == "mixed"
    assert columns[0]["plot_label"] == "q1 (+1)"
    assert [row["column_order"] for row in query_map] == [1, 1, 2]
    assert [row["marker_source"] for row in query_map] == ["direct", "best", "best"]
    assert [row["merged_query_count"] for row in query_map] == [2, 2, 1]
    shared = [row for row in glyphs if row["species"] == "Other_species"]
    assert len(shared) == 1
    assert shared[0]["relation"] == "shared_ancestral"
    assert shared[0]["anchor_query_ids"] == "q1;q2;q3"
    assert {row["reference_species"] for row in tree_nodes} == {"query_gene"}

    output_columns = mod.query_columns_for_output(columns)
    output_glyphs = mod.query_glyphs_for_output(glyphs)
    output_tree = mod.query_tree_for_output(tree_nodes)
    assert set(output_columns[0]) == set(mod.QUERY_COLUMN_FIELDS)
    assert set(output_glyphs[0]) == set(mod.QUERY_GLYPH_FIELDS)
    assert set(output_tree[0]) == set(mod.QUERY_TREE_FIELDS)
    assert output_columns[0]["anchor_cds_fasta_id"] == "Species_A_GENE_A"
    assert output_tree[0]["basis"] == "query_gene"


def test_gene_summary_wires_query_gene_ortholog_tables_and_plot():
    text = CORE_SCRIPT.read_text(encoding="utf-8")
    entrypoint_text = ENTRYPOINT_SCRIPT.read_text(encoding="utf-8")
    registry_text = CONFIG_REGISTRY.read_text(encoding="utf-8")

    assert 'query_gene_orthologs.py"' in text
    assert "query2family_reference_gene_orthologs.columns.tsv" in text
    assert "query2family_reference_gene_orthologs.glyphs.tsv" in text
    assert "query2family_reference_gene_orthologs.tree.tsv" in text
    assert "query2family_reference_gene_orthologs.synteny.tsv" in text
    assert "query2family_reference_gene_orthologs.ufboot.tsv" in text
    assert "query2family_query_gene_orthologs.columns.tsv" in text
    assert "query2family_query_gene_orthologs.glyphs.tsv" in text
    assert "query2family_query_gene_orthologs.tree.tsv" in text
    assert "query2family_query_gene_orthologs.synteny.tsv" in text
    assert "query2family_query_gene_orthologs.ufboot.tsv" in text
    assert "query2family_query_gene_orthologs.query_map.tsv" in text
    assert '--ortholog_column_table="${file_query_columns}"' in text
    assert '--ortholog_glyph_table="${file_query_glyphs}"' in text
    assert '--ortholog_tree_table="${file_query_tree}"' in text
    assert '--ortholog_synteny_table="${file_query_synteny}"' in text
    assert '--ortholog_ufboot_table="${file_query_ufboot}"' in text
    assert '--species_mapping_tree="${file_species_mapping_tree}"' in text
    assert '--evidence_layout="${presence_absence_evidence_layout}"' in text
    assert 'presence_absence_evidence_layout must be "band", "rail", "glyph", or "off"' in text
    assert '--ortholog_basis=query_gene' in text
    assert 'presence_absence_ortholog_basis must be "reference_species", "query_gene", or "both"' in text
    assert "resolve_presence_absence_species_mapping_tree" in text
    assert '--reference_species "${reference_species_resolved}"' in text
    assert '--reference_species="${reference_species_resolved}"' in text
    assert 'GG_COMMON_REFERENCE_SPECIES:-auto' in text
    assert "query2family_reference_gene_orthologs.pdf" in text
    assert 'presence_absence_ortholog_basis="${presence_absence_ortholog_basis:-reference_species}"' in entrypoint_text
    assert "presence_absence_ortholog_basis" in registry_text

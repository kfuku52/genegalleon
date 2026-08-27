from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas
import pytest

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
CORE_SCRIPT = Path(__file__).resolve().parents[1] / "core" / "gg_gene_summary_core.sh"
ENTRYPOINT_SCRIPT = Path(__file__).resolve().parents[1] / "gg_gene_summary_entrypoint.sh"
CONFIG_REGISTRY = SUPPORT_DIR / "gg_entrypoint_config_vars.sh"


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

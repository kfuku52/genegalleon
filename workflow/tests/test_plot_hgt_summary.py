import subprocess
import sys
from pathlib import Path

import pandas
import pytest

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT_DIR))

import plot_hgt_summary as plotter  # noqa: E402

SCRIPT_PATH = SUPPORT_DIR / "plot_hgt_summary.py"


def test_transfer_tree_preserves_numeric_internal_branch_labels(tmp_path):
    path = tmp_path / "tree.nwk"
    path.write_text("((A:1,B:1)42:1,C:2)99;")
    tree, _, _, labels = plotter.load_species_tree_layout(str(path))
    assert labels["42"] == "42"
    assert labels["99"] == "99"
    edges = plotter.build_transfer_edge_table(pandas.DataFrame({"generax_transfer": ["Y@C@42"]}), labels)
    measured = plotter.add_transfer_distances(edges, tree, labels, 0)
    assert measured.iloc[0].mapped_to_species_tree == 1
    assert measured.iloc[0].phylogenetic_distance == 3


def test_plot_hgt_summary_generates_overview_and_taxonomy_flow_pdfs(tmp_path: Path):
    branch_tsv = tmp_path / "hgt_branch_candidates.tsv"
    gene_tsv = tmp_path / "hgt_gene_candidates.tsv"
    overview_pdf = tmp_path / "hgt_branch_overview.pdf"
    flow_pdf = tmp_path / "hgt_taxonomy_flow.pdf"
    readme_md = tmp_path / "README.md"

    pandas.DataFrame(
        [
            {
                "orthogroup": "OG0001",
                "branch_id": 3,
                "candidate_gene_count": 2,
                "matched_leaf_count": 2,
                "besthit_gene_count": 2,
                "besthit_taxid_count": 0,
                "besthit_taxonomy_method": "name_heuristic",
                "besthit_same_superkingdom_fraction": 0.0,
                "besthit_lca_rank_mode": "genus_mismatch",
                "intron_support_fraction": 1.0,
                "expression_measured_fraction": 1.0,
                "clade_min_expression_pearsoncor": 0.4,
                "synteny_support_fraction": 0.5,
                "synteny_mean_support_score": 0.6,
                "contamination_incompatible_fraction": 0.0,
                "contamination_top_lca_sciname": "",
            },
            {
                "orthogroup": "OG0002",
                "branch_id": 7,
                "candidate_gene_count": 1,
                "matched_leaf_count": 1,
                "besthit_gene_count": 1,
                "besthit_taxid_count": 0,
                "besthit_taxonomy_method": "name_heuristic",
                "besthit_same_superkingdom_fraction": 1.0,
                "besthit_lca_rank_mode": "species",
                "intron_support_fraction": 0.0,
                "expression_measured_fraction": 1.0,
                "clade_min_expression_pearsoncor": 0.0,
                "synteny_support_fraction": 0.0,
                "synteny_mean_support_score": 0.0,
                "contamination_incompatible_fraction": 1.0,
                "contamination_top_lca_sciname": "Escherichia coli",
            },
        ]
    ).to_csv(branch_tsv, sep="\t", index=False)
    pandas.DataFrame(
        [
            {
                "orthogroup": "OG0001",
                "gene_id": "geneA",
                "gene_taxon": "Arabidopsis thaliana",
                "candidate_branch_count": 1,
                "candidate_branch_ids": "3",
                "besthit_accession": "P1",
                "besthit_organism": "Escherichia coli",
                "besthit_taxid": "",
                "besthit_taxonomy_method": "name_heuristic",
                "besthit_lca_rank": "genus_mismatch",
                "besthit_same_superkingdom": 0,
                "intron_supported": True,
                "expression_measured": True,
                "synteny_support_score": 0.8,
                "contamination_lca_taxid": "",
                "contamination_lca_sciname": "",
                "contamination_is_compatible_lineage": True,
            },
            {
                "orthogroup": "OG0002",
                "gene_id": "geneB",
                "gene_taxon": "Arabidopsis thaliana",
                "candidate_branch_count": 1,
                "candidate_branch_ids": "7",
                "besthit_accession": "P2",
                "besthit_organism": "Bacillus subtilis",
                "besthit_taxid": "",
                "besthit_taxonomy_method": "name_heuristic",
                "besthit_lca_rank": "genus_mismatch",
                "besthit_same_superkingdom": 0,
                "intron_supported": False,
                "expression_measured": True,
                "synteny_support_score": 0.0,
                "contamination_lca_taxid": "562",
                "contamination_lca_sciname": "Escherichia coli",
                "contamination_is_compatible_lineage": False,
            },
        ]
    ).to_csv(gene_tsv, sep="\t", index=False)

    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--branch_tsv",
            str(branch_tsv),
            "--gene_tsv",
            str(gene_tsv),
            "--overview_pdf",
            str(overview_pdf),
            "--taxonomy_flow_pdf",
            str(flow_pdf),
            "--flow_rank",
            "phylum",
            "--flow_max_categories",
            "8",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert overview_pdf.exists() and overview_pdf.stat().st_size > 0
    assert flow_pdf.exists() and flow_pdf.stat().st_size > 0
    assert readme_md.exists()
    readme_text = readme_md.read_text(encoding="utf-8")
    assert "`Cand`" in readme_text
    assert "`TopContam`" in readme_text


def test_plot_taxonomy_flow_uses_precomputed_phylum_columns_without_taxonomy_db(tmp_path, monkeypatch):
    captured = {}

    def capture_flow_frame(count_df, left_col, right_col, max_categories):
        captured["frame"] = count_df.copy()
        return count_df

    monkeypatch.setattr(plotter, "collapse_to_top_categories", capture_flow_frame)

    gene_df = pandas.DataFrame(
        [
            {
                "gene_taxon": "Arabidopsis thaliana",
                "besthit_organism": "Escherichia coli",
                "besthit_taxid": "562",
                "recipient_phylum": "Streptophyta",
                "donor_phylum": "Pseudomonadota",
            },
            {
                "gene_taxon": "Arabidopsis thaliana",
                "besthit_organism": "Bacillus subtilis",
                "besthit_taxid": "1423",
                "recipient_phylum": "Streptophyta",
                "donor_phylum": "Bacillota",
            },
        ]
    )
    flow_pdf = tmp_path / "flow.pdf"
    plotter.plot_taxonomy_flow(
        gene_df=gene_df,
        out_pdf=str(flow_pdf),
        resolver=plotter.TaxonomyResolver(""),
        preferred_rank="phylum",
        max_categories=8,
    )

    assert flow_pdf.exists() and flow_pdf.stat().st_size > 0
    assert set(captured["frame"]["recipient_label"]) == {"Streptophyta"}
    assert set(captured["frame"]["besthit_label"]) == {"Pseudomonadota", "Bacillota"}


def test_plot_transfer_tree_aggregates_directed_counts_and_limits_display_edges(tmp_path: Path):
    species_tree = tmp_path / "species_tree.nwk"
    species_tree.write_text("((A:1,B:1)n1:1,(C:1,D:1)n2:1)n0;\n", encoding="utf-8")
    branch_df = pandas.DataFrame(
        [
            {"orthogroup": "OG0001", "generax_transfer": "Y@A@B"},
            {"orthogroup": "OG0002", "generax_transfer": "Y@A@B"},
            {"orthogroup": "OG0003", "generax_transfer": "Y@B@C"},
            {"orthogroup": "OG0004", "generax_transfer": "Y@n1@D"},
            {"orthogroup": "OG0005", "generax_transfer": "Y@missing@A"},
            {"orthogroup": "OG0006", "generax_transfer": "N@A@C"},
        ]
    )
    transfer_pdf = tmp_path / "plots" / "hgt_transfer_tree.pdf"
    edges_tsv = tmp_path / "plots" / "hgt_transfer_edges.tsv"

    plotter.plot_transfer_tree(
        branch_df=branch_df,
        out_pdf=str(transfer_pdf),
        species_tree_path=str(species_tree),
        edges_tsv=str(edges_tsv),
        max_edges=2,
    )

    assert transfer_pdf.exists() and transfer_pdf.stat().st_size > 0
    edges = pandas.read_csv(edges_tsv, sep="\t")
    assert edges.shape[0] == 4
    assert set(edges.columns) == set(plotter.TRANSFER_EDGE_COLUMNS)
    top_edge = edges.loc[(edges["donor_node"] == "A") & (edges["recipient_node"] == "B")].iloc[0]
    assert int(top_edge["hgt_event_count"]) == 2
    assert int(top_edge["orthogroup_count"]) == 2
    assert int(top_edge["mapped_to_species_tree"]) == 1
    assert int(edges["displayed"].sum()) == 2
    missing_edge = edges.loc[edges["donor_node"] == "missing"].iloc[0]
    assert int(missing_edge["mapped_to_species_tree"]) == 0
    assert int(missing_edge["displayed"]) == 0
    assert top_edge["phylogenetic_distance"] == 2
    assert top_edge["distance_metric"] == "branch_length"
    assert pandas.isna(missing_edge["phylogenetic_distance"])
    distant = edges.loc[edges.donor_node.eq("B")].iloc[0]
    assert distant.phylogenetic_distance == 4
    assert distant.displayed == 1
    assert distant.selection_reason == "distance"


def test_transfer_distance_topology_fallback_and_ancestor(tmp_path):
    path = tmp_path / "tree.nwk"
    path.write_text("((A,B)n1,C)n0;")
    tree, _, _, labels = plotter.load_species_tree_layout(str(path))
    frame = pandas.DataFrame({"generax_transfer": ["Y@n1@A", "Y@A@C", "Y@C@A", "Y@A@A"]})
    edges = plotter.add_transfer_distances(
        plotter.build_transfer_edge_table(frame, labels), tree, labels, 0)
    assert set(edges.distance_metric) == {"topology_edges"}
    values = {(r.donor_node, r.recipient_node): r.phylogenetic_distance for r in edges.itertuples()}
    assert values == {("n1", "A"): 1, ("A", "C"): 3, ("C", "A"): 3, ("A", "A"): 0}
    assert edges.displayed.sum() == 4


def test_species_branch_anchors_use_incoming_horizontal_segment(tmp_path):
    path = tmp_path / "tree.nwk"
    path.write_text("((A:2,B:0)n1:1,C:3)n0;")
    tree, x, y, _ = plotter.load_species_tree_layout(str(path))
    anchors = plotter.species_branch_anchors(tree, x, y)
    nodes = {c.name: c for c in tree.find_clades()}
    for parent, child in [("n0", "n1"), ("n1", "A"), ("n1", "B")]:
        p, c = id(nodes[parent]), id(nodes[child])
        assert anchors[c] == ((x[p] + x[c]) / 2, y[c])
    assert anchors[id(tree.root)] == (x[id(tree.root)] / 2, y[id(tree.root)])


def test_reciprocal_completion_and_shared_curve(tmp_path):
    path = tmp_path / "tree.nwk"
    path.write_text("(A:1,B:1,C:1)n0;")
    tree, x, y, labels = plotter.load_species_tree_layout(str(path))
    frame = pandas.DataFrame({"generax_transfer": ["Y@A@B"] * 8 + ["Y@B@A", "Y@A@C"]})
    edges = plotter.add_transfer_distances(plotter.build_transfer_edge_table(frame, labels), tree, labels, 1)
    shown = edges.loc[edges.displayed.eq(1)]
    assert len(shown) == 2
    assert shown.loc[shown.donor_node.eq("B"), "selection_reason"].iloc[0] == "reciprocal"
    groups = plotter.transfer_connections(shown)
    assert len(groups) == 1
    assert {r.recipient_node: r.hgt_event_count for r in groups[0][1]} == {"B": 8, "A": 1}
    plt, *_ = plotter.get_pyplot()
    fig, ax = plt.subplots()
    a, b = (0.1, 0.2), (0.8, 0.9)
    left, right = plotter.transfer_half_paths(ax, a, b)
    import numpy as np
    np.testing.assert_allclose(left.vertices[0], right.vertices[0])
    np.testing.assert_allclose(left.vertices[-1], a)
    np.testing.assert_allclose(right.vertices[-1], b)
    plotter.draw_species_tree(ax, tree, x, y, set())
    assert {t.get_text() for t in ax.texts} == {"A", "B", "C", "n0"}
    anchors = plotter.species_branch_anchors(tree, x, y)
    nodes = {c.name: c for c in tree.find_clades()}
    for text in ax.texts:
        node = nodes[text.get_text()]
        assert text.xy == ((x[id(node)], y[id(node)]) if node.is_terminal() else anchors[id(node)])
        assert text.get_position() == ((4, 0) if node.is_terminal() else (0, 3))
        assert text.get_fontsize() == 8
    plt.close(fig)


def test_transfer_traits_preserve_zero_unknown_and_tip_order(tmp_path):
    path = tmp_path / "traits.tsv"
    path.write_text("species\tgall\nB\t0\nA species\t1\nC\tNA\n")
    traits = plotter.read_transfer_traits(str(path))
    assert traits.loc["A_species", "gall"] == 1
    assert traits.loc["B", "gall"] == 0
    assert pandas.isna(traits.loc["C", "gall"])
    tree_path = tmp_path / "tree.nwk"
    tree_path.write_text("(A_species:1,B:1,C:1,D:1)n0;")
    tree, _, y, _ = plotter.load_species_tree_layout(str(tree_path))
    plt, *_ = plotter.get_pyplot()
    fig, ax = plt.subplots()
    plotter.draw_transfer_traits(ax, tree, y, traits)
    values = [t.get_text() for t in ax.texts][1:-1]
    assert values == ["1", "0", "NA", "NA"]
    assert all(t.get_fontsize() == 8 for t in ax.texts)
    plt.close(fig)
    plotter.plot_transfer_tree(pandas.DataFrame({"generax_transfer": ["Y@A_species@B"]}),
                              str(tmp_path / "trait_plot.pdf"), str(tree_path),
                              species_trait_path=str(path))
    assert (tmp_path / "trait_plot.pdf").stat().st_size > 0
    path.write_text("species\tgall\nA species\t1\nA_species\t0\n")
    with pytest.raises(ValueError, match="unique"):
        plotter.read_transfer_traits(str(path))

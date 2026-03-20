from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
import types


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "orthogroup_statistics.py"


def load_module():
    if "kftools" not in sys.modules:
        pkg = types.ModuleType("kftools")
        pkg.__path__ = []
        sys.modules["kftools"] = pkg
    if "kftools.kfog" not in sys.modules:
        sys.modules["kftools.kfog"] = types.ModuleType("kftools.kfog")
    if "kftools.kfphylo" not in sys.modules:
        sys.modules["kftools.kfphylo"] = types.ModuleType("kftools.kfphylo")
    spec = spec_from_file_location("orthogroup_statistics", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_load_fimo_hits_parses_modern_fimo_tsv(tmp_path):
    mod = load_module()
    infile = tmp_path / "fimo.tsv"
    infile.write_text(
        (
            "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n"
            "M1\tALT1\tgeneA\t10\t15\t+\t11.5\t1e-4\t1e-3\tACGTAC\n"
            "M2\t.\tgeneB\t20\t24\t-\t9.3\t1e-3\t2e-2\tTTGGCC\n"
        ),
        encoding="utf-8",
    )

    df = mod.load_fimo_hits(str(infile))
    assert df.shape[0] == 2
    assert df.columns.tolist() == [
        "motif_id",
        "motif_alt_id",
        "sequence_name",
        "start",
        "stop",
        "strand",
        "q-value",
    ]
    assert df.loc[0, "motif_alt_id"] == "ALT1"
    # "." should fall back to motif_id for consistent downstream grouping.
    assert df.loc[1, "motif_alt_id"] == "M2"


def test_orthogroup_statistics_preserves_qualified_species_prefixes():
    mod = load_module()

    assert mod.extract_species_label("Dictyostelium_discoideum_cf_geneA") == "Dictyostelium_discoideum_cf"
    assert mod.extract_species_label("Bacillus_subtilis_subsp_subtilis_geneB") == "Bacillus_subtilis_subsp_subtilis"
    assert mod.scientific_name_from_label("Dictyostelium_discoideum_cf_geneA") == "Dictyostelium cf. discoideum"


def test_load_fimo_hits_parses_legacy_fimo_txt_header(tmp_path):
    mod = load_module()
    infile = tmp_path / "fimo.txt"
    infile.write_text(
        (
            "#pattern name\tsequence name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched sequence\n"
            "MA0123.1\tgeneC\t30\t39\t+\t7.0\t1e-3\t4e-2\tACGTTGCAAA\n"
            "MA0567.1\tgeneD\t40\t47\t-\t6.1\t2e-3\t5e-2\tTTAACCGG\n"
        ),
        encoding="utf-8",
    )

    df = mod.load_fimo_hits(str(infile))
    assert df.shape[0] == 2
    assert set(df["sequence_name"].tolist()) == {"geneC", "geneD"}
    # legacy format has no alt id; should be mirrored from motif_id.
    assert df.loc[df["sequence_name"] == "geneC", "motif_alt_id"].iloc[0] == "MA0123.1"


def test_load_scm_intron_branch_table_maps_dated_tree_labels_to_branch_ids(tmp_path):
    mod = load_module()
    rooted_tree = tmp_path / "rooted.nwk"
    dated_tree = tmp_path / "dated.nwk"
    scm_intron = tmp_path / "scm.tsv"

    rooted_tree.write_text("((A:1,B:1)n1:1,(C:1,D:1)n2:1)n0;\n", encoding="utf-8")
    dated_tree.write_text("((A:1,B:1)n11:1,(C:1,D:1)n12:1)n99;\n", encoding="utf-8")
    scm_intron.write_text(
        (
            "leaf\tnum_intron\tintron_present\tintron_absent\n"
            "n99\t2\t1\t0\n"
            "n11\t1\t0.7\t0.3\n"
            "n12\t0\t0.2\t0.8\n"
            "A\t3\t1\t0\n"
            "B\t0\t0\t1\n"
            "C\t4\t1\t0\n"
            "D\t0\t0\t1\n"
            "n404\t9\t1\t0\n"
        ),
        encoding="utf-8",
    )

    df_scm = mod.load_scm_intron_branch_table(str(scm_intron), str(dated_tree))
    assert "branch_id" in df_scm.columns
    assert df_scm["branch_id"].isna().sum() == 0
    assert "n404" not in set(df_scm["node_name"].tolist())

    rooted = mod._ensure_branch_ids(mod.new_tree(str(rooted_tree), format=1))
    df_branch = mod.pandas.DataFrame(
        {
            "branch_id": [mod._get_node_label(node) for node in rooted.traverse()],
            "node_name": [node.name for node in rooted.traverse()],
        }
    )
    merged = mod.pandas.merge(
        df_branch,
        df_scm.drop(columns=["node_name"]),
        on="branch_id",
        how="outer",
    )

    assert merged["branch_id"].isna().sum() == 0
    assert merged.loc[merged["node_name"] == "n0", "num_intron"].iloc[0] == 2
    assert merged.loc[merged["node_name"] == "n1", "num_intron"].iloc[0] == 1
    assert merged.loc[merged["node_name"] == "n2", "num_intron"].iloc[0] == 0


def test_flatten_trait_variable_stats_builds_tree_info_keys():
    mod = load_module()
    df = mod.pandas.DataFrame(
        [
            {"trait": "salt", "variable": "omega", "coef": 1.5, "pvalue": 0.01},
            {"trait": "cold", "variable": "omega", "coef": -0.5, "pvalue": 0.20},
        ]
    )

    out = mod.flatten_trait_variable_stats(df, "pgls_geneTree_")

    assert out == {
        "pgls_geneTree_coef_salt_omega": 1.5,
        "pgls_geneTree_coef_cold_omega": -0.5,
        "pgls_geneTree_pvalue_salt_omega": 0.01,
        "pgls_geneTree_pvalue_cold_omega": 0.20,
    }


def test_load_synteny_summary_aggregates_per_node_metrics(tmp_path):
    mod = load_module()
    infile = tmp_path / "synteny.tsv"
    infile.write_text(
        (
            "node_name\tspecies\tdirection\toffset\tneighbor_gene\tgroup_id\tgroup_size\n"
            "geneA\tsp1\tupstream\t-1\tn1\tG1\t5\n"
            "geneA\tsp1\tdownstream\t2\tn2\tG2\t2\n"
            "geneB\tsp2\tupstream\t-3\tn3\tG1\t5\n"
            "geneB\tsp2\tdownstream\t1\tn4\tG3\t1\n"
            # duplicate row should not inflate counts
            "geneA\tsp1\tupstream\t-1\tn1\tG1\t5\n"
        ),
        encoding="utf-8",
    )

    out = mod.load_synteny_summary(str(infile))
    assert set(out["node_name"].tolist()) == {"geneA", "geneB"}
    row_a = out.loc[out["node_name"] == "geneA"].iloc[0]
    row_b = out.loc[out["node_name"] == "geneB"].iloc[0]

    assert row_a["synteny_edge_count"] == 2
    assert row_a["synteny_group_count"] == 2
    assert row_a["synteny_shared_group_count"] == 1
    assert row_a["synteny_shared_edge_count"] == 1
    assert row_a["synteny_shared_upstream_edge_count"] == 1
    assert row_a["synteny_shared_downstream_edge_count"] == 0
    assert row_a["synteny_min_abs_offset_shared"] == 1
    assert row_a["synteny_max_group_size"] == 5
    assert row_a["synteny_max_shared_group_tip_count"] == 2
    assert row_a["synteny_support_score"] == 0.5

    assert row_b["synteny_edge_count"] == 2
    assert row_b["synteny_group_count"] == 2
    assert row_b["synteny_shared_group_count"] == 1
    assert row_b["synteny_shared_edge_count"] == 1
    assert row_b["synteny_shared_upstream_edge_count"] == 1
    assert row_b["synteny_shared_downstream_edge_count"] == 0
    assert row_b["synteny_min_abs_offset_shared"] == 3
    assert row_b["synteny_max_group_size"] == 5
    assert row_b["synteny_max_shared_group_tip_count"] == 2
    assert row_b["synteny_support_score"] == 0.5


def test_load_synteny_summary_returns_empty_for_missing_required_columns(tmp_path):
    mod = load_module()
    infile = tmp_path / "synteny_missing.tsv"
    infile.write_text(
        (
            "node_name\toffset\tneighbor_gene\n"
            "geneA\t-1\tn1\n"
        ),
        encoding="utf-8",
    )
    out = mod.load_synteny_summary(str(infile))
    assert out.empty

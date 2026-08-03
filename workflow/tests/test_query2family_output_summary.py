from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace

import pandas

from workflow.support.gene_family_output_store import archive_completed_outputs, family_context

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "query2family_output_summary.py"


def load_module():
    spec = spec_from_file_location("query2family_output_summary", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_amas(path: Path, no_of_taxa: int):
    pandas.DataFrame(
        [
            {
                "No_of_taxa": no_of_taxa,
                "Alignment_length": 10,
                "Total_matrix_cells": 20,
                "Undetermined_characters": 1,
                "Missing_percent": 5.0,
                "No_variable_sites": 3,
                "Parsimony_informative_sites": 2,
                "GC_content": 0.4,
            }
        ]
    ).to_csv(path, sep="\t", index=False)


def test_extract_query_id_prefers_longest_valid_prefix():
    mod = load_module()
    matchers = mod._query_id_matchers(["A", "A_B", "A_B.extra"])

    assert mod._extract_query_id("A_B_tree_plot.pdf", matchers) == "A_B"
    assert mod._extract_query_id("A_B.extra.csubst_cb_3.tsv", matchers) == "A_B.extra"
    assert mod._extract_query_id("A_uniprot.tsv", matchers) == "A"
    assert mod._extract_query_id("AHA_tree_plot.pdf", matchers) is None


def test_run_writes_query_completion_summary_from_input_order(tmp_path: Path):
    mod = load_module()
    query_gene = tmp_path / "input" / "query_gene"
    query2family = tmp_path / "output" / "query2family"
    query_gene.mkdir(parents=True)
    query2family.mkdir(parents=True)

    (query_gene / "AHA").write_text("gene1\n", encoding="utf-8")
    (query_gene / "2_WOX").write_text("gene2\n", encoding="utf-8")
    (query_gene / ".hidden").write_text("hidden\n", encoding="utf-8")
    (query2family / "tree_plot").mkdir()
    (query2family / "stat_branch").mkdir()
    (query2family / "tmp" / "1_AHA").mkdir(parents=True)
    (query2family / ".hidden_stage").mkdir()
    (query2family / "parameters").mkdir()
    (query2family / "parameters" / "species_tree.nwk").write_text("(A,B);\n", encoding="utf-8")
    (query2family / "tree_plot" / "2_WOX_tree_plot.pdf").write_text("", encoding="utf-8")
    (query2family / "tree_plot" / "UNKNOWN_tree_plot.pdf").write_text("", encoding="utf-8")
    (query2family / "stat_branch" / "AHA_stat.branch.tsv").write_text("", encoding="utf-8")
    (query2family / "stat_branch" / ".DS_Store").write_text("", encoding="utf-8")

    out_tsv = tmp_path / "query2family_summary.tsv"
    mod.run(
        SimpleNamespace(
            dir_query2family=str(query2family),
            dir_query_gene=str(query_gene),
            out=str(out_tsv),
            ncpu=1,
        )
    )

    out = pandas.read_csv(out_tsv, sep="\t", index_col=0)
    assert out.index.tolist() == ["2_WOX", "AHA"]
    assert out.loc["2_WOX", "GG_ARRAY_TASK_ID"] == 1
    assert out.loc["AHA", "GG_ARRAY_TASK_ID"] == 2
    assert ".hidden_stage" not in out.columns
    assert "tmp" not in out.columns
    assert "parameters" not in out.columns
    assert out.loc["2_WOX", "tree_plot"] == 1
    assert out.loc["AHA", "tree_plot"] == 0
    assert out.loc["2_WOX", "stat_branch"] == 0
    assert out.loc["AHA", "stat_branch"] == 1


def test_run_appends_optional_amas_stats_without_extra_rows(tmp_path: Path):
    mod = load_module()
    query_gene = tmp_path / "input" / "query_gene"
    query2family = tmp_path / "output" / "query2family"
    query_gene.mkdir(parents=True)
    (query_gene / "AHA").write_text("gene1\n", encoding="utf-8")
    (query_gene / "WOX").write_text("gene2\n", encoding="utf-8")
    amas_original = query2family / "amas_original"
    amas_original.mkdir(parents=True)
    _write_amas(amas_original / "AHA_amas.original.tsv", 7)
    _write_amas(amas_original / "UNKNOWN_amas.original.tsv", 9)

    out_tsv = tmp_path / "query2family_summary.tsv"
    mod.run(
        SimpleNamespace(
            dir_query2family=str(query2family),
            dir_query_gene=str(query_gene),
            out=str(out_tsv),
            ncpu=1,
        )
    )

    out = pandas.read_csv(out_tsv, sep="\t", index_col=0)
    assert set(out.index.tolist()) == {"AHA", "WOX"}
    assert out.loc["AHA", "No_of_taxa_original"] == 7
    assert pandas.isna(out.loc["WOX", "No_of_taxa_original"])
    assert out.loc["AHA", "amas_original"] == 1


def test_run_reads_completion_and_amas_members_from_zip_shards(tmp_path: Path):
    mod = load_module()
    query_gene = tmp_path / "input" / "query_gene"
    query2family = tmp_path / "output" / "query2family"
    query_gene.mkdir(parents=True)
    (query_gene / "AHA").write_text("gene1\n", encoding="utf-8")
    for subdir, suffix in (
        ("stat_branch", "_stat.branch.tsv"),
        ("stat_tree", "_stat.tree.tsv"),
        ("tree_plot", "_tree_plot.pdf"),
    ):
        path = query2family / subdir / f"AHA{suffix}"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"{subdir}\n", encoding="utf-8")
    amas_path = query2family / "amas_original" / "AHA_amas.original.tsv"
    amas_path.parent.mkdir(parents=True)
    _write_amas(amas_path, 11)
    family_ids, family_from_name = family_context("query2family", query_dir=query_gene)
    archive_completed_outputs(query2family, "query2family", family_ids, family_from_name)

    out_tsv = tmp_path / "query2family_summary.tsv"
    mod.run(
        SimpleNamespace(
            dir_query2family=str(query2family),
            dir_query_gene=str(query_gene),
            out=str(out_tsv),
            ncpu=1,
        )
    )

    out = pandas.read_csv(out_tsv, sep="\t", index_col=0)
    assert out.loc["AHA", "tree_plot"] == 1
    assert out.loc["AHA", "stat_branch"] == 1
    assert out.loc["AHA", "No_of_taxa_original"] == 11
    assert not (query2family / "amas_original").exists()

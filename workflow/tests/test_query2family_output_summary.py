import gzip
import json
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace

import pandas
import pytest

from workflow.support.gene_family_output_store import archive_completed_outputs, family_context

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "query2family_output_summary.py"


def load_module():
    spec = spec_from_file_location("query2family_output_summary", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_alignment_stats(path: Path, no_of_taxa: int):
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


def test_run_appends_optional_alignment_stats_stats_without_extra_rows(tmp_path: Path):
    mod = load_module()
    query_gene = tmp_path / "input" / "query_gene"
    query2family = tmp_path / "output" / "query2family"
    query_gene.mkdir(parents=True)
    (query_gene / "AHA").write_text("gene1\n", encoding="utf-8")
    (query_gene / "WOX").write_text("gene2\n", encoding="utf-8")
    alignment_stats_original = query2family / "alignment_stats_original"
    alignment_stats_original.mkdir(parents=True)
    _write_alignment_stats(alignment_stats_original / "AHA_alignment_stats.original.tsv", 7)
    _write_alignment_stats(alignment_stats_original / "UNKNOWN_alignment_stats.original.tsv", 9)

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
    assert out.loc["AHA", "alignment_stats_original"] == 1


def test_run_reads_completion_and_alignment_stats_members_from_zip_shards(tmp_path: Path):
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
    alignment_stats_path = query2family / "alignment_stats_original" / "AHA_alignment_stats.original.tsv"
    alignment_stats_path.parent.mkdir(parents=True)
    _write_alignment_stats(alignment_stats_path, 11)
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
    assert not (query2family / "alignment_stats_original").exists()


@pytest.mark.parametrize("archive", [False, True])
def test_query_search_settings_report_uses_recorded_and_verified_legacy_values(tmp_path, archive):
    mod = load_module()
    query_gene = tmp_path / "input" / "query_gene"
    root = tmp_path / "output" / "query2family"
    query_gene.mkdir(parents=True)
    (root / "artifact_provenance").mkdir(parents=True)
    (root / "query_aa_fasta").mkdir()
    families = ["AUTO", "CCDC15", "FIXED", "MISSING", "STALE"]
    for query in families:
        (query_gene / query).write_text("gene1\n")
        for subdir, extension in (
            ("stat_branch", "stat.branch.tsv"), ("stat_tree", "stat.tree.tsv"), ("tree_plot", "tree_plot.pdf"),
        ):
            path = root / subdir / f"{query}_{extension}"
            path.parent.mkdir(exist_ok=True)
            path.write_text("complete\n")
        if query == "MISSING":
            continue

        query_path = root / "query_aa_fasta" / f"{query}_query.aa.fa.gz"
        raw = gzip.compress(b">short\n" + b"A" * 40 + b"\n>long\n" + b"L" * 150 + b"\n")
        query_path.write_bytes(raw)
        query_input = dict(label="query_fasta", scope="logical",
                           path=query_path.relative_to(root).as_posix(), sha256=sha256(raw).hexdigest())
        search_output = dict(label="query_blast", scope="logical", path=f"query_blast/{query}_query_blast.tsv",
                             sha256=sha256(query.encode()).hexdigest())
        search_args = dict(method="diamond", diamond_sensitivity="ultra-sensitive", evalue="auto",
                           auto_evalue_cutoffs="80:100,150:10,300:1,inf:0.01")
        diagnostics = {}
        coverage = "0.25"
        if query == "CCDC15":
            diagnostics = dict(effective_evalue="0.01", max_query_aa_length="992")
            coverage = "0.15"
        elif query == "FIXED":
            search_args.update(method="tblastn", evalue="1e-5")
        elif query == "STALE":
            query_path.write_bytes(gzip.compress(b">different\nAA\n"))
        primary_input = search_output.copy()
        if query == "STALE":
            primary_input["sha256"] = "0" * 64
        search = dict(inputs=[query_input], outputs=[search_output], parameters=search_args, diagnostics=diagnostics)
        primary = dict(inputs=[primary_input], parameters=dict(query_blast_coverage=coverage, max_query_hits="5000"))
        for step, manifest in (("query_blast", search), ("primary_fasta", primary)):
            manifest.update(schema_version=1, step=step, family_id=query)
            (root / "artifact_provenance" / f"{query}.{step}.json").write_text(json.dumps(manifest))
    if archive:
        ids, from_name = family_context("query2family", query_dir=query_gene)
        archive_completed_outputs(root, "query2family", ids, from_name)
        assert not (root / "query_aa_fasta" / "AUTO_query.aa.fa.gz").exists()

    path = tmp_path / "query2family_search_settings.tsv"
    mod.run(SimpleNamespace(dir_query2family=str(root), dir_query_gene=str(query_gene), ncpu=1,
                            out=tmp_path / "query2family_summary.tsv", out_search_settings=path))
    out = pandas.read_csv(path, sep="\t", dtype=str, keep_default_na=False, index_col="query")
    assert out.index.tolist() == families
    assert out.loc["AUTO", ["GG_ARRAY_TASK_ID", "effective_evalue", "evalue_source", "max_query_aa_length",
                            "min_query_coverage", "max_query_hits"]].tolist() == [
        "1", "10", "derived_from_query_fasta", "150", "0.25", "5000"]
    assert out.loc["CCDC15", ["search_method", "diamond_sensitivity", "requested_evalue", "effective_evalue",
                              "evalue_source", "max_query_aa_length", "min_query_coverage"]].tolist() == [
        "diamond", "ultra-sensitive", "auto", "0.01", "recorded", "992", "0.15"]
    assert out.loc["FIXED", ["search_method", "diamond_sensitivity", "effective_evalue", "evalue_source"]].tolist() == [
        "tblastn", "NA", "1e-5", "configured"]
    assert out.loc["STALE", ["effective_evalue", "min_query_coverage"]].tolist() == ["NA", "NA"]
    assert out.loc["MISSING", "search_method"] == "NA"

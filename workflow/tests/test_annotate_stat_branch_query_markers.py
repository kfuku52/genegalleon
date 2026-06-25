from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "annotate_stat_branch_query_markers.py"


def load_module():
    spec = spec_from_file_location("annotate_stat_branch_query_markers", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_stat_branch(path: Path):
    pandas.DataFrame(
        [
            {"branch_id": "0", "node_name": "Root", "so_event": "S"},
            {"branch_id": "1", "node_name": "Sp_geneA", "so_event": "L"},
            {"branch_id": "2", "node_name": "Sp_geneB", "so_event": "L"},
            {"branch_id": "3", "node_name": "Sp_geneC", "so_event": "L"},
        ]
    ).to_csv(path, sep="\t", index=False)


def test_gene_list_query_marks_direct_query_and_best_blast_hit(tmp_path):
    mod = load_module()
    stat_branch = tmp_path / "stat.branch.tsv"
    query_gene = tmp_path / "query_gene.txt"
    query_aa = tmp_path / "query.aa.fa"
    query_blast = tmp_path / "query_blast.tsv"
    outfile = tmp_path / "annotated.tsv"

    write_stat_branch(stat_branch)
    query_gene.write_text("geneA\ngeneB\n", encoding="utf-8")
    query_aa.write_text(">Sp_geneA\nMAAA\n", encoding="utf-8")
    pandas.DataFrame(
        [
            {"qacc": "Sp_geneA", "sacc": "Sp_geneA", "qjointcov": 1.0, "min_evalue": 0.0, "max_bitscore": 500},
            {"qacc": "Sp_geneA", "sacc": "Sp_geneB", "qjointcov": 0.8, "min_evalue": 1e-20, "max_bitscore": 300},
        ]
    ).to_csv(query_blast, sep="\t", index=False)

    mod.annotate_stat_branch(
        stat_branch_path=stat_branch,
        query_gene_path=query_gene,
        query_aa_fasta_path=query_aa,
        query_blast_path=query_blast,
        outfile=outfile,
        min_query_blast_coverage=0.25,
    )

    out = pandas.read_csv(outfile, sep="\t", dtype=str, keep_default_na=False).set_index("node_name")
    assert out.loc["Sp_geneA", "query_marker"] == "Best hit"
    assert out.loc["Sp_geneA", "query_marker_source"] == "direct:Sp_geneA;geneA|best:Sp_geneA"
    assert out.loc["Sp_geneB", "query_marker"] == "-"
    assert out.loc["Sp_geneB", "query_marker_source"] == "direct:geneB"
    assert out.loc["Root", "query_marker"] == "-"


def test_fasta_query_marks_best_retained_tip_per_query(tmp_path):
    mod = load_module()
    stat_branch = tmp_path / "stat.branch.tsv"
    query_gene = tmp_path / "query.fa"
    query_blast = tmp_path / "query_blast.tsv"
    outfile = tmp_path / "annotated.tsv"

    write_stat_branch(stat_branch)
    query_gene.write_text(">external_query\nMAAA\n", encoding="utf-8")
    pandas.DataFrame(
        [
            {"qacc": "external_query", "sacc": "Sp_geneB", "qjointcov": 0.7, "min_evalue": 1e-30, "max_bitscore": 250},
            {"qacc": "external_query", "sacc": "Sp_geneC", "qjointcov": 0.9, "min_evalue": 1e-5, "max_bitscore": 200},
            {"qacc": "external_query", "sacc": "not_in_tree", "qjointcov": 1.0, "min_evalue": 0.0, "max_bitscore": 999},
        ]
    ).to_csv(query_blast, sep="\t", index=False)

    mod.annotate_stat_branch(
        stat_branch_path=stat_branch,
        query_gene_path=query_gene,
        query_aa_fasta_path="",
        query_blast_path=query_blast,
        outfile=outfile,
        min_query_blast_coverage=0.25,
    )

    out = pandas.read_csv(outfile, sep="\t", dtype=str, keep_default_na=False).set_index("node_name")
    assert out.loc["Sp_geneC", "query_marker"] == "Best hit"
    assert out.loc["Sp_geneC", "query_marker_source"] == "best:external_query"
    assert out.loc["Sp_geneC", "query_marker_best_qjointcov"] == "0.9"
    assert out.loc["Sp_geneB", "query_marker"] == "-"
    assert out.loc["Sp_geneA", "query_marker"] == "-"


def test_existing_query_marker_columns_are_replaced(tmp_path):
    mod = load_module()
    stat_branch = tmp_path / "stat.branch.tsv"
    query_gene = tmp_path / "query.fa"
    query_blast = tmp_path / "query_blast.tsv"
    outfile = tmp_path / "annotated.tsv"

    write_stat_branch(stat_branch)
    df = pandas.read_csv(stat_branch, sep="\t", dtype=str)
    df["query_marker"] = "stale"
    df.to_csv(stat_branch, sep="\t", index=False)
    query_gene.write_text(">external_query\nMAAA\n", encoding="utf-8")
    pandas.DataFrame(
        [{"qacc": "external_query", "sacc": "Sp_geneB", "qjointcov": 0.9, "min_evalue": 1e-10, "max_bitscore": 100}]
    ).to_csv(query_blast, sep="\t", index=False)

    mod.annotate_stat_branch(
        stat_branch_path=stat_branch,
        query_gene_path=query_gene,
        query_aa_fasta_path="",
        query_blast_path=query_blast,
        outfile=outfile,
        min_query_blast_coverage=0.25,
    )

    out = pandas.read_csv(outfile, sep="\t", dtype=str, keep_default_na=False)
    assert list(out.columns).count("query_marker") == 1
    assert "stale" not in set(out["query_marker"])

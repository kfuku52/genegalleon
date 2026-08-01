from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace

import pandas

from workflow.support.gene_family_output_store import archive_completed_outputs, family_context

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "orthogroup_output_summary.py"


def load_module():
    spec = spec_from_file_location("orthogroup_output_summary", SCRIPT_PATH)
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


def test_extract_orthogroup_id_accepts_prefixed_filenames():
    mod = load_module()
    assert mod._extract_orthogroup_id("HOG0000010_cds.fa.gz") == "HOG0000010"
    assert mod._extract_orthogroup_id("OG0001234_amas.original.tsv") == "OG0001234"
    assert mod._extract_orthogroup_id("SP0000001.treefile") == "SP0000001"
    assert mod._extract_orthogroup_id("README.txt") is None


def test_extract_orthogroup_id_requires_boundary_after_numeric_id():
    mod = load_module()
    assert mod._extract_orthogroup_id("OG1abc.tsv") is None
    assert mod._extract_orthogroup_id("HOG0000010extra.tsv") is None
    assert mod._extract_orthogroup_id("SP1-extra.tsv") == "SP1"


def test_get_amas_stats_uses_index_ids_without_creating_extra_rows(tmp_path):
    mod = load_module()
    df = pandas.DataFrame({"Total": [1, 2]}, index=["HOG0000010", "HOG0000011"])
    amas_dir = tmp_path / "amas_original"
    amas_dir.mkdir()
    _write_amas(amas_dir / "HOG0000010_amas.original.tsv", 7)
    _write_amas(amas_dir / "HOG9999999_amas.original.tsv", 9)

    out = mod.get_amas_stats(df.copy(), str(amas_dir), "original", ncpu=1)

    assert set(out.index.tolist()) == {"HOG0000010", "HOG0000011"}
    assert out.loc["HOG0000010", "No_of_taxa_original"] == 7
    assert pandas.isna(out.loc["HOG0000011", "No_of_taxa_original"])


def test_run_handles_subdir_files_with_suffixes(tmp_path):
    mod = load_module()
    dir_og = tmp_path / "orthogroup"
    dir_og.mkdir()
    (dir_og / "cds_fasta").mkdir()
    (dir_og / "cds_fasta" / "HOG0000010_cds.fa.gz").write_text("", encoding="utf-8")
    (dir_og / "cds_fasta" / "HOG9999999_cds.fa.gz").write_text("", encoding="utf-8")

    genecount = tmp_path / "Orthogroups.GeneCount.selected.tsv"
    pandas.DataFrame(
        [
            {"Orthogroup": "HOG0000010", "Total": 1},
            {"Orthogroup": "HOG0000011", "Total": 1},
        ]
    ).to_csv(genecount, sep="\t", index=False)

    out_tsv = tmp_path / "orthogroup_summary.tsv"
    mod.run(
        SimpleNamespace(
            dir_og=str(dir_og),
            genecount=str(genecount),
            out=str(out_tsv),
            ncpu=1,
        )
    )

    out = pandas.read_csv(out_tsv, sep="\t", index_col=0)
    assert "GG_ARRAY_TASK_ID" in out.columns
    assert out.loc["HOG0000010", "cds_fasta"] == 1
    assert out.loc["HOG0000011", "cds_fasta"] == 0


def test_run_realigns_existing_amas_table_to_original_index(tmp_path):
    mod = load_module()
    dir_og = tmp_path / "orthogroup"
    dir_og.mkdir()

    genecount = tmp_path / "Orthogroups.GeneCount.selected.tsv"
    pandas.DataFrame(
        [
            {"Orthogroup": "HOG0000010", "Total": 1},
            {"Orthogroup": "HOG0000011", "Total": 1},
        ]
    ).to_csv(genecount, sep="\t", index=False)

    existing = tmp_path / "Orthogroups.GeneCount.selected.amas.tsv"
    pandas.DataFrame(
        [
            {"Orthogroup": "HOG0000010", "Total": 1, "No_of_taxa_original": 7},
            {"Orthogroup": "HOG0000011", "Total": 1, "No_of_taxa_original": 8},
            {"Orthogroup": "HOG9999999_cds", "Total": 1, "No_of_taxa_original": 9},
        ]
    ).to_csv(existing, sep="\t", index=False)

    out_tsv = tmp_path / "orthogroup_summary.tsv"
    mod.run(
        SimpleNamespace(
            dir_og=str(dir_og),
            genecount=str(genecount),
            out=str(out_tsv),
            ncpu=1,
        )
    )

    out = pandas.read_csv(out_tsv, sep="\t", index_col=0)
    assert set(out.index.tolist()) == {"HOG0000010", "HOG0000011"}


def test_run_ignores_hidden_dirs_and_hidden_files(tmp_path):
    mod = load_module()
    dir_og = tmp_path / "orthogroup"
    dir_og.mkdir()
    (dir_og / "cds_fasta").mkdir()
    (dir_og / ".hidden_stage").mkdir()
    (dir_og / "cds_fasta" / "HOG0000010_cds.fa.gz").write_text("", encoding="utf-8")
    (dir_og / "cds_fasta" / ".DS_Store").write_text("", encoding="utf-8")
    (dir_og / ".hidden_stage" / "HOG0000010.tmp").write_text("", encoding="utf-8")

    genecount = tmp_path / "Orthogroups.GeneCount.selected.tsv"
    pandas.DataFrame(
        [
            {"Orthogroup": "HOG0000010", "Total": 1},
            {"Orthogroup": "HOG0000011", "Total": 1},
        ]
    ).to_csv(genecount, sep="\t", index=False)

    out_tsv = tmp_path / "orthogroup_summary.tsv"
    mod.run(
        SimpleNamespace(
            dir_og=str(dir_og),
            genecount=str(genecount),
            out=str(out_tsv),
            ncpu=1,
        )
    )

    out = pandas.read_csv(out_tsv, sep="\t", index_col=0)
    assert "cds_fasta" in out.columns
    assert ".hidden_stage" not in out.columns
    assert out.loc["HOG0000010", "cds_fasta"] == 1


def test_get_amas_stats_ignores_non_file_entries_in_amas_dir(tmp_path):
    mod = load_module()
    df = pandas.DataFrame({"Total": [1]}, index=["HOG0000010"])
    amas_dir = tmp_path / "amas_original"
    amas_dir.mkdir()
    _write_amas(amas_dir / "HOG0000010_amas.original.tsv", 7)
    (amas_dir / "HOG0000011_placeholder").mkdir()

    out = mod.get_amas_stats(df.copy(), str(amas_dir), "original", ncpu=1)
    assert out.loc["HOG0000010", "No_of_taxa_original"] == 7


def test_run_reads_orthogroup_artifacts_from_zip_shards(tmp_path):
    mod = load_module()
    dir_og = tmp_path / "orthogroup"
    family_id = "HOG0000010"
    genecount = tmp_path / "Orthogroups.GeneCount.selected.tsv"
    pandas.DataFrame(
        [{"Orthogroup": family_id, "Total": 1}]
    ).to_csv(genecount, sep="\t", index=False)
    for subdir, suffix in (
        ("stat_branch", "_stat.branch.tsv"),
        ("stat_tree", "_stat.tree.tsv"),
        ("tree_plot", "_tree_plot.pdf"),
        ("cds_fasta", "_cds.fa.gz"),
    ):
        path = dir_og / subdir / f"{family_id}{suffix}"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"{subdir}\n", encoding="utf-8")
    amas_path = dir_og / "amas_original" / f"{family_id}_amas.original.tsv"
    amas_path.parent.mkdir(parents=True)
    _write_amas(amas_path, 13)
    family_ids, family_from_name = family_context("orthogroup", genecount=genecount)
    archive_completed_outputs(dir_og, "orthogroup", family_ids, family_from_name)

    out_tsv = tmp_path / "orthogroup_summary.tsv"
    mod.run(
        SimpleNamespace(
            dir_og=str(dir_og),
            genecount=str(genecount),
            out=str(out_tsv),
            ncpu=1,
        )
    )

    out = pandas.read_csv(out_tsv, sep="\t", index_col=0)
    assert out.loc[family_id, "cds_fasta"] == 1
    assert out.loc[family_id, "tree_plot"] == 1
    assert out.loc[family_id, "No_of_taxa_original"] == 13

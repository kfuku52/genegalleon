from pathlib import Path
import subprocess
import sys

import pandas


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "collect_common_BUSCO_genes.py"


def write_busco_full_table(path, rows):
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# BUSCO full table mock\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")


def test_collect_common_busco_filters_all_missing_and_merges_metadata_by_busco_id(tmp_path):
    busco_dir = tmp_path / "busco_full"
    busco_dir.mkdir()
    outfile = tmp_path / "busco_summary.tsv"

    write_busco_full_table(
        busco_dir / "speciesA.tsv",
        [
            ("b1", "Complete", "a_gene1:1-10", "100", "10", "-", "-"),
            ("b2", "Missing", "-", "-", "-", "-", "-"),
            ("b4", "Complete", "a_gene4_1:1-9", "90", "9", "url4", "desc4"),
            ("b4", "Complete", "a_gene4_2:1-9", "90", "9", "url4", "desc4"),
        ],
    )
    write_busco_full_table(
        busco_dir / "speciesB.tsv",
        [
            ("b1", "Complete", "b_gene1:5-20", "100", "16", "url1", "desc1"),
            ("b2", "Missing", "-", "-", "-", "-", "-"),
            ("b3", "Missing", "-", "-", "-", "-", "-"),
            ("b4", "Missing", "-", "-", "-", "-", "-"),
        ],
    )

    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--busco_outdir",
            str(busco_dir),
            "--outfile",
            str(outfile),
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr

    out = pandas.read_csv(outfile, sep="\t")
    assert out.columns.tolist() == ["busco_id", "orthodb_url", "description", "speciesA", "speciesB"]
    assert out["busco_id"].tolist() == ["b1", "b4"]

    row_b1 = out.loc[out["busco_id"] == "b1", :].iloc[0]
    assert row_b1["orthodb_url"] == "url1"
    assert row_b1["description"] == "desc1"
    assert row_b1["speciesA"] == "a_gene1"
    assert row_b1["speciesB"] == "b_gene1"

    row_b4 = out.loc[out["busco_id"] == "b4", :].iloc[0]
    assert row_b4["speciesA"] == "a_gene4_1,a_gene4_2"
    assert row_b4["speciesB"] == "-"


def test_collect_common_busco_ignores_hidden_and_non_file_tsv_entries(tmp_path):
    busco_dir = tmp_path / "busco_full"
    busco_dir.mkdir()
    outfile = tmp_path / "busco_summary.tsv"

    write_busco_full_table(
        busco_dir / "speciesA.tsv",
        [
            ("b1", "Complete", "a_gene1:1-10", "100", "10", "url1", "desc1"),
        ],
    )
    write_busco_full_table(
        busco_dir / ".hidden.tsv",
        [
            ("b1", "Complete", "hidden_gene:1-10", "100", "10", "urlH", "descH"),
        ],
    )
    (busco_dir / "fake.tsv").mkdir()

    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--busco_outdir",
            str(busco_dir),
            "--outfile",
            str(outfile),
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr

    out = pandas.read_csv(outfile, sep="\t")
    assert out.columns.tolist() == ["busco_id", "orthodb_url", "description", "speciesA"]
    assert out.loc[0, "speciesA"] == "a_gene1"


def test_collect_common_busco_returns_header_only_when_busco_dir_is_missing(tmp_path):
    outfile = tmp_path / "busco_summary.tsv"
    missing_dir = tmp_path / "missing_busco_full"

    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--busco_outdir",
            str(missing_dir),
            "--outfile",
            str(outfile),
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr

    out = pandas.read_csv(outfile, sep="\t")
    assert out.columns.tolist() == ["busco_id", "orthodb_url", "description"]
    assert out.empty


def test_collect_common_busco_accepts_short_missing_rows_with_parallel_workers(tmp_path):
    busco_dir = tmp_path / "busco_full"
    busco_dir.mkdir()
    outfile = tmp_path / "busco_summary.tsv"

    write_busco_full_table(
        busco_dir / "speciesA.tsv",
        [
            ("b1", "Complete", "a_gene1:1-10", "100", "10", "url1", "desc1"),
            ("b2", "Missing"),
            ("b3", "Missing"),
        ],
    )
    write_busco_full_table(
        busco_dir / "speciesB.tsv",
        [
            ("b1", "Missing"),
            ("b2", "Missing"),
            ("b3", "Complete", "b_gene3:1-7", "90", "7", "url3", "desc3"),
        ],
    )

    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--busco_outdir",
            str(busco_dir),
            "--outfile",
            str(outfile),
            "--ncpu",
            "2",
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr

    out = pandas.read_csv(outfile, sep="\t")
    assert out.columns.tolist() == ["busco_id", "orthodb_url", "description", "speciesA", "speciesB"]
    assert out["busco_id"].tolist() == ["b1", "b3"]

    row_b1 = out.loc[out["busco_id"] == "b1", :].iloc[0]
    assert row_b1["orthodb_url"] == "url1"
    assert row_b1["description"] == "desc1"
    assert row_b1["speciesA"] == "a_gene1"
    assert row_b1["speciesB"] == "-"

    row_b3 = out.loc[out["busco_id"] == "b3", :].iloc[0]
    assert row_b3["orthodb_url"] == "url3"
    assert row_b3["description"] == "desc3"
    assert row_b3["speciesA"] == "-"
    assert row_b3["speciesB"] == "b_gene3"

from pathlib import Path
import subprocess
import sys

import pandas


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "script" / "collect_common_BUSCO_genes.py"


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
    assert out.columns.tolist() == ["busco_id", "orthodb_url", "description", "speciesA.tsv", "speciesB.tsv"]
    assert out["busco_id"].tolist() == ["b1", "b4"]

    row_b1 = out.loc[out["busco_id"] == "b1", :].iloc[0]
    assert row_b1["orthodb_url"] == "url1"
    assert row_b1["description"] == "desc1"
    assert row_b1["speciesA.tsv"] == "a_gene1"
    assert row_b1["speciesB.tsv"] == "b_gene1"

    row_b4 = out.loc[out["busco_id"] == "b4", :].iloc[0]
    assert row_b4["speciesA.tsv"] == "a_gene4_1,a_gene4_2"
    assert row_b4["speciesB.tsv"] == "-"

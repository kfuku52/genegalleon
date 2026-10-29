import csv
import subprocess
import sys
from pathlib import Path

import pandas
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "collect_common_BUSCO_genes.py"


def test_collector_preserves_colon_numeric_gene_identifiers(tmp_path):
    source = tmp_path / "full"
    source.mkdir()
    write_busco_full_table(source / "A.tsv", [
        ("b1", "Complete", "gene:001", "100", "3", "-", "-"),
        ("b2", "Complete", "gene:1.2", "100", "3", "-", "-"),
        ("b3", "Complete", "gene:1-9", "100", "3", "-", "-"),
        ("b4", "Complete", "gene:1..9", "100", "3", "-", "-"),
        ("b5", "Complete", "gene:001:1-9", "100", "3", "-", "-"),
        ("b6", "Complete", "gene:", "100", "3", "-", "-"),
    ])
    output = tmp_path / "summary.tsv"
    subprocess.run([sys.executable, str(SCRIPT_PATH), "--busco_outdir", str(source), "--outfile", str(output)], check=True)
    with output.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["A"] for row in rows] == ["gene:001", "gene:1.2", "gene", "gene", "gene:001", "gene:"]


def test_collector_rejects_ambiguous_species_tables_before_replacing_output(tmp_path):
    source = tmp_path / "full"
    source.mkdir()
    for name in ("A.tsv", "A.busco.full.tsv"):
        write_busco_full_table(source / name, [("b1", "Complete", "g1", "100", "3", "-", "-")])
    output = tmp_path / "summary.tsv"
    output.write_text("previous valid result\n")
    result = subprocess.run([sys.executable, str(SCRIPT_PATH), "--busco_outdir", str(source), "--outfile", str(output)],
                            text=True, capture_output=True)
    assert result.returncode != 0
    assert "Ambiguous BUSCO full tables" in result.stderr
    assert output.read_text() == "previous valid result\n"


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


def run_collector(busco_dir, outfile, workers):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "--busco_outdir", str(busco_dir),
         "--outfile", str(outfile), "--ncpu", str(workers)],
        capture_output=True, text=True, check=False,
    )


def test_collect_common_busco_metadata_fields_have_independent_first_nonmissing_values(tmp_path):
    busco_dir = tmp_path / "tables"
    busco_dir.mkdir()
    descriptions = ["説明 α\tquoted field\nsecond line", "nan", "None"]
    rows = {
        "speciesA": [
            ("b1", "Complete", "a1:1-9", "100", "9", "-", descriptions[0]),
            ("b1", "Duplicated", "a1_copy:1-9", "100", "9", "ignored-url", "ignored-description"),
            ("b2", "Complete", "a2:1-9", "100", "9", "-", "-"),
            ("b2", "Duplicated", "a2_copy:1-9", "100", "9", "ignored-url", "ignored-description"),
            ("001", "Complete", "a3:1-9", "100", "9", "NA", descriptions[1]),
            ("b4", "Complete", "a4:1-9", "100", "9", "null", descriptions[2]),
        ],
        "speciesB": [
            ("b1", "Complete", "b1:1-9", "100", "9", "url-from-B", "-"),
            ("b2", "Complete", "b2:1-9", "100", "9", "-", "-"),
            ("001", "Complete", "b3:1-9", "100", "9", "later-url", "later-description"),
            ("b4", "Missing"),
        ],
    }
    for species, records in rows.items():
        with (busco_dir / f"{species}.tsv").open("w", encoding="utf-8", newline="") as handle:
            csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(records)
    outputs = []
    for workers in (1, 2):
        output = tmp_path / f"summary-{workers}.tsv"
        completed = run_collector(busco_dir, output, workers)
        assert completed.returncode == 0, completed.stderr
        assert completed.stderr == ""
        assert set(completed.stdout.splitlines()) == {
            "Starting collect_common_BUSCO_gene.py", "Working on speciesA.tsv", "Working on speciesB.tsv",
            "Ending collect_common_BUSCO_gene.py",
        }
        outputs.append(output.read_bytes())
    assert outputs[0] == outputs[1]
    out = pandas.read_csv(tmp_path / "summary-1.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert out.columns.tolist() == ["busco_id", "orthodb_url", "description", "speciesA", "speciesB"]
    assert out.to_dict("records") == [
        {"busco_id": "001", "orthodb_url": "NA", "description": "nan", "speciesA": "a3", "speciesB": "b3"},
        {"busco_id": "b1", "orthodb_url": "url-from-B", "description": descriptions[0],
         "speciesA": "a1,a1_copy", "speciesB": "b1"},
        {"busco_id": "b2", "orthodb_url": "-", "description": "-", "speciesA": "a2,a2_copy", "speciesB": "b2"},
        {"busco_id": "b4", "orthodb_url": "null", "description": "None", "speciesA": "a4", "speciesB": "-"},
    ]


@pytest.mark.parametrize("workers", [1, 2])
@pytest.mark.parametrize("with_gene", [False, True])
def test_collect_common_busco_preserves_empty_species_columns(tmp_path, workers, with_gene):
    busco_dir = tmp_path / "tables"
    busco_dir.mkdir()
    (busco_dir / "speciesA.tsv").write_text("# only a comment\n", encoding="utf-8")
    (busco_dir / "speciesB.tsv").write_text(
        "b1\tComplete\tb_gene:1-9\t100\t9\t-\t-\n" if with_gene else "", encoding="utf-8",
    )
    output = tmp_path / "summary.tsv"
    completed = run_collector(busco_dir, output, workers)
    assert completed.returncode == 0, completed.stderr
    assert completed.stderr == ""
    expected = "busco_id\torthodb_url\tdescription\tspeciesA\tspeciesB\n"
    if with_gene:
        expected += "b1\t-\t-\t-\tb_gene\n"
    assert output.read_text() == expected


def test_collect_common_busco_rejects_invalid_utf8_before_replacing_output(tmp_path):
    busco_dir = tmp_path / "tables"
    busco_dir.mkdir()
    (busco_dir / "speciesA.tsv").write_bytes(b"b1\tComplete\ta1\t100\t9\turl\tdescription\n\xff\n")
    output = tmp_path / "summary.tsv"
    output.write_text("existing output\n")
    completed = run_collector(busco_dir, output, 1)
    assert completed.returncode != 0
    assert "UnicodeDecodeError" in completed.stderr
    assert output.read_text() == "existing output\n"

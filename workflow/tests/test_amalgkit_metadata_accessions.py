import csv
import sys
from pathlib import Path

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
if str(SUPPORT_DIR) not in sys.path:
    sys.path.insert(0, str(SUPPORT_DIR))

import amalgkit_metadata_accessions as ama  # noqa: E402


def write_tsv(path, fieldnames, rows):
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_missing_accessions_preserves_requested_order(tmp_path, capsys):
    metadata = tmp_path / "metadata.tsv"
    accessions = tmp_path / "accessions.txt"
    write_tsv(metadata, ["run", "lib_source"], [{"run": "SRR1", "lib_source": "TRANSCRIPTOMIC"}])
    accessions.write_text("SRR1\nSRR2\nSRR2\nSRR3\n", encoding="utf-8")

    ama.print_missing_accessions(metadata, accessions)

    assert capsys.readouterr().out.splitlines() == ["SRR2", "SRR3"]


def test_extract_transcriptomic_rows_filters_requested_runs(tmp_path):
    metadata = tmp_path / "metadata.tsv"
    accessions = tmp_path / "accessions.txt"
    output = tmp_path / "relaxed.tsv"
    write_tsv(
        metadata,
        ["run", "lib_source", "lib_strategy"],
        [
            {"run": "SRR1", "lib_source": "TRANSCRIPTOMIC", "lib_strategy": "RNA-Seq"},
            {"run": "SRR2", "lib_source": "GENOMIC", "lib_strategy": "WGS"},
            {"run": "SRR3", "lib_source": "", "lib_strategy": "CLONE"},
        ],
    )
    accessions.write_text("SRR1\nSRR2\nSRR3\n", encoding="utf-8")

    ama.extract_transcriptomic_rows(metadata, accessions, output)

    assert [row["run"] for row in read_tsv(output)] == ["SRR1", "SRR3"]


def test_merge_metadata_tables_by_run_keeps_primary_rows_first(tmp_path):
    primary = tmp_path / "primary.tsv"
    extra = tmp_path / "extra.tsv"
    output = tmp_path / "merged.tsv"
    write_tsv(primary, ["run", "species"], [{"run": "SRR1", "species": "A"}])
    write_tsv(
        extra,
        ["run", "species", "extra"],
        [
            {"run": "SRR1", "species": "A2", "extra": "x"},
            {"run": "SRR2", "species": "B", "extra": "y"},
        ],
    )

    ama.merge_metadata_tables_by_run(primary, extra, output)

    rows = read_tsv(output)
    assert [row["run"] for row in rows] == ["SRR1", "SRR2"]
    assert rows[0]["species"] == "A"
    assert rows[1]["extra"] == "y"

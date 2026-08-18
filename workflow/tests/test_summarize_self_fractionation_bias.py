import csv
import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "workflow" / "support" / "summarize_self_fractionation_bias.py"


def test_summarizer_separates_interchromosomal_best_retention(tmp_path: Path) -> None:
    table = tmp_path / "fractionation_bias_pairs.tsv"
    table.write_text(
        "analysis_id\ttarget_species\tquery_species\tquota\tmode\n"
        "pair\tTarget\tQuery\t1:2\tcompare\n"
        "self1\tSpecies_one\tSpecies_one\t1:1\tself\n",
        encoding="utf-8",
    )
    result_dir = tmp_path / "results" / "self1"
    result_dir.mkdir(parents=True)
    (result_dir / "self1.summary.json").write_text(
        json.dumps(
            {
                "analysis_mode": "self_synteny_retention",
                "parameters": {"window_size": 2, "denominator": "all"},
                "counts": {
                    "target_gene_count": 10,
                    "synteny_pair_count": 4,
                    "interchromosomal_pair_count": 3,
                    "intrachromosomal_pair_count": 1,
                },
            }
        ),
        encoding="utf-8",
    )
    (result_dir / "self1.windows.tsv").write_text(
        "target_seqid\tquery_seqid\twindow_index\tstart_rank\tend_rank\tretention_percent\n"
        "chr1\tchr1\t1\t1\t2\t100\n"
        "chr1\tchr2\t1\t1\t2\t50\n"
        "chr1\tchr3\t1\t1\t2\t75\n"
        "chr1\tchr2\t2\t2\t3\t25\n"
        "chr1\tchr3\t2\t2\t3\t0\n",
        encoding="utf-8",
    )
    output_tsv = tmp_path / "summary.tsv"
    output_pdf = tmp_path / "summary.pdf"
    output_png = tmp_path / "summary.png"
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--table",
            str(table),
            "--results-dir",
            str(tmp_path / "results"),
            "--output-tsv",
            str(output_tsv),
            "--output-pdf",
            str(output_pdf),
            "--output-png",
            str(output_png),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    with output_tsv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["analysis_id"] == "self1"
    assert rows[0]["analyzed_window_count"] == "2"
    assert rows[0]["mean_best_interchromosomal_retention_percent"] == "50"
    assert rows[0]["median_best_interchromosomal_retention_percent"] == "50"
    assert rows[0]["maximum_best_interchromosomal_retention_percent"] == "75"
    assert rows[0]["mean_intrachromosomal_retention_percent"] == "100"
    assert output_pdf.stat().st_size > 0
    assert output_png.stat().st_size > 0

from pathlib import Path
import subprocess
import sys

import pandas


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "extract_gene_id_from_blast_table.py"


def run_extract(tmp_path, rows, min_cov=0.25, max_hits=5000):
    infile = tmp_path / "query_blast.tsv"
    outfile = tmp_path / "gene_ids.txt"
    pandas.DataFrame(rows).to_csv(infile, sep="\t", index=False)

    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--infile",
            str(infile),
            "--outfile",
            str(outfile),
            "--min_query_blast_coverage",
            str(min_cov),
            "--max_num_gene_blast_hit_retrieval",
            str(max_hits),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    return completed, outfile


def test_extract_parses_semicolon_evalue_and_sorts_by_best_evalue(tmp_path):
    rows = [
        {"sacc": "gene_A", "qjointcov": 0.5, "evalue": "1e-2;1e-20"},
        {"sacc": "gene_B", "qjointcov": 0.5, "evalue": "1e-5"},
        {"sacc": "gene_C", "qjointcov": 0.2, "evalue": "1e-30"},
    ]
    completed, outfile = run_extract(tmp_path, rows, min_cov=0.25, max_hits=5)
    assert completed.returncode == 0, completed.stderr

    genes = [line.strip() for line in outfile.read_text().splitlines() if line.strip()]
    assert genes == ["gene_A", "gene_B"]


def test_extract_uses_min_evalue_column_when_present(tmp_path):
    rows = [
        {"sacc": "gene_X", "qjointcov": 0.5, "evalue": "1e-2;1e-20", "min_evalue": 1e-2},
        {"sacc": "gene_Y", "qjointcov": 0.5, "evalue": "1e-5;1e-7", "min_evalue": 1e-7},
    ]
    completed, outfile = run_extract(tmp_path, rows, min_cov=0.25, max_hits=5)
    assert completed.returncode == 0, completed.stderr

    genes = [line.strip() for line in outfile.read_text().splitlines() if line.strip()]
    assert genes == ["gene_Y", "gene_X"]


def test_extract_keeps_gene_when_best_evalue_is_nan(tmp_path):
    rows = [
        {"sacc": "gene_A", "qjointcov": 0.9, "evalue": ""},
        {"sacc": "gene_A", "qjointcov": 0.7, "evalue": "1e-20"},
        {"sacc": "gene_B", "qjointcov": 0.85, "evalue": "1e-30"},
    ]
    completed, outfile = run_extract(tmp_path, rows, min_cov=0.25, max_hits=5)
    assert completed.returncode == 0, completed.stderr

    genes = [line.strip() for line in outfile.read_text().splitlines() if line.strip()]
    assert genes == ["gene_A", "gene_B"]


def test_extract_writes_empty_file_when_all_hits_fail_coverage_threshold(tmp_path):
    rows = [
        {"sacc": "gene_A", "qjointcov": 0.1, "evalue": "1e-10"},
        {"sacc": "gene_B", "qjointcov": 0.2, "evalue": "1e-20"},
    ]
    completed, outfile = run_extract(tmp_path, rows, min_cov=0.25, max_hits=5)
    assert completed.returncode == 0, completed.stderr
    assert outfile.read_text() == ""

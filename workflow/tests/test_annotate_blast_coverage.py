from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import subprocess
import sys

import numpy
import pandas


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "script" / "annotate_blast_coverage.py"


def load_module():
    spec = spec_from_file_location("annotate_blast_coverage", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_multi_query_sparse_pairs_do_not_crash_and_compute_qjointcov():
    mod = load_module()
    df = pandas.DataFrame(
        [
            {
                "qacc": "q1",
                "sacc": "s1",
                "pident": "99.0",
                "length": 10,
                "mismatch": 0,
                "gapopen": 0,
                "qstart": 1,
                "qend": 10,
                "sstart": 101,
                "send": 110,
                "evalue": "1e-20",
                "bitscore": "300",
                "frames": "0/1",
                "qlen": 100,
                "slen": 500,
            },
            {
                "qacc": "q1",
                "sacc": "s1",
                "pident": "95.0",
                "length": 11,
                "mismatch": 1,
                "gapopen": 2,
                "qstart": 20,
                "qend": 30,
                "sstart": 120,
                "send": 130,
                "evalue": "3e-10",
                "bitscore": "200",
                "frames": "0/1",
                "qlen": 100,
                "slen": 500,
            },
            {
                "qacc": "q2",
                "sacc": "s2",
                "pident": "90.0",
                "length": 11,
                "mismatch": 1,
                "gapopen": 0,
                "qstart": 5,
                "qend": 15,
                "sstart": 205,
                "send": 215,
                "evalue": "5e-8",
                "bitscore": "180",
                "frames": "0/1",
                "qlen": 100,
                "slen": 600,
            },
        ]
    )

    out = mod.annotate_blast_coverage(df)
    assert len(out) == 2

    q1_row = out.loc[(out["qacc"] == "q1") & (out["sacc"] == "s1"), :].iloc[0]
    q2_row = out.loc[(out["qacc"] == "q2") & (out["sacc"] == "s2"), :].iloc[0]
    q1_cov = q1_row["qjointcov"]
    q2_cov = q2_row["qjointcov"]
    assert q1_cov == 0.21
    assert q2_cov == 0.11
    assert q1_row["num_hits"] == 2
    assert q1_row["qstart"] == "1;20"
    assert q1_row["gapopen"] == "0;2"
    assert q1_row["evalue"] == "1e-20;3e-10"
    assert q1_row["qhitcov"] == "0.1;0.11"
    assert q1_row["qhitcov"] != str(q1_row["qjointcov"])


def test_reverse_coordinates_are_handled_for_per_hsp_and_joint_coverage():
    mod = load_module()
    df = pandas.DataFrame(
        [
            {
                "qacc": "q1",
                "sacc": "s1",
                "pident": "50.0",
                "length": 41,
                "mismatch": 20,
                "gapopen": 1,
                "qstart": 50,
                "qend": 10,
                "sstart": 900,
                "send": 700,
                "evalue": "1e-5",
                "bitscore": "120",
                "frames": "0/-1",
                "qlen": 100,
                "slen": 1000,
            }
        ]
    )

    out = mod.annotate_blast_coverage(df)
    assert out["qhitcov"].iloc[0] == "0.41"
    assert out["shitcov"].iloc[0] == "0.201"
    assert out["qjointcov"].iloc[0] == 0.41


def test_parallel_matches_serial_output():
    mod = load_module()
    df = pandas.DataFrame(
        [
            {
                "qacc": "q1",
                "sacc": "s1",
                "pident": "99.0",
                "length": 10,
                "mismatch": 0,
                "gapopen": 0,
                "qstart": 1,
                "qend": 10,
                "sstart": 101,
                "send": 110,
                "evalue": "1e-20",
                "bitscore": "300",
                "frames": "0/1",
                "qlen": 100,
                "slen": 500,
            },
            {
                "qacc": "q1",
                "sacc": "s1",
                "pident": "95.0",
                "length": 11,
                "mismatch": 1,
                "gapopen": 2,
                "qstart": 20,
                "qend": 30,
                "sstart": 120,
                "send": 130,
                "evalue": "3e-10",
                "bitscore": "200",
                "frames": "0/1",
                "qlen": 100,
                "slen": 500,
            },
            {
                "qacc": "q2",
                "sacc": "s2",
                "pident": "90.0",
                "length": 11,
                "mismatch": 1,
                "gapopen": 0,
                "qstart": 5,
                "qend": 15,
                "sstart": 205,
                "send": 215,
                "evalue": "5e-8",
                "bitscore": "180",
                "frames": "0/1",
                "qlen": 100,
                "slen": 600,
            },
        ]
    )
    serial = mod.annotate_blast_coverage(df, ncpu=1).reset_index(drop=True)
    parallel = mod.annotate_blast_coverage(df, ncpu=2).reset_index(drop=True)
    pandas.testing.assert_frame_equal(serial, parallel)


def test_nan_values_in_concatenated_columns_are_rendered_as_empty_tokens():
    mod = load_module()
    df = pandas.DataFrame(
        [
            {
                "qacc": "q1",
                "sacc": "s1",
                "pident": 90.0,
                "length": 10,
                "mismatch": 0,
                "gapopen": 0,
                "qstart": 1,
                "qend": 10,
                "sstart": 1,
                "send": 10,
                "evalue": "1e-20",
                "bitscore": 100,
                "frames": "0/1",
                "qlen": 100,
                "slen": 200,
            },
            {
                "qacc": "q1",
                "sacc": "s1",
                "pident": numpy.nan,
                "length": 11,
                "mismatch": 1,
                "gapopen": numpy.nan,
                "qstart": 20,
                "qend": 30,
                "sstart": 20,
                "send": 30,
                "evalue": "1e-10",
                "bitscore": 90,
                "frames": "0/1",
                "qlen": 100,
                "slen": 200,
            },
        ]
    )

    out = mod.annotate_blast_coverage(df, ncpu=1)
    row = out.iloc[0]
    assert row["pident"] == "90.0;"
    assert row["gapopen"] == "0.0;"


def test_cli_writes_expected_columns(tmp_path):
    infile = tmp_path / "in.tsv"
    outfile = tmp_path / "out.tsv"

    df = pandas.DataFrame(
        [
            {
                "qacc": "q1",
                "sacc": "s1",
                "pident": "80.0",
                "length": 10,
                "mismatch": 1,
                "gapopen": 0,
                "qstart": 1,
                "qend": 10,
                "sstart": 1,
                "send": 10,
                "evalue": "1e-6",
                "bitscore": "50",
                "frames": "0/1",
                "qlen": 100,
                "slen": 400,
            }
        ]
    )
    df.to_csv(infile, sep="\t", index=False)

    completed = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "--in", str(infile), "--out", str(outfile)],
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    out = pandas.read_csv(outfile, sep="\t")
    for col in ("qhitcov", "shitcov", "qjointcov", "num_hits", "min_evalue"):
        assert col in out.columns
    assert out["qjointcov"].iloc[0] == 0.1
    assert out["num_hits"].iloc[0] == 1


def test_cli_accepts_ncpu(tmp_path):
    infile = tmp_path / "in.tsv"
    outfile = tmp_path / "out.tsv"

    df = pandas.DataFrame(
        [
            {
                "qacc": "q1",
                "sacc": "s1",
                "pident": "80.0",
                "length": 10,
                "mismatch": 1,
                "gapopen": 0,
                "qstart": 1,
                "qend": 10,
                "sstart": 1,
                "send": 10,
                "evalue": "1e-6",
                "bitscore": "50",
                "frames": "0/1",
                "qlen": 100,
                "slen": 400,
            }
        ]
    )
    df.to_csv(infile, sep="\t", index=False)

    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--in",
            str(infile),
            "--out",
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
    assert len(out) == 1
    assert out["num_hits"].iloc[0] == 1


def test_empty_input_includes_joint_coverage_columns():
    mod = load_module()
    df = pandas.DataFrame(
        columns=["qacc", "sacc", "qstart", "qend", "sstart", "send", "qlen", "slen"]
    )
    out = mod.annotate_blast_coverage(df)
    assert "qjointcov" in out.columns
    assert "sjointcov" in out.columns


def test_cli_headerless_outfmt_and_frame_filter(tmp_path):
    infile = tmp_path / "in.tsv"
    outfile = tmp_path / "out.tsv"

    # Headerless rows (BLAST outfmt 6 style)
    infile.write_text(
        "\n".join(
            [
                "q1\ts1\t80\t10\t1\t0\t1\t10\t1\t10\t1e-6\t50\t0/1\t100\t400",
                "q1\ts1\t85\t10\t0\t0\t11\t20\t11\t20\t1e-7\t55\t0/-1\t100\t400",
            ]
        )
        + "\n"
    )

    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--in",
            str(infile),
            "--out",
            str(outfile),
            "--outfmt-columns",
            "qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore frames qlen slen",
            "--frame-filter",
            "0/1",
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    out = pandas.read_csv(outfile, sep="\t")
    assert len(out) == 1
    assert out["num_hits"].iloc[0] == 1
    assert out["frames"].iloc[0] == "0/1"


def test_cli_frame_filter_requires_frames_column(tmp_path):
    infile = tmp_path / "in.tsv"
    outfile = tmp_path / "out.tsv"

    infile.write_text("q1\ts1\t1\t10\t1\t10\t100\t400\n")

    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--in",
            str(infile),
            "--out",
            str(outfile),
            "--outfmt-columns",
            "qacc sacc qstart qend sstart send qlen slen",
            "--frame-filter",
            "0/1",
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode != 0
    assert "'frames' column was not found" in completed.stderr

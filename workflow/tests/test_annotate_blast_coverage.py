import gzip
import subprocess
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy
import pandas
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "annotate_blast_coverage.py"


def load_module():
    spec = spec_from_file_location("annotate_blast_coverage", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize('headerless', [False, True])
@pytest.mark.parametrize('compressed', [False, True])
def test_blast_coverage_preserves_literal_identifiers(tmp_path, headerless, compressed):
    mod = load_module()
    columns = ['qacc', 'sacc', 'qstart', 'qend', 'sstart', 'send', 'qlen', 'slen']
    genes = ['NA', 'NULL', 'nan', '001']
    subjects = genes[::-1]
    text = '' if headerless else '\t'.join(columns) + '\n'
    text += ''.join(f'{gene}\t{subject}\t1\t10\t1\t10\t100\t200\n'
                    for gene, subject in zip(genes, subjects, strict=True))
    text += 'NA\t001\t11\t20\t20\t11\t100\t200\n'
    path = tmp_path / ('input.tsv.gz' if compressed else 'input.tsv')
    if compressed:
        with gzip.open(path, 'wt') as handle:
            handle.write(text)
    else:
        path.write_text(text)
    frame = mod.read_blast_table(path, columns if headerless else None)
    assert frame.qacc.tolist() == genes + ['NA']
    assert frame.sacc.tolist() == subjects + ['001']
    assert str(frame.qstart.dtype) == 'int64'
    output = tmp_path / 'coverage.tsv'
    command = [sys.executable, str(SCRIPT_PATH), '--in', str(path), '--out', str(output)]
    if headerless:
        command += ['--outfmt-columns', ' '.join(columns)]
    completed = subprocess.run(command, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    actual = pandas.read_csv(output, sep='\t', converters={'qacc': str, 'sacc': str})
    assert actual.qacc.tolist() == genes
    assert actual.sacc.tolist() == subjects
    assert actual.num_hits.tolist() == [2, 1, 1, 1]
    assert actual.qjointcov.tolist() == [0.2, 0.1, 0.1, 0.1]
    assert actual.sjointcov.tolist() == [0.1, 0.05, 0.05, 0.05]


def test_blast_reader_retains_missing_empty_accessions_and_numeric_values(tmp_path):
    mod = load_module()
    path = tmp_path / 'missing.tsv'
    path.write_text('qacc\tsacc\tpident\n\tNA\tNULL\n001\t\tNaN\n')
    frame = mod.read_blast_table(path)
    assert pandas.isna(frame.qacc.iloc[0]) and frame.qacc.iloc[1] == '001'
    assert frame.sacc.iloc[0] == 'NA' and pandas.isna(frame.sacc.iloc[1])
    assert frame.pident.isna().all()


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
        [sys.executable, str(SCRIPT_PATH), "--in", str(infile), "--out", str(outfile), "--ncpu", "2"],
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


@pytest.mark.parametrize("optional_columns", [False, True])
def test_single_hsp_matches_identical_repeated_interval_contract(optional_columns):
    mod = load_module()
    lengths = [100.0, 0.0, -100.0, numpy.nan, numpy.inf, -numpy.inf, 3.0]
    df = pandas.DataFrame({
        "qacc": ["NA", "q2", "NULL", "q4", "nan", "001", "wide"],
        "sacc": ["s1", "s2", "s1", "s2", "s1", "s2", "s3"],
        "qstart": [10, 2, 12, 5, 7, 9, -(2**63)],
        "qend": [1, 2, 3, 15, 1, 20, 2**63 - 1],
        "sstart": [1, 12, 2, 15, 5, 20, -(2**63)],
        "send": [10, 2, 3, 5, 7, 9, 2**63 - 1],
        "qlen": lengths,
        "slen": lengths[::-1],
    })
    if optional_columns:
        df["pident"] = [90.0, numpy.nan, 0.0, 0.1, 99.0, 20.0, 100.0]
        df["length"] = [10, 1, 10, 11, 7, 12, 1]
        df["gapopen"] = [0, 1, 2, 3, 4, 5, 0]
        df["evalue"] = ["1e-20", "invalid", None, "", "0.0", "-0.0", "nan"]
        df["bitscore"] = [100.0, numpy.nan, -0.0, numpy.inf, -numpy.inf, 0.0, 20.5]
        df["frames"] = ["0/1", "0/-1", None, "0/1", "0/1", "0/-1", "0/1"]
    missing = df.iloc[[0, 1]].copy()
    missing.loc[missing.index[0], "qacc"] = None
    missing.loc[missing.index[1], "sacc"] = None
    df = pandas.concat([df, missing], ignore_index=True)

    # Repeating an HSP keeps its interval union and numeric extrema unchanged.
    # It also provides a reference using the general multi-HSP aggregation path.
    with numpy.errstate(divide="ignore", invalid="ignore", over="ignore"):
        single = mod.annotate_blast_coverage(df)
        repeated = mod.annotate_blast_coverage(pandas.concat([df, df], ignore_index=True))
    assert single["qacc"].tolist() == ["NA", "q2", "NULL", "q4", "nan", "001", "wide"]
    assert single.dtypes.equals(repeated.dtypes)
    repeated["num_hits"] = 1
    for column in (*mod.HIT_COLUMNS_TO_CONCAT, "qhitcov", "shitcov"):
        if column in repeated:
            dtype = repeated[column].dtype
            repeated[column] = repeated[column].str.split(";", n=1).str[0].astype(dtype)
    pandas.testing.assert_frame_equal(single, repeated, check_exact=True)
    assert single.to_csv(sep="\t", index=False) == repeated.to_csv(sep="\t", index=False)


def test_mixed_hsp_groups_keep_input_tokens_union_and_group_order():
    mod = load_module()
    df = pandas.DataFrame({
        "qacc": ["q2", "q1", "q2", "q3", "q2", "q1"],
        "sacc": ["s1", "s2", "s1", "s3", "s1", "s4"],
        "qstart": [1, 30, 10, 20, 50, 5],
        "qend": [10, 21, 20, 11, 41, 14],
        "sstart": [50, 1, 41, 40, 100, 10],
        "send": [41, 10, 31, 31, 91, 1],
        "qlen": [100, 100, 200, 100, 100, 100],
        "slen": [200, 200, 200, 200, 200, 200],
        "evalue": ["1e-5", None, "1e-10", "invalid", "1e-8", "1e-7"],
        "frames": ["0/1", "0/-1", "0/1", "0/-1", "0/1", "0/1"],
    })
    out = mod.annotate_blast_coverage(df)
    assert list(zip(out.qacc, out.sacc, strict=True)) == [("q2", "s1"), ("q1", "s2"), ("q1", "s4"), ("q3", "s3")]
    assert out.num_hits.tolist() == [3, 1, 1, 1]
    assert out.qstart.tolist() == ["1;10;50", "30", "5", "20"]
    assert out.evalue.tolist() == ["1e-5;1e-10;1e-8", "", "1e-7", "invalid"]
    assert out.frames.tolist() == ["0/1;0/1;0/1", "0/-1", "0/1", "0/-1"]
    assert out.qjointcov.tolist() == [0.15, 0.1, 0.1, 0.1]
    assert out.sjointcov.tolist() == [0.15, 0.05, 0.05, 0.05]
    assert out.qhitcov.tolist() == ["0.1;0.055;0.1", "0.1", "0.1", "0.1"]

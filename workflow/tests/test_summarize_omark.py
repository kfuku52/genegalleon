from pathlib import Path
import subprocess
import sys

import pandas


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "summarize_omark.py"


def run_script(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def test_summarize_omark_parses_sum_files(tmp_path):
    omark_outdir = tmp_path / "omark"
    species_dir = omark_outdir / "Arabidopsis_thaliana"
    species_dir.mkdir(parents=True)
    outfile = tmp_path / "omark_summary.tsv"

    (species_dir / "Arabidopsis_thaliana.sum").write_text(
        (
            "#The selected clade was Brassicales\n"
            "#Number of conserved HOGs is: 1234\n"
            "#Results on conserved HOGs is:\n"
            "#S:Single:S, D:Duplicated[U:Unexpected,E:Expected],M:Missing\n"
            "S:1200,D:20[U:5,E:15],M:14\n"
            "S:97.24%,D:1.62%[U:0.41%,E:1.21%],M:1.13%\n"
            "#On the whole proteome, there are 25000 proteins\n"
            "#Of which:\n"
            "#A:Consistent (taxonomically)[P:Partial hits,F:Fragmented], I: Inconsistent (taxonomically)[P:Partial hits,F:Fragmented], C: Likely Contamination[P:Partial hits,F:Fragmented], U: Unknown \n"
            "A:24000[P:100,F:50],I:500[P:30,F:20],C:50[P:10,F:5],U:450\n"
            "A:96.00%[P:0.40%,F:0.20%],I:2.00%[P:0.12%,F:0.08%],C:0.20%[P:0.04%,F:0.02%],U:1.80%\n"
            "#From HOG placement, the detected species are:\n"
            "#Clade\tNCBI taxid\tNumber of associated proteins\tPercentage of proteome's total\n"
            "Arabidopsis thaliana\t3702\t22000\t88.0%\n"
            "Brassicales\t3699\t500\t2.0%\n"
        ),
        encoding="utf-8",
    )

    completed = run_script("--omark_outdir", str(omark_outdir), "--outfile", str(outfile))

    assert completed.returncode == 0, completed.stderr
    df = pandas.read_csv(outfile, sep="\t")
    assert df.shape[0] == 1
    row = df.iloc[0]
    assert row["species"] == "Arabidopsis_thaliana"
    assert row["selected_clade"] == "Brassicales"
    assert row["conserved_hogs"] == 1234
    assert row["single_count"] == 1200
    assert row["duplicated_unexpected_count"] == 5
    assert row["protein_count"] == 25000
    assert row["consistent_pct"] == 96.0
    assert row["top_detected_clade"] == "Arabidopsis thaliana"
    assert row["top_detected_taxid"] == 3702
    assert row["detected_clade_count"] == 2


def test_summarize_omark_writes_header_only_when_no_sum_files_exist(tmp_path):
    omark_outdir = tmp_path / "omark"
    omark_outdir.mkdir()
    outfile = tmp_path / "omark_summary.tsv"

    completed = run_script("--omark_outdir", str(omark_outdir), "--outfile", str(outfile))

    assert completed.returncode == 0, completed.stderr
    df = pandas.read_csv(outfile, sep="\t")
    assert df.empty
    assert "species" in df.columns

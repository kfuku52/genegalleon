from pathlib import Path
import subprocess
import sys


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "summarize_inputprep_runs.py"


def run_script(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def write_summary_tsv(path):
    path.write_text(
        (
            "started_utc\tended_utc\tduration_sec\texit_code\tprovider\trun_build_manifest\trun_format_inputs\trun_validate_inputs\tstrict\toverwrite\tdownload_only\tdry_run\tdownload_timeout\tdataset_root\tinput_dir\tdownload_manifest\tdownload_dir\tmanifest_output\tspecies_cds_dir\tspecies_gff_dir\tnum_species_cds\tnum_species_gff\tstage_manifest_status\tstage_format_status\tstage_validate_status\tconfig_file\n"
            "2026-02-19T10:00:00Z\t2026-02-19T10:10:00Z\t600\t0\tall\t1\t1\t1\t0\t0\t0\t0\t120\t/a\t\t/m1\t/d1\t/o1\t/c1\t/g1\t100\t100\tok\tok\tok\t/p1\n"
            "2026-02-19T12:00:00Z\t2026-02-19T12:01:00Z\t60\t1\tphycocosm\t1\t1\t1\t1\t0\t0\t0\t120\t/a\t\t/m2\t/d2\t/o2\t/c2\t/g2\t50\t49\tok\tfailed\tskipped\t/p2\n"
            "2026-02-19T13:00:00Z\t2026-02-19T13:05:00Z\t300\t0\tensemblplants\t1\t1\t1\t0\t1\t0\t0\t120\t/a\t\t/m3\t/d3\t/o3\t/c3\t/g3\t120\t120\tok\tok\tok\t/p3\n"
        ),
        encoding="utf-8",
    )


def test_summarize_inputprep_runs_stdout(tmp_path):
    infile = tmp_path / "inputprep_runs.tsv"
    write_summary_tsv(infile)

    completed = run_script("--infile", str(infile), "--last-n", "2")
    assert completed.returncode == 0, completed.stderr
    text = completed.stdout
    assert "InputPrep Run Summary" in text
    assert "total_runs\t3" in text
    assert "successful_runs\t2" in text
    assert "failed_runs\t1" in text
    assert "latest_provider\tensemblplants" in text
    assert "Runs By Provider" in text
    assert "all\t1" in text
    assert "phycocosm\t1" in text
    assert "ensemblplants\t1" in text
    assert "Recent Runs (n=2)" in text
    assert "2026-02-19T12:00:00Z\t1\tphycocosm" in text
    assert "2026-02-19T13:00:00Z\t0\tensemblplants" in text


def test_summarize_inputprep_runs_outfile(tmp_path):
    infile = tmp_path / "inputprep_runs.tsv"
    write_summary_tsv(infile)
    outfile = tmp_path / "summary.txt"

    completed = run_script("--infile", str(infile), "--outfile", str(outfile))
    assert completed.returncode == 0, completed.stderr
    assert outfile.exists()
    out = outfile.read_text(encoding="utf-8")
    assert "InputPrep Run Summary" in out


def test_summarize_inputprep_runs_missing_file_fails(tmp_path):
    missing = tmp_path / "no_file.tsv"
    completed = run_script("--infile", str(missing))
    assert completed.returncode != 0
    assert "Input file not found" in completed.stderr

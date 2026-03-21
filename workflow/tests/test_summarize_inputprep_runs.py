from pathlib import Path
import subprocess
import sys


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "summarize_gg_input_generation_runs.py"


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
            "started_utc\tended_utc\tduration_sec\texit_code\tprovider\tinput_generation_mode\trun_format_inputs\trun_validate_inputs\trun_cds_fx2tab\trun_species_busco\trun_multispecies_summary\trun_generate_species_trait\tbusco_lineage\tbusco_lineage_resolved\tstrict\toverwrite\tdownload_only\tdry_run\tdownload_timeout\tinput_dir\tdownload_manifest\tdownload_dir\ttask_plan_output\tspecies_cds_dir\tspecies_cds_fx2tab_dir\tspecies_busco_full_dir\tspecies_busco_short_dir\tspecies_gff_dir\tspecies_genome_dir\tspecies_trait_output\tfile_multispecies_summary\tnum_species_cds\tnum_species_cds_fx2tab\tnum_species_gff\tnum_species_genome\tnum_species_busco_full\tnum_species_busco_short\tnum_species_trait\tnum_trait_columns\tcds_sequences_before\tcds_sequences_after\tcds_first_sequence_name\tstage_format_status\tstage_validate_status\tstage_cds_fx2tab_status\tstage_species_busco_status\tstage_multispecies_summary_status\tstage_trait_status\tconfig_file\n"
            "2026-02-19T10:00:00Z\t2026-02-19T10:10:00Z\t600\t0\tall\tsingle\t1\t1\t1\t1\t1\t0\tauto\teukaryota_odb12\t0\t0\t0\t0\t120\t/i1\t/m1\t/d1\t/t1\t/c1\t/cf1\t/bf1\t/bs1\t/g1\t/gn1\t/st1\t/ms1\t100\t100\t100\t100\t100\t100\t0\t0\t3000\t2900\tGenus_species_A\tok\tok\tok\tok\tok\tskipped\t/p1\n"
            "2026-02-19T12:00:00Z\t2026-02-19T12:01:00Z\t60\t1\tphycocosm\tsingle\t1\t1\t1\t1\t1\t1\tauto\teukaryota_odb12\t1\t0\t0\t0\t120\t/i2\t/m2\t/d2\t/t2\t/c2\t/cf2\t/bf2\t/bs2\t/g2\t/gn2\t/st2\t/ms2\t50\t49\t49\t40\t50\t50\t0\t0\t2000\t1800\tMicroglena_spYARC_gene1\tfailed\tskipped\tskipped\tskipped\tskipped\tfailed\t/p2\n"
            "2026-02-19T13:00:00Z\t2026-02-19T13:05:00Z\t300\t0\tensemblplants\tarray_finalize\t1\t1\t1\t1\t1\t1\tauto\teukaryota_odb12\t0\t1\t0\t0\t120\t/i3\t/m3\t/d3\t/t3\t/c3\t/cf3\t/bf3\t/bs3\t/g3\t/gn3\t/st3\t/ms3\t120\t120\t120\t120\t120\t120\t118\t4\t4000\t3950\tOstreococcus_lucimarinus_OSTLU_25062\tok\tok\tok\tok\tok\tok\t/p3\n"
        ),
        encoding="utf-8",
    )


def test_summarize_inputprep_runs_stdout(tmp_path):
    infile = tmp_path / "inputprep_runs.tsv"
    write_summary_tsv(infile)

    completed = run_script("--infile", str(infile), "--last-n", "2")
    assert completed.returncode == 0, completed.stderr
    text = completed.stdout
    assert "GG Input Generation Run Summary" in text
    assert "total_runs\t3" in text
    assert "successful_runs\t2" in text
    assert "failed_runs\t1" in text
    assert "latest_provider\tensemblplants" in text
    assert "latest_stage_status\tformat=ok;validate=ok;cds_fx2tab=ok;species_busco=ok;multispecies_summary=ok;trait=ok" in text
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
    assert "GG Input Generation Run Summary" in out


def test_summarize_inputprep_runs_missing_file_fails(tmp_path):
    missing = tmp_path / "no_file.tsv"
    completed = run_script("--infile", str(missing))
    assert completed.returncode != 0
    assert "Input file not found" in completed.stderr

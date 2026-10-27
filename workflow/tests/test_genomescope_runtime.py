import hashlib
import json
import os
import shutil
import subprocess
import zipfile
from pathlib import Path

from test_smudgeplot_runtime import make_diploid_reads

ROOT = Path(__file__).resolve().parents[2]
CORE = ROOT / "workflow/core/gg_genome_annotation_core.sh"


def run_stage(workspace, *, initial="12", policy="stop", histogram="", program=None, guided=False, source_mode="smudgeplot", k="21", smudge_k="21", enabled="1"):
    source = CORE.read_text()
    stage = 'task="GenomeScope"' + source.split('task="GenomeScope"', 1)[1].split('task="JCVI synteny dotplot"', 1)[0]
    if guided:
        stage = '# Read k-mer profiling:' + source.split('# Read k-mer profiling:', 1)[1].split('task="JCVI synteny dotplot"', 1)[0]
    task_tmp = workspace / "output/tmp/task"
    task_tmp.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ, gg_support_dir=str(ROOT / "workflow/support"),
                       gg_workspace_dir=str(workspace), gg_workspace_output_dir=str(workspace / "output"),
                       dir_sp_dnaseq=str(workspace / "input/species_dnaseq"), sp_ub="Test_species",
                       dir_sp_tmp=str(task_tmp), run_genomescope=enabled,
                       file_sp_genomescope=str(workspace / "output/species_dnaseq_genomescope/Test_species_genomescope.zip"),
                       annotation_provenance_dir=str(workspace / "output/artifact_provenance/genome_annotation"),
                       run_smudgeplot="0", smudgeplot_kmer_length=smudge_k, smudgeplot_lower_count="auto",
                       smudgeplot_aggregation_distance="2",
                       file_sp_smudgeplot=str(workspace / "output/species_dnaseq_smudgeplot/Test_species_smudgeplot.zip"),
                       genomescope_kmer_length=k, genomescope_kmer_lower="1", genomescope_kmer_upper="10000",
                       genomescope_guided="0", genomescope_coverage_source=source_mode,
                       genomescope_ploidy="2", genomescope_initial_coverage=initial, genomescope_max_kmercov="-1",
                       genomescope_histogram=histogram,
                       genomescope_program=program or os.environ.get("GG_TEST_GENOMESCOPE_PROGRAM", "genomescope.R"),
                       GG_TASK_CPUS="2", GG_MEM_TOOL_GB="2", artifact_stale_policy=policy)
    return subprocess.run(["bash", "-c", 'set -euo pipefail\nsource "$gg_support_dir/gg_util.sh"\n' + stage],
                          cwd=task_tmp, env=environment, capture_output=True, text=True, timeout=300)


def test_native_count_fit_refit_and_failed_publication(tmp_path):
    reads, compressed = make_diploid_reads(tmp_path, 12000)
    workspace = tmp_path / "workspace"
    inputs = workspace / "input/species_dnaseq/Test_species"
    inputs.mkdir(parents=True)
    (inputs / "linked.fastq.gz").symlink_to(compressed)
    first = run_stage(workspace)
    assert first.returncode == 0, first.stdout + first.stderr
    archive = workspace / "output/species_dnaseq_genomescope/Test_species_genomescope.zip"
    with zipfile.ZipFile(archive) as handle:
        prefix = "Test_species.genomescope/"
        metadata = json.loads(handle.read(prefix + "histogram_metadata.json"))
        assert metadata["ceiling_bin_distinct_kmers"] == 0
        assert metadata["total_kmer_observations"] == 24000 * (1000 - 21 + 1)
        assert handle.read(prefix + "transformed_linear_plot.png").startswith(b"\x89PNG")
        assert "Genome Haploid Length" in handle.read(prefix + "summary.txt").decode()
        old_histogram = tmp_path / "existing.histo"
        old_histogram.write_bytes(handle.read(prefix + "kmc.histo.tsv"))
    unchanged = run_stage(workspace)
    assert unchanged.returncode == 0, unchanged.stdout + unchanged.stderr
    assert "Skipped: GenomeScope" in unchanged.stdout
    original = hashlib.sha256(archive.read_bytes()).hexdigest()
    stale = run_stage(workspace, initial="11")
    assert stale.returncode != 0
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == original
    shutil.rmtree(inputs)
    refit = run_stage(workspace, initial="11", histogram=str(old_histogram), policy="rebuild")
    assert refit.returncode == 0, refit.stdout + refit.stderr
    with zipfile.ZipFile(archive) as handle:
        assert "-l 11" in handle.read(prefix + "command.txt").decode()
    final = archive.read_bytes()
    failed = run_stage(workspace, initial="10", histogram=str(old_histogram),
                       program="/missing/GenomeScope", policy="rebuild")
    assert failed.returncode != 0
    assert archive.read_bytes() == final
    provenance = workspace / "output/artifact_provenance/genome_annotation/Test_species.genomescope.json"
    final_manifest = provenance.read_bytes()
    no_fit = tmp_path / "no_fit"
    no_fit.write_text("#!/usr/bin/env bash\nexit 0\n")
    no_fit.chmod(0o755)
    unfit = run_stage(workspace, initial="10", histogram=str(old_histogram), program=str(no_fit), policy="rebuild")
    assert unfit.returncode != 0 and "primary fit failed" in unfit.stderr
    assert archive.read_bytes() == final and provenance.read_bytes() == final_manifest


def test_counter_boundary_retains_instead_of_dropping_high_copy_kmers(tmp_path):
    reads = tmp_path / "repeat.fastq"
    reads.write_text("@repeat\n" + "A" * 1000 + "\n+\n" + "I" * 1000 + "\n")
    database = tmp_path / "reads"
    subprocess.run(["kmc", "-k21", "-m2", "-t2", "-ci1", "-cx4294967295", "-cs100",
                    str(reads), str(database), str(tmp_path)], check=True, capture_output=True, timeout=60)
    output = tmp_path / "histogram"
    subprocess.run(["kmc_tools", "transform", str(database), "histogram", str(output), "-cx100"],
                   check=True, capture_output=True, timeout=60)
    nonzero = [line.split() for line in output.read_text().splitlines() if int(line.split()[1])]
    assert nonzero == [["100", "1"]]


def test_guided_native_default_sensitivity_cache_and_source_provenance(tmp_path):
    _, compressed = make_diploid_reads(tmp_path, 24000)
    workspace = tmp_path / "workspace"
    inputs = workspace / "input/species_dnaseq/Test_species"
    inputs.mkdir(parents=True)
    (inputs / "linked.fastq.gz").symlink_to(compressed)
    # Reject incompatible counting before running either native counter.
    mismatch = run_stage(workspace, initial="-1", guided=True, smudge_k="31")
    assert mismatch.returncode != 0 and "same reads and k" in mismatch.stderr
    assert not list((workspace / "output").rglob("*.zip"))
    first = run_stage(workspace, initial="-1", guided=True)
    assert first.returncode == 0, first.stdout + first.stderr
    archive = workspace / "output/species_dnaseq_genomescope/Test_species_genomescope.zip"
    smudge_archive = workspace / "output/species_dnaseq_smudgeplot/Test_species_smudgeplot.zip"
    with zipfile.ZipFile(archive) as handle:
        prefix = "Test_species.genomescope/"
        selection = json.loads(handle.read(prefix + "coverage_selection.json"))
        assert selection["selected_source"] == "smudgeplot"
        assert 20 <= selection["initial_coverage"] <= 28
        assert selection["ploidy"] == 2
        assert selection["smudgeplot"]["archive_sha256"] == hashlib.sha256(smudge_archive.read_bytes()).hexdigest()
        assert [fit["label"] for fit in selection["results"]] == ["primary", "half_H", "twice_H", "auto"]
        assert selection["results"][0]["status"] == "success"
        for result in selection["results"]:
            relative = "" if result["label"] == "primary" else "sensitivity/" + result["label"] + "/"
            assert handle.read(prefix + relative + "command.txt")
            assert prefix + relative + "fit.log" in handle.namelist()
        report = json.loads(handle.read(prefix + "report.json"))
        assert report["input_parameters"]["initial_kmer_coverage"] == selection["initial_coverage"]
    original = archive.read_bytes()
    provenance = workspace / "output/artifact_provenance/genome_annotation/Test_species.genomescope.json"
    manifest = provenance.read_bytes()
    second = run_stage(workspace, initial="-1", guided=True)
    assert second.returncode == 0, second.stdout + second.stderr
    assert "Skipped: Smudgeplot" in second.stdout and "Skipped: GenomeScope" in second.stdout
    assert archive.read_bytes() == original and provenance.read_bytes() == manifest
    # A missing ZIP with retained provenance must also restore its dependency
    # when the generic stale policy enables rebuilding a disabled stage.
    archive.unlink()
    smudge_archive.unlink()
    restored = run_stage(workspace, initial="-1", guided=True, policy="rebuild", enabled="0")
    assert restored.returncode == 0, restored.stdout + restored.stderr
    assert archive.is_file() and smudge_archive.is_file()
    # Manual coverage takes precedence even with a configured Smudgeplot k mismatch.
    manual = run_stage(workspace, initial="24", guided=True, smudge_k="31", policy="rebuild")
    assert manual.returncode == 0, manual.stdout + manual.stderr
    with zipfile.ZipFile(archive) as handle:
        selection = json.loads(handle.read(prefix + "coverage_selection.json"))
        assert selection["selected_source"] == "manual" and len(selection["results"]) == 1


def test_legacy_automatic_initialization_without_smudgeplot(tmp_path):
    _, compressed = make_diploid_reads(tmp_path, 12000)
    workspace = tmp_path / "workspace"
    inputs = workspace / "input/species_dnaseq/Test_species"
    inputs.mkdir(parents=True)
    (inputs / "linked.fastq.gz").symlink_to(compressed)
    result = run_stage(workspace, initial="-1", source_mode="auto")
    assert result.returncode == 0, result.stdout + result.stderr
    archive = workspace / "output/species_dnaseq_genomescope/Test_species_genomescope.zip"
    with zipfile.ZipFile(archive) as handle:
        selection = json.loads(handle.read("Test_species.genomescope/coverage_selection.json"))
        assert selection["selected_source"] == "auto" and selection["initial_coverage"] == -1
    assert not list((workspace / "output/species_dnaseq_smudgeplot").glob("*.zip"))

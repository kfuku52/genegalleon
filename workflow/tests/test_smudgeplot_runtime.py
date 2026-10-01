import gzip
import hashlib
import json
import os
import random
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "workflow/support/run_smudgeplot.py"
CORE = ROOT / "workflow/core/gg_genome_annotation_core.sh"


def make_diploid_reads(directory, reads_per_haplotype):
    rng = random.Random(19)
    genome = "".join(rng.choices("ACGT", k=1000000))
    alternative = list(genome)
    for index in range(41, len(alternative), 83):
        alternative[index] = rng.choice([base for base in "ACGT" if base != alternative[index]])
    reads = directory / "diploid.fastq"
    with reads.open("w") as handle:
        for haplotype_index, haplotype in enumerate((genome, "".join(alternative))):
            for index in range(reads_per_haplotype):
                start = rng.randrange(len(haplotype) - 1000)
                sequence = haplotype[start:start + 1000]
                handle.write(f"@{haplotype_index}_{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")
    compressed = directory / "diploid.fastq.gz"
    with reads.open("rb") as source, gzip.open(compressed, "wb") as destination:
        shutil.copyfileobj(source, destination)
    return reads, compressed


@pytest.fixture(scope="module")
def diploid_reads(tmp_path_factory):
    return make_diploid_reads(tmp_path_factory.mktemp("diploid"), 12000)


@pytest.fixture(scope="module")
def high_coverage_reads(tmp_path_factory):
    return make_diploid_reads(tmp_path_factory.mktemp("diploid_high_coverage"), 24000)


@pytest.mark.parametrize(("compressed", "lower"), [(False, "4"), (True, "4"), (True, "auto")])
def test_real_fastk_smudgeplot_pipeline(tmp_path, diploid_reads, high_coverage_reads, compressed, lower):
    dataset = high_coverage_reads if lower == "auto" else diploid_reads
    reads = dataset[int(compressed)]
    alias = tmp_path / ("linked.fastq.gz" if compressed else "linked.fastq")
    alias.symlink_to(reads)
    output = tmp_path / "result"
    arguments = [sys.executable, str(SCRIPT), "--reads", str(alias), "--output-dir", str(output),
                 "--lower-count", lower, "--threads", "2", "--memory-gb", "2"]
    if compressed:
        arguments.extend(["--database-dir", str(tmp_path / "database")])
    completed = subprocess.run(arguments, text=True, capture_output=True, timeout=300)
    log = (output / "commands.log").read_text() if output.exists() else ""
    assert completed.returncode == 0, completed.stderr + "\n" + log[-8000:]
    report = json.loads((output / "smudgeplot_smudgeplot_report.json").read_text())
    if lower == "auto":
        assert 20 < report["haploid_coverage"] < 28
    else:
        assert 8 < report["haploid_coverage"] < 16
    assert max(report["smudges"], key=lambda row: row["fraction"])["structure"] == "AB"
    assert (output / "smudgeplot_smudgeplot.pdf").read_bytes().startswith(b"%PDF-")
    assert (output / "smudgeplot_smudgeplot_log10.pdf").read_bytes().startswith(b"%PDF-")
    assert not list(output.glob("fastk-*"))
    metadata = json.loads((output / "run_metadata.json").read_text())
    expected_lower = 4
    if lower == "auto":
        cutoff = subprocess.run(["smudgeplot", "cutoff", str(output / "fastk.histo.tsv"), "L"],
                                capture_output=True, text=True, check=True, timeout=30)
        expected_lower = int(cutoff.stdout.strip())
    assert metadata["lower_count_used"] == expected_lower
    assert metadata["reads"][0]["path"] == str(reads.resolve())
    if compressed:
        assert (tmp_path / "database" / "reads.ktab").is_file()


def run_stage(workspace, *, lower="4", policy="stop", enabled="1"):
    # Execute the actual ordered stage in a temporary workspace, with the real
    # provenance and publication helpers. Other annotation stages need not run.
    source = CORE.read_text()
    stage = 'task="Smudgeplot"' + source.split('task="Smudgeplot"', 1)[1].split('task="JCVI synteny dotplot"', 1)[0]
    task_tmp = workspace / "output/tmp/task"
    task_tmp.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ, gg_support_dir=str(SCRIPT.parent),
                       gg_workspace_dir=str(workspace), gg_workspace_output_dir=str(workspace / "output"),
                       dir_sp_dnaseq=str(workspace / "input/species_dnaseq"), sp_ub="Test_species",
                       dir_sp_tmp=str(task_tmp), run_smudgeplot=enabled,
                       file_sp_smudgeplot=str(workspace / "output/species_dnaseq_smudgeplot/Test_species_smudgeplot.zip"),
                       annotation_provenance_dir=str(workspace / "output/artifact_provenance/genome_annotation"),
                       smudgeplot_kmer_length="21", smudgeplot_lower_count=lower,
                       smudgeplot_aggregation_distance="2", GG_TASK_CPUS="2", GG_MEM_TOOL_GB="2",
                       artifact_stale_policy=policy)
    return subprocess.run(["bash", "-c", 'set -euo pipefail\nsource "$gg_support_dir/gg_util.sh"\n' + stage],
                          cwd=task_tmp, env=environment, capture_output=True, text=True, timeout=300)


def test_stage_provenance_reruns_and_failed_rebuild_preserves_archive(tmp_path, diploid_reads):
    workspace = tmp_path / "workspace"
    external = tmp_path / "actual_reads"
    external.mkdir()
    reads = external / "diploid.fastq.gz"
    shutil.copyfile(diploid_reads[1], reads)
    (external / "notes.txt").write_text("ignored notes")
    inputs = workspace / "input/species_dnaseq"
    inputs.mkdir(parents=True)
    (inputs / "Test_species").symlink_to(external, target_is_directory=True)
    archive = workspace / "output/species_dnaseq_smudgeplot/Test_species_smudgeplot.zip"
    manifest = workspace / "output/artifact_provenance/genome_annotation/Test_species.smudgeplot.json"
    first = run_stage(workspace)
    assert first.returncode == 0, first.stdout + first.stderr
    with zipfile.ZipFile(archive) as handle:
        metadata = json.loads(handle.read("Test_species.smudgeplot/run_metadata.json"))
        assert metadata["lower_count_used"] == 4
        assert metadata["reads"][0]["path"] == str(reads.resolve())
        assert not any(".ktab" in name or "fastk-" in name for name in handle.namelist())
    assert manifest.is_file()
    original_digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    (external / "notes.txt").write_text("changed unrelated notes")
    unchanged = run_stage(workspace)
    assert unchanged.returncode == 0, unchanged.stdout + unchanged.stderr
    assert "Skipped: Smudgeplot" in unchanged.stdout
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == original_digest
    stale = run_stage(workspace, lower="5")
    assert stale.returncode != 0
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == original_digest
    rebuilt = run_stage(workspace, lower="5", policy="rebuild", enabled="0")
    assert rebuilt.returncode == 0, rebuilt.stdout + rebuilt.stderr
    with zipfile.ZipFile(archive) as handle:
        assert json.loads(handle.read("Test_species.smudgeplot/run_metadata.json"))["lower_count_used"] == 5
    successful_archive = archive.read_bytes()
    successful_manifest = manifest.read_bytes()
    # Broken discovery must stop before publication, retaining the last result.
    (external / "broken.fastq").symlink_to(external / "missing.fastq")
    failed = run_stage(workspace, policy="rebuild")
    assert failed.returncode != 0
    assert archive.read_bytes() == successful_archive
    assert manifest.read_bytes() == successful_manifest
    (external / "broken.fastq").unlink()
    # A real native-tool failure after a changed input must also preserve it.
    reads.write_bytes(gzip.compress(b"not FASTQ\n"))
    failed_count = run_stage(workspace, policy="rebuild")
    assert failed_count.returncode != 0
    assert archive.read_bytes() == successful_archive
    assert manifest.read_bytes() == successful_manifest
    # A directory with no FASTQ must not silently skip a requested rebuild.
    reads.rename(external / "reads.backup")
    missing_reads = run_stage(workspace, policy="rebuild")
    assert missing_reads.returncode != 0
    assert archive.read_bytes() == successful_archive
    assert manifest.read_bytes() == successful_manifest


@pytest.mark.parametrize("directory_exists", [False, True])
def test_stage_without_reads_skips_new_analysis(tmp_path, directory_exists):
    workspace = tmp_path / "workspace"
    if directory_exists:
        (workspace / "input/species_dnaseq/Test_species").mkdir(parents=True)
    result = run_stage(workspace)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Skipped: Smudgeplot" in result.stdout
    assert not list((workspace / "output").rglob("*.zip"))

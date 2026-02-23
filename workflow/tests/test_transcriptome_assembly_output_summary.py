from pathlib import Path
import subprocess
import sys

import pandas


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "workflow" / "support" / "transcriptome_assembly_output_summary.py"


def run_summary(args, cwd: Path):
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        cwd=str(cwd),
        capture_output=True,
        text=True,
        check=False,
    )
    return proc


def test_summary_uses_new_fastq_input_path(tmp_path: Path):
    workspace = tmp_path / "workspace"
    input_root = workspace / "input"
    output_dir = workspace / "output" / "transcriptome_assembly"

    (input_root / "species_rnaseq" / "Arabidopsis_thaliana").mkdir(parents=True)
    (output_dir / "longestcds").mkdir(parents=True)
    (output_dir / "longestcds" / "Arabidopsis_thaliana.fa.gz").write_text("", encoding="utf-8")

    out_tsv = tmp_path / "summary.tsv"
    proc = run_summary(
        [
            "--dir_transcriptome_assembly",
            str(output_dir),
            "--dir_pg_input",
            str(input_root),
            "--mode",
            "fastq",
            "--out",
            str(out_tsv),
        ],
        cwd=tmp_path,
    )

    assert proc.returncode == 0, f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
    df = pandas.read_csv(out_tsv, sep="\t", index_col=0)
    assert "Arabidopsis_thaliana" in df.index
    assert int(df.loc["Arabidopsis_thaliana", "longestcds"]) == 1


def test_summary_uses_new_metadata_file_suffix(tmp_path: Path):
    workspace = tmp_path / "workspace"
    input_root = workspace / "input"
    output_dir = workspace / "output" / "transcriptome_assembly"

    metadata_dir = input_root / "transcriptome_assembly" / "amalgkit_metadata"
    metadata_dir.mkdir(parents=True)
    (metadata_dir / "Homo_sapiens_metadata.tsv").write_text("run\tscientific_name\n", encoding="utf-8")
    (output_dir / "isoform").mkdir(parents=True)
    (output_dir / "isoform" / "Homo_sapiens.isoform.fa.gz").write_text("", encoding="utf-8")

    out_tsv = tmp_path / "summary.tsv"
    proc = run_summary(
        [
            "--dir_transcriptome_assembly",
            str(output_dir),
            "--dir_pg_input",
            str(input_root),
            "--mode",
            "metadata",
            "--out",
            str(out_tsv),
        ],
        cwd=tmp_path,
    )

    assert proc.returncode == 0, f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
    df = pandas.read_csv(out_tsv, sep="\t", index_col=0)
    assert "Homo_sapiens" in df.index
    assert int(df.loc["Homo_sapiens", "isoform"]) == 1


def test_summary_falls_back_to_output_derived_species_when_input_is_missing(tmp_path: Path):
    output_dir = tmp_path / "workspace" / "output" / "transcriptome_assembly"
    (output_dir / "longestcds").mkdir(parents=True)
    (output_dir / "longestcds" / "Mus_musculus.fa.gz").write_text("", encoding="utf-8")

    out_tsv = tmp_path / "summary.tsv"
    proc = run_summary(
        [
            "--dir_transcriptome_assembly",
            str(output_dir),
            "--mode",
            "auto",
            "--out",
            str(out_tsv),
        ],
        cwd=tmp_path,
    )

    assert proc.returncode == 0, f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
    df = pandas.read_csv(out_tsv, sep="\t", index_col=0)
    assert "Mus_musculus" in df.index
    assert int(df.loc["Mus_musculus", "longestcds"]) == 1


def test_summary_ignores_hidden_entries_in_input_and_output_dirs(tmp_path: Path):
    workspace = tmp_path / "workspace"
    input_root = workspace / "input"
    output_dir = workspace / "output" / "transcriptome_assembly"

    (input_root / "species_rnaseq" / "Arabidopsis_thaliana").mkdir(parents=True)
    (input_root / "species_rnaseq" / ".hidden_species").mkdir(parents=True)
    (output_dir / "longestcds").mkdir(parents=True)
    (output_dir / "longestcds" / "Arabidopsis_thaliana.fa.gz").write_text("", encoding="utf-8")
    (output_dir / "longestcds" / ".DS_Store").write_text("", encoding="utf-8")

    out_tsv = tmp_path / "summary.tsv"
    proc = run_summary(
        [
            "--dir_transcriptome_assembly",
            str(output_dir),
            "--dir_pg_input",
            str(input_root),
            "--mode",
            "fastq",
            "--out",
            str(out_tsv),
        ],
        cwd=tmp_path,
    )

    assert proc.returncode == 0, f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
    df = pandas.read_csv(out_tsv, sep="\t", index_col=0)
    assert "Arabidopsis_thaliana" in df.index
    assert ".hidden_species" not in df.index
    assert ".DS_Store" not in df.index

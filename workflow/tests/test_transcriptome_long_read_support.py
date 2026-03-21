from pathlib import Path
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
DETECT_SCRIPT = REPO_ROOT / "workflow" / "support" / "detect_amalgkit_read_technology.py"
RENAME_SCRIPT = REPO_ROOT / "workflow" / "support" / "rename_rnabloom_transcripts.py"


def _run_python(script: Path, *args: str, cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(script), *args],
        cwd=str(cwd),
        capture_output=True,
        text=True,
        check=False,
    )


def test_detect_amalgkit_read_technology_classifies_runs_and_summary(tmp_path: Path):
    metadata_path = tmp_path / "metadata.tsv"
    classification_path = tmp_path / "classification.tsv"
    summary_path = tmp_path / "summary.sh"

    metadata_path.write_text(
        "\n".join(
            [
                "run\tinstrument\tlib_layout\texp_title",
                "SRR1\tIllumina NovaSeq 6000\tsingle\tbulk RNA-seq",
                "SRR2\tPacBio Sequel II\tsingle\tIso-Seq cDNA",
                "SRR3\tOxford Nanopore MinION\tsingle\tcDNA transcriptome",
                "SRR4\tOxford Nanopore PromethION\tsingle\tdirect RNA transcriptome",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    completed = _run_python(
        DETECT_SCRIPT,
        "--metadata",
        str(metadata_path),
        "--classification-out",
        str(classification_path),
        "--summary-sh",
        str(summary_path),
        cwd=tmp_path,
    )

    assert completed.returncode == 0, f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"

    classification = classification_path.read_text(encoding="utf-8")
    assert "SRR1\tshort_read\tillumina\t0\t0\tsingle\tIllumina NovaSeq 6000" in classification
    assert "SRR2\tpacbio\tpacbio\t1\t0\tsingle\tPacBio Sequel II" in classification
    assert "SRR3\tont_cdna\tont\t1\t0\tsingle\tOxford Nanopore MinION" in classification
    assert "SRR4\tont_direct_rna\tont\t1\t1\tsingle\tOxford Nanopore PromethION" in classification

    summary = summary_path.read_text(encoding="utf-8")
    assert "detected_metadata_run_count=4" in summary
    assert "detected_short_read_run_count=1" in summary
    assert "detected_pacbio_run_count=1" in summary
    assert "detected_ont_cdna_run_count=1" in summary
    assert "detected_ont_direct_rna_run_count=1" in summary
    assert "detected_has_long_reads=1" in summary
    assert "detected_has_short_reads=1" in summary
    assert "detected_has_pacbio=1" in summary
    assert "detected_has_ont=1" in summary
    assert "detected_input_class=mixed_long_platforms" in summary
    assert "detected_metadata_instrument_field=instrument" in summary


def test_rename_rnabloom_transcripts_uses_corset_clusters(tmp_path: Path):
    input_fasta = tmp_path / "rnabloom.transcripts.fa"
    clusters = tmp_path / "clusters.txt"
    output_fasta = tmp_path / "renamed.fa"

    input_fasta.write_text(
        "\n".join(
            [
                ">transcriptA",
                "ATGAAATAG",
                ">transcriptB",
                "ATGCCCTAG",
                ">transcriptC",
                "ATGGGGTAG",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    clusters.write_text(
        "\n".join(
            [
                "transcriptA\tcluster_10",
                "transcriptB\tcluster_20",
                "transcriptC\tcluster_10",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    completed = _run_python(
        RENAME_SCRIPT,
        "--input-fasta",
        str(input_fasta),
        "--clusters",
        str(clusters),
        "--species-prefix",
        "Test_species",
        "--output-fasta",
        str(output_fasta),
        cwd=tmp_path,
    )

    assert completed.returncode == 0, f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"

    output = output_fasta.read_text(encoding="utf-8")
    assert ">Test_species_c000001-i000001\nATGAAATAG\n" in output
    assert ">Test_species_c000002-i000001\nATGCCCTAG\n" in output
    assert ">Test_species_c000001-i000002\nATGGGGTAG\n" in output

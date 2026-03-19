from pathlib import Path
import csv
import gzip
import json
import subprocess
import sys


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "validate_longest_cds_selection.py"


def run_script(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def write_gzip_text(path, text):
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def write_species_summary(path, rows):
    fieldnames = [
        "provider",
        "species_key",
        "species_prefix",
        "cds_input_path",
        "gff_input_path",
        "genome_input_path",
        "cds_output_path",
        "cds_sequences_before",
        "cds_sequences_after",
        "aggregated_cds_removed",
        "gene_grouping_mode",
    ]
    with open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def test_validate_longest_cds_selection_passes_for_raw_cds_inputs(tmp_path):
    raw_dir = tmp_path / "raw"
    cds_dir = tmp_path / "species_cds"
    raw_dir.mkdir()
    cds_dir.mkdir()

    raw_cds = raw_dir / "Arabidopsis_thaliana.raw.fa"
    raw_cds.write_text(
        (
            ">gene1.t1 gene=gene1\nATGAAA\n"
            ">gene1.t2 gene=gene1\nATGAAATTT\n"
            ">gene2.t1 gene=gene2\nATGCCC\n"
        ),
        encoding="utf-8",
    )

    formatted_cds = cds_dir / "Arabidopsis_thaliana_demo.fa.gz"
    write_gzip_text(
        formatted_cds,
        (
            ">Arabidopsis_thaliana_gene1\nATGAAATTT\n"
            ">Arabidopsis_thaliana_gene2\nATGCCC\n"
        ),
    )

    summary_path = tmp_path / "gg_input_generation_species.tsv"
    write_species_summary(
        summary_path,
        [
            {
                "provider": "direct",
                "species_key": "Arabidopsis_thaliana",
                "species_prefix": "Arabidopsis_thaliana",
                "cds_input_path": str(raw_cds),
                "gff_input_path": "",
                "genome_input_path": "",
                "cds_output_path": str(formatted_cds),
                "cds_sequences_before": "3",
                "cds_sequences_after": "2",
                "aggregated_cds_removed": "1",
            }
        ],
    )

    stats_path = tmp_path / "stats.json"
    completed = run_script(
        "--species-cds-dir",
        str(cds_dir),
        "--species-summary",
        str(summary_path),
        "--nthreads",
        "2",
        "--stats-output",
        str(stats_path),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert "[Arabidopsis_thaliana] Longest CDS validation OK:" in completed.stdout
    stats = json.loads(stats_path.read_text(encoding="utf-8"))
    assert stats["nthreads"] == 2
    assert stats["species_checked"] == 1
    assert stats["species_passed"] == 1
    assert stats["transcripts_total"] == 3
    assert stats["genes_total"] == 2
    assert stats["multi_isoform_genes_total"] == 1


def test_validate_longest_cds_selection_fails_when_shorter_isoform_is_kept(tmp_path):
    raw_dir = tmp_path / "raw"
    cds_dir = tmp_path / "species_cds"
    raw_dir.mkdir()
    cds_dir.mkdir()

    raw_cds = raw_dir / "Arabidopsis_thaliana.raw.fa"
    raw_cds.write_text(
        (
            ">gene1.t1 gene=gene1\nATGAAA\n"
            ">gene1.t2 gene=gene1\nATGAAATTT\n"
        ),
        encoding="utf-8",
    )

    formatted_cds = cds_dir / "Arabidopsis_thaliana_demo.fa.gz"
    write_gzip_text(
        formatted_cds,
        ">Arabidopsis_thaliana_gene1\nATGAAA\n",
    )

    summary_path = tmp_path / "gg_input_generation_species.tsv"
    write_species_summary(
        summary_path,
        [
            {
                "provider": "direct",
                "species_key": "Arabidopsis_thaliana",
                "species_prefix": "Arabidopsis_thaliana",
                "cds_input_path": str(raw_cds),
                "gff_input_path": "",
                "genome_input_path": "",
                "cds_output_path": str(formatted_cds),
                "cds_sequences_before": "2",
                "cds_sequences_after": "1",
                "aggregated_cds_removed": "1",
            }
        ],
    )

    completed = run_script(
        "--species-cds-dir",
        str(cds_dir),
        "--species-summary",
        str(summary_path),
    )
    assert completed.returncode == 1
    assert "length_mismatch=1 sample=Arabidopsis_thaliana_gene1" in completed.stderr
    assert "sequence_mismatch=1 sample=Arabidopsis_thaliana_gene1" in completed.stderr


def test_validate_longest_cds_selection_passes_for_gff_derived_inputs(tmp_path):
    raw_dir = tmp_path / "raw"
    cds_dir = tmp_path / "species_cds"
    raw_dir.mkdir()
    cds_dir.mkdir()

    gff_path = raw_dir / "Arabidopsis_thaliana.annotation.gff3"
    genome_path = raw_dir / "Arabidopsis_thaliana.genome.fa"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tID=gene1",
                "chr1\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=gene1.t1;Parent=gene1",
                "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cds1;Parent=gene1.t1",
                "chr1\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=gene1.t2;Parent=gene1",
                "chr1\tsrc\tCDS\t1\t12\t.\t+\t0\tID=cds2;Parent=gene1.t2",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(">chr1\nATGAAATTTCCC\n", encoding="utf-8")

    formatted_cds = cds_dir / "Arabidopsis_thaliana_demo.fa.gz"
    write_gzip_text(
        formatted_cds,
        ">Arabidopsis_thaliana_gene1\nATGAAATTTCCC\n",
    )

    summary_path = tmp_path / "gg_input_generation_species.tsv"
    write_species_summary(
        summary_path,
        [
            {
                "provider": "direct",
                "species_key": "Arabidopsis_thaliana",
                "species_prefix": "Arabidopsis_thaliana",
                "cds_input_path": "{} + {} (derived CDS)".format(gff_path, genome_path),
                "gff_input_path": str(gff_path),
                "genome_input_path": str(genome_path),
                "cds_output_path": str(formatted_cds),
                "cds_sequences_before": "2",
                "cds_sequences_after": "1",
                "aggregated_cds_removed": "1",
            }
        ],
    )

    completed = run_script(
        "--species-cds-dir",
        str(cds_dir),
        "--species-summary",
        str(summary_path),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert "[Arabidopsis_thaliana] Longest CDS validation OK:" in completed.stdout


def test_validate_longest_cds_selection_respects_rescue_overlap_mode(tmp_path):
    raw_dir = tmp_path / "raw"
    cds_dir = tmp_path / "species_cds"
    raw_dir.mkdir()
    cds_dir.mkdir()

    gff_path = raw_dir / "Arabidopsis_thaliana.annotation.gff3"
    genome_path = raw_dir / "Arabidopsis_thaliana.genome.fa"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t18\t.\t+\t.\tID=badGeneA",
                "chr1\tsrc\tmRNA\t1\t18\t.\t+\t.\tID=locusX.t1;Parent=badGeneA",
                "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cds1;Parent=locusX.t1",
                "chr1\tsrc\tCDS\t13\t18\t.\t+\t0\tID=cds2;Parent=locusX.t1",
                "chr1\tsrc\tgene\t1\t18\t.\t+\t.\tID=badGeneB",
                "chr1\tsrc\tmRNA\t1\t18\t.\t+\t.\tID=locusX.t2;Parent=badGeneB",
                "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=cds3;Parent=locusX.t2",
                "chr1\tsrc\tCDS\t13\t18\t.\t+\t0\tID=cds4;Parent=locusX.t2",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(">chr1\nATGAAACCCGGGTTTCCC\n", encoding="utf-8")

    formatted_cds = cds_dir / "Arabidopsis_thaliana_demo.fa.gz"
    write_gzip_text(
        formatted_cds,
        ">Arabidopsis_thaliana_locusX\nATGAAATTTCCC\n",
    )

    summary_path = tmp_path / "gg_input_generation_species.tsv"
    write_species_summary(
        summary_path,
        [
            {
                "provider": "direct",
                "species_key": "Arabidopsis_thaliana",
                "species_prefix": "Arabidopsis_thaliana",
                "cds_input_path": "{} + {} (derived CDS)".format(gff_path, genome_path),
                "gff_input_path": str(gff_path),
                "genome_input_path": str(genome_path),
                "cds_output_path": str(formatted_cds),
                "cds_sequences_before": "2",
                "cds_sequences_after": "1",
                "aggregated_cds_removed": "1",
                "gene_grouping_mode": "rescue_overlap",
            }
        ],
    )

    completed = run_script(
        "--species-cds-dir",
        str(cds_dir),
        "--species-summary",
        str(summary_path),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert "[Arabidopsis_thaliana] Longest CDS validation OK:" in completed.stdout

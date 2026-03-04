from pathlib import Path
import csv
import shutil
import subprocess
import sys


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "build_download_manifest.py"
SMALL_DATASET_ROOT = Path(__file__).resolve().parent / "data" / "small_gfe_dataset"


def run_script(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def read_manifest(path):
    with open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_build_download_manifest_from_small_fixture_all(tmp_path):
    out = tmp_path / "download_manifest.tsv"
    completed = run_script(
        "--provider",
        "all",
        "--dataset-root",
        str(SMALL_DATASET_ROOT),
        "--output",
        str(out),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert out.exists()

    rows = read_manifest(out)
    assert len(rows) == 3
    providers = sorted([row["provider"] for row in rows])
    assert providers == ["ensemblplants", "phycocosm", "phytozome"]

    for row in rows:
        assert row["cds_url"].startswith("file://")
        assert row["gff_url"].startswith("file://")
        assert "genome_url" in row
        assert row["cds_filename"] != ""
        assert row["gff_filename"] != ""
        assert "genome_filename" in row


def test_build_download_manifest_strict_fails_on_missing_pair(tmp_path):
    dataset_copy = tmp_path / "small_dataset_copy"
    shutil.copytree(SMALL_DATASET_ROOT, dataset_copy)
    missing_gff = dataset_copy / "PhycoCosm" / "species_wise_original" / "Microglena_spYARC_MicrYARC1" / "MicrYARC1_GeneCatalog_genes_20220803.gff3"
    missing_gff.unlink()

    out = tmp_path / "download_manifest.tsv"
    completed = run_script(
        "--provider",
        "phycocosm",
        "--dataset-root",
        str(dataset_copy),
        "--output",
        str(out),
        "--strict",
    )
    assert completed.returncode == 1, completed.stderr + "\n" + completed.stdout
    rows = read_manifest(out)
    assert len(rows) == 0


def test_build_download_manifest_ncbi_provider_from_species_dir_fixture(tmp_path):
    dataset_root = tmp_path / "dataset"
    ncbi_species_dir = dataset_root / "NCBI_Genome" / "species_wise_original" / "GCF_000001405.40"
    ncbi_species_dir.mkdir(parents=True)
    (ncbi_species_dir / "GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz").write_bytes(b"dummy")
    (ncbi_species_dir / "GCF_000001405.40_GRCh38.p14_genomic.gff.gz").write_bytes(b"dummy")
    (ncbi_species_dir / "GCF_000001405.40_GRCh38.p14_genomic.fna.gz").write_bytes(b"dummy")

    out = tmp_path / "download_manifest.tsv"
    completed = run_script(
        "--provider",
        "ncbi",
        "--dataset-root",
        str(dataset_root),
        "--output",
        str(out),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    rows = read_manifest(out)
    assert len(rows) == 1
    assert rows[0]["provider"] == "ncbi"
    assert rows[0]["id"] == "GCF_000001405.40"
    assert rows[0]["species_key"] == "GCF_000001405.40"
    assert rows[0]["genome_url"].startswith("file://")
    assert rows[0]["genome_filename"] == "GCF_000001405.40_GRCh38.p14_genomic.fna.gz"


def test_build_download_manifest_coge_and_cngb_provider_from_species_dir_fixture(tmp_path):
    dataset_root = tmp_path / "dataset"
    species_key = "Arabidopsis_thaliana"
    provider_dir_name = {"coge": "CoGe", "cngb": "CNGB"}

    for provider in ("coge", "cngb"):
        species_dir = dataset_root / provider_dir_name[provider] / "species_wise_original" / species_key
        species_dir.mkdir(parents=True, exist_ok=True)
        (species_dir / (species_key + ".cds.fa")).write_text(">AT1G01010_t1\nATGAA\n", encoding="utf-8")
        (species_dir / (species_key + ".gene.gff3")).write_text(
            "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=AT1G01010\n",
            encoding="utf-8",
        )
        (species_dir / (species_key + ".genome.fa")).write_text(">chr1\nATGCATGC\n", encoding="utf-8")

        out = tmp_path / ("download_manifest_{}.tsv".format(provider))
        completed = run_script(
            "--provider",
            provider,
            "--dataset-root",
            str(dataset_root),
            "--output",
            str(out),
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
        rows = read_manifest(out)
        assert len(rows) == 1
        assert rows[0]["provider"] == provider
        assert rows[0]["id"] == species_key
        assert rows[0]["species_key"] == species_key
        assert rows[0]["cds_filename"] == species_key + ".cds.fa"
        assert rows[0]["gff_filename"] == species_key + ".gene.gff3"
        assert rows[0]["genome_filename"] == species_key + ".genome.fa"

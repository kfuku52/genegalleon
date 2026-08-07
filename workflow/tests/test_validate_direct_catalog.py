import csv
import gzip
import subprocess
import sys
from pathlib import Path

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "validate_direct_catalog.py"
CURATED_CATALOG = Path(__file__).resolve().parents[2] / "workspace" / "input" / "input_generation" / "direct_catalog_curated.tsv"
FIELDS = ("provider", "id", "species_key", "cds_url", "gff_url", "gbff_url", "genome_url", "gff_filename")


def run_script(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def write_manifest(path, rows):
    with open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_curated_direct_catalog_passes_static_validation():
    completed = run_script("--manifest", str(CURATED_CATALOG))
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout


def test_catalog_rejects_transcriptome_as_cds_only(tmp_path):
    manifest = tmp_path / "catalog.tsv"
    write_manifest(
        manifest,
        [
            {
                "provider": "direct",
                "id": "onekp",
                "species_key": "Croton_tiglium",
                "cds_url": "https://example.test/SOAPdenovo-Trans-assembly.fa.gz",
            }
        ],
    )

    completed = run_script("--manifest", str(manifest))

    assert completed.returncode == 1
    assert "transcriptome assembly is not an accepted CDS-only bundle" in completed.stderr


def test_remote_catalog_check_rejects_complete_gff_without_cds(tmp_path):
    gff_path = tmp_path / "genes.gff.gz"
    with gzip.open(gff_path, "wt", encoding="utf-8") as handle:
        handle.write("chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=G1\n")
    genome_path = tmp_path / "genome.fa"
    genome_path.write_text(">chr1\nATG\n", encoding="utf-8")
    manifest = tmp_path / "catalog.tsv"
    write_manifest(
        manifest,
        [
            {
                "provider": "direct",
                "id": "fixture",
                "species_key": "Test_species",
                "gff_url": gff_path.as_uri(),
                "genome_url": genome_path.as_uri(),
                "gff_filename": gff_path.name,
            }
        ],
    )

    completed = run_script("--manifest", str(manifest), "--check-remote-gff")

    assert completed.returncode == 1
    assert "complete GFF contained no CDS feature" in completed.stderr

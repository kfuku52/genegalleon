from pathlib import Path
import csv
import shutil
import subprocess
import sys


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "script" / "build_download_manifest.py"
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
        assert row["cds_filename"] != ""
        assert row["gff_filename"] != ""


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

from pathlib import Path
import csv
import os
import subprocess
import sys


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "format_species_inputs.py"
SMALL_DATASET_ROOT = Path(__file__).resolve().parent / "data" / "small_gfe_dataset"


def run_script(*args, env=None):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), *args],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )


def make_manifest(path, rows):
    with open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "provider",
                "id",
                "species_key",
                "cds_url",
                "gff_url",
                "genome_url",
                "cds_filename",
                "gff_filename",
                "genome_filename",
                "local_cds_path",
                "local_gff_path",
                "local_genome_path",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            payload = dict(row)
            if str(payload.get("id", "") or "").strip() == "":
                payload["id"] = str(payload.get("species_key", "") or "").strip()
            writer.writerow(payload)


def to_file_url(path):
    return path.resolve().as_uri()


def test_format_dry_run_does_not_write_outputs(tmp_path):
    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    completed = run_script(
        "--provider",
        "all",
        "--input-dir",
        str(SMALL_DATASET_ROOT),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--dry-run",
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert out_cds.exists()
    assert out_gff.exists()
    assert list(out_cds.iterdir()) == []
    assert list(out_gff.iterdir()) == []
    assert "dry-run" in completed.stdout


def test_download_dry_run_does_not_fetch_files(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    cds_source = source_dir / "cds.fa"
    gff_source = source_dir / "gene.gff3"
    cds_source.write_text(">a\nATG\n", encoding="utf-8")
    gff_source.write_text("scaf\tsrc\tgene\t1\t3\t.\t+\t.\tID=a\n", encoding="utf-8")

    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "coge",
                "id": "24739",
                "species_key": "Arabidopsis_thaliana",
                "cds_url": to_file_url(cds_source),
                "gff_url": to_file_url(gff_source),
                "cds_filename": "Arabidopsis_thaliana.cds.fa",
                "gff_filename": "Arabidopsis_thaliana.gff3",
            }
        ],
    )

    download_dir = tmp_path / "download_cache"
    completed = run_script(
        "--provider",
        "coge",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(download_dir),
        "--download-only",
        "--dry-run",
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    raw_dir = download_dir / "CoGe" / "species_wise_original" / "Arabidopsis_thaliana"
    assert raw_dir.exists()
    assert (raw_dir / "Arabidopsis_thaliana.cds.fa").exists() is False
    assert (raw_dir / "Arabidopsis_thaliana.gff3").exists() is False
    assert "planned downloads=2" in completed.stdout


def test_invalid_http_header_fails_fast(tmp_path):
    completed = run_script(
        "--provider",
        "ensemblplants",
        "--input-dir",
        str(SMALL_DATASET_ROOT / "20230216_EnsemblPlants" / "original_files"),
        "--http-header",
        "INVALID_HEADER_WITHOUT_COLON",
    )
    assert completed.returncode != 0
    assert "Invalid --http-header value" in completed.stderr


def test_auth_bearer_token_env_missing_fails_fast(tmp_path):
    completed = run_script(
        "--provider",
        "ensemblplants",
        "--input-dir",
        str(SMALL_DATASET_ROOT / "20230216_EnsemblPlants" / "original_files"),
        "--auth-bearer-token-env",
        "THIS_ENV_SHOULD_NOT_EXIST_FOR_TEST",
    )
    assert completed.returncode != 0
    assert "bearer token is empty or undefined" in completed.stderr


def test_auth_bearer_token_env_allows_download(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    cds_source = source_dir / "cds.fa"
    gff_source = source_dir / "gene.gff3"
    cds_source.write_text(">a\nATG\n", encoding="utf-8")
    gff_source.write_text("scaf\tsrc\tgene\t1\t3\t.\t+\t.\tID=a\n", encoding="utf-8")

    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "coge",
                "id": "24739",
                "species_key": "Arabidopsis_thaliana",
                "cds_url": to_file_url(cds_source),
                "gff_url": to_file_url(gff_source),
                "cds_filename": "Arabidopsis_thaliana.cds.fa",
                "gff_filename": "Arabidopsis_thaliana.gff3",
            }
        ],
    )

    env = dict(os.environ)
    env["TEST_BEARER_TOKEN"] = "dummy-token"
    download_dir = tmp_path / "download_cache"
    completed = run_script(
        "--provider",
        "coge",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(download_dir),
        "--download-only",
        "--auth-bearer-token-env",
        "TEST_BEARER_TOKEN",
        env=env,
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    raw_dir = download_dir / "CoGe" / "species_wise_original" / "Arabidopsis_thaliana"
    assert (raw_dir / "Arabidopsis_thaliana.cds.fa").exists()
    assert (raw_dir / "Arabidopsis_thaliana.gff3").exists()

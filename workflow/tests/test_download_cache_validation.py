import gzip
import json
import sys
from pathlib import Path

import pytest

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
if str(SUPPORT_DIR) not in sys.path:
    sys.path.insert(0, str(SUPPORT_DIR))

from format_species_download.cache_validation import (  # noqa: E402
    GzipValidationCache,
    gzip_validation_key_for_target,
)
from format_species_download.local import (  # noqa: E402
    quarantine_corrupt_gzip,
    validate_gzip_with_cache,
)
from format_species_download.manifest import download_from_manifest  # noqa: E402


def _write_gzip(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def test_gzip_validation_receipt_skips_unchanged_full_read(tmp_path, monkeypatch):
    target_root = tmp_path / "staged" / "plan-hash"
    target = target_root / "Direct" / "species_wise_original" / "Good_species" / "Good_species.cds.fa.gz"
    _write_gzip(target, ">gene1\nATG\n")
    cache = GzipValidationCache(tmp_path / "staged" / ".gg-gzip-validation")
    key = gzip_validation_key_for_target(target, target_root, "https://example.test/good.cds.gz")

    assert validate_gzip_with_cache(target, validation_cache=cache, validation_key=key) is None
    receipt = tmp_path / "staged" / ".gg-gzip-validation" / (key + ".json")
    assert receipt.exists()
    assert json.loads(receipt.read_text(encoding="utf-8"))["validation_key"] == key

    def fail_if_read(*_args, **_kwargs):
        raise AssertionError("unchanged gzip should use its validation receipt")

    monkeypatch.setattr("format_species_download.local.gzip.open", fail_if_read)
    assert validate_gzip_with_cache(target, validation_cache=cache, validation_key=key) is None
    assert cache.diagnostics()["validation_cache_hits"] == 1


def test_gzip_validation_receipt_invalidates_on_file_change_and_corruption(tmp_path):
    target_root = tmp_path / "staged" / "plan-hash"
    target = target_root / "Direct" / "species_wise_original" / "Good_species" / "Good_species.cds.fa.gz"
    _write_gzip(target, ">gene1\nATG\n")
    cache = GzipValidationCache(tmp_path / "staged" / ".gg-gzip-validation")
    key = gzip_validation_key_for_target(target, target_root, "https://example.test/good.cds.gz")
    assert validate_gzip_with_cache(target, validation_cache=cache, validation_key=key) is None

    _write_gzip(target, ">gene1\nATGAAATTT\n")
    assert validate_gzip_with_cache(target, validation_cache=cache, validation_key=key) is None
    target.write_bytes(b"\x1f\x8b\x08\x00truncated")
    warnings = []
    assert quarantine_corrupt_gzip(
        target,
        warnings,
        "[test]",
        validation_cache=cache,
        validation_key=key,
    )
    assert not target.exists()
    assert list(target.parent.glob(target.name + ".corrupt.*"))


def test_download_manifest_reuses_validation_receipt(tmp_path):
    source = tmp_path / "source.cds.fa.gz"
    _write_gzip(source, ">gene1\nATG\n")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "provider\tid\tspecies_key\tcds_url\tcds_filename\n"
        "direct\tfixture\tGood_species\t{}\tGood_species.cds.fa.gz\n".format(source.as_uri()),
        encoding="utf-8",
    )
    download_root = tmp_path / "downloads"
    kwargs = {
        "manifest_path": manifest,
        "download_root": download_root,
        "provider_filter": "direct",
        "overwrite": False,
        "headers": {},
        "timeout": 10,
        "dry_run": False,
        "jobs": 1,
    }
    first = download_from_manifest(**kwargs)
    assert first["errors"] == []
    assert first["download_diagnostics"]["validation_cache_records"] == 1
    second = download_from_manifest(**kwargs)
    assert second["errors"] == []
    assert second["download_diagnostics"]["validation_cache_hits"] == 1
    assert second["download_diagnostics"]["download_jobs"] == 0


def test_staging_preserves_successful_task_receipts_on_partial_failure(tmp_path, monkeypatch):
    import format_species_inputs as fsi
    import stage_input_generation_downloads as stage_module

    download_dir = tmp_path / "downloads"
    plan = tmp_path / "task_plan.json"
    tasks = []
    for species in ("Good_species", "Bad_species"):
        tasks.append(
            {
                "provider": "direct",
                "species_key": species,
                "species_prefix": species,
                "manifest_row": {
                    "provider": "direct",
                    "id": species,
                    "species_key": species,
                },
                "download_dir": str(download_dir),
                "gene_grouping_mode": "rescue_overlap",
                "gff_repair_mode": "safe",
                "format_strict": False,
                "input_sha256": {},
            }
        )
    plan.write_text(
        json.dumps(
            {
                "version": 2,
                "provider": "direct",
                "task_count": len(tasks),
                "species": [task["species_prefix"] for task in tasks],
                "tasks": tasks,
                "download_mode": "staged",
            }
        ),
        encoding="utf-8",
    )

    def fake_download_from_manifest(**kwargs):
        root = Path(kwargs["download_root"])
        species_dir = root / "Direct" / "species_wise_original" / "Good_species"
        species_dir.mkdir(parents=True)
        paths = {
            "cds_path": species_dir / "Good_species.cds.fa",
            "gff_path": species_dir / "Good_species.gff",
            "genome_path": species_dir / "Good_species.genome.fa",
        }
        for path in paths.values():
            path.write_text(">gene1\nATG\n", encoding="utf-8")
        return {
            "warnings": [],
            "errors": ["failed Bad_species"],
            "downloaded": 3,
            "resolved_rows": [{"provider": "direct", "species_key": "Good_species"}],
        }

    def fake_discover_tasks(_provider, input_dir, **_kwargs):
        species_dir = Path(input_dir) / "Good_species"
        return [
            {
                "provider": "direct",
                "species_key": "Good_species",
                "species_prefix": "Good_species",
                "cds_path": species_dir / "Good_species.cds.fa",
                "gff_path": species_dir / "Good_species.gff",
                "genome_path": species_dir / "Good_species.genome.fa",
            }
        ], [], []

    monkeypatch.setattr(fsi, "download_from_manifest", fake_download_from_manifest)
    monkeypatch.setattr(fsi, "discover_tasks", fake_discover_tasks)
    with pytest.raises(ValueError, match=r"Staged 1 of 2 pending species"):
        stage_module.stage_downloads(plan, jobs=2, headers={})

    assert (Path(str(plan) + ".tasks") / "1.json").exists()
    assert not (Path(str(plan) + ".tasks") / "2.json").exists()

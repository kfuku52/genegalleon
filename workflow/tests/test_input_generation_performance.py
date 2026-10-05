import argparse
import hashlib
import json
import os
import struct
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
import input_generation_stage_resume as resume
import performance_metrics as metrics
import version_report_identity as identity
from shared_namespace_lock import namespace_lock

from workflow.benchmarks.benchmark_stage_resume import fixture


@pytest.fixture
def native_import(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["fixture"])
    _, (source, target) = fixture(Path(resume.__file__).parent, tmp_path, 1)
    return argparse.Namespace(task_plan=target[0], root=target[1], source_plan=source[0],
                              source_root=source[1], source_plan_sha256=resume.digest(source[0]),
                              format_contract_version="24", task_index=1, source_only=False)


def test_current_format_proof_does_not_read_completion_receipt(native_import, monkeypatch):
    monkeypatch.setattr(resume, "verify_receipt", lambda *a, **k: pytest.fail("Unneeded completion read"))
    resume.import_stages(native_import)
    for stage in ("format", "validate"):
        assert resume.valid(native_import.task_plan, 1, native_import.root, stage, "24")


def test_prepare_checks_only_source_contract(native_import, monkeypatch):
    native_import.source_only = True
    monkeypatch.setattr(resume, "snapshot", lambda *a, **k: pytest.fail("Prepare imported a species"))
    resume.import_stages(native_import)
    assert not (native_import.root / "tmp/task_meta_shards/1.json").exists()


def test_worker_import_selects_only_its_task(native_import, capsys):
    plan = json.loads(native_import.task_plan.read_text())
    plan["tasks"].append({**plan["tasks"][0], "species_prefix": "Other_species", "species_key": "Other_species"})
    plan["task_count"] = 2
    resume.atomic_json(native_import.task_plan, plan)
    owner_path = native_import.root / ".array-plan.json"
    owner = json.loads(owner_path.read_text())
    owner["plan_sha256"] = resume.digest(native_import.task_plan)
    resume.atomic_json(owner_path, owner)
    resume.import_stages(native_import)
    report = json.loads(capsys.readouterr().out.splitlines()[-1])
    assert report["without_verified_format"] == []
    assert [item["species"] for item in report["imported"]] == ["Example_species"]
    assert not (native_import.root / "tmp/task_meta_shards/2.json").exists()


def test_import_allows_other_readers_but_excludes_same_species_writer(native_import):
    donor_lock = Path(str(native_import.source_plan) + ".locks/1.lock")
    with namespace_lock(native_import.source_root / ".array-phase.lock", exclusive=False):
        with namespace_lock(donor_lock, exclusive=True):
            with pytest.raises(ValueError, match="active workers"):
                resume.import_stages(native_import)
        resume.import_stages(native_import)


def test_target_duplicate_and_foreign_owner_are_rejected(native_import):
    lock = Path(str(native_import.task_plan) + ".locks/1.lock")
    with namespace_lock(lock, exclusive=True):
        with pytest.raises(ValueError, match="active owner"):
            resume.import_stages(native_import)
        native_import.target_lock_token = "0" * 32
        with pytest.raises(ValueError, match="ownership changed"):
            resume.import_stages(native_import)


@pytest.mark.parametrize("mutation", ["source", "destination"])
def test_corruption_during_copy_cannot_publish_checkpoint(native_import, monkeypatch, mutation):
    original = resume.copy_atomic
    def corrupt(source, destination, **kwargs):
        if mutation == "source":
            Path(source).write_text("corrupted bytes\n")
        original(source, destination, **kwargs)
        if mutation == "destination":
            Path(destination).write_text("corrupted bytes\n")
    monkeypatch.setattr(resume, "copy_atomic", corrupt)
    with pytest.raises(ValueError, match="differs from its verified source"):
        resume.import_stages(native_import)
    assert not resume.checkpoint_path(native_import.root, "Example_species", "format").exists()


def test_mutable_image_is_fully_hashed_after_same_size_change(tmp_path):
    image = tmp_path / "image.sif"
    image.write_bytes(b"original")
    before = image.stat()
    first = identity.image_identity(image)
    image.write_bytes(b"modified")
    os.utime(image, ns=(before.st_atime_ns, before.st_mtime_ns))
    assert identity.image_identity(image) != first
    assert identity.image_identity(image) == "sha256:" + hashlib.sha256(b"modified").hexdigest()


def test_kernel_verified_digest_avoids_full_read(tmp_path, monkeypatch):
    image = tmp_path / "image.sif"
    image.write_bytes(b"original")
    def measured(fd, operation, buffer, mutate):
        assert operation == 0xC0046686
        buffer[:36] = struct.pack("=HH", 1, 32) + bytes.fromhex("ab" * 32)
    monkeypatch.setattr(identity.fcntl, "ioctl", measured)
    monkeypatch.setattr(identity.sys, "platform", "linux")
    monkeypatch.setattr(identity, "digest", lambda *a: pytest.fail("Read immutable image"))
    assert identity.image_identity(image) == "fsverity:1:" + "ab" * 32


def test_cache_context_tracks_tool_environment_without_exposing_values(monkeypatch):
    before = identity.environment_identity()
    monkeypatch.setenv("GG_JOB_ID", "diagnostic-job")
    assert identity.environment_identity() == before
    monkeypatch.setenv("PYTHONPATH", "fixture-private-value")
    changed = identity.environment_identity()
    assert changed != before and "fixture-private-value" not in changed


def test_cache_rejects_tampering_and_symlink(tmp_path):
    log = tmp_path / "inventory.log"
    log.write_text("verified inventory")
    receipt = {"schema_version": 1, "key": "fixture", "sha256": identity.digest(log)}
    log.with_suffix(".json").write_text(json.dumps(receipt))
    assert identity.cache_ready(log, "fixture")
    log.write_text("corrupt inventory")
    assert not identity.cache_ready(log, "fixture")
    alias = tmp_path / "alias.log"
    alias.symlink_to(log)
    assert not identity.cache_ready(alias, "fixture")


def test_performance_records_are_optional_and_advisory(tmp_path, monkeypatch):
    monkeypatch.setenv("GG_PERFORMANCE_DIR", str(tmp_path))
    with metrics.measure("fixture"):
        metrics.count("sha256_bytes", 123)
    rows = [json.loads(line) for line in (tmp_path / f"{os.getpid()}.jsonl").read_text().splitlines()]
    assert rows[-1]["phase"] == "fixture" and rows[-1]["counters"]["sha256_bytes"] == 123
    assert rows[-1]["status"] == "ok" and rows[-1]["process_peak_rss_bytes"] > 0

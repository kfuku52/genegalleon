import argparse
import hashlib
import io
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


@pytest.mark.parametrize("label", ["cds", "gff", "stats", "summary", "ownership_qc", "mapping_qc"])
def test_late_destination_change_cannot_certify_imported_validation(native_import, monkeypatch, label):
    original = resume.copy_atomic
    def corrupt_after_qc(source, destination, **kwargs):
        original(source, destination, **kwargs)
        if str(destination).endswith(".mapping.json"):
            path = resume.context(native_import.task_plan, 1, native_import.root, "validate")[3][label]
            Path(path).write_text("unvalidated bytes\n")
    monkeypatch.setattr(resume, "copy_atomic", corrupt_after_qc)
    with pytest.raises(ValueError, match="differs from its verified source"):
        resume.import_stages(native_import)
    assert not resume.checkpoint_path(native_import.root, "Example_species", "validate").exists()


def test_donor_output_alias_cannot_be_written_under_a_reader_lock(native_import):
    settings_path = Path(str(native_import.task_plan) + ".settings.json")
    settings = json.loads(settings_path.read_text())
    donor_settings = json.loads(Path(str(native_import.source_plan) + ".settings.json").read_text())
    alias = native_import.root.parent.parent / "aliased-cds"
    alias.symlink_to(donor_settings["species_cds_dir"], target_is_directory=True)
    settings["species_cds_dir"] = str(alias)
    resume.atomic_json(settings_path, settings)
    cds = Path(resume.context(native_import.source_plan, 1, native_import.source_root, "format")[3]["cds"])
    before = cds.stat()
    with pytest.raises(ValueError, match="overlap"):
        resume.import_stages(native_import)
    assert cds.stat().st_ino == before.st_ino
    assert not resume.checkpoint_path(native_import.root, "Example_species", "format").exists()


def test_hard_link_to_raw_input_cannot_be_published(native_import):
    settings = json.loads(Path(str(native_import.task_plan) + ".settings.json").read_text())
    source_paths = resume.context(native_import.source_plan, 1, native_import.source_root, "format")[3]
    target = Path(settings["species_cds_dir"]) / Path(source_paths["cds"]).name
    target.parent.mkdir(parents=True)
    os.link(source_paths["cds_path"], target)
    before = target.stat()
    with pytest.raises(ValueError, match="overlap"):
        resume.import_stages(native_import)
    assert target.stat().st_ino == before.st_ino


def test_unwritten_directory_setting_does_not_block_resume(native_import):
    source_settings = json.loads(Path(str(native_import.source_plan) + ".settings.json").read_text())
    settings_path = Path(str(native_import.task_plan) + ".settings.json")
    target_settings = json.loads(settings_path.read_text())
    shared_unused = str(native_import.root.parent / "unused-fx2tab")
    for path, settings in ((Path(str(native_import.source_plan) + ".settings.json"), source_settings),
                           (settings_path, target_settings)):
        settings["species_cds_fx2tab_dir"] = shared_unused
        resume.atomic_json(path, settings)
    prepared_path = Path(str(native_import.source_plan) + ".prepared.json")
    prepared = json.loads(prepared_path.read_text())
    prepared["settings_sha256"] = resume.digest(Path(str(native_import.source_plan) + ".settings.json"))
    resume.atomic_json(prepared_path, prepared)
    resume.import_stages(native_import)
    assert resume.valid(native_import.task_plan, 1, native_import.root, "validate", "24")


def test_task_metadata_alias_cannot_overwrite_donor(native_import):
    source_dir = native_import.source_root / "tmp/task_meta_shards"
    (native_import.root / "tmp/task_meta_shards").symlink_to(source_dir, target_is_directory=True)
    before = (source_dir / "1.json").read_bytes()
    with pytest.raises(ValueError, match="overlap"):
        resume.import_stages(native_import)
    assert (source_dir / "1.json").read_bytes() == before


def test_summary_parser_uses_the_verified_bytes_even_if_source_is_restored(native_import, monkeypatch):
    summary = Path(resume.context(native_import.source_plan, 1, native_import.source_root, "format")[3]["summary"])
    original_bytes = summary.read_bytes()
    original_reader = resume.csv.DictReader
    def change_then_restore(handle, **kwargs):
        summary.write_text("species_prefix\tcds_output_path\nExample_species\tunverified-path\n")
        captured = handle.read()
        summary.write_bytes(original_bytes)
        return original_reader(io.StringIO(captured), **kwargs)
    monkeypatch.setattr(resume.csv, "DictReader", change_then_restore)
    resume.import_stages(native_import)
    paths = resume.context(native_import.task_plan, 1, native_import.root, "format")[3]
    with Path(paths["summary"]).open() as handle:
        row = next(original_reader(handle, delimiter="\t"))
    assert row["cds_output_path"] == paths["cds"]
    assert summary.read_bytes() == original_bytes


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


def test_fresh_batch_deduplicates_hardlinks_and_rejects_late_change(tmp_path, monkeypatch):
    import input_generation_array_state as state
    first, alias, other = (tmp_path / name for name in ("first", "alias", "other"))
    first.write_bytes(b"content")
    os.link(first, alias)
    other.write_bytes(b"other")
    original = state.digest
    reads = []
    def counted(path):
        reads.append(str(path))
        return original(path)
    monkeypatch.setattr(state, "digest", counted)
    batch = state.FreshDigestBatch()
    assert batch.read([first, alias])[str(first)] == original(first)
    batch.read([alias, other])
    assert len(reads) == 2
    info = first.stat()
    first.write_bytes(b"CONTENT")
    os.utime(first, ns=(info.st_atime_ns, info.st_mtime_ns + 2_000_000_000))
    with pytest.raises(OSError, match="changed while hashing"):
        batch.read([other])
    # A separate boundary reads the modified bytes instead of an old digest.
    assert state.digest_paths([first])[str(first)] == original(first)


@pytest.mark.parametrize("changed", ["own", "other", "shared", "unknown", "plan", "settings"])
def test_task_prepared_scope_still_checks_shared_and_sealed_files(tmp_path, changed):
    import input_generation_array_state as state
    plan = tmp_path / "plan.json"
    state.atomic_json(plan, {"task_count": 3, "tasks": [{"species_prefix": f"Species_{i}"} for i in range(3)]})
    settings = Path(str(plan) + ".settings.json")
    settings.write_text("{}")
    staging = Path(str(plan) + ".tasks")
    staging.mkdir()
    files = [staging / f"{index}{suffix}" for index in (1, 2, 3) for suffix in (".json", ".resolved.tsv")]
    shared, unknown = tmp_path / "lineage", staging / "shared.json"
    files.extend([shared, unknown])
    for path in files:
        path.write_text("sealed")
    state.atomic_json(str(plan) + ".prepared.json", {
        "plan_sha256": state.digest(plan), "settings_sha256": state.digest(settings),
        "files": {str(path): state.digest(path) for path in files}})
    assert state.prepared(plan, task_index=1)
    assert not state.prepared(plan, task_index=4)
    path = {"own": staging / "1.json", "other": staging / "2.resolved.tsv", "shared": shared,
            "unknown": unknown, "plan": plan, "settings": settings}[changed]
    path.write_text("changed")
    assert state.prepared(plan, task_index=1) == (changed == "other")
    assert not state.prepared(plan)


@pytest.mark.parametrize("bad", ["none", "json", "nul"])
def test_shell_metadata_batch_preserves_values_and_propagates_errors(tmp_path, bad):
    import subprocess
    core = (Path(resume.__file__).parents[1] / "core/gg_input_generation_core.sh").read_text()
    helper = core.split("read_stats_json_fields() {", 1)[1].split("\ninput_generation_effective_input_dir_path()", 1)[0]
    payload = {"first": None, "second": "space ' quote \\ tab\tline\nend\n", "third": True,
               "fourth": "$(touch injected);`false`"}
    if bad == "nul":
        payload["second"] = "bad\0value"
    path = tmp_path / "stats.json"
    path.write_text("invalid" if bad == "json" else json.dumps(payload))
    script = "read_stats_json_fields() {" + helper + '\ndownload_tmp_root="$2"\n'
    script += 'read_stats_json_fields "$1" a=first b=second c=third d=fourth e=missing || exit $?\n'
    script += 'printf "%s\\0" "$a" "$b" "$c" "$d" "$e"\n'
    result = subprocess.run(["bash", "-c", script, "metadata-test", str(path), str(tmp_path)], capture_output=True)
    if bad != "none":
        assert result.returncode != 0
    else:
        assert result.returncode == 0, result.stderr
        assert result.stdout.split(b"\0") == [b"", payload["second"].rstrip("\n").encode(), b"True",
                                               payload["fourth"].encode(), b"", b""]
    assert not (tmp_path / "injected").exists()
    assert not list(tmp_path.glob("json-fields.*"))


def test_import_binds_generated_summary_bytes(native_import, monkeypatch):
    original = Path.write_bytes
    def change_summary(path, data):
        result = original(path, data)
        if path == native_import.root / "tmp/species_summary_shards/1.tsv":
            original(path, data.replace(b"Example_species", b"Changed_species"))
        return result
    monkeypatch.setattr(Path, "write_bytes", change_summary)
    with pytest.raises(ValueError, match="differs from its verified source"):
        resume.import_stages(native_import)
    assert not resume.checkpoint_path(native_import.root, "Example_species", "format").exists()


def test_final_qc_session_rejects_change_between_validators(tmp_path, monkeypatch):
    import input_generation_array_state as state
    import input_validation_reuse as reuse
    import validate_cds_gff_mapping as mapping
    import validate_longest_cds_selection as ownership
    source = tmp_path / "genome.fa"
    source.write_text("original")
    output_paths = [tmp_path / "mapping.json", tmp_path / "ownership.json"]
    def first(argv, *, verification_session):
        batch = state.FreshDigestBatch()
        batch.read([source])
        verification_session.batches["Species_one"] = batch
        output_paths[0].write_text('{"species_passed": 1}')
        return 0
    def second(argv, *, verification_session):
        info = source.stat()
        source.write_text("modified")
        os.utime(source, ns=(info.st_atime_ns, info.st_mtime_ns + 2_000_000_000))
        output_paths[1].write_text('{"species_passed": 1}')
        return 0
    monkeypatch.setattr(mapping, "main", first)
    monkeypatch.setattr(ownership, "main", second)
    argv = ["final-qc"]
    for key in ("species-cds-dir", "species-gff-dir", "species-genome-dir", "species-summary"):
        argv.extend(["--" + key, str(tmp_path)])
    for key, path in zip(("mapping-stats-output", "ownership-stats-output"), output_paths, strict=True):
        argv.extend(["--" + key, str(path)])
    monkeypatch.setattr(sys, "argv", argv)
    assert reuse.main() == 1
    assert all(not path.exists() for path in output_paths)

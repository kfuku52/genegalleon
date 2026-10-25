import json
import os
import stat
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

from workflow.support import artifact_provenance as PROVENANCE
from workflow.support import orthofinder_output_store as STORE
from workflow.support import workspace_output_storage as WORKSPACE
from workflow.support.shared_namespace_lock import namespace_lock


def result(tmp_path):
    root = tmp_path / "workspace/output/orthofinder"
    root.mkdir(parents=True)
    (root / "Log.txt").write_text("2026-08-03: OrthoFinder run completed\n")
    for name in STORE.MANAGED_DIRECTORIES:
        raw = root / name
        (raw / "empty").mkdir(parents=True)
        for og in ("OG0000001", "OG9999999"):
            (raw / (og + ".fa")).write_text(">Species_gene\nACGT--ACGT\n")
        os.chmod(raw / "OG9999999.fa", 0o640)
        os.utime(raw / "OG9999999.fa", ns=(1700000000123456789, 1700000000123456789))
    tables = root / "Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv"
    tables.parent.mkdir()
    tables.write_text("Orthogroup\tTotal\nOG0000001\t1\n")
    return root


@pytest.mark.parametrize("workers", [1, 2, 4])
def test_lossless_all_og_roundtrip_and_logical_provenance(tmp_path, workers):
    root = result(tmp_path)
    with STORE.storage_lock(root):
        pass
    before = PROVENANCE.sha256_path(root)
    core_before = PROVENANCE.sha256_path(root / "core")
    output_before = PROVENANCE.sha256_path(root.parent)
    original = {
        str(p.relative_to(root)): (p.read_bytes(), p.stat().st_mode & 0o777, p.stat().st_mtime_ns)
        for p in root.rglob("*")
        if p.is_file()
    }
    STORE.convert(root, workers=workers)
    assert PROVENANCE.sha256_path(root) == before
    assert PROVENANCE.sha256_path(root / "core") == core_before
    assert PROVENANCE.sha256_path(root.parent) == output_before
    for name in STORE.MANAGED_DIRECTORIES:
        assert not (root / name).exists()
        with STORE.open_binary(root, name + "/OG9999999.fa") as handle:
            assert handle.read() == b">Species_gene\nACGT--ACGT\n"
        data = STORE.archive_records(root, name, deep=True)
        assert len(data["files"]) == 2
    with STORE.materialized_directory(root, "core/Orthogroup_Sequences", tmp_path) as native:
        assert (native / "OG9999999.fa").read_bytes() == b">Species_gene\nACGT--ACGT\n"
        assert (native / "empty").is_dir()
    assert not native.exists()
    STORE.convert(root, to="raw", workers=workers)
    assert PROVENANCE.sha256_path(root) == before
    for relative, metadata in original.items():
        p = root / relative
        assert (p.read_bytes(), p.stat().st_mode & 0o777, p.stat().st_mtime_ns) == metadata
    assert not list(root.rglob("*.zip"))


def test_workspace_staged_conversion_dry_run_and_deep_audit(tmp_path):
    root = result(tmp_path)
    base = [
        "convert",
        "--workspace",
        str(root.parents[1]),
        "--target",
        "orthofinder",
        "--to",
        "zip",
        "--orthofinder-directory",
        "Orthogroup_Sequences",
        "--workers",
        "2",
    ]
    code, report, _, _ = WORKSPACE.run(WORKSPACE.build_parser().parse_args(base + ["--dry-run"]))
    assert code == 0 and report["targets"][0]["status"] == "dry-run"
    assert not list(root.rglob("*.zip"))
    code, report, _, _ = WORKSPACE.run(WORKSPACE.build_parser().parse_args(base))
    assert code == 0 and report["totals"]["file_reduction"] == 1
    assert (root / "MultipleSequenceAlignments/OG9999999.fa").is_file()
    code, report, _, _ = WORKSPACE.run(
        WORKSPACE.build_parser().parse_args(
            ["audit", "--workspace", str(root.parents[1]), "--target", "orthofinder", "--verification", "deep"]
        )
    )
    assert code == 0 and report["targets"][0]["verification"]["zip_files"] == 1


def test_interrupted_deletion_resumes_only_verified_remaining_members(tmp_path, monkeypatch):
    root = result(tmp_path)
    before = PROVENANCE.sha256_path(root)
    unlink = Path.unlink

    def interrupt(path, *args, **kwargs):
        if path.name == "OG9999999.fa":
            raise OSError("interrupted deletion")
        return unlink(path, *args, **kwargs)

    with monkeypatch.context() as scoped:
        scoped.setattr(Path, "unlink", interrupt)
        with pytest.raises(OSError, match="interrupted"):
            STORE.convert(root, directories=["Orthogroup_Sequences"])
    assert (root / "Orthogroup_Sequences.zip").is_file()
    assert not (root / "Orthogroup_Sequences/OG0000001.fa").exists()
    with pytest.raises(STORE.OrthoFinderStorageError, match="Incomplete conversion"):
        STORE.logical_records(root)
    STORE.convert(root)
    assert PROVENANCE.sha256_path(root) == before


def test_changed_source_after_zip_publication_preserves_raw_and_zip(tmp_path, monkeypatch):
    root = result(tmp_path)
    retire = STORE._retire

    def changed(raw, data, **kwargs):
        (raw / "OG9999999.fa").write_bytes(b"changed")
        return retire(raw, data, **kwargs)

    with monkeypatch.context() as scoped:
        scoped.setattr(STORE, "_retire", changed)
        with pytest.raises(STORE.OrthoFinderStorageError, match="differs"):
            STORE.convert(root, directories=["Orthogroup_Sequences"])
    assert (root / "Orthogroup_Sequences/OG0000001.fa").exists()
    with pytest.raises(STORE.OrthoFinderStorageError, match="differs"):
        STORE.convert(root, directories=["Orthogroup_Sequences"])


@pytest.mark.parametrize("problem", ["incomplete", "symlink", "unmarked", "extra-directory"])
def test_unsafe_sources_and_conflicts_never_delete_input(tmp_path, problem):
    root = result(tmp_path)
    source = root / "Orthogroup_Sequences/OG9999999.fa"
    if problem == "incomplete":
        (root / "Log.txt").write_text("OrthoFinder running\n")
    elif problem == "symlink":
        source.unlink()
        source.symlink_to(root / "Log.txt")
    elif problem == "extra-directory":
        STORE.convert(root, directories=["Orthogroup_Sequences"])
        source.parent.mkdir()
        (source.parent / "new-empty-directory").mkdir()
    else:
        with zipfile.ZipFile(root / "Orthogroup_Sequences.zip", "w") as archive:
            archive.writestr("foo", "bar")
    with pytest.raises(STORE.OrthoFinderStorageError):
        STORE.convert(root, directories=["Orthogroup_Sequences"])
    assert source.exists() or source.is_symlink() or (source.parent / "new-empty-directory").exists()


@pytest.mark.parametrize("problem", ["sha", "traversal", "duplicate"])
def test_zip_payload_and_manifest_validation(tmp_path, problem):
    root = result(tmp_path)
    STORE.convert(root, directories=["Orthogroup_Sequences"])
    path = root / "Orthogroup_Sequences.zip"
    with zipfile.ZipFile(path) as archive:
        values = {info.filename: archive.read(info) for info in archive.infolist()}
    data = json.loads(values[STORE.MANIFEST])
    if problem == "sha":
        data["files"][0]["sha256"] = "0" * 64
    elif problem == "traversal":
        data["files"][0]["path"] = "../escape"
    else:
        data["files"].append(data["files"][0])
    values[STORE.MANIFEST] = json.dumps(data).encode()
    with zipfile.ZipFile(path, "w") as archive:
        for name, value in values.items():
            archive.writestr(name, value)
        archive.comment = STORE.MARKER + b"Orthogroup_Sequences"
    with pytest.raises(STORE.OrthoFinderStorageError):
        STORE.archive_records(root, "Orthogroup_Sequences", deep=True)
    with pytest.raises(STORE.OrthoFinderStorageError):
        STORE.convert(root, to="raw", directories=["Orthogroup_Sequences"])
    assert path.exists() and not (root / "Orthogroup_Sequences").exists()


def test_lock_is_project_local_and_blocks_a_second_writer(tmp_path):
    root = result(tmp_path)
    other = result(tmp_path / "other")
    with STORE.storage_lock(root):
        with namespace_lock(STORE.lock_path(root), exclusive=True, nonblocking=True) as acquired:
            assert not acquired
        STORE.convert(other, directories=["Orthogroup_Sequences"])


def test_quota_preflight_includes_concurrent_directories(tmp_path):
    root = result(tmp_path)
    args = WORKSPACE.build_parser().parse_args(
        [
            "convert",
            "--workspace",
            str(root.parents[1]),
            "--target",
            "orthofinder",
            "--to",
            "zip",
            "--workers",
            "4",
            "--available-bytes",
            "2000000",
            "--dry-run",
        ]
    )
    code, report, _, _ = WORKSPACE.run(args)
    assert code == 1 and not report["quota_plan"]["sufficient"]
    assert not list(root.rglob("*.zip"))


def test_existing_directory_provenance_remains_current_after_conversion(tmp_path):
    root = result(tmp_path)
    workspace = root.parents[1]
    output = workspace / "summary.tsv"
    output.write_text("unchanged summary\n")
    args = [
        "--manifest",
        str(workspace / "summary.json"),
        "--step",
        "csubst_site",
        "--family-id",
        "orthogroup",
        "--logical-root",
        str(workspace / "output/.gg_global_artifacts"),
        "--workspace-root",
        str(workspace),
        "--input",
        f"orthofinder={root}",
        "--output",
        f"summary={output}",
    ]
    cli = Path(PROVENANCE.__file__)
    record = subprocess.run([sys.executable, str(cli), "record", *args], capture_output=True, text=True)
    assert record.returncode == 0, record.stderr
    manifest = (workspace / "summary.json").read_bytes()
    STORE.convert(root)
    current = subprocess.run([sys.executable, str(cli), "needs-run", *args], capture_output=True, text=True)
    assert current.returncode == 1, current.stdout + current.stderr
    assert (workspace / "summary.json").read_bytes() == manifest


def test_zip64_member_count_above_classic_zip_limit(tmp_path):
    root = tmp_path / "orthofinder"
    raw = root / "Orthogroup_Sequences"
    raw.mkdir(parents=True)
    (root / "Log.txt").write_text("OrthoFinder run completed\n")
    for index in range(65536):
        (raw / f"OG{index:07d}.fa").touch()
    STORE.convert(root, directories=["Orthogroup_Sequences"])
    data = STORE.archive_records(root, "Orthogroup_Sequences", deep=True)
    assert len(data["files"]) == 65536 and not raw.exists()
    with zipfile.ZipFile(root / "Orthogroup_Sequences.zip") as archive:
        assert len(archive.infolist()) == 65537


def test_logical_digest_keeps_path_component_order_for_prefix_collisions(tmp_path):
    root = result(tmp_path)
    (root / "Orthogroups").mkdir()
    (root / "Orthogroups/table.tsv").write_text("table\n")
    (root / "Orthogroups.tsv").write_text("sibling with shared prefix\n")
    before = PROVENANCE.sha256_path(root)
    STORE.convert(root)
    assert PROVENANCE.sha256_path(root) == before


def test_orthofinder_conversion_does_not_require_an_unrelated_selected_catalog(tmp_path):
    root = result(tmp_path)
    (root / "Orthogroups.GeneCount.selected.tsv").write_text("a different selected catalog\n")
    args = WORKSPACE.build_parser().parse_args([
        "convert", "--workspace", str(root.parents[1]), "--target", "orthofinder", "--to", "zip", "--dry-run",
    ])
    code, report, _, _ = WORKSPACE.run(args)
    assert code == 0 and report["targets"][0]["status"] == "dry-run"


@pytest.mark.parametrize("to", ["zip", "raw"])
def test_directory_sync_failure_preserves_recovery_data(tmp_path, monkeypatch, to):
    root = result(tmp_path)
    raw = root / "Orthogroup_Sequences"
    archive = raw.with_suffix(".zip")
    if to == "raw":
        STORE.convert(root, directories=[raw.name])
    fsync = os.fsync

    def fail_directory_sync(fd):
        if stat.S_ISDIR(os.fstat(fd).st_mode):
            raise OSError("directory sync failed")
        return fsync(fd)

    monkeypatch.setattr(os, "fsync", fail_directory_sync)
    with pytest.raises(OSError, match="directory sync failed"):
        STORE.convert(root, to=to, directories=[raw.name])
    assert archive.is_file() and (raw / "OG9999999.fa").is_file()


def test_restore_syncs_nested_directory_entries_before_removing_zip(tmp_path, monkeypatch):
    root = result(tmp_path)
    raw = root / "Orthogroup_Sequences"
    (raw / "nested").mkdir()
    (raw / "nested/member.fa").write_text(">gene\nACGT\n")
    STORE.convert(root, directories=[raw.name])
    synced = set()
    fsync, unlink = os.fsync, Path.unlink

    def record_sync(fd):
        if stat.S_ISDIR(os.fstat(fd).st_mode):
            synced.add(os.fstat(fd).st_ino)
        return fsync(fd)

    def check_commit(path, *args, **kwargs):
        if path == raw.with_suffix(".zip"):
            assert {p.stat().st_ino for p in [raw, raw / "nested", raw / "empty", raw.parent]} <= synced
        return unlink(path, *args, **kwargs)

    monkeypatch.setattr(os, "fsync", record_sync)
    monkeypatch.setattr(Path, "unlink", check_commit)
    STORE.convert(root, to="raw", directories=[raw.name])


def test_changed_directory_mode_is_not_discarded_by_retirement(tmp_path, monkeypatch):
    root = result(tmp_path)
    retire = STORE._retire

    def changed(raw, data, **kwargs):
        (raw / "empty").chmod(0o750)
        return retire(raw, data, **kwargs)

    monkeypatch.setattr(STORE, "_retire", changed)
    with pytest.raises(STORE.OrthoFinderStorageError, match="metadata"):
        STORE.convert(root, directories=["Orthogroup_Sequences"])
    assert (root / "Orthogroup_Sequences/empty").stat().st_mode & 0o777 == 0o750
    assert (root / "Orthogroup_Sequences/OG9999999.fa").is_file()


@pytest.mark.parametrize("name", STORE.MANAGED_DIRECTORIES)
def test_managed_logical_directory_digest_accepts_versioned_manifest(tmp_path, name):
    root = result(tmp_path)
    before = PROVENANCE.raw_or_zip_directory_digest(root / name)
    STORE.convert(root, directories=[name])
    assert PROVENANCE.raw_or_zip_directory_digest(root / name) == before


@pytest.mark.parametrize("bad_data", [[], {"schema_version": 1, "directory": "Orthogroup_Sequences", "files": None}])
def test_invalid_manifest_types_fail_with_storage_error(tmp_path, bad_data):
    root = result(tmp_path)
    path = root / "Orthogroup_Sequences.zip"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(STORE.MANIFEST, json.dumps(bad_data))
        archive.comment = STORE.MARKER + b"Orthogroup_Sequences"
    with pytest.raises(STORE.OrthoFinderStorageError, match="Invalid"):
        STORE.archive_records(root, "Orthogroup_Sequences")
    assert (root / "Orthogroup_Sequences/OG9999999.fa").is_file()


def test_parallel_backend_checks_combined_temporary_space(tmp_path):
    root = result(tmp_path)
    with pytest.raises(STORE.OrthoFinderStorageError, match="Insufficient temporary space"):
        STORE.convert(root, workers=4, available_bytes=2_000_000)
    assert not list(root.rglob("*.zip"))


def test_empty_raw_directory_still_requires_completed_run_in_dry_run(tmp_path):
    root = tmp_path / "workspace/output/orthofinder"
    (root / "Orthogroup_Sequences").mkdir(parents=True)
    args = WORKSPACE.build_parser().parse_args([
        "convert", "--workspace", str(root.parents[1]), "--target", "orthofinder", "--to", "zip", "--dry-run",
    ])
    code, report, _, _ = WORKSPACE.run(args)
    assert code == 1 and "completed" in report["targets"][0]["error"]


def test_native_child_can_convert_without_releasing_its_parent_fence(tmp_path):
    root = result(tmp_path)
    with STORE.storage_lock(root):
        owner = STORE.inspect_lock(STORE.lock_path(root))["exclusive"]
        command = [sys.executable, STORE.__file__, "pack", "--root", str(root),
                   "--directory", "Orthogroup_Sequences", "--producer-lock-token", owner["token"]]
        child = subprocess.run(command, capture_output=True, text=True, timeout=15)
        assert child.returncode == 0, child.stderr
        assert STORE.inspect_lock(STORE.lock_path(root))["exclusive"] == owner
        command[-1] = "f" * 32
        child = subprocess.run(command, capture_output=True, text=True, timeout=15)
        assert child.returncode != 0 and "ownership" in child.stderr
        assert STORE.inspect_lock(STORE.lock_path(root))["exclusive"] == owner
        with pytest.raises(STORE.OrthoFinderStorageError, match="ownership"):
            STORE.convert(root, producer_lock_token=owner["token"])


def test_interrupted_restore_resumes_and_keeps_zip_until_full_verification(tmp_path, monkeypatch):
    root = result(tmp_path)
    raw = root / "Orthogroup_Sequences"
    STORE.convert(root, directories=[raw.name])
    link = os.link

    def interrupt(source, target):
        if Path(target).name == "OG9999999.fa":
            raise OSError("restore interrupted")
        return link(source, target)

    with monkeypatch.context() as scoped:
        scoped.setattr(os, "link", interrupt)
        with pytest.raises(OSError, match="restore interrupted"):
            STORE.convert(root, to="raw", directories=[raw.name])
    assert (raw / "OG0000001.fa").is_file() and raw.with_suffix(".zip").is_file()
    STORE.convert(root, to="raw", directories=[raw.name])
    assert (raw / "OG9999999.fa").is_file() and not raw.with_suffix(".zip").exists()


def test_recovery_zip_replacement_during_retirement_never_deletes_source(tmp_path, monkeypatch):
    root = result(tmp_path)
    retire = STORE._retire

    def replaced(raw, data, **kwargs):
        path = raw.with_suffix(".zip")
        replacement = path.with_suffix(".replacement")
        replacement.write_bytes(path.read_bytes())
        os.replace(replacement, path)
        return retire(raw, data, **kwargs)

    monkeypatch.setattr(STORE, "_retire", replaced)
    with pytest.raises(STORE.OrthoFinderStorageError, match="Recovery archive changed"):
        STORE.convert(root, directories=["Orthogroup_Sequences"])
    assert (root / "Orthogroup_Sequences/OG0000001.fa").is_file()


def test_raw_change_during_restore_commit_preserves_recovery_zip(tmp_path, monkeypatch):
    root = result(tmp_path)
    raw = root / "Orthogroup_Sequences"
    STORE.convert(root, directories=[raw.name])
    sync = STORE._sync_directory

    def changed(directory):
        if directory == raw:
            (raw / "OG9999999.fa").write_bytes(b"concurrent writer")
        return sync(directory)

    monkeypatch.setattr(STORE, "_sync_directory", changed)
    with pytest.raises(STORE.OrthoFinderStorageError, match="inventory changed"):
        STORE.convert(root, to="raw", directories=[raw.name])
    assert raw.with_suffix(".zip").is_file()


def test_audit_reports_leftover_restore_files(tmp_path):
    root = result(tmp_path)
    partial = root / ".Orthogroup_Sequences.restore.interrupted"
    partial.write_bytes(b"partial restored member")
    args = WORKSPACE.build_parser().parse_args([
        "audit", "--workspace", str(root.parents[1]), "--target", "orthofinder",
    ])
    code, report, _, _ = WORKSPACE.run(args)
    assert code == 1 and str(partial) in report["targets"][0]["preflight"]["partial_paths"]


def test_resumed_retirement_requires_a_durable_recovery_zip(tmp_path, monkeypatch):
    root = result(tmp_path)
    raw = root / "Orthogroup_Sequences"
    with monkeypatch.context() as scoped:
        def interrupt(*args, **kwargs):
            raise OSError("retirement interrupted")
        scoped.setattr(STORE, "_retire", interrupt)
        with pytest.raises(OSError, match="retirement interrupted"):
            STORE.convert(root, directories=[raw.name])

    def fail_sync(path):
        raise OSError("recovery publication not synced")

    monkeypatch.setattr(STORE, "_sync_directory", fail_sync)
    with pytest.raises(OSError, match="not synced"):
        STORE.convert(root, directories=[raw.name])
    assert (raw / "OG9999999.fa").is_file() and raw.with_suffix(".zip").is_file()


@pytest.mark.parametrize("member", ["OG0000001.fa", "OG9999999.fa"])
def test_member_reader_rejects_mixed_raw_and_zip_sources(tmp_path, member):
    root = result(tmp_path)
    STORE.convert(root, directories=["Orthogroup_Sequences"])
    raw = root / "Orthogroup_Sequences"
    raw.mkdir()
    (raw / "OG9999999.fa").write_bytes(b"unverified raw member")
    with pytest.raises(STORE.OrthoFinderStorageError, match="Incomplete conversion"):
        with STORE.open_binary(root, raw.name + "/" + member):
            pytest.fail("mixed member was exposed")


def test_unmarked_zip_named_directory_remains_an_ordinary_provenance_tree(tmp_path):
    root = result(tmp_path)
    unrelated = root / "unrelated/Orthogroup_Sequences.zip"
    unrelated.mkdir(parents=True)
    (unrelated / "ordinary.txt").write_bytes(b"ordinary directory contents")
    before = PROVENANCE.sha256_path(root)
    STORE.convert(root)
    assert PROVENANCE.sha256_path(root) == before


def test_cached_logical_zip_digest_cannot_hide_a_new_raw_zip_conflict(tmp_path):
    root = result(tmp_path)
    raw = root / "Orthogroup_Sequences"
    contents = {p.name: p.read_bytes() for p in raw.glob("*.fa")}
    PROVENANCE.configure_digest_cache(tmp_path / "digests.sqlite")
    try:
        before = PROVENANCE.raw_or_zip_directory_digest(raw)
        STORE.convert(root, directories=[raw.name])
        assert PROVENANCE.raw_or_zip_directory_digest(raw) == before  # Populate the archive cache.
        raw.mkdir()
        for name, value in contents.items():
            (raw / name).write_bytes(value)
        with pytest.raises(PROVENANCE.ProvenanceError, match="Incomplete conversion"):
            PROVENANCE.raw_or_zip_directory_digest(raw)
    finally:
        PROVENANCE.configure_digest_cache(None)

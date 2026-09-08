"""Storage lifecycle: conversion, deletion history, and generation recovery."""
import os

import pytest

from workflow.support import gene_family_output_store as store

FAMILY = "OG0000001"


def archived(root):
    path = root / "mafft" / f"{FAMILY}_alignment.fa"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"original")
    store.archive_completed_outputs(root, "orthogroup", [FAMILY], store.orthogroup_id_from_name,
                                    include_incomplete=True)
    inventory = store.family_inventory_path(root, FAMILY)
    inventory.mkdir(parents=True)
    (inventory / "worker.paths").write_bytes(os.fsencode(path.relative_to(root)) + b"\0")
    store.enqueue_family_archive(root, "orthogroup", FAMILY)
    return path


def drain(root):
    return store.drain_archive_queue(root, "orthogroup", store.orthogroup_id_from_name)


@pytest.mark.parametrize("pure_raw", [False, True])
def test_raw_conversion_cancels_pending_zip_requests(tmp_path, pure_raw):
    path = archived(tmp_path)
    store.convert_storage_to_raw(tmp_path, "orthogroup", pure_raw=pure_raw)
    assert store.archive_queue_status(tmp_path)["pending_families"] == 0
    assert drain(tmp_path)["archived_files"] == 0
    assert path.read_bytes() == b"original"


@pytest.mark.parametrize("target", ["raw", "zip"])
def test_collector_defers_interrupted_conversion(tmp_path, target):
    path = archived(tmp_path)
    store._write_storage_conversion_marker(tmp_path, "orthogroup", target)
    path.parent.mkdir(exist_ok=True)
    path.write_bytes(b"new result")
    assert drain(tmp_path)["status"] == "conversion-pending"
    assert path.read_bytes() == b"new result"


def test_conversion_cannot_begin_during_collection(tmp_path):
    archived(tmp_path)
    with store._bucket_lock(store._archive_state_root(tmp_path) / "collector.lock", exclusive=True):
        with pytest.raises(store.ArchiveStoreError, match="collector"):
            with store.storage_conversion_session(tmp_path, "orthogroup", "raw"):
                pytest.fail("Conversion overlapped active collector")
    assert store._read_storage_conversion_marker(tmp_path) is None


def test_missing_counter_never_reuses_a_deletion_generation(tmp_path):
    path = archived(tmp_path)
    logical = store.GeneFamilyOutputStore(tmp_path)
    logical.delete(f"mafft/{path.name}")
    (logical.archive_root / store.GENERATION_FILE).unlink()
    path.parent.mkdir(exist_ok=True)
    path.write_bytes(b"new result")
    assert drain(tmp_path)["archived_files"] == 1
    with store.GeneFamilyOutputStore(tmp_path).open_binary("mafft", path.name) as handle:
        assert handle.read() == b"new result"


def test_repair_advances_stale_counter_past_physical_archive_generation(tmp_path):
    path = archived(tmp_path)
    logical = store.GeneFamilyOutputStore(tmp_path)
    artifact = logical.artifact("mafft", path.name)
    store._compact_artifact_chunk(tmp_path / "archives", "mafft", [artifact], "orthogroup", 100)
    store.repair_archive_index(tmp_path)
    recovered = store.GeneFamilyOutputStore(tmp_path)
    with store.archive_lock(recovered.archive_root):
        assert recovered._next_generation() > 100


@pytest.mark.parametrize("operation", ["restore", "undelete-purge"])
@pytest.mark.parametrize("physical_current", [False, True])
def test_restore_through_current_name_overrides_older_alias_deletion(tmp_path, operation, physical_current):
    source = (tmp_path / "stat_branch" / f"{FAMILY}_stat.branch.tsv" if physical_current
              else tmp_path / "stat.branch" / f"{FAMILY}.stat.branch.tsv")
    source.parent.mkdir()
    source.write_bytes(b"historical result")
    store.archive_completed_outputs(tmp_path, "orthogroup", [FAMILY], store.orthogroup_id_from_name,
                                    include_incomplete=True)
    logical = store.GeneFamilyOutputStore(tmp_path)
    logical.delete(f"stat.branch/{FAMILY}.stat.branch.tsv")
    current = f"stat_branch/{FAMILY}_stat.branch.tsv"
    if operation == "restore":
        assert logical.restore(current).read_bytes() == b"historical result"
    else:
        logical.undelete(current)
        store.purge_archives(tmp_path, "orthogroup")
        with store.GeneFamilyOutputStore(tmp_path).open_binary(*current.split("/")) as handle:
            assert handle.read() == b"historical result"


def test_missing_counter_uses_compaction_generation(tmp_path):
    path = archived(tmp_path)
    logical = store.GeneFamilyOutputStore(tmp_path)
    store._compact_artifact_chunk(tmp_path / "archives", "mafft", [logical.artifact("mafft", path.name)],
                                  "orthogroup", 100)
    (logical.archive_root / store.GENERATION_FILE).unlink()
    with store.archive_lock(logical.archive_root):
        assert logical._next_generation() > 100



def test_enqueue_rejects_pending_conversion(tmp_path):
    archived(tmp_path)
    request = store._queue_path(tmp_path, FAMILY)
    request.unlink()
    store._write_storage_conversion_marker(tmp_path, "orthogroup", "raw")
    with pytest.raises(store.ArchiveStoreError, match="conversion"):
        store.enqueue_family_archive(tmp_path, "orthogroup", FAMILY)
    assert not request.exists()


def test_conversion_excludes_collection_between_phases(tmp_path):
    path = archived(tmp_path)
    seen = []

    def between_phases(**fields):
        if fields.get("phase") == "removing-archives":
            seen.append(drain(tmp_path)["status"])

    store.convert_storage_to_raw(tmp_path, "orthogroup", progress_callback=between_phases)
    assert seen == ["collector-busy"]
    assert path.read_bytes() == b"original"


def test_restore_failure_keeps_alias_deletion(tmp_path, monkeypatch):
    source = tmp_path / "stat.branch" / f"{FAMILY}.stat.branch.tsv"
    source.parent.mkdir()
    source.write_bytes(b"historical result")
    store.archive_completed_outputs(tmp_path, "orthogroup", [FAMILY], store.orthogroup_id_from_name,
                                    include_incomplete=True)
    logical = store.GeneFamilyOutputStore(tmp_path)
    logical.delete(f"stat.branch/{source.name}")

    def no_space(*args, **kwargs):
        raise OSError("No space left on device")

    monkeypatch.setattr(logical, "materialize", no_space)
    current = f"stat_branch/{FAMILY}_stat.branch.tsv"
    with pytest.raises(OSError, match="No space"):
        logical.restore(current)
    assert not store.GeneFamilyOutputStore(tmp_path).logical_exists(current)



def test_raw_conversion_resumes_after_failure_without_rearchiving(tmp_path):
    path = archived(tmp_path)

    def interrupted(**fields):
        if fields.get("phase") == "removing-archives":
            raise OSError("simulated interruption")

    with pytest.raises(OSError, match="simulated interruption"):
        store.convert_storage_to_raw(tmp_path, "orthogroup", progress_callback=interrupted)
    assert drain(tmp_path)["status"] == "conversion-pending"
    result = store.convert_storage_to_raw(tmp_path, "orthogroup", require_resume=True)
    assert result["conversion_resumed"]
    assert store.archive_queue_status(tmp_path)["pending_families"] == 0
    assert path.read_bytes() == b"original"


def test_repair_preserves_higher_reserved_generation(tmp_path):
    archived(tmp_path)
    counter = store._archive_state_root(tmp_path) / store.GENERATION_FILE
    counter.write_text("1000\n")
    store.repair_archive_index(tmp_path)
    logical = store.GeneFamilyOutputStore(tmp_path)
    with store.archive_lock(logical.archive_root):
        assert logical._next_generation() > 1000



@pytest.mark.parametrize("corrupt_counter", ["", "0", "-1", "broken"])
def test_invalid_counter_blocks_writes_until_explicit_repair(tmp_path, corrupt_counter):
    path = archived(tmp_path)
    counter = store._archive_state_root(tmp_path) / store.GENERATION_FILE
    counter.write_text(corrupt_counter)
    path.parent.mkdir(exist_ok=True)
    path.write_bytes(b"new result")
    with pytest.raises(store.ArchiveStoreError, match="generation counter"):
        drain(tmp_path)
    assert path.read_bytes() == b"new result"
    store.repair_archive_index(tmp_path)
    assert drain(tmp_path)["archived_files"] == 1
    with store.GeneFamilyOutputStore(tmp_path).open_binary("mafft", path.name) as handle:
        assert handle.read() == b"new result"

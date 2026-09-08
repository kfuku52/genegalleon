"""Regression coverage for ZIP encoding, destinations, and materialization."""
import os
import zipfile
from dataclasses import replace

import pytest

from workflow.support import gene_family_output_store as store

CONTENT = b"ACTGACTG\tgene-family-result\n" * 4000


def raw_archive(root, *, level=6, mtime_ns=None, destination=None):
    source = root / "mafft" / "OG0000001_alignment.fa"
    source.parent.mkdir(parents=True, exist_ok=True)
    source.write_bytes(CONTENT)
    if mtime_ns is not None:
        os.utime(source, ns=(mtime_ns, mtime_ns))
    return store._archive_chunk(
        root, root / "archives", "mafft", [source], "orthogroup", 1,
        store.orthogroup_id_from_name, compression="deflate", compression_level=level,
        destination_path=destination,
    )


@pytest.mark.parametrize("compact", [False, True])
def test_compression_level_controls_streamed_member_bytes(tmp_path, compact):
    sizes = []
    for level in (0, 9):
        root = tmp_path / str(level)
        archive_path, artifacts, _ = raw_archive(root, level=level)
        if compact:
            archive_path, _ = store._compact_artifact_chunk(
                root / "archives", "mafft", artifacts, "orthogroup", 2,
                compression="deflate", compression_level=level,
            )
        with zipfile.ZipFile(archive_path) as archive:
            assert archive.read(artifacts[0].member_name) == CONTENT
            sizes.append(archive.getinfo(artifacts[0].member_name).compress_size)
    assert sizes[0] >= len(CONTENT)  # DEFLATE level zero stores uncompressed blocks.
    assert sizes[1] < len(CONTENT) // 10


@pytest.mark.parametrize("mtime_ns,year", [(0, 1980), (7258118400_000000000, 2107)])
@pytest.mark.parametrize("compact", [False, True])
def test_zip_timestamp_limits_preserve_precise_manifest_mtime(tmp_path, mtime_ns, year, compact):
    archive_path, artifacts, _ = raw_archive(tmp_path, mtime_ns=None if compact else mtime_ns)
    if compact:
        artifacts = [replace(artifacts[0], mtime_ns=mtime_ns)]
        archive_path, artifacts = store._compact_artifact_chunk(
            tmp_path / "archives", "mafft", artifacts, "orthogroup", 2,
        )
    with zipfile.ZipFile(archive_path) as archive:
        assert archive.read(artifacts[0].member_name) == CONTENT
        assert archive.getinfo(artifacts[0].member_name).date_time[0] == year
    assert artifacts[0].mtime_ns == mtime_ns
    store.repair_archive_index(tmp_path)
    source = tmp_path / artifacts[0].logical_path
    source.unlink()
    restored = store.GeneFamilyOutputStore(tmp_path).materialize("mafft", source.name)
    assert restored.read_bytes() == CONTENT
    assert restored.stat().st_mtime_ns == mtime_ns


@pytest.mark.parametrize("compact", [False, True])
@pytest.mark.parametrize("link_parent", [False, True, "ancestor"])
def test_zip_writers_reject_symlink_destinations(tmp_path, compact, link_parent):
    root = tmp_path / "output"
    _, artifacts, _ = raw_archive(root)
    outside = tmp_path / "outside"
    outside.mkdir()
    if link_parent:
        alias = root / "alias"
        alias.symlink_to(outside, target_is_directory=True)
        if link_parent == "ancestor":
            (outside / "nested").mkdir()
            destination = alias / "nested" / "mafft.zip"
        else:
            destination = alias / "mafft.zip"
    else:
        destination = root / "mafft.zip"
        destination.symlink_to(outside / "mafft.zip")
    with pytest.raises(store.ArchiveStoreError, match="Symlinked"):
        if compact:
            store._compact_artifact_chunk(root / "archives", "mafft", artifacts, "orthogroup", 2,
                                          destination_path=destination)
        else:
            raw_archive(root, destination=destination)
    assert not list(outside.rglob("*.zip"))


def test_materialize_validates_paths_before_existing_file_shortcut(tmp_path):
    root = tmp_path / "output"
    (root / "mafft").mkdir(parents=True)
    outside = tmp_path / "outside.txt"
    outside.write_bytes(b"keep")
    with pytest.raises(store.ArchiveStoreError, match="Unsafe logical"):
        store.GeneFamilyOutputStore(root).materialize("mafft", "../../outside.txt")
    assert outside.read_bytes() == b"keep"


def test_materialization_publishes_before_releasing_reader_lock(tmp_path, monkeypatch):
    import contextlib

    _, artifacts, _ = raw_archive(tmp_path)
    store.repair_archive_index(tmp_path)
    destination = tmp_path / artifacts[0].logical_path
    destination.unlink()
    original_lock = store.producer_read_lock

    @contextlib.contextmanager
    def writer_after_reader_unlock(archive_root, **kwargs):
        with original_lock(archive_root, **kwargs) as acquired:
            yield acquired
        # A maintenance writer can take over as soon as materialization releases
        # the reader lock; its newer output must not be overwritten afterwards.
        with store.producer_quiescence_lock(archive_root) as acquired:
            assert acquired
            destination.write_bytes(b"new generation")

    monkeypatch.setattr(store, "producer_read_lock", writer_after_reader_unlock)
    store.GeneFamilyOutputStore(tmp_path).materialize("mafft", destination.name)
    assert destination.read_bytes() == b"new generation"


@pytest.mark.parametrize("damage", ["corrupt", "missing"])
@pytest.mark.parametrize("lost_copy", [None, "index", "subdir_index"])
def test_repair_does_not_erase_index_of_corrupted_final_zip(tmp_path, damage, lost_copy):
    final_zip = tmp_path / "mafft.zip"
    raw_archive(tmp_path, destination=final_zip)
    store.repair_archive_index(tmp_path)
    if lost_copy is not None:
        directory = tmp_path / ".gg_store" / (store.INDEX_DIR_NAME if lost_copy == "index" else store.SUBDIR_INDEX_DIR_NAME)
        for path in directory.glob("*.json"):
            path.unlink()
    index_dir = tmp_path / ".gg_store"
    before = {path.relative_to(index_dir): path.read_bytes() for path in index_dir.rglob("*.json")}
    if damage == "corrupt":
        final_zip.write_bytes(b"broken central directory")
    else:
        final_zip.unlink()
    with pytest.raises(store.ArchiveStoreError):
        store.repair_archive_index(tmp_path)
    assert {path.relative_to(index_dir): path.read_bytes() for path in index_dir.rglob("*.json")} == before



def test_repair_rebuilds_damaged_json_when_final_zip_is_healthy(tmp_path):
    raw_archive(tmp_path, destination=tmp_path / "mafft.zip")
    store.repair_archive_index(tmp_path)
    for name in (store.INDEX_DIR_NAME, store.SUBDIR_INDEX_DIR_NAME):
        for path in (tmp_path / ".gg_store" / name).glob("*.json"):
            path.write_text("{incomplete JSON")
    store.repair_archive_index(tmp_path)
    assert store.GeneFamilyOutputStore(tmp_path).verify()


def test_live_directory_symlink_cannot_override_archived_content(tmp_path):
    _, artifacts, _ = raw_archive(tmp_path)
    store.repair_archive_index(tmp_path)
    source = tmp_path / artifacts[0].logical_path
    outside = tmp_path / "unrelated"
    source.parent.rename(outside)
    source.parent.symlink_to(outside, target_is_directory=True)
    (outside / source.name).write_bytes(b"unrelated external result")
    logical = store.GeneFamilyOutputStore(tmp_path)
    with logical.open_binary("mafft", source.name) as handle:
        assert handle.read() == CONTENT

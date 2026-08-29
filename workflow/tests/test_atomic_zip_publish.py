from __future__ import annotations

import os
import zipfile
from pathlib import Path

import pytest

from workflow.support import atomic_zip_publish as PUBLISH


def _write_zip(path: Path, payload: bytes = b"payload\n") -> None:
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("result/value.txt", payload)


def test_publish_zip_verifies_and_atomically_replaces_destination(tmp_path: Path):
    source = tmp_path / "source.zip"
    destination = tmp_path / "published" / "result.zip"
    _write_zip(source)
    destination.parent.mkdir()
    destination.write_bytes(b"old")

    result = PUBLISH.publish_zip(source, destination)

    assert result == destination
    assert source.is_file()
    with zipfile.ZipFile(destination) as archive:
        assert archive.read("result/value.txt") == b"payload\n"
    assert not list(destination.parent.glob(".result.zip.partial.*"))


def test_publish_zip_remove_source_happens_after_success(tmp_path: Path):
    source = tmp_path / "source.zip"
    destination = tmp_path / "result.zip"
    _write_zip(source)

    PUBLISH.publish_zip(source, destination, remove_source=True)

    assert destination.is_file()
    assert not source.exists()


def test_remove_source_preserves_destination_reached_through_symlinked_parent(
    tmp_path: Path,
):
    real_parent = tmp_path / "real"
    real_parent.mkdir()
    alias_parent = tmp_path / "alias"
    alias_parent.symlink_to(real_parent, target_is_directory=True)
    source = real_parent / "result.zip"
    destination = alias_parent / "result.zip"
    _write_zip(source)

    PUBLISH.publish_zip(source, destination, remove_source=True)

    assert source.is_file()
    assert destination.is_file()
    with zipfile.ZipFile(destination) as archive:
        assert archive.read("result/value.txt") == b"payload\n"


def test_corrupt_source_does_not_replace_destination(tmp_path: Path):
    source = tmp_path / "source.zip"
    destination = tmp_path / "result.zip"
    _write_zip(source)
    source.write_bytes(source.read_bytes()[:-12])
    destination.write_bytes(b"old")

    with pytest.raises(PUBLISH.ZipPublishError, match="Failed to verify"):
        PUBLISH.publish_zip(source, destination, remove_source=True)

    assert destination.read_bytes() == b"old"
    assert source.is_file()
    assert not list(tmp_path.glob(".result.zip.partial.*"))


def test_copy_failure_does_not_replace_destination(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    source = tmp_path / "source.zip"
    destination = tmp_path / "result.zip"
    _write_zip(source)
    destination.write_bytes(b"old")

    def fail_copy(input_handle, output_handle, length=0):
        output_handle.write(input_handle.read(8))
        raise OSError("injected copy failure")

    monkeypatch.setattr(PUBLISH.shutil, "copyfileobj", fail_copy)

    with pytest.raises(OSError, match="injected copy failure"):
        PUBLISH.publish_zip(source, destination)

    assert destination.read_bytes() == b"old"
    assert not list(tmp_path.glob(".result.zip.partial.*"))


def test_same_size_source_change_with_restored_mtime_is_detected(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    source = tmp_path / "source.zip"
    destination = tmp_path / "result.zip"
    _write_zip(source)
    destination.write_bytes(b"old")
    original_copy = PUBLISH.shutil.copyfileobj

    def copy_then_mutate(input_handle, output_handle, length=0):
        original_copy(input_handle, output_handle, length=length)
        original_mtime_ns = source.stat().st_mtime_ns
        source.write_bytes(b"X" * source.stat().st_size)
        os.utime(source, ns=(original_mtime_ns, original_mtime_ns))

    monkeypatch.setattr(PUBLISH.shutil, "copyfileobj", copy_then_mutate)

    with pytest.raises(PUBLISH.ZipPublishError, match="changed while being published"):
        PUBLISH.publish_zip(source, destination)

    assert destination.read_bytes() == b"old"
    assert not list(tmp_path.glob(".result.zip.partial.*"))


def test_expected_prefix_rejects_unexpected_members_without_publication(tmp_path: Path):
    source = tmp_path / "source.zip"
    destination = tmp_path / "result.zip"
    with zipfile.ZipFile(source, "w") as archive:
        archive.writestr("../escaped.txt", b"unsafe")
    destination.write_bytes(b"old")

    with pytest.raises(PUBLISH.ZipPublishError, match="Failed to verify"):
        PUBLISH.publish_zip(source, destination, expected_prefix="result")

    assert destination.read_bytes() == b"old"

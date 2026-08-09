from __future__ import annotations

import os
import stat
import zipfile
from pathlib import Path

import pytest

from workflow.support import safe_zip_extract as SAFE


def _write_zip(path: Path, members: list[tuple[str | zipfile.ZipInfo, bytes]]) -> None:
    with zipfile.ZipFile(path, "w") as archive:
        for name, payload in members:
            archive.writestr(name, payload)


def test_extract_expected_prefix_replaces_directory_atomically(tmp_path: Path):
    archive_path = tmp_path / "result.zip"
    _write_zip(
        archive_path,
        [
            ("result/", b""),
            ("result/tree.nwk", b"(A,B);\n"),
            ("result/nested/value.txt", b"value\n"),
        ],
    )
    destination = tmp_path / "output" / "result"
    destination.mkdir(parents=True)
    (destination / "old.txt").write_text("old\n", encoding="utf-8")

    extracted = SAFE.extract_expected_prefix(
        archive_path,
        destination.parent,
        "result",
    )

    assert extracted == destination.resolve()
    assert (destination / "tree.nwk").read_text(encoding="utf-8") == "(A,B);\n"
    assert (destination / "nested" / "value.txt").read_text(encoding="utf-8") == "value\n"
    assert not (destination / "old.txt").exists()
    assert not list(destination.parent.glob(".result.extract.*"))


@pytest.mark.parametrize(
    "member_name",
    [
        "../victim/payload",
        "/tmp/payload",
        "other/payload",
        "result/../victim",
        "result\\payload",
        "result//payload",
    ],
)
def test_unsafe_members_are_rejected_without_replacing_existing_directory(
    tmp_path: Path,
    member_name: str,
):
    archive_path = tmp_path / "unsafe.zip"
    _write_zip(archive_path, [(member_name, b"malicious")])
    destination = tmp_path / "output" / "result"
    destination.mkdir(parents=True)
    sentinel = destination / "sentinel.txt"
    sentinel.write_text("keep\n", encoding="utf-8")

    with pytest.raises(SAFE.SafeZipError):
        SAFE.extract_expected_prefix(archive_path, destination.parent, "result")

    assert sentinel.read_text(encoding="utf-8") == "keep\n"
    assert not list(destination.parent.glob(".result.extract.*"))


def test_symlink_member_is_rejected(tmp_path: Path):
    archive_path = tmp_path / "symlink.zip"
    symlink = zipfile.ZipInfo("result/link")
    symlink.create_system = 3
    symlink.external_attr = (stat.S_IFLNK | 0o777) << 16
    _write_zip(archive_path, [(symlink, b"../../victim")])

    with pytest.raises(SAFE.SafeZipError, match="Symlinked ZIP members"):
        SAFE.extract_expected_prefix(archive_path, tmp_path / "output", "result")


def test_special_member_is_rejected(tmp_path: Path):
    archive_path = tmp_path / "fifo.zip"
    fifo = zipfile.ZipInfo("result/fifo")
    fifo.create_system = 3
    fifo.external_attr = (stat.S_IFIFO | 0o600) << 16
    _write_zip(archive_path, [(fifo, b"")])

    with pytest.raises(SAFE.SafeZipError, match="Special ZIP members"):
        SAFE.extract_expected_prefix(archive_path, tmp_path / "output", "result")


def test_special_permission_bits_are_not_restored(tmp_path: Path):
    archive_path = tmp_path / "setuid.zip"
    executable = zipfile.ZipInfo("result/tool")
    executable.create_system = 3
    executable.external_attr = (stat.S_IFREG | stat.S_ISUID | 0o755) << 16
    _write_zip(archive_path, [(executable, b"#!/bin/sh\n")])

    extracted = SAFE.extract_expected_prefix(
        archive_path,
        tmp_path / "output",
        "result",
    )

    assert stat.S_IMODE((extracted / "tool").stat().st_mode) == 0o755


def test_inode_preflight_treats_zero_as_exhausted(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    archive_path = tmp_path / "result.zip"
    _write_zip(archive_path, [("result/value.txt", b"value")])
    real_statvfs = SAFE.os.statvfs

    def exhausted(path):
        values = list(real_statvfs(path))
        values[7] = 0
        return os.statvfs_result(values)

    monkeypatch.setattr(SAFE.os, "statvfs", exhausted)

    with pytest.raises(SAFE.SafeZipError, match="available=0"):
        SAFE.extract_expected_prefix(archive_path, tmp_path / "output", "result")


def test_corrupt_archive_preserves_existing_directory(tmp_path: Path):
    archive_path = tmp_path / "truncated.zip"
    _write_zip(archive_path, [("result/tree.nwk", b"(A,B);\n")])
    archive_path.write_bytes(archive_path.read_bytes()[:-12])
    destination = tmp_path / "output" / "result"
    destination.mkdir(parents=True)
    sentinel = destination / "sentinel.txt"
    sentinel.write_text("keep\n", encoding="utf-8")

    with pytest.raises(SAFE.SafeZipError):
        SAFE.extract_expected_prefix(archive_path, destination.parent, "result")

    assert sentinel.read_text(encoding="utf-8") == "keep\n"
    assert not list(destination.parent.glob(".result.extract.*"))


def test_failed_destination_replace_restores_previous_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    archive_path = tmp_path / "result.zip"
    _write_zip(archive_path, [("result/new.txt", b"new\n")])
    destination = tmp_path / "output" / "result"
    destination.mkdir(parents=True)
    sentinel = destination / "sentinel.txt"
    sentinel.write_text("keep\n", encoding="utf-8")
    real_replace = SAFE.os.replace

    def fail_new_directory_replace(source: os.PathLike[str], target: os.PathLike[str]):
        source_path = Path(source)
        if Path(target) == destination and ".extract.partial." in source_path.parent.name:
            raise OSError("injected replace failure")
        return real_replace(source, target)

    monkeypatch.setattr(SAFE.os, "replace", fail_new_directory_replace)

    with pytest.raises(OSError, match="injected replace failure"):
        SAFE.extract_expected_prefix(archive_path, destination.parent, "result")

    assert sentinel.read_text(encoding="utf-8") == "keep\n"
    assert not (destination / "new.txt").exists()
    assert not list(destination.parent.glob(".result.extract.*"))

#!/usr/bin/env python3
"""Validate and atomically extract one expected top-level ZIP directory."""

from __future__ import annotations

import argparse
import os
import shutil
import stat
import sys
import time
import uuid
import zipfile
from pathlib import Path, PurePosixPath
from typing import Sequence


class SafeZipError(RuntimeError):
    """Raised when a ZIP cannot be extracted without leaving its expected prefix."""


def _fsync_directory(path: Path) -> None:
    try:
        descriptor = os.open(path, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(descriptor)
    except OSError:
        pass
    finally:
        os.close(descriptor)


def _preflight_extraction(
    destination_root: Path,
    files: Sequence[zipfile.ZipInfo],
    directories: Sequence[zipfile.ZipInfo],
) -> None:
    filesystem_stats = os.statvfs(destination_root)
    block_size = int(
        filesystem_stats.f_frsize or filesystem_stats.f_bsize or 4096
    )
    required_bytes = block_size
    implied_directories: set[PurePosixPath] = set()
    for info in files:
        required_bytes += max(
            block_size,
            ((int(info.file_size) + block_size - 1) // block_size) * block_size,
        )
        parent = PurePosixPath(info.filename).parent
        while parent.as_posix() != ".":
            implied_directories.add(parent)
            parent = parent.parent
    for info in directories:
        implied_directories.add(PurePosixPath(info.filename))
    required_bytes += len(implied_directories) * block_size
    available_bytes = int(shutil.disk_usage(destination_root).free)
    if required_bytes > available_bytes:
        raise SafeZipError(
            "Insufficient filesystem space for ZIP extraction: "
            f"required={required_bytes}, available={available_bytes}"
        )
    required_inodes = len(files) + len(implied_directories) + 1
    total_inodes = int(filesystem_stats.f_files)
    available_inodes = int(filesystem_stats.f_favail)
    if (
        total_inodes > 0
        and available_inodes >= 0
        and required_inodes > available_inodes
    ):
        raise SafeZipError(
            "Insufficient filesystem inodes for ZIP extraction: "
            f"required={required_inodes}, available={available_inodes}"
        )


def _validated_prefix(value: str) -> str:
    prefix = PurePosixPath(value)
    if (
        not value
        or prefix.is_absolute()
        or len(prefix.parts) != 1
        or prefix.as_posix() != value
        or value in {".", ".."}
        or "\\" in value
    ):
        raise SafeZipError(
            f"Expected ZIP prefix must be one safe top-level directory name: {value!r}"
        )
    return value


def validated_members(
    archive: zipfile.ZipFile,
    archive_path: Path | str,
    expected_prefix: str,
) -> tuple[list[zipfile.ZipInfo], list[zipfile.ZipInfo]]:
    """Return validated file and directory members below ``expected_prefix``."""

    expected_prefix = _validated_prefix(expected_prefix)
    archive_path = Path(archive_path)
    prefix = f"{expected_prefix}/"
    seen: set[str] = set()
    file_names: set[str] = set()
    files: list[zipfile.ZipInfo] = []
    directories: list[zipfile.ZipInfo] = []
    for info in archive.infolist():
        member = PurePosixPath(info.filename)
        canonical_name = member.as_posix()
        canonical_archive_name = (
            f"{canonical_name}/" if info.is_dir() else canonical_name
        )
        root_directory = info.is_dir() and canonical_name == expected_prefix
        if (
            not canonical_name
            or canonical_name in seen
            or info.filename.startswith("/")
            or "\\" in info.filename
            or canonical_archive_name != info.filename
            or "." in member.parts
            or ".." in member.parts
            or (
                not root_directory
                and (
                    not info.filename.startswith(prefix)
                    or len(member.parts) < 2
                )
            )
        ):
            raise SafeZipError(
                f"Unsafe or unexpected ZIP member in {archive_path}: {info.filename!r}"
            )
        seen.add(canonical_name)
        mode = info.external_attr >> 16
        if stat.S_ISLNK(mode):
            raise SafeZipError(
                f"Symlinked ZIP members are unsupported in {archive_path}: "
                f"{info.filename}"
            )
        file_type = stat.S_IFMT(mode)
        expected_types = {0, stat.S_IFDIR} if info.is_dir() else {0, stat.S_IFREG}
        if file_type not in expected_types:
            raise SafeZipError(
                f"Special ZIP members are unsupported in {archive_path}: "
                f"{info.filename}"
            )
        if info.is_dir():
            directories.append(info)
        else:
            files.append(info)
            file_names.add(canonical_name)

    for canonical_name in seen:
        parent = PurePosixPath(canonical_name).parent
        while parent.as_posix() not in {".", expected_prefix}:
            if parent.as_posix() in file_names:
                raise SafeZipError(
                    f"ZIP member is nested below a file in {archive_path}: "
                    f"{canonical_name!r}"
                )
            parent = parent.parent
    return files, directories


def extract_expected_prefix(
    archive_path: Path | str,
    destination_root: Path | str,
    expected_prefix: str,
    *,
    replace: bool = True,
) -> Path:
    """Extract a validated prefix through a sibling temporary directory."""

    expected_prefix = _validated_prefix(expected_prefix)
    archive_path = Path(archive_path)
    if archive_path.is_symlink() or not archive_path.is_file():
        raise SafeZipError(f"ZIP archive is not a regular file: {archive_path}")
    destination_root = Path(destination_root).resolve()
    destination_root.mkdir(parents=True, exist_ok=True)
    destination = destination_root / expected_prefix
    if destination.is_symlink():
        raise SafeZipError(f"Symlinked ZIP destination is unsupported: {destination}")
    if destination.exists() and not destination.is_dir():
        raise SafeZipError(f"ZIP destination is not a directory: {destination}")
    if destination.exists() and not replace:
        raise SafeZipError(f"ZIP destination already exists: {destination}")

    temporary_root = destination_root / (
        f".{expected_prefix}.extract.partial.{os.getpid()}.{uuid.uuid4().hex}"
    )
    temporary_root.mkdir()
    extracted_prefix = temporary_root / expected_prefix
    try:
        try:
            with zipfile.ZipFile(archive_path, "r") as archive:
                files, directories = validated_members(
                    archive,
                    archive_path,
                    expected_prefix,
                )
                _preflight_extraction(destination_root, files, directories)
                extracted_prefix.mkdir()
                directory_metadata: list[tuple[Path, int, tuple[int, ...]]] = []
                for info in sorted(
                    directories,
                    key=lambda item: len(PurePosixPath(item.filename).parts),
                ):
                    relative = PurePosixPath(info.filename).relative_to(expected_prefix)
                    target = extracted_prefix.joinpath(*relative.parts)
                    target.mkdir(parents=True, exist_ok=True)
                    directory_metadata.append(
                        (target, info.external_attr >> 16, info.date_time)
                    )
                for info in files:
                    relative = PurePosixPath(info.filename).relative_to(expected_prefix)
                    target = extracted_prefix.joinpath(*relative.parts)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with archive.open(info, "r") as source, target.open("xb") as output:
                        shutil.copyfileobj(source, output, length=1024 * 1024)
                        output.flush()
                        os.fsync(output.fileno())
                    mode = info.external_attr >> 16
                    if mode:
                        target.chmod(stat.S_IMODE(mode) & 0o777)
                    modified = time.mktime((*info.date_time, 0, 0, -1))
                    os.utime(target, (modified, modified))
                for target, mode, date_time in reversed(directory_metadata):
                    if mode:
                        target.chmod(stat.S_IMODE(mode) & 0o777)
                    modified = time.mktime((*date_time, 0, 0, -1))
                    os.utime(target, (modified, modified))
                extracted_directories = [
                    path
                    for path in temporary_root.rglob("*")
                    if path.is_dir() and not path.is_symlink()
                ]
                for path in sorted(
                    extracted_directories,
                    key=lambda item: len(item.parts),
                    reverse=True,
                ):
                    _fsync_directory(path)
                _fsync_directory(temporary_root)
        except (
            OSError,
            RuntimeError,
            ValueError,
            OverflowError,
            zipfile.BadZipFile,
        ) as exc:
            raise SafeZipError(f"Failed to extract {archive_path}: {exc}") from exc

        backup = destination_root / (
            f".{expected_prefix}.extract.backup.{os.getpid()}.{uuid.uuid4().hex}"
        )
        destination_was_moved = False
        destination_replaced = False
        try:
            if destination.exists():
                os.replace(destination, backup)
                destination_was_moved = True
            os.replace(extracted_prefix, destination)
            destination_replaced = True
        except BaseException:
            if destination_was_moved and not destination.exists() and backup.exists():
                os.replace(backup, destination)
            raise
        finally:
            if destination_replaced and backup.exists():
                shutil.rmtree(backup)
        _fsync_directory(destination_root)
        return destination
    finally:
        shutil.rmtree(temporary_root, ignore_errors=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--destination-root", required=True, type=Path)
    parser.add_argument("--expected-prefix", required=True)
    parser.add_argument(
        "--keep-existing",
        action="store_true",
        help="Fail instead of replacing an existing expected-prefix directory",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        path = extract_expected_prefix(
            args.archive,
            args.destination_root,
            args.expected_prefix,
            replace=not args.keep_existing,
        )
    except SafeZipError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

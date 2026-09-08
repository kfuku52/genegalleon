#!/usr/bin/env python3
"""Manage high-file-count species-tree stage directories as visible ZIP files."""

from __future__ import annotations

import argparse
import contextlib
import fcntl
import fnmatch
import hashlib
import json
import os
import shutil
import stat
import sys
import time
import uuid
import zipfile
from pathlib import Path, PurePosixPath
from typing import Callable, Iterable, Iterator, Sequence

SCHEMA_VERSION = 1
COMMENT_PREFIX = b"GeneGalleon species-tree directory archive\n"
MANAGED_DIRECTORIES = (
    "single_copy_cds_fasta",
    "single_copy_mafft",
    "single_copy_trimal",
    "single_copy_iqtree_pep",
    "single_copy_iqtree_dna",
)
ALREADY_COMPRESSED_SUFFIXES = {
    ".7z",
    ".bz2",
    ".gz",
    ".jpeg",
    ".jpg",
    ".pdf",
    ".png",
    ".rar",
    ".tgz",
    ".xz",
    ".zip",
    ".zst",
}


class SpeciesTreeArchiveError(RuntimeError):
    """Raised when a managed species-tree archive is unsafe or inconsistent."""


class ProgressReporter:
    def __init__(self, interval: float) -> None:
        self.interval = max(0.0, float(interval))
        self.last = 0.0

    def update(self, *, force: bool = False, **fields: object) -> None:
        now = time.monotonic()
        if not force and self.interval > 0 and now - self.last < self.interval:
            return
        self.last = now
        print(
            "progress\t" + json.dumps(fields, sort_keys=True),
            file=sys.stderr,
            flush=True,
        )


def _directory_name(value: str) -> str:
    if value not in MANAGED_DIRECTORIES:
        expected = ", ".join(MANAGED_DIRECTORIES)
        raise SpeciesTreeArchiveError(f"Unsupported species-tree directory {value!r}; expected one of: {expected}")
    return value


def _selected_directories(values: Sequence[str] | None) -> tuple[str, ...]:
    if not values:
        return MANAGED_DIRECTORIES
    selected: list[str] = []
    for value in values:
        name = _directory_name(value)
        if name not in selected:
            selected.append(name)
    return tuple(selected)


def _archive_path(root: Path, name: str) -> Path:
    return root / f"{name}.zip"


def _raw_path(root: Path, name: str) -> Path:
    return root / name


def _assert_regular_path(path: Path, description: str) -> None:
    if path.is_symlink():
        raise SpeciesTreeArchiveError(f"Symlinked {description} is not supported: {path}")


@contextlib.contextmanager
def _store_lock(root: Path) -> Iterator[None]:
    root.mkdir(parents=True, exist_ok=True)
    directory_fd = os.open(root, os.O_RDONLY)
    fcntl.flock(directory_fd, fcntl.LOCK_EX)
    try:
        yield
    finally:
        fcntl.flock(directory_fd, fcntl.LOCK_UN)
        os.close(directory_fd)


def _cleanup_stale_partials(root: Path, name: str) -> None:
    patterns = (
        f".{name}.zip.partial.*",
        f".{name}.materialize.partial.*",
    )
    for pattern in patterns:
        for path in root.glob(pattern):
            _assert_regular_path(path, "partial archive path")
            if path.is_dir():
                shutil.rmtree(path)
            elif path.exists():
                path.unlink()


def _fsync_directory(path: Path) -> None:
    try:
        directory_fd = os.open(path, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(directory_fd)
    except OSError:
        pass
    finally:
        os.close(directory_fd)


def _compression_for(path: Path, compression: str) -> int:
    if compression == "store":
        return zipfile.ZIP_STORED
    if compression == "deflate":
        return zipfile.ZIP_DEFLATED
    if path.suffix.lower() in ALREADY_COMPRESSED_SUFFIXES:
        return zipfile.ZIP_STORED
    return zipfile.ZIP_DEFLATED


def _source_inventory(raw_path: Path) -> tuple[list[Path], list[Path]]:
    _assert_regular_path(raw_path, "managed directory")
    directories: list[Path] = [raw_path]
    files: list[Path] = []

    def raise_walk_error(exc: OSError) -> None:
        raise SpeciesTreeArchiveError(f"Failed to scan {raw_path}: {exc}") from exc

    for current_root, directory_names, file_names in os.walk(raw_path, onerror=raise_walk_error):
        current = Path(current_root)
        for directory_name in directory_names:
            directory = current / directory_name
            _assert_regular_path(directory, "source directory")
            directories.append(directory)
        for file_name in file_names:
            path = current / file_name
            _assert_regular_path(path, "source file")
            if not path.is_file():
                raise SpeciesTreeArchiveError(f"Non-regular source file is not supported: {path}")
            files.append(path)
    return (
        sorted(
            directories,
            key=lambda path: (
                len(path.relative_to(raw_path).parts),
                path.relative_to(raw_path).as_posix(),
            ),
        ),
        sorted(files, key=lambda path: path.relative_to(raw_path).as_posix()),
    )


def _iter_source_files(raw_path: Path) -> list[Path]:
    return _source_inventory(raw_path)[1]


def _count_source_files_best_effort(raw_path: Path) -> int:
    """Count a possibly active stage without requiring a stable inventory."""
    _assert_regular_path(raw_path, "managed directory")
    count = 0

    def raise_walk_error(exc: OSError) -> None:
        raise SpeciesTreeArchiveError(f"Failed to scan {raw_path}: {exc}") from exc

    for current_root, directory_names, file_names in os.walk(raw_path, onerror=raise_walk_error):
        current = Path(current_root)
        for directory_name in directory_names:
            _assert_regular_path(current / directory_name, "source directory")
        for file_name in file_names:
            path = current / file_name
            _assert_regular_path(path, "source file")
            if path.exists() and not path.is_file():
                raise SpeciesTreeArchiveError(f"Non-regular source file is not supported: {path}")
            count += 1
    return count


def _raw_directory_stats(raw_path: Path) -> tuple[int, int]:
    files = _iter_source_files(raw_path)
    return len(files), sum(int(path.stat().st_size) for path in files)


def _raw_to_zip_peak_bytes(files: Sequence[Path]) -> int:
    logical_bytes = sum(int(path.stat().st_size) for path in files)
    return logical_bytes + max(64 * 1024, (logical_bytes + 99) // 100) + len(files) * 512 + 1024 * 1024


def _zip_to_raw_requirements(
    root: Path,
    file_members: Sequence[zipfile.ZipInfo],
    directory_members: Sequence[zipfile.ZipInfo],
) -> tuple[int, int]:
    filesystem_stats = os.statvfs(root)
    block_size = int(filesystem_stats.f_frsize or filesystem_stats.f_bsize or 4096)
    allocated_file_bytes = sum(
        ((int(info.file_size) + block_size - 1) // block_size) * block_size for info in file_members
    )
    directory_names = {PurePosixPath(info.filename).as_posix().rstrip("/") for info in directory_members}
    for info in file_members:
        parent = PurePosixPath(info.filename).parent
        while parent.as_posix() != ".":
            directory_names.add(parent.as_posix())
            parent = parent.parent
    # The temporary extraction root and every explicit or implicit archived
    # directory can consume a block and an inode while the source ZIP exists.
    directory_count = len(directory_names) + 1
    return (
        allocated_file_bytes + directory_count * block_size,
        len(file_members) + directory_count,
    )


def _signature(path: Path) -> tuple[int, int, int, int, int]:
    value = path.stat()
    return (
        value.st_dev,
        value.st_ino,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def _archive_comment(name: str, file_count: int, directory_count: int) -> bytes:
    payload = {
        "schema_version": SCHEMA_VERSION,
        "directory": name,
        "file_count": file_count,
        "directory_count": directory_count,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    return COMMENT_PREFIX + json.dumps(payload, sort_keys=True).encode("utf-8")


def _read_archive_metadata(archive: zipfile.ZipFile, path: Path) -> dict:
    if not archive.comment.startswith(COMMENT_PREFIX):
        raise SpeciesTreeArchiveError(f"Refusing to manage a ZIP without a GeneGalleon species-tree marker: {path}")
    try:
        payload = json.loads(archive.comment[len(COMMENT_PREFIX) :].decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SpeciesTreeArchiveError(f"Invalid archive metadata in {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise SpeciesTreeArchiveError(f"Archive metadata must be a JSON object in {path}")
    if type(payload.get("schema_version")) is not int or payload.get("schema_version") != SCHEMA_VERSION:
        raise SpeciesTreeArchiveError(
            f"Unsupported species-tree archive schema in {path}: {payload.get('schema_version')!r}"
        )
    return payload


def _validated_members(
    archive: zipfile.ZipFile,
    path: Path,
    expected_name: str,
    *,
    check_crc: bool = True,
) -> tuple[list[zipfile.ZipInfo], list[zipfile.ZipInfo]]:
    payload = _read_archive_metadata(archive, path)
    if payload.get("directory") != expected_name:
        raise SpeciesTreeArchiveError(f"Archive directory mismatch in {path}: {payload.get('directory')!r}")
    file_members: list[zipfile.ZipInfo] = []
    directory_members: list[zipfile.ZipInfo] = []
    seen: set[str] = set()
    file_names: set[str] = set()
    prefix = f"{expected_name}/"
    for info in archive.infolist():
        member = PurePosixPath(info.filename)
        canonical_name = member.as_posix()
        canonical_archive_name = f"{canonical_name}/" if info.is_dir() else canonical_name
        root_directory_member = info.is_dir() and canonical_name == expected_name
        if (
            canonical_name in seen
            or info.filename.startswith("/")
            or "\\" in info.filename
            or canonical_archive_name != info.filename
            or ".." in member.parts
            or (not root_directory_member and (not info.filename.startswith(prefix) or len(member.parts) < 2))
        ):
            raise SpeciesTreeArchiveError(f"Unsafe ZIP member in {path}: {info.filename!r}")
        seen.add(canonical_name)
        mode = info.external_attr >> 16
        if stat.S_ISLNK(mode):
            raise SpeciesTreeArchiveError(f"Symlinked ZIP members are not supported in {path}: {info.filename}")
        if info.is_dir():
            directory_members.append(info)
        else:
            file_members.append(info)
            file_names.add(canonical_name)
    for canonical_name in seen:
        parent = PurePosixPath(canonical_name).parent
        while parent.as_posix() not in {".", expected_name}:
            if parent.as_posix() in file_names:
                raise SpeciesTreeArchiveError(f"ZIP member is nested below a file in {path}: {canonical_name!r}")
            parent = parent.parent
    expected_count = payload.get("file_count")
    if type(expected_count) is not int or expected_count != len(file_members):
        raise SpeciesTreeArchiveError(
            f"Archive file count differs from its metadata in {path}: "
            f"metadata={expected_count!r}, actual={len(file_members)}"
        )
    expected_directory_count = payload.get("directory_count")
    if expected_directory_count is not None and (
        type(expected_directory_count) is not int or expected_directory_count != len(directory_members)
    ):
        raise SpeciesTreeArchiveError(
            f"Archive directory count differs from its metadata in {path}: "
            f"metadata={expected_directory_count!r}, actual={len(directory_members)}"
        )
    if check_crc:
        bad_member = archive.testzip()
        if bad_member is not None:
            raise SpeciesTreeArchiveError(f"CRC verification failed in {path}: {bad_member}")
    return file_members, directory_members


def _verify_member_payloads(
    archive: zipfile.ZipFile,
    path: Path,
    members: Sequence[zipfile.ZipInfo],
    *,
    progress_callback: Callable[..., None] | None = None,
) -> int:
    total_bytes = sum(int(info.file_size) for info in members)
    completed_bytes = 0
    completed_files = 0
    if progress_callback is not None:
        progress_callback(
            force=True,
            phase="verifying",
            current_zip=path.name,
            verify_files_completed=0,
            verify_files_total=len(members),
            verify_bytes_completed=0,
            verify_bytes_total=total_bytes,
        )
    for info in members:
        try:
            with archive.open(info, "r") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    completed_bytes += len(chunk)
                    if progress_callback is not None:
                        progress_callback(
                            phase="verifying",
                            current_zip=path.name,
                            verify_files_completed=completed_files,
                            verify_files_total=len(members),
                            verify_bytes_completed=completed_bytes,
                            verify_bytes_total=total_bytes,
                        )
        except zipfile.BadZipFile as exc:
            raise SpeciesTreeArchiveError(f"CRC verification failed in {path}: {info.filename}") from exc
        completed_files += 1
    if progress_callback is not None:
        progress_callback(
            force=True,
            phase="verifying",
            current_zip=path.name,
            verify_files_completed=completed_files,
            verify_files_total=len(members),
            verify_bytes_completed=completed_bytes,
            verify_bytes_total=total_bytes,
        )
    return completed_bytes


def verify_archive(
    root: Path,
    name: str,
    *,
    check_crc: bool = True,
    progress_callback: Callable[..., None] | None = None,
) -> dict:
    path = _archive_path(root, name)
    _assert_regular_path(path, "ZIP archive")
    if not path.exists():
        return {
            "directory": name,
            "state": "absent",
            "files": 0,
            "verification_mode": "deep" if check_crc else "quick",
            "verified_payload_bytes": 0,
        }
    if not path.is_file():
        raise SpeciesTreeArchiveError(f"Archive path is not a file: {path}")
    try:
        with zipfile.ZipFile(path, "r") as archive:
            members, directory_members = _validated_members(
                archive,
                path,
                name,
                check_crc=False,
            )
            verified_payload_bytes = (
                _verify_member_payloads(
                    archive,
                    path,
                    members,
                    progress_callback=progress_callback,
                )
                if check_crc
                else 0
            )
    except (OSError, RuntimeError, zipfile.BadZipFile) as exc:
        raise SpeciesTreeArchiveError(f"Failed to verify {path}: {exc}") from exc
    required_bytes, required_inodes = _zip_to_raw_requirements(
        root,
        members,
        directory_members,
    )
    return {
        "directory": name,
        "state": "archived",
        "files": len(members),
        "logical_bytes": sum(int(info.file_size) for info in members),
        "zip_bytes": path.stat().st_size,
        "raw_materialize_allocated_bytes": required_bytes,
        "raw_materialize_required_inodes": required_inodes,
        "verification_mode": "deep" if check_crc else "quick",
        "verified_payload_bytes": verified_payload_bytes,
    }


def pack_directory(
    root: Path,
    name: str,
    *,
    compression: str = "adaptive",
    compression_level: int = 6,
    available_bytes: int | None = None,
    progress_callback: Callable[..., None] | None = None,
    include_metrics: bool = False,
) -> dict:
    name = _directory_name(name)
    _validate_options(compression, compression_level)
    raw_path = _raw_path(root, name)
    archive_path = _archive_path(root, name)
    with _store_lock(root):
        _cleanup_stale_partials(root, name)
        _assert_regular_path(archive_path, "ZIP archive")
        if not raw_path.exists():
            if archive_path.exists():
                return verify_archive(root, name, check_crc=False)
            return {"directory": name, "state": "absent", "files": 0}
        if not raw_path.is_dir():
            raise SpeciesTreeArchiveError(f"Managed path is not a directory: {raw_path}")
        if archive_path.exists():
            verify_archive(root, name, check_crc=False)
            raise SpeciesTreeArchiveError(
                f"Both raw and ZIP forms exist for {name}. Run materialize first to "
                "merge them safely before packing again."
            )

        directories, files = _source_inventory(raw_path)
        estimated_peak_new_bytes = _raw_to_zip_peak_bytes(files)
        filesystem_free_bytes = int(shutil.disk_usage(root).free)
        effective_available_bytes = filesystem_free_bytes
        if available_bytes is not None:
            if available_bytes < 0:
                raise SpeciesTreeArchiveError("available_bytes must be non-negative")
            effective_available_bytes = min(
                effective_available_bytes,
                int(available_bytes),
            )
        if estimated_peak_new_bytes > effective_available_bytes:
            raise SpeciesTreeArchiveError(
                "Insufficient temporary space for raw-to-ZIP conversion: "
                f"estimated_peak_new_bytes={estimated_peak_new_bytes}, "
                f"effective_available_bytes={effective_available_bytes}. "
                "Pass quota remaining as --available-bytes for a quota-aware preflight."
            )
        signatures = {path: _signature(path) for path in directories + files}
        archived_digests: dict[Path, str] = {}
        partial_path = root / f".{name}.zip.partial.{os.getpid()}.{uuid.uuid4().hex}"
        total_source_bytes = sum(int(path.stat().st_size) for path in files)
        completed_source_bytes = 0
        completed_files = 0
        if progress_callback is not None:
            progress_callback(
                force=True,
                phase="packing",
                directory=name,
                files_completed=0,
                files_total=len(files),
                bytes_completed=0,
                bytes_total=total_source_bytes,
            )
        try:
            with zipfile.ZipFile(
                partial_path,
                "w",
                compression=zipfile.ZIP_DEFLATED,
                compresslevel=compression_level,
                allowZip64=True,
            ) as archive:
                for path in directories:
                    relative = path.relative_to(raw_path).as_posix()
                    archive_name = f"{name}/" if relative == "." else f"{name}/{relative}/"
                    info = zipfile.ZipInfo.from_file(path, arcname=archive_name, strict_timestamps=False)
                    info.compress_type = zipfile.ZIP_STORED
                    archive.writestr(info, b"")
                for path in files:
                    digest = hashlib.sha256()
                    relative = path.relative_to(raw_path).as_posix()
                    info = zipfile.ZipInfo.from_file(
                        path,
                        arcname=f"{name}/{relative}",
                        strict_timestamps=False,
                    )
                    info.compress_type = _compression_for(path, compression)
                    with path.open("rb") as source, archive.open(info, "w", force_zip64=True) as destination:

                        class ProgressDestination:
                            def write(
                                self,
                                chunk: bytes,
                                *,
                                _digest: "hashlib._Hash" = digest,
                                _completed_files: int = completed_files,
                            ) -> int:
                                nonlocal completed_source_bytes
                                written = destination.write(chunk)
                                _digest.update(memoryview(chunk)[:written])
                                completed_source_bytes += int(written)
                                if progress_callback is not None:
                                    progress_callback(
                                        phase="packing",
                                        directory=name,
                                        files_completed=_completed_files,
                                        files_total=len(files),
                                        bytes_completed=completed_source_bytes,
                                        bytes_total=total_source_bytes,
                                    )
                                return int(written)

                        shutil.copyfileobj(
                            source,
                            ProgressDestination(),
                            length=1024 * 1024,
                        )
                    archived_digests[path] = digest.hexdigest()
                    completed_files += 1
                archive.comment = _archive_comment(name, len(files), len(directories))

            with zipfile.ZipFile(partial_path, "r") as archive:
                members, _ = _validated_members(
                    archive,
                    partial_path,
                    name,
                    check_crc=False,
                )
                _verify_member_payloads(
                    archive,
                    partial_path,
                    members,
                    progress_callback=progress_callback,
                )
            for path, signature in signatures.items():
                if not path.exists() or path.is_symlink() or _signature(path) != signature:
                    raise SpeciesTreeArchiveError(f"Source changed while {name} was being archived: {path}")
            for path, archived_digest in archived_digests.items():
                current_digest = hashlib.sha256()
                with path.open("rb") as source:
                    for chunk in iter(lambda: source.read(1024 * 1024), b""):
                        current_digest.update(chunk)
                if (
                    current_digest.hexdigest() != archived_digest
                    or path.is_symlink()
                    or _signature(path) != signatures[path]
                ):
                    raise SpeciesTreeArchiveError(f"Source changed while {name} was being archived: {path}")
            current_directories, current_files = _source_inventory(raw_path)
            if current_directories != directories or current_files != files:
                raise SpeciesTreeArchiveError(f"Source inventory changed while {name} was being archived")

            with partial_path.open("rb") as handle:
                os.fsync(handle.fileno())
            os.replace(partial_path, archive_path)
            _fsync_directory(root)
            shutil.rmtree(raw_path)
            _fsync_directory(root)
        finally:
            partial_path.unlink(missing_ok=True)
    result = {
        "directory": name,
        "state": "archived",
        "files": len(files),
    }
    if include_metrics:
        result.update(
            {
                "zip_bytes": archive_path.stat().st_size,
                "estimated_peak_new_bytes": estimated_peak_new_bytes,
                "filesystem_free_bytes": filesystem_free_bytes,
                "effective_available_bytes": effective_available_bytes,
            }
        )
    return result


def _safe_destination(base: Path, member_name: str) -> Path:
    destination = base.joinpath(*PurePosixPath(member_name).parts)
    try:
        destination.resolve().relative_to(base.resolve())
    except ValueError as exc:
        raise SpeciesTreeArchiveError(f"ZIP member escapes its extraction root: {member_name!r}") from exc
    return destination


def _extract_archive_to(
    archive_path: Path,
    name: str,
    destination: Path,
    *,
    progress_callback: Callable[..., None] | None = None,
) -> int:
    try:
        with zipfile.ZipFile(archive_path, "r") as archive:
            members, directory_members = _validated_members(archive, archive_path, name, check_crc=False)
            directory_metadata: list[tuple[Path, int, tuple[int, ...]]] = []
            total_bytes = sum(int(info.file_size) for info in members)
            completed_bytes = 0
            completed_files = 0
            if progress_callback is not None:
                progress_callback(
                    force=True,
                    phase="materializing",
                    directory=name,
                    files_completed=0,
                    files_total=len(members),
                    bytes_completed=0,
                    bytes_total=total_bytes,
                )
            for info in sorted(
                directory_members,
                key=lambda item: len(PurePosixPath(item.filename).parts),
            ):
                target = _safe_destination(destination, info.filename)
                target.mkdir(parents=True, exist_ok=True)
                directory_metadata.append((target, info.external_attr >> 16, info.date_time))
            for info in members:
                target = _safe_destination(destination, info.filename)
                target.parent.mkdir(parents=True, exist_ok=True)
                with archive.open(info, "r") as source, target.open("wb") as output:

                    class ProgressOutput:
                        def write(
                            self,
                            chunk: bytes,
                            *,
                            _completed_files: int = completed_files,
                        ) -> int:
                            nonlocal completed_bytes
                            written = output.write(chunk)
                            completed_bytes += int(written)
                            if progress_callback is not None:
                                progress_callback(
                                    phase="materializing",
                                    directory=name,
                                    files_completed=_completed_files,
                                    files_total=len(members),
                                    bytes_completed=completed_bytes,
                                    bytes_total=total_bytes,
                                )
                            return int(written)

                    shutil.copyfileobj(
                        source,
                        ProgressOutput(),
                        length=1024 * 1024,
                    )
                    output.flush()
                    os.fsync(output.fileno())
                mode = info.external_attr >> 16
                if mode:
                    target.chmod(stat.S_IMODE(mode))
                modified = time.mktime((*info.date_time, 0, 0, -1))
                os.utime(target, (modified, modified))
                completed_files += 1
            for target, mode, date_time in reversed(directory_metadata):
                if mode:
                    target.chmod(stat.S_IMODE(mode))
                modified = time.mktime((*date_time, 0, 0, -1))
                os.utime(target, (modified, modified))
            return len(members)
    except (OSError, RuntimeError, zipfile.BadZipFile) as exc:
        raise SpeciesTreeArchiveError(f"Failed to extract {archive_path}: {exc}") from exc


def materialize_directory(
    root: Path,
    name: str,
    *,
    available_bytes: int | None = None,
    progress_callback: Callable[..., None] | None = None,
    include_metrics: bool = False,
) -> dict:
    name = _directory_name(name)
    raw_path = _raw_path(root, name)
    archive_path = _archive_path(root, name)
    with _store_lock(root):
        _cleanup_stale_partials(root, name)
        _assert_regular_path(raw_path, "managed directory")
        _assert_regular_path(archive_path, "ZIP archive")
        if not archive_path.exists():
            if not raw_path.exists():
                return {"directory": name, "state": "absent", "files": 0}
            if not raw_path.is_dir():
                raise SpeciesTreeArchiveError(f"Managed path is not a directory: {raw_path}")
            return {
                "directory": name,
                "state": "raw",
                "files": len(_iter_source_files(raw_path)),
            }
        if not archive_path.is_file():
            raise SpeciesTreeArchiveError(f"Archive path is not a file: {archive_path}")

        with zipfile.ZipFile(archive_path, "r") as archive:
            members, directory_members = _validated_members(
                archive,
                archive_path,
                name,
                check_crc=False,
            )
        required_bytes, required_inodes = _zip_to_raw_requirements(
            root,
            members,
            directory_members,
        )
        filesystem_free_bytes = int(shutil.disk_usage(root).free)
        effective_available_bytes = filesystem_free_bytes
        if available_bytes is not None:
            if int(available_bytes) < 0:
                raise SpeciesTreeArchiveError("available_bytes must be non-negative")
            effective_available_bytes = min(
                effective_available_bytes,
                int(available_bytes),
            )
        if required_bytes > effective_available_bytes:
            raise SpeciesTreeArchiveError(
                "Insufficient temporary space for ZIP-to-raw conversion: "
                f"allocated_required={required_bytes}, "
                f"effective_available_bytes={effective_available_bytes}. "
                "Pass quota remaining as --available-bytes for a quota-aware "
                "preflight."
            )
        filesystem_stats = os.statvfs(root)
        total_inodes = int(filesystem_stats.f_files)
        available_inodes = int(filesystem_stats.f_favail)
        if (
            total_inodes > 0
            and available_inodes >= 0
            and required_inodes > available_inodes
        ):
            raise SpeciesTreeArchiveError(
                "Insufficient filesystem inodes for ZIP-to-raw conversion: "
                f"required={required_inodes}, available={available_inodes}"
            )

        temporary = root / f".{name}.materialize.partial.{os.getpid()}.{uuid.uuid4().hex}"
        temporary.mkdir()
        try:
            (temporary / name).mkdir()
            file_count = _extract_archive_to(
                archive_path,
                name,
                temporary,
                progress_callback=progress_callback,
            )
            extracted = temporary / name
            if raw_path.exists():
                if not raw_path.is_dir():
                    raise SpeciesTreeArchiveError(f"Managed path is not a directory: {raw_path}")
                _source_inventory(raw_path)
                extracted_directories, extracted_files = _source_inventory(extracted)
                new_directories: list[tuple[Path, Path]] = []
                for source_directory in extracted_directories[1:]:
                    relative = source_directory.relative_to(extracted)
                    destination_directory = raw_path / relative
                    if destination_directory.exists():
                        if destination_directory.is_symlink() or not destination_directory.is_dir():
                            raise SpeciesTreeArchiveError(
                                f"Live path conflicts with an archived directory: {destination_directory}"
                            )
                    else:
                        destination_directory.mkdir(parents=True)
                        new_directories.append((source_directory, destination_directory))
                for source in extracted_files:
                    relative = source.relative_to(extracted)
                    destination = raw_path / relative
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    if destination.is_symlink() or (destination.exists() and not destination.is_file()):
                        raise SpeciesTreeArchiveError(f"Live path conflicts with an archived file: {destination}")
                    if not destination.exists():
                        os.replace(source, destination)
                for source_directory, destination_directory in reversed(new_directories):
                    shutil.copystat(source_directory, destination_directory)
            else:
                os.replace(extracted, raw_path)
            _fsync_directory(raw_path)
            archive_path.unlink()
            _fsync_directory(root)
        except OSError as exc:
            raise SpeciesTreeArchiveError(f"Failed to materialize {archive_path}: {exc}") from exc
        finally:
            shutil.rmtree(temporary, ignore_errors=True)
    result = {"directory": name, "state": "raw", "files": file_count}
    if include_metrics:
        result.update(
            {
                "required_peak_new_bytes": required_bytes,
                "required_inodes": required_inodes,
                "filesystem_free_bytes": filesystem_free_bytes,
                "effective_available_bytes": effective_available_bytes,
            }
        )
    return result


def status(root: Path, names: Iterable[str]) -> list[dict]:
    records: list[dict] = []
    with _store_lock(root):
        for name in names:
            raw_path = _raw_path(root, name)
            archive_path = _archive_path(root, name)
            _assert_regular_path(raw_path, "managed directory")
            _assert_regular_path(archive_path, "ZIP archive")
            if raw_path.exists() and not raw_path.is_dir():
                raise SpeciesTreeArchiveError(f"Managed path is not a directory: {raw_path}")
            if raw_path.is_dir():
                raw_inventory = _iter_source_files(raw_path)
                raw_files = len(raw_inventory)
                raw_bytes = sum(int(path.stat().st_size) for path in raw_inventory)
                estimated_peak_new_bytes = _raw_to_zip_peak_bytes(raw_inventory)
            else:
                raw_files = 0
                raw_bytes = 0
                estimated_peak_new_bytes = 0
            archived_files = 0
            archived_logical_bytes = 0
            zip_bytes = 0
            raw_materialize_allocated_bytes = 0
            raw_materialize_required_inodes = 0
            if archive_path.exists():
                archive_status = verify_archive(root, name, check_crc=False)
                archived_files = int(archive_status["files"])
                archived_logical_bytes = int(archive_status.get("logical_bytes", 0))
                zip_bytes = int(archive_status.get("zip_bytes", 0))
                raw_materialize_allocated_bytes = int(archive_status.get("raw_materialize_allocated_bytes", 0))
                raw_materialize_required_inodes = int(archive_status.get("raw_materialize_required_inodes", 0))
            state = (
                "mixed"
                if raw_path.exists() and archive_path.exists()
                else "raw"
                if raw_path.exists()
                else "archived"
                if archive_path.exists()
                else "absent"
            )
            records.append(
                {
                    "directory": name,
                    "state": state,
                    "raw_files": raw_files,
                    "raw_bytes": raw_bytes,
                    "archived_files": archived_files,
                    "archived_logical_bytes": archived_logical_bytes,
                    "zip_bytes": zip_bytes,
                    "estimated_peak_new_bytes": estimated_peak_new_bytes,
                    "raw_materialize_allocated_bytes": raw_materialize_allocated_bytes,
                    "raw_materialize_required_inodes": raw_materialize_required_inodes,
                }
            )
    return records


def count_matching_files(root: Path, name: str, pattern: str) -> int:
    """Count non-hidden, top-level files across raw and archived forms."""
    name = _directory_name(name)
    raw_path = _raw_path(root, name)
    archive_path = _archive_path(root, name)
    matching_names: set[str] = set()
    with _store_lock(root):
        _assert_regular_path(raw_path, "managed directory")
        _assert_regular_path(archive_path, "ZIP archive")
        if raw_path.exists():
            if not raw_path.is_dir():
                raise SpeciesTreeArchiveError(f"Managed path is not a directory: {raw_path}")
            _, raw_files = _source_inventory(raw_path)
            for path in raw_files:
                relative = path.relative_to(raw_path)
                if (
                    len(relative.parts) == 1
                    and not relative.name.startswith(".")
                    and fnmatch.fnmatchcase(relative.name, pattern)
                ):
                    matching_names.add(relative.name)
        if archive_path.exists():
            if not archive_path.is_file():
                raise SpeciesTreeArchiveError(f"Archive path is not a file: {archive_path}")
            try:
                with zipfile.ZipFile(archive_path, "r") as archive:
                    members, _ = _validated_members(archive, archive_path, name, check_crc=False)
                    for info in members:
                        relative = PurePosixPath(info.filename).relative_to(name)
                        if (
                            len(relative.parts) == 1
                            and not relative.name.startswith(".")
                            and fnmatch.fnmatchcase(relative.name, pattern)
                        ):
                            matching_names.add(relative.name)
            except (OSError, RuntimeError, zipfile.BadZipFile) as exc:
                raise SpeciesTreeArchiveError(f"Failed to inspect {archive_path}: {exc}") from exc
    return len(matching_names)


def _validate_options(compression: str, compression_level: int) -> None:
    if compression not in {"adaptive", "deflate", "store"}:
        raise SpeciesTreeArchiveError(f"Invalid compression {compression!r}; expected adaptive, deflate, or store")
    if compression_level < 0 or compression_level > 9:
        raise SpeciesTreeArchiveError("Compression level must be from 0 through 9")


def _emit(records: object) -> None:
    print(json.dumps(records, indent=2, sort_keys=True))


def _run_each(names: Iterable[str], operation: Callable[[str], dict]) -> tuple[list[dict], list[str]]:
    records: list[dict] = []
    errors: list[str] = []
    for name in names:
        try:
            records.append(operation(name))
        except SpeciesTreeArchiveError as exc:
            errors.append(f"{name}: {exc}")
    return records, errors


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add_common(subparser: argparse.ArgumentParser) -> None:
        subparser.add_argument("--root", required=True, type=Path)
        subparser.add_argument("--directory", action="append", choices=MANAGED_DIRECTORIES)
        subparser.add_argument(
            "--json",
            action="store_true",
            help="Emit JSON (species-tree manager output is JSON by default).",
        )

    pack = subparsers.add_parser("pack", help="Atomically replace raw directories with ZIPs")
    add_common(pack)
    pack.add_argument("--compression", choices=("adaptive", "deflate", "store"), default="adaptive")
    pack.add_argument("--compression-level", type=int, default=6)
    pack.add_argument(
        "--available-bytes",
        type=int,
        help="Known quota headroom, combined with filesystem free space",
    )
    pack.add_argument("--progress-interval", type=float, default=30.0)

    materialize = subparsers.add_parser("materialize", help="Restore ZIPs as raw directories and consume the ZIPs")
    add_common(materialize)
    materialize.add_argument(
        "--available-bytes",
        type=int,
        help="Known quota headroom, combined with filesystem free space",
    )
    materialize.add_argument("--progress-interval", type=float, default=30.0)

    convert = subparsers.add_parser("convert-storage", help="Convert managed directories")
    add_common(convert)
    convert.add_argument("--to", required=True, choices=("zip", "raw", "files"))
    convert.add_argument("--compression", choices=("adaptive", "deflate", "store"), default="adaptive")
    convert.add_argument("--compression-level", type=int, default=6)
    convert.add_argument(
        "--available-bytes",
        type=int,
        help="Known quota headroom, combined with filesystem free space",
    )
    convert.add_argument("--dry-run", action="store_true")
    convert.add_argument("--progress-interval", type=float, default=30.0)

    verify = subparsers.add_parser("verify", help="Validate managed archives")
    add_common(verify)
    verify.add_argument("--progress-interval", type=float, default=30.0)
    verify_depth = verify.add_mutually_exclusive_group()
    verify_depth.add_argument(
        "--quick",
        action="store_true",
        help="Validate metadata and ZIP member inventories without reading payload bytes.",
    )
    verify_depth.add_argument(
        "--deep",
        action="store_true",
        help="Read every member and validate CRC (the default).",
    )

    status_parser = subparsers.add_parser("status", help="Show raw/ZIP state")
    add_common(status_parser)
    count = subparsers.add_parser("count", help="Count matching top-level files without materializing a ZIP")
    count.add_argument("--root", required=True, type=Path)
    count.add_argument("--directory", required=True, choices=MANAGED_DIRECTORIES)
    count.add_argument("--pattern", required=True)
    return parser


def run_cli(args: argparse.Namespace) -> int:
    root = args.root.resolve()
    if args.command == "count":
        print(count_matching_files(root, args.directory, args.pattern))
        return 0
    names = _selected_directories(args.directory)
    if args.command in {"pack", "convert-storage", "materialize"}:
        if args.command != "materialize":
            _validate_options(args.compression, args.compression_level)
        if args.available_bytes is not None and args.available_bytes < 0:
            raise SpeciesTreeArchiveError("--available-bytes must be non-negative")
        if args.progress_interval < 0:
            raise SpeciesTreeArchiveError("--progress-interval must be non-negative")
    if args.command == "verify" and args.progress_interval < 0:
        raise SpeciesTreeArchiveError("--progress-interval must be non-negative")
    if args.command == "materialize" and args.progress_interval < 0:
        raise SpeciesTreeArchiveError("--progress-interval must be non-negative")
    if args.command == "convert-storage" and args.dry_run:
        records = status(root, names)
        filesystem_free_bytes = int(shutil.disk_usage(root).free)
        effective_available_bytes = filesystem_free_bytes
        if args.available_bytes is not None:
            effective_available_bytes = min(
                effective_available_bytes,
                int(args.available_bytes),
            )
        for record in records:
            requested_target = "raw" if args.to == "files" else args.to
            record["requested_target"] = requested_target
            record["filesystem_free_bytes"] = filesystem_free_bytes
            record["effective_available_bytes"] = effective_available_bytes
            required_peak_new_bytes = int(
                record["estimated_peak_new_bytes" if requested_target == "zip" else "raw_materialize_allocated_bytes"]
            )
            record["required_peak_new_bytes"] = required_peak_new_bytes
            record["temporary_space_sufficient"] = required_peak_new_bytes <= effective_available_bytes
        _emit(records)
        return 0
    if args.command == "pack" or (args.command == "convert-storage" and args.to == "zip"):
        reporter = ProgressReporter(args.progress_interval)
        records, errors = _run_each(
            names,
            lambda name: pack_directory(
                root,
                name,
                compression=args.compression,
                compression_level=args.compression_level,
                available_bytes=args.available_bytes,
                progress_callback=reporter.update,
                include_metrics=True,
            ),
        )
        _emit(records)
        for error in errors:
            print(f"Error: {error}", file=sys.stderr)
        return 1 if errors else 0
    if args.command == "materialize" or (args.command == "convert-storage" and args.to in {"raw", "files"}):
        reporter = ProgressReporter(args.progress_interval)
        records, errors = _run_each(
            names,
            lambda name: materialize_directory(
                root,
                name,
                available_bytes=args.available_bytes,
                progress_callback=reporter.update,
                include_metrics=True,
            ),
        )
        _emit(records)
        for error in errors:
            print(f"Error: {error}", file=sys.stderr)
        return 1 if errors else 0
    if args.command == "verify":
        reporter = ProgressReporter(args.progress_interval)
        records, errors = _run_each(
            names,
            lambda name: verify_archive(
                root,
                name,
                check_crc=not args.quick,
                progress_callback=(None if args.quick else reporter.update),
            ),
        )
        _emit(records)
        for error in errors:
            print(f"Error: {error}", file=sys.stderr)
        return 1 if errors else 0
    if args.command == "status":
        _emit(status(root, names))
        return 0
    raise SpeciesTreeArchiveError(f"Unsupported command: {args.command}")


def main(argv: Sequence[str] | None = None) -> int:
    try:
        return run_cli(build_parser().parse_args(argv))
    except SpeciesTreeArchiveError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    except (OSError, RuntimeError, zipfile.BadZipFile) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

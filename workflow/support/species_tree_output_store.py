#!/usr/bin/env python3
"""Manage high-file-count species-tree stage directories as visible ZIP files."""

from __future__ import annotations

import argparse
import contextlib
import fcntl
import fnmatch
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


def _directory_name(value: str) -> str:
    if value not in MANAGED_DIRECTORIES:
        expected = ", ".join(MANAGED_DIRECTORIES)
        raise SpeciesTreeArchiveError(
            f"Unsupported species-tree directory {value!r}; expected one of: {expected}"
        )
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

    for current_root, directory_names, file_names in os.walk(
        raw_path, onerror=raise_walk_error
    ):
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

    for current_root, directory_names, file_names in os.walk(
        raw_path, onerror=raise_walk_error
    ):
        current = Path(current_root)
        for directory_name in directory_names:
            _assert_regular_path(current / directory_name, "source directory")
        for file_name in file_names:
            path = current / file_name
            _assert_regular_path(path, "source file")
            if path.exists() and not path.is_file():
                raise SpeciesTreeArchiveError(
                    f"Non-regular source file is not supported: {path}"
                )
            count += 1
    return count


def _signature(path: Path) -> tuple[int, int, int, int]:
    value = path.stat()
    return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns


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
        raise SpeciesTreeArchiveError(
            f"Refusing to manage a ZIP without a GeneGalleon species-tree marker: {path}"
        )
    try:
        payload = json.loads(archive.comment[len(COMMENT_PREFIX) :].decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SpeciesTreeArchiveError(f"Invalid archive metadata in {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise SpeciesTreeArchiveError(
            f"Archive metadata must be a JSON object in {path}"
        )
    if (
        type(payload.get("schema_version")) is not int
        or payload.get("schema_version") != SCHEMA_VERSION
    ):
        raise SpeciesTreeArchiveError(
            f"Unsupported species-tree archive schema in {path}: "
            f"{payload.get('schema_version')!r}"
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
        raise SpeciesTreeArchiveError(
            f"Archive directory mismatch in {path}: {payload.get('directory')!r}"
        )
    file_members: list[zipfile.ZipInfo] = []
    directory_members: list[zipfile.ZipInfo] = []
    seen: set[str] = set()
    file_names: set[str] = set()
    prefix = f"{expected_name}/"
    for info in archive.infolist():
        member = PurePosixPath(info.filename)
        canonical_name = member.as_posix()
        canonical_archive_name = (
            f"{canonical_name}/" if info.is_dir() else canonical_name
        )
        root_directory_member = info.is_dir() and canonical_name == expected_name
        if (
            canonical_name in seen
            or info.filename.startswith("/")
            or "\\" in info.filename
            or canonical_archive_name != info.filename
            or ".." in member.parts
            or (
                not root_directory_member
                and (
                    not info.filename.startswith(prefix)
                    or len(member.parts) < 2
                )
            )
        ):
            raise SpeciesTreeArchiveError(f"Unsafe ZIP member in {path}: {info.filename!r}")
        seen.add(canonical_name)
        mode = info.external_attr >> 16
        if stat.S_ISLNK(mode):
            raise SpeciesTreeArchiveError(
                f"Symlinked ZIP members are not supported in {path}: {info.filename}"
            )
        if info.is_dir():
            directory_members.append(info)
        else:
            file_members.append(info)
            file_names.add(canonical_name)
    for canonical_name in seen:
        parent = PurePosixPath(canonical_name).parent
        while parent.as_posix() not in {".", expected_name}:
            if parent.as_posix() in file_names:
                raise SpeciesTreeArchiveError(
                    f"ZIP member is nested below a file in {path}: {canonical_name!r}"
                )
            parent = parent.parent
    expected_count = payload.get("file_count")
    if type(expected_count) is not int or expected_count != len(file_members):
        raise SpeciesTreeArchiveError(
            f"Archive file count differs from its metadata in {path}: "
            f"metadata={expected_count!r}, actual={len(file_members)}"
        )
    expected_directory_count = payload.get("directory_count")
    if expected_directory_count is not None and (
        type(expected_directory_count) is not int
        or expected_directory_count != len(directory_members)
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


def verify_archive(root: Path, name: str, *, check_crc: bool = True) -> dict:
    path = _archive_path(root, name)
    _assert_regular_path(path, "ZIP archive")
    if not path.exists():
        return {"directory": name, "state": "absent", "files": 0}
    if not path.is_file():
        raise SpeciesTreeArchiveError(f"Archive path is not a file: {path}")
    try:
        with zipfile.ZipFile(path, "r") as archive:
            members, _ = _validated_members(archive, path, name, check_crc=check_crc)
    except (OSError, RuntimeError, zipfile.BadZipFile) as exc:
        raise SpeciesTreeArchiveError(f"Failed to verify {path}: {exc}") from exc
    return {"directory": name, "state": "archived", "files": len(members)}


def pack_directory(
    root: Path,
    name: str,
    *,
    compression: str = "adaptive",
    compression_level: int = 6,
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
        signatures = {path: _signature(path) for path in directories + files}
        partial_path = root / f".{name}.zip.partial.{os.getpid()}.{uuid.uuid4().hex}"
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
                    info = zipfile.ZipInfo.from_file(
                        path, arcname=archive_name, strict_timestamps=False
                    )
                    info.compress_type = zipfile.ZIP_STORED
                    archive.writestr(info, b"")
                for path in files:
                    relative = path.relative_to(raw_path).as_posix()
                    info = zipfile.ZipInfo.from_file(
                        path,
                        arcname=f"{name}/{relative}",
                        strict_timestamps=False,
                    )
                    info.compress_type = _compression_for(path, compression)
                    with path.open("rb") as source, archive.open(
                        info, "w", force_zip64=True
                    ) as destination:
                        shutil.copyfileobj(source, destination, length=1024 * 1024)
                archive.comment = _archive_comment(name, len(files), len(directories))

            with zipfile.ZipFile(partial_path, "r") as archive:
                _validated_members(archive, partial_path, name)
            for path, signature in signatures.items():
                if (
                    not path.exists()
                    or path.is_symlink()
                    or _signature(path) != signature
                ):
                    raise SpeciesTreeArchiveError(
                        f"Source changed while {name} was being archived: {path}"
                    )
            current_directories, current_files = _source_inventory(raw_path)
            if current_directories != directories or current_files != files:
                raise SpeciesTreeArchiveError(
                    f"Source inventory changed while {name} was being archived"
                )

            with partial_path.open("rb") as handle:
                os.fsync(handle.fileno())
            os.replace(partial_path, archive_path)
            _fsync_directory(root)
            shutil.rmtree(raw_path)
            _fsync_directory(root)
        finally:
            partial_path.unlink(missing_ok=True)
    return {"directory": name, "state": "archived", "files": len(files)}


def _safe_destination(base: Path, member_name: str) -> Path:
    destination = base.joinpath(*PurePosixPath(member_name).parts)
    try:
        destination.resolve().relative_to(base.resolve())
    except ValueError as exc:
        raise SpeciesTreeArchiveError(
            f"ZIP member escapes its extraction root: {member_name!r}"
        ) from exc
    return destination


def _extract_archive_to(archive_path: Path, name: str, destination: Path) -> int:
    try:
        with zipfile.ZipFile(archive_path, "r") as archive:
            members, directory_members = _validated_members(
                archive, archive_path, name, check_crc=False
            )
            directory_metadata: list[tuple[Path, int, tuple[int, ...]]] = []
            for info in sorted(
                directory_members,
                key=lambda item: len(PurePosixPath(item.filename).parts),
            ):
                target = _safe_destination(destination, info.filename)
                target.mkdir(parents=True, exist_ok=True)
                directory_metadata.append(
                    (target, info.external_attr >> 16, info.date_time)
                )
            for info in members:
                target = _safe_destination(destination, info.filename)
                target.parent.mkdir(parents=True, exist_ok=True)
                with archive.open(info, "r") as source, target.open("wb") as output:
                    shutil.copyfileobj(source, output, length=1024 * 1024)
                    output.flush()
                    os.fsync(output.fileno())
                mode = info.external_attr >> 16
                if mode:
                    target.chmod(stat.S_IMODE(mode))
                modified = time.mktime((*info.date_time, 0, 0, -1))
                os.utime(target, (modified, modified))
            for target, mode, date_time in reversed(directory_metadata):
                if mode:
                    target.chmod(stat.S_IMODE(mode))
                modified = time.mktime((*date_time, 0, 0, -1))
                os.utime(target, (modified, modified))
            return len(members)
    except (OSError, RuntimeError, zipfile.BadZipFile) as exc:
        raise SpeciesTreeArchiveError(f"Failed to extract {archive_path}: {exc}") from exc


def materialize_directory(root: Path, name: str) -> dict:
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

        temporary = root / f".{name}.materialize.partial.{os.getpid()}.{uuid.uuid4().hex}"
        temporary.mkdir()
        try:
            (temporary / name).mkdir()
            file_count = _extract_archive_to(archive_path, name, temporary)
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
                        if (
                            destination_directory.is_symlink()
                            or not destination_directory.is_dir()
                        ):
                            raise SpeciesTreeArchiveError(
                                "Live path conflicts with an archived directory: "
                                f"{destination_directory}"
                            )
                    else:
                        destination_directory.mkdir(parents=True)
                        new_directories.append(
                            (source_directory, destination_directory)
                        )
                for source in extracted_files:
                    relative = source.relative_to(extracted)
                    destination = raw_path / relative
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    if destination.is_symlink() or (
                        destination.exists() and not destination.is_file()
                    ):
                        raise SpeciesTreeArchiveError(
                            f"Live path conflicts with an archived file: {destination}"
                        )
                    if not destination.exists():
                        os.replace(source, destination)
                for source_directory, destination_directory in reversed(
                    new_directories
                ):
                    shutil.copystat(source_directory, destination_directory)
            else:
                os.replace(extracted, raw_path)
            _fsync_directory(raw_path)
            archive_path.unlink()
            _fsync_directory(root)
        except OSError as exc:
            raise SpeciesTreeArchiveError(
                f"Failed to materialize {archive_path}: {exc}"
            ) from exc
        finally:
            shutil.rmtree(temporary, ignore_errors=True)
    return {"directory": name, "state": "raw", "files": file_count}


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
            raw_files = _count_source_files_best_effort(raw_path) if raw_path.is_dir() else 0
            archived_files = 0
            if archive_path.exists():
                archived_files = int(
                    verify_archive(root, name, check_crc=False)["files"]
                )
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
                    "archived_files": archived_files,
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
                    members, _ = _validated_members(
                        archive, archive_path, name, check_crc=False
                    )
                    for info in members:
                        relative = PurePosixPath(info.filename).relative_to(name)
                        if (
                            len(relative.parts) == 1
                            and not relative.name.startswith(".")
                            and fnmatch.fnmatchcase(relative.name, pattern)
                        ):
                            matching_names.add(relative.name)
            except (OSError, RuntimeError, zipfile.BadZipFile) as exc:
                raise SpeciesTreeArchiveError(
                    f"Failed to inspect {archive_path}: {exc}"
                ) from exc
    return len(matching_names)


def _validate_options(compression: str, compression_level: int) -> None:
    if compression not in {"adaptive", "deflate", "store"}:
        raise SpeciesTreeArchiveError(
            f"Invalid compression {compression!r}; expected adaptive, deflate, or store"
        )
    if compression_level < 0 or compression_level > 9:
        raise SpeciesTreeArchiveError("Compression level must be from 0 through 9")


def _emit(records: object) -> None:
    print(json.dumps(records, indent=2, sort_keys=True))


def _run_each(
    names: Iterable[str], operation: Callable[[str], dict]
) -> tuple[list[dict], list[str]]:
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

    pack = subparsers.add_parser("pack", help="Atomically replace raw directories with ZIPs")
    add_common(pack)
    pack.add_argument("--compression", choices=("adaptive", "deflate", "store"), default="adaptive")
    pack.add_argument("--compression-level", type=int, default=6)

    materialize = subparsers.add_parser(
        "materialize", help="Restore ZIPs as raw directories and consume the ZIPs"
    )
    add_common(materialize)

    convert = subparsers.add_parser("convert-storage", help="Convert managed directories")
    add_common(convert)
    convert.add_argument("--to", required=True, choices=("zip", "raw", "files"))
    convert.add_argument("--compression", choices=("adaptive", "deflate", "store"), default="adaptive")
    convert.add_argument("--compression-level", type=int, default=6)

    verify = subparsers.add_parser("verify", help="CRC-check managed archives")
    add_common(verify)

    status_parser = subparsers.add_parser("status", help="Show raw/ZIP state")
    add_common(status_parser)
    count = subparsers.add_parser(
        "count", help="Count matching top-level files without materializing a ZIP"
    )
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
    if args.command in {"pack", "convert-storage"}:
        _validate_options(args.compression, args.compression_level)
    if args.command == "pack" or (
        args.command == "convert-storage" and args.to == "zip"
    ):
        records, errors = _run_each(
            names,
            lambda name: pack_directory(
                root,
                name,
                compression=args.compression,
                compression_level=args.compression_level,
            ),
        )
        _emit(records)
        for error in errors:
            print(f"Error: {error}", file=sys.stderr)
        return 1 if errors else 0
    if args.command == "materialize" or (
        args.command == "convert-storage" and args.to in {"raw", "files"}
    ):
        records, errors = _run_each(
            names, lambda name: materialize_directory(root, name)
        )
        _emit(records)
        for error in errors:
            print(f"Error: {error}", file=sys.stderr)
        return 1 if errors else 0
    if args.command == "verify":
        records, errors = _run_each(names, lambda name: verify_archive(root, name))
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

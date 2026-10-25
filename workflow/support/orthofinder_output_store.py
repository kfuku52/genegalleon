#!/usr/bin/env python3
"""Lossless ZIP64 storage for completed OrthoFinder result directories.

Shared summary tables stay raw. Every OG (including unselected OGs) belongs to this store.
No WorkingDirectory or native restart checkpoint is inferred or manufactured.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import os
import shutil
import socket
import stat
import tempfile
import time
import uuid
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path, PurePosixPath

try:
    from . import species_tree_output_store as DIRECTORY
    from .directory_storage_policy import ALL_DIRECTORIES, DIRECTORY_SCOPES, ORTHOFINDER_DIRECTORIES
    from .shared_namespace_lock import inspect_lock, namespace_lock
except ImportError:
    import species_tree_output_store as DIRECTORY
    from directory_storage_policy import ALL_DIRECTORIES, DIRECTORY_SCOPES, ORTHOFINDER_DIRECTORIES
    from shared_namespace_lock import inspect_lock, namespace_lock

MANAGED_DIRECTORIES = ORTHOFINDER_DIRECTORIES
MANIFEST = "__gg_orthofinder_manifest__.json"
MARKER = b"GeneGalleon orthofinder directory archive v1\n"
DIRECTORY_MARKER = b"GeneGalleon retained directory archive v1\n"
DIRECTORY_MANIFEST = "__gg_directory_manifest__.json"
ARCHIVE_BASENAMES = {Path(name).name + ".zip" for name in ALL_DIRECTORIES}
CHUNK = 1024 * 1024


def _format(name):
    return (MARKER, MANIFEST) if name in MANAGED_DIRECTORIES else (DIRECTORY_MARKER, DIRECTORY_MANIFEST)


def _namespace_directories(namespace):
    _require(namespace == "orthofinder" or namespace in DIRECTORY_SCOPES, "Invalid storage namespace")
    return MANAGED_DIRECTORIES if namespace == "orthofinder" else DIRECTORY_SCOPES[namespace]


class OrthoFinderStorageError(RuntimeError):
    pass


def _require(condition, message):
    if not condition:
        raise OrthoFinderStorageError(message)


def _physical(path):
    _require(path.resolve() == path.absolute(), f"Redirected storage path: {path}")
    _require(not path.is_symlink(), f"Symlinked storage path: {path}")


def lock_path(root, namespace="orthofinder"):
    if namespace == "orthofinder":
        return root.parent / ".gg_orthofinder_storage.lock"
    _require(namespace in DIRECTORY_SCOPES, "Invalid storage namespace")
    return root.parent / (".gg_" + root.name + "_directory_storage.lock")


@contextlib.contextmanager
def storage_lock(root, *, exclusive=True, producer_lock_token=None, namespace="orthofinder"):
    _physical(root)
    if producer_lock_token is not None:
        # A native producer holds the fence through its child converter. The
        # child may use, but must never release, its live parent's exact lock.
        owner = inspect_lock(lock_path(root, namespace))["exclusive"]
        _require(
            exclusive and owner is not None and owner.get("token") == producer_lock_token
            and owner.get("pid") == os.getppid() and owner.get("host") == socket.gethostname(),
            "Invalid parent producer lock ownership",
        )
        yield
        return
    with namespace_lock(lock_path(root, namespace), exclusive=exclusive, timeout=300):
        yield


def _sync_directory(path):
    """A failed durability barrier must prevent deletion of the recovery form."""
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


class _StreamingWriter:
    """ZIP data descriptors avoid a flush/seek/header rewrite per small file."""

    def __init__(self, handle):
        self.handle = handle

    def write(self, data):
        return self.handle.write(data)

    def flush(self):
        return self.handle.flush()


def _paths(root, name):
    _require(name in ALL_DIRECTORIES, f"Unmanaged output directory: {name}")
    raw = root / name
    archive = raw.with_name(raw.name + ".zip")
    _physical(raw)
    _physical(archive)
    return raw, archive


def _signature(value):
    return (value.st_dev, value.st_ino, value.st_mode, value.st_size, value.st_mtime_ns, value.st_ctime_ns)


def _inventory(raw):
    _physical(raw)
    directories, files = [], []

    def scan(path):
        value = path.lstat()
        _require(stat.S_ISDIR(value.st_mode), f"Non-directory source: {path}")
        directories.append((path, value))
        with os.scandir(path) as entries:
            for entry in entries:
                child, metadata = Path(entry.path), entry.stat(follow_symlinks=False)
                if stat.S_ISDIR(metadata.st_mode):
                    scan(child)
                else:
                    _require(stat.S_ISREG(metadata.st_mode), f"Non-regular source: {child}")
                    files.append((child, metadata))

    scan(raw)
    return (sorted(directories, key=lambda row: row[0].as_posix()), sorted(files, key=lambda row: row[0].as_posix()))


def _hash(handle):
    digest, size = hashlib.sha256(), 0
    for chunk in iter(lambda: handle.read(CHUNK), b""):
        digest.update(chunk)
        size += len(chunk)
    return digest.hexdigest(), size


@contextlib.contextmanager
def _source(path, expected):
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    with os.fdopen(fd, "rb") as handle:
        _require(_signature(os.fstat(handle.fileno())) == _signature(expected), f"Source changed: {path}")
        yield handle
        _require(_signature(os.fstat(handle.fileno())) == _signature(expected), f"Source changed: {path}")


def _metadata(path, value, raw):
    return dict(path=path.relative_to(raw).as_posix(), mode=stat.S_IMODE(value.st_mode), mtime_ns=value.st_mtime_ns)


def require_completed(root):
    _physical(root)
    log = root / "Log.txt"
    _physical(log)
    _require(
        log.is_file()
        and any(line.rstrip().endswith("OrthoFinder run completed") for line in log.read_text().splitlines()),
        f"A completed all-species OrthoFinder result is required: {log}",
    )


def _manifest(archive, path, name):
    marker, manifest_name = _format(name)
    _require(archive.comment == marker + name.encode(), f"Unmarked directory archive: {path}")
    infos = archive.infolist()
    _require(len({info.filename for info in infos}) == len(infos), f"Duplicate ZIP members: {path}")
    member = archive.getinfo(manifest_name)
    _require(member.file_size <= 256 * 1024 * 1024, f"Oversized manifest: {path}")
    data = json.loads(archive.read(member))
    _require(
        isinstance(data, dict) and type(data.get("schema_version")) is int
        and data["schema_version"] == 1 and data.get("directory") == name
        and isinstance(data.get("files"), list) and isinstance(data.get("directories"), list)
        and all(isinstance(row, dict) for row in data["files"] + data["directories"]),
        f"Invalid OrthoFinder archive manifest: {path}",
    )
    records = data["files"]
    for row in records + data["directories"]:
        _require(
            type(row.get("mode")) is int and 0 <= row["mode"] <= 0o7777 and type(row.get("mtime_ns")) is int,
            "Invalid member metadata",
        )
    expected = {manifest_name}
    prefix = Path(name).name
    for record in records:
        relative = PurePosixPath(record["path"])
        _require(
            relative.parts
            and not relative.is_absolute()
            and ".." not in relative.parts
            and relative.as_posix() == record["path"]
            and "\\" not in record["path"],
            f"Unsafe manifest path: {record['path']}",
        )
        _require(
            type(record["size"]) is int
            and record["size"] >= 0
            and isinstance(record["sha256"], str)
            and len(record["sha256"]) == 64,
            "Invalid manifest content",
        )
        bytes.fromhex(record["sha256"])
        full = prefix + "/" + record["path"]
        _require(full not in expected, f"Duplicate manifest member: {full}")
        expected.add(full)
        info = archive.getinfo(full)
        _require(
            not info.is_dir()
            and stat.S_IFMT(info.external_attr >> 16) in (0, stat.S_IFREG)
            and info.file_size == record["size"],
            f"Invalid ZIP member: {full}",
        )
    directory_names = set()
    for record in data["directories"]:
        relative = PurePosixPath(record["path"])
        _require(
            not relative.is_absolute()
            and ".." not in relative.parts
            and relative.as_posix() == record["path"]
            and "\\" not in record["path"],
            "Unsafe directory metadata",
        )
        _require(record["path"] not in directory_names, "Duplicate directory metadata")
        directory_names.add(record["path"])
    _require("." in directory_names, "Missing root directory metadata")
    files = {record["path"] for record in records}
    for relative in [PurePosixPath(value) for value in files | directory_names]:
        _require(
            relative.as_posix() not in files or relative.as_posix() not in directory_names, "File/directory collision"
        )
        for parent in relative.parents:
            _require(parent.as_posix() in directory_names, "Missing parent directory metadata")
    _require({info.filename for info in infos} == expected, f"Unexpected ZIP members: {path}")
    return data


def archive_records(root, name, *, deep=False, progress_callback=None):
    _, path = _paths(root, name)
    _require(path.is_file(), f"Archive is not a regular file: {path}")
    before = _signature(path.stat())
    try:
        with zipfile.ZipFile(path) as archive:
            data = _manifest(archive, path, name)
            for index, record in enumerate(data["files"]):
                if deep:
                    with archive.open(Path(name).name + "/" + record["path"]) as handle:
                        actual = _hash(handle)
                    _require(
                        actual == (record["sha256"], record["size"]), f"SHA256 mismatch in {path}: {record['path']}"
                    )
                if progress_callback is not None:
                    progress_callback(
                        phase="verifying", directory=name, files_completed=index + 1, files_total=len(data["files"])
                    )
    except (KeyError, ValueError, TypeError, OSError, zipfile.BadZipFile) as exc:
        raise OrthoFinderStorageError(f"Invalid OrthoFinder archive {path}: {exc}") from exc
    _require(before == _signature(path.stat()), f"Archive changed during read: {path}")
    return data


def _retire(raw, data, *, recovery=None):
    """Resume a partial deletion only if every remaining file matches the ZIP."""
    directories, files = _inventory(raw)
    records = {row["path"]: row for row in data["files"]}
    known_directories = {row["path"]: row for row in data["directories"]}
    for path, value in directories:
        record = known_directories.get(path.relative_to(raw).as_posix())
        _require(record is not None, f"New source directory: {path}")
        # Directory mtimes change during partial retirement; modes do not.
        _require(stat.S_IMODE(value.st_mode) == record["mode"], f"Source directory metadata differs: {path}")
    for path, value in files:
        record = records.get(path.relative_to(raw).as_posix())
        _require(record is not None, f"Unarchived source file: {path}")
        _require(
            stat.S_IMODE(value.st_mode) == record["mode"] and value.st_mtime_ns == record["mtime_ns"],
            f"Source metadata differs from ZIP: {path}",
        )
        with _source(path, value) as handle:
            _require(_hash(handle) == (record["sha256"], record["size"]), f"Source differs from ZIP: {path}")
    if recovery is not None:
        archive, expected_signature = recovery
        _require(_signature(archive.lstat()) == expected_signature, f"Recovery archive changed: {archive}")
    for path, value in files:
        _require(_signature(path.lstat()) == _signature(value), f"Source changed before deletion: {path}")
        path.unlink()
    for path, value in sorted(directories, key=lambda row: len(row[0].parts), reverse=True):
        _require(stat.S_IMODE(path.lstat().st_mode) == stat.S_IMODE(value.st_mode), f"Source metadata changed: {path}")
        path.rmdir()  # A new file is never recursively deleted.


def _pack(
    root, name, *, compression="adaptive", compression_level=6, available_bytes=None, progress_callback=None, **_unused
):
    raw, path = _paths(root, name)
    if path.exists():
        archive_signature = _signature(path.stat())
        data = archive_records(root, name, deep=raw.exists(), progress_callback=progress_callback)
        if raw.exists():
            with path.open("rb") as handle:
                os.fsync(handle.fileno())
            _sync_directory(path.parent)
            _retire(raw, data, recovery=(path, archive_signature))
            _sync_directory(path.parent)
        return dict(directory=name, state="archived", files=len(data["files"]))
    if not raw.exists():
        return dict(directory=name, state="absent", files=0)
    directories, files = _inventory(raw)
    peak = sum(value.st_size + value.st_size // 1000 + 2048 for _, value in files) + 1024 * 1024
    _space(root, peak, available_bytes)
    partial = path.with_name("." + path.name + ".partial." + uuid.uuid4().hex)
    data = dict(schema_version=1, directory=name, files=[], directories=[_metadata(p, s, raw) for p, s in directories])
    DIRECTORY._validate_options(compression, compression_level)
    try:
        with partial.open("xb", buffering=CHUNK) as output, zipfile.ZipFile(
            _StreamingWriter(output), "w", allowZip64=True
        ) as archive:
            for index, (source, value) in enumerate(files):
                info = zipfile.ZipInfo.from_file(
                    source, arcname=raw.name + "/" + source.relative_to(raw).as_posix(), strict_timestamps=False
                )
                info.compress_type = DIRECTORY._compression_for(source, compression)
                info._compresslevel = compression_level
                digest = hashlib.sha256()
                with _source(source, value) as handle, archive.open(info, "w", force_zip64=True) as target:
                    for chunk in iter(lambda: handle.read(CHUNK), b""):
                        target.write(chunk)
                        digest.update(chunk)
                data["files"].append(dict(_metadata(source, value, raw), size=value.st_size, sha256=digest.hexdigest()))
                if progress_callback is not None:
                    progress_callback(
                        phase="packing", directory=name, files_completed=index + 1, files_total=len(files)
                    )
            marker, manifest_name = _format(name)
            archive.writestr(manifest_name, json.dumps(data, sort_keys=True).encode(), compress_type=zipfile.ZIP_DEFLATED)
            archive.comment = marker + name.encode()
        # Reuse the exact public verifier before publication, without changing its scope.
        with zipfile.ZipFile(partial) as archive:
            _manifest(archive, partial, name)
            for record in data["files"]:
                with archive.open(raw.name + "/" + record["path"]) as handle:
                    _require(_hash(handle) == (record["sha256"], record["size"]), "New ZIP payload mismatch")
        current_directories, current_files = _inventory(raw)
        _require(
            [(p, _signature(s)) for p, s in current_directories + current_files]
            == [(p, _signature(s)) for p, s in directories + files],
            f"Source inventory changed: {raw}",
        )
        with partial.open("rb") as handle:
            os.fsync(handle.fileno())
        _require(not path.exists(), f"Archive appeared during conversion: {path}")
        os.replace(partial, path)
        _sync_directory(path.parent)
        _retire(raw, data, recovery=(path, _signature(path.stat())))
        _sync_directory(path.parent)
    finally:
        partial.unlink(missing_ok=True)
    return dict(directory=name, state="archived", files=len(files), zip_bytes=path.stat().st_size)


def _space(root, required, available):
    free = shutil.disk_usage(root).free
    if available is not None:
        _require(available >= 0, "available_bytes must be non-negative")
        free = min(free, available)
    _require(required <= free, f"Insufficient temporary space: required={required}, available={free}")


@contextlib.contextmanager
def _restore_failure_cleanup(raw, data, *, was_raw):
    try:
        yield
    except BaseException:
        if not was_raw:
            # A rejected first payload must not leave a newly visible empty
            # result tree. rmdir never removes an existing or concurrent file.
            for row in sorted(data["directories"], key=lambda row: len(PurePosixPath(row["path"]).parts), reverse=True):
                try:
                    (raw / row["path"]).rmdir()
                except OSError:
                    pass
        raise


def _materialize(root, name, *, available_bytes=None, progress_callback=None, **_unused):
    raw, path = _paths(root, name)
    if not path.exists():
        return dict(directory=name, state="raw" if raw.exists() else "absent")
    archive_signature = _signature(path.stat())
    data = archive_records(root, name, progress_callback=progress_callback)
    _require(_signature(path.stat()) == archive_signature, "Archive changed during restoration validation")
    _space(
        root,
        sum(max(4096, ((row["size"] + 4095) // 4096) * 4096) for row in data["files"])
        + len(data["directories"]) * 4096,
        available_bytes,
    )
    was_raw = raw.exists()
    if was_raw:
        # An interrupted pack may leave an equal subset; reject changed or extra files.
        existing_directories, existing = _inventory(raw)
        known_directories = {row["path"] for row in data["directories"]}
        _require(
            all(p.relative_to(raw).as_posix() in known_directories for p, _ in existing_directories),
            "Unarchived directory in raw/ZIP conflict",
        )
        records = {row["path"]: row for row in data["files"]}
        for p, value in existing:
            record = records.get(p.relative_to(raw).as_posix())
            _require(record is not None, f"Unarchived file in raw/ZIP conflict: {p}")
            _require(
                stat.S_IMODE(value.st_mode) == record["mode"] and value.st_mtime_ns == record["mtime_ns"],
                f"Raw/ZIP metadata conflict: {p}",
            )
            with _source(p, value) as handle:
                _require(_hash(handle) == (record["sha256"], record["size"]), f"Raw/ZIP conflict: {p}")
    else:
        raw.mkdir(parents=True)
    with _restore_failure_cleanup(raw, data, was_raw=was_raw):
        with zipfile.ZipFile(path) as archive:
            for record in data["directories"]:
                (raw / record["path"]).mkdir(parents=True, exist_ok=True)
            for record in data["files"]:
                target = raw / record["path"]
                if target.exists():
                    continue
                partial = raw.parent / ("." + raw.name + ".restore." + uuid.uuid4().hex)
                try:
                    with archive.open(raw.name + "/" + record["path"]) as source, partial.open("xb") as destination:
                        digest = hashlib.sha256()
                        copied = 0
                        while chunk := source.read(CHUNK):
                            destination.write(chunk)
                            digest.update(chunk)
                            copied += len(chunk)
                        _require((digest.hexdigest(), copied) == (record["sha256"], record["size"]),
                                 f"ZIP payload mismatch during restoration: {target}")
                        destination.flush()
                        os.fchmod(destination.fileno(), record["mode"])
                        os.utime(destination.fileno(), ns=(record["mtime_ns"], record["mtime_ns"]))
                    # Unlike replace(), linking a verified temporary file cannot
                    # overwrite a file another producer created after our inventory.
                    try:
                        os.link(partial, target)
                    except FileExistsError as exc:
                        raise OrthoFinderStorageError(f"Raw/ZIP publication conflict; target appeared: {target}") from exc
                except zipfile.BadZipFile as exc:
                    raise OrthoFinderStorageError(f"Invalid ZIP payload during restoration: {path}: {exc}") from exc
                finally:
                    partial.unlink(missing_ok=True)
            for record in reversed(data["directories"]):
                target = raw / record["path"]
                os.chmod(target, record["mode"])
                os.utime(target, ns=(record["mtime_ns"], record["mtime_ns"]))
        # Read every restored byte before removing its recovery ZIP.
        restored_directories, files = _inventory(raw)
        _require(
            {p.relative_to(raw).as_posix() for p, _ in restored_directories}
            == {row["path"] for row in data["directories"]},
            "Restored directory membership differs",
        )
        expected = {row["path"]: row for row in data["files"]}
        _require(len(files) == len(expected), "Restored membership differs")
        for p, value in files:
            record = expected[p.relative_to(raw).as_posix()]
            _require(
                stat.S_IMODE(value.st_mode) == record["mode"] and value.st_mtime_ns == record["mtime_ns"],
                f"Restored metadata differs: {p}",
            )
            with _source(p, value) as handle:
                _require(_hash(handle) == (record["sha256"], record["size"]), f"Restored content differs: {p}")
                # Also sync equal pre-existing files from an interrupted conversion.
                os.fsync(handle.fileno())
        for directory, _ in sorted(restored_directories, key=lambda row: len(row[0].parts), reverse=True):
            _sync_directory(directory)
        _sync_directory(raw.parent)
        final_directories, final_files = _inventory(raw)
        _require(
            [(p, _signature(s)) for p, s in final_directories + final_files]
            == [(p, _signature(s)) for p, s in restored_directories + files],
            "Restored inventory changed before commit",
        )
        _require(_signature(path.stat()) == archive_signature, "Archive changed before restoration commit")
    path.unlink()
    _sync_directory(path.parent)
    return dict(directory=name, state="raw", files=len(files))


def convert(root, *, to="zip", directories=None, workers=1, producer_lock_token=None, namespace="orthofinder", **options):
    root = Path(root).absolute()
    allowed = _namespace_directories(namespace)
    names = tuple(directories or allowed)
    _require(
        all(name in allowed for name in names) and len(names) == len(set(names)),
        "Invalid directory selection",
    )
    _require(1 <= workers <= 4 and to in ("zip", "raw"), "Invalid conversion settings")
    with storage_lock(root, producer_lock_token=producer_lock_token, namespace=namespace):
        if namespace == "orthofinder" and to == "zip" and any((root / name).exists() for name in names):
            require_completed(root)
        if workers > 1:
            # The standalone backend needs the same aggregate quota guard as
            # the workspace CLI, before any parallel worker publishes/deletes.
            rows = status(root, names)
            key = "estimated_peak_new_bytes" if to == "zip" else "raw_materialize_allocated_bytes"
            peak = sum(sorted((row[key] for row in rows), reverse=True)[:workers])
            growth = sum(max(0, row[key] - row["raw_bytes" if to == "zip" else "zip_bytes"]) for row in rows)
            _space(root, peak + growth, options.get("available_bytes"))
        operation = _pack if to == "zip" else _materialize
        with ThreadPoolExecutor(max_workers=workers) as pool:
            records = list(pool.map(lambda name: operation(root, name, **options), names))
    return {"directories": records}


def status(root, directories=None):
    records = []
    for name in directories or MANAGED_DIRECTORIES:
        raw, path = _paths(root, name)
        data = archive_records(root, name) if path.exists() else {"files": [], "directories": []}
        _, files = _inventory(raw) if raw.exists() else ([], [])
        raw_bytes = sum(value.st_size for _, value in files)
        logical_bytes = sum(row["size"] for row in data["files"])
        records.append(
            dict(
                directory=name,
                state="mixed"
                if raw.exists() and path.exists()
                else "raw"
                if raw.exists()
                else "archived"
                if path.exists()
                else "absent",
                raw_files=len(files),
                raw_bytes=raw_bytes,
                archived_files=len(data["files"]),
                archived_logical_bytes=logical_bytes,
                zip_bytes=path.stat().st_size if path.exists() else 0,
                estimated_peak_new_bytes=(sum(value.st_size + value.st_size // 1000 + 2048 for _, value in files)
                + 1024 * 1024) if raw.exists() and not path.exists() else 0,
                raw_materialize_allocated_bytes=sum(
                    max(4096, ((row["size"] + 4095) // 4096) * 4096) for row in data["files"]
                )
                + len(data["directories"]) * 4096,
            )
        )
    return records


def logical_directory_exists(raw, *, observe=lambda _path: None):
    raw = Path(raw)
    return _archive_location(raw.with_name(raw.name + ".zip"), observe=observe) is not None


def logical_file_location(path, *, observe=lambda _path: None):
    """Resolve a missing physical file inside a recognized directory archive."""
    path = Path(path).absolute()
    for parent in path.parents:
        if parent.name + ".zip" not in ARCHIVE_BASENAMES:
            continue
        location = _archive_location(parent.with_name(parent.name + ".zip"), observe=observe)
        if location is None:
            continue
        _require(not parent.exists(), f"Incomplete conversion; raw and ZIP both exist: {parent}")
        root, name = location
        relative = path.relative_to(parent).as_posix()
        records = archive_records(root, name)["files"]
        if not any(row["path"] == relative for row in records):
            return None
        namespace = "orthofinder" if name in MANAGED_DIRECTORIES else next(
            scope for scope, names in DIRECTORY_SCOPES.items() if name in names)
        return root, name, relative, namespace
    return None


def logical_directory_records(raw):
    """Read logical files, including managed ZIPs nested inside an archived directory."""
    raw = Path(raw)
    path = raw.with_name(raw.name + ".zip")
    location = _archive_location(path)
    if location is None:
        return None
    _require(not raw.exists(), f"Incomplete conversion; raw and ZIP both exist: {raw}")
    return _logical_archive_records(*location)


def _logical_members(archive, path, name):
    data = _manifest(archive, path, name)
    rows = []
    directories = {record["path"] for record in data["directories"]}
    for record in data["files"]:
        relative = PurePosixPath(record["path"])
        member = Path(name).name + "/" + record["path"]
        if relative.name not in ARCHIVE_BASENAMES:
            with archive.open(member) as handle:
                actual = _hash(handle)
            _require(actual == (record["sha256"], record["size"]), f"SHA256 mismatch in {path}: {relative}")
            rows.append(record)
            continue
        # Nested ZIPs need seeking. Spool only these members, with a bounded
        # in-memory buffer; large nested archives spill into caller-local scratch.
        with tempfile.SpooledTemporaryFile(max_size=8 * CHUNK) as spool:
            digest = hashlib.sha256()
            size = 0
            with archive.open(member) as source:
                while chunk := source.read(CHUNK):
                    spool.write(chunk)
                    digest.update(chunk)
                    size += len(chunk)
            _require((digest.hexdigest(), size) == (record["sha256"], record["size"]),
                     f"SHA256 mismatch in {path}: {relative}")
            spool.seek(0)
            try:
                child = zipfile.ZipFile(spool)
            except zipfile.BadZipFile:
                rows.append(record)  # An ordinary, unmarked .zip is a physical file.
                continue
            with child:
                child_name = _marked_name(child)
                if child_name is None:
                    rows.append(record)
                    continue
                suffix = PurePosixPath(child_name + ".zip").parts
                _require(relative.parts[-len(suffix):] == suffix, f"Misplaced nested archive: {path}:{relative}")
                prefix = relative.with_suffix("").as_posix()
                _require(prefix not in directories, f"Incomplete conversion; nested raw and ZIP both exist: {prefix}")
                rows.extend(dict(row, path=prefix + "/" + row["path"])
                            for row in _logical_members(child, f"{path}:{relative}", child_name))
    return rows


def _logical_archive_records(root, name):
    _, path = _paths(root, name)
    before = _signature(path.stat())
    try:
        with zipfile.ZipFile(path) as archive:
            rows = _logical_members(archive, path, name)
    except (KeyError, ValueError, TypeError, OSError, zipfile.BadZipFile) as exc:
        raise OrthoFinderStorageError(f"Invalid logical directory archive {path}: {exc}") from exc
    _require(_signature(path.stat()) == before, f"Archive changed during read: {path}")
    return rows


def _marked_name(archive):
    marker = next((prefix for prefix in (MARKER, DIRECTORY_MARKER) if archive.comment.startswith(prefix)), None)
    if marker is None:
        return None
    name = archive.comment[len(marker):].decode()
    _require(name in ALL_DIRECTORIES and _format(name)[0] == marker, "Invalid archive directory marker")
    return name


def _archive_location(path, *, observe=lambda _path: None):
    if path.name not in ARCHIVE_BASENAMES:
        return None
    observe(path)
    _physical(path)
    if not path.is_file():
        return None
    try:
        with zipfile.ZipFile(path) as archive:
            name = _marked_name(archive)
            if name is None:
                return None
    except UnicodeDecodeError as exc:
        raise OrthoFinderStorageError(f"Invalid archive directory marker: {path}") from exc
    except zipfile.BadZipFile:
        return None  # An unmarked physical ZIP is ordinary provenance input.
    base = path.parents[len(PurePosixPath(name).parts) - 1]
    _, expected = _paths(base, name)
    _require(path == expected, f"Misplaced OrthoFinder archive: {path}")
    return base, name


def logical_records(root, *, observe=lambda _path: None, children=None):
    """Normalize managed archives anywhere below a directory without state writes.

    Callers can share their one physical inventory. Actual ZIP bytes, rather than
    manifest claims alone, prove logical contents. Unmarked ZIPs remain ordinary
    files. This also preserves digests of core/ and ancestor output directories.
    """
    inventory = children if children is not None else sorted(root.rglob("*"))
    archives = {}
    for path in inventory:
        if path.name not in ARCHIVE_BASENAMES:
            continue
        observe(path)
        location = _archive_location(path)
        if location is None:
            continue
        base, name = location
        raw, expected = _paths(base, name)
        _require(path == expected and not raw.exists(), f"Incomplete conversion; raw and ZIP both exist: {raw}")
        archives[path] = (base, name, raw.relative_to(root).as_posix())
    if not archives:
        return None
    rows = []
    for path in inventory:
        observe(path)
        _require(not path.is_symlink(), f"Symlinked provenance member: {path}")
        if path in archives:
            base, name, prefix = archives[path]
            records = _logical_archive_records(base, name)
            rows.extend((prefix + "/" + row["path"], row["sha256"], row["size"]) for row in records)
        elif path.is_file():
            before = path.stat()
            with _source(path, before) as handle:
                digest, size = _hash(handle)
            rows.append((path.relative_to(root).as_posix(), digest, size))
    # Match Path ordering used by legacy directory digests, including a/ vs a.tsv.
    return sorted(rows, key=lambda row: PurePosixPath(row[0]).parts)


@contextlib.contextmanager
def materialized_directory(root, name, scratch, *, namespace="orthofinder"):
    """Give native tools a raw directory, expanding only into caller-owned scratch."""
    root = Path(root).absolute()
    allowed = _namespace_directories(namespace)
    _require(name in allowed, "Directory is outside the selected storage namespace")
    with storage_lock(root, exclusive=False, namespace=namespace):
        raw, path = _paths(root, name)
        if raw.exists() and not path.exists():
            yield raw
            return
        _require(not raw.exists(), f"Incomplete conversion: {raw}")
        archive_records(root, name)
        with tempfile.TemporaryDirectory(prefix="gg-orthofinder-", dir=scratch) as temporary:
            stage = Path(temporary)
            target = stage / (name + ".zip")
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, target)
            _materialize(stage, name)
            yield stage / name


@contextlib.contextmanager
def open_binary(root, relative, *, namespace="orthofinder"):
    """Read a result member without expanding its directory onto NFS."""
    root, relative = Path(root).absolute(), PurePosixPath(relative)
    _require(not relative.is_absolute() and ".." not in relative.parts, "Unsafe logical path")
    allowed = _namespace_directories(namespace)
    with storage_lock(root, exclusive=False, namespace=namespace):
        path = root / relative
        _physical(path)
        for name in allowed:
            if PurePosixPath(name) in relative.parents:
                raw, archive_path = _paths(root, name)
                _require(not (raw.exists() and archive_path.exists()), f"Incomplete conversion: {raw}")
        if path.is_file():
            with _source(path, path.stat()) as handle:
                yield handle
            return
        for name in allowed:
            prefix = PurePosixPath(name)
            if prefix in relative.parents:
                _, archive_path = _paths(root, name)
                before = _signature(archive_path.stat())
                with zipfile.ZipFile(archive_path) as archive:
                    _manifest(archive, archive_path, name)
                    with archive.open(Path(name).name + "/" + relative.relative_to(prefix).as_posix()) as handle:
                        yield handle
                _require(before == _signature(archive_path.stat()), f"Archive changed during read: {archive_path}")
                return
        raise FileNotFoundError(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=("pack", "materialize", "verify", "status"))
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--directory", choices=MANAGED_DIRECTORIES, action="append")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--progress-interval", type=float, default=30)
    parser.add_argument("--producer-lock-token", help=argparse.SUPPRESS)
    args = parser.parse_args()
    reporter = DIRECTORY.ProgressReporter(args.progress_interval)
    started = time.monotonic()
    if args.operation in ("pack", "materialize"):
        result = convert(
            args.root,
            to="zip" if args.operation == "pack" else "raw",
            directories=args.directory,
            workers=args.workers,
            producer_lock_token=args.producer_lock_token,
            progress_callback=reporter.update,
        )
    elif args.operation == "verify":
        result = {
            "directories": [
                archive_records(args.root, name, deep=True, progress_callback=reporter.update)["directory"]
                for name in args.directory or MANAGED_DIRECTORIES
                if (args.root / (name + ".zip")).exists()
            ]
        }
    else:
        result = {"directories": status(args.root, args.directory)}
    print(json.dumps(dict(result, elapsed_seconds=time.monotonic() - started)))


if __name__ == "__main__":
    main()

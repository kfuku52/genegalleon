"""Import verified local copies of explicit public downloads."""

import hashlib
import os
import re
import shutil
import stat
import tempfile
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlparse

from .local import validate_gzip_with_cache
from .locking import _fsync_directory, acquire_download_lock, release_download_lock


@dataclass(frozen=True)
class DownloadImport:
    path: Path
    size: int
    algorithm: str
    checksum: str


def manifest_download_imports(row):
    imports = {}
    for role in ("cds", "gff", "gbff", "genome"):
        path, origin, size, checksum = (
            (row.get(f"{role}_import_{field}") or "").strip()
            for field in ("path", "url", "size", "checksum")
        )
        if not any((path, origin, size, checksum)):
            continue
        if not all((path, origin, size, checksum)):
            raise ValueError(f"{role} import requires path, URL, size and checksum")
        url = (row.get(f"{role}_url") or "").strip()
        parsed = urlparse(url)
        if ((row.get("provider") or "").strip().lower() != "direct"
                or parsed.scheme not in ("http", "https", "ftp") or not parsed.hostname
                or parsed.username is not None or row.get(f"{role}_archive_member")):
            raise ValueError(f"{role} import requires a direct public file URL without an archive member")
        if origin != url:
            raise ValueError(f"{role} import URL does not match the manifest source URL")
        if not Path(path).is_absolute():
            raise ValueError(f"{role} import path must be absolute")
        if not re.fullmatch(r"[1-9][0-9]*", size):
            raise ValueError(f"{role} import size must be a positive byte count")
        algorithm, _, digest = checksum.lower().partition(":")
        if not ((algorithm == "sha256" and re.fullmatch(r"[a-f0-9]{64}", digest))
                or (algorithm == "md5" and re.fullmatch(r"[a-f0-9]{32}", digest))):
            raise ValueError(f"{role} import checksum must be sha256:<hex> or md5:<hex>")
        imports[role.upper()] = DownloadImport(Path(path), int(size), algorithm, digest)
    return imports


def verify_import(path, spec, destination):
    # Open without following symlinks or blocking on named pipes.
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(descriptor, "rb") as source:
        before = os.fstat(source.fileno())
        if not stat.S_ISREG(before.st_mode) or before.st_size != spec.size:
            raise ValueError("import file is not a regular file of the expected size")
        digest = hashlib.file_digest(source, lambda: hashlib.new(spec.algorithm, usedforsecurity=False))
        if digest.hexdigest() != spec.checksum:
            raise ValueError(f"import file {spec.algorithm} checksum mismatch")
        error = validate_gzip_with_cache(path, expected_path=destination)
        if error is not None:
            raise ValueError(f"import file failed format validation: {error}")
        current = path.stat(follow_symlinks=False)
        if file_signature(before) != file_signature(current):
            raise ValueError("import file changed during validation")


def file_signature(info):
    return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns


def materialize_download_import(spec, destination, overwrite, lock_stale_seconds, warnings, context):
    destination.parent.mkdir(parents=True, exist_ok=True)
    lock = destination.with_name(destination.name + ".lock")
    ownership = acquire_download_lock(lock, lock_stale_seconds, warnings, context)
    try:
        if destination.is_symlink():
            raise ValueError("import destination is a symlink")
        if destination.exists() and not overwrite:
            verify_import(destination, spec, destination)
            return False
        descriptor = os.open(spec.path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
        with os.fdopen(descriptor, "rb") as source:
            info = os.fstat(source.fileno())
            if not stat.S_ISREG(info.st_mode) or info.st_size != spec.size:
                raise ValueError("import source is not a regular file of the expected size")
            descriptor, temporary = tempfile.mkstemp(dir=destination.parent, prefix=".gg-import-")
            try:
                with os.fdopen(descriptor, "wb") as staged:
                    shutil.copyfileobj(source, staged, length=1024 * 1024)
                    staged.flush()
                    verify_import(Path(temporary), spec, destination)
                    os.fsync(staged.fileno())
                if overwrite:
                    os.replace(temporary, destination)
                else:
                    os.link(temporary, destination)
                _fsync_directory(destination.parent)
            finally:
                Path(temporary).unlink(missing_ok=True)
    finally:
        release_download_lock(lock, ownership)
    return True

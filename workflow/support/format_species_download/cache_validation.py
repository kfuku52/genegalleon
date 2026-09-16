"""Persistent, fail-closed validation receipts for downloaded gzip files."""

from __future__ import annotations

import hashlib
import json
import os
import stat
import tempfile
import threading
import time
from pathlib import Path

GZIP_VALIDATION_CACHE_SCHEMA = 1


def build_gzip_validation_key(relative_target, source_url="", archive_member=""):
    """Return a stable identity for one downloaded manifest target.

    The target is deliberately relative to the download root.  Array plans
    use a hash-specific staging directory, so an absolute target path would
    prevent a validated file from being reused by a later plan.
    """
    payload = {
        "archive_member": str(archive_member or ""),
        "relative_target": str(relative_target).replace(os.sep, "/"),
        "source_url": str(source_url or ""),
    }
    encoded = json.dumps(payload, ensure_ascii=True, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def gzip_validation_key_for_target(target, download_root, source_url="", archive_member=""):
    """Build a stable validation key from a target path and its download root."""
    target = Path(target)
    download_root = Path(download_root)
    try:
        relative_target = target.relative_to(download_root)
    except ValueError:
        relative_target = Path(os.path.relpath(str(target), str(download_root)))
    return build_gzip_validation_key(relative_target, source_url, archive_member)


def _file_identity(path):
    stat_result = Path(path).stat()
    if not stat.S_ISREG(stat_result.st_mode):
        raise OSError("gzip validation requires a regular file: {}".format(path))
    # Do not include ctime: creating a hardlink changes ctime on some filesystems
    # even though the bytes and inode are unchanged.  Hardlinked plan caches are
    # an intended reuse case.
    return {
        "st_dev": int(stat_result.st_dev),
        "st_ino": int(stat_result.st_ino),
        "st_mtime_ns": int(stat_result.st_mtime_ns),
        "st_size": int(stat_result.st_size),
    }


def _atomic_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True) + "\n"
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(data)
        handle.flush()
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


class GzipValidationCache:
    """A sidecar receipt store shared by plan-specific download directories.

    Each receipt is keyed by source identity and relative target, rather than
    by an absolute plan directory.  A receipt is usable only when the current
    file has the same device, inode, size, and mtime.  Any malformed, missing,
    or stale receipt falls back to a complete gzip read.
    """

    def __init__(self, cache_dir):
        self.cache_dir = Path(cache_dir).expanduser().resolve()
        self._lock = threading.Lock()
        self._hits = 0
        self._misses = 0
        self._records = 0

    def _receipt_path(self, validation_key):
        key = str(validation_key or "").strip().lower()
        if len(key) != 64 or any(character not in "0123456789abcdef" for character in key):
            return None
        return self.cache_dir / (key + ".json")

    def _count(self, field):
        with self._lock:
            setattr(self, field, getattr(self, field) + 1)

    def is_valid(self, path, validation_key):
        receipt_path = self._receipt_path(validation_key)
        if receipt_path is None:
            self._count("_misses")
            return False
        try:
            payload = json.loads(receipt_path.read_text(encoding="utf-8"))
            if payload.get("schema_version") != GZIP_VALIDATION_CACHE_SCHEMA:
                raise ValueError("unsupported gzip validation receipt schema")
            if payload.get("validation_key") != str(validation_key):
                raise ValueError("gzip validation receipt key mismatch")
            expected = payload.get("file_identity")
            if not isinstance(expected, dict) or _file_identity(path) != expected:
                raise ValueError("gzip validation receipt does not match file identity")
        except (OSError, ValueError, TypeError, json.JSONDecodeError):
            self._count("_misses")
            return False
        self._count("_hits")
        return True

    def record(self, path, validation_key, *, source_url="", archive_member="", relative_target=""):
        receipt_path = self._receipt_path(validation_key)
        if receipt_path is None:
            return False
        try:
            identity = _file_identity(path)
        except OSError:
            return False
        payload = {
            "archive_member": str(archive_member or ""),
            "file_identity": identity,
            "relative_target": str(relative_target),
            "schema_version": GZIP_VALIDATION_CACHE_SCHEMA,
            "source_url": str(source_url or ""),
            "validated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "validation_key": str(validation_key),
        }
        try:
            _atomic_json(receipt_path, payload)
        except OSError:
            # A read-only or temporarily unavailable sidecar directory must not
            # turn a successfully validated download into a failed download.
            return False
        self._count("_records")
        return True

    def validate(self, path, validation_key, *, source_url="", archive_member="", relative_target="", validator):
        """Validate a gzip, using a receipt when possible.

        ``validator`` is injected by the download runtime to keep this module
        independent of gzip I/O and to make the cache behavior easy to test.
        """
        if self.is_valid(path, validation_key):
            return None
        before = None
        try:
            before = _file_identity(path)
        except OSError:
            before = None
        error = validator(path)
        if error is None:
            try:
                after = _file_identity(path)
            except OSError as exc:
                error = exc
            else:
                if before is None or before != after:
                    error = OSError("file changed during gzip validation: {}".format(path))
        if error is None:
            self.record(
                path,
                validation_key,
                source_url=source_url,
                archive_member=archive_member,
                relative_target=relative_target,
            )
        return error

    def diagnostics(self):
        with self._lock:
            return {
                "validation_cache_hits": self._hits,
                "validation_cache_misses": self._misses,
                "validation_cache_records": self._records,
            }

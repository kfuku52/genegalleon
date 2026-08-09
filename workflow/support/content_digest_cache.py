#!/usr/bin/env python3
"""Best-effort persistent cache for content digests.

Cache hits require an exact match of device, inode, size, nanosecond mtime, and
nanosecond ctime.  A cache failure never weakens provenance checks: callers
fall back to reading and hashing the content.
"""

from __future__ import annotations

import hashlib
import os
import re
import sqlite3
import threading
from pathlib import Path
from typing import BinaryIO

CHUNK_SIZE = 1024 * 1024
_CACHE: "ContentDigestCache | None" = None


def _sha256_stream(handle: BinaryIO) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    while chunk := handle.read(CHUNK_SIZE):
        digest.update(chunk)
        size += len(chunk)
    return digest.hexdigest(), size


class ContentDigestCache:
    """Small SQLite cache shared by independent workflow processes."""

    def __init__(self, database: Path):
        self.database = database.absolute()
        self.enabled = True
        self._connection: sqlite3.Connection | None = None
        self._lock = threading.RLock()
        try:
            self.database.parent.mkdir(parents=True, exist_ok=True)
            self._connection = sqlite3.connect(self.database, timeout=30, check_same_thread=False)
            self._connection.execute("PRAGMA busy_timeout=30000")
            self._connection.execute("PRAGMA synchronous=NORMAL")
            with self._lock:
                with self._connect() as connection:
                    connection.execute(
                        """
                        CREATE TABLE IF NOT EXISTS digest_cache (
                            namespace TEXT NOT NULL,
                            path TEXT NOT NULL,
                            device INTEGER NOT NULL,
                            inode INTEGER NOT NULL,
                            stat_size INTEGER NOT NULL,
                            mtime_ns INTEGER NOT NULL,
                            ctime_ns INTEGER NOT NULL,
                            digest TEXT NOT NULL,
                            logical_size INTEGER NOT NULL,
                            member_count INTEGER,
                            PRIMARY KEY (namespace, path)
                        )
                        """
                    )
        except (OSError, sqlite3.Error):
            self.enabled = False
            self.close()

    def _connect(self) -> sqlite3.Connection:
        if self._connection is None:
            raise sqlite3.ProgrammingError("Digest cache connection is closed")
        return self._connection

    def close(self) -> None:
        with self._lock:
            if self._connection is not None:
                self._connection.close()
                self._connection = None

    @staticmethod
    def signature(path: Path) -> tuple[int, int, int, int, int]:
        stat = path.stat()
        return (
            int(stat.st_dev),
            int(stat.st_ino),
            int(stat.st_size),
            int(stat.st_mtime_ns),
            int(stat.st_ctime_ns),
        )

    def get(self, namespace: str, path: Path) -> tuple[str, int, int | None] | None:
        if not self.enabled:
            return None
        try:
            signature_before = self.signature(path)
            with self._lock:
                with self._connect() as connection:
                    row = connection.execute(
                        """
                        SELECT device, inode, stat_size, mtime_ns, ctime_ns,
                               digest, logical_size, member_count
                        FROM digest_cache WHERE namespace = ? AND path = ?
                        """,
                        (namespace, str(path.absolute())),
                    ).fetchone()
            signature_after = self.signature(path)
            if signature_before != signature_after:
                return None
            if row is None or tuple(int(value) for value in row[:5]) != signature_after:
                return None
            digest = str(row[5])
            if re.fullmatch(r"[0-9a-f]{64}", digest) is None:
                return None
            return digest, int(row[6]), None if row[7] is None else int(row[7])
        except (OSError, sqlite3.Error, ValueError, TypeError):
            return None

    def put(
        self,
        namespace: str,
        path: Path,
        digest: str,
        logical_size: int,
        member_count: int | None = None,
        expected_signature: tuple[int, int, int, int, int] | None = None,
    ) -> None:
        if not self.enabled:
            return
        try:
            signature = self.signature(path)
            if expected_signature is not None and signature != expected_signature:
                return
            with self._lock:
                with self._connect() as connection:
                    connection.execute(
                        """
                        INSERT INTO digest_cache (
                            namespace, path, device, inode, stat_size, mtime_ns,
                            ctime_ns, digest, logical_size, member_count
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                        ON CONFLICT(namespace, path) DO UPDATE SET
                            device=excluded.device,
                            inode=excluded.inode,
                            stat_size=excluded.stat_size,
                            mtime_ns=excluded.mtime_ns,
                            ctime_ns=excluded.ctime_ns,
                            digest=excluded.digest,
                            logical_size=excluded.logical_size,
                            member_count=excluded.member_count
                        """,
                        (
                            namespace,
                            str(path.absolute()),
                            *signature,
                            digest,
                            int(logical_size),
                            member_count,
                        ),
                    )
        except (OSError, sqlite3.Error, ValueError):
            return


def configure(database: Path | None) -> None:
    """Configure the process-wide cache; ``None`` disables it."""

    global _CACHE
    resolved_database = None if database is None else database.absolute()
    if _CACHE is not None and resolved_database is not None and _CACHE.database == resolved_database and _CACHE.enabled:
        return
    if _CACHE is not None:
        _CACHE.close()
    _CACHE = None if resolved_database is None else ContentDigestCache(resolved_database)


def configure_from_environment() -> None:
    raw_path = os.environ.get("GG_CONTENT_DIGEST_CACHE", "").strip()
    if raw_path:
        configure(Path(raw_path))


def cache() -> ContentDigestCache | None:
    return _CACHE


def cached_sha256_file(path: Path) -> tuple[str, int]:
    current_cache = cache()
    if current_cache is not None:
        cached = current_cache.get("file_sha256", path)
        if cached is not None:
            return cached[0], cached[1]
    signature_before = None
    if current_cache is not None:
        try:
            signature_before = current_cache.signature(path)
        except OSError:
            pass
    with path.open("rb") as handle:
        digest, size = _sha256_stream(handle)
    if current_cache is not None and signature_before is not None:
        current_cache.put(
            "file_sha256",
            path,
            digest,
            size,
            expected_signature=signature_before,
        )
    return digest, size


configure_from_environment()

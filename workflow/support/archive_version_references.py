"""References to exact published members of the gene-family archive."""
from __future__ import annotations

import contextlib
import hashlib
import json
import re
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Mapping, Sequence


class ArchiveReferenceError(RuntimeError):
    """A requested member version cannot be retained or resolved."""


def _canonical(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


@dataclass(frozen=True)
class MemberVersion:
    logical_path: str
    generation: int
    sha256: str

    def __post_init__(self):
        if not isinstance(self.logical_path, str) or not isinstance(self.sha256, str):
            raise ArchiveReferenceError("Invalid archive member version")
        pieces = self.logical_path.split("/")
        if (len(pieces) != 2 or any(not item or item.startswith(".") or chr(92) in item or chr(0) in item
                                   for item in pieces)
                or type(self.generation) is not int or self.generation < 1
                or re.fullmatch(r"[0-9a-f]{64}", self.sha256) is None):
            raise ArchiveReferenceError("Invalid archive member version")

    @property
    def key(self) -> str:
        value = [self.logical_path, self.generation, self.sha256]
        return hashlib.sha256(_canonical(value).encode()).hexdigest()


@dataclass(frozen=True)
class VersionLocation:
    version: MemberVersion
    zip_path: str
    record: Mapping

    def __post_init__(self):
        if not Path(self.zip_path).is_absolute():
            raise ArchiveReferenceError("Archive version location must be absolute")
        if (self.record.get("member_name") != self.version.logical_path
                or self.record.get("generation") != self.version.generation
                or self.record.get("sha256") != self.version.sha256):
            raise ArchiveReferenceError("Archive member record does not identify the requested version")


class ArchiveVersionReferences:
    """Native archive metadata; callers hold the archive publication lock for writes."""

    def __init__(self, archive_root: Path):
        self.path = Path(archive_root).absolute() / "version-references.sqlite3"
        self.retained = self.path.parent / "retained-versions"

    @contextlib.contextmanager
    def _connection(self, *, write: bool = False, create: bool = False) -> Iterator[sqlite3.Connection | None]:
        if self.path.is_symlink() or self.path.parent.is_symlink():
            raise ArchiveReferenceError("Symlinked archive reference metadata")
        if self.path.exists() and not self.path.is_file():
            raise ArchiveReferenceError("Archive reference metadata is not a regular file")
        if not self.path.exists() and not create:
            if write:
                raise ArchiveReferenceError("No archive reference catalogue has been created")
            yield None
            return
        mode = "rwc" if create else "rw" if write else "ro"
        with contextlib.closing(sqlite3.connect(self.path.as_uri() + "?mode=" + mode, uri=True, timeout=15)) as db:
            db.execute("PRAGMA foreign_keys=ON")
            if create:
                # Keep the schema even when a first reservation exceeds quota.
                with db:
                    db.execute("BEGIN IMMEDIATE")
                    db.execute("CREATE TABLE IF NOT EXISTS consumers (id TEXT PRIMARY KEY, selection TEXT NOT NULL,"
                               "released INTEGER NOT NULL DEFAULT 0)")
                    db.execute("CREATE TABLE IF NOT EXISTS versions (id TEXT PRIMARY KEY, logical_path TEXT NOT NULL,"
                               "generation INTEGER NOT NULL, sha256 TEXT NOT NULL, zip_path TEXT NOT NULL, record TEXT NOT NULL)")
                    db.execute("CREATE INDEX IF NOT EXISTS archive_version_zip ON versions(zip_path)")
                    db.execute("CREATE TABLE IF NOT EXISTS references_v1 (consumer TEXT NOT NULL REFERENCES consumers(id)"
                               " ON DELETE CASCADE, version TEXT NOT NULL REFERENCES versions(id), PRIMARY KEY (consumer,version))")
                    db.execute("CREATE INDEX IF NOT EXISTS archive_version_readers ON references_v1(version)")
                    db.execute("CREATE TABLE IF NOT EXISTS pending_removal (zip_path TEXT PRIMARY KEY,"
                               " consumer TEXT NOT NULL REFERENCES consumers(id))")
            with db:
                if write:
                    db.execute("BEGIN IMMEDIATE")
                yield db

    @staticmethod
    def _require_consumer(consumer: str):
        if re.fullmatch(r"[a-z][a-z0-9-]{11,95}", consumer) is None:
            raise ArchiveReferenceError("Invalid archive reference owner")

    def reserve(self, consumer: str, locations: Sequence[VersionLocation], *, max_bytes: int) -> str:
        """Record a canonical selection once. Repeating it adds no copies."""
        self._require_consumer(consumer)
        if type(max_bytes) is not int or max_bytes < 1:
            raise ArchiveReferenceError("Archive reference capacity must be positive")
        keys = sorted(location.version.key for location in locations)
        if not keys or len(keys) != len(set(keys)):
            raise ArchiveReferenceError("An archive selection must contain unique versions")
        selection = hashlib.sha256(_canonical(keys).encode()).hexdigest()
        with self._connection(write=True, create=True) as db:
            before = db.execute("SELECT selection,released FROM consumers WHERE id=?", (consumer,)).fetchone()
            if before is not None:
                if before != (selection, 0):
                    raise ArchiveReferenceError("Consumer already registered a different or released archive selection")
                return selection
            db.execute("INSERT INTO consumers(id,selection) VALUES (?,?)", (consumer, selection))
            for source in locations:
                version = source.version
                db.execute("INSERT OR IGNORE INTO versions VALUES (?,?,?,?,?,?)",
                           (version.key, version.logical_path, version.generation, version.sha256,
                            source.zip_path, _canonical(dict(source.record))))
                db.execute("INSERT INTO references_v1 VALUES (?,?)", (consumer, version.key))
            total, = db.execute("SELECT COALESCE(SUM(json_extract(record, '$.size')),0) FROM versions").fetchone()
            if total > max_bytes:
                raise ArchiveReferenceError(f"Archive reference capacity exceeded: {total} > {max_bytes}")
        return selection

    def member(self, consumer: str, selection: str, version: MemberVersion) -> VersionLocation:
        self._require_consumer(consumer)
        with self._connection() as db:
            row = None if db is None else db.execute("""
                SELECT v.zip_path,v.record FROM versions AS v
                JOIN references_v1 AS r ON r.version=v.id
                JOIN consumers AS c ON c.id=r.consumer
                WHERE c.id=? AND c.selection=? AND v.id=?
                """, (consumer, selection, version.key)).fetchone()
        if row is None:
            raise ArchiveReferenceError("No held reference for the exact archive version")
        return VersionLocation(version, row[0], json.loads(row[1]))

    def locations_for_archive(self, archive: Path) -> list[VersionLocation]:
        with self._connection() as db:
            rows = [] if db is None else db.execute("""
                SELECT logical_path,generation,sha256,record FROM versions WHERE zip_path=?
                ORDER BY logical_path,generation,sha256
                """, (str(archive),)).fetchall()
        return [VersionLocation(MemberVersion(path, generation, sha256), str(archive), json.loads(record))
                for path, generation, sha256, record in rows]

    def relocate(self, locations: Sequence[VersionLocation]) -> None:
        """Update storage locations only. The consumer selections cannot change."""
        if not locations:
            return
        with self._connection(write=True) as db:
            for item in locations:
                changed = db.execute("UPDATE versions SET zip_path=?, record=? WHERE id=?",
                                     (item.zip_path, _canonical(dict(item.record)), item.version.key))
                if changed.rowcount != 1:
                    raise ArchiveReferenceError("Archive version disappeared during relocation")

    def release(self, consumer: str, selection: str) -> list[Path]:
        """Release exact references after the owning scheduler verifies termination."""
        self._require_consumer(consumer)
        with self._connection(write=True) as db:
            current = db.execute("SELECT selection,released FROM consumers WHERE id=?", (consumer,)).fetchone()
            if current is None or current[0] != selection:
                raise ArchiveReferenceError("Archive selection changed")
            if not current[1]:
                db.execute("UPDATE consumers SET released=1 WHERE id=?", (consumer,))
                db.execute("DELETE FROM references_v1 WHERE consumer=?", (consumer,))
                orphaned = "NOT EXISTS (SELECT 1 FROM references_v1 WHERE version=versions.id)"
                dropped = [row[0] for row in db.execute("SELECT DISTINCT zip_path FROM versions WHERE " + orphaned)]
                db.execute("DELETE FROM versions WHERE " + orphaned)
                live = {row[0] for row in db.execute("SELECT DISTINCT zip_path FROM versions WHERE zip_path IN ("
                                                     "SELECT value FROM json_each(?))", (_canonical(dropped),))}
                removable = [(path, consumer) for path in dropped if path not in live
                             and Path(path).parent.parent == self.retained]
                db.executemany("INSERT OR IGNORE INTO pending_removal VALUES (?,?)", removable)
            return [Path(path) for path, in db.execute(
                "SELECT zip_path FROM pending_removal WHERE consumer=?", (consumer,))]

    def record_removal(self, consumer: str, path: Path) -> None:
        with self._connection(write=True) as db:
            db.execute("DELETE FROM pending_removal WHERE zip_path=? AND consumer=?", (str(path), consumer))

    def has_references(self) -> bool:
        with self._connection() as db:
            return db is not None and db.execute("SELECT 1 FROM references_v1 LIMIT 1").fetchone() is not None

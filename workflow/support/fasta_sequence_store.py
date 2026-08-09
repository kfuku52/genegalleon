#!/usr/bin/env python3
"""Build and query a reusable SQLite index for species FASTA inputs."""

from __future__ import annotations

import argparse
import contextlib
import fcntl
import gzip
import hashlib
import json
import os
import shutil
import sqlite3
import sys
import tempfile
import textwrap
import zlib
from pathlib import Path
from typing import Iterator, TextIO

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from content_digest_cache import cached_sha256_file
from content_digest_cache import configure as configure_digest_cache

SCHEMA_VERSION = 2
DEFAULT_MAX_DATABASE_BYTES = 20 * 1024 * 1024 * 1024
DEFAULT_MINIMUM_FREE_BYTES = 1024 * 1024 * 1024


class SequenceStoreError(RuntimeError):
    pass


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def fasta_records(path: Path) -> Iterator[tuple[str, str, str]]:
    header = ""
    sequence: list[str] = []
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line.startswith(">"):
                if header:
                    identifier = header.split(None, 1)[0]
                    yield identifier, header, "".join(sequence)
                header = line[1:]
                sequence = []
            elif header:
                sequence.append(line.strip())
            elif line.strip():
                raise SequenceStoreError(f"FASTA sequence data precedes the first header: {path}")
    if header:
        identifier = header.split(None, 1)[0]
        yield identifier, header, "".join(sequence)


def read_sources(path: Path) -> list[tuple[Path, str]]:
    sources: list[tuple[Path, str]] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t", 1)
            source = Path(fields[0]).absolute()
            species = fields[1] if len(fields) == 2 else ""
            if not source.is_file() or source.is_symlink():
                raise SequenceStoreError(f"Invalid FASTA source at line {line_number}: {source}")
            sources.append((source, species))
    if not sources:
        raise SequenceStoreError(f"No FASTA sources were listed: {path}")
    return sources


def source_signature(path: Path, species: str) -> dict[str, object]:
    stat = path.stat()
    return {
        "path": str(path),
        "species": species,
        "device": int(stat.st_dev),
        "inode": int(stat.st_ino),
        "size": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
        "ctime_ns": int(stat.st_ctime_ns),
    }


def file_signature(path: Path) -> dict[str, int]:
    stat = path.stat()
    return {
        "device": int(stat.st_dev),
        "inode": int(stat.st_ino),
        "size": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
        "ctime_ns": int(stat.st_ctime_ns),
    }


def database_state_path(database: Path) -> Path:
    return database.with_suffix(database.suffix + ".state.json")


def sha256_file_uncached(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def record_database_state(database: Path) -> None:
    signature_before = file_signature(database)
    digest = sha256_file_uncached(database)
    signature_after = file_signature(database)
    if signature_before != signature_after:
        raise SequenceStoreError(f"FASTA sequence database changed while it was hashed: {database}")
    write_json_atomic(
        database_state_path(database),
        {
            "schema_version": SCHEMA_VERSION,
            "database_signature": signature_after,
            "sha256": digest,
        },
    )


def database_content_current(database: Path) -> bool:
    state_path = database_state_path(database)
    if not database.is_file() or database.is_symlink() or not state_path.is_file() or state_path.is_symlink():
        return False
    try:
        state = json.loads(state_path.read_text(encoding="utf-8"))
        if state.get("schema_version") != SCHEMA_VERSION:
            return False
        signature_before = file_signature(database)
        if state.get("database_signature") == signature_before:
            return True
        digest = sha256_file_uncached(database)
        signature_after = file_signature(database)
        if signature_before != signature_after or state.get("sha256") != digest:
            return False
        # Copies and atomic moves legitimately change inode/stat identity. Once
        # content is verified, refresh the private fast-path state atomically.
        write_json_atomic(
            state_path,
            {
                "schema_version": SCHEMA_VERSION,
                "database_signature": signature_after,
                "sha256": digest,
            },
        )
        return True
    except (OSError, TypeError, ValueError, json.JSONDecodeError):
        return False


def manifest_current(database: Path, manifest: Path, sources: list[tuple[Path, str]]) -> bool:
    if not database.is_file() or database.is_symlink() or not manifest.is_file() or manifest.is_symlink():
        return False
    try:
        if not database_content_current(database):
            return False
        payload = json.loads(manifest.read_text(encoding="utf-8"))
        if payload.get("schema_version") != SCHEMA_VERSION:
            return False
        connection = sqlite3.connect(f"file:{database}?mode=ro", uri=True)
        try:
            metadata = dict(connection.execute("SELECT key, value FROM metadata"))
        finally:
            connection.close()
        if int(metadata.get("schema_version", "-1")) != SCHEMA_VERSION:
            return False
        expected_signatures = [source_signature(path, species) for path, species in sources]
        if json.loads(metadata.get("source_signatures", "null")) != expected_signatures:
            return False
        return json.loads(metadata.get("source_identity", "null")) == payload.get("sources")
    except (OSError, ValueError, TypeError, json.JSONDecodeError, sqlite3.Error):
        return False


def refresh_signatures_if_content_current(
    database: Path,
    manifest: Path,
    sources: list[tuple[Path, str]],
) -> bool:
    """Refresh private path/stat state when the portable content identity matches."""

    if not database.is_file() or database.is_symlink() or not manifest.is_file() or manifest.is_symlink():
        return False
    try:
        if not database_content_current(database):
            return False
        payload = json.loads(manifest.read_text(encoding="utf-8"))
        if payload.get("schema_version") != SCHEMA_VERSION:
            return False
        recorded = payload.get("sources")
        if not isinstance(recorded, list) or len(recorded) != len(sources):
            return False
        current_identity: list[dict[str, object]] = []
        signatures: list[dict[str, object]] = []
        for (source, species), entry in zip(sources, recorded):
            if not isinstance(entry, dict) or entry.get("species") != species:
                return False
            signature_before = source_signature(source, species)
            digest, _size = cached_sha256_file(source)
            signature_after = source_signature(source, species)
            if signature_before != signature_after:
                return False
            if entry.get("sha256") != digest:
                return False
            current_identity.append(entry)
            signatures.append(signature_after)
        with sqlite3.connect(database) as connection:
            metadata = dict(connection.execute("SELECT key, value FROM metadata"))
            if int(metadata.get("schema_version", "-1")) != SCHEMA_VERSION:
                return False
            if json.loads(metadata.get("source_identity", "null")) != current_identity:
                return False
            connection.execute(
                "UPDATE metadata SET value=? WHERE key='source_signatures'",
                (json.dumps(signatures, sort_keys=True),),
            )
        record_database_state(database)
        return True
    except (OSError, ValueError, TypeError, json.JSONDecodeError, sqlite3.Error):
        return False


@contextlib.contextmanager
def exclusive_lock(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def write_json_atomic(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def enforce_storage_budget(
    connection: sqlite3.Connection,
    directory: Path,
    max_database_bytes: int,
    minimum_free_bytes: int,
) -> None:
    page_size = int(connection.execute("PRAGMA page_size").fetchone()[0])
    page_count = int(connection.execute("PRAGMA page_count").fetchone()[0])
    database_bytes = page_size * page_count
    if max_database_bytes and database_bytes > max_database_bytes:
        raise SequenceStoreError(
            f"FASTA sequence database exceeded its {max_database_bytes}-byte size limit "
            f"(current allocation: {database_bytes} bytes)"
        )
    free_bytes = shutil.disk_usage(directory).free
    if minimum_free_bytes and free_bytes < minimum_free_bytes:
        raise SequenceStoreError(
            f"Insufficient free space while building the FASTA sequence database: "
            f"{free_bytes} bytes free, {minimum_free_bytes} bytes required"
        )


def build_store(
    database: Path,
    manifest: Path,
    sources: list[tuple[Path, str]],
    max_database_bytes: int,
    minimum_free_bytes: int,
) -> None:
    database.parent.mkdir(parents=True, exist_ok=True)
    free_bytes = shutil.disk_usage(database.parent).free
    if minimum_free_bytes and free_bytes < minimum_free_bytes:
        raise SequenceStoreError(
            f"Insufficient free space before building the FASTA sequence database: "
            f"{free_bytes} bytes free, {minimum_free_bytes} bytes required"
        )
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{database.name}.", dir=database.parent)
    os.close(descriptor)
    temporary = Path(temporary_name)
    temporary.unlink()
    source_rows: list[dict[str, object]] = []
    try:
        connection = sqlite3.connect(temporary)
        try:
            connection.executescript(
                """
                PRAGMA journal_mode=OFF;
                PRAGMA synchronous=OFF;
                PRAGMA temp_store=MEMORY;
                CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
                CREATE TABLE sequences (
                    identifier TEXT PRIMARY KEY,
                    identifier_lower TEXT NOT NULL,
                    header TEXT NOT NULL,
                    sequence BLOB NOT NULL,
                    species TEXT NOT NULL,
                    source_order INTEGER NOT NULL,
                    record_order INTEGER NOT NULL
                );
                CREATE INDEX sequences_lower_idx ON sequences(identifier_lower);
                """
            )
            connection.execute(
                "INSERT INTO metadata(key, value) VALUES('schema_version', ?)",
                (str(SCHEMA_VERSION),),
            )
            for source_order, (source, species) in enumerate(sources):
                signature_before = source_signature(source, species)
                digest, _size = cached_sha256_file(source)
                batch: list[tuple[str, str, str, bytes, str, int, int]] = []
                record_count = 0
                for record_order, (identifier, header, sequence) in enumerate(fasta_records(source)):
                    record_count += 1
                    if not identifier:
                        raise SequenceStoreError(f"Empty FASTA identifier: {source}")
                    batch.append(
                        (
                            identifier,
                            identifier.lower(),
                            header,
                            sqlite3.Binary(zlib.compress(sequence.encode("utf-8"), level=6)),
                            species,
                            source_order,
                            record_order,
                        )
                    )
                    if len(batch) >= 2000:
                        try:
                            connection.executemany("INSERT INTO sequences VALUES (?, ?, ?, ?, ?, ?, ?)", batch)
                        except sqlite3.IntegrityError as exc:
                            raise SequenceStoreError(
                                f"Duplicate FASTA identifier while indexing {source}: {exc}"
                            ) from exc
                        batch.clear()
                        enforce_storage_budget(
                            connection,
                            database.parent,
                            max_database_bytes,
                            minimum_free_bytes,
                        )
                if batch:
                    try:
                        connection.executemany("INSERT INTO sequences VALUES (?, ?, ?, ?, ?, ?, ?)", batch)
                    except sqlite3.IntegrityError as exc:
                        raise SequenceStoreError(f"Duplicate FASTA identifier while indexing {source}: {exc}") from exc
                    enforce_storage_budget(
                        connection,
                        database.parent,
                        max_database_bytes,
                        minimum_free_bytes,
                    )
                signature_after = source_signature(source, species)
                if signature_after != signature_before:
                    raise SequenceStoreError(f"FASTA source changed while it was indexed: {source}")
                signature_after["sha256"] = digest
                signature_after["record_count"] = record_count
                source_rows.append(signature_after)
            source_signatures = [
                {key: value for key, value in row.items() if key not in {"sha256", "record_count"}}
                for row in source_rows
            ]
            source_identity = [
                {
                    "species": row["species"],
                    "sha256": row["sha256"],
                    "record_count": row["record_count"],
                }
                for row in source_rows
            ]
            connection.executemany(
                "INSERT INTO metadata(key, value) VALUES(?, ?)",
                (
                    ("source_signatures", json.dumps(source_signatures, sort_keys=True)),
                    ("source_identity", json.dumps(source_identity, sort_keys=True)),
                ),
            )
            connection.commit()
            enforce_storage_budget(
                connection,
                database.parent,
                max_database_bytes,
                minimum_free_bytes,
            )
        finally:
            connection.close()
        os.replace(temporary, database)
        write_json_atomic(
            manifest,
            {
                "schema_version": SCHEMA_VERSION,
                # This public identity is deliberately independent of absolute
                # paths and timestamps so artifact provenance remains portable.
                # Fast stat signatures stay private inside the SQLite store.
                "sources": source_identity,
            },
        )
        record_database_state(database)
    finally:
        temporary.unlink(missing_ok=True)


def ensure(args: argparse.Namespace) -> int:
    if args.max_database_bytes < 0 or args.minimum_free_bytes < 0:
        raise SequenceStoreError("FASTA sequence store byte limits must be nonnegative integers")
    sources = read_sources(args.source_list)
    if manifest_current(args.database, args.manifest, sources):
        print(args.manifest)
        return 0
    with exclusive_lock(args.database.with_suffix(args.database.suffix + ".lock")):
        if not manifest_current(args.database, args.manifest, sources):
            if not refresh_signatures_if_content_current(
                args.database,
                args.manifest,
                sources,
            ):
                build_store(
                    args.database,
                    args.manifest,
                    sources,
                    args.max_database_bytes,
                    args.minimum_free_bytes,
                )
    print(args.manifest)
    return 0


def read_patterns(path: Path) -> list[str]:
    patterns: list[str] = []
    seen: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            pattern = raw_line.strip()
            if pattern and pattern not in seen:
                patterns.append(pattern)
                seen.add(pattern)
    return patterns


def write_record(handle: TextIO, header: str, sequence: str) -> None:
    handle.write(f">{header}\n")
    for line in textwrap.wrap(sequence, width=60) or [""]:
        handle.write(f"{line}\n")


def decode_sequence(value: object) -> str:
    if not isinstance(value, (bytes, bytearray, memoryview)):
        raise SequenceStoreError("Indexed FASTA sequence has an invalid storage encoding")
    try:
        return zlib.decompress(bytes(value)).decode("utf-8")
    except (UnicodeError, zlib.error) as exc:
        raise SequenceStoreError("Indexed FASTA sequence is corrupt") from exc


def extract(args: argparse.Namespace) -> int:
    patterns = read_patterns(args.pattern_file)
    if not patterns:
        args.output.write_text("", encoding="utf-8")
        return 0
    connection = sqlite3.connect(f"file:{args.database}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    rows: list[sqlite3.Row] = []
    try:
        connection.execute("CREATE TEMP TABLE requested(pattern TEXT PRIMARY KEY)")
        requested = patterns
        if args.query_variants:
            species = [str(row[0]) for row in connection.execute("SELECT DISTINCT species FROM sequences")]
            requested = list(patterns)
            requested.extend(f"{sp}_{pattern}" for sp in species for pattern in patterns if sp)
        normalize = (lambda value: value.lower()) if args.ignore_case else (lambda value: value)
        connection.executemany(
            "INSERT OR IGNORE INTO requested(pattern) VALUES (?)",
            ((normalize(pattern),) for pattern in requested),
        )
        identifier_column = "identifier_lower" if args.ignore_case else "identifier"
        rows = list(
            connection.execute(
                f"""
                SELECT s.identifier, s.header, s.sequence, s.species
                FROM sequences AS s JOIN requested AS r ON s.{identifier_column}=r.pattern
                ORDER BY s.source_order, s.record_order
                """
            )
        )
    finally:
        connection.close()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for row in rows:
            header = str(row["header"])
            if args.prefix_species:
                species = str(row["species"])
                identifier = str(row["identifier"])
                for separator in ("_", "-", "."):
                    prefix = f"{species}{separator}"
                    if species and identifier.startswith(prefix):
                        identifier = identifier[len(prefix) :]
                        break
                header = f"{species}_{identifier}" if species else identifier
            write_record(handle, header, decode_sequence(row["sequence"]))
    if args.require_all and len(rows) != len(patterns):
        found = {str(row["identifier"]).lower() if args.ignore_case else str(row["identifier"]) for row in rows}
        missing = [pattern for pattern in patterns if normalize(pattern) not in found]
        print(
            f"Expected {len(patterns)} indexed FASTA records but extracted {len(rows)}; missing: "
            + ", ".join(missing[:20]),
            file=sys.stderr,
        )
        return 3
    print(f"Extracted {len(rows)} indexed FASTA record(s) to {args.output}")
    return 0


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    subparsers = root.add_subparsers(dest="command", required=True)
    ensure_parser = subparsers.add_parser("ensure")
    ensure_parser.add_argument("--database", required=True, type=Path)
    ensure_parser.add_argument("--manifest", required=True, type=Path)
    ensure_parser.add_argument("--source-list", required=True, type=Path)
    ensure_parser.add_argument("--digest-cache", type=Path)
    ensure_parser.add_argument(
        "--max-database-bytes",
        type=int,
        default=os.environ.get("GG_FASTA_SEQUENCE_STORE_MAX_BYTES", str(DEFAULT_MAX_DATABASE_BYTES)),
        help="Maximum allocated SQLite size (0 disables; default: 20 GiB)",
    )
    ensure_parser.add_argument(
        "--minimum-free-bytes",
        type=int,
        default=os.environ.get("GG_FASTA_SEQUENCE_STORE_MIN_FREE_BYTES", str(DEFAULT_MINIMUM_FREE_BYTES)),
        help="Minimum filesystem free space retained during builds (0 disables; default: 1 GiB)",
    )
    extract_parser = subparsers.add_parser("extract")
    extract_parser.add_argument("--database", required=True, type=Path)
    extract_parser.add_argument("--pattern-file", required=True, type=Path)
    extract_parser.add_argument("--output", required=True, type=Path)
    extract_parser.add_argument("--ignore-case", action="store_true")
    extract_parser.add_argument("--query-variants", action="store_true")
    extract_parser.add_argument("--prefix-species", action="store_true")
    extract_parser.add_argument("--require-all", action="store_true")
    return root


def main() -> int:
    args = parser().parse_args()
    try:
        if args.command == "ensure":
            if args.digest_cache:
                configure_digest_cache(args.digest_cache)
            return ensure(args)
        if args.command == "extract":
            return extract(args)
    except (OSError, sqlite3.Error, SequenceStoreError, UnicodeError) as exc:
        print(f"FASTA sequence store error: {exc}", file=sys.stderr)
        return 2
    return 2


if __name__ == "__main__":
    raise SystemExit(main())

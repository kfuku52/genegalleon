"""Lossless SQLite packs for completed workflow metadata, never execution authority.

Writers publish ordinary files. An offline collector builds locally, verifies
all bodies, atomically publishes an immutable database, then retires unchanged
sources. Readers merge live files with archived records; live files take priority.
Scientific TSV/FASTA/tree payloads and distributed control state are excluded.
"""
from __future__ import annotations

import argparse
import contextlib
import contextvars
import hashlib
import json
import os
import re
import shutil
import sqlite3
import stat
import tempfile
import threading
import uuid
import zlib
from collections import OrderedDict, deque
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path, PurePosixPath

if __package__:
    from .shared_namespace_lock import namespace_lock
else:
    from shared_namespace_lock import namespace_lock

DATABASE = ".gg-metadata.sqlite3"
APPLICATION_ID = 0x47474D44
SCHEMA_VERSION = 1
READONLY_SOURCE_LINK_CHURN_V1 = True
MAX_BODY = 8 * 1024 * 1024
INDEX_ROOTS = {"index", "index-by-subdir"}
ROOTS = {"performance", "runtime", "observations", "resource_metrics", "artifact_provenance", "versions", "version_inventory", "output-inventory"}
_LOCAL = threading.local()
_SNAPSHOT = contextvars.ContextVar("gg_metadata_snapshot", default=None)


def signature(path):
    value = Path(path).lstat()
    return (value.st_dev, value.st_ino, value.st_mode, value.st_size,
            value.st_mtime_ns, value.st_ctime_ns, value.st_uid, value.st_nlink, value.st_gid)


def scope(path):
    path = Path(os.path.abspath(path))
    for root in (path, *path.parents):
        if ((root.name in ROOTS and (root.name != "performance" or root.parent.name == "versions"))
                or (root.name in INDEX_ROOTS and root.parent.name == ".gg_store")):
            if root.is_symlink():
                raise ValueError(f"Metadata root must be a physical directory: {root}")
            if not root.exists():
                return None
            if not root.is_dir():
                raise ValueError(f"Metadata root must be a physical directory: {root}")
            if path != root:
                for parent in path.parents:
                    if parent == root:
                        break
                    if parent.is_symlink():
                        raise ValueError(f"Symlinked metadata ancestor: {parent}")
            return root
    return None


def _key(root, path):
    raw = Path(path).relative_to(root).as_posix()
    relative = PurePosixPath(raw)
    if relative.is_absolute() or not relative.parts or any(p in {".", "..", ""} for p in relative.parts):
        raise ValueError("Unsafe metadata key")
    return raw


def _db_signature(path):
    value = signature(path)
    if not stat.S_ISREG(value[2]) or value[6] != os.getuid() or value[7] != 1 or value[2] & 0o022:
        raise ValueError(f"Metadata database must be an owned regular file: {path}")
    if any(Path(str(path) + suffix).exists() for suffix in ("-wal", "-shm", "-journal")):
        raise ValueError("Published metadata must be a closed, standalone SQLite database")
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(descriptor, "rb") as handle:
        opened = os.fstat(handle.fileno())
        if (opened.st_dev, opened.st_ino) != value[:2]:
            raise ValueError("Metadata database changed while opening its header")
        header = handle.read(100)
    if signature(path) != value:
        raise ValueError("Metadata database changed while reading its header")
    if len(header) != 100 or header[:16] != b"SQLite format 3\0":
        raise ValueError("Unsupported metadata database header")
    # SQLite commits increment this counter even when file size and coarse
    # filesystem timestamps remain unchanged. No full-database hash is needed.
    return (*value, int.from_bytes(header[24:28], "big"), int.from_bytes(header[92:96], "big"))


@contextlib.contextmanager
def read_snapshot():
    if _SNAPSHOT.get() is not None:
        yield
        return
    opened = {}
    token = _SNAPSHOT.set(opened)
    try:
        yield
        for path, connection in opened.items():
            current = None if not path.exists() else _db_signature(path)
            if current != (None if connection is None else connection[2]):
                raise ValueError("Metadata database changed during observation")
    finally:
        _SNAPSHOT.reset(token)


def _connection(root, observe=None):
    path = root / DATABASE
    snapshot = _SNAPSHOT.get()
    if snapshot is not None and path in snapshot:
        return snapshot[path]
    if observe is not None:
        observe(path)
    try:
        identity = _db_signature(path)
    except FileNotFoundError:
        if snapshot is not None:
            snapshot[path] = None
        return None
    cache = getattr(_LOCAL, "connections", None)
    if cache is None:
        cache = _LOCAL.connections = OrderedDict()
    key = (str(path), identity)
    if key not in cache:
        db = sqlite3.connect(path.absolute().as_uri() + "?mode=ro", uri=True)
        try:
            db.execute("PRAGMA query_only=ON")
            db.execute("PRAGMA trusted_schema=OFF")
            if (db.execute("PRAGMA application_id").fetchone()[0] != APPLICATION_ID
                    or db.execute("PRAGMA user_version").fetchone()[0] != SCHEMA_VERSION
                    or db.execute("SELECT value FROM metadata WHERE key='kind'").fetchone() != (root.name,)):
                raise ValueError("Unsupported or misplaced metadata database")
            if _db_signature(path) != identity:
                raise ValueError("Metadata database changed while opening")
        except BaseException:
            db.close()
            raise
        cache[key] = db
        while len(cache) > 8:
            cache.popitem(last=False)[1].close()
    cache.move_to_end(key)
    opened = cache[key], path, identity
    if snapshot is not None:
        snapshot[path] = opened
    return opened


def _query(root, sql, values=(), observe=None):
    opened = _connection(root, observe)
    if opened is None:
        return []
    db, path, before = opened
    rows = db.execute(sql, values).fetchall()
    if _SNAPSHOT.get() is None and _db_signature(path) != before:
        raise ValueError("Metadata database changed during read")
    return rows


def body_limit(root):
    return 128 * 1024 * 1024 if root is not None and root.name in (INDEX_ROOTS | {"performance"}) else MAX_BODY


def archived_bytes(path, *, limit=None, observe=None):
    root = scope(path)
    if root is None:
        return None
    limit = body_limit(root) if limit is None else limit
    query = "SELECT CASE WHEN length(body)<=? THEN body ELSE NULL END,sha256,bytes FROM files WHERE path=?"
    if root.name == "runtime":
        query = ("SELECT CASE WHEN length(b.body)<=? THEN b.body ELSE NULL END,f.sha256,f.bytes "
                 "FROM files f LEFT JOIN runtime_bodies b ON b.sha256=f.sha256 WHERE f.path=?")
    rows = _query(root, query, (limit * 2, _key(root, path)), observe)
    if not rows:
        return None
    body, expected, size = rows[0]
    if type(size) is not int or not 0 <= size <= limit or not isinstance(body, bytes):
        raise ValueError("Archived metadata exceeds read bound")
    decoder = zlib.decompressobj()
    raw = decoder.decompress(body, limit + 1)
    if (len(raw) != size or not decoder.eof or decoder.unused_data
            or hashlib.sha256(raw).hexdigest() != expected):
        raise ValueError("Archived metadata body is corrupt")
    return raw


def archived_exists(path, *, observe=None):
    root = scope(path)
    return bool(root is not None and _query(root, "SELECT 1 FROM files WHERE path=?", (_key(root, path),), observe))


def archived_children(directory, *, directories=False, observe=None):
    root = scope(directory)
    if root is None:
        return []
    prefix = "" if directory == root else _key(root, directory) + "/"
    rows = _query(root, "SELECT path FROM files WHERE path>=? AND path<? ORDER BY path",
                  (prefix, prefix + chr(0x10ffff)), observe)
    names = set()
    for (key,) in rows:
        _key(root, root / key)
        tail = key[len(prefix):]
        if directories:
            if "/" in tail:
                names.add(tail.split("/", 1)[0])
        elif "/" not in tail:
            names.add(tail)
    return [directory / name for name in sorted(names)]


@contextlib.contextmanager
def writer(path):
    root = scope(path)
    if root is None or not (root / DATABASE).exists():
        yield
    else:
        with namespace_lock(root / ".gg-metadata.guard", exclusive=False):
            yield


def _read_source(path, *, root=None):
    before = signature(path)
    runtime = root is not None and Path(root).name == "runtime"
    if (not stat.S_ISREG(before[2]) or before[6] != os.getuid()
            or (before[7] != 1 and not (runtime and not before[2] & 0o222))):
        raise ValueError(f"Metadata source must be an owned, unlinked regular file: {path}")
    limit = body_limit(scope(path) if root is None else root)
    if before[3] > limit:
        raise ValueError(f"Metadata source exceeds bound: {path}")
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(fd, "rb") as handle:
        if (os.fstat(handle.fileno()).st_dev, os.fstat(handle.fileno()).st_ino) != before[:2]:
            raise ValueError("Metadata source replaced while opening")
        raw = handle.read(limit + 1)
    after = signature(path)
    if after != before:
        immutable = runtime and not before[2] & 0o222
        protected = before[:5] + before[6:7] + before[8:]
        if not immutable or after[:5] + after[6:7] + after[8:] != protected:
            raise ValueError("Metadata source changed while reading")
        # A different package can link this immutable code inode. Require a
        # second identical body; ordinary sources keep their full signature.
        fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
        with os.fdopen(fd, "rb") as handle:
            if (os.fstat(handle.fileno()).st_dev, os.fstat(handle.fileno()).st_ino) != before[:2]:
                raise ValueError("Metadata source replaced while rereading")
            repeated = handle.read(limit + 1)
        after = signature(path)
        if repeated != raw or after[:5] + after[6:7] + after[8:] != protected:
            raise ValueError("Metadata source changed while rereading")
    if len(raw) != before[3]:
        raise ValueError("Metadata source changed while reading")
    return raw, after


@contextlib.contextmanager
def _source_reads(root, keys):
    """Bound small-record read-ahead; keep SQLite writes on the caller thread."""
    def read(key):
        raw, identity = _read_source(root / key, root=root)
        return key, raw, identity
    if root.name in INDEX_ROOTS or len(keys) < 2:
        # Native index bodies may be 128 MiB; never queue several of them.
        yield map(read, keys)
        return
    def ordered(pool):
        pending = deque()
        for key in keys:
            pending.append(pool.submit(read, key))
            if len(pending) >= 16:
                yield pending.popleft().result()
        while pending:
            yield pending.popleft().result()
    # At most 16 reads/results (8 MiB each) ahead of the current record.
    # Join every reader on both source and caller errors before cleanup.
    with ThreadPoolExecutor(max_workers=16) as pool:
        yield ordered(pool)


def _eligible(root, key, closed):
    if any(part.startswith(".") for part in PurePosixPath(key).parts):
        return False
    if root.name == "runtime":
        from runtime_snapshot_storage import eligible
        return eligible(key)
    if root.name == "observations":
        parts = PurePosixPath(key).parts
        if len(parts) != 2 or re.fullmatch(r"[0-9a-f]{32}", parts[0]) is None or not parts[1].endswith(".json"):
            return False
        if parts[0] not in closed:
            run = root / parts[0] / "run.json"
            raw = _read_source(run, root=root)[0] if run.exists() else archived_bytes(run)
            value = json.loads(raw) if raw is not None else {}
            closed[parts[0]] = (value.get("execution_state") == "exited"
                                and type(value.get("finished_at_ns")) is int)
        return closed[parts[0]]
    if root.name in INDEX_ROOTS:
        return "/" not in key and key.endswith(".json")
    if root.name == "performance":
        return re.fullmatch(r"(?:[1-9][0-9]*|trace-(?:initial|[a-f0-9]{32})-[0-3][a-f0-9])\.jsonl", key) is not None
    if root.name == "resource_metrics":
        return key.endswith((".json", ".jsonl"))
    if root.name == "artifact_provenance":
        return "/" not in key and key.endswith(".json")
    if root.name in {"versions", "version_inventory"}:
        return "/" not in key and key.endswith((".log", ".json"))
    return key.endswith(".paths")  # output-inventory; state and queues excluded


def _closed_observations(root):
    def load(entry):
        if entry.is_symlink():
            raise ValueError("Symlinked metadata directory")
        return entry.name, (_eligible(root, entry.name + "/run.json", {})
                            if entry.is_dir(follow_symlinks=False) else False)
    closed = {}
    # Run bodies are independent, bounded reads; retain a small submission
    # window and join every worker before publication or an exception escapes.
    with os.scandir(root) as entries, ThreadPoolExecutor(max_workers=16) as pool:
        pending = deque()
        for entry in entries:
            if re.fullmatch(r"[0-9a-f]{32}", entry.name) is None:
                continue
            pending.append(pool.submit(load, entry))
            if len(pending) >= 32:
                key, value = pending.popleft().result()
                closed[key] = value
        while pending:
            key, value = pending.popleft().result()
            closed[key] = value
    return closed


def _inventory(root, max_files):
    selected = []
    closed = _closed_observations(root) if root.name == "observations" else {}
    for directory, dirs, files in os.walk(root, followlinks=False):
        dirs[:] = sorted(d for d in dirs if not d.startswith("."))
        for name in dirs:
            if (Path(directory) / name).is_symlink():
                raise ValueError("Symlinked metadata directory")
        if root.name in {"performance", "versions", "version_inventory", "artifact_provenance", *INDEX_ROOTS}:
            dirs.clear()  # Flat record roots never consume nested coordination namespaces.
        for name in sorted(files):
            path = Path(directory) / name
            key = _key(root, path)
            if _eligible(root, key, closed):
                selected.append(key)
                if len(selected) > max_files:
                    raise ValueError("Metadata inventory exceeds bound; narrow the source")
    return sorted(selected)


def _fields(raw):
    try:
        value = json.loads(raw)
        if not isinstance(value, dict):
            return {}
        names = ("workflow", "attempt_id", "family_id", "step", "exit_code", "started_at", "finished_at",
                 "started_at_ns", "finished_at_ns", "runtime_id", "cpus", "memory_gb", "wall_seconds",
                 "cpu_seconds", "peak_process_rss_gb")
        fields = {key: value[key] for key in names if key in value}
        runtime = value.get("runtime", {})
        scheduler = runtime.get("scheduler", {}) if isinstance(runtime, dict) else {}
        if isinstance(scheduler, dict):
            fields["job_id"] = scheduler.get("GG_ARRAY_JOB_ID") or scheduler.get("GG_JOB_ID")
            fields["array_task_id"] = scheduler.get("GG_ARRAY_TASK_ID")
        if isinstance(runtime, dict):
            fields.setdefault("runtime_id", runtime.get("runtime_id"))
        for name in ("started_at", "finished_at"):
            nanos = fields.get(name + "_ns")
            if name not in fields and type(nanos) is int:
                fields[name] = nanos / 1e9
        return fields
    except (ValueError, UnicodeError):
        return {}


def _file_sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def archived_info(path, *, observe=None):
    root = scope(path)
    if root is None:
        return None
    rows = _query(root, "SELECT sha256,bytes,mtime_ns,generation FROM files WHERE path=?", (_key(root, path),), observe)
    if not rows:
        return None
    sha, size, mtime, generation = rows[0]
    if (not isinstance(sha, str) or re.fullmatch(r"[0-9a-f]{64}", sha) is None
            or type(size) is not int or not 0 <= size <= body_limit(root) or type(mtime) is not int
            or type(generation) is not int or generation < 0):
        raise ValueError("Invalid archived metadata identity")
    return {"database": root / DATABASE, "sha256": sha, "size": size, "mtime_ns": mtime,
            "generation": generation}


def maximum_generation(root):
    root = Path(root)
    rows = _query(root, "SELECT coalesce(max(generation),0) FROM files") if root.is_dir() else []
    return int(rows[0][0]) if rows else 0


def cache_policy(root):
    rows = _query(Path(root), "SELECT value FROM metadata WHERE key='cache_limit'")
    limit = int(rows[0][0]) if rows else None
    if limit is not None and not 1 <= limit <= 1024:
        raise ValueError("Invalid version cache retention policy")
    return limit


def read_bytes(path, *, observe=None):
    return _read_source(path)[0] if Path(path).exists() else archived_bytes(path, observe=observe)


def journal_paths(directory):
    archived = archived_children(directory)
    live = list(directory.glob("*.paths"))
    canonical = directory / "inventory.paths"
    # A compacted inventory supersedes old archived append journals.
    if canonical in live or canonical in archived:
        archived = [canonical]
    return sorted(set(live + [p for p in archived if p.suffix == ".paths"]))


def _sync_dir(path):
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def index_members(root):
    raw, _ = _read_source(root.parent / "index-catalog.json")
    value = json.loads(raw)
    names = ([str(bucket) + ".json" for bucket in value["buckets"]] if root.name == "index"
             else list(value["subdir_indexes"].values()))
    if len(names) != len(set(names)) or any(_key(root, root / name) != name or "/" in name for name in names):
        raise ValueError("Invalid native index catalog membership")
    return set(names)


def _retire(root, expected, operation_id):
    if root.name == "runtime":
        from runtime_snapshot_storage import retire
        return retire(root, expected)
    # Replays may remove only the exact files bound before this operation's publication.
    removed = 0
    touched = set()
    try:
        # The exclusive namespace guard excludes participating publications.
        # Pin one immutable DB; retain the before/after identity fence.
        with read_snapshot():
            for key, (before, sha) in expected.items():
                path = root / key
                _key(root, path)
                if not os.path.lexists(path):
                    continue
                raw, current = _read_source(path, root=root)
                if list(current) != list(before) or hashlib.sha256(raw).hexdigest() != sha or archived_bytes(path) != raw:
                    raise ValueError("Source changed or archive cannot restore it; originals retained")
                if signature(path) != current:
                    raise ValueError("Source identity changed before retirement")
                path.unlink()
                touched.add(path.parent)
                removed += 1
    finally:
        # The database is already durable. Flush every changed directory once,
        # including successful removals before a later interruption.
        for directory in sorted(touched):
            _sync_dir(directory)
    return removed


def pack(root, *, scratch=None, delete_originals=False, writers_stopped=False, max_files=500000,
         operation_id="", cache_limit=None):
    if not isinstance(operation_id, str) or len(operation_id) > 128 or (operation_id and re.fullmatch(r"[A-Za-z0-9_.:-]+", operation_id) is None):
        raise ValueError("Invalid metadata operation identity")
    root = Path(root).absolute()
    if scope(root) != root or root.stat().st_uid != os.getuid():
        raise ValueError("Only an owned, registered metadata root can be packed")
    if not writers_stopped:
        raise ValueError("Stop all workspace writers, including older runtimes, before packing")
    with namespace_lock(root / ".gg-metadata.guard", exclusive=True), contextlib.ExitStack() as locks:
        keys = _inventory(root, max_files)
        keep = index_members(root) if root.name in INDEX_ROOTS else None
        if keep is not None and not set(keys).issubset(keep):
            raise ValueError("Native index has unmatched source files")
        generation = 0
        original = _connection(root)
        before_db = None if original is None else original[2]
        if original is not None and root.name == "runtime":
            from runtime_snapshot_storage import verify_database
            verify_database(original[0])
        if original is not None and operation_id:
            saved = original[0].execute("SELECT value FROM metadata WHERE key='operation_id'").fetchone()
            if saved == (operation_id,):
                expected = {key: (json.loads(identity), sha) for key, identity, sha in original[0].execute(
                    "SELECT path,source_identity,sha256 FROM operation_files")}
                if not set(keys).issubset(expected):
                    raise ValueError("Resumed metadata operation has new source members")
                removed = _retire(root, expected, operation_id) if delete_originals else 0
                return {"kind": root.name, "selected_files": len(expected),
                        "archived_files": original[0].execute("SELECT count(*) FROM files").fetchone()[0],
                        "removed_files": len(expected) if delete_originals else 0,
                        "removed_this_call": removed, "database_bytes": (root / DATABASE).stat().st_size,
                        "operation_id": operation_id, "idempotent_replay": True}
        if root.name == "version_inventory" and keys:
            generation = 1 + (original[0].execute("SELECT coalesce(max(generation),0) FROM files").fetchone()[0]
                              if original is not None else 0)
        if cache_limit is not None:
            if root.name != "version_inventory" or type(cache_limit) is not int or not 1 <= cache_limit <= 1024 or keys or original is None:
                raise ValueError("Cache expiry requires a collected version inventory and a bounded limit")
            rows = original[0].execute("SELECT path FROM files ORDER BY generation DESC,mtime_ns DESC,path").fetchall()
            contexts = [key[:-4] for (key,) in rows if re.fullmatch(r"[0-9a-f]{64}\.log", key)]
            expired = set(contexts[cache_limit:])
            keep = {key for (key,) in rows if not (key.endswith((".log", ".json")) and key.rsplit(".", 1)[0] in expired)}
        if original is not None and not keys:
            existing = {row[0] for row in original[0].execute("SELECT path FROM files")}
            prior_limit = original[0].execute("SELECT value FROM metadata WHERE key='cache_limit'").fetchone()
            if (keep is None or existing == keep) and (cache_limit is None or prior_limit == (str(cache_limit),)):
                return {"kind": root.name, "selected_files": 0, "archived_files": len(existing),
                        "removed_files": 0, "database_bytes": (root / DATABASE).stat().st_size}
        if root.name == "artifact_provenance" and keys:
            from gene_family_output_store import (
                GeneFamilyOutputStore,
                all_family_bucket_locks,
                archive_lock,
                producer_quiescence_lock,
            )
            logical = GeneFamilyOutputStore(root.parent)
            locks.enter_context(all_family_bucket_locks(logical.archive_root))
            locks.enter_context(producer_quiescence_lock(logical.archive_root))
            locks.enter_context(archive_lock(logical.archive_root))
            logical.archive_root.mkdir(parents=True, exist_ok=True)
            generation = logical._next_generation()
        with tempfile.TemporaryDirectory(prefix="gg-metadata-", dir=scratch) as work:
            local = Path(work) / "metadata.sqlite3"
            db = sqlite3.connect(local)
            pending = root / (DATABASE + ".partial." + uuid.uuid4().hex)
            snapshots = {}
            try:
                if original is not None:
                    original[0].backup(db)
                else:
                    db.executescript("""
                        PRAGMA journal_mode=DELETE;
                        CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL);
                        CREATE TABLE files(path TEXT PRIMARY KEY,body BLOB NOT NULL,sha256 TEXT NOT NULL,
                            bytes INTEGER NOT NULL,mode INTEGER NOT NULL,mtime_ns INTEGER NOT NULL,
                            fields TEXT NOT NULL,generation INTEGER NOT NULL);
                        CREATE TABLE operation_files(path TEXT PRIMARY KEY,source_identity TEXT NOT NULL,sha256 TEXT NOT NULL);
                        CREATE TABLE records(path TEXT PRIMARY KEY,workflow TEXT,attempt_id TEXT,
                            family_id TEXT,step TEXT,runtime_id TEXT,exit_code INTEGER,
                            started_at REAL,finished_at REAL,cpus INTEGER,memory_gb REAL,
                            wall_seconds REAL,cpu_seconds REAL,peak_process_rss_gb REAL,
                            job_id TEXT,array_task_id TEXT);
                        CREATE INDEX records_family ON records(family_id,step);
                        CREATE INDEX records_workflow ON records(workflow,runtime_id,exit_code);
                        CREATE INDEX records_attempt ON records(attempt_id);
                        CREATE INDEX records_job ON records(job_id,array_task_id);
                        CREATE INDEX records_started ON records(workflow,started_at);
                        CREATE TABLE archive_entries(logical_path TEXT PRIMARY KEY,bucket_path TEXT NOT NULL,
                            subdir TEXT NOT NULL,family_id TEXT,generation INTEGER NOT NULL,
                            zip_path TEXT NOT NULL,zip_location TEXT NOT NULL,member_name TEXT NOT NULL,
                            sha256 TEXT,size INTEGER NOT NULL);
                        CREATE INDEX archive_entries_family ON archive_entries(family_id,subdir);
                    """)
                    db.execute(f"PRAGMA application_id={APPLICATION_ID}")
                    db.execute(f"PRAGMA user_version={SCHEMA_VERSION}")
                    db.execute("INSERT INTO metadata VALUES('kind',?)", (root.name,))
                    if root.name == "runtime":
                        from runtime_snapshot_storage import initialize
                        initialize(db)
                runtime_bodies = ({r[0] for r in db.execute("SELECT sha256 FROM runtime_bodies")}
                                  if root.name == "runtime" else set())
                db.execute("DELETE FROM operation_files")
                with _source_reads(root, keys) as sources:
                    for key, raw, before in sources:
                        snapshots[key] = (before, hashlib.sha256(raw).hexdigest())
                        db.execute("INSERT INTO operation_files VALUES(?,?,?)", (key, json.dumps(before), snapshots[key][1]))
                        if root.name == "runtime":
                            if snapshots[key][1] not in runtime_bodies:
                                db.execute("INSERT INTO runtime_bodies VALUES(?,?,?)",
                                           (snapshots[key][1], zlib.compress(raw, 1), len(raw)))
                                runtime_bodies.add(snapshots[key][1])
                            body = b""
                        else:
                            body = zlib.compress(raw, 1)
                        db.execute("INSERT OR REPLACE INTO files VALUES(?,?,?,?,?,?,?,?)",
                                   (key, body, snapshots[key][1], len(raw),
                                    stat.S_IMODE(before[2]), before[4], json.dumps(_fields(raw), sort_keys=True), generation))
                        if root.name == "index":
                            entries = json.loads(raw)["artifacts"]
                            db.execute("DELETE FROM archive_entries WHERE bucket_path=?", (key,))
                            for logical_path, record in entries.items():
                                db.execute("INSERT INTO archive_entries VALUES(?,?,?,?,?,?,?,?,?,?)",
                                           (logical_path, key, logical_path.split("/", 1)[0], record.get("family_id"),
                                            record["generation"], record["zip_path"], record["zip_location"],
                                            record["member_name"], record.get("sha256"), record["size"]))
                        fields = _fields(raw)
                        names = ("workflow", "attempt_id", "family_id", "step", "runtime_id", "exit_code",
                                 "started_at", "finished_at", "cpus", "memory_gb", "wall_seconds", "cpu_seconds",
                                 "peak_process_rss_gb", "job_id", "array_task_id")
                        values = [fields.get(name) for name in names]
                        if all(v is None or type(v) in (str, int, float) for v in values):
                            db.execute("INSERT OR REPLACE INTO records VALUES(" + ",".join(["?"] * (len(names) + 1)) + ")", (key, *values))
                if keep is not None:
                    for (key,) in db.execute("SELECT path FROM files").fetchall():
                        if key not in keep:
                            db.execute("DELETE FROM files WHERE path=?", (key,))
                            db.execute("DELETE FROM records WHERE path=?", (key,))
                            db.execute("DELETE FROM archive_entries WHERE bucket_path=?", (key,))
                    if {row[0] for row in db.execute("SELECT path FROM files")} != keep:
                        raise ValueError("Native index database membership is incomplete")
                db.execute("INSERT OR REPLACE INTO metadata VALUES('operation_id',?)", (operation_id or uuid.uuid4().hex,))
                if cache_limit is not None:
                    db.execute("INSERT OR REPLACE INTO metadata VALUES('cache_limit',?)", (str(cache_limit),))
                db.execute("INSERT OR REPLACE INTO metadata VALUES('generation',?)", (uuid.uuid4().hex,))
                db.commit()
                if cache_limit is not None:
                    db.execute("VACUUM")
                if db.execute("PRAGMA integrity_check").fetchall() != [("ok",)]:
                    raise ValueError("Metadata database integrity check failed")
                # Full readback, including prior retained generations' active records.
                if root.name == "runtime":
                    from runtime_snapshot_storage import readback_rows
                    bodies = readback_rows(db)
                else:
                    bodies = db.execute("SELECT body,sha256,bytes FROM files")
                for body, sha, size in bodies:
                    decoder = zlib.decompressobj()
                    raw = decoder.decompress(body, body_limit(root) + 1)
                    if (len(raw) != size or not decoder.eof or decoder.unused_data
                            or hashlib.sha256(raw).hexdigest() != sha):
                        raise ValueError("Metadata database readback differs")
                total = db.execute("SELECT count(*) FROM files").fetchone()[0]
                db.close()
                if _inventory(root, max_files) != keys:
                    raise ValueError("Metadata source membership changed")
                with _source_reads(root, snapshots) as sources:
                    for key, raw, current in sources:
                        before, sha = snapshots[key]
                        if current != before or hashlib.sha256(raw).hexdigest() != sha:
                            raise ValueError("Metadata source changed before publication")
                if (None if not (root / DATABASE).exists() else _db_signature(root / DATABASE)) != before_db:
                    raise ValueError("Metadata database changed before publication")
                with local.open("rb") as source, pending.open("xb") as target:
                    os.chmod(pending, 0o600)
                    shutil.copyfileobj(source, target, 1024 * 1024)
                    target.flush()
                    os.fsync(target.fileno())
                if _file_sha(local) != _file_sha(pending):
                    raise ValueError("Published database copy differs")
                os.replace(pending, root / DATABASE)
                _sync_dir(root)
                removed = _retire(root, snapshots, operation_id) if delete_originals else 0
                return {"kind": root.name, "selected_files": len(keys), "archived_files": total,
                        "removed_files": removed, "database_bytes": (root / DATABASE).stat().st_size}
            finally:
                db.close()
                pending.unlink(missing_ok=True)


def restore(root, *, writers_stopped=False):
    root = Path(root).absolute()
    if not writers_stopped:
        raise ValueError("Stop all workspace writers before restoration")
    with namespace_lock(root / ".gg-metadata.guard", exclusive=True), contextlib.ExitStack() as locks:
        logical = None
        if root.name == "artifact_provenance":
            from gene_family_output_store import (
                GeneFamilyOutputStore,
                all_family_bucket_locks,
                archive_lock,
                producer_quiescence_lock,
                read_only_observation,
            )
            logical = GeneFamilyOutputStore(root.parent)
            locks.enter_context(all_family_bucket_locks(logical.archive_root))
            locks.enter_context(producer_quiescence_lock(logical.archive_root))
            locks.enter_context(archive_lock(logical.archive_root))
        rows = _query(root, "SELECT path,mode,mtime_ns FROM files ORDER BY path")
        restored = 0
        for key, _mode, _mtime in rows:
            _key(root, root / key)
        for key, mode, mtime in rows:
            path = root / key
            scope(path)  # Reject redirected ancestors before creating anything.
            if logical is not None and not path.exists():
                with read_only_observation(), logical.read_snapshot():
                    artifact = logical.source_artifact(root.name, path.name)
                if artifact is None or artifact.metadata_path is None:
                    continue  # Preserve native deletion or a newer ZIP generation.
            raw = archived_bytes(path)
            if path.exists() or path.is_symlink():
                if _read_source(path)[0] != raw:
                    raise ValueError(f"Restoration refuses to overwrite changed metadata: {path}")
                continue
            path.parent.mkdir(parents=True, exist_ok=True)
            temp = path.with_name(".restore." + uuid.uuid4().hex)
            try:
                with temp.open("xb") as handle:
                    handle.write(raw)
                    handle.flush()
                    os.fsync(handle.fileno())
                os.chmod(temp, mode)
                os.utime(temp, ns=(mtime, mtime))
                os.link(temp, path, follow_symlinks=False)
                temp.unlink()
                _sync_dir(path.parent)
                if _read_source(path)[0] != raw:
                    raise ValueError("Restored metadata differs")
                restored += 1
            finally:
                temp.unlink(missing_ok=True)
        return {"restored_files": restored, "database_retained": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("pack")
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--scratch", type=Path)
    p.add_argument("--writers-stopped", action="store_true")
    p.add_argument("--delete-originals", action="store_true")
    p.add_argument("--max-files", type=int, default=500000)
    p.add_argument("--operation-id", default="")
    p.add_argument("--cache-limit", type=int)
    p = sub.add_parser("restore")
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--writers-stopped", action="store_true")
    p = sub.add_parser("read")
    p.add_argument("--path", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "pack":
        result = pack(args.root, scratch=args.scratch, delete_originals=args.delete_originals,
                      writers_stopped=args.writers_stopped, max_files=args.max_files,
                      operation_id=args.operation_id, cache_limit=args.cache_limit)
    elif args.command == "restore":
        result = restore(args.root, writers_stopped=args.writers_stopped)
    else:
        raw = _read_source(args.path)[0] if args.path.exists() else archived_bytes(args.path)
        if raw is None:
            raise FileNotFoundError(args.path)
        import sys
        sys.stdout.buffer.write(raw)
        return
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

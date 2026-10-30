"""Journaled offline retirement of redundant metadata directories/markers.

Only zero-byte legacy cache locks and empty archived observation directories
are eligible. Native callers must independently establish stopped writers.
"""

from __future__ import annotations

import argparse
import contextlib
import gzip
import hashlib
import json
import os
import re
import stat
import tempfile
import time
import zlib
from pathlib import Path

import metadata_store as metadata
from shared_namespace_lock import namespace_lock

MAX_ENTRIES = 100000


def _sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _check(path, *, directory=False):
    info = path.lstat()
    if (
        not (stat.S_ISDIR(info.st_mode) if directory else stat.S_ISREG(info.st_mode))
        or info.st_uid != os.getuid()
        or info.st_mode & 0o022
        or not directory
        and info.st_nlink != 1
    ):
        raise ValueError("Unsafe maintenance path: " + str(path))
    return info


def _identity(path):
    info = path.lstat()
    return [
        info.st_dev,
        info.st_ino,
        info.st_mode,
        info.st_uid,
        info.st_gid,
        info.st_nlink,
        info.st_size,
        info.st_mtime_ns,
        info.st_ctime_ns,
    ]


def _stable_directory(path, original):
    info = _check(path, directory=True)
    if [info.st_dev, info.st_ino, info.st_mode, info.st_uid, info.st_gid] != original[:5]:
        raise ValueError("Maintenance directory replaced")
    return info


def _save(path, value):
    raw = (json.dumps(value, sort_keys=True) + "\n").encode()
    if path.exists() or path.is_symlink():
        _check(path)
        if path.read_bytes() != raw:
            raise ValueError("Maintenance journal identity changed")
        return
    descriptor, temporary = tempfile.mkstemp(prefix=".maintenance-", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, path)
        metadata._sync_dir(path.parent)
    finally:
        os.unlink(temporary)


def _root(root, kind):
    root = Path(root).absolute()
    if root.is_symlink() or kind not in {"prune-observations", "retire-version-locks"}:
        raise ValueError("Invalid maintenance root/kind")
    root = root.resolve()
    if (
        root.name != ("observations" if kind == "prune-observations" else "version_inventory")
        or metadata.scope(root) != root
    ):
        raise ValueError("Invalid maintenance root")
    _check(root, directory=True)
    _check(root / metadata.DATABASE)
    metadata._query(root, 'SELECT value FROM metadata WHERE key="kind"')
    return root


def _decoded(body, expected, size):
    if not isinstance(body, bytes) or type(size) is not int or not 0 <= size <= metadata.MAX_BODY:
        raise ValueError("Invalid archived observation body")
    decoder = zlib.decompressobj()
    raw = decoder.decompress(body, size + 1)
    if len(raw) != size or not decoder.eof or decoder.unused_data or hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError("Corrupt archived observation body")
    return json.loads(raw)


def _observation_items(root):
    selected, protected = [], 0
    archived = set()
    physical = {p.name for p in root.iterdir() if re.fullmatch("[a-f0-9]{32}", p.name)}
    rows = metadata._query(root, "SELECT path,body,sha256,bytes FROM files WHERE path LIKE '%/run.json'")
    for key, body, sha, size in rows:
        if re.fullmatch("[a-f0-9]{32}/run[.]json", key) is None:
            raise ValueError("Unexpected archived observation identity")
        name = key.split("/")[0]
        archived.add(name)
        path = root / name
        if not os.path.lexists(path):
            continue
        value = _decoded(body, sha, size)
        closed = (
            isinstance(value, dict)
            and value.get("schema") == "genegalleon-observation-v1"
            and value.get("attempt_id") == name
            and value.get("execution_state") == "exited"
            and type(value.get("finished_at_ns")) is int
            and 0 < value["finished_at_ns"] <= time.time_ns() - 86400 * 10**9
        )
        try:
            _check(path, directory=True)
        except ValueError:
            protected += 1
            continue
        if not closed or any(path.iterdir()):
            protected += 1
            continue
        selected.append({"path": name, "kind": "directory", "identity": _identity(path), "run_sha256": sha})
    return selected, protected + len(physical - archived)


def _cache_items(root):
    names = set()
    for path in root.iterdir():
        match = re.fullmatch(r"([a-f0-9]{64})[.]log[.]lock(?:[.]namespace-v1)?", path.name)
        if match:
            names.add(match[1])
    selected, protected = [], 0
    for name in sorted(names):
        marker = root / (name + ".log.lock")
        namespace = root / (marker.name + ".namespace-v1")
        paths = []
        try:
            if os.path.lexists(marker):
                if _check(marker).st_size != 0:
                    raise ValueError("Legacy cache marker is nonempty")
                paths.append((marker, "file"))
            if os.path.lexists(namespace):
                _check(namespace, directory=True)
                children = list(namespace.iterdir())
                if any(p.name != "readers" for p in children):
                    raise ValueError("Legacy cache owner or unknown entry is present")
                readers = namespace / "readers"
                if os.path.lexists(readers):
                    _check(readers, directory=True)
                    if any(readers.iterdir()):
                        raise ValueError("Legacy cache readers are still present")
                    paths.append((readers, "directory"))
                paths.append((namespace, "directory"))
        except (ValueError, OSError):
            protected += 1
            continue
        selected.extend(
            {"path": str(path.relative_to(root)), "kind": kind, "identity": _identity(path)} for path, kind in paths
        )
    return selected, protected


def _manifest(work):
    path = work / "source.jsonl.gz"
    _check(path)
    before, rows, size = _identity(path), [], 0
    with gzip.open(path, "rb") as stream:
        while line := stream.readline(8193):
            size += len(line)
            if len(line) > 8192 or size > 64 * 1024**2 or len(rows) > MAX_ENTRIES:
                raise ValueError("Maintenance manifest exceeds its bound")
            value = json.loads(line)
            if not isinstance(value, dict):
                raise ValueError("Invalid maintenance manifest row")
            rows.append(value)
    if not rows or _identity(path) != before:
        raise ValueError("Maintenance manifest changed")
    return rows[0], rows[1:]


def _validate_manifest(root, work, kind=None, operation_id=None):
    header, rows = _manifest(work)
    _check(work / "intent.json")
    intent = json.loads((work / "intent.json").read_bytes())
    if (
        intent != {"header": header, "manifest_sha256": _sha(work / "source.jsonl.gz")}
        or header.get("schema_version") != 1
        or header.get("root") != str(root)
        or kind is not None
        and header.get("kind") != kind
        or operation_id is not None
        and header.get("operation_id") != operation_id
    ):
        raise ValueError("Maintenance manifest/operation changed")
    if (
        set(header)
        != {
            "schema_version",
            "root",
            "root_identity",
            "kind",
            "operation_id",
            "database_sha256",
            "selected_entries",
            "protected_items",
        }
        or type(header["schema_version"]) is not int
        or not isinstance(header["root_identity"], list)
        or len(header["root_identity"]) != 9
        or any(type(x) is not int or x < 0 for x in header["root_identity"])
        or header["root_identity"][3] != os.getuid()
        or not isinstance(header["operation_id"], str)
        or re.fullmatch("[A-Za-z0-9_.:-]{1,128}", header["operation_id"]) is None
        or type(header["selected_entries"]) is not int
        or header["selected_entries"] != len(rows)
        or type(header["protected_items"]) is not int
        or header["protected_items"] < 0
        or not isinstance(header["database_sha256"], str)
        or re.fullmatch("[a-f0-9]{64}", header["database_sha256"]) is None
    ):
        raise ValueError("Invalid maintenance journal header")
    _stable_directory(root, header["root_identity"])
    if _sha(root / metadata.DATABASE) != header["database_sha256"]:
        raise ValueError("Maintenance database changed")
    seen = set()
    for row in rows:
        key = row["path"]
        if (
            key in seen
            or row["kind"] not in {"file", "directory"}
            or not isinstance(row["identity"], list)
            or len(row["identity"]) != 9
            or any(type(x) is not int or x < 0 for x in row["identity"])
            or row["identity"][3] != os.getuid()
            or row["identity"][2] & 0o022
            or not (
                stat.S_ISDIR(row["identity"][2])
                if row["kind"] == "directory"
                else stat.S_ISREG(row["identity"][2]) and row["identity"][5] == 1 and row["identity"][6] == 0
            )
        ):
            raise ValueError("Invalid maintenance membership")
        seen.add(key)
        metadata._key(root, root / key)
        if header["kind"] == "prune-observations":
            if row["kind"] != "directory" or re.fullmatch("[a-f0-9]{32}", key) is None:
                raise ValueError("Observation maintenance escaped its namespace")
        elif header["kind"] == "retire-version-locks":
            if re.fullmatch(r"[a-f0-9]{64}[.]log[.]lock(?:[.]namespace-v1(?:/readers)?)?", key) is None or (
                row["kind"] == "file"
            ) != key.endswith(".log.lock"):
                raise ValueError("Cache maintenance escaped its namespace")
        else:
            raise ValueError("Unknown maintenance journal kind")
    return header, rows, intent["manifest_sha256"]


def _prepare(root, work, kind, operation_id):
    if (
        not isinstance(operation_id, str)
        or re.fullmatch("[A-Za-z0-9_.:-]{1,128}", operation_id) is None
        or work == root
        or root in work.parents
    ):
        raise ValueError("Invalid maintenance journal location/identity")
    _check(work.parent, directory=True)
    work.mkdir(mode=0o700, exist_ok=True)
    if _check(work, directory=True).st_mode & 0o077:
        raise ValueError("Maintenance journal is not private")
    if (work / "intent.json").exists():
        return _validate_manifest(root, work, kind, operation_id)
    if os.path.lexists(work / "source.jsonl.gz"):
        raise ValueError("Unsealed maintenance manifest must be preserved for review")
    before = metadata._db_signature(root / metadata.DATABASE)
    with metadata.read_snapshot():
        rows, protected = _observation_items(root) if kind == "prune-observations" else _cache_items(root)
    if len(rows) > MAX_ENTRIES:
        raise ValueError("Maintenance selection exceeds its bound")
    header = {
        "schema_version": 1,
        "root": str(root),
        "root_identity": _identity(root),
        "kind": kind,
        "operation_id": operation_id,
        "database_sha256": _sha(root / metadata.DATABASE),
        "selected_entries": len(rows),
        "protected_items": protected,
    }
    if metadata._db_signature(root / metadata.DATABASE) != before:
        raise ValueError("Maintenance database changed while selecting")
    data = b"".join((json.dumps(row, sort_keys=True) + "\n").encode() for row in [header, *rows])
    with (work / "source.jsonl.gz").open("xb") as stream:
        stream.write(gzip.compress(data, mtime=0))
        stream.flush()
        os.fsync(stream.fileno())
    metadata._sync_dir(work)
    _save(work / "intent.json", {"header": header, "manifest_sha256": _sha(work / "source.jsonl.gz")})
    return _validate_manifest(root, work, kind, operation_id)


def maintain(root, work, *, kind, operation_id, writers_stopped=False):
    if writers_stopped is not True:
        raise ValueError("Stop all workspace writers before maintenance")
    root = _root(root, kind)
    work = Path(work).absolute()
    if work.is_symlink():
        raise ValueError("Symlinked maintenance journal")
    work = work.resolve()
    removed = 0
    with contextlib.ExitStack() as locks:
        if kind == "retire-version-locks":
            locks.enter_context(namespace_lock(root / ".versions-cache.lock", exclusive=True))
        locks.enter_context(namespace_lock(root / ".gg-metadata.guard", exclusive=True))
        header, rows, digest = _prepare(root, work, kind, operation_id)
        try:
            for row in sorted(rows, key=lambda item: (-item["path"].count("/"), item["path"])):
                path = root / row["path"]
                metadata._key(root, path)
                if not os.path.lexists(path):
                    continue  # Only this same durable intent admits a missing source.
                if row["kind"] == "directory":
                    _stable_directory(path, row["identity"])
                    if header["kind"] == "prune-observations" and _identity(path) != row["identity"]:
                        raise ValueError("Observation directory metadata changed")
                    path.rmdir()
                else:
                    _check(path)
                    if _identity(path) != row["identity"] or path.stat().st_size != 0:
                        raise ValueError("Legacy cache marker changed")
                    path.unlink()
                removed += 1
        finally:
            metadata._sync_dir(root)
        result = {
            "operation": kind,
            "operation_id": operation_id,
            "manifest_sha256": digest,
            "database_sha256": header["database_sha256"],
            "removed_entries": len(rows),
            "protected_items": header["protected_items"],
        }
        if not verify(root, work):
            raise ValueError("Maintenance verification failed")
        _save(work / "complete.json", result)
        return {**result, "removed_this_call": removed}


def verify(root, work):
    root, work = Path(root).resolve(), Path(work).resolve()
    _header, rows, _digest = _validate_manifest(root, work)
    return all(not os.path.lexists(root / row["path"]) for row in rows)


def restore(root, work, *, writers_stopped=False):
    if writers_stopped is not True:
        raise ValueError("Stop all workspace writers before restoration")
    root, work = Path(root).resolve(), Path(work).resolve()
    with contextlib.ExitStack() as locks:
        if root.name == "version_inventory":
            locks.enter_context(namespace_lock(root / ".versions-cache.lock", exclusive=True))
        locks.enter_context(namespace_lock(root / ".gg-metadata.guard", exclusive=True))
        header, rows, _digest = _validate_manifest(root, work)
        _root(root, header["kind"])
        _save(work / "restore-intent.json", {"manifest_sha256": _sha(work / "source.jsonl.gz")})
        children = {}
        for row in rows:
            relative = Path(row["path"])
            children.setdefault(relative.parent, set()).add(relative.name)
        for row in sorted(rows, key=lambda item: (item["path"].count("/"), item["path"])):
            path = root / row["path"]
            metadata._key(root, path)
            if row["kind"] == "directory":
                path.mkdir(mode=0o700, exist_ok=True)
                _check(path, directory=True)
                allowed = children.get(Path(row["path"]), set())
                if any(p.name not in allowed for p in path.iterdir()):
                    raise ValueError("Restoration directory contains new data")
            else:
                if not os.path.lexists(path):
                    descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY | os.O_NOFOLLOW, 0o600)
                    os.fsync(descriptor)
                    os.close(descriptor)
                if _check(path).st_size != 0:
                    raise ValueError("Restoration refuses changed marker")
        for row in sorted(rows, key=lambda item: (-item["path"].count("/"), item["path"])):
            path, original = root / row["path"], row["identity"]
            descriptor = os.open(
                path, os.O_RDONLY | os.O_NOFOLLOW | (os.O_DIRECTORY if row["kind"] == "directory" else 0)
            )
            try:
                if os.fstat(descriptor).st_gid != original[4]:
                    os.fchown(descriptor, -1, original[4])
                os.fchmod(descriptor, stat.S_IMODE(original[2]))
                os.utime(descriptor, ns=(original[7], original[7]))
                os.fsync(descriptor)
            finally:
                os.close(descriptor)
        metadata._sync_dir(root)
        return {"restored_entries": len(rows), "database_retained": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prune-observations", "retire-version-locks", "restore", "verify"))
    parser.add_argument("root", type=Path)
    parser.add_argument("--journal", type=Path, required=True)
    parser.add_argument("--operation-id")
    parser.add_argument("--writers-stopped", action="store_true")
    args = parser.parse_args()
    if args.action == "restore":
        result = restore(args.root, args.journal, writers_stopped=args.writers_stopped)
    elif args.action == "verify":
        result = {"verified": verify(args.root, args.journal)}
    else:
        result = maintain(
            args.root,
            args.journal,
            kind=args.action,
            operation_id=args.operation_id,
            writers_stopped=args.writers_stopped,
        )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

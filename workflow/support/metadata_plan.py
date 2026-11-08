"""Read-only, bounded storage plans; a plan never authorizes stopped writers.

Format membership stays with the same implementation used by collection and
maintenance. No journal, lock, directory or SQLite connection with write access
is created. Consumers must revalidate identities and scheduler ownership at the
mutation boundary; a plan is an observation, not a lease.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import stat
import time
from pathlib import Path

import metadata_maintenance as maintenance
import metadata_store as metadata
from shared_namespace_lock import inspect_lock

PUBLIC_STORAGE_PLAN_V1 = True
OPERATIONS = {"collect", "prune-observations", "retire-version-locks"}


def _sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _directory(path):
    info = path.lstat()
    if not stat.S_ISDIR(info.st_mode) or info.st_uid != os.getuid() or info.st_mode & 0o022:
        raise ValueError("Storage plan requires owned physical directories")
    return maintenance._identity(path)


def _inventory(root, limit):
    """Fence every visible source directory, including ineligible live attempts."""
    identities = {}
    files = directories = 0
    for directory, dirs, names in os.walk(root, followlinks=False):
        parent = Path(directory)
        identities[parent] = _directory(parent)
        dirs[:] = sorted(name for name in dirs if not name.startswith("."))
        names = sorted(name for name in names if not name.startswith("."))
        directories += len(dirs)
        files += len(names)
        if files + directories > limit:
            raise ValueError("Storage plan inventory exceeds bound")
        for name in dirs:
            _directory(parent / name)
        for name in names:
            path = parent / name
            info = path.lstat()
            if not stat.S_ISREG(info.st_mode) or info.st_uid != os.getuid() or info.st_mode & 0o022:
                raise ValueError("Storage plan source must be an owned regular file")
            identities[path] = maintenance._identity(path)
    return identities, files, directories


def plan(root, *, operation="collect", max_entries=500000,
         minimum_files=1, minimum_idle_seconds=0, now_ns=None):
    """Return schema v1 counts, generation and exact source evidence.

    Protected counts cover visible source entries outside format membership;
    maintenance's protected_items additionally counts recognized unsafe groups.
    They are directory entries, not deduplicated or quota inode counts.
    """
    if operation not in OPERATIONS:
        raise ValueError("Unknown storage plan operation")
    for value, low, high in ((max_entries, 1, 1000000), (minimum_files, 1, 500000),
                             (minimum_idle_seconds, 0, 604800)):
        if type(value) is not int or not low <= value <= high:
            raise ValueError("Invalid storage plan bound")
    if now_ns is None:
        now_ns = time.time_ns()
    if type(now_ns) is not int or now_ns <= 0:
        raise ValueError("Invalid storage plan time")
    root = Path(os.path.abspath(root))
    if metadata.scope(root) != root or any(p.is_symlink() for p in (root, *root.parents)):
        raise ValueError("Storage plan requires an exact physical metadata root")
    original = _directory(root)
    locks = [root / ".gg-metadata.guard"]
    if root.name == "version_inventory":
        locks.append(root / ".versions-cache.lock")
    states = [inspect_lock(path) for path in locks]
    reasons = []
    if any(s["exclusive"] is not None or s["gate_exists"] or s["shared"] for s in states):
        reasons.append("namespace-owned")
    database = root / metadata.DATABASE
    with metadata.read_snapshot():
        identity = metadata._db_signature(database) if database.exists() else None
        generation = metadata.maximum_generation(root)
        database_sha = _sha(database) if identity is not None else ""
        identities, physical_files, physical_directories = _inventory(root, max_entries)
        if operation == "collect":
            keys = metadata._inventory(root, max_entries)
            selected = [(root / key, "file") for key in keys]
            protected_items = physical_files - len(selected)
        else:
            maintenance._root(root, operation)
            rows, protected_items = (maintenance._observation_items(root) if operation == "prune-observations"
                                     else maintenance._cache_items(root))
            selected = [(root / row["path"], row["kind"]) for row in rows]
        if len(selected) > max_entries:
            raise ValueError("Storage plan selection exceeds bound")
        # Include directory identities for maintenance and file identities for pack.
        evidence = [(str(path.relative_to(root)), kind, maintenance._identity(path)) for path, kind in selected]
        latest = max((max(row[2][7:9]) for row in evidence), default=0)
        if not selected:
            reasons.append("no-eligible-entries")
        elif len(selected) < minimum_files:
            reasons.append("below-minimum-entries")
        if latest and now_ns - latest < minimum_idle_seconds * 10**9:
            reasons.append("sources-too-recent")
        if any(maintenance._identity(path) != before for path, before in identities.items()):
            raise ValueError("Storage plan source changed during observation")
        if (maintenance._identity(root) != original or states != [inspect_lock(path) for path in locks]
                or identity != (metadata._db_signature(database) if database.exists() else None)):
            raise ValueError("Storage plan generation or ownership changed during observation")
    return {
        "schema_version": 1, "root": str(root), "operation": operation,
        "root_identity": original, "database_sha256": database_sha,
        "database_generation": generation, "observed_at_ns": now_ns,
        "latest_source_change_ns": latest,
        "source_identity_sha256": hashlib.sha256(json.dumps(evidence, separators=(",", ":")).encode()).hexdigest(),
        "physical_files": physical_files, "physical_directories": physical_directories,
        "selected_files": sum(kind == "file" for _, kind in selected),
        "selected_directories": sum(kind == "directory" for _, kind in selected),
        "protected_items": protected_items, "eligible": not reasons, "deferred_reasons": reasons,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--operation", choices=sorted(OPERATIONS), default="collect")
    parser.add_argument("--max-entries", type=int, default=500000)
    parser.add_argument("--minimum-files", type=int, default=1)
    parser.add_argument("--minimum-idle-seconds", type=int, default=0)
    args = parser.parse_args()
    print(json.dumps(plan(args.root, operation=args.operation, max_entries=args.max_entries,
                          minimum_files=args.minimum_files, minimum_idle_seconds=args.minimum_idle_seconds), sort_keys=True))


if __name__ == "__main__":
    main()

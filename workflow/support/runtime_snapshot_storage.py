"""Lossless offline runtime-source collection using the existing metadata store.

These records preserve captured scripts for inspection/restoration. They never
select a scientific runtime or confer execution authority.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import stat
import uuid
import zlib
from collections import Counter
from pathlib import Path

import metadata_store as metadata
from shared_namespace_lock import namespace_lock

FILE = re.compile(r"VERSION|gg_common_params\.sh|gg_[a-z0-9_]+_(?:core|entrypoint)\.sh")
RUN = re.compile(r"[A-Za-z0-9_.-]{1,192}")


def eligible(key):
    parts = Path(key).parts
    return (len(parts) == 3 and re.fullmatch(r"gg_[a-z0-9_]+_entrypoint", parts[0]) is not None
            and parts[1] not in {".", ".."} and RUN.fullmatch(parts[1]) is not None
            and FILE.fullmatch(parts[2]) is not None)


def initialize(database):
    database.execute("CREATE TABLE runtime_bodies(sha256 TEXT PRIMARY KEY,body BLOB NOT NULL,bytes INTEGER NOT NULL)")


def readback_rows(database):
    for key, in database.execute("SELECT path FROM files"):
        if not eligible(key):
            raise ValueError("Unexpected archived runtime key")
    if database.execute("SELECT count(*) FROM files f LEFT JOIN runtime_bodies b ON b.sha256=f.sha256 "
                        "WHERE b.sha256 IS NULL OR b.bytes!=f.bytes OR length(f.body)!=0").fetchone()[0]:
        raise ValueError("Runtime body reference is missing or changed")
    return database.execute("SELECT body,sha256,bytes FROM runtime_bodies")


def verify_database(database):
    for body, sha, size in readback_rows(database):
        decoder = zlib.decompressobj()
        raw = decoder.decompress(body, metadata.MAX_BODY + 1)
        if (len(raw) != size or not decoder.eof or decoder.unused_data
                or hashlib.sha256(raw).hexdigest() != sha):
            raise ValueError("Runtime captured body is corrupt")


def retire(root, expected):
    """Retire readonly links while accounting for this operation's own unlinks."""
    missing = Counter()
    for key, (before, _sha) in expected.items():
        if not eligible(key):
            raise ValueError("Runtime snapshot key escaped its registered scope")
        if not os.path.lexists(root / key):
            missing[tuple(before[:2])] += 1
    removed, touched = 0, set()
    try:
        with metadata.read_snapshot():
            for key, (before, sha) in expected.items():
                path = root / key
                metadata.scope(path)
                if not os.path.lexists(path):
                    continue
                raw, current = metadata._read_source(path, root=root)
                stable = (0, 1, 2, 3, 4, 6, 8)
                inode = tuple(before[:2])
                if (any(current[i] != before[i] for i in stable)
                        or current[7] != before[7] - missing[inode]
                        or (not missing[inode] and current[5] != before[5])
                        or hashlib.sha256(raw).hexdigest() != sha or metadata.archived_bytes(path) != raw):
                    raise ValueError("Runtime snapshot changed or archive cannot restore it")
                if metadata.signature(path) != current:
                    raise ValueError("Runtime snapshot changed before unlink")
                path.unlink()
                missing[inode] += 1
                touched.add(path.parent)
                removed += 1
    finally:
        for directory in sorted(touched):
            metadata._sync_dir(directory)
    return removed


def prune(root, *, writers_stopped=False):
    """Remove only empty run directories represented in the verified database."""
    if not writers_stopped:
        raise ValueError("Stop all runtime writers before pruning")
    root = Path(root).absolute()
    if root.name != "runtime" or metadata.scope(root) != root or root.stat().st_uid != os.getuid():
        raise ValueError("Invalid runtime snapshot root")
    removed = 0
    for key, in metadata._query(root, "SELECT path FROM files ORDER BY path"):
        if not eligible(key):
            raise ValueError("Unexpected runtime snapshot record")
    parents = {str(Path(key).parent) for key, in metadata._query(root, "SELECT path FROM files")}
    for relative in sorted(parents):
        path = root / relative
        metadata.scope(path)
        try:
            info = path.lstat()
        except FileNotFoundError:
            continue
        if not stat.S_ISDIR(info.st_mode) or info.st_uid != os.getuid():
            raise ValueError("Runtime snapshot directory changed")
        if any(path.iterdir()):
            continue
        path.rmdir()
        metadata._sync_dir(path.parent)
        removed += 1
    return removed


def prepare(directory, sources):
    """Restore just this captured run after checking the current source bytes."""
    directory = Path(directory).absolute()
    root = directory.parent.parent
    if root.name != "runtime" or not eligible(str(directory.relative_to(root) / "VERSION")):
        raise ValueError("Runtime snapshot directory escaped its scope")
    if metadata.scope(root) != root:
        raise ValueError("Invalid runtime snapshot root")
    with namespace_lock(root / ".gg-metadata.guard", exclusive=True):
        metadata.scope(directory / "VERSION")
        archived = list(metadata.archived_children(directory))
        if not archived:
            return 0
        source_map = {Path(p).name: Path(p) for p in sources if Path(p).is_file()}
        if {p.name for p in archived} != set(source_map) or not all(FILE.fullmatch(n) for n in source_map):
            raise ValueError("Captured runtime source membership changed")
        rows = []
        for path in archived:
            source = source_map[path.name]
            raw = metadata._read_source(source, root=root)[0]
            if raw != metadata.archived_bytes(path):
                raise ValueError("Captured runtime source content changed")
            info = metadata.archived_info(path)
            mode, = metadata._query(root, "SELECT mode FROM files WHERE path=?", (metadata._key(root, path),))[0]
            if type(mode) is not int or not 0 <= mode <= 0o777:
                raise ValueError("Invalid captured runtime mode")
            info["mode"] = mode
            rows.append((path, raw, info))
        restored = 0
        for path, raw, info in rows:
            if os.path.lexists(path):
                if metadata._read_source(path, root=root)[0] != raw:
                    raise ValueError("Live runtime snapshot differs from its archive")
                continue
            path.parent.mkdir(parents=True, exist_ok=True)
            temp = path.with_name(".restore." + uuid.uuid4().hex)
            try:
                with temp.open("xb") as stream:
                    stream.write(raw)
                    stream.flush()
                    os.fsync(stream.fileno())
                os.chmod(temp, info["mode"])
                os.utime(temp, ns=(info["mtime_ns"], info["mtime_ns"]))
                os.link(temp, path, follow_symlinks=False)
                metadata._sync_dir(path.parent)
                restored += 1
            finally:
                temp.unlink(missing_ok=True)
        return restored


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("sources", type=Path, nargs="*")
    parser.add_argument("--collect", action="store_true")
    parser.add_argument("--writers-stopped", action="store_true")
    parser.add_argument("--operation-id", default="")
    args = parser.parse_args()
    if args.collect:
        if args.sources or not args.writers_stopped or not args.operation_id:
            parser.error("Collection requires stopped writers and one exact operation identity")
        value = metadata.pack(args.directory, delete_originals=True, writers_stopped=True,
                              operation_id=args.operation_id)
        value["removed_directories"] = prune(args.directory, writers_stopped=True)
        print(json.dumps(value))
    else:
        prepare(args.directory, args.sources)


if __name__ == "__main__":
    main()

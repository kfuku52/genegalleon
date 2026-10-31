"""Bounded live JSONL stripes and lossless archived performance readers.

The immutable metadata database changes only during an offline collection.
Its generation gives the next live stripes distinct names, preserving history.
"""

from __future__ import annotations

import json
import os
import re
import stat
from pathlib import Path

if __package__:
    from . import metadata_store as metadata
    from .shared_namespace_lock import namespace_lock
else:
    import metadata_store as metadata
    from shared_namespace_lock import namespace_lock

STRIPES = 64


def append(directory, pid, raw):
    if type(pid) is not int or pid < 1 or not isinstance(raw, bytes) or not raw.endswith(b"\n"):
        raise ValueError("Invalid performance log record")
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=True)
    if root.name != "performance" or root.parent.name != "versions":
        # Preserve the documented override's historical per-process layout.
        with (root / f"{pid}.jsonl").open("ab") as stream:
            stream.write(raw)
        return
    info = root.lstat()
    if not stat.S_ISDIR(info.st_mode) or info.st_uid != os.getuid() or info.st_mode & 0o022:
        raise ValueError("Unsafe performance log directory")
    with metadata.writer(root / "performance.jsonl"):
        rows = metadata._query(root, "SELECT value FROM metadata WHERE key='generation'")
        generation = rows[0][0] if rows else "initial"
        if generation != "initial" and (
            not isinstance(generation, str) or re.fullmatch("[a-f0-9]{32}", generation) is None
        ):
            raise ValueError("Invalid performance log generation")
        stripe = pid % STRIPES
        name = f"trace-{generation}-{stripe:02x}.jsonl"
        with namespace_lock(root / f".gg-performance-{stripe:02x}.lock", exclusive=True):
            path = root / name
            descriptor = os.open(path, os.O_WRONLY | os.O_APPEND | os.O_CREAT | os.O_NOFOLLOW, 0o600)
            try:
                before = os.fstat(descriptor)
                if (
                    not stat.S_ISREG(before.st_mode)
                    or before.st_uid != os.getuid()
                    or before.st_nlink != 1
                    or before.st_mode & 0o022
                    or before.st_size + len(raw) > metadata.body_limit(root)
                ):
                    raise ValueError("Unsafe or oversized performance stripe")
                if os.write(descriptor, raw) != len(raw):
                    raise OSError("Incomplete performance log write")
            finally:
                os.close(descriptor)


def records(directory):
    """Yield every original JSONL record across collected and live generations."""
    root = Path(directory)
    with metadata.read_snapshot():
        paths = set(root.glob("*.jsonl")) | set(metadata.archived_children(root))
        for path in sorted(paths):
            if path.suffix != ".jsonl":
                continue
            raw = metadata.read_bytes(path)
            if raw is None:
                raise ValueError("Performance source disappeared")
            for line in raw.splitlines():
                if line.strip():
                    value = json.loads(line)
                    if not isinstance(value, dict):
                        raise ValueError("Performance record must be an object")
                    yield value

"""Retire obsolete per-output publication guards during fenced offline maintenance.

The caller must exclude every old writer before invoking this API. Empty guard
namespaces alone never establish that their primary publication lease is absent.
Current striped guards and scientific outputs are outside this operation.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import stat
import tempfile
from pathlib import Path


class PublicationLockError(RuntimeError):
    pass


_PATTERN = re.compile(r"^\.\.(.+)\.gg-bundle\.lock\.guard(?:\.namespace-v1)?$")


def _identity(path, directory=False):
    info = path.lstat()
    kind = stat.S_ISDIR if directory else stat.S_ISREG
    if not kind(info.st_mode) or info.st_uid != os.getuid():
        raise PublicationLockError(f"Unsafe publication metadata: {path}")
    return [getattr(info, name) for name in (
        "st_dev", "st_ino", "st_mode", "st_uid", "st_gid", "st_nlink",
        "st_size", "st_mtime_ns", "st_ctime_ns")]


def _root(root):
    root = Path(root).absolute()
    if root.resolve(strict=True) != root:
        raise PublicationLockError("Publication root must be physical")
    _identity(root, True)
    return root


def _names(root):
    for directory in sorted(root.iterdir()):
        if directory.name.startswith(".") or not directory.is_dir():
            continue
        _identity(directory, True)
        names = set()
        for path in directory.iterdir():
            if path.name.startswith(".") and path.name.endswith(".gg-bundle.lock"):
                raise PublicationLockError(f"Primary publication lease still exists: {path}")
            if _PATTERN.fullmatch(path.name):
                names.add(path.name.removesuffix(".namespace-v1"))
        for name in sorted(names):
            yield directory.name, name


def _inspect(root, directory, name, *, partial=False):
    if not _PATTERN.fullmatch(name) or name.endswith(".namespace-v1"):
        raise PublicationLockError("Invalid legacy publication guard")
    parent = root / directory
    if partial and not os.path.lexists(parent):
        return [None, None, None]
    _identity(parent, True)
    marker = parent / name
    lease = parent / name[1:-len(".guard")]
    if os.path.lexists(lease):
        raise PublicationLockError(f"Primary publication lease still exists: {lease}")
    namespace = parent / (name + ".namespace-v1")
    paths = (marker, namespace, namespace / "readers")
    identities = []
    for index, path in enumerate(paths):
        if not os.path.lexists(path):
            if not partial:
                raise PublicationLockError(f"Incomplete publication guard: {path}")
            identities.append(None)
            continue
        identity = _identity(path, index > 0)
        if identity[2] & 0o077 or (index == 0 and (identity[5] != 1 or identity[6] != 0)):
            raise PublicationLockError(f"Nonempty or nonprivate publication guard: {path}")
        identities.append(identity)
    if identities[1] is not None and set(p.name for p in namespace.iterdir()) - {"readers"}:
        raise PublicationLockError(f"Publication namespace is occupied: {namespace}")
    if identities[2] is not None and any(paths[2].iterdir()):
        raise PublicationLockError(f"Publication readers still exist: {paths[2]}")
    return identities


def _encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def _records(root):
    yield {"format": "legacy-publication-guards-v1", "root": str(root),
           "identity": _identity(root, True)[:5]}
    for directory, name in _names(root):
        yield {"directory": directory, "name": name,
               "identities": _inspect(root, directory, name)}


def snapshot(root, output=None):
    """Stream exact empty metadata identities; never infer ownership from age."""
    root = _root(root)
    digest, count = hashlib.sha256(), -1
    for row in _records(root):
        data = _encoded(row)
        digest.update(data)
        count += 1
        if output is not None:
            output.write(data)
    return {"groups": count, "entries": count * 3, "sha256": digest.hexdigest()}


def _rows(root, manifest, expected):
    if _identity(manifest)[2] & 0o077:
        raise PublicationLockError("Publication retirement journal must be private")
    digest, count = hashlib.sha256(), 0
    with gzip.open(manifest, "rb") as stream:
        header = next(stream)
        digest.update(header)
        if json.loads(header) != next(_records(root)):
            raise PublicationLockError("Publication root changed")
        for line in stream:
            digest.update(line)
            row = json.loads(line)
            if set(row) != {"directory", "name", "identities"} or row["directory"] in {"", ".", ".."} or Path(row["directory"]).name != row["directory"]:
                raise PublicationLockError("Invalid publication journal row")
            count += 1
            yield row
    if expected != {"groups": count, "entries": count * 3, "sha256": digest.hexdigest()}:
        raise PublicationLockError("Publication retirement journal changed")


def _remaining(root, manifest, expected):
    names = iter(_names(root))
    current = next(names, None)
    for row in _rows(root, manifest, expected):
        key = row["directory"], row["name"]
        if current is not None and current < key:
            raise PublicationLockError("Unrecorded publication guard appeared")
        observed = _inspect(root, *key, partial=True)
        for index, (before, after) in enumerate(zip(row["identities"], observed, strict=True)):
            # Removing readers updates their parent namespace metadata. Its
            # stable object identity still has to match the original journal.
            if after is not None and (after[:5] != before[:5] or (index != 1 and after != before)):
                raise PublicationLockError("Publication guard identity changed")
        if current == key:
            current = next(names, None)
    if current is not None:
        raise PublicationLockError("Unrecorded publication guard appeared")


def retire(root, journal, expected, *, writers_stopped=False):
    """Delete only journaled empty guards after the caller fences all old users.

    The immutable journal permits resuming a partially removed *same* operation.
    It is recovery evidence, never permission to steal a lease or a namespace.
    """
    if writers_stopped is not True:
        raise PublicationLockError("All old publication writers must be fenced")
    root, journal = _root(root), Path(journal)
    _identity(journal, True)
    manifest = journal / "publication-guards-original.jsonl.gz"
    if not manifest.exists():
        fd, temporary = tempfile.mkstemp(prefix=".publication-guards-", dir=journal)
        try:
            with os.fdopen(fd, "wb") as raw:
                with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as output:
                    observed = snapshot(root, output)
                raw.flush()
                os.fsync(raw.fileno())
            if observed != expected:
                raise PublicationLockError("Publication snapshot changed before retirement")
            os.link(temporary, manifest)
        finally:
            Path(temporary).unlink()
    # Validate the whole journal and every remaining group before any deletion.
    _remaining(root, manifest, expected)
    for row in _rows(root, manifest, expected):
        directory, name = row["directory"], row["name"]
        observed = _inspect(root, directory, name, partial=True)
        marker = root / directory / name
        paths = (marker, Path(str(marker) + ".namespace-v1"), Path(str(marker) + ".namespace-v1") / "readers")
        for index in (2, 1, 0):
            if observed[index] is not None:
                current = _identity(paths[index], index > 0)
                expected_identity = row["identities"][index]
                if current[:5] != expected_identity[:5] or (index != 1 and current != expected_identity):
                    raise PublicationLockError("Publication guard changed during retirement")
                paths[index].rmdir() if index else paths[index].unlink()
    return verify_retired(root, journal, expected)


def verify_retired(root, journal, expected):
    root = _root(root)
    manifest = Path(journal) / "publication-guards-original.jsonl.gz"
    for row in _rows(root, manifest, expected):
        if any(_inspect(root, row["directory"], row["name"], partial=True)):
            raise PublicationLockError("Publication retirement is incomplete")
    if next(_names(root), None) is not None:
        raise PublicationLockError("Unrecorded publication guard appeared")
    return {**expected, "retired_entries": expected["entries"], "verified": True}

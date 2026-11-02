"""Retire idle run/family coordination guards during fenced offline maintenance.

The caller must exclude every old writer before invoking this API. Empty guard
namespaces alone never establish that their primary coordination lease is absent.
Stable global gates, publication guards and scientific outputs are outside scope.
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


class CoordinationLockError(RuntimeError):
    pass


def _identity(path, directory=False):
    info = path.lstat()
    kind = stat.S_ISDIR if directory else stat.S_ISREG
    if not kind(info.st_mode) or info.st_uid != os.getuid():
        raise CoordinationLockError(f"Unsafe coordination metadata: {path}")
    return [getattr(info, name) for name in (
        "st_dev", "st_ino", "st_mode", "st_uid", "st_gid", "st_nlink",
        "st_size", "st_mtime_ns", "st_ctime_ns")]


def _root(root):
    root = Path(root).absolute()
    if root.resolve(strict=True) != root:
        raise CoordinationLockError("Coordination root must be physical")
    _identity(root, True)
    return root


def _valid(directory, name):
    if directory == '.gg_run_locks':
        return re.fullmatch(r'[.]family[.][0-9a-f]{64}[.]lock[.]guard', name) is not None
    bucket = re.fullmatch(r'[.]gg_store_locks/family-locks-v2/([0-9a-f]{2})', directory)
    return bool(bucket and re.fullmatch(r'[0-9a-f]{64}[.]lock', name) and name[:2] == bucket[1])


def _names(root):
    run = root / '.gg_run_locks'
    if run.exists():
        _identity(run, True)
        names = set()
        for path in run.iterdir():
            if re.fullmatch(r'family[.][0-9a-f]{64}[.]lock', path.name):
                raise CoordinationLockError(f'Primary producer lease still exists: {path}')
            name = path.name.removesuffix('.namespace-v1')
            if _valid('.gg_run_locks', name):
                names.add(name)
        for name in sorted(names):
            yield '.gg_run_locks', name
    base = root / '.gg_store_locks/family-locks-v2'
    if base.exists():
        _identity(base, True)
        for bucket in sorted(base.iterdir()):
            if not re.fullmatch('[0-9a-f]{2}', bucket.name):
                continue
            _identity(bucket, True)
            directory = '.gg_store_locks/family-locks-v2/' + bucket.name
            names = {p.name.removesuffix('.namespace-v1') for p in bucket.iterdir()}
            for name in sorted(names):
                if _valid(directory, name):
                    yield directory, name


def _inspect(root, directory, name, *, partial=False):
    if not _valid(directory, name):
        raise CoordinationLockError("Invalid legacy coordination guard")
    parent = root / directory
    if partial and not os.path.lexists(parent):
        return [None, None, None]
    _identity(parent, True)
    marker = parent / name
    lease = parent / name[1:-len(".guard")] if directory == ".gg_run_locks" else None
    if lease is not None and os.path.lexists(lease):
        raise CoordinationLockError(f"Primary coordination lease still exists: {lease}")
    namespace = parent / (name + ".namespace-v1")
    paths = (marker, namespace, namespace / "readers")
    identities = []
    for index, path in enumerate(paths):
        if not os.path.lexists(path):
            if not partial:
                raise CoordinationLockError(f"Incomplete coordination guard: {path}")
            identities.append(None)
            continue
        identity = _identity(path, index > 0)
        if identity[2] & 0o077 or (index == 0 and (identity[5] != 1 or identity[6] != 0)):
            raise CoordinationLockError(f"Nonempty or nonprivate coordination guard: {path}")
        identities.append(identity)
    if identities[1] is not None and set(p.name for p in namespace.iterdir()) - {"readers"}:
        raise CoordinationLockError(f"Coordination namespace is occupied: {namespace}")
    if identities[2] is not None and any(paths[2].iterdir()):
        raise CoordinationLockError(f"Coordination readers still exist: {paths[2]}")
    return identities


def _encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def _records(root):
    yield {"format": "legacy-coordination-guards-v1", "root": str(root),
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
        raise CoordinationLockError("Coordination retirement journal must be private")
    digest, count = hashlib.sha256(), 0
    with gzip.open(manifest, "rb") as stream:
        header = next(stream)
        digest.update(header)
        if json.loads(header) != next(_records(root)):
            raise CoordinationLockError("Coordination root changed")
        for line in stream:
            digest.update(line)
            row = json.loads(line)
            if set(row) != {"directory", "name", "identities"} or not _valid(row["directory"], row["name"]):
                raise CoordinationLockError("Invalid coordination journal row")
            count += 1
            yield row
    if expected != {"groups": count, "entries": count * 3, "sha256": digest.hexdigest()}:
        raise CoordinationLockError("Coordination retirement journal changed")


def _remaining(root, manifest, expected):
    names = iter(_names(root))
    current = next(names, None)
    for row in _rows(root, manifest, expected):
        key = row["directory"], row["name"]
        if current is not None and current < key:
            raise CoordinationLockError("Unrecorded coordination guard appeared")
        observed = _inspect(root, *key, partial=True)
        for index, (before, after) in enumerate(zip(row["identities"], observed, strict=True)):
            # Removing readers updates their parent namespace metadata. Its
            # stable object identity still has to match the original journal.
            if after is not None and (after[:5] != before[:5] or (index != 1 and after != before)):
                raise CoordinationLockError("Coordination guard identity changed")
        if current == key:
            current = next(names, None)
    if current is not None:
        raise CoordinationLockError("Unrecorded coordination guard appeared")


def _retire(root, journal, expected, *, writers_stopped=False):
    """Delete only journaled empty guards after the caller fences all old users.

    The immutable journal permits resuming a partially removed *same* operation.
    It is recovery evidence, never permission to steal a lease or a namespace.
    """
    if writers_stopped is not True:
        raise CoordinationLockError("All old coordination writers must be fenced")
    root, journal = _root(root), Path(journal)
    _identity(journal, True)
    manifest = journal / "coordination-guards-original.jsonl.gz"
    if not manifest.exists():
        fd, temporary = tempfile.mkstemp(prefix=".coordination-guards-", dir=journal)
        try:
            with os.fdopen(fd, "wb") as raw:
                with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as output:
                    observed = snapshot(root, output)
                raw.flush()
                os.fsync(raw.fileno())
            if observed != expected:
                raise CoordinationLockError("Coordination snapshot changed before retirement")
            os.link(temporary, manifest)
        finally:
            Path(temporary).unlink()
    _sync(journal)
    # Validate the whole journal and every remaining group before any deletion.
    _remaining(root, manifest, expected)
    touched = set()
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
                    raise CoordinationLockError("Coordination guard changed during retirement")
                paths[index].rmdir() if index else paths[index].unlink()
        touched.add(root / directory)
    for parent in sorted(touched):
        _sync(parent)
    return verify_retired(root, journal, expected)


def verify_retired(root, journal, expected):
    root = _root(root)
    manifest = Path(journal) / "coordination-guards-original.jsonl.gz"
    for row in _rows(root, manifest, expected):
        if any(_inspect(root, row["directory"], row["name"], partial=True)):
            raise CoordinationLockError("Coordination retirement is incomplete")
    if next(_names(root), None) is not None:
        raise CoordinationLockError("Unrecorded coordination guard appeared")
    return {**expected, "retired_entries": expected["entries"], "verified": True}


def _sync(path):
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def retire(root, journal, expected, *, writers_stopped=False):
    """Retire only recorded idle guards; callers must fence all old producer starts.

    Current archive readers/maintenance and every family operation are additionally
    excluded through the original store gates. Live/unknown/stopped owner records
    are never removed, and scientific data and stable global gates stay unchanged.
    """
    if writers_stopped is not True:
        raise CoordinationLockError('All old producer starts must be fenced')
    try:
        from . import gene_family_output_store as store
    except ImportError:
        import gene_family_output_store as store
    root = _root(root)
    archive_root = store._archive_state_root(root)
    if store._read_storage_conversion_marker(root) is not None:
        raise CoordinationLockError('Storage conversion needs to finish first')
    with store.producer_quiescence_lock(archive_root, nonblocking=True) as readers_idle:
        if not readers_idle:
            raise CoordinationLockError('Archive readers are active')
        with store.archive_lock(archive_root, nonblocking=True) as archive_idle:
            if not archive_idle:
                raise CoordinationLockError('Archive maintenance is active')
            with store._bucket_lock(store.family_gate_path(archive_root), exclusive=True,
                                    nonblocking=True) as families_idle:
                if not families_idle:
                    raise CoordinationLockError('Family operations are active')
                return _retire(root, journal, expected, writers_stopped=True)

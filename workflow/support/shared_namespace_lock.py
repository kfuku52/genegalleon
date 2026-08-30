#!/usr/bin/env python3
"""Shared/exclusive locks using atomic shared-filesystem namespace operations.

Unlike flock, mkdir/O_EXCL are coordinated across Lustre localflock clients.
No age-based stealing: a crashed owner leaves a fail-closed lock for explicit
reconciliation after all users of the workspace have stopped. Never mix this
protocol with an older flock-only runtime in the same workspace.
"""

from __future__ import annotations

import argparse
import contextlib
import json
import os
import re
import socket
import time
import uuid
from pathlib import Path


class NamespaceLockError(RuntimeError):
    pass


def _directory(path: Path) -> None:
    if path.is_symlink():
        raise NamespaceLockError(f"Symlinked lock directory: {path}")
    path.mkdir(parents=True, exist_ok=True, mode=0o700)
    if not path.is_dir() or path.is_symlink():
        raise NamespaceLockError(f"Invalid lock directory: {path}")


def _root(path: Path) -> Path:
    _directory(path.parent)
    # Keep the old regular marker for metadata tooling, not for synchronization.
    fd = os.open(path, os.O_CREAT | os.O_WRONLY | os.O_NOFOLLOW, 0o600)
    os.close(fd)
    root = Path(str(path) + ".namespace-v1")
    _directory(root)
    _directory(root / "readers")
    return root


def _write_owner(path: Path, token: str, owner_pid: int) -> None:
    fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY | os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, "w") as handle:
        json.dump({"token": token, "pid": owner_pid, "host": socket.gethostname(),
                   "created_ns": time.time_ns()}, handle)


def _remove_owned(path: Path, token: str) -> None:
    if path.is_symlink() or json.loads(path.read_text())["token"] != token:
        raise NamespaceLockError(f"Lock ownership changed: {path}")
    path.unlink()


def acquire(path: Path, *, exclusive: bool, nonblocking: bool = False,
            timeout: float = 300, owner_pid: int | None = None) -> str | None:
    if timeout < 0:
        raise ValueError("Lock timeout must be nonnegative")
    root = _root(Path(path))
    gate = root / "gate"
    token = uuid.uuid4().hex
    deadline = time.monotonic() + timeout
    while True:
        try:
            gate.mkdir(mode=0o700)
        except FileExistsError as exc:
            if gate.is_symlink() or not gate.is_dir():
                raise NamespaceLockError(f"Invalid lock gate: {gate}") from exc
        else:
            keep_gate = False
            owner_written = False
            try:
                _write_owner(gate / "owner.json", token, owner_pid or os.getpid())
                owner_written = True
                if exclusive:
                    # Do not wait while holding the gate: an existing reader
                    # may need a nested shared lock before it can finish.
                    if not any((root / "readers").iterdir()):
                        keep_gate = True
                        return token
                else:
                    _write_owner(root / "readers" / token, token, owner_pid or os.getpid())
                    return token
            finally:
                if not keep_gate:
                    if owner_written:
                        _remove_owned(gate / "owner.json", token)
                    gate.rmdir()
        if nonblocking:
            return None
        if time.monotonic() >= deadline:
            raise NamespaceLockError(
                f"Timed out waiting for shared-filesystem lock: {path}; "
                "do not remove owner records until all workspace users are stopped")
        time.sleep(min(0.05, max(0, deadline - time.monotonic())))


def release(path: Path, token: str, *, exclusive: bool) -> None:
    if not re.fullmatch(r"[a-f0-9]{32}", token):
        raise NamespaceLockError("Invalid lock ownership token")
    root = Path(str(path) + ".namespace-v1")
    if root.is_symlink() or (root / "readers").is_symlink():
        raise NamespaceLockError(f"Symlinked lock root: {root}")
    if exclusive:
        gate = root / "gate"
        if gate.is_symlink():
            raise NamespaceLockError(f"Symlinked lock gate: {gate}")
        _remove_owned(gate / "owner.json", token)
        gate.rmdir()
    else:
        _remove_owned(root / "readers" / token, token)


@contextlib.contextmanager
def namespace_lock(path: Path, *, exclusive: bool, nonblocking: bool = False,
                   timeout: float = 300):
    token = acquire(path, exclusive=exclusive, nonblocking=nonblocking, timeout=timeout)
    try:
        yield token is not None
    finally:
        if token is not None:
            release(path, token, exclusive=exclusive)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=["acquire-shared", "release-shared"])
    parser.add_argument("path", type=Path)
    parser.add_argument("--owner-pid", type=int)
    parser.add_argument("--token")
    parser.add_argument("--timeout", type=float, default=300)
    args = parser.parse_args()
    if args.operation == "acquire-shared":
        print(acquire(args.path, exclusive=False, owner_pid=args.owner_pid, timeout=args.timeout))
    else:
        release(args.path, args.token or "", exclusive=False)


if __name__ == "__main__":
    main()

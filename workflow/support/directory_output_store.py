#!/usr/bin/env python3
"""Manage retained results and closed scheduler-log batches as verified ZIP64 files."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path, PurePosixPath

try:
    from . import orthofinder_output_store as STORE
    from .directory_storage_policy import DIRECTORY_SCOPES
except ImportError:
    import orthofinder_output_store as STORE
    from directory_storage_policy import DIRECTORY_SCOPES


def convert(root, *, kind, **options):
    return STORE.convert(root, namespace=kind, **options)


def open_binary(root, relative, *, kind):
    return STORE.open_binary(root, relative, namespace=kind)


def materialized_directory(root, name, scratch, *, kind):
    return STORE.materialized_directory(root, name, scratch, namespace=kind)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=("pack", "materialize", "verify", "status", "read"))
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--kind", required=True, choices=DIRECTORY_SCOPES)
    parser.add_argument("--directory", action="append")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--compression", choices=("adaptive", "deflate", "store"), default="adaptive")
    parser.add_argument("--compression-level", type=int, default=6)
    parser.add_argument("--available-bytes", type=int)
    parser.add_argument("--writers-stopped", action="store_true",
                        help="Scheduler log pack/materialize only: all jobs and output delivery for this batch have stopped.")
    parser.add_argument("--relative-path", help="Logical member to stream for read, such as scheduler_output/task.err.")
    parser.add_argument("--producer-lock-token", help=argparse.SUPPRESS)
    args = parser.parse_args()
    names = args.directory or DIRECTORY_SCOPES[args.kind]
    if not all(name in DIRECTORY_SCOPES[args.kind] for name in names):
        parser.error("directory is outside the selected output kind")
    if args.operation == "read":
        relative = PurePosixPath(args.relative_path or "")
        if (not args.relative_path or relative.is_absolute() or ".." in relative.parts
                or not any(PurePosixPath(name) in relative.parents for name in names)):
            parser.error("read requires --relative-path inside the selected output kind")
        with open_binary(args.root, relative, kind=args.kind) as handle:
            shutil.copyfileobj(handle, sys.stdout.buffer, length=1024 * 1024)
        return
    if args.relative_path is not None:
        parser.error("--relative-path is valid only for read")
    if args.operation in ("pack", "materialize"):
        result = convert(
            args.root, kind=args.kind, to="zip" if args.operation == "pack" else "raw", directories=names,
            workers=args.workers, compression=args.compression, compression_level=args.compression_level,
            available_bytes=args.available_bytes, producer_lock_token=args.producer_lock_token,
            writers_stopped=args.writers_stopped,
        )
    elif args.operation == "verify":
        result = {"directories": [STORE.archive_records(args.root, name, deep=True) for name in names
                                  if (args.root / (name + ".zip")).exists()]}
    else:
        result = {"directories": STORE.status(args.root, names)}
    print(json.dumps(result))


if __name__ == "__main__":
    main()

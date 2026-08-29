#!/usr/bin/env python3
"""Verify a ZIP and atomically publish it at a final destination."""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import sys
import uuid
import zipfile
from pathlib import Path
from typing import Sequence

try:
    from .safe_zip_extract import validated_members
except ImportError:  # pragma: no cover - direct script execution
    from safe_zip_extract import validated_members


class ZipPublishError(RuntimeError):
    """Raised when a source ZIP cannot be safely published."""


def _fsync_directory(path: Path) -> None:
    try:
        descriptor = os.open(path, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(descriptor)
    except OSError:
        pass
    finally:
        os.close(descriptor)


def _signature(path: Path) -> tuple[int, int, int, int, int]:
    value = path.stat()
    return (
        int(value.st_dev),
        int(value.st_ino),
        int(value.st_size),
        int(value.st_mtime_ns),
        int(value.st_ctime_ns),
    )


def _matches_source(
    path: Path,
    signature: tuple[int, int, int, int, int],
    sha256: str,
) -> bool:
    if path.is_symlink():
        return False
    try:
        if _signature(path) != signature:
            return False
    except FileNotFoundError:
        return False
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest() == sha256 and _signature(path) == signature
    except FileNotFoundError:
        return False


def publish_zip(
    source: Path | str,
    destination: Path | str,
    *,
    remove_source: bool = False,
    expected_prefix: str | None = None,
) -> Path:
    source = Path(source)
    destination = Path(destination)
    if source.is_symlink() or not source.is_file():
        raise ZipPublishError(f"Source ZIP is not a regular file: {source}")
    if destination.is_symlink():
        raise ZipPublishError(f"Symlinked ZIP destinations are unsupported: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    source_location = source.resolve(strict=True)
    destination_location = destination.resolve(strict=False)
    source_signature = _signature(source)
    temporary = destination.parent / (
        f".{destination.name}.partial.{os.getpid()}.{uuid.uuid4().hex}"
    )
    try:
        copied_digest = hashlib.sha256()
        with source.open("rb") as input_handle, temporary.open("xb") as output_handle:
            class HashingDestination:
                def write(self, chunk: bytes) -> int:
                    written = output_handle.write(chunk)
                    copied_digest.update(memoryview(chunk)[:written])
                    return int(written)

            shutil.copyfileobj(
                input_handle,
                HashingDestination(),
                length=1024 * 1024,
            )
            output_handle.flush()
            os.fsync(output_handle.fileno())
        try:
            source_metadata_changed = (
                source.is_symlink() or _signature(source) != source_signature
            )
        except FileNotFoundError:
            source_metadata_changed = True
        if source_metadata_changed:
            raise ZipPublishError(f"Source ZIP changed while being published: {source}")
        try:
            with zipfile.ZipFile(temporary, "r") as archive:
                if expected_prefix is not None:
                    validated_members(archive, source, expected_prefix)
                bad_member = archive.testzip()
        except (OSError, RuntimeError, zipfile.BadZipFile) as exc:
            raise ZipPublishError(f"Failed to verify source ZIP {source}: {exc}") from exc
        if bad_member is not None:
            raise ZipPublishError(
                f"CRC verification failed in source ZIP {source}: {bad_member}"
            )
        if not _matches_source(
            source,
            source_signature,
            copied_digest.hexdigest(),
        ):
            raise ZipPublishError(f"Source ZIP changed while being published: {source}")
        os.replace(temporary, destination)
        _fsync_directory(destination.parent)
        if remove_source and source_location != destination_location:
            source.unlink()
            _fsync_directory(source.parent)
        return destination
    finally:
        temporary.unlink(missing_ok=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--destination", required=True, type=Path)
    parser.add_argument("--remove-source", action="store_true")
    parser.add_argument("--expected-prefix")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        path = publish_zip(
            args.source,
            args.destination,
            remove_source=args.remove_source,
            expected_prefix=args.expected_prefix,
        )
    except ZipPublishError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

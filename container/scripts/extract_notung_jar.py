#!/usr/bin/env python3
"""Safely select and atomically extract the Notung JAR from its release ZIP."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import stat
import sys
import uuid
import zipfile
from pathlib import Path, PurePosixPath

MAX_NOTUNG_JAR_BYTES = 256 * 1024 * 1024


class NotungArchiveError(RuntimeError):
    """Raised when a Notung release archive is unsafe or incomplete."""


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


def _natural_key(value: str) -> tuple[object, ...]:
    return tuple(
        int(token) if token.isdigit() else token.lower()
        for token in re.split(r"([0-9]+)", value)
    )


def _safe_member(info: zipfile.ZipInfo) -> bool:
    member = PurePosixPath(info.filename)
    canonical = member.as_posix()
    mode = info.external_attr >> 16
    file_type = stat.S_IFMT(mode)
    return (
        not info.is_dir()
        and not member.is_absolute()
        and canonical == info.filename
        and "." not in member.parts
        and ".." not in member.parts
        and "\\" not in info.filename
        and not stat.S_ISLNK(mode)
        and file_type in {0, stat.S_IFREG}
    )


def extract_notung_jar(archive_path: Path | str, destination: Path | str) -> Path:
    archive_path = Path(archive_path)
    destination = Path(destination)
    if archive_path.is_symlink() or not archive_path.is_file():
        raise NotungArchiveError(f"Notung ZIP is not a regular file: {archive_path}")
    if destination.is_symlink():
        raise NotungArchiveError(f"Symlinked Notung destination is unsupported: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.parent / (
        f".{destination.name}.partial.{os.getpid()}.{uuid.uuid4().hex}"
    )
    try:
        try:
            with zipfile.ZipFile(archive_path, "r") as archive:
                candidates = [
                    info
                    for info in archive.infolist()
                    if _safe_member(info)
                    and re.fullmatch(r"Notung-2[.]9.*[.]jar", PurePosixPath(info.filename).name)
                ]
                if not candidates:
                    raise NotungArchiveError(
                        f"Notung-2.9*.jar was not found in {archive_path}"
                    )
                selected = sorted(
                    candidates,
                    key=lambda info: _natural_key(PurePosixPath(info.filename).name),
                )[-1]
                if selected.file_size <= 0 or selected.file_size > MAX_NOTUNG_JAR_BYTES:
                    raise NotungArchiveError(
                        f"Notung JAR has an implausible size ({selected.file_size} bytes): "
                        f"{selected.filename}"
                    )
                with archive.open(selected, "r") as source, temporary.open("xb") as output:
                    shutil.copyfileobj(source, output, length=1024 * 1024)
                    output.flush()
                    os.fsync(output.fileno())
        except (OSError, RuntimeError, zipfile.BadZipFile) as exc:
            if isinstance(exc, NotungArchiveError):
                raise
            raise NotungArchiveError(f"Failed to extract {archive_path}: {exc}") from exc
        os.chmod(temporary, 0o755)
        os.replace(temporary, destination)
        _fsync_directory(destination.parent)
        return destination
    finally:
        temporary.unlink(missing_ok=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--destination", required=True, type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        extract_notung_jar(args.archive, args.destination)
    except NotungArchiveError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Read or update GeneGalleon's semantic VERSION marker."""

from __future__ import annotations

import argparse
import os
import re
import tempfile
from pathlib import Path

SEMVER_RE = re.compile(r"^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)$")
DEFAULT_VERSION_FILE = Path(__file__).resolve().parents[2] / "VERSION"


def parse_version(value: str) -> tuple[int, int, int]:
    match = SEMVER_RE.fullmatch(value.strip())
    if match is None:
        raise ValueError(f"VERSION must be MAJOR.MINOR.PATCH, found: {value!r}")
    return tuple(int(part) for part in match.groups())


def next_version(current: tuple[int, int, int], component: str) -> tuple[int, int, int]:
    major, minor, patch = current
    if component == "major":
        return major + 1, 0, 0
    if component == "minor":
        return major, minor + 1, 0
    return major, minor, patch + 1


def write_atomic(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(value + "\n")
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("component", choices=("major", "minor", "patch"))
    parser.add_argument("--dry-run", action="store_true", help="Print the next version without changing VERSION.")
    parser.add_argument("--version-file", type=Path, default=DEFAULT_VERSION_FILE, help=argparse.SUPPRESS)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    current = parse_version(args.version_file.read_text(encoding="utf-8"))
    updated = ".".join(str(part) for part in next_version(current, args.component))
    if not args.dry_run:
        write_atomic(args.version_file, updated)
    print(updated)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

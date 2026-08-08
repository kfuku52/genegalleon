#!/usr/bin/env python3
"""Validate that all CSUBST and stat_branch IDs identify the same clades."""

from __future__ import annotations

import argparse
import tempfile
import zipfile
from pathlib import Path, PurePosixPath

from csubst_site_wrapper import validate_all_csubst_stat_branch_identity


def safe_extract_iqtree_archive(archive_path: Path, destination: Path) -> Path:
    with zipfile.ZipFile(archive_path, "r") as archive:
        bad_member = archive.testzip()
        if bad_member is not None:
            raise ValueError(f"Corrupt IQ-TREE ancestral-state ZIP member: {bad_member}")
        for member in archive.infolist():
            pure_name = PurePosixPath(member.filename)
            if pure_name.is_absolute() or ".." in pure_name.parts:
                raise ValueError(f"Unsafe IQ-TREE ancestral-state ZIP member: {member.filename}")
        archive.extractall(destination)
    candidates = sorted(path.parent for path in destination.rglob("csubst.nwk"))
    if len(candidates) != 1:
        raise ValueError(f"Expected one csubst.nwk in {archive_path}, found {len(candidates)}")
    return candidates[0]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iqtree-anc", required=True, type=Path)
    parser.add_argument("--stat-branch", required=True, type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    with tempfile.TemporaryDirectory(prefix="gg-csubst-branch-check-") as temporary:
        iqtree_anc_dir = safe_extract_iqtree_archive(args.iqtree_anc, Path(temporary) / "iqtree_anc")
        validate_all_csubst_stat_branch_identity(str(args.stat_branch), str(iqtree_anc_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

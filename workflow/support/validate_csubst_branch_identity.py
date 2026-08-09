#!/usr/bin/env python3
"""Validate that all CSUBST and stat_branch IDs identify the same clades."""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

from csubst_site_wrapper import validate_all_csubst_stat_branch_identity
from safe_zip_extract import extract_expected_prefix


def safe_extract_iqtree_archive(archive_path: Path, destination: Path) -> Path:
    archive_path = Path(archive_path)
    stem = archive_path.stem
    if stem.endswith("_iqtree.anc"):
        expected_prefix = f"{stem.removesuffix('_iqtree.anc')}.iqtree.anc"
    elif stem.endswith(".iqtree.anc"):
        expected_prefix = stem
    else:
        raise ValueError(f"Unexpected IQ-TREE ancestral-state ZIP name: {archive_path}")
    extracted = extract_expected_prefix(
        archive_path,
        destination,
        expected_prefix,
    )
    csubst_tree = extracted / "csubst.nwk"
    if not csubst_tree.is_file():
        raise ValueError(f"Expected {csubst_tree} in {archive_path}")
    return extracted


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

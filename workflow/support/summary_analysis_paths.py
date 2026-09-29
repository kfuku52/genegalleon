#!/usr/bin/env python3
"""Recover analysis input paths selected when a gene-family summary was recorded."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path, PurePosixPath

SLOTS = {
    "input_2": "untrimmed_aln",
    "input_3": "trimmed_aln",
    "input_4": "unrooted_tree",
    "input_5": "rooted_tree",
    "input_9": "dated_tree",
}


def recorded_analysis_paths(manifest: Path, logical_root: Path) -> dict[str, Path]:
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    if payload.get("step") != "summary_statistics" or not isinstance(payload.get("inputs"), list):
        raise ValueError("summary provenance manifest is invalid")
    root = logical_root.resolve(strict=True)
    result: dict[str, Path] = {}
    for item in payload["inputs"]:
        label = item.get("label") if isinstance(item, dict) else None
        if label not in SLOTS:
            continue
        if label in result or item.get("scope") != "logical" or item.get("artifact_type") != "file":
            raise ValueError(f"invalid recorded analysis input: {label}")
        relative = PurePosixPath(item.get("path", ""))
        if relative.is_absolute() or not relative.parts or any(part in {"", ".", ".."} for part in relative.parts):
            raise ValueError(f"unsafe recorded analysis path: {label}")
        path = (root / Path(*relative.parts)).resolve(strict=True)
        if not path.is_relative_to(root) or not path.is_file():
            raise ValueError(f"recorded analysis input is unavailable: {label}")
        result[label] = path
    if "input_2" not in result or "input_3" not in result:
        raise ValueError("summary provenance lacks its analysis alignments")
    return {SLOTS[label]: path for label, path in result.items()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--logical-root", type=Path, required=True)
    args = parser.parse_args()
    try:
        for slot, path in recorded_analysis_paths(args.manifest, args.logical_root).items():
            print(f"{slot}\t{path}")
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"Cannot recover recorded summary analysis inputs: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

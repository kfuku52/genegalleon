#!/usr/bin/env python3
"""Safely promote provenance-verified input-generation BUSCO outputs."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from content_digest_cache import cached_sha256_file
from content_digest_cache import configure as configure_digest_cache


class ReuseError(RuntimeError):
    pass


def entry_by_label(payload: dict[str, object], collection: str, label: str) -> dict[str, object]:
    entries = payload.get(collection)
    if not isinstance(entries, list):
        raise ReuseError(f"Manifest {collection} is not a list")
    matches = [entry for entry in entries if isinstance(entry, dict) and entry.get("label") == label]
    if len(matches) != 1:
        raise ReuseError(f"Manifest does not contain exactly one {collection} entry named {label}")
    return matches[0]


def verify_file(path: Path, entry: dict[str, object]) -> None:
    if not path.is_file() or path.is_symlink():
        raise ReuseError(f"BUSCO reuse source is missing or unsafe: {path}")
    digest, size = cached_sha256_file(path)
    if digest != entry.get("sha256") or size != entry.get("size_bytes"):
        raise ReuseError(f"BUSCO reuse source differs from its provenance manifest: {path}")


def copy_atomic(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{destination.name}.", dir=destination.parent)
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        shutil.copyfile(source, temporary)
        shutil.copymode(source, temporary)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def run(args: argparse.Namespace) -> int:
    payload = json.loads(args.manifest.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1 or payload.get("step") != "input_generation_species_busco":
        raise ReuseError("Unsupported input-generation BUSCO provenance manifest")
    if str(payload.get("family_id")) != args.species:
        raise ReuseError("BUSCO manifest species does not match the requested species")
    parameters = payload.get("parameters")
    if not isinstance(parameters, dict):
        raise ReuseError("BUSCO manifest parameters are malformed")
    required_parameters = {
        "busco_mode": args.mode,
        "busco_lineage_resolved": args.lineage,
        "evalue": "1e-03",
        "limit": "20",
    }
    for key, expected in required_parameters.items():
        if str(parameters.get(key, "")) != expected:
            raise ReuseError(
                f"BUSCO manifest parameter mismatch for {key}: expected {expected!r}, found {parameters.get(key)!r}"
            )
    verify_file(args.input, entry_by_label(payload, "inputs", "species_cds"))
    verify_file(args.full_source, entry_by_label(payload, "outputs", "busco_full"))
    verify_file(args.short_source, entry_by_label(payload, "outputs", "busco_short"))
    copy_atomic(args.full_source, args.full_destination)
    copy_atomic(args.short_source, args.short_destination)
    print(f"Reused input-generation BUSCO outputs for {args.species}")
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--manifest", required=True, type=Path)
    result.add_argument("--species", required=True)
    result.add_argument("--input", required=True, type=Path)
    result.add_argument("--lineage", required=True)
    result.add_argument("--mode", default="transcriptome")
    result.add_argument("--full-source", required=True, type=Path)
    result.add_argument("--short-source", required=True, type=Path)
    result.add_argument("--full-destination", required=True, type=Path)
    result.add_argument("--short-destination", required=True, type=Path)
    result.add_argument("--digest-cache", type=Path)
    return result


def main() -> int:
    args = parser().parse_args()
    try:
        if args.digest_cache:
            configure_digest_cache(args.digest_cache)
        return run(args)
    except (OSError, ValueError, TypeError, json.JSONDecodeError, ReuseError) as exc:
        print(f"BUSCO output is not reusable: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

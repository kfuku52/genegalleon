#!/usr/bin/env python3
"""Publish verified per-species quantification ZIPs and stage merge inputs."""

from __future__ import annotations

import argparse
import os
import shutil
import stat
import sys
import tempfile
import zipfile
from pathlib import Path

from artifact_provenance import (
    SCHEMA_VERSION,
    entries_by_label,
    load_manifest,
    raw_or_zip_directory_digest,
    zip_directory_digest,
)
from directory_storage_policy import is_transcriptome_quant_species
from filesystem_io import fsync_directory
from safe_zip_extract import validated_members
from validate_transcriptome_quant_outputs import selected_runs


def inventory(source: Path) -> dict[Path, os.stat_result]:
    if source.is_symlink() or not source.is_dir():
        raise ValueError(f"Quantification source must be a regular directory: {source}")
    members = {}
    for item in sorted(source.rglob("*")):
        status = item.lstat()
        if stat.S_ISREG(status.st_mode):
            members[item] = status
        elif not stat.S_ISDIR(status.st_mode):
            raise ValueError(f"Unexpected quantification member type: {item}")
    if not members:
        raise ValueError(f"No quantification files found: {source}")
    return members


def verify_inventory(source: Path, before: dict[Path, os.stat_result]) -> None:
    def signatures(members):
        return {item: (value.st_dev, value.st_ino, value.st_mode, value.st_size,
                       value.st_mtime_ns, value.st_ctime_ns) for item, value in members.items()}

    if signatures(inventory(source)) != signatures(before):
        raise ValueError(f"Quantification files changed during packing: {source}")


def recorded_identity(manifest: Path, raw: Path) -> tuple[str, int]:
    recorded = load_manifest(manifest)
    entry = entries_by_label(recorded, "outputs").get("quant_directory", {})
    if (recorded["schema_version"] != SCHEMA_VERSION or recorded.get("step") != "transcriptome_quant"
            or recorded.get("family_id") != raw.name or entry.get("artifact_type") != "directory"
            or Path(entry.get("path", "")).parts[-4:] != raw.parts[-4:]):
        raise ValueError(f"Manifest does not describe this quantification directory: {manifest}")
    return entry["sha256"], entry["size_bytes"]


def archive_path(raw: Path) -> Path:
    if not is_transcriptome_quant_species(raw) or raw.is_symlink() or raw.parent.is_symlink():
        raise ValueError(f"Unexpected quantification species path: {raw}")
    archive = raw.with_name(raw.name + ".zip")
    if archive.is_symlink() or (archive.exists() and not archive.is_file()):
        raise ValueError(f"Invalid quantification ZIP path: {archive}")
    return archive


def pack_quant(raw: Path, source: Path, manifest: Path | None = None) -> Path:
    """Caller holds the species lock for publication and retirement."""
    archive = archive_path(raw)
    retire = source == raw
    if retire != (manifest is not None) or (not retire and raw.exists()):
        raise ValueError("Retiring a published directory requires its original quantification manifest")
    expected = recorded_identity(manifest, raw) if retire else None
    before = inventory(source)
    raw.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{raw.name}.packing-", dir=raw.parent) as work:
        if retire and archive.exists():
            identity = raw_or_zip_directory_digest(raw)
            if identity[:2] != expected:
                raise ValueError(f"Existing quantification ZIP differs from recorded outputs: {archive}")
            verify_inventory(source, before)
        else:
            partial = Path(work) / archive.name
            with zipfile.ZipFile(partial, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=1) as output:
                for item in before:
                    output.write(item, raw.name + "/" + item.relative_to(source).as_posix(),
                                 compress_type=zipfile.ZIP_STORED if item.suffix == ".h5" else None)
            identity = zip_directory_digest(partial, raw.name, by_parts=True)
            verify_inventory(source, before)
            if expected is not None and identity[:2] != expected:
                raise ValueError(f"Quantification ZIP differs from recorded outputs: {raw}")
            with partial.open("rb") as handle:
                os.fsync(handle.fileno())
            os.replace(partial, archive)
            fsync_directory(raw.parent)
        if retire:
            # A rename avoids exposing a partially removed directory after interruption.
            os.replace(raw, Path(work) / "retired")
            fsync_directory(raw.parent)
    return archive


def stage_merge(raw: Path, metadata: Path, destination: Path) -> None:
    """Extract the two files per selected run that amalgkit merge reads."""
    archive_pathname = archive_path(raw)
    runs = selected_runs(metadata)
    with zipfile.ZipFile(archive_pathname) as archive:
        members, _ = validated_members(archive, archive_pathname, raw.name)
        present = {info.filename: info for info in members}
        with tempfile.TemporaryDirectory(prefix=".quant-merge-", dir=destination.parent) as work:
            stage = Path(work) / "quant"
            stage.mkdir()
            for run in runs:
                run_dir = stage / run
                run_dir.mkdir()
                for suffix in ("abundance.tsv", "run_info.json"):
                    name = f"{run}_{suffix}"
                    member = f"{raw.name}/{run}/{name}"
                    if member not in present or present[member].file_size == 0:
                        raise ValueError(f"Missing or empty merge input in {archive_pathname}: {member}")
                    with archive.open(present[member]) as source, (run_dir / name).open("xb") as output:
                        shutil.copyfileobj(source, output)
            os.rename(stage, destination)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    pack = commands.add_parser("pack")
    pack.add_argument("--raw", type=Path, required=True)
    pack.add_argument("--source", type=Path, required=True)
    pack.add_argument("--manifest", type=Path)
    merge = commands.add_parser("stage-merge")
    merge.add_argument("--raw", type=Path, required=True)
    merge.add_argument("--metadata", type=Path, required=True)
    merge.add_argument("--destination", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "pack":
            print(f"Published verified quantification ZIP: {pack_quant(args.raw, args.source, args.manifest)}")
        else:
            stage_merge(args.raw, args.metadata, args.destination)
    except (OSError, RuntimeError, ValueError, KeyError, zipfile.BadZipFile) as exc:
        print(f"Quantification ZIP error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Materialize generated workspace inputs without unsafe directory symlinks."""

from __future__ import annotations

import argparse
import datetime as dt
import errno
import hashlib
import json
import os
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path

SUPPORTED_INPUTS = ("species_cds", "species_cds_fx2tab", "species_gff", "species_genome")


class MaterializeError(RuntimeError):
    """Raised when an input cannot be materialized without changing meaning."""


@dataclass(frozen=True)
class TreeSummary:
    sha256: str
    files: int
    bytes: int


@dataclass
class PreparedInput:
    name: str
    source: Path
    destination: Path
    original_link: str
    stage: Path
    summary: TreeSummary
    hardlinks: int
    copies: int


def _is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def summarize_tree(root: Path) -> TreeSummary:
    if root.is_symlink() or not root.is_dir():
        raise MaterializeError(f"Expected a real directory: {root}")
    digest = hashlib.sha256()
    file_count = 0
    total_bytes = 0
    for current, directories, files in os.walk(root):
        directories.sort()
        files.sort()
        current_path = Path(current)
        for directory in directories:
            child = current_path / directory
            if child.is_symlink():
                raise MaterializeError(f"Symlinked source directory member is unsupported: {child}")
        for filename in files:
            child = current_path / filename
            if child.is_symlink() or not child.is_file():
                raise MaterializeError(f"Non-regular source member is unsupported: {child}")
            relative = child.relative_to(root).as_posix().encode("utf-8")
            child_digest = hashlib.sha256()
            size = 0
            with child.open("rb") as handle:
                while chunk := handle.read(1024 * 1024):
                    child_digest.update(chunk)
                    size += len(chunk)
            digest.update(len(relative).to_bytes(8, "big"))
            digest.update(relative)
            digest.update(child_digest.digest())
            digest.update(size.to_bytes(8, "big"))
            file_count += 1
            total_bytes += size
    if not file_count:
        raise MaterializeError(f"Source directory contains no regular files: {root}")
    return TreeSummary(digest.hexdigest(), file_count, total_bytes)


def _clone_tree(source: Path, stage: Path) -> tuple[int, int]:
    hardlinks = 0
    copies = 0
    for current, directories, files in os.walk(source):
        directories.sort()
        files.sort()
        current_path = Path(current)
        relative_dir = current_path.relative_to(source)
        output_dir = stage / relative_dir
        output_dir.mkdir(parents=True, exist_ok=True)
        for directory in directories:
            child = current_path / directory
            if child.is_symlink():
                raise MaterializeError(f"Symlinked source directory member is unsupported: {child}")
        for filename in files:
            source_file = current_path / filename
            if source_file.is_symlink() or not source_file.is_file():
                raise MaterializeError(f"Non-regular source member is unsupported: {source_file}")
            destination_file = output_dir / filename
            try:
                os.link(source_file, destination_file)
                hardlinks += 1
            except OSError as exc:
                if exc.errno not in {errno.EXDEV, errno.EPERM, errno.EACCES, errno.EMLINK}:
                    raise
                shutil.copy2(source_file, destination_file)
                copies += 1
    return hardlinks, copies


def inspect_input(workspace: Path, name: str) -> tuple[Path, Path, str, TreeSummary]:
    workspace_resolved = workspace.resolve(strict=True)
    source = workspace / "output" / "input_generation" / name
    destination = workspace / "input" / name
    source_resolved = source.resolve(strict=True)
    if not _is_relative_to(source_resolved, workspace_resolved):
        raise MaterializeError(f"Generated input source escapes the workspace: {source}")
    if not destination.is_symlink():
        raise MaterializeError(f"Expected the current workspace input to be a symlink: {destination}")
    original_link = os.readlink(destination)
    destination_resolved = destination.resolve(strict=True)
    if destination_resolved != source_resolved:
        raise MaterializeError(
            f"Workspace input does not point to its generated source: {destination} -> {original_link}"
        )
    return source_resolved, destination, original_link, summarize_tree(source_resolved)


def prepare_input(workspace: Path, name: str) -> PreparedInput:
    source, destination, original_link, source_summary = inspect_input(workspace, name)
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{name}.materialize-", dir=destination.parent))
    try:
        hardlinks, copies = _clone_tree(source, stage)
        staged_summary = summarize_tree(stage)
        if staged_summary != source_summary:
            raise MaterializeError(f"Staged input digest does not match its source: {name}")
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return PreparedInput(
        name=name,
        source=source,
        destination=destination,
        original_link=original_link,
        stage=stage,
        summary=source_summary,
        hardlinks=hardlinks,
        copies=copies,
    )


def _write_json_atomic(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def apply(workspace: Path, backup_root: Path, names: list[str]) -> list[dict[str, object]]:
    if backup_root.exists() and any(backup_root.iterdir()):
        raise MaterializeError(f"Backup root must be absent or empty: {backup_root}")
    prepared: list[PreparedInput] = []
    switched: list[tuple[PreparedInput, Path]] = []
    try:
        prepared = [prepare_input(workspace, name) for name in names]
        backup_root.mkdir(parents=True, exist_ok=True)
        for item in prepared:
            backup = backup_root / f"{item.name}.original-symlink"
            os.replace(item.destination, backup)
            try:
                os.replace(item.stage, item.destination)
            except Exception:
                os.replace(backup, item.destination)
                raise
            switched.append((item, backup))

        receipts: list[dict[str, object]] = []
        for item, backup in switched:
            current = summarize_tree(item.destination)
            if current != item.summary:
                raise MaterializeError(f"Post-switch digest does not match the source: {item.name}")
            receipt: dict[str, object] = {
                "schema_version": 1,
                "applied_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                "name": item.name,
                "source": str(item.source),
                "destination": str(item.destination),
                "original_symlink": item.original_link,
                "backup_symlink": str(backup),
                "sha256": item.summary.sha256,
                "files": item.summary.files,
                "bytes": item.summary.bytes,
                "hardlinks": item.hardlinks,
                "copies": item.copies,
            }
            _write_json_atomic(backup_root / f"{item.name}.receipt.json", receipt)
            receipts.append(receipt)
        return receipts
    except Exception:
        for item, backup in reversed(switched):
            failed_materialized = backup_root / f"{item.name}.failed-materialized"
            if item.destination.exists() and not failed_materialized.exists():
                os.replace(item.destination, failed_materialized)
            if backup.is_symlink() and not item.destination.exists():
                os.replace(backup, item.destination)
        raise
    finally:
        for item in prepared:
            if item.stage.exists():
                shutil.rmtree(item.stage, ignore_errors=True)


def verify(backup_root: Path, names: list[str]) -> list[dict[str, object]]:
    verified: list[dict[str, object]] = []
    for name in names:
        receipt_path = backup_root / f"{name}.receipt.json"
        try:
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise MaterializeError(f"Failed to read materialization receipt {receipt_path}: {exc}") from exc
        destination = Path(str(receipt.get("destination", "")))
        backup = Path(str(receipt.get("backup_symlink", "")))
        if destination.is_symlink() or not destination.is_dir():
            raise MaterializeError(f"Materialized destination is not a real directory: {destination}")
        if not backup.is_symlink() or os.readlink(backup) != receipt.get("original_symlink"):
            raise MaterializeError(f"Original recovery symlink is missing or changed: {backup}")
        summary = summarize_tree(destination)
        expected = TreeSummary(
            str(receipt.get("sha256", "")),
            int(receipt.get("files", -1)),
            int(receipt.get("bytes", -1)),
        )
        if summary != expected:
            raise MaterializeError(f"Materialized destination no longer matches its receipt: {destination}")
        verified.append(receipt)
    return verified


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("plan", "apply", "verify"))
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--backup-root", required=True, type=Path)
    parser.add_argument("--name", action="append", choices=SUPPORTED_INPUTS, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    names = list(dict.fromkeys(args.name))
    try:
        if args.command == "plan":
            rows = []
            for name in names:
                source, destination, original_link, summary = inspect_input(args.workspace, name)
                rows.append(
                    {
                        "name": name,
                        "source": str(source),
                        "destination": str(destination),
                        "original_symlink": original_link,
                        "sha256": summary.sha256,
                        "files": summary.files,
                        "bytes": summary.bytes,
                    }
                )
        elif args.command == "apply":
            rows = apply(args.workspace, args.backup_root, names)
        else:
            rows = verify(args.backup_root, names)
    except (OSError, ValueError, MaterializeError) as exc:
        print(json.dumps({"status": "error", "error": str(exc)}, sort_keys=True))
        return 1
    print(json.dumps({"status": "ok", "command": args.command, "inputs": rows}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

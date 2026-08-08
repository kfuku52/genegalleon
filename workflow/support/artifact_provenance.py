#!/usr/bin/env python3
"""Content-based provenance contracts for GeneGalleon workflow artifacts.

Only declared inputs, outputs, and output-affecting parameters participate in
freshness decisions. Producer, tool, container, and GeneGalleon versions are
retained as diagnostics and intentionally do not invalidate artifacts.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import datetime as dt
import hashlib
import io
import json
import os
import re
import sys
import tempfile
import zipfile
from pathlib import Path, PurePosixPath
from typing import BinaryIO, Iterable

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from gene_family_output_store import (
    GeneFamilyOutputStore,
    query_id_from_name,
    query_id_matchers,
)

SCHEMA_VERSION = 1
CURRENT = 1
NEEDS_RUN = 0
ERROR = 2
CHUNK_SIZE = 1024 * 1024
MANIFEST_SUBDIR = "artifact_provenance"
DEFAULT_REQUIRED_STEP_SUBDIRS = {
    "iqtree_anc": "iqtree_anc",
    "csubst": "csubst_b",
    "csubst_scan": "csubst_scan",
    "summary_statistics": "stat_branch",
    "tree_plot": "tree_plot",
}


class ProvenanceError(RuntimeError):
    """Raised for malformed or unsafe artifact provenance data."""


def parse_key_value(raw: str, option: str) -> tuple[str, str]:
    if "=" not in raw:
        raise ProvenanceError(f"{option} must use KEY=VALUE syntax: {raw!r}")
    key, value = raw.split("=", 1)
    key = key.strip()
    if not key:
        raise ProvenanceError(f"{option} key must not be empty: {raw!r}")
    return key, value


def parse_unique_pairs(values: Iterable[str], option: str) -> list[tuple[str, str]]:
    pairs = [parse_key_value(value, option) for value in values]
    keys = [key for key, _ in pairs]
    duplicates = sorted({key for key in keys if keys.count(key) > 1})
    if duplicates:
        raise ProvenanceError(f"Duplicate {option} key(s): {', '.join(duplicates)}")
    return pairs


def sha256_stream(handle: BinaryIO) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    while chunk := handle.read(CHUNK_SIZE):
        digest.update(chunk)
        size += len(chunk)
    return digest.hexdigest(), size


def sha256_path(path: Path) -> tuple[str, int, str]:
    if path.is_symlink():
        raise ProvenanceError(f"Symlinked provenance inputs and outputs are unsupported: {path}")
    if path.is_file():
        with path.open("rb") as handle:
            digest, size = sha256_stream(handle)
        return digest, size, "file"
    if not path.is_dir():
        raise FileNotFoundError(path)

    digest = hashlib.sha256()
    total_size = 0
    for child in sorted(path.rglob("*")):
        if child.is_symlink():
            raise ProvenanceError(f"Symlinked directory members are unsupported: {child}")
        if not child.is_file():
            continue
        relative = child.relative_to(path).as_posix().encode("utf-8")
        child_digest, child_size, _ = sha256_path(child)
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(bytes.fromhex(child_digest))
        digest.update(child_size.to_bytes(8, "big"))
        total_size += child_size
    return digest.hexdigest(), total_size, "directory"


def is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def path_reference(path: Path, logical_root: Path, workspace_root: Path) -> dict[str, str]:
    absolute = path.absolute()
    logical_absolute = logical_root.absolute()
    workspace_absolute = workspace_root.absolute()
    if is_relative_to(absolute, logical_absolute):
        relative = absolute.relative_to(logical_absolute)
        if len(relative.parts) != 2:
            raise ProvenanceError(f"Gene-family logical artifacts must use SUBDIR/FILENAME paths: {path}")
        return {"scope": "logical", "path": relative.as_posix()}
    if is_relative_to(absolute, workspace_absolute):
        return {
            "scope": "workspace",
            "path": absolute.relative_to(workspace_absolute).as_posix(),
        }
    return {"scope": "absolute", "path": str(absolute)}


def resolve_reference(reference: dict[str, str], logical_root: Path, workspace_root: Path) -> Path:
    scope = reference.get("scope")
    raw_path = reference.get("path", "")
    if scope == "logical":
        pure_path = PurePosixPath(raw_path)
        if len(pure_path.parts) != 2 or pure_path.is_absolute() or ".." in pure_path.parts:
            raise ProvenanceError(f"Unsafe logical provenance path: {raw_path!r}")
        return logical_root.joinpath(*pure_path.parts)
    if scope == "workspace":
        pure_path = PurePosixPath(raw_path)
        if pure_path.is_absolute() or ".." in pure_path.parts:
            raise ProvenanceError(f"Unsafe workspace provenance path: {raw_path!r}")
        return workspace_root.joinpath(*pure_path.parts)
    if scope == "absolute":
        path = Path(raw_path)
        if not path.is_absolute():
            raise ProvenanceError(f"Expected an absolute provenance path: {raw_path!r}")
        return path
    raise ProvenanceError(f"Unsupported provenance path scope: {scope!r}")


def describe_paths(pairs: list[tuple[str, str]], logical_root: Path, workspace_root: Path) -> list[dict[str, object]]:
    entries: list[dict[str, object]] = []
    for label, raw_path in pairs:
        path = Path(raw_path)
        digest, size, artifact_type = sha256_path(path)
        entry: dict[str, object] = {
            "label": label,
            **path_reference(path, logical_root, workspace_root),
            "artifact_type": artifact_type,
            "sha256": digest,
            "size_bytes": size,
        }
        entries.append(entry)
    return entries


def normalized_parameters(values: Iterable[str]) -> dict[str, str]:
    return dict(sorted(parse_unique_pairs(values, "--parameter")))


def normalized_diagnostics(values: Iterable[str]) -> dict[str, str]:
    return dict(sorted(parse_unique_pairs(values, "--diagnostic")))


def build_contract(args: argparse.Namespace, include_diagnostics: bool) -> dict[str, object]:
    logical_root = args.logical_root.absolute()
    workspace_root = args.workspace_root.absolute()
    contract: dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "step": args.step,
        "family_id": args.family_id,
        "inputs": describe_paths(parse_unique_pairs(args.input, "--input"), logical_root, workspace_root),
        "outputs": describe_paths(parse_unique_pairs(args.output, "--output"), logical_root, workspace_root),
        "parameters": normalized_parameters(args.parameter),
    }
    if include_diagnostics:
        contract["created_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        contract["diagnostics"] = normalized_diagnostics(args.diagnostic)
    return contract


def contract_comparison_payload(contract: dict[str, object]) -> dict[str, object]:
    return {
        "schema_version": contract.get("schema_version"),
        "step": contract.get("step"),
        "family_id": contract.get("family_id"),
        "inputs": contract.get("inputs"),
        "outputs": contract.get("outputs"),
        "parameters": contract.get("parameters"),
    }


def load_manifest(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ProvenanceError(f"Failed to read provenance manifest {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ProvenanceError(f"Provenance manifest must contain a JSON object: {path}")
    return payload


def write_manifest_atomic(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    encoded = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")
    try:
        with temporary.open("wb") as handle:
            handle.write(encoded)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def needs_run(args: argparse.Namespace) -> int:
    if not args.manifest.is_file():
        output_pairs = parse_unique_pairs(args.output, "--output")
        if not output_pairs:
            print(f"Artifact will be regenerated because no outputs were declared: {args.manifest}")
            return NEEDS_RUN
        try:
            for _label, raw_path in output_pairs:
                sha256_path(Path(raw_path))
        except FileNotFoundError as exc:
            print(f"Artifact will be regenerated because a declared output is missing: {exc.filename}")
            return NEEDS_RUN

        try:
            adopted = build_contract(args, include_diagnostics=False)
        except FileNotFoundError as exc:
            print(
                "Reusing legacy artifact without a provenance manifest; the manifest could not "
                f"be backfilled because a declared input is missing: {exc.filename}"
            )
            return CURRENT
        adopted["created_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        adopted["diagnostics"] = {
            "provenance_state": "adopted_legacy_output_without_rebuild",
        }
        write_manifest_atomic(args.manifest, adopted)
        print(f"Reusing legacy artifact and backfilled its provenance manifest: {args.manifest}")
        return CURRENT
    try:
        recorded = load_manifest(args.manifest)
        current = build_contract(args, include_diagnostics=False)
    except FileNotFoundError as exc:
        print(f"Artifact will be regenerated because a declared path is missing: {exc.filename}")
        return NEEDS_RUN
    if contract_comparison_payload(recorded) != contract_comparison_payload(current):
        print(
            "Artifact will be regenerated because an input, output, or output-affecting "
            f"parameter changed: {args.manifest}"
        )
        return NEEDS_RUN
    print(f"Artifact provenance is current: {args.manifest}")
    return CURRENT


def record(args: argparse.Namespace) -> int:
    payload = build_contract(args, include_diagnostics=True)
    write_manifest_atomic(args.manifest, payload)
    print(f"Recorded artifact provenance: {args.manifest}")
    return 0


def read_manifest_from_store(store: GeneFamilyOutputStore, name: str) -> dict[str, object]:
    try:
        with store.open_binary(MANIFEST_SUBDIR, name) as handle:
            payload = json.load(io.TextIOWrapper(handle, encoding="utf-8"))
    except (OSError, json.JSONDecodeError, UnicodeError, zipfile.BadZipFile) as exc:
        raise ProvenanceError(f"Failed to read logical provenance manifest {name}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ProvenanceError(f"Logical provenance manifest must contain a JSON object: {name}")
    return payload


@contextlib.contextmanager
def open_reference(
    reference: dict[str, object],
    store: GeneFamilyOutputStore,
    logical_root: Path,
    workspace_root: Path,
):
    scope = reference.get("scope")
    raw_path = str(reference.get("path", ""))
    if scope == "logical":
        pure_path = PurePosixPath(raw_path)
        if len(pure_path.parts) != 2 or pure_path.is_absolute() or ".." in pure_path.parts:
            raise ProvenanceError(f"Unsafe logical provenance path: {raw_path!r}")
        with store.open_binary(pure_path.parts[0], pure_path.parts[1]) as handle:
            yield handle
        return
    path = resolve_reference({"scope": str(scope), "path": raw_path}, logical_root, workspace_root)
    if not path.is_file() or path.is_symlink():
        raise FileNotFoundError(path)
    with path.open("rb") as handle:
        yield handle


def audit_manifest(
    payload: dict[str, object],
    store: GeneFamilyOutputStore,
    logical_root: Path,
    workspace_root: Path,
) -> tuple[str, str]:
    if payload.get("schema_version") != SCHEMA_VERSION:
        return "invalid_manifest", f"unsupported schema_version={payload.get('schema_version')!r}"
    for collection_name in ("inputs", "outputs"):
        entries = payload.get(collection_name)
        if not isinstance(entries, list):
            return "invalid_manifest", f"{collection_name} must be a list"
        for entry in entries:
            if not isinstance(entry, dict):
                return "invalid_manifest", f"{collection_name} contains a non-object entry"
            label = str(entry.get("label", ""))
            try:
                with open_reference(entry, store, logical_root, workspace_root) as handle:
                    digest, size = sha256_stream(handle)
            except FileNotFoundError:
                return f"missing_{collection_name[:-1]}", label
            except Exception as exc:
                return "invalid_manifest", f"{label}: {exc}"
            if digest != entry.get("sha256") or size != entry.get("size_bytes"):
                return f"changed_{collection_name[:-1]}", label
    return "current", ""


def infer_orthogroup_id(name: str) -> str | None:
    match = re.match(r"^(OG[0-9]+)(?:[_.]|$)", name)
    return match.group(1) if match else None


def parse_required_step_subdirs(values: Iterable[str]) -> dict[str, str]:
    if not values:
        return dict(DEFAULT_REQUIRED_STEP_SUBDIRS)
    return dict(parse_unique_pairs(values, "--require-step-for-subdir"))


def artifact_family_id(artifact, mode: str, query_matcher_list=None) -> str | None:
    if artifact.family_id:
        return str(artifact.family_id)
    if mode == "orthogroup":
        return infer_orthogroup_id(artifact.name)
    if mode == "query2family" and query_matcher_list:
        return query_id_from_name(artifact.name, query_matcher_list)
    return None


def safe_extract_zip_bytes(data: bytes, destination: Path) -> Path:
    with zipfile.ZipFile(io.BytesIO(data), "r") as archive:
        for member in archive.infolist():
            pure_name = PurePosixPath(member.filename)
            if pure_name.is_absolute() or ".." in pure_name.parts:
                raise ProvenanceError(f"Unsafe ZIP member in IQ-TREE archive: {member.filename}")
        archive.extractall(destination)
    candidates = sorted(path.parent for path in destination.rglob("csubst.nwk"))
    if len(candidates) != 1:
        raise ProvenanceError(f"Expected one extracted csubst.nwk, found {len(candidates)} in IQ-TREE archive")
    return candidates[0]


def branch_identity_rows(
    store: GeneFamilyOutputStore,
    mode: str,
    query_matcher_list=None,
) -> list[dict[str, str]]:
    iqtree_by_family = {
        family_id: artifact
        for artifact in store.artifacts("iqtree_anc")
        if (family_id := artifact_family_id(artifact, mode, query_matcher_list)) is not None
    }
    stat_by_family = {
        family_id: artifact
        for artifact in store.artifacts("stat_branch")
        if (family_id := artifact_family_id(artifact, mode, query_matcher_list)) is not None
    }
    rows: list[dict[str, str]] = []
    try:
        from csubst_site_wrapper import validate_all_csubst_stat_branch_identity
    except ImportError as exc:
        return [
            {
                "family_id": "-",
                "step": "branch_identity",
                "status": "audit_error",
                "reason": f"failed to import CSUBST branch validator: {exc}",
                "manifest": "",
            }
        ]
    for family_id in sorted(set(iqtree_by_family).intersection(stat_by_family)):
        iqtree_artifact = iqtree_by_family[family_id]
        stat_artifact = stat_by_family[family_id]
        try:
            with store.open_binary(iqtree_artifact.subdir, iqtree_artifact.name) as handle:
                iqtree_data = handle.read()
            with store.open_binary(stat_artifact.subdir, stat_artifact.name) as handle:
                stat_data = handle.read()
            with tempfile.TemporaryDirectory(prefix="gg-branch-audit-") as temporary:
                temporary_path = Path(temporary)
                iqtree_dir = safe_extract_zip_bytes(iqtree_data, temporary_path / "iqtree")
                stat_path = temporary_path / "stat_branch.tsv"
                stat_path.write_bytes(stat_data)
                validate_all_csubst_stat_branch_identity(str(stat_path), str(iqtree_dir))
            status, reason = "current", ""
        except Exception as exc:
            status, reason = "semantic_mismatch", str(exc)
        rows.append(
            {
                "family_id": family_id,
                "step": "branch_identity",
                "status": status,
                "reason": reason,
                "manifest": "",
            }
        )
    return rows


def audit(args: argparse.Namespace) -> int:
    logical_root = args.logical_root.absolute()
    workspace_root = args.workspace_root.absolute()
    store = GeneFamilyOutputStore(logical_root)
    query_matcher_list = None
    if args.mode == "query2family" and args.query_dir is not None:
        if not args.query_dir.is_dir():
            raise ProvenanceError(f"Query directory not found: {args.query_dir}")
        query_matcher_list = query_id_matchers(
            sorted(path.name for path in args.query_dir.iterdir() if path.is_file() and not path.name.startswith("."))
        )
    rows: list[dict[str, str]] = []
    manifested_steps: set[tuple[str, str]] = set()
    if MANIFEST_SUBDIR in store.logical_subdirs():
        for artifact in store.artifacts(MANIFEST_SUBDIR):
            try:
                payload = read_manifest_from_store(store, artifact.name)
                family_id = str(payload.get("family_id", ""))
                step = str(payload.get("step", ""))
                if not family_id or not step:
                    raise ProvenanceError("family_id and step are required")
                status, reason = audit_manifest(payload, store, logical_root, workspace_root)
                manifested_steps.add((family_id, step))
            except Exception as exc:
                family_id = artifact_family_id(artifact, args.mode, query_matcher_list) or "-"
                step = "unknown"
                status, reason = "invalid_manifest", str(exc)
            rows.append(
                {
                    "family_id": family_id,
                    "step": step,
                    "status": status,
                    "reason": reason,
                    "manifest": artifact.logical_path,
                }
            )

    for step, subdir in parse_required_step_subdirs(args.require_step_for_subdir).items():
        if subdir not in store.logical_subdirs():
            continue
        for artifact in store.artifacts(subdir):
            family_id = artifact_family_id(artifact, args.mode, query_matcher_list)
            if family_id is None or (family_id, step) in manifested_steps:
                continue
            rows.append(
                {
                    "family_id": family_id or "-",
                    "step": step,
                    "status": "legacy_untracked",
                    "reason": artifact.logical_path,
                    "manifest": "",
                }
            )

    if args.check_csubst_branches:
        rows.extend(branch_identity_rows(store, args.mode, query_matcher_list))

    rows.sort(key=lambda row: (row["family_id"], row["step"], row["status"]))
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output_tsv.with_name(f".{args.output_tsv.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=("family_id", "step", "status", "reason", "manifest"),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, args.output_tsv)
    finally:
        temporary.unlink(missing_ok=True)

    nonfailure_statuses = {"current", "legacy_untracked"}
    failures = [row for row in rows if row["status"] not in nonfailure_statuses]
    print(f"Artifact provenance audit: checked={len(rows)}, failures={len(failures)}, report={args.output_tsv}")
    for row in failures[:20]:
        print(
            f"{row['family_id']}\t{row['step']}\t{row['status']}\t{row['reason']}",
            file=sys.stderr,
        )
    return 1 if failures else 0


def add_contract_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--step", required=True)
    parser.add_argument("--family-id", required=True)
    parser.add_argument("--logical-root", required=True, type=Path)
    parser.add_argument("--workspace-root", required=True, type=Path)
    parser.add_argument("--input", action="append", default=[], metavar="LABEL=PATH")
    parser.add_argument("--output", action="append", default=[], metavar="LABEL=PATH")
    parser.add_argument("--parameter", action="append", default=[], metavar="KEY=VALUE")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    needs_parser = subparsers.add_parser(
        "needs-run", help="Return 0 when a declared artifact must be regenerated and 1 when current"
    )
    add_contract_arguments(needs_parser)

    record_parser = subparsers.add_parser("record", help="Write a completed artifact contract")
    add_contract_arguments(record_parser)
    record_parser.add_argument("--diagnostic", action="append", default=[], metavar="KEY=VALUE")

    audit_parser = subparsers.add_parser("audit", help="Audit raw and ZIP-backed provenance manifests")
    audit_parser.add_argument("--logical-root", required=True, type=Path)
    audit_parser.add_argument("--workspace-root", required=True, type=Path)
    audit_parser.add_argument("--output-tsv", required=True, type=Path)
    audit_parser.add_argument("--mode", choices=("orthogroup", "query2family"), required=True)
    audit_parser.add_argument("--query-dir", type=Path)
    audit_parser.add_argument(
        "--require-step-for-subdir",
        action="append",
        default=[],
        metavar="STEP=SUBDIR",
        help="Require a manifest for every family represented in the logical subdirectory",
    )
    audit_parser.add_argument("--check-csubst-branches", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "needs-run":
            return needs_run(args)
        if args.command == "record":
            return record(args)
        if args.command == "audit":
            return audit(args)
        parser.error(f"Unsupported command: {args.command}")
    except (OSError, ProvenanceError, ValueError) as exc:
        print(f"Artifact provenance error: {exc}", file=sys.stderr)
        return ERROR
    return ERROR


if __name__ == "__main__":
    raise SystemExit(main())

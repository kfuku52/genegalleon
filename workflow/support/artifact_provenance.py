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
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path, PurePosixPath
from typing import BinaryIO, Iterable

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from content_digest_cache import cache as digest_cache
from content_digest_cache import cached_sha256_file
from content_digest_cache import configure as configure_digest_cache
from gene_family_output_store import (
    GeneFamilyOutputStore,
    query_id_from_name,
    query_id_matchers,
)
from safe_zip_extract import extract_expected_prefix, validated_members

SCHEMA_VERSION = 1
CURRENT = 1
NEEDS_RUN = 0
ERROR = 2
STALE_STOP = 3
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


def parse_path_pairs(values: Iterable[str], option: str) -> list[tuple[str, str]]:
    pairs = parse_unique_pairs(values, option)
    empty = [label for label, path in pairs if not path.strip()]
    if empty:
        raise ProvenanceError(f"{option} path must not be empty for key(s): {', '.join(empty)}")
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
        digest, size = cached_sha256_file(path)
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


def describe_optional_paths(
    pairs: list[tuple[str, str]], logical_root: Path, workspace_root: Path
) -> list[dict[str, object]]:
    """Describe managed outputs whose absence can be a valid result."""

    entries: list[dict[str, object]] = []
    for label, raw_path in pairs:
        path = Path(raw_path)
        reference = path_reference(path, logical_root, workspace_root)
        try:
            digest, size, artifact_type = sha256_path(path)
        except FileNotFoundError:
            entries.append(
                {
                    "label": label,
                    **reference,
                    "state": "absent",
                }
            )
            continue
        entries.append(
            {
                "label": label,
                **reference,
                "state": "present",
                "artifact_type": artifact_type,
                "sha256": digest,
                "size_bytes": size,
            }
        )
    return entries


def gene_family_store_digest(root: Path) -> tuple[str, int, int]:
    store = GeneFamilyOutputStore(root)
    digest = hashlib.sha256()
    total_size = 0
    member_count = 0
    for subdir in store.logical_subdirs():
        if subdir == MANIFEST_SUBDIR:
            continue
        for artifact in store.artifacts(subdir):
            artifact_digest = artifact.sha256
            artifact_size = artifact.size
            if not artifact_digest or artifact_size is None:
                with store.open_binary(artifact.subdir, artifact.name) as handle:
                    artifact_digest, artifact_size = sha256_stream(handle)
            logical_path = artifact.logical_path.encode("utf-8")
            digest.update(len(logical_path).to_bytes(8, "big"))
            digest.update(logical_path)
            digest.update(bytes.fromhex(str(artifact_digest)))
            digest.update(int(artifact_size).to_bytes(8, "big"))
            total_size += int(artifact_size)
            member_count += 1
    return digest.hexdigest(), total_size, member_count


def gene_family_subdir_digest(root: Path, subdir: str) -> tuple[str, int, int]:
    """Hash one logical gene-family subdirectory independent of raw/ZIP layout."""

    pure_subdir = PurePosixPath(subdir)
    if (
        pure_subdir.is_absolute()
        or len(pure_subdir.parts) != 1
        or pure_subdir.parts[0] in {"", ".", ".."}
    ):
        raise ProvenanceError(f"Unsafe gene-family logical subdirectory: {subdir!r}")
    if not root.is_dir() or root.is_symlink():
        raise FileNotFoundError(root)
    store = GeneFamilyOutputStore(root)
    digest = hashlib.sha256()
    total_size = 0
    member_count = 0
    for artifact in store.artifacts(subdir):
        artifact_digest = artifact.sha256
        artifact_size = artifact.size
        if not artifact_digest or artifact_size is None:
            with store.open_binary(artifact.subdir, artifact.name) as handle:
                artifact_digest, artifact_size = sha256_stream(handle)
        logical_path = artifact.logical_path.encode("utf-8")
        digest.update(len(logical_path).to_bytes(8, "big"))
        digest.update(logical_path)
        digest.update(bytes.fromhex(str(artifact_digest)))
        digest.update(int(artifact_size).to_bytes(8, "big"))
        total_size += int(artifact_size)
        member_count += 1
    return digest.hexdigest(), total_size, member_count


def gene_family_artifact_digest(root: Path, subdir: str, name: str) -> tuple[str, int]:
    """Hash one logical artifact independent of raw/ZIP storage layout."""

    pure_subdir = PurePosixPath(subdir)
    pure_name = PurePosixPath(name)
    if (
        pure_subdir.is_absolute()
        or len(pure_subdir.parts) != 1
        or pure_subdir.parts[0] in {"", ".", ".."}
        or pure_name.is_absolute()
        or len(pure_name.parts) != 1
        or pure_name.parts[0] in {"", ".", ".."}
    ):
        raise ProvenanceError(
            f"Unsafe gene-family logical artifact: subdir={subdir!r}, name={name!r}"
        )
    if not root.is_dir() or root.is_symlink():
        raise FileNotFoundError(root)
    store = GeneFamilyOutputStore(root)
    artifact = store.artifact(subdir, name)
    if artifact is None:
        raise FileNotFoundError(root / subdir / name)
    artifact_digest = artifact.sha256
    artifact_size = artifact.size
    if not artifact_digest or artifact_size is None:
        with store.open_binary(artifact.subdir, artifact.name) as handle:
            artifact_digest, artifact_size = sha256_stream(handle)
    return str(artifact_digest), int(artifact_size)


def parse_gene_family_subdir_pairs(
    values: Iterable[str], option: str
) -> list[tuple[str, str, str]]:
    pairs = parse_path_pairs(values, option)
    parsed: list[tuple[str, str, str]] = []
    for label, raw_value in pairs:
        if "::" not in raw_value:
            raise ProvenanceError(f"{option} must use LABEL=ROOT::SUBDIR syntax: {raw_value!r}")
        root, subdir = raw_value.rsplit("::", 1)
        if not root.strip() or not subdir.strip():
            raise ProvenanceError(f"{option} must use LABEL=ROOT::SUBDIR syntax: {raw_value!r}")
        parsed.append((label, root, subdir))
    return parsed


def parse_gene_family_artifact_pairs(
    values: Iterable[str], option: str
) -> list[tuple[str, str, str, str]]:
    pairs = parse_path_pairs(values, option)
    parsed: list[tuple[str, str, str, str]] = []
    for label, raw_value in pairs:
        parts = raw_value.rsplit("::", 2)
        if len(parts) != 3 or any(not part.strip() for part in parts):
            raise ProvenanceError(
                f"{option} must use LABEL=ROOT::SUBDIR::FILENAME syntax: {raw_value!r}"
            )
        root, subdir, name = parts
        parsed.append((label, root, subdir, name))
    return parsed


def describe_gene_family_stores(
    pairs: list[tuple[str, str]], logical_root: Path, workspace_root: Path
) -> list[dict[str, object]]:
    entries: list[dict[str, object]] = []
    for label, raw_path in pairs:
        path = Path(raw_path)
        if not path.is_dir() or path.is_symlink():
            raise FileNotFoundError(path)
        digest, size, member_count = gene_family_store_digest(path)
        entries.append(
            {
                "label": label,
                **path_reference(path, logical_root, workspace_root),
                "artifact_type": "gene_family_store",
                "sha256": digest,
                "size_bytes": size,
                "member_count": member_count,
            }
        )
    return entries


def describe_gene_family_subdirs(
    pairs: list[tuple[str, str, str]], logical_root: Path, workspace_root: Path
) -> list[dict[str, object]]:
    entries: list[dict[str, object]] = []
    for label, raw_root, subdir in pairs:
        root = Path(raw_root)
        digest, size, member_count = gene_family_subdir_digest(root, subdir)
        entries.append(
            {
                "label": label,
                **path_reference(root, logical_root, workspace_root),
                "artifact_type": "gene_family_subdir",
                "subdir": subdir,
                "sha256": digest,
                "size_bytes": size,
                "member_count": member_count,
            }
        )
    return entries


def describe_gene_family_artifacts(
    pairs: list[tuple[str, str, str, str]], logical_root: Path, workspace_root: Path
) -> list[dict[str, object]]:
    entries: list[dict[str, object]] = []
    for label, raw_root, subdir, name in pairs:
        root = Path(raw_root)
        digest, size = gene_family_artifact_digest(root, subdir, name)
        entries.append(
            {
                "label": label,
                **path_reference(root, logical_root, workspace_root),
                "artifact_type": "gene_family_artifact",
                "subdir": subdir,
                "name": name,
                "sha256": digest,
                "size_bytes": size,
            }
        )
    return entries


def raw_or_zip_directory_digest(path: Path) -> tuple[str, int, int]:
    archive_path = path.with_name(f"{path.name}.zip")
    raw_result: tuple[str, int, int] | None = None
    archive_result: tuple[str, int, int] | None = None
    if path.exists():
        digest, size, artifact_type = sha256_path(path)
        if artifact_type != "directory":
            raise ProvenanceError(f"Managed logical directory path is not a directory: {path}")
        member_count = sum(1 for child in path.rglob("*") if child.is_file())
        raw_result = (digest, size, member_count)
    if archive_path.exists():
        if archive_path.is_symlink() or not archive_path.is_file():
            raise ProvenanceError(f"Managed logical directory archive is not a regular file: {archive_path}")
        current_cache = digest_cache()
        cached_archive = (
            current_cache.get("logical_zip_directory", archive_path)
            if current_cache is not None
            else None
        )
        if cached_archive is not None and cached_archive[2] is not None:
            archive_result = (cached_archive[0], cached_archive[1], cached_archive[2])
        else:
            archive_signature = None
            if current_cache is not None:
                try:
                    archive_signature = current_cache.signature(archive_path)
                except OSError:
                    pass
            digest = hashlib.sha256()
            total_size = 0
            member_count = 0
            try:
                with zipfile.ZipFile(archive_path, "r") as archive:
                    infos, _directories = validated_members(
                        archive,
                        archive_path,
                        path.name,
                    )
                    prefix = PurePosixPath(path.name)
                    members: list[tuple[PurePosixPath, zipfile.ZipInfo]] = []
                    for info in infos:
                        member = PurePosixPath(info.filename)
                        if member.is_absolute() or ".." in member.parts:
                            raise ProvenanceError(f"Unsafe logical directory ZIP member: {info.filename!r}")
                        try:
                            relative = member.relative_to(prefix)
                        except ValueError as exc:
                            raise ProvenanceError(
                                f"Logical directory ZIP member is outside {path.name}/: {info.filename!r}"
                            ) from exc
                        if not relative.parts:
                            continue
                        members.append((relative, info))
                    for relative, info in sorted(members, key=lambda item: item[0].as_posix()):
                        relative_bytes = relative.as_posix().encode("utf-8")
                        with archive.open(info, "r") as handle:
                            member_digest, member_size = sha256_stream(handle)
                        digest.update(len(relative_bytes).to_bytes(8, "big"))
                        digest.update(relative_bytes)
                        digest.update(bytes.fromhex(member_digest))
                        digest.update(member_size.to_bytes(8, "big"))
                        total_size += member_size
                        member_count += 1
            except (OSError, RuntimeError, zipfile.BadZipFile) as exc:
                raise ProvenanceError(f"Failed to read logical directory ZIP {archive_path}: {exc}") from exc
            archive_result = (digest.hexdigest(), total_size, member_count)
            if current_cache is not None and archive_signature is not None:
                current_cache.put(
                    "logical_zip_directory",
                    archive_path,
                    archive_result[0],
                    archive_result[1],
                    archive_result[2],
                    expected_signature=archive_signature,
                )
    if raw_result is None and archive_result is None:
        raise FileNotFoundError(path)
    if raw_result is not None and archive_result is not None and raw_result != archive_result:
        raise ProvenanceError(
            f"Raw and ZIP forms of a managed logical directory differ: {path} and {archive_path}"
        )
    return raw_result if raw_result is not None else archive_result  # type: ignore[return-value]


def describe_logical_directories(
    pairs: list[tuple[str, str]], logical_root: Path, workspace_root: Path
) -> list[dict[str, object]]:
    entries: list[dict[str, object]] = []
    for label, raw_path in pairs:
        path = Path(raw_path)
        digest, size, member_count = raw_or_zip_directory_digest(path)
        entries.append(
            {
                "label": label,
                **path_reference(path, logical_root, workspace_root),
                "artifact_type": "logical_directory",
                "sha256": digest,
                "size_bytes": size,
                "member_count": member_count,
            }
        )
    return entries


def normalized_parameters(values: Iterable[str]) -> dict[str, str]:
    return dict(sorted(parse_unique_pairs(values, "--parameter")))


def normalized_diagnostics(values: Iterable[str]) -> dict[str, str]:
    return dict(sorted(parse_unique_pairs(values, "--diagnostic")))


def build_contract(args: argparse.Namespace, include_diagnostics: bool) -> dict[str, object]:
    logical_root = args.logical_root.absolute()
    workspace_root = args.workspace_root.absolute()
    input_pairs = parse_path_pairs(args.input, "--input")
    store_pairs = parse_path_pairs(args.input_gene_family_store, "--input-gene-family-store")
    store_subdir_pairs = parse_gene_family_subdir_pairs(
        args.input_gene_family_subdir, "--input-gene-family-subdir"
    )
    store_artifact_pairs = parse_gene_family_artifact_pairs(
        args.input_gene_family_artifact, "--input-gene-family-artifact"
    )
    logical_input_pairs = parse_path_pairs(args.input_logical_directory, "--input-logical-directory")
    output_pairs = parse_path_pairs(args.output, "--output")
    logical_output_pairs = parse_path_pairs(args.output_logical_directory, "--output-logical-directory")
    input_labels = [label for label, _path in input_pairs + store_pairs + logical_input_pairs]
    input_labels.extend(label for label, _root, _subdir in store_subdir_pairs)
    input_labels.extend(label for label, _root, _subdir, _name in store_artifact_pairs)
    duplicate_input_labels = sorted({label for label in input_labels if input_labels.count(label) > 1})
    if duplicate_input_labels:
        raise ProvenanceError(f"Duplicate input key(s): {', '.join(duplicate_input_labels)}")
    output_labels = [label for label, _path in output_pairs + logical_output_pairs]
    duplicate_output_labels = sorted({label for label in output_labels if output_labels.count(label) > 1})
    if duplicate_output_labels:
        raise ProvenanceError(f"Duplicate output key(s): {', '.join(duplicate_output_labels)}")
    contract: dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "step": args.step,
        "family_id": args.family_id,
        "inputs": describe_paths(input_pairs, logical_root, workspace_root)
        + describe_gene_family_stores(store_pairs, logical_root, workspace_root)
        + describe_gene_family_subdirs(store_subdir_pairs, logical_root, workspace_root)
        + describe_gene_family_artifacts(store_artifact_pairs, logical_root, workspace_root)
        + describe_logical_directories(logical_input_pairs, logical_root, workspace_root),
        "outputs": describe_paths(output_pairs, logical_root, workspace_root)
        + describe_logical_directories(logical_output_pairs, logical_root, workspace_root),
        "optional_outputs": describe_optional_paths(
            parse_path_pairs(args.optional_output, "--optional-output"), logical_root, workspace_root
        ),
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
        "optional_outputs": contract.get("optional_outputs", []),
        "parameters": contract.get("parameters"),
    }


def entries_by_label(contract: dict[str, object], collection: str) -> dict[str, object]:
    entries = contract.get(collection, [])
    if not isinstance(entries, list):
        return {}
    return {
        str(entry.get("label", "")): entry
        for entry in entries
        if isinstance(entry, dict) and entry.get("label")
    }


def describe_contract_difference(recorded: dict[str, object], current: dict[str, object]) -> str:
    for key in ("schema_version", "step", "family_id"):
        if recorded.get(key) != current.get(key):
            return f"{key} changed: recorded={recorded.get(key)!r}, current={current.get(key)!r}"

    recorded_parameters = recorded.get("parameters", {})
    current_parameters = current.get("parameters", {})
    if recorded_parameters != current_parameters:
        recorded_parameters = recorded_parameters if isinstance(recorded_parameters, dict) else {}
        current_parameters = current_parameters if isinstance(current_parameters, dict) else {}
        for key in sorted(set(recorded_parameters) | set(current_parameters)):
            if recorded_parameters.get(key) != current_parameters.get(key):
                return (
                    f"output-affecting parameter {key!r} changed: "
                    f"recorded={recorded_parameters.get(key)!r}, current={current_parameters.get(key)!r}"
                )
        return "output-affecting parameters changed"

    for collection, description in (
        ("inputs", "input"),
        ("outputs", "output"),
        ("optional_outputs", "managed optional output"),
    ):
        recorded_entries = entries_by_label(recorded, collection)
        current_entries = entries_by_label(current, collection)
        if recorded_entries.keys() != current_entries.keys():
            return f"declared {description} set changed"
        for label in sorted(recorded_entries):
            if recorded_entries[label] != current_entries[label]:
                return f"{description} {label!r} changed"
    return "artifact contract changed"


def print_stale_message(args: argparse.Namespace, reason: str, action: str) -> None:
    print("Stale artifact detected.", file=sys.stderr)
    print(f"Family: {args.family_id}", file=sys.stderr)
    print(f"Step: {args.step}", file=sys.stderr)
    print(f"Reason: {reason}", file=sys.stderr)
    print(f"Manifest: {args.manifest}", file=sys.stderr)
    print(f"Policy action: {action}", file=sys.stderr)


def stale_policy_result(
    args: argparse.Namespace, reason: str, *, reusable: bool = True
) -> int:
    policy = args.stale_policy
    if policy == "rebuild":
        print_stale_message(args, reason, "regenerate without confirmation")
        return NEEDS_RUN
    if policy == "reuse" and reusable:
        print_stale_message(args, reason, "reuse stale output without changing its manifest")
        return CURRENT

    print_stale_message(args, reason, "stop before modifying outputs")
    print("No artifact files were modified.", file=sys.stderr)
    print("Choose one of:", file=sys.stderr)
    print("  artifact_stale_policy=rebuild  # regenerate automatically", file=sys.stderr)
    if reusable:
        print("  artifact_stale_policy=reuse    # continue with the stale output", file=sys.stderr)
    return STALE_STOP


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
        output_pairs = parse_path_pairs(args.output, "--output")
        logical_output_pairs = parse_path_pairs(
            args.output_logical_directory, "--output-logical-directory"
        )
        if not output_pairs and not logical_output_pairs:
            optional_pairs = parse_path_pairs(args.optional_output, "--optional-output")
            present_optional = [
                label for label, raw_path in optional_pairs if Path(raw_path).is_file() or Path(raw_path).is_dir()
            ]
            if not present_optional:
                print(
                    "Artifact will be generated because no provenance manifest or declared optional outputs exist: "
                    f"{args.manifest}"
                )
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
                "provenance_state": "adopted_legacy_optional_output_without_rebuild",
            }
            write_manifest_atomic(args.manifest, adopted)
            print(f"Reusing legacy optional artifact and backfilled its provenance manifest: {args.manifest}")
            return CURRENT
        present_outputs = []
        missing_outputs = []
        for label, raw_path in output_pairs:
            try:
                sha256_path(Path(raw_path))
                present_outputs.append(label)
            except FileNotFoundError:
                missing_outputs.append(label)
        for label, raw_path in logical_output_pairs:
            try:
                raw_or_zip_directory_digest(Path(raw_path))
                present_outputs.append(label)
            except FileNotFoundError:
                missing_outputs.append(label)
        if missing_outputs and not present_outputs:
            print(f"Artifact will be generated because no declared outputs exist: {args.manifest}")
            return NEEDS_RUN
        if missing_outputs:
            return stale_policy_result(
                args,
                "legacy output set is incomplete; missing required output(s): "
                + ", ".join(missing_outputs),
                reusable=False,
            )

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
    except FileNotFoundError as exc:
        print(f"Artifact provenance error: declared path is missing: {exc.filename}", file=sys.stderr)
        return ERROR
    try:
        current = build_contract(args, include_diagnostics=False)
    except FileNotFoundError as exc:
        missing_path = str(exc.filename or exc)
        input_paths = {raw_path for _label, raw_path in parse_path_pairs(args.input, "--input")}
        input_paths.update(
            raw_path
            for _label, raw_path in parse_path_pairs(
                args.input_gene_family_store, "--input-gene-family-store"
            )
        )
        input_paths.update(
            raw_root
            for _label, raw_root, _subdir in parse_gene_family_subdir_pairs(
                args.input_gene_family_subdir, "--input-gene-family-subdir"
            )
        )
        for _label, raw_root, subdir, name in parse_gene_family_artifact_pairs(
            args.input_gene_family_artifact, "--input-gene-family-artifact"
        ):
            input_paths.add(raw_root)
            input_paths.add(str(Path(raw_root) / subdir / name))
        input_paths.update(
            raw_path
            for _label, raw_path in parse_path_pairs(
                args.input_logical_directory, "--input-logical-directory"
            )
        )
        if missing_path in input_paths:
            print(f"Artifact provenance error: declared input is missing: {missing_path}", file=sys.stderr)
            return ERROR
        return stale_policy_result(
            args,
            f"required output is missing: {missing_path}",
            reusable=False,
        )

    # Optional outputs were introduced after the first manifest schema. Adopt
    # their current presence/absence once so existing tracked artifacts remain
    # compatible with the legacy-output policy.
    if "optional_outputs" not in recorded and current.get("optional_outputs"):
        comparison_without_optional = contract_comparison_payload(current)
        comparison_without_optional["optional_outputs"] = []
        if contract_comparison_payload(recorded) == comparison_without_optional:
            recorded["optional_outputs"] = current["optional_outputs"]
            diagnostics = recorded.setdefault("diagnostics", {})
            if isinstance(diagnostics, dict):
                diagnostics["optional_output_provenance_state"] = "adopted_legacy_optional_outputs"
            write_manifest_atomic(args.manifest, recorded)
            print(f"Adopted previously untracked optional outputs: {args.manifest}")
            return CURRENT
    if contract_comparison_payload(recorded) != contract_comparison_payload(current):
        return stale_policy_result(
            args,
            describe_contract_difference(recorded, current),
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
            if collection_name == "inputs" and entry.get("artifact_type") in {
                "gene_family_store",
                "gene_family_subdir",
                "gene_family_artifact",
            }:
                try:
                    store_path = resolve_reference(entry, logical_root, workspace_root)
                    if entry.get("artifact_type") == "gene_family_artifact":
                        digest, size = gene_family_artifact_digest(
                            store_path,
                            str(entry.get("subdir", "")),
                            str(entry.get("name", "")),
                        )
                        member_count = None
                    elif entry.get("artifact_type") == "gene_family_subdir":
                        digest, size, member_count = gene_family_subdir_digest(
                            store_path, str(entry.get("subdir", ""))
                        )
                    else:
                        digest, size, member_count = gene_family_store_digest(store_path)
                except FileNotFoundError:
                    return "missing_input", label
                except Exception as exc:
                    return "invalid_manifest", f"{label}: {exc}"
                if (
                    digest != entry.get("sha256")
                    or size != entry.get("size_bytes")
                    or (
                        entry.get("artifact_type") != "gene_family_artifact"
                        and member_count != entry.get("member_count")
                    )
                ):
                    return "changed_input", label
                continue
            try:
                with open_reference(entry, store, logical_root, workspace_root) as handle:
                    digest, size = sha256_stream(handle)
            except FileNotFoundError:
                return f"missing_{collection_name[:-1]}", label
            except Exception as exc:
                return "invalid_manifest", f"{label}: {exc}"
            if digest != entry.get("sha256") or size != entry.get("size_bytes"):
                return f"changed_{collection_name[:-1]}", label

    optional_entries = payload.get("optional_outputs", [])
    if not isinstance(optional_entries, list):
        return "invalid_manifest", "optional_outputs must be a list"
    for entry in optional_entries:
        if not isinstance(entry, dict):
            return "invalid_manifest", "optional_outputs contains a non-object entry"
        label = str(entry.get("label", ""))
        state = entry.get("state")
        if state not in {"present", "absent"}:
            return "invalid_manifest", f"optional output {label!r} has invalid state={state!r}"
        try:
            with open_reference(entry, store, logical_root, workspace_root) as handle:
                digest, size = sha256_stream(handle)
        except FileNotFoundError:
            if state == "absent":
                continue
            return "missing_optional_output", label
        except Exception as exc:
            return "invalid_manifest", f"{label}: {exc}"
        if state == "absent":
            return "unexpected_optional_output", label
        if digest != entry.get("sha256") or size != entry.get("size_bytes"):
            return "changed_optional_output", label
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


def safe_extract_iqtree_zip(
    archive_path: Path,
    destination: Path,
    family_id: str,
) -> Path:
    expected_prefix = f"{family_id}.iqtree.anc"
    extracted = extract_expected_prefix(
        archive_path,
        destination,
        expected_prefix,
    )
    csubst_tree = extracted / "csubst.nwk"
    if not csubst_tree.is_file():
        raise ProvenanceError(f"Expected {csubst_tree} in IQ-TREE archive")
    return extracted


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
            with tempfile.TemporaryDirectory(prefix="gg-branch-audit-") as temporary:
                temporary_path = Path(temporary)
                iqtree_path = temporary_path / "iqtree.zip"
                with store.open_binary(
                    iqtree_artifact.subdir,
                    iqtree_artifact.name,
                ) as source, iqtree_path.open("wb") as output:
                    shutil.copyfileobj(source, output, length=CHUNK_SIZE)
                    output.flush()
                    os.fsync(output.fileno())
                iqtree_dir = safe_extract_iqtree_zip(
                    iqtree_path,
                    temporary_path / "iqtree",
                    family_id,
                )
                stat_path = temporary_path / "stat_branch.tsv"
                with store.open_binary(
                    stat_artifact.subdir,
                    stat_artifact.name,
                ) as source, stat_path.open("wb") as output:
                    shutil.copyfileobj(source, output, length=CHUNK_SIZE)
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
    if args.stale_policy == "reuse":
        nonfailure_statuses.update(
            {
                "changed_input",
                "changed_output",
                "changed_optional_output",
                "missing_optional_output",
                "unexpected_optional_output",
            }
        )
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
    parser.add_argument(
        "--input-gene-family-store",
        action="append",
        default=[],
        metavar="LABEL=ROOT",
        help="Hash logical gene-family artifacts while ignoring ZIP/raw storage layout and provenance manifests",
    )
    parser.add_argument(
        "--input-gene-family-subdir",
        action="append",
        default=[],
        metavar="LABEL=ROOT::SUBDIR",
        help="Hash one logical gene-family subdirectory while ignoring ZIP/raw storage layout",
    )
    parser.add_argument(
        "--input-gene-family-artifact",
        action="append",
        default=[],
        metavar="LABEL=ROOT::SUBDIR::FILENAME",
        help="Hash one logical gene-family artifact while ignoring ZIP/raw storage layout",
    )
    parser.add_argument(
        "--input-logical-directory",
        action="append",
        default=[],
        metavar="LABEL=PATH",
        help="Hash a managed directory identically in raw PATH and sibling PATH.zip forms",
    )
    parser.add_argument("--output", action="append", default=[], metavar="LABEL=PATH")
    parser.add_argument(
        "--output-logical-directory",
        action="append",
        default=[],
        metavar="LABEL=PATH",
        help="Manage a directory output identically in raw PATH and sibling PATH.zip forms",
    )
    parser.add_argument(
        "--optional-output",
        action="append",
        default=[],
        metavar="LABEL=PATH",
        help="Manage an output whose recorded absence is a valid completed result",
    )
    parser.add_argument("--parameter", action="append", default=[], metavar="KEY=VALUE")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    needs_parser = subparsers.add_parser(
        "needs-run", help="Return 0 when a declared artifact must be regenerated and 1 when current"
    )
    add_contract_arguments(needs_parser)
    needs_parser.add_argument(
        "--stale-policy",
        choices=("stop", "reuse", "rebuild"),
        default="stop",
        help="Action when an existing tracked artifact is stale",
    )

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
    audit_parser.add_argument(
        "--stale-policy",
        choices=("stop", "reuse", "rebuild"),
        default="stop",
        help="Allow explicitly requested reuse of readable stale artifacts; structural errors remain fatal",
    )
    return parser


def dispatch(argv: list[str]) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    workspace_root = getattr(args, "workspace_root", None)
    if workspace_root is not None:
        cache_path = Path(workspace_root).absolute() / ".gg_cache" / "content_digests.sqlite3"
        os.environ["GG_CONTENT_DIGEST_CACHE"] = str(cache_path)
        configure_digest_cache(cache_path)
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


def serve() -> int:
    """Serve NUL-delimited argv requests for one workflow shell process."""

    pending = bytearray()
    arguments: list[str] = []
    while chunk := sys.stdin.buffer.read1(65536):
        pending.extend(chunk)
        while True:
            try:
                end = pending.index(0)
            except ValueError:
                break
            token = bytes(pending[:end]).decode("utf-8")
            del pending[: end + 1]
            if token:
                arguments.append(token)
                continue
            if not arguments:
                continue
            # Stdout is reserved for the one-line status protocol. Normal
            # provenance diagnostics remain visible through stderr.
            with contextlib.redirect_stdout(sys.stderr):
                status = dispatch(arguments)
            print(status, flush=True)
            arguments = []
    return 0


def main(argv: list[str] | None = None) -> int:
    resolved_argv = list(sys.argv[1:] if argv is None else argv)
    if resolved_argv == ["serve"]:
        return serve()
    return dispatch(resolved_argv)


if __name__ == "__main__":
    raise SystemExit(main())

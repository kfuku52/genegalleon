#!/usr/bin/env python3
"""Audit or convert all managed ZIP-backed outputs in one GeneGalleon workspace."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
import uuid
from pathlib import Path
from typing import Callable, Optional, Sequence

try:
    from . import gene_family_output_store as FAMILY
    from . import species_tree_output_store as SPECIES
except ImportError:
    import gene_family_output_store as FAMILY
    import species_tree_output_store as SPECIES


REPORT_SCHEMA_VERSION = 1
TARGET_NAMES = ("query2family", "orthogroup", "species_tree")


def _environment_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if raw in {None, ""}:
        return default
    try:
        value = int(raw)
    except ValueError as exc:
        raise WorkspaceStorageError(f"{name} must be an integer") from exc
    if value < 0:
        raise WorkspaceStorageError(f"{name} must be non-negative")
    return value


class WorkspaceStorageError(RuntimeError):
    """Raised when workspace-level storage management cannot proceed safely."""


def _utc_timestamp() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _filename_timestamp() -> str:
    return time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())


def _git_commit() -> str:
    repository = Path(__file__).resolve().parents[2]
    try:
        result = subprocess.run(
            ["git", "-C", str(repository), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    return result.stdout.strip() or "unknown"


def _git_dirty() -> Optional[bool]:
    repository = Path(__file__).resolve().parents[2]
    try:
        result = subprocess.run(
            ["git", "-C", str(repository), "status", "--porcelain"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return bool(result.stdout.strip())


def _managed_children(path: Path) -> bool:
    return any((path / name).exists() for name in TARGET_NAMES + ("orthofinder",))


def resolve_workspace_layout(workspace: Path) -> tuple[Path, Path, Path]:
    workspace = workspace.resolve()
    if not workspace.is_dir():
        raise WorkspaceStorageError(f"Workspace directory was not found: {workspace}")
    if (workspace / "output").is_dir():
        return workspace, (workspace / "output").resolve(), (workspace / "input").resolve()
    nested_workspace = workspace / "workspace"
    if (nested_workspace / "output").is_dir():
        return (
            nested_workspace.resolve(),
            (nested_workspace / "output").resolve(),
            (nested_workspace / "input").resolve(),
        )
    if (workspace / "gfe_data").is_dir() and _managed_children(workspace / "gfe_data"):
        return workspace, (workspace / "gfe_data").resolve(), (workspace / "input").resolve()
    if _managed_children(workspace):
        input_root = (
            workspace.parent / "input"
            if (workspace.parent / "input").is_dir()
            else workspace / "input"
        )
        return workspace, workspace, input_root.resolve()
    raise WorkspaceStorageError(
        "Could not find output/{query2family,orthogroup,species_tree}, gfe_data, "
        f"or direct managed output roots below {workspace}"
    )


def _find_unique(candidates: Sequence[Path], description: str) -> Optional[Path]:
    existing = []
    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved.is_file() and resolved not in existing:
            existing.append(resolved)
    if len(existing) > 1:
        raise WorkspaceStorageError(
            f"Multiple {description} candidates were found; select one explicitly: "
            + ", ".join(str(path) for path in existing)
        )
    return existing[0] if existing else None


def discover_genecount(output_root: Path, override: Optional[Path]) -> Optional[Path]:
    if override is not None:
        path = override.resolve()
        if not path.is_file():
            raise WorkspaceStorageError(f"Gene-count table was not found: {path}")
        return path
    orthofinder = output_root / "orthofinder"
    preferred = [
        orthofinder
        / "Orthogroups_filtered"
        / "Orthogroups.GeneCount.selected.tsv",
        orthofinder / "Orthogroups.GeneCount.selected.tsv",
    ]
    selected = _find_unique(preferred, "selected orthogroup gene-count table")
    if selected is not None:
        return selected
    fallback = sorted(orthofinder.glob("**/Orthogroups.GeneCount.selected.tsv"))
    return _find_unique(fallback, "selected orthogroup gene-count table")


def discover_query_dir(input_root: Path, override: Optional[Path]) -> Optional[Path]:
    if override is not None:
        path = override.resolve()
        if not path.is_dir():
            raise WorkspaceStorageError(f"Query-gene directory was not found: {path}")
        return path
    candidate = input_root / "query_gene"
    return candidate.resolve() if candidate.is_dir() else None


def _target_specs(output_root: Path, selected: set[str]) -> list[dict]:
    specs = []
    for name in TARGET_NAMES:
        if name not in selected:
            continue
        root = output_root / name
        if root.is_dir():
            specs.append(
                {
                    "name": name,
                    "kind": "species-tree" if name == "species_tree" else "gene-family",
                    "mode": name if name != "species_tree" else None,
                    "root": root.resolve(),
                }
            )
    return specs


def _excluded_paths(output_root: Path) -> list[str]:
    paths = []
    orthofinder = output_root / "orthofinder"
    if orthofinder.exists():
        paths.append(str(orthofinder.resolve()))
    genome_evolution = output_root / "genome_evolution"
    if genome_evolution.is_dir():
        paths.extend(
            str(path.resolve())
            for path in sorted(genome_evolution.glob("busco_*"))
        )
    return paths


def _gene_context(
    spec: dict,
    *,
    query_dir: Optional[Path],
    genecount: Optional[Path],
    query_family_id_file: Optional[Path],
    orthogroup_family_id_file: Optional[Path],
    require_catalog: bool,
) -> tuple[list[str], Callable[[str], Optional[str]], list[str]]:
    mode = str(spec["mode"])
    primary = query_dir if mode == "query2family" else genecount
    supplement = (
        query_family_id_file
        if mode == "query2family"
        else orthogroup_family_id_file
    )
    if primary is not None or supplement is not None:
        family_ids, matcher = FAMILY.family_context_with_supplement(
            mode,
            query_dir=query_dir if mode == "query2family" else None,
            genecount=genecount if mode == "orthogroup" else None,
            family_id_file=supplement,
        )
        sources = [
            str(path.resolve())
            for path in (primary, supplement)
            if path is not None
        ]
        return family_ids, matcher, sources

    root = Path(spec["root"])
    if FAMILY._physical_archive_paths(root):
        store = FAMILY.GeneFamilyOutputStore(root)
        store._load_archives()
        family_ids = sorted(
            {
                artifact.family_id
                for artifact in store._archived.values()
                if artifact.family_id is not None
            }
        )
        if mode == "orthogroup":
            matcher = FAMILY.orthogroup_id_from_name
        else:
            matchers = FAMILY.query_id_matchers(family_ids)

            def matcher(name: str) -> Optional[str]:
                return FAMILY.query_id_from_name(name, matchers)

        if not require_catalog or family_ids:
            return family_ids, matcher, []

    raise WorkspaceStorageError(
        f"{mode} requires "
        + (
            "workspace/input/query_gene or --query-dir"
            if mode == "query2family"
            else "Orthogroups.GeneCount.selected.tsv or --genecount"
        )
        + " for raw-file ownership detection"
    )


def _available_space_fields(root: Path, requested: Optional[int]) -> dict:
    filesystem_free = int(shutil.disk_usage(root).free)
    effective = filesystem_free
    if requested is not None:
        effective = min(effective, int(requested))
    return {
        "filesystem_free_bytes": filesystem_free,
        "effective_available_bytes": effective,
        "quota_available_bytes_argument": requested,
    }


def _gene_preflight(
    spec: dict,
    args: argparse.Namespace,
    context: tuple[list[str], Callable[[str], Optional[str]], list[str]],
) -> dict:
    family_ids, matcher, catalog_sources = context
    summary, unmatched = FAMILY.storage_conversion_summary(
        Path(spec["root"]),
        set(family_ids),
        matcher,
        workers=args.workers,
        large_zip_warning_bytes=args.large_zip_warning_bytes,
        max_final_zip_bytes=args.max_final_zip_bytes,
    )
    summary.update(_available_space_fields(Path(spec["root"]), args.available_bytes))
    required_peak_new_bytes = int(
        summary[
            "raw_zip_peak_new_bytes"
            if args.command != "convert" or args.to == "zip"
            else "raw_materialize_allocated_bytes"
        ]
    )
    summary["required_peak_new_bytes"] = required_peak_new_bytes
    summary["temporary_space_sufficient"] = (
        required_peak_new_bytes <= int(summary["effective_available_bytes"])
    )
    summary["catalog_sources"] = catalog_sources
    summary["unmatched_examples"] = [str(path) for path in unmatched[:20]]
    summary["conversion_status"] = FAMILY.storage_conversion_status(
        Path(spec["root"])
    )
    issues = []
    if int(summary["unsupported_symlinks"]):
        issues.append("unsupported-symlinks")
    if int(summary["conflicting_final_zip_files"]):
        issues.append("conflicting-final-zips")
    if args.command == "convert":
        if not summary["temporary_space_sufficient"]:
            issues.append("insufficient-temporary-space")
        if (
            args.to == "zip"
            and args.strict_unmatched
            and int(summary["unmatched_live_files"])
        ):
            issues.append("unmatched-live-files")
    conversion_status = summary["conversion_status"]
    if args.command == "audit":
        if conversion_status["conversion"] != "idle":
            issues.append("conversion-in-progress")
        if int(conversion_status["partial_zip_files"]):
            issues.append("partial-zip-files")
    summary["issues"] = issues
    return summary


def _species_preflight(spec: dict, args: argparse.Namespace) -> dict:
    root = Path(spec["root"])
    records = SPECIES.status(root, SPECIES.MANAGED_DIRECTORIES)
    space = _available_space_fields(root, args.available_bytes)
    for record in records:
        record.update(space)
        required_peak_new_bytes = int(
            record[
                "estimated_peak_new_bytes"
                if args.command != "convert" or args.to == "zip"
                else "raw_materialize_allocated_bytes"
            ]
        )
        record["required_peak_new_bytes"] = required_peak_new_bytes
        record["temporary_space_sufficient"] = (
            required_peak_new_bytes <= int(record["effective_available_bytes"])
        )
    partials = sorted(
        str(path)
        for pattern in (".*.zip.partial.*", ".*.materialize.partial.*")
        for path in root.glob(pattern)
    )
    issues = []
    if any(record["state"] == "mixed" for record in records):
        issues.append("mixed-raw-and-zip")
    if partials:
        issues.append("partial-paths")
    if (
        args.command == "convert"
        and any(not record["temporary_space_sufficient"] for record in records)
    ):
        issues.append("insufficient-temporary-space")
    return {"directories": records, "partial_paths": partials, "issues": issues}


def _workspace_quota_plan(
    report: dict,
    args: argparse.Namespace,
    output_root: Path,
) -> dict:
    """Estimate the initial headroom required across sequential conversions."""

    space = _available_space_fields(output_root, args.available_bytes)
    operations = []
    cumulative_growth = 0
    required_from_initial = 0
    for target in report["targets"]:
        if target.get("status") != "preflight-ok":
            continue
        preflight = target["preflight"]
        if target["kind"] == "gene-family":
            if args.to == "zip":
                temporary_bytes = int(preflight["raw_zip_peak_new_bytes"])
                maximum_net_growth = int(
                    preflight["raw_zip_max_net_growth_bytes"]
                )
            else:
                temporary_bytes = int(
                    preflight["raw_materialize_allocated_bytes"]
                )
                maximum_net_growth = max(
                    0,
                    temporary_bytes - int(preflight["zip_bytes"]),
                )
            operation_specs = [
                (str(target["name"]), temporary_bytes, maximum_net_growth)
            ]
        else:
            operation_specs = []
            for record in preflight["directories"]:
                if args.to == "zip":
                    temporary_bytes = int(record["estimated_peak_new_bytes"])
                    maximum_net_growth = max(
                        0,
                        temporary_bytes - int(record["raw_bytes"]),
                    )
                else:
                    temporary_bytes = int(
                        record["raw_materialize_allocated_bytes"]
                    )
                    maximum_net_growth = max(
                        0,
                        temporary_bytes - int(record["zip_bytes"]),
                    )
                operation_specs.append(
                    (
                        f"{target['name']}/{record['directory']}",
                        temporary_bytes,
                        maximum_net_growth,
                    )
                )
        for name, temporary_bytes, maximum_net_growth in operation_specs:
            operation_required = cumulative_growth + temporary_bytes
            required_from_initial = max(required_from_initial, operation_required)
            operations.append(
                {
                    "operation": name,
                    "temporary_bytes": temporary_bytes,
                    "maximum_prior_net_growth_bytes": cumulative_growth,
                    "required_from_initial_bytes": operation_required,
                    "maximum_net_growth_bytes": maximum_net_growth,
                }
            )
            cumulative_growth += maximum_net_growth
    return {
        **space,
        "required_from_initial_bytes": required_from_initial,
        "maximum_final_net_growth_bytes": cumulative_growth,
        "sufficient": required_from_initial <= int(space["effective_available_bytes"]),
        "operations": operations,
    }


def _verify_gene(
    root: Path,
    verification: str,
    progress_interval: float,
) -> dict:
    if verification == "none":
        return {"mode": "none", "zip_files": 0}
    reporter = FAMILY.ProgressReporter(progress_interval)
    reporter.start()
    try:
        paths = FAMILY.GeneFamilyOutputStore(root).verify(
            progress_callback=reporter.update,
            deep=(verification == "deep"),
        )
    finally:
        reporter.close()
    return {"mode": verification, "zip_files": len(paths)}


def _verify_species(
    root: Path,
    verification: str,
    progress_interval: float,
) -> dict:
    if verification == "none":
        return {"mode": "none", "zip_files": 0, "directories": []}
    reporter = SPECIES.ProgressReporter(progress_interval)
    records, errors = SPECIES._run_each(
        SPECIES.MANAGED_DIRECTORIES,
        lambda name: SPECIES.verify_archive(
            root,
            name,
            check_crc=(verification == "deep"),
            progress_callback=(
                reporter.update if verification == "deep" else None
            ),
        ),
    )
    if errors:
        raise WorkspaceStorageError("; ".join(errors))
    return {
        "mode": verification,
        "zip_files": sum(record["state"] == "archived" for record in records),
        "directories": records,
    }


def _convert_gene(
    spec: dict,
    args: argparse.Namespace,
    context: tuple[list[str], Callable[[str], Optional[str]], list[str]],
) -> dict:
    root = Path(spec["root"])
    family_ids, matcher, catalog_source_strings = context
    reporter = FAMILY.ProgressReporter(args.progress_interval)
    reporter.start()
    try:
        if args.to == "zip":
            return FAMILY.convert_storage_to_zip(
                root,
                str(spec["mode"]),
                family_ids,
                matcher,
                max_files_per_shard=args.max_files_per_shard,
                strict_unmatched=args.strict_unmatched,
                compression=args.compression,
                compression_level=args.compression_level,
                workers=args.workers,
                catalog_sources=[Path(value) for value in catalog_source_strings],
                available_bytes=args.available_bytes,
                large_zip_warning_bytes=args.large_zip_warning_bytes,
                max_final_zip_bytes=args.max_final_zip_bytes,
                progress_callback=reporter.update,
            )
        return FAMILY.convert_storage_to_raw(
            root,
            str(spec["mode"]),
            max_files_per_shard=args.max_files_per_shard,
            pure_raw=args.pure_raw,
            available_bytes=args.available_bytes,
            progress_callback=reporter.update,
        )
    finally:
        reporter.close()


def _convert_species(spec: dict, args: argparse.Namespace) -> dict:
    root = Path(spec["root"])
    reporter = SPECIES.ProgressReporter(args.progress_interval)
    if args.to == "zip":
        records, errors = SPECIES._run_each(
            SPECIES.MANAGED_DIRECTORIES,
            lambda name: SPECIES.pack_directory(
                root,
                name,
                compression=args.compression,
                compression_level=args.compression_level,
                available_bytes=args.available_bytes,
                progress_callback=reporter.update,
                include_metrics=True,
            ),
        )
    else:
        records, errors = SPECIES._run_each(
            SPECIES.MANAGED_DIRECTORIES,
            lambda name: SPECIES.materialize_directory(
                root,
                name,
                available_bytes=args.available_bytes,
                progress_callback=reporter.update,
                include_metrics=True,
            ),
        )
    if errors:
        raise WorkspaceStorageError("; ".join(errors))
    return {"directories": records}


def _atomic_write(path: Path, contents: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial.{os.getpid()}.{uuid.uuid4().hex}")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            handle.write(contents)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _metric(payload: object, key: str, default: object = "") -> object:
    return payload.get(key, default) if isinstance(payload, dict) else default


def _tsv_rows(report: dict) -> list[dict]:
    rows = []
    for target in report["targets"]:
        before = target.get("preflight", {})
        after = target.get("postflight", {})
        if target["kind"] == "species-tree":
            before_dirs = before.get("directories", [])
            after_dirs = after.get("directories", [])
            before_logical_files = sum(
                max(int(row.get("raw_files", 0)), int(row.get("archived_files", 0)))
                for row in before_dirs
            )
            after_physical_files = sum(
                int(row.get("raw_files", 0))
                + int(row.get("state") == "archived")
                for row in after_dirs
            )
            before_logical_bytes = sum(
                int(row.get("raw_bytes", 0))
                or int(row.get("archived_logical_bytes", 0))
                for row in before_dirs
            )
            after_physical_bytes = sum(
                int(row.get("raw_bytes", 0)) + int(row.get("zip_bytes", 0))
                for row in after_dirs
            )
            issues = after.get("issues", [])
        else:
            before_logical_files = int(before.get("managed_logical_files", 0))
            before_logical_bytes = int(before.get("managed_logical_bytes", 0))
            after_physical_files = int(after.get("physical_managed_files", 0))
            after_physical_bytes = int(after.get("physical_managed_bytes", 0))
            issues = after.get("issues", [])
        rows.append(
            {
                "target": target["name"],
                "kind": target["kind"],
                "root": target["root"],
                "status": target.get("status", "unknown"),
                "before_managed_logical_files": before_logical_files,
                "after_managed_physical_files": after_physical_files,
                "before_managed_logical_bytes": before_logical_bytes,
                "after_managed_physical_bytes": after_physical_bytes,
                "elapsed_seconds": target.get("elapsed_seconds", 0),
                "verification": _metric(target.get("verification", {}), "mode", "none"),
                "issues": ",".join(str(value) for value in issues),
                "error": target.get("error", ""),
            }
        )
    totals = report.get("totals")
    if isinstance(totals, dict):
        rows.append(
            {
                "target": "TOTAL",
                "kind": "workspace",
                "root": report["output_root"],
                "status": report["status"],
                "before_managed_logical_files": totals[
                    "before_managed_logical_files"
                ],
                "after_managed_physical_files": totals[
                    "after_managed_physical_files"
                ],
                "before_managed_logical_bytes": totals[
                    "before_managed_logical_bytes"
                ],
                "after_managed_physical_bytes": totals[
                    "after_managed_physical_bytes"
                ],
                "elapsed_seconds": report.get("elapsed_seconds", 0),
                "verification": report.get("verification", "none"),
                "issues": "",
                "error": "",
            }
        )
    return rows


def _report_totals(report: dict) -> dict:
    report_without_totals = dict(report)
    report_without_totals.pop("totals", None)
    rows = _tsv_rows(report_without_totals)
    before_files = sum(int(row["before_managed_logical_files"]) for row in rows)
    after_files = sum(int(row["after_managed_physical_files"]) for row in rows)
    before_bytes = sum(int(row["before_managed_logical_bytes"]) for row in rows)
    after_bytes = sum(int(row["after_managed_physical_bytes"]) for row in rows)
    return {
        "before_managed_logical_files": before_files,
        "after_managed_physical_files": after_files,
        "file_reduction": before_files - after_files,
        "file_reduction_percent": (
            round((before_files - after_files) * 100 / before_files, 4)
            if before_files
            else 0.0
        ),
        "before_managed_logical_bytes": before_bytes,
        "after_managed_physical_bytes": after_bytes,
        "byte_reduction": before_bytes - after_bytes,
        "byte_reduction_percent": (
            round((before_bytes - after_bytes) * 100 / before_bytes, 4)
            if before_bytes
            else 0.0
        ),
    }


def write_report(report: dict, report_dir: Path, prefix: Optional[str]) -> tuple[Path, Path]:
    report_dir = report_dir.resolve()
    stem = prefix or f"gg_storage_{_filename_timestamp()}_{os.getpid()}"
    if "/" in stem or stem in {"", ".", ".."}:
        raise WorkspaceStorageError("--report-prefix must be a filename stem")
    json_path = report_dir / f"{stem}.json"
    tsv_path = report_dir / f"{stem}.tsv"
    _atomic_write(json_path, json.dumps(report, indent=2, sort_keys=True) + "\n")
    rows = _tsv_rows(report)
    columns = list(rows[0]) if rows else [
        "target",
        "kind",
        "root",
        "status",
        "before_managed_logical_files",
        "after_managed_physical_files",
        "before_managed_logical_bytes",
        "after_managed_physical_bytes",
        "elapsed_seconds",
        "verification",
        "issues",
        "error",
    ]
    temporary = tsv_path.with_name(
        f".{tsv_path.name}.partial.{os.getpid()}.{uuid.uuid4().hex}"
    )
    try:
        temporary.parent.mkdir(parents=True, exist_ok=True)
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, tsv_path)
    finally:
        temporary.unlink(missing_ok=True)
    return json_path, tsv_path


def build_parser() -> argparse.ArgumentParser:
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=formatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add_common(subparser: argparse.ArgumentParser) -> None:
        subparser.add_argument(
            "--workspace",
            required=True,
            type=Path,
            help=(
                "Project root, workspace root, output root, or legacy gfe_data "
                "root to discover and manage"
            ),
        )
        subparser.add_argument(
            "--target",
            action="append",
            choices=TARGET_NAMES,
            help="Manage only this target; repeat the option to select several",
        )
        subparser.add_argument(
            "--query-dir",
            type=Path,
            help="Override the discovered input/query_gene family catalog",
        )
        subparser.add_argument(
            "--genecount",
            type=Path,
            help="Override the discovered Orthogroups.GeneCount.selected.tsv",
        )
        subparser.add_argument(
            "--query-family-id-file",
            type=Path,
            help="Add query family IDs listed one per line to the current catalog",
        )
        subparser.add_argument(
            "--orthogroup-family-id-file",
            type=Path,
            help="Add orthogroup IDs listed one per line to the current catalog",
        )
        subparser.add_argument(
            "--verification",
            choices=("quick", "deep", "none"),
            default="quick",
            help=(
                "Postflight ZIP verification depth; deep reads every payload, "
                "quick checks ZIP/index metadata, and none skips verification"
            ),
        )
        subparser.add_argument(
            "--progress-interval",
            type=float,
            default=30.0,
            help="Minimum seconds between periodic progress messages; 0 is unthrottled",
        )
        subparser.add_argument(
            "--workers",
            type=int,
            default=1,
            help="Concurrent gene-family ZIP writers (1-4)",
        )
        subparser.add_argument(
            "--large-zip-warning-bytes",
            type=int,
            default=_environment_int(
                "GG_COMMON_GENE_FAMILY_LARGE_ZIP_WARNING_BYTES",
                FAMILY.DEFAULT_LARGE_ZIP_WARNING_BYTES,
            ),
            help="Report gene-family subdirectories above this logical size; 0 disables",
        )
        subparser.add_argument(
            "--max-final-zip-bytes",
            type=int,
            default=_environment_int(
                "GG_COMMON_GENE_FAMILY_FINAL_ZIP_MAX_BYTES",
                0,
            ),
            help=(
                "Retain readable subdir.part-N.zip shards above this logical "
                "size; 0 permits one final ZIP of any size"
            ),
        )
        subparser.add_argument(
            "--available-bytes",
            type=int,
            help=(
                "Known remaining user/project quota; effective space is the "
                "smaller of this value and filesystem free space"
            ),
        )
        subparser.add_argument(
            "--report-dir",
            type=Path,
            help="Directory for the atomic JSON and TSV reports",
        )
        subparser.add_argument(
            "--report-prefix",
            help="Filename stem for both reports (without .json or .tsv)",
        )

    audit = subparsers.add_parser(
        "audit",
        help="Inspect all managed stores without converting them",
        formatter_class=formatter,
    )
    add_common(audit)

    convert = subparsers.add_parser(
        "convert",
        help="Preflight, convert, verify, and report all managed stores",
        formatter_class=formatter,
    )
    add_common(convert)
    convert.add_argument(
        "--to",
        required=True,
        choices=("zip", "raw", "files"),
        help="Destination storage mode; files is an alias for raw",
    )
    convert.add_argument(
        "--dry-run",
        action="store_true",
        help="Run all preflights and write reports without changing managed storage",
    )
    convert.add_argument(
        "--compression",
        choices=("adaptive", "deflate", "store"),
        default="adaptive",
        help="Compression used for newly created ZIP members",
    )
    convert.add_argument(
        "--compression-level",
        type=int,
        default=6,
        help="Deflate level for newly created ZIP members (0-9)",
    )
    convert.add_argument(
        "--max-files-per-shard",
        type=int,
        default=5000,
        help="Maximum member count in each intermediate or retained ZIP shard",
    )
    convert.add_argument(
        "--strict-unmatched",
        action="store_true",
        help="Block ZIP conversion instead of preserving unmatched gene-family files",
    )
    convert.add_argument(
        "--pure-raw",
        action="store_true",
        help="With --to raw, also remove generated gene-family archive metadata/history",
    )
    convert.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Convert targets whose preflight passed even if another target failed",
    )
    return parser


def _validate_args(args: argparse.Namespace) -> None:
    if args.progress_interval < 0:
        raise WorkspaceStorageError("--progress-interval must be non-negative")
    if not 1 <= args.workers <= 4:
        raise WorkspaceStorageError("--workers must be between 1 and 4")
    for name in ("large_zip_warning_bytes", "max_final_zip_bytes", "available_bytes"):
        value = getattr(args, name)
        if value is not None and value < 0:
            raise WorkspaceStorageError(f"--{name.replace('_', '-')} must be non-negative")
    if args.command == "convert":
        args.to = "raw" if args.to == "files" else args.to
        if not 0 <= args.compression_level <= 9:
            raise WorkspaceStorageError("--compression-level must be between 0 and 9")
        if args.max_files_per_shard < 1:
            raise WorkspaceStorageError("--max-files-per-shard must be positive")
        if args.to == "zip" and args.pure_raw:
            raise WorkspaceStorageError("--pure-raw is valid only with --to raw")
        if args.to == "raw" and args.strict_unmatched:
            raise WorkspaceStorageError("--strict-unmatched is valid only with --to zip")


def run(args: argparse.Namespace) -> tuple[int, dict, Path, Path]:
    _validate_args(args)
    workspace, output_root, input_root = resolve_workspace_layout(args.workspace)
    selected = set(args.target or TARGET_NAMES)
    specs = _target_specs(output_root, selected)
    if not specs:
        raise WorkspaceStorageError(
            f"No selected managed output roots were found below {output_root}"
        )
    query_dir = discover_query_dir(input_root, args.query_dir)
    genecount = discover_genecount(output_root, args.genecount)
    started = time.monotonic()
    report = {
        "schema_version": REPORT_SCHEMA_VERSION,
        "command": args.command,
        "requested_target_storage": (
            args.to if args.command == "convert" else None
        ),
        "dry_run": bool(getattr(args, "dry_run", False)),
        "status": "running",
        "started_utc": _utc_timestamp(),
        "workspace": str(workspace),
        "output_root": str(output_root),
        "input_root": str(input_root),
        "genegalleon_commit": _git_commit(),
        "genegalleon_dirty": _git_dirty(),
        "verification": args.verification,
        "excluded_paths": _excluded_paths(output_root),
        "options": {
            "workers": args.workers,
            "available_bytes": args.available_bytes,
            "large_zip_warning_bytes": args.large_zip_warning_bytes,
            "max_final_zip_bytes": args.max_final_zip_bytes,
            "strict_unmatched": bool(getattr(args, "strict_unmatched", False)),
            "pure_raw": bool(getattr(args, "pure_raw", False)),
        },
        "targets": [],
    }
    contexts: dict[str, tuple[list[str], Callable[[str], Optional[str]], list[str]]] = {}
    preflight_failed = False
    for spec in specs:
        target = {
            "name": spec["name"],
            "kind": spec["kind"],
            "mode": spec["mode"],
            "root": str(spec["root"]),
            "status": "preflight",
        }
        report["targets"].append(target)
        try:
            if spec["kind"] == "gene-family":
                context = _gene_context(
                    spec,
                    query_dir=query_dir,
                    genecount=genecount,
                    query_family_id_file=args.query_family_id_file,
                    orthogroup_family_id_file=args.orthogroup_family_id_file,
                    require_catalog=(
                        args.command == "convert"
                        and getattr(args, "to", None) == "zip"
                    ),
                )
                contexts[str(spec["name"])] = context
                target["preflight"] = _gene_preflight(spec, args, context)
            else:
                target["preflight"] = _species_preflight(spec, args)
            issues = target["preflight"].get("issues", [])
            blocking_issues = set(issues)
            if blocking_issues:
                raise WorkspaceStorageError(
                    "preflight issues: " + ", ".join(sorted(blocking_issues))
                )
            target["status"] = "preflight-ok"
        except Exception as exc:
            target["status"] = "failed"
            target["error"] = str(exc)
            preflight_failed = True

    if args.command == "convert":
        quota_plan = _workspace_quota_plan(report, args, output_root)
        report["quota_plan"] = quota_plan
        if not quota_plan["sufficient"]:
            preflight_failed = True
            message = (
                "workspace-wide temporary-space estimate exceeds available "
                f"headroom: required_from_initial_bytes="
                f"{quota_plan['required_from_initial_bytes']}, "
                f"effective_available_bytes="
                f"{quota_plan['effective_available_bytes']}"
            )
            for target in report["targets"]:
                if target["status"] != "preflight-ok":
                    continue
                target["preflight"].setdefault("issues", []).append(
                    "insufficient-workspace-temporary-space"
                )
                target["status"] = "failed"
                target["error"] = message

    should_convert = (
        args.command == "convert"
        and not args.dry_run
        and (not preflight_failed or args.continue_on_error)
    )
    if args.command == "convert" and preflight_failed and not args.continue_on_error:
        for target in report["targets"]:
            if target["status"] == "preflight-ok":
                target["status"] = "not-run"

    for spec, target in zip(specs, report["targets"]):
        if target["status"] in {"failed", "not-run"}:
            continue
        if args.command == "convert" and not args.dry_run and not should_convert:
            target["status"] = "not-run"
            target["elapsed_seconds"] = 0
            continue
        target_started = time.monotonic()
        try:
            if should_convert:
                if spec["kind"] == "gene-family":
                    target["conversion"] = _convert_gene(
                        spec,
                        args,
                        contexts[str(spec["name"])],
                    )
                else:
                    target["conversion"] = _convert_species(spec, args)
            if args.command == "audit" or should_convert:
                target["verification"] = (
                    _verify_gene(
                        Path(spec["root"]),
                        args.verification,
                        args.progress_interval,
                    )
                    if spec["kind"] == "gene-family"
                    else _verify_species(
                        Path(spec["root"]),
                        args.verification,
                        args.progress_interval,
                    )
                )
                if spec["kind"] == "gene-family":
                    target["postflight"] = _gene_preflight(
                        spec,
                        args,
                        contexts[str(spec["name"])],
                    )
                else:
                    target["postflight"] = _species_preflight(spec, args)
            else:
                target["postflight"] = target["preflight"]
            target["status"] = (
                "dry-run" if args.command == "convert" and args.dry_run else "complete"
            )
        except Exception as exc:
            target["status"] = "failed"
            target["error"] = str(exc)
            if args.command == "convert" and not args.continue_on_error:
                should_convert = False
        finally:
            target["elapsed_seconds"] = round(time.monotonic() - target_started, 3)

    report["finished_utc"] = _utc_timestamp()
    report["elapsed_seconds"] = round(time.monotonic() - started, 3)
    failed = any(target["status"] in {"failed", "not-run"} for target in report["targets"])
    report["status"] = "failed" if failed else "complete"
    report["totals"] = _report_totals(report)
    report_dir = (
        args.report_dir.resolve()
        if args.report_dir is not None
        else (output_root.parent / "storage_reports").resolve()
    )
    json_path, tsv_path = write_report(report, report_dir, args.report_prefix)
    return (1 if failed else 0), report, json_path, tsv_path


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        exit_code, report, json_path, tsv_path = run(args)
    except (WorkspaceStorageError, FAMILY.ArchiveStoreError, SPECIES.SpeciesTreeArchiveError, OSError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    print(f"status\t{report['status']}")
    print(f"report_json\t{json_path}")
    print(f"report_tsv\t{tsv_path}")
    for target in report["targets"]:
        print(
            "target\t"
            f"{target['name']}\t{target['status']}\t"
            f"{target.get('error', '')}"
        )
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())

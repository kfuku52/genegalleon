#!/usr/bin/env python3
"""Download an array's inputs together and freeze local tasks before workers start."""

import argparse
import json
import sys
from pathlib import Path

import format_species_inputs as fsi
from format_species_manifest import resolved_manifest_fieldnames, write_resolved_manifest_tsv
from format_species_provider_config import DEFAULT_INPUT_RELATIVE_DIRS
from input_generation_array_state import atomic_json, digest, export_manifest, load_plan


def has_genome_input(task):
    return any(
        path and Path(path).is_file() and Path(path).stat().st_size > 0
        for path in (task.get("genome_path"), task.get("gbff_path"))
    )


def stage_downloads(plan_path, *, jobs=4, timeout=120, headers=None, require_genome=False):
    plan_path = Path(plan_path).resolve()
    plan = load_plan(plan_path)
    if plan.get("download_mode") != "staged" or any("manifest_row" not in task for task in plan["tasks"]):
        raise ValueError("Staging requires a manifest plan created with --stage-downloads")
    plan_hash = digest(plan_path)
    task_root = Path(str(plan_path) + ".tasks")
    task_root.mkdir(parents=True, exist_ok=True)
    pending = []
    for index, task in enumerate(plan["tasks"], 1):
        for path, expected in task.get("input_sha256", {}).items():
            if digest(path) != expected:
                raise ValueError("Local manifest input changed after planning: " + path)
        cached_path = task_root / f"{index}.json"
        if cached_path.exists():
            cached = json.loads(cached_path.read_text())
            if cached.get("plan_sha256") != plan_hash or cached.get("task_index") != index:
                raise ValueError("Staged input belongs to another plan/task")
            for path, expected in cached["task"]["input_sha256"].items():
                if digest(path) != expected:
                    raise ValueError("Staged raw input changed; use a new workspace: " + path)
            if not (task_root / f"{index}.resolved.tsv").is_file():
                raise ValueError("Staged resolved manifest is missing; use a new workspace")
            if require_genome and not has_genome_input(cached["task"]):
                raise ValueError("Required genome input is missing for " + task["species_prefix"])
        else:
            pending.append((index, task))
    if not pending:
        print("All staged inputs verified; no downloads needed.")
        return

    manifest = task_root / "download_manifest.tsv"
    pending_tasks = [task for _index, task in pending]
    export_manifest(plan, manifest, tasks=pending_tasks)
    # A dedicated immutable-plan cache avoids overwriting another plan's inputs.
    download_root = Path(plan["tasks"][0]["download_dir"]) / "staged" / plan_hash
    validation_cache_dir = Path(plan["tasks"][0]["download_dir"]) / "staged" / ".gg-gzip-validation"
    report = fsi.download_from_manifest(
        manifest_path=manifest, download_root=download_root, provider_filter=plan["provider"],
        overwrite=False, headers=headers or {}, timeout=timeout, dry_run=False, jobs=jobs,
        validation_cache_dir=validation_cache_dir)
    for warning in report["warnings"]:
        print("Warning: " + warning, file=sys.stderr)
    # Unscoped validation/merge errors cannot certify any pending species.
    # Keep downloads for retry, but never freeze inputs that failed validation.
    if report["errors"]:
        raise ValueError("Staged 0 of {} pending species before failure: {}".format(
            len(pending), "; ".join(report["errors"])))
    resolved_rows = {(row["provider"], row["species_key"]): row for row in report["resolved_rows"]}
    discovered = {}
    discovery_errors = []
    failed_providers = set()
    for provider in sorted({task["provider"] for task in plan["tasks"]}):
        allowed_species_keys = {
            task["species_key"] for task in plan["tasks"] if task["provider"] == provider
        }
        tasks, warnings, errors = fsi.discover_tasks(
            provider,
            download_root / DEFAULT_INPUT_RELATIVE_DIRS[provider],
            allowed_species_keys=allowed_species_keys,
        )
        for warning in warnings:
            print("Warning: " + warning, file=sys.stderr)
        discovery_errors.extend(errors)
        if errors:
            failed_providers.add(provider)
            print("Warning: partial staged discovery for {}: {}".format(provider, "; ".join(errors)), file=sys.stderr)
        for task in tasks:
            key = (provider, task["species_prefix"])
            if key in discovered:
                raise ValueError("Duplicate staged species: " + repr(key))
            discovered[key] = task

    staged_count = 0
    staging_errors = []
    for index, task in pending:
        key = (task["provider"], task["species_prefix"])
        if key not in discovered or key not in resolved_rows:
            staging_errors.append("Missing staged species: " + repr(key))
            continue
        if task["provider"] in failed_providers:
            continue
        actual = discovered[key]
        if require_genome and not has_genome_input(actual):
            staging_errors.append("Required genome input is missing for " + task["species_prefix"])
            continue
        actual["input_sha256"] = {
            **task.get("input_sha256", {}),
            **{str(actual[k]): digest(actual[k]) for k in ("cds_path", "gff_path", "gbff_path", "genome_path") if actual.get(k)},
        }
        for path, expected in task.get("input_sha256", {}).items():
            if digest(path) != expected:
                raise ValueError("Local input changed during staging: " + path)
        for setting in ("gene_grouping_mode", "gff_repair_mode", "format_strict"):
            actual[setting] = task[setting]
        rows = [resolved_rows[key]]
        write_resolved_manifest_tsv(task_root / f"{index}.resolved.tsv", resolved_manifest_fieldnames(rows), rows)
        atomic_json(task_root / f"{index}.json", {
            "plan_sha256": plan_hash, "task_index": index,
            "task": {k: str(v) if isinstance(v, Path) else v for k, v in actual.items()},
        }, immutable=True)
        staged_count += 1

    failures = list(report["errors"]) + discovery_errors + staging_errors
    if failures:
        raise ValueError(
            "Staged {} of {} pending species before failure: {}".format(
                staged_count, len(pending), "; ".join(failures)
            )
        )
    print(f"Staged {staged_count} species; downloaded {report['downloaded']} files. Workers use verified local inputs.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-plan", required=True, type=Path)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--download-timeout", type=float, default=120)
    parser.add_argument("--http-header", action="append", default=[])
    parser.add_argument("--auth-bearer-token-env", default="")
    parser.add_argument("--require-genome", action="store_true")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    stage_downloads(args.task_plan, jobs=args.jobs, timeout=args.download_timeout, require_genome=args.require_genome,
                    headers=fsi.parse_http_headers(args.http_header, args.auth_bearer_token_env))


if __name__ == "__main__":
    main()

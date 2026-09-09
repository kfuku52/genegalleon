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


def stage_downloads(plan_path, *, jobs=4, timeout=120, headers=None):
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
        else:
            pending.append((index, task))
    if not pending:
        print("All staged inputs verified; no downloads needed.")
        return

    manifest = task_root / "download_manifest.tsv"
    export_manifest(plan, manifest)
    # A dedicated immutable-plan cache avoids overwriting another plan's inputs.
    download_root = Path(plan["tasks"][0]["download_dir"]) / "staged" / plan_hash
    report = fsi.download_from_manifest(
        manifest_path=manifest, download_root=download_root, provider_filter=plan["provider"],
        overwrite=False, headers=headers or {}, timeout=timeout, dry_run=False, jobs=jobs)
    for warning in report["warnings"]:
        print("Warning: " + warning, file=sys.stderr)
    if report["errors"]:
        raise ValueError("; ".join(report["errors"]))
    resolved_rows = {(row["provider"], row["species_key"]): row for row in report["resolved_rows"]}
    discovered = {}
    for provider in sorted({task["provider"] for task in plan["tasks"]}):
        tasks, warnings, errors = fsi.discover_tasks(provider, download_root / DEFAULT_INPUT_RELATIVE_DIRS[provider])
        if errors:
            raise ValueError("; ".join(errors))
        for task in tasks:
            key = (provider, task["species_prefix"])
            if key in discovered:
                raise ValueError("Duplicate staged species: " + repr(key))
            discovered[key] = task

    for index, task in pending:
        key = (task["provider"], task["species_prefix"])
        if key not in discovered or key not in resolved_rows:
            raise ValueError("Missing staged species: " + repr(key))
        actual = discovered[key]
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
    print(f"Staged {len(pending)} species; downloaded {report['downloaded']} files. Workers use verified local inputs.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-plan", required=True, type=Path)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--download-timeout", type=float, default=120)
    parser.add_argument("--http-header", action="append", default=[])
    parser.add_argument("--auth-bearer-token-env", default="")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    stage_downloads(args.task_plan, jobs=args.jobs, timeout=args.download_timeout,
                    headers=fsi.parse_http_headers(args.http_header, args.auth_bearer_token_env))


if __name__ == "__main__":
    main()

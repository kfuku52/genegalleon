#!/usr/bin/env python3
"""Download an array's inputs together and freeze local tasks before workers start."""

import argparse
import json
import sys
from pathlib import Path
from urllib.parse import unquote, urlparse

import format_species_inputs as fsi
from format_species_annotation.tasks import task_missing_annotation_label
from format_species_common import is_fasta_filename, is_gbff_filename, is_gff_filename
from format_species_download.local import validate_gzip_with_cache
from format_species_manifest import resolved_manifest_fieldnames, write_resolved_manifest_tsv
from format_species_provider_config import DEFAULT_INPUT_RELATIVE_DIRS
from format_species_provider_resolvers import provider_raw_dir
from format_species_providers.catalogs import validate_coge_export_gff_file
from input_generation_array_state import atomic_json, digest, digest_paths, export_manifest, load_plan
from input_generation_staging_reuse import FreshReadFence, StagedProofReader


def has_required_source(task, keys):
    return any(
        path and Path(path).is_file() and Path(path).stat().st_size > 0
        for path in (task.get(key) for key in keys)
    )


def explicit_manifest_task(task, row, download_root):
    """Bind resolved manifest roles without inferring species from filenames."""
    provider = task["provider"]
    roles = ("cds", "gff", "gbff", "genome")
    if not any((row.get(role + "_url") or "").strip() for role in roles):
        return None
    species_key = task["species_key"]
    if species_key in ("", ".", "..") or Path(species_key).name != species_key:
        raise ValueError("Unsafe manifest species key: " + species_key)
    raw_dir = provider_raw_dir(provider, download_root, species_key)
    actual = {
        "provider": provider, "species_key": species_key,
        "species_prefix": task["species_prefix"],
        "gff_auto_selected_from_multiple": False,
        "gff_selection_candidates": (),
    }
    validators = {"cds": is_fasta_filename, "gff": is_gff_filename,
                  "gbff": is_gbff_filename, "genome": is_fasta_filename}
    for role in roles:
        url = (row.get(role + "_url") or "").strip()
        filename = (row.get(role + "_filename") or "").strip()
        path = None
        if url:
            if not filename or filename in (".", "..") or Path(filename).name != filename or not validators[role](filename):
                raise ValueError("Invalid explicit {} filename for {}".format(role, species_key))
            path = raw_dir / filename
            if path.is_symlink():
                raise ValueError("Explicit {} input is a symlink for {}".format(role, species_key))
            if not path.is_file() or path.stat().st_size == 0:
                path = None
        actual[role + "_path"] = path
    missing = task_missing_annotation_label(
        actual["cds_path"], actual["gff_path"], actual["gbff_path"], actual["genome_path"])
    if missing:
        raise ValueError("[{}] {}: missing {}".format(provider, species_key, missing))
    if provider == "coge" and actual["gff_path"] is not None:
        validate_coge_export_gff_file(actual["gff_path"], gid=row.get("id", ""))
    actual["gff_selection_candidates"] = (
        (actual["gff_path"].name,) if actual["gff_path"] is not None else ()
    )
    return actual


def bound_local_manifest_task(task, verified_inputs=None):
    """Explicitly bind frozen file sources without creating another raw copy."""
    row = task["manifest_row"]
    mode = str(row.get("bind_local_sources", "") or "").strip()
    if mode in ("", "0"):
        return None
    if mode != "1":
        raise ValueError("bind_local_sources must be 0 or 1")
    actual = {"provider": task["provider"], "species_key": task["species_key"],
              "species_prefix": task["species_prefix"], "gff_auto_selected_from_multiple": False}
    for role in ("cds", "gff", "gbff", "genome"):
        url = row.get(role + "_url", "") or ""
        path = None
        if row.get(role + "_archive_member"):
            raise ValueError("Bound local sources cannot select archive members")
        if url:
            parsed = urlparse(url)
            if parsed.scheme != "file" or parsed.netloc or parsed.query or parsed.fragment:
                raise ValueError("Bound local sources require absolute file URLs for every supplied role")
            source = Path(unquote(parsed.path))
            if not source.is_absolute() or source.is_symlink():
                raise ValueError("Unsafe bound local source: " + str(source))
            path = source.resolve()
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError("Missing or empty bound local source: " + str(path))
            expected = task.get("input_sha256", {}).get(str(path))
            if not expected:
                raise ValueError("Bound local source does not match the frozen plan: " + str(path))
        actual[role + "_path"] = path
    paths = [actual[key] for key in ("cds_path", "gff_path", "gbff_path", "genome_path") if actual[key]]
    observed_hashes = verified_inputs.verified(paths) if verified_inputs else digest_paths(paths)
    for path, observed in observed_hashes.items():
        if observed != task["input_sha256"][path]:
            raise ValueError("Bound local source does not match the frozen plan: " + path)
        error = None if verified_inputs else validate_gzip_with_cache(Path(path))
        if error is not None:
            raise ValueError("Invalid bound local source: {} ({})".format(path, error))
    missing = task_missing_annotation_label(actual["cds_path"], actual["gff_path"],
                                            actual["gbff_path"], actual["genome_path"])
    if missing:
        raise ValueError("Bound local task is missing " + missing)
    if task["provider"] == "coge" and actual["gff_path"] is not None:
        validate_coge_export_gff_file(actual["gff_path"], gid=row.get("id", ""))
    actual["gff_selection_candidates"] = (actual["gff_path"].name,) if actual["gff_path"] else ()
    return actual


def stage_downloads(plan_path, *, jobs=4, timeout=120, headers=None, require_gff=False, require_genome=False):
    plan_path = Path(plan_path).resolve()
    plan = load_plan(plan_path)
    if plan.get("download_mode") != "staged" or any("manifest_row" not in task for task in plan["tasks"]):
        raise ValueError("Staging requires a manifest plan created with --stage-downloads")
    plan_hash = digest(plan_path)
    task_root = Path(str(plan_path) + ".tasks")
    task_root.mkdir(parents=True, exist_ok=True)
    pending = []
    reuse_reader = StagedProofReader()
    reuse_fences = {}
    for index, task in enumerate(plan["tasks"], 1):
        original_hashes = task.get("input_sha256", {})
        cached_path = task_root / f"{index}.json"
        cached = None
        if cached_path.exists():
            cached = json.loads(cached_path.read_text())
            if cached.get("plan_sha256") != plan_hash or cached.get("task_index") != index:
                raise ValueError("Staged input belongs to another plan/task")
        # A resumed receipt often names the original bound sources again.
        # Share this full read only within preflight; binding/publication still
        # perform their independent content checks after intervening work.
        cached_hashes = cached["task"]["input_sha256"] if cached is not None else {}
        reuse = reuse_reader.resolve(task)
        if reuse != task.get("staged_input_reuse"):
            raise ValueError("Staging reuse proof differs from the frozen plan")
        if reuse:
            if original_hashes != reuse["input_sha256"]:
                raise ValueError("Staging reuse hashes differ from the frozen plan")
            reuse_fences[index] = FreshReadFence([*original_hashes, *cached_hashes])
        observed_hashes = digest_paths([*original_hashes, *cached_hashes])
        if reuse:
            if cached_hashes and cached_hashes != original_hashes:
                raise ValueError("Staging reuse cache hashes differ from the frozen plan")
            reuse_fences[index].certify(observed_hashes, original_hashes)
        for path, expected in original_hashes.items():
            if observed_hashes[path] != expected:
                raise ValueError("Local manifest input changed after planning: " + path)
        if cached is not None:
            for path, expected in cached["task"]["input_sha256"].items():
                if observed_hashes[path] != expected:
                    raise ValueError("Staged raw input changed; use a new workspace: " + path)
            if not (task_root / f"{index}.resolved.tsv").is_file():
                raise ValueError("Staged resolved manifest is missing; use a new workspace")
            if require_gff and not has_required_source(cached["task"], ("gff_path", "gbff_path")):
                raise ValueError("Required GFF input is missing for " + task["species_prefix"])
            if require_genome and not has_required_source(cached["task"], ("genome_path", "gbff_path")):
                raise ValueError("Required genome input is missing for " + task["species_prefix"])
        else:
            pending.append((index, task))
    if not pending:
        reuse_reader.check()
        for fence in reuse_fences.values():
            fence.check()
        print("All staged inputs verified; no downloads needed.")
        return

    manifest = task_root / "download_manifest.tsv"
    discovered = {}
    bound = {}
    pending_tasks = []
    for _index, task in pending:
        fence = reuse_fences.get(_index)
        actual = bound_local_manifest_task(task, verified_inputs=fence) if fence else bound_local_manifest_task(task)
        key = (task["provider"], task["species_prefix"])
        if actual is None:
            pending_tasks.append(task)
        else:
            discovered[key] = actual
            bound[key] = task["manifest_row"]
    export_manifest(plan, manifest, tasks=pending_tasks)
    # A dedicated immutable-plan cache avoids overwriting another plan's inputs.
    download_root = Path(plan["tasks"][0]["download_dir"]) / "staged" / plan_hash
    validation_cache_dir = Path(plan["tasks"][0]["download_dir"]) / "staged" / ".gg-gzip-validation"
    report = fsi.download_from_manifest(
        manifest_path=manifest, download_root=download_root, provider_filter=plan["provider"],
        overwrite=False, headers=headers or {}, timeout=timeout, dry_run=False, jobs=jobs,
        validation_cache_dir=validation_cache_dir) if pending_tasks else {
            "warnings": [], "errors": [], "resolved_rows": [], "downloaded": 0}
    for warning in report["warnings"]:
        print("Warning: " + warning, file=sys.stderr)
    # Unscoped validation/merge errors cannot certify any pending species.
    # Keep downloads for retry, but never freeze inputs that failed validation.
    if report["errors"]:
        raise ValueError("Staged 0 of {} pending species before failure: {}".format(
            len(pending), "; ".join(report["errors"])))
    resolved_rows = {(row["provider"], row["species_key"]): row for row in report["resolved_rows"]}
    resolved_rows.update(bound)
    discovery_errors = []
    failed_providers = set()
    explicit_errors = {}
    explicit_keys = set()
    fallback_species = {}
    for _index, task in pending:
        provider = task["provider"]
        key = (provider, task["species_prefix"])
        if key in bound:
            continue
        row = resolved_rows.get((provider, task["species_key"]))
        try:
            actual = explicit_manifest_task(task, row or {}, download_root)
        except ValueError as exc:
            explicit_errors[key] = str(exc)
            continue
        if actual is None:
            fallback_species.setdefault(provider, set()).add(task["species_key"])
        else:
            discovered[key] = actual
            explicit_keys.add(key)
    for provider, allowed_species_keys in sorted(fallback_species.items()):
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
        if key in explicit_errors:
            staging_errors.append(explicit_errors[key])
            continue
        if key not in discovered or key not in resolved_rows:
            staging_errors.append("Missing staged species: " + repr(key))
            continue
        if task["provider"] in failed_providers and key not in explicit_keys:
            continue
        actual = discovered[key]
        if require_gff and not has_required_source(actual, ("gff_path", "gbff_path")):
            staging_errors.append("Required GFF input is missing for " + task["species_prefix"])
            continue
        if require_genome and not has_required_source(actual, ("genome_path", "gbff_path")):
            staging_errors.append("Required genome input is missing for " + task["species_prefix"])
            continue
        actual_paths = [str(actual[k]) for k in ("cds_path", "gff_path", "gbff_path", "genome_path") if actual.get(k)]
        # Check the union before merging: a bound path serves both as original
        # input and actual output, and must match the frozen hash we publish.
        fence = reuse_fences.get(index)
        observed_hashes = (fence.verified([*actual_paths, *task.get("input_sha256", {})]) if fence
                           else digest_paths([*actual_paths, *task.get("input_sha256", {})]))
        for path, expected in task.get("input_sha256", {}).items():
            if observed_hashes[path] != expected:
                raise ValueError("Local input changed during staging: " + path)
        actual["input_sha256"] = {**task.get("input_sha256", {}),
                                  **{path: observed_hashes[path] for path in actual_paths}}
        for setting in ("gene_grouping_mode", "gff_repair_mode", "format_strict"):
            actual[setting] = task[setting]
        if fence:
            actual["staged_input_reuse"] = task["staged_input_reuse"]
        rows = [resolved_rows[key]]
        write_resolved_manifest_tsv(task_root / f"{index}.resolved.tsv", resolved_manifest_fieldnames(rows), rows)
        atomic_json(task_root / f"{index}.json", {
            "plan_sha256": plan_hash, "task_index": index,
            "task": {k: str(v) if isinstance(v, Path) else v for k, v in actual.items()},
        }, immutable=True)
        staged_count += 1

    reuse_reader.check()
    for fence in reuse_fences.values():
        fence.check()
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
    parser.add_argument("--require-gff", action="store_true")
    parser.add_argument("--require-genome", action="store_true")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    stage_downloads(args.task_plan, jobs=args.jobs, timeout=args.download_timeout,
                    require_gff=args.require_gff, require_genome=args.require_genome,
                    headers=fsi.parse_http_headers(args.http_header, args.auth_bearer_token_env))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Read-only inventory of scratch, quarantines and retained project caches."""

import argparse
import datetime as dt
import json
import os
import re
import stat
import time
from pathlib import Path


def identity(path):
    info = path.lstat()
    return {"device": info.st_dev, "inode": info.st_ino, "uid": info.st_uid,
            "mtime_ns": info.st_mtime_ns, "ctime_ns": info.st_ctime_ns}


def measure(path):
    result = {"files": 0, "bytes": 0, "directories": 0, "symlinks": 0,
              "foreign_entries": 0, "errors": [], "newest_mtime": 0}
    stack = [path]
    while stack:
        current = stack.pop()
        try:
            info = current.lstat()
            result["newest_mtime"] = max(result["newest_mtime"], info.st_mtime)
            result["foreign_entries"] += info.st_uid != os.getuid()
            if stat.S_ISLNK(info.st_mode):
                result["symlinks"] += 1
            elif stat.S_ISREG(info.st_mode):
                result["files"] += 1
                result["bytes"] += info.st_size
            elif stat.S_ISDIR(info.st_mode):
                result["directories"] += 1
                with os.scandir(current) as entries:
                    stack.extend(Path(entry.path) for entry in entries)
            else:
                result["errors"].append({"path": str(current), "error": "special entry"})
        except OSError as exc:
            result["errors"].append({"path": str(current), "error": str(exc)})
    result["age_days"] = round((time.time() - result.pop("newest_mtime")) / 86400, 2)
    return result


def inventory(workspace, *, legacy=False):
    workspace = Path(workspace).expanduser().resolve(strict=True)
    output = workspace if legacy else workspace / "output"
    records, errors, seen = [], [], set()

    def blocked_parent(path):
        for parent in path.parents:
            if parent == workspace:
                break
            if parent.is_symlink():
                errors.append({"path": str(path), "error": "symlinked parent: not traversed"})
                return True
        return False

    def add(path, kind, policy):
        if path in seen or blocked_parent(path) or not os.path.lexists(path):
            return
        seen.add(path)
        try:
            records.append({"path": str(path), "kind": kind, "policy": policy,
                            "identity": identity(path), **measure(path)})
        except OSError as exc:
            errors.append({"path": str(path), "error": str(exc)})

    add(output / "species_tree/tmp", "failed_scratch", "verify_idle_and_retired")
    for root in (output / "tmp", output / "transcriptome_assembly/tmp",
                 output / "query2family/tmp", output / "orthogroup/tmp"):
        if blocked_parent(root):
            continue
        if root.is_symlink():
            add(root, "symlink", "preserve")
        elif root.is_dir():
            try:
                for path in root.iterdir():
                    if re.match(r"^[0-9]+_", path.name):
                        add(path, "failed_scratch", "verify_idle_and_retired")
                    elif root == output / "tmp" and path.name == "kffractbias":
                        for task in path.iterdir() if not path.is_symlink() else []:
                            add(task, "failed_scratch", "verify_idle_and_retired")
            except OSError as exc:
                errors.append({"path": str(root), "error": str(exc)})
    add(output / "orthofinder/core", "core_results", "retain_unless_explicitly_retired")
    if not legacy:
        add(workspace / "downloads/tmp/species_genetic_code.resolved.tsv",
            "derived_genetic_code", "verify_idle_and_retained_inputs")
        add(output / "input_generation/tmp/input_download_cache", "provider_cache",
            "retain_for_plan_retry")

    # Directory discovery does not follow symlinks. Counts describe filename
    # entries and logical bytes, including hardlinked names and backup copies.
    stack = [workspace]
    while stack:
        directory = stack.pop()
        try:
            if directory.is_symlink():
                errors.append({"path": str(directory), "error": "symlink: not traversed"})
                continue
            with os.scandir(directory) as entries:
                for entry in entries:
                    path = Path(entry.path)
                    if entry.name == ".archive_cache":
                        add(path, "archive_cache", "retain_unless_explicitly_retired")
                        continue
                    if re.search(r"\.corrupt\.\d{14}\.\d+(?:\.\d+)?$", entry.name):
                        add(path, "download_quarantine", "verify_good_replacement_and_idle")
                    elif (entry.name == directory.name + ".query.fa"
                          and directory.parent == output / "genome_evolution/omark"):
                        add(path, "legacy_omamer_query", "verify_results_and_migrate_summary_provenance")
                    if entry.is_dir(follow_symlinks=False) and (
                        path not in seen or path == output / "input_generation/tmp/input_download_cache"
                    ):
                        stack.append(path)
        except OSError as exc:
            errors.append({"path": str(directory), "error": str(exc)})
    return {"schema_version": 1, "read_only": True, "workspace": str(workspace),
            "observed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
            "counts": "regular filename entries and logical bytes; not unique inodes or quota",
            "records": sorted(records, key=lambda row: row["path"]), "errors": errors}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--legacy-output-root", action="store_true",
                        help="Treat the supplied gfe_data path as the legacy output root")
    args = parser.parse_args()
    print(json.dumps(inventory(args.workspace, legacy=args.legacy_output_root), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

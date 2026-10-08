#!/usr/bin/env python3
"""List only explicitly declared local source files for contained native arrays."""

import argparse
import csv
import json
import re
from pathlib import Path
from urllib.parse import unquote, urlparse


def local_path(value, parent, workspace):
    parsed = urlparse(value)
    if parsed.scheme and parsed.scheme.lower() != "file":
        return None
    if parsed.scheme:
        if parsed.netloc not in ("", "localhost") or parsed.query or parsed.fragment:
            raise ValueError("Local source URL must name a file on this host")
        value = unquote(parsed.path)
    if value.startswith("/workspace/"):
        value = str(workspace / value[len("/workspace/"):])
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = parent / path
    path = path.resolve(strict=True)
    if any(char in str(path) for char in (":", ",")) or any(ord(char) < 32 or ord(char) == 127 for char in str(path)):
        raise ValueError("Local source path cannot be represented safely as a container bind")
    if not path.is_file():
        raise ValueError("Local source bind requires an explicitly declared regular file: " + str(path))
    return path


def source_files(source, workspace, plan=False):
    source, workspace = Path(source), Path(workspace)
    if plan:
        data = json.loads(source.read_text())
        rows = []
        for task in data["tasks"]:
            rows.append((task.get("manifest_row", {}), Path(task.get("manifest_parent", source.parent))))
            rows.append(({key: task.get(key) for key in ("cds_path", "gff_path", "gbff_path", "genome_path")}, source.parent))
            rows.extend(({"local_source_path": name}, source.parent) for name in task.get("input_sha256", {}))
    else:
        with source.open(newline="") as handle:
            first = handle.readline()
            handle.seek(0)
            rows = [(row, source.parent) for row in csv.DictReader(handle, delimiter="\t" if "\t" in first else ",")]
    paths = set()
    for row, parent in rows:
        if str(parent) == "/workspace" or str(parent).startswith("/workspace/"):
            parent = workspace / str(parent).removeprefix("/workspace").lstrip("/")
        for key, value in row.items():
            if not value or key not in {"cds_url", "gff_url", "gbff_url", "genome_url", "cds_path", "gff_path", "gbff_path", "genome_path", "local_source_path", "local_cds_path", "local_gff_path", "local_gbff_path", "local_genome_path"}:
                continue
            path = local_path(str(value), parent, workspace)
            if path is not None:
                paths.add(path)
    return sorted(paths)


def coalesce_project_sources(paths, project, workspace, destinations=()):
    """Restrict already-visible raw input trees without growing container scope.

    External sources and the writable workspace always retain exact file binds.
    Only native staged-cache generations or flat input directories are grouped;
    donor task plans and their lock namespaces must remain writable.
    """
    project, workspace = Path(project).resolve(strict=True), Path(workspace).resolve(strict=True)
    paths = set(paths)
    reserved = {Path(name).resolve() for name in destinations} - {project}

    def allowed(directory):
        return (directory != project and directory.is_relative_to(project)
                and not directory.is_relative_to(workspace) and not workspace.is_relative_to(directory)
                and not any(directory.is_relative_to(name) or name.is_relative_to(directory) for name in reserved))

    groups = set()
    for path in paths:
        for parent in path.parents:
            if (parent.parent.name == "staged" and parent.parent.parent.name == "input_download_cache"
                    and re.fullmatch(r"[0-9a-f]{64}", parent.name) and allowed(parent)):
                groups.add(parent)
                break
    remaining = {path for path in paths if not any(path.is_relative_to(group) for group in groups)}
    for directory in {path.parent for path in remaining}:
        if allowed(directory) and "input" in directory.relative_to(project).parts:
            members = {path for path in remaining if path.parent == directory}
            if len(members) > 1 and all(item.is_file() and not item.is_symlink() and not item.name.startswith(".")
                                        and not item.name.endswith(".lock") for item in directory.iterdir()):
                groups.add(directory)
    return sorted(groups | {path for path in remaining if not any(path.is_relative_to(group) for group in groups)})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--manifest", type=Path)
    group.add_argument("--plan", type=Path)
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--coalesce-project-sources", type=Path)
    parser.add_argument("--existing-bind-destination", action="append", default=[])
    args = parser.parse_args()
    try:
        paths = source_files(args.plan or args.manifest, args.workspace, plan=args.plan is not None)
        if args.coalesce_project_sources is not None:
            paths = coalesce_project_sources(paths, args.coalesce_project_sources, args.workspace,
                                             args.existing_bind_destination)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.error(str(exc))
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()

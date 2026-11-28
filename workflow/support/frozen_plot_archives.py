"""Bind archived dependencies of an existing frozen-summary provenance contract."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from dataclasses import asdict
from pathlib import Path, PurePosixPath

from archive_version_references import ArchiveReferenceError, ArchiveVersionReferences, MemberVersion
from gene_family_output_store import GeneFamilyOutputStore, read_only_observation, strict_json_loads


def logical_files(payload, manifest_name):
    paths = {"artifact_provenance/" + manifest_name: None}
    for collection in ("inputs", "outputs", "optional_outputs"):
        for entry in payload.get(collection, ()):
            if (entry.get("scope") != "logical" or entry.get("artifact_type", "file") != "file"
                    or entry.get("state") == "absent"):
                continue
            path = PurePosixPath(entry["path"])
            if path.is_absolute() or len(path.parts) != 2 or ".." in path.parts:
                raise ArchiveReferenceError("Invalid logical frozen plot dependency")
            digest = entry.get("sha256")
            if str(path) in paths and paths[str(path)] != digest:
                raise ArchiveReferenceError("Conflicting logical frozen plot dependency")
            paths[str(path)] = digest
    return paths


def audited_family_versions(store, manifest_name, payload):
    if payload.get("family_id") != store.family_filter or payload.get("step") != "summary_statistics":
        raise ArchiveReferenceError("Frozen plot archive selection does not match the recorded family")
    for path, digest in logical_files(payload, manifest_name).items():
        artifact = store.source_artifact(*path.split("/"))
        if artifact is None:
            raise ArchiveReferenceError("Frozen plot dependency is unavailable: " + path)
        if artifact.zip_path is None or artifact.sha256 is None:
            # Raw inputs and untracked archives still receive the normal full content audit.
            continue
        if artifact.member_name != path or artifact.family_id != store.family_filter:
            raise ArchiveReferenceError("Frozen plot archive member is outside the exact family: " + path)
        if digest is not None and artifact.sha256 != digest:
            raise ArchiveReferenceError("Frozen plot archive member changed: " + path)
        yield MemberVersion(path, artifact.generation, artifact.sha256)


def verify_reserved_family(store, consumer, manifest_name, payload, memo):
    references = ArchiveVersionReferences(store.archive_root)
    keys = logical_files(payload, manifest_name)
    versions = references.consumer_family_versions(consumer, store.family_filter)
    memo.guard(("frozen-archive-reference", consumer, store.family_filter), versions,
               lambda: references.consumer_family_versions(consumer, store.family_filter))
    for member in versions:
        path = member.logical_path
        if path not in keys:
            raise ArchiveReferenceError("Reserved plot dependency is outside the recorded manifest: " + path)
        subdir, name = path.split("/")
        artifact = store.source_artifact(subdir, name)
        if artifact is None:
            raise ArchiveReferenceError("Reserved plot dependency is unavailable: " + path)
        if artifact.zip_path is not None:
            if (artifact.member_name, artifact.generation, artifact.sha256) != (
                    path, member.generation, member.sha256):
                raise ArchiveReferenceError("Reserved plot dependency version changed: " + path)
        elif keys[path] is None:
            # A materialized provenance manifest has no enclosing manifest digest.
            with store.open_binary(subdir, name) as stream:
                if hashlib.file_digest(stream, "sha256").hexdigest() != member.sha256:
                    raise ArchiveReferenceError("Materialized reserved manifest changed: " + path)
        if keys[path] is not None and keys[path] != member.sha256:
            raise ArchiveReferenceError("Reserved plot dependency digest changed: " + path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--catalog", required=True, type=Path)
    parser.add_argument("--tasks", required=True)
    parser.add_argument("--consumer")
    parser.add_argument("--max-bytes", type=int)
    args = parser.parse_args()
    tasks = json.load(sys.stdin) if args.tasks == "-" else json.loads(args.tasks)
    if (not isinstance(tasks, list) or not tasks or any(type(i) is not int or i < 1 for i in tasks)
            or len(set(tasks)) != len(tasks) or bool(args.consumer) != (args.max_bytes is not None)):
        parser.error("provide unique one-based tasks and pair consumer with max bytes")
    wanted = set(tasks)
    with args.catalog.open() as source:
        families = [row["Orthogroup"] for index, row in enumerate(csv.DictReader(source, delimiter="\t"), 1)
                    if index in wanted]
    if len(families) != len(tasks) or len(set(families)) != len(families):
        raise ArchiveReferenceError("Frozen plot catalogue task identity changed")
    versions = []
    with read_only_observation():
        for family in families:
            store = GeneFamilyOutputStore(args.root, family_filter=family)
            name = family + ".summary_statistics.json"
            with store.open_binary("artifact_provenance", name) as stream:
                payload = strict_json_loads(stream.read())
            versions.extend(audited_family_versions(store, name, payload))
    if args.consumer:
        selection = GeneFamilyOutputStore(args.root).reserve_archive_versions(
            args.consumer, versions, max_bytes=args.max_bytes, verify_legacy=True)
        print(json.dumps({"consumer": args.consumer, "selection": selection, "count": len(versions)}))
    else:
        for version in versions:
            print(json.dumps(asdict(version), sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, ArchiveReferenceError) as exc:
        sys.exit(str(exc))

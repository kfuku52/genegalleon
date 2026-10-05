#!/usr/bin/env python3
"""Plan a BUSCO-filtered CDS/GFF/genome export without copying research inputs."""

import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
import shutil
import tempfile
from collections import defaultdict
from pathlib import Path

from busco_quality_metadata import parse_short_summary
from input_generation_array_state import digest_paths, safe_component


def table(path, raw=None):
    with (io.StringIO(raw.decode("utf-8")) if raw is not None else Path(path).open(encoding="utf-8", newline="")) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or any(not name or not name.isascii() for name in reader.fieldnames):
            raise ValueError("Metadata requires nonempty English column names")
        rows = list(reader)
        if any(None in row or any(value is None for value in row.values()) for row in rows):
            raise ValueError("Malformed metadata table: " + str(path))
        return reader.fieldnames, rows


def complete_groups(path, raw=None):
    groups = defaultdict(set)
    with (io.StringIO(raw.decode("utf-8")) if raw is not None else Path(path).open(encoding="utf-8")) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 2 or not fields[0] or fields[1] not in {"Complete", "Duplicated", "Fragmented", "Missing"}:
                raise ValueError("Invalid BUSCO full table: " + str(path))
            groups[fields[0]].add(fields[1])
    if not groups or any(len(statuses) != 1 for statuses in groups.values()):
        raise ValueError("Empty or contradictory BUSCO groups: " + str(path))
    return sum(next(iter(statuses)) in {"Complete", "Duplicated"} for statuses in groups.values()), len(groups)


def write_table(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def resolve_provenance_path(entry, workspace):
    if entry.get("scope") == "workspace":
        path = Path(entry["path"])
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("Unsafe workspace provenance path")
        return workspace / path
    if entry.get("scope") == "absolute" and Path(entry["path"]).is_absolute():
        return Path(entry["path"])
    raise ValueError("BUSCO provenance requires a workspace or absolute path")


def build(args):
    if (not math.isfinite(args.minimum_complete_busco) or not 0 <= args.minimum_complete_busco <= 100
            or not math.isfinite(args.target_part_gb) or args.target_part_gb <= 0):
        raise ValueError("BUSCO cutoff must be 0–100 and target part size must be positive")
    expected = {}
    source_map = {}
    identities = {}

    def remember(path):
        fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
        value = tuple(getattr(path.stat(), field) for field in fields)
        if str(path) in identities and identities[str(path)] != value:
            raise ValueError("Input changed while planning the export: " + str(path))
        identities[str(path)] = value

    def capture(path):
        path = Path(path).resolve(strict=True)
        remember(path)
        raw = path.read_bytes()
        remember(path)
        expected[str(path)] = hashlib.sha256(raw).hexdigest()
        return path, raw

    summary, raw = capture(args.species_summary)
    _fields, rows = table(summary, raw)
    seen, chosen, excluded = set(), [], []
    lineage_mode = set()
    for row in rows:
        species = row.get("species_prefix", "")
        if not safe_component(species) or species in seen:
            raise ValueError("Missing, unsafe, or duplicate species in summary: " + species)
        seen.add(species)
        short, short_raw = capture(args.busco_short_dir / f"{species}.busco.short.txt")
        full, _full_raw = capture(args.busco_full_dir / f"{species}.busco.full.tsv")
        provenance_path, provenance_raw = capture(args.busco_provenance_dir / f"busco.{species}.json")
        quality, provenance = parse_short_summary(short_raw.decode("utf-8"), short), json.loads(provenance_raw)
        complete, total = complete_groups(full, _full_raw)
        declared_total = re.search(r"(?:^|[,\s])n:(\d+)(?:\s|$)", short_raw.decode("utf-8"))
        if declared_total is None or int(declared_total[1]) != total:
            raise ValueError("BUSCO short/full group totals disagree for " + species)
        percentage = 100 * complete / total
        displayed = re.search(r"(?:^|[,\s])C:([0-9.]+)%\[", short_raw.decode("utf-8"))[1]
        precision = len(displayed.partition(".")[2])
        if abs(percentage - quality["busco_complete_pct"]) > 0.5 * 10 ** (-precision) + 1e-9:
            raise ValueError("BUSCO short/full complete counts disagree for " + species)
        lineage_mode.add((quality["lineage"], quality["mode"]))
        if not quality["lineage"] or quality["mode"] != "transcriptome":
            raise ValueError("Expected a named lineage and CDS transcriptome BUSCO for " + species)
        if (provenance.get("schema_version") != 1 or provenance.get("step") != "input_generation_species_busco"
                or provenance.get("family_id") != species
                or provenance.get("parameters", {}).get("busco_mode") != quality["mode"]
                or provenance.get("parameters", {}).get("busco_lineage_resolved") != quality["lineage"]):
            raise ValueError("BUSCO provenance contract differs for " + species)
        inputs, outputs = provenance.get("inputs", []), provenance.get("outputs", [])
        for entries, label, path in ((inputs, "species_cds", row["cds_output_path"]),
                                     (outputs, "busco_full", full), (outputs, "busco_short", short)):
            matches = [item for item in entries if item.get("label") == label]
            path = Path(path).resolve(strict=True)
            if len(matches) != 1 or resolve_provenance_path(matches[0], args.workspace_root).resolve(strict=True) != path:
                raise ValueError("BUSCO provenance source path differs for " + species)
            item = matches[0]
            remember(path)
            if str(path) in expected and expected[str(path)] != item["sha256"]:
                raise ValueError("BUSCO provenance content differs for " + species)
            expected[str(path)] = item["sha256"]
            if path.stat().st_size != item["size_bytes"]:
                raise ValueError("BUSCO provenance source size differs for " + species)
        if percentage < args.minimum_complete_busco:
            excluded.append({"species": species, "complete_percent": percentage})
            continue
        item = {"species": species, "taxid": row.get("taxid", ""), "busco_complete_percent": percentage,
                "busco_complete_groups": complete, "busco_total_groups": total,
                "busco_lineage": quality["lineage"], "busco_version": quality["busco_version"], "files": []}
        for role in ("cds", "gff", "genome"):
            source = Path(row[role + "_output_path"]).resolve(strict=True)
            remember(source)
            if (not safe_component(source.name) or not source.name.endswith(".gz")
                    or not source.is_file() or not source.stat().st_size):
                raise ValueError("Export requires a nonempty gzip triplet: " + str(source))
            relative = f"species/{species}/{source.name}"
            if relative in source_map or str(source) in source_map.values():
                raise ValueError("Duplicate export file or aliased source: " + relative)
            source_map[relative] = str(source)
            item["files"].append({"species": species, "role": role, "export_path": relative,
                                  "size_bytes": source.stat().st_size, "source_path": str(source)})
        item["triplet_size_bytes"] = sum(file["size_bytes"] for file in item["files"])
        chosen.append(item)
    if len(lineage_mode) != 1:
        raise ValueError("A cohort must use one BUSCO lineage and mode")
    if not chosen:
        raise ValueError("No species pass the requested BUSCO cutoff")
    metadata = {}
    for specification in args.metadata_table:
        name, separator, path = specification.partition("=")
        if not separator or not safe_component(name) or name in metadata or name in {"cohort", "files"}:
            raise ValueError("Metadata must be a unique safe NAME=PATH")
        source, metadata_raw = capture(path)
        fields, records = table(source, metadata_raw)
        keys = [key for key in ("species", "species_prefix") if key in fields]
        if len(keys) != 1:
            raise ValueError("Metadata needs exactly one species or species_prefix column")
        by_species = {}
        for record in records:
            species = record[keys[0]]
            if species in by_species:
                raise ValueError("Duplicate metadata species: " + species)
            by_species[species] = record
        if any(item["species"] not in by_species for item in chosen):
            raise ValueError("Metadata is missing selected species: " + name)
        metadata[name] = (fields, [by_species[item["species"]] for item in sorted(chosen, key=lambda x: x["species"])])
    # All selected triplets and selection evidence receive fresh hashes. Large
    # files are streamed; already parsed metadata must match those fresh reads.
    observed = digest_paths([*expected, *source_map.values()])
    if any(observed[path] != sha for path, sha in expected.items()):
        raise ValueError("Input or BUSCO evidence changed while planning the export")
    for path in identities:
        remember(Path(path))
    parts = []
    target = args.target_part_gb * 10 ** 9
    for item in sorted(chosen, key=lambda x: (-x["triplet_size_bytes"], x["species"])):
        suitable = [part for part in parts if part["gzip_bytes"] + item["triplet_size_bytes"] <= target]
        if suitable:
            part = min(suitable, key=lambda x: (x["gzip_bytes"], x["part"]))
        else:
            part = {"part": len(parts) + 1, "gzip_bytes": 0, "species": []}
            parts.append(part)
        item["part"] = part["part"]
        part["gzip_bytes"] += item["triplet_size_bytes"]
        part["species"].append(item["species"])
        for file in item["files"]:
            file.update(part=part["part"], sha256=observed[file["source_path"]])
    return chosen, excluded, parts, metadata, observed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-summary", type=Path, required=True)
    parser.add_argument("--workspace-root", type=Path, required=True, help="Resolve native BUSCO provenance paths.")
    parser.add_argument("--busco-short-dir", type=Path, required=True)
    parser.add_argument("--busco-full-dir", type=Path, required=True)
    parser.add_argument("--busco-provenance-dir", type=Path, required=True)
    parser.add_argument("--minimum-complete-busco", type=float, required=True)
    parser.add_argument("--target-part-gb", type=float, default=20)
    parser.add_argument("--metadata-table", action="append", default=[], help="Optional NAME=TSV with species labels.")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Export output must be a new directory")
    try:
        chosen, excluded, parts, metadata, evidence = build(args)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.error(str(exc))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".cohort-export-", dir=args.output.parent))
    try:
        files, cohort = [], []
        for item in sorted(chosen, key=lambda x: x["species"]):
            files.extend(item["files"])
            cohort.append({**{key: value for key, value in item.items() if key != "files"},
                           **{file["role"] + "_path": file["export_path"] for file in item["files"]}})
        write_table(temporary / "cohort.tsv", list(cohort[0]), cohort)
        shared_files = [{key: value for key, value in file.items() if key != "source_path"} for file in files]
        write_table(temporary / "files.tsv", list(shared_files[0]), shared_files)
        for name, (fields, rows) in metadata.items():
            write_table(temporary / f"{name}.tsv", fields, rows)
        payload = {"schema_version": 1, "minimum_complete_busco": args.minimum_complete_busco,
                   "selected_species": len(chosen), "triplet_file_count": len(files),
                   "gzip_bytes": sum(item["triplet_size_bytes"] for item in chosen),
                   "target_part_bytes": int(args.target_part_gb * 10 ** 9), "parts": parts,
                   "excluded": excluded, "source_files": files, "verified_source_sha256": evidence,
                   "status": "metadata_only; archives and uploads have not been created"}
        (temporary / "export_plan.json").write_text(json.dumps(payload, indent=2) + "\n")
        (temporary / "README.txt").write_text(
            "CDS/GFF/genome cohort export plan\n\n"
            "BUSCO filtering uses exact Complete + Duplicated group counts from CDS transcriptome BUSCO.\n"
            "cohort.tsv and files.tsv use recipient-relative paths. Metadata tables contain selected species only.\n"
            "export_plan.json contains internal source paths and fresh content hashes for future packaging.\n"
            "Part sizes count existing gzip payloads; tar headers/padding are additional.\n"
            "A species larger than the target receives its own part.\n"
            "This directory contains metadata only; input files were not copied, changed, or uploaded.\n")
        checksums = [hashlib.sha256(path.read_bytes()).hexdigest() + "  " + path.name
                     for path in sorted(temporary.iterdir())]
        (temporary / "SHA256SUMS").write_text("\n".join(checksums) + "\n")
        # Reserve the destination exclusively; directory rename can otherwise
        # replace an existing empty directory on POSIX systems.
        args.output.mkdir()
        for path in temporary.iterdir():
            os.rename(path, args.output / path.name)
        print(json.dumps({key: payload[key] for key in ("selected_species", "triplet_file_count", "gzip_bytes", "status")}))
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


if __name__ == "__main__":
    main()

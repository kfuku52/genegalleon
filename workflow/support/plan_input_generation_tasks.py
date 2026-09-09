#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path
from urllib.parse import unquote, urlparse

import format_species_inputs as fsi
from format_species_manifest import read_download_manifest
from format_species_provider_config import DOWNLOAD_MANIFEST_SUPPORTED_PROVIDERS
from format_species_taxonomy import invalid_species_key_error, normalize_species_key_for_runtime
from input_generation_array_state import atomic_json, digest, safe_component


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Discover gg_input_generation species tasks and write a reusable task-plan JSON."
    )
    parser.add_argument(
        "--provider",
        choices=("all",) + fsi.PROVIDERS,
        required=True,
        help="Input provider type. Use 'all' with --input-dir pointing to a provider-root directory.",
    )
    parser.add_argument("--download-manifest", default="")
    parser.add_argument("--download-dir", default="")
    parser.add_argument("--stage-downloads", action="store_true", help="Require prepare to stage manifest inputs before workers run.")
    parser.add_argument(
        "--input-dir",
        default="",
        help="Provider input directory to scan for species tasks.",
    )
    parser.add_argument(
        "--outfile",
        required=True,
        help="Output JSON path for the discovered task plan.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit with error if task discovery reports any species-level errors.",
    )
    parser.add_argument(
        "--gene-grouping-mode",
        choices=fsi.GENE_GROUPING_MODES,
        default="rescue_overlap",
        help="Gene grouping mode to embed in each planned task.",
    )
    parser.add_argument(
        "--gff-repair-mode",
        choices=fsi.GFF_REPAIR_MODES,
        default="safe",
        help="GFF repair mode to embed in each planned task.",
    )
    return parser


def serialize_task(task):
    payload = {}
    for key, value in task.items():
        if isinstance(value, Path):
            payload[key] = str(value)
        elif value is None:
            payload[key] = ""
        else:
            payload[key] = value
    return payload


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    try:
        provider_inputs = [] if args.download_manifest else fsi.resolve_provider_inputs(args)
    except ValueError as exc:
        parser.error(str(exc))

    all_tasks = []
    all_warnings = []
    all_errors = []
    resolved_inputs = []

    for provider, input_dir in provider_inputs:
        resolved_inputs.append({"provider": provider, "input_dir": str(input_dir)})
        if not input_dir.exists() or not input_dir.is_dir():
            message = "[{}] input directory not found: {}".format(provider, input_dir)
            if args.provider == "all":
                all_warnings.append(message)
            else:
                all_errors.append(message)
            continue
        tasks, warnings, errors = fsi.discover_tasks(provider, input_dir)
        for task in tasks:
            task["gene_grouping_mode"] = args.gene_grouping_mode
            task["gff_repair_mode"] = args.gff_repair_mode
            task["format_strict"] = bool(args.strict)
            task["input_sha256"] = {str(task[key]): digest(task[key]) for key in ("cds_path", "gff_path", "gbff_path", "genome_path") if task.get(key)}
        all_tasks.extend(tasks)
        all_warnings.extend(warnings)
        all_errors.extend(errors)

    if args.download_manifest:
        if not args.download_dir:
            parser.error("--download-dir is required with --download-manifest")
        manifest = Path(args.download_manifest).expanduser().resolve()
        for row in read_download_manifest(manifest):
            provider = row.get("provider", "").strip().lower()
            if args.provider != "all" and args.provider != provider:
                continue
            species = normalize_species_key_for_runtime(row.get("species_key", ""))
            if not safe_component(species) or any(row.get(key) and not safe_component(row[key]) for key in ("cds_filename", "gff_filename", "gbff_filename", "genome_filename")):
                parser.error("Species keys and download filenames must be safe, non-hidden filename components")
            if not row.get("id", "").strip() or provider not in DOWNLOAD_MANIFEST_SUPPORTED_PROVIDERS or not species or invalid_species_key_error(provider, species):
                parser.error("Array manifests require a supported provider and explicit valid species_key on every selected row")
            if provider == "local":
                from format_species_download.local import resolve_local_manifest_row
                row = {**row, **resolve_local_manifest_row(provider, row.get("id", ""), species, row,
                          manifest.parent, [], 0)}
            row = {**row, "provider": provider, "species_key": species}
            source_hashes = {}
            for key in ("cds_url", "gff_url", "gbff_url", "genome_url"):
                parsed = urlparse(row.get(key, ""))
                if parsed.scheme == "file":
                    source = str(Path(unquote(parsed.path)).resolve())
                    source_hashes[source] = digest(source)
            all_tasks.append({"input_sha256": source_hashes, "provider": provider, "species_key": species, "species_prefix": species,
                              "manifest_row": row, "manifest_parent": str(manifest.parent),
                              "download_dir": str(Path(args.download_dir).expanduser().resolve()),
                              "gene_grouping_mode": args.gene_grouping_mode,
                              "gff_repair_mode": args.gff_repair_mode, "format_strict": bool(args.strict)})
    species = [task["species_prefix"] for task in all_tasks]
    if any(not safe_component(name) for name in species):
        parser.error("Species prefixes must be safe, non-hidden filename components")
    if len(species) != len(set(species)):
        parser.error("Duplicate species prefixes would overwrite worker outputs")

    for warning in all_warnings:
        sys.stderr.write("Warning: {}\n".format(warning))
    for error in all_errors:
        sys.stderr.write("Error: {}\n".format(error))

    if args.strict and all_errors:
        return 1
    if not all_tasks:
        sys.stderr.write("No species tasks were discovered.\n")
        return 1

    outfile = Path(args.outfile).expanduser().resolve()
    outfile.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "version": 2,
        "provider": args.provider,
        "input_dir": str(Path(args.input_dir).expanduser().resolve()),
        "provider_inputs": resolved_inputs,
        "task_count": len(all_tasks),
        "species": [task["species_prefix"] for task in all_tasks],
        "tasks": [serialize_task(task) for task in all_tasks],
    }
    if args.stage_downloads:
        payload["download_mode"] = "staged"
    atomic_json(outfile, payload, immutable=True)

    print("Discovered {} species tasks -> {}".format(len(all_tasks), outfile))
    return 0


if __name__ == "__main__":
    sys.exit(main())

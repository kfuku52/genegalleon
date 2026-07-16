#!/usr/bin/env python3
"""Public API and CLI facade for GeneGalleon species input formatting."""
# ruff: noqa: F401,F403,F405

import sys
from pathlib import Path

SUPPORT_DIR = Path(__file__).resolve().parent
if str(SUPPORT_DIR) not in sys.path:
    sys.path.insert(0, str(SUPPORT_DIR))

import json

from format_species_annotations import *
from format_species_cli import *
from format_species_common import *
from format_species_constants import *
from format_species_discovery import *
from format_species_download_runtime import *
from format_species_download_stage import apply_download_input_dir, run_download_stage
from format_species_manifest import *
from format_species_provider_inputs import (
    manifest_declared_providers,
    manifest_declared_species_keys,
    resolve_provider_inputs,
)
from format_species_provider_resolvers import *
from format_species_provider_urls import *
from format_species_summary import *
from format_species_taxonomy import *
from format_species_writers import *


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    run_started_utc = utc_now_iso()

    if args.download_only and args.download_manifest == "":
        parser.error("--download-only requires --download-manifest.")
    if args.download_manifest != "" and args.provider in ("phycocosm", "phytozome"):
        parser.error(
            "--download-manifest does not support provider '{}'. Use --input-dir for local formatting.".format(
                args.provider
            )
        )

    try:
        http_headers = parse_http_headers(
            args.http_header,
            args.auth_bearer_token_env,
            auth_cookie_env=args.auth_cookie_env,
            jgi_login_env=args.jgi_login_env,
            jgi_password_env=args.jgi_password_env,
            timeout=args.download_timeout,
        )
    except ValueError as exc:
        parser.error(str(exc))

    download_report = None
    if args.download_manifest != "":
        download_report = run_download_stage(
            args,
            http_headers,
            download_from_manifest=download_from_manifest,
            format_download_diagnostics_line=format_download_diagnostics_line,
            resolve_parallel_jobs=resolve_parallel_jobs,
            stderr=sys.stderr,
        )
        if args.download_only:
            return 0 if len(download_report["errors"]) == 0 else 1
        if len(download_report["errors"]) > 0 and args.strict:
            return 1

    apply_download_input_dir(args)

    manifest_rows_for_provider_inputs = None
    if download_report is not None:
        manifest_rows_for_provider_inputs = list(download_report.get("resolved_rows", ()))

    try:
        provider_inputs = resolve_provider_inputs(args, manifest_rows=manifest_rows_for_provider_inputs)
    except ValueError as exc:
        parser.error(str(exc))

    output_cds_dir = Path(args.species_cds_dir).expanduser().resolve()
    output_gff_dir = Path(args.species_gff_dir).expanduser().resolve()
    output_genome_dir = Path(args.species_genome_dir).expanduser().resolve()
    species_summary_path = Path(args.species_summary_output).expanduser().resolve()
    output_cds_dir.mkdir(parents=True, exist_ok=True)
    output_gff_dir.mkdir(parents=True, exist_ok=True)
    output_genome_dir.mkdir(parents=True, exist_ok=True)
    species_summary_rows = retain_existing_species_summary_rows(read_species_summary_rows(species_summary_path))
    write_species_summary_rows(species_summary_path, species_summary_rows)
    taxonomy_resolver = SpeciesTaxonomyMetadataResolver.from_environment()

    all_tasks = []
    all_warnings = []
    all_errors = []

    for provider, input_dir in provider_inputs:
        if not input_dir.exists() or not input_dir.is_dir():
            message = "[{}] input directory not found: {}".format(provider, input_dir)
            if args.provider == "all":
                all_warnings.append(message)
            else:
                all_errors.append(message)
            continue
        allowed_species_keys = None
        if manifest_rows_for_provider_inputs is not None:
            allowed_species_keys = manifest_declared_species_keys(manifest_rows_for_provider_inputs, provider)
        tasks, warnings, errors = discover_tasks(
            provider,
            input_dir,
            allowed_species_keys=allowed_species_keys,
        )
        for task in tasks:
            task["gene_grouping_mode"] = args.gene_grouping_mode
        all_tasks.extend(tasks)
        all_warnings.extend(warnings)
        all_errors.extend(errors)

    for warning in all_warnings:
        sys.stderr.write("Warning: {}\n".format(warning))
    for error in all_errors:
        sys.stderr.write("Error: {}\n".format(error))

    if args.strict and len(all_errors) > 0:
        return 1
    if len(all_tasks) == 0:
        sys.stderr.write("No species inputs with CDS, GBFF, or GFF plus genome were discovered.\n")
        return 1

    processed = 0
    total_duplicates = 0
    total_cds_before = 0
    total_cds_after = 0
    first_cds_sequence_name = ""
    species_with_genome = 0
    failed_format_tasks = 0
    for task in all_tasks:
        cds_result = format_cds(task, output_cds_dir, args.overwrite, args.dry_run)
        gff_result = format_gff(task, output_gff_dir, args.overwrite, args.dry_run)
        genome_result = format_genome(task, output_genome_dir, args.overwrite, args.dry_run)
        stale_gff_outputs = []
        if gff_result["status"] in ("write", "skip"):
            stale_gff_outputs = remove_stale_ensembl_like_partial_gff_outputs(
                task["provider"],
                task["species_prefix"],
                output_gff_dir,
                gff_result.get("output_path"),
            )
            for stale_name in stale_gff_outputs:
                sys.stderr.write(
                    "Warning: [{}] {}: removed stale partial GFF output '{}'.\n".format(
                        task["provider"],
                        task["species_prefix"],
                        stale_name,
                    )
                )
        if cds_result["status"] == "empty":
            failed_format_tasks += 1
            sys.stderr.write(
                "Warning: [{}] {}: derived CDS contained no records; no CDS output was written.\n".format(
                    task["provider"],
                    task["species_prefix"],
                )
            )
        if gff_result["status"] == "empty":
            failed_format_tasks += 1
            sys.stderr.write(
                "Warning: [{}] {}: derived GFF contained no feature rows; no GFF output was written.\n".format(
                    task["provider"],
                    task["species_prefix"],
                )
            )
        if format_task_succeeded(cds_result, gff_result, genome_result, args.dry_run):
            key = species_row_key(task["provider"], task["species_key"], task["species_prefix"])
            taxonomy_metadata = taxonomy_resolver.resolve(task["species_prefix"])
            species_summary_rows[key] = build_species_summary_row(
                task,
                cds_result,
                gff_result,
                genome_result,
                run_started_utc=run_started_utc,
                overwrite=args.overwrite,
                dry_run=args.dry_run,
                taxonomy_metadata=taxonomy_metadata,
            )
            # Persist each successful species incrementally.
            write_species_summary_rows(species_summary_path, species_summary_rows)
        processed += 1
        total_duplicates += cds_result["duplicates"]
        total_cds_before += cds_result["before_count"]
        total_cds_after += cds_result["after_count"]
        if first_cds_sequence_name == "":
            first_cds_sequence_name = cds_result["first_sequence_name"]
        if genome_result["status"] != "missing":
            species_with_genome += 1
        print(
            "[{}] {}: CDS={} ({}, {}, aggregated_away={}, before={}, after={}), GFF={} ({}, lines={}), GENOME={} ({})".format(
                task["provider"],
                task["species_prefix"],
                task["cds_path"].name
                if task.get("cds_path") is not None
                else Path(str(describe_task_cds_input(task) or "derived_cds")).name,
                cds_result["status"],
                result_output_name(cds_result),
                cds_result["duplicates"],
                cds_result["before_count"],
                cds_result["after_count"],
                Path(str(describe_task_gff_input(task) or "derived_gff")).name,
                gff_result["status"],
                gff_result["lines"],
                Path(str(describe_task_genome_input(task) or "NA")).name,
                genome_result["status"],
            )
        )

    stats = {
        "species_processed": processed,
        "num_species_cds_files": processed,
        "num_species_gff_files": processed,
        "num_species_genome_files": species_with_genome,
        "cds_sequences_before": total_cds_before,
        "cds_sequences_after": total_cds_after,
        "cds_first_sequence_name": first_cds_sequence_name,
        "aggregated_cds_removed": total_duplicates,
        "duplicate_cds_ids_skipped": total_duplicates,
        "dry_run": int(args.dry_run),
    }
    if args.stats_output != "":
        stats_path = Path(args.stats_output).expanduser().resolve()
        stats_path.parent.mkdir(parents=True, exist_ok=True)
        with open(stats_path, "wt", encoding="utf-8") as handle:
            json.dump(stats, handle, ensure_ascii=True, indent=2, sort_keys=True)

    print(
        "Finished. species processed: {}, CDS aggregated away: {}, CDS before/after: {}/{}, first CDS sequence: {}, output CDS dir: {}, output GFF dir: {}, output genome dir: {}, species summary: {}, dry_run={}".format(
            processed,
            total_duplicates,
            total_cds_before,
            total_cds_after,
            first_cds_sequence_name if first_cds_sequence_name != "" else "NA",
            output_cds_dir,
            output_gff_dir,
            output_genome_dir,
            species_summary_path,
            int(args.dry_run),
        )
    )
    if len(all_errors) > 0:
        return 2
    if failed_format_tasks > 0:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())

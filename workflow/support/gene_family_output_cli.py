"""Argument parser for the gene-family output-store command line."""

from __future__ import annotations

import argparse
from pathlib import Path


def build_parser(
    *,
    description: str,
    default_large_zip_warning_bytes: int,
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=description)
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add_root(subparser: argparse.ArgumentParser) -> None:
        subparser.add_argument("--root", required=True, type=Path)

    def add_family_context(subparser: argparse.ArgumentParser) -> None:
        subparser.add_argument(
            "--mode",
            choices=["query2family", "orthogroup"],
            help="Omit when an existing ZIP store or one unambiguous catalog argument identifies the mode.",
        )
        subparser.add_argument("--query-dir", type=Path)
        subparser.add_argument("--genecount", type=Path)

    def add_zip_write_options(subparser: argparse.ArgumentParser) -> None:
        subparser.add_argument(
            "--compression",
            choices=["adaptive", "deflate", "store"],
            default="adaptive",
        )
        subparser.add_argument("--compression-level", type=int, default=6)
        subparser.add_argument("--workers", type=int, default=1)

    archive_parser = subparsers.add_parser("archive-completed")
    add_root(archive_parser)
    add_family_context(archive_parser)
    archive_parser.add_argument("--min-files", type=int, default=1)
    archive_parser.add_argument("--max-files-per-shard", type=int, default=5000)
    archive_parser.add_argument("--nonblocking", action="store_true")
    archive_parser.add_argument("--progress-interval", type=float, default=30.0)
    archive_parser.add_argument("--max-final-zip-bytes", type=int, default=0)
    add_zip_write_options(archive_parser)

    archive_family_parser = subparsers.add_parser(
        "archive-family",
        help="Archive live artifacts for one family even when it is incomplete.",
    )
    add_root(archive_family_parser)
    add_family_context(archive_family_parser)
    archive_family_parser.add_argument("--family-id", required=True)
    archive_family_parser.add_argument("--max-files-per-shard", type=int, default=5000)
    archive_family_parser.add_argument("--nonblocking", action="store_true")
    archive_family_parser.add_argument("--max-final-zip-bytes", type=int, default=0)
    add_zip_write_options(archive_family_parser)

    convert_parser = subparsers.add_parser(
        "convert-storage",
        help="Convert a gene-family output root between raw files and ZIP storage.",
    )
    add_root(convert_parser)
    add_family_context(convert_parser)
    convert_parser.add_argument("--to", required=True, choices=["zip", "raw", "files"])
    convert_parser.add_argument("--family-id-file", type=Path)
    convert_parser.add_argument("--max-files-per-shard", type=int, default=5000)
    convert_parser.add_argument("--strict-unmatched", action="store_true")
    convert_parser.add_argument("--pure-raw", action="store_true")
    convert_parser.add_argument("--dry-run", action="store_true")
    convert_parser.add_argument("--json", action="store_true")
    convert_parser.add_argument(
        "--available-bytes",
        type=int,
        help="Quota-aware bytes available for temporary ZIP creation or raw materialization; combined with filesystem free space.",
    )
    convert_parser.add_argument(
        "--large-zip-warning-bytes",
        type=int,
        default=default_large_zip_warning_bytes,
        help="Warn when a logical subdirectory exceeds this many bytes; 0 disables warnings.",
    )
    convert_parser.add_argument(
        "--max-final-zip-bytes",
        type=int,
        default=0,
        help="Keep human-readable archives/<subdir>/<subdir>.part-N.zip files above this logical size; 0 allows one final ZIP of any size.",
    )
    convert_parser.add_argument(
        "--resume",
        action="store_true",
        help="Require and resume a matching interrupted conversion.",
    )
    convert_parser.add_argument("--progress-interval", type=float, default=30.0)
    add_zip_write_options(convert_parser)

    conversion_status_parser = subparsers.add_parser(
        "conversion-status",
        help="Report resumable storage-conversion state without requiring a catalog.",
    )
    add_root(conversion_status_parser)
    conversion_status_parser.add_argument("--json", action="store_true")

    optimize_metadata_parser = subparsers.add_parser(
        "optimize-metadata",
        help="Prune obsolete pre-stripe lock files after all project jobs are stopped.",
    )
    add_root(optimize_metadata_parser)

    storage_status_parser = subparsers.add_parser(
        "storage-status",
        help="Report logical, raw, ZIP, and conversion counts.",
    )
    add_root(storage_status_parser)
    add_family_context(storage_status_parser)
    storage_status_parser.add_argument("--family-id-file", type=Path)
    storage_status_parser.add_argument("--json", action="store_true")

    materialize_parser = subparsers.add_parser("materialize-family")
    add_root(materialize_parser)
    add_family_context(materialize_parser)
    materialize_parser.add_argument("--family-id", required=True)
    materialize_parser.add_argument("--destination-root", "--destination", dest="destination_root", type=Path)
    materialize_parser.add_argument("--subdirs")
    materialize_parser.add_argument("--receipt", type=Path)
    materialize_parser.add_argument("--run-token", default="")

    materialize_many_parser = subparsers.add_parser("materialize-families")
    add_root(materialize_many_parser)
    add_family_context(materialize_many_parser)
    materialize_many_parser.add_argument("--family-id-file", required=True, type=Path)
    materialize_many_parser.add_argument(
        "--destination-root",
        "--destination",
        dest="destination_root",
        required=True,
        type=Path,
    )
    materialize_many_parser.add_argument("--subdirs")

    export_parser = subparsers.add_parser("export-current")
    add_root(export_parser)
    export_parser.add_argument("--destination-root", required=True, type=Path)

    complete_parser = subparsers.add_parser("is-complete")
    add_root(complete_parser)
    complete_parser.add_argument("--family-id", required=True)

    delete_parser = subparsers.add_parser("delete")
    add_root(delete_parser)
    delete_parser.add_argument("--path", required=True)
    delete_parser.add_argument("--keep-live", action="store_true")
    delete_parser.add_argument("--family-id")

    undelete_parser = subparsers.add_parser("undelete")
    add_root(undelete_parser)
    undelete_parser.add_argument("--path", required=True)
    undelete_parser.add_argument("--family-id")

    restore_parser = subparsers.add_parser("restore")
    add_root(restore_parser)
    restore_parser.add_argument("--path", required=True)
    restore_parser.add_argument("--overwrite", action="store_true")
    restore_parser.add_argument("--family-id")

    verify_parser = subparsers.add_parser("verify")
    add_root(verify_parser)
    verify_parser.add_argument("--progress-interval", type=float, default=30.0)
    verify_parser.add_argument("--json", action="store_true")
    verify_depth = verify_parser.add_mutually_exclusive_group()
    verify_depth.add_argument(
        "--quick",
        action="store_true",
        help="Validate ZIP inventories, manifests, and indexes without reading every payload byte.",
    )
    verify_depth.add_argument(
        "--deep",
        action="store_true",
        help="Read every member and validate CRC and SHA256 (the default).",
    )

    list_parser = subparsers.add_parser("list")
    add_root(list_parser)
    list_parser.add_argument("--subdir")

    has_files_parser = subparsers.add_parser("has-files")
    add_root(has_files_parser)
    has_files_parser.add_argument("--subdir", required=True)
    has_files_parser.add_argument("--suffix", default="")

    lock_path_parser = subparsers.add_parser("lock-path")
    add_root(lock_path_parser)
    lock_path_parser.add_argument("--family-id", required=True)

    status_parser = subparsers.add_parser("status")
    add_root(status_parser)
    add_family_context(status_parser)

    refresh_status_parser = subparsers.add_parser(
        "refresh-status",
        help="Refresh ARCHIVE_STATUS.tsv after manual live-file changes.",
    )
    add_root(refresh_status_parser)

    compact_parser = subparsers.add_parser("compact")
    add_root(compact_parser)
    compact_parser.add_argument("--mode", choices=["query2family", "orthogroup"])
    compact_parser.add_argument("--max-files-per-shard", type=int, default=5000)
    compact_parser.add_argument("--max-final-zip-bytes", type=int, default=0)
    compact_parser.add_argument("--nonblocking", action="store_true")
    add_zip_write_options(compact_parser)

    finalize_parser = subparsers.add_parser(
        "finalize",
        help="Consolidate each selected logical subdirectory into <subdirectory>.zip.",
    )
    add_root(finalize_parser)
    add_family_context(finalize_parser)
    finalize_parser.add_argument("--subdirs")
    finalize_parser.add_argument("--nonblocking", action="store_true")
    finalize_parser.add_argument("--progress-interval", type=float, default=30.0)
    finalize_parser.add_argument("--max-final-zip-bytes", type=int, default=0)
    add_zip_write_options(finalize_parser)

    migrate_layout_parser = subparsers.add_parser(
        "migrate-layout",
        help="Move a legacy .gg_archives store into the visible ZIP layout.",
    )
    add_root(migrate_layout_parser)

    repair_parser = subparsers.add_parser("repair")
    add_root(repair_parser)
    repair_parser.add_argument("--remove-orphans", action="store_true")
    repair_parser.add_argument("--progress-interval", type=float, default=30.0)

    purge_parser = subparsers.add_parser("purge")
    add_root(purge_parser)
    add_family_context(purge_parser)
    purge_parser.add_argument("--max-files-per-shard", type=int, default=5000)
    purge_parser.add_argument("--drop-unlisted", action="store_true")
    add_zip_write_options(purge_parser)

    cleanup_tmp_parser = subparsers.add_parser("cleanup-tmp")
    add_root(cleanup_tmp_parser)
    cleanup_tmp_parser.add_argument("--older-than-days", type=float, required=True)
    cleanup_tmp_parser.add_argument("--max-directories", type=int, default=0)
    cleanup_tmp_parser.add_argument("--max-bytes", type=int, default=0)
    cleanup_tmp_parser.add_argument("--max-files", type=int, default=0)
    cleanup_tmp_parser.add_argument("--nonblocking", action="store_true")

    cleanup_materialized_parser = subparsers.add_parser("cleanup-materialized")
    cleanup_materialized_parser.add_argument("--receipt", required=True, type=Path)
    cleanup_materialized_parser.add_argument("--nonblocking", action="store_true")

    for command, status in (
        ("mark-running", "running"),
        ("mark-complete", "complete"),
        ("mark-failed", "failed"),
    ):
        state_parser = subparsers.add_parser(command)
        add_root(state_parser)
        state_parser.add_argument("--family-id", required=True)
        state_parser.add_argument("--run-token", default="")
        state_parser.set_defaults(family_status=status)
    return parser

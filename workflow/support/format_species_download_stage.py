from pathlib import Path

from format_species_provider_config import DEFAULT_INPUT_RELATIVE_DIRS


def run_download_stage(
    args,
    http_headers,
    *,
    download_from_manifest,
    format_download_diagnostics_line,
    resolve_parallel_jobs,
    stderr,
):
    parallel_jobs = resolve_parallel_jobs(args.jobs)
    download_root = Path(args.download_dir).expanduser().resolve()
    download_root.mkdir(parents=True, exist_ok=True)
    resolved_manifest_output_path = None
    if args.resolved_manifest_output != "":
        resolved_manifest_output_path = Path(args.resolved_manifest_output).expanduser().resolve()

    report = download_from_manifest(
        manifest_path=Path(args.download_manifest).expanduser().resolve(),
        download_root=download_root,
        provider_filter=args.provider,
        overwrite=args.overwrite,
        headers=http_headers,
        timeout=float(args.download_timeout),
        dry_run=args.dry_run,
        jobs=parallel_jobs,
        resolved_manifest_output_path=resolved_manifest_output_path,
    )
    for warning in report["warnings"]:
        stderr.write("Warning: {}\n".format(warning))
    for error in report["errors"]:
        stderr.write("Error: {}\n".format(error))
    message_template = (
        "Download stage: rows processed={}, files downloaded={}, planned downloads={}, "
        "download root={}, resolved manifest={}, dry_run={}, jobs={}"
    )
    print(
        message_template.format(
            report["processed"],
            report["downloaded"],
            report["planned"],
            download_root,
            report.get("resolved_manifest_output", ""),
            int(args.dry_run),
            parallel_jobs,
        )
    )
    print(format_download_diagnostics_line(report.get("download_diagnostics", {})))
    return report


def apply_download_input_dir(args):
    if args.input_dir != "" or args.download_manifest == "":
        return
    download_input_root = Path(args.download_dir).expanduser().resolve()
    if args.provider == "all":
        args.input_dir = str(download_input_root)
        return
    args.input_dir = str((download_input_root / DEFAULT_INPUT_RELATIVE_DIRS[args.provider]).resolve())

"""Download runtime implementation: manifest."""

import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import urlparse

from format_species_constants import (
    FERNBASE_CONFIDENCE_MODE_FIELD,
    FERNBASE_CONFIDENCE_MODE_HIGH_LOW_COMBINED,
)
from format_species_manifest import (
    build_resolved_manifest_row,
    manifest_has_usable_source_bundle,
    read_download_manifest,
    resolved_manifest_fieldnames,
    write_resolved_manifest_tsv,
)
from format_species_provider_config import (
    DOWNLOAD_MANIFEST_SUPPORTED_PROVIDERS,
    ENSEMBL_LIKE_PROVIDERS,
    ORYZA_MINUTA_PROVIDER,
    PROVIDERS,
)
from format_species_provider_resolvers import (
    build_fernbase_combined_filename,
    effective_fernbase_confidence_mode,
    extract_coge_gid_candidate,
    merge_fernbase_confidence_bundle,
    merge_oryza_minuta_multisource_bundle,
    normalize_manifest_source_id,
    parse_fernbase_confidence_mode,
    provider_raw_dir,
    resolve_fernbase_split_bundle,
    resolve_ncbi_download_urls_from_id,
    resolve_non_ncbi_download_urls_from_id,
    resolve_oryza_minuta_download_bundle_from_id,
    resolve_preferred_ensembl_like_gff_url,
    validate_coge_export_gff_file,
)
from format_species_provider_urls import (
    extract_ncbi_accession_from_source_id,
)
from format_species_taxonomy import invalid_species_key_error, normalize_species_key_for_runtime

from .local import (
    quarantine_corrupt_gzip,
    resolve_local_manifest_row,
)
from .locking import (
    download_url_to_file,
    resolve_download_lock_stale_seconds,
    scan_download_cache_diagnostics,
    summarize_download_diagnostics,
)
from .targets import (
    default_download_filename,
    download_ncbi_datasets_file_from_id,
    resolve_download_urls_from_templates,
    resolve_provider_download_limits,
)


def execute_download_target_job(
    job,
    headers,
    timeout,
    overwrite,
    lock_stale_seconds,
    provider_semaphores,
):
    provider = job["provider"]
    source_id = job["source_id"]
    species_key = job["species_key"]
    label = job["label"]
    url = job["url"]
    target = job["target"]
    archive_member = str(job.get("archive_member") or "").strip()
    local_warnings = []
    local_errors = []
    downloaded = 0
    failed = []

    sem = provider_semaphores.get(provider)
    if sem is None:
        sem = threading.Semaphore(1)

    with sem:
        if target.exists() and target.stat().st_size > 0 and not overwrite:
            if quarantine_corrupt_gzip(
                target,
                local_warnings,
                "[download:{}] {} {}".format(provider, species_key, label),
            ):
                pass
            else:
                local_warnings.append(
                    "[download:{}] {} {} already exists. Skipping: {}".format(provider, species_key, label, target)
                )
                return {"warnings": local_warnings, "errors": local_errors, "downloaded": downloaded, "failed": failed}

        try:
            did_download = download_url_to_file(
                url,
                target,
                headers=headers,
                timeout=timeout,
                dry_run=False,
                overwrite=overwrite,
                lock_stale_seconds=lock_stale_seconds,
                warnings=local_warnings,
                lock_context="[download:{}] {} {}".format(provider, species_key, label),
                archive_member=archive_member,
            )
            if did_download:
                downloaded += 1
            elif target.exists() and target.stat().st_size > 0 and not overwrite:
                local_warnings.append(
                    "[download:{}] {} {} already exists after lock. Skipping: {}".format(
                        provider, species_key, label, target
                    )
                )
        except Exception as exc:
            fallback_exc = None
            if provider in ("ncbi", "refseq", "genbank") and source_id != "":
                try:
                    did_download = download_ncbi_datasets_file_from_id(
                        source_id=source_id,
                        label=label,
                        destination=target,
                        headers=headers,
                        timeout=timeout,
                        dry_run=False,
                        overwrite=overwrite,
                        lock_stale_seconds=lock_stale_seconds,
                        warnings=local_warnings,
                        lock_context="[download:{}] {} {} datasets".format(provider, species_key, label),
                    )
                    if did_download:
                        downloaded += 1
                        local_warnings.append(
                            "[download:{}] {} {} fallback via NCBI Datasets API for id '{}'".format(
                                provider, species_key, label, source_id
                            )
                        )
                        return {
                            "warnings": local_warnings,
                            "errors": local_errors,
                            "downloaded": downloaded,
                            "failed": failed,
                        }
                except Exception as fallback_error:
                    fallback_exc = fallback_error
            if fallback_exc is None:
                failed.append(
                    {
                        "row_id": job.get("row_id"),
                        "message": "[download:{}] failed {} {} from {} -> {} ({})".format(
                            provider, species_key, label, url, target, exc
                        ),
                    }
                )
            else:
                failed.append(
                    {
                        "row_id": job.get("row_id"),
                        "message": "[download:{}] failed {} {} from {} -> {} ({}) ; fallback datasets failed ({})".format(
                            provider, species_key, label, url, target, exc, fallback_exc
                        ),
                    }
                )
    return {"warnings": local_warnings, "errors": local_errors, "downloaded": downloaded, "failed": failed}


def download_from_manifest(
    manifest_path,
    download_root,
    provider_filter,
    overwrite,
    headers,
    timeout,
    dry_run,
    jobs,
    resolved_manifest_output_path=None,
):
    rows = read_download_manifest(manifest_path)
    warnings = []
    errors = []
    processed = 0
    downloaded = 0
    planned = 0
    resolved_rows = []
    resolved_fieldnames = resolved_manifest_fieldnames(rows)
    lock_stale_seconds = resolve_download_lock_stale_seconds()
    download_jobs = []
    merge_jobs = []
    cleanup_paths = []
    failed_downloads = []
    row_target_paths = {}
    preexisting_cache_diagnostics = scan_download_cache_diagnostics(download_root)

    if len(rows) == 0:
        errors.append("Download manifest is empty: {}".format(manifest_path))
        final_cache_diagnostics = scan_download_cache_diagnostics(download_root)
        return {
            "warnings": warnings,
            "errors": errors,
            "processed": processed,
            "downloaded": downloaded,
            "planned": planned,
            "resolved_rows": resolved_rows,
            "download_diagnostics": summarize_download_diagnostics(
                preexisting_cache_diagnostics,
                final_cache_diagnostics,
                warnings,
                len(download_jobs),
                len(failed_downloads),
            ),
        }
    manifest_parent_dir = manifest_path.parent
    header_cols = set(rows[0].keys())
    if "provider" not in header_cols or "id" not in header_cols:
        errors.append("Download manifest must contain required columns provider,id: {}".format(manifest_path))
        final_cache_diagnostics = scan_download_cache_diagnostics(download_root)
        return {
            "warnings": warnings,
            "errors": errors,
            "processed": processed,
            "downloaded": downloaded,
            "planned": planned,
            "resolved_rows": resolved_rows,
            "download_diagnostics": summarize_download_diagnostics(
                preexisting_cache_diagnostics,
                final_cache_diagnostics,
                warnings,
                len(download_jobs),
                len(failed_downloads),
            ),
        }

    for i, row in enumerate(rows, start=2):
        provider = (row.get("provider") or "").strip().lower()
        source_id_raw = (row.get("id") or "").strip()
        source_id = normalize_manifest_source_id(provider, source_id_raw)
        species_key = (row.get("species_key") or "").strip()
        cds_url = (row.get("cds_url") or "").strip()
        gff_url = (row.get("gff_url") or "").strip()
        gbff_url = (row.get("gbff_url") or "").strip()
        genome_url = (row.get("genome_url") or "").strip()
        cds_archive_member = (row.get("cds_archive_member") or "").strip()
        gff_archive_member = (row.get("gff_archive_member") or "").strip()
        gbff_archive_member = (row.get("gbff_archive_member") or "").strip()
        genome_archive_member = (row.get("genome_archive_member") or "").strip()
        cds_filename = (row.get("cds_filename") or "").strip()
        gff_filename = (row.get("gff_filename") or "").strip()
        gbff_filename = (row.get("gbff_filename") or "").strip()
        genome_filename = (row.get("genome_filename") or "").strip()
        gff_validation = ""
        fernbase_confidence_mode_raw = (row.get(FERNBASE_CONFIDENCE_MODE_FIELD) or "").strip()
        resolved_ncbi = None
        oryza_minuta_bundle = None

        if provider_filter != "all" and provider != provider_filter:
            continue
        processed += 1

        if provider == "":
            errors.append("Manifest line {}: provider is required".format(i))
            continue
        if source_id == "":
            errors.append("Manifest line {}: id is required".format(i))
            continue
        if provider not in PROVIDERS:
            errors.append("Manifest line {}: unsupported provider '{}'".format(i, provider))
            continue
        try:
            fernbase_confidence_mode = parse_fernbase_confidence_mode(fernbase_confidence_mode_raw)
        except ValueError as exc:
            errors.append("Manifest line {}: {}".format(i, exc))
            continue
        if provider != "fernbase" and fernbase_confidence_mode != "":
            warnings.append(
                "Manifest line {}: {} is ignored for provider '{}'".format(
                    i,
                    FERNBASE_CONFIDENCE_MODE_FIELD,
                    provider,
                )
            )
        if provider == "coge":
            coge_gid = extract_coge_gid_candidate(source_id)
            if coge_gid == "":
                errors.append("Manifest line {}: provider=coge requires numeric genome_id (gid) in id column".format(i))
                continue
            source_id = coge_gid
        if provider not in DOWNLOAD_MANIFEST_SUPPORTED_PROVIDERS:
            errors.append(
                "Manifest line {}: provider '{}' is not supported for --download-manifest "
                "(use --input-dir for local formatting)".format(i, provider)
            )
            continue

        if provider == ORYZA_MINUTA_PROVIDER:
            try:
                oryza_minuta_bundle = resolve_oryza_minuta_download_bundle_from_id(source_id, species_key)
                source_id = str(oryza_minuta_bundle.get("canonical_id") or source_id).strip()
                species_key = str(oryza_minuta_bundle.get("species_key") or species_key).strip()
                if cds_filename == "":
                    cds_filename = str(oryza_minuta_bundle.get("cds_filename") or "").strip()
                if gff_filename == "":
                    gff_filename = str(oryza_minuta_bundle.get("gff_filename") or "").strip()
                if genome_filename == "":
                    genome_filename = str(oryza_minuta_bundle.get("genome_filename") or "").strip()
            except Exception as exc:
                errors.append(
                    "Manifest line {}: failed to resolve id '{}' (provider={}): {}".format(i, source_id, provider, exc)
                )
                continue
        elif provider == "local":
            try:
                resolved_local = resolve_local_manifest_row(
                    provider=provider,
                    source_id=source_id,
                    species_key=species_key,
                    row=row,
                    manifest_parent_dir=manifest_parent_dir,
                    warnings=warnings,
                    line_number=i,
                )
                species_key = str(resolved_local.get("species_key") or "").strip()
                cds_url = str(resolved_local.get("cds_url") or "").strip()
                gff_url = str(resolved_local.get("gff_url") or "").strip()
                gbff_url = str(resolved_local.get("gbff_url") or "").strip()
                genome_url = str(resolved_local.get("genome_url") or "").strip()
                if cds_filename == "":
                    cds_filename = str(resolved_local.get("cds_filename") or "").strip()
                if gff_filename == "":
                    gff_filename = str(resolved_local.get("gff_filename") or "").strip()
                if gbff_filename == "":
                    gbff_filename = str(resolved_local.get("gbff_filename") or "").strip()
                if genome_filename == "":
                    genome_filename = str(resolved_local.get("genome_filename") or "").strip()
            except Exception as exc:
                errors.append(
                    "Manifest line {}: failed to resolve local paths from id '{}' (provider={}): {}".format(
                        i, source_id, provider, exc
                    )
                )
                continue
        elif provider in ("ncbi", "refseq", "genbank"):
            accession_candidate = extract_ncbi_accession_from_source_id(source_id)
            requires_ncbi_resolution = not manifest_has_usable_source_bundle(cds_url, gff_url, gbff_url, genome_url)
            should_enrich_from_ncbi = accession_candidate != "" and (species_key == "" or genome_url == "")
            if requires_ncbi_resolution or should_enrich_from_ncbi:
                try:
                    preferred_ncbi_source = "auto"
                    if provider in ("refseq", "genbank"):
                        preferred_ncbi_source = provider
                    resolved_ncbi = resolve_ncbi_download_urls_from_id(
                        source_id,
                        timeout=float(timeout),
                        ncbi_source=preferred_ncbi_source,
                    )
                except Exception as exc:
                    if requires_ncbi_resolution:
                        errors.append(
                            "Manifest line {}: failed to resolve id '{}' (provider={}): {}".format(
                                i, source_id, provider, exc
                            )
                        )
                        continue
                    warnings.append(
                        "Manifest line {}: skipped optional NCBI enrichment from id '{}' (provider={}): {}".format(
                            i, source_id, provider, exc
                        )
                    )
                    resolved_ncbi = None

        if species_key == "" and resolved_ncbi is not None:
            species_key = (resolved_ncbi.get("species_key") or "").strip()

        if provider not in ("local", ORYZA_MINUTA_PROVIDER) and (
            not manifest_has_usable_source_bundle(cds_url, gff_url, gbff_url, genome_url)
        ):
            try:
                resolved_urls = None
                if provider in ("ncbi", "refseq", "genbank"):
                    resolved_urls = resolved_ncbi
                else:
                    resolved_urls = resolve_download_urls_from_templates(provider, source_id, species_key, row)
                    if resolved_urls is None and (
                        not manifest_has_usable_source_bundle(cds_url, gff_url, gbff_url, genome_url)
                    ):
                        resolved_urls = resolve_non_ncbi_download_urls_from_id(
                            provider=provider,
                            source_id=source_id,
                            species_key=species_key,
                            timeout=float(timeout),
                            headers=headers,
                            row=row,
                        )
                if resolved_urls is not None:
                    gff_validation = str(resolved_urls.get("gff_validation") or "").strip()
                    if cds_url == "":
                        cds_url = (resolved_urls.get("cds_url") or "").strip()
                    if gff_url == "":
                        gff_url = (resolved_urls.get("gff_url") or "").strip()
                    if gbff_url == "":
                        gbff_url = (resolved_urls.get("gbff_url") or "").strip()
                    if genome_url == "":
                        genome_url = (resolved_urls.get("genome_url") or "").strip()
                    if cds_archive_member == "":
                        cds_archive_member = (resolved_urls.get("cds_archive_member") or "").strip()
                    if gff_archive_member == "":
                        gff_archive_member = (resolved_urls.get("gff_archive_member") or "").strip()
                    if gbff_archive_member == "":
                        gbff_archive_member = (resolved_urls.get("gbff_archive_member") or "").strip()
                    if genome_archive_member == "":
                        genome_archive_member = (resolved_urls.get("genome_archive_member") or "").strip()
                    resolved_cds_filename = (resolved_urls.get("cds_filename") or "").strip()
                    resolved_gff_filename = (resolved_urls.get("gff_filename") or "").strip()
                    resolved_gbff_filename = (resolved_urls.get("gbff_filename") or "").strip()
                    resolved_genome_filename = (resolved_urls.get("genome_filename") or "").strip()
                    if cds_filename == "" or (
                        cds_archive_member != "" and cds_filename == (row.get("cds_filename") or "").strip()
                    ):
                        cds_filename = resolved_cds_filename or cds_filename
                    if gff_filename == "" or (
                        gff_archive_member != "" and gff_filename == (row.get("gff_filename") or "").strip()
                    ):
                        gff_filename = resolved_gff_filename or gff_filename
                    if gbff_filename == "" or (
                        gbff_archive_member != "" and gbff_filename == (row.get("gbff_filename") or "").strip()
                    ):
                        gbff_filename = resolved_gbff_filename or gbff_filename
                    if genome_filename == "" or (
                        genome_archive_member != "" and genome_filename == (row.get("genome_filename") or "").strip()
                    ):
                        genome_filename = resolved_genome_filename or genome_filename
                    if species_key == "":
                        species_key = (resolved_urls.get("species_key") or "").strip()
            except Exception as exc:
                errors.append(
                    "Manifest line {}: failed to resolve urls from id '{}' (provider={}): {}".format(
                        i, source_id, provider, exc
                    )
                )
                continue

        if species_key == "":
            species_key = source_id
        if species_key == "":
            errors.append("Manifest line {}: species_key/id is empty".format(i))
            continue
        species_key = normalize_species_key_for_runtime(species_key)
        invalid_species_key = invalid_species_key_error(provider, species_key)
        if invalid_species_key != "":
            errors.append("Manifest line {}: {}".format(i, invalid_species_key))
            continue

        if provider != ORYZA_MINUTA_PROVIDER and not manifest_has_usable_source_bundle(
            cds_url, gff_url, gbff_url, genome_url
        ):
            errors.append(
                "Manifest line {}: CDS, GBFF, or GFF plus genome are required for {} "
                "(set urls directly or provide resolvable id)".format(i, species_key)
            )
            continue

        if provider in ENSEMBL_LIKE_PROVIDERS and gff_archive_member == "":
            try:
                preferred_gff_url = resolve_preferred_ensembl_like_gff_url(
                    provider,
                    gff_url,
                    float(timeout),
                    headers,
                )
                if preferred_gff_url != gff_url:
                    warnings.append(
                        "Manifest line {}: replaced partial Ensembl-like GFF '{}' with '{}'".format(
                            i,
                            Path(urlparse(gff_url).path).name,
                            Path(urlparse(preferred_gff_url).path).name,
                        )
                    )
                    gff_url = preferred_gff_url
                    gff_filename = ""
            except Exception as exc:
                warnings.append(
                    "Manifest line {}: could not check Ensembl-like GFF alternatives for '{}': {}".format(
                        i,
                        Path(urlparse(gff_url).path).name if gff_url != "" else "",
                        exc,
                    )
                )

        raw_dir = provider_raw_dir(provider, download_root, species_key)
        raw_dir.mkdir(parents=True, exist_ok=True)

        if cds_filename == "":
            cds_filename = default_download_filename(provider, species_key, "cds", cds_url, cds_archive_member)
        if gff_url != "" and gff_filename == "":
            gff_filename = default_download_filename(provider, species_key, "gff", gff_url, gff_archive_member)
        if gbff_url != "" and gbff_filename == "":
            gbff_filename = default_download_filename(provider, species_key, "gbff", gbff_url, gbff_archive_member)
        if genome_url != "" and genome_filename == "":
            genome_filename = default_download_filename(
                provider,
                species_key,
                "genome",
                genome_url,
                genome_archive_member,
            )

        resolved_cds_filename = cds_filename
        resolved_gff_filename = gff_filename
        download_targets = []
        if provider == ORYZA_MINUTA_PROVIDER:
            representative_urls = (
                str(oryza_minuta_bundle.get("cds_url") or "").strip(),
                str(oryza_minuta_bundle.get("gff_url") or "").strip(),
                str(oryza_minuta_bundle.get("genome_url") or "").strip(),
            )
            if cds_url == "":
                cds_url = representative_urls[0]
            if gff_url == "":
                gff_url = representative_urls[1]
            if genome_url == "":
                genome_url = representative_urls[2]
            download_targets = [
                (
                    target["label"],
                    target["url"],
                    raw_dir / target["filename"],
                    "",
                )
                for target in oryza_minuta_bundle.get("download_targets", ())
            ]
            merge_jobs.append(
                {
                    "provider": ORYZA_MINUTA_PROVIDER,
                    "row_id": i,
                    "species_key": species_key,
                    "cds_paths": tuple(
                        raw_dir / target["filename"]
                        for target in oryza_minuta_bundle.get("download_targets", ())
                        if str(target.get("label") or "").startswith("CDS_")
                    ),
                    "gff_paths": tuple(
                        raw_dir / target["filename"]
                        for target in oryza_minuta_bundle.get("download_targets", ())
                        if str(target.get("label") or "").startswith("GFF_")
                    ),
                    "genome_paths": tuple(
                        raw_dir / target["filename"]
                        for target in oryza_minuta_bundle.get("download_targets", ())
                        if str(target.get("label") or "").startswith("GENOME_")
                    ),
                    "combined_cds_path": raw_dir / cds_filename,
                    "combined_gff_path": raw_dir / gff_filename,
                    "combined_genome_path": raw_dir / genome_filename,
                }
            )
        else:
            if cds_url != "":
                download_targets.append(("CDS", cds_url, raw_dir / cds_filename, cds_archive_member))
            if gff_url != "":
                download_targets.append(("GFF", gff_url, raw_dir / gff_filename, gff_archive_member))
            if gbff_url != "":
                download_targets.append(("GBFF", gbff_url, raw_dir / gbff_filename, gbff_archive_member))
            if genome_url != "":
                download_targets.append(("GENOME", genome_url, raw_dir / genome_filename, genome_archive_member))

        if provider == "fernbase":
            effective_mode = effective_fernbase_confidence_mode(fernbase_confidence_mode)
            if effective_mode == FERNBASE_CONFIDENCE_MODE_HIGH_LOW_COMBINED:
                try:
                    split_bundle = resolve_fernbase_split_bundle(
                        cds_url=cds_url,
                        gff_url=gff_url,
                        genome_url=genome_url,
                        timeout=float(timeout),
                        headers=headers,
                    )
                except Exception as exc:
                    errors.append(
                        "Manifest line {}: failed to resolve FernBase high/low confidence bundle for '{}': {}".format(
                            i,
                            species_key,
                            exc,
                        )
                    )
                    continue
                if split_bundle.get("is_split"):
                    high_cds_filename = str(split_bundle.get("high_cds_filename") or cds_filename).strip()
                    low_cds_filename = str(split_bundle.get("low_cds_filename") or "").strip()
                    high_gff_filename = str(split_bundle.get("high_gff_filename") or gff_filename).strip()
                    low_gff_filename = str(split_bundle.get("low_gff_filename") or "").strip()
                    if high_cds_filename != "":
                        cds_filename = high_cds_filename
                    if high_gff_filename != "":
                        gff_filename = high_gff_filename
                    resolved_cds_filename = build_fernbase_combined_filename(cds_filename)
                    resolved_gff_filename = build_fernbase_combined_filename(gff_filename)
                    download_targets = []
                    if split_bundle.get("high_cds_url", "") != "":
                        download_targets.append(
                            ("CDS", split_bundle["high_cds_url"], raw_dir / cds_filename, cds_archive_member)
                        )
                    if split_bundle.get("high_gff_url", "") != "":
                        download_targets.append(
                            ("GFF", split_bundle["high_gff_url"], raw_dir / gff_filename, gff_archive_member)
                        )
                    if split_bundle.get("low_cds_url", "") != "":
                        download_targets.append(
                            ("CDS_LOW", split_bundle["low_cds_url"], raw_dir / low_cds_filename, "")
                        )
                    if split_bundle.get("low_gff_url", "") != "":
                        download_targets.append(
                            ("GFF_LOW", split_bundle["low_gff_url"], raw_dir / low_gff_filename, "")
                        )
                    if genome_url != "":
                        download_targets.append(
                            ("GENOME", genome_url, raw_dir / genome_filename, genome_archive_member)
                        )
                    merge_jobs.append(
                        {
                            "row_id": i,
                            "species_key": species_key,
                            "high_cds_path": raw_dir / cds_filename,
                            "low_cds_path": raw_dir / low_cds_filename,
                            "high_gff_path": raw_dir / gff_filename,
                            "low_gff_path": raw_dir / low_gff_filename,
                            "combined_cds_path": raw_dir / resolved_cds_filename,
                            "combined_gff_path": raw_dir / resolved_gff_filename,
                        }
                    )
                else:
                    warnings.append(
                        "Manifest line {}: {}='{}' ignored for FernBase row '{}' because no legacy split high/low bundle was found".format(
                            i,
                            FERNBASE_CONFIDENCE_MODE_FIELD,
                            effective_mode,
                            species_key,
                        )
                    )
            elif ("highconfidence" in cds_filename.lower() or "lowconfidence" in cds_filename.lower()) and (
                "highconfidence" in gff_filename.lower() or "lowconfidence" in gff_filename.lower()
            ):
                stale_combined_cds = raw_dir / build_fernbase_combined_filename(cds_filename)
                stale_combined_gff = raw_dir / build_fernbase_combined_filename(gff_filename)
                if stale_combined_cds != raw_dir / cds_filename:
                    cleanup_paths.append(stale_combined_cds)
                if stale_combined_gff != raw_dir / gff_filename:
                    cleanup_paths.append(stale_combined_gff)

        resolved_row = build_resolved_manifest_row(
            raw_row=row,
            fieldnames=resolved_fieldnames,
            provider=provider,
            source_id=source_id,
            species_key=species_key,
            cds_url=cds_url,
            gff_url=gff_url,
            gbff_url=gbff_url,
            genome_url=genome_url,
            cds_archive_member=cds_archive_member,
            gff_archive_member=gff_archive_member,
            gbff_archive_member=gbff_archive_member,
            genome_archive_member=genome_archive_member,
            cds_filename=resolved_cds_filename,
            gff_filename=resolved_gff_filename,
            gbff_filename=gbff_filename,
            genome_filename=genome_filename,
        )
        resolved_row[FERNBASE_CONFIDENCE_MODE_FIELD] = fernbase_confidence_mode
        resolved_rows.append(resolved_row)
        row_target_paths[i] = {
            "provider": provider,
            "source_id": source_id,
            "species_key": species_key,
            "gff_validation": gff_validation,
            "paths": {label: target for label, _url, target, _archive_member in download_targets},
        }
        for label, url, target, archive_member in download_targets:
            if target.exists() and target.stat().st_size > 0 and not overwrite:
                if (not dry_run) and quarantine_corrupt_gzip(
                    target,
                    warnings,
                    "[download:{}] {} {}".format(provider, species_key, label),
                ):
                    pass
                else:
                    warnings.append(
                        "[download:{}] {} {} already exists. Skipping: {}".format(provider, species_key, label, target)
                    )
                    continue
            if dry_run:
                planned += 1
                continue
            download_jobs.append(
                {
                    "row_id": i,
                    "provider": provider,
                    "source_id": source_id,
                    "species_key": species_key,
                    "label": label,
                    "url": url,
                    "target": target,
                    "archive_member": archive_member,
                }
            )

    if len(download_jobs) > 0:
        max_workers = max(1, int(jobs))
        provider_limits = resolve_provider_download_limits(max_workers)
        provider_semaphores = {
            provider_name: threading.Semaphore(provider_limits.get(provider_name, 1)) for provider_name in PROVIDERS
        }

        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = [
                pool.submit(
                    execute_download_target_job,
                    job,
                    headers,
                    timeout,
                    overwrite,
                    lock_stale_seconds,
                    provider_semaphores,
                )
                for job in download_jobs
            ]
            for future in as_completed(futures):
                try:
                    result = future.result()
                except Exception as exc:
                    errors.append("Unhandled download worker error: {}".format(exc))
                    continue
                warnings.extend(result.get("warnings", []))
                errors.extend(result.get("errors", []))
                downloaded += int(result.get("downloaded", 0))
                failed_downloads.extend(result.get("failed", []))

    if not dry_run:
        for row_info in row_target_paths.values():
            if row_info.get("gff_validation") != "coge_export_cds":
                continue
            gff_path = row_info.get("paths", {}).get("GFF")
            if gff_path is None or not gff_path.exists() or gff_path.stat().st_size == 0:
                continue
            try:
                validate_coge_export_gff_file(gff_path, gid=row_info.get("source_id", ""))
            except Exception as exc:
                errors.append(
                    "[download:coge] {} invalid GFF export at {} ({})".format(
                        row_info.get("species_key", ""),
                        gff_path,
                        exc,
                    )
                )

        for cleanup_path in cleanup_paths:
            try:
                if cleanup_path.exists() and cleanup_path.is_file():
                    cleanup_path.unlink()
            except OSError as exc:
                warnings.append("Failed to remove stale FernBase combined file '{}': {}".format(cleanup_path, exc))

        for merge_job in merge_jobs:
            try:
                merge_provider = str(merge_job.get("provider") or "fernbase").strip().lower()
                if merge_provider == "fernbase":
                    result = merge_fernbase_confidence_bundle(
                        species_key=merge_job["species_key"],
                        high_cds_path=merge_job["high_cds_path"],
                        low_cds_path=merge_job["low_cds_path"],
                        high_gff_path=merge_job["high_gff_path"],
                        low_gff_path=merge_job["low_gff_path"],
                        combined_cds_path=merge_job["combined_cds_path"],
                        combined_gff_path=merge_job["combined_gff_path"],
                        overwrite=overwrite,
                    )
                elif merge_provider == ORYZA_MINUTA_PROVIDER:
                    result = merge_oryza_minuta_multisource_bundle(
                        species_key=merge_job["species_key"],
                        cds_paths=merge_job["cds_paths"],
                        gff_paths=merge_job["gff_paths"],
                        genome_paths=merge_job["genome_paths"],
                        combined_cds_path=merge_job["combined_cds_path"],
                        combined_gff_path=merge_job["combined_gff_path"],
                        combined_genome_path=merge_job["combined_genome_path"],
                        overwrite=overwrite,
                    )
                else:
                    warnings.append("Unknown merge job provider '{}'. Skipping.".format(merge_provider))
                    continue
            except Exception as exc:
                if str(merge_job.get("provider") or "") == ORYZA_MINUTA_PROVIDER:
                    errors.append(
                        "[download:{}] {} failed to merge public multi-source bundle ({})".format(
                            ORYZA_MINUTA_PROVIDER,
                            merge_job["species_key"],
                            exc,
                        )
                    )
                else:
                    errors.append(
                        "[download:fernbase] {} failed to merge high/low confidence bundle ({})".format(
                            merge_job["species_key"],
                            exc,
                        )
                    )
                continue
            warnings.extend(result.get("warnings", []))
            errors.extend(result.get("errors", []))

    for failure in failed_downloads:
        row_info = row_target_paths.get(failure.get("row_id"), {})
        paths = row_info.get("paths", {})
        cds_ok = paths.get("CDS") is not None and paths["CDS"].exists() and paths["CDS"].stat().st_size > 0
        gff_ok = paths.get("GFF") is not None and paths["GFF"].exists() and paths["GFF"].stat().st_size > 0
        gbff_ok = paths.get("GBFF") is not None and paths["GBFF"].exists() and paths["GBFF"].stat().st_size > 0
        genome_ok = paths.get("GENOME") is not None and paths["GENOME"].exists() and paths["GENOME"].stat().st_size > 0
        if cds_ok or gbff_ok or (gff_ok and genome_ok):
            warnings.append(
                "{} ; continuing because a usable source bundle is available".format(failure.get("message", ""))
            )
        else:
            errors.append(failure.get("message", "download failed"))

    if processed == 0:
        warnings.append("No manifest rows matched --provider {}.".format(provider_filter))

    if resolved_manifest_output_path is not None:
        try:
            write_resolved_manifest_tsv(resolved_manifest_output_path, resolved_fieldnames, resolved_rows)
        except Exception as exc:
            errors.append("Failed to write resolved manifest TSV '{}': {}".format(resolved_manifest_output_path, exc))

    final_cache_diagnostics = scan_download_cache_diagnostics(download_root)
    return {
        "warnings": warnings,
        "errors": errors,
        "processed": processed,
        "downloaded": downloaded,
        "planned": planned,
        "resolved_rows": resolved_rows,
        "download_diagnostics": summarize_download_diagnostics(
            preexisting_cache_diagnostics,
            final_cache_diagnostics,
            warnings,
            len(download_jobs),
            len(failed_downloads),
        ),
        "resolved_manifest_output": str(resolved_manifest_output_path)
        if resolved_manifest_output_path is not None
        else "",
    }

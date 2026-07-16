"""Download runtime implementation: local."""

import gzip
import os
import re
import time
from pathlib import Path
from urllib.parse import unquote, urlparse

from format_species_common import (
    is_fasta_filename,
    is_gbff_filename,
    is_gff_filename,
    is_probable_genome_filename,
    pick_single_file,
)
from format_species_manifest import (
    manifest_has_usable_source_bundle,
)
from format_species_provider_urls import (
    is_url_like,
)


def is_gzip_path(path):
    return str(path).lower().endswith(".gz")


def gzip_integrity_error(path):
    if not is_gzip_path(path):
        return None
    try:
        with gzip.open(path, "rb") as handle:
            while handle.read(1024 * 1024):
                pass
    except Exception as exc:
        return exc
    return None


def quarantine_existing_file(path, warnings, context, reason):
    suffix = ".corrupt.{}.{}".format(time.strftime("%Y%m%d%H%M%S"), os.getpid())
    quarantine_path = Path(str(path) + suffix)
    counter = 1
    while quarantine_path.exists():
        counter += 1
        quarantine_path = Path(str(path) + suffix + ".{}".format(counter))
    path.replace(quarantine_path)
    warnings.append("{} found corrupt gzip cache {}; moved to {} ({})".format(context, path, quarantine_path, reason))
    return quarantine_path


def quarantine_corrupt_gzip(path, warnings, context):
    validation_error = gzip_integrity_error(path)
    if validation_error is None:
        return False
    quarantine_existing_file(path, warnings, context, validation_error)
    return True


def resolve_local_reference_path(reference, manifest_parent_dir):
    text = str(reference or "").strip()
    if text == "":
        return None
    if is_url_like(text):
        parsed = urlparse(text)
        if parsed.scheme.lower() != "file":
            raise ValueError("expected local file path or file:// URL, got '{}'".format(reference))
        path = Path(unquote(parsed.path)).expanduser()
        return path.resolve()
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = (manifest_parent_dir / path).resolve()
    else:
        path = path.resolve()
    return path


def resolve_local_source_id_with_label_fallback(source_id, manifest_parent_dir):
    text = str(source_id or "").strip()
    if text == "":
        return text
    try:
        direct_path = resolve_local_reference_path(text, manifest_parent_dir)
    except Exception:
        direct_path = None
    if direct_path is not None and direct_path.exists():
        return text

    # Accept spreadsheet label format like:
    #   /path/to/species_dir (Species name)
    # only when the stripped path actually exists.
    match = re.match(r"^(.*\S)\s+\([^()]*\)$", text)
    if match is None:
        return text
    stripped = match.group(1).strip()
    if stripped == "":
        return text
    try:
        stripped_path = resolve_local_reference_path(stripped, manifest_parent_dir)
    except Exception:
        return text
    if stripped_path is not None and stripped_path.exists():
        return stripped
    return text


def local_reference_to_file_url(reference, manifest_parent_dir):
    path = resolve_local_reference_path(reference, manifest_parent_dir)
    if path is None:
        return ""
    return path.as_uri()


def resolve_local_manifest_row(provider, source_id, species_key, row, manifest_parent_dir, warnings, line_number):
    cds_url = local_reference_to_file_url(row.get("cds_url", ""), manifest_parent_dir)
    gff_url = local_reference_to_file_url(row.get("gff_url", ""), manifest_parent_dir)
    gbff_url = local_reference_to_file_url(row.get("gbff_url", ""), manifest_parent_dir)
    genome_url = local_reference_to_file_url(row.get("genome_url", ""), manifest_parent_dir)
    cds_filename = str(row.get("cds_filename") or "").strip()
    gff_filename = str(row.get("gff_filename") or "").strip()
    gbff_filename = str(row.get("gbff_filename") or "").strip()
    genome_filename = str(row.get("genome_filename") or "").strip()

    local_cds_path = str(row.get("local_cds_path") or "").strip()
    local_gff_path = str(row.get("local_gff_path") or "").strip()
    local_gbff_path = str(row.get("local_gbff_path") or "").strip()
    local_genome_path = str(row.get("local_genome_path") or "").strip()
    if cds_url == "" and local_cds_path != "":
        cds_url = local_reference_to_file_url(local_cds_path, manifest_parent_dir)
    if gff_url == "" and local_gff_path != "":
        gff_url = local_reference_to_file_url(local_gff_path, manifest_parent_dir)
    if gbff_url == "" and local_gbff_path != "":
        gbff_url = local_reference_to_file_url(local_gbff_path, manifest_parent_dir)
    if genome_url == "" and local_genome_path != "":
        genome_url = local_reference_to_file_url(local_genome_path, manifest_parent_dir)

    source_id_text = str(source_id or "").strip()
    source_id_resolved = resolve_local_source_id_with_label_fallback(source_id_text, manifest_parent_dir)
    if source_id_resolved != source_id_text:
        warnings.append(
            "Manifest line {}: provider=local id looked like a labeled choice. Using '{}'".format(
                line_number, source_id_resolved
            )
        )
    source_dir = resolve_local_reference_path(source_id_resolved, manifest_parent_dir)
    if species_key == "":
        if source_dir is not None and source_dir.exists() and source_dir.is_dir():
            species_key = source_dir.name
        else:
            species_key = source_id_resolved

    if (
        (not manifest_has_usable_source_bundle(cds_url, gff_url, gbff_url, genome_url))
        and source_dir is not None
        and source_dir.exists()
        and source_dir.is_dir()
    ):
        files = [path for path in source_dir.iterdir() if path.is_file()]
        cds_matches = [
            path
            for path in files
            if is_fasta_filename(path.name)
            and any(marker in path.name.lower() for marker in ("cds", "transcript", "mrna", "cdna"))
        ]
        if len(cds_matches) == 0:
            cds_matches = [
                path
                for path in files
                if is_fasta_filename(path.name) and not is_probable_genome_filename(provider, path.name)
            ]
        gff_matches = [path for path in files if is_gff_filename(path.name)]
        gbff_matches = [path for path in files if is_gbff_filename(path.name)]
        genome_matches = [path for path in files if is_probable_genome_filename(provider, path.name)]

        cds_path = pick_single_file(cds_matches, provider, species_key, "CDS", warnings)
        gff_path = pick_single_file(gff_matches, provider, species_key, "GFF", warnings)
        gbff_path = pick_single_file(gbff_matches, provider, species_key, "GBFF", warnings)
        genome_path = pick_single_file(genome_matches, provider, species_key, "genome", warnings)

        if cds_url == "" and cds_path is not None:
            cds_url = cds_path.resolve().as_uri()
            if cds_filename == "":
                cds_filename = cds_path.name
        if gff_url == "" and gff_path is not None:
            gff_url = gff_path.resolve().as_uri()
            if gff_filename == "":
                gff_filename = gff_path.name
        if gbff_url == "" and gbff_path is not None:
            gbff_url = gbff_path.resolve().as_uri()
            if gbff_filename == "":
                gbff_filename = gbff_path.name
        if genome_url == "" and genome_path is not None:
            genome_url = genome_path.resolve().as_uri()
            if genome_filename == "":
                genome_filename = genome_path.name

    if not manifest_has_usable_source_bundle(cds_url, gff_url, gbff_url, genome_url):
        raise ValueError(
            "provider=local requires CDS, GBFF, or GFF plus genome paths/URLs, or a species directory id. "
            "line {} id='{}'".format(line_number, source_id_resolved)
        )

    return {
        "species_key": species_key,
        "cds_url": cds_url,
        "gff_url": gff_url,
        "gbff_url": gbff_url,
        "genome_url": genome_url,
        "cds_filename": cds_filename,
        "gff_filename": gff_filename,
        "gbff_filename": gbff_filename,
        "genome_filename": genome_filename,
    }

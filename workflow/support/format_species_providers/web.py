"""Provider resolver implementation: web."""

import html
import json
import re
import threading
from pathlib import Path
from urllib.parse import quote, unquote, urljoin, urlparse

from format_species_annotations import (
    sanitize_identifier,
)
from format_species_common import (
    is_fasta_filename,
    is_gff_filename,
    is_probable_cds_filename,
    is_probable_genome_filename,
    provider_candidate_sort_key,
)
from format_species_constants import (
    INSECTBASE_IBG_ID_PATTERN,
    PLANTGARDEN_GENOME_ID_PATTERN,
    PLANTGARDEN_TAXID_PATTERN,
)
from format_species_manifest import manifest_has_usable_source_bundle
from format_species_provider_urls import (
    fetch_text_with_headers,
    is_url_like,
    parse_links_from_document,
    resolve_figshare_api_base_url,
    resolve_insectbase_api_base_url,
    resolve_plantgarden_web_base_url,
    resolve_veupathdb_service_base_url,
    strip_provider_prefix,
)

from .common import (
    extract_figshare_article_id_candidate,
    fetch_json_with_headers,
    is_probable_gff_url,
    normalize_lookup_text,
    parse_species_key_candidate,
    source_id_candidates,
)


def extract_insectbase_ibg_id_candidate(source_id):
    text = str(source_id or "").strip()
    if text == "":
        return ""
    stripped = strip_provider_prefix(text, "insectbase")
    for candidate in (text, stripped):
        if INSECTBASE_IBG_ID_PATTERN.match(candidate):
            return candidate.upper()
        match = re.search(r"(IBG_[0-9]+)", candidate, flags=re.IGNORECASE)
        if match is not None:
            return match.group(1).upper()
    if is_url_like(stripped):
        parts = [unquote(part) for part in urlparse(stripped).path.split("/") if part]
        for part in reversed(parts):
            if INSECTBASE_IBG_ID_PATTERN.match(part):
                return part.upper()
    return ""


def strip_html_tags(text):
    return re.sub(r"<[^>]+>", " ", str(text or ""))


def infer_figshare_species_key(article_payload, fallback):
    fallback_text = str(fallback or "").strip()
    if fallback_text != "":
        return fallback_text
    title = strip_html_tags(html.unescape(str((article_payload or {}).get("title", "") or "")))
    match = re.fullmatch(r"\s*([A-Z][a-z]+)\s+([a-z][a-z-]+)\s*", title)
    if match is None:
        return ""
    return "{}_{}".format(match.group(1), match.group(2))


def figshare_file_download_url(file_record):
    direct = str((file_record or {}).get("download_url", "") or "").strip()
    if direct != "":
        return direct
    file_id = str((file_record or {}).get("id", "") or "").strip()
    if file_id == "":
        return ""
    return "https://ndownloader.figshare.com/files/{}".format(file_id)


def figshare_candidate_name(file_record):
    return Path(str((file_record or {}).get("name", "") or "").strip()).name


def normalize_figshare_filename_hint(value):
    return Path(str(value or "").strip()).name.lower()


def filter_figshare_files_by_requested_name(files, requested_name):
    hint = normalize_figshare_filename_hint(requested_name)
    if hint == "":
        return list(files)
    return [file_record for file_record in files if figshare_candidate_name(file_record).lower() == hint]


def filter_figshare_files_by_species_key(files, species_key):
    text = str(species_key or "").strip().lower()
    if text == "":
        return list(files)
    tokens = [token for token in re.split(r"[_\s]+", text) if token]
    patterns = []
    if len(tokens) >= 2:
        patterns.extend(
            (
                "_".join(tokens),
                "-".join(tokens),
                "".join(tokens),
                "{} {}".format(tokens[0], tokens[1]),
                "{}_{}".format(tokens[0], tokens[-1]),
                "{} {}".format(tokens[0], tokens[-1]),
            )
        )
    patterns.extend(tokens)
    unique_patterns = []
    for pattern in patterns:
        if pattern and pattern not in unique_patterns:
            unique_patterns.append(pattern)
    out = []
    for file_record in files:
        lowered = figshare_candidate_name(file_record).lower()
        if any(pattern in lowered for pattern in unique_patterns):
            out.append(file_record)
    return out


def select_figshare_file_record(files, provider_label, requested_name, species_key, article_id):
    if str(requested_name or "").strip() != "":
        requested_candidates = filter_figshare_files_by_requested_name(files, requested_name)
        if len(requested_candidates) == 0:
            raise ValueError(
                "figshare article '{}' did not contain requested {} file '{}'".format(
                    article_id,
                    provider_label.lower(),
                    requested_name,
                )
            )
        ordered = sorted(
            requested_candidates,
            key=lambda record: provider_candidate_sort_key("figshare", provider_label, figshare_candidate_name(record)),
        )
        return ordered[0]

    filtered = []
    for file_record in files:
        name = figshare_candidate_name(file_record)
        if provider_label == "CDS" and is_probable_cds_filename("figshare", name):
            filtered.append(file_record)
        elif provider_label == "GFF" and is_probable_gff_url(name):
            filtered.append(file_record)
        elif provider_label == "GENOME" and is_probable_genome_filename("figshare", name):
            filtered.append(file_record)
    if len(filtered) == 0:
        return None
    species_filtered = filter_figshare_files_by_species_key(filtered, species_key)
    if len(species_filtered) == 1:
        return species_filtered[0]
    if len(filtered) == 1:
        return filtered[0]
    raise ValueError(
        "figshare article '{}' exposed multiple {} candidates; set {}_filename to disambiguate".format(
            article_id,
            provider_label.lower(),
            provider_label.lower(),
        )
    )


def resolve_figshare_download_urls_from_id(source_id, species_key, timeout, headers, row=None):
    article_id = extract_figshare_article_id_candidate(source_id)
    if article_id == "":
        raise ValueError("figshare id '{}' did not contain a public article id".format(source_id))
    article_url = "{}/articles/{}".format(resolve_figshare_api_base_url(), article_id)
    payload = fetch_json_with_headers(article_url, timeout, headers)
    files = payload.get("files", [])
    if not isinstance(files, list) or len(files) == 0:
        files = fetch_json_with_headers(article_url + "/files?page_size=1000", timeout, headers)
    if not isinstance(files, list) or len(files) == 0:
        raise ValueError("figshare article '{}' did not expose downloadable files".format(article_id))
    row = dict(row or {})
    resolved = {}
    for label, field_prefix in (("CDS", "cds"), ("GFF", "gff"), ("GENOME", "genome")):
        selected = select_figshare_file_record(
            files=files,
            provider_label=label,
            requested_name=row.get("{}_filename".format(field_prefix), ""),
            species_key=species_key,
            article_id=article_id,
        )
        if selected is None:
            continue
        resolved["{}_url".format(field_prefix)] = figshare_file_download_url(selected)
        archive_member = str(row.get("{}_archive_member".format(field_prefix), "") or "").strip()
        if archive_member != "":
            resolved["{}_archive_member".format(field_prefix)] = archive_member
            resolved["{}_filename".format(field_prefix)] = Path(archive_member).name
        else:
            resolved["{}_filename".format(field_prefix)] = figshare_candidate_name(selected)
    inferred_species_key = infer_figshare_species_key(payload, species_key)
    if inferred_species_key != "":
        resolved["species_key"] = inferred_species_key
    if not manifest_has_usable_source_bundle(
        resolved.get("cds_url", ""),
        resolved.get("gff_url", ""),
        "",
        resolved.get("genome_url", ""),
    ):
        raise ValueError("figshare article '{}' did not expose CDS, GBFF, or GFF plus genome".format(article_id))
    return resolved


def normalize_plantgarden_source_page_url(source_id):
    source_clean = strip_provider_prefix(source_id, "plantgarden").strip()
    if source_clean == "":
        raise ValueError("PlantGARDEN id is empty")
    base = resolve_plantgarden_web_base_url().rstrip("/") + "/"
    if is_url_like(source_clean):
        return source_clean
    if source_clean.startswith("/"):
        return urljoin(base, source_clean.lstrip("/"))
    genome_match = PLANTGARDEN_GENOME_ID_PATTERN.match(source_clean)
    if genome_match is not None:
        return urljoin(base, "en/list/{}/genome/{}".format(genome_match.group(1), source_clean))
    if PLANTGARDEN_TAXID_PATTERN.match(source_clean):
        return urljoin(base, "en/list/{}/genome".format(source_clean))
    lowered = source_clean.lower()
    if lowered.startswith(("en/", "ja/", "list/", "download/")):
        return urljoin(base, source_clean)
    return urljoin(base, source_clean)


def extract_plantgarden_page_metadata(page_text):
    metadata = {}
    for key, value in re.findall(r'data-pgtag-([a-z_]+)="([^"]*)"', str(page_text or "")):
        metadata[str(key or "").strip()] = str(value or "").strip()
    return metadata


def infer_plantgarden_species_key(metadata, fallback):
    species_name = str((metadata or {}).get("species_name", "") or "").strip()
    inferred = parse_species_key_candidate(species_name)
    if inferred != "":
        return inferred
    return str(fallback or "").strip()


def select_best_plantgarden_candidate(links, label):
    best_url = ""
    best_rank = None
    for link in links:
        parsed = urlparse(link)
        filename = Path(parsed.path).name
        lower = filename.lower()
        if lower in ("", ".", ".."):
            continue
        if label == "CDS":
            if not is_fasta_filename(lower):
                continue
            if any(marker in lower for marker in ("protein", "pep")):
                continue
            if "cds" in lower:
                rank = (0, 0, lower)
            elif any(marker in lower for marker in ("transcript", "mrna", "cdna")):
                rank = (0, 1, lower)
            else:
                continue
        elif label == "GFF":
            if not is_gff_filename(lower):
                continue
            if any(marker in lower for marker in ("repeat", "repeatmask")):
                continue
            rank = (0, lower)
        elif label == "GENOME":
            if not is_fasta_filename(lower):
                continue
            if any(marker in lower for marker in ("cds", "transcript", "mrna", "cdna", "protein", "pep", "repeat")):
                continue
            preferred = any(
                marker in lower
                for marker in ("genome", "assembly", ".chr", "chr.", "chrom", "scaffold", "pmol", "pseudomol")
            )
            rank = (0 if preferred else 1, lower)
        else:
            continue
        if best_url == "" or rank < best_rank:
            best_url = link
            best_rank = rank
    return best_url


def resolve_plantgarden_download_page_url(source_id, timeout, headers):
    page_url = normalize_plantgarden_source_page_url(source_id)
    parsed = urlparse(page_url)
    locale = "ja" if "/ja/" in parsed.path.lower() else "en"
    direct_match = re.search(r"/download/(t[0-9]+)/(t[0-9]+[.]G[0-9]+)/?$", parsed.path, flags=re.IGNORECASE)
    if direct_match is not None:
        return page_url.rstrip("/") + "/"

    page_text = fetch_text_with_headers(page_url, timeout, headers)
    metadata = extract_plantgarden_page_metadata(page_text)
    species_id = str(metadata.get("species_id", "") or "").strip()
    genome_id = str(metadata.get("genome_assembly_id", "") or "").strip()
    if species_id != "" and genome_id != "":
        return urljoin(
            resolve_plantgarden_web_base_url().rstrip("/") + "/",
            "{}/download/{}/{}/".format(locale, species_id, genome_id),
        )

    links = parse_links_from_document(page_url, page_text)
    for link in links:
        link_path = urlparse(link).path
        if re.search(r"/download/(t[0-9]+)/(t[0-9]+[.]G[0-9]+)/?$", link_path, flags=re.IGNORECASE):
            return link.rstrip("/") + "/"
    for link in links:
        link_path = urlparse(link).path
        genome_match = re.search(r"/list/(t[0-9]+)/genome/(t[0-9]+[.]G[0-9]+)/?$", link_path, flags=re.IGNORECASE)
        if genome_match is not None:
            return urljoin(
                resolve_plantgarden_web_base_url().rstrip("/") + "/",
                "{}/download/{}/{}/".format(locale, genome_match.group(1), genome_match.group(2)),
            )
    for link in links:
        link_path = urlparse(link).path
        if re.search(r"/list/t[0-9]+/genome/?$", link_path, flags=re.IGNORECASE) and link.rstrip(
            "/"
        ) != page_url.rstrip("/"):
            return resolve_plantgarden_download_page_url(link, timeout, headers)
    raise ValueError("PlantGARDEN page did not expose a downloadable assembly page: {}".format(page_url))


def resolve_plantgarden_download_urls_from_id(source_id, species_key, timeout, headers):
    page_url = normalize_plantgarden_source_page_url(source_id)
    page_text = fetch_text_with_headers(page_url, timeout, headers)
    metadata = extract_plantgarden_page_metadata(page_text)
    download_page_url = resolve_plantgarden_download_page_url(page_url, timeout, headers)
    download_text = fetch_text_with_headers(download_page_url, timeout, headers)
    links = parse_links_from_document(download_page_url, download_text)
    cds_url = select_best_plantgarden_candidate(links, "CDS")
    gff_url = select_best_plantgarden_candidate(links, "GFF")
    genome_url = select_best_plantgarden_candidate(links, "GENOME")
    if not manifest_has_usable_source_bundle(cds_url, gff_url, "", genome_url):
        raise ValueError("PlantGARDEN id '{}' did not resolve to a usable CDS-or-(GFF+genome) bundle".format(source_id))
    inferred_species_key = str(species_key or "").strip() or infer_plantgarden_species_key(metadata, species_key)
    return {
        "species_key": inferred_species_key,
        "cds_url": cds_url,
        "gff_url": gff_url,
        "genome_url": genome_url,
        "cds_filename": Path(urlparse(cds_url).path).name if cds_url != "" else "",
        "gff_filename": Path(urlparse(gff_url).path).name if gff_url != "" else "",
        "genome_filename": Path(urlparse(genome_url).path).name if genome_url != "" else "",
    }


_veupathdb_records_cache = {}


_veupathdb_records_lock = threading.Lock()


def veupathdb_source_id_from_url(url):
    path = urlparse(str(url or "")).path
    parts = [unquote(part) for part in path.split("/") if part]
    for idx, token in enumerate(parts[:-1]):
        if token == "Current_Release" and idx + 1 < len(parts):
            return parts[idx + 1]
    return ""


def infer_veupathdb_cds_url(attributes):
    protein_url = str(attributes.get("URLproteinFasta", "") or "").strip()
    if protein_url != "":
        return re.sub(
            r"_AnnotatedProteins[.]fasta$",
            "_AnnotatedCDSs.fasta",
            protein_url,
            flags=re.IGNORECASE,
        )
    genome_url = str(attributes.get("URLGenomeFasta", "") or "").strip()
    if genome_url == "":
        return ""
    parsed = urlparse(genome_url)
    directory = parsed.path.rsplit("/", 1)[0]
    basename = parsed.path.rsplit("/", 1)[-1]
    stem = re.sub(r"_Genome[.]fasta$", "", basename, flags=re.IGNORECASE)
    if stem == basename:
        return ""
    return parsed._replace(path="{}/{}_AnnotatedCDSs.fasta".format(directory, stem)).geturl()


def fetch_veupathdb_records(timeout, headers):
    service_base = resolve_veupathdb_service_base_url()
    with _veupathdb_records_lock:
        cached = _veupathdb_records_cache.get(service_base)
    if cached is not None:
        return cached
    config = {
        "attributes": ["primary_key", "species", "URLGenomeFasta", "URLgff", "URLproteinFasta", "project_id"],
        "tables": [],
        "attributeFormat": "text",
    }
    service_url = "{}/record-types/organism/searches/GenomeDataTypes/reports/standard?reportConfig={}".format(
        service_base,
        quote(json.dumps(config, separators=(",", ":")), safe=""),
    )
    payload = fetch_json_with_headers(service_url, timeout, headers)
    records = list(payload.get("records", []))
    with _veupathdb_records_lock:
        _veupathdb_records_cache[service_base] = records
    return records


def resolve_veupathdb_download_urls_from_id(source_id, species_key, timeout, headers):
    candidates = source_id_candidates("veupathdb", source_id, species_key)
    if len(candidates) == 0:
        candidates = [str(source_id or "").strip()]
    candidate_keys = [normalize_lookup_text(value) for value in candidates if normalize_lookup_text(value) != ""]
    if len(candidate_keys) == 0:
        raise ValueError("VEuPathDB id is empty")

    best_record = None
    best_rank = None
    for record in fetch_veupathdb_records(timeout, headers):
        attributes = record.get("attributes", {})
        genome_url = str(attributes.get("URLGenomeFasta", "") or "").strip()
        gff_url = str(attributes.get("URLgff", "") or "").strip()
        cds_url = infer_veupathdb_cds_url(attributes)
        derived_source_id = veupathdb_source_id_from_url(genome_url) or veupathdb_source_id_from_url(gff_url)
        normalized_source_id = normalize_lookup_text(derived_source_id)
        normalized_primary = normalize_lookup_text(attributes.get("primary_key", ""))
        normalized_species = normalize_lookup_text(attributes.get("species", ""))
        normalized_project = normalize_lookup_text(attributes.get("project_id", ""))
        source_id_exact = any(
            candidate == normalized_source_id for candidate in candidate_keys if normalized_source_id != ""
        )
        primary_exact = any(candidate == normalized_primary for candidate in candidate_keys if normalized_primary != "")
        species_exact = any(candidate == normalized_species for candidate in candidate_keys if normalized_species != "")
        project_exact = any(candidate == normalized_project for candidate in candidate_keys if normalized_project != "")
        source_id_contains = any(
            candidate in normalized_source_id or normalized_source_id in candidate
            for candidate in candidate_keys
            if normalized_source_id != ""
        )
        primary_contains = any(
            candidate in normalized_primary or normalized_primary in candidate
            for candidate in candidate_keys
            if normalized_primary != ""
        )
        species_contains = any(
            candidate in normalized_species or normalized_species in candidate
            for candidate in candidate_keys
            if normalized_species != ""
        )
        rank = (
            1 if genome_url != "" else 0,
            1 if gff_url != "" else 0,
            1 if cds_url != "" else 0,
            1 if source_id_exact else 0,
            1 if primary_exact else 0,
            1 if species_exact else 0,
            1 if project_exact else 0,
            1 if source_id_contains else 0,
            1 if primary_contains else 0,
            1 if species_contains else 0,
            -len(normalized_source_id),
            normalized_source_id,
        )
        if best_record is None or rank > best_rank:
            best_record = {
                "record": record,
                "derived_source_id": derived_source_id,
                "cds_url": cds_url,
            }
            best_rank = rank

    if best_record is None or max(best_rank[3:10]) == 0:
        raise ValueError("VEuPathDB id '{}' was not found in GenomeDataTypes".format(source_id))

    attributes = best_record["record"].get("attributes", {})
    genome_url = str(attributes.get("URLGenomeFasta", "") or "").strip()
    gff_url = str(attributes.get("URLgff", "") or "").strip()
    cds_url = str(best_record.get("cds_url", "") or "").strip()
    if gff_url == "" or genome_url == "":
        raise ValueError("VEuPathDB id '{}' did not resolve to GFF/genome URLs".format(source_id))

    inferred_species_key = parse_species_key_candidate(attributes.get("species", ""))
    if inferred_species_key == "":
        inferred_species_key = str(species_key or "").strip()
    if inferred_species_key == "":
        inferred_species_key = sanitize_identifier(best_record.get("derived_source_id", "") or str(source_id or ""))
    prefix = sanitize_identifier(inferred_species_key)
    source_token = best_record.get("derived_source_id", "") or strip_provider_prefix(source_id, "veupathdb")
    source_token = sanitize_identifier(source_token)

    return {
        "species_key": inferred_species_key,
        "cds_url": cds_url,
        "gff_url": gff_url,
        "genome_url": genome_url,
        "cds_filename": "{}.veupathdb.{}.cds.fa".format(prefix, source_token),
        "gff_filename": "{}.veupathdb.{}.gene.gff3".format(prefix, source_token),
        "genome_filename": "{}.veupathdb.{}.genome.fa".format(prefix, source_token),
    }


def resolve_insectbase_download_urls_from_id(source_id, species_key, timeout, headers):
    ibg_id = extract_insectbase_ibg_id_candidate(source_id)
    if ibg_id == "":
        raise ValueError("InsectBase id '{}' did not contain an IBG_* genome identifier".format(source_id))

    api_base = resolve_insectbase_api_base_url()
    detail_url = "{}/genomes/{}/".format(api_base, quote(ibg_id, safe=""))
    record = fetch_json_with_headers(detail_url, timeout, headers)
    species_name = str(record.get("species", "") or "").strip()
    if species_name == "":
        raise ValueError("InsectBase id '{}' did not resolve to a species name".format(source_id))

    species_token = re.sub(r"\s+", "_", species_name)
    inferred_species_key = str(species_key or "").strip()
    if inferred_species_key == "":
        inferred_species_key = parse_species_key_candidate(species_name)
    if inferred_species_key == "":
        inferred_species_key = species_token

    api_parts = urlparse(api_base)
    site_root = "{}://{}".format(api_parts.scheme, api_parts.netloc)
    base_data_url = "{}/data/genome/{}/{}".format(
        site_root.rstrip("/"),
        quote(species_token, safe="._-"),
        quote(species_token, safe="._-"),
    )

    return {
        "species_key": inferred_species_key,
        "cds_url": base_data_url + ".cds.fa",
        "gff_url": base_data_url + ".gff3",
        "genome_url": base_data_url + ".genome.fa.tar.bz2",
        "cds_filename": species_token + ".cds.fa",
        "gff_filename": species_token + ".gff3",
        "genome_filename": species_token + ".genome.fa.tar.bz2",
    }

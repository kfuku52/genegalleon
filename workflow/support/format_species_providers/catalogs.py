"""Provider resolver implementation: catalogs."""

import re
from pathlib import Path
from urllib.parse import quote, urljoin, urlparse

from format_species_annotations import (
    sanitize_identifier,
)
from format_species_common import (
    is_gff_filename,
    is_transient_network_error,
)
from format_species_constants import (
    CITRUSGENOMEDB_BINOMIAL_PATTERN,
    COGE_GID_PATTERN,
)
from format_species_manifest import manifest_has_usable_source_bundle
from format_species_provider_urls import (
    extract_ncbi_accession_from_source_id,
    fetch_text_with_headers,
    is_url_like,
    parse_links_from_document,
    resolve_citrusgenomedb_web_base_url,
    resolve_cngb_cnsa_base_url,
    resolve_coge_api_base_url,
    resolve_coge_web_base_url,
    resolve_plantaedb_web_base_url,
    strip_provider_prefix,
)

from .common import (
    fetch_json_with_headers,
    normalize_lookup_text,
    parse_species_key_candidate,
    remote_resource_exists,
    resolve_urls_from_index_url,
    source_id_candidates,
)
from .ncbi import (
    resolve_ncbi_download_urls_from_id,
)


def extract_coge_gid_candidate(source_id):
    candidates = source_id_candidates("coge", source_id, species_key="")
    for candidate in candidates:
        stripped = strip_provider_prefix(candidate, "coge")
        if COGE_GID_PATTERN.match(stripped):
            return stripped
    return ""


def pick_best_coge_genome(genomes, source_id):
    query = strip_provider_prefix(source_id, "coge")
    query_normalized = normalize_lookup_text(query)
    selected = None
    selected_rank = None
    for genome in genomes:
        if not isinstance(genome, dict):
            continue
        gid = str(genome.get("id", "") or "").strip()
        if not gid.isdigit():
            continue
        name = str(genome.get("name", "") or "")
        info = str(genome.get("info", "") or "")
        organism = str((genome.get("organism") or {}).get("name", "") or "")
        searchable = [name, info, organism]
        exact = any(query_normalized != "" and normalize_lookup_text(text) == query_normalized for text in searchable)
        contains = any(
            query_normalized != "" and query_normalized in normalize_lookup_text(text) for text in searchable
        )
        deleted = bool(genome.get("deleted"))
        certified = bool(genome.get("certified"))
        rank = (
            0 if deleted else 1,
            1 if exact else 0,
            1 if contains else 0,
            1 if certified else 0,
            int(gid),
        )
        if selected is None or rank > selected_rank:
            selected = genome
            selected_rank = rank
    return selected


def resolve_coge_genome_from_source_id(source_id, timeout, headers):
    gid = extract_coge_gid_candidate(source_id)
    if gid == "":
        raise ValueError("CoGe id must be genome_id (numeric gid). got '{}'".format(source_id))
    return gid, {"id": int(gid)}


def resolve_coge_download_urls_from_id(source_id, species_key, timeout, headers):
    gid, genome_info = resolve_coge_genome_from_source_id(source_id, timeout, headers)
    web_base = resolve_coge_web_base_url()
    api_base = resolve_coge_api_base_url()

    gff_meta_url = "{}/GenomeInfo.pl?fname=get_gff&gid={}&id_type=name&cds=0&annos=0&nu=0&upa=0&chr=".format(
        web_base, gid
    )
    gff_meta = fetch_json_with_headers(gff_meta_url, timeout, headers)
    gff_candidates = [str(item).strip() for item in gff_meta.get("files", []) if str(item).strip() != ""]
    gff_url = gff_candidates[0] if len(gff_candidates) > 0 else ""
    gff_filename = str(gff_meta.get("file", "") or "").strip()
    if gff_url == "" and gff_filename != "":
        gff_url = "{}/api/v1/downloads/?gid={}&filename={}".format(web_base, gid, quote(gff_filename, safe=""))
    if gff_url == "":
        raise ValueError("CoGe did not return a downloadable GFF URL for gid {}".format(gid))

    organism_name = str((genome_info.get("organism") or {}).get("name", "") or "")
    inferred_species = parse_species_key_candidate(organism_name)
    if inferred_species == "" and gff_filename != "":
        prefix = re.sub(r"[.]gid[0-9]+.*$", "", gff_filename, flags=re.IGNORECASE)
        inferred_species = parse_species_key_candidate(prefix.replace("_", " "))
    if inferred_species == "":
        inferred_species = species_key
    if inferred_species == "":
        inferred_species = "CoGe_gid{}".format(gid)
    prefix = sanitize_identifier(inferred_species)

    if gff_filename == "":
        gff_filename = "{}.coge.gid{}.gff".format(prefix, gid)
    elif not is_gff_filename(gff_filename):
        gff_filename = "{}.coge.gid{}.gff".format(prefix, gid)

    return {
        "species_key": inferred_species,
        "cds_url": "{}/get_seqs_for_feattype_for_genome.pl?ftid=3;dsgid={};".format(web_base, gid),
        "gff_url": gff_url,
        "genome_url": "{}/genomes/{}/sequence".format(api_base, gid),
        "cds_filename": "{}.coge.gid{}.cds.fa".format(prefix, gid),
        "gff_filename": gff_filename,
        "genome_filename": "{}.coge.gid{}.genome.fa".format(prefix, gid),
    }


def resolve_cngb_summary_from_id(cngb_id, timeout, headers):
    base = resolve_cngb_cnsa_base_url()
    summary_url = "{}/assembly/public_view/?q={}".format(base, quote(cngb_id, safe=""))
    payload = fetch_json_with_headers(summary_url, timeout, headers)
    if payload.get("code") != 0:
        error_message = str((payload.get("error") or {}).get("content", "") or "").strip()
        if error_message == "":
            error_message = "unknown CNGB API error"
        raise ValueError("CNGB summary lookup failed for '{}': {}".format(cngb_id, error_message))
    summary = payload.get("data", {}).get("summary_data", {})
    if not isinstance(summary, dict):
        raise ValueError("CNGB summary response was malformed for '{}'".format(cngb_id))
    return summary


def resolve_cngb_download_urls_from_id(source_id, timeout, headers):
    raw_candidates = source_id_candidates("cngb", source_id, species_key="")
    cngb_id = strip_provider_prefix(source_id, "cngb")
    if cngb_id == "":
        cngb_id = source_id

    accession = ""
    for candidate in raw_candidates:
        accession = extract_ncbi_accession_from_source_id(candidate)
        if accession != "":
            break

    summary = {}
    if accession == "":
        summary = resolve_cngb_summary_from_id(cngb_id, timeout, headers)
        for key in (
            "refseq_assembly_accession",
            "genbank_assembly_accession",
            "accession_id",
            "external_accession_id",
        ):
            accession = extract_ncbi_accession_from_source_id(summary.get(key, ""))
            if accession != "":
                break

    if accession == "":
        raise ValueError("CNGB id '{}' did not map to a downloadable NCBI assembly accession".format(source_id))

    resolved = resolve_ncbi_download_urls_from_id(accession, timeout, ncbi_source="auto")
    if resolved.get("species_key", "") == "":
        organism_name = str((summary.get("organism") or {}).get("name", "") or "")
        species_candidate = parse_species_key_candidate(organism_name)
        if species_candidate != "":
            resolved["species_key"] = species_candidate
    return resolved


def resolve_plantaedb_source_page_url(source_id):
    source_clean = str(source_id or "").strip()
    stripped = strip_provider_prefix(source_clean, "plantaedb").strip()
    if stripped == "":
        raise ValueError("PlantaeDB id is empty")
    if is_url_like(stripped):
        return stripped
    base = resolve_plantaedb_web_base_url().rstrip("/") + "/"
    if stripped.startswith("/"):
        stripped = stripped[1:]
    return urljoin(base, stripped)


def extract_plantaedb_ncbi_accession(page_url, page_text):
    links = parse_links_from_document(page_url, page_text)
    preferred_links = []
    fallback_links = []
    for link in links:
        if "ncbi.nlm.nih.gov" not in link.lower():
            continue
        accession = extract_ncbi_accession_from_source_id(link)
        if accession == "":
            continue
        if "/data-hub/genome/" in link or "/datasets/genome/" in link:
            preferred_links.append(accession)
        else:
            fallback_links.append(accession)
    for accession in preferred_links + fallback_links:
        if accession != "":
            return accession
    return extract_ncbi_accession_from_source_id(page_text)


def resolve_plantaedb_download_urls_from_id(source_id, species_key, timeout, headers):
    page_url = resolve_plantaedb_source_page_url(source_id)
    page_text = fetch_text_with_headers(page_url, timeout, headers)
    accession = extract_plantaedb_ncbi_accession(page_url, page_text)
    if accession == "":
        raise ValueError("PlantaeDB page did not expose an NCBI assembly accession: {}".format(page_url))
    resolved = resolve_ncbi_download_urls_from_id(accession, timeout, ncbi_source="auto")
    if str(resolved.get("species_key", "") or "").strip() == "" and species_key != "":
        resolved["species_key"] = species_key
    return resolved


def normalize_citrusgenomedb_source_page_url(source_id):
    source_clean = strip_provider_prefix(source_id, "citrusgenomedb").strip()
    if source_clean == "":
        raise ValueError("Citrus Genome Database id is empty")
    base = resolve_citrusgenomedb_web_base_url().rstrip("/") + "/"
    if is_url_like(source_clean):
        return source_clean
    lowered = source_clean.lower()
    if source_clean.startswith("/"):
        return urljoin(base, source_clean.lstrip("/"))
    if lowered.startswith("analysis/") or lowered.startswith("organism/"):
        return urljoin(base, source_clean)
    if source_clean.isdigit():
        return urljoin(base, "Analysis/{}".format(source_clean))
    return urljoin(base, source_clean)


def infer_citrusgenomedb_species_key(page_text, links, fallback):
    text = str(page_text or "")
    match = CITRUSGENOMEDB_BINOMIAL_PATTERN.search(text)
    if match is not None:
        genus = match.group(1).strip()
        epithet = re.sub(r"\s+", "_", match.group(2).strip())
        candidate = "{}_{}".format(genus, epithet)
        if candidate != "":
            return candidate
    for link in links:
        parts = [part for part in urlparse(link).path.split("/") if part]
        if "citrus_downloads" not in parts:
            continue
        try:
            species_part = parts[parts.index("citrus_downloads") + 1]
        except Exception:
            continue
        species_key = sanitize_identifier(species_part)
        if species_key != "":
            return species_key
    return str(fallback or "").strip()


def merge_resolved_manifest_bundle_fields(base_resolved, candidate_resolved):
    merged = dict(base_resolved or {})
    candidate = dict(candidate_resolved or {})
    for key in ("cds_url", "gff_url", "genome_url"):
        if str(merged.get(key, "") or "").strip() == "" and str(candidate.get(key, "") or "").strip() != "":
            merged[key] = candidate[key]
    return merged


def resolved_manifest_bundle_urls_available(resolved, timeout, headers):
    urls = [
        str(resolved.get(key, "") or "").strip()
        for key in ("cds_url", "gff_url", "gbff_url", "genome_url")
        if str(resolved.get(key, "") or "").strip() != ""
    ]
    for url in urls:
        try:
            if not remote_resource_exists(url, timeout, headers):
                return False
        except Exception as exc:
            if is_transient_network_error(exc):
                continue
            raise
    return True


def citrusgenomedb_repository_request_headers(headers):
    req_headers = dict(headers or {})
    if "User-Agent" not in req_headers:
        req_headers["User-Agent"] = "Mozilla/5.0 (compatible; genegalleon-input-generation)"
    return req_headers


def resolve_citrusgenomedb_repository_bundle(index_url, timeout, headers):
    normalized_index_url = str(index_url).rstrip("/") + "/"
    repo_headers = citrusgenomedb_repository_request_headers(headers)
    resolved = resolve_urls_from_index_url("citrusgenomedb", index_url, timeout, repo_headers)
    if manifest_has_usable_source_bundle(
        resolved.get("cds_url", ""),
        resolved.get("gff_url", ""),
        "",
        resolved.get("genome_url", ""),
    ):
        return resolved
    page_text = fetch_text_with_headers(index_url, timeout, repo_headers)
    links = parse_links_from_document(index_url, page_text)
    merged = dict(resolved)
    for link in links:
        lower = link.lower()
        if "/citrus_downloads/" not in lower:
            continue
        if link.rstrip("/") == str(index_url).rstrip("/"):
            continue
        if not link.startswith(normalized_index_url):
            continue
        candidate = resolve_urls_from_index_url("citrusgenomedb", link, timeout, repo_headers)
        merged = merge_resolved_manifest_bundle_fields(merged, candidate)
    return merged


def resolve_citrusgenomedb_bundle_from_page(page_url, species_key, timeout, headers):
    page_text = fetch_text_with_headers(page_url, timeout, headers)
    if "Bulk download of the assembly files will become available once the data is published" in page_text:
        return None
    links = parse_links_from_document(page_url, page_text)
    resolved = resolve_urls_from_index_url("citrusgenomedb", page_url, timeout, headers)
    if not manifest_has_usable_source_bundle(
        resolved.get("cds_url", ""),
        resolved.get("gff_url", ""),
        "",
        resolved.get("genome_url", ""),
    ):
        for link in links:
            if "/citrus_downloads/" not in link.lower():
                continue
            candidate = resolve_citrusgenomedb_repository_bundle(link, timeout, headers)
            resolved = merge_resolved_manifest_bundle_fields(resolved, candidate)
    if not manifest_has_usable_source_bundle(
        resolved.get("cds_url", ""),
        resolved.get("gff_url", ""),
        "",
        resolved.get("genome_url", ""),
    ):
        return None
    repo_headers = citrusgenomedb_repository_request_headers(headers)
    if not resolved_manifest_bundle_urls_available(resolved, timeout, repo_headers):
        return None
    resolved["species_key"] = str(species_key or "").strip() or infer_citrusgenomedb_species_key(
        page_text, links, species_key
    )
    resolved["cds_filename"] = Path(urlparse(resolved["cds_url"]).path).name if resolved.get("cds_url", "") else ""
    resolved["gff_filename"] = Path(urlparse(resolved["gff_url"]).path).name if resolved.get("gff_url", "") else ""
    resolved["genome_filename"] = (
        Path(urlparse(resolved["genome_url"]).path).name if resolved.get("genome_url", "") else ""
    )
    return resolved


def resolve_citrusgenomedb_download_urls_from_id(source_id, species_key, timeout, headers):
    page_url = normalize_citrusgenomedb_source_page_url(source_id)
    if "/analysis/" in urlparse(page_url).path.lower():
        resolved = resolve_citrusgenomedb_bundle_from_page(page_url, species_key, timeout, headers)
        if resolved is None:
            raise ValueError(
                "Citrus Genome Database analysis page '{}' did not expose downloadable genome/GFF bundle".format(
                    page_url
                )
            )
        return resolved

    page_text = fetch_text_with_headers(page_url, timeout, headers)
    links = parse_links_from_document(page_url, page_text)

    best = None
    best_rank = None
    direct_bundle = resolve_citrusgenomedb_bundle_from_page(page_url, species_key, timeout, headers)
    if direct_bundle is not None:
        genome_url = str(direct_bundle.get("genome_url", "") or "").strip().lower()
        best = direct_bundle
        best_rank = (
            1 if "/assembly/" in genome_url else 0,
            1 if direct_bundle.get("gff_url", "") and direct_bundle.get("genome_url", "") else 0,
            1 if direct_bundle.get("cds_url", "") else 0,
            0,
        )

    seen = set()
    for order, link in enumerate(links, start=1):
        if "/analysis/" not in link.lower():
            continue
        if link in seen:
            continue
        seen.add(link)
        try:
            resolved = resolve_citrusgenomedb_bundle_from_page(link, species_key, timeout, headers)
        except Exception:
            continue
        if resolved is None:
            continue
        genome_url = str(resolved.get("genome_url", "") or "").strip().lower()
        rank = (
            1 if "/assembly/" in genome_url else 0,
            1 if resolved.get("gff_url", "") and resolved.get("genome_url", "") else 0,
            1 if resolved.get("cds_url", "") else 0,
            -order,
        )
        if best is None or rank > best_rank:
            best = resolved
            best_rank = rank
    if best is None:
        raise ValueError(
            "Citrus Genome Database id '{}' did not resolve to downloadable GFF/genome URLs".format(source_id)
        )
    return best

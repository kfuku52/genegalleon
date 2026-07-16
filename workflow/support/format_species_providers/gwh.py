"""Provider resolver implementation: gwh."""

import html
import re
from pathlib import Path
from urllib.parse import unquote, urlencode, urljoin, urlparse

from format_species_annotations import (
    sanitize_identifier,
)
from format_species_constants import (
    GWH_ASSEMBLY_ACCESSION_PATTERN,
    GWH_ASSEMBLY_ACCESSION_SEARCH_PATTERN,
)
from format_species_manifest import manifest_has_usable_source_bundle
from format_species_provider_urls import (
    fetch_text_with_headers,
    is_url_like,
    parse_links_from_document,
    resolve_gwh_download_base_url,
    resolve_gwh_search_api_url,
    resolve_gwh_web_base_url,
    strip_provider_prefix,
)

from .common import (
    fetch_json_with_headers,
    parse_species_key_candidate,
    resolve_urls_from_index_url,
    source_id_candidates,
)


def extract_gwh_accession_candidate(source_id):
    candidates = source_id_candidates("gwh", source_id, species_key="")
    for candidate in candidates:
        stripped = strip_provider_prefix(candidate, "gwh")
        match = GWH_ASSEMBLY_ACCESSION_SEARCH_PATTERN.search(stripped)
        if match is None:
            continue
        accession = match.group(1).upper()
        if GWH_ASSEMBLY_ACCESSION_PATTERN.match(accession):
            return accession
    return ""


def gwh_accession_stem(accession):
    return str(accession or "").strip().upper().split(".", 1)[0]


def parse_last_path_token(url):
    path = urlparse(str(url or "")).path.rstrip("/")
    if path == "":
        return ""
    return unquote(path.split("/")[-1])


def strip_html_to_text(text):
    unescaped = html.unescape(str(text or ""))
    unescaped = re.sub(r"<[^>]+>", " ", unescaped)
    return re.sub(r"\s+", " ", unescaped).strip()


def infer_gwh_species_key_from_show_page(page_text, fallback):
    clean = strip_html_to_text(page_text)
    match = re.search(r"Scientific Name\s+(.+?)\s+Common Names", clean, flags=re.IGNORECASE)
    if match is not None:
        inferred = parse_species_key_candidate(match.group(1))
        if inferred != "":
            return inferred
    return str(fallback or "").strip()


def resolve_gwh_bundle_from_show_page(page_url, species_key, timeout, headers):
    page_text = fetch_text_with_headers(page_url, timeout, headers)
    resolved = resolve_urls_from_index_url("gwh", page_url, timeout, headers)
    if not manifest_has_usable_source_bundle(
        resolved.get("cds_url", ""),
        resolved.get("gff_url", ""),
        "",
        resolved.get("genome_url", ""),
    ):
        raise ValueError("GWH show page '{}' did not expose a usable bundle".format(page_url))
    resolved["species_key"] = str(species_key or "").strip() or infer_gwh_species_key_from_show_page(
        page_text, species_key
    )
    resolved["cds_filename"] = Path(urlparse(resolved["cds_url"]).path).name if resolved.get("cds_url", "") else ""
    resolved["gff_filename"] = Path(urlparse(resolved["gff_url"]).path).name if resolved.get("gff_url", "") else ""
    resolved["genome_filename"] = (
        Path(urlparse(resolved["genome_url"]).path).name if resolved.get("genome_url", "") else ""
    )
    return resolved


def resolve_gwh_show_page_url_from_id(source_id, timeout, headers):
    source_clean = str(source_id or "").strip()
    accession = extract_gwh_accession_candidate(source_clean)
    if accession == "":
        raise ValueError("GWH id '{}' did not contain a GWH accession".format(source_id))
    if (
        is_url_like(source_clean)
        and "/gwh/" in source_clean.lower()
        and source_clean.rstrip("/").lower().endswith("/show")
    ):
        return source_clean, accession

    search_url = resolve_gwh_search_api_url()
    query_url = search_url + "?" + urlencode({"term": gwh_accession_stem(accession)})
    payload = fetch_json_with_headers(query_url, timeout, headers)
    hits = payload.get("data", [])
    if not isinstance(hits, list) or len(hits) == 0:
        raise ValueError("GWH accession '{}' returned no search hits from {}".format(accession, search_url))

    normalized_accession = accession.lower()
    accession_stem = gwh_accession_stem(accession).lower()
    candidates = []
    for hit in hits:
        url = str(hit.get("url", "") or "").strip()
        if url == "":
            continue
        if not is_url_like(url):
            url = urljoin(resolve_gwh_web_base_url().rstrip("/") + "/", url.lstrip("/"))
        attrs = hit.get("attrs", {}) or {}
        url_lower = url.lower()
        score = (
            1 if str(attrs.get("has_annotation", "")).lower() == "yes" else 0,
            1 if str(attrs.get("source", "")).lower() == "direct submission" else 0,
            1 if "/assembly/" in url_lower else 0,
            1 if "/ncbi_assembly/" not in url_lower else 0,
            url_lower,
        )
        candidates.append((score, url))
    if len(candidates) == 0:
        raise ValueError("GWH accession '{}' returned search hits without usable show page URLs".format(accession))

    last_error = None
    for _, page_url in sorted(candidates, reverse=True):
        try:
            resolved = resolve_gwh_bundle_from_show_page(page_url, species_key="", timeout=timeout, headers=headers)
        except Exception as exc:
            last_error = exc
            continue
        bundle_urls = " ".join(
            [
                str(resolved.get("cds_url", "") or ""),
                str(resolved.get("gff_url", "") or ""),
                str(resolved.get("genome_url", "") or ""),
            ]
        ).lower()
        if normalized_accession in bundle_urls or accession_stem in bundle_urls:
            return page_url, accession
    if last_error is not None:
        raise ValueError("GWH accession '{}' did not resolve via show-page fallback: {}".format(accession, last_error))
    raise ValueError("GWH accession '{}' did not match any show-page download bundle".format(accession))


def resolve_gwh_folder_url_from_id(source_id, timeout, headers):
    source_clean = str(source_id or "").strip()
    accession = extract_gwh_accession_candidate(source_clean)
    if accession == "":
        raise ValueError("GWH id '{}' did not contain a GWH accession".format(source_id))

    if is_url_like(source_clean) and "/gwh/" in source_clean.lower():
        parsed = urlparse(source_clean)
        if parsed.path.endswith("/"):
            return source_clean, accession
        return parsed._replace(path=parsed.path.rsplit("/", 1)[0] + "/").geturl(), accession

    root_url = resolve_gwh_download_base_url().rstrip("/") + "/"
    root_text = fetch_text_with_headers(root_url, timeout, headers)
    category_urls = []
    for link in parse_links_from_document(root_url, root_text):
        if not urlparse(link).path.endswith("/"):
            continue
        token = parse_last_path_token(link)
        if token == "":
            continue
        category_urls.append(link)
    if len(category_urls) == 0:
        raise ValueError("GWH root '{}' did not expose category directories".format(root_url))

    accession_stem = gwh_accession_stem(accession)
    best_url = ""
    best_rank = None
    last_error = None
    for category_url in sorted(
        set(category_urls),
        key=lambda url: (
            0 if parse_last_path_token(url).lower() == "plants" else 1,
            parse_last_path_token(url).lower(),
        ),
    ):
        try:
            category_text = fetch_text_with_headers(category_url, timeout, headers)
        except Exception as exc:
            last_error = exc
            continue
        for link in parse_links_from_document(category_url, category_text):
            if not urlparse(link).path.endswith("/"):
                continue
            folder_name = parse_last_path_token(link)
            if folder_name == "":
                continue
            folder_lower = folder_name.lower()
            has_full = accession.lower() in folder_lower
            has_stem = accession_stem.lower() in folder_lower
            if not has_full and not has_stem:
                continue
            rank = (
                1 if has_full else 0,
                1 if folder_lower.endswith(accession.lower()) else 0,
                1 if parse_last_path_token(category_url).lower() == "plants" else 0,
                -len(folder_name),
                folder_lower,
            )
            if best_url == "" or rank > best_rank:
                best_url = link
                best_rank = rank
    if best_url != "":
        return best_url, accession
    if last_error is not None:
        raise ValueError("GWH accession '{}' did not resolve: {}".format(accession, last_error))
    raise ValueError("GWH accession '{}' was not found under {}".format(accession, root_url))


def resolve_gwh_download_urls_from_id(source_id, species_key, timeout, headers):
    source_clean = str(source_id or "").strip()
    accession = extract_gwh_accession_candidate(source_clean)
    errors = []

    if (
        is_url_like(source_clean)
        and "/gwh/" in source_clean.lower()
        and source_clean.rstrip("/").lower().endswith("/show")
    ):
        try:
            return resolve_gwh_bundle_from_show_page(source_clean, species_key, timeout, headers)
        except Exception as exc:
            errors.append("show page: {}".format(exc))

    try:
        folder_url, accession = resolve_gwh_folder_url_from_id(source_id, timeout, headers)
        resolved = resolve_urls_from_index_url("gwh", folder_url, timeout, headers)
        if resolved.get("gff_url", "") == "" or resolved.get("genome_url", "") == "":
            raise ValueError(
                "GWH accession '{}' did not resolve to downloadable GFF/genome URLs from {}".format(
                    accession, folder_url
                )
            )

        folder_name = parse_last_path_token(folder_url)
        inferred_species_key = str(species_key or "").strip()
        if inferred_species_key == "":
            inferred_species_key = parse_species_key_candidate(folder_name.replace("_", " "))
        if inferred_species_key == "":
            inferred_species_key = sanitize_identifier(accession)
        return {
            "species_key": inferred_species_key,
            "cds_url": resolved["cds_url"],
            "gff_url": resolved["gff_url"],
            "genome_url": resolved["genome_url"],
            "cds_filename": Path(urlparse(resolved["cds_url"]).path).name,
            "gff_filename": Path(urlparse(resolved["gff_url"]).path).name,
            "genome_filename": Path(urlparse(resolved["genome_url"]).path).name,
        }
    except Exception as exc:
        errors.append("public index: {}".format(exc))

    try:
        page_url, accession = resolve_gwh_show_page_url_from_id(source_id, timeout, headers)
        resolved = resolve_gwh_bundle_from_show_page(page_url, species_key, timeout, headers)
        if str(resolved.get("species_key", "")).strip() == "":
            resolved["species_key"] = sanitize_identifier(accession)
        return resolved
    except Exception as exc:
        errors.append("show fallback: {}".format(exc))

    raise ValueError(
        "GWH accession '{}' did not resolve via public index or show-page fallback ({})".format(
            accession or source_id,
            "; ".join(errors),
        )
    )

"""Provider resolver implementation: ddbj."""

from urllib.parse import quote

from format_species_provider_urls import (
    resolve_ddbj_public_wgs_base_url,
    resolve_ddbj_search_api_base_url,
)

from .common import (
    extract_ddbj_bioproject_id,
    extract_ddbj_wgs_prefix,
    fetch_json_with_headers,
    parse_species_key_candidate,
    remote_resource_exists,
)


def ddbj_entry_organism_name(entry):
    organism = entry.get("organism")
    if isinstance(organism, dict):
        return str(organism.get("name") or "").strip()
    return str(organism or "").strip()


def ddbj_wgs_candidate_urls(prefix):
    normalized_prefix = str(prefix or "").strip().upper()
    if normalized_prefix == "":
        return []
    base = resolve_ddbj_public_wgs_base_url().rstrip("/")
    candidates = []

    def add(url):
        if url not in candidates:
            candidates.append(url)

    if len(normalized_prefix) >= 4:
        add(
            "{}/{}/{}/{}.gz".format(
                base,
                quote(normalized_prefix[:2], safe=""),
                quote(normalized_prefix[2:4], safe=""),
                quote(normalized_prefix, safe=""),
            )
        )
    if len(normalized_prefix) >= 2:
        add(
            "{}/{}/{}.gz".format(
                base,
                quote(normalized_prefix[:2], safe=""),
                quote(normalized_prefix, safe=""),
            )
        )
    return candidates


def collect_ddbj_wgs_prefixes_from_bioproject_entry(entry):
    prefixes = []
    for item in entry.get("dbXrefs", ()):
        if not isinstance(item, dict):
            continue
        if str(item.get("type") or "").strip().lower() != "insdc-master":
            continue
        prefix = extract_ddbj_wgs_prefix(str(item.get("identifier") or "").strip())
        if prefix == "" or prefix in prefixes:
            continue
        prefixes.append(prefix)
    return prefixes


def resolve_ddbj_download_urls_from_id(source_id, species_key, timeout, headers):
    inferred_species_key = str(species_key or "").strip()
    prefix_candidates = []
    bioproject_id = extract_ddbj_bioproject_id(source_id)
    entry_url = ""
    if bioproject_id != "":
        entry_url = "{}/entries/bioproject/{}".format(
            resolve_ddbj_search_api_base_url(),
            quote(bioproject_id, safe=""),
        )
        entry = fetch_json_with_headers(entry_url, timeout, headers)
        if inferred_species_key == "":
            inferred_species_key = parse_species_key_candidate(ddbj_entry_organism_name(entry))
        prefix_candidates.extend(collect_ddbj_wgs_prefixes_from_bioproject_entry(entry))
        if len(prefix_candidates) == 0:
            dbxrefs_entry = fetch_json_with_headers(entry_url + "/dbxrefs.json", timeout, headers)
            prefix_candidates.extend(collect_ddbj_wgs_prefixes_from_bioproject_entry(dbxrefs_entry))

    direct_prefix = extract_ddbj_wgs_prefix(source_id)
    if direct_prefix != "" and direct_prefix not in prefix_candidates:
        prefix_candidates.append(direct_prefix)
    if len(prefix_candidates) == 0:
        raise ValueError(
            "DDBJ id '{}' did not resolve to a public WGS BioProject or master accession".format(source_id)
        )

    attempted_urls = []
    for prefix in prefix_candidates:
        for candidate_url in ddbj_wgs_candidate_urls(prefix):
            attempted_urls.append(candidate_url)
            if remote_resource_exists(candidate_url, timeout, headers):
                return {
                    "species_key": inferred_species_key,
                    "gbff_url": candidate_url,
                    "gbff_filename": "{}.gbff.gz".format(prefix),
                }

    raise ValueError(
        "DDBJ id '{}' resolved to WGS prefix(es) {} but no public GBFF was found at {}".format(
            source_id,
            ",".join(prefix_candidates),
            ", ".join(attempted_urls),
        )
    )

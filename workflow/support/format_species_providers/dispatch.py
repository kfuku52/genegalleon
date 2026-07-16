"""Provider resolver implementation: dispatch."""

import os

from format_species_constants import (
    DOWNLOAD_LABELS,
)
from format_species_manifest import manifest_has_usable_source_bundle
from format_species_provider_config import (
    DOWNLOAD_MANIFEST_SUPPORTED_PROVIDERS,
    ORYZA_MINUTA_PROVIDER,
)
from format_species_provider_urls import (
    ENSEMBLGENOMES_DEFAULT_ID_URL_TEMPLATES,
    PROVIDER_DEFAULT_ID_PAGE_URL_TEMPLATES,
    ensembl_like_default_id_url_templates,
    is_url_like,
    provider_env_prefix,
    render_id_url_template,
)

from .catalogs import (
    resolve_citrusgenomedb_download_urls_from_id,
    resolve_cngb_download_urls_from_id,
    resolve_coge_download_urls_from_id,
    resolve_plantaedb_download_urls_from_id,
)
from .common import (
    expand_ensemblgenomes_id_candidates,
    resolve_urls_from_index_url,
    source_id_candidates,
)
from .ddbj import (
    resolve_ddbj_download_urls_from_id,
)
from .fernbase_oryza import (
    resolve_fernbase_download_urls_from_id,
    resolve_oryza_minuta_download_urls_from_id,
)
from .gwh import (
    resolve_gwh_download_urls_from_id,
)
from .web import (
    resolve_figshare_download_urls_from_id,
    resolve_insectbase_download_urls_from_id,
    resolve_plantgarden_download_urls_from_id,
    resolve_veupathdb_download_urls_from_id,
)


class ProviderIdResolver:
    """Adapt provider-specific functions to one manifest-resolution interface."""

    def __init__(self, handler, *, accepts_species_key=True, accepts_row=False):
        self.handler = handler
        self.accepts_species_key = accepts_species_key
        self.accepts_row = accepts_row

    def resolve(self, *, source_id, species_key, timeout, headers, row=None):
        if self.accepts_row:
            return self.handler(source_id, species_key, timeout, headers, row=row)
        if self.accepts_species_key:
            return self.handler(source_id, species_key, timeout, headers)
        return self.handler(source_id, timeout, headers)


PROVIDER_ID_RESOLVERS = {
    "coge": ProviderIdResolver(resolve_coge_download_urls_from_id),
    "cngb": ProviderIdResolver(resolve_cngb_download_urls_from_id, accepts_species_key=False),
    "ddbj": ProviderIdResolver(resolve_ddbj_download_urls_from_id),
    "gwh": ProviderIdResolver(resolve_gwh_download_urls_from_id),
    "citrusgenomedb": ProviderIdResolver(resolve_citrusgenomedb_download_urls_from_id),
    "figshare": ProviderIdResolver(resolve_figshare_download_urls_from_id, accepts_row=True),
    "plantgarden": ProviderIdResolver(resolve_plantgarden_download_urls_from_id),
    "plantaedb": ProviderIdResolver(resolve_plantaedb_download_urls_from_id),
    "fernbase": ProviderIdResolver(resolve_fernbase_download_urls_from_id),
    "veupathdb": ProviderIdResolver(resolve_veupathdb_download_urls_from_id),
    "insectbase": ProviderIdResolver(resolve_insectbase_download_urls_from_id),
    ORYZA_MINUTA_PROVIDER: ProviderIdResolver(resolve_oryza_minuta_download_urls_from_id),
}


def resolve_provider_specific_download_urls_from_id(provider, source_id, species_key, timeout, headers, row=None):
    resolver = PROVIDER_ID_RESOLVERS.get(provider)
    if resolver is None:
        return None
    return resolver.resolve(
        source_id=source_id,
        species_key=species_key,
        timeout=timeout,
        headers=headers,
        row=row,
    )


def resolve_non_ncbi_download_urls_from_id(provider, source_id, species_key, timeout, headers, row=None):
    if provider not in DOWNLOAD_MANIFEST_SUPPORTED_PROVIDERS:
        raise ValueError(
            "provider '{}' is not supported for --download-manifest (use --input-dir for local formatting)".format(
                provider
            )
        )

    resolved = {}
    provider_specific_error = None

    env_prefix = provider_env_prefix(provider)
    label_to_key = {"CDS": "cds_url", "GFF": "gff_url", "GENOME": "genome_url"}
    id_candidates = source_id_candidates(provider, source_id, species_key)
    if len(id_candidates) == 0:
        id_candidates = [str(source_id or "").strip()]
    effective_candidates = id_candidates
    if provider in ENSEMBLGENOMES_DEFAULT_ID_URL_TEMPLATES:
        try:
            effective_candidates = expand_ensemblgenomes_id_candidates(provider, id_candidates, timeout, headers)
        except Exception:
            effective_candidates = id_candidates

    for label in DOWNLOAD_LABELS:
        template_env = "{}_{}_URL_TEMPLATE".format(env_prefix, label)
        template = os.environ.get(template_env, "").strip()
        if template == "":
            template = ensembl_like_default_id_url_templates(provider).get(label, "")
        if template == "":
            continue
        for candidate in effective_candidates:
            url = render_id_url_template(template, provider, candidate, species_key)
            if url.endswith("/"):
                try:
                    discovered = resolve_urls_from_index_url(provider, url, timeout, headers)
                    inferred_url = discovered.get(label_to_key[label], "")
                    if inferred_url != "":
                        resolved[label_to_key[label]] = inferred_url
                        break
                except Exception:
                    continue
            else:
                resolved[label_to_key[label]] = url
                break

    if manifest_has_usable_source_bundle(
        resolved.get("cds_url", ""),
        resolved.get("gff_url", ""),
        resolved.get("gbff_url", ""),
        resolved.get("genome_url", ""),
    ):
        return resolved

    try:
        provider_specific = resolve_provider_specific_download_urls_from_id(
            provider=provider,
            source_id=source_id,
            species_key=species_key,
            timeout=timeout,
            headers=headers,
            row=row,
        )
    except Exception as exc:
        provider_specific = None
        provider_specific_error = exc
    if provider_specific is not None:
        for key in (
            "cds_url",
            "gff_url",
            "gbff_url",
            "genome_url",
            "cds_archive_member",
            "gff_archive_member",
            "gbff_archive_member",
            "genome_archive_member",
            "cds_filename",
            "gff_filename",
            "gbff_filename",
            "genome_filename",
            "species_key",
        ):
            if resolved.get(key, "") == "" and str(provider_specific.get(key, "") or "").strip() != "":
                resolved[key] = str(provider_specific[key]).strip()
        if manifest_has_usable_source_bundle(
            resolved.get("cds_url", ""),
            resolved.get("gff_url", ""),
            resolved.get("gbff_url", ""),
            resolved.get("genome_url", ""),
        ):
            return resolved

    source_clean = str(source_id or "").strip()
    if is_url_like(source_clean):
        discovered = resolve_urls_from_index_url(provider, source_clean, timeout, headers)
        for key in ("cds_url", "gff_url", "gbff_url", "genome_url"):
            if resolved.get(key, "") == "" and discovered.get(key, "") != "":
                resolved[key] = discovered[key]

    id_page_template_env = "{}_ID_URL_TEMPLATE".format(env_prefix)
    id_page_template = os.environ.get(id_page_template_env, "").strip()
    page_templates = []
    if id_page_template != "":
        page_templates.append(id_page_template)
    page_templates.extend(PROVIDER_DEFAULT_ID_PAGE_URL_TEMPLATES.get(provider, ()))

    for page_template in page_templates:
        for candidate in id_candidates:
            page_url = render_id_url_template(page_template, provider, candidate, species_key)
            try:
                discovered = resolve_urls_from_index_url(provider, page_url, timeout, headers)
            except Exception:
                continue
            for key in ("cds_url", "gff_url", "gbff_url", "genome_url"):
                if resolved.get(key, "") == "" and discovered.get(key, "") != "":
                    resolved[key] = discovered[key]
            if manifest_has_usable_source_bundle(
                resolved.get("cds_url", ""),
                resolved.get("gff_url", ""),
                resolved.get("gbff_url", ""),
                resolved.get("genome_url", ""),
            ):
                break
        if manifest_has_usable_source_bundle(
            resolved.get("cds_url", ""),
            resolved.get("gff_url", ""),
            resolved.get("gbff_url", ""),
            resolved.get("genome_url", ""),
        ):
            break

    if not manifest_has_usable_source_bundle(
        resolved.get("cds_url", ""),
        resolved.get("gff_url", ""),
        resolved.get("gbff_url", ""),
        resolved.get("genome_url", ""),
    ):
        detail = ""
        if provider_specific_error is not None:
            detail = " (provider-specific resolver: {})".format(provider_specific_error)
        raise ValueError(
            "could not infer a usable source bundle (CDS, GBFF, or GFF plus genome) from id '{}' for provider '{}'{}".format(
                source_id, provider, detail
            )
        )
    return resolved

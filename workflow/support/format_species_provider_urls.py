#!/usr/bin/env python3
"""Provider URL, HTTP, and source-id helpers for format_species_inputs.py."""

import os
import re
from http.cookiejar import CookieJar
from urllib.parse import urlencode, urljoin, urlparse
from urllib.request import HTTPCookieProcessor, Request, build_opener

from format_species_network import guarded_urlopen as urlopen

NCBI_ASSEMBLY_ACCESSION_PATTERN = re.compile(r"^GC[AF]_[0-9]+\.[0-9]+$", re.IGNORECASE)

DEFAULT_NCBI_EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DEFAULT_NCBI_FTP_BASE_URL = "https://ftp.ncbi.nlm.nih.gov"
DEFAULT_NCBI_DATASETS_BASE_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2"
DEFAULT_DDBJ_SEARCH_API_BASE_URL = "https://ddbj.nig.ac.jp/search/api"
DEFAULT_DDBJ_PUBLIC_WGS_BASE_URL = "https://ddbj.nig.ac.jp/public/ddbj_database/wgs"
DEFAULT_COGE_API_BASE_URL = "https://genomevolution.org/coge/api/v1"
DEFAULT_COGE_WEB_BASE_URL = "https://genomevolution.org/coge"
DEFAULT_CNGB_CNSA_BASE_URL = "https://db.cngb.org/cnsa/ajax"
DEFAULT_GWH_DOWNLOAD_BASE_URL = "https://download.cncb.ac.cn/gwh"
DEFAULT_GWH_WEB_BASE_URL = "https://ngdc.cncb.ac.cn/gwh"
DEFAULT_CITRUSGENOMEDB_WEB_BASE_URL = "https://www.citrusgenomedb.org"
DEFAULT_FIGSHARE_API_BASE_URL = "https://api.figshare.com/v2"
DEFAULT_PLANTGARDEN_WEB_BASE_URL = "https://plantgarden.jp"
DEFAULT_JGI_SIGNON_BASE_URL = "https://signon.jgi.doe.gov"
DEFAULT_PLANTAEDB_WEB_BASE_URL = "https://plantaedb.com"
DEFAULT_VEUPATHDB_SERVICE_BASE_URL = "https://veupathdb.org/veupathdb/service"
DEFAULT_INSECTBASE_API_BASE_URL = "https://www.insect-genome.com/api/genome"
DEFAULT_ORYZA_MINUTA_GRAMENE_BASE_URL = "https://ftp.gramene.org/oryza/tetraploids"

ENSEMBL_DEFAULT_ID_URL_TEMPLATES = {
    "CDS": "https://ftp.ensembl.org/pub/current_fasta/{id_lower}/cds/",
    "GFF": "https://ftp.ensembl.org/pub/current_gff3/{id_lower}/",
    "GENOME": "https://ftp.ensembl.org/pub/current_fasta/{id_lower}/dna/",
}
ENSEMBLPLANTS_DEFAULT_ID_URL_TEMPLATES = {
    "CDS": "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/{id_lower}/cds/",
    "GFF": "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/gff3/{id_lower}/",
    "GENOME": "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/{id_lower}/dna/",
}
ENSEMBLMETAZOA_DEFAULT_ID_URL_TEMPLATES = {
    "CDS": "https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/metazoa/fasta/{id_lower}/cds/",
    "GFF": "https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/metazoa/gff3/{id_lower}/",
    "GENOME": "https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/metazoa/fasta/{id_lower}/dna/",
}
ENSEMBLPROTISTS_DEFAULT_ID_URL_TEMPLATES = {
    "CDS": "https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/protists/fasta/{id_lower}/cds/",
    "GFF": "https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/protists/gff3/{id_lower}/",
    "GENOME": "https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/protists/fasta/{id_lower}/dna/",
}
ENSEMBLGENOMES_DEFAULT_ID_URL_TEMPLATES = {
    "ensemblplants": ENSEMBLPLANTS_DEFAULT_ID_URL_TEMPLATES,
    "ensemblmetazoa": ENSEMBLMETAZOA_DEFAULT_ID_URL_TEMPLATES,
    "ensemblprotists": ENSEMBLPROTISTS_DEFAULT_ID_URL_TEMPLATES,
}
ENSEMBLGENOMES_FASTA_INDEX_URLS = {
    "ensemblplants": "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/",
    "ensemblmetazoa": "https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/metazoa/fasta/",
    "ensemblprotists": "https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/protists/fasta/",
}
PROVIDER_DEFAULT_ID_PAGE_URL_TEMPLATES = {
    "phycocosm": (
        "https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism={id}",
    ),
    "phytozome": (
        "https://phytozome-next.jgi.doe.gov/download/file/{id}",
        "https://phytozome-next.jgi.doe.gov/info/{id}",
        "https://data.jgi.doe.gov/refine-download/phytozome?genome_id={id}",
    ),
    "coge": (
        "https://genomevolution.org/coge/GenomeInfo.pl?gid={id}",
    ),
    "cngb": (
        "https://db.cngb.org/data_resources/project/{id}",
    ),
    "gwh": (
        "https://download.cncb.ac.cn/gwh/",
    ),
    "citrusgenomedb": (
        "https://www.citrusgenomedb.org/Analysis/{id}",
        "https://www.citrusgenomedb.org/organism/{id}",
    ),
    "figshare": (
        "https://figshare.com/articles/dataset/{id}",
    ),
    "plantgarden": (
        "https://plantgarden.jp/en/list/{id}/genome",
        "https://plantgarden.jp/en/list/{id}/genome/{id}",
    ),
    "plantaedb": (
        "https://plantaedb.com/",
    ),
    "fernbase": (
        "https://fernbase.org/ftp/{id}/",
    ),
    "veupathdb": (
        "https://veupathdb.org/",
    ),
    "dictybase": (
        "https://dictybase.org/Downloads/",
        "https://dictybase.org/download/",
    ),
    "insectbase": (
        "https://www.insect-genome.com/genome/{id}",
    ),
}


def strip_provider_prefix(source_id, provider):
    text = str(source_id or "").strip()
    prefix = provider.lower() + ":"
    if text.lower().startswith(prefix):
        return text[len(prefix) :]
    return text


def extract_ncbi_accession_from_source_id(source_id):
    text = str(source_id or "").strip()
    if NCBI_ASSEMBLY_ACCESSION_PATTERN.match(text):
        return text
    match = re.search(r"(GC[AF]_[0-9]+\.[0-9]+)", text, flags=re.IGNORECASE)
    if match is None:
        return ""
    return match.group(1)


def provider_env_prefix(provider):
    normalized = re.sub(r"[^A-Za-z0-9]+", "_", provider.upper())
    return "GG_{}".format(normalized)


def is_url_like(text):
    parsed = urlparse(str(text or ""))
    return parsed.scheme in ("http", "https", "ftp", "file")


def render_id_url_template(template, provider, source_id, species_key):
    source_clean = str(source_id or "").strip()
    id_value = strip_provider_prefix(source_clean, provider)
    if id_value == "":
        id_value = source_clean
    mapping = {
        "id": id_value,
        "id_raw": source_clean,
        "id_lower": id_value.lower(),
        "species_key": species_key,
        "provider": provider,
    }
    return template.format(**mapping)


def fetch_text_with_headers(url, timeout, headers):
    req_headers = dict(headers)
    if "User-Agent" not in req_headers:
        req_headers["User-Agent"] = "genegalleon-input-generation"
    request = Request(url, headers=req_headers)
    with urlopen(request, timeout=timeout) as response:
        payload = response.read()
    try:
        return payload.decode("utf-8")
    except UnicodeDecodeError:
        return payload.decode("utf-8", errors="replace")


def append_cookie_header(existing_cookie_header, added_cookie_header):
    existing = str(existing_cookie_header or "").strip()
    added = str(added_cookie_header or "").strip()
    if existing == "":
        return added
    if added == "":
        return existing
    return "{}; {}".format(existing, added)


def extract_hidden_input_value(html_text, input_name):
    pattern = (
        r"""<input[^>]*name=["']{name}["'][^>]*value=["']([^"']*)["'][^>]*>"""
    ).format(name=re.escape(str(input_name or "").strip()))
    match = re.search(pattern, html_text, flags=re.IGNORECASE)
    if match is None:
        return ""
    return str(match.group(1) or "").strip()


def fetch_jgi_session_cookie(login, password, timeout, headers):
    login_text = str(login or "").strip()
    password_text = str(password or "").strip()
    if login_text == "" or password_text == "":
        raise ValueError("JGI credentials are empty")
    signon_url = urljoin(resolve_jgi_signon_base_url().rstrip("/") + "/", "signon")
    cookie_jar = CookieJar()
    opener = build_opener(HTTPCookieProcessor(cookie_jar))
    req_headers = dict(headers)
    if "User-Agent" not in req_headers:
        req_headers["User-Agent"] = "genegalleon-input-generation"
    with opener.open(Request(signon_url, headers=req_headers), timeout=timeout) as response:
        login_html = response.read().decode("utf-8", errors="replace")
    form_match = re.search(
        r"""<form[^>]*action=["']([^"']+)["'][^>]*>(.*?)</form>""",
        login_html,
        flags=re.IGNORECASE | re.DOTALL,
    )
    if form_match is None:
        raise ValueError("JGI sign-on page did not expose a login form")
    action_url = urljoin(signon_url, str(form_match.group(1) or "").strip())
    form_html = str(form_match.group(2) or "")
    authenticity_token = extract_hidden_input_value(form_html, "authenticity_token")
    if authenticity_token == "":
        raise ValueError("JGI sign-on page did not expose authenticity_token")
    utf8_value = extract_hidden_input_value(form_html, "utf8") or "\u2713"
    commit_value = extract_hidden_input_value(form_html, "commit") or "Sign In"
    post_body = urlencode(
        {
            "utf8": utf8_value,
            "authenticity_token": authenticity_token,
            "login": login_text,
            "password": password_text,
            "commit": commit_value,
        }
    ).encode("utf-8")
    post_headers = dict(req_headers)
    post_headers["Content-Type"] = "application/x-www-form-urlencoded"
    with opener.open(Request(action_url, data=post_body, headers=post_headers), timeout=timeout) as response:
        final_url = response.geturl()
        final_html = response.read().decode("utf-8", errors="replace")
    cookie_header = "; ".join(
        "{}={}".format(cookie.name, cookie.value)
        for cookie in cookie_jar
        if str(cookie.name or "").strip() != ""
    )
    if cookie_header == "":
        raise ValueError("JGI sign-on did not yield session cookies")
    final_url_lower = str(final_url or "").lower()
    if "/signon" in final_url_lower and re.search(r"""name=["']password["']""", final_html, flags=re.IGNORECASE):
        raise ValueError("JGI sign-on did not complete; credentials may be invalid")
    return cookie_header


def parse_links_from_document(base_url, text):
    links = []
    seen = set()
    for href in re.findall(r"""href\s*=\s*["']([^"']+)["']""", text, flags=re.IGNORECASE):
        candidate = href.strip()
        if candidate == "" or candidate.startswith("javascript:") or candidate.startswith("#"):
            continue
        absolute = urljoin(base_url, candidate)
        if absolute in seen:
            continue
        seen.add(absolute)
        links.append(absolute)
    for raw in re.findall(r"""https?://[^\s"'<>]+""", text):
        candidate = raw.strip()
        if candidate in seen:
            continue
        seen.add(candidate)
        links.append(candidate)
    return links


def parse_dirname_from_index_link(url):
    path = urlparse(url).path
    if path.endswith("/"):
        return path.rstrip("/").split("/")[-1]
    if path.endswith("/index.html"):
        return path[: -len("/index.html")].rstrip("/").split("/")[-1]
    return ""


def ensembl_like_default_id_url_templates(provider):
    if provider == "ensembl":
        return ENSEMBL_DEFAULT_ID_URL_TEMPLATES
    return ENSEMBLGENOMES_DEFAULT_ID_URL_TEMPLATES.get(provider, {})


def fetch_ensemblgenomes_dir_ids(provider, timeout, headers):
    index_url = ENSEMBLGENOMES_FASTA_INDEX_URLS.get(provider, "")
    if index_url == "":
        return []
    html_text = fetch_text_with_headers(index_url, timeout, headers)
    out = []
    for link in parse_links_from_document(index_url, html_text):
        source_id = parse_dirname_from_index_link(link)
        if source_id in ("", "current", "fasta", "pub", "plants", "metazoa", "protists"):
            continue
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", source_id):
            continue
        out.append(source_id)
    return sorted(set(out), key=lambda value: (value.lower(), len(value)))


def expand_ensemblgenomes_id_candidates_from_available(id_candidates, available_ids):
    if not available_ids:
        return list(id_candidates)

    available_lookup = {value.lower(): value for value in available_ids}
    expanded = []
    seen = set()

    def add(value):
        text = str(value or "").strip()
        if text == "":
            return
        lowered = text.lower()
        if lowered in seen:
            return
        seen.add(lowered)
        expanded.append(text)

    for candidate in id_candidates:
        stripped = str(candidate or "").strip()
        if stripped == "":
            continue
        exact = available_lookup.get(stripped.lower())
        if exact is not None:
            add(exact)
        else:
            prefix_matches = [
                value
                for value in available_ids
                if value.lower().startswith(stripped.lower() + "_")
            ]
            if prefix_matches:
                add(sorted(prefix_matches, key=lambda value: (len(value), value.lower()))[0])
        add(stripped)
    return expanded


def expand_ensemblgenomes_id_candidates(provider, id_candidates, timeout, headers):
    available_ids = fetch_ensemblgenomes_dir_ids(provider, timeout, headers)
    return expand_ensemblgenomes_id_candidates_from_available(id_candidates, available_ids)


def resolve_oryza_minuta_gramene_base_url():
    return os.environ.get("GG_ORYZA_MINUTA_GRAMENE_BASE_URL", DEFAULT_ORYZA_MINUTA_GRAMENE_BASE_URL).rstrip("/")


def resolve_ddbj_search_api_base_url():
    return os.environ.get("GG_DDBJ_SEARCH_API_BASE_URL", DEFAULT_DDBJ_SEARCH_API_BASE_URL).rstrip("/")


def resolve_ddbj_public_wgs_base_url():
    return os.environ.get("GG_DDBJ_PUBLIC_WGS_BASE_URL", DEFAULT_DDBJ_PUBLIC_WGS_BASE_URL).rstrip("/")


def resolve_coge_api_base_url():
    return os.environ.get("GG_COGE_API_BASE_URL", DEFAULT_COGE_API_BASE_URL).rstrip("/")


def resolve_coge_web_base_url():
    return os.environ.get("GG_COGE_WEB_BASE_URL", DEFAULT_COGE_WEB_BASE_URL).rstrip("/")


def resolve_jgi_signon_base_url():
    return os.environ.get("GG_JGI_SIGNON_BASE_URL", DEFAULT_JGI_SIGNON_BASE_URL).rstrip("/")


def resolve_plantaedb_web_base_url():
    return os.environ.get("GG_PLANTAEDB_WEB_BASE_URL", DEFAULT_PLANTAEDB_WEB_BASE_URL).rstrip("/")


def resolve_plantgarden_web_base_url():
    return os.environ.get("GG_PLANTGARDEN_WEB_BASE_URL", DEFAULT_PLANTGARDEN_WEB_BASE_URL).rstrip("/")


def resolve_cngb_cnsa_base_url():
    return os.environ.get("GG_CNGB_CNSA_BASE_URL", DEFAULT_CNGB_CNSA_BASE_URL).rstrip("/")


def resolve_gwh_download_base_url():
    return os.environ.get("GG_GWH_DOWNLOAD_BASE_URL", DEFAULT_GWH_DOWNLOAD_BASE_URL).rstrip("/")


def resolve_gwh_web_base_url():
    return os.environ.get("GG_GWH_WEB_BASE_URL", DEFAULT_GWH_WEB_BASE_URL).rstrip("/")


def resolve_gwh_search_api_url():
    default = resolve_gwh_web_base_url().rstrip("/") + "/gwhSearch/api"
    return os.environ.get("GG_GWH_SEARCH_API_URL", default).rstrip("/")


def resolve_citrusgenomedb_web_base_url():
    return os.environ.get("GG_CITRUSGENOMEDB_WEB_BASE_URL", DEFAULT_CITRUSGENOMEDB_WEB_BASE_URL).rstrip("/")


def resolve_figshare_api_base_url():
    return os.environ.get("GG_FIGSHARE_API_BASE_URL", DEFAULT_FIGSHARE_API_BASE_URL).rstrip("/")


def resolve_veupathdb_service_base_url():
    return os.environ.get("GG_VEUPATHDB_SERVICE_BASE_URL", DEFAULT_VEUPATHDB_SERVICE_BASE_URL).rstrip("/")


def resolve_insectbase_api_base_url():
    return os.environ.get("GG_INSECTBASE_API_BASE_URL", DEFAULT_INSECTBASE_API_BASE_URL).rstrip("/")


def resolve_ncbi_eutils_base_url():
    return os.environ.get("GG_NCBI_EUTILS_BASE_URL", DEFAULT_NCBI_EUTILS_BASE_URL).rstrip("/")


def resolve_ncbi_ftp_base_url():
    return os.environ.get("GG_NCBI_FTP_BASE_URL", DEFAULT_NCBI_FTP_BASE_URL).rstrip("/")


def resolve_ncbi_datasets_base_url():
    return os.environ.get("GG_NCBI_DATASETS_BASE_URL", DEFAULT_NCBI_DATASETS_BASE_URL).rstrip("/")

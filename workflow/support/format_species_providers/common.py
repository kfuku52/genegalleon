"""Provider resolver implementation: common."""

import json
import re
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import parse_qs, unquote, urlparse
from urllib.request import Request

from format_species_annotations import (
    species_prefix_from_value,
)
from format_species_common import (
    is_fasta_filename,
    is_gff_filename,
    is_probable_cds_filename,
    is_probable_genome_filename,
    provider_candidate_sort_key,
)
from format_species_constants import (
    CITRUSGENOMEDB_ID_HINT_PATTERN,
    CNGB_ASSEMBLY_ACCESSION_PATTERN,
    CNGB_ID_HINT_PATTERN,
    COGE_ID_HINT_PATTERN,
    DDBJ_BIOPROJECT_PATTERN,
    DDBJ_ID_HINT_PATTERN,
    DDBJ_WGS_ACCESSION_PATTERN,
    DDBJ_WGS_MASTER_ACCESSION_PATTERN,
    DDBJ_WGS_PROJECT_PREFIX_PATTERN,
    DICTYBASE_ID_HINT_PATTERN,
    ENSEMBL_ID_HINT_PATTERN,
    ENSEMBLMETAZOA_ID_HINT_PATTERN,
    ENSEMBLPROTISTS_ID_HINT_PATTERN,
    FERNBASE_ID_HINT_PATTERN,
    FIGSHARE_ARTICLE_ID_PATTERN,
    FIGSHARE_ID_HINT_PATTERN,
    GWH_ASSEMBLY_ACCESSION_PATTERN,
    GWH_ASSEMBLY_ACCESSION_SEARCH_PATTERN,
    GWH_ID_HINT_PATTERN,
    INSECTBASE_IBG_ID_PATTERN,
    INSECTBASE_ID_HINT_PATTERN,
    ORYZA_MINUTA_ID_HINT_PATTERN,
    PLANTAEDB_ID_HINT_PATTERN,
    PLANTGARDEN_ID_HINT_PATTERN,
    TAXONOMIC_INFRASPECIFIC_RANK_ALIASES,
    VEUPATHDB_ID_HINT_PATTERN,
)
from format_species_network import guarded_urlopen as urlopen
from format_species_provider_config import (
    DEFAULT_INPUT_RELATIVE_DIRS,
    ENSEMBL_LIKE_PROVIDERS,
    ORYZA_MINUTA_PROVIDER,
)
from format_species_provider_urls import (
    NCBI_ASSEMBLY_ACCESSION_PATTERN,
    expand_ensemblgenomes_id_candidates_from_available,
    fetch_text_with_headers,
    is_url_like,
    parse_links_from_document,
    strip_provider_prefix,
)
from format_species_provider_urls import fetch_ensemblgenomes_dir_ids as _fetch_ensemblgenomes_dir_ids
from format_species_taxonomy import (
    build_species_key_from_tokens,
    tokenize_taxonomic_name,
)


def provider_raw_dir(provider, download_root, species_key):
    if provider in ENSEMBL_LIKE_PROVIDERS:
        return download_root / DEFAULT_INPUT_RELATIVE_DIRS[provider]
    if provider in (
        "phycocosm",
        "phytozome",
        "ncbi",
        "ddbj",
        "refseq",
        "genbank",
        "coge",
        "cngb",
        "gwh",
        "citrusgenomedb",
        "figshare",
        "plantgarden",
        "plantaedb",
        "flybase",
        "wormbase",
        "vectorbase",
        "fernbase",
        "veupathdb",
        "dictybase",
        "insectbase",
        ORYZA_MINUTA_PROVIDER,
        "direct",
        "local",
    ):
        return download_root / DEFAULT_INPUT_RELATIVE_DIRS[provider] / species_key
    raise ValueError("Unknown provider: {}".format(provider))


def infer_provider_from_id(source_id):
    if source_id == "":
        return ""
    lowered = source_id.lower()
    if NCBI_ASSEMBLY_ACCESSION_PATTERN.match(source_id):
        return "ncbi"
    if "ncbi.nlm.nih.gov/datasets/genome/" in lowered:
        return "ncbi"
    if DDBJ_ID_HINT_PATTERN.match(source_id):
        return "ddbj"
    if ENSEMBL_ID_HINT_PATTERN.match(source_id):
        return "ensembl"
    if ENSEMBLMETAZOA_ID_HINT_PATTERN.match(source_id):
        return "ensemblmetazoa"
    if ENSEMBLPROTISTS_ID_HINT_PATTERN.match(source_id):
        return "ensemblprotists"
    if FERNBASE_ID_HINT_PATTERN.match(source_id):
        return "fernbase"
    if VEUPATHDB_ID_HINT_PATTERN.match(source_id):
        return "veupathdb"
    if DICTYBASE_ID_HINT_PATTERN.match(source_id):
        return "dictybase"
    if INSECTBASE_ID_HINT_PATTERN.match(source_id):
        return "insectbase"
    if ORYZA_MINUTA_ID_HINT_PATTERN.match(source_id):
        return ORYZA_MINUTA_PROVIDER
    if INSECTBASE_IBG_ID_PATTERN.match(source_id):
        return "insectbase"
    if "ftp.ensembl.org" in lowered or "ensembl.org/pub/current_" in lowered:
        return "ensembl"
    if (
        "ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/" in lowered
        or "ftp.ensemblgenomes.ebi.ac.uk/pub/current/plants/" in lowered
    ):
        return "ensemblplants"
    if "ftp.ensemblgenomes.ebi.ac.uk/pub/current/metazoa/" in lowered:
        return "ensemblmetazoa"
    if "ftp.ensemblgenomes.ebi.ac.uk/pub/current/protists/" in lowered:
        return "ensemblprotists"
    if COGE_ID_HINT_PATTERN.match(source_id):
        return "coge"
    if GWH_ID_HINT_PATTERN.match(source_id):
        return "gwh"
    if CITRUSGENOMEDB_ID_HINT_PATTERN.match(source_id):
        return "citrusgenomedb"
    if FIGSHARE_ID_HINT_PATTERN.match(source_id):
        return "figshare"
    if PLANTGARDEN_ID_HINT_PATTERN.match(source_id):
        return "plantgarden"
    if PLANTAEDB_ID_HINT_PATTERN.match(source_id):
        return "plantaedb"
    if CNGB_ID_HINT_PATTERN.match(source_id):
        return "cngb"
    if extract_ddbj_bioproject_id(source_id) != "":
        return "ddbj"
    if extract_ddbj_wgs_prefix(source_id) != "":
        return "ddbj"
    if GWH_ASSEMBLY_ACCESSION_PATTERN.match(source_id):
        return "gwh"
    if CNGB_ASSEMBLY_ACCESSION_PATTERN.match(source_id):
        return "cngb"
    if "genomevolution.org" in lowered:
        return "coge"
    if "/gwh/" in lowered:
        return "gwh"
    if "citrusgenomedb.org/" in lowered:
        return "citrusgenomedb"
    if "figshare.com/articles/" in lowered or "api.figshare.com/" in lowered:
        return "figshare"
    if "plantgarden.jp/" in lowered:
        return "plantgarden"
    if "plantaedb.com/" in lowered:
        return "plantaedb"
    if "cngb.org" in lowered or "cncb.ac.cn" in lowered:
        return "cngb"
    if "ddbj.nig.ac.jp/" in lowered or "getentry.ddbj.nig.ac.jp/" in lowered:
        return "ddbj"
    if "flybase.org" in lowered:
        return "flybase"
    if "wormbase.org" in lowered:
        return "wormbase"
    if "vectorbase.org" in lowered:
        return "vectorbase"
    if "fernbase.org" in lowered:
        return "fernbase"
    if "veupathdb.org" in lowered:
        return "veupathdb"
    if "dictybase.org" in lowered:
        return "dictybase"
    if "insect-genome.com" in lowered:
        return "insectbase"
    return ""


def fetch_ensemblgenomes_dir_ids(provider, timeout, headers):
    return _fetch_ensemblgenomes_dir_ids(provider, timeout, headers)


def expand_ensemblgenomes_id_candidates(provider, id_candidates, timeout, headers):
    available_ids = fetch_ensemblgenomes_dir_ids(provider, timeout, headers)
    return expand_ensemblgenomes_id_candidates_from_available(id_candidates, available_ids)


def is_probable_cds_url(provider, url):
    lower = url.lower()
    return is_probable_cds_filename(provider, lower)


def is_probable_gff_url(url):
    return is_gff_filename(url.lower())


def is_probable_genome_url(provider, url):
    path_lower = urlparse(url).path.lower()
    if provider == "citrusgenomedb":
        filename_lower = Path(path_lower).name
        if "/assembly/" in path_lower and is_fasta_filename(path_lower):
            return True
        if (
            any(marker in path_lower for marker in ("/genes/", "/annotation/"))
            and is_fasta_filename(path_lower)
            and any(
                marker in filename_lower
                for marker in ("gene", "genes", "mrna", "rna", "cdna", "cds", "protein", "pep", "transcript")
            )
        ):
            return False
    return is_probable_genome_filename(provider, path_lower)


def select_best_url_for_label(provider, label, candidates):
    if len(candidates) == 0:
        return ""
    filtered = []
    for candidate in candidates:
        path_lower = urlparse(candidate).path.lower()
        if label == "CDS":
            if is_probable_cds_url(provider, path_lower):
                filtered.append(candidate)
            continue
        if label == "GFF":
            if is_probable_gff_url(path_lower):
                filtered.append(candidate)
            continue
        if label == "GENOME":
            if is_probable_genome_url(provider, candidate):
                filtered.append(candidate)
            continue
    if len(filtered) == 0:
        return ""
    preferred = sorted(
        filtered,
        key=lambda x: provider_candidate_sort_key(provider, label, urlparse(x).path.split("/")[-1]),
    )
    return preferred[0]


def ensembl_like_gff_url_is_disfavored(url):
    name = urlparse(str(url or "")).path.split("/")[-1].lower()
    if name == "":
        return False
    return (
        "abinitio" in name
        or re.search(r"(?:^|[._-])(?:chr|chromosome)(?:[._-])", name) is not None
        or ".primary_assembly." in name
    )


def resolve_preferred_ensembl_like_gff_url(provider, gff_url, timeout, headers):
    if provider not in ENSEMBL_LIKE_PROVIDERS or gff_url == "":
        return gff_url
    if not ensembl_like_gff_url_is_disfavored(gff_url):
        return gff_url
    parsed = urlparse(gff_url)
    if parsed.scheme not in ("http", "https", "ftp"):
        return gff_url
    parent_url = gff_url.rsplit("/", 1)[0] + "/"
    discovered = resolve_urls_from_index_url(provider, parent_url, timeout, headers)
    preferred = str(discovered.get("gff_url") or "").strip()
    if preferred != "" and preferred != gff_url:
        return preferred
    return gff_url


def remove_stale_ensembl_like_partial_gff_outputs(provider, species_prefix, output_dir, keep_path):
    if provider not in ENSEMBL_LIKE_PROVIDERS or keep_path is None:
        return []
    keep_path = Path(keep_path)
    if ensembl_like_gff_url_is_disfavored(keep_path.name):
        return []
    removed = []
    for path in sorted(Path(output_dir).iterdir()):
        if path == keep_path or not path.is_file():
            continue
        if species_prefix_from_value(path.name) != species_prefix:
            continue
        if not is_gff_filename(path.name):
            continue
        if not ensembl_like_gff_url_is_disfavored(path.name):
            continue
        path.unlink()
        audit_path = Path(str(path) + ".repair.json")
        if audit_path.is_file():
            audit_path.unlink()
        removed.append(path.name)
    return removed


def resolve_urls_from_index_url(provider, index_url, timeout, headers):
    text = fetch_text_with_headers(index_url, timeout, headers)
    links = parse_links_from_document(index_url, text)
    return {
        "cds_url": select_best_url_for_label(provider, "CDS", links),
        "gff_url": select_best_url_for_label(provider, "GFF", links),
        "genome_url": select_best_url_for_label(provider, "GENOME", links),
    }


def fetch_json_with_headers(url, timeout, headers):
    text = fetch_text_with_headers(url, timeout, headers)
    return json.loads(text)


def extract_ddbj_bioproject_id(source_id):
    text = str(source_id or "").strip()
    if text == "":
        return ""
    match = DDBJ_BIOPROJECT_PATTERN.search(text)
    if match is None:
        return ""
    return str(match.group(1) or "").upper()


def extract_ddbj_wgs_prefix(source_id):
    text = str(source_id or "").strip()
    if text == "":
        return ""
    stripped = strip_provider_prefix(text, "ddbj").strip()
    if stripped == "":
        stripped = text
    if is_url_like(stripped):
        parsed = urlparse(stripped)
        for key in ("accession_number", "id"):
            for value in parse_qs(parsed.query).get(key, []):
                prefix = extract_ddbj_wgs_prefix(value)
                if prefix != "":
                    return prefix
        for part in parsed.path.split("/"):
            prefix = extract_ddbj_wgs_prefix(unquote(part))
            if prefix != "":
                return prefix
        return ""

    candidate = re.sub(r"(?:[.]gz|[.]gbff|[.]gbk|[.]gb|[.]genbank)+$", "", stripped, flags=re.IGNORECASE)
    candidate = candidate.split()[0].strip().upper()
    if candidate == "":
        return ""
    if DDBJ_WGS_PROJECT_PREFIX_PATTERN.fullmatch(candidate):
        return candidate
    for pattern in (DDBJ_WGS_MASTER_ACCESSION_PATTERN, DDBJ_WGS_ACCESSION_PATTERN):
        match = pattern.match(candidate)
        if match is not None:
            return str(match.group(1) or "").upper()
    search = re.search(r"\b([A-Z]{4,6})(?:0{8,9}|[0-9]{8,9})(?:\.[0-9]+)?\b", candidate, flags=re.IGNORECASE)
    if search is None:
        return ""
    return str(search.group(1) or "").upper()


def remote_resource_exists(url, timeout, headers):
    req_headers = dict(headers)
    if "User-Agent" not in req_headers:
        req_headers["User-Agent"] = "genegalleon-input-generation"
    try:
        request = Request(url, headers=req_headers, method="HEAD")
        with urlopen(request, timeout=timeout):
            return True
    except HTTPError as exc:
        if exc.code == 404:
            return False
        if exc.code not in (405, 501):
            raise
    fallback_headers = dict(req_headers)
    if "Range" not in fallback_headers:
        fallback_headers["Range"] = "bytes=0-0"
    try:
        request = Request(url, headers=fallback_headers)
        with urlopen(request, timeout=timeout):
            return True
    except HTTPError as exc:
        if exc.code == 404:
            return False
        raise


def parse_species_key_candidate(text):
    return build_species_key_from_tokens(tokenize_taxonomic_name(text))


def source_id_candidates(provider, source_id, species_key):
    candidates = []

    def add(value):
        text = str(value or "").strip()
        if text == "":
            return
        if text in candidates:
            return
        candidates.append(text)

    source_clean = str(source_id or "").strip()
    stripped = strip_provider_prefix(source_clean, provider)
    add(source_clean)
    add(stripped)

    if "_" in stripped:
        add(stripped.rsplit("_", 1)[-1])
    if species_key != "":
        add(species_key)
        if "_" in species_key:
            add(species_key.rsplit("_", 1)[-1])

    if is_url_like(stripped):
        parsed = urlparse(stripped)
        query = parse_qs(parsed.query)
        if provider == "coge":
            for gid in query.get("gid", []):
                add(gid)
        if provider == "cngb":
            for qval in query.get("q", []):
                add(qval)
            parts = [token for token in parsed.path.split("/") if token != ""]
            for idx, token in enumerate(parts[:-1]):
                if token.lower() == "assembly":
                    add(parts[idx + 1])
        if provider == "gwh":
            for qvals in query.values():
                for qval in qvals:
                    for match in GWH_ASSEMBLY_ACCESSION_SEARCH_PATTERN.findall(qval):
                        add(match)
            for part in parsed.path.split("/"):
                for match in GWH_ASSEMBLY_ACCESSION_SEARCH_PATTERN.findall(unquote(part)):
                    add(match)

    if provider == "phytozome":
        match = re.search(r"_([0-9]+)_", stripped)
        if match is not None:
            add(match.group(1))
    if provider in ENSEMBL_LIKE_PROVIDERS:
        species_tokens = [token for token in str(species_key or "").split("_") if token != ""]
        if len(species_tokens) >= 2:
            add("{}_{}".format(species_tokens[0], species_tokens[1]))
        if len(species_tokens) >= 3 and species_tokens[1].lower() == "x":
            add("{}_{}_{}".format(species_tokens[0], species_tokens[1], species_tokens[2]))
        if len(species_tokens) >= 4 and species_tokens[2].lower() in TAXONOMIC_INFRASPECIFIC_RANK_ALIASES:
            add("{}_{}".format(species_tokens[0], species_tokens[1]))
        if len(species_tokens) >= 4 and species_tokens[3].lower() in ("group", "subgroup"):
            add("{}_{}_{}".format(species_tokens[0], species_tokens[1], species_tokens[2]))
        if NCBI_ASSEMBLY_ACCESSION_PATTERN.match(source_clean):
            accession_candidate = source_clean.lower().replace("_", "").replace(".", "v")
            for base_candidate in list(candidates):
                if re.search(r"[A-Za-z]", base_candidate) is None:
                    continue
                add("{}_{}".format(base_candidate, accession_candidate))
                add("{}_{}cm".format(base_candidate, accession_candidate))
    if provider == "gwh":
        for match in GWH_ASSEMBLY_ACCESSION_SEARCH_PATTERN.findall(stripped):
            add(match)
    if provider == "figshare":
        article_id = extract_figshare_article_id_candidate(stripped)
        if article_id != "":
            add(article_id)

    return candidates


def normalize_manifest_source_id(provider, source_id):
    text = str(source_id or "").strip()
    if text == "":
        return ""
    if provider in ("local", "direct"):
        return text
    if is_url_like(text):
        return text
    match = re.match(r"^(\S+)\s+\(.*\)$", text)
    if match is not None:
        return match.group(1)
    return text


def normalize_lookup_text(text):
    return re.sub(r"[^a-z0-9]+", "", str(text or "").lower())


def extract_figshare_article_id_candidate(source_id):
    text = str(source_id or "").strip()
    if text == "":
        return ""
    stripped = strip_provider_prefix(text, "figshare").strip()
    for candidate in (text, stripped):
        if str(candidate).isdigit():
            return str(candidate)
        match = FIGSHARE_ARTICLE_ID_PATTERN.search(str(candidate))
        if match is not None:
            return str(match.group(1) or "").strip()
        if is_url_like(candidate):
            parts = [unquote(part) for part in urlparse(candidate).path.split("/") if part]
            lowered_parts = [part.lower() for part in parts]
            if "articles" in lowered_parts:
                for part in reversed(parts):
                    if str(part).isdigit():
                        return str(part)
    return ""

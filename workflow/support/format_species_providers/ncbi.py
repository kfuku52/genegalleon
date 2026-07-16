"""Provider resolver implementation: ncbi."""

import json
import os
import time
from pathlib import Path
from urllib.parse import quote, urlparse
from urllib.request import Request

from format_species_common import (
    is_transient_network_error,
    parse_positive_int,
)
from format_species_constants import (
    DEFAULT_NCBI_EUTILS_RPS_NO_API_KEY,
    DEFAULT_NCBI_EUTILS_RPS_WITH_API_KEY,
    _ncbi_eutils_rate_lock,
    _ncbi_eutils_request_times,
)
from format_species_network import guarded_urlopen as urlopen
from format_species_provider_urls import (
    extract_ncbi_accession_from_source_id,
    resolve_ncbi_eutils_base_url,
    resolve_ncbi_ftp_base_url,
)

from .common import (
    parse_species_key_candidate,
)


def resolve_ncbi_api_key():
    return os.environ.get("GG_NCBI_API_KEY", "").strip()


def resolve_ncbi_eutils_max_rps():
    fallback = DEFAULT_NCBI_EUTILS_RPS_NO_API_KEY
    if resolve_ncbi_api_key() != "":
        fallback = DEFAULT_NCBI_EUTILS_RPS_WITH_API_KEY
    env_value = os.environ.get("GG_NCBI_EUTILS_MAX_RPS", "").strip()
    return parse_positive_int(env_value, fallback)


def throttle_ncbi_eutils_request():
    max_rps = resolve_ncbi_eutils_max_rps()
    if max_rps < 1:
        return
    window_seconds = 1.0
    while True:
        with _ncbi_eutils_rate_lock:
            now = time.time()
            while len(_ncbi_eutils_request_times) > 0 and (now - _ncbi_eutils_request_times[0]) >= window_seconds:
                _ncbi_eutils_request_times.popleft()
            if len(_ncbi_eutils_request_times) < max_rps:
                _ncbi_eutils_request_times.append(now)
                return
            oldest = _ncbi_eutils_request_times[0]
            sleep_seconds = window_seconds - (now - oldest)
        time.sleep(max(0.001, sleep_seconds))


def fetch_json(url, timeout):
    request = Request(url, headers={"User-Agent": "genegalleon-input-generation"})
    last_error = None
    for attempt in range(1, 4):
        try:
            with urlopen(request, timeout=timeout) as response:
                payload = response.read()
            return json.loads(payload.decode("utf-8"))
        except Exception as exc:
            last_error = exc
            if attempt >= 3 or not is_transient_network_error(exc):
                raise
            time.sleep(float(attempt))
    raise last_error


def normalize_ncbi_ftp_path(ftppath):
    if ftppath == "":
        return ""
    ftp_base = resolve_ncbi_ftp_base_url()
    if ftppath.startswith("ftp://ftp.ncbi.nlm.nih.gov"):
        return ftp_base + ftppath[len("ftp://ftp.ncbi.nlm.nih.gov") :]
    if ftppath.startswith("https://ftp.ncbi.nlm.nih.gov"):
        return ftp_base + ftppath[len("https://ftp.ncbi.nlm.nih.gov") :]
    if ftppath.startswith("http://ftp.ncbi.nlm.nih.gov"):
        return ftp_base + ftppath[len("http://ftp.ncbi.nlm.nih.gov") :]
    if ftppath.startswith("ftp://"):
        return "https://" + ftppath[len("ftp://") :]
    return ftppath


def infer_ncbi_species_key_from_doc(doc, fallback):
    candidates = []
    speciesname = str(doc.get("speciesname", "") or "").strip()
    organism = str(doc.get("organism", "") or "").strip()
    if speciesname != "":
        candidates.append(speciesname)
    if organism != "":
        candidates.append(organism)
    for raw in candidates:
        inferred = parse_species_key_candidate(raw)
        if inferred != "":
            return inferred
    return fallback


def resolve_ncbi_download_urls_from_id(source_id, timeout, ncbi_source="auto"):
    accession = extract_ncbi_accession_from_source_id(source_id)
    if accession == "":
        raise ValueError(
            "ncbi id must include an assembly accession like GCF_000001405.40 or GCA_000001405.29: {}".format(source_id)
        )

    eutils_base = resolve_ncbi_eutils_base_url()
    term = quote("{}[Assembly Accession]".format(accession), safe="")
    api_key = resolve_ncbi_api_key()
    api_key_query = ""
    if api_key != "":
        api_key_query = "&api_key={}".format(quote(api_key, safe=""))
    esearch_url = "{}/esearch.fcgi?db=assembly&term={}&retmode=json{}".format(eutils_base, term, api_key_query)
    throttle_ncbi_eutils_request()
    esearch_data = fetch_json(esearch_url, timeout=timeout)
    uid_list = esearch_data.get("esearchresult", {}).get("idlist", [])
    if len(uid_list) == 0:
        raise ValueError("NCBI assembly was not found for id: {}".format(source_id))
    uid = uid_list[0]

    esummary_url = "{}/esummary.fcgi?db=assembly&id={}&retmode=json{}".format(
        eutils_base,
        quote(uid, safe=""),
        api_key_query,
    )
    throttle_ncbi_eutils_request()
    esummary_data = fetch_json(esummary_url, timeout=timeout)
    doc = esummary_data.get("result", {}).get(uid, {})

    ftppath_refseq = str(doc.get("ftppath_refseq", "") or "").strip()
    ftppath_genbank = str(doc.get("ftppath_genbank", "") or "").strip()
    selected_source = "auto"
    if ncbi_source == "refseq":
        selected_source = "refseq"
        ftp_dir = ftppath_refseq if ftppath_refseq != "" else ftppath_genbank
        if ftppath_refseq == "" and ftppath_genbank != "":
            selected_source = "genbank"
    elif ncbi_source == "genbank":
        selected_source = "genbank"
        ftp_dir = ftppath_genbank if ftppath_genbank != "" else ftppath_refseq
        if ftppath_genbank == "" and ftppath_refseq != "":
            selected_source = "refseq"
    else:
        ftp_dir = ftppath_refseq if ftppath_refseq != "" else ftppath_genbank
        if ftppath_refseq != "":
            selected_source = "refseq"
        elif ftppath_genbank != "":
            selected_source = "genbank"
    if ftp_dir == "":
        raise ValueError("NCBI FTP path was not found in assembly summary for id: {}".format(source_id))

    normalized_ftp_dir = normalize_ncbi_ftp_path(ftp_dir).rstrip("/")
    assembly_dir_name = Path(urlparse(ftp_dir).path.rstrip("/")).name
    if assembly_dir_name == "":
        raise ValueError("Could not parse assembly directory name from: {}".format(ftp_dir))

    cds_filename = "{}_cds_from_genomic.fna.gz".format(assembly_dir_name)
    gff_filename = "{}_genomic.gff.gz".format(assembly_dir_name)
    gbff_filename = "{}_genomic.gbff.gz".format(assembly_dir_name)
    genome_filename = "{}_genomic.fna.gz".format(assembly_dir_name)
    species_key = infer_ncbi_species_key_from_doc(doc, accession)
    return {
        "species_key": species_key,
        "cds_url": "{}/{}".format(normalized_ftp_dir, cds_filename),
        "gff_url": "{}/{}".format(normalized_ftp_dir, gff_filename),
        "gbff_url": "{}/{}".format(normalized_ftp_dir, gbff_filename),
        "genome_url": "{}/{}".format(normalized_ftp_dir, genome_filename),
        "cds_filename": cds_filename,
        "gff_filename": gff_filename,
        "gbff_filename": gbff_filename,
        "genome_filename": genome_filename,
        "ncbi_source_db": selected_source,
    }

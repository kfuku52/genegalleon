"""Resolve published 1KP coding nucleotides by the original four-letter sample code."""
import re
from functools import lru_cache
from urllib.parse import urljoin, urlparse

from .common import fetch_json_with_headers

SAMPLE_INDEX = "https://db.cngb.org/onekp/api/sample/get_samples/?per_page=2000&page=1"


@lru_cache(maxsize=8)
def sample_index(timeout, headers):
    data = fetch_json_with_headers(SAMPLE_INDEX, timeout, dict(headers))["data"]
    samples = data["samples"]
    if len(samples) != data["total"]:
        raise ValueError("1KP sample index is incomplete; refusing a partial code lookup")
    index = {}
    for sample in samples:
        code = sample["code"]
        if not re.fullmatch(r"[A-Z]{4}", code) or code in index:
            raise ValueError("1KP sample index contains an invalid or duplicate code")
        index[code] = sample
    return index


def resolve_onekp_download_urls_from_id(source_id, species_key, timeout, headers):
    code = str(source_id).strip().upper()
    if not re.fullmatch(r"[A-Z]{4}", code):
        raise ValueError("provider=onekp requires the four-letter sample code, not a species name")
    sample = sample_index(float(timeout), tuple(sorted(headers.items()))).get(code)
    if sample is None:
        raise ValueError("Unknown 1KP sample code: " + code)
    if sample.get("is_mislabeled") is not False:
        raise ValueError("1KP sample is mislabeled or its identity is unverified: " + code)
    base = sample["assembly_url"]
    parsed = urlparse(base)
    if parsed.scheme != "https" or parsed.netloc != "ftp.cngb.org" or not parsed.path.startswith("/pub/SciRAID/onekp/assemblies/" + code + "-") or not base.endswith("/"):
        raise ValueError("Unexpected 1KP sample source location")
    # These are the published coding nucleotide sequences. SOAP assemblies and
    # translated proteins are different biological roles and must not enter CDS.
    return {"cds_url": urljoin(base, code + "-translated-nucleotides.fa.gz"),
            "cds_filename": code + ".cds.fa.gz",
            "species_key": species_key or sample["species"].replace(" ", "_"),
            "onekp_reported_species": sample["species"], "onekp_sample_code": code}

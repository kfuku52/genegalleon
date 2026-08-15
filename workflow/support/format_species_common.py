"""Shared filename classification, candidate ordering, and network predicates."""

import bz2
import gzip
import re
import socket
from http.client import IncompleteRead, RemoteDisconnected
from urllib.error import HTTPError, URLError

from format_species_constants import (
    FASTA_EXTENSIONS,
    FERNBASE_COMBINED_FILENAME_MARKER,
    GENBANK_EXTENSIONS,
    GFF_EXTENSIONS,
    TRANSIENT_HTTP_STATUS_CODES,
)
from format_species_provider_config import ENSEMBL_LIKE_PROVIDERS, ORYZA_MINUTA_PROVIDER

FERNBASE_GFF_EXCLUDE_PATTERN = re.compile(
    r"(?:^|[._-])(te|teanno|repeat|repeats|transpos(?:on)?)(?:[._-]|$)",
    re.IGNORECASE,
)
GFF_EXCLUDE_PATTERN = re.compile(
    r"(?:^|[._-])(?:edta|red|repeat|repeats|te|teanno|transpos(?:on)?)(?:[._-]|$)",
    re.IGNORECASE,
)
GFF_REDUCED_MODEL_PATTERN = re.compile(
    r"(?:only[._-]?long(?:est)?[._-]?transcripts?|primary[._-]?transcripts?|longest)",
    re.IGNORECASE,
)


def parse_positive_int(value, fallback):
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return int(fallback)
    if parsed < 1:
        return int(fallback)
    return parsed


def is_probable_cds_filename(provider, name):
    lower = name.lower()
    if provider == "fernbase" and (lower.endswith(".cds") or lower.endswith(".cds.gz")):
        return True
    if not is_fasta_filename(name):
        return False
    return any(marker in lower for marker in ("cds", "transcript", "mrna", "cdna", "genecatalog_cds"))


def provider_candidate_sort_key(provider, label, name):
    lower = str(name or "").lower()
    label_upper = str(label or "").upper()
    if provider == ORYZA_MINUTA_PROVIDER:
        if label_upper == "CDS":
            return (
                0 if ".merged." in lower else 1,
                0 if "cds" in lower else 1,
                lower,
            )
        if label_upper == "GFF":
            return (
                0 if ".merged." in lower else 1,
                lower,
            )
        if label_upper == "GENOME":
            return (
                0 if ".merged." in lower else 1,
                0 if any(marker in lower for marker in ("dna", "genome", "toplevel")) else 1,
                lower,
            )
    if provider == "citrusgenomedb":
        if label_upper == "CDS":
            return (
                0 if "cds" in lower else 1,
                1 if any(marker in lower for marker in ("cdna", "mrna", "transcript")) else 0,
                1 if any(marker in lower for marker in ("protein", "pep")) else 0,
                lower,
            )
        if label_upper == "GFF":
            return (
                1 if any(marker in lower for marker in ("repeat", "interpro", "go", "ipr")) else 0,
                lower,
            )
        if label_upper == "GENOME":
            return (
                0 if any(marker in lower for marker in ("genome", "assembly", "hap", "chr", "chrom")) else 1,
                lower,
            )
    if provider == "figshare":
        if label_upper == "CDS":
            return (
                0 if "cds" in lower else 1,
                1 if any(marker in lower for marker in ("cdna", "mrna", "transcript")) else 0,
                1 if any(marker in lower for marker in ("protein", "pep")) else 0,
                lower,
            )
        if label_upper == "GFF":
            return (
                1 if GFF_EXCLUDE_PATTERN.search(lower) else 0,
                1 if GFF_REDUCED_MODEL_PATTERN.search(lower) else 0,
                lower,
            )
        if label_upper == "GENOME":
            return (
                0
                if any(
                    marker in lower for marker in ("genome", "assembly", "genomic", "chromosome", "_chr", ".chr", "hap")
                )
                else 1,
                lower,
            )
    if provider != "fernbase":
        if provider in ENSEMBL_LIKE_PROVIDERS and label_upper == "GFF":
            return (
                1 if "abinitio" in lower else 0,
                1 if re.search(r"(?:^|[._-])(?:chr|chromosome)(?:[._-])", lower) else 0,
                1 if ".primary_assembly." in lower else 0,
                lower,
            )
        if label_upper == "GENOME":
            return (".chromosome." in lower, lower)
        if label_upper == "CDS":
            return (
                0 if re.search(r"(?:^|[._-])cds(?:[._-]|$)", lower) else 1,
                1 if any(marker in lower for marker in ("cdna", "mrna", "transcript")) else 0,
                1 if any(marker in lower for marker in ("protein", "pep")) else 0,
                lower,
            )
        if label_upper == "GFF":
            return (
                1 if GFF_EXCLUDE_PATTERN.search(lower) else 0,
                1 if GFF_REDUCED_MODEL_PATTERN.search(lower) else 0,
                lower,
            )
        return (lower,)
    if label_upper == "CDS":
        return (
            0 if FERNBASE_COMBINED_FILENAME_MARKER in lower else 1,
            0 if "highconfidence" in lower else 1,
            1 if "lowconfidence" in lower else 0,
            0 if "cds" in lower else 1,
            1 if any(marker in lower for marker in ("transcript", "mrna", "cdna")) else 0,
            lower,
        )
    if label_upper == "GFF":
        return (
            1 if FERNBASE_GFF_EXCLUDE_PATTERN.search(lower) else 0,
            0 if FERNBASE_COMBINED_FILENAME_MARKER in lower else 1,
            0 if "highconfidence" in lower else 1,
            1 if "lowconfidence" in lower else 0,
            lower,
        )
    if label_upper == "GENOME":
        return (
            1 if "chloroplast" in lower else 0,
            0 if any(marker in lower for marker in ("genome", "assembly", "chr", "chromosome")) else 1,
            lower,
        )
    return (lower,)


def gff_contains_cds_feature(path):
    """Return whether a local GFF candidate contains at least one CDS row."""
    name = str(path.name or "").lower()
    opener = gzip.open if name.endswith(".gz") else bz2.open if name.endswith(".bz2") else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            for raw_line in handle:
                if raw_line == "" or raw_line.startswith("#"):
                    continue
                parts = raw_line.rstrip("\n\r").split("\t", 3)
                if len(parts) >= 3 and parts[2].strip().lower() == "cds":
                    return True
    except (OSError, UnicodeError):
        return False
    return False


def gff_candidate_sort_key(provider, path):
    """Prefer full gene-model GFFs over repeat and reduced annotations."""
    name = str(path.name or "")
    excluded_by_name = GFF_EXCLUDE_PATTERN.search(name) is not None
    reduced_by_name = GFF_REDUCED_MODEL_PATTERN.search(name) is not None
    contains_cds = False if excluded_by_name else gff_contains_cds_feature(path)
    return (
        1 if excluded_by_name else 0,
        1 if reduced_by_name else 0,
        0 if contains_cds else 1,
        provider_candidate_sort_key(provider, "GFF", name),
    )


def figshare_gff_candidate_sort_key(path):
    """Backward-compatible wrapper for the shared GFF candidate policy."""
    return gff_candidate_sort_key("figshare", path)


def is_transient_network_error(exc):
    pending = [exc]
    seen = set()
    while len(pending) > 0:
        current = pending.pop()
        if current is None:
            continue
        current_id = id(current)
        if current_id in seen:
            continue
        seen.add(current_id)

        if isinstance(current, HTTPError):
            if int(getattr(current, "code", 0)) in TRANSIENT_HTTP_STATUS_CODES:
                return True
            continue

        if isinstance(
            current,
            (
                RemoteDisconnected,
                IncompleteRead,
                TimeoutError,
                socket.timeout,
                ConnectionResetError,
                ConnectionAbortedError,
                BrokenPipeError,
                EOFError,
            ),
        ):
            return True

        if isinstance(current, OSError):
            if getattr(current, "errno", None) in (104, 110, 111, 113):
                return True
            lowered_message = str(current).lower()
            if (
                "unexpected_eof_while_reading" in lowered_message
                or "eof occurred in violation of protocol" in lowered_message
            ):
                return True

        if isinstance(current, URLError):
            reason = getattr(current, "reason", None)
            if isinstance(reason, str):
                lowered = reason.lower()
                if (
                    "timed out" in lowered
                    or "temporarily unavailable" in lowered
                    or "connection reset" in lowered
                    or "remote end closed connection" in lowered
                    or "unexpected_eof_while_reading" in lowered
                    or "eof occurred in violation of protocol" in lowered
                ):
                    return True
            elif isinstance(reason, BaseException):
                pending.append(reason)

        cause = getattr(current, "__cause__", None)
        context = getattr(current, "__context__", None)
        if isinstance(cause, BaseException):
            pending.append(cause)
        if isinstance(context, BaseException):
            pending.append(context)
    return False


def is_fasta_filename(name):
    lower = name.lower()
    return any(lower.endswith(ext) for ext in FASTA_EXTENSIONS)


def is_gff_filename(name):
    lower = name.lower()
    return any(lower.endswith(ext) for ext in GFF_EXTENSIONS)


def is_gbff_filename(name):
    lower = name.lower()
    return any(lower.endswith(ext) for ext in GENBANK_EXTENSIONS)


def has_haplotype_assembly_marker(name):
    lower = str(name or "").lower()
    return re.search(r"(?:^|[._-])(?:hap[0-9][0-9a-z]*|haplotype[0-9a-z]*)(?:[._-]|$)", lower) is not None


def is_probable_genome_filename(provider, name):
    if not is_fasta_filename(name):
        return False
    lower = name.lower()
    non_genome_markers = (
        "cds",
        "mrna",
        "trna",
        "rrna",
        "ncrna",
        "rna_from_genomic",
        "transcript",
        "cdna",
        "pep",
        "protein",
    )
    if any(marker in lower for marker in non_genome_markers):
        return False
    if provider == "fernbase":
        # FernBase often exposes the assembly as a plain ".fa"/".fasta" filename.
        return True
    if provider == "plantgarden":
        return any(
            marker in lower
            for marker in ("genomic", "genome", "assembly", "pseudomolecule", "pmol")
        )
    if provider in ("ncbi", "refseq", "genbank", "plantaedb"):
        return "genomic" in lower
    if provider == "figshare":
        return any(
            marker in lower for marker in ("genome", "assembly", "genomic", "dna", "chromosome", "_chr", ".chr", "hap")
        )
    if provider in ENSEMBL_LIKE_PROVIDERS:
        return any(marker in lower for marker in ("dna", "genome", "toplevel", "primary_assembly", "chromosome"))
    return any(
        marker in lower
        for marker in (
            "genome",
            "assembly",
            "genomic",
            "dna",
            "chromosome",
            "_chr",
            ".chr",
            "chr.",
            "chrm",
            "scaffold",
            "asm",
            "softmasked",
            "hardmasked",
            "masked",
        )
    ) or has_haplotype_assembly_marker(lower)


def pick_single_file(matches, provider, species_key, label, warnings):
    if len(matches) == 0:
        return None
    if str(label or "").upper() == "GFF":
        ordered = sorted(matches, key=lambda path: gff_candidate_sort_key(provider, path))
    else:
        ordered = sorted(matches, key=lambda path: provider_candidate_sort_key(provider, label, path.name))
    if len(ordered) > 1:
        warnings.append(
            "[{}] {}: multiple {} files found. Using '{}'".format(provider, species_key, label, ordered[0].name)
        )
    return ordered[0]

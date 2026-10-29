"""Shared filename classification, candidate ordering, and network predicates."""

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
from format_species_file_roles import (
    ASSEMBLY_MARKER,
    CDS_MARKER,
    HAPLOTYPE_MARKER,
    PROTEIN_MARKER,
    REPEAT_MARKER,
    TRANSCRIPT_MARKER,
    filename_role_text,
    has_filename_marker,
    has_non_cds_role,
)
from format_species_provider_config import DEFAULT_INPUT_RELATIVE_DIRS, ENSEMBL_LIKE_PROVIDERS, ORYZA_MINUTA_PROVIDER
from format_species_writers import open_text

FERNBASE_GFF_EXCLUDE_PATTERN = re.compile(
    r"(?:^|[._-])(?:" + REPEAT_MARKER + r")(?=[._-]|$)",
    re.IGNORECASE,
)
GFF_EXCLUDE_PATTERN = re.compile(
    r"(?:^|[._-])(?:edta|red|" + REPEAT_MARKER + r")(?=[._-]|$)",
    re.IGNORECASE,
)
GFF_REDUCED_MODEL_PATTERN = re.compile(
    r"(?:^|[._-])(?:only[._-]?long(?:est)?[._-]?transcripts?|primary[._-]?transcripts?|longest)(?=[._-]|$)",
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


def is_probable_cds_filename(provider, name, *, species_key=""):
    if not is_fasta_filename(name) or has_non_cds_role(name, species_key):
        return False
    return has_filename_marker(name, CDS_MARKER + "|" + TRANSCRIPT_MARKER, species_key)


def provider_candidate_sort_key(provider, label, name, *, species_key=""):
    lower = str(name or "").lower()
    role_text = filename_role_text(name, species_key)
    label_upper = str(label or "").upper()

    def marked(pattern):
        return has_filename_marker(role_text, pattern)

    if provider == ORYZA_MINUTA_PROVIDER:
        if label_upper == "CDS":
            return (0 if marked("merged") else 1, 0 if marked(CDS_MARKER) else 1, lower)
        if label_upper == "GFF":
            return (0 if marked("merged") else 1, lower)
        if label_upper == "GENOME":
            return (0 if marked("merged") else 1, 0 if marked(ASSEMBLY_MARKER) else 1, lower)
    if provider in ("citrusgenomedb", "figshare"):
        if label_upper == "CDS":
            return (0 if marked(CDS_MARKER) else 1, int(marked(TRANSCRIPT_MARKER)),
                    int(marked(PROTEIN_MARKER)), lower)
        if label_upper == "GFF":
            if provider == "citrusgenomedb":
                return (int(marked(REPEAT_MARKER + r"|interpro|go|ipr")), lower)
            return (int(GFF_EXCLUDE_PATTERN.search(role_text) is not None),
                    int(GFF_REDUCED_MODEL_PATTERN.search(role_text) is not None), lower)
        if label_upper == "GENOME":
            return (0 if marked(ASSEMBLY_MARKER + "|" + HAPLOTYPE_MARKER) else 1, lower)
    if provider != "fernbase":
        if provider in ENSEMBL_LIKE_PROVIDERS and label_upper == "GFF":
            return (int(marked("abinitio")), int(marked(r"chr|chromosome")),
                    int(marked(r"primary[._-]assembly")), lower)
        if label_upper == "GENOME":
            return (marked("chromosome"), lower)
        if label_upper == "CDS":
            return (0 if marked(CDS_MARKER) else 1, int(marked(TRANSCRIPT_MARKER)),
                    int(marked(PROTEIN_MARKER)), lower)
        if label_upper == "GFF":
            return (int(GFF_EXCLUDE_PATTERN.search(role_text) is not None),
                    int(GFF_REDUCED_MODEL_PATTERN.search(role_text) is not None), lower)
        return (lower,)
    if label_upper == "CDS":
        return (0 if marked(FERNBASE_COMBINED_FILENAME_MARKER) else 1,
                0 if marked("highconfidence") else 1, int(marked("lowconfidence")),
                0 if marked(CDS_MARKER) else 1, int(marked(TRANSCRIPT_MARKER)), lower)
    if label_upper == "GFF":
        return (int(FERNBASE_GFF_EXCLUDE_PATTERN.search(role_text) is not None),
                0 if marked(FERNBASE_COMBINED_FILENAME_MARKER) else 1,
                0 if marked("highconfidence") else 1, int(marked("lowconfidence")), lower)
    if label_upper == "GENOME":
        return (int(marked("chloroplast")), 0 if marked(ASSEMBLY_MARKER) else 1, lower)
    return (lower,)


def gff_contains_cds_feature(path):
    """Return whether a local GFF candidate contains at least one CDS row."""
    try:
        with open_text(path, "rt", errors="replace") as handle:
            for raw_line in handle:
                if raw_line == "" or raw_line.startswith("#"):
                    continue
                parts = raw_line.rstrip("\n\r").split("\t", 3)
                if len(parts) >= 3 and parts[2].strip().lower() == "cds":
                    return True
    except (OSError, UnicodeError):
        return False
    return False


def gff_candidate_sort_key(provider, path, *, species_key=""):
    """Prefer full gene-model GFFs over repeat and reduced annotations."""
    name = str(path.name or "")
    role_text = filename_role_text(name, species_key)
    excluded_by_name = GFF_EXCLUDE_PATTERN.search(role_text) is not None
    reduced_by_name = GFF_REDUCED_MODEL_PATTERN.search(role_text) is not None
    contains_cds = False if excluded_by_name else gff_contains_cds_feature(path)
    return (
        1 if excluded_by_name else 0,
        1 if reduced_by_name else 0,
        0 if contains_cds else 1,
        provider_candidate_sort_key(provider, "GFF", name, species_key=species_key),
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
    return has_filename_marker(name, HAPLOTYPE_MARKER)


def is_probable_genome_filename(provider, name, *, species_key=""):
    if (not is_fasta_filename(name) or has_non_cds_role(name, species_key)
            or is_probable_cds_filename(provider, name, species_key=species_key)):
        return False
    if provider == "fernbase":
        # FernBase often exposes the assembly as a plain ".fa"/".fasta" filename.
        return True
    if provider == "plantgarden":
        return has_filename_marker(name, r"genomic|genomes?|assembly|pseudomolecules?|pmol", species_key)
    if provider in ("ncbi", "refseq", "genbank", "plantaedb"):
        return has_filename_marker(name, "genomic", species_key)
    if provider == "figshare":
        return has_filename_marker(name, ASSEMBLY_MARKER + "|" + HAPLOTYPE_MARKER, species_key)
    if provider in ENSEMBL_LIKE_PROVIDERS:
        return has_filename_marker(name, ASSEMBLY_MARKER, species_key)
    return has_filename_marker(name, ASSEMBLY_MARKER + "|" + HAPLOTYPE_MARKER, species_key)


def pick_single_file(matches, provider, species_key, label, warnings):
    if len(matches) == 0:
        return None
    if str(label or "").upper() == "GFF":
        ordered = sorted(matches, key=lambda path: gff_candidate_sort_key(provider, path, species_key=species_key))
    else:
        ordered = sorted(matches, key=lambda path: provider_candidate_sort_key(provider, label, path.name, species_key=species_key))
    if len(ordered) > 1:
        warnings.append(
            "[{}] {}: multiple {} files found. Using '{}'".format(provider, species_key, label, ordered[0].name)
        )
    return ordered[0]

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
        "onekp",
        "direct",
        "local",
    ):
        return download_root / DEFAULT_INPUT_RELATIVE_DIRS[provider] / species_key
    raise ValueError("Unknown provider: {}".format(provider))

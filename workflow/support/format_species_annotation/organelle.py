"""Identify organellar records that must not enter a nuclear genome input.

NCBI GFF3 files commonly bundle mitochondrial and plastid annotations with the
nuclear assembly.  The records are valid GenBank annotations, but they are not
part of the nuclear gene set used by GeneGalleon.  This module keeps the
filtering rule in one place so CDS, GFF and genome outputs make the same
decision.
"""

import hashlib
import re
from functools import lru_cache
from pathlib import Path

from format_species_writers import apply_common_replacements, open_text
from gff_attribute_syntax import normalise_line

from .common import first_token, parse_gff_attributes
from .source_identity import source_annotation_path

ORGANELLE_SCAN_KEYS = frozenset({"genome", "genome_type", "organelle",
                                 "chromosome", "molecule", "mol_type"})

ORGANELLE_GENOME_TERMS = frozenset(
    {
        "mitochondrion",
        "mitochondria",
        "mitochondrial",
        "chloroplast",
        "chloroplastid",
        "plastid",
        "apicoplast",
        "kinetoplast",
        "chromoplast",
        "cyanelle",
    }
)
ORGANELLE_CHROMOSOME_TERMS = frozenset(
    {
        "m",
        "mt",
        "mtDNA".lower(),
        "mitochondrion",
        "mitochondrial",
        "cp",
        "pt",
        "pltd",
        "plastid",
        "chloroplast",
        "chloroplastid",
    }
)

ORGANELLE_GENOME_LABEL = re.compile(
    r"(?:" + "|".join(sorted(ORGANELLE_GENOME_TERMS)) + r")(?:[ _-]+(?:dna|genome|chromosome))?$")


def _normalized_values(attrs, keys):
    values = []
    for key in keys:
        for raw_value in attrs.get(key, ()):
            value = str(raw_value or "").strip().lower()
            if value != "":
                values.append(value)
    return values


def is_organelle_annotation(feature_type, attrs):
    """Return true only for strong GFF evidence of an organellar record."""
    genome_values = _normalized_values(attrs, ("genome", "genome_type", "organelle"))
    if any(
        ORGANELLE_GENOME_LABEL.fullmatch(value) is not None
        for value in genome_values
    ):
        return True
    if str(feature_type or "").strip().lower() == "region":
        chromosome_values = _normalized_values(attrs, ("chromosome", "molecule", "mol_type"))
        if any(value in ORGANELLE_CHROMOSOME_TERMS for value in chromosome_values):
            return True
    return False


def gff_organelle_seqids(path):
    """Return sequence IDs marked as mitochondrial/plastid in a GFF3 file."""
    source = Path(source_annotation_path(path)).resolve()
    metadata = source.stat()
    return _cached_organelle_seqids(str(source), metadata.st_dev, metadata.st_ino,
                                   metadata.st_size, metadata.st_mtime_ns, metadata.st_ctime_ns,
                                   _annotation_digest(source))


def _annotation_digest(path):
    with open(path, "rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


@lru_cache(maxsize=8)
def _cached_organelle_seqids(path, device, inode, size, mtime, ctime, content_digest):
    seqids = set()
    with open_text(path, "rt", errors="replace") as handle:
        for raw_line in handle:
            if raw_line.startswith("##FASTA"):
                break
            line = apply_common_replacements(raw_line.rstrip("\n\r"))
            if line == "" or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            seqid = str(parts[0] or "").strip()
            if seqid == "":
                continue
            attrs = parse_gff_attributes(parts[8], keys=ORGANELLE_SCAN_KEYS)
            if is_organelle_annotation(parts[2], attrs):
                seqids.add(seqid)
    metadata = Path(path).stat()
    if (metadata.st_dev, metadata.st_ino, metadata.st_size, metadata.st_mtime_ns,
            metadata.st_ctime_ns) != (device, inode, size, mtime, ctime) or (
                _annotation_digest(path) != content_digest):
        raise OSError("GFF changed during organelle detection: " + path)
    return frozenset(sorted(seqids))


def gff_organelle_aliases(path, seqids=None):
    """Return stable feature aliases on organellar sequence IDs."""
    excluded = gff_organelle_seqids(path) if seqids is None else frozenset(seqids)
    aliases = set()
    with open_text(source_annotation_path(path), "rt", errors="replace") as handle:
        for raw_line in handle:
            if raw_line.startswith("##FASTA"):
                break
            line = apply_common_replacements(raw_line.rstrip("\n\r"))
            if not gff_data_line_is_organelle(line, excluded):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            attrs = parse_gff_attributes(parts[8])
            for key in ("ID", "Parent", "transcript_id", "protein_id", "orig_transcript_id", "orig_protein_id"):
                aliases.update(str(value).strip() for value in attrs.get(key, ()) if str(value).strip() != "")
    return frozenset(sorted(aliases))


def gff_data_line_is_organelle(line, seqids):
    if not seqids:
        return False
    parts = str(line).rstrip("\n\r").split("\t")
    return len(parts) >= 9 and str(parts[0] or "").strip() in seqids


def iter_non_organelle_gff_lines(path, seqids=None, *, attribute_changes=None):
    """Yield source GFF lines after removing complete organellar records."""
    excluded = gff_organelle_seqids(path) if seqids is None else frozenset(seqids)
    with open_text(source_annotation_path(path), "rt", errors="replace") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            if raw_line.startswith("##FASTA"):
                break
            line = apply_common_replacements(raw_line)
            if not gff_data_line_is_organelle(line, excluded):
                if attribute_changes is not None or raw_line.rstrip("\r\n").endswith("\t"):
                    raw_line = normalise_line(raw_line, path, line_number,
                                              attribute_changes if attribute_changes is not None else [])
                yield apply_common_replacements(raw_line)


def count_organelle_gff_features(path, seqids=None):
    excluded = gff_organelle_seqids(path) if seqids is None else frozenset(seqids)
    count = 0
    with open_text(source_annotation_path(path), "rt", errors="replace") as handle:
        for raw_line in handle:
            if raw_line.startswith("##FASTA"):
                break
            if gff_data_line_is_organelle(raw_line, excluded):
                count += 1
    return count


def cds_header_seqid(header):
    """Extract an NCBI-style source sequence ID from a CDS FASTA header."""
    token = first_token(str(header or "")).lstrip(">")
    if token.startswith("lcl|"):
        token = token[len("lcl|") :]
    match = re.match(r"^(.+)_cds_[^\s]+(?:_[0-9]+)?$", token)
    return str(match.group(1) or "").strip() if match is not None else ""


def fasta_header_is_organelle(header, seqids, aliases=()):
    seqid = cds_header_seqid(header)
    if seqid != "" and seqid in seqids:
        return True
    if not aliases:
        return False
    token = first_token(str(header or "")).lstrip(">").removeprefix("lcl|")
    if token in aliases:
        return True
    for tag in ("protein_id", "transcript_id", "orig_protein_id", "orig_transcript_id"):
        match = re.search(r"\[{}=([^\]]+)\]".format(re.escape(tag)), str(header or ""))
        if match is not None and match.group(1).strip() in aliases:
            return True
    return False

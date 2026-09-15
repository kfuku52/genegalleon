"""Identify organellar records that must not enter a nuclear genome input.

NCBI GFF3 files commonly bundle mitochondrial and plastid annotations with the
nuclear assembly.  The records are valid GenBank annotations, but they are not
part of the nuclear gene set used by GeneGalleon.  This module keeps the
filtering rule in one place so CDS, GFF and genome outputs make the same
decision.
"""

import re

from format_species_writers import apply_common_replacements, open_text

from .common import first_token, parse_gff_attributes

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
        value in ORGANELLE_GENOME_TERMS or "mitochond" in value or "chloroplast" in value
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
    seqids = set()
    with open_text(path, "rt", errors="replace") as handle:
        for raw_line in handle:
            line = apply_common_replacements(raw_line.rstrip("\n\r"))
            if line == "" or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            seqid = str(parts[0] or "").strip()
            if seqid == "":
                continue
            attrs = parse_gff_attributes(parts[8])
            if is_organelle_annotation(parts[2], attrs):
                seqids.add(seqid)
    return frozenset(sorted(seqids))


def gff_organelle_aliases(path, seqids=None):
    """Return stable feature aliases on organellar sequence IDs."""
    excluded = gff_organelle_seqids(path) if seqids is None else frozenset(seqids)
    aliases = set()
    with open_text(path, "rt", errors="replace") as handle:
        for raw_line in handle:
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


def iter_non_organelle_gff_lines(path, seqids=None):
    """Yield source GFF lines after removing complete organellar records."""
    excluded = gff_organelle_seqids(path) if seqids is None else frozenset(seqids)
    with open_text(path, "rt", errors="replace") as handle:
        for raw_line in handle:
            line = apply_common_replacements(raw_line)
            if not gff_data_line_is_organelle(line, excluded):
                yield line


def count_organelle_gff_features(path, seqids=None):
    excluded = gff_organelle_seqids(path) if seqids is None else frozenset(seqids)
    count = 0
    with open_text(path, "rt", errors="replace") as handle:
        for raw_line in handle:
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

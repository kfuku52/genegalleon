#!/usr/bin/env python3

from pathlib import Path
import re


TAXONOMIC_PROXIMITY_QUALIFIERS = frozenset(("cf", "aff", "nr"))
TAXONOMIC_GENUS_ONLY_PLACEHOLDERS = frozenset(("sp", "spp"))
TAXONOMIC_INFRASPECIFIC_RANKS = frozenset((
    "subsp",
    "var",
    "forma",
    "strain",
    "substrain",
    "serovar",
    "serotype",
    "serogroup",
    "pathovar",
    "biovar",
    "biotype",
    "chemovar",
    "morphovar",
    "cultivar",
    "isolate",
    "group",
    "subgroup",
    "complex",
    "clade",
    "lineage",
    "section",
    "series",
    "ecotype",
    "breed",
))
TAXONOMIC_INFRASPECIFIC_RANK_ALIASES = {
    "subsp": "subsp",
    "ssp": "subsp",
    "subspecies": "subsp",
    "var": "var",
    "variety": "var",
    "forma": "forma",
    "form": "forma",
    "f": "forma",
    "strain": "strain",
    "substrain": "substrain",
    "serovar": "serovar",
    "serotype": "serotype",
    "serogroup": "serogroup",
    "pathovar": "pathovar",
    "pv": "pathovar",
    "biovar": "biovar",
    "biotype": "biotype",
    "chemovar": "chemovar",
    "morphovar": "morphovar",
    "cultivar": "cultivar",
    "cv": "cultivar",
    "isolate": "isolate",
    "group": "group",
    "subgroup": "subgroup",
    "complex": "complex",
    "clade": "clade",
    "lineage": "lineage",
    "section": "section",
    "series": "series",
    "ecotype": "ecotype",
    "breed": "breed",
}
TAXONOMIC_DISPLAY_RANKS = {
    "subsp": "subsp.",
    "var": "var.",
    "forma": "f.",
    "strain": "strain",
    "substrain": "substrain",
    "serovar": "serovar",
    "serotype": "serotype",
    "serogroup": "serogroup",
    "pathovar": "pathovar",
    "biovar": "biovar",
    "biotype": "biotype",
    "chemovar": "chemovar",
    "morphovar": "morphovar",
    "cultivar": "cultivar",
    "isolate": "isolate",
    "group": "group",
    "subgroup": "subgroup",
    "complex": "complex",
    "clade": "clade",
    "lineage": "lineage",
    "section": "section",
    "series": "series",
    "ecotype": "ecotype",
    "breed": "breed",
}

SPECIES_LABEL_TERMINAL_SUFFIXES = (
    ".fasta.busco.full.tsv",
    ".fa.busco.full.tsv",
    ".faa.busco.full.tsv",
    ".fna.busco.full.tsv",
    ".ffn.busco.full.tsv",
    ".fasta.busco.short.txt",
    ".fa.busco.short.txt",
    ".faa.busco.short.txt",
    ".fna.busco.short.txt",
    ".ffn.busco.short.txt",
    ".busco.full.tsv",
    "_busco.full.tsv",
    ".busco.short.txt",
    "_busco.short.txt",
    ".busco.full",
    "_busco.full",
    ".busco.short",
    "_busco.short",
    ".busco",
    "_busco",
    ".longestcds.transcript.fa.gz",
    ".longestcds_contamination_removal.fa.gz",
    ".longestcds.fx2tab_cds.tsv",
    ".longestcds.mmseqs2taxonomy.tsv",
    ".longestcds.fa.gz",
    ".isoform.fa.gz",
    ".fastq.gz",
    ".fq.gz",
    ".fasta.gz",
    ".fa.gz",
    ".faa.gz",
    ".fna.gz",
    ".ffn.gz",
    ".gff3.gz",
    ".gff.gz",
    ".gtf.gz",
    ".tsv.gz",
    ".txt.gz",
    ".csv.gz",
    ".fastq",
    ".fq",
    ".fasta",
    ".fa",
    ".faa",
    ".fna",
    ".ffn",
    ".gff3",
    ".gff",
    ".gtf",
    ".tsv",
    ".txt",
    ".csv",
)


def normalize_taxonomic_name_text(text):
    raw = str(text or "").strip()
    raw = raw.replace("−", "-")
    raw = re.sub(r"\([^)]*\)", " ", raw)
    raw = raw.replace("_", " ")
    return re.sub(r"\s+", " ", raw).strip()


def tokenize_taxonomic_name(text):
    return [token for token in normalize_taxonomic_name_text(text).split(" ") if token]


def taxonomic_token_key(token):
    lowered = str(token or "").strip().rstrip(".").lower()
    if lowered in TAXONOMIC_GENUS_ONLY_PLACEHOLDERS:
        return "sp"
    if lowered in TAXONOMIC_PROXIMITY_QUALIFIERS:
        return lowered
    if lowered in TAXONOMIC_INFRASPECIFIC_RANK_ALIASES:
        return TAXONOMIC_INFRASPECIFIC_RANK_ALIASES[lowered]
    return lowered


def canonical_taxonomic_token(token):
    cleaned = str(token or "").strip()
    key = taxonomic_token_key(cleaned)
    if key == "sp":
        return "sp"
    if key in TAXONOMIC_PROXIMITY_QUALIFIERS:
        return key
    if key in TAXONOMIC_INFRASPECIFIC_RANKS:
        return key
    return cleaned


def species_key_token(token):
    text = str(token or "").strip()
    text = text.replace("/", "_").replace("\\", "_")
    return re.sub(r"\s+", "_", text)


def display_rank_token(rank):
    if rank == "sp":
        return "sp."
    if rank in TAXONOMIC_PROXIMITY_QUALIFIERS:
        return "{}.".format(rank)
    return TAXONOMIC_DISPLAY_RANKS.get(rank, rank)


def species_label_from_taxonomic_text(text):
    tokens = [token for token in tokenize_taxonomic_name(text) if str(token or "").strip() != ""]
    if len(tokens) < 2:
        return ""

    genus = tokens[0][:1].upper() + tokens[0][1:].lower()
    second_key = taxonomic_token_key(tokens[1])
    if second_key in TAXONOMIC_PROXIMITY_QUALIFIERS:
        if len(tokens) >= 3:
            return "{}_{}_{}".format(genus, display_rank_token(second_key), species_key_token(tokens[2].lower()))
        return "{}_{}".format(genus, display_rank_token(second_key))
    if second_key == "sp":
        label = "".join(species_key_token(token) for token in tokens[2:]).strip()
        if label != "":
            return "{}_{}_{}".format(genus, display_rank_token("sp"), label)
        return "{}_{}".format(genus, display_rank_token("sp"))

    species = species_key_token(tokens[1].lower())
    if len(tokens) >= 3:
        third_key = taxonomic_token_key(tokens[2])
        if third_key in TAXONOMIC_PROXIMITY_QUALIFIERS:
            return "{}_{}_{}".format(genus, species, display_rank_token(third_key))
        if third_key in TAXONOMIC_INFRASPECIFIC_RANKS:
            rank = third_key
            value = "".join(species_key_token(token) for token in tokens[3:]).strip()
            if value != "":
                return "{}_{}_{}_{}".format(genus, species, display_rank_token(rank), value)
            return "{}_{}_{}".format(genus, species, display_rank_token(rank))
    return "{}_{}".format(genus, species)


def species_prefix_token_count(parts):
    normalized = [str(part or "").strip() for part in parts if str(part or "").strip() != ""]
    if len(normalized) < 2:
        return 0
    second = taxonomic_token_key(normalized[1])
    third = taxonomic_token_key(normalized[2]) if len(normalized) >= 3 else ""
    if second == "sp":
        return 3 if len(normalized) >= 3 else 2
    if second in TAXONOMIC_PROXIMITY_QUALIFIERS:
        return 3 if len(normalized) >= 3 else 2
    if third in TAXONOMIC_PROXIMITY_QUALIFIERS:
        return 3
    if third in TAXONOMIC_INFRASPECIFIC_RANKS:
        return 4 if len(normalized) >= 4 else 3
    return 2


def strip_species_label_terminal_suffixes(name):
    stem = str(name or "")
    stem_lower = stem.lower()
    for suffix in SPECIES_LABEL_TERMINAL_SUFFIXES:
        if stem_lower.endswith(suffix):
            return stem[:-len(suffix)]
    return stem


def extract_species_label(value, strip_extension=False):
    name = Path(str(value or "")).name
    if strip_extension:
        name = strip_species_label_terminal_suffixes(name)
    parts = [part for part in name.split("_") if part != ""]
    count = species_prefix_token_count(parts)
    if count == 0:
        return ""
    return "_".join(parts[:count])


def strip_species_label(value):
    text = str(value or "")
    species_label = extract_species_label(text, strip_extension=False)
    if species_label != "" and text.startswith(species_label + "_"):
        return text[len(species_label) + 1 :]
    return text


def scientific_name_from_label(value):
    species_label = extract_species_label(value, strip_extension=False)
    if species_label == "":
        species_label = str(value or "").strip()
    parts = [part for part in species_label.split("_") if part != ""]
    if len(parts) >= 3 and taxonomic_token_key(parts[1]) in TAXONOMIC_PROXIMITY_QUALIFIERS:
        return "{} {}. {}".format(parts[0], taxonomic_token_key(parts[1]), parts[2])
    if len(parts) >= 3 and taxonomic_token_key(parts[2]) in TAXONOMIC_PROXIMITY_QUALIFIERS:
        return "{} {}. {}".format(parts[0], taxonomic_token_key(parts[2]), parts[1])
    if len(parts) >= 3 and taxonomic_token_key(parts[1]) == "sp":
        return "{} sp. {}".format(parts[0], parts[2])
    if len(parts) >= 4 and taxonomic_token_key(parts[2]) in TAXONOMIC_INFRASPECIFIC_RANKS:
        rank = taxonomic_token_key(parts[2])
        return "{} {} {} {}".format(parts[0], parts[1], TAXONOMIC_DISPLAY_RANKS.get(rank, rank), parts[3])
    return species_label.replace("_", " ")


def base_species_label(value):
    species_label = extract_species_label(value, strip_extension=False)
    if species_label == "":
        species_label = str(value or "").strip()
    parts = [part for part in species_label.split("_") if part != ""]
    if len(parts) >= 3 and taxonomic_token_key(parts[1]) == "sp":
        return parts[0]
    if len(parts) >= 3 and taxonomic_token_key(parts[1]) in TAXONOMIC_PROXIMITY_QUALIFIERS:
        return "{}_{}".format(parts[0], parts[2])
    if len(parts) >= 3 and taxonomic_token_key(parts[2]) in TAXONOMIC_PROXIMITY_QUALIFIERS:
        return "{}_{}".format(parts[0], parts[1])
    if len(parts) >= 2:
        return "{}_{}".format(parts[0], parts[1])
    return species_label

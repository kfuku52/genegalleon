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


def normalize_taxonomic_name_text(text):
    raw = str(text or "").strip()
    raw = raw.replace("−", "-")
    raw = re.sub(r"\([^)]*\)", " ", raw)
    raw = raw.replace("_", " ")
    return re.sub(r"\s+", " ", raw).strip()


def tokenize_taxonomic_name(text):
    return re.findall(r"[A-Za-z0-9]+", normalize_taxonomic_name_text(text))


def canonical_taxonomic_token(token):
    cleaned = str(token or "").strip()
    lowered = cleaned.lower()
    if lowered in TAXONOMIC_GENUS_ONLY_PLACEHOLDERS:
        return "sp"
    if lowered in TAXONOMIC_PROXIMITY_QUALIFIERS:
        return lowered
    if lowered in TAXONOMIC_INFRASPECIFIC_RANK_ALIASES:
        return TAXONOMIC_INFRASPECIFIC_RANK_ALIASES[lowered]
    return cleaned


def species_label_from_taxonomic_text(text):
    normalized = [canonical_taxonomic_token(token) for token in tokenize_taxonomic_name(text) if str(token or "").strip() != ""]
    if len(normalized) < 2:
        return ""

    genus = normalized[0][:1].upper() + normalized[0][1:].lower()
    second = normalized[1].lower()
    if second in TAXONOMIC_PROXIMITY_QUALIFIERS:
        if len(normalized) >= 3:
            return "{}_{}_{}".format(genus, second, normalized[2].lower())
        return "{}_{}".format(genus, second)
    if second == "sp":
        label = "".join(normalized[2:]).strip()
        if label != "":
            return "{}_sp_{}".format(genus, label)
        return "{}_sp".format(genus)

    species = normalized[1].lower()
    if len(normalized) >= 3:
        third = normalized[2].lower()
        if third in TAXONOMIC_PROXIMITY_QUALIFIERS:
            return "{}_{}_{}".format(genus, species, third)
        if third in TAXONOMIC_INFRASPECIFIC_RANK_ALIASES:
            rank = TAXONOMIC_INFRASPECIFIC_RANK_ALIASES[third]
            value = "".join(normalized[3:]).strip()
            if value != "":
                return "{}_{}_{}_{}".format(genus, species, rank, value)
            return "{}_{}_{}".format(genus, species, rank)
    return "{}_{}".format(genus, species)


def species_prefix_token_count(parts):
    normalized = [str(part or "").strip() for part in parts if str(part or "").strip() != ""]
    if len(normalized) < 2:
        return 0
    second = normalized[1].lower()
    third = normalized[2].lower() if len(normalized) >= 3 else ""
    if second == "sp":
        return 3 if len(normalized) >= 3 else 2
    if second in TAXONOMIC_PROXIMITY_QUALIFIERS:
        return 3 if len(normalized) >= 3 else 2
    if third in TAXONOMIC_PROXIMITY_QUALIFIERS:
        return 3
    if third in TAXONOMIC_INFRASPECIFIC_RANKS:
        return 4 if len(normalized) >= 4 else 3
    return 2


def extract_species_label(value, strip_extension=False):
    name = Path(str(value or "")).name
    if strip_extension and "." in name:
        name = name.split(".", 1)[0]
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
    if len(parts) >= 3 and parts[1].lower() in TAXONOMIC_PROXIMITY_QUALIFIERS:
        return "{} {}. {}".format(parts[0], parts[1].lower(), parts[2])
    if len(parts) >= 3 and parts[2].lower() in TAXONOMIC_PROXIMITY_QUALIFIERS:
        return "{} {}. {}".format(parts[0], parts[2].lower(), parts[1])
    if len(parts) >= 3 and parts[1].lower() == "sp":
        return "{} sp. {}".format(parts[0], parts[2])
    if len(parts) >= 4 and parts[2].lower() in TAXONOMIC_INFRASPECIFIC_RANKS:
        rank = parts[2].lower()
        return "{} {} {} {}".format(parts[0], parts[1], TAXONOMIC_DISPLAY_RANKS.get(rank, rank), parts[3])
    return species_label.replace("_", " ")


def base_species_label(value):
    species_label = extract_species_label(value, strip_extension=False)
    if species_label == "":
        species_label = str(value or "").strip()
    parts = [part for part in species_label.split("_") if part != ""]
    if len(parts) >= 3 and parts[1].lower() == "sp":
        return parts[0]
    if len(parts) >= 3 and parts[1].lower() in TAXONOMIC_PROXIMITY_QUALIFIERS:
        return "{}_{}".format(parts[0], parts[2])
    if len(parts) >= 3 and parts[2].lower() in TAXONOMIC_PROXIMITY_QUALIFIERS:
        return "{}_{}".format(parts[0], parts[1])
    if len(parts) >= 2:
        return "{}_{}".format(parts[0], parts[1])
    return species_label

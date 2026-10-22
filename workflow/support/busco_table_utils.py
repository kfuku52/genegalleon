"""Lightweight BUSCO identifier conventions shared by summary and extraction."""

import re

BUSCO_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
SOURCE_COORDINATE_PATTERN = re.compile(r":[-\.0-9]*$")


def normalize_species_label(path_label):
    normalized = str(path_label)
    if normalized.endswith(".tsv"):
        normalized = normalized[:-4]
    for suffix in (".busco.full", "_busco.full", ".busco", "_busco", ".full", "_full"):
        if normalized.endswith(suffix):
            normalized = normalized[:-len(suffix)]
            break
    normalized = re.sub(r"\s+", "_", normalized.strip())
    return normalized or str(path_label)

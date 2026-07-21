"""Strict taxonomy-backed species aliases for transcriptome quantification."""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from species_labeling import (
    base_species_label,
    normalize_taxonomic_name_text,
    scientific_name_from_label,
    species_label_from_taxonomic_text,
)


class SpeciesAliasError(ValueError):
    """Raised when a species alias cannot be validated safely."""


@dataclass(frozen=True)
class SpeciesAliasResolution:
    canonical_prefix: str
    canonical_scientific_name: str
    metadata_prefix: str
    metadata_scientific_name: str
    declared_taxid: str
    resolved_taxid: str
    method: str


def taxonomic_name_lookup_candidates(value: object) -> list[str]:
    """Return conservative exact-name and normalized species-name candidates."""

    raw_name = normalize_taxonomic_name_text(value)
    species_prefix = species_label_from_taxonomic_text(raw_name)
    candidates = []

    def add(candidate: object) -> None:
        text = str(candidate or "").strip()
        if text and text not in candidates:
            candidates.append(text)

    add(raw_name)
    add(scientific_name_from_label(species_prefix))
    base_prefix = base_species_label(species_prefix)
    if base_prefix and base_prefix != species_prefix:
        add(scientific_name_from_label(base_prefix))
    return candidates


def normalize_species_prefix(value: object) -> str:
    text = str(value or "").strip()
    if text == "":
        return ""
    label = species_label_from_taxonomic_text(text)
    if label != "":
        return label
    return text.replace("/", "_").replace("\\", "_").replace(" ", "_")


def metadata_declared_species_taxid(row: Mapping[str, object]) -> str:
    """Return taxid_species when present, otherwise taxid, as an integer string."""

    raw_value = row.get("taxid_species", "") or row.get("taxid", "")
    value = str(raw_value or "").strip()
    if value == "":
        return ""
    try:
        integer_value = int(value)
    except ValueError as exc:
        raise SpeciesAliasError("Metadata species taxid is not an integer: {!r}".format(value)) from exc
    if integer_value <= 0:
        raise SpeciesAliasError("Metadata species taxid must be positive: {!r}".format(value))
    return str(integer_value)


class SpeciesTaxidResolver:
    """Resolve a taxonomic name uniquely to a species-rank NCBI taxid."""

    def __init__(self, taxonomy_dbfile: object):
        self.db_path = Path(str(taxonomy_dbfile or "")).expanduser()

    def _connect(self) -> sqlite3.Connection:
        if str(self.db_path) == "." or not self.db_path.is_file():
            raise SpeciesAliasError(
                "Taxonomy DB is required to validate a cross-name species alias but was not found: {}".format(
                    self.db_path
                )
            )
        try:
            return sqlite3.connect("file:{}?mode=ro".format(self.db_path), uri=True)
        except sqlite3.Error as exc:
            raise SpeciesAliasError("Could not open taxonomy DB {}: {}".format(self.db_path, exc)) from exc

    def resolve(self, scientific_name: object) -> str:
        name = str(scientific_name or "").strip().replace("_", " ")
        candidates = taxonomic_name_lookup_candidates(name)
        if not candidates:
            raise SpeciesAliasError("Could not construct taxonomy lookup candidates for: {!r}".format(name))

        taxids = set()
        try:
            with self._connect() as conn:
                for candidate in candidates:
                    rows = conn.execute(
                        """
                        SELECT taxid
                        FROM species
                        WHERE spname = ? COLLATE NOCASE AND rank = 'species'
                        UNION
                        SELECT s.taxid
                        FROM synonym AS sy
                        JOIN species AS s ON s.taxid = sy.taxid
                        WHERE sy.spname = ? COLLATE NOCASE AND s.rank = 'species'
                        """,
                        (candidate, candidate),
                    ).fetchall()
                    taxids.update(str(row[0]) for row in rows)
        except sqlite3.Error as exc:
            raise SpeciesAliasError("Failed to query taxonomy DB {}: {}".format(self.db_path, exc)) from exc

        if not taxids:
            raise SpeciesAliasError(
                "Taxonomic name did not resolve to a species-rank taxid: {!r} (candidates: {})".format(
                    name, ", ".join(candidates)
                )
            )
        if len(taxids) != 1:
            raise SpeciesAliasError(
                "Taxonomic name resolved ambiguously to multiple species taxids: {!r} -> {}".format(
                    name, ", ".join(sorted(taxids, key=int))
                )
            )
        return next(iter(taxids))


def validate_species_alias(
    canonical_prefix: object,
    metadata_scientific_name: object,
    taxonomy_dbfile: object = "",
    declared_taxid: object = "",
) -> SpeciesAliasResolution:
    """Validate a metadata name against a GeneGalleon species key.

    Exact and infraspecific aliases retain the legacy label-based behavior.
    Cross-name aliases must resolve uniquely to the same NCBI species taxid.
    """

    canonical_prefix_normalized = normalize_species_prefix(canonical_prefix)
    metadata_name = str(metadata_scientific_name or "").strip()
    metadata_prefix = normalize_species_prefix(metadata_name)
    canonical_name = scientific_name_from_label(canonical_prefix_normalized)
    if canonical_name == "":
        canonical_name = canonical_prefix_normalized.replace("_", " ")
    declared = str(declared_taxid or "").strip()

    if canonical_prefix_normalized == "" or metadata_prefix == "":
        raise SpeciesAliasError(
            "Could not determine species prefixes from canonical={!r}, metadata={!r}".format(
                canonical_prefix, metadata_scientific_name
            )
        )

    canonical_base = base_species_label(canonical_prefix_normalized) or canonical_prefix_normalized
    metadata_base = base_species_label(metadata_prefix) or metadata_prefix
    if metadata_prefix == canonical_prefix_normalized:
        return SpeciesAliasResolution(
            canonical_prefix_normalized,
            canonical_name,
            metadata_prefix,
            metadata_name,
            declared,
            "",
            "exact_species_key",
        )
    if canonical_prefix_normalized == canonical_base and metadata_base == canonical_prefix_normalized:
        return SpeciesAliasResolution(
            canonical_prefix_normalized,
            canonical_name,
            metadata_prefix,
            metadata_name,
            declared,
            "",
            "same_base_species_key",
        )

    resolver = SpeciesTaxidResolver(taxonomy_dbfile)
    canonical_taxid = resolver.resolve(canonical_name)
    metadata_taxid = resolver.resolve(metadata_name)
    if canonical_taxid != metadata_taxid:
        raise SpeciesAliasError(
            "Metadata scientific_name is not an alias of GeneGalleon species_key {!r}: {!r} "
            "resolved to taxid {}, canonical resolved to taxid {}".format(
                canonical_prefix_normalized,
                metadata_name,
                metadata_taxid,
                canonical_taxid,
            )
        )
    if declared and declared != metadata_taxid:
        raise SpeciesAliasError(
            "Metadata species taxid {} does not match resolved taxid {} for {!r}".format(
                declared, metadata_taxid, metadata_name
            )
        )
    return SpeciesAliasResolution(
        canonical_prefix_normalized,
        canonical_name,
        metadata_prefix,
        metadata_name,
        declared,
        metadata_taxid,
        "shared_species_taxid",
    )

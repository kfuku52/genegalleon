"""Taxonomic name normalization and metadata resolution."""

import os
import re
import sqlite3
import tarfile
from pathlib import Path

from format_species_constants import (
    PLASTID_GENETIC_CODE_LINEAGE_DEFAULTS,
    TAXONOMIC_GENUS_ONLY_PLACEHOLDERS,
    TAXONOMIC_HYBRID_CONNECTORS,
    TAXONOMIC_INFRASPECIFIC_RANK_ALIASES,
    TAXONOMIC_PROXIMITY_QUALIFIERS,
)


def blank_species_taxonomy_metadata():
    return {
        "taxid": "",
        "nuclear_genetic_code_id": "",
        "nuclear_genetic_code_name": "",
        "mitochondrial_genetic_code_id": "",
        "mitochondrial_genetic_code_name": "",
        "plastid_genetic_code_id": "",
        "plastid_genetic_code_name": "",
    }


def normalize_taxonomic_name_text(text):
    normalized = str(text or "").replace("\u00d7", " x ")
    return re.sub(r"\s*\(.*?\)\s*", " ", normalized).strip()


def tokenize_taxonomic_name(text):
    return re.findall(r"[A-Za-z0-9]+", normalize_taxonomic_name_text(text))


def canonical_taxonomic_token(token):
    cleaned = str(token or "").strip()
    stripped = cleaned.rstrip(".")
    lowered = stripped.lower()
    if lowered in TAXONOMIC_GENUS_ONLY_PLACEHOLDERS:
        return "sp"
    if lowered in TAXONOMIC_PROXIMITY_QUALIFIERS:
        return lowered
    if lowered in TAXONOMIC_INFRASPECIFIC_RANK_ALIASES:
        return TAXONOMIC_INFRASPECIFIC_RANK_ALIASES[lowered]
    return cleaned


def taxonomic_token_key(token):
    return str(canonical_taxonomic_token(token)).lower()


def is_hybrid_connector_token(token):
    return str(token or "").strip().rstrip(".").lower() in TAXONOMIC_HYBRID_CONNECTORS


def is_hybrid_binomial_connector(tokens, index=2):
    if len(tokens) < index + 3:
        return False
    if not is_hybrid_connector_token(tokens[index]):
        return False
    next_genus = str(tokens[index + 1] or "").strip()
    return bool(next_genus) and next_genus[:1].isupper()


def build_species_key_from_tokens(tokens):
    normalized = [canonical_taxonomic_token(token) for token in tokens if str(token or "").strip() != ""]
    if len(normalized) == 0:
        return ""
    if len(normalized) == 1:
        return normalized[0]

    genus = normalized[0]
    second = normalized[1].lower()
    if is_hybrid_connector_token(normalized[1]):
        if len(normalized) >= 3:
            return "{}_x_{}".format(genus, normalized[2])
        return "{}_x".format(genus)
    if second in TAXONOMIC_PROXIMITY_QUALIFIERS:
        if len(normalized) >= 3:
            return "{}_{}_{}".format(genus, second, normalized[2])
        return "{}_{}".format(genus, normalized[1])
    if second == "sp":
        label = "".join(normalized[2:]).strip()
        if label != "":
            return "{}_sp_{}".format(genus, label)
        return "{}_sp_unknown".format(genus)

    species = normalized[1]
    if len(normalized) >= 3:
        if is_hybrid_binomial_connector(normalized, 2):
            hybrid_genus = normalized[3]
            hybrid_species = normalized[4]
            return "{}_{}_x_{}_{}".format(genus, species, hybrid_genus, hybrid_species)
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


def taxonomic_name_rank_variants(rank):
    if rank == "subsp":
        return ("subsp.", "subsp", "ssp.", "ssp", "subspecies")
    if rank == "var":
        return ("var.", "var", "variety")
    if rank == "forma":
        return ("forma", "forma.", "f.", "f", "form")
    if rank == "pathovar":
        return ("pathovar", "pv.", "pv")
    if rank == "cultivar":
        return ("cultivar", "cv.", "cv")
    if rank == "f":
        return ("f.", "f", "forma", "forma.")
    return (rank,)


def taxonomic_name_lookup_candidates(text):
    species_key = build_species_key_from_tokens(tokenize_taxonomic_name(str(text or "").replace("_", " ")))
    if species_key == "":
        return []

    parts = [part for part in species_key.split("_") if part != ""]
    candidates = []

    def add(value):
        candidate = str(value or "").strip()
        if candidate == "":
            return
        if candidate in candidates:
            return
        candidates.append(candidate)

    add(" ".join(parts))
    if len(parts) >= 3 and is_hybrid_connector_token(parts[1]):
        add("{} x {}".format(parts[0], parts[2]))
        add("{} \u00d7 {}".format(parts[0], parts[2]))
        return candidates
    if len(parts) >= 5 and is_hybrid_connector_token(parts[2]):
        add("{} {} x {} {}".format(parts[0], parts[1], parts[3], parts[4]))
        add("{} {} \u00d7 {} {}".format(parts[0], parts[1], parts[3], parts[4]))
        return candidates
    if len(parts) >= 3 and parts[1].lower() in TAXONOMIC_PROXIMITY_QUALIFIERS:
        add("{} {}. {}".format(parts[0], parts[1], parts[2]))
        add("{} {} {}".format(parts[0], parts[1], parts[2]))
        add("{} {}".format(parts[0], parts[2]))
        return candidates
    if len(parts) >= 3 and parts[2].lower() in TAXONOMIC_PROXIMITY_QUALIFIERS:
        add("{} {}. {}".format(parts[0], parts[2], parts[1]))
        add("{} {} {}".format(parts[0], parts[2], parts[1]))
        add("{} {}".format(parts[0], parts[1]))
        return candidates
    if len(parts) >= 3 and parts[1].lower() in TAXONOMIC_PROXIMITY_QUALIFIERS:
        add("{} {}. {}".format(parts[0], parts[1], parts[2]))
        add("{} {} {}".format(parts[0], parts[1], parts[2]))
        add("{} {}".format(parts[0], parts[2]))
        return candidates
    if len(parts) >= 3 and parts[1].lower() == "sp":
        add("{} sp. {}".format(parts[0], parts[2]))
        add("{} sp {}".format(parts[0], parts[2]))
        add(parts[0])
        return candidates
    if len(parts) >= 4 and parts[2].lower() in TAXONOMIC_INFRASPECIFIC_RANK_ALIASES:
        for rank_variant in taxonomic_name_rank_variants(parts[2].lower()):
            add("{} {} {} {}".format(parts[0], parts[1], rank_variant, parts[3]))
        add("{} {}".format(parts[0], parts[1]))
        return candidates
    if len(parts) >= 2:
        add("{} {}".format(parts[0], parts[1]))
    return candidates


def species_prefix_token_count(tokens):
    if len(tokens) < 2:
        return 0
    second = taxonomic_token_key(tokens[1])
    third = taxonomic_token_key(tokens[2]) if len(tokens) >= 3 else ""
    if is_hybrid_connector_token(tokens[1]):
        return 3 if len(tokens) >= 3 else 2
    if second == "sp":
        return 3 if len(tokens) >= 3 else 2
    if second in TAXONOMIC_PROXIMITY_QUALIFIERS:
        return 3 if len(tokens) >= 3 else 2
    if is_hybrid_binomial_connector(tokens, 2):
        return 5
    if third in TAXONOMIC_PROXIMITY_QUALIFIERS:
        return 3
    if third in TAXONOMIC_INFRASPECIFIC_RANK_ALIASES:
        return 4 if len(tokens) >= 4 else 3
    return 2


def invalid_species_key_reason(value):
    name = Path(str(value or "")).name
    tokens = [token for token in name.split("_") if token != ""]
    if len(tokens) < 2:
        return ""

    second = taxonomic_token_key(tokens[1])
    third = taxonomic_token_key(tokens[2]) if len(tokens) >= 3 else ""
    if second in TAXONOMIC_PROXIMITY_QUALIFIERS and len(tokens) < 3:
        return "taxonomic qualifier '{}' requires a following epithet in species key '{}'".format(second, value)
    if third in TAXONOMIC_INFRASPECIFIC_RANK_ALIASES and len(tokens) < 4:
        return "infraspecific rank '{}' requires a following value in species key '{}'".format(third, value)
    return ""


def normalize_species_key_for_runtime(value):
    species_key = str(value or "").strip()
    if species_key == "":
        return ""
    tokens = [token for token in Path(species_key).name.split("_") if token != ""]
    if len(tokens) >= 2 and tokens[1].lower() == "sp" and len(tokens) < 3:
        return "{}_sp_unknown".format(tokens[0])
    return species_key


def invalid_species_key_error(provider, species_key):
    reason = invalid_species_key_reason(species_key)
    if reason == "":
        return ""
    return (
        "[{}] {}. Incomplete species labels break downstream file naming; "
        "use a concrete key such as 'Dictyostelium_cf_discoideum' or 'Amoeba_sp_sample1'."
    ).format(provider, reason)


class SpeciesTaxonomyMetadataResolver:
    def __init__(self, taxonomy_dbfile="", taxonomy_taxdumpfile=""):
        self.taxonomy_dbfile = str(taxonomy_dbfile or "").strip()
        self.taxonomy_taxdumpfile = str(taxonomy_taxdumpfile or "").strip()
        self._conn = None
        self._nodes = None
        self._genetic_codes = None
        self._cache = {}

    @classmethod
    def from_environment(cls):
        return cls(
            taxonomy_dbfile=os.environ.get("GG_TAXONOMY_DBFILE", ""),
            taxonomy_taxdumpfile=os.environ.get("GG_TAXONOMY_TAXDUMPFILE", ""),
        )

    def _ensure_conn(self):
        if self._conn is not None:
            return self._conn
        if self.taxonomy_dbfile == "":
            return None
        db_path = Path(self.taxonomy_dbfile).expanduser()
        if not db_path.exists():
            return None
        self._conn = sqlite3.connect(str(db_path))
        return self._conn

    def _ensure_taxdump(self):
        if self._nodes is not None and self._genetic_codes is not None:
            return
        self._nodes = {}
        self._genetic_codes = {}
        if self.taxonomy_taxdumpfile == "":
            return
        taxdump_path = Path(self.taxonomy_taxdumpfile).expanduser()
        if not taxdump_path.exists():
            return
        with tarfile.open(taxdump_path, "r:gz") as archive:
            with archive.extractfile("gencode.dmp") as handle:
                for raw in handle:
                    parts = [part.strip() for part in raw.decode("utf-8").split("|")]
                    if len(parts) < 3 or parts[0] == "":
                        continue
                    self._genetic_codes[parts[0]] = {
                        "abbr": parts[1],
                        "name": parts[2],
                    }
            with archive.extractfile("nodes.dmp") as handle:
                for raw in handle:
                    parts = [part.strip() for part in raw.decode("utf-8").split("|")]
                    if len(parts) < 9 or parts[0] == "":
                        continue
                    self._nodes[parts[0]] = {
                        "genetic_code_id": parts[6],
                        "mitochondrial_genetic_code_id": parts[8],
                    }

    def _lookup_species_row(self, species_name):
        conn = self._ensure_conn()
        if conn is None:
            return None
        cur = conn.cursor()
        for candidate in taxonomic_name_lookup_candidates(species_name):
            row = cur.execute(
                "SELECT taxid, spname, rank, track FROM species WHERE spname = ? COLLATE NOCASE",
                (candidate,),
            ).fetchone()
            if row is not None:
                return row
            row = cur.execute(
                """
                SELECT s.taxid, s.spname, s.rank, s.track
                FROM synonym sy
                JOIN species s ON sy.taxid = s.taxid
                WHERE sy.spname = ? COLLATE NOCASE
                """,
                (candidate,),
            ).fetchone()
            if row is not None:
                return row
        return None

    def _code_name(self, code_id):
        if code_id == "":
            return ""
        self._ensure_taxdump()
        return str(self._genetic_codes.get(str(code_id), {}).get("name", ""))

    def _infer_plastid_code_id(self, lineage_taxids):
        lineage_taxid_set = {str(taxid).strip() for taxid in lineage_taxids if str(taxid).strip() != ""}
        for lineage_taxid, code_id in PLASTID_GENETIC_CODE_LINEAGE_DEFAULTS.items():
            if lineage_taxid in lineage_taxid_set:
                return code_id
        return ""

    def resolve(self, species_name):
        normalized_name = str(species_name or "").strip().replace("_", " ")
        if normalized_name == "":
            return blank_species_taxonomy_metadata()
        if normalized_name in self._cache:
            return dict(self._cache[normalized_name])

        metadata = blank_species_taxonomy_metadata()
        row = self._lookup_species_row(normalized_name)
        if row is None:
            self._cache[normalized_name] = dict(metadata)
            return metadata

        taxid, _spname, _rank, track = row
        taxid_str = str(taxid)
        metadata["taxid"] = taxid_str

        self._ensure_taxdump()
        node = self._nodes.get(taxid_str, {})
        nuclear_code_id = str(node.get("genetic_code_id", "") or "")
        mitochondrial_code_id = str(node.get("mitochondrial_genetic_code_id", "") or "")
        metadata["nuclear_genetic_code_id"] = nuclear_code_id
        metadata["nuclear_genetic_code_name"] = self._code_name(nuclear_code_id)
        metadata["mitochondrial_genetic_code_id"] = mitochondrial_code_id
        metadata["mitochondrial_genetic_code_name"] = self._code_name(mitochondrial_code_id)

        lineage_taxids = [part.strip() for part in str(track or "").split(",") if part.strip() != ""]
        plastid_code_id = self._infer_plastid_code_id(lineage_taxids)
        metadata["plastid_genetic_code_id"] = plastid_code_id
        metadata["plastid_genetic_code_name"] = self._code_name(plastid_code_id)

        self._cache[normalized_name] = dict(metadata)
        return metadata

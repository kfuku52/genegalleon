"""Stable gene-identity provenance for GFF-backed CDS grouping."""

import re
from collections import defaultdict

from format_species_constants import ENSEMBL_GENE_ID_PATTERN

from .common import choose_first_gff_attribute

KNOWN_ALIAS_PREFIXES = (
    "cds-",
    "cds:",
    "gene-",
    "gene:",
    "mrna-",
    "mrna:",
    "protein-",
    "protein:",
    "rna-",
    "rna:",
    "transcript-",
    "transcript:",
)


def gff_dbxref_gene_token(attrs):
    dbxrefs = []
    for key in ("Dbxref", "db_xref", "dbxref"):
        dbxrefs.extend(attrs.get(key, ()))
    for raw_value in dbxrefs:
        value = str(raw_value or "").strip()
        if value.lower().startswith("ensembl:"):
            candidate = value.split(":", 1)[1].strip()
            if ENSEMBL_GENE_ID_PATTERN.match(candidate):
                return candidate
    for namespace in ("Araport", "TAIR"):
        prefix = namespace.lower() + ":"
        for raw_value in dbxrefs:
            value = str(raw_value or "").strip()
            if not value.lower().startswith(prefix):
                continue
            candidate = value[len(prefix) :].strip()
            match = re.match(
                r"^(AT(?:[1-5]|C|M)G[0-9]+)(?:\.[0-9]+)?$",
                candidate,
                flags=re.IGNORECASE,
            )
            if match is not None:
                return match.group(1)
    for raw_value in dbxrefs:
        value = str(raw_value or "").strip()
        if value.lower().startswith("geneid:"):
            candidate = value.split(":", 1)[1].strip()
            if candidate != "":
                return "GeneID{}".format(candidate)
    return ""


def strip_gff_feature_prefix(value):
    current = str(value or "").strip()
    while current != "":
        lowered = current.lower()
        stripped = ""
        for prefix in KNOWN_ALIAS_PREFIXES:
            if lowered.startswith(prefix):
                stripped = current[len(prefix) :].strip()
                break
        if stripped == "" or stripped == current:
            break
        current = stripped
    return current


def gff_authoritative_gene_token(attrs):
    direct = choose_first_gff_attribute(attrs, ("locus_tag", "gene_id"))
    if direct != "":
        return direct
    dbxref_token = gff_dbxref_gene_token(attrs)
    if dbxref_token != "":
        return dbxref_token
    return choose_first_gff_attribute(attrs, ("Parent_Accession", "Accession"))


def gff_stable_gene_token(attrs, feature_type, include_symbol=True):
    authoritative = gff_authoritative_gene_token(attrs)
    if authoritative != "":
        return authoritative
    if str(feature_type or "").strip().lower() == "gene":
        direct = choose_first_gff_attribute(attrs, ("ID", "Name"))
        if direct != "":
            return strip_gff_feature_prefix(direct)
    if include_symbol:
        return choose_first_gff_attribute(attrs, ("gene", "geneName"))
    return ""


def resolve_grouping_feature_authoritative_gene_tokens(feature_id, feature_records, cache, active):
    feature_text = str(feature_id or "").strip()
    if feature_text == "" or feature_text in active:
        return ()
    if feature_text in cache:
        return cache[feature_text]
    record = feature_records.get(feature_text)
    if record is None:
        cache[feature_text] = ()
        return cache[feature_text]

    active.add(feature_text)
    tokens = set(record.get("authoritative_gene_tokens", ()))
    for parent_id in record.get("parents", ()):
        tokens.update(
            resolve_grouping_feature_authoritative_gene_tokens(parent_id, feature_records, cache, active)
        )
    active.remove(feature_text)
    cache[feature_text] = tuple(sorted(token for token in tokens if str(token or "").strip() != ""))
    return cache[feature_text]


def resolve_grouping_feature_gene_feature_ids(feature_id, feature_records, cache, active):
    feature_text = str(feature_id or "").strip()
    if feature_text == "" or feature_text in active:
        return ()
    if feature_text in cache:
        return cache[feature_text]
    record = feature_records.get(feature_text)
    if record is None:
        cache[feature_text] = ()
        return cache[feature_text]
    if record.get("feature_type") == "gene":
        cache[feature_text] = (feature_text,)
        return cache[feature_text]

    active.add(feature_text)
    gene_feature_ids = set()
    for parent_id in record.get("parents", ()):
        gene_feature_ids.update(
            resolve_grouping_feature_gene_feature_ids(parent_id, feature_records, cache, active)
        )
    active.remove(feature_text)
    cache[feature_text] = tuple(sorted(gene_feature_ids))
    return cache[feature_text]


def reject_gff_gene_prefix_normalization_collisions(
    task,
    original_gene_tokens,
    gene_feature_ids_by_transcript,
    authoritative_gene_tokens_by_transcript,
):
    gene_feature_ids_by_token = defaultdict(set)
    authoritatively_confirmed_by_token = defaultdict(lambda: True)
    for transcript_id, gene_token in original_gene_tokens.items():
        matching_gene_feature_ids = tuple(
            gene_feature_id
            for gene_feature_id in gene_feature_ids_by_transcript.get(transcript_id, ())
            if strip_gff_feature_prefix(gene_feature_id) == gene_token
        )
        if len(matching_gene_feature_ids) == 0:
            continue
        gene_feature_ids_by_token[gene_token].update(matching_gene_feature_ids)
        authoritative_tokens = set(authoritative_gene_tokens_by_transcript.get(transcript_id, ()))
        if gene_token not in authoritative_tokens:
            authoritatively_confirmed_by_token[gene_token] = False

    prefix_collisions = {
        gene_token: tuple(sorted(gene_feature_ids))
        for gene_token, gene_feature_ids in gene_feature_ids_by_token.items()
        if len(gene_feature_ids) > 1 and not authoritatively_confirmed_by_token[gene_token]
    }
    if len(prefix_collisions) == 0:
        return

    examples = []
    for gene_token in sorted(prefix_collisions)[:5]:
        examples.append("{} <- {}".format(gene_token, ",".join(prefix_collisions[gene_token])))
    raise ValueError(
        "GFF-backed CDS grouping for {} has distinct gene feature IDs that collapse after "
        "prefix normalization: {}".format(task.get("species_prefix", ""), "; ".join(examples))
    )

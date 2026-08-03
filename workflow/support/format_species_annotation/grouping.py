"""Annotation implementation: grouping."""

import re
from collections import defaultdict

try:
    from Bio import SeqIO
except Exception:  # pragma: no cover - runtime without biopython
    SeqIO = None

from format_species_constants import (
    ENSEMBL_GENE_ID_PATTERN,
    RESCUE_SAME_TERMINAL_MIN_LONGER_OVERLAP,
    RESCUE_SAME_TERMINAL_MIN_SHORTER_OVERLAP,
    RESCUE_SHARED_JUNCTION_MIN_LONGER_OVERLAP,
    RESCUE_SHARED_JUNCTION_MIN_SHORTER_OVERLAP,
)
from format_species_writers import open_text

from .common import (
    choose_first_gff_attribute,
    collapse_transcript_suffix,
    compute_interval_overlap_bases,
    extract_header_tag_value,
    extract_provider_id,
    extract_provider_transcript_id,
    first_token,
    gene_grouping_mode_for_task,
    merge_coordinate_intervals,
    parse_gff_attributes,
    transcript_feature_gene_token,
)

GFF_CDS_ALIAS_KEYS = frozenset(
    (
        "accession",
        "alias",
        "cds",
        "dbxref",
        "db_xref",
        "gene",
        "gene_id",
        "genename",
        "id",
        "locus_tag",
        "mrna",
        "name",
        "orig_protein_id",
        "orig_transcript_id",
        "parent",
        "parent_accession",
        "protein_id",
        "transcript_id",
    )
)

CDS_HEADER_PRIMARY_TAGS = (
    "protein_id",
    "transcript_id",
    "orig_protein_id",
    "orig_transcript_id",
)

CDS_HEADER_GENE_TAGS = (
    "locus_tag",
    "gene_id",
    "gene",
)

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

GFF_GROUPING_PARENT_FEATURE_TYPES = frozenset(
    (
        "gene",
        "lncrna",
        "mrna",
        "ncrna",
        "polypeptide",
        "primary_transcript",
        "protein",
        "pseudogene",
        "rna",
        "rrna",
        "transcript",
        "trna",
    )
)


def gff_alias_variants(value):
    text = str(value or "").strip()
    if text == "":
        return ()
    variants = []

    def add(candidate, normalize_separators=True):
        normalized = str(candidate or "").strip()
        if normalized != "" and normalized not in variants:
            variants.append(normalized)
        if normalize_separators:
            separator_normalized = re.sub(r"[-_]+", "_", normalized)
            if separator_normalized != "" and separator_normalized not in variants:
                variants.append(separator_normalized)

    has_known_prefix = any(text.lower().startswith(prefix) for prefix in KNOWN_ALIAS_PREFIXES)
    add(text, normalize_separators=not has_known_prefix)
    current = text
    while True:
        lowered = current.lower()
        stripped = ""
        for prefix in KNOWN_ALIAS_PREFIXES:
            if lowered.startswith(prefix):
                stripped = current[len(prefix) :].strip()
                break
        if stripped == "" or stripped == current:
            break
        add(stripped)
        current = stripped
    if "|" in text:
        add(text.rsplit("|", 1)[-1])
    if ":" in text:
        namespace, candidate = text.split(":", 1)
        if namespace.lower() in (
            "ensembl",
            "genbank",
            "geneid",
            "ncbi_gene",
            "ncbi_gp",
            "refseq",
        ):
            add(candidate)
    return tuple(variants)


def gff_alias_values_from_attributes(attrs):
    aliases = []
    for key, values in attrs.items():
        if str(key or "").strip().lower() not in GFF_CDS_ALIAS_KEYS:
            continue
        for value in values:
            alias = str(value or "").strip()
            if alias != "" and alias not in aliases:
                aliases.append(alias)
    return tuple(aliases)


def gff_aliases_from_attributes(attrs):
    aliases = []
    for value in gff_alias_values_from_attributes(attrs):
        for alias in gff_alias_variants(value):
            if alias not in aliases:
                aliases.append(alias)
    return tuple(aliases)


def gff_grouping_parent_feature_type(feature_type):
    feature_type_lower = str(feature_type or "").strip().lower()
    return (
        feature_type_lower in GFF_GROUPING_PARENT_FEATURE_TYPES
        or "transcript" in feature_type_lower
        or feature_type_lower.endswith("gene")
        or feature_type_lower.endswith("rna")
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


def gff_stable_gene_token(attrs, feature_type, include_symbol=True):
    direct = choose_first_gff_attribute(attrs, ("locus_tag", "gene_id"))
    if direct != "":
        return direct
    dbxref_token = gff_dbxref_gene_token(attrs)
    if dbxref_token != "":
        return dbxref_token
    direct = choose_first_gff_attribute(attrs, ("Parent_Accession", "Accession"))
    if direct != "":
        return direct
    if str(feature_type or "").strip().lower() == "gene":
        direct = choose_first_gff_attribute(attrs, ("ID", "Name"))
        if direct != "":
            return strip_gff_feature_prefix(direct)
    if include_symbol:
        return choose_first_gff_attribute(attrs, ("gene", "geneName"))
    return ""


def resolve_grouping_feature_gene_token(feature_id, feature_records, provider, cache, active):
    feature_text = str(feature_id or "").strip()
    if feature_text == "":
        return ""
    if feature_text in cache:
        return cache[feature_text]
    if feature_text in active:
        collapsed = collapse_transcript_suffix(provider, feature_text)
        return collapsed if collapsed != "" else feature_text

    record = feature_records.get(feature_text)
    if record is None:
        collapsed = collapse_transcript_suffix(provider, feature_text)
        cache[feature_text] = collapsed if collapsed != "" else strip_gff_feature_prefix(feature_text)
        return cache[feature_text]

    active.add(feature_text)
    resolved = ""
    if record.get("feature_type") != "gene":
        for parent_id in record.get("parents", ()):
            resolved = resolve_grouping_feature_gene_token(parent_id, feature_records, provider, cache, active)
            if resolved != "":
                break
    if resolved == "":
        resolved = str(record.get("gene_token", "") or "").strip()
    if resolved == "" and record.get("feature_type") == "gene":
        for parent_id in record.get("parents", ()):
            resolved = resolve_grouping_feature_gene_token(parent_id, feature_records, provider, cache, active)
            if resolved != "":
                break
    if resolved == "":
        resolved = collapse_transcript_suffix(provider, feature_text)
    active.remove(feature_text)
    cache[feature_text] = str(resolved or "").strip()
    return cache[feature_text]


def build_gff_cds_grouping_index(task):
    """Build a GFF-backed alias-to-gene index for provider CDS FASTA records."""
    gff_path = task.get("gff_path")
    if gff_path is None:
        return None

    feature_records = {}
    child_features_by_parent = defaultdict(set)
    aliases_by_transcript = defaultdict(set)
    fallback_gene_tokens = {}
    use_coordinate_rescue = gene_grouping_mode_for_task(task) == "rescue_overlap"
    cds_features_by_transcript = defaultdict(list) if use_coordinate_rescue else None

    with open_text(gff_path, "rt") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            if line == "" or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            seqid, _source, feature_type, start_text, end_text, _score, strand, _phase, attr_text = parts[:9]
            feature_type_lower = str(feature_type or "").strip().lower()
            attrs = parse_gff_attributes(attr_text)
            feature_id = choose_first_gff_attribute(attrs, ("ID", "transcript_id", "protein_id", "Name"))
            direct_parents = tuple(
                str(value or "").strip()
                for value in attrs.get("Parent", ())
                if str(value or "").strip() != ""
            )
            parents = list(direct_parents)
            for relationship_key in ("Parent", "Derives_from", "derives_from"):
                for value in attrs.get(relationship_key, ()):
                    parent_text = str(value or "").strip()
                    if parent_text != "" and parent_text not in parents:
                        parents.append(parent_text)
            parents = tuple(parents)
            if feature_id != "" and gff_grouping_parent_feature_type(feature_type_lower):
                feature_records[feature_id] = {
                    "feature_type": feature_type_lower,
                    "parents": parents,
                    "gene_token": gff_stable_gene_token(attrs, feature_type_lower),
                    "aliases": gff_alias_values_from_attributes(attrs),
                }
                if feature_type_lower in ("polypeptide", "protein") or any(
                    len(attrs.get(key, ())) > 0 for key in ("Derives_from", "derives_from")
                ):
                    for parent_id in parents:
                        child_features_by_parent[parent_id].add(feature_id)
            if feature_type_lower != "cds":
                continue
            try:
                start = int(start_text)
                end = int(end_text)
            except Exception:
                continue
            if start > end:
                start, end = end, start
            transcript_ids = list(direct_parents if len(direct_parents) > 0 else parents)
            if len(transcript_ids) == 0:
                transcript_id = choose_first_gff_attribute(attrs, ("transcript_id", "protein_id", "ID", "Name"))
                if transcript_id == "":
                    transcript_id = "{}:{}-{}".format(seqid, start, end)
                transcript_ids = [transcript_id]
            fallback_gene = gff_stable_gene_token(attrs, feature_type_lower)
            attribute_aliases = gff_alias_values_from_attributes(attrs)
            for transcript_id in transcript_ids:
                transcript_text = str(transcript_id or "").strip()
                if transcript_text == "":
                    continue
                for alias in gff_alias_variants(transcript_text):
                    aliases_by_transcript[transcript_text].add(alias)
                for raw_alias in attribute_aliases:
                    aliases_by_transcript[transcript_text].update(gff_alias_variants(raw_alias))
                if fallback_gene != "":
                    previous_fallback = fallback_gene_tokens.get(transcript_text, "")
                    if previous_fallback == "" or fallback_gene < previous_fallback:
                        fallback_gene_tokens[transcript_text] = fallback_gene
                if use_coordinate_rescue:
                    cds_features_by_transcript[transcript_text].append(
                        {
                            "seqid": str(seqid or "").strip(),
                            "start": start,
                            "end": end,
                            "strand": str(strand or "").strip() or "+",
                            "gene_token": fallback_gene,
                        }
                    )

    gene_cache = {}
    if use_coordinate_rescue:
        for transcript_id, features in cds_features_by_transcript.items():
            resolved_gene = ""
            if transcript_id in feature_records:
                resolved_gene = resolve_grouping_feature_gene_token(
                    transcript_id,
                    feature_records,
                    task["provider"],
                    gene_cache,
                    set(),
                )
            if resolved_gene == "":
                resolved_gene = fallback_gene_tokens.get(transcript_id, "")
            for feature in features:
                feature["gene_token"] = resolved_gene
        original_gene_tokens = {
            transcript_id: transcript_feature_gene_token(features)
            for transcript_id, features in cds_features_by_transcript.items()
        }
        resolved_gene_tokens = build_rescued_gene_tokens_for_transcripts(task, cds_features_by_transcript)
    else:
        original_gene_tokens = {}
        for transcript_id in aliases_by_transcript:
            resolved_gene = ""
            if transcript_id in feature_records:
                resolved_gene = resolve_grouping_feature_gene_token(
                    transcript_id,
                    feature_records,
                    task["provider"],
                    gene_cache,
                    set(),
                )
            if resolved_gene == "":
                resolved_gene = fallback_gene_tokens.get(transcript_id, "")
            if resolved_gene == "":
                resolved_gene = collapse_transcript_suffix(task["provider"], transcript_id) or transcript_id
            original_gene_tokens[transcript_id] = resolved_gene
        resolved_gene_tokens = dict(original_gene_tokens)

    for transcript_id in aliases_by_transcript:
        pending = [transcript_id]
        visited = set()
        while len(pending) > 0:
            feature_id = pending.pop()
            if feature_id in visited:
                continue
            visited.add(feature_id)
            for alias in gff_alias_variants(feature_id):
                aliases_by_transcript[transcript_id].add(alias)
            record = feature_records.get(feature_id)
            if record is None:
                continue
            for raw_alias in record.get("aliases", ()):
                aliases_by_transcript[transcript_id].update(gff_alias_variants(raw_alias))
            for parent_id in record.get("parents", ()):
                parent_text = str(parent_id or "").strip()
                if parent_text != "" and parent_text not in visited:
                    pending.append(parent_text)

        pending = list(child_features_by_parent.get(transcript_id, ()))
        visited = set()
        while len(pending) > 0:
            feature_id = pending.pop()
            if feature_id in visited:
                continue
            visited.add(feature_id)
            for alias in gff_alias_variants(feature_id):
                aliases_by_transcript[transcript_id].add(alias)
            record = feature_records.get(feature_id)
            if record is not None:
                for raw_alias in record.get("aliases", ()):
                    aliases_by_transcript[transcript_id].update(gff_alias_variants(raw_alias))
            for child_id in child_features_by_parent.get(feature_id, ()):
                if child_id not in visited:
                    pending.append(child_id)

    feature_records.clear()
    child_features_by_parent.clear()
    gene_cache.clear()
    fallback_gene_tokens.clear()
    if cds_features_by_transcript is not None:
        cds_features_by_transcript.clear()

    alias_to_gene_tokens = defaultdict(set)
    for transcript_id, aliases in aliases_by_transcript.items():
        gene_token = str(resolved_gene_tokens.get(transcript_id, "") or "").strip()
        if gene_token == "":
            gene_token = collapse_transcript_suffix(task["provider"], transcript_id) or transcript_id
        for alias in aliases:
            alias_to_gene_tokens[alias].add(gene_token)

    rescued_transcripts = sorted(
        transcript_id
        for transcript_id in resolved_gene_tokens
        if str(original_gene_tokens.get(transcript_id, "") or "").strip()
        != str(resolved_gene_tokens.get(transcript_id, "") or "").strip()
    )
    return {
        "gff_path": str(gff_path),
        "alias_to_gene_tokens": {
            alias: tuple(sorted(gene_tokens)) for alias, gene_tokens in alias_to_gene_tokens.items()
        },
        "transcript_gene_tokens": dict(resolved_gene_tokens),
        "transcripts_total": len(aliases_by_transcript),
        "coordinate_rescued_transcripts": len(rescued_transcripts),
        "coordinate_rescued_groups": len(
            {
                resolved_gene_tokens[transcript_id]
                for transcript_id in rescued_transcripts
                if str(resolved_gene_tokens.get(transcript_id, "") or "").strip() != ""
            }
        ),
    }


def extract_cds_header_alias_tiers(task, header):
    primary = []
    gene_level = []

    def add(target, value):
        for alias in gff_alias_variants(value):
            if alias not in target:
                target.append(alias)
            species_prefix = str(task.get("species_prefix", "") or "").strip()
            for separator in (".", "_", "-"):
                prefixed = species_prefix + separator
                if species_prefix != "" and alias.startswith(prefixed):
                    stripped = alias[len(prefixed) :].strip()
                    if stripped != "" and stripped not in target:
                        target.append(stripped)

    add(primary, first_token(header))
    add(primary, extract_provider_transcript_id(task["provider"], header))
    add(primary, extract_provider_id(task["provider"], header))
    for tag in CDS_HEADER_PRIMARY_TAGS:
        add(primary, extract_header_tag_value(header, tag))
    for tag in CDS_HEADER_GENE_TAGS:
        add(gene_level, extract_header_tag_value(header, tag))

    for match in re.finditer(
        r"(?:^|\s)(transcript|mrna|rna|protein|gene|gene_id|locus_tag):([^\s\]]+)",
        str(header or ""),
        flags=re.IGNORECASE,
    ):
        label = match.group(1).lower()
        target = gene_level if label in ("gene", "gene_id", "locus_tag") else primary
        add(target, match.group(2))

    for dbxref_tag in ("db_xref", "Dbxref"):
        raw_dbxref = extract_header_tag_value(header, dbxref_tag)
        for token in str(raw_dbxref or "").split(","):
            add(primary, token)

    raw_token = first_token(header)
    ncbi_cds_match = re.search(r"_cds_([^\s]+?)(?:_[0-9]+)?$", raw_token)
    if ncbi_cds_match is not None:
        add(primary, ncbi_cds_match.group(1))
    return (tuple(primary), tuple(gene_level))


def resolve_cds_header_gff_gene(task, header, grouping_index=None):
    index = grouping_index if grouping_index is not None else task.get("_gff_cds_grouping_index")
    if index is None:
        return {
            "status": "not_applicable",
            "gene_token": "",
            "matched_aliases": (),
            "candidate_gene_tokens": (),
        }
    alias_index = index.get("alias_to_gene_tokens", {})
    evidence = []
    seen_aliases = set()
    unique_gene_tokens = set()
    candidate_gene_tokens = set()
    for aliases in extract_cds_header_alias_tiers(task, header):
        for alias in aliases:
            if alias in seen_aliases:
                continue
            seen_aliases.add(alias)
            candidates = tuple(alias_index.get(alias, ()))
            if len(candidates) > 0:
                evidence.append((alias, candidates))
                candidate_gene_tokens.update(candidates)
            if len(candidates) == 1:
                unique_gene_tokens.add(candidates[0])
    if len(unique_gene_tokens) == 1:
        selected = next(iter(unique_gene_tokens))
        if all(selected in candidates for _alias, candidates in evidence):
            return {
                "status": "mapped",
                "gene_token": selected,
                "matched_aliases": tuple(alias for alias, candidates in evidence if selected in candidates),
                "candidate_gene_tokens": (selected,),
            }
    if len(candidate_gene_tokens) > 0:
        return {
            "status": "ambiguous",
            "gene_token": "",
            "matched_aliases": tuple(alias for alias, _candidates in evidence),
            "candidate_gene_tokens": tuple(sorted(candidate_gene_tokens)),
        }
    return {
        "status": "unmapped",
        "gene_token": "",
        "matched_aliases": (),
        "candidate_gene_tokens": (),
    }


def build_transcript_grouping_entry(provider, transcript_id, features):
    if len(features) == 0:
        return {
            "transcript_id": transcript_id,
            "gene_token": "",
            "seqid": "",
            "strand": "+",
            "merged_parts": (),
            "junctions": frozenset(),
            "cds_length": 0,
            "span_start": 0,
            "span_end": 0,
            "terminal_coord": 0,
            "rescue_ineligible": True,
            "collapsed_transcript_id": collapse_transcript_suffix(provider, transcript_id),
        }
    seqid = str(features[0].get("seqid", "") or "").strip()
    strand = str(features[0].get("strand", "") or "").strip() or "+"
    rescue_ineligible = False
    intervals = []
    for feature in features:
        feature_seqid = str(feature.get("seqid", "") or "").strip()
        feature_strand = str(feature.get("strand", "") or "").strip() or "+"
        if feature_seqid != seqid or feature_strand != strand:
            rescue_ineligible = True
        intervals.append((int(feature["start"]), int(feature["end"])))
    merged_parts = merge_coordinate_intervals(intervals)
    if len(merged_parts) == 0:
        rescue_ineligible = True
        span_start = 0
        span_end = 0
        terminal_coord = 0
        cds_length = 0
        junctions = frozenset()
    else:
        span_start = merged_parts[0][0]
        span_end = merged_parts[-1][1]
        terminal_coord = span_end if strand == "+" else span_start
        cds_length = sum(end - start + 1 for start, end in merged_parts)
        junctions = frozenset(
            (merged_parts[index][1], merged_parts[index + 1][0]) for index in range(len(merged_parts) - 1)
        )
    return {
        "transcript_id": transcript_id,
        "gene_token": transcript_feature_gene_token(features),
        "seqid": seqid,
        "strand": strand,
        "merged_parts": tuple(merged_parts),
        "junctions": junctions,
        "cds_length": cds_length,
        "span_start": span_start,
        "span_end": span_end,
        "terminal_coord": terminal_coord,
        "rescue_ineligible": rescue_ineligible,
        "collapsed_transcript_id": collapse_transcript_suffix(provider, transcript_id),
    }


def should_rescue_overlapping_transcripts(left, right):
    if left["rescue_ineligible"] or right["rescue_ineligible"]:
        return False
    if left["seqid"] != right["seqid"] or left["strand"] != right["strand"]:
        return False
    if left["cds_length"] <= 0 or right["cds_length"] <= 0:
        return False
    if left["span_end"] < right["span_start"] or right["span_end"] < left["span_start"]:
        return False
    overlap_bases = compute_interval_overlap_bases(left["merged_parts"], right["merged_parts"])
    if overlap_bases <= 0:
        return False
    shorter_overlap = overlap_bases / float(min(left["cds_length"], right["cds_length"]))
    longer_overlap = overlap_bases / float(max(left["cds_length"], right["cds_length"]))
    shared_junctions = len(left["junctions"].intersection(right["junctions"]))
    if (
        shared_junctions > 0
        and shorter_overlap >= RESCUE_SHARED_JUNCTION_MIN_SHORTER_OVERLAP
        and longer_overlap >= RESCUE_SHARED_JUNCTION_MIN_LONGER_OVERLAP
    ):
        return True
    if (
        left["terminal_coord"] == right["terminal_coord"]
        and shorter_overlap >= RESCUE_SAME_TERMINAL_MIN_SHORTER_OVERLAP
        and longer_overlap >= RESCUE_SAME_TERMINAL_MIN_LONGER_OVERLAP
    ):
        return True
    return False


def choose_rescue_cluster_gene_token(provider, entries):
    collapsed_transcript_ids = sorted(
        {
            str(entry.get("collapsed_transcript_id", "") or "").strip()
            for entry in entries
            if str(entry.get("collapsed_transcript_id", "") or "").strip() != ""
        }
    )
    if len(collapsed_transcript_ids) == 1:
        return collapsed_transcript_ids[0]
    gene_tokens = sorted(
        {
            str(entry.get("gene_token", "") or "").strip()
            for entry in entries
            if str(entry.get("gene_token", "") or "").strip() != ""
        }
    )
    if len(gene_tokens) > 0:
        return gene_tokens[0]
    if len(collapsed_transcript_ids) > 0:
        return collapsed_transcript_ids[0]
    transcript_ids = sorted(
        str(entry.get("transcript_id", "") or "").strip()
        for entry in entries
        if str(entry.get("transcript_id", "") or "").strip() != ""
    )
    if len(transcript_ids) > 0:
        fallback = collapse_transcript_suffix(provider, transcript_ids[0])
        return fallback if fallback != "" else transcript_ids[0]
    return ""


def build_rescued_gene_tokens_for_transcripts(task, cds_features_by_transcript):
    resolved = {
        transcript_id: transcript_feature_gene_token(features)
        for transcript_id, features in cds_features_by_transcript.items()
    }
    if gene_grouping_mode_for_task(task) != "rescue_overlap":
        return resolved

    provider = task["provider"]
    entries_by_id = {}
    buckets = defaultdict(list)
    for transcript_id, features in cds_features_by_transcript.items():
        entry = build_transcript_grouping_entry(provider, transcript_id, features)
        entries_by_id[transcript_id] = entry
        buckets[(entry["seqid"], entry["strand"])].append(entry)

    adjacency = defaultdict(set)
    for bucket_entries in buckets.values():
        ordered = sorted(
            bucket_entries,
            key=lambda item: (item["span_start"], item["span_end"], item["transcript_id"]),
        )
        for index, left in enumerate(ordered):
            for right in ordered[index + 1 :]:
                if right["span_start"] > left["span_end"]:
                    break
                if left["gene_token"] == right["gene_token"]:
                    continue
                if should_rescue_overlapping_transcripts(left, right):
                    adjacency[left["transcript_id"]].add(right["transcript_id"])
                    adjacency[right["transcript_id"]].add(left["transcript_id"])

    visited = set()
    for transcript_id in sorted(entries_by_id.keys()):
        if transcript_id in visited:
            continue
        stack = [transcript_id]
        component_ids = []
        while len(stack) > 0:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)
            component_ids.append(current)
            for neighbor in sorted(adjacency.get(current, ())):
                if neighbor not in visited:
                    stack.append(neighbor)
        if len(component_ids) <= 1:
            continue
        component_entries = [entries_by_id[current_id] for current_id in sorted(component_ids)]
        cluster_gene_token = choose_rescue_cluster_gene_token(provider, component_entries)
        if cluster_gene_token == "":
            continue
        if not any(str(entry.get("gene_token", "") or "").strip() != cluster_gene_token for entry in component_entries):
            continue
        for entry in component_entries:
            resolved[entry["transcript_id"]] = cluster_gene_token
    return resolved

"""Annotation implementation: grouping."""

from collections import defaultdict

try:
    from Bio import SeqIO
except Exception:  # pragma: no cover - runtime without biopython
    SeqIO = None

from format_species_constants import (
    RESCUE_SAME_TERMINAL_MIN_LONGER_OVERLAP,
    RESCUE_SAME_TERMINAL_MIN_SHORTER_OVERLAP,
    RESCUE_SHARED_JUNCTION_MIN_LONGER_OVERLAP,
    RESCUE_SHARED_JUNCTION_MIN_SHORTER_OVERLAP,
)

from .common import (
    collapse_transcript_suffix,
    compute_interval_overlap_bases,
    gene_grouping_mode_for_task,
    merge_coordinate_intervals,
    transcript_feature_gene_token,
)


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

"""Conservative, inspectable duplication-origin evidence, not posterior inference."""

import csv
import math
from collections import defaultdict
from pathlib import Path


def read_table(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        fields = reader.fieldnames or []
        if not fields or len(set(fields)) != len(fields) or any(not field for field in fields):
            raise ValueError(f"Invalid TSV header: {path}")
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"Invalid TSV row width: {path}")
    return rows


def number(value):
    try:
        result = float(value)
    except (ValueError, TypeError):
        return None
    return result if math.isfinite(result) and result >= 0 else None


def branch_for_ks(focal, ks, tree, bounds):
    """Require the whole divergence interval to lie between corrected boundaries."""
    if not isinstance(ks, (tuple, list)):
        ks = (ks, ks)
    ks = tuple(number(value) for value in ks)
    if len(ks) != 2 or any(value is None for value in ks) or ks[0] > ks[1]:
        return None, "missing_or_invalid_ks"
    leaf = next((node for node in tree.leaves() if node.name == focal), None)
    if leaf is None:
        return None, "unknown_species"
    from nwkit.clade_index import CladeIndex

    index = CladeIndex(tree)
    child = leaf
    lower = 0.0
    while child.up is not None:
        parent = child.up
        row = bounds.get((focal, index.clade_id_for_node(parent)))
        if (row is None or row.get("status") != "ok"
                or row.get("monotone_from_younger_node") not in {"yes", "not_comparable"}):
            return None, "unresolved_divergence_boundary"
        # Use the primary boundary interval, not optional bootstrap diagnostics.
        lo, hi = number(row.get("ci_lower")), number(row.get("ci_upper"))
        if row.get("interval_status") != "ok" or lo is None or hi is None or hi < lo or lo < lower:
            return None, "unresolved_boundary_uncertainty"
        if ks[0] > lower and ks[1] < lo:
            return index.clade_id_for_node(child), "interval_supported"
        if ks[0] <= hi:
            return None, "boundary_overlap"
        lower = hi
        child = parent
    return None, "older_than_root_unresolved"


def valid_position(row):
    values = [number(row.get(field)) for field in ("start", "end", "rank")]
    return (all(value is not None and value == int(value) for value in values)
            and values[0] < values[1] and values[2] >= 1
            and all(row.get(field) for field in ("species", "locus_id", "seqid")))


def positional_feature(left, right, positions, proximal_distance=10):
    a, b = positions.get(left), positions.get(right)
    if a is None or b is None:
        return "unmapped", None
    if a["species"] != b["species"]:
        return "between_species", None
    if a["locus_id"] and a["locus_id"] == b["locus_id"]:
        return "same_locus", 0
    if not valid_position(a) or not valid_position(b):
        return "missing_or_invalid_coordinates", None
    if a["seqid"] != b["seqid"]:
        return "different_chromosomes", None
    distance = abs(int(a["rank"]) - int(b["rank"]))
    if max(number(a["start"]), number(b["start"])) < min(number(a["end"]), number(b["end"])):
        return "overlapping_loci", distance
    if distance == 0:
        return "ambiguous_locus_rank", distance
    if distance == 1:
        return "tandem", distance
    if 1 < distance <= proximal_distance:
        return "proximal", distance
    return "distant_same_chromosome", distance


def valid_anchor(row, positions):
    left, right = row.get("gene_a"), row.get("gene_b")
    ks = number(row.get("ks"))
    return (row.get("placement_status") == "interval_supported" and row.get("ks_status") == "ok"
            and ks is not None and ks > 0 and bool(row.get("species_event_id")) and bool(row.get("block_id"))
            and left in positions and right in positions
            and positions[left]["species"] == positions[right]["species"] == row.get("species")
            and positional_feature(left, right, positions)[0] in {
                "different_chromosomes", "distant_same_chromosome", "proximal"})


def combine_node(mapped_branch, pairs, positions, anchors, events, proximal_distance=10, terminal_cherry=True):
    """Only direct positive evidence supports an origin; missing evidence is not SSD."""
    features = []
    anchored = []
    for left, right in pairs:
        feature, distance = positional_feature(left, right, positions, proximal_distance)
        features.append((left, right, feature, distance))
        anchored.extend(anchors.get(tuple(sorted((left, right))), []))
    if any(feature[2] == "same_locus" for feature in features):
        return "unresolved", "same_locus_annotation_ambiguity", features, anchored
    if any(feature[2] in {"overlapping_loci", "ambiguous_locus_rank"} for feature in features):
        return "unresolved", "overlapping_or_ambiguous_loci", features, anchored
    if any(feature == "tandem" and anchors.get(tuple(sorted((left, right))))
           for left, right, feature, _ in features):
        return "unresolved", "tandem_and_collinearity_conflict", features, anchored
    pair_ids = {tuple(sorted(pair)) for pair in pairs}
    supported = [row for row in anchored
                 if row.get("species_event_id") == mapped_branch
                 and tuple(sorted((row.get("gene_a", ""), row.get("gene_b", "")))) in pair_ids
                 and valid_anchor(row, positions)
                 and events.get(mapped_branch, {}).get("event_support") == "WGD-supported"]
    if supported:
        return "WGD-supported", "calibrated_count_and_branch_matched_anchor_ks", features, anchored
    # A single terminal cherry has direct adjacency evidence. Larger ancestral
    # nodes and distant copies are not assigned SSD from absence of anchors.
    if terminal_cherry and len(pairs) == 1 and features[0][2] == "tandem":
        if anchored:
            return "unresolved", "tandem_and_collinearity_conflict", features, anchored
        return "SSD-supported", "terminal_tandem_adjacency", features, anchored
    return "unresolved", "insufficient_or_conflicting_origin_evidence", features, anchored


def summarize_events(candidates, anchors, species_summaries, min_coverage=0.2, min_blocks=3, positions=None):
    positions = positions or {}
    evidence = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    genes = defaultdict(lambda: defaultdict(set))
    for row in anchors:
        if valid_anchor(row, positions):
            loci = {positions[row[gene]]["locus_id"] for gene in ("gene_a", "gene_b")}
            evidence[row["species_event_id"]][row["species"]][row["block_id"]].update(loci)
            genes[row["species_event_id"]][row["species"]].update(loci)
    result = []
    for candidate in candidates:
        branch = candidate["species_event_id"]
        taxa = next(csv.reader([candidate["descendant_taxa"]]))
        species = []
        for name, blocks in sorted(evidence[branch].items()):
            summary = species_summaries.get(name, {})
            denominator = number(summary.get("num_annotated_loci"))
            # Shared locus arms join redundant/overlapping blocks into one
            # evidence component; fragmented block IDs are not replication.
            parents = {block: block for block in blocks}

            def root(block, parents=parents):
                while parents[block] != block:
                    parents[block] = parents[parents[block]]
                    block = parents[block]
                return block

            owner = {}
            for block, loci in blocks.items():
                for locus in loci:
                    if locus in owner:
                        parents[root(block)] = root(owner[locus])
                    owner[locus] = block
            independent = len({root(block) for block in blocks})
            if (name in taxa and denominator is not None and denominator > 0 and denominator == int(denominator)
                    and independent >= min_blocks
                    and min_coverage <= len(genes[branch][name]) / denominator <= 1):
                species.append(name)
        required = min(2, len(taxa))
        supported = (bool(taxa) and len(set(taxa)) == len(taxa)
                     and candidate["count_support"] == "count_supported_conditional" and len(species) >= required
                     and all(str(candidate.get(field, "")).lower() in {"false", "0"} for field in (
                         "nuisance_bound_reached", "background_nuisance_bound_reached",
                         "branch_burst_nuisance_bound_reached")))
        result.append({**candidate, "event_support": "WGD-supported" if supported else "unresolved",
                       "num_branch_matched_synteny_species": len(species),
                       "branch_matched_synteny_species": ",".join(species),
                       "minimum_branch_matched_blocks": min_blocks, "minimum_synteny_coverage": min_coverage,
                       "synteny_coverage_definition": "unique_branch_matched_anchor_loci_over_annotated_loci",
                       "synteny_block_definition": "locus_disjoint_block_components",
                       "support_meaning": "experimental_evidence_rule_not_posterior"})
    return result

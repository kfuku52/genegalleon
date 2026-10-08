"""Bounded genome-only exploration and conservative annotation placement.

Homology nominates a protein search; it never assigns orthology or a WGD copy.
This module consumes already prepared, receipt-verified pairwise comparisons.
"""
import hashlib
import json
import math
from collections import Counter, defaultdict, deque
from pathlib import Path

try:
    from fasta_sequence_store import fasta_records
    from rescue_coding_paths import coding_shape, compatible_paths, overlap_components
except ImportError:
    from .fasta_sequence_store import fasta_records
    from .rescue_coding_paths import coding_shape, compatible_paths, overlap_components


def _proteins(path):
    result = {}
    for identifier, _, sequence in fasta_records(path):
        if identifier in result:
            raise ValueError("Duplicate prepared protein: " + identifier)
        if not sequence or set(sequence.upper()) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO"):
            raise ValueError("Invalid prepared donor protein: " + identifier)
        result[identifier] = sequence
    return result


def _representation(path, target_lengths, donor_lengths, reverse, params):
    """Keep per-target matches; never union disjoint paralogs into a full hit."""
    matches = defaultdict(dict)
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 12:
                raise ValueError(f"Invalid pairwise protein alignment at {path}:{number}")
            target, donor = (fields[1], fields[0]) if reverse else fields[:2]
            if target not in target_lengths or donor not in donor_lengths:
                raise ValueError("Protein comparison references an unknown prepared protein")
            identity, bitscore = float(fields[2]) / 100, float(fields[11])
            coordinates = fields[6:8] if reverse else fields[8:10]
            start, end = sorted(map(int, coordinates))
            if (not 0 < start <= end <= donor_lengths[donor] or not math.isfinite(identity)
                    or not 0 <= identity <= 1 or not math.isfinite(bitscore) or bitscore <= 0):
                raise ValueError("Invalid pairwise protein alignment coordinates or score")
            # Alignment length includes gaps; donor coordinate coverage avoids
            # counting a long insertion in the target as donor representation.
            coverage = (end - start + 1) / donor_lengths[donor]
            record = {"target": target, "coverage": coverage, "identity": identity, "bitscore": bitscore}
            old = matches[donor].get(target)
            if old is None or (coverage, bitscore) > (old["coverage"], old["bitscore"]):
                matches[donor][target] = record
    result = {}
    unique_target_owners = defaultdict(set)
    for donor, values in matches.items():
        full = [m for m in values.values() if m["coverage"] >= params["minimum_coverage"]
                and m["identity"] >= params["minimum_identity"]]
        best = max((m["bitscore"] for m in full), default=0)
        contenders = [m for m in full if m["bitscore"] >= best * params.get("cscore", .7)]
        if len(contenders) == 1:
            unique_target_owners[contenders[0]["target"]].add(donor)
        strongest = max(values.values(), key=lambda m: (m["identity"] >= params["minimum_identity"], m["coverage"], m["bitscore"]))
        result[donor] = {"reason": "represented" if len(contenders) == 1 else
                         "copy_ambiguous" if contenders else "partial_target_match",
                         "target_matches": len(values), "full_target_matches": len(contenders),
                         "best_target_coverage": strongest["coverage"], "best_target_identity": strongest["identity"],
                         "best_target_bitscore": strongest["bitscore"],
                         "full_target_ids": sorted(m["target"] for m in contenders)}
    for state in result.values():
        if state["reason"] == "represented" and len(unique_target_owners[state["full_target_ids"][0]]) > 1:
            state["reason"] = "copy_ambiguous"
            state["ambiguity"] = "multiple_donor_genes_share_one_full_target"
    return result


def _balanced_queue(rows, nearest, minimum_identity):
    """Weighted donor/bucket round robin; input labels never define priority."""
    queues = defaultdict(lambda: defaultdict(list))
    for row in rows:
        state = row["nomination"]
        reason = state["reason"]
        bucket = ("strong_partial" if state.get("best_target_identity", 0) >= minimum_identity
                  else "weak_partial") if reason == "partial_target_match" else reason
        row["screen_priority"] = bucket
        queues[row["donor"]][bucket].append(row)
    tracks = {}
    for donor, buckets in queues.items():
        tracks[donor] = {key: deque(sorted(values, key=lambda r: (-r["nomination"].get("best_target_identity", 0),
                                                               -r["nomination"].get("best_target_coverage", 0), r["query"])))
                         for key, values in buckets.items()}
    groups = {"nearest": deque(sorted(set(tracks) & set(nearest))),
              "common": deque(sorted(set(tracks) - set(nearest)))}
    buckets = ("strong_partial", "strong_partial", "no_target_match", "copy_ambiguous")
    cursors = defaultdict(int)
    output = []
    schedule = ("nearest", "nearest", "common")
    while any(groups.values()):
        for group in schedule:
            donors = groups[group] or groups["common" if group == "nearest" else "nearest"]
            if not donors:
                continue
            donor = donors.popleft()
            track = tracks[donor]
            available = [key for key in buckets if track.get(key)] or [key for key in track if track[key]]
            if not available:
                continue
            chosen = None
            for _ in buckets:
                key = buckets[cursors[donor] % len(buckets)]
                cursors[donor] += 1
                if track.get(key):
                    chosen = key
                    break
            chosen = chosen or available[0]
            output.append(track[chosen].popleft())
            if any(track.values()):
                donors.append(donor)
    return output


def nominate_genome_only_candidates(root, plan, name, regions, params):
    """Return missing/partial/ambiguous donor queries not covered by two anchors.

    Only frozen nearest/common donor comparisons are read. Existing interval
    candidates are untouched, including their independent WGD-copy windows.
    Rows have no invented genomic coordinates and must bypass interval search.
    """
    root = Path(root)
    limit = params.get("max_genome_queries", 20_000)
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 0:
        raise ValueError("max_genome_queries must be a nonnegative integer")
    if not params.get("genome_fallback", True):
        return [], {"policy": "disabled_genome_fallback", "nominated": 0, "max_genome_queries": limit,
                    "deferred_by_limit": 0, "deferred_queries": [], "counts": {}, "by_donor": {}}
    target_lengths = {key: len(seq) for key, seq in _proteins(root / "prepared" / name / "genes.pep").items()}
    nominated = {(r["donor"], r["query"]) for r in regions}
    donors = set(plan["donors"][name]) - {name}
    rows, counts, by_donor = [], Counter(), {}
    seen = set()
    for job in sorted(plan["synteny_jobs"], key=lambda j: j["id"]):
        if job["a"] == job["b"] or name not in {job["a"], job["b"]}:
            continue
        reverse = job["b"] == name
        donor = job["a"] if reverse else job["b"]
        if donor not in donors:
            continue
        if donor in seen:
            raise ValueError("Duplicate selected donor comparison: " + donor)
        seen.add(donor)
        proteins = _proteins(root / "prepared" / donor / "genes.pep")
        representation = _representation(root / "synteny" / job["id"] / "target.query.last",
                                         target_lengths, {key: len(seq) for key, seq in proteins.items()}, reverse, params)
        donor_counts = Counter()
        for query in sorted(proteins):
            state = representation.get(query, {"reason": "no_target_match", "target_matches": 0,
                                               "full_target_matches": 0, "best_target_coverage": 0,
                                               "best_target_identity": 0, "best_target_bitscore": 0,
                                               "full_target_ids": []})
            reason = "already_nominated" if (donor, query) in nominated else state["reason"]
            counts[reason] += 1
            donor_counts[reason] += 1
            if reason in {"represented", "already_nominated"}:
                continue
            row = {"target": name, "donor": donor, "query": query, "comparison": job["id"],
                   "genome_only": True, "placement": "unanchored", "nomination": state,
                   "orthology": "unassigned", "expected_copy": "unassigned"}
            row["id"] = "genome_query_" + hashlib.sha256(json.dumps(row, sort_keys=True).encode()).hexdigest()[:20]
            rows.append(row)
        by_donor[donor] = dict(donor_counts)
    if seen != donors:
        raise ValueError("Missing selected donor protein comparison: " + ", ".join(sorted(donors - seen)))
    nearest = set(plan.get("nearest_references", {}).get(name, ()))
    rows = _balanced_queue(rows, nearest, params["minimum_identity"])
    deferred = rows[limit:]
    return rows[:limit], {"policy": "selected_donor_protein_representation", "counts": dict(counts),
                          "by_donor": by_donor, "max_genome_queries": limit, "nominated": min(len(rows), limit),
                          "selection_policy": "nearest_common_weight_2_to_1_donor_and_evidence_bucket_round_robin",
                          "selected_by_donor": dict(Counter(r["donor"] for r in rows[:limit])),
                          "selected_by_priority": dict(Counter(r["screen_priority"] for r in rows[:limit])),
                          "deferred_by_limit": len(deferred),
                          "deferred_queries": [{"donor": r["donor"], "query": r["query"],
                                                "reason": r["nomination"]["reason"], "priority": r["screen_priority"]} for r in deferred]}


def reassess_unanchored_models(models, minimum_species=2, *, ownership=None):
    """Admit only intact, unique paths corroborated by independent species.

    Compatible coding paths may corroborate annotation placement. A donor
    protein mapping to multiple unannotated loci cannot corroborate uniqueness.
    Existing annotated placements retain the original all-loci assessment for
    revision nominations. Missing-gene annotation uses unannotated loci only;
    neither assessment assigns orthology or an expected copy.
    """
    if isinstance(minimum_species, bool) or not isinstance(minimum_species, int) or minimum_species < 2:
        raise ValueError("Unanchored placement needs at least two donor species")
    placement_problems = {"outside_expected_synteny_interval", "unanchored_genome_search"}
    eligible = [m for m in models if m.get("sequence") and m.get("cds")
                and not (set(m.get("problems", ())) - placement_problems)
                and all(not isinstance(m.get(key), bool) and isinstance(m.get(key), (int, float))
                        and math.isfinite(m[key]) and 0 <= m[key] <= 1 for key in ("coverage", "identity"))]
    shapes = {}
    for model in eligible:
        key = coding_shape(model)
        if key not in shapes:
            shapes[key] = {"model": model, "records": []}
        shapes[key]["records"].append(model)
    unique = [row["model"] for row in shapes.values()]
    for row in unique:
        row["start"] = min(e[0] for e in row["cds"])
        row["end"] = max(e[1] for e in row["cds"])
    def assessments(paths):
        components = overlap_components(paths)
        query_loci = defaultdict(set)
        for locus, indexes in enumerate(components):
            for index in indexes:
                for row in shapes[coding_shape(paths[index])]["records"]:
                    evidence = row.get("evidence", {})
                    if evidence.get("donor") and evidence.get("query"):
                        query_loci[evidence["donor"], evidence["query"]].add(locus)
        result = {}
        for indexes in components:
            component = [paths[i] for i in indexes]
            # Every pair must be frame-compatible; a bridge cannot fuse loci.
            # Only placement is reassessed, using verified primary alignment
            # metrics. Per-path support and scores remain unchanged.
            clean = [{**path, "problems": [], "support": []} for path in component]
            coherent = all(compatible_paths(a, b) for i, a in enumerate(clean) for b in clean[i + 1:])
            records = (m for path in component for m in shapes[coding_shape(path)]["records"])
            support = set()
            for model in records:
                evidence = model.get("evidence", {})
                key = evidence.get("donor"), evidence.get("query")
                if (key[0] and key[1] and key[0] != evidence.get("target")
                        and len(query_loci[key]) == 1):
                    support.add(key[0])
            reason = "supported_unanchored_annotation" if coherent and len(support) >= minimum_species else \
                     "incompatible_unanchored_paths" if not coherent else "insufficient_unique_donor_support"
            for path in component:
                result[coding_shape(path)] = {"reason": reason, "independent_donor_species": sorted(support)}
        return result
    original = assessments(unique)
    owners = ({coding_shape(path): ownership.overlapping(path) for path in unique} if ownership is not None else {})
    unannotated = [path for path in unique if not owners.get(coding_shape(path))]
    missing = assessments(unannotated) if ownership is not None else original
    counts = Counter()
    newly_supported = 0
    for path in unique:
        shape = coding_shape(path)
        owned = bool(owners.get(shape))
        decision = original[shape] if owned else missing[shape]
        reason = decision["reason"]
        for model in shapes[shape]["records"]:
            if not (set(model.get("problems", ())) & placement_problems):
                continue
            counts[reason] += 1
            model["placement_evidence"] = {"classification": "unanchored", "reason": reason,
                                           "independent_donor_species": decision["independent_donor_species"],
                                           "minimum_species": minimum_species, "orthology": "unassigned",
                                           "expected_copy": "unassigned"}
            if ownership is not None:
                model["placement_evidence"].update({
                    "uniqueness_scope": "all_coding_loci_for_existing_model_revision" if owned else "unannotated_coding_loci",
                    "original_all_loci_assessment": original[shape],
                    "counting_unit": "independent_donor_species"})
                if owned:
                    model["placement_evidence"]["original_annotation_owner_ids"] = sorted({
                        row["gene_id"] for row in owners[shape] if row.get("gene_id")})
                elif reason == "supported_unanchored_annotation" and original[shape]["reason"] != reason:
                    model["placement_evidence"]["uniqueness_reason"] = "annotated_alternatives_excluded_from_missing_gene_uniqueness"
                    newly_supported += 1
            if reason == "supported_unanchored_annotation":
                model["problems"] = [p for p in model["problems"] if p not in placement_problems]
    result = {"minimum_species": minimum_species, "counts": dict(counts), "intact_unique_paths": len(unique)}
    if ownership is not None:
        result["annotation_aware_uniqueness"] = {
            "annotated_intact_coding_paths": len(unique) - len(unannotated),
            "unannotated_intact_coding_paths": len(unannotated),
            "newly_supported_missing_annotation_records": newly_supported,
            "existing_model_revision_scoring": "original_all_coding_loci",
            "orthology": "unassigned", "expected_copy": "unassigned"}
    return result

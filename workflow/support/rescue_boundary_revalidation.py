"""Recheck nearby coding boundaries using source DNA and complete donor alignments.

This is a second prediction step, not a relaxation of the strict CDS validator.
Original predictor records and failed or ambiguous alternatives remain auditable.
"""
from __future__ import annotations

import copy
import hashlib
import itertools
import math

from Bio.Align import PairwiseAligner
from Bio.Data import CodonTable
from Bio.Seq import Seq

try:
    from rescue_homology_policy import homology_supported
except ImportError:
    from .rescue_homology_policy import homology_supported

_CANONICAL = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
_RECHECKABLE = {"missing_start", "missing_stop", "noncanonical_splice"}
_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
_SPLICE_RADIUS = 30
_PATH_BUDGET = 128


def _sha(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _dna(model, genome):
    pieces = [genome.fetch(model["seqid"], a, b).upper() for a, b, _ in model["cds"]]
    if model["strand"] == "-":
        pieces = [str(Seq(piece).reverse_complement()) for piece in pieces]
    return "".join(pieces)


def _phases(blocks):
    result, cumulative = [], 0
    for a, b, _ in blocks:
        result.append([a, b, (3 - cumulative % 3) % 3])
        cumulative += b - a
    return result


def _junction(model, previous, following, genome):
    a, b = (previous[1], following[0]) if model["strand"] == "+" else (following[1], previous[0])
    if a < 0 or b > genome.get_reference_length(model["seqid"]) or b - a < 4:
        return None
    donor = genome.fetch(model["seqid"], a, a + 2).upper()
    acceptor = genome.fetch(model["seqid"], b - 2, b).upper()
    if model["strand"] == "-":
        donor, acceptor = str(Seq(acceptor).reverse_complement()), str(Seq(donor).reverse_complement())
    return donor, acceptor


def _splice_choices(model, genome, maximum_intron):
    choices = []
    for i, (previous, following) in enumerate(zip(model["cds"], model["cds"][1:], strict=False)):
        if _junction(model, previous, following, genome) in _CANONICAL:
            choices.append([(0, 0)])
            continue
        found = []
        left_field, right_field = (1, 0) if model["strand"] == "+" else (0, 1)
        for x in range(-_SPLICE_RADIUS, _SPLICE_RADIUS + 1):
            for y in range(-_SPLICE_RADIUS, _SPLICE_RADIUS + 1):
                p, c = previous.copy(), following.copy()
                p[left_field] += x
                c[right_field] += y
                if not (0 <= p[0] < p[1] <= genome.get_reference_length(model["seqid"])
                        and 0 <= c[0] < c[1] <= genome.get_reference_length(model["seqid"])):
                    continue
                intron = c[0] - p[1] if model["strand"] == "+" else p[0] - c[1]
                if 4 <= intron <= maximum_intron and _junction(model, p, c, genome) in _CANONICAL:
                    found.append((x, y))
        if not found:
            return [], "no_nearby_canonical_splice"
        choices.append(found)
        if math.prod(map(len, choices)) > _PATH_BUDGET:
            return [], "splice_candidate_budget_exceeded"
    return choices, None


def _starts(model, genome, maximum, allowed, needed):
    if not needed:
        return [model["cds"]]
    first = model["cds"][0]
    length = genome.get_reference_length(model["seqid"])
    result = []
    # A prediction boundary may be off by one base. Examine all positions,
    # then require a complete, frame-consistent CDS and a full donor alignment.
    for offset in range(-min(_SPLICE_RADIUS, first[1] - first[0] - 1), maximum + 1):
        blocks = copy.deepcopy(model["cds"])
        field = 0 if model["strand"] == "+" else 1
        blocks[0][field] += -offset if model["strand"] == "+" else offset
        a, b = blocks[0][:2]
        if not 0 <= a < b <= length:
            continue
        pieces, needed_bases = [], 3
        for exon_start, exon_end, _ in blocks:
            extent = min(needed_bases, exon_end - exon_start)
            piece = (genome.fetch(model["seqid"], exon_start, exon_start + extent) if model["strand"] == "+"
                     else str(Seq(genome.fetch(model["seqid"], exon_end - extent, exon_end)).reverse_complement()))
            pieces.append(piece.upper())
            needed_bases -= extent
            if not needed_bases:
                break
        start = "".join(pieces)
        if start in allowed:
            result.append(blocks)
    return result


def _stop(model, genome, maximum, stop_codons):
    sequence = _dna(model, genome)
    if len(sequence) % 3 == 0 and sequence[-3:] in stop_codons:
        return model
    blocks = copy.deepcopy(model["cds"])
    last = blocks[-1]
    length = genome.get_reference_length(model["seqid"])
    room = length - last[1] if model["strand"] == "+" else last[0]
    extent = min(maximum, room)
    downstream = (genome.fetch(model["seqid"], last[1], last[1] + extent) if model["strand"] == "+"
                  else str(Seq(genome.fetch(model["seqid"], last[0] - extent, last[0])).reverse_complement())).upper()
    dna = sequence + downstream
    # The first genomic in-frame stop is the only contiguous terminal choice.
    for end in range(len(sequence) + 1, len(dna) + 1):
        if end % 3 == 0 and dna[end - 3:end] in stop_codons:
            if model["strand"] == "+":
                last[1] += end - len(sequence)
            else:
                last[0] -= end - len(sequence)
            return {**model, "cds": _phases(blocks)}
    return None


def _alignment(protein, donor, params, overhang):
    if len(protein) * len(donor) > params.get("terminal_max_alignment_cells", 25_000_000):
        return None, "alignment_cell_budget_exceeded"
    aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                              open_gap_score=-5, extend_gap_score=-.5)
    alignments = aligner.align(protein, donor)
    limit = params.get("terminal_max_optimal_alignments", 8)
    if type(limit) is not int or not 0 < limit <= 32:
        raise ValueError("Invalid donor alignment traceback budget")
    evidence = []
    for index in range(limit):
        try:
            alignment = alignments[index]
        except IndexError:
            break
        pairs = {}
        for (a, b), (c, d) in zip(*alignment.aligned, strict=True):
            pairs.update((int(a + k), int(c + k)) for k in range(int(b - a)))
        aligned = len(pairs)
        identity = sum(protein[a] == donor[b] for a, b in pairs.items()) / max(1, aligned)
        coverage, recipient = aligned / len(donor), aligned / len(protein)
        terminal = (pairs.get(0) == 0 and any(b == len(donor) - 1 and a >= len(protein) - 1 - overhang
                                           for a, b in pairs.items()))
        evidence.append({"coverage": coverage, "identity": identity, "recipient_coverage": recipient,
                         "score": identity * min(coverage, recipient), "donor_termini_supported": terminal,
                         "query_start": min(pairs.values(), default=0),
                         "query_end": max(pairs.values(), default=-1) + 1})
    truncated = False
    try:
        alignments[limit]
    except IndexError:
        pass
    else:
        truncated = True
    valid = [e for e in evidence if e["donor_termini_supported"]
             and homology_supported(min(e["coverage"], e["recipient_coverage"]), e["identity"], params)]
    if not valid:
        return {"optimal_alignments": evidence}, ("donor_alignment_budget_exceeded" if truncated
                                                  else "no_complete_donor_alignment")
    best = max(valid, key=lambda e: e["score"])
    best["optimal_alignments_reviewed"] = len(evidence)
    best["optimal_alignments_truncated"] = truncated
    # Multiple tracebacks for one identical genomic CDS do not create multiple
    # genomic paths. Distinct CDS coordinates are compared separately below.
    return best, "donor_alignment_budget_exceeded" if truncated else None


def revalidate_boundaries(model, genome, genetic_code, params, donor_protein, validate_model, *, copy_model=None):
    problems = set(model.get("problems", ()))
    if not problems & _RECHECKABLE or problems - _RECHECKABLE:
        return model
    if model.get("frameshift") or not homology_supported(model["coverage"], model["identity"], params):
        return model
    copier = copy.deepcopy if copy_model is None else copy_model
    result = copier(model)
    maximum = params.get("terminal_max_extension", 300)
    overhang = params.get("terminal_max_unaligned_c_overhang", 2)
    margin = params.get("terminal_minimum_margin", .02)
    if type(maximum) is not int or maximum < 1 or type(overhang) is not int or not 0 <= overhang <= 2:
        raise ValueError("Invalid coding-boundary search limits")
    if type(margin) not in (int, float) or not math.isfinite(margin) or not 0 <= margin <= 1:
        raise ValueError("Invalid coding-boundary score margin")
    donor = str(donor_protein or "").upper().removesuffix("*")
    audit = {"schema": 1, "status": "proposal", "splice_radius_nt": _SPLICE_RADIUS,
             "maximum_candidate_paths": _PATH_BUDGET, "original_cds": copy.deepcopy(model["cds"]),
             "original_sequence_sha256": _sha(model.get("sequence", "")), "original_problems": sorted(problems),
             "original_alignment": {k: model.get(k) for k in ("coverage", "identity", "query_start", "query_end", "query_length", "paf")},
             "attempts": [], "reasons": []}
    result["boundary_revalidation"] = audit

    def reject(reason):
        audit.update(status="rejected", reasons=[reason])
        return result

    if not donor or set(donor) - _AMINO_ACIDS:
        return reject("missing_or_ambiguous_donor_protein")
    if _dna(model, genome) != model.get("sequence"):
        return reject("source_sequence_mismatch")
    baseline = validate_model(copier(model), genome, genetic_code, params)
    if baseline["cds"] != model["cds"] or baseline["sequence"] != model["sequence"] or set(baseline["problems"]) != problems:
        return reject("input_not_already_strictly_checked")
    choices, error = _splice_choices(model, genome, params.get("max_intron", math.inf))
    if error:
        return reject(error)
    codons = CodonTable.unambiguous_dna_by_id[genetic_code]
    candidates, examined = {}, set()
    audit["genomic_rejection_counts"], audit["genomic_rejection_samples"] = {}, []
    left_field, right_field = (1, 0) if model["strand"] == "+" else (0, 1)
    for shifts in itertools.product(*choices):
        blocks = copy.deepcopy(model["cds"])
        for i, (x, y) in enumerate(shifts):
            blocks[i][left_field] += x
            blocks[i + 1][right_field] += y
        # Adjacent junction proposals can otherwise invert a short middle
        # exon even though each junction was valid in isolation.
        if any(not 0 <= a < b <= genome.get_reference_length(model["seqid"]) for a, b, _ in blocks):
            continue
        proposed = {**model, "cds": _phases(blocks)}
        for start_blocks in _starts(proposed, genome, maximum, codons.start_codons, "missing_start" in problems):
            proposed = {**model, "cds": _phases(start_blocks)}
            c_bound = min(maximum, 3 * (overhang + 1)) if model.get("query_end") == len(donor) else maximum
            stopped = _stop(proposed, genome, c_bound, codons.stop_codons)
            if stopped is None:
                continue
            key = tuple(tuple(b) for b in stopped["cds"])
            if key in examined:
                continue
            examined.add(key)
            audit["examined_genomic_paths"] = len(examined)
            if len(examined) > _PATH_BUDGET * 32:
                return reject("examined_genomic_path_budget_exceeded")
            checked = validate_model(copier(stopped), genome, genetic_code, params)
            if checked["problems"]:
                for problem in checked["problems"]:
                    counts = audit["genomic_rejection_counts"]
                    counts[problem] = counts.get(problem, 0) + 1
                if len(audit["genomic_rejection_samples"]) < 16:
                    audit["genomic_rejection_samples"].append({"cds": checked["cds"], "problems": checked["problems"]})
                continue
            candidates[key] = checked
            if len(candidates) > _PATH_BUDGET:
                return reject("strict_genomic_candidate_budget_exceeded")
    if not candidates:
        return reject("no_complete_genomic_boundary_path_within_bounds")
    valid, uncertain = [], []
    for candidate in candidates.values():
        checked = validate_model(copier(candidate), genome, genetic_code, params)
        attempt = {"cds": copy.deepcopy(checked["cds"]), "problems": checked["problems"]}
        audit["attempts"].append(attempt)
        if checked["problems"]:
            continue
        protein = str(Seq(checked["sequence"]).translate(table=genetic_code)).removesuffix("*")
        if checked["sequence"][:3] in codons.start_codons and donor.startswith("M"):
            protein = "M" + protein[1:]
        alignment, error = _alignment(protein, donor, params, overhang)
        attempt.update(alignment=alignment, sequence_sha256=_sha(checked["sequence"]))
        if error:
            attempt["reasons"] = [error]
            if error in {"alignment_cell_budget_exceeded", "donor_alignment_budget_exceeded"}:
                uncertain.append(attempt)
            if not alignment or not alignment.get("donor_termini_supported"):
                continue
        checked.update(coverage=alignment["coverage"], identity=alignment["identity"],
                       query_start=alignment["query_start"], query_end=alignment["query_end"],
                       query_length=len(donor), query_span_coverage=(alignment["query_end"] - alignment["query_start"]) / len(donor),
                       alignment_source="boundary_revalidation_global")
        valid.append((alignment["score"], checked, attempt))
    audit["supports_intact_homolog"] = bool(valid)
    audit["intact_donor_supported_paths"] = [
        {"cds": checked["cds"], "sequence_sha256": attempt["sequence_sha256"],
         "coverage": checked["coverage"], "identity": checked["identity"],
         "recipient_coverage": attempt["alignment"]["recipient_coverage"]}
        for _, checked, attempt in valid]
    audit["orthology"], audit["expected_copy"] = "unassigned", "unassigned"
    if uncertain:
        audit.update(status="ambiguous", reasons=["unresolved_donor_alignment_budget"])
        return result
    if not valid:
        return reject("no_strict_donor_supported_boundary_path")
    valid.sort(key=lambda row: row[0], reverse=True)
    if len(valid) > 1 and valid[0][0] - valid[1][0] <= margin:
        audit.update(status="ambiguous", reasons=["alternative_genomic_boundary_paths_within_score_margin"])
        return result
    _, selected, attempt = valid[0]
    selected = validate_model(selected, genome, genetic_code, params)
    if selected["problems"]:
        return reject("selected_boundary_path_failed_final_validation")
    audit.update(status="completed", selected_cds=attempt["cds"],
                 selected_sequence_sha256=attempt["sequence_sha256"],
                 score_margin=valid[0][0] - valid[1][0] if len(valid) > 1 else None)
    selected["boundary_revalidation"] = audit
    selected["partial_evidence"] = {"partial": False, "completed_from_partial": bool(problems & {"missing_start", "missing_stop"}),
                                    "original_missing_termini": sorted(problems & {"missing_start", "missing_stop"})}
    return selected

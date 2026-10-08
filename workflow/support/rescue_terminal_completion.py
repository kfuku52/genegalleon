"""Bounded, donor-supported extension of genuinely partial coding termini.

Only the outer ends of the first/last coding exon can move.  Existing CDS
residues, splice junctions, phases and the genetic code stay unchanged.  Failed
and ambiguous extensions remain evidence, never intact representative models.
"""

from __future__ import annotations

import copy
import hashlib
import math

from Bio.Align import PairwiseAligner
from Bio.Data import CodonTable
from Bio.Seq import Seq

_TERMINAL_PROBLEMS = {"missing_start", "missing_stop"}
_RECHECKABLE = _TERMINAL_PROBLEMS | {"low_coverage", "low_identity"}
_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def _sha(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _chain_sequence(model, genome):
    pieces = [genome.fetch(model["seqid"], start, end).upper() for start, end, _ in model["cds"]]
    if model["strand"] == "-":
        pieces = [str(Seq(piece).reverse_complement()) for piece in pieces]
    return "".join(pieces)


def _flank(model, genome, maximum, terminus):
    """Fetch each contiguous window once, in transcription order and frame."""
    exon = model["cds"][0 if terminus == "n" else -1]
    length = genome.get_reference_length(model["seqid"])
    left = (model["strand"] == "+") == (terminus == "n")
    boundary = exon[0] if left else exon[1]
    room = boundary if left else length - boundary
    extent = min(maximum, room) // 3 * 3
    start, end = (boundary - extent, boundary) if left else (boundary, boundary + extent)
    sequence = genome.fetch(model["seqid"], start, end).upper() if extent else ""
    if model["strand"] == "-":
        sequence = str(Seq(sequence).reverse_complement())
    return sequence, extent < maximum


def _alignment_evidence(protein, donor, n_added, c_added, original_begin, original_end, coordinates, c_overhang):
    """Project newly added residues onto the donor's missing terminal segments."""
    pairs, matches = {}, 0
    for index in range(coordinates.shape[1] - 1):
        a, b = int(coordinates[0, index]), int(coordinates[1, index])
        aa, bb = int(coordinates[0, index + 1]), int(coordinates[1, index + 1])
        if aa > a and bb > b:
            for offset in range(aa - a):
                pairs[a + offset] = b + offset
                matches += protein[a + offset] == donor[b + offset]
    aligned = len(pairs)
    matched_donor = set(pairs.values())
    identity = matches / max(1, aligned)
    coverage = aligned / len(donor)
    recipient_coverage = aligned / len(protein)
    reasons = []
    last_supported = len(protein) - c_added - 1 if c_overhang else len(protein) - 1
    if pairs.get(0) != 0 or pairs.get(last_supported) != len(donor) - 1:
        reasons.append("donor_termini_not_aligned")
    if any(index not in pairs or pairs[index] >= original_begin for index in range(n_added)):
        reasons.append("n_extension_not_supported_by_missing_donor_terminus")
    if c_overhang and any(index in pairs for index in range(len(protein) - c_added, len(protein))):
        reasons.append("short_c_overhang_did_not_follow_donor_terminal")
    elif not c_overhang and any(index not in pairs or pairs[index] < original_end
                               for index in range(len(protein) - c_added, len(protein))):
        reasons.append("c_extension_not_supported_by_missing_donor_terminus")
    if any(index not in matched_donor for index in range(original_begin, original_end)):
        reasons.append("internal_donor_gap")
    # Both previously missing terminal segments and the unchanged core must
    # contribute independent matching residues, not only aligned mismatches.
    for side, positions in (("n", range(n_added)), ("c", range(len(protein) - c_added, len(protein)))):
        if positions and not (side == "c" and c_overhang) and not any(
                index in pairs and protein[index] == donor[pairs[index]] for index in positions):
            reasons.append(side + "_extension_has_no_matching_residue")
    positions = sorted({0, last_supported, *range(n_added), *range(len(protein) - c_added, len(protein))})
    return {"identity": identity, "coverage": coverage, "recipient_coverage": recipient_coverage,
            "score": identity * min(coverage, recipient_coverage), "aligned_residues": aligned,
            "query_start": min(matched_donor, default=0), "query_end": max(matched_donor, default=-1) + 1,
            "c_terminal_support": "short_genomic_overhang_after_aligned_donor" if c_overhang else "donor_alignment",
            "unaligned_c_overhang_residues": c_added if c_overhang else 0,
            "terminal_projection": [[index, pairs.get(index)] for index in positions], "reasons": reasons}


def _protein_evidence(protein, donor, n_added, c_added, original_begin, original_end, aligner, c_overhang, limit):
    """Do not turn an arbitrary optimal traceback into unique terminal support."""
    alignments = aligner.align(protein, donor)
    examined = []
    for index in range(limit):
        try:
            coordinates = alignments[index].coordinates
        except IndexError:
            break
        examined.append(_alignment_evidence(protein, donor, n_added, c_added, original_begin, original_end,
                                             coordinates, c_overhang))
    best = min(examined, key=lambda row: (len(row["reasons"]), -row["score"]))
    signatures = {tuple(tuple(pair) for pair in row["terminal_projection"]) for row in examined}
    best["optimal_alignments_reviewed"] = len(examined)
    best["terminal_projection_count"] = len(signatures)
    if len(signatures) > 1:
        best["reasons"].append("ambiguous_donor_terminal_alignment")
    try:
        alignments[limit]
    except IndexError:
        pass
    else:
        best["reasons"].append("donor_alignment_budget_exceeded")
    return best


def complete_terminals(model, genome, genetic_code, params, donor_protein, validate_model):
    """Return a checked model or the unchanged partial model with an audit.

    ``validate_model(model, genome, genetic_code, params)`` is the original
    non-recursive strict validator.  Input CDS blocks are 0-based half-open in
    transcription order, including phases.  Alignment metrics updated after
    completion are explicitly distinguished from the original miniprot PAF.
    """
    problems = set(model.get("problems", []))
    needed = problems & _TERMINAL_PROBLEMS
    if not needed:
        # Most predictions never need this helper.  Nested provenance can be
        # large; it is read-only here, while the caller's problem list is private.
        unchanged = {**model, "problems": list(model.get("problems", []))}
        if problems & {"incomplete_frame", "no_cds"}:
            unchanged["partial_evidence"] = {"partial": True, "missing_termini": [],
                                              "partial_reasons": sorted(problems), "representative_eligible": False}
        return unchanged
    result = copy.deepcopy(model)
    maximum = params.get("terminal_max_extension", 300)
    candidate_limit = params.get("terminal_max_candidates", 32)
    cell_limit = params.get("terminal_max_alignment_cells", 25_000_000)
    alignment_limit = params.get("terminal_max_optimal_alignments", 8)
    max_overhang = params.get("terminal_max_unaligned_c_overhang", 2)
    margin = params.get("terminal_minimum_margin", .02)
    for value, name in ((maximum, "terminal_max_extension"), (candidate_limit, "terminal_max_candidates"),
                        (cell_limit, "terminal_max_alignment_cells")):
        if type(value) is not int or value <= 0:
            raise ValueError(name + " must be a positive integer")
    if type(alignment_limit) is not int or not 0 < alignment_limit <= 32:
        raise ValueError("terminal_max_optimal_alignments must be an integer between one and 32")
    if type(max_overhang) is not int or not 0 <= max_overhang <= 2:
        raise ValueError("terminal_max_unaligned_c_overhang must be an integer between zero and two")
    if not isinstance(margin, (int, float)) or isinstance(margin, bool) or not math.isfinite(margin) or not 0 <= margin <= 1:
        raise ValueError("terminal_minimum_margin must be between zero and one")
    audit = {"schema": 1, "status": "proposal", "maximum_extension_nt": maximum,
             "maximum_unaligned_c_overhang_residues": max_overhang,
             "maximum_optimal_alignments_reviewed": alignment_limit,
             "original_cds": copy.deepcopy(result["cds"]),
             "original_sequence_sha256": _sha(result.get("sequence", "")),
             "original_problems": sorted(problems),
             "original_alignment": {key: result.get(key) for key in
                                    ("coverage", "identity", "query_start", "query_end", "query_length", "paf")},
             "attempts": [], "search_limits": [], "reasons": []}
    result["terminal_completion"] = audit
    result["partial_evidence"] = {"partial": True, "missing_termini": sorted(needed),
                                  "partial_reasons": sorted(problems), "representative_eligible": False}

    def reject(reason):
        audit["status"] = "rejected"
        audit["reasons"].append(reason)
        return result

    if problems - _RECHECKABLE or result.get("frameshift", False):
        return reject("non_terminal_failure_preserved")
    if result["strand"] not in {"+", "-"} or not result["cds"]:
        return reject("invalid_coding_chain")
    # Recompute from the actual genome before considering any extension.
    if _chain_sequence(result, genome) != result.get("sequence", ""):
        return reject("source_sequence_mismatch")
    baseline = validate_model(copy.deepcopy(result), genome, genetic_code, params)
    if set(baseline.get("problems", [])) - _RECHECKABLE:
        return reject("strict_source_validation_failed")
    if baseline["cds"] != result["cds"] or baseline["sequence"] != result["sequence"]:
        return reject("input_not_already_strictly_checked")
    donor = str(donor_protein or "").upper().removesuffix("*")
    if not donor or set(donor) - _AMINO_ACIDS:
        return reject("missing_or_ambiguous_donor_protein")
    begin, end, length = (result.get(key) for key in ("query_start", "query_end", "query_length"))
    if any(type(value) is not int for value in (begin, end, length)) or not 0 <= begin < end <= length or length != len(donor):
        return reject("invalid_donor_query_coordinates")
    internal_unaligned = (end - begin) / length - result["coverage"]
    if internal_unaligned > 1e-10:
        return reject("internal_donor_gap")
    if "missing_start" in needed and begin == 0:
        return reject("donor_n_terminus_already_aligned")
    c_overhang = "missing_stop" in needed and end == length
    if c_overhang and not max_overhang:
        return reject("donor_c_terminus_already_aligned")
    codons = CodonTable.unambiguous_dna_by_id[genetic_code]
    n_choices, c_sequence = [""], ""
    if "missing_start" in needed:
        upstream, bounded = _flank(result, genome, maximum, "n")
        n_choices = []
        for offset in range(len(upstream) - 3, -1, -3):
            codon = upstream[offset:offset + 3]
            if set(codon) - set("ACGT") or codon in codons.stop_codons:
                audit["search_limits"].append("upstream_ambiguity_or_in_frame_stop")
                break
            if codon in codons.start_codons:
                n_choices.append(upstream[offset:])
        if not n_choices:
            return reject("no_in_frame_start_before_contig_boundary" if bounded else "no_in_frame_start_within_bound")
    if len(n_choices) > candidate_limit:
        return reject("start_candidate_budget_exceeded")
    if "missing_stop" in needed:
        c_bound = min(maximum, 3 * (max_overhang + 1)) if c_overhang else maximum
        downstream, bounded = _flank(result, genome, c_bound, "c")
        for offset in range(0, len(downstream), 3):
            codon = downstream[offset:offset + 3]
            if set(codon) - set("ACGT"):
                return reject("downstream_assembly_ambiguity")
            if codon in codons.stop_codons:
                c_sequence = downstream[:offset + 3]
                break
        if not c_sequence:
            return reject("no_in_frame_stop_before_contig_boundary" if bounded else "no_in_frame_stop_within_bound")
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score, aligner.mismatch_score = 2., -1.
    aligner.open_gap_score, aligner.extend_gap_score = -5., -.5
    eligible, uncertain = [], []
    for n_sequence in n_choices:
        candidate = copy.deepcopy(model)
        blocks = candidate["cds"] = [list(block) for block in candidate["cds"]]
        if candidate["strand"] == "+":
            blocks[0][0] -= len(n_sequence)
            blocks[-1][1] += len(c_sequence)
        else:
            blocks[0][1] += len(n_sequence)
            blocks[-1][0] -= len(c_sequence)
        dna = n_sequence + model["sequence"] + c_sequence
        attempt = {"cds": copy.deepcopy(blocks), "n_extension_dna": n_sequence, "c_extension_dna": c_sequence,
                   "sequence_sha256": _sha(dna), "reasons": []}
        audit["attempts"].append(attempt)
        protein = str(Seq(dna).translate(table=genetic_code)).removesuffix("*")
        # A code-permitted alternative initiation codon supplies methionine
        # at initiation; this changes only the alignment's amino-acid view.
        if dna[:3] in codons.start_codons and donor.startswith("M") and protein:
            attempt["alternative_initiation_normalized_for_alignment"] = protein[0] != "M"
            protein = "M" + protein[1:]
        if "*" in protein or set(protein) - _AMINO_ACIDS:
            attempt["reasons"].append("disrupted_extension_translation")
            continue
        n_added, c_added = len(n_sequence) // 3, max(0, len(c_sequence) // 3 - 1)
        if n_added > begin or (not c_overhang and c_added > length - end):
            attempt["reasons"].append("extension_longer_than_missing_donor_terminus")
            continue
        if min(len(protein), len(donor)) / max(len(protein), len(donor)) < params["minimum_coverage"]:
            attempt["reasons"].append("low_reciprocal_coverage_after_extension")
            continue
        if len(protein) * len(donor) > cell_limit:
            attempt["reasons"].append("alignment_cell_budget_exceeded")
            continue
        evidence = _protein_evidence(protein, donor, n_added, c_added, begin, end, aligner, c_overhang, alignment_limit)
        attempt["alignment"] = evidence
        attempt["reasons"].extend(evidence["reasons"])
        if evidence["identity"] < params["minimum_identity"]:
            attempt["reasons"].append("low_identity_after_extension")
        if min(evidence["coverage"], evidence["recipient_coverage"]) < params["minimum_coverage"]:
            attempt["reasons"].append("low_reciprocal_coverage_after_extension")
        uncertain_reasons = {"ambiguous_donor_terminal_alignment", "donor_alignment_budget_exceeded"}
        if "donor_alignment_budget_exceeded" in attempt["reasons"]:
            # Unexamined optimal paths may change these terminal projections.
            # Keep the reported failures but never infer a uniquely better
            # start from a traceback budget that was insufficient to decide.
            uncertain_reasons.update({"donor_termini_not_aligned", "internal_donor_gap",
                                      "n_extension_not_supported_by_missing_donor_terminus",
                                      "c_extension_not_supported_by_missing_donor_terminus",
                                      "short_c_overhang_did_not_follow_donor_terminal",
                                      "n_extension_has_no_matching_residue", "c_extension_has_no_matching_residue"})
        if set(attempt["reasons"]) - uncertain_reasons:
            continue
        candidate.update(coverage=evidence["coverage"], identity=evidence["identity"],
                         query_start=evidence["query_start"], query_end=evidence["query_end"],
                         query_span_coverage=(evidence["query_end"] - evidence["query_start"]) / len(donor),
                         alignment_source="terminal_completion_global")
        checked = validate_model(candidate, genome, genetic_code, params)
        if checked["problems"] or checked["cds"] != blocks or checked["sequence"] != dna:
            attempt["reasons"].extend(checked["problems"] or ["strict_extension_chain_mismatch"])
            continue
        (uncertain if attempt["reasons"] else eligible).append((evidence["score"], checked, attempt))
    if not eligible:
        if uncertain:
            audit["status"] = "ambiguous"
            audit["reasons"].append("donor_terminal_alignment_is_ambiguous")
            return result
        audit["reasons"].append("no_strict_donor_supported_extension")
        return result
    eligible.sort(key=lambda row: row[0], reverse=True)
    if uncertain and max(row[0] for row in uncertain) >= eligible[0][0] - margin:
        audit["status"] = "ambiguous"
        audit["reasons"].append("alternative_starts_within_score_margin")
        return result
    if len(eligible) > 1 and eligible[0][0] - eligible[1][0] <= margin:
        audit["status"] = "ambiguous"
        audit["reasons"].append("alternative_starts_within_score_margin")
        return result
    _, completed, selected = eligible[0]
    audit.update(status="completed", selected_cds=selected["cds"], selected_sequence_sha256=selected["sequence_sha256"],
                 score_margin=eligible[0][0] - eligible[1][0] if len(eligible) > 1 else None)
    completed["terminal_completion"] = audit
    completed["partial_evidence"] = {"partial": False, "completed_from_partial": True,
                                     "original_missing_termini": sorted(needed)}
    return completed

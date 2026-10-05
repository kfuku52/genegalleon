"""Advisory evidence axes for predicted CDS; these do not change admission."""
import math

from Bio.Data import CodonTable


def model_quality(model, genetic_code=1):
    sequence = model.get("sequence", "").upper()
    start = sequence[:3] if len(sequence) >= 3 else None
    flags = []
    allowed = CodonTable.unambiguous_dna_by_id[genetic_code].start_codons
    if start in allowed and start != "ATG":
        flags.append("alternative_start")
    terminal = {"n_aligned": None, "c_aligned": None, "query_length": model.get("query_length")}
    alignment = {"aligned_query_fraction": None, "query_span_fraction": None,
                 "internal_unaligned_query_fraction": None}
    if all(model.get(key) is not None for key in ("query_start", "query_end", "query_length")):
        begin, end, length = (model[key] for key in ("query_start", "query_end", "query_length"))
        if any(type(value) is not int for value in (begin, end, length)) or not 0 <= begin < end <= length:
            raise ValueError("Invalid terminal alignment coordinates")
        terminal.update(n_aligned=begin == 0, c_aligned=end == length)
        alignment["query_span_fraction"] = (end - begin) / length
        if begin:
            flags.append("donor_n_terminus_unaligned")
        if end != length:
            flags.append("donor_c_terminus_unaligned")
    else:
        flags.append("terminal_alignment_not_assessed")
    coverage = model.get("coverage")
    if coverage is not None:
        if isinstance(coverage, bool) or not isinstance(coverage, (int, float)) or not math.isfinite(coverage) or not 0 <= coverage <= 1:
            raise ValueError("Invalid aligned query coverage")
        alignment["aligned_query_fraction"] = float(coverage)
        if alignment["query_span_fraction"] is not None:
            internal = alignment["query_span_fraction"] - coverage
            if internal < -1e-12:
                raise ValueError("Aligned query coverage exceeds query span")
            alignment["internal_unaligned_query_fraction"] = max(0.0, internal)
            if internal > 1e-12:
                flags.append("donor_internal_unaligned_query")
    evidence = model.get("support") or ([model["evidence"]] if model.get("evidence") else [])
    donor_species = sorted({row["donor"] for row in evidence if row.get("donor")})
    if len(donor_species) == 1:
        flags.append("single_donor_species")
    inconsistent = sum(row.get("expected_strand") in {"+", "-"}
                       and row["expected_strand"] != model.get("strand") for row in evidence)
    if inconsistent:
        flags.append("expected_strand_conflict")
    return {"start_codon": start, "code_allows_start": start in allowed,
            "translation_initiation": "not_established", "terminal_alignment": terminal,
            "query_alignment": alignment,
            "native_terminal_completeness": "not_established", "donor_species": donor_species,
            "strand_conflicting_support": inconsistent, "flags": sorted(flags)}

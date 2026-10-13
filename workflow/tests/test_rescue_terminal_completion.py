"""Terminal repair must follow source DNA, donor termini, and strict validation."""

import copy
import hashlib

import pytest
from Bio.Seq import Seq

from workflow.support.rescue_gene_models import validate_model
from workflow.support.rescue_terminal_completion import complete_terminals

PARAMS = {"minimum_coverage": .95, "minimum_identity": .6, "max_intron": 1000}


class Genome:
    def __init__(self, dna):
        self.dna = dna
        self.calls = []

    def fetch(self, seqid, start, end):
        assert seqid == "chr" and 0 <= start <= end <= len(self.dna)
        self.calls.append((start, end))
        return self.dna[start:end]

    def get_reference_length(self, seqid):
        assert seqid == "chr"
        return len(self.dna)


def fixture(upstream, core, downstream, donor, begin, end, strand="+", code=1, **changes):
    dna = upstream + core + downstream
    start, stop = len(upstream), len(upstream) + len(core)
    if strand == "-":
        start, stop = len(dna) - stop, len(dna) - start
        dna = str(Seq(dna).reverse_complement())
    genome = Genome(dna)
    model = {"id": "candidate", "seqid": "chr", "strand": strand,
             "cds": [[start, stop, 0]], "frameshift": False,
             "coverage": (end - begin) / len(donor), "identity": 1.,
             "query_start": begin, "query_end": end, "query_length": len(donor),
             "paf": "unaltered source miniprot PAF"}
    model.update(changes)
    return validate_model(model, genome, code, PARAMS), genome


def complete(model, genome, donor, code=1, **params):
    return complete_terminals(model, genome, code, {**PARAMS, **params}, donor, validate_model)


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("terminus", ["n", "c", "both"])
def test_recover_missing_donor_termini_without_changing_source(strand, terminus):
    upstream = "TAAATG" if terminus in {"n", "both"} else "TAA"
    core = ("" if terminus in {"n", "both"} else "ATG") + "AAAGGA"
    core += "" if terminus in {"c", "both"} else "CCCTAA"
    downstream = "CCCTAA" if terminus in {"c", "both"} else "AAA"
    model, genome = fixture(upstream, core, downstream, "MKGP", int(terminus in {"n", "both"}),
                            3 if terminus in {"c", "both"} else 4, strand)
    original, dna = copy.deepcopy(model), genome.dna
    checked = complete(model, genome, "MKGP")
    assert checked["sequence"] == "ATGAAAGGACCCTAA"
    assert checked["problems"] == [] and checked["coverage"] == 1.
    assert checked["terminal_completion"]["status"] == "completed"
    assert checked["terminal_completion"]["original_alignment"]["paf"] == model["paf"]
    assert checked["alignment_source"] == "terminal_completion_global"
    assert checked["partial_evidence"]["partial"] is False
    assert checked["terminal_completion"]["selected_sequence_sha256"] == hashlib.sha256(
        checked["sequence"].encode()).hexdigest()
    assert model == original and genome.dna == dna


@pytest.mark.parametrize("strand", ["+", "-"])
def test_preserve_split_codon_phases_and_splice_junctions(strand):
    first, intron, second = "AAAGG", "GT" + "CCC" * 4 + "AG", "ACCCTAA"
    dna = "TAAATG" + first + intron + second
    blocks = [[6, 11, 0], [11 + len(intron), len(dna), 1]]
    if strand == "-":
        blocks = [[len(dna) - end, len(dna) - start, phase] for start, end, phase in blocks]
        dna = str(Seq(dna).reverse_complement())
    genome = Genome(dna)
    model = validate_model({"id": "spliced", "seqid": "chr", "strand": strand, "cds": blocks,
                            "frameshift": False, "coverage": .75, "identity": 1.,
                            "query_start": 1, "query_end": 4, "query_length": 4}, genome, 1, PARAMS)
    checked = complete(model, genome, "MKGP")
    assert checked["problems"] == [] and checked["sequence"] == "ATGAAAGGACCCTAA"
    assert checked["cds"][1] == model["cds"][1]
    assert [block[2] for block in checked["cds"]] == [0, 1]
    assert checked["cds"][0][1 if strand == "+" else 0] == model["cds"][0][1 if strand == "+" else 0]


@pytest.mark.parametrize("problem", ["frameshift", "invalid_phase", "invalid_splice", "noncanonical_splice",
                                     "internal_stop", "assembly_gap_or_ambiguity", "assembly_gap_within_model_span",
                                     "sequence_mismatch", "incomplete_frame", "intron_too_long"])
def test_other_failures_cannot_be_repaired_or_masked(problem):
    model, genome = fixture("TAAATG", "AAAGGACCCTAA", "", "MKGP", 1, 4)
    model["problems"].append(problem)
    original = copy.deepcopy(model)
    checked = complete(model, genome, "MKGP")
    assert checked["sequence"] == model["sequence"] and checked["cds"] == model["cds"]
    assert problem in checked["problems"]
    assert checked["terminal_completion"]["reasons"] == ["non_terminal_failure_preserved"]
    assert checked["partial_evidence"]["representative_eligible"] is False and model == original


def test_source_dna_mismatch_is_not_overwritten_by_extension():
    model, genome = fixture("TAAATG", "AAAGGACCCTAA", "", "MKGP", 1, 4)
    model["sequence"] = "CCCGGACCCTAA"
    checked = complete(model, genome, "MKGP")
    assert checked["sequence"] == model["sequence"]
    assert checked["terminal_completion"]["reasons"] == ["source_sequence_mismatch"]


def test_internal_donor_insertion_cannot_be_called_terminal_recovery():
    model, genome = fixture("TAAATG", "AAAGGACCCTAA", "", "MKRG P".replace(" ", ""), 1, 5)
    model["coverage"] = 3 / 5
    checked = complete(model, genome, "MKRGP")
    assert checked["cds"] == model["cds"] and "missing_start" in checked["problems"]
    assert checked["terminal_completion"]["reasons"] == ["internal_donor_gap"]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_upstream_in_frame_stop_blocks_distant_start(strand):
    model, genome = fixture("ATGTAA", "AAAGGACCCTAA", "", "MKGP", 1, 4, strand)
    checked = complete(model, genome, "MKGP")
    assert checked["cds"] == model["cds"]
    assert "upstream_ambiguity_or_in_frame_stop" in checked["terminal_completion"]["search_limits"]


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("terminus", ["n", "c"])
def test_contig_boundary_preserves_partial_without_invalid_fetch(strand, terminus):
    core = "AAAGGACCCTAA" if terminus == "n" else "ATGAAAGGA"
    begin, end = (1, 4) if terminus == "n" else (0, 3)
    model, genome = fixture("", core, "", "MKGP", begin, end, strand)
    checked = complete(model, genome, "MKGP")
    assert checked["cds"] == model["cds"]
    assert any("contig_boundary" in reason for reason in checked["terminal_completion"]["reasons"])


@pytest.mark.parametrize("strand", ["+", "-"])
def test_nonstandard_stop_code_and_source_internal_stop(strand):
    model, genome = fixture("TAA", "ATGTGAAAA", "CCCTAA", "MWKP", 0, 3, strand, 4)
    checked = complete(model, genome, "MWKP", 4)
    assert checked["problems"] == [] and checked["sequence"] == "ATGTGAAAACCCTAA"
    standard, _ = fixture("TAA", "ATGTGAAAA", "CCCTAA", "MWKP", 0, 3, strand, 1)
    unchanged = complete(standard, genome, "MWKP", 1)
    assert "internal_stop" in unchanged["problems"]
    assert unchanged["cds"] == standard["cds"]


def test_code_permitted_alternative_start_is_actual_genomic_dna():
    model, genome = fixture("TAAGTG", "AAAGGACCCTAA", "", "MKGP", 1, 4, code=11)
    checked = complete(model, genome, "MKGP", 11)
    assert checked["problems"] == [] and checked["sequence"] == "GTGAAAGGACCCTAA"
    assert checked["terminal_completion"]["attempts"][0]["alternative_initiation_normalized_for_alignment"]


@pytest.mark.parametrize("extra", [1, 2])
@pytest.mark.parametrize("strand", ["+", "-"])
def test_short_target_specific_tail_requires_unchanged_donor_terminal(extra, strand):
    core, donor = "ATG" + "AAA" * 50, "M" + "K" * 50
    model, genome = fixture("TAA", core, "GCT" * extra + "TAA", donor, 0, len(donor), strand)
    checked = complete(model, genome, donor)
    assert checked["problems"] == [] and checked["sequence"] == core + "GCT" * extra + "TAA"
    evidence = checked["terminal_completion"]["attempts"][0]["alignment"]
    assert evidence["c_terminal_support"] == "short_genomic_overhang_after_aligned_donor"
    assert evidence["unaligned_c_overhang_residues"] == extra


def test_short_overhang_cannot_pass_reciprocal_coverage_on_tiny_donor():
    model, genome = fixture("TAA", "ATGAAA", "GCTGCTTAA", "MK", 0, 2)
    checked = complete(model, genome, "MK")
    assert checked["cds"] == model["cds"]
    assert "low_reciprocal_coverage_after_extension" in checked["terminal_completion"]["attempts"][0]["reasons"]


def test_long_target_specific_tail_or_disabled_overhang_remains_partial():
    core, donor = "ATG" + "AAA" * 50, "M" + "K" * 50
    model, genome = fixture("TAA", core, "GCT" * 3 + "TAA", donor, 0, len(donor))
    checked = complete(model, genome, donor)
    assert checked["cds"] == model["cds"] and "missing_stop" in checked["problems"]
    disabled = complete(model, genome, donor, terminal_max_unaligned_c_overhang=0)
    assert disabled["terminal_completion"]["reasons"] == ["donor_c_terminus_already_aligned"]


@pytest.mark.parametrize("terminus", ["n", "c"])
def test_assembly_ambiguity_cannot_be_skipped(terminus):
    if terminus == "n":
        model, genome = fixture("ATGNNN", "AAAGGACCCTAA", "", "MKGP", 1, 4)
    else:
        model, genome = fixture("TAA", "ATGAAAGGA", "NNNCCCTAA", "MKGP", 0, 3)
    checked = complete(model, genome, "MKGP")
    assert checked["cds"] == model["cds"] and checked["terminal_completion"]["status"] != "completed"


def test_unrelated_downstream_orf_has_no_donor_terminal_support():
    model, genome = fixture("TAA", "ATGAAA", "AAAAAATAA", "MKGP", 0, 2)
    checked = complete(model, genome, "MKGP")
    assert checked["cds"] == model["cds"]
    assert "low_identity_after_extension" in checked["terminal_completion"]["attempts"][0]["reasons"]


def test_alternative_start_candidates_without_decisive_donor_score_stay_proposals():
    donor = "MMM" + "A" * 100
    model, genome = fixture("TAA" + "ATG" * 3, "GCT" * 100 + "TAA", "", donor, 3, len(donor))
    checked = complete(model, genome, donor)
    assert checked["terminal_completion"]["status"] == "ambiguous"
    assert checked["cds"] == model["cds"] and "missing_start" in checked["problems"]
    assert len(checked["terminal_completion"]["attempts"]) == 3


def test_optimal_alignment_budget_cannot_manufacture_unique_start_support():
    donor = "MMM" + "A" * 100
    model, genome = fixture("TAA" + "ATG" * 3, "GCT" * 100 + "TAA", "", donor, 3, len(donor))
    checked = complete(model, genome, donor, terminal_max_optimal_alignments=1)
    assert checked["terminal_completion"]["status"] == "ambiguous"
    assert checked["cds"] == model["cds"] and checked["sequence"] == model["sequence"]
    assert any("donor_alignment_budget_exceeded" in attempt["reasons"]
               for attempt in checked["terminal_completion"]["attempts"])


def test_extension_bound_and_alignment_budget_do_not_relax_guards():
    model, genome = fixture("TAAATG" + "AAA" * 100, "GGACCCTAA", "", "M" + "K" * 100 + "GP", 101, 103)
    checked = complete(model, genome, "M" + "K" * 100 + "GP")
    assert checked["cds"] == model["cds"] and checked["terminal_completion"]["status"] != "completed"
    model, genome = fixture("TAAATG", "AAAGGACCCTAA", "", "MKGP", 1, 4)
    limited = complete(model, genome, "MKGP", terminal_max_alignment_cells=1)
    assert limited["terminal_completion"]["attempts"][0]["reasons"] == ["alignment_cell_budget_exceeded"]


def test_candidate_budget_is_bounded_and_audited():
    donor = "M" * 33 + "GP"
    model, genome = fixture("TAA" + "ATG" * 33, "GGACCCTAA", "", donor, 33, 35)
    checked = complete(model, genome, donor)
    assert checked["terminal_completion"]["reasons"] == ["start_candidate_budget_exceeded"]
    assert checked["terminal_completion"]["attempts"] == []


def test_intact_model_is_unchanged_without_extra_alignment_work():
    model, genome = fixture("TAA", "ATGAAAGGACCCTAA", "", "MKGP", 0, 4)
    genome.calls.clear()
    checked = complete(model, genome, "MKGP")
    assert checked == model and genome.calls == []
    checked["problems"].append("outside_expected_synteny_interval")
    assert model["problems"] == []


def test_incomplete_frame_stays_explicit_partial_without_inventing_phase():
    model, genome = fixture("TAA", "ATGAAAA", "TAA", "MKGP", 0, 4)
    original = copy.deepcopy(model)
    genome.calls.clear()
    checked = complete(model, genome, "MKGP")
    assert checked["cds"] == model["cds"] and checked["sequence"] == model["sequence"]
    assert checked["problems"] == ["incomplete_frame"]
    assert checked["partial_evidence"]["partial"] and not checked["partial_evidence"]["representative_eligible"]
    assert model == original and genome.calls == []


@pytest.mark.parametrize("params", [{"terminal_max_extension": 0}, {"terminal_max_candidates": True},
                                    {"terminal_max_alignment_cells": -1}, {"terminal_minimum_margin": float("nan")},
                                    {"terminal_max_unaligned_c_overhang": 3}, {"terminal_max_optimal_alignments": 0},
                                    {"terminal_max_optimal_alignments": 33}])
def test_invalid_bounded_completion_parameters_raise(params):
    model, genome = fixture("TAAATG", "AAAGGACCCTAA", "", "MKGP", 1, 4)
    with pytest.raises(ValueError):
        complete(model, genome, "MKGP", **params)

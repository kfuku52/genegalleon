"""Nearby boundaries must recover actual genomic ORFs, not hide defects."""
import copy

import pytest
from Bio.Seq import Seq

from workflow.support import rescue_boundary_revalidation as boundaries
from workflow.support.rescue_boundary_revalidation import revalidate_boundaries, validated_boundary_candidates
from workflow.support.rescue_gene_models import validate_model
from workflow.support.rescue_terminal_completion import complete_terminals

PARAMS = {"minimum_coverage": .95, "minimum_identity": .5, "max_intron": 1000}


class Genome:
    def __init__(self, dna):
        self.dna = dna

    def fetch(self, seqid, start, end):
        assert seqid == "chr" and 0 <= start <= end <= len(self.dna)
        return self.dna[start:end]

    def get_reference_length(self, seqid):
        assert seqid == "chr"
        return len(self.dna)


def checked(dna, blocks, donor, strand="+", begin=0):
    if strand == "-":
        blocks = [[len(dna) - b, len(dna) - a, p] for a, b, p in blocks]
        dna = str(Seq(dna).reverse_complement())
    genome = Genome(dna)
    model = {"seqid": "chr", "strand": strand, "cds": blocks, "frameshift": False,
             "coverage": (len(donor) - begin) / len(donor), "identity": 1.,
             "query_start": begin, "query_end": len(donor), "query_length": len(donor),
             "paf": "original predictor alignment",
             "evidence": {"donor": "Relative", "query": "gene", "target": "Target"}}
    return validate_model(model, genome, 1, PARAMS), genome


def recheck(model, genome, donor, **params):
    return revalidate_boundaries(model, genome, 1, {**PARAMS, **params}, donor, validate_model)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_one_genomic_tail_path_survives_ambiguous_donor_tracebacks(strand):
    donor = "M" + "A" * 90 + "NN"
    core = "ATG" + "GCT" * 90 + "AACAAC"
    model, genome = checked("TAA" + core + "AACAATTAA", [[3, 3 + len(core), 0]], donor, strand)
    prior = complete_terminals(model, genome, 1, PARAMS, donor, validate_model)
    assert prior["terminal_completion"]["status"] == "ambiguous"
    original = copy.deepcopy(prior)
    result = recheck(prior, genome, donor)
    assert result["problems"] == [] and result["sequence"] == core + "AACAATTAA"
    assert result["boundary_revalidation"]["status"] == "completed"
    assert result["boundary_revalidation"]["attempts"][0]["alignment"]["optimal_alignments_reviewed"] > 1
    assert result["paf"] == model["paf"] and prior == original


@pytest.mark.parametrize("strand", ["+", "-"])
def test_recover_canonical_nearby_splice_and_recompute_split_phases(strand):
    first, second = "ATG" + "GAC" * 9, "GAA" * 90 + "TAA"
    intron = "GT" + "C" * 20 + "AG"
    dna = "TAA" + first + intron + second
    blocks = [[3, 3 + len(first) + 3, 0], [3 + len(first) + len(intron), len(dna), 0]]
    donor = "M" + "D" * 9 + "E" * 90
    model, genome = checked(dna, blocks, donor, strand)
    assert model["problems"] == ["noncanonical_splice"]
    original = copy.deepcopy(model)
    result = recheck(model, genome, donor)
    assert result["problems"] == [] and result["sequence"] == first + second
    assert result["boundary_revalidation"]["status"] == "completed" and result["identity"] == 1
    assert model == original


@pytest.mark.parametrize("strand", ["+", "-"])
def test_recover_donor_supported_start_from_all_nearby_positions(strand):
    donor = "M" + "A" * 90
    core = "GCT" * 90 + "TAA"
    model, genome = checked("TAAATG" + core, [[6, 6 + len(core), 0]], donor, strand, begin=1)
    result = recheck(model, genome, donor)
    assert result["problems"] == [] and result["sequence"] == "ATG" + core
    assert result["query_start"] == 0 and result["paf"] == model["paf"]


@pytest.mark.parametrize("problem", ["frameshift", "invalid_phase", "incomplete_frame",
                                    "internal_stop", "assembly_gap_or_ambiguity",
                                    "assembly_gap_within_model_span", "intron_too_long"])
def test_real_disruptions_are_preserved(problem):
    donor = "M" + "A" * 90
    core = "GCT" * 90 + "TAA"
    model, genome = checked("TAAATG" + core, [[6, 6 + len(core), 0]], donor, begin=1)
    model["problems"].append(problem)
    original = copy.deepcopy(model)
    result = recheck(model, genome, donor)
    assert result == original and "boundary_revalidation" not in result


def test_long_unaligned_tail_remains_partial():
    donor = "M" + "A" * 90
    core = "ATG" + "GCT" * 90
    model, genome = checked("TAA" + core + "GCT" * 21 + "TAA", [[3, 3 + len(core), 0]], donor)
    result = recheck(model, genome, donor)
    assert result["cds"] == model["cds"] and "missing_stop" in result["problems"]
    assert result["boundary_revalidation"]["status"] == "rejected"


def test_source_sequence_mismatch_is_not_overwritten():
    donor = "M" + "A" * 90
    core = "GCT" * 90 + "TAA"
    model, genome = checked("TAAATG" + core, [[6, 6 + len(core), 0]], donor, begin=1)
    model["sequence"] = "AAA" + model["sequence"][3:]
    result = recheck(model, genome, donor)
    assert result["boundary_revalidation"]["reasons"] == ["source_sequence_mismatch"]


def test_intact_prediction_does_not_trigger_a_second_search():
    donor = "M" + "A" * 90
    core = "ATG" + "GCT" * 90 + "TAA"
    model, genome = checked(core, [[0, len(core), 0]], donor)
    assert recheck(model, genome, donor) is model

def test_distinct_complete_start_paths_remain_ambiguous():
    donor = "M" + "A" * 90
    core = "GCT" * 90 + "TAA"
    model, genome = checked("TAAATGATGATG" + core, [[12, 12 + len(core), 0]], donor, begin=1)
    result = recheck(model, genome, donor)
    assert result["cds"] == model["cds"] and "missing_start" in result["problems"]
    assert result["boundary_revalidation"]["status"] == "ambiguous"
    assert result["boundary_revalidation"]["reasons"] == ["alternative_genomic_boundary_paths_within_score_margin"]
    assert result["boundary_revalidation"]["supports_intact_homolog"] is True
    assert len(result["boundary_revalidation"]["intact_donor_supported_paths"]) > 1


@pytest.mark.parametrize('strand', ['+', '-'])
def test_complete_boundary_alternatives_enter_final_selection_without_rewriting_source(strand):
    from workflow.support.rescue_gene_models import consolidate
    donor = 'M' + 'A' * 90
    core = 'GCT' * 90 + 'TAA'
    model, genome = checked('TAAATGATGATG' + core, [[12, 12 + len(core), 0]], donor, strand, begin=1)
    source = recheck(model, genome, donor)
    original = copy.deepcopy(source)
    source['query'] = 'region'
    candidates = validated_boundary_candidates(source, genome, 1, PARAMS, donor, validate_model)
    assert len(candidates) > 1 and all(not row['problems'] for row in candidates)
    consolidate(candidates, [], 'Target')
    primary, = [row for row in candidates if row['status'] == 'accepted']
    assert primary['path_selection']['representative_status'] == 'ambiguous'
    assert len(primary['alternative_coding_paths']) == len(candidates) - 1
    assert source == {**original, 'query': 'region'}
    assert all(row['boundary_candidate']['representative_status'] == 'ambiguous' for row in candidates)


@pytest.mark.parametrize('failure', ['budget', 'changed_source', 'disruption', 'missing_donor'])
def test_boundary_presence_flag_alone_cannot_admit_candidates(failure):
    donor = 'M' + 'A' * 90
    core = 'GCT' * 90 + 'TAA'
    model, genome = checked('TAAATGATGATG' + core, [[12, 12 + len(core), 0]], donor, begin=1)
    source = recheck(model, genome, donor)
    assert source['boundary_revalidation']['supports_intact_homolog']
    if failure == 'budget':
        source['boundary_revalidation']['reasons'] = ['unresolved_donor_alignment_budget']
    elif failure == 'changed_source':
        genome.dna = genome.dna.replace('GCT', 'GAT')
    elif failure == 'disruption':
        source['problems'].append('frameshift')
    elif failure == 'missing_donor':
        donor = ''
    assert validated_boundary_candidates(source, genome, 1, PARAMS, donor, validate_model) == []


def test_unexamined_optimal_alignments_cannot_select_a_path():
    donor = "M" + "A" * 90 + "NN"
    core = "ATG" + "GCT" * 90 + "AACAAC"
    model, genome = checked("TAA" + core + "AACAATTAA", [[3, 3 + len(core), 0]], donor)
    result = recheck(model, genome, donor, terminal_max_optimal_alignments=1)
    assert result["cds"] == model["cds"] and "missing_stop" in result["problems"]
    assert result["boundary_revalidation"]["reasons"] == ["unresolved_donor_alignment_budget"]
    assert result["boundary_revalidation"]["supports_intact_homolog"] is True


@pytest.mark.parametrize("strand", ["+", "-"])
def test_recover_start_codon_split_across_short_first_exon(strand):
    first, intron, second = "CCCA", "GT" + "C" * 10 + "AG", "TG" + "GCT" * 90 + "TAA"
    dna = "TAA" + first + intron + second
    blocks = [[3, 7, 0], [7 + len(intron), len(dna), 2]]
    donor = "M" + "A" * 90
    model, genome = checked(dna, blocks, donor, strand)
    assert model["problems"] == ["missing_start"]
    result = recheck(model, genome, donor)
    assert result["problems"] == [] and result["sequence"] == "ATG" + "GCT" * 90 + "TAA"
    assert result["cds"][0][1] - result["cds"][0][0] == 1


def test_adjacent_boundary_choices_cannot_invert_a_middle_exon(monkeypatch):
    first, middle, last = "ATG" + "GCT" * 10, "AAA", "GCT" * 10 + "TAA"
    intron = "AA" + "C" * 10 + "GG"
    dna = first + intron + middle + intron + last
    a, b = len(first) + len(intron), len(first) + len(intron) + len(middle)
    blocks = [[0, len(first), 0], [a, b, 0], [b + len(intron), len(dna), 0]]
    donor = "M" + "A" * 10 + "K" + "A" * 10
    model, genome = checked(dna, blocks, donor)
    assert model["problems"] == ["noncanonical_splice"]
    monkeypatch.setattr(boundaries, "_splice_choices", lambda *args: ([[(0, 2)], [(-2, 0)]], None))
    result = recheck(model, genome, donor)
    assert result["cds"] == model["cds"]
    assert result["boundary_revalidation"]["status"] == "rejected"

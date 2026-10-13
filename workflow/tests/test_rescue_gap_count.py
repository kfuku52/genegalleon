"""Preserve complete genomic QC results while counting span ambiguities in C."""

import copy
import inspect

import pytest

from workflow.support import rescue_gene_models as rescue


class Genome:
    def __init__(self, sequence):
        self.sequence = sequence

    def fetch(self, seqid, start, end):
        assert seqid == "chr"
        return self.sequence[start:end]

    def get_reference_length(self, seqid):
        assert seqid == "chr"
        return len(self.sequence)


@pytest.fixture
def legacy_validator():
    # Freeze the former ambiguity count while exercising the complete validator,
    # including query metrics, raw phases, stop extension and splice checks.
    source = inspect.getsource(rescue.validate_model)
    current = 'len(span) - span.count("A") - span.count("C") - span.count("G") - span.count("T")'
    assert source.count(current) == 1
    namespace = dict(rescue.__dict__)
    exec(compile(source.replace(current, 'sum(base not in "ACGT" for base in span)'),
                 "legacy_rescue_validator", "exec"), namespace)
    return namespace["validate_model"]


PARAMS = {"minimum_coverage": .95, "minimum_identity": .5, "max_intron": 20000}


def spliced_model(payload, strand, *, first="ATGAAA", last="CCCTAA"):
    sequence = first + "GT" + payload + "AG" + last
    cds = [[0, 6, 0], [len(sequence) - 6, len(sequence), 0]]
    if strand == "-":
        # Preserve arbitrary Unicode intronic code points; coding and splice
        # residues stay ordinary DNA in both transcription directions.
        sequence = sequence.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]
        cds = [[len(sequence) - 6, len(sequence), 0], [0, 6, 0]]
    model = {"id": "prediction", "seqid": "chr", "strand": strand, "cds": cds,
             "frameshift": False, "coverage": 1., "identity": .8,
             "query": "query_A", "evidence": {"donor": "relative_A", "query": "protein_A",
                                               "owner": "original_gene", "rna_supported": True}}
    return model, Genome(sequence)


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("payload,ambiguities", [
    ("", 0), ("ACGTacgt" * 16, 0), ("NNnnry", 6), ("- .\t\r\n", 6), ("éβ🙂𝔸ß", 6),
])
def test_complete_qc_parity_for_intronic_ambiguities(legacy_validator, strand, payload, ambiguities):
    model, genome = spliced_model(payload, strand)
    before = copy.deepcopy(model)
    checked = rescue.validate_model(model, genome, 1, PARAMS)
    assert checked == legacy_validator(model, genome, 1, PARAMS)
    assert checked["assembly_ambiguous_bases"] == ambiguities
    assert checked["problems"] == (["assembly_gap_within_model_span"] if ambiguities else [])
    assert checked["sequence"] == "ATGAAACCCTAA"
    assert checked["evidence"] == before["evidence"]
    assert model == before


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("case", ["exonic_ambiguity", "internal_stop", "missing_stop", "raw_phase", "query_flags"])
def test_disrupted_genomic_and_query_qc_parity(legacy_validator, strand, case):
    model, genome = spliced_model("NN", strand,
                                   first="ATGNAA" if case == "exonic_ambiguity" else
                                   "ATGTAA" if case == "internal_stop" else "ATGAAA",
                                   last="CCCAAA" if case == "missing_stop" else "CCCTAA")
    if case == "raw_phase":
        model["cds"][1][2] = 1
    if case == "query_flags":
        model.update(coverage=.4, identity=.2, frameshift=True)
    before = copy.deepcopy(model)
    checked = rescue.validate_model(model, genome, 1, PARAMS)
    assert checked == legacy_validator(model, genome, 1, PARAMS)
    expected = {"exonic_ambiguity": "assembly_gap_or_ambiguity", "internal_stop": "internal_stop",
                "missing_stop": "missing_stop", "raw_phase": "invalid_phase", "query_flags": "frameshift"}
    assert expected[case] in checked["problems"]
    if case == "query_flags":
        assert {"frameshift", "low_coverage", "low_identity"} <= set(checked["problems"])
    assert model == before


@pytest.mark.parametrize("code", [1, 4])
def test_genetic_code_qc_parity(legacy_validator, code):
    model = {"seqid": "chr", "strand": "+", "cds": [[0, 12, 0]],
             "frameshift": False, "coverage": 1., "identity": 1.}
    genome = Genome("ATGTGACCCTAA")
    checked = rescue.validate_model(model, genome, code, PARAMS)
    assert checked == legacy_validator(model, genome, code, PARAMS)
    assert ("internal_stop" in checked["problems"]) == (code == 1)


def test_empty_chain_and_out_of_bounds_keep_original_behavior(legacy_validator):
    genome = Genome("ATGAAATAA")
    empty = {"seqid": "chr", "strand": "+", "cds": []}
    assert rescue.validate_model(empty, genome, 1, PARAMS) == legacy_validator(empty, genome, 1, PARAMS)
    model = {**empty, "cds": [[0, 12, 0]], "frameshift": False, "coverage": 1., "identity": 1.}
    for validator in (legacy_validator, rescue.validate_model):
        with pytest.raises(ValueError, match="Predicted CDS outside genome"):
            validator(model, genome, 1, PARAMS)

import copy
import math
from collections import OrderedDict

import pytest
from Bio.Seq import Seq

from workflow.support.rescue_gene_models import validate_model
from workflow.support.rescue_homology_policy import homology_problems
from workflow.support.rescue_raw_validation import RawValidationMemo, copy_validation_data


class Genome:
    def __init__(self, sequence="ATGAAATAA"):
        self.sequence = sequence
        self.fetches = 0

    def fetch(self, _seqid, start, end):
        self.fetches += 1
        return self.sequence[start:end]

    def get_reference_length(self, _seqid):
        return len(self.sequence)


PARAMS = {"minimum_coverage": .95, "minimum_identity": .5, "max_intron": 20000}


def model(**changes):
    return {"query": "first", "evidence": {"donor": "A", "query": "q"},
            "seqid": "chr", "strand": "+", "cds": [[0, 9, 0]],
            "frameshift": False, "coverage": 1.0, "identity": .8, **changes}


def signature(call):
    try:
        value = call()
        return "result", value, list(value)
    except Exception as error:
        return "exception", type(error), str(error)


@pytest.mark.parametrize("genomic_only", [False, True])
@pytest.mark.parametrize("changes", [
    {}, {"cds": [[0, 6, 0]]}, {"cds": [[0, 7, 0]]},
    {"frameshift": True}, {"coverage": .8}, {"identity": .3},
    {"cds": [[0, 9, 1]]}, {"strand": "-"}, {"cds": []},
    {"coverage": math.nan}, {"identity": math.inf}, {"coverage": None},
    {"strand": "?"}, {"strand": []}, {"cds": [[0, 3]]}, {"cds": [[-1, 9, 0]]},
    {"cds": [[0, 10, 0]]}, {"frameshift": 1},
    {"sequence": "old", "problems": ["old"], "assembly_ambiguous_bases": 123},
    {"genetic_code": 2},
])
def test_full_validator_output_order_and_exceptions(changes, genomic_only):
    genome = Genome()
    memo = RawValidationMemo(validate_model, genome, 1, PARAMS,
                             candidate_homology=homology_problems if genomic_only else None)
    for _ in range(2):
        row = model(**changes)
        expected = signature(lambda row=row: validate_model(copy.deepcopy(row), genome, 1, PARAMS))
        assert signature(lambda row=row: memo.validate(copy.deepcopy(row))) == expected


@pytest.mark.parametrize("row", [
    {"strand": "+", "cds": []}, {"strand": "?", "cds": []},
    {"cds": []}, {"strand": "+"}, {"strand": "+", "cds": [[0, 3, 0]]},
])
def test_key_does_not_preempt_original_errors_or_no_cds(row):
    genome = Genome()
    memo = RawValidationMemo(validate_model, genome, 1, PARAMS)
    expected = signature(lambda row=row: validate_model(copy.deepcopy(row), genome, 1, PARAMS))
    assert signature(lambda row=row: memo.validate(copy.deepcopy(row))) == expected


def test_alias_hit_preserves_current_unknowns_and_isolates_dna():
    genome = Genome()
    memo = RawValidationMemo(validate_model, genome, 1, PARAMS)
    first = model(unknown={"nested": [1]}, raw_prediction={"query": "first"})
    a = memo.validate(first)
    fetches = genome.fetches
    a["cds"][0][1] = -99
    a["problems"].append("poison")
    second = model(query="second", evidence={"donor": "B"},
                   unknown={"nested": [2]}, raw_prediction={"query": "second"})
    b = memo.validate(second)
    assert genome.fetches == fetches
    assert b == validate_model(second, genome, 1, PARAMS)
    assert b["unknown"] is second["unknown"] and b["evidence"] is second["evidence"]
    assert b["raw_prediction"] is second["raw_prediction"]
    b["cds"][0][1] = -42
    b["problems"].append("poison")
    assert memo.validate(second)["cds"] == [[0, 9, 0]]
    assert memo.validate(second)["problems"] == []


def test_genome_code_and_quality_scopes_are_independent():
    row = model()
    a = RawValidationMemo(validate_model, Genome(), 1, PARAMS)
    b = RawValidationMemo(validate_model, Genome("ATGAAATGA"), 2, PARAMS)
    assert a.validate(row)["sequence"] != b.validate(row)["sequence"]
    params = dict(PARAMS)
    c = RawValidationMemo(validate_model, Genome(), 1, params)
    assert c.validate(row)["problems"] == []
    params["minimum_identity"] = .9
    assert c.validate(row)["problems"] == ["low_identity"]
    assert c.diagnostics()["scope_invalidations"] == 1


@pytest.mark.parametrize("intron", ["GTAAAAG", "GCAAAAG", "ATAAAAC", "ATAAAAG", "GTNNNAG"])
@pytest.mark.parametrize("strand", ["+", "-"])
def test_splice_gap_and_transcription_order_are_not_relaxed(intron, strand):
    sequence = "ATG" + intron + "AAATAA"
    blocks = [[0, 3, 0], [10, 16, 0]]
    if strand == "-":
        sequence = str(Seq(sequence).reverse_complement())
        blocks = [[0, 6, 0], [13, 16, 0]]
    genome = Genome(sequence)
    row = model(cds=blocks, strand=strand)
    memo = RawValidationMemo(validate_model, genome, 1, PARAMS)
    expected = validate_model(row, genome, 1, PARAMS)
    for query in ("first", "another"):
        current = {**row, "query": query}
        assert memo.validate(current) == {**expected, "query": query}
    if intron == "ATAAAAG":
        assert "noncanonical_splice" in expected["problems"]
    if "N" in intron:
        assert "assembly_gap_within_model_span" in expected["problems"]


@pytest.mark.parametrize("field,value", [
    ("identity", .7), ("coverage", .9), ("frameshift", True),
    ("strand", "-"), ("cds", [[0, 6, 0]]), ("seqid", "other"),
])
def test_changed_validator_input_never_borrows_qc(field, value):
    genome = Genome()
    memo = RawValidationMemo(validate_model, genome, 1, PARAMS)
    memo.validate(model())
    memo.validate(model(**{field: value}))
    assert memo.diagnostics()["hits"] == 0


def test_bounded_entries_and_bytes_and_lru_eviction():
    memo = RawValidationMemo(validate_model, Genome(), 1, PARAMS, max_entries=2)
    for seqid in ("a", "b", "a", "c", "b"):
        memo.validate(model(seqid=seqid))
    stats = memo.diagnostics()
    assert stats["entries"] == stats["peak_entries"] == 2
    assert stats["hits"] == 1 and stats["evictions"] == 2
    small = RawValidationMemo(validate_model, Genome(), 1, PARAMS, max_bytes=1)
    small.validate(model())
    assert small.diagnostics()["entries"] == 0
    assert small.diagnostics()["oversized"] == 1
    assert stats["peak_accounted_bytes"] <= stats["max_bytes"]


def test_disabled_cache_and_future_validator_output_are_safe():
    def future(row, genome, code, params):
        return {**validate_model(row, genome, code, params), "new_derived": row["query"]}
    memo = RawValidationMemo(future, Genome(), 1, PARAMS)
    assert memo.validate(model())["new_derived"] == "first"
    assert memo.validate(model(query="second"))["new_derived"] == "second"
    assert not memo.diagnostics()["entries"]
    disabled = RawValidationMemo(validate_model, Genome(), 1, PARAMS, max_entries=0)
    disabled.validate(model())
    assert disabled.diagnostics()["bypasses"] == 1


def test_mapping_subclass_and_preexisting_field_order():
    row = OrderedDict(model(sequence="old", assembly_ambiguous_bases=12, problems=["old"]))
    genome = Genome()
    memo = RawValidationMemo(validate_model, genome, 1, PARAMS)
    assert signature(lambda: memo.validate(row)) == signature(lambda: validate_model(row, genome, 1, PARAMS))


@pytest.mark.parametrize("kwargs", [{"max_entries": -1}, {"max_bytes": -1}, {"max_entries": True}])
def test_invalid_memory_limits(kwargs):
    with pytest.raises(ValueError):
        RawValidationMemo(validate_model, Genome(), 1, PARAMS, **kwargs)


def test_overlay_copy_preserves_order_aliases_cycles_and_independent_mutables():
    blocks = [[0, 9, 0]]
    source = {"sequence": "ATG" * 10000, "cds": blocks, "also_cds": blocks,
              "tuple": (blocks,), "completion": OrderedDict([("end", [9]), ("start", [0])])}
    source["cycle"] = [source]
    result = copy_validation_data(source)
    assert list(result) == list(source)
    assert result["cds"] is result["also_cds"] is result["tuple"][0]
    assert result["cycle"][0] is result
    assert result["sequence"] is source["sequence"]
    assert type(result["completion"]) is OrderedDict
    assert list(result["completion"]) == ["end", "start"]
    result["cds"][0][1] = -1
    result["completion"]["end"].append(12)
    assert blocks == [[0, 9, 0]]
    assert source["completion"]["end"] == [9]
    second = copy_validation_data(source)
    assert second["cds"] == [[0, 9, 0]]
    assert second["cycle"][0] is second


def test_overlay_copy_delegates_subclasses_and_custom_copy_protocol():
    class Custom:
        def __deepcopy__(self, memo):
            return {"blocks": copy.deepcopy(self.blocks, memo)}

    class Blocks(list):
        pass

    blocks = Blocks([[0, 9, 0]])
    custom = Custom()
    custom.blocks = blocks
    source = {"blocks": blocks, "custom": custom, ("tuple", 2): b"bytes"}
    result = copy_validation_data(source)
    assert result == copy.deepcopy(source)
    assert type(result["blocks"]) is Blocks
    assert result["custom"]["blocks"] is result["blocks"]
    result["blocks"][0].append(3)
    assert blocks == [[0, 9, 0]]


def test_overlay_copy_preserves_float_objects_and_copy_errors():
    source = [math.nan, math.inf, -0.0, 2 ** 100, "α"]
    result = copy_validation_data(source)
    assert all(a is b for a, b in zip(result, source, strict=True))
    assert math.copysign(1, result[2]) == -1

    class Broken:
        def __deepcopy__(self, memo):
            raise ValueError("original copy failure")

    with pytest.raises(ValueError, match="original copy failure"):
        copy_validation_data({"nested": [Broken()]})


@pytest.mark.parametrize("sequence", ["ATGAAATAA", "ATGNNNTAA", "ATGAAATGA", "ATGTAATAA", "AAACCCTAA"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("code", [1, 2])
@pytest.mark.parametrize("strand", ["+", "-"])
def test_genomic_memo_rechecks_every_donor_metric_without_repeating_dna(sequence, reverse, code, strand):
    params = {**PARAMS, "minimum_coverage": .7, "minimum_broad_coverage": .9,
              "minimum_lower_coverage_identity": .7}
    pairs = [(1., .8), (.7, .7), (.85, .65), (.9, .5), (.69, .9), (.95, .49)]
    sequence = sequence if strand == "+" else str(Seq(sequence).reverse_complement())
    genome = Genome(sequence)
    memo = RawValidationMemo(validate_model, genome, code, params, candidate_homology=homology_problems)
    fetched = None
    for i, (coverage, identity) in enumerate(reversed(pairs) if reverse else pairs):
        row = model(query=f"donor{i}", coverage=coverage, identity=identity, strand=strand)
        expected = validate_model(row, Genome(sequence), code, params)
        actual = memo.validate(row)
        assert actual == expected and list(actual) == list(expected)
        assert actual["evidence"] is row["evidence"]
        if fetched is None:
            fetched = genome.fetches
        assert genome.fetches == fetched
        actual["cds"][0][1] = -99
        actual["problems"].append("poison")
    assert memo.diagnostics()["hits"] == len(pairs) - 1
    assert memo.diagnostics()["homology_evaluated_per_candidate"] is True
    params["minimum_identity"] = .95
    row = model(identity=.8, strand=strand)
    assert memo.validate(row) == validate_model(row, Genome(sequence), code, params)
    assert memo.diagnostics()["scope_invalidations"] == 1


@pytest.mark.parametrize("field,value", [(field, value) for field in ("coverage", "identity")
                         for value in (-.1, 1.1, math.nan, math.inf, True, "0.8", None)])
def test_genomic_memo_preserves_metric_errors_on_a_populated_cache(field, value):
    memo = RawValidationMemo(validate_model, Genome(), 1, PARAMS, candidate_homology=homology_problems)
    memo.validate(model())
    row = model(**{field: value})
    assert signature(lambda: memo.validate(row)) == signature(lambda: validate_model(row, Genome(), 1, PARAMS))
    assert memo.diagnostics()["hits"] == 0

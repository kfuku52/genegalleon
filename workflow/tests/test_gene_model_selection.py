import copy
import itertools
import json
import random

import pytest
from Bio.Align import PairwiseAligner
from Bio.Data import CodonTable

from workflow.support import gene_model_selection as selection
from workflow.support.gene_model_selection import pair_score, select_representatives


def protein(seed=1, length=90):
    randomizer = random.Random(seed)
    return "M" + "".join(randomizer.choices("ACDEFGHIKLNPQRSTVWY", k=length - 1))


def candidate(identifier, sequence, *, blocks=None, quality=None, **fields):
    return {"candidate_id": identifier, "source_transcript_id": identifier, "protein": sequence,
            "cds": "ATG" * len(sequence), "blocks": blocks or [[0, len(sequence) * 3, 0]],
            "quality": {"usable": True, "valid_orf": True} | (quality or {}), "origin": "original"} | fields


def catalog(species, *loci):
    return {"schema": 1, "species": species, "loci": [
        {"species": species, "gene_id": gene, "candidates": candidates} for gene, candidates in loci]}


def edge(a, b, gene_a="g", gene_b="g", **fields):
    return {"species_a": a, "gene_a": gene_a, "species_b": b, "gene_b": gene_b,
            "weight": 1, "evidence": "two_flank_synteny"} | fields


def extension_fixture():
    sequence = protein()
    return [catalog("A", ("g", [candidate("A_long", sequence + protein(5, 35)),
                                candidate("A_conserved", sequence)])),
            catalog("B", ("g", [candidate("B_main", sequence)])),
            catalog("C", ("g", [candidate("C_main", sequence)]))], [edge("A", "B"), edge("A", "C")]


def coding_candidate(identifier, sequence, **fields):
    codons = {}
    for triplet, residue in sorted(CodonTable.unambiguous_dna_by_id[1].forward_table.items()):
        codons.setdefault(residue, triplet)
    cds = "".join(codons[residue] for residue in sequence) + "TAA"
    return candidate(identifier, sequence, cds=cds, blocks=[[0, len(cds), 0]], **fields)


def legacy_pair_score(a, b, *, max_alignment_cells=25_000_000):
    """Pre-optimization native alignment and full-residue mapping oracle."""
    pa, pb = (str(item.get("protein", "")).upper().removesuffix("*") for item in (a, b))
    if not pa or not pb or len(pa) * len(pb) > max_alignment_cells:
        return selection.PairScore(0.0, 0.0, 0.0, 0.0, 0.0, 0, False, False, True)
    aligner = PairwiseAligner()
    aligner.mode, aligner.match_score, aligner.mismatch_score = "global", 2.0, -1.0
    aligner.open_gap_score, aligner.extend_gap_score = -5.0, -.5
    coordinates = aligner.align(pa, pb)[0].coordinates
    ma, mb = {}, {}
    aligned = matches = columns = 0
    loss_a = loss_b = False
    for index in range(coordinates.shape[1] - 1):
        sa, ea = int(coordinates[0, index]), int(coordinates[0, index + 1])
        sb, eb = int(coordinates[1, index]), int(coordinates[1, index + 1])
        wa, wb = ea - sa, eb - sb
        if wa and wb:
            aligned += wa
            matches += sum(x == y and x not in "XBZJ" for x, y in zip(pa[sa:ea], pb[sb:eb], strict=True))
        elif wa >= 9 and 0 < sb < len(pb):
            loss_a = True
        elif wb >= 9 and 0 < sa < len(pa):
            loss_b = True
        ma.update((sa + offset, columns + offset) for offset in range(wa))
        mb.update((sb + offset, columns + offset) for offset in range(wb))
        columns += max(wa, wb)

    def projected(item, mapping, size):
        result = set()
        for residue, phase in selection._junction_positions(item):
            if residue == size:
                result.add((columns, phase))
            elif residue in mapping:
                result.add((mapping[residue], phase))
        return result

    ja, jb = projected(a, ma, len(pa)), projected(b, mb, len(pb))
    union = ja | jb
    similarity = len(ja & jb) / len(union) if union else 1.0
    identity = matches / max(1, aligned)
    coverage_a, coverage_b = aligned / len(pa), aligned / len(pb)
    score = min(coverage_a, coverage_b) * (.85 * identity + .15 * similarity)
    return selection.PairScore(score, identity, coverage_a, coverage_b, similarity, aligned, loss_a, loss_b)


def pair_json(value):
    return json.dumps(value.as_dict(), sort_keys=True, separators=(",", ":"))


def test_pair_score_exhaustive_short_sequences_junctions_and_caps_match_legacy_json():
    sequences = ["".join(value) for length in range(5) for value in itertools.product("AX", repeat=length)]
    for pa, pb in itertools.product(sequences, repeat=2):
        a = candidate("a", pa.lower() + "*", protein_junctions=[[r, phase] for r in range(-1, len(pa) + 2) for phase in range(3)])
        b = candidate("b", pb, protein_junctions=[[r, phase] for r in range(-1, len(pb) + 2) for phase in (0, 2)])
        for cap in (9, 25_000_000):
            assert pair_json(pair_score(a, b, max_alignment_cells=cap)) == pair_json(
                legacy_pair_score(a, b, max_alignment_cells=cap))


def test_pair_score_random_gaps_ambiguous_residues_and_splice_phases_match_legacy_json():
    rng = random.Random(93829)
    for index in range(250):
        pa = "M" + "".join(rng.choices("ACDEFGHIKLNPQRSTVWYXBZJ", k=rng.randrange(20, 80)))
        cut = rng.randrange(1, len(pa) - 10)
        if index % 5 == 0:
            pb = pa
        elif index % 5 == 1:
            pb = pa[:cut] + "W" * 12 + pa[cut:]
        elif index % 5 == 2:
            pb = pa[:cut] + pa[cut + 10:]
        elif index % 5 == 3:
            pb = "Q" * 10 + pa + "W" * 9
        else:
            pb = "".join(rng.choices("ACDEFGHIKLNPQRSTVWYXBZJ", k=rng.randrange(20, 80)))
        a, b = candidate("a", pa), candidate("b", pb)
        for item, sequence in ((a, pa), (b, pb)):
            item["protein_junctions"] = [[rng.randrange(-2, len(sequence) + 3), rng.randrange(3)] for _ in range(7)]
            item["protein_junctions"] += [[len(sequence), 2], [0, 1]]
        assert pair_json(pair_score(a, b)) == pair_json(legacy_pair_score(a, b))
        assert pair_json(pair_score(b, a)) == pair_json(legacy_pair_score(b, a))
    for strand in ("+", "-"):
        a = candidate("a", "MABCXBZJ**", strand=strand, blocks=[[1, 8, 1], [12, 21, 2], [30, 43, 2]])
        b = candidate("b", "MABCXBZJ**", strand=strand, blocks=[[5, 13, 2], [24, 33, 0], [40, 53, 0]])
        assert pair_json(pair_score(a, b)) == pair_json(legacy_pair_score(a, b))


def test_pair_score_public_calls_observe_candidate_mutation_and_preserve_pre_junction_cap():
    a, b = candidate("a", "MABXJ"), candidate("b", "MABXJ")
    before = pair_score(a, b)
    a.update(protein="MABCQJ", protein_junctions=[[3, 1], [6, 2]])
    after = pair_score(a, b)
    assert pair_json(after) == pair_json(legacy_pair_score(a, b))
    assert before != after
    a["protein_junctions"] = [["not_an_integer", 0]]
    assert pair_json(pair_score(a, b, max_alignment_cells=1)) == pair_json(
        legacy_pair_score(a, b, max_alignment_cells=1))
    a["protein"] = ""
    assert pair_score(a, b).bounded


def test_identical_ascii_inputs_keep_native_results_or_errors():
    for value in range(128):
        sequence = "M" + chr(value) + "A"
        a = candidate("a", sequence, protein_junctions=[[1, 0], [3, 2], [4, 1]])
        b = candidate("b", sequence, protein_junctions=[[1, 1], [3, 2], [-1, 0]])
        try:
            expected = legacy_pair_score(a, b)
        except Exception as error:
            with pytest.raises(type(error)) as actual:
                pair_score(a, b)
            assert str(actual.value) == str(error)
        else:
            assert pair_json(pair_score(a, b)) == pair_json(expected)


def test_pair_cache_reuses_one_lazy_aligner_and_preserves_canonical_direction(monkeypatch):
    created = []
    original = selection._aligner

    def counted():
        value = original()
        created.append(value)
        return value

    monkeypatch.setattr(selection, "_aligner", counted)
    bounded = selection._PairScores(1, False, 2)
    bounded_inputs = [(candidate("one", "M"), candidate("two", "M")),
                      (candidate("one", "MM"), candidate("two", "MM"))]
    assert bounded.get(*bounded_inputs[0]).score == 1.0
    assert bounded.get(*bounded_inputs[1]).bounded
    assert not created and bounded.aligner is None
    cache = selection._PairScores(25_000_000, False, 2)
    # The cache belongs to one frozen selection invocation, which retains its
    # catalog candidates throughout the call. Keep that same ownership here.
    inputs = []
    for index in range(25):
        a = candidate("a", protein(index, 30), protein_junctions=[[10, 1], [30, 2]])
        b = candidate("b", a["protein"][:8] + "W" * 10 + a["protein"][8:], protein_junctions=[[10, 2], [40, 2]])
        inputs.append((a, b))
    for a, b in inputs:
        reverse = cache.signature(a) > cache.signature(b)
        expected = legacy_pair_score(b, a) if reverse else legacy_pair_score(a, b)
        if reverse:
            expected = selection.PairScore(expected.score, expected.identity, expected.coverage_b, expected.coverage_a,
                                           expected.junction_similarity, expected.aligned_residues,
                                           expected.internal_loss_b, expected.internal_loss_a, expected.bounded)
        assert pair_json(cache.get(a, b)) == pair_json(expected)
    assert len(created) == 1 and cache.aligner is created[0]


def invalid_donor_fixture(quality, origin="original"):
    source = protein(71, 100)
    alternative = source[:30] + "W" * 40 + source[70:]
    catalogs = [catalog("A", ("g", [coding_candidate("A_source", source),
                                    coding_candidate("A_alt", alternative)])),
                catalog("B", ("g", [coding_candidate("B_source", source)])),
                catalog("C", ("g", [coding_candidate("C_source", source)]))]
    catalogs[0]["loci"][0]["source_baseline_candidate_id"] = "A_source"
    edges = [edge("A", "B"), edge("A", "C")]
    donors = [catalog(f"P{index}", ("g", [coding_candidate(f"P{index}_excluded", alternative,
                                                        quality=quality, origin=origin)])) for index in range(5)]
    donor_edges = [edge("A", donor["species"]) for donor in donors]
    return catalogs, edges, donors, donor_edges


def invalid_length_normalization_fixture(quality=None):
    sequence = protein(71, 100)
    catalogs = [catalog("A", ("g", [coding_candidate("A_source", sequence + protein(72, 12)),
                                    coding_candidate("A_core", sequence)])),
                catalog("B", ("g", [coding_candidate("B_source", sequence)])),
                catalog("C", ("g", [coding_candidate("C_source", sequence)]))]
    catalogs[0]["loci"][0]["source_baseline_candidate_id"] = "A_source"
    if quality is not None:
        catalogs[0]["loci"][0]["candidates"].append(coding_candidate("A_invalid", sequence * 10, quality=quality))
    return catalogs, [edge("A", "B"), edge("A", "C")]


@pytest.mark.parametrize("quality", [
    {"usable": False, "annotated_pseudogene": True}, {"sequence_mismatch": True},
    {"phase_conflict": True}, {"translation_uncertain": True},
])
def test_ineligible_long_original_cannot_inflate_conserved_margin_or_change_quality_scale(quality):
    catalogs, edges = invalid_length_normalization_fixture()
    expected = select_representatives(catalogs, edges, exact_limit=2)
    assert expected["selections"][0]["candidate_id"] == "A_source"
    assert 0 < expected["selections"][0]["proposed_margin"] < .10
    modified, edges = invalid_length_normalization_fixture(quality)
    for cached, reordered in [(True, False), (False, False), (True, True)]:
        inputs = copy.deepcopy(modified)
        if reordered:
            inputs.reverse()
            for item in inputs:
                item["loci"][0]["candidates"].reverse()
        actual = select_representatives(inputs, list(reversed(edges)) if reordered else edges,
                                        cache_pair_scores=cached, exact_limit=2)
        assert actual["selections"] == expected["selections"]
        assert actual["audit"] == expected["audit"]
        assert [row for row in actual["scores"] if row["candidate_id"] != "A_invalid"] == expected["scores"]


def weighted_irregular_fixture():
    randomizer = random.Random(0)
    sequence = protein(71, 40)
    species = list("ABCDEFGH")
    catalogs = [catalog(name, ("g", [coding_candidate(name + "_source", sequence[:20] + protein(index + 90, 20),
                                                      quality={"rna_supported": index == 1, "partial": index == 2})]))
                for index, name in enumerate(species)]
    catalogs[0]["loci"][0]["candidates"].append(coding_candidate("A_other", sequence[:10] + "W" * 20 + sequence[30:]))
    edges = [edge(a, b, weight=randomizer.uniform(.01, 20))
             for index, a in enumerate(species) for b in species[index + 1:]]
    randomizer.shuffle(edges)
    return catalogs, edges


def test_weighted_irregular_graph_is_exactly_invariant_to_edge_and_catalog_order():
    catalogs, edges = weighted_irregular_fixture()
    result = select_representatives(catalogs, edges, exact_limit=2)
    reversed_result = select_representatives(list(reversed(catalogs)), list(reversed(edges)), exact_limit=2)
    oracle = select_representatives(catalogs, edges, cache_pair_scores=False, exact_limit=2)
    for key in ["selections", "scores", "audit", "omitted_edges"]:
        assert result[key] == reversed_result[key] == oracle[key]


@pytest.mark.parametrize("cached", [True, False])
def test_selector_irregular_graph_matches_legacy_pair_kernel_exact_json(monkeypatch, cached):
    catalogs, edges = weighted_irregular_fixture()
    optimized = select_representatives(catalogs, edges, exact_limit=2, cache_pair_scores=cached)

    def native_oracle(signature_a, signature_b, *, max_alignment_cells, aligner_factory):
        a, b = ({"protein": signature[0] + "*", "protein_junctions": signature[1]}
                for signature in (signature_a, signature_b))
        return legacy_pair_score(a, b, max_alignment_cells=max_alignment_cells)

    monkeypatch.setattr(selection, "_score_signatures", native_oracle)
    oracle = select_representatives(list(reversed(catalogs)), list(reversed(edges)),
                                   exact_limit=2, cache_pair_scores=cached)
    assert json.dumps(optimized, sort_keys=True) == json.dumps(oracle, sort_keys=True)


@pytest.mark.parametrize("quality,origin", [
    ({"usable": False, "valid_orf": False, "annotated_pseudogene": True}, "original"),
    ({"sequence_mismatch": True}, "original"),
    ({"phase_conflict": True}, "original"),
    ({"translation_uncertain": True, "rna_supported": True}, "original"),
    ({"representative_eligible": False}, "homology_predicted"),
])
def test_excluded_donors_cannot_change_selection_margin_or_global_objective(quality, origin):
    catalogs, edges, donors, donor_edges = invalid_donor_fixture(quality, origin)
    expected = select_representatives(catalogs, edges, exact_limit=10)
    result = select_representatives(catalogs + donors, edges + donor_edges, exact_limit=10)
    assert expected["selections"][0]["candidate_id"] == "A_source"
    assert result["selections"][:3] == expected["selections"]
    assert result["scores"] == expected["scores"]
    assert result["audit"] == expected["audit"]
    assert result["metrics"]["trusted_edges"] == expected["metrics"]["trusted_edges"] == 2
    assert result["metrics"]["pair_requests"] == expected["metrics"]["pair_requests"]
    assert all(row["reason"] == "no_eligible_candidate" for row in result["selections"][3:])
    oracle = select_representatives(catalogs + donors, edges + donor_edges, cache_pair_scores=False, exact_limit=10)
    assert oracle["selections"] == result["selections"]
    assert oracle["scores"] == result["scores"]
    assert oracle["audit"] == result["audit"]


@pytest.mark.parametrize("cached", [True, False])
def test_eligibility_is_compiled_once_per_invocation_and_rebuilt_for_changed_candidates(monkeypatch, cached):
    catalogs, edges, donors, donor_edges = invalid_donor_fixture({"translation_uncertain": True})
    inputs, graph = catalogs + donors, edges + donor_edges
    original = selection._eligible
    calls = {}

    def counted(candidate):
        identifier = id(candidate)
        calls[identifier] = calls.get(identifier, 0) + 1
        return original(candidate)

    monkeypatch.setattr(selection, "_eligible", counted)
    before = select_representatives(inputs, graph, cache_pair_scores=cached)
    records = [candidate for item in inputs for locus in item["loci"] for candidate in locus["candidates"]]
    assert calls == {id(candidate): 1 for candidate in records}
    assert before["selections"][0]["candidate_id"] == "A_source"
    for donor in donors:
        donor["loci"][0]["candidates"][0]["quality"]["translation_uncertain"] = False
    calls.clear()
    after = select_representatives(inputs, graph, cache_pair_scores=cached)
    assert calls == {id(candidate): 1 for candidate in records}
    assert after["selections"][0]["candidate_id"] == "A_alt"


@pytest.mark.parametrize("cached", [True, False])
def test_disconnected_component_trials_match_separate_calls_with_exact_objectives(cached):
    fixtures = [extension_fixture(), invalid_length_normalization_fixture({"sequence_mismatch": True})]
    catalogs, edges, parts = [], [], []
    for index, (items, graph) in enumerate(fixtures):
        prefix = "family" + str(index) + "_"
        for item in items:
            for locus in item["loci"]:
                locus["gene_id"] = prefix + locus["gene_id"]
                for candidate in locus["candidates"]:
                    candidate["candidate_id"] = prefix + candidate["candidate_id"]
                    candidate["source_transcript_id"] = prefix + candidate["source_transcript_id"]
                for field in ("baseline_candidate_id", "source_baseline_candidate_id"):
                    if locus.get(field):
                        locus[field] = prefix + locus[field]
        for edge in graph:
            edge["gene_a"], edge["gene_b"] = prefix + edge["gene_a"], prefix + edge["gene_b"]
        parts.append(select_representatives(items, graph, cache_pair_scores=cached, exact_limit=2))
        catalogs.extend(items)
        edges.extend(graph)
    assert parts[0]["selections"][0]["candidate_id"] == "family0_A_conserved"
    assert parts[1]["selections"][0]["candidate_id"] == "family1_A_source"
    combined = select_representatives(list(reversed(catalogs)), list(reversed(edges)),
                                      cache_pair_scores=cached, exact_limit=2)
    keys = {"selections": lambda row: (row["species"], row["gene_id"]),
            "scores": lambda row: (row["species"], row["gene_id"], row["candidate_id"]),
            "audit": lambda row: row["loci"], "omitted_edges": lambda row: (row["a"], row["b"])}
    for field, key in keys.items():
        assert sorted(combined[field], key=key) == sorted([row for part in parts for row in part[field]], key=key)


def test_donor_returned_to_ineligible_source_is_removed_from_final_degree_normalization():
    catalogs, edges, donors, donor_edges = invalid_donor_fixture({"sequence_mismatch": True})
    expected = select_representatives(catalogs, edges)
    donor = donors[0]
    unusable = donor["loci"][0]["candidates"][0]
    donor["loci"][0]["source_baseline_candidate_id"] = unusable["candidate_id"]
    donor["loci"][0]["candidates"].append(coding_candidate("P0_proposal", unusable["protein"]))
    result = select_representatives(catalogs + [donor], edges + donor_edges[:1])
    assert result["selections"][0] == expected["selections"][0]
    assert result["scores"][:len(expected["scores"])] == expected["scores"]
    returned = result["selections"][-1]
    assert returned["candidate_id"] == unusable["candidate_id"]
    assert returned["proposed_candidate_id"] == "P0_proposal"
    source_score = next(row for row in result["scores"] if row["candidate_id"] == unusable["candidate_id"])
    assert returned["score"] == pytest.approx(source_score["quality_score"])


def test_excluded_donor_cannot_supply_the_second_independent_species_for_adoption():
    catalogs, edges = extension_fixture()
    catalogs[2]["loci"][0]["candidates"][0]["quality"]["sequence_mismatch"] = True
    result = select_representatives(catalogs, edges)
    assert result["selections"][0]["candidate_id"] == "A_long"
    assert result["selections"][0]["reason"] == "insufficient_independent_donor_species"


def test_conserved_choice_changes_longest_wrong_terminal_extension():
    catalogs, edges = extension_fixture()
    result = select_representatives(catalogs, edges)
    selected = result["selections"][0]
    assert selected["candidate_id"] == "A_conserved"
    assert selected["status"] == "conserved"
    assert selected["margin"] >= 0.10
    assert result["metrics"]["changed_representatives"] == 1
    assert result["parameters"]["cutoffs_calibrated"] is False


def test_longest_policy_uses_corrected_unpadded_length_and_stable_tie():
    catalogs = [catalog("A", ("g", [candidate("z", protein(), corrected_cds_length=268),
                                   candidate("a", protein(), corrected_cds_length=269)]))]
    assert select_representatives(catalogs, [], "longest")["selections"][0]["candidate_id"] == "a"
    catalogs[0]["loci"][0]["candidates"][0]["corrected_cds_length"] = 269
    assert select_representatives(catalogs, [], "longest")["selections"][0]["candidate_id"] == "a"


@pytest.mark.parametrize("quality,expected", [
    ({"representative_eligible": False}, "A_source"),
    ({"representative_eligible": True, "rna_supported": True}, "A_prediction"),
    ({"representative_eligible": True}, "A_prediction"),
    ({"representative_eligible": True, "usable": False}, "A_source"),
    ({"representative_eligible": True, "sequence_mismatch": True}, "A_source"),
])
def test_longest_predictions_require_usable_and_representative_permission(quality, expected):
    sequence = protein(71, 90)
    source = coding_candidate("A_source", sequence)
    predicted = coding_candidate("A_prediction", sequence + protein(72, 50), quality=quality, origin="predicted")
    catalogs = [catalog("A", ("g", [source, predicted]))]
    assert select_representatives(catalogs, [], "longest")["selections"][0]["candidate_id"] == expected
    # Quality exclusions on originals do not change the legacy longest rule.
    predicted["origin"] = "original"
    assert select_representatives(catalogs, [], "longest")["selections"][0]["candidate_id"] == "A_prediction"


@pytest.mark.parametrize("origin", ["original", "predicted"])
def test_single_eligible_repair_compares_with_ineligible_source_and_requires_donors(origin):
    sequence = protein(73, 100)
    source = coding_candidate("A_source", sequence[:35], quality={"usable": False, "valid_orf": False,
                                                                 "frameshift": True})
    repair = coding_candidate("A_repair", sequence, quality={"representative_eligible": True}, origin=origin)
    catalogs = [catalog("A", ("g", [source, repair])),
                catalog("B", ("g", [coding_candidate("B_source", sequence)])),
                catalog("C", ("g", [coding_candidate("C_source", sequence)]))]
    catalogs[0]["loci"][0]["source_baseline_candidate_id"] = "A_source"
    edges = [edge("A", "B"), edge("A", "C")]
    result = select_representatives(catalogs, edges, exact_limit=1)
    chosen = result["selections"][0]
    assert chosen["candidate_id"] == "A_repair"
    assert chosen["status"] == "conserved"
    source_row = next(row for row in result["scores"] if row["candidate_id"] == "A_source")
    assert chosen["margin"] == pytest.approx(chosen["score"] - source_row["quality_score"])
    assert chosen["margin"] > result["parameters"]["min_margin"]
    assert result["audit"][0]["heuristic_matches_exact"] is True
    catalogs[2]["loci"][0]["candidates"][0]["quality"]["sequence_mismatch"] = True
    unsupported = select_representatives(catalogs, edges)["selections"][0]
    assert unsupported["candidate_id"] == "A_source"
    assert unsupported["reason"] == "insufficient_independent_donor_species"
    if origin == "predicted":
        repair["quality"]["representative_eligible"] = False
        unaccepted = select_representatives(catalogs, edges)["selections"][0]
        assert unaccepted["candidate_id"] == "A_source"
        assert unaccepted["reason"] == "no_eligible_candidate"


def test_single_eligible_repair_needs_measured_baseline_improvement():
    sequence = protein(74, 100)
    source = coding_candidate("A_source", sequence, quality={"sequence_mismatch": True,
                                                             "independent_support": 1})
    repair = coding_candidate("A_repair", sequence[:30], quality={"representative_eligible": True,
                                                                "partial": True}, origin="predicted")
    catalogs = [catalog("A", ("g", [source, repair])),
                catalog("B", ("g", [coding_candidate("B_source", sequence)])),
                catalog("C", ("g", [coding_candidate("C_source", sequence)]))]
    catalogs[0]["loci"][0]["source_baseline_candidate_id"] = "A_source"
    result = select_representatives(catalogs, [edge("A", "B"), edge("A", "C")])
    assert result["selections"][0]["candidate_id"] == "A_source"
    assert result["selections"][0]["reason"] == "candidate_score_margin"
    assert result["selections"][0]["proposed_margin"] < 0


def test_species_specific_internal_region_is_preserved_without_target_evidence():
    sequence = protein()
    species_specific = sequence[:40] + protein(6, 20) + sequence[40:]
    catalogs = [catalog("A", ("g", [candidate("A_long", species_specific), candidate("A_shared", sequence)])),
                catalog("B", ("g", [candidate("B", sequence)])),
                catalog("C", ("g", [candidate("C", sequence)]))]
    result = select_representatives(catalogs, [edge("A", "B"), edge("A", "C")])
    record = result["selections"][0]
    assert record["candidate_id"] == "A_long"
    assert record["status"] == "insufficient_evidence"
    assert record["reason"] == "possible_species_specific_internal_region"
    catalogs[0]["loci"][0]["candidates"][1]["quality"]["full_length_supported"] = True
    result = select_representatives(catalogs, [edge("A", "B"), edge("A", "C")])
    assert result["selections"][0]["candidate_id"] == "A_shared"


def test_common_short_fragment_cannot_replace_major_coding_span():
    full = protein(length=150)
    fragment = full[:60]
    catalogs = [catalog("A", ("g", [candidate("A_full", full), candidate("A_short", fragment)])),
                catalog("B", ("g", [candidate("B", fragment)])),
                catalog("C", ("g", [candidate("C", fragment)]))]
    result = select_representatives(catalogs, [edge("A", "B"), edge("A", "C")])
    assert result["selections"][0]["candidate_id"] == "A_full"
    assert result["selections"][0]["reason"] == "shorter_than_source_major_coding_span"


def test_wgd_copies_are_never_collapsed_and_unresolved_copies_abstain():
    seq1, seq2 = protein(1), protein(2)
    catalogs = [catalog(species, ("copy1", [candidate(f"{species}_1", seq1)]),
                        ("copy2", [candidate(f"{species}_2", seq2)])) for species in ["A", "B", "C"]]
    edges = [edge(a, b, gene_a=gene, gene_b=gene) for a, b in [("A", "B"), ("A", "C")]
             for gene in ["copy1", "copy2"]]
    result = select_representatives(catalogs, edges)
    assert len(result["selections"]) == 6
    assert len({(record["species"], record["gene_id"]) for record in result["selections"]}) == 6
    result = select_representatives(catalogs, edges + [edge("A", "B", "copy1", "copy2")])
    assert result["selections"][0]["status"] == "ambiguous_correspondence"
    assert len(result["selections"]) == 6


def test_ties_keep_source_baseline_and_low_margin_abstains():
    catalogs, edges = extension_fixture()
    records = catalogs[0]["loci"][0]["candidates"]
    records[0] = candidate("A_long", records[1]["protein"])
    catalogs[0]["loci"][0]["source_baseline_candidate_id"] = "A_long"
    result = select_representatives(catalogs, edges)
    assert result["selections"][0]["candidate_id"] == "A_long"
    assert result["selections"][0]["margin"] == pytest.approx(0)
    catalogs, edges = extension_fixture()
    result = select_representatives(catalogs, edges, min_margin=10)
    selected = result["selections"][0]
    assert selected["candidate_id"] == "A_long"
    assert selected["status"] == "insufficient_evidence"
    assert selected["reason"] == "candidate_score_margin"
    actual = next(score for score in result["scores"] if score["candidate_id"] == "A_long")
    assert selected["score"] == pytest.approx(actual["score"])
    assert selected["proposed_candidate_id"] == "A_conserved"


def test_common_reference_star_does_not_inflate_conservation_score_or_margin():
    sequence = protein()
    selections = []
    for donor_count in [2, 24]:
        catalogs = [catalog("A", ("g", [candidate("A_long", sequence + protein(5, 12)),
                                      candidate("A_core", sequence)]))]
        donors = [f"D{index:03d}" for index in range(donor_count)]
        catalogs += [catalog(species, ("g", [candidate(species + "_core", sequence)])) for species in donors]
        result = select_representatives(catalogs, [edge("A", species) for species in donors], exact_limit=2)
        selected = result["selections"][0]
        assert selected["candidate_id"] == "A_long"
        assert selected["reason"] == "candidate_score_margin"
        assert 0 < selected["proposed_margin"] < result["parameters"]["min_margin"]
        assert all(row["score"] <= row["quality_score"] + 1 + 1e-12 for row in result["scores"])
        assert result["audit"][0]["heuristic_matches_exact"] is True
        selections.append(selected)
    for field in ["score", "margin", "proposed_score", "proposed_margin"]:
        assert selections[0][field] == pytest.approx(selections[1][field])


def test_weighted_edges_change_mean_conservation_and_symmetric_objective_consistently():
    sequence = protein()
    alternative = sequence[:35] + "A" * 20 + sequence[55:]
    a = candidate("A_a", sequence, quality={"full_length_supported": True})
    b = candidate("A_b", alternative, quality={"full_length_supported": True})
    catalogs = [catalog("A", ("g", [a, b])),
                catalog("B", ("g", [candidate("B_a", sequence)])),
                catalog("C", ("g", [candidate("C_b", alternative)]))]
    catalogs[0]["loci"][0]["source_baseline_candidate_id"] = "A_a"
    overlap = pair_score(a, b).score
    for weight_b, weight_c, expected in [(9, 1, "A_a"), (1, 9, "A_b")]:
        result = select_representatives(catalogs, [edge("A", "B", weight=weight_b),
                                                  edge("A", "C", weight=weight_c)],
                                        min_support=1, min_margin=0, exact_limit=2)
        selected = result["selections"][0]
        assert selected["candidate_id"] == expected
        # Symmetric coefficients are .95 and .55, giving I_A=1.5.
        assert selected["score"] == pytest.approx(.75 + (.95 + .55 * overlap) / 1.5)
        assert selected["margin"] == pytest.approx(.40 * (1 - overlap) / 1.5)
        assert result["audit"][0]["objective"] == pytest.approx(1.875 + .95 + .55 * overlap)
        assert result["audit"][0]["heuristic_matches_exact"] is True
        scaled = select_representatives(catalogs, [edge("A", "B", weight=10 * weight_b),
                                                  edge("A", "C", weight=10 * weight_c)],
                                        min_support=1, min_margin=0, exact_limit=2)
        for scaled_row, row in zip(scaled["selections"], result["selections"], strict=True):
            assert {key: value for key, value in scaled_row.items() if key not in {"score", "margin"}} == {
                key: value for key, value in row.items() if key not in {"score", "margin"}}
            assert scaled_row["score"] == pytest.approx(row["score"])
            assert scaled_row["margin"] == pytest.approx(row["margin"])


def test_projected_junctions_use_alignment_positions_and_phase_not_genomic_coordinates():
    sequence = protein(length=30)
    a = candidate("A", sequence, blocks=[[100, 130, 0], [500, 560, 0]])
    b = candidate("B", sequence, blocks=[[300, 330, 0], [900, 960, 0]])
    assert pair_score(a, b).junction_similarity == 1
    b["blocks"] = [[300, 331, 0], [900, 959, 2]]
    assert pair_score(a, b).junction_similarity == 0
    b["strand"] = "-"
    b["blocks"] = [[300, 360, 0], [900, 930, 0]]
    assert pair_score(a, b).junction_similarity == 1


def test_invalid_models_and_homology_only_predictions_cannot_win():
    catalogs, edges = extension_fixture()
    proposed = catalogs[0]["loci"][0]["candidates"][1]
    proposed["quality"]["phase_conflict"] = True
    result = select_representatives(catalogs, edges)
    assert result["selections"][0]["candidate_id"] == "A_long"
    proposed["quality"].pop("phase_conflict")
    proposed["origin"] = "homology_predicted"
    result = select_representatives(catalogs, edges)
    assert result["selections"][0]["candidate_id"] == "A_long"
    proposed["quality"]["representative_eligible"] = True
    result = select_representatives(catalogs, edges)
    assert result["selections"][0]["candidate_id"] == "A_conserved"


def test_sparse_graph_cache_is_equivalent_to_uncached_oracle_and_deterministic():
    catalogs, edges = extension_fixture()
    cached = select_representatives(catalogs, edges, exact_limit=100)
    uncached = select_representatives(catalogs, edges, cache_pair_scores=False, exact_limit=100)
    assert cached["selections"] == uncached["selections"]
    assert cached["scores"] == uncached["scores"]
    assert cached["audit"] == uncached["audit"]
    assert cached["metrics"]["pair_evaluations"] < uncached["metrics"]["pair_evaluations"]
    assert all(component["heuristic_matches_exact"] for component in cached["audit"])
    shuffled = copy.deepcopy(catalogs)
    shuffled.reverse()
    for item in shuffled:
        item["loci"].reverse()
        for locus in item["loci"]:
            locus["candidates"].reverse()
    reversed_result = select_representatives(shuffled, list(reversed(edges)), exact_limit=100)
    assert cached == reversed_result
    evicting = select_representatives(catalogs, edges, pair_cache_limit=1, exact_limit=100)
    assert evicting["selections"] == cached["selections"]
    assert evicting["scores"] == cached["scores"]
    assert evicting["audit"] == cached["audit"]
    assert evicting["metrics"]["pair_cache_peak_entries"] == 1
    assert evicting["metrics"]["pair_cache_evictions"] > 0


def test_joint_small_component_reports_comparison_to_exact_state_space():
    sequence = protein(length=60)
    catalogs = [catalog(species, ("g", [candidate(f"{species}_long", sequence + protein(seed, 20)),
                                      candidate(f"{species}_core", sequence),
                                      candidate(f"{species}_short", sequence[:30])]))
                for species, seed in [("A", 5), ("B", 7), ("C", 9)]]
    result = select_representatives(catalogs, [edge("A", "B"), edge("B", "C"), edge("A", "C")],
                                    exact_limit=100)
    assert len(result["audit"]) == 1
    assert result["audit"][0]["state_count"] == 27
    assert result["audit"][0]["heuristic_matches_exact"] is True
    assert result["audit"][0]["objective"] == pytest.approx(result["audit"][0]["exact_objective"])


def test_ambiguous_edges_do_not_supply_votes_and_duplicate_edges_do_not_inflate_support():
    catalogs, edges = extension_fixture()
    result = select_representatives(catalogs, [edges[0], edges[0]])
    assert result["selections"][0]["candidate_id"] == "A_long"
    assert result["selections"][0]["reason"] == "insufficient_independent_donor_species"
    result = select_representatives(catalogs, [edges[0] | {"ambiguous": True}, edges[1]])
    assert result["selections"][0]["status"] == "ambiguous_correspondence"


def test_ambiguous_locus_coordinates_abstain_and_cannot_vote_as_resolved_donor():
    catalogs, edges = extension_fixture()
    catalogs[0]["loci"][0]["ambiguous_coordinates"] = True
    result = select_representatives(catalogs, edges)
    assert result["selections"][0]["candidate_id"] == "A_long"
    assert result["selections"][0]["reason"] == "ambiguous_locus_coordinates"
    assert result["metrics"]["trusted_edges"] == 0


def test_large_candidate_sets_and_alignments_abstain_with_recorded_limits():
    catalogs, edges = extension_fixture()
    result = select_representatives(catalogs, edges, candidate_limit=1)
    assert result["selections"][0]["reason"] == "candidate_limit_exceeded"
    result = select_representatives(catalogs, edges, max_alignment_cells=1)
    assert result["selections"][0]["candidate_id"] == "A_long"
    assert result["metrics"]["bounded_alignments"] > 0
    assert pair_score(candidate("a", "M"), candidate("b", "M"), max_alignment_cells=0).bounded


@pytest.mark.parametrize("mutate", [
    lambda catalogs: catalogs.append(copy.deepcopy(catalogs[0])),
    lambda catalogs: catalogs[0]["loci"][0]["candidates"].append(copy.deepcopy(
        catalogs[0]["loci"][0]["candidates"][0])),
    lambda catalogs: catalogs[0]["loci"][0].update(baseline_candidate_id="absent"),
    lambda catalogs: catalogs[0].update(schema=999),
    lambda catalogs: catalogs[0]["loci"].append(
        dict(copy.deepcopy(catalogs[0]["loci"][0]), gene_id="another_owner")),
])
def test_invalid_catalog_contract_fails(mutate):
    catalogs, edges = extension_fixture()
    mutate(catalogs)
    with pytest.raises(ValueError):
        select_representatives(catalogs, edges)


def test_empty_catalog_is_well_defined_and_invalid_edges_fail():
    assert select_representatives([], [])["selections"] == []
    catalogs, _ = extension_fixture()
    with pytest.raises(ValueError, match="absent locus"):
        select_representatives(catalogs, [edge("A", "D")])
    with pytest.raises(ValueError, match="different species"):
        select_representatives(catalogs, [edge("A", "A")])
    with pytest.raises(ValueError, match="positive"):
        select_representatives(catalogs, [edge("A", "B", weight=0)])

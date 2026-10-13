"""Unanchored annotation does not imply orthology or an expected WGD copy."""
import copy
import hashlib
import json
import os
import shlex
import subprocess
from pathlib import Path

import pytest

from workflow.support.rescue_additional_candidates import (
    nominate_genome_only_candidates,
    reassess_unanchored_models,
)
from workflow.support.rescue_prediction_cache import (
    LEGACY_DIRECT_IMPLEMENTATION,
    LEGACY_INHERITED_IMPLEMENTATION,
    LEGACY_INHERITED_VERIFIER,
    LEGACY_MINIPROT_SHA256,
    SEARCH_CONTRACT_FILE,
    frozen_prediction_cache_key,
    prediction_search_contract,
    stream_json_array,
    verify_prediction_cache,
)

PARAMS = {"minimum_coverage": .8, "minimum_identity": .5, "cscore": .7, "max_intron": 20000,
          "genome_fallback": True, "max_interval": 200000, "padding": 1000,
          "min_anchors": 3, "distance": 20, "diagonal_bound": 10}


def proteins(root, species, ids):
    directory = root / "prepared" / species
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "genes.pep").write_text("".join(f">{key}\n{'M' * length}\n" for key, length in ids.items()))


def comparison(root, rows, reverse=False):
    proteins(root, "T", {"t1": 100, "t2": 100})
    proteins(root, "D", {"absent": 100, "partial": 100, "full": 100, "ambiguous": 100, "anchored": 100})
    directory = root / "synteny" / "pair1"
    directory.mkdir(parents=True)
    (directory / "target.query.last").write_text("".join("\t".join(map(str, row)) + "\n" for row in rows))
    return {"donors": {"T": ["D"]}, "synteny_jobs": [{"id": "pair1", "a": "D" if reverse else "T", "b": "T" if reverse else "D"}]}


def hit(target, donor, length=100, identity=90, score=100):
    return [target, donor, identity, length, 0, 0, 1, length, 1, length, 1e-10, score]


def test_nomination_preserves_partial_and_copy_ambiguous_queries(tmp_path):
    plan = comparison(tmp_path, [hit("t1", "partial", 60), hit("t1", "full"),
                                 hit("t1", "ambiguous"), hit("t2", "ambiguous", score=90)])
    regions = [{"donor": "D", "query": "anchored", "id": "r1"}]
    rows, diagnostic = nominate_genome_only_candidates(tmp_path, plan, "T", regions, PARAMS)
    assert [row["query"] for row in rows] == ["partial", "absent", "ambiguous"]
    assert all(row["genome_only"] and "seqid" not in row and row["orthology"] == "unassigned" for row in rows)
    assert diagnostic["counts"]["represented"] == 1
    assert diagnostic["counts"]["already_nominated"] == 1
    assert regions == [{"donor": "D", "query": "anchored", "id": "r1"}]


def test_nomination_reverse_uses_donor_query_coordinates(tmp_path):
    row = hit("full", "t1")
    row[8:10] = [1, 20]  # target partial does not erase full donor coverage
    plan = comparison(tmp_path, [row], reverse=True)
    rows, diagnostic = nominate_genome_only_candidates(tmp_path, plan, "T", [], PARAMS)
    assert "full" not in {r["query"] for r in rows}
    assert diagnostic["counts"]["represented"] == 1


def test_nomination_threshold_uses_single_target_match_not_union(tmp_path):
    a, b = hit("t1", "partial", 60), hit("t2", "partial", 60)
    b[8:10] = [41, 100]
    plan = comparison(tmp_path, [a, b, hit("t1", "full", identity=20)])
    rows, _ = nominate_genome_only_candidates(tmp_path, plan, "T", [], PARAMS)
    assert {"partial", "full"} <= {row["query"] for row in rows}


def test_nomination_cap_is_deterministic_and_reported(tmp_path):
    plan = comparison(tmp_path, [])
    rows, diagnostic = nominate_genome_only_candidates(tmp_path, plan, "T", [], {**PARAMS, "max_genome_queries": 2})
    assert len(rows) == 2 and diagnostic["deferred_by_limit"] == 3
    assert len(diagnostic["deferred_queries"]) == 3


def test_many_to_one_donor_paralogs_are_not_claimed_represented(tmp_path):
    plan = comparison(tmp_path, [hit("t1", "full"), hit("t1", "partial")])
    rows, diagnostic = nominate_genome_only_candidates(tmp_path, plan, "T", [], PARAMS)
    assert diagnostic["counts"]["copy_ambiguous"] == 2
    assert all(row["nomination"].get("ambiguity") == "multiple_donor_genes_share_one_full_target"
               for row in rows if row["query"] in {"full", "partial"})


def test_cap_balances_donors_and_evidence_buckets():
    from workflow.support.rescue_additional_candidates import _balanced_queue
    rows = [{"donor": donor, "query": f"{reason}{index:03}", "nomination": {
        "reason": reason, "best_target_identity": .9, "best_target_coverage": .7}}
        for donor in ("A", "B", "C", "D") for reason in ("partial_target_match", "no_target_match", "copy_ambiguous")
        for index in range(12)]
    selected = _balanced_queue(rows, {"C", "D"}, .5)[:36]
    counts = {donor: sum(r["donor"] == donor for r in selected) for donor in "ABCD"}
    assert counts == {"A": 6, "B": 6, "C": 12, "D": 12}
    assert {r["screen_priority"] for r in selected} == {"strong_partial", "no_target_match", "copy_ambiguous"}


def test_nomination_unknown_alignment_fails(tmp_path):
    plan = comparison(tmp_path, [hit("unknown", "full")])
    with pytest.raises(ValueError, match="unknown prepared"):
        nominate_genome_only_candidates(tmp_path, plan, "T", [], PARAMS)


def test_nomination_disabled_fallback_does_not_search(tmp_path):
    rows, diagnostic = nominate_genome_only_candidates(tmp_path, {}, "T", [], {"genome_fallback": False})
    assert not rows and diagnostic["policy"] == "disabled_genome_fallback"


def model(donor="D", gene="q", start=0, seqid="chr1", problems=None):
    return {"seqid": seqid, "strand": "+", "cds": [[start, start + 99, 0]],
            "sequence": "ATG" + "AAA" * 31 + "TAA", "coverage": 1, "identity": .9,
            "problems": ["outside_expected_synteny_interval"] if problems is None else problems,
            "evidence": {"target": "T", "donor": donor, "query": gene}}


def test_off_block_requires_two_independent_species():
    rows = [model(), model("E")]
    diagnostic = reassess_unanchored_models(rows)
    assert rows[0]["problems"] == rows[1]["problems"] == []
    assert rows[0]["placement_evidence"]["orthology"] == "unassigned"
    assert diagnostic["counts"]["supported_unanchored_annotation"] == 2


def test_paralogs_and_self_do_not_inflate_support():
    rows = [model("D", "q1"), model("D", "q2"), model("T", "self")]
    reassess_unanchored_models(rows)
    assert all("outside_expected_synteny_interval" in row["problems"] for row in rows)
    assert rows[0]["placement_evidence"]["independent_donor_species"] == ["D"]


def test_multiple_genomic_copies_remain_proposals():
    rows = [model(), model("E"), model(start=200), model("E", start=200)]
    reassess_unanchored_models(rows)
    assert all(row["problems"] for row in rows)
    assert not rows[0]["placement_evidence"]["independent_donor_species"]


def test_hard_disruption_does_not_supply_support():
    rows = [model(), model("E", problems=["frameshift", "outside_expected_synteny_interval"])]
    reassess_unanchored_models(rows)
    assert all(row["problems"] for row in rows)
    assert "placement_evidence" not in rows[1]


def test_unchecked_coordinates_cannot_supply_unanchored_support():
    rows = [model(), model("E")]
    for row in rows:
        del row["identity"]
    reassess_unanchored_models(rows)
    assert all(row["problems"] for row in rows)


def test_compatible_isoforms_can_corroborate_annotation():
    rows = [model(), model("E", start=3)]
    reassess_unanchored_models(rows)
    assert not rows[0]["problems"] and not rows[1]["problems"]


def test_incompatible_paths_do_not_bridge_annotations():
    rows = [model(), model("E", start=1)]
    reassess_unanchored_models(rows)
    assert all(row["problems"] for row in rows)
    assert rows[0]["placement_evidence"]["reason"] == "incompatible_unanchored_paths"


@pytest.mark.parametrize("text", ['[{}] trailing', '[{},]', '[{}', '{"not":"array"}', '[NaN]', '[{"key":1,"key":2}]'])
def test_stream_reader_rejects_malformed_or_nonfinite_records(tmp_path, text):
    path = tmp_path / "array.json"
    path.write_text(text)
    with pytest.raises(ValueError):
        list(stream_json_array(path, chunk_size=2))


def test_stream_reader_boundaries_and_content_hash(tmp_path):
    path = tmp_path / "array.json"
    content = json.dumps([{"name": "µRNA", "values": [1, 2]}, {}, "a"])
    path.write_text(content)
    hasher = hashlib.sha256()
    assert list(stream_json_array(path, chunk_size=2, hasher=hasher)) == json.loads(content)
    assert hasher.hexdigest() == hashlib.sha256(content.encode()).hexdigest()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_receipt(directory, key):
    files = {path.name: sha(path) for path in directory.iterdir() if path.is_file() and path.name != "receipt.json"}
    (directory / "receipt.json").write_text(json.dumps({"key": key, "files": files}))


def cache_fixture(tmp_path):
    old, new = tmp_path / "old", tmp_path / "new"
    old.mkdir()
    new.mkdir()
    files, sources = {}, {}
    for name in ("T", "D", "E"):
        source = {"genetic_code": 1}
        for field in ("genome", "fasta", "gff"):
            path = tmp_path / f"{name}.{field}"
            path.write_text(name + field)
            source[field] = str(path)
            files[str(path)] = sha(path)
        sources[name] = source
        proteins(old, name, {name + "_q": 33})
        proteins(new, name, {name + "_q": 33})
    plan = {"species": ["T"], "request": {"sources": sources, "files": files, "parameters": PARAMS,
                                              "tools": {"miniprot": "1", "miniprot_sha256": "a" * 64,
                                                        "prediction_search_contract": prediction_search_contract()}}}
    (old / "plan.json").write_text(json.dumps(plan))
    plan_hash = sha(old / "plan.json")
    prepared = {}
    for name in ("T", "D", "E"):
        write_receipt(old / "prepared" / name, {"plan": plan_hash, "species": name})
        prepared[name] = sha(old / "prepared" / name / "receipt.json")
    region = {"id": "region1", "target": "T", "donor": "D", "query": "D_q", "seqid": "chr1",
              "start": 0, "end": 99, "expected_start": 0, "expected_end": 99}
    prediction = {**model(), "evidence": region, "query": "region1", "search": "synteny_interval",
                  "frameshift": False, "status": "accepted", "model_id": "old_decision", "support": [region]}
    directory = old / "rescued" / "T"
    directory.mkdir(parents=True)
    (directory / "models.json").write_text(json.dumps([prediction]))
    (directory / "candidates.json").write_text(json.dumps([region]))
    (directory / "genome_query_mapping.tsv").write_text("candidate\trepresentative\nregion1\tregion1\n")
    (directory / SEARCH_CONTRACT_FILE).write_text(json.dumps(prediction_search_contract()))
    write_receipt(directory, {"plan": plan_hash, "species": "T", "prepared": prepared})
    return old, new, plan, region


def rebind_fixture_plan(root, update):
    """Controlled synthetic producer change; recompute all fixture receipts."""
    plan = json.loads((root / "plan.json").read_text())
    update(plan)
    (root / "plan.json").write_text(json.dumps(plan))
    plan_hash = sha(root / "plan.json")
    prepared = {}
    for directory in (root / "prepared").iterdir():
        write_receipt(directory, {"plan": plan_hash, "species": directory.name})
        prepared[directory.name] = sha(directory / "receipt.json")
    for directory in (root / "rescued").iterdir():
        write_receipt(directory, {"plan": plan_hash, "species": directory.name, "prepared": prepared})
    return plan


def legacy_fixture(root, implementation=LEGACY_DIRECT_IMPLEMENTATION, parent=None):
    (root / "rescued/T" / SEARCH_CONTRACT_FILE).unlink(missing_ok=True)
    def update(plan):
        tools = plan["request"]["tools"]
        tools.pop("prediction_search_contract", None)
        tools.update(implementation=implementation, miniprot="0.18-r281", miniprot_sha256=LEGACY_MINIPROT_SHA256)
        if parent:
            plan["request"]["prediction_cache"] = parent
            tools["verify_prediction_cache_implementation"] = LEGACY_INHERITED_VERIFIER
    return rebind_fixture_plan(root, update)


def current_fixture_plan(plan):
    result = copy.deepcopy(plan)
    result["request"]["tools"]["prediction_search_contract"] = prediction_search_contract()
    return result


def test_search_contract_is_frozen_with_modern_predictions(tmp_path):
    old, new, plan, region = cache_fixture(tmp_path)
    frozen = frozen_prediction_cache_key(old)
    assert frozen["species"]["T"]["files"][SEARCH_CONTRACT_FILE] == sha(old / "rescued/T" / SEARCH_CONTRACT_FILE)
    cache = verify_prediction_cache(old, new, plan, "T", [region], PARAMS, frozen)
    assert cache.local_search_compatible and cache.genome_search_compatible
    assert cache.search_contract == prediction_search_contract()


@pytest.mark.parametrize("legacy", [False, True])
def test_narrow_genome_search_does_not_reuse_positive_or_empty_coverage(tmp_path, legacy):
    old, new, plan, region = cache_fixture(tmp_path)
    directory = old / "rescued/T"
    predictions = json.loads((directory / "models.json").read_text())
    predictions.append({**predictions[0], "search": "genome_fallback"})
    (directory / "models.json").write_text(json.dumps(predictions))
    if legacy:
        plan = current_fixture_plan(legacy_fixture(old))
    else:
        contract = prediction_search_contract()
        contract["genome"]["output_score_ratio"] = .99
        (directory / SEARCH_CONTRACT_FILE).write_text(json.dumps(contract))
        rebind_fixture_plan(old, lambda p: p["request"]["tools"].update(prediction_search_contract=contract))
    cache = verify_prediction_cache(old, new, plan, "T", [region], PARAMS)
    assert cache.candidate_ids() == {"region1"}  # Positive/empty local searches are still reusable.
    assert not cache.genome_candidate_ids()  # Even a completed empty genome search must be rerun.
    assert not cache.genome_query_mapping()
    rows = list(cache.iter_models())
    assert [row["search"] for row in rows] == ["synteny_interval"]
    assert {"status", "model_id", "support", "problems"}.isdisjoint(rows[0])
    assert cache.genome_reuse["policy"] == "search_contract_mismatch"


def test_legacy_empty_local_result_reuses_only_local_coverage_without_logs(tmp_path):
    old, new, plan, region = cache_fixture(tmp_path)
    (old / "rescued/T/models.json").write_text("[]")
    plan = current_fixture_plan(legacy_fixture(old))
    cache = verify_prediction_cache(old, new, plan, "T", [region], PARAMS)
    assert not (old / "rescued/T/logs").exists()
    assert cache.candidate_ids() == {"region1"} and not cache.genome_candidate_ids()
    assert list(cache.iter_models()) == []


def inherited_fixture(tmp_path, *, parent_implementation=LEGACY_DIRECT_IMPLEMENTATION):
    import shutil
    old, new, _, region = cache_fixture(tmp_path)
    ancestor = tmp_path / "ancestor"
    shutil.copytree(old, ancestor)
    legacy_fixture(ancestor, parent_implementation)
    parent = frozen_prediction_cache_key(ancestor)
    plan = current_fixture_plan(legacy_fixture(old, LEGACY_INHERITED_IMPLEMENTATION, parent))
    return old, new, plan, region, ancestor


def test_transitive_legacy_only_local_proof_needs_no_command_logs(tmp_path):
    old, new, plan, region, ancestor = inherited_fixture(tmp_path)
    cache = verify_prediction_cache(old, new, plan, "T", [region], PARAMS)
    assert cache.local_search_compatible and not cache.genome_search_compatible
    assert cache.candidate_ids() == {"region1"}
    assert not cache.genome_candidate_ids()
    assert len(list(cache.iter_models())) == 1
    assert not (ancestor / "rescued/T/logs").exists() and not (old / "rescued/T/logs").exists()
    with (ancestor / "plan.json").open("a") as handle:
        handle.write("\n")
    with pytest.raises(OSError, match="changed"):
        cache.check()


@pytest.mark.parametrize("field", ["plan_sha256", "receipt_sha256", "models.json"])
def test_forged_frozen_ancestor_is_an_integrity_error(tmp_path, field):
    old, new, plan, region, _ = inherited_fixture(tmp_path)
    def forge(producer):
        parent = producer["request"]["prediction_cache"]
        if field == "plan_sha256":
            parent[field] = "0" * 64
        elif field == "receipt_sha256":
            parent["species"]["T"][field] = "0" * 64
        else:
            parent["species"]["T"]["files"][field] = "0" * 64
    rebind_fixture_plan(old, forge)
    with pytest.raises(ValueError, match="changed|differs"):
        verify_prediction_cache(old, new, plan, "T", [region], PARAMS)


@pytest.mark.parametrize("unknown", ["implementation", "miniprot_sha256", "ancestor"])
def test_unknown_legacy_scope_abstains_from_prediction_and_coverage_reuse(tmp_path, unknown):
    if unknown == "ancestor":
        old, new, plan, region, _ = inherited_fixture(tmp_path, parent_implementation="unknown")
    else:
        old, new, _, region = cache_fixture(tmp_path)
        legacy_fixture(old)
        legacy = rebind_fixture_plan(old, lambda p: p["request"]["tools"].update({unknown: "unknown"}))
        plan = current_fixture_plan(legacy)
    cache = verify_prediction_cache(old, new, plan, "T", [region], PARAMS)
    assert not cache.local_search_compatible and not cache.genome_search_compatible
    assert not cache.candidate_ids() and not cache.genome_candidate_ids()
    assert list(cache.iter_models()) == []


def test_mixed_inherited_and_current_legacy_local_bounds_are_not_complete_coverage(tmp_path):
    import shutil
    old, new, _, region = cache_fixture(tmp_path)
    ancestor = tmp_path / "ancestor"
    shutil.copytree(old, ancestor)
    contract = prediction_search_contract()
    contract["local"]["output_score_ratio"] = .25
    (ancestor / "rescued/T" / SEARCH_CONTRACT_FILE).write_text(json.dumps(contract))
    rebind_fixture_plan(ancestor, lambda p: p["request"]["tools"].update(prediction_search_contract=contract))
    plan = current_fixture_plan(legacy_fixture(old, LEGACY_INHERITED_IMPLEMENTATION, frozen_prediction_cache_key(ancestor)))
    cache = verify_prediction_cache(old, new, plan, "T", [region], PARAMS)
    assert cache.search_contract["local"] is None and not cache.candidate_ids()
    assert list(cache.iter_models()) == []


@pytest.mark.parametrize("problem", ["missing", "unbound", "invalid_ratio", "invalid_secondary", "bool_schema"])
def test_malformed_or_unbound_modern_contract_is_rejected(tmp_path, problem):
    old, new, plan, region = cache_fixture(tmp_path)
    path = old / "rescued/T" / SEARCH_CONTRACT_FILE
    contract = prediction_search_contract()
    if problem == "missing":
        path.unlink()
    else:
        if problem == "unbound":
            contract["local"]["output_score_ratio"] = .25
        elif problem == "invalid_ratio":
            contract["local"]["output_score_ratio"] = float("nan")
        elif problem == "invalid_secondary":
            contract["genome"]["max_secondary"] = True
        else:
            contract["schema"] = True
        path.write_text(json.dumps(contract))
    write_receipt(path.parent, json.loads((path.parent / "receipt.json").read_text())["key"])
    with pytest.raises(ValueError, match="contract"):
        verify_prediction_cache(old, new, plan, "T", [region], PARAMS)


def test_verified_cache_reuses_predictions_not_acceptance(tmp_path):
    old, new, plan, region = cache_fixture(tmp_path)
    key = frozen_prediction_cache_key(old)
    cache = verify_prediction_cache(old, new, plan, "T", [region], PARAMS, key)
    rows = list(cache.iter_models())
    assert cache.candidate_ids() == cache.genome_candidate_ids() == {"region1"}
    assert len(rows) == 1 and rows[0]["cds"] == [[0, 99, 0]]
    assert {"sequence", "status", "support", "quality_evidence", "problems", "model_id"}.isdisjoint(rows[0])


def test_frozen_cache_metadata_reads_plan_and_selected_receipt_once(tmp_path, monkeypatch):
    from workflow.support import rescue_prediction_cache as implementation
    old, _, _, _ = cache_fixture(tmp_path)
    original = implementation._json_snapshot
    reads = []
    def track(path):
        reads.append(Path(path))
        return original(path)
    monkeypatch.setattr(implementation, "_json_snapshot", track)
    frozen = frozen_prediction_cache_key(old, ["T"])
    assert reads == [old / "plan.json", old / "rescued/T/receipt.json"]
    assert frozen["plan_sha256"] == sha(old / "plan.json")
    assert frozen["species"]["T"]["receipt_sha256"] == sha(old / "rescued/T/receipt.json")


@pytest.mark.parametrize("surface", ["models.json", "candidates.json", "receipt.json"])
def test_cache_content_mutation_is_rejected(tmp_path, surface):
    old, new, plan, region = cache_fixture(tmp_path)
    key = frozen_prediction_cache_key(old)
    path = old / "rescued" / "T" / surface
    path.write_text(path.read_text() + " ")
    with pytest.raises(ValueError, match="changed"):
        verify_prediction_cache(old, new, plan, "T", [region], PARAMS, key)


@pytest.mark.parametrize("change", ["miniprot", "search", "source", "prepared"])
def test_cache_incompatible_predictor_inputs_are_rejected(tmp_path, change):
    old, new, plan, region = cache_fixture(tmp_path)
    key = frozen_prediction_cache_key(old)
    params = dict(PARAMS)
    plan = copy.deepcopy(plan)
    if change == "miniprot":
        plan["request"]["tools"]["miniprot"] = "2"
    elif change == "search":
        params["max_intron"] += 1
    elif change == "source":
        plan["request"]["files"][plan["request"]["sources"]["T"]["genome"]] = "b" * 64
    else:
        (new / "prepared" / "D" / "genes.pep").write_text(">D_q\nMAD\n")
    with pytest.raises(ValueError):
        verify_prediction_cache(old, new, plan, "T", [region], params, key)


def test_cache_candidate_change_is_rejected(tmp_path):
    old, new, plan, region = cache_fixture(tmp_path)
    region["end"] = 150
    with pytest.raises(ValueError, match="definition changed"):
        verify_prediction_cache(old, new, plan, "T", [region], PARAMS)


def test_cache_mutation_after_verification_is_rejected(tmp_path):
    old, new, plan, region = cache_fixture(tmp_path)
    cache = verify_prediction_cache(old, new, plan, "T", [region], PARAMS)
    path = old / "rescued" / "T" / "models.json"
    path.write_text(path.read_text() + " ")
    with pytest.raises(OSError, match="changed"):
        list(cache.iter_models())


def test_genome_cache_rebinds_exact_protein_from_different_donor(tmp_path):
    old, new, plan, region = cache_fixture(tmp_path)
    path = old / "rescued" / "T" / "models.json"
    rows = json.loads(path.read_text())
    rows[0]["search"] = "genome_fallback"
    rows[0]["paf_query"] = "region1"
    rows[0]["paf"] = "##PAF\tregion1\t33\t0\t33\t+\tchr1\t99\t0\t99\t33\t33\t60"
    path.write_text(json.dumps(rows))
    write_receipt(path.parent, json.loads((path.parent / "receipt.json").read_text())["key"])
    additional = {"id": "new_genome_query", "target": "T", "donor": "E", "query": "E_q", "genome_only": True}
    cache = verify_prediction_cache(old, new, plan, "T", [region, additional], PARAMS)
    predictions = list(cache.iter_models())
    assert cache.genome_candidate_ids() == {"region1", "new_genome_query"}
    assert cache.candidate_ids() == {"region1"}
    assert {row["query"] for row in predictions} == {"region1", "new_genome_query"}
    added = next(row for row in predictions if row["query"] == "new_genome_query")
    assert added["evidence"] == additional
    assert added["paf"].startswith("##PAF\tnew_genome_query\t")
    assert added["cds"] == rows[0]["cds"] and "sequence" not in added
    assert cache.genome_reuse["new_genome_only_queries"] == 1
    assert cache.genome_query_mapping() == {"region1": "new_genome_query", "new_genome_query": "new_genome_query"}


def test_exact_genome_empty_search_is_reused_without_interval_results(tmp_path):
    old, new, plan, region = cache_fixture(tmp_path)
    directory = old / "rescued" / "T"
    (directory / "models.json").write_text("[]")
    write_receipt(directory, json.loads((directory / "receipt.json").read_text())["key"])
    additional = {"id": "new_genome_query", "target": "T", "donor": "E", "query": "E_q", "genome_only": True}
    cache = verify_prediction_cache(old, new, plan, "T", [additional], PARAMS)
    assert not list(cache.iter_models())
    assert not cache.candidate_ids() and cache.genome_candidate_ids() == {"new_genome_query"}


def test_interval_search_cannot_authorise_sequence_keyed_genome_reuse(tmp_path):
    old, new, plan, region = cache_fixture(tmp_path)
    directory = old / "rescued" / "T"
    (directory / "genome_query_mapping.tsv").unlink()
    write_receipt(directory, json.loads((directory / "receipt.json").read_text())["key"])
    additional = {"id": "new_genome_query", "target": "T", "donor": "E", "query": "E_q", "genome_only": True}
    cache = verify_prediction_cache(old, new, plan, "T", [region, additional], PARAMS)
    assert not cache.genome_candidate_ids()
    assert [row["query"] for row in cache.iter_models()] == ["region1"]


def test_cache_uses_pristine_coordinates_before_terminal_completion(tmp_path):
    old, new, plan, region = cache_fixture(tmp_path)
    directory = old / "rescued" / "T"
    path = directory / "models.json"
    models = json.loads(path.read_text())
    raw = {key: value for key, value in models[0].items() if key not in {"sequence", "status", "problems", "model_id", "support"}}
    models[0]["raw_prediction"] = raw
    models[0]["cds"] = [[0, 105, 0]]
    path.write_text(json.dumps(models))
    write_receipt(directory, json.loads((directory / "receipt.json").read_text())["key"])
    cache = verify_prediction_cache(old, new, plan, "T", [region], PARAMS)
    assert list(cache.iter_models())[0]["cds"] == [[0, 99, 0]]


def test_completed_terminal_model_without_raw_prediction_is_rejected(tmp_path):
    old, new, plan, region = cache_fixture(tmp_path)
    directory = old / "rescued" / "T"
    path = directory / "models.json"
    models = json.loads(path.read_text())
    models[0]["terminal_completion"] = {"status": "completed", "selected_cds": [[0, 105, 0]]}
    path.write_text(json.dumps(models))
    write_receipt(directory, json.loads((directory / "receipt.json").read_text())["key"])
    cache = verify_prediction_cache(old, new, plan, "T", [region], PARAMS)
    with pytest.raises(ValueError, match="pristine cached prediction"):
        list(cache.iter_models())


def test_cache_unsafe_receipt_path_is_rejected(tmp_path):
    old, _, _, _ = cache_fixture(tmp_path)
    path = old / "rescued" / "T" / "receipt.json"
    receipt = json.loads(path.read_text())
    receipt["files"]["../../escape"] = "a" * 64
    path.write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match="Unsafe"):
        frozen_prediction_cache_key(old)


def test_rescue_config_overrides_are_exported_and_forwarded_to_plan(tmp_path):
    from workflow.support.entrypoint_config_schema import parse_entrypoint
    root = Path(__file__).resolve().parents[2]
    new_parameters = {"prediction_cache": ("", "old verified predictions"), "max_genome_queries": ("20000", "17"),
                      "unanchored_min_species": ("2", "3"), "terminal_max_extension": ("300", "180"),
                      "terminal_max_unaligned_c_overhang": ("2", "1")}
    metadata = {parameter.name: parameter for parameter in parse_entrypoint("gg_input_generation_entrypoint.sh")}
    overrides = {}
    for suffix, (default, value) in new_parameters.items():
        name = "gene_model_rescue_" + suffix
        assert metadata[name].default == default
        assert metadata[name].environment == "GG_INPUT_" + name.upper()
        overrides[metadata[name].environment] = value
    core = (root / "workflow/core/gg_input_generation_core.sh").read_text()
    defaults = core[core.index('gene_model_rescue_tree='):core.index('require_cds=')]
    function = core[core.index('prepare_gene_model_rescue() {'):core.index('\ngene_model_rescue_busco_species() {')]
    tree = tmp_path / "tree.nwk"
    tree.write_text("(T,D);\n")
    capture = tmp_path / "arguments"
    commands = ["set -euo pipefail", defaults,
                "source " + shlex.quote(str(root / "workflow/support/gg_entrypoint_config_vars.sh")),
                "source " + shlex.quote(str(root / "workflow/support/gg_util/01_runtime_config.sh")),
                "source " + shlex.quote(str(root / "workflow/support/gg_util/02_container_scheduler.sh")),
                "gg_apply_registered_env_overrides gg_input_generation_entrypoint.sh",
                "forward_config_vars_to_container_env gg_input_generation_entrypoint.sh", function,
                'python() { printf "%s\\n" "$@" > ' + shlex.quote(str(capture)) + '; }',
                "gene_model_rescue_tree=" + shlex.quote(str(tree)), "gg_support_dir=/support",
                "species_cds_dir=/cds", "species_gff_dir=/gff", "species_genome_dir=/genome",
                "species_busco_short_dir=/busco", "gg_workspace_input_dir=" + shlex.quote(str(tmp_path)),
                "prepare_gene_model_rescue"]
    for suffix, (_, value) in new_parameters.items():
        name = "gene_model_rescue_" + suffix
        commands.append('[[ "${SINGULARITYENV_' + name + '}" == ' + shlex.quote(value) + ' ]]')
        commands.append('[[ "${APPTAINERENV_' + name + '}" == ' + shlex.quote(value) + ' ]]')
    result = subprocess.run(["bash", "-c", "\n".join(commands)], env={**os.environ, **overrides},
                            capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    arguments = capture.read_text().splitlines()
    for suffix, (_, value) in new_parameters.items():
        flag = "--" + suffix.replace("_", "-")
        assert arguments[arguments.index(flag) + 1] == value

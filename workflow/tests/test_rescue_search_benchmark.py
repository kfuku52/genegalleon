"""Benchmark provenance and retention-free replay, without expensive searches."""
import gzip
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "benchmarks/benchmark_rescue_search.py"
SPEC = importlib.util.spec_from_file_location("rescue_search_benchmark_contract", SCRIPT)
benchmark = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(benchmark)


def publish(root, plan, source):
    (root / "plan.json").write_text(json.dumps(plan))
    files = {p.relative_to(source).as_posix(): benchmark.rescue.digest(p)
             for p in source.rglob("*") if p.is_file() and p.name != "receipt.json"}
    receipt = {"key": {"plan": benchmark.rescue.digest(root / "plan.json"), "species": "Target", "prepared": {}}, "files": files}
    (source / "receipt.json").write_text(json.dumps(receipt))
    return receipt


def fixture(tmp_path, *, compact=False, legacy=False, ratio=0.5, secondary=30):
    root = tmp_path / "evidence"
    source = root / "rescued/Target"
    (source / "intervals/1").mkdir(parents=True)
    (source / "intervals/1/models.gff").write_text("# original local alignment\n")
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nATGAAATAA\n")
    contract = benchmark.prediction_cache.prediction_search_contract()
    contract["genome"] = {"output_score_ratio": ratio, "max_secondary": secondary}
    tools = {"prediction_search_contract": contract}
    if legacy:
        tools = {"implementation": benchmark.prediction_cache.LEGACY_DIRECT_IMPLEMENTATION,
                 "miniprot_sha256": benchmark.prediction_cache.LEGACY_MINIPROT_SHA256}
    else:
        (source / benchmark.prediction_cache.SEARCH_CONTRACT_FILE).write_text(json.dumps(contract))
    plan = {"species": ["Target"], "request": {"tools": tools, "parameters": {"max_intron": 20000},
                "sources": {"Target": {"genome": str(genome), "genetic_code": 1}},
                "files": {str(genome): benchmark.rescue.digest(genome)}}}
    if compact:
        for name in ["local_search_inputs.jsonl.gz", "genome_search_inputs.jsonl.gz"]:
            with gzip.open(source / name, "wt") as handle:
                handle.write(json.dumps({"schema": 1}) + "\n")
        (source / "genome_prediction_query_mapping.tsv").write_text("candidate\trepresentative\nquery\tquery\n")
        (source / "genome.unique.gff").write_text("# original expanded alignment\n")
        (root / "plan.json").write_text(json.dumps(plan))
        (source / "search_inputs.json").write_text(json.dumps({"schema": 1,
            "plan_sha256": benchmark.rescue.digest(root / "plan.json"), "retained_search_inputs": False,
            "source_files": {str(genome): benchmark.rescue.digest(genome)},
            "expanded_genome_gff_sha256": benchmark.rescue.digest(source / "genome.unique.gff")}))
    else:
        (source / "intervals/1/region.fa").write_text(">interval\nATGAAATAA\n")
        (source / "intervals/1/queries.fa").write_text(">query\nMK\n")
        (source / "unresolved.fa").write_text(">query\nMK\n")
        (source / "genome.gff").write_text("# original expanded alignment\n")
    publish(root, plan, source)
    return root, source, plan


def fake_export(monkeypatch, source, *, mutate=None):
    calls = []
    def export(root, plan, name, destination, *, combined, cpus):
        calls.append((destination, combined, cpus))
        (destination / "intervals/1").mkdir(parents=True)
        (destination / "intervals/1/region.fa").write_text(">interval\nATGAAATAA\n")
        (destination / "intervals/1/queries.fa").write_text(">query\nMK\n")
        (destination / "unresolved.fa").write_text(">query\nMK\n")
        (destination / "genome.gff").write_bytes((source / "genome.unique.gff").read_bytes())
        files = {p.relative_to(destination).as_posix(): benchmark.rescue.digest(p) for p in destination.rglob("*") if p.is_file()}
        proof = {"origin_receipt_sha256": benchmark.rescue.digest(source / "receipt.json"),
                 "plan_sha256": benchmark.rescue.digest(root / "plan.json"), "prepared_receipts": {},
                 "local_windows": 1, "combined_diagnostics": True, "files": files}
        if mutate:
            mutate(proof, destination)
        (destination / "export_receipt.json").write_text(json.dumps(proof))
        return destination
    monkeypatch.setattr(benchmark.rescue, "export_worker_inputs", export)
    return calls


def test_retention_free_replay_is_bound_and_never_writes_source(tmp_path, monkeypatch):
    root, source, _ = fixture(tmp_path, compact=True)
    before = {p.relative_to(source): p.read_bytes() for p in source.rglob("*") if p.is_file()}
    calls = fake_export(monkeypatch, source)
    output = tmp_path / "benchmark/inputs"
    evidence = benchmark.FrozenEvidence(root, "Target", export_directory=output, cpus=4)
    assert calls == [(output, True, 4)]
    assert evidence.interval_ids() == [1]
    assert evidence.inputs == output and evidence.fallback_recorded
    assert evidence.checked(source / "intervals/1/models.gff").read_text() == "# original local alignment\n"
    evidence.recheck()
    assert before == {p.relative_to(source): p.read_bytes() for p in source.rglob("*") if p.is_file()}


@pytest.mark.parametrize("field", ["origin_receipt_sha256", "plan_sha256", "local_windows", "combined_diagnostics", "prepared_receipts"])
def test_foreign_export_proof_rejected(tmp_path, monkeypatch, field):
    root, source, _ = fixture(tmp_path, compact=True)
    fake_export(monkeypatch, source, mutate=lambda proof, dest: proof.__setitem__(field, "foreign"))
    with pytest.raises(ValueError, match="foreign export"):
        benchmark.FrozenEvidence(root, "Target", export_directory=tmp_path / "output/inputs")


@pytest.mark.parametrize("member", ["intervals/1/region.fa", "genome.gff"])
def test_export_member_mutation_after_receipt_rejected(tmp_path, monkeypatch, member):
    root, source, _ = fixture(tmp_path, compact=True)
    fake_export(monkeypatch, source, mutate=lambda proof, dest: (dest / member).write_text("changed"))
    with pytest.raises(ValueError, match="source changed"):
        benchmark.FrozenEvidence(root, "Target", export_directory=tmp_path / "output/inputs")


def test_expanded_genome_checksum_is_original_authority(tmp_path, monkeypatch):
    root, source, plan = fixture(tmp_path, compact=True)
    metadata = json.loads((source / "search_inputs.json").read_text())
    metadata["expanded_genome_gff_sha256"] = "0" * 64
    (source / "search_inputs.json").write_text(json.dumps(metadata))
    publish(root, plan, source)
    fake_export(monkeypatch, source)
    with pytest.raises(ValueError, match="original evidence"):
        benchmark.FrozenEvidence(root, "Target", export_directory=tmp_path / "output/inputs")


@pytest.mark.parametrize("missing", ["unresolved.fa", "genome.gff", "intervals/1/queries.fa"])
def test_retained_members_missing_from_receipt_are_not_valid_omissions(tmp_path, missing):
    root, source, _ = fixture(tmp_path)
    receipt = json.loads((source / "receipt.json").read_text())
    del receipt["files"][missing]
    (source / "receipt.json").write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match="Incomplete"):
        benchmark.FrozenEvidence(root, "Target", export_directory=tmp_path / "output/inputs")


@pytest.mark.parametrize("missing", ["genome_prediction_query_mapping.tsv", "genome.unique.gff", "genome_search_inputs.jsonl.gz"])
def test_retention_free_fallback_must_have_complete_frozen_members(tmp_path, missing):
    root, source, _ = fixture(tmp_path, compact=True)
    receipt = json.loads((source / "receipt.json").read_text())
    del receipt["files"][missing]
    (source / "receipt.json").write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match="Incomplete fallback"):
        benchmark.FrozenEvidence(root, "Target", export_directory=tmp_path / "output/inputs")


def test_retention_flag_and_explicit_export_directory_required(tmp_path):
    root, source, plan = fixture(tmp_path, compact=True)
    with pytest.raises(ValueError, match="valid omission"):
        benchmark.FrozenEvidence(root, "Target")
    metadata = json.loads((source / "search_inputs.json").read_text())
    metadata["retained_search_inputs"] = True
    (source / "search_inputs.json").write_text(json.dumps(metadata))
    publish(root, plan, source)
    with pytest.raises(ValueError, match="Missing original interval"):
        benchmark.FrozenEvidence(root, "Target", export_directory=tmp_path / "output/inputs")


def test_genome_runner_uses_recorded_legacy_or_modern_bounds(tmp_path, monkeypatch):
    root, _, _ = fixture(tmp_path, legacy=True)
    evidence = benchmark.FrozenEvidence(root, "Target")
    assert evidence.contract["genome"] == {"output_score_ratio": 0.99, "max_secondary": 30}
    commands = []
    def run(command, directory, label, stdout=None):
        commands.append(command)
        if label == "map":
            stdout.write_text("# original expanded alignment\n")
            (directory / "logs").mkdir()
            (directory / "logs/map.log").write_text("Peak RSS: 0.01 GB\n")
    monkeypatch.setattr(benchmark.rescue, "run", run)
    monkeypatch.setattr(benchmark.rescue, "expand_miniprot_queries", lambda src, dst, mapping: dst.write_bytes(src.read_bytes()))
    output = tmp_path / "benchmark"
    output.mkdir()
    result = {"samples": {"genome_unique": []}}
    args = SimpleNamespace(output=output, species="Target", cpus=2, query_count=1, check_existing=True, repeats=1)
    benchmark.benchmark_genome(args, evidence, result)
    mapping = next(c for c in commands if "--gff" in c)
    assert "--outs=0.99" in mapping and mapping[mapping.index("-N") + 1] == 30
    evidence.recheck()


def test_custom_modern_genome_contract_is_not_replaced_by_defaults(tmp_path):
    root, source, plan = fixture(tmp_path, ratio=0.7, secondary=17)
    plan["request"]["species_profiles"] = {"Target": {"max_intron": 777}}
    publish(root, plan, source)
    evidence = benchmark.FrozenEvidence(root, "Target")
    assert evidence.contract["genome"] == {"output_score_ratio": 0.7, "max_secondary": 17}
    assert evidence.parameters["max_intron"] == 777


def test_unknown_legacy_and_changed_contract_fail(tmp_path):
    root, source, plan = fixture(tmp_path, legacy=True)
    plan["request"]["tools"]["implementation"] = "unknown"
    publish(root, plan, source)
    with pytest.raises(ValueError, match="Unknown producer"):
        benchmark.FrozenEvidence(root, "Target")


@pytest.mark.parametrize("forged", [False, True])
def test_inherited_legacy_contract_binds_parent_plan_and_receipt(tmp_path, forged):
    root, source, plan = fixture(tmp_path, legacy=True)
    parent = tmp_path / "ancestor"
    producer = parent / "rescued/Target"
    producer.mkdir(parents=True)
    parent_plan = {"species": ["Target"], "request": {"tools": dict(plan["request"]["tools"])}}
    (parent / "plan.json").write_text(json.dumps(parent_plan))
    (producer / "models.json").write_text("[]")
    (producer / "candidates.json").write_text("[]")
    parent_sha = benchmark.rescue.digest(parent / "plan.json")
    files = {name: benchmark.rescue.digest(producer / name) for name in ("models.json", "candidates.json")}
    (producer / "receipt.json").write_text(json.dumps({"key": {"plan": parent_sha, "species": "Target"}, "files": files}))
    plan["request"]["tools"].update({"implementation": benchmark.prediction_cache.LEGACY_INHERITED_IMPLEMENTATION,
        "verify_prediction_cache_implementation": benchmark.prediction_cache.LEGACY_INHERITED_VERIFIER})
    plan["request"]["prediction_cache"] = {"schema": 1, "root": str(parent),
        "plan_sha256": "f" * 64 if forged else parent_sha,
        "species": {"Target": {"receipt_sha256": benchmark.rescue.digest(producer / "receipt.json"), "files": files}}}
    publish(root, plan, source)
    if forged:
        with pytest.raises(ValueError, match="source changed"):
            benchmark.FrozenEvidence(root, "Target")
    else:
        evidence = benchmark.FrozenEvidence(root, "Target")
        assert evidence.contract["genome"]["output_score_ratio"] == .99
        assert evidence.hashes[str(parent / "plan.json")] == parent_sha
        (producer / "receipt.json").write_text((producer / "receipt.json").read_text() + " ")
        with pytest.raises(ValueError, match="source changed"):
            evidence.recheck()


def test_export_and_original_hashes_are_rechecked(tmp_path, monkeypatch):
    root, source, _ = fixture(tmp_path, compact=True)
    fake_export(monkeypatch, source)
    evidence = benchmark.FrozenEvidence(root, "Target", export_directory=tmp_path / "output/inputs")
    (evidence.inputs / "unresolved.fa").write_text(">changed\nMA\n")
    with pytest.raises(ValueError, match="source changed"):
        evidence.recheck()


def test_export_destination_cannot_be_inside_source(tmp_path):
    root, _, _ = fixture(tmp_path, compact=True)
    with pytest.raises(ValueError, match="outside original"):
        benchmark.FrozenEvidence(root, "Target", export_directory=root / "benchmark-export")

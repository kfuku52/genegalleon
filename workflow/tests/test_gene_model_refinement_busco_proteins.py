"""Protein admission and historical-result reuse with synthetic BUSCO outputs."""
import ast
import copy
import hashlib
import json
import os
import shutil
import subprocess
from importlib import import_module
from pathlib import Path
from types import SimpleNamespace

import pytest

from workflow.tests.test_gene_model_refinement_busco import summary
from workflow.tests.test_plot_gene_model_refinement import completed

busco = import_module("gene_model_refinement_busco")
reuse_api = import_module("refinement_busco_reuse")
inputs_api = import_module("refinement_busco_inputs")


def test_version_probe_isolates_import_time_writes_and_preserves_caller(tmp_path, monkeypatch):
    """An actual child reproduces BUSCO's logger on a read-only launch cwd."""
    executable = tmp_path / "busco"
    executable.write_text("#!/bin/sh\nmkdir -p \"$MPLCONFIGDIR/cache\"\nprintf log > busco.log\nprintf 'BUSCO fixture\\n'\n")
    executable.chmod(0o700)
    caller = tmp_path / "sealed-package"
    caller.mkdir(mode=0o500)
    monkeypatch.chdir(caller)
    monkeypatch.setenv("PATH", str(tmp_path) + os.pathsep + os.environ["PATH"])
    monkeypatch.setenv("MPLCONFIGDIR", str(caller / "no-cache"))
    original = dict(os.environ)
    real = subprocess.check_output
    scratch = []

    def observe(command, **options):
        scratch.append(Path(options["cwd"]))
        assert options["timeout"] == 30
        assert options["env"]["MPLCONFIGDIR"] == options["cwd"]
        assert Path(options["cwd"]).stat().st_mode & 0o077 == 0
        return real(command, **options)

    monkeypatch.setattr(busco.subprocess, "check_output", observe)
    assert busco.busco_version() == "BUSCO fixture"
    assert Path.cwd() == caller and not list(caller.iterdir())
    assert dict(os.environ) == original and all(not path.exists() for path in scratch)


@pytest.mark.parametrize("failure", [subprocess.CalledProcessError(1, ["busco"]),
                                    subprocess.TimeoutExpired(["busco"], 30)])
def test_version_probe_propagates_failure_without_retry_and_cleans_scratch(monkeypatch, failure):
    scratch = []

    def fail(command, **options):
        scratch.append(Path(options["cwd"]))
        raise failure

    monkeypatch.setattr(busco.subprocess, "check_output", fail)
    with pytest.raises(type(failure)):
        busco.busco_version()
    assert len(scratch) == 1 and not scratch[0].exists()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.mark.parametrize("sequence,reason", [("ATGAAAA", "partial_terminal_codon"),
    ("ATGTGATAA", "internal_stop"), ("ATG???TAA", "invalid_base")])
def test_strict_protein_translation_preserves_rejection(sequence, reason, tmp_path):
    source = tmp_path / "source.fa"
    source.write_text(">gene description\n" + sequence + "\n")
    audit = busco.translate_cds(source, tmp_path / "protein.fa", 1)
    assert audit["status"] == "rejected"
    assert reason in audit["invalid_records"][0]["reasons"]


def test_genetic_code_ambiguity_and_duplicate_identifiers_are_not_guessed(tmp_path):
    source = tmp_path / "source.fa"
    source.write_text(">gene\nATGTGATAA\n")
    audit = busco.translate_cds(source, tmp_path / "protein.fa", 2)
    assert audit["status"] == "translated"
    assert (tmp_path / "protein.fa").read_text() == ">gene\nMW\n"
    source.write_text(">gene\nATGAAATAA\n>gene alias\nATGAAATAA\n")
    audit = busco.translate_cds(source, tmp_path / "protein.fa", 1)
    assert audit["status"] == "rejected"
    assert "duplicate_identifier" in audit["invalid_records"][0]["reasons"]
    with pytest.raises(ValueError, match="Context-dependent"):
        busco.protein_genetic_code(27)


def test_native_admitted_view_uses_completed_publications_and_rejects_tampering(tmp_path):
    root = completed(tmp_path)
    pairs = busco.input_pairs(root)
    codes = busco.read_genetic_codes(pairs, root)
    destination = tmp_path / "admitted"
    descriptors = inputs_api.prepare_admitted_protein_inputs(pairs, root, destination, codes)
    assert inputs_api.prepare_admitted_protein_inputs(pairs, root, destination, codes) == descriptors
    for pair in pairs:
        name = pair["species"]
        descriptor = descriptors[name]["after"]
        assert Path(descriptor["protein_path"]).read_bytes() == (root / "effective/species_protein" / (name + ".fa")).read_bytes()
        audit = busco._validate_admitted_input(descriptor, sha(pair["after"]), codes[name], species=name, phase="after")
        assert audit["counts"]["records"] == audit["counts"]["included"] + audit["counts"]["excluded"]
    target = Path(next(iter(descriptors.values()))["before"]["protein_path"])
    target.write_text(">tampered\nMK\n")
    with pytest.raises(ValueError):
        inputs_api.prepare_admitted_protein_inputs(pairs, root, destination, codes)


def synthetic_descriptor(root, dna, peptide="MK"):
    root.mkdir(parents=True)
    source, protein, audit_path = root / "source.fa", root / "protein.fa", root / "audit.json"
    source.write_text(">gene\n" + dna + "\n")
    protein.write_text(">gene\n" + peptide + "\n")
    audit = {"schema":1, "status":"translated", "species":"Species_target", "phase":"before",
             "policy":inputs_api.ADMITTED_PROTEIN_POLICY, "genetic_code":1, "source_sha256":sha(source),
             "protein_sha256":sha(protein), "counts":{"records":1,"included":1,"excluded":0}, "excluded_records":[]}
    audit_path.write_text(json.dumps(audit))
    descriptor = {"protein_path":str(protein), "audit_path":str(audit_path), "protein_sha256":sha(protein),
                  "audit_sha256":sha(audit_path), "source_sha256":sha(source), "policy":inputs_api.ADMITTED_PROTEIN_POLICY}
    pair = {"species":"Species_target", "before":str(source), "after":str(source),
            "protein_inputs":{"before":descriptor}, "refinement_status":"analysed", "reason":""}
    return pair


def archived_producer(directory):
    """Synthetic prior renderer with the exact computational function bodies."""
    directory.mkdir()
    current = Path(busco.__file__)
    raw = current.read_text()
    tree = ast.parse(raw)
    sections = []
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and (node.name in reuse_api.KERNELS or node.name == "_run_one"):
            if node.name == "run_one":
                continue
            source = ast.get_source_segment(raw, node)
            if node.name == "_run_one":
                source = source.replace("def _run_one(", "def run_one(", 1)
            sections.append(source)
        elif isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "PROTEIN_TRANSLATION_POLICY" for t in node.targets):
            sections.append(ast.get_source_segment(raw, node))
    producer = directory / "producer.py"
    producer.write_text("\n\n".join(sections) + "\n# An older rendering implementation.\n")
    for name in reuse_api.DEPENDENCIES:
        shutil.copyfile(current.parent / name, directory / name)
    return producer


@pytest.fixture
def reuse_case(tmp_path, monkeypatch):
    previous = tmp_path / "previous"
    previous.mkdir()
    old_pair = synthetic_descriptor(tmp_path / "old_inputs", "ATGAAATAA")
    # Different DNA and input/audit paths can still describe identical proteins.
    pair = synthetic_descriptor(tmp_path / "new_inputs", "ATGAAGTAA")
    producer = archived_producer(tmp_path / "producer")
    contract = {"schema":1, "mode":"proteins", "busco_version":"fixture", "lineage_sha256":"a"*64,
                "lineage_path":"/lineage", "download_path":"/downloads", "evalue":"1e-03", "limit":20,
                "cpus_per_job":4, "tool_sha256":{"busco":"b"*64,"metaeuk":"c"*64,"hmmsearch":"d"*64},
                "implementation":sha(busco.__file__), "wrapper_sha256":"e"*64,"hmmsearch_wrapper_sha256":"f"*64,
                "translation":{"input_kind":"admitted","policy":inputs_api.ADMITTED_PROTEIN_POLICY,
                               "genetic_codes":{"Species_target":1}}}
    old_contract = copy.deepcopy(contract)
    old_contract.update(implementation=sha(producer), lineage_path="/old/lineage", download_path="/old/downloads")
    (previous / "contract.json").write_text(json.dumps({"contract":old_contract,"pairs":[old_pair]}))
    calls = []

    def fake_run(command, **kwargs):
        assert command[command.index("--mode") + 1] == "proteins"
        calls.append(Path(command[command.index("--in") + 1]).read_text())
        scratch = Path(command[command.index("--out_path") + 1])
        result = scratch / "busco"
        (result / "run_fixture").mkdir(parents=True)
        (result / "short_summary.specific.fixture.txt").write_text(summary().replace("euk_tran", "proteins"))
        (result / "run_fixture/full_table.tsv").write_text("# fixture\ngroup\tComplete\tgene\n")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(busco.subprocess, "run", fake_run)
    busco._run_one(old_pair, "before", previous, old_contract, 4)
    return SimpleNamespace(previous=previous, producer=producer, contract=contract, old_contract=old_contract,
                           pair=pair, old_pair=old_pair, report=tmp_path / "new_report", calls=calls)


def context(case):
    return reuse_api.prepare_reuse(case.previous, case.producer, case.contract, Path(busco.__file__))


def test_identical_proteins_reuse_counts_but_keep_new_dna_admission_and_audit(reuse_case):
    c = reuse_case
    old_directory = c.previous / "runs/Species_target/before"
    original = {p.name:sha(p) for p in old_directory.iterdir() if p.is_file()}
    result = busco.run_one(c.pair, "before", c.report, c.contract, 4, reuse=context(c))
    assert result["mode"] == "proteins" and len(c.calls) == 1
    directory = c.report / "runs/Species_target/before"
    saved = json.loads((directory / "receipt.json").read_text())
    assert saved["key"]["source_sha256"] == sha(c.pair["before"]) != sha(c.old_pair["before"])
    assert saved["translation_audit_sha256"] == c.pair["protein_inputs"]["before"]["audit_sha256"]
    assert saved["reused_verified_previous"]["producer_sha256"] == sha(c.producer)
    assert sha(directory / "full_table.tsv") == sha(old_directory / "full_table.tsv")
    assert {p.name:sha(p) for p in old_directory.iterdir() if p.is_file()} == original
    busco.run_one(c.pair, "before", c.report, c.contract, 4, reuse=context(c))
    assert len(c.calls) == 1


@pytest.mark.parametrize("fault", ["producer", "kernel", "dependency", "runtime", "lineage", "translation"])
def test_reuse_rejects_unverified_computation_identity(reuse_case, fault):
    c = reuse_case
    if fault == "producer":
        c.producer.write_text(c.producer.read_text() + "\n")
    if fault == "kernel":
        c.producer.write_text(c.producer.read_text().replace('"--limit", "20"', '"--limit", "21"'))
        c.old_contract["implementation"] = sha(c.producer)
        (c.previous / "contract.json").write_text(json.dumps({"contract":c.old_contract,"pairs":[c.old_pair]}))
    if fault == "dependency":
        (c.producer.parent / "fasta_sequence_store.py").write_text("# changed parser\n")
    if fault == "runtime":
        c.contract["tool_sha256"]["hmmsearch"] = "0"*64
    if fault == "lineage":
        c.contract["lineage_sha256"] = "0"*64
    if fault == "translation":
        c.contract["translation"]["genetic_codes"]["Species_target"] = 2
    with pytest.raises(ValueError):
        context(c)
    assert not c.report.exists()


@pytest.mark.parametrize("member", ["summary.txt", "full_table.tsv", "protein.faa", "translation_audit.json"])
def test_reuse_rejects_changed_old_leaf_before_publishing(reuse_case, member):
    c = reuse_case
    (c.previous / "runs/Species_target/before" / member).write_text("tampered\n")
    with pytest.raises(ValueError, match="leaf changed"):
        busco.run_one(c.pair, "before", c.report, c.contract, 4, reuse=context(c))
    assert not (c.report / "runs/Species_target/before/receipt.json").exists()
    assert len(c.calls) == 1


def test_different_valid_protein_runs_new_busco_instead_of_importing_old_result(reuse_case):
    c = reuse_case
    pair = synthetic_descriptor(c.report.parent / "different_inputs", "ATGTTTTAA", "MF")
    busco.run_one(pair, "before", c.report, c.contract, 4, reuse=context(c))
    assert c.calls == [">gene\nMK\n", ">gene\nMF\n"]
    saved = json.loads((c.report / "runs/Species_target/before/receipt.json").read_text())
    assert "reused_verified_previous" not in saved


def test_changed_current_dna_cannot_use_a_stale_admitted_descriptor(reuse_case):
    c = reuse_case
    Path(c.pair["before"]).write_text(">gene\nATGTTTTAA\n")
    with pytest.raises(ValueError, match="DNA collection"):
        busco.run_one(c.pair, "before", c.report, c.contract, 4, reuse=context(c))
    assert not (c.report / "runs/Species_target/before/receipt.json").exists()


def test_cpu_count_must_match_the_recorded_protein_computation(reuse_case):
    c = reuse_case
    with pytest.raises(ValueError, match="CPU count"):
        busco.run_one(c.pair, "before", c.report, c.contract, 8, reuse=context(c))
    assert not c.report.exists()


def test_evaluate_builds_explicit_protein_contract_and_assesses_both_phases(reuse_case, monkeypatch):
    c = reuse_case
    pair = copy.deepcopy(c.pair)
    before = pair.pop("protein_inputs")["before"]
    after = copy.deepcopy(before)
    audit_path = Path(after["audit_path"]).with_name("after_audit.json")
    audit = json.loads(Path(after["audit_path"]).read_text())
    audit["phase"] = "after"
    audit_path.write_text(json.dumps(audit))
    after.update(audit_path=str(audit_path), audit_sha256=sha(audit_path))
    lineage = c.report.parent / "lineage"
    lineage.mkdir()
    (lineage / "dataset.cfg").write_text("fixture\n")
    monkeypatch.setattr(busco.subprocess, "check_output", lambda *a, **k: "BUSCO fixture")
    monkeypatch.setattr(busco.shutil, "which", lambda name: busco.__file__)
    monkeypatch.setattr(busco, "plot_comparison", lambda *a: None)
    rows = busco.evaluate([pair], c.report, lineage, c.report.parent / "downloads", cpus=4, jobs=1,
                         mode="proteins", genetic_codes={pair["species"]:1},
                         protein_inputs={pair["species"]:{"before":before,"after":after}})
    assert len(c.calls) == 3 and len(rows) == 1
    frozen = json.loads((c.report / "contract.json").read_text())["contract"]
    assert frozen["mode"] == "proteins" and frozen["translation"]["input_kind"] == "admitted"
    assert rows[0]["before_result"]["mode"] == rows[0]["after_result"]["mode"] == "proteins"
    assert rows[0]["delta_complete"] == 0


def test_cds_only_species_requires_explicit_code_and_native_code_cannot_change(tmp_path):
    root = completed(tmp_path)
    pairs = busco.input_pairs(root)
    extra = {"species":"Drosophyllum_lusitanicum"}
    with pytest.raises(ValueError, match="CDS-only"):
        busco.read_genetic_codes([*pairs, extra], root)
    tsv = tmp_path / "codes.tsv"
    tsv.write_text("species\tgenetic_code\nDrosophyllum_lusitanicum\t1\n")
    assert busco.read_genetic_codes([*pairs, extra], root, tsv)[extra["species"]] == 1
    name = pairs[0]["species"]
    tsv.write_text("species\tgenetic_code\n" + name + "\t2\n")
    with pytest.raises(ValueError, match="frozen refinement plan"):
        busco.read_genetic_codes(pairs, root, tsv)

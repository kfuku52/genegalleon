"""Opt-in six-clade qualification in an allocated GPU container; no quality claim."""
import argparse
import importlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))
selection = importlib.import_module("annevo_model_selection")
cache = importlib.import_module("annevo_runtime_cache")
runtime = importlib.import_module("annevo_runtime")


@pytest.fixture(scope="module")
def qualification(tmp_path_factory):
    locator = os.environ.get("GG_ANNEVO_GPU_FIXTURES", "")
    if not locator:
        pytest.skip("requires an allocated CUDA GPU and a sealed GG_ANNEVO_GPU_FIXTURES manifest")
    document = json.loads(Path(locator).read_text())
    assert set(document["fixtures"]) == set(selection.PROFILES)
    for profile, fixture in document["fixtures"].items():
        assert fixture["profile"] == profile
        assert selection.file_sha(fixture["genome"]) == fixture["sha256"]
        assert 0 < Path(fixture["genome"]).stat().st_size <= 1100000
    assert selection.file_sha(document["taxonomy_db"]) == document["taxonomy_sha256"]
    root = tmp_path_factory.mktemp("annevo-six-clades")
    mapping = root / "taxids.tsv"
    mapping.write_text("species\ttaxid\n" + "".join(f"{f['species']}\t{f['taxid']}\n" for f in document["fixtures"].values()))
    return document, root, mapping, {}


def invoke(root, name, script, *args):
    with (root / (name + ".stdout.txt")).open("w") as stdout, (root / (name + ".stderr.txt")).open("w") as stderr:
        result = subprocess.run([sys.executable, str(SUPPORT / script), *map(str, args)],
                                stdout=stdout, stderr=stderr, timeout=600)
    assert result.returncode == 0, (root / (name + ".stderr.txt")).read_text()[-4000:]


@pytest.mark.parametrize("profile", selection.PROFILES)
def test_cold_model_warm_offline_gff_and_cds(profile, qualification, monkeypatch):
    # Only this explicitly provisioned, bounded GPU integration downloads assets.
    monkeypatch.setenv("GG_TEST_ALLOW_ONLY_LOOPBACK_HTTP", "0")
    document, root, mapping, seen = qualification
    fixture = document["fixtures"][profile]
    species = fixture["species"]
    work = root / species
    work.mkdir()
    choice_path = work / "model_selection.json"
    choice = selection.select(argparse.Namespace(species=species, genome=fixture["genome"],
        output=str(choice_path), model_profile="auto", taxid_map=str(mapping), genome_metadata="",
        taxonomy_db=document["taxonomy_db"]))
    assert choice["model_profile"] == profile
    prepared_path = work / "runtime.json"
    prepared = cache.prepare(root / "cache", profile, document["runtime_locator"], receipt_path=prepared_path)
    if seen:
        assert prepared["generation"] == next(iter(seen.values()))
        assert prepared["reused"] is True
    else:
        assert prepared["reused"] is False
    seen[profile] = prepared["generation"]
    args = argparse.Namespace(receipt=str(prepared_path), genome=fixture["genome"],
        output=str(work / "cold.gff3"), output_receipt=str(work / "cold.json"), model_profile=profile,
        threads=4, batch_size=32, sequence_threshold=25, selection_receipt=str(choice_path), species=species)
    cold = runtime.run(args)
    before = prepared_path.read_bytes(), prepared_path.stat().st_mtime_ns
    def no_network(*a, **k):
        pytest.fail("warm offline reuse accessed the network")
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(cache, "resolve_manifest", no_network)
        transport = importlib.import_module("annevo_runtime_download")
        patch.setattr(transport, "download_url_to_file", no_network)
        patch.setattr(transport, "limited_urlopen", no_network)
        warm = cache.prepare(root / "cache", profile, document["runtime_locator"],
                             receipt_path=prepared_path, no_download=True)
        assert warm["generation"] == prepared["generation"]
        assert (prepared_path.read_bytes(), prepared_path.stat().st_mtime_ns) == before
        args.output, args.output_receipt = str(work / "warm.gff3"), str(work / "warm.json")
        hot = runtime.run(args)
    assert cold["GFF_sha256"] == hot["GFF_sha256"]
    assert cold["gpu_probe"]["gpu_verified"] and cold["model_selection"] == choice
    assert cold["inference_settings"]["lineage"] == profile
    assert cold["inference_settings"]["overlap_pred"] == (profile in {"Mammalia", "Actinopteri"})
    assert set(p.name for p in (root / "cache/models").glob("*/*.pt")) == {"ANNEVO_" + p + ".pt" for p in seen}
    raw = work / "input" / species
    raw.mkdir(parents=True)
    shutil.copyfile(fixture["genome"], raw / "genome.fa")
    shutil.copyfile(work / "cold.gff3", raw / "annotation.gff3")
    common = ["--species-cds-dir", work / "cds", "--species-gff-dir", work / "gff", "--species-genome-dir", work / "genome"]
    summary = work / "species.tsv"
    invoke(work, "format", "format_species_inputs.py", "--provider", "local", "--input-dir", work / "input",
           *common, "--species-summary-output", summary, "--gene-grouping-mode", "strict",
           "--gff-repair-mode", "safe", "--stats-output", work / "format.json")
    invoke(work, "mapping", "validate_cds_gff_mapping.py", *common, "--species-summary", summary,
           "--nthreads", 1, "--stats-output", work / "mapping.json")
    invoke(work, "longest", "validate_longest_cds_selection.py", "--species-cds-dir", work / "cds",
           "--species-summary", summary, "--summary-outputs-only", "--nthreads", 1, "--stats-output", work / "longest.json")
    proof = json.loads((work / "mapping.json").read_text())
    assert proof["phase_conflicts_total"] == proof["fallback_ids_total"] == 0
    assert proof["mapped_ids_total"] == proof["cds_ids_total"] > 0
    assert proof["species_passed"] == proof["species_checked"] == 1
    longest = json.loads((work / "longest.json").read_text())
    assert longest["species_passed"] == longest["species_checked"] == 1 and longest["genes_total"] > 0

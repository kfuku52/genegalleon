"""Explicit profiles and genome-index reuse preserve source/receipt contracts."""
import json
from pathlib import Path

import pytest

from workflow.tests.test_gene_model_refinement import refinement, tiny_inputs

profiles = __import__("gene_model_species_profiles")
catalog = __import__("gene_model_catalog")


def test_profile_overrides_only_named_species_and_freezes_file(tmp_path):
    inputs, edges, _ = tiny_inputs(tmp_path)
    path = tmp_path / "profiles.tsv"
    path.write_text("species\tmax_intron\tminimum_coverage\tmin_support\nSpecies_target\t100000\t0.98\t3\n")
    root = tmp_path / "run"
    plan = refinement.plan(root, inputs=inputs, edges=edges, mode="off", species_profiles=path)
    target = profiles.parameters_for(plan["request"], "Species_target")
    assert (target["max_intron"], target["minimum_coverage"], target["min_support"]) == (100000, .98, 3)
    assert profiles.parameters_for(plan["request"], "Species_donor1")["max_intron"] == 20000
    path.write_text(path.read_text().replace("100000", "200000"))
    with pytest.raises(ValueError, match="changed"):
        refinement.load(root)


@pytest.mark.parametrize("contents", [
    "species\tunsupported\nSpecies_target\t1\n",
    "species\tmax_intron\nUnknown_species\t10\n",
    "species\tmax_intron\nSpecies_target\t0\n",
    "species\tminimum_identity\nSpecies_target\tnan\n",
    "species\tminimum_coverage\nSpecies_target\t1.1\n",
    "species\tpadding\nSpecies_target\t-1\n",
    "species\tmin_support\nSpecies_target\t1.5\n",
    "species\tmax_intron\nSpecies_target\t10\nSpecies_target\t20\n",
])
def test_invalid_profile_is_rejected(tmp_path, contents):
    path = tmp_path / "profile.tsv"
    path.write_text(contents)
    with pytest.raises(ValueError):
        profiles.read_profiles(path, {"Species_target"})


def test_cached_compressed_genome_corruption_is_rebuilt_and_source_is_unchanged(tmp_path, monkeypatch):
    import gzip
    genome = tmp_path / "source.fa.gz"
    with gzip.open(genome, "wt") as handle:
        handle.write(">chr1\nATGAAATAA\n")
    original = genome.read_bytes()
    cache = tmp_path / "cache"
    monkeypatch.setenv("GG_GENOME_INDEX_CACHE", str(cache))
    with catalog.indexed_genome(genome) as indexed:
        assert indexed.fetch("chr1") == "ATGAAATAA"
    staged = next(cache.glob("*/indexed.fa"))
    receipt = json.loads((staged.parent / "receipt.json").read_text())
    assert "indexed.fa" in receipt["files"]
    staged.write_text(">chr1\nATGCCCTAA\n")
    with catalog.indexed_genome(genome) as indexed:
        assert indexed.fetch("chr1") == "ATGAAATAA"
    assert genome.read_bytes() == original
    assert not Path(str(genome) + ".fai").exists()


def test_zero_prediction_windows_skip_predictor_genome_index(tmp_path, monkeypatch):
    inputs, edges, _ = tiny_inputs(tmp_path)
    root = tmp_path / "run"
    # Plan provenance requires a predictor identity even when no work is queued.
    # This unit test must not depend on the host installing a real predictor.
    predictor = tmp_path / "miniprot"
    predictor.write_text("#!/bin/sh\nexit 99\n")
    predictor.chmod(0o755)
    original_which = refinement.shutil.which
    monkeypatch.setattr(refinement.shutil, "which",
                        lambda name: str(predictor) if name == "miniprot" else original_which(name))
    value = refinement.plan(root, inputs=inputs, edges=edges)
    def forbidden(_path):
        raise AssertionError("Prediction indexed a genome without a candidate window")
    monkeypatch.setattr(refinement, "indexed_genome", forbidden)
    refinement.predict_species(root, value, "Species_target")
    assert json.loads((root / "predictions/Species_target/summary.json").read_text())["windows"] == 0


def test_cached_genome_mutation_during_use_is_rejected(tmp_path, monkeypatch):
    import gzip
    source = tmp_path / "source.fa.gz"
    with gzip.open(source, "wt") as handle:
        handle.write(">chr\nATGAAATAA\n")
    monkeypatch.setenv("GG_GENOME_INDEX_CACHE", str(tmp_path / "cache"))
    with pytest.raises(OSError, match="cache changed"):
        with catalog.indexed_genome(source) as genome:
            Path(genome.filename.decode()).write_text(">chr\nATGCCCTAA\n")

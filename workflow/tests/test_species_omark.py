import json
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
import species_omark as qc
from test_genome_evolution_protein_mode import _prepare_stub_binaries


def test_omark_cache_uses_query_db_and_tool_identity(tmp_path, monkeypatch):
    tools = _prepare_stub_binaries(tmp_path)
    monkeypatch.setenv("PATH", str(tools) + os.pathsep + os.environ["PATH"])
    database = tmp_path / "LUCA.h5"
    database.write_bytes(b"test database")
    # Ensure the fixture filesystem can distinguish two same-size writes.
    os.utime(database, ns=(1_000_000_100, 1_000_000_100))
    name = "Plant_species"
    proteins = tmp_path / "input.fa"
    proteins.write_text(">g2 first description\nMI*\n>g1\nMA\n")
    captured = tmp_path / "capture/omamer_queries.txt"
    def runner():
        return qc.SpeciesOMArk(tmp_path / "cache", database, tmp_path / "scratch")
    first = runner()
    key = first.protein(name, proteins)
    sums = first.result(name, key) / "omark" / (name + ".sum")
    with monkeypatch.context() as patch:
        original_open = Path.open
        def no_db_read(path, *args, **kwargs):
            if path == database:
                pytest.fail("Reread the unchanged database while reusing a result")
            return original_open(path, *args, **kwargs)
        patch.setattr(Path, "open", no_db_read)
        proteins.write_text(">g1 new description\nm\na\n>g2\nMI\n")
        assert runner().protein(name, proteins) == key
        assert len(captured.read_text().splitlines()) == 1
    proteins.write_text(">g1\nMA\n>renamed\nMI\n")
    assert runner().protein(name, proteins) != key
    database.write_bytes(b"next database")
    next_db = runner()
    next_key = next_db.protein(name, proteins)
    with (tools / "omark").open("a") as tool:
        tool.write("# changed wrapper\n")
    next_tool = runner()
    assert next_tool.protein(name, proteins) == next_key
    assert next_tool.result(name, next_key) != next_db.result(name, next_key)
    assert len(captured.read_text().splitlines()) == 4
    assert sums.is_file()


def test_cached_baseline_does_not_reopen_catalog_and_handles_no_proteins(tmp_path, monkeypatch):
    tools = _prepare_stub_binaries(tmp_path)
    monkeypatch.setenv("PATH", str(tools) + os.pathsep + os.environ["PATH"])
    database = tmp_path / "db"
    database.write_text("db")
    catalog = tmp_path / "loci.jsonl"
    name = "Plant_assembly"
    candidate = dict(candidate_id="t", protein="ML", cds="ATGCTT", quality={"usable": False})
    catalog.write_text(json.dumps(dict(gene_id="g", candidates=[candidate])) + "\n")
    runner = qc.SpeciesOMArk(tmp_path / "cache", database, tmp_path / "scratch")
    assert qc.baseline_keys(runner, name, catalog, None, True) == {"input": qc.EMPTY_QUERY}
    open_path = Path.open
    def no_catalog_read(path, *args, **kwargs):
        if path == catalog:
            pytest.fail("Reopened a cached catalogue")
        return open_path(path, *args, **kwargs)
    monkeypatch.setattr(Path, "open", no_catalog_read)
    assert qc.baseline_keys(runner, name, catalog, None, True) == {"input": qc.EMPTY_QUERY}
    assert qc.comparison_fields(None)["missing_hogs"] == ""
    assert not (tmp_path / "capture/omamer_queries.txt").exists()



def test_omark_reuses_first_taxonomy_generation_and_detects_replacement(tmp_path, monkeypatch):
    tools = _prepare_stub_binaries(tmp_path)
    monkeypatch.setenv("PATH", str(tools) + os.pathsep + os.environ["PATH"])
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    with (tools / "omark").open("a") as script:
        script.write("\nfrom pathlib import Path\np = Path.home() / '.etetoolkit/taxa.sqlite'\n"
                     "if not p.exists():\n    p.parent.mkdir(parents=True)\n    p.write_text('first')\n")
    database = tmp_path / "db"
    database.write_text("fixture")
    def runner():
        return qc.SpeciesOMArk(tmp_path / "cache", database, tmp_path / "scratch")
    cold = runner()
    waiting = runner()
    data = b">g\nML\n"
    key = cold.evaluate("Plant_target", data)
    assert cold.taxonomy_sha256
    assert waiting.evaluate("Plant_target", data) == key and waiting.searches == 0
    warm = runner()
    assert warm.contract == cold.contract
    assert warm.evaluate("Plant_target", data) == key and warm.searches == 0
    warm.taxonomy.write_text("next!")
    updated = runner()
    assert updated.contract != cold.contract
    assert updated.evaluate("Plant_target", data) == key and updated.searches == 1


def test_excluded_rescue_locus_adds_no_protein(tmp_path, monkeypatch):
    tools = _prepare_stub_binaries(tmp_path)
    monkeypatch.setenv("PATH", str(tools) + os.pathsep + os.environ["PATH"])
    database = tmp_path / "db"
    database.write_text("fixture")
    runner = qc.SpeciesOMArk(tmp_path / "cache", database, tmp_path / "scratch")
    cds = tmp_path / "augmented.fa"
    cds.write_text(">original\nATGAAA\n>added\nATGTAA\n")
    augmentation = {"species": "Plant_target", "cds": str(cds), "original_cds": "1", "rescued_models": "1"}
    catalog = tmp_path / "loci.jsonl"
    candidate = dict(candidate_id="t", protein="MK", cds="ATGAAA", blocks=[[1, 6, 0]], quality={"usable": True})
    catalog.write_text(json.dumps(dict(gene_id="g", candidates=[candidate])) + "\n")
    with pytest.raises(ValueError, match="Unmapped newly rescued CDS"):
        qc.baseline_keys(runner, "Plant_target", catalog, augmentation, False)
    (tmp_path / "excluded_loci.jsonl").write_text(json.dumps(dict(source_gene_id="added", candidates=[])) + "\n")
    keys = qc.baseline_keys(runner, "Plant_target", catalog, augmentation, False)
    assert keys["input"] == keys["rescued"] and runner.searches == 1



def test_shared_ete_database_is_part_of_the_run_identity(tmp_path, monkeypatch):
    binaries = _prepare_stub_binaries(tmp_path)
    monkeypatch.setenv("PATH", str(binaries) + os.pathsep + os.environ["PATH"])
    db = tmp_path / "db"
    db.write_bytes(b"OMAmer")
    taxonomy = tmp_path / "shared.sqlite"
    taxonomy.write_bytes(b"taxonomy")
    monkeypatch.setenv("GG_TAXONOMY_DBFILE", str(taxonomy))
    runner = qc.SpeciesOMArk(tmp_path / "cache", db, tmp_path / "scratch")
    assert runner.taxonomy == taxonomy and runner.taxonomy_sha256
    assert runner.evaluate("Plant_species", b">g\nMA\n") != qc.EMPTY_QUERY

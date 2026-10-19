"""Taxonomy-backed model choices, rejected evidence and per-species provenance."""
import argparse
import importlib
import json
import sqlite3
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))
selection = importlib.import_module("annevo_model_selection")
runtime = importlib.import_module("annevo_runtime")


def taxonomy(path):
    conn = sqlite3.connect(path)
    conn.executescript("CREATE TABLE species(taxid INTEGER PRIMARY KEY,spname TEXT,rank TEXT,track TEXT);"
                       "CREATE TABLE synonym(taxid INTEGER,spname TEXT);"
                       "CREATE TABLE merged(taxid_old INTEGER,taxid_new INTEGER);")
    conn.execute("INSERT INTO species VALUES(1,'root','no rank','1')")
    for index, (profile, anchor) in enumerate(selection.PROFILES.items()):
        conn.execute("INSERT INTO species VALUES(?,?,?,?)", (anchor, profile, "class", f"{anchor},1"))
        conn.execute("INSERT INTO species VALUES(?,?,?,?)", (1000000 + index, f"{profile} species", "species", f"{1000000 + index},{anchor},1"))
    conn.execute("INSERT INTO species VALUES(2000000,'Oryza sativa','species','2000000,3398,1')")
    conn.execute("INSERT INTO species VALUES(2000001,'Other species','species','2000001,1')")
    conn.execute("INSERT INTO merged VALUES(999999,1000000)")
    conn.commit()
    conn.close()


def fixture(tmp_path, taxid="1000000"):
    db, genome, mapping = [tmp_path / n for n in ("taxa.sqlite", "genome.fa", "taxids.tsv")]
    taxonomy(db)
    genome.write_text(">s\nATGAAATAA\n")
    mapping.write_text("species\ttaxid\nTest_species\t" + taxid + "\n")
    return argparse.Namespace(species="Test_species", genome=str(genome), model_profile="auto",
                              taxonomy_db=str(db), taxid_map=str(mapping), genome_metadata="",
                              output=str(tmp_path / "selection.json"))


@pytest.mark.parametrize("profile,taxid", [(p, str(1000000 + i)) for i, p in enumerate(selection.PROFILES)] + [("Magnoliopsida", "2000000")])
def test_all_clades_and_monocots_follow_ncbi_ancestors(profile, taxid, tmp_path):
    args = fixture(tmp_path, taxid)
    record = selection.select(args)
    assert record["model_profile"] == profile
    assert record["anchor_taxid"] == selection.PROFILES[profile]
    assert record["overlap_pred"] == (profile in {"Mammalia", "Actinopteri"})
    assert record["taxonomy"]["sha256"] == selection.file_sha(args.taxonomy_db)
    assert record["evidence"][0]["sha256"] == selection.file_sha(args.taxid_map)
    before = Path(args.output).read_bytes(), Path(args.output).stat().st_mtime_ns
    selection.select(args)
    assert (Path(args.output).read_bytes(), Path(args.output).stat().st_mtime_ns) == before


@pytest.mark.parametrize("taxid,message", [("2000001", "No supported"), ("123456", "Unresolved"), ("40674", "does not identify a species")])
def test_unsupported_and_unresolved_fail_with_a_record(taxid, message, tmp_path):
    args = fixture(tmp_path, taxid)
    with pytest.raises(ValueError, match=message):
        selection.select(args)
    assert not Path(args.output).exists()
    failed = json.loads(Path(args.output + ".failed.json").read_text())
    assert failed["status"] == "error" and message in failed["reason"]


def test_merged_taxid_is_recorded_and_equivalent_evidence_agrees(tmp_path):
    args = fixture(tmp_path, "999999")
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("species\ttaxid\tgenome_sha256\nTest_species\t1000000\t" + selection.file_sha(args.genome) + "\n")
    args.genome_metadata = str(metadata)
    record = selection.select(args)
    assert record["taxid"] == 1000000 and record["merged_taxids"] == [999999]


@pytest.mark.parametrize("fault,message", [("conflict", "Conflicting"), ("hash", "checksum"), ("name", "contradictory scientific_name"), ("ambiguous", "Ambiguous"), ("duplicate", "duplicate species"), ("track", "contradictory NCBI")])
def test_bad_metadata_cannot_replace_previous_selection(fault, message, tmp_path):
    args = fixture(tmp_path)
    selection.select(args)
    previous = Path(args.output).read_bytes()
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("species\ttaxid\tgenome_sha256\tscientific_name\nTest_species\t" +
                        ("1000001" if fault == "conflict" else "1000000") + "\t" +
                        ("a" * 64 if fault == "hash" else selection.file_sha(args.genome)) + "\t" +
                        ("Other species" if fault == "name" else "Mammalia species") + "\n")
    args.genome_metadata = str(metadata)
    if fault == "duplicate":
        Path(args.taxid_map).write_text("species\ttaxid\nTest_species\t1000000\nTest_species\t1000000\n")
    if fault in {"ambiguous", "track"}:
        with sqlite3.connect(args.taxonomy_db) as conn:
            if fault == "ambiguous":
                conn.execute("INSERT INTO synonym VALUES(2000001,'Mammalia species')")
            else:
                conn.execute("UPDATE species SET track='40674,8782,1' WHERE taxid=40674")
    with pytest.raises(ValueError, match=message):
        selection.select(args)
    assert Path(args.output).read_bytes() == previous


def test_unmapped_and_qualified_labels_are_not_guessed(tmp_path):
    args = fixture(tmp_path)
    args.species = "Mammalia_cf_species"
    with pytest.raises(ValueError, match="species labels are not guessed"):
        selection.select(args)


def test_different_species_selections_do_not_collide(tmp_path):
    args = fixture(tmp_path)
    Path(args.taxid_map).write_text("species\ttaxid\nOne_species\t1000000\nTwo_species\t1000004\n")
    def choose(species):
        local = argparse.Namespace(**{**vars(args), "species": species, "output": str(tmp_path / species / "selection.json")})
        return selection.select(local)
    with ThreadPoolExecutor(2) as pool:
        records = list(pool.map(choose, ["One_species", "Two_species"]))
    assert [r["model_profile"] for r in records] == ["Mammalia", "Fungi"]
    assert len(list(tmp_path.glob("*/selection.json"))) == 2


@pytest.mark.parametrize("profile", selection.PROFILES)
def test_native_settings_are_read_from_the_selected_source(profile, tmp_path):
    source = tmp_path / "prediction.py"
    source.write_text("def main():\n length_config=" + repr({p: [30720, 5120, 32] for p in selection.PROFILES}) + "\n")
    settings = runtime.inference_settings({"native_entrypoint": str(tmp_path / "annotation.py")}, profile)
    assert settings["lineage"] == profile and settings["window_size"] == 30720
    assert settings["prediction_source_sha256"] == selection.file_sha(source)


def test_native_selection_requires_same_species_genome_and_model(tmp_path):
    args = fixture(tmp_path)
    record = selection.select(args)
    run_args = argparse.Namespace(selection_receipt=args.output, species=args.species, model_profile="Mammalia")
    assert runtime.selection_record(run_args, args.genome)[0] == record
    run_args.species = "Other_species"
    with pytest.raises(ValueError, match="differs"):
        runtime.selection_record(run_args, args.genome)


@pytest.mark.parametrize("boundary", ["hash", "lineage"])
def test_wal_writer_cannot_publish_a_selection_from_unhashed_bytes(tmp_path, monkeypatch, boundary):
    args = fixture(tmp_path)
    db = Path(args.taxonomy_db)
    conn = sqlite3.connect(db)
    assert conn.execute("PRAGMA journal_mode=WAL").fetchone()[0] == "wal"
    conn.close()
    selection.select(args)
    previous = Path(args.output).read_bytes()
    before = selection.stable_signature(db)
    writers = []

    def update():
        if not writers:
            writer = sqlite3.connect(db)
            writer.execute("UPDATE species SET track='1000000,4751,1' WHERE taxid=1000000")
            writer.commit()
            writers.append(writer)

    if boundary == "hash":
        real_sha = selection.file_sha
        def hash_then_update(path):
            digest = real_sha(path)
            if Path(path) == db:
                update()
            return digest
        monkeypatch.setattr(selection, "file_sha", hash_then_update)
    else:
        real_node = selection.Taxonomy.node
        def node_then_update(self, taxid):
            node = real_node(self, taxid)
            update()
            return node
        monkeypatch.setattr(selection.Taxonomy, "node", node_then_update)
    try:
        with pytest.raises(ValueError, match="active transaction"):
            selection.select(args)
        assert selection.stable_signature(db) == before
        assert Path(args.output).read_bytes() == previous
        assert json.loads(Path(args.output + ".failed.json").read_text())["status"] == "error"
    finally:
        for writer in writers:
            writer.close()


def test_wal_change_after_selection_is_rejected_by_prediction(tmp_path):
    args = fixture(tmp_path)
    db = Path(args.taxonomy_db)
    writer = sqlite3.connect(db)
    writer.execute("PRAGMA journal_mode=WAL")
    writer.close()
    selection.select(args)
    writer = sqlite3.connect(db)
    try:
        writer.execute("INSERT INTO synonym VALUES(1000000,'late synonym')")
        writer.commit()
        run_args = argparse.Namespace(selection_receipt=args.output, species=args.species, model_profile="Mammalia")
        with pytest.raises(ValueError, match="active transaction"):
            runtime.selection_record(run_args, args.genome)
    finally:
        writer.close()


@pytest.mark.parametrize("role", ["genome", "taxid_map", "genome_metadata", "taxonomy_db"])
@pytest.mark.parametrize("publication", ["selected", "failed"])
def test_selector_output_cannot_replace_any_input(tmp_path, role, publication):
    args = fixture(tmp_path)
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("species\ttaxid\tgenome_sha256\nTest_species\t1000000\t" + selection.file_sha(args.genome) + "\n")
    args.genome_metadata = str(metadata)
    if publication == "failed":
        renamed = tmp_path / "protected.failed.json"
        Path(getattr(args, role)).rename(renamed)
        setattr(args, role, str(renamed))
        args.output = str(renamed).removesuffix(".failed.json")
    else:
        args.output = getattr(args, role)
    inputs = [Path(getattr(args, name)) for name in ("genome", "taxid_map", "genome_metadata", "taxonomy_db")]
    before = {p: p.read_bytes() for p in inputs}
    with pytest.raises(ValueError, match="outputs must not overwrite inputs"):
        selection.select(args)
    assert {p: p.read_bytes() for p in inputs} == before

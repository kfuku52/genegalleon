import copy
import json
import sqlite3
import weakref

import pytest

from workflow.support import gene_model_store as store
from workflow.support.gene_model_selection import select_representatives
from workflow.support.gene_model_store import build_store, iter_loci, iter_locus_keys, load_locus, select_from_store
from workflow.tests.test_gene_model_selection import (
    candidate,
    catalog,
    edge,
    extension_fixture,
    invalid_donor_fixture,
    invalid_length_normalization_fixture,
    protein,
    weighted_irregular_fixture,
)


def write_catalogs(tmp_path, catalogs):
    directories = []
    for data in catalogs:
        directory = tmp_path / data["species"]
        directory.mkdir()
        metadata = {key: value for key, value in data.items() if key != "loci"}
        (directory / "catalog_metadata.json").write_text(json.dumps(metadata) + "\n")
        (directory / "loci.jsonl").write_text("".join(json.dumps(locus) + "\n" for locus in data["loci"]))
        directories.append(directory)
    return directories


def scientific_output(result):
    return {key: value for key, value in result.items() if key != "metrics"}


def test_bulky_metadata_preserves_catalog_and_foreign_key_integrity(tmp_path):
    catalogs, _ = extension_fixture()
    catalogs[0]["fasta_mapping"] = [{"source_cds": "A" * 16384} for _ in range(128)]
    database = tmp_path / "models.sqlite"
    build_store(write_catalogs(tmp_path, catalogs), database)
    with sqlite3.connect(database) as connection:
        connection.execute("PRAGMA foreign_keys=ON")
        metadata = json.loads(connection.execute("SELECT json FROM catalogs WHERE species='A'").fetchone()[0])
        assert metadata["fasta_mapping"] == catalogs[0]["fasta_mapping"]
        assert connection.execute("PRAGMA foreign_key_check").fetchall() == []
        with pytest.raises(sqlite3.IntegrityError):
            connection.execute("INSERT INTO loci VALUES ('absent','g','chr1','+','{}')")
        with pytest.raises(sqlite3.IntegrityError):
            connection.execute("INSERT INTO candidate_owners VALUES ('A','absent','missing_locus')")
    assert list(iter_loci(database, "A")) == catalogs[0]["loci"]


def test_store_streams_complete_loci_and_species_keys_in_stable_order(tmp_path):
    catalogs, _ = extension_fixture()
    catalogs.append(catalog("Z", ("g1", [candidate("Z1", protein())]),
                            ("g2", [candidate("Z2", protein(2))])))
    directories = write_catalogs(tmp_path, catalogs)
    database = tmp_path / "models.sqlite"
    receipt = build_store(reversed(directories), database)
    assert receipt["species_count"] == 4
    assert receipt["loci_count"] == 5
    assert receipt["candidate_count"] == 6
    assert len(receipt["sources"][0]["loci_sha256"]) == 64
    assert list(iter_locus_keys(database, "Z")) == [("Z", "g1"), ("Z", "g2")]
    assert list(iter_locus_keys(database)) == [("A", "g"), ("B", "g"), ("C", "g"), ("Z", "g1"), ("Z", "g2")]
    assert load_locus(database, "A", "g") == catalogs[0]["loci"][0]
    assert list(iter_loci(database, "Z")) == catalogs[3]["loci"]
    assert len(list(iter_loci(database))) == 5
    with pytest.raises(KeyError):
        load_locus(database, "Z", "absent")


@pytest.mark.parametrize("policy", ["conserved", "longest"])
def test_component_selection_equivalent_to_memory_engine_with_wgd_and_isolates(tmp_path, policy):
    catalogs, edges = extension_fixture()
    sequence = protein(2)
    for item in catalogs:
        species = item["species"]
        item["loci"].append({"species": species, "gene_id": "wgd2", "candidates": [candidate(species + "_copy2", sequence)]})
    catalogs.append(catalog("D", ("unlinked", [candidate("D_long", sequence + "AAAA"),
                                             candidate("D_source", sequence)])))
    catalogs[-1]["loci"][0]["source_baseline_candidate_id"] = "D_source"
    edges += [edge("A", "B", "wgd2", "wgd2"), edge("A", "C", "wgd2", "wgd2")]
    expected = select_representatives(catalogs, edges, policy=policy, exact_limit=100)
    database = tmp_path / "models.sqlite"
    build_store(write_catalogs(tmp_path, catalogs), database)
    result = select_from_store(database, iter(edges), policy=policy, exact_limit=100)
    assert scientific_output(result) == scientific_output(expected)
    assert result["metrics"]["peak_loaded_loci"] == 3
    assert result["metrics"]["store_loci"] == 7
    assert result["metrics"]["limited_components"] == 0


def test_ambiguity_remains_abstention_and_giant_components_are_bounded(tmp_path):
    catalogs, edges = extension_fixture()
    catalogs[1]["loci"].append({"species": "B", "gene_id": "copy2", "candidates": [candidate("B_copy2", protein(2))]})
    edges.append(edge("A", "B", "g", "copy2"))
    database = tmp_path / "models.sqlite"
    build_store(write_catalogs(tmp_path, catalogs), database)
    result = select_from_store(database, edges)
    expected = select_representatives(catalogs, edges)
    assert scientific_output(result) == scientific_output(expected)
    bounded = select_from_store(database, edges, max_component_loci=2)
    # Global one-to-many ambiguity removes the connecting votes before the
    # memory cap, without turning the target into an apparently unique copy.
    assert scientific_output(bounded) == scientific_output(expected)
    assert bounded["metrics"]["limited_components"] == 0
    assert bounded["metrics"]["peak_loaded_loci"] == 1
    assert len(bounded["selections"]) == 4
    assert bounded["selections"][0]["reason"] == "unresolved_paralog_copy"
    assert bounded["selections"][0]["candidate_id"] == "A_long"


def test_giant_usable_component_still_abstains_without_loading_its_sequences_together(tmp_path):
    catalogs, edges = extension_fixture()
    database = tmp_path / "models.sqlite"
    build_store(write_catalogs(tmp_path, catalogs), database)
    result = select_from_store(database, edges, max_component_loci=2)
    assert result["metrics"]["limited_components"] == 1
    assert result["metrics"]["peak_loaded_loci"] == 1
    assert len(result["selections"]) == 3
    assert all(row["reason"] == "component_limit_exceeded" for row in result["selections"])
    assert result["selections"][0]["candidate_id"] == "A_long"


def test_sqlite_overlay_visible_without_reloading_other_species(tmp_path):
    catalogs, edges = extension_fixture()
    database = tmp_path / "models.sqlite"
    build_store(write_catalogs(tmp_path, catalogs), database)
    locus = load_locus(database, "A", "g")
    locus["candidates"][1]["quality"]["full_length_supported"] = True
    with sqlite3.connect(database) as connection:
        connection.execute("UPDATE loci SET json=? WHERE species=? AND gene_id=?", (json.dumps(locus), "A", "g"))
        assert load_locus(connection, "A", "g") == locus
        assert list(iter_locus_keys(connection, "B")) == [("B", "g")]
    assert select_from_store(database, edges)["selections"][0]["candidate_id"] == "A_conserved"


def test_owned_reader_keeps_a_consistent_snapshot_and_releases_it(tmp_path):
    catalogs, _ = extension_fixture()
    database = tmp_path / "models.sqlite"
    build_store(write_catalogs(tmp_path, catalogs), database)
    original = load_locus(database, "A", "g")
    changed = copy.deepcopy(original)
    changed["candidates"][1]["quality"]["full_length_supported"] = True
    with sqlite3.connect(database) as writer:
        writer.execute("PRAGMA journal_mode=WAL")
        with store._connection(database) as reader:
            assert load_locus(reader, "A", "g") == original
            writer.execute("UPDATE loci SET json=? WHERE species=? AND gene_id=?", (json.dumps(changed), "A", "g"))
            writer.commit()
            assert load_locus(reader, "A", "g") == original
            with pytest.raises(sqlite3.OperationalError, match="readonly"):
                reader.execute("DELETE FROM loci")
    assert load_locus(database, "A", "g") == changed


@pytest.mark.parametrize("invalid", ["duplicate_locus", "duplicate_candidate", "wrong_species", "schema", "empty"])
def test_atomic_import_rejects_invalid_contract_and_preserves_previous_database(tmp_path, invalid):
    catalogs, _ = extension_fixture()
    directories = write_catalogs(tmp_path, catalogs)
    database = tmp_path / "models.sqlite"
    build_store(directories, database)
    previous = database.read_bytes()
    path = directories[0] / "loci.jsonl"
    locus = copy.deepcopy(catalogs[0]["loci"][0])
    if invalid == "duplicate_locus":
        path.write_text(json.dumps(locus) + "\n" + json.dumps(locus) + "\n")
    elif invalid == "duplicate_candidate":
        other = dict(locus, gene_id="other")
        path.write_text(json.dumps(locus) + "\n" + json.dumps(other) + "\n")
    elif invalid == "wrong_species":
        locus["species"] = "wrong"
        path.write_text(json.dumps(locus) + "\n")
    elif invalid == "schema":
        (directories[0] / "catalog_metadata.json").write_text(json.dumps({"schema": 2, "species": "A"}))
    else:
        path.write_text("\n")
    with pytest.raises(ValueError):
        build_store(directories, database)
    assert database.read_bytes() == previous
    assert not list(tmp_path.glob(".locus-store-*"))


def test_store_requires_streaming_artifacts_and_refuses_incompatible_store(tmp_path):
    directory = tmp_path / "A"
    directory.mkdir()
    (directory / "catalog.json").write_text(json.dumps(extension_fixture()[0][0]))
    with pytest.raises(FileNotFoundError):
        build_store([directory], tmp_path / "models.sqlite")
    database = tmp_path / "unknown.sqlite"
    sqlite3.connect(database).close()
    with pytest.raises(ValueError, match="schema"):
        list(iter_loci(database))


def test_metadata_counts_and_explicit_candidate_ownership_are_checked(tmp_path):
    catalogs, _ = extension_fixture()
    directories = write_catalogs(tmp_path, catalogs)
    metadata_path = directories[0] / "catalog_metadata.json"
    metadata_path.write_text(json.dumps({"schema_version": 1, "species": "A", "summary": {"loci": 99}}))
    with pytest.raises(ValueError, match="summary"):
        build_store(directories, tmp_path / "models.sqlite")
    metadata_path.write_text(json.dumps({"schema_version": 1, "species": "A"}))
    locus = catalogs[0]["loci"][0]
    locus["candidates"][0]["gene_id"] = "different_owning_gene"
    (directories[0] / "loci.jsonl").write_text(json.dumps(locus) + "\n")
    with pytest.raises(ValueError, match="ownership"):
        build_store(directories, tmp_path / "models.sqlite")


def test_empty_store_and_absent_edge_are_explicit(tmp_path):
    database = tmp_path / "models.sqlite"
    build_store([], database)
    assert select_from_store(database, [])["selections"] == []
    assert select_from_store(database, [])["metrics"]["peak_loaded_loci"] == 0
    with pytest.raises(ValueError, match="absent"):
        select_from_store(database, [edge("A", "B")])


def test_store_excluded_donors_cannot_change_selection_or_margin(tmp_path):
    catalogs, edges, donors, donor_edges = invalid_donor_fixture({"sequence_mismatch": True})
    expected = select_representatives(catalogs, edges, exact_limit=10)
    database = tmp_path / "models.sqlite"
    build_store(write_catalogs(tmp_path, catalogs + donors), database)
    result = select_from_store(database, edges + donor_edges, exact_limit=10)
    memory = select_representatives(catalogs + donors, edges + donor_edges, exact_limit=10)
    assert result["selections"][:3] == expected["selections"]
    assert result["scores"] == expected["scores"]
    assert result["audit"] == expected["audit"]
    assert scientific_output(result) == scientific_output(memory)


def test_store_ineligible_long_original_uses_same_eligible_quality_scale(tmp_path):
    catalogs, edges = invalid_length_normalization_fixture({"usable": False, "sequence_mismatch": True})
    expected_catalogs, _ = invalid_length_normalization_fixture()
    expected = select_representatives(expected_catalogs, edges, exact_limit=2)
    for item in catalogs:
        item["loci"][0]["candidates"].reverse()
    database = tmp_path / "models.sqlite"
    build_store(write_catalogs(tmp_path, list(reversed(catalogs))), database)
    memory = select_representatives(catalogs, edges, exact_limit=2)
    stored = select_from_store(database, list(reversed(edges)), exact_limit=2)
    assert scientific_output(stored) == scientific_output(memory)
    assert stored["selections"] == expected["selections"]
    assert stored["audit"] == expected["audit"]
    assert [row for row in stored["scores"] if row["candidate_id"] != "A_invalid"] == expected["scores"]


def test_store_weighted_irregular_graph_matches_exact_memory_audit_bytes(tmp_path):
    catalogs, edges = weighted_irregular_fixture()
    database = tmp_path / "models.sqlite"
    build_store(write_catalogs(tmp_path, catalogs), database)
    direct = select_representatives(catalogs, edges, exact_limit=2)
    stored = select_from_store(database, list(reversed(edges)), exact_limit=2)
    assert scientific_output(stored) == scientific_output(direct)
    assert json.dumps(scientific_output(stored), sort_keys=True) == json.dumps(scientific_output(direct), sort_keys=True)


@pytest.mark.parametrize("exclusion", ["unusable", "candidate_limit", "explicit_ambiguity", "coordinates"])
def test_excluded_edges_do_not_push_valid_component_over_cap(tmp_path, exclusion):
    catalogs, edges = extension_fixture()
    expected = select_representatives(catalogs, edges)
    sequence = catalogs[0]["loci"][0]["candidates"][1]["protein"]
    for index in range(7):
        species = f"P{index}"
        records = [candidate(species + "_source", sequence)]
        if exclusion == "unusable":
            records[0]["quality"].update(usable=False, annotated_pseudogene=True)
        if exclusion == "candidate_limit":
            records += [candidate(species + "_alternate" + str(i), sequence) for i in range(32)]
        donor = catalog(species, ("g", records))
        donor_edge = edge("A", species)
        if exclusion == "explicit_ambiguity":
            # An ambiguous edge is attached to an isolated target, avoiding
            # any change to the valid A/B/C correspondence while testing cap
            # formation from a long chain of explicitly excluded edges.
            donor_edge = edge("P0", species) if index else None
            if donor_edge:
                donor_edge["ambiguous"] = True
        if exclusion == "coordinates":
            donor["loci"][0]["ambiguous_coordinates"] = True
            donor_edge = edge("P0", species) if index else None
        catalogs.append(donor)
        if donor_edge:
            edges.append(donor_edge)
    database = tmp_path / "models.sqlite"
    build_store(write_catalogs(tmp_path, catalogs), database)
    result = select_from_store(database, edges, max_component_loci=3)
    memory = select_representatives(catalogs, edges)
    assert result["selections"][:3] == expected["selections"]
    assert scientific_output(result) == scientific_output(memory)
    assert result["metrics"]["limited_components"] == 0
    assert result["metrics"]["peak_loaded_loci"] == 3
    assert len(result["selections"]) == len(catalogs)
    assert len(result["omitted_edges"]) == (7 if exclusion in {"unusable", "candidate_limit"} else 6)


def test_streaming_admission_keeps_only_one_locus_before_loading_a_bounded_component(tmp_path, monkeypatch):
    catalogs, edges = extension_fixture()
    for index in range(40):
        species = f"P{index:02d}"
        catalogs.append(catalog(species, ("g", [candidate(species + "_excluded", protein(index + 100),
                                                         quality={"sequence_mismatch": True})])))
        edges.append(edge("A", species))
    database = tmp_path / "models.sqlite"
    build_store(write_catalogs(tmp_path, catalogs), database)
    alive = weakref.WeakSet()
    peak = 0
    class TrackedLocus(dict):
        __hash__ = object.__hash__
    def tracked(locus):
        nonlocal peak
        result = TrackedLocus(locus)
        alive.add(result)
        peak = max(peak, len(alive))
        return result
    original_iter, original_load = store.iter_loci, store.load_locus
    def tracked_iter(*args, **kwargs):
        for locus in original_iter(*args, **kwargs):
            yield tracked(locus)
    def tracked_load(*args, **kwargs):
        return tracked(original_load(*args, **kwargs))
    monkeypatch.setattr(store, "iter_loci", tracked_iter)
    monkeypatch.setattr(store, "load_locus", tracked_load)
    result = store.select_from_store(database, edges, max_component_loci=3)
    assert result["metrics"]["peak_loaded_loci"] == 3
    assert peak <= 3
    assert len(result["selections"]) == 43

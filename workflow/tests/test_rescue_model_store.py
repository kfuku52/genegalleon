"""Lossless compact publications, selective integrity, and legacy interoperability."""
import copy
import gzip
import hashlib
import json
import sqlite3
from pathlib import Path

import pytest

from workflow.support import rescue_model_store as store


def receipt(directory):
    files = {p.relative_to(directory).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
             for p in directory.rglob("*") if p.is_file() and p.name != "receipt.json"}
    (directory / "receipt.json").write_text(json.dumps({"key": {"species": "test"}, "files": files}))


def model(i=0, status="proposal"):
    paf = f"#PAF\tquery{i}\t10\t0\t10\t+\tchr1\t100\t3\t33\t10\t10\t60\n"
    raw = {"seqid": "chr1", "strand": "+", "cds": [[3, 33, 0]], "frameshift": False,
           "coverage": 1.0, "identity": 0.8, "query": f"query{i}", "id": f"old{i}",
           "paf_query": f"query{i}", "paf": paf, "evidence": {"donor": f"species{i}", "query": f"gene{i}"}}
    return {"query": f"query{i}", "seqid": "chr1", "cds": [[3, 36, 0]], "strand": "+",
            "sequence": "ATG" + "AAA" * 9 + "TAA", "coverage": 1.0, "identity": 0.8,
            "id": f"new{i}", "paf": paf, "paf_query": f"query{i}", "raw_prediction": raw,
            "evidence": copy.deepcopy(raw["evidence"]), "status": status, "support": [raw["evidence"]],
            "problems": [], "terminal_completion": {"status": "completed", "query": f"gene{i}"},
            "quality_evidence": {"coding": True}, "partial_evidence": {"partial": i % 2 == 1},
            "unknown": {"null": None, "empty": [], "unicode": "β植物", "number": -0.0}}


def published(tmp_path, records, *, partial=None, revisions=None, codec="gzip"):
    directory = tmp_path / "worker"
    manifest = store.write_model_store(directory, iter(records), partial_models=partial,
                                       revisions=revisions, shard_bytes=900, codec=codec)
    receipt(directory)
    return directory, manifest


@pytest.mark.parametrize("codec", ["gzip", "auto", "zstd"])
def test_exact_all_records_fields_order_and_raw_binding(tmp_path, codec):
    if codec == "zstd":
        pytest.importorskip("zstandard")
    original = [model(0, "accepted"), model(1, "accepted_alternative_path"), model(2), model(3)]
    revisions = [{**model(2), "support": [model(0)["evidence"], model(3)["evidence"]]}]
    directory, manifest = published(tmp_path, original, revisions=revisions, codec=codec)
    restored = list(store.iter_models(directory))
    assert restored == original
    assert json.dumps(restored, ensure_ascii=False) == json.dumps(original, ensure_ascii=False)
    assert list(store.iter_accepted_models(directory)) == [original[0]]
    assert list(store.iter_models(directory, statuses={"accepted_alternative_path"})) == [original[1]]
    assert list(store.iter_partial_models(directory)) == [original[1], original[3]]
    assert list(store.iter_revision_models(directory)) == revisions
    assert manifest["counts"]["bodies"] == 1
    assert manifest["codec"] in {"gzip", "zstd"}
    assert list(store.iter_models(directory / "model_store")) == original


def test_mutating_restored_candidate_does_not_change_other_candidates(tmp_path):
    directory, _ = published(tmp_path, [model(0), model(1)])
    iterator = store.iter_models(directory)
    first = next(iterator)
    first["cds"][0][0] = 999
    first["raw_prediction"]["cds"][0][0] = 888
    first["problems"].append("bad")
    second = next(iterator)
    assert second == model(1)
    assert list(iterator) == []


@pytest.mark.parametrize("prefix", ["##PAF\t", "#PAF\t", ""])
@pytest.mark.parametrize("newline", ["\n", "\r\n", ""])
def test_real_miniprot_paf_aliases_share_body_with_exact_nested_binding(tmp_path, prefix, newline):
    records = [model(0), model(1)]
    for record in records:
        record["paf"] = prefix + record["paf"].removeprefix("#PAF\t").rstrip("\r\n") + newline
        record["raw_prediction"]["paf"] = record["paf"]
    body0, envelope0 = store._pack(records[0])
    body1, envelope1 = store._pack(records[1])
    assert body0 == body1
    assert envelope0 != envelope1
    assert b"query0" not in body0 and b"query1" not in body1
    assert store._unpack(envelope0, body0) == records[0]
    directory, manifest = published(tmp_path, records)
    assert manifest["counts"]["bodies"] == 1
    assert list(store.iter_models(directory)) == records


@pytest.mark.parametrize("kind", ["models", "accepted", "partial", "revision"])
def test_frozen_key_verify_and_atomic_legacy_export(tmp_path, kind):
    records = [model(0, "accepted"), model(1), model(2)]
    revisions = [model(9)]
    directory, _ = published(tmp_path, records, revisions=revisions)
    key = store.frozen_model_store_key(directory, kind=kind)
    assert key["kind"] == kind
    assert "model_store/manifest.json" in key["files"]
    assert store.verify_model_store_key(directory, key) is True
    output = tmp_path / f"{kind}.json"
    store.export_legacy_models(directory, output, kind=kind, frozen_key=key)
    expected = {"models": records, "accepted": records[:1], "partial": records[1:2], "revision": revisions}[kind]
    assert json.loads(output.read_text()) == expected


def test_accepted_and_revision_only_verify_their_own_dependencies(tmp_path, monkeypatch):
    directory, manifest = published(tmp_path, [model(0, "accepted"), model(1)], revisions=[model(2)])
    calls = []
    original = store._digest
    monkeypatch.setattr(store, "_digest", lambda p: (calls.append(Path(p).name), original(p))[1])
    key = store.frozen_model_store_key(directory, kind="accepted")
    assert calls == []  # Capture only reads small manifest/receipt snapshots.
    assert set(key["files"]) == {"model_store/manifest.json", "model_store/accepted.sqlite"}
    (directory / "model_store" / manifest["shards"][1]["path"]).write_bytes(b"broken rejected shard")
    assert list(store.iter_accepted_models(directory)) == [model(0, "accepted")]
    assert list(store.iter_revision_models(directory)) == [model(2)]
    assert "bodies.sqlite" not in calls
    with pytest.raises((ValueError, gzip.BadGzipFile)):
        list(store.iter_models(directory))


def test_partial_scope_excludes_unrelated_shards_and_preserves_explicit_order(tmp_path):
    records = [model(0), model(1), model(2)]
    directory, manifest = published(tmp_path, records, partial=[records[2], records[0], records[2]])
    key = store.frozen_model_store_key(directory, kind="partial")
    unused = manifest["shards"][1]["path"]
    assert f"model_store/{unused}" not in key["files"]
    (directory / "model_store" / unused).write_bytes(b"unused bad shard")
    assert list(store.iter_partial_models(directory, frozen_key=key)) == [records[2], records[0], records[2]]


def test_partial_sparse_repeated_reordered_records_keep_full_qc_and_independence(tmp_path):
    records = [model(i, "accepted_alternative_path" if i == 2 else "proposal") for i in range(13)]
    for i, record in enumerate(records):
        record["problems"] = ["missing_stop", "internal_stop"] if i % 3 == 0 else []
        record["unknown"]["lineage"] = {"source_rescue": [i, "植物🙂"]}
    partial = [records[i] for i in (12, 2, 2, 7, 0, 12)]
    directory = tmp_path / "worker"
    manifest = store.write_model_store(directory, records, partial_models=partial,
                                       shard_bytes=5500, codec="gzip")
    assert len(manifest["shards"]) > 1
    receipt(directory)
    restored = list(store.iter_partial_models(directory))
    assert json.dumps(restored, ensure_ascii=False) == json.dumps(partial, ensure_ascii=False)
    iterator = store.iter_partial_models(directory)
    first = next(iterator)
    first["unknown"]["lineage"]["source_rescue"].append("changed")
    first["raw_prediction"]["cds"][0][0] = -1
    second = next(iterator)
    second["problems"].append("changed")
    assert next(iterator) == records[2]
    assert list(iterator) == partial[3:]


def test_partial_lookup_sql_is_bounded_by_shards_instead_of_envelope_count(tmp_path, monkeypatch):
    records = [model(i) for i in range(120)]
    partial = [records[i] for i in (117, 2, 2, 79, 119)]
    directory = tmp_path / "worker"
    manifest = store.write_model_store(directory, records, partial_models=partial,
                                       shard_bytes=1024 * 1024, codec="gzip")
    assert len(manifest["shards"]) == 1
    receipt(directory)
    statements = []
    original = store._readonly
    def traced(path):
        db = original(path)
        if Path(path).name == "partial.sqlite":
            db.set_trace_callback(statements.append)
        return db
    monkeypatch.setattr(store, "_readonly", traced)
    assert list(store.iter_partial_models(directory)) == partial
    reference_reads = [sql for sql in statements if sql.lstrip().upper().startswith("SELECT")]
    assert len(reference_reads) <= 2 + len(manifest["shards"])


@pytest.mark.parametrize("line", [-1, 9999, None, 0.5])
def test_partial_invalid_line_reference_rejected_after_authenticating_receipt(tmp_path, line):
    records = [model(i) for i in range(3)]
    directory, _ = published(tmp_path, records, partial=[records[1]])
    with sqlite3.connect(directory / "model_store/partial.sqlite") as db:
        db.execute("UPDATE refs SET line=?", (line,))
    path = directory / "model_store/manifest.json"
    manifest = json.loads(path.read_text())
    manifest["files"]["partial.sqlite"] = hashlib.sha256(
        (directory / "model_store/partial.sqlite").read_bytes()).hexdigest()
    path.write_text(json.dumps(manifest))
    receipt(directory)
    with pytest.raises(ValueError, match="partial.*(line|reference)"):
        list(store.iter_partial_models(directory))


def test_partial_consumes_trailing_unreferenced_bytes_and_checks_hash(tmp_path):
    records = [model(i) for i in range(3)]
    directory = tmp_path / "worker"
    manifest = store.write_model_store(directory, records, partial_models=[records[0]],
                                       shard_bytes=1024 * 1024, codec="gzip")
    receipt(directory)
    key = store.frozen_model_store_key(directory, kind="partial")
    path = directory / "model_store" / manifest["shards"][0]["path"]
    with path.open("ab") as handle:
        handle.write(gzip.compress(b""))
    iterator = store.iter_partial_models(directory, frozen_key=key)
    assert next(iterator) == records[0]
    with pytest.raises(ValueError, match="shard checksum"):
        list(iterator)


def test_partial_authentic_shard_wrong_count_rejected_at_end(tmp_path):
    records = [model(i) for i in range(3)]
    directory = tmp_path / "worker"
    store.write_model_store(directory, records, partial_models=[records[0]],
                            shard_bytes=1024 * 1024, codec="gzip")
    path = directory / "model_store/manifest.json"
    manifest = json.loads(path.read_text())
    manifest["shards"][0]["count"] += 1
    manifest["counts"]["models"] += 1
    path.write_text(json.dumps(manifest))
    receipt(directory)
    iterator = store.iter_partial_models(directory)
    assert next(iterator) == records[0]
    with pytest.raises(ValueError, match="shard count"):
        list(iterator)


@pytest.mark.parametrize("member", ["accepted.sqlite", "bodies.sqlite", "models-000000.jsonl.gz"])
def test_mutated_or_missing_consumed_member_rejected(tmp_path, member):
    directory, _ = published(tmp_path, [model(0, "accepted")])
    kind = "accepted" if member == "accepted.sqlite" else "models"
    key = store.frozen_model_store_key(directory, kind=kind)
    path = directory / "model_store" / member
    with path.open("ab") as handle:
        handle.write(b"mutation")
    with pytest.raises(ValueError, match="checksum"):
        store.verify_model_store_key(directory, key)
    path.unlink()
    with pytest.raises(ValueError, match="Missing"):
        store.verify_model_store_key(directory, key)


def test_generation_and_end_fences(tmp_path):
    directory, _ = published(tmp_path, [model(0), model(1)])
    key = store.frozen_model_store_key(directory)
    iterator = store.iter_models(directory, frozen_key=key)
    assert next(iterator) == model(0)
    receipt_path = directory / "receipt.json"
    receipt_path.write_bytes(receipt_path.read_bytes() + b" ")
    with pytest.raises(ValueError, match="generation"):
        list(iterator)
    with pytest.raises(ValueError, match="generation"):
        store.verify_model_store_key(directory, key)


def test_member_mutation_after_initial_verification_rejected_at_end(tmp_path):
    directory, _ = published(tmp_path, [model(0), model(1)])
    iterator = store.iter_models(directory)
    next(iterator)
    path = directory / "model_store/bodies.sqlite"
    with path.open("ab") as handle:
        handle.write(b"changed while reading")
    with pytest.raises(ValueError, match="during consumption"):
        list(iterator)


@pytest.mark.parametrize("bad", ["../outside", "/absolute", "a/../b", "a\\b", "a//b", "./a", "bad\nname"])
def test_unsafe_receipt_members_rejected(tmp_path, bad):
    directory, _ = published(tmp_path, [])
    value = json.loads((directory / "receipt.json").read_text())
    value["files"][bad] = "0" * 64
    (directory / "receipt.json").write_text(json.dumps(value))
    with pytest.raises(ValueError, match="Unsafe"):
        store.frozen_model_store_key(directory)


def test_symlink_members_and_unreceipted_manifest_rejected(tmp_path):
    directory, _ = published(tmp_path, [model(0)])
    member = directory / "model_store/bodies.sqlite"
    saved = tmp_path / "saved.sqlite"
    member.rename(saved)
    member.symlink_to(saved)
    key = store.frozen_model_store_key(directory)
    with pytest.raises(ValueError, match="Symlink"):
        store.verify_model_store_key(directory, key)
    value = json.loads((directory / "receipt.json").read_text())
    del value["files"]["model_store/manifest.json"]
    (directory / "receipt.json").write_text(json.dumps(value))
    with pytest.raises(ValueError, match="not bound"):
        store.frozen_model_store_key(directory)


@pytest.mark.parametrize("kind", ["models", "accepted", "partial", "revision"])
def test_legacy_reader_export_and_missing_optional_files(tmp_path, kind):
    records = [model(0, "accepted"), model(1, "accepted_alternative_path")]
    (tmp_path / "models.json").write_text(json.dumps(records))
    receipt(tmp_path)
    expected = {"models": records, "accepted": records[:1], "partial": [], "revision": []}[kind]
    key = store.frozen_model_store_key(tmp_path, kind=kind)
    assert store.verify_model_store_key(tmp_path, key)
    output = tmp_path.parent / f"export-{tmp_path.name}-{kind}.json"
    store.export_legacy_models(tmp_path, output, kind=kind, frozen_key=key)
    assert json.loads(output.read_text()) == expected
    if kind in {"partial", "revision"}:
        optional = "partial_models.json" if kind == "partial" else "revision_candidates.json"
        (tmp_path / optional).write_text("[]")
        with pytest.raises(ValueError):
            store.verify_model_store_key(tmp_path, key)


def test_legacy_malformed_or_mutated_array_rejected(tmp_path):
    path = tmp_path / "models.json"
    path.write_text('[{"query":"x"},]')
    receipt(tmp_path)
    with pytest.raises(ValueError):
        list(store.iter_models(tmp_path))
    path.write_text("[]")
    with pytest.raises(ValueError, match="checksum"):
        list(store.iter_models(tmp_path))


def test_failed_write_cleans_and_does_not_replace_existing_publication(tmp_path):
    def bad():
        yield model(0)
        raise RuntimeError("producer failed")
    with pytest.raises(RuntimeError, match="producer failed"):
        store.write_model_store(tmp_path, bad(), codec="gzip")
    assert list(tmp_path.iterdir()) == []
    store.write_model_store(tmp_path, [model(0)], codec="gzip")
    before = (tmp_path / "model_store/manifest.json").read_bytes()
    with pytest.raises(ValueError, match="already exists"):
        store.write_model_store(tmp_path, [], codec="gzip")
    assert (tmp_path / "model_store/manifest.json").read_bytes() == before


def test_partial_nonmember_and_nonfinite_record_rejected(tmp_path):
    with pytest.raises(ValueError, match="exact"):
        store.write_model_store(tmp_path, [model(0)], partial_models=[model(1)], codec="gzip")
    assert list(tmp_path.iterdir()) == []
    with pytest.raises(ValueError):
        store.write_model_store(tmp_path, [{**model(0), "identity": float("nan")}], codec="gzip")
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("records", [[], [model(0)]])
def test_producer_can_read_before_receipt_only_with_explicit_unverified_mode(tmp_path, records):
    store.write_model_store(tmp_path, records, codec="gzip")
    assert list(store.iter_models(tmp_path, verify=False)) == records
    with pytest.raises(ValueError, match="receipt"):
        list(store.iter_models(tmp_path))
    with pytest.raises(ValueError, match="immutable"):
        store.export_legacy_models(tmp_path, tmp_path / "model_store/legacy.json", verify=False)


def test_missing_null_empty_fields_and_unknown_raw_fields_preserved(tmp_path):
    records = [{"x": 1}, {"sequence": None, "raw_prediction": None},
               {"raw_prediction": {"custom": {}, "cds": [], "paf": "not a PAF"}, "cds": []}]
    directory, _ = published(tmp_path, records)
    assert list(store.iter_models(directory)) == records


def test_actual_body_dedup_and_no_large_partial_duplicate(tmp_path):
    records = [model(i) for i in range(50)]
    directory, manifest = published(tmp_path, records)
    with sqlite3.connect(directory / "model_store/bodies.sqlite") as db:
        assert db.execute("SELECT COUNT(*) FROM bodies").fetchone()[0] == 1
        assert {x[0] for x in db.execute("SELECT name FROM sqlite_master WHERE type='table'")} == {"bodies"}
    with sqlite3.connect(directory / "model_store/partial.sqlite") as db:
        assert db.execute("SELECT COUNT(*) FROM refs").fetchone()[0] == 25
        assert db.execute("SELECT DISTINCT typeof(shard) FROM refs").fetchall() == [("integer",)]
        assert {x[0] for x in db.execute("SELECT name FROM sqlite_master WHERE type='table'")} == {"refs"}
    assert manifest["counts"]["models"] == 50


@pytest.mark.parametrize("explicit", [False, True])
def test_all_model_lookup_index_is_ephemeral_or_never_created(tmp_path, explicit, monkeypatch):
    records = [model(i) for i in range(6)]
    original = store.sqlite3.connect
    lookup_paths = []
    def track(path, *args, **kwargs):
        if Path(path).name == ".lookup.sqlite":
            lookup_paths.append(Path(path))
        return original(path, *args, **kwargs)
    monkeypatch.setattr(store.sqlite3, "connect", track)
    partial = [records[5], records[1]] if explicit else None
    directory, manifest = published(tmp_path, records, partial=partial)
    assert bool(lookup_paths) == explicit
    assert not any(p.exists() for p in lookup_paths)
    assert all("lookup" not in member for member in manifest["files"])
    with original(directory / "model_store/bodies.sqlite") as db:
        assert db.execute("SELECT COUNT(*) FROM sqlite_master WHERE name='model_refs'").fetchone()[0] == 0
    assert list(store.iter_partial_models(directory)) == (partial if explicit else records[1::2])


@pytest.mark.parametrize("explicit", [False, True])
def test_model_order_partial_reader_does_not_materialize_context_spool(tmp_path, explicit, monkeypatch):
    records = [model(i) for i in range(6)]
    partial = [records[1], records[1], records[5]] if explicit else None
    directory, manifest = published(tmp_path, records, partial=partial)
    assert manifest["partial_ordered"] is True
    def forbid(*args, **kwargs):
        raise AssertionError("Monotonic partial references must stream without spool")
    monkeypatch.setattr(store.tempfile, "TemporaryDirectory", forbid)
    assert list(store.iter_partial_models(directory)) == (partial if explicit else records[1::2])


def test_atomic_directory_export_with_bound_receipt(tmp_path):
    records = [model(0, "accepted"), model(1)]
    revisions = [model(2)]
    source, _ = published(tmp_path, records, revisions=revisions)
    destination = tmp_path / "legacy-export"
    proof = store.export_legacy_models(source, destination, directory=True)
    assert json.loads((destination / "models.json").read_text()) == records
    assert json.loads((destination / "partial_models.json").read_text()) == records[1:]
    assert json.loads((destination / "revision_candidates.json").read_text()) == revisions
    assert proof == json.loads((destination / "export_receipt.json").read_text())
    assert proof["source_keys"]["models"] == store.frozen_model_store_key(source)
    with pytest.raises(ValueError, match="already exists"):
        store.export_legacy_models(source, destination, directory=True)


def test_directory_export_source_generation_change_discards_partial_output(tmp_path, monkeypatch):
    source, _ = published(tmp_path, [model(0)])
    original = store.export_legacy_models
    def changing(*args, **kwargs):
        result = original(*args, **kwargs)
        if not kwargs.get("directory"):
            (source / "receipt.json").write_text((source / "receipt.json").read_text() + " ")
        return result
    monkeypatch.setattr(store, "export_legacy_models", changing)
    destination = tmp_path / "must-not-publish"
    with pytest.raises(ValueError, match="generation"):
        original(source, destination, directory=True)
    assert not destination.exists()
    assert not list(tmp_path.glob(".model-export-*"))


@pytest.mark.parametrize("change", ["schema", "codec", "count", "scope", "unsafe_shard"])
def test_unknown_or_malformed_manifest_rejected_even_with_new_receipt(tmp_path, change):
    directory, manifest = published(tmp_path, [model(0)])
    if change == "schema":
        manifest["schema"] = True
    elif change == "codec":
        manifest["codec"] = manifest["body_codec"] = "unknown"
    elif change == "count":
        manifest["counts"]["models"] = 9
    elif change == "scope":
        manifest["scopes"]["accepted"] = ["bodies.sqlite"]
    else:
        manifest["files"]["../outside"] = "0" * 64
    (directory / "model_store/manifest.json").write_text(json.dumps(manifest))
    receipt(directory)
    with pytest.raises(ValueError):
        store.frozen_model_store_key(directory)


def test_bounded_body_and_envelope_fail_without_partial_publication(tmp_path, monkeypatch):
    monkeypatch.setattr(store, "MAX_RECORD_BYTES", 128)
    with pytest.raises(ValueError, match="bounded"):
        store.write_model_store(tmp_path, [model(0)], codec="gzip")
    assert not list(tmp_path.iterdir())


def test_legacy_export_cannot_overwrite_its_source(tmp_path):
    records = [model(0)]
    (tmp_path / "models.json").write_text(json.dumps(records))
    receipt(tmp_path)
    with pytest.raises(ValueError, match="consumed source"):
        store.export_legacy_models(tmp_path, tmp_path / "models.json")
    assert json.loads((tmp_path / "models.json").read_text()) == records

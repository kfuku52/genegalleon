"""Writer validation/publication boundaries and selective-copy independence."""
import copy
import hashlib

import pytest

from workflow.support import rescue_model_store as store
from workflow.tests.test_rescue_model_store import model, receipt

BAD_CONTEXT = [
    pytest.param(float("nan"), ValueError, "Out of range float values", id="nonfinite"),
    pytest.param({"unencodable"}, TypeError, "not JSON serializable", id="nonjson"),
]


def assert_unpublished(directory):
    assert not (directory / "model_store").exists()
    assert not list(directory.iterdir())


def seeded_collision(monkeypatch):
    """Exercise the actual collision check against unequal compressed body bytes."""
    original = store._insert_body

    def insert(db, body, level, codec, cache=None):
        digest = hashlib.sha256(body).hexdigest()
        db.execute("INSERT OR REPLACE INTO bodies VALUES (?,?)",
                   (digest, store._compress(b"{}", codec, level)))
        return original(db, body, level, codec, cache)

    monkeypatch.setattr(store, "_insert_body", insert)


@pytest.mark.parametrize("stage", ["models", "revision"])
@pytest.mark.parametrize("nested", [False, True], ids=["candidate", "raw_prediction"])
@pytest.mark.parametrize("invalid,error,message", BAD_CONTEXT)
def test_invalid_context_is_rejected_before_real_body_collision(
        tmp_path, monkeypatch, stage, nested, invalid, error, message):
    record = model(0, "accepted")
    target = record["raw_prediction"] if nested else record
    target["evidence"]["invalid"] = invalid
    seeded_collision(monkeypatch)
    directory = tmp_path / "worker"
    with pytest.raises(error, match=message):
        store.write_model_store(directory, [record] if stage == "models" else [],
                                revisions=[record] if stage == "revision" else [],
                                codec="gzip")
    assert_unpublished(directory)


def test_collision_fixture_rejects_valid_json_without_publication(tmp_path, monkeypatch):
    seeded_collision(monkeypatch)
    directory = tmp_path / "worker"
    with pytest.raises(ValueError, match="body hash collision"):
        store.write_model_store(directory, [model()], codec="gzip")
    assert_unpublished(directory)


@pytest.mark.parametrize("invalid,error,message", BAD_CONTEXT)
def test_default_pack_still_validates_candidate_context(invalid, error, message):
    record = model()
    record["evidence"]["invalid"] = invalid
    with pytest.raises(error, match=message):
        store._pack(record)


def test_body_encoding_failure_precedes_invalid_candidate_context(tmp_path):
    record = model()
    record["sequence"] = {"not_json"}
    record["evidence"]["invalid"] = float("nan")
    directory = tmp_path / "worker"
    with pytest.raises(TypeError, match="not JSON serializable"):
        store.write_model_store(directory, [record], codec="gzip")
    assert_unpublished(directory)


def test_default_overlay_inserts_body_before_context_encoding(tmp_path):
    body, envelope = store._pack(model())
    envelope["context"]["nonfinite"] = float("nan")
    db = store._overlay(tmp_path / "direct.sqlite")
    try:
        with pytest.raises(ValueError, match="Out of range float values"):
            store._add_overlay(db, 0, body, envelope, 1, "gzip")
        assert db.execute("SELECT data FROM bodies").fetchone() is not None
        assert db.execute("SELECT ordinal FROM entries").fetchone() is None
    finally:
        db.close()


def test_default_overlay_preserves_collision_before_invalid_context(tmp_path, monkeypatch):
    body, envelope = store._pack(model())
    envelope["context"]["nonfinite"] = float("nan")
    db = store._overlay(tmp_path / "direct.sqlite")
    seeded_collision(monkeypatch)
    try:
        with pytest.raises(ValueError, match="body hash collision"):
            store._add_overlay(db, 0, body, envelope, 1, "gzip")
    finally:
        db.close()


@pytest.mark.parametrize("codec", ["gzip", "auto"])
@pytest.mark.parametrize("status", ["proposal", "accepted"])
def test_record_bound_includes_jsonl_newline(tmp_path, monkeypatch, codec, status):
    record = model(0, status)
    body, envelope = store._pack(record)
    envelope["ordinal"] = 0
    final_bytes = store._json(envelope)
    assert len(body) < len(final_bytes)

    monkeypatch.setattr(store, "MAX_RECORD_BYTES", len(final_bytes))
    rejected = tmp_path / "without_newline_capacity"
    with pytest.raises(ValueError, match="record exceeds bounded buffer"):
        store.write_model_store(rejected, [record], codec=codec)
    assert_unpublished(rejected)

    monkeypatch.setattr(store, "MAX_RECORD_BYTES", len(final_bytes) + 1)
    accepted = tmp_path / "with_newline_capacity"
    store.write_model_store(accepted, [record], codec=codec)
    receipt(accepted)
    assert list(store.iter_models(accepted)) == [record]
    assert list(store.iter_accepted_models(accepted)) == ([record] if status == "accepted" else [])


@pytest.mark.parametrize("codec", ["gzip", "auto"])
def test_shared_bodies_have_independent_accepted_partial_and_revision_contexts(tmp_path, codec):
    records = [model(0, "accepted"), model(1), model(2, "accepted")]
    revision = copy.deepcopy(records[0])
    revision["revision_owner_ids"] = ["original_gene"]
    revision["evidence"]["revision_note"] = "retained"
    partials = [records[2], records[0], records[2]]
    directory = tmp_path / "worker"
    store.write_model_store(directory, records, partial_models=partials,
                            revisions=[revision], codec=codec, shard_bytes=900)
    receipt(directory)

    accepted = list(store.iter_accepted_models(directory))
    assert accepted == [records[0], records[2]]
    assert list(accepted[0]) == list(records[0])
    assert list(accepted[0]["raw_prediction"]) == list(records[0]["raw_prediction"])
    accepted[0]["cds"][0][0] = 999
    accepted[0]["raw_prediction"]["evidence"]["donor"] = "changed"
    accepted[0]["unknown"]["empty"].append("changed")
    assert accepted[1] == records[2]
    assert list(store.iter_models(directory)) == records
    assert list(store.iter_revision_models(directory)) == [revision]
    assert list(store.iter_accepted_models(directory)) == [records[0], records[2]]

    restored_partial = list(store.iter_partial_models(directory))
    assert restored_partial == partials
    restored_partial[0]["raw_prediction"]["cds"][0][0] = 888
    restored_partial[0]["support"].append({"donor": "changed"})
    assert restored_partial[2] == records[2]
    assert list(store.iter_partial_models(directory)) == partials

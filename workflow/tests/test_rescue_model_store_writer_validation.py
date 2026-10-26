"""Writer input validation, atomic publication and selective-copy independence."""
import copy

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


@pytest.mark.parametrize("stage", ["models", "revision"])
@pytest.mark.parametrize("nested", [False, True], ids=["candidate", "raw_prediction"])
@pytest.mark.parametrize("invalid,error,message", BAD_CONTEXT)
def test_writer_rejects_invalid_context_without_publication(
        tmp_path, stage, nested, invalid, error, message):
    record = model(0, "accepted")
    target = record["raw_prediction"] if nested else record
    target["evidence"]["invalid"] = invalid
    directory = tmp_path / "worker"
    with pytest.raises(error, match=message):
        store.write_model_store(directory, [record] if stage == "models" else [],
                                revisions=[record] if stage == "revision" else [],
                                codec="gzip")
    assert_unpublished(directory)


def test_writer_rejects_non_json_body_without_publication(tmp_path):
    record = model()
    record["sequence"] = {"not_json"}
    directory = tmp_path / "worker"
    with pytest.raises(TypeError, match="not JSON serializable"):
        store.write_model_store(directory, [record], codec="gzip")
    assert_unpublished(directory)


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

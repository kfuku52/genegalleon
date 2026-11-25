import base64
import hashlib
import json
import zipfile

import pytest

from workflow.support import archive_rescue_model_store as migration
from workflow.support import rescue_model_store as store
from workflow.tests.test_rescue_model_store import model, published


def test_archives_and_resumes_after_only_the_manifest_was_published(tmp_path, monkeypatch):
    models = [model(0, "accepted"), model(1)]
    worker, initial = published(tmp_path, models)
    directory = worker / "model_store"
    receipt = worker / "receipt.json"
    original_receipt = receipt.read_bytes()
    original_manifest = (directory / "manifest.json").read_bytes()
    original_shards = {shard["path"]: (directory / shard["path"]).read_bytes() for shard in initial["shards"]}
    original_replace = migration._replace

    def stop_before_receipt(path, raw):
        if path == receipt:
            raise OSError("simulated interrupt")
        return original_replace(path, raw)

    with monkeypatch.context() as patch:
        patch.setattr(migration, "_replace", stop_before_receipt)
        with pytest.raises(OSError, match="simulated interrupt"):
            migration.archive_closed_store(worker)
    assert receipt.read_bytes() == original_receipt
    assert all((directory / name).read_bytes() == data for name, data in original_shards.items())

    converted = migration.archive_closed_store(worker)
    assert converted["removed"] == len(initial["shards"])
    assert converted["original_receipt_sha256"] == hashlib.sha256(original_receipt).hexdigest()
    assert list(store.iter_models(worker)) == models
    assert list(store.iter_accepted_models(worker)) == models[:1]
    assert list(store.iter_partial_models(worker)) == models[1:]
    assert migration.archive_closed_store(worker)["removed"] == 0
    saved = json.loads((directory / migration.RECOVERY).read_bytes())
    assert base64.b64decode(saved["receipt"]) == original_receipt
    assert base64.b64decode(saved["manifest"]) == original_manifest
    with zipfile.ZipFile(directory / store.SHARD_ARCHIVE) as archive:
        assert {name: archive.read(name) for name in archive.namelist()} == original_shards


def test_corrupt_source_never_loses_loose_shards(tmp_path):
    worker, initial = published(tmp_path, [model(0), model(1)])
    directory = worker / "model_store"
    shard = initial["shards"][0]["path"]
    (directory / shard).write_bytes(b"bad source")
    with pytest.raises(ValueError, match="source shard changed"):
        migration.archive_closed_store(worker)
    assert (directory / shard).read_bytes() == b"bad source"
    assert not (directory / store.SHARD_ARCHIVE).exists()
    assert not (directory / migration.RECOVERY).exists()

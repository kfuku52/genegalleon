"""Differential publication checks and real content-mutation rejection."""
import hashlib
import sys

import pytest

from workflow.support import rescue_gene_models as rescue

state = sys.modules[rescue.digest.__module__]
hash_outputs = rescue.hash_outputs


@pytest.mark.parametrize("workers", [1, 4])
def test_manifest_matches_regular_outputs_and_retains_legacy_exclusions(tmp_path, workers):
    expected = {}
    for i in range(12):
        directory = tmp_path / "intervals" / str(i) / "logs"
        directory.mkdir(parents=True)
        for name, payload in (("command.json", b"[]\n"), ("models.gff", bytes([i]) * 2000), ("empty", b"")):
            path = directory / name
            path.write_bytes(payload)
            expected[str(path.relative_to(tmp_path))] = hashlib.sha256(payload).hexdigest()
        for excluded in ("genome.fa", "genome.mpi"):
            (directory / excluded).write_text("execution scratch")
    (tmp_path / "file_link").symlink_to(tmp_path / "intervals/0/logs/models.gff")
    (tmp_path / "directory_link").symlink_to(tmp_path / "intervals/0", target_is_directory=True)
    assert hash_outputs(tmp_path, workers=workers) == expected


@pytest.mark.parametrize("workers", [1, 4])
def test_receipt_verification_rejects_changed_and_missing_bytes(tmp_path, workers):
    key = {"plan": "unchanged"}
    rescue.stage(tmp_path, "published", key, lambda p: (p / "model.json").write_text("original"), hash_workers=workers)
    directory = tmp_path / "published"
    assert rescue.verified(directory, key, hash_workers=workers)
    (directory / "model.json").write_text("modified")
    assert not rescue.verified(directory, key, hash_workers=workers)
    (directory / "model.json").unlink()
    assert not rescue.verified(directory, key, hash_workers=workers)


def test_parallel_hash_detects_mutation_and_preserves_previous_publication(tmp_path, monkeypatch):
    key = {"plan": "old"}
    rescue.stage(tmp_path, "published", key, lambda p: (p / "stable").write_text("previous"))
    original_count = state.count
    changed = False
    target = None

    def mutate_during_read(name, value=1):
        nonlocal changed
        if name == "sha256_bytes" and not changed:
            changed = True
            with target.open("ab") as handle:
                handle.write(b"changed while hashing")
        original_count(name, value)

    def build(directory):
        nonlocal target
        target = directory / "large"
        target.write_bytes(b"A" * (4 * 1024 * 1024))

    with monkeypatch.context() as patch:
        patch.setattr(state, "count", mutate_during_read)
        with pytest.raises(OSError, match="File changed while hashing"):
            rescue.stage(tmp_path, "published", {"plan": "new"}, build, hash_workers=4)
    assert changed
    assert rescue.verified(tmp_path / "published", key, hash_workers=4)
    assert (tmp_path / "published.failed/large").is_file()

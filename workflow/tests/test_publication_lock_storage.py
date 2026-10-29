"""Offline publication metadata retirement preserves real data and live owners."""
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
import publication_lock_storage as storage


@pytest.fixture
def fixture(tmp_path):
    root = tmp_path / "orthogroup"
    directory = root / "mafft"
    directory.mkdir(parents=True)
    output = directory / "OG0000001_cds.aln.fa.gz"
    output.write_bytes(b"scientific output remains exact")
    marker = directory / "..OG0000001_cds.aln.fa.gz.gg-bundle.lock.guard"
    marker.touch(mode=0o600)
    namespace = Path(str(marker) + ".namespace-v1")
    namespace.mkdir(mode=0o700)
    (namespace / "readers").mkdir(mode=0o700)
    journal = tmp_path / "journal"
    journal.mkdir(mode=0o700)
    stripes = directory / ".gg_publish_lock_stripes-v1"
    stripes.mkdir()
    (stripes / "00.lock").write_bytes(b"current stripe unchanged")
    return root, output, marker, namespace, journal, stripes


def test_exact_retirement_preserves_outputs_current_stripes_and_resume(fixture):
    root, output, marker, namespace, journal, stripes = fixture
    before = storage._identity(output)
    expected = storage.snapshot(root)
    assert expected["groups"] == 1 and expected["entries"] == 3
    result = storage.retire(root, journal, expected, writers_stopped=True)
    assert result["verified"] and result["retired_entries"] == 3
    assert not marker.exists() and not namespace.exists()
    assert output.read_bytes() == b"scientific output remains exact"
    assert storage._identity(output) == before
    assert (stripes / "00.lock").read_bytes() == b"current stripe unchanged"
    assert storage.retire(root, journal, expected, writers_stopped=True) == result


def test_quiescence_is_explicit_and_empty_guard_does_not_clear_primary_lease(fixture):
    root, _, marker, _, journal, _ = fixture
    expected = storage.snapshot(root)
    with pytest.raises(storage.PublicationLockError, match="must be fenced"):
        storage.retire(root, journal, expected)
    lease = marker.with_name(marker.name[1:-len(".guard")])
    lease.write_text("live primary owner")
    with pytest.raises(storage.PublicationLockError, match="Primary publication lease"):
        storage.retire(root, journal, expected, writers_stopped=True)
    assert marker.exists() and lease.read_text() == "live primary owner"
    assert not list(journal.iterdir())


@pytest.mark.parametrize("member", ["gate", "readers/live-reader", "unexpected"])
def test_busy_namespace_fails_closed_without_removing_anything(fixture, member):
    root, _, marker, namespace, journal, _ = fixture
    expected = storage.snapshot(root)
    target = namespace / member
    target.mkdir() if member == "gate" else target.write_text("owned")
    with pytest.raises(storage.PublicationLockError):
        storage.retire(root, journal, expected, writers_stopped=True)
    assert marker.exists() and target.exists() and not list(journal.iterdir())


def test_symlink_and_nonempty_marker_are_rejected(fixture, tmp_path):
    root, _, marker, _, _, _ = fixture
    marker.write_text("not an empty marker")
    with pytest.raises(storage.PublicationLockError, match="Nonempty"):
        storage.snapshot(root)
    marker.unlink()
    target = tmp_path / "foreign"
    target.write_text("retain")
    marker.symlink_to(target)
    with pytest.raises(storage.PublicationLockError, match="Unsafe"):
        storage.snapshot(root)
    assert target.read_text() == "retain"


def test_changed_inventory_and_added_guard_rejected_before_deletion(fixture):
    root, _, marker, namespace, journal, _ = fixture
    expected = storage.snapshot(root)
    before = marker.stat()
    os.utime(marker, ns=(before.st_atime_ns, before.st_mtime_ns + 1_000_000_000))
    with pytest.raises(storage.PublicationLockError, match="snapshot changed"):
        storage.retire(root, journal, expected, writers_stopped=True)
    assert marker.exists() and namespace.exists()
    assert not list(journal.iterdir())


def test_partial_same_operation_can_resume_but_changed_identity_cannot(fixture, monkeypatch):
    root, output, marker, namespace, journal, _ = fixture
    expected = storage.snapshot(root)
    original = Path.rmdir
    def fail_after_reader(path):
        if path == namespace:
            raise OSError("interrupted after reader directory retirement")
        return original(path)
    with monkeypatch.context() as patch:
        patch.setattr(Path, "rmdir", fail_after_reader)
        with pytest.raises(OSError, match="interrupted"):
            storage.retire(root, journal, expected, writers_stopped=True)
    assert not (namespace / "readers").exists() and namespace.exists() and marker.exists()
    assert storage.retire(root, journal, expected, writers_stopped=True)["verified"]
    assert output.exists()


def test_original_journal_corruption_rejected_before_deletion(fixture, monkeypatch):
    root, _, marker, namespace, journal, _ = fixture
    expected = storage.snapshot(root)
    with monkeypatch.context() as patch:
        patch.setattr(storage, "_remaining", lambda *args: (_ for _ in ()).throw(OSError("stop")))
        with pytest.raises(OSError):
            storage.retire(root, journal, expected, writers_stopped=True)
    (journal / "publication-guards-original.jsonl.gz").write_bytes(b"corrupt")
    with pytest.raises(OSError):
        storage.retire(root, journal, expected, writers_stopped=True)
    assert marker.exists() and namespace.exists()


def test_new_guard_after_journal_does_not_authorize_new_deletions(fixture, monkeypatch):
    root, _, marker, namespace, journal, _ = fixture
    expected = storage.snapshot(root)
    with monkeypatch.context() as patch:
        patch.setattr(storage, "_remaining", lambda *args: (_ for _ in ()).throw(OSError("stop")))
        with pytest.raises(OSError):
            storage.retire(root, journal, expected, writers_stopped=True)
    new = marker.with_name("..AAA.gg-bundle.lock.guard")
    new.touch(mode=0o600)
    ns = Path(str(new) + ".namespace-v1")
    ns.mkdir(mode=0o700)
    (ns / "readers").mkdir(mode=0o700)
    with pytest.raises(storage.PublicationLockError, match="Unrecorded"):
        storage.retire(root, journal, expected, writers_stopped=True)
    assert marker.exists() and namespace.exists() and new.exists()


def test_incomplete_unjournaled_guard_cannot_be_retired(fixture):
    root, _, marker, namespace, journal, _ = fixture
    expected = storage.snapshot(root)
    (namespace / "readers").rmdir()
    with pytest.raises(storage.PublicationLockError, match="Incomplete"):
        storage.retire(root, journal, expected, writers_stopped=True)
    assert marker.exists() and namespace.exists() and not list(journal.iterdir())


def test_new_inode_after_original_journal_is_rejected(fixture, monkeypatch):
    root, _, marker, namespace, journal, _ = fixture
    expected = storage.snapshot(root)
    with monkeypatch.context() as patch:
        patch.setattr(storage, "_remaining", lambda *args: (_ for _ in ()).throw(OSError("stop")))
        with pytest.raises(OSError):
            storage.retire(root, journal, expected, writers_stopped=True)
    original = marker.with_name("retained-original-marker")
    marker.rename(original)
    marker.touch(mode=0o600)
    with pytest.raises(storage.PublicationLockError, match="identity changed"):
        storage.retire(root, journal, expected, writers_stopped=True)
    assert marker.exists() and namespace.exists() and original.exists()


def test_native_zip_materialization_preserves_complete_and_failed_states(fixture, tmp_path):
    import gene_family_output_store as store_module
    root, output, _, _, journal, _ = fixture
    other = output.with_name("OG0000002_cds.aln.fa.gz")
    other.write_bytes(b"failed family partial output")
    store = store_module.GeneFamilyOutputStore(root)
    store.mark_family_state("OG0000001", "complete", "original-complete")
    store.mark_family_state("OG0000002", "failed", "original-failed")
    families = ["OG0000001", "OG0000002"]
    before = {family: store.family_observation(family) for family in families}
    contents = {path.name: path.read_bytes() for path in [output, other]}
    expected = storage.snapshot(root)
    storage.retire(root, journal, expected, writers_stopped=True)
    store_module.convert_storage_to_zip(root, "orthogroup", families,
        store_module.orthogroup_id_from_name, available_bytes=1024**3)
    assert store.verify(deep=True)
    for family in families:
        assert store.family_observation(family) == before[family]
        restored = store.materialize_family(family, store_module.orthogroup_id_from_name,
            destination_root=tmp_path / ("resume-" + family))
        assert restored and all(path.read_bytes() == contents[path.name] for path in restored)
    assert storage.verify_retired(root, journal, expected)["retired_entries"] == 3


def test_unguarded_primary_lease_is_also_rejected(fixture):
    root, _, marker, _, journal, _ = fixture
    expected = storage.snapshot(root)
    lease = marker.parent / ".another-output.gg-bundle.lock"
    lease.write_text("unreconciled primary owner")
    with pytest.raises(storage.PublicationLockError, match="Primary publication lease"):
        storage.retire(root, journal, expected, writers_stopped=True)
    assert marker.exists() and lease.exists() and not list(journal.iterdir())

"""Queue/collector invariants, including failures at publication boundaries."""
import concurrent.futures
import os
from pathlib import Path

import pytest

from workflow.support import gene_family_output_store as store


def read_bytes(logical, subdir, name):
    with logical.open_binary(subdir, name) as handle:
        return handle.read()


def queued(root, family="OG0000001", token="run-1", content=b"alignment\n"):
    path = root / "mafft" / f"{family}_cds.aln.fa.gz"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    inventory = store.family_inventory_path(root, family)
    inventory.mkdir(parents=True, exist_ok=True)
    (inventory / "worker.paths").write_bytes(os.fsencode(path) + b"\0")
    logical = store.GeneFamilyOutputStore(root)
    logical.mark_family_state(family, "running", token)
    logical.mark_family_state(family, "complete", token)
    request = store.enqueue_family_archive(root, "orthogroup", family, token)
    return path, request


def drain(root, **kwargs):
    return store.drain_archive_queue(root, "orthogroup", store.orthogroup_id_from_name, **kwargs)


def test_collector_batches_without_scanning_outputs_or_compacting(tmp_path, monkeypatch):
    paths = [queued(tmp_path, f"OG{i:07d}")[0] for i in range(7)]
    original_iterdir = Path.iterdir

    def no_output_scan(path):
        if path == tmp_path / "mafft":
            raise AssertionError("Collector scanned the output directory")
        return original_iterdir(path)

    monkeypatch.setattr(Path, "iterdir", no_output_scan)
    monkeypatch.setattr(store, "_compact_artifact_chunk", lambda *a, **k: pytest.fail("Inline compaction"))
    # Source cleanup may check whether a directory is empty; do not enumerate
    # it on the queue hot path, which deliberately retains live directories.
    for expected in (3, 3, 1):
        result = drain(tmp_path, batch_families=3)
        assert result["selected_families"] == expected
        assert result["archived_files"] == expected
    assert store.archive_queue_status(tmp_path)["pending_families"] == 0
    assert all(not path.exists() for path in paths)
    logical = store.GeneFamilyOutputStore(tmp_path)
    for path in paths:
        assert read_bytes(logical, "mafft", path.name) == b"alignment\n"
    assert len(list((tmp_path / "archives" / "mafft").glob("*.zip"))) == 3


def test_compression_does_not_hold_global_locks(tmp_path, monkeypatch):
    queued(tmp_path)
    original = store._archive_chunk

    def compress(*args, **kwargs):
        with store.producer_quiescence_lock(tmp_path / ".gg_store", nonblocking=True) as acquired:
            assert acquired
            with store.archive_lock(tmp_path / ".gg_store", nonblocking=True) as acquired:
                assert acquired
        with store.family_bucket_lock(tmp_path / ".gg_store", "OG0000002", exclusive=True, nonblocking=True) as acquired:
            assert acquired
        return original(*args, **kwargs)

    monkeypatch.setattr(store, "_archive_chunk", compress)
    assert drain(tmp_path)["status"] == "committed"


def test_busy_family_does_not_block_unrelated_family(tmp_path):
    first, _ = queued(tmp_path, "OG0000001")
    second, _ = queued(tmp_path, "OG0000002")
    with store.family_bucket_lock(tmp_path / ".gg_store", "OG0000001", exclusive=False):
        result = drain(tmp_path)
    assert result["archived_files"] == 1
    assert result["deferred_families"] == 1
    assert first.exists() and not second.exists()
    assert drain(tmp_path)["archived_files"] == 1


@pytest.mark.parametrize("failure", ["compress", "index", "remove", "ack"])
def test_failure_preserves_results_and_request_for_retry(tmp_path, monkeypatch, failure):
    path, request = queued(tmp_path)
    original_unlink = Path.unlink
    with monkeypatch.context() as patch:
        def fail(*args, **kwargs):
            raise OSError("injected failure")
        if failure == "compress":
            patch.setattr(store, "_archive_chunk", fail)
        elif failure == "index":
            # Crash after the durable pending marker but before index commit.
            patch.setattr(store.GeneFamilyOutputStore, "_merge_index_subdirs_uncommitted", fail)
        elif failure == "remove":
            patch.setattr(store, "_remove_archived_sources", fail)
        else:
            def unlink(target, *args, **kwargs):
                if target == request:
                    fail()
                return original_unlink(target, *args, **kwargs)
            patch.setattr(Path, "unlink", unlink)
        with pytest.raises(OSError, match="injected"):
            drain(tmp_path)
    assert request.exists()
    if failure == "index":
        assert path.read_bytes() == b"alignment\n"
        store.repair_archive_index(tmp_path)
    assert drain(tmp_path)["status"] == "committed"
    assert read_bytes(store.GeneFamilyOutputStore(tmp_path), "mafft", path.name) == b"alignment\n"
    assert not path.exists() and not request.exists()


def test_changed_source_is_not_committed_or_deleted(tmp_path, monkeypatch):
    path, request = queued(tmp_path)
    original = store._archive_chunk

    def compress(*args, **kwargs):
        result = original(*args, **kwargs)
        path.write_bytes(b"new result\n")
        return result

    monkeypatch.setattr(store, "_archive_chunk", compress)
    with pytest.raises(store.ArchiveStoreError, match="changed"):
        drain(tmp_path)
    assert path.read_bytes() == b"new result\n" and request.exists()
    assert not list((tmp_path / "archives").glob("*/*.zip"))


def test_stale_request_cannot_archive_new_run(tmp_path):
    path, request = queued(tmp_path)
    logical = store.GeneFamilyOutputStore(tmp_path)
    logical.mark_family_state("OG0000001", "running", "run-2")
    assert drain(tmp_path)["deferred_families"] == 1
    logical.mark_family_state("OG0000001", "complete", "run-2")
    assert drain(tmp_path)["stale_requests"] == 1
    assert path.exists() and request.exists()
    with pytest.raises(store.ArchiveStoreError, match="Stale"):
        store.enqueue_family_archive(tmp_path, "orthogroup", "OG0000001", "run-1")
    store.enqueue_family_archive(tmp_path, "orthogroup", "OG0000001", "run-2")
    assert drain(tmp_path)["archived_files"] == 1


def test_byte_limit_and_oversize_single_family(tmp_path):
    for i in range(3):
        queued(tmp_path, f"OG{i:07d}", content=b"x" * 100)
    assert drain(tmp_path, batch_bytes=150)["selected_families"] == 1
    assert drain(tmp_path, batch_bytes=1)["selected_families"] == 1


def test_maintenance_gate_excludes_new_families(tmp_path):
    with store.all_family_bucket_locks(tmp_path / ".gg_store"):
        with store.family_bucket_lock(tmp_path / ".gg_store", "OG-new", exclusive=False, nonblocking=True) as acquired:
            assert not acquired
    with store.family_bucket_lock(tmp_path / ".gg_store", "OG-new", exclusive=False):
        with store._bucket_lock(store.family_gate_path(tmp_path / ".gg_store"), exclusive=True, nonblocking=True) as acquired:
            assert not acquired


def _publish_from_process(args):
    root, index = args
    return queued(root, f"OG{index:07d}")


def test_concurrent_queue_publication(tmp_path):
    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
        paths = list(executor.map(_publish_from_process, [(tmp_path, i) for i in range(128)]))
    assert store.archive_queue_status(tmp_path)["pending_families"] == 128
    while drain(tmp_path)["status"] == "committed":
        pass
    assert store.archive_queue_status(tmp_path)["pending_families"] == 0
    logical = store.GeneFamilyOutputStore(tmp_path)
    assert all(read_bytes(logical, "mafft", path.name) == b"alignment\n" for path, _ in paths)


def test_incomplete_inventory_is_not_acknowledged(tmp_path):
    path, request = queued(tmp_path)
    (store.family_inventory_path(tmp_path, "OG0000001") / "worker.paths").write_bytes(os.fsencode(path))
    with pytest.raises(store.ArchiveStoreError, match="Incomplete"):
        drain(tmp_path)
    assert path.exists() and request.exists()


@pytest.mark.parametrize("boundary", ["compressed", "published", "indexed", "removed"])
def test_sigkill_requires_explicit_owner_reconciliation_then_retries(tmp_path, boundary):
    import signal
    import subprocess
    import sys

    from workflow.support import shared_namespace_lock as locks

    path, request = queued(tmp_path)
    script = '''
import os, signal, sys
from pathlib import Path
from workflow.support import gene_family_output_store as m
root, boundary = Path(sys.argv[1]), sys.argv[2]
def killed():
    os.kill(os.getpid(), signal.SIGKILL)
if boundary == "compressed":
    original = m._archive_chunk
    def compress(*a, **k):
        result = original(*a, **k)
        killed()
    m._archive_chunk = compress
elif boundary == "published":
    original = os.replace
    def publish(source, destination):
        result = original(source, destination)
        if Path(destination).parent.parent == root / "archives":
            killed()
        return result
    os.replace = publish
elif boundary == "indexed":
    original = m.GeneFamilyOutputStore._merge_index_subdirs_uncommitted
    def index(*a, **k):
        result = original(*a, **k)
        killed()
    m.GeneFamilyOutputStore._merge_index_subdirs_uncommitted = index
else:
    original = m._remove_archived_sources
    def remove(*a, **k):
        result = original(*a, **k)
        killed()
    m._remove_archived_sources = remove
m.drain_archive_queue(root, "orthogroup", m.orthogroup_id_from_name)
'''
    process = subprocess.Popen([sys.executable, "-c", script, str(tmp_path), boundary],
                               env={**os.environ, "GG_JOB_ID": "test-killed-collector"})
    assert process.wait(timeout=30) == -signal.SIGKILL
    assert request.exists()
    assert path.exists() == (boundary != "removed")
    assert drain(tmp_path)["status"] == "collector-busy"
    released = 0
    # Simulate an operator only after wait() confirms this test process died.
    for namespace in (tmp_path / ".gg_store_locks").rglob("*.namespace-v1"):
        lock_path = Path(str(namespace).removesuffix(".namespace-v1"))
        evidence = locks.inspect_lock(lock_path)
        for owner in evidence["shared"]:
            assert owner["pid"] == process.pid
            assert owner["job_id"] == "test-killed-collector"
            locks.release(lock_path, owner["token"], exclusive=False)
            released += 1
        if evidence["exclusive"]:
            owner = evidence["exclusive"]
            assert owner["pid"] == process.pid
            locks.release(lock_path, owner["token"], exclusive=True)
            released += 1
    assert released >= 3
    if boundary in {"published", "indexed"}:
        assert (tmp_path / ".gg_store" / store.INDEX_UPDATE_FILE).exists()
        store.repair_archive_index(tmp_path)
    assert drain(tmp_path)["archived_files"] == (0 if boundary == "removed" else 1)
    assert read_bytes(store.GeneFamilyOutputStore(tmp_path), "mafft", path.name) == b"alignment\n"
    assert not list((tmp_path / ".gg_store" / "queue-staging").iterdir())


def test_commit_contention_keeps_request_and_live_files(tmp_path, monkeypatch):
    from workflow.support.shared_namespace_lock import acquire, release

    path, request = queued(tmp_path)
    original = store._archive_chunk
    producer_lock = store._store_lock_path(tmp_path / ".gg_store" / store.PRODUCER_LOCK_FILE)
    tokens = []

    def compress(*args, **kwargs):
        result = original(*args, **kwargs)
        tokens.append(acquire(producer_lock, exclusive=False))
        return result

    with monkeypatch.context() as patch:
        patch.setattr(store, "_archive_chunk", compress)
        result = drain(tmp_path, nonblocking=True)
    assert result["status"] == "readers-busy"
    assert request.exists() and path.exists()
    assert not list((tmp_path / "archives").glob("*/*.zip"))
    for token in tokens:
        release(producer_lock, token, exclusive=False)
    assert drain(tmp_path)["archived_files"] == 1


@pytest.mark.parametrize("symlink_root", [False, True])
def test_publisher_journals_are_compacted_and_cli_settings_are_preserved(tmp_path, symlink_root):
    import subprocess
    import zipfile

    repo = Path(__file__).resolve().parents[2]
    family = "OG0000001"
    inventory = store.family_inventory_path(tmp_path, family)
    inventory.mkdir(parents=True)
    script = r'''
set -euo pipefail
source "$1/workflow/support/gg_util.sh"
export GG_FAMILY_OUTPUT_ROOT="$2" GG_FAMILY_OUTPUT_INVENTORY="$3"
GG_FAMILY_OUTPUT_CANONICAL_ROOT=$(cd "$2" && pwd -P)
export GG_FAMILY_OUTPUT_CANONICAL_ROOT
mkdir -p "$2/mafft" "$2/scratch"
printf 'first\n' > "$2/scratch/a"
printf 'second\n' > "$2/scratch/b"
mv_out "$2/scratch/a" "$2/mafft/OG0000001_first.fa"
cp_out "$2/scratch/b" "$2/mafft/OG0000001_second.fa"
printf 'third\n' | cp_out "$2/mafft/OG0000001_third.fa"
'''
    published_root = tmp_path
    if symlink_root:
        published_root = tmp_path / "alias"
        published_root.symlink_to(tmp_path, target_is_directory=True)
    subprocess.run(["bash", "-c", script, "bash", str(repo), str(published_root), str(inventory)], check=True)
    args = store.build_parser().parse_args([
        "enqueue-family", "--root", str(tmp_path), "--mode", "orthogroup", "--family-id", family,
        "--compression", "store", "--compression-level", "2", "--workers", "4",
        "--max-final-zip-bytes", "20",
    ])
    assert store.run_cli(args) == 0
    assert drain(tmp_path)["archived_files"] == 3
    assert [path.name for path in inventory.glob("*.paths")] == ["inventory.paths"]
    for archive_path in (tmp_path / "archives" / "mafft").glob("*.zip"):
        with zipfile.ZipFile(archive_path) as archive:
            for info in archive.infolist():
                if info.filename != store.MANIFEST_MEMBER:
                    assert info.compress_type == zipfile.ZIP_STORED


def test_inventory_tracks_existing_live_and_restored_outputs(tmp_path, monkeypatch):
    path, _ = queued(tmp_path)
    assert drain(tmp_path)["archived_files"] == 1
    inventory = store.family_inventory_path(tmp_path, "OG0000001")
    for journal in inventory.glob("*.paths"):
        journal.unlink()
    monkeypatch.setenv("GG_FAMILY_OUTPUT_INVENTORY", str(inventory))
    logical = store.GeneFamilyOutputStore(tmp_path)
    logical.materialize_family("OG0000001", store.orthogroup_id_from_name)
    assert path.exists()
    assert store._inventory_paths(tmp_path, "OG0000001", store.orthogroup_id_from_name) == [path]
    for journal in inventory.glob("*.paths"):
        journal.unlink()
    logical.materialize_family("OG0000001", store.orthogroup_id_from_name)
    assert store._inventory_paths(tmp_path, "OG0000001", store.orthogroup_id_from_name) == [path]


def test_query_families_with_overlapping_names_stay_separate(tmp_path):
    query_dir = tmp_path / "queries"
    query_dir.mkdir()
    for family in ("A", "A_long"):
        (query_dir / family).write_text("gene\n")
        path = tmp_path / "mafft" / f"{family}_cds.aln.fa.gz"
        path.parent.mkdir(exist_ok=True)
        path.write_bytes(family.encode())
        inventory = store.family_inventory_path(tmp_path, family)
        inventory.mkdir(parents=True)
        (inventory / "worker.paths").write_bytes(os.fsencode(path) + b"\0")
        store.enqueue_family_archive(tmp_path, "query2family", family)
    _, matcher = store.family_context("query2family", query_dir=query_dir)
    result = store.drain_archive_queue(tmp_path, "query2family", matcher)
    assert result["archived_files"] == 2
    logical = store.GeneFamilyOutputStore(tmp_path)
    assert read_bytes(logical, "mafft", "A_cds.aln.fa.gz") == b"A"
    assert read_bytes(logical, "mafft", "A_long_cds.aln.fa.gz") == b"A_long"


def test_state_file_has_one_lock_for_all_of_its_families():
    locks_by_state_file = {}
    for i in range(10000):
        family = f"OG{i:07d}"
        locks_by_state_file.setdefault(store._family_index_bucket(family), set()).add(store._state_lock_bucket(family))
    assert len(locks_by_state_file) == 256
    assert all(len(lock_names) == 1 for lock_names in locks_by_state_file.values())


def test_queue_commit_preserves_unrelated_family_index_bucket(tmp_path):
    first, _ = queued(tmp_path, "OG0000001")
    assert drain(tmp_path)["archived_files"] == 1
    first_bucket = store._family_index_bucket("OG0000001")
    bucket_path = tmp_path / ".gg_store" / store.INDEX_DIR_NAME / f"{first_bucket}.json"
    before = bucket_path.stat().st_mtime_ns, bucket_path.read_bytes()
    family = next(f"OG{i:07d}" for i in range(2, 100) if store._family_index_bucket(f"OG{i:07d}") != first_bucket)
    second, _ = queued(tmp_path, family)
    assert drain(tmp_path)["archived_files"] == 1
    assert (bucket_path.stat().st_mtime_ns, bucket_path.read_bytes()) == before
    logical = store.GeneFamilyOutputStore(tmp_path)
    assert read_bytes(logical, "mafft", first.name) == b"alignment\n"
    assert read_bytes(logical, "mafft", second.name) == b"alignment\n"


def test_relocated_inventory_still_archives_published_outputs(tmp_path, monkeypatch):
    root = tmp_path / 'original'
    root.mkdir()
    path, _ = queued(root)
    inventory = store.family_inventory_path(root, 'OG0000001')
    for journal in inventory.glob('*.paths'):
        journal.unlink()
    monkeypatch.setenv('GG_FAMILY_OUTPUT_INVENTORY', str(inventory))
    store.record_family_output_inventory(root, 'OG0000001', path)
    store.enqueue_family_archive(root, 'orthogroup', 'OG0000001', 'run-1')
    moved = tmp_path / 'moved'
    root.rename(moved)
    assert drain(moved)['archived_files'] == 1
    assert read_bytes(store.GeneFamilyOutputStore(moved), 'mafft', path.name) == b'alignment\n'


@pytest.mark.parametrize('missing', ['directory', 'journals'])
def test_missing_inventory_never_acknowledges_request(tmp_path, missing):
    path, request = queued(tmp_path)
    inventory = store.family_inventory_path(tmp_path, 'OG0000001')
    for journal in inventory.glob('*.paths'):
        journal.unlink()
    if missing == 'directory':
        inventory.rmdir()
    with pytest.raises(store.ArchiveStoreError, match='inventory'):
        drain(tmp_path)
    assert request.exists() and path.exists()


def test_enqueue_does_not_touch_active_publisher_journal(tmp_path):
    path, _ = queued(tmp_path)
    journal = store.family_inventory_path(tmp_path, "OG0000001") / "worker.paths"
    # An open append handle must remain attached to the collected journal.
    with journal.open("ab") as publisher:
        store.enqueue_family_archive(tmp_path, "orthogroup", "OG0000001", "run-1")
        second = path.with_name("OG0000001_second.fa")
        publisher.write(os.fsencode(second.relative_to(tmp_path)) + b"\0")
        publisher.flush()
        second.write_bytes(b"second")
    assert drain(tmp_path)["archived_files"] == 2


def test_cancel_request_preserves_live_outputs(tmp_path):
    path, request = queued(tmp_path)
    args = store.build_parser().parse_args([
        "cancel-family-archive", "--root", str(tmp_path), "--family-id", "OG0000001",
    ])
    assert store.run_cli(args) == 0
    assert store.run_cli(args) == 0
    assert not request.exists()
    assert drain(tmp_path)["archived_files"] == 0
    assert path.exists()


def test_collector_defers_during_offline_maintenance(tmp_path):
    path, request = queued(tmp_path)
    with store.all_family_bucket_locks(store._archive_state_root(tmp_path)):
        assert drain(tmp_path)["status"] == "maintenance-busy"
    assert path.exists() and request.exists()


@pytest.mark.parametrize("level", [0, 1, 2])
def test_inventory_symlink_ancestors_are_rejected(tmp_path, level):
    path, request = queued(tmp_path)
    directory = store.family_inventory_path(tmp_path, "OG0000001")
    for _ in range(level):
        directory = directory.parent
    moved = tmp_path / "outside-inventory"
    directory.rename(moved)
    directory.symlink_to(moved, target_is_directory=True)
    with pytest.raises(store.ArchiveStoreError, match="Symlinked output inventory"):
        drain(tmp_path)
    assert path.exists() and request.exists()


@pytest.mark.parametrize("storage,debug,expected", [
    ("files", 0, False), ("zip", 1, False), ("zip", 0, True),
])
def test_core_startup_only_queues_normal_zip_runs(tmp_path, storage, debug, expected):
    import subprocess

    repo = Path(__file__).resolve().parents[2]
    core = (repo / "workflow/core/gg_gene_evolution_core.sh").read_text()
    block = core.split("# The run lock and shared family lock", 1)[1].split('\ndir_sp_genome=', 1)[0]
    block = "# The run lock and shared family lock" + block
    script = r'''
set -eu
root=$1
log="$root/commands"
gene_family_output_storage=$2
gg_debug_mode=$3
dir_output_active=$root
gene_family_store_script=store.py
og_id=OG0000001
mode_gene_evolution=orthogroup
gene_family_archive_write_args=()
python() {
  printf '%s\n' "$2" >> "$log"
  if [[ "$2" == inventory-path ]]; then printf '%s\n' "$root/inventory"; fi
}
''' + block
    subprocess.run(["bash", "-c", script, "bash", str(tmp_path), storage, str(debug)], check=True)
    commands = (tmp_path / "commands").read_text().splitlines()
    assert ("enqueue-family" in commands) is expected
    assert ("cancel-family-archive" in commands) is (not expected)


@pytest.mark.parametrize("entry", ["../mafft/OG0000001.fa", "/old-workspace/mafft/OG0000001.fa"])
def test_unsafe_inventory_entry_is_not_acknowledged(tmp_path, entry):
    path, request = queued(tmp_path)
    inventory = store.family_inventory_path(tmp_path, "OG0000001")
    (inventory / "worker.paths").write_bytes(os.fsencode(entry) + b"\0")
    with pytest.raises(store.ArchiveStoreError, match="inventory"):
        drain(tmp_path)
    assert path.exists() and request.exists()

import fcntl
import json
import os
import shutil
import threading
import time
import zipfile
from pathlib import Path
from types import SimpleNamespace

import pytest

from workflow.support.gene_family_output_store import (
    MAX_REFERENCED_SHARDS_PER_SUBDIR,
    ArchiveStoreError,
    Artifact,
    GeneFamilyOutputStore,
    _zip_to_raw_requirements,
    archive_completed_outputs,
    build_parser,
    cleanup_materialization_receipt,
    cleanup_stale_tmp,
    compact_archives,
    convert_storage_to_raw,
    convert_storage_to_zip,
    family_context,
    family_context_with_supplement,
    family_lock_path,
    orthogroup_id_from_name,
    purge_archives,
    repair_archive_index,
    run_cli,
    storage_conversion_summary,
)


def _write_family_outputs(root: Path, family_id: str, complete: bool = True):
    paths = {
        "mafft": root / "mafft" / f"{family_id}_cds.aln.fa.gz",
        "stat_branch": root / "stat_branch" / f"{family_id}_stat.branch.tsv",
        "stat_tree": root / "stat_tree" / f"{family_id}_stat.tree.tsv",
        "tree_plot": root / "tree_plot" / f"{family_id}_tree_plot.pdf",
    }
    for subdir, path in paths.items():
        if not complete and subdir != "mafft":
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"{subdir}:{family_id}\n".encode("utf-8"))
    return paths


def test_archives_only_completed_families_from_partially_complete_subdirs(tmp_path: Path):
    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "A").write_text("geneA\n", encoding="utf-8")
    (query_dir / "B").write_text("geneB\n", encoding="utf-8")
    paths_a = _write_family_outputs(root, "A", complete=True)
    paths_b = _write_family_outputs(root, "B", complete=False)
    family_ids, family_from_name = family_context("query2family", query_dir=query_dir)

    archived = archive_completed_outputs(
        root=root,
        mode="query2family",
        family_ids=family_ids,
        family_from_name=family_from_name,
        min_files=1,
    )

    assert archived
    assert not paths_a["mafft"].exists()
    assert paths_b["mafft"].exists()
    store = GeneFamilyOutputStore(root)
    assert store.logical_exists("mafft/A_cds.aln.fa.gz")
    assert store.logical_exists("mafft/B_cds.aln.fa.gz")
    assert store.logical_exists("tree_plot/A_tree_plot.pdf")
    assert not store.logical_exists("tree_plot/B_tree_plot.pdf")


def test_empty_terminal_artifact_does_not_mark_family_complete(tmp_path: Path):
    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "A").write_text("geneA\n", encoding="utf-8")
    paths = _write_family_outputs(root, "A", complete=True)
    paths["tree_plot"].write_bytes(b"")
    family_ids, family_from_name = family_context("query2family", query_dir=query_dir)

    archived = archive_completed_outputs(
        root,
        "query2family",
        family_ids,
        family_from_name,
    )

    assert archived == []
    assert paths["stat_branch"].is_file()


def test_query_addition_creates_a_later_shard_without_rewriting_old_shards(tmp_path: Path):
    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "A").write_text("geneA\n", encoding="utf-8")
    _write_family_outputs(root, "A", complete=True)
    family_ids, family_from_name = family_context("query2family", query_dir=query_dir)
    archive_completed_outputs(root, "query2family", family_ids, family_from_name)
    old_shards = set((root / ".gg_archives").glob("*/*.zip"))

    (query_dir / "B").write_text("geneB\n", encoding="utf-8")
    _write_family_outputs(root, "B", complete=True)
    family_ids, family_from_name = family_context("query2family", query_dir=query_dir)
    archive_completed_outputs(root, "query2family", family_ids, family_from_name)
    new_shards = set((root / ".gg_archives").glob("*/*.zip"))

    assert old_shards
    assert old_shards < new_shards
    store = GeneFamilyOutputStore(root)
    assert store.logical_exists("stat_branch/A_stat.branch.tsv")
    assert store.logical_exists("stat_branch/B_stat.branch.tsv")


def test_live_file_overrides_archive_and_delete_uses_tombstone(tmp_path: Path):
    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "A").write_text("geneA\n", encoding="utf-8")
    _write_family_outputs(root, "A", complete=True)
    family_ids, family_from_name = family_context("query2family", query_dir=query_dir)
    archive_completed_outputs(root, "query2family", family_ids, family_from_name)

    live_path = root / "stat_branch" / "A_stat.branch.tsv"
    live_path.parent.mkdir(parents=True, exist_ok=True)
    live_path.write_text("new-live-version\n", encoding="utf-8")
    store = GeneFamilyOutputStore(root)
    with store.open_binary("stat_branch", live_path.name) as handle:
        assert handle.read() == b"new-live-version\n"

    store.delete("stat_branch/A_stat.branch.tsv")
    refreshed = GeneFamilyOutputStore(root)
    assert not refreshed.logical_exists("stat_branch/A_stat.branch.tsv")

    refreshed.undelete("stat_branch/A_stat.branch.tsv")
    restored_store = GeneFamilyOutputStore(root)
    assert restored_store.logical_exists("stat_branch/A_stat.branch.tsv")


def test_tombstone_log_is_compacted_to_latest_records(tmp_path: Path, monkeypatch):
    import workflow.support.gene_family_output_store as output_store_module

    root = tmp_path / "orthogroup"
    monkeypatch.setattr(output_store_module, "TOMBSTONE_LOG_COMPACT_BYTES", 1)
    store = GeneFamilyOutputStore(root)
    logical_path = "stat_branch/OG0000001_stat.branch.tsv"
    store.delete(logical_path)
    store.undelete(logical_path)
    store.delete(logical_path)

    lines = (
        root / ".gg_archives" / "tombstones.jsonl"
    ).read_text(encoding="utf-8").splitlines()
    assert len(lines) == 1
    assert '"operation":"delete"' in lines[0]


def test_materialize_family_restores_only_requested_family(tmp_path: Path):
    root = tmp_path / "orthogroup"
    genecount = tmp_path / "Orthogroups.GeneCount.selected.tsv"
    genecount.write_text(
        "Orthogroup\tTotal\nHOG0000001\t1\nHOG0000002\t1\n",
        encoding="utf-8",
    )
    _write_family_outputs(root, "HOG0000001", complete=True)
    _write_family_outputs(root, "HOG0000002", complete=True)
    family_ids, family_from_name = family_context("orthogroup", genecount=genecount)
    archive_completed_outputs(root, "orthogroup", family_ids, family_from_name)

    store = GeneFamilyOutputStore(root)
    restored = store.materialize_family("HOG0000001", family_from_name)

    assert restored
    assert (root / "mafft" / "HOG0000001_cds.aln.fa.gz").is_file()
    assert not (root / "mafft" / "HOG0000002_cds.aln.fa.gz").exists()


def test_materialize_family_cli_does_not_require_full_family_catalog(tmp_path: Path):
    root = tmp_path / "orthogroup"
    genecount = tmp_path / "Orthogroups.GeneCount.selected.tsv"
    genecount.write_text("Orthogroup\tTotal\nHOG0000001\t1\n", encoding="utf-8")
    _write_family_outputs(root, "HOG0000001", complete=True)
    family_ids, family_from_name = family_context("orthogroup", genecount=genecount)
    archive_completed_outputs(root, "orthogroup", family_ids, family_from_name)
    destination = tmp_path / "selected"

    assert run_cli(
        SimpleNamespace(
            command="materialize-family",
            root=root,
            mode="orthogroup",
            family_id="HOG0000001",
            destination_root=destination,
            subdirs="stat_branch",
        )
    ) == 0

    assert (destination / "stat_branch" / "HOG0000001_stat.branch.tsv").is_file()


def test_materialize_families_cli_restores_selected_subdir_in_one_pass(tmp_path: Path):
    root = tmp_path / "orthogroup"
    for family_id in ("OG0000001", "OG0000002", "OG0000003"):
        _write_family_outputs(root, family_id, complete=True)
    archive_completed_outputs(
        root,
        "orthogroup",
        ["OG0000001", "OG0000002", "OG0000003"],
        lambda name: name.split("_", 1)[0],
    )
    family_id_file = tmp_path / "selected.txt"
    family_id_file.write_text("OG0000001\nOG0000003\n", encoding="utf-8")
    destination = tmp_path / "selected"

    assert run_cli(
        SimpleNamespace(
            command="materialize-families",
            root=root,
            mode="orthogroup",
            query_dir=None,
            family_id_file=family_id_file,
            destination_root=destination,
            subdirs="mafft",
        )
    ) == 0

    assert (destination / "mafft" / "OG0000001_cds.aln.fa.gz").is_file()
    assert not (destination / "mafft" / "OG0000002_cds.aln.fa.gz").exists()
    assert (destination / "mafft" / "OG0000003_cds.aln.fa.gz").is_file()


def test_verify_rejects_directly_modified_zip(tmp_path: Path):
    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "A").write_text("geneA\n", encoding="utf-8")
    _write_family_outputs(root, "A", complete=True)
    family_ids, family_from_name = family_context("query2family", query_dir=query_dir)
    archive_completed_outputs(root, "query2family", family_ids, family_from_name)
    zip_path = next((root / ".gg_archives").glob("*/*.zip"))

    with zipfile.ZipFile(zip_path, "a") as archive:
        archive.writestr("manually-added.txt", "unsupported")

    with pytest.raises(ArchiveStoreError, match="differs from its manifest"):
        GeneFamilyOutputStore(root).verify()


def test_verify_streams_each_manifest_once_without_testzip(
    tmp_path: Path,
    monkeypatch,
):
    root = tmp_path / "orthogroup"
    _write_family_outputs(root, "OG0000001", complete=True)
    archive_completed_outputs(
        root,
        "orthogroup",
        ["OG0000001"],
        lambda name: "OG0000001" if name.startswith("OG0000001_") else None,
    )
    store = GeneFamilyOutputStore(root)
    original_read_manifest = store._read_manifest

    class SinglePassMembers(list):
        iterated = False

        def __iter__(self):
            if self.iterated:
                raise AssertionError("manifest members were iterated more than once")
            self.iterated = True
            return super().__iter__()

    def read_manifest_once(zip_path):
        manifest = original_read_manifest(zip_path)
        manifest["members"] = SinglePassMembers(manifest["members"])
        return manifest

    monkeypatch.setattr(store, "_read_manifest", read_manifest_once)
    monkeypatch.setattr(
        zipfile.ZipFile,
        "testzip",
        lambda self: (_ for _ in ()).throw(AssertionError("testzip was called")),
    )

    assert store.verify()


def test_list_cli_reports_live_override_and_archived_members(tmp_path: Path, capsys):
    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "A").write_text("geneA\n", encoding="utf-8")
    _write_family_outputs(root, "A", complete=True)
    family_ids, family_from_name = family_context("query2family", query_dir=query_dir)
    archive_completed_outputs(root, "query2family", family_ids, family_from_name)
    live = root / "stat_branch" / "A_stat.branch.tsv"
    live.parent.mkdir(parents=True)
    live.write_text("replacement\n", encoding="utf-8")

    assert run_cli(SimpleNamespace(command="list", root=root, subdir="stat_branch")) == 0

    fields = capsys.readouterr().out.strip().split("\t")
    assert fields[0] == "stat_branch/A_stat.branch.tsv"
    assert fields[1] == str(len("replacement\n"))
    assert fields[2] == "live"

    assert run_cli(
        SimpleNamespace(
            command="has-files",
            root=root,
            subdir="stat_branch",
            suffix="_stat.branch.tsv",
        )
    ) == 0


def test_materialize_preserves_source_mtime_for_staleness_checks(tmp_path: Path):
    root = tmp_path / "orthogroup"
    paths = _write_family_outputs(root, "OG0000001", complete=True)
    source_mtime_ns = paths["stat_branch"].stat().st_mtime_ns

    archive_completed_outputs(
        root,
        "orthogroup",
        ["OG0000001"],
        lambda name: "OG0000001" if name.startswith("OG0000001_") else None,
    )
    time.sleep(0.01)
    restored = GeneFamilyOutputStore(
        root,
        family_filter="OG0000001",
    ).materialize("stat_branch", "OG0000001_stat.branch.tsv")

    assert restored.stat().st_mtime_ns == source_mtime_ns


def test_explicit_complete_state_allows_archive_without_optional_tree_outputs(tmp_path: Path):
    root = tmp_path / "orthogroup"
    mafft = root / "mafft" / "OG0000001_cds.aln.fa.gz"
    mafft.parent.mkdir(parents=True)
    mafft.write_bytes(b"alignment")
    store = GeneFamilyOutputStore(root, family_filter="OG0000001")
    store.mark_family_state("OG0000001", "running", "run-1")

    assert archive_completed_outputs(
        root,
        "orthogroup",
        ["OG0000001"],
        lambda name: "OG0000001" if name.startswith("OG0000001_") else None,
    ) == []
    assert mafft.is_file()

    GeneFamilyOutputStore(root, family_filter="OG0000001").mark_family_state(
        "OG0000001",
        "complete",
        "run-1",
    )
    assert archive_completed_outputs(
        root,
        "orthogroup",
        ["OG0000001"],
        lambda name: "OG0000001" if name.startswith("OG0000001_") else None,
    )
    assert not mafft.exists()


def test_stale_run_cannot_overwrite_newer_family_state(tmp_path: Path):
    root = tmp_path / "orthogroup"
    store = GeneFamilyOutputStore(root, family_filter="OG0000001")
    assert store.mark_family_state("OG0000001", "running", "older")
    assert store.mark_family_state("OG0000001", "running", "newer")
    assert not store.mark_family_state("OG0000001", "failed", "older")
    assert store.family_state("OG0000001") == "running"
    assert store.mark_family_state("OG0000001", "complete", "newer")
    assert not store.mark_family_state("OG0000001", "complete", "older")
    assert store.family_state("OG0000001") == "complete"


def test_family_filtered_store_refreshes_current_state(tmp_path: Path):
    root = tmp_path / "orthogroup"
    family_id = "OG0000001"
    existing = GeneFamilyOutputStore(root, family_filter=family_id)
    assert existing.mark_family_state(family_id, "complete", "first")
    assert existing.family_state(family_id) == "complete"

    assert GeneFamilyOutputStore(
        root,
        family_filter=family_id,
    ).mark_family_state(family_id, "running", "second")

    assert existing.family_state(family_id) == "running"


def test_active_family_shared_lock_prevents_nonblocking_archive(tmp_path: Path):
    root = tmp_path / "orthogroup"
    paths = _write_family_outputs(root, "OG0000001", complete=True)
    archive_root = root / ".gg_archives"
    lock_path = family_lock_path(archive_root, "OG0000001")
    lock_path.parent.mkdir(parents=True)
    lock_handle = lock_path.open("a+b")
    fcntl.flock(lock_handle.fileno(), fcntl.LOCK_SH)
    try:
        assert archive_completed_outputs(
            root,
            "orthogroup",
            ["OG0000001"],
            lambda name: "OG0000001" if name.startswith("OG0000001_") else None,
            nonblocking=True,
        ) == []
        assert paths["stat_branch"].is_file()
    finally:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
        lock_handle.close()


def test_active_family_does_not_block_archiving_an_unrelated_family(tmp_path: Path):
    root = tmp_path / "orthogroup"
    active_family = "OG0000001"
    completed_family = "OG0000002"
    while family_lock_path(
        root / ".gg_archives",
        active_family,
    ) == family_lock_path(root / ".gg_archives", completed_family):
        completed_family = f"OG{int(completed_family[2:]) + 1:07d}"
    active_paths = _write_family_outputs(root, active_family, complete=True)
    completed_paths = _write_family_outputs(root, completed_family, complete=True)
    archive_root = root / ".gg_archives"
    lock_path = family_lock_path(archive_root, active_family)
    lock_path.parent.mkdir(parents=True)
    lock_handle = lock_path.open("a+b")
    fcntl.flock(lock_handle.fileno(), fcntl.LOCK_SH)
    try:
        archived = archive_completed_outputs(
            root,
            "orthogroup",
            [active_family, completed_family],
            lambda name: name.split("_", 1)[0],
            nonblocking=True,
        )
    finally:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
        lock_handle.close()

    assert archived
    assert active_paths["stat_branch"].is_file()
    assert not completed_paths["stat_branch"].exists()
    store = GeneFamilyOutputStore(root)
    assert store.logical_exists(
        f"stat_branch/{completed_family}_stat.branch.tsv"
    )


def test_managed_delete_waits_for_the_active_family_lock(tmp_path: Path):
    root = tmp_path / "orthogroup"
    family_id = "OG0000001"
    live_path = root / "stat_branch" / f"{family_id}_stat.branch.tsv"
    live_path.parent.mkdir(parents=True)
    live_path.write_text("active output\n", encoding="utf-8")
    lock_path = family_lock_path(root / ".gg_archives", family_id)
    lock_path.parent.mkdir(parents=True)
    lock_handle = lock_path.open("a+b")
    fcntl.flock(lock_handle.fileno(), fcntl.LOCK_SH)
    started = threading.Event()
    finished = threading.Event()
    failures = []

    def delete_in_thread():
        started.set()
        try:
            GeneFamilyOutputStore(root).delete(
                f"stat_branch/{live_path.name}",
                family_id=family_id,
            )
        except Exception as exc:  # pragma: no cover - asserted below
            failures.append(exc)
        finally:
            finished.set()

    worker = threading.Thread(target=delete_in_thread)
    worker.start()
    assert started.wait(timeout=1)
    time.sleep(0.05)
    assert not finished.is_set()
    assert live_path.is_file()
    fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
    lock_handle.close()
    worker.join(timeout=2)

    assert finished.is_set()
    assert failures == []
    assert not live_path.exists()


@pytest.mark.parametrize("operation", ["delete", "undelete", "restore"])
def test_managed_query_operation_requires_unknown_family_id(
    tmp_path: Path,
    operation: str,
):
    root = tmp_path / "query2family"
    logical_path = "stat_branch/query_name_stat.branch.tsv"
    live_path = root / logical_path
    live_path.parent.mkdir(parents=True)
    live_path.write_text("manual output\n", encoding="utf-8")
    store = GeneFamilyOutputStore(root)

    with pytest.raises(ArchiveStoreError, match="supply --family-id explicitly"):
        getattr(store, operation)(logical_path)

    assert live_path.is_file()


def test_managed_query_delete_accepts_explicit_unknown_family_id(tmp_path: Path):
    root = tmp_path / "query2family"
    logical_path = "stat_branch/query_name_stat.branch.tsv"
    live_path = root / logical_path
    live_path.parent.mkdir(parents=True)
    live_path.write_text("manual output\n", encoding="utf-8")

    GeneFamilyOutputStore(root).delete(
        logical_path,
        family_id="query_name",
    )

    assert not live_path.exists()
    assert not GeneFamilyOutputStore(root).logical_exists(logical_path)


def test_managed_query_delete_rejects_explicit_id_that_does_not_match_filename(
    tmp_path: Path,
):
    root = tmp_path / "query2family"
    logical_path = "stat_branch/query_name_stat.branch.tsv"
    live_path = root / logical_path
    live_path.parent.mkdir(parents=True)
    live_path.write_text("manual output\n", encoding="utf-8")

    with pytest.raises(ArchiveStoreError, match="does not match the logical filename"):
        GeneFamilyOutputStore(root).delete(
            logical_path,
            family_id="another_family",
        )

    assert live_path.is_file()


@pytest.mark.parametrize("operation", ["delete", "undelete", "restore"])
def test_managed_operation_rejects_family_id_different_from_archive_metadata(
    tmp_path: Path,
    operation: str,
):
    root = tmp_path / "query2family"
    family_id = "query_name"
    logical_path = f"stat_branch/{family_id}_stat.branch.tsv"
    _write_family_outputs(root, family_id, complete=True)
    archive_completed_outputs(
        root,
        "query2family",
        [family_id],
        lambda name: family_id if name.startswith(family_id + "_") else None,
        include_incomplete=True,
    )
    store = GeneFamilyOutputStore(root)

    with pytest.raises(ArchiveStoreError, match="differs from the inferred family ID"):
        getattr(store, operation)(logical_path, family_id="another_family")

    assert GeneFamilyOutputStore(root).logical_exists(logical_path)


def test_active_store_reader_prevents_nonblocking_archive(tmp_path: Path):
    root = tmp_path / "orthogroup"
    _write_family_outputs(root, "OG0000001", complete=True)
    archive_completed_outputs(
        root,
        "orthogroup",
        ["OG0000001"],
        lambda name: "OG0000001" if name.startswith("OG0000001_") else None,
    )
    next_path = root / "mafft" / "OG0000002_cds.aln.fa.gz"
    next_path.parent.mkdir(parents=True, exist_ok=True)
    next_path.write_bytes(b"next")
    GeneFamilyOutputStore(root, family_filter="OG0000002").mark_family_state(
        "OG0000002",
        "complete",
        "",
    )

    store = GeneFamilyOutputStore(root, family_filter="OG0000001")
    with store.open_binary("mafft", "OG0000001_cds.aln.fa.gz") as handle:
        assert handle.read()
        assert archive_completed_outputs(
            root,
            "orthogroup",
            ["OG0000002"],
            lambda name: "OG0000002" if name.startswith("OG0000002_") else None,
            nonblocking=True,
        ) == []
        assert next_path.is_file()

    assert archive_completed_outputs(
        root,
        "orthogroup",
        ["OG0000002"],
        lambda name: "OG0000002" if name.startswith("OG0000002_") else None,
        nonblocking=True,
    )


def test_query_prefix_family_metadata_prevents_cross_family_materialization(tmp_path: Path):
    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "q").write_text("short\n", encoding="utf-8")
    (query_dir / "q_long").write_text("long\n", encoding="utf-8")
    _write_family_outputs(root, "q_long", complete=True)
    family_ids, family_from_name = family_context("query2family", query_dir=query_dir)
    archive_completed_outputs(
        root,
        "query2family",
        family_ids,
        family_from_name,
    )

    destination = tmp_path / "selected"
    restored = GeneFamilyOutputStore(root, family_filter="q").materialize_family(
        "q",
        lambda name: "q" if name == "q" or name.startswith(("q_", "q.")) else None,
        destination_root=destination,
    )

    assert restored == []
    assert not destination.exists()


def test_incremental_archives_are_compacted_to_bounded_shard_count(tmp_path: Path):
    root = tmp_path / "orthogroup"
    subdirs = {
        "stat_branch": "_stat.branch.tsv",
        "stat_tree": "_stat.tree.tsv",
        "tree_plot": "_tree_plot.pdf",
        "mafft": "_cds.aln.fa.gz",
    }
    for index in range(20):
        family_id = f"OG{index:07d}"
        for subdir, suffix in subdirs.items():
            path = root / subdir / f"{family_id}{suffix}"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"x")
        archive_completed_outputs(
            root,
            "orthogroup",
            [family_id],
            lambda name, family_id=family_id: (
                family_id if name.startswith(f"{family_id}_") else None
            ),
        )

    shard_count = len(list((root / ".gg_archives").glob("*/*.zip")))
    assert shard_count <= len(subdirs) * (MAX_REFERENCED_SHARDS_PER_SUBDIR - 1)
    assert sum(
        len(GeneFamilyOutputStore(root).file_names(subdir))
        for subdir in subdirs
    ) == 20 * len(subdirs)


def test_required_multi_shard_store_allows_only_bounded_fragmentation(tmp_path: Path):
    root = tmp_path / "orthogroup"
    for index in range(30):
        family_id = f"OG{index:07d}"
        path = root / "mafft" / f"{family_id}_cds.aln.fa.gz"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(family_id.encode("utf-8"))
        store = GeneFamilyOutputStore(root, family_filter=family_id)
        store.mark_family_state(family_id, "complete", "")
        archive_completed_outputs(
            root,
            "orthogroup",
            [family_id],
            lambda name, family_id=family_id: (
                family_id if name.startswith(f"{family_id}_") else None
            ),
            max_files_per_shard=2,
        )

    shard_count = len(list((root / ".gg_archives" / "mafft").glob("*.zip")))
    minimum_shards = 15
    assert minimum_shards <= shard_count <= (
        minimum_shards + MAX_REFERENCED_SHARDS_PER_SUBDIR - 1
    )
    assert len(GeneFamilyOutputStore(root).file_names("mafft")) == 30


def test_family_scoped_parameters_are_archived(tmp_path: Path):
    root = tmp_path / "orthogroup"
    _write_family_outputs(root, "OG0000001", complete=True)
    parameter = root / "parameters" / "OG0000001_species_genetic_code.resolved.tsv"
    parameter.parent.mkdir(parents=True)
    parameter.write_text("species\tcode\n", encoding="utf-8")

    archive_completed_outputs(
        root,
        "orthogroup",
        ["OG0000001"],
        lambda name: "OG0000001" if name.startswith("OG0000001_") else None,
    )

    assert not parameter.exists()
    assert GeneFamilyOutputStore(root).logical_exists(
        "parameters/OG0000001_species_genetic_code.resolved.tsv"
    )


def test_export_current_applies_live_overrides_and_tombstones(tmp_path: Path):
    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "A").write_text("geneA\n", encoding="utf-8")
    _write_family_outputs(root, "A", complete=True)
    family_ids, family_from_name = family_context("query2family", query_dir=query_dir)
    archive_completed_outputs(root, "query2family", family_ids, family_from_name)
    store = GeneFamilyOutputStore(root)
    store.delete("tree_plot/A_tree_plot.pdf")
    live = root / "stat_branch" / "A_stat.branch.tsv"
    live.parent.mkdir(parents=True, exist_ok=True)
    live.write_text("replacement\n", encoding="utf-8")

    destination = tmp_path / "export"
    store = GeneFamilyOutputStore(root)
    store.materialize_all(destination)

    assert (destination / "stat_branch" / "A_stat.branch.tsv").read_text(
        encoding="utf-8"
    ) == "replacement\n"
    assert not (destination / "tree_plot" / "A_tree_plot.pdf").exists()


def test_export_current_cli_requires_a_clean_destination(tmp_path: Path):
    root = tmp_path / "orthogroup"
    _write_family_outputs(root, "OG0000001", complete=True)
    destination = tmp_path / "export"
    destination.mkdir()
    (destination / "stale.txt").write_text("stale", encoding="utf-8")

    with pytest.raises(ValueError, match="absent or empty"):
        run_cli(
            SimpleNamespace(
                command="export-current",
                root=root,
                destination_root=destination,
            )
        )


def test_verify_rejects_index_reference_to_missing_zip(tmp_path: Path):
    root = tmp_path / "orthogroup"
    _write_family_outputs(root, "OG0000001", complete=True)
    archive_completed_outputs(
        root,
        "orthogroup",
        ["OG0000001"],
        lambda name: "OG0000001" if name.startswith("OG0000001_") else None,
    )
    zip_path = next((root / ".gg_archives").glob("*/*.zip"))
    zip_path.unlink()

    with pytest.raises(ArchiveStoreError, match="references a missing ZIP"):
        GeneFamilyOutputStore(root).verify()


def test_verify_rejects_manifest_members_missing_from_index(tmp_path: Path):
    root = tmp_path / "orthogroup"
    _write_family_outputs(root, "OG0000001", complete=True)
    archive_completed_outputs(
        root,
        "orthogroup",
        ["OG0000001"],
        lambda name: "OG0000001" if name.startswith("OG0000001_") else None,
    )
    for index_path in (root / ".gg_archives" / "index").glob("*.json"):
        index_path.unlink()

    with pytest.raises(ArchiveStoreError, match="absent from the archive index"):
        GeneFamilyOutputStore(root).verify()


def test_failed_restore_preserves_prior_logical_deletion(tmp_path: Path):
    root = tmp_path / "orthogroup"
    _write_family_outputs(root, "OG0000001", complete=True)
    archive_completed_outputs(
        root,
        "orthogroup",
        ["OG0000001"],
        lambda name: "OG0000001" if name.startswith("OG0000001_") else None,
    )
    logical_path = "stat_branch/OG0000001_stat.branch.tsv"
    store = GeneFamilyOutputStore(root)
    store.delete(logical_path)
    zip_path = next((root / ".gg_archives").glob("stat_branch/*.zip"))
    zip_path.unlink()

    with pytest.raises(FileNotFoundError):
        GeneFamilyOutputStore(root).restore(logical_path)
    assert not GeneFamilyOutputStore(root).logical_exists(logical_path)


def test_stale_tmp_cleanup_removes_only_old_task_directories(tmp_path: Path):
    root = tmp_path / "orthogroup"
    old_dir = root / "tmp" / "1_OG0000001"
    recent_dir = root / "tmp" / "2_OG0000002"
    old_file = old_dir / "work.txt"
    recent_file = recent_dir / "work.txt"
    old_file.parent.mkdir(parents=True)
    recent_file.parent.mkdir(parents=True)
    old_file.write_text("old", encoding="utf-8")
    recent_file.write_text("recent", encoding="utf-8")
    old_ns = time.time_ns() - 10 * 86400 * 1_000_000_000
    os.utime(old_file, ns=(old_ns, old_ns))
    os.utime(old_dir, ns=(old_ns, old_ns))

    removed = cleanup_stale_tmp(root, older_than_days=7)

    assert removed == [old_dir]
    assert not old_dir.exists()
    assert recent_dir.is_dir()


def test_stale_tmp_cleanup_skips_only_the_active_family(tmp_path: Path):
    root = tmp_path / "orthogroup"
    unrelated_family = "OG0000002"
    while family_lock_path(
        root / ".gg_archives",
        "OG0000001",
    ) == family_lock_path(root / ".gg_archives", unrelated_family):
        unrelated_family = f"OG{int(unrelated_family[2:]) + 1:07d}"
    old_dir = root / "tmp" / "1_OG0000001"
    unrelated_dir = root / "tmp" / f"2_{unrelated_family}"
    old_dir.mkdir(parents=True)
    unrelated_dir.mkdir(parents=True)
    old_ns = time.time_ns() - 10 * 86400 * 1_000_000_000
    os.utime(old_dir, ns=(old_ns, old_ns))
    os.utime(unrelated_dir, ns=(old_ns, old_ns))
    archive_root = root / ".gg_archives"
    lock_path = family_lock_path(archive_root, "OG0000001")
    lock_path.parent.mkdir(parents=True)
    lock_handle = lock_path.open("a+b")
    fcntl.flock(lock_handle.fileno(), fcntl.LOCK_SH)
    try:
        removed = cleanup_stale_tmp(
            root,
            older_than_days=7,
            nonblocking=True,
        )
        assert removed == [unrelated_dir]
        assert old_dir.is_dir()
        assert not unrelated_dir.exists()
    finally:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
        lock_handle.close()


def test_tmp_cleanup_caps_failed_directories_and_ignores_manual_directories(
    tmp_path: Path,
):
    root = tmp_path / "orthogroup"
    for index in range(5):
        candidate = root / "tmp" / f"{index + 1}_OG{index:07d}"
        candidate.mkdir(parents=True)
        mtime_ns = time.time_ns() - (10 - index) * 1_000_000_000
        os.utime(candidate, ns=(mtime_ns, mtime_ns))
    manual = root / "tmp" / "manual-notes"
    manual.mkdir()
    old_ns = time.time_ns() - 30 * 86400 * 1_000_000_000
    os.utime(manual, ns=(old_ns, old_ns))

    removed = cleanup_stale_tmp(
        root,
        older_than_days=0,
        max_directories=2,
    )

    assert {path.name for path in removed} == {
        "1_OG0000000",
        "2_OG0000001",
        "3_OG0000002",
    }
    assert (root / "tmp" / "4_OG0000003").is_dir()
    assert (root / "tmp" / "5_OG0000004").is_dir()
    assert manual.is_dir()


def test_materialization_receipt_removes_only_unchanged_archived_inputs(
    tmp_path: Path,
):
    root = tmp_path / "orthogroup"
    family_id = "OG0000001"
    _write_family_outputs(root, family_id, complete=True)
    archive_completed_outputs(
        root,
        "orthogroup",
        [family_id],
        lambda name: family_id if name.startswith(f"{family_id}_") else None,
    )
    receipt = root / "tmp" / f"1_{family_id}" / ".gg_materialized.jsonl"

    restored = GeneFamilyOutputStore(root, family_filter=family_id).materialize_family(
        family_id,
        lambda name: family_id if name.startswith(f"{family_id}_") else None,
        receipt_path=receipt,
        run_token="run-1",
    )
    modified = root / "stat_branch" / f"{family_id}_stat.branch.tsv"
    modified.write_text("new failed-run output\n", encoding="utf-8")

    removed = cleanup_materialization_receipt(receipt)

    assert modified.is_file()
    assert modified not in removed
    assert not receipt.exists()
    assert all(
        not path.exists()
        for path in restored
        if path != modified
    )


def test_materialization_receipt_cleans_all_unchanged_inputs(tmp_path: Path):
    root = tmp_path / "orthogroup"
    family_id = "OG0000001"
    _write_family_outputs(root, family_id, complete=True)
    archive_completed_outputs(
        root,
        "orthogroup",
        [family_id],
        lambda name: family_id if name.startswith(f"{family_id}_") else None,
    )
    receipt = root / "tmp" / f"1_{family_id}" / ".gg_materialized.jsonl"
    restored = GeneFamilyOutputStore(root, family_filter=family_id).materialize_family(
        family_id,
        lambda name: family_id if name.startswith(f"{family_id}_") else None,
        receipt_path=receipt,
    )

    assert cleanup_materialization_receipt(receipt) == restored
    assert not receipt.exists()
    assert all(not path.exists() for path in restored)


def test_receipt_is_durable_before_materialization_begins(
    tmp_path: Path,
    monkeypatch,
):
    root = tmp_path / "orthogroup"
    family_id = "OG0000001"
    _write_family_outputs(root, family_id, complete=True)
    archive_completed_outputs(
        root,
        "orthogroup",
        [family_id],
        lambda name: family_id if name.startswith(f"{family_id}_") else None,
    )
    receipt = root / "tmp" / f"1_{family_id}" / ".gg_materialized.jsonl"
    store = GeneFamilyOutputStore(root, family_filter=family_id)

    def fail_materialize(*args, **kwargs):
        raise OSError("simulated materialization failure")

    monkeypatch.setattr(store, "materialize", fail_materialize)
    with pytest.raises(OSError, match="simulated materialization failure"):
        store.materialize_family(
            family_id,
            lambda name: family_id if name.startswith(f"{family_id}_") else None,
            receipt_path=receipt,
        )

    assert receipt.is_file()
    assert cleanup_materialization_receipt(receipt) == []
    assert not receipt.exists()


def test_existing_store_refreshes_after_archive_replacement(tmp_path: Path):
    root = tmp_path / "orthogroup"
    family_id = "OG0000001"
    paths = _write_family_outputs(root, family_id, complete=True)
    archive_completed_outputs(
        root,
        "orthogroup",
        [family_id],
        lambda name: family_id if name.startswith(f"{family_id}_") else None,
    )
    store = GeneFamilyOutputStore(root, family_filter=family_id)
    with store.open_binary("stat_branch", paths["stat_branch"].name) as handle:
        assert handle.read() == f"stat_branch:{family_id}\n".encode()

    replacement = root / "stat_branch" / paths["stat_branch"].name
    replacement.parent.mkdir(parents=True, exist_ok=True)
    replacement.write_bytes(b"replacement\n")
    GeneFamilyOutputStore(root, family_filter=family_id).mark_family_state(
        family_id,
        "complete",
    )
    archive_completed_outputs(
        root,
        "orthogroup",
        [family_id],
        lambda name: family_id if name.startswith(f"{family_id}_") else None,
    )

    with store.open_binary("stat_branch", replacement.name) as handle:
        assert handle.read() == b"replacement\n"


def test_existing_store_refreshes_after_tombstone_change(tmp_path: Path):
    root = tmp_path / "orthogroup"
    family_id = "OG0000001"
    _write_family_outputs(root, family_id, complete=True)
    archive_completed_outputs(
        root,
        "orthogroup",
        [family_id],
        lambda name: family_id if name.startswith(f"{family_id}_") else None,
    )
    logical_path = f"stat_branch/{family_id}_stat.branch.tsv"
    existing = GeneFamilyOutputStore(root)
    assert existing.logical_exists(logical_path)

    GeneFamilyOutputStore(root).delete(logical_path)

    assert not existing.logical_exists(logical_path)


def test_verify_rejects_an_unreferenced_duplicate_zip(tmp_path: Path):
    root = tmp_path / "orthogroup"
    family_id = "OG0000001"
    _write_family_outputs(root, family_id, complete=True)
    archive_completed_outputs(
        root,
        "orthogroup",
        [family_id],
        lambda name: family_id if name.startswith(f"{family_id}_") else None,
    )
    source = next((root / ".gg_archives" / "mafft").glob("*.zip"))
    duplicate = source.with_name(f"orphan-{source.name}")
    shutil.copy2(source, duplicate)

    with pytest.raises(ArchiveStoreError, match="Unreferenced ZIP shards"):
        GeneFamilyOutputStore(root).verify()


def test_repair_rebuilds_missing_index_and_can_remove_orphans(tmp_path: Path):
    root = tmp_path / "orthogroup"
    family_id = "OG0000001"
    _write_family_outputs(root, family_id, complete=True)
    archive_completed_outputs(
        root,
        "orthogroup",
        [family_id],
        lambda name: family_id if name.startswith(f"{family_id}_") else None,
    )
    source = next((root / ".gg_archives" / "mafft").glob("*.zip"))
    duplicate = source.with_name(f"orphan-{source.name}")
    shutil.copy2(source, duplicate)
    for index_path in (root / ".gg_archives" / "index").glob("*.json"):
        index_path.unlink()

    orphaned = repair_archive_index(root, remove_orphans=True)

    assert len(orphaned) == 1
    assert not orphaned[0].exists()
    assert GeneFamilyOutputStore(root).verify()
    assert GeneFamilyOutputStore(root).logical_exists(
        f"mafft/{family_id}_cds.aln.fa.gz"
    )


def test_purge_physically_applies_tombstones_live_overrides_and_family_filter(
    tmp_path: Path,
):
    root = tmp_path / "orthogroup"
    retained_family = "OG0000001"
    dropped_family = "OG0000002"
    for family_id in (retained_family, dropped_family):
        _write_family_outputs(root, family_id, complete=True)
        GeneFamilyOutputStore(root, family_filter=family_id).mark_family_state(
            family_id,
            "complete",
        )
    archive_completed_outputs(
        root,
        "orthogroup",
        [retained_family, dropped_family],
        lambda name: name.split("_", 1)[0],
    )
    store = GeneFamilyOutputStore(root)
    deleted_path = f"stat_branch/{retained_family}_stat.branch.tsv"
    store.delete(deleted_path)
    live_override = root / "mafft" / f"{retained_family}_cds.aln.fa.gz"
    live_override.parent.mkdir(parents=True, exist_ok=True)
    live_override.write_bytes(b"manual replacement\n")

    purge_archives(
        root,
        "orthogroup",
        valid_family_ids={retained_family},
        family_from_name=lambda name: name.split("_", 1)[0],
    )

    refreshed = GeneFamilyOutputStore(root)
    assert not refreshed.logical_exists(deleted_path)
    with refreshed.open_binary("mafft", live_override.name) as handle:
        assert handle.read() == b"manual replacement\n"
    assert not refreshed.logical_exists(
        f"tree_plot/{dropped_family}_tree_plot.pdf"
    )
    assert GeneFamilyOutputStore(
        root,
        family_filter=retained_family,
    ).family_state(retained_family) == "complete"
    assert GeneFamilyOutputStore(
        root,
        family_filter=dropped_family,
    ).family_state(dropped_family) is None
    assert not (root / ".gg_archives" / "tombstones.jsonl").exists()
    refreshed.verify()


def test_tmp_cleanup_enforces_total_byte_and_file_limits(tmp_path: Path):
    root = tmp_path / "orthogroup"
    candidates = []
    now_ns = time.time_ns()
    for index, size in enumerate((3, 4, 5), start=1):
        candidate = root / "tmp" / f"{index}_OG{index:07d}"
        candidate.mkdir(parents=True)
        (candidate / "work.bin").write_bytes(b"x" * size)
        mtime_ns = now_ns - (4 - index) * 1_000_000_000
        os.utime(candidate / "work.bin", ns=(mtime_ns, mtime_ns))
        os.utime(candidate, ns=(mtime_ns, mtime_ns))
        candidates.append(candidate)

    removed = cleanup_stale_tmp(
        root,
        older_than_days=0,
        max_bytes=8,
    )

    assert removed == [candidates[1]]
    assert candidates[0].is_dir()
    assert candidates[2].is_dir()

    removed = cleanup_stale_tmp(
        root,
        older_than_days=0,
        max_files=1,
    )
    assert removed == [candidates[0]]
    assert candidates[2].is_dir()


def test_family_states_use_bounded_bucketed_current_state_files(tmp_path: Path):
    root = tmp_path / "orthogroup"
    for index in range(300):
        family_id = f"OG{index:07d}"
        GeneFamilyOutputStore(root, family_filter=family_id).mark_family_state(
            family_id,
            "complete",
            f"run-{index}",
        )

    state_root = root / ".gg_archives" / "family-state"
    assert (state_root / ".complete").is_file()
    assert 1 <= len(list(state_root.glob("*.json"))) <= 256
    assert not (root / ".gg_archives" / "family-state.jsonl").exists()
    assert GeneFamilyOutputStore(
        root,
        family_filter="OG0000299",
    ).family_state("OG0000299") == "complete"


def test_query_catalog_addition_can_rebucket_existing_logical_members(
    tmp_path: Path,
):
    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "q").write_text("short\n", encoding="utf-8")
    _write_family_outputs(root, "q_long", complete=True)
    GeneFamilyOutputStore(root, family_filter="q").mark_family_state(
        "q",
        "complete",
    )
    family_ids, family_from_name = family_context("query2family", query_dir=query_dir)
    archive_completed_outputs(
        root,
        "query2family",
        family_ids,
        family_from_name,
    )

    (query_dir / "q_long").write_text("long\n", encoding="utf-8")
    _write_family_outputs(root, "q_long", complete=True)
    GeneFamilyOutputStore(root, family_filter="q_long").mark_family_state(
        "q_long",
        "complete",
    )
    family_ids, family_from_name = family_context("query2family", query_dir=query_dir)
    archive_completed_outputs(
        root,
        "query2family",
        family_ids,
        family_from_name,
    )

    GeneFamilyOutputStore(root).verify()
    assert not GeneFamilyOutputStore(
        root,
        family_filter="q",
    ).logical_exists("stat_branch/q_long_stat.branch.tsv")
    assert GeneFamilyOutputStore(
        root,
        family_filter="q_long",
    ).logical_exists("stat_branch/q_long_stat.branch.tsv")


def test_interrupted_index_update_fails_closed_and_repair_recovers(
    tmp_path: Path,
    monkeypatch,
):
    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "A").write_text("geneA\n", encoding="utf-8")
    paths = _write_family_outputs(root, "A", complete=True)
    family_ids, family_from_name = family_context("query2family", query_dir=query_dir)

    original = GeneFamilyOutputStore._write_json_index_payload
    calls = 0

    def fail_during_denormalized_index_update(self, path, payload):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("simulated index update interruption")
        return original(self, path, payload)

    monkeypatch.setattr(
        GeneFamilyOutputStore,
        "_write_json_index_payload",
        fail_during_denormalized_index_update,
    )
    with pytest.raises(OSError, match="simulated index update interruption"):
        archive_completed_outputs(
            root,
            "query2family",
            family_ids,
            family_from_name,
        )

    marker = root / ".gg_archives" / "index-update.pending"
    assert marker.is_file()
    assert paths["mafft"].is_file()
    with pytest.raises(ArchiveStoreError, match="interrupted archive index update"):
        GeneFamilyOutputStore(root).logical_exists("mafft/A_cds.aln.fa.gz")

    monkeypatch.setattr(
        GeneFamilyOutputStore,
        "_write_json_index_payload",
        original,
    )
    repair_archive_index(root)
    assert not marker.exists()
    repaired = GeneFamilyOutputStore(root)
    repaired.verify()
    assert repaired.logical_exists("mafft/A_cds.aln.fa.gz")


def test_verify_compares_family_and_subdirectory_index_views(tmp_path: Path):
    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "A").write_text("geneA\n", encoding="utf-8")
    _write_family_outputs(root, "A", complete=True)
    family_ids, family_from_name = family_context("query2family", query_dir=query_dir)
    archive_completed_outputs(root, "query2family", family_ids, family_from_name)

    archive_root = root / ".gg_archives"
    catalog = json.loads(
        (archive_root / "index-catalog.json").read_text(encoding="utf-8")
    )
    subdir_index = (
        archive_root
        / "index-by-subdir"
        / catalog["subdir_indexes"]["stat_branch"]
    )
    payload = json.loads(subdir_index.read_text(encoding="utf-8"))
    payload["artifacts"].clear()
    subdir_index.write_text(json.dumps(payload) + "\n", encoding="utf-8")

    with pytest.raises(ArchiveStoreError, match="Subdirectory index count differs"):
        GeneFamilyOutputStore(root).verify()

    repair_archive_index(root)
    GeneFamilyOutputStore(root).verify()


def test_invalid_family_bucket_record_is_not_silently_ignored(tmp_path: Path):
    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "A").write_text("geneA\n", encoding="utf-8")
    _write_family_outputs(root, "A", complete=True)
    family_ids, family_from_name = family_context("query2family", query_dir=query_dir)
    archive_completed_outputs(root, "query2family", family_ids, family_from_name)

    bucket_path = next((root / ".gg_archives" / "index").glob("*.json"))
    payload = json.loads(bucket_path.read_text(encoding="utf-8"))
    logical_path = next(iter(payload["artifacts"]))
    payload["artifacts"][logical_path] = "invalid"
    bucket_path.write_text(json.dumps(payload) + "\n", encoding="utf-8")

    with pytest.raises(ArchiveStoreError, match="invalid artifact record"):
        GeneFamilyOutputStore(root).verify()


def test_long_lived_unfiltered_store_refreshes_family_state(tmp_path: Path):
    root = tmp_path / "orthogroup"
    observer = GeneFamilyOutputStore(root)
    assert observer.family_state("OG0000001") is None

    writer = GeneFamilyOutputStore(root, family_filter="OG0000001")
    writer.mark_family_state("OG0000001", "complete", "run-1")

    assert observer.family_state("OG0000001") == "complete"


def test_symlinked_archive_root_is_rejected(tmp_path: Path):
    root = tmp_path / "orthogroup"
    root.mkdir()
    external = tmp_path / "external-archive"
    external.mkdir()
    (root / ".gg_archives").symlink_to(external, target_is_directory=True)

    with pytest.raises(ArchiveStoreError, match="Symlinked GeneGalleon archive roots"):
        GeneFamilyOutputStore(root)


def test_repair_prefers_newer_compaction_when_member_generations_tie(
    tmp_path: Path,
):
    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "A").write_text("geneA\n", encoding="utf-8")
    _write_family_outputs(root, "A", complete=True)
    family_ids, family_from_name = family_context("query2family", query_dir=query_dir)
    archive_completed_outputs(root, "query2family", family_ids, family_from_name)

    (query_dir / "B").write_text("geneB\n", encoding="utf-8")
    _write_family_outputs(root, "B", complete=True)
    family_ids, family_from_name = family_context("query2family", query_dir=query_dir)
    archive_completed_outputs(root, "query2family", family_ids, family_from_name)
    archive_root = root / ".gg_archives"
    old_shards = {
        path.relative_to(archive_root): path.read_bytes()
        for path in archive_root.glob("*/*.zip")
    }

    compact_archives(root, "query2family")
    compacted = GeneFamilyOutputStore(root).artifact(
        "stat_branch",
        "A_stat.branch.tsv",
    )
    assert compacted is not None
    assert compacted.zip_path is not None
    assert compacted.zip_path.name.startswith("pack-")

    # Recreate the crash state after the new index commit but before obsolete
    # shards were removed.
    for relative_path, contents in old_shards.items():
        path = archive_root / relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(contents)

    repair_archive_index(root, remove_orphans=True)
    repaired = GeneFamilyOutputStore(root)
    repaired_artifact = repaired.artifact("stat_branch", "A_stat.branch.tsv")
    assert repaired_artifact is not None
    assert repaired_artifact.zip_path is not None
    assert repaired_artifact.zip_path.name.startswith("pack-")
    repaired.verify()


def test_archive_keeps_source_changed_with_same_size_and_mtime(
    tmp_path: Path,
    monkeypatch,
):
    import workflow.support.gene_family_output_store as output_store_module

    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "A").write_text("geneA\n", encoding="utf-8")
    paths = _write_family_outputs(root, "A", complete=True)
    target = paths["stat_branch"]
    original_contents = target.read_bytes()
    replacement = b"X" * len(original_contents)
    original_archive_chunk = output_store_module._archive_chunk
    changed = False

    def archive_then_change_source(*args, **kwargs):
        nonlocal changed
        result = original_archive_chunk(*args, **kwargs)
        archived_paths = args[3]
        if target in archived_paths and not changed:
            original_mtime_ns = target.stat().st_mtime_ns
            target.write_bytes(replacement)
            os.utime(target, ns=(original_mtime_ns, original_mtime_ns))
            changed = True
        return result

    monkeypatch.setattr(
        output_store_module,
        "_archive_chunk",
        archive_then_change_source,
    )
    family_ids, family_from_name = family_context("query2family", query_dir=query_dir)
    archive_completed_outputs(root, "query2family", family_ids, family_from_name)

    assert changed
    assert target.read_bytes() == replacement
    with GeneFamilyOutputStore(root).open_binary("stat_branch", target.name) as handle:
        assert handle.read() == replacement


def test_manual_live_replacement_and_addition_are_archived_as_new_versions(
    tmp_path: Path,
):
    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "A").write_text("geneA\n", encoding="utf-8")
    _write_family_outputs(root, "A", complete=True)
    family_ids, family_from_name = family_context("query2family", query_dir=query_dir)
    archive_completed_outputs(root, "query2family", family_ids, family_from_name)

    replacement = root / "stat_branch" / "A_stat.branch.tsv"
    replacement.parent.mkdir(parents=True, exist_ok=True)
    replacement.write_bytes(b"manual replacement\n")
    added = root / "custom_result" / "A_manual.tsv"
    added.parent.mkdir(parents=True, exist_ok=True)
    added.write_bytes(b"manual addition\n")

    archive_completed_outputs(root, "query2family", family_ids, family_from_name)

    assert not replacement.exists()
    assert not added.exists()
    store = GeneFamilyOutputStore(root)
    with store.open_binary("stat_branch", replacement.name) as handle:
        assert handle.read() == b"manual replacement\n"
    with store.open_binary("custom_result", added.name) as handle:
        assert handle.read() == b"manual addition\n"
    store.verify()


def test_cleanup_materialized_cli_does_not_require_root(tmp_path: Path):
    receipt = tmp_path / ".gg_materialized.jsonl"
    receipt.write_text("", encoding="utf-8")
    args = build_parser().parse_args(
        ["cleanup-materialized", "--receipt", str(receipt)]
    )

    assert run_cli(args) == 0
    assert not receipt.exists()


def test_materialize_rejects_symlinked_destination_subdirectory(tmp_path: Path):
    root = tmp_path / "orthogroup"
    family_id = "OG0000001"
    _write_family_outputs(root, family_id, complete=True)
    archive_completed_outputs(
        root,
        "orthogroup",
        [family_id],
        lambda name: name.split("_", 1)[0],
    )
    external = tmp_path / "external"
    external.mkdir()
    (root / "stat_branch").symlink_to(external, target_is_directory=True)

    with pytest.raises(ArchiveStoreError, match="Symlinked materialization destinations"):
        GeneFamilyOutputStore(root).materialize(
            "stat_branch",
            f"{family_id}_stat.branch.tsv",
        )

    assert not any(external.iterdir())


def test_malformed_materialization_receipt_fails_closed(tmp_path: Path):
    receipt = tmp_path / ".gg_materialized.jsonl"
    receipt.write_text("[]\n", encoding="utf-8")

    with pytest.raises(ArchiveStoreError, match="not a JSON object"):
        cleanup_materialization_receipt(receipt)

    assert receipt.is_file()


def test_raw_to_zip_conversion_archives_partial_family_without_marking_complete(
    tmp_path: Path,
):
    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "A").write_text("geneA\n", encoding="utf-8")
    paths = _write_family_outputs(root, "A", complete=False)
    shared = root / "parameters" / "common.tsv"
    shared.parent.mkdir(parents=True)
    shared.write_text("shared\n", encoding="utf-8")
    family_ids, family_from_name = family_context(
        "query2family",
        query_dir=query_dir,
    )

    result = convert_storage_to_zip(
        root,
        "query2family",
        family_ids,
        family_from_name,
    )

    assert result["storage"] == "zip"
    assert not paths["mafft"].exists()
    assert shared.is_file()
    store = GeneFamilyOutputStore(root, family_filter="A")
    assert store.logical_exists("mafft/A_cds.aln.fa.gz")
    assert store.family_state("A") is None
    assert not (root / ".gg_archives" / "storage-conversion.pending").exists()


def test_partial_family_can_be_materialized_updated_and_rearchived(tmp_path: Path):
    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "A").write_text("geneA\n", encoding="utf-8")
    paths = _write_family_outputs(root, "A", complete=False)
    original_mtime_ns = paths["mafft"].stat().st_mtime_ns
    family_ids, family_from_name = family_context(
        "query2family",
        query_dir=query_dir,
    )
    convert_storage_to_zip(
        root,
        "query2family",
        family_ids,
        family_from_name,
    )

    restored = GeneFamilyOutputStore(root, family_filter="A").materialize_family(
        "A",
        family_from_name,
    )
    assert restored == [paths["mafft"]]
    assert paths["mafft"].stat().st_mtime_ns == original_mtime_ns
    paths["mafft"].write_bytes(b"updated partial output\n")
    GeneFamilyOutputStore(root, family_filter="A").mark_family_state(
        "A",
        "failed",
        "run-1",
    )
    archived = archive_completed_outputs(
        root,
        "query2family",
        ["A"],
        family_from_name,
        include_incomplete=True,
    )

    assert archived
    assert not paths["mafft"].exists()
    with GeneFamilyOutputStore(root, family_filter="A").open_binary(
        "mafft",
        paths["mafft"].name,
    ) as handle:
        assert handle.read() == b"updated partial output\n"
    assert GeneFamilyOutputStore(root, family_filter="A").family_state("A") == "failed"


def test_zip_to_raw_conversion_applies_tombstones_and_live_overrides(tmp_path: Path):
    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "A").write_text("geneA\n", encoding="utf-8")
    paths = _write_family_outputs(root, "A", complete=True)
    family_ids, family_from_name = family_context(
        "query2family",
        query_dir=query_dir,
    )
    convert_storage_to_zip(
        root,
        "query2family",
        family_ids,
        family_from_name,
    )
    store = GeneFamilyOutputStore(root)
    store.delete("tree_plot/A_tree_plot.pdf")
    replacement = root / "stat_branch" / "A_stat.branch.tsv"
    replacement.parent.mkdir(parents=True, exist_ok=True)
    replacement.write_bytes(b"manual live override\n")

    result = convert_storage_to_raw(root, "query2family")

    assert result["storage"] == "raw"
    assert replacement.read_bytes() == b"manual live override\n"
    assert not paths["tree_plot"].exists()
    assert paths["mafft"].is_file()
    assert not list((root / ".gg_archives").glob("*/*.zip"))
    GeneFamilyOutputStore(root).verify()

    convert_storage_to_zip(
        root,
        "query2family",
        family_ids,
        family_from_name,
    )
    assert not paths["mafft"].exists()
    assert GeneFamilyOutputStore(root).logical_exists(
        "mafft/A_cds.aln.fa.gz"
    )
    assert not GeneFamilyOutputStore(root).logical_exists(
        "tree_plot/A_tree_plot.pdf"
    )


def test_pure_raw_conversion_removes_archive_control_metadata(tmp_path: Path):
    root = tmp_path / "orthogroup"
    family_id = "OG0000001"
    paths = _write_family_outputs(root, family_id, complete=False)
    convert_storage_to_zip(
        root,
        "orthogroup",
        [family_id],
        lambda name: name.split("_", 1)[0],
    )
    interrupted_partial = (
        root / "mafft" / ".OG0000001_cds.aln.fa.gz.materialize.dead"
    )
    interrupted_partial.parent.mkdir(parents=True, exist_ok=True)
    interrupted_partial.write_bytes(b"partial")

    result = convert_storage_to_raw(root, "orthogroup", pure_raw=True)

    assert result["pure_raw"] is True
    assert paths["mafft"].is_file()
    assert not interrupted_partial.exists()
    assert not (root / ".gg_archives").exists()
    assert not (root / ".gg_archives.pure-raw-retired").exists()


def test_conversion_summary_reports_mixed_store_and_unmatched_files(tmp_path: Path):
    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "A").write_text("geneA\n", encoding="utf-8")
    _write_family_outputs(root, "A", complete=False)
    family_ids, family_from_name = family_context(
        "query2family",
        query_dir=query_dir,
    )
    convert_storage_to_zip(
        root,
        "query2family",
        family_ids,
        family_from_name,
    )
    live_override = root / "mafft" / "A_cds.aln.fa.gz"
    live_override.parent.mkdir(parents=True, exist_ok=True)
    live_override.write_bytes(b"override\n")
    unmatched = root / "parameters" / "common.tsv"
    unmatched.parent.mkdir(parents=True, exist_ok=True)
    unmatched.write_bytes(b"shared\n")

    summary, unmatched_paths = storage_conversion_summary(
        root,
        set(family_ids),
        family_from_name,
    )

    assert summary["storage"] == "mixed"
    assert summary["owned_live_files"] == 1
    assert unmatched_paths == [unmatched]


def test_raw_conversion_summary_combines_size_and_symlink_inventory(tmp_path: Path):
    root = tmp_path / "orthogroup"
    owned = root / "mafft" / "OG0000001_alignment.tsv"
    owned.parent.mkdir(parents=True)
    owned.write_bytes(b"owned")
    unmatched = root / "parameters" / "common.tsv"
    unmatched.parent.mkdir(parents=True)
    unmatched.write_bytes(b"shared")
    linked = root / "mafft" / "OG0000002_alignment.tsv"
    linked.symlink_to(owned)

    summary, unmatched_paths = storage_conversion_summary(
        root,
        {"OG0000001", "OG0000002"},
        orthogroup_id_from_name,
    )

    assert summary["storage"] == "raw"
    assert summary["logical_files"] == 2
    assert summary["owned_live_files"] == 1
    assert summary["owned_live_bytes"] == len(b"owned")
    assert summary["unmatched_live_files"] == 1
    assert summary["unsupported_symlinks"] == 1
    assert unmatched_paths == [unmatched]


def test_convert_storage_cli_accepts_files_as_raw_alias(tmp_path: Path):
    root = tmp_path / "orthogroup"
    family_id = "OG0000001"
    _write_family_outputs(root, family_id, complete=False)
    convert_storage_to_zip(
        root,
        "orthogroup",
        [family_id],
        lambda name: name.split("_", 1)[0],
    )
    args = build_parser().parse_args(
        [
            "convert-storage",
            "--root",
            str(root),
            "--mode",
            "orthogroup",
            "--to",
            "files",
        ]
    )

    assert run_cli(args) == 0
    assert (root / "mafft" / f"{family_id}_cds.aln.fa.gz").is_file()
    assert not list((root / ".gg_archives").glob("*/*.zip"))


def test_raw_to_zip_conversion_resumes_after_source_cleanup_interruption(
    tmp_path: Path,
    monkeypatch,
):
    import workflow.support.gene_family_output_store as output_store_module

    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "A").write_text("geneA\n", encoding="utf-8")
    paths = _write_family_outputs(root, "A", complete=False)
    family_ids, family_from_name = family_context(
        "query2family",
        query_dir=query_dir,
    )
    original_remove = output_store_module._remove_archived_sources

    def interrupt_source_cleanup(root, signatures):
        raise OSError("simulated conversion interruption")

    monkeypatch.setattr(
        output_store_module,
        "_remove_archived_sources",
        interrupt_source_cleanup,
    )
    with pytest.raises(OSError, match="simulated conversion interruption"):
        convert_storage_to_zip(
            root,
            "query2family",
            family_ids,
            family_from_name,
        )

    marker = root / ".gg_archives" / "storage-conversion.pending"
    assert marker.is_file()
    assert paths["mafft"].is_file()
    monkeypatch.setattr(
        output_store_module,
        "_remove_archived_sources",
        original_remove,
    )

    result = convert_storage_to_zip(
        root,
        "query2family",
        family_ids,
        family_from_name,
    )

    assert result["storage"] == "zip"
    assert not marker.exists()
    assert not paths["mafft"].exists()
    GeneFamilyOutputStore(root).verify()


def test_removed_query_family_can_be_added_to_conversion_catalog(tmp_path: Path):
    root = tmp_path / "query2family"
    paths = _write_family_outputs(root, "removed_query", complete=False)
    family_id_file = tmp_path / "legacy-family-ids.txt"
    family_id_file.write_text("removed_query\n", encoding="utf-8")
    family_ids, family_from_name = family_context_with_supplement(
        "query2family",
        family_id_file=family_id_file,
    )

    convert_storage_to_zip(
        root,
        "query2family",
        family_ids,
        family_from_name,
    )

    assert not paths["mafft"].exists()
    assert GeneFamilyOutputStore(
        root,
        family_filter="removed_query",
    ).logical_exists("mafft/removed_query_cds.aln.fa.gz")


def test_strict_unmatched_conversion_fails_before_creating_marker(tmp_path: Path):
    root = tmp_path / "orthogroup"
    root.mkdir()
    shared = root / "parameters" / "common.tsv"
    shared.parent.mkdir()
    shared.write_text("shared\n", encoding="utf-8")

    with pytest.raises(ArchiveStoreError, match="Unmatched live output files"):
        convert_storage_to_zip(
            root,
            "orthogroup",
            ["OG0000001"],
            lambda name: name.split("_", 1)[0] if "_" in name else None,
            strict_unmatched=True,
        )

    assert shared.is_file()
    assert not (root / ".gg_archives" / "storage-conversion.pending").exists()


def test_archive_family_cli_archives_failed_partial_outputs(tmp_path: Path):
    root = tmp_path / "query2family"
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "A").write_text("geneA\n", encoding="utf-8")
    paths = _write_family_outputs(root, "A", complete=False)
    GeneFamilyOutputStore(root, family_filter="A").mark_family_state(
        "A",
        "failed",
        "run-1",
    )
    args = build_parser().parse_args(
        [
            "archive-family",
            "--root",
            str(root),
            "--mode",
            "query2family",
            "--query-dir",
            str(query_dir),
            "--family-id",
            "A",
        ]
    )

    assert run_cli(args) == 0
    assert not paths["mafft"].exists()
    assert GeneFamilyOutputStore(root, family_filter="A").logical_exists(
        "mafft/A_cds.aln.fa.gz"
    )
    assert GeneFamilyOutputStore(root, family_filter="A").family_state("A") == "failed"


def test_raw_to_zip_conversion_rejects_symlinked_family_artifact(tmp_path: Path):
    root = tmp_path / "orthogroup"
    external = tmp_path / "external.tsv"
    external.write_text("external\n", encoding="utf-8")
    linked = root / "mafft" / "OG0000001_cds.aln.fa.gz"
    linked.parent.mkdir(parents=True)
    linked.symlink_to(external)

    with pytest.raises(ArchiveStoreError, match="Symlinked output paths"):
        convert_storage_to_zip(
            root,
            "orthogroup",
            ["OG0000001"],
            lambda name: name.split("_", 1)[0],
        )

    assert linked.is_symlink()
    assert external.read_text(encoding="utf-8") == "external\n"


def test_archive_rejects_mode_different_from_existing_shards(tmp_path: Path):
    root = tmp_path / "query2family"
    family_id = "OG0000001"
    _write_family_outputs(root, family_id, complete=False)
    matcher = lambda name: family_id if name.startswith(family_id + "_") else None
    archive_completed_outputs(
        root,
        "query2family",
        [family_id],
        matcher,
        include_incomplete=True,
    )
    replacement = root / "rpsblast" / f"{family_id}_rpsblast.tsv"
    replacement.parent.mkdir(parents=True)
    replacement.write_text("replacement\n", encoding="utf-8")

    for maintenance in (
        lambda: archive_completed_outputs(
            root,
            "orthogroup",
            [family_id],
            matcher,
            include_incomplete=True,
        ),
        lambda: compact_archives(root, "orthogroup"),
        lambda: purge_archives(root, "orthogroup"),
    ):
        with pytest.raises(ArchiveStoreError, match="different gene-family mode"):
            maintenance()

    assert replacement.is_file()
    GeneFamilyOutputStore(root).verify()


def test_archive_mode_mismatch_is_rejected_without_archive_candidates(
    tmp_path: Path,
):
    root = tmp_path / "query2family"
    family_id = "OG0000001"
    _write_family_outputs(root, family_id, complete=False)
    matcher = lambda name: family_id if name.startswith(family_id + "_") else None
    archive_completed_outputs(
        root,
        "query2family",
        [family_id],
        matcher,
        include_incomplete=True,
    )

    with pytest.raises(ArchiveStoreError, match="different gene-family mode"):
        archive_completed_outputs(
            root,
            "orthogroup",
            [family_id],
            matcher,
            min_files=100,
            include_incomplete=True,
        )


def test_verify_rejects_referenced_shards_with_mixed_modes(
    tmp_path: Path,
    monkeypatch,
):
    import workflow.support.gene_family_output_store as output_store_module

    root = tmp_path / "query2family"
    family_id = "OG0000001"
    _write_family_outputs(root, family_id, complete=False)
    matcher = lambda name: family_id if name.startswith(family_id + "_") else None
    archive_completed_outputs(
        root,
        "query2family",
        [family_id],
        matcher,
        include_incomplete=True,
    )
    second = root / "rpsblast" / f"{family_id}_rpsblast.tsv"
    second.parent.mkdir(parents=True)
    second.write_text("second mode\n", encoding="utf-8")
    monkeypatch.setattr(output_store_module, "_assert_archive_mode", lambda *_: None)
    archive_completed_outputs(
        root,
        "orthogroup",
        [family_id],
        matcher,
        include_incomplete=True,
    )

    with pytest.raises(ArchiveStoreError, match="mixed gene-family modes"):
        repair_archive_index(root)
    with pytest.raises(ArchiveStoreError, match="mixed gene-family modes"):
        GeneFamilyOutputStore(root).verify()


def test_zip_to_raw_requirements_include_blocks_and_missing_subdirectories(
    tmp_path: Path,
):
    root = tmp_path / "orthogroup"
    (root / "existing").mkdir(parents=True)
    artifacts = [
        Artifact(
            logical_path="existing/OG0000001_a.tsv",
            subdir="existing",
            name="OG0000001_a.tsv",
            generation=1,
            size=1,
        ),
        Artifact(
            logical_path="missing/OG0000001_b.tsv",
            subdir="missing",
            name="OG0000001_b.tsv",
            generation=1,
            size=4097,
        ),
    ]

    required_bytes, required_inodes = _zip_to_raw_requirements(
        root,
        artifacts,
        block_size=4096,
    )

    assert required_bytes == 4 * 4096
    assert required_inodes == 4


def test_storage_status_cli_can_inspect_physical_store_without_catalog(
    tmp_path: Path,
    capsys,
):
    root = tmp_path / "orthogroup"
    root.mkdir()
    args = build_parser().parse_args(
        [
            "storage-status",
            "--root",
            str(root),
            "--mode",
            "orthogroup",
        ]
    )

    assert run_cli(args) == 0
    assert "storage\traw" in capsys.readouterr().out


def test_convert_storage_dry_run_does_not_create_archive_metadata(
    tmp_path: Path,
    capsys,
):
    root = tmp_path / "orthogroup"
    genecount = tmp_path / "Orthogroups.GeneCount.selected.tsv"
    genecount.write_text(
        "Orthogroup\tTotal\nOG0000001\t1\n",
        encoding="utf-8",
    )
    paths = _write_family_outputs(root, "OG0000001", complete=False)
    args = build_parser().parse_args(
        [
            "convert-storage",
            "--root",
            str(root),
            "--mode",
            "orthogroup",
            "--genecount",
            str(genecount),
            "--to",
            "zip",
            "--dry-run",
        ]
    )

    assert run_cli(args) == 0
    output = capsys.readouterr().out
    assert "requested_target\tzip" in output
    assert paths["mafft"].is_file()
    assert not (root / ".gg_archives").exists()

from __future__ import annotations

import importlib.util
import json
import os
import re
import stat
import struct
import subprocess
import sys
import time
import zipfile
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO_ROOT / "workflow" / "support" / "species_tree_output_store.py"
WRAPPER_PATH = REPO_ROOT / "workflow" / "gg_species_tree_archive.sh"
CORE_PATH = REPO_ROOT / "workflow" / "core" / "gg_genome_evolution_core.sh"
PROGRESS_CORE_PATH = REPO_ROOT / "workflow" / "core" / "gg_progress_summary_core.sh"

SPEC = importlib.util.spec_from_file_location("species_tree_output_store", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
STORE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(STORE)


def _write_marked_archive(
    root: Path,
    name: str,
    members: list[tuple[zipfile.ZipInfo | str, bytes]],
) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    archive_path = root / f"{name}.zip"
    file_count = sum(
        not item.is_dir() if isinstance(item, zipfile.ZipInfo) else not item.endswith("/")
        for item, _ in members
    )
    directory_count = len(members) - file_count
    with zipfile.ZipFile(archive_path, "w") as archive:
        for member, data in members:
            archive.writestr(member, data)
        archive.comment = STORE._archive_comment(name, file_count, directory_count)
    return archive_path


def test_pack_and_materialize_visible_named_archive(tmp_path: Path):
    root = tmp_path / "species_tree"
    raw = root / "single_copy_iqtree_dna"
    raw.mkdir(parents=True)
    (raw / "BUSCO1.dna.nwk").write_text("(A,B);\n", encoding="utf-8")
    (raw / "nested").mkdir()
    (raw / "nested" / "BUSCO2.dna.nwk").write_text("(A,C);\n", encoding="utf-8")

    result = STORE.pack_directory(root, "single_copy_iqtree_dna")

    archive_path = root / "single_copy_iqtree_dna.zip"
    assert result == {
        "directory": "single_copy_iqtree_dna",
        "state": "archived",
        "files": 2,
    }
    assert archive_path.is_file()
    assert not raw.exists()
    with zipfile.ZipFile(archive_path) as archive:
        assert sorted(archive.namelist()) == [
            "single_copy_iqtree_dna/",
            "single_copy_iqtree_dna/BUSCO1.dna.nwk",
            "single_copy_iqtree_dna/nested/",
            "single_copy_iqtree_dna/nested/BUSCO2.dna.nwk",
        ]
        assert archive.comment.startswith(STORE.COMMENT_PREFIX)

    restored = STORE.materialize_directory(root, "single_copy_iqtree_dna")

    assert restored["state"] == "raw"
    assert not archive_path.exists()
    assert (raw / "BUSCO1.dna.nwk").read_text(encoding="utf-8") == "(A,B);\n"
    assert (raw / "nested" / "BUSCO2.dna.nwk").read_text(encoding="utf-8") == "(A,C);\n"


def test_materialize_quota_preflight_preserves_archive(tmp_path: Path):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_dna"
    raw = root / name
    raw.mkdir(parents=True)
    (raw / "BUSCO1.dna.nwk").write_text("(A,B);\n", encoding="utf-8")
    STORE.pack_directory(root, name)

    with pytest.raises(STORE.SpeciesTreeArchiveError, match="ZIP-to-raw"):
        STORE.materialize_directory(root, name, available_bytes=0)

    assert not raw.exists()
    assert (root / f"{name}.zip").is_file()
    assert STORE.verify_archive(root, name)["state"] == "archived"


def test_materialize_inode_preflight_treats_zero_as_exhausted(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_dna"
    raw = root / name
    raw.mkdir(parents=True)
    (raw / "BUSCO1.dna.nwk").write_text("(A,B);\n", encoding="utf-8")
    STORE.pack_directory(root, name)
    real_stats = os.statvfs(root)
    exhausted = SimpleNamespace(
        **{
            field: (0 if field == "f_favail" else getattr(real_stats, field))
            for field in dir(real_stats)
            if field.startswith("f_")
        }
    )
    monkeypatch.setattr(STORE.os, "statvfs", lambda _path: exhausted)

    with pytest.raises(STORE.SpeciesTreeArchiveError, match="available=0"):
        STORE.materialize_directory(root, name)

    assert not raw.exists()
    assert (root / f"{name}.zip").is_file()


def test_materialize_accepts_filesystem_without_inode_accounting(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_dna"
    raw = root / name
    raw.mkdir(parents=True)
    expected = raw / "BUSCO1.dna.nwk"
    expected.write_text("(A,B);\n", encoding="utf-8")
    STORE.pack_directory(root, name)
    real_stats = os.statvfs(root)
    unavailable = SimpleNamespace(
        **{
            field: (
                0
                if field in {"f_files", "f_favail"}
                else getattr(real_stats, field)
            )
            for field in dir(real_stats)
            if field.startswith("f_")
        }
    )
    monkeypatch.setattr(STORE.os, "statvfs", lambda _path: unavailable)

    STORE.materialize_directory(root, name)

    assert expected.read_text(encoding="utf-8") == "(A,B);\n"


def test_manual_delete_and_add_are_preserved_after_unpack_and_repack(tmp_path: Path):
    root = tmp_path / "species_tree"
    raw = root / "single_copy_trimal"
    raw.mkdir(parents=True)
    removed = raw / "BUSCO1.trimal.fa.gz"
    removed.write_bytes(b"old")
    STORE.pack_directory(root, "single_copy_trimal")
    STORE.materialize_directory(root, "single_copy_trimal")

    removed.unlink()
    added = raw / "BUSCO2.trimal.fa.gz"
    added.write_bytes(b"new")
    STORE.pack_directory(root, "single_copy_trimal")

    with zipfile.ZipFile(root / "single_copy_trimal.zip") as archive:
        files = [info.filename for info in archive.infolist() if not info.is_dir()]
        assert files == ["single_copy_trimal/BUSCO2.trimal.fa.gz"]
        assert archive.read(files[0]) == b"new"


def test_empty_directories_and_directory_modes_survive_round_trip(tmp_path: Path):
    root = tmp_path / "species_tree"
    raw = root / "single_copy_mafft"
    empty = raw / "nested" / "empty"
    empty.mkdir(parents=True)
    raw.chmod(0o750)
    empty.chmod(0o700)

    STORE.pack_directory(root, "single_copy_mafft")
    STORE.materialize_directory(root, "single_copy_mafft")

    assert empty.is_dir()
    assert raw.stat().st_mode & 0o777 == 0o750
    assert empty.stat().st_mode & 0o777 == 0o700


def test_pack_refuses_ambiguous_raw_and_zip_forms(tmp_path: Path):
    root = tmp_path / "species_tree"
    raw = root / "single_copy_mafft"
    raw.mkdir(parents=True)
    (raw / "BUSCO1.aln.fa.gz").write_bytes(b"one")
    STORE.pack_directory(root, "single_copy_mafft")
    raw.mkdir()
    (raw / "BUSCO2.aln.fa.gz").write_bytes(b"two")

    with pytest.raises(STORE.SpeciesTreeArchiveError, match="Both raw and ZIP forms exist"):
        STORE.pack_directory(root, "single_copy_mafft")

    STORE.materialize_directory(root, "single_copy_mafft")
    assert (raw / "BUSCO1.aln.fa.gz").read_bytes() == b"one"
    assert (raw / "BUSCO2.aln.fa.gz").read_bytes() == b"two"


def test_materialize_keeps_live_file_when_recovering_mixed_state(tmp_path: Path):
    root = tmp_path / "species_tree"
    raw = root / "single_copy_cds_fasta"
    raw.mkdir(parents=True)
    path = raw / "BUSCO1.fa.gz"
    path.write_bytes(b"archived")
    STORE.pack_directory(root, "single_copy_cds_fasta")
    raw.mkdir()
    path.write_bytes(b"newer-live")

    STORE.materialize_directory(root, "single_copy_cds_fasta")

    assert path.read_bytes() == b"newer-live"
    assert not (root / "single_copy_cds_fasta.zip").exists()


def test_verify_rejects_unmarked_zip(tmp_path: Path):
    root = tmp_path / "species_tree"
    root.mkdir()
    archive_path = root / "single_copy_iqtree_pep.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("single_copy_iqtree_pep/BUSCO1.pep.nwk", "(A,B);\n")

    with pytest.raises(STORE.SpeciesTreeArchiveError, match="without a GeneGalleon"):
        STORE.verify_archive(root, "single_copy_iqtree_pep")


def test_verify_rejects_non_object_archive_metadata(tmp_path: Path):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_pep"
    root.mkdir()
    archive_path = root / f"{name}.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.comment = STORE.COMMENT_PREFIX + b"[]"

    with pytest.raises(STORE.SpeciesTreeArchiveError, match="JSON object"):
        STORE.verify_archive(root, name)


@pytest.mark.parametrize(
    "member_name",
    [
        "single_copy_iqtree_pep/../../escaped.txt",
        "/single_copy_iqtree_pep/absolute.txt",
        "single_copy_iqtree_pep//noncanonical.txt",
        "single_copy_iqtree_pep/..\\escaped.txt",
        "different_directory/file.txt",
    ],
)
def test_materialize_rejects_unsafe_member_paths(tmp_path: Path, member_name: str):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_pep"
    archive_path = _write_marked_archive(
        root,
        name,
        [(member_name, b"unsafe")],
    )

    with pytest.raises(STORE.SpeciesTreeArchiveError, match="Unsafe ZIP member"):
        STORE.materialize_directory(root, name)

    assert archive_path.is_file()
    assert not (tmp_path / "escaped.txt").exists()
    assert not (root / name).exists()


def test_materialize_rejects_symlink_member(tmp_path: Path):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_pep"
    member = zipfile.ZipInfo(f"{name}/link")
    member.create_system = 3
    member.external_attr = (stat.S_IFLNK | 0o777) << 16
    archive_path = _write_marked_archive(root, name, [(member, b"../../outside")])

    with pytest.raises(STORE.SpeciesTreeArchiveError, match="Symlinked ZIP members"):
        STORE.materialize_directory(root, name)

    assert archive_path.is_file()
    assert not (root / name).exists()


def test_materialize_rejects_member_nested_below_file(tmp_path: Path):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_pep"
    archive_path = _write_marked_archive(
        root,
        name,
        [
            (f"{name}/parent", b"file"),
            (f"{name}/parent/child", b"nested"),
        ],
    )

    with pytest.raises(STORE.SpeciesTreeArchiveError, match="nested below a file"):
        STORE.materialize_directory(root, name)

    assert archive_path.is_file()
    assert not (root / name).exists()


def test_verify_rejects_duplicate_member_names(tmp_path: Path):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_pep"
    archive_path = root / f"{name}.zip"
    root.mkdir()
    with pytest.warns(UserWarning, match="Duplicate name"):
        with zipfile.ZipFile(archive_path, "w") as archive:
            archive.writestr(f"{name}/duplicate", b"first")
            archive.writestr(f"{name}/duplicate", b"second")
            archive.comment = STORE._archive_comment(name, 2, 0)

    with pytest.raises(STORE.SpeciesTreeArchiveError, match="Unsafe ZIP member"):
        STORE.verify_archive(root, name)


def test_verify_rejects_directory_count_mismatch(tmp_path: Path):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_pep"
    archive_path = _write_marked_archive(
        root,
        name,
        [(f"{name}/", b""), (f"{name}/empty/", b"")],
    )
    with zipfile.ZipFile(archive_path, "a") as archive:
        archive.comment = STORE._archive_comment(name, 0, 1)

    with pytest.raises(STORE.SpeciesTreeArchiveError, match="directory count"):
        STORE.verify_archive(root, name)


def test_verify_rejects_archive_directory_path(tmp_path: Path):
    root = tmp_path / "species_tree"
    (root / "single_copy_iqtree_pep.zip").mkdir(parents=True)

    with pytest.raises(STORE.SpeciesTreeArchiveError, match="Archive path is not a file"):
        STORE.verify_archive(root, "single_copy_iqtree_pep")


def test_corrupt_archive_is_preserved_when_verify_and_materialize_fail(tmp_path: Path):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_pep"
    raw = root / name
    raw.mkdir(parents=True)
    (raw / "BUSCO1.pep.nwk").write_bytes(b"0123456789" * 100)
    STORE.pack_directory(root, name, compression="store")
    archive_path = root / f"{name}.zip"
    with zipfile.ZipFile(archive_path) as archive:
        info = archive.getinfo(f"{name}/BUSCO1.pep.nwk")
    with archive_path.open("r+b") as handle:
        handle.seek(info.header_offset + 26)
        name_length, extra_length = struct.unpack("<HH", handle.read(4))
        data_offset = info.header_offset + 30 + name_length + extra_length
        handle.seek(data_offset + 10)
        original = handle.read(1)
        handle.seek(data_offset + 10)
        handle.write(bytes([original[0] ^ 0xFF]))

    status = STORE.status(root, [name])
    assert status[0]["state"] == "archived"
    assert STORE.verify_archive(root, name, check_crc=False)["verification_mode"] == "quick"
    with pytest.raises(STORE.SpeciesTreeArchiveError, match="CRC|Bad CRC"):
        STORE.verify_archive(root, name)
    with pytest.raises(STORE.SpeciesTreeArchiveError, match="extract|CRC|Bad CRC"):
        STORE.materialize_directory(root, name)

    assert archive_path.is_file()
    assert not raw.exists()
    assert not list(root.glob(f".{name}.materialize.partial.*"))


def test_convert_storage_dry_run_reports_space_without_writing(tmp_path: Path):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_dna"
    raw = root / name
    raw.mkdir(parents=True)
    (raw / "BUSCO1.dna.nwk").write_text("(A,B);\n", encoding="utf-8")
    args = STORE.build_parser().parse_args(
        [
            "convert-storage",
            "--root",
            str(root),
            "--directory",
            name,
            "--to",
            "zip",
            "--dry-run",
            "--available-bytes",
            "0",
        ]
    )

    assert STORE.run_cli(args) == 0
    assert raw.is_dir()
    assert not (root / f"{name}.zip").exists()


def test_pack_rejects_source_symlinks_without_deleting_raw(tmp_path: Path):
    root = tmp_path / "species_tree"
    name = "single_copy_trimal"
    raw = root / name
    raw.mkdir(parents=True)
    target = tmp_path / "outside.txt"
    target.write_text("outside", encoding="utf-8")
    (raw / "link").symlink_to(target)

    with pytest.raises(STORE.SpeciesTreeArchiveError, match="Symlinked source file"):
        STORE.pack_directory(root, name)

    assert raw.is_dir()
    assert target.read_text(encoding="utf-8") == "outside"
    assert not (root / f"{name}.zip").exists()


def test_materialize_rejects_live_symlink_in_mixed_state(tmp_path: Path):
    root = tmp_path / "species_tree"
    name = "single_copy_trimal"
    raw = root / name
    raw.mkdir(parents=True)
    (raw / "BUSCO1.fa.gz").write_bytes(b"archived")
    STORE.pack_directory(root, name)
    raw.mkdir()
    target = tmp_path / "outside.txt"
    target.write_text("outside", encoding="utf-8")
    (raw / "BUSCO1.fa.gz").symlink_to(target)

    with pytest.raises(STORE.SpeciesTreeArchiveError, match="Symlinked source file"):
        STORE.materialize_directory(root, name)

    assert (root / f"{name}.zip").is_file()
    assert (raw / "BUSCO1.fa.gz").is_symlink()
    assert target.read_text(encoding="utf-8") == "outside"


def test_next_operation_removes_stale_same_directory_partials(tmp_path: Path):
    root = tmp_path / "species_tree"
    raw = root / "single_copy_iqtree_dna"
    raw.mkdir(parents=True)
    (raw / "BUSCO1.dna.nwk").write_text("(A,B);\n", encoding="utf-8")
    stale_zip = root / ".single_copy_iqtree_dna.zip.partial.100.dead"
    stale_zip.write_bytes(b"partial")
    stale_extract = root / ".single_copy_iqtree_dna.materialize.partial.100.dead"
    stale_extract.mkdir()
    (stale_extract / "leftover").write_text("partial", encoding="utf-8")

    STORE.pack_directory(root, "single_copy_iqtree_dna")

    assert not stale_zip.exists()
    assert not stale_extract.exists()
    assert not (root / ".gg_species_tree_archive.lock").exists()


def test_pack_keeps_raw_when_atomic_replace_fails(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_dna"
    raw = root / name
    raw.mkdir(parents=True)
    source = raw / "BUSCO1.dna.nwk"
    source.write_text("(A,B);\n", encoding="utf-8")
    archive_path = root / f"{name}.zip"
    real_replace = STORE.os.replace

    def fail_archive_replace(source_path, destination_path):
        if Path(destination_path) == archive_path:
            raise OSError("injected replace failure")
        return real_replace(source_path, destination_path)

    monkeypatch.setattr(STORE.os, "replace", fail_archive_replace)

    with pytest.raises(OSError, match="injected replace failure"):
        STORE.pack_directory(root, name)

    assert source.read_text(encoding="utf-8") == "(A,B);\n"
    assert not archive_path.exists()
    assert not list(root.glob(f".{name}.zip.partial.*"))


def test_pack_keeps_raw_and_removes_partial_when_archive_write_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_dna"
    raw = root / name
    raw.mkdir(parents=True)
    source = raw / "BUSCO1.dna.nwk"
    content = b"(A,B);\n" * 1000
    source.write_bytes(content)
    real_copyfileobj = STORE.shutil.copyfileobj

    def fail_after_first_chunk(source_handle, destination_handle, length=0):
        destination_handle.write(source_handle.read(128))
        raise OSError("injected archive write failure")

    monkeypatch.setattr(STORE.shutil, "copyfileobj", fail_after_first_chunk)

    with pytest.raises(OSError, match="injected archive write failure"):
        STORE.pack_directory(root, name)

    monkeypatch.setattr(STORE.shutil, "copyfileobj", real_copyfileobj)
    assert source.read_bytes() == content
    assert not (root / f"{name}.zip").exists()
    assert not list(root.glob(f".{name}.zip.partial.*"))


def test_materialize_keeps_zip_and_removes_partial_when_extraction_write_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_dna"
    raw = root / name
    raw.mkdir(parents=True)
    content = b"(A,B);\n" * 1000
    (raw / "BUSCO1.dna.nwk").write_bytes(content)
    STORE.pack_directory(root, name)
    archive_path = root / f"{name}.zip"

    def fail_after_first_chunk(source_handle, destination_handle, length=0):
        destination_handle.write(source_handle.read(128))
        raise OSError("injected extraction write failure")

    monkeypatch.setattr(STORE.shutil, "copyfileobj", fail_after_first_chunk)

    with pytest.raises(STORE.SpeciesTreeArchiveError, match="extraction write failure"):
        STORE.materialize_directory(root, name)

    assert archive_path.is_file()
    assert not raw.exists()
    assert not list(root.glob(f".{name}.materialize.partial.*"))


def test_failed_raw_removal_leaves_recoverable_mixed_state(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_dna"
    raw = root / name
    raw.mkdir(parents=True)
    source = raw / "BUSCO1.dna.nwk"
    source.write_text("(A,B);\n", encoding="utf-8")
    real_rmtree = STORE.shutil.rmtree

    def fail_raw_removal(path, *args, **kwargs):
        if Path(path) == raw:
            raise OSError("injected removal failure")
        return real_rmtree(path, *args, **kwargs)

    monkeypatch.setattr(STORE.shutil, "rmtree", fail_raw_removal)
    with pytest.raises(OSError, match="injected removal failure"):
        STORE.pack_directory(root, name)

    archive_path = root / f"{name}.zip"
    assert archive_path.is_file()
    assert source.is_file()

    monkeypatch.setattr(STORE.shutil, "rmtree", real_rmtree)
    STORE.materialize_directory(root, name)
    assert source.read_text(encoding="utf-8") == "(A,B);\n"
    assert not archive_path.exists()


def test_failed_archive_unlink_leaves_recoverable_mixed_state(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_pep"
    raw = root / name
    raw.mkdir(parents=True)
    source = raw / "BUSCO1.pep.nwk"
    source.write_text("(A,B);\n", encoding="utf-8")
    STORE.pack_directory(root, name)
    archive_path = root / f"{name}.zip"
    real_unlink = Path.unlink

    def fail_archive_unlink(path: Path, *args, **kwargs):
        if path == archive_path:
            raise OSError("injected unlink failure")
        return real_unlink(path, *args, **kwargs)

    monkeypatch.setattr(Path, "unlink", fail_archive_unlink)
    with pytest.raises(STORE.SpeciesTreeArchiveError, match="injected unlink failure"):
        STORE.materialize_directory(root, name)

    assert archive_path.is_file()
    assert source.read_text(encoding="utf-8") == "(A,B);\n"

    monkeypatch.setattr(Path, "unlink", real_unlink)
    STORE.materialize_directory(root, name)
    assert source.is_file()
    assert not archive_path.exists()


def test_source_change_during_pack_never_replaces_raw(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_dna"
    raw = root / name
    raw.mkdir(parents=True)
    source = raw / "BUSCO1.dna.nwk"
    source.write_text("before\n", encoding="utf-8")
    real_validate = STORE._validated_members
    mutated = False

    def mutate_after_zip_validation(*args, **kwargs):
        nonlocal mutated
        result = real_validate(*args, **kwargs)
        archive_path = Path(args[1])
        if ".zip.partial." in archive_path.name and not mutated:
            source.write_text("after\n", encoding="utf-8")
            mutated = True
        return result

    monkeypatch.setattr(STORE, "_validated_members", mutate_after_zip_validation)

    with pytest.raises(STORE.SpeciesTreeArchiveError, match="Source changed"):
        STORE.pack_directory(root, name)

    assert source.read_text(encoding="utf-8") == "after\n"
    assert not (root / f"{name}.zip").exists()
    assert not list(root.glob(f".{name}.zip.partial.*"))


def test_same_size_source_change_with_restored_mtime_is_detected_during_pack(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_dna"
    raw = root / name
    raw.mkdir(parents=True)
    source = raw / "BUSCO1.dna.nwk"
    source.write_bytes(b"before\n")
    real_validate = STORE._validated_members
    mutated = False

    def mutate_after_zip_validation(*args, **kwargs):
        nonlocal mutated
        result = real_validate(*args, **kwargs)
        archive_path = Path(args[1])
        if ".zip.partial." in archive_path.name and not mutated:
            original_mtime_ns = source.stat().st_mtime_ns
            source.write_bytes(b"AFTER!\n")
            os.utime(source, ns=(original_mtime_ns, original_mtime_ns))
            mutated = True
        return result

    monkeypatch.setattr(STORE, "_validated_members", mutate_after_zip_validation)

    with pytest.raises(STORE.SpeciesTreeArchiveError, match="Source changed"):
        STORE.pack_directory(root, name)

    assert source.read_bytes() == b"AFTER!\n"
    assert not (root / f"{name}.zip").exists()
    assert not list(root.glob(f".{name}.zip.partial.*"))


@pytest.mark.parametrize(
    ("compression", "expected_text_type", "expected_gzip_type"),
    [
        ("adaptive", zipfile.ZIP_DEFLATED, zipfile.ZIP_STORED),
        ("deflate", zipfile.ZIP_DEFLATED, zipfile.ZIP_DEFLATED),
        ("store", zipfile.ZIP_STORED, zipfile.ZIP_STORED),
    ],
)
def test_compression_modes(
    tmp_path: Path,
    compression: str,
    expected_text_type: int,
    expected_gzip_type: int,
):
    root = tmp_path / compression
    name = "single_copy_cds_fasta"
    raw = root / name
    raw.mkdir(parents=True)
    (raw / "marker.nwk").write_text("(A,B);\n" * 100, encoding="utf-8")
    (raw / "marker.fa.gz").write_bytes(os.urandom(1024))

    STORE.pack_directory(root, name, compression=compression)

    with zipfile.ZipFile(root / f"{name}.zip") as archive:
        assert archive.getinfo(f"{name}/marker.nwk").compress_type == expected_text_type
        assert archive.getinfo(f"{name}/marker.fa.gz").compress_type == expected_gzip_type


def test_unicode_binary_and_file_mode_round_trip(tmp_path: Path):
    root = tmp_path / "species_tree"
    name = "single_copy_cds_fasta"
    raw = root / name
    raw.mkdir(parents=True)
    path = raw / "遺伝子 α.bin"
    content = bytes(range(256)) * 4
    path.write_bytes(content)
    path.chmod(0o640)

    STORE.pack_directory(root, name)
    STORE.materialize_directory(root, name)

    assert path.read_bytes() == content
    assert path.stat().st_mode & 0o777 == 0o640


def test_pre_1980_timestamp_is_clamped_to_standard_zip_range(tmp_path: Path):
    root = tmp_path / "species_tree"
    name = "single_copy_cds_fasta"
    raw = root / name
    raw.mkdir(parents=True)
    path = raw / "old.pep.fa.gz"
    path.write_bytes(b"old")
    os.utime(path, (0, 0))

    STORE.pack_directory(root, name)
    STORE.materialize_directory(root, name)

    assert path.read_bytes() == b"old"
    assert time.localtime(path.stat().st_mtime).tm_year == 1980


def test_concurrent_pack_and_materialize_are_serialized(tmp_path: Path):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_dna"
    raw = root / name
    raw.mkdir(parents=True)
    for index in range(200):
        (raw / f"BUSCO{index:04d}.dna.nwk").write_text(
            f"(A:{index},B:{index});\n", encoding="utf-8"
        )
    base_command = [
        sys.executable,
        str(MODULE_PATH),
        "pack",
        "--root",
        str(root),
        "--directory",
        name,
    ]
    processes = [
        subprocess.Popen(base_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for _ in range(2)
    ]
    results = [process.communicate(timeout=30) + (process.returncode,) for process in processes]
    assert all(returncode == 0 for _, _, returncode in results), results
    assert not raw.exists()
    assert STORE.verify_archive(root, name)["files"] == 200

    materialize_command = base_command.copy()
    materialize_command[2] = "materialize"
    processes = [
        subprocess.Popen(
            materialize_command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        for _ in range(2)
    ]
    results = [process.communicate(timeout=30) + (process.returncode,) for process in processes]
    assert all(returncode == 0 for _, _, returncode in results), results
    assert len(list(raw.glob("*.nwk"))) == 200
    assert not (root / f"{name}.zip").exists()


def test_standard_unzip_accepts_archive(tmp_path: Path):
    unzip = "/usr/bin/unzip"
    if not Path(unzip).is_file():
        pytest.skip("system unzip is unavailable")
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_dna"
    raw = root / name
    raw.mkdir(parents=True)
    (raw / "BUSCO1.dna.nwk").write_text("(A,B);\n", encoding="utf-8")
    STORE.pack_directory(root, name)

    completed = subprocess.run(
        [unzip, "-t", str(root / f"{name}.zip")],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "No errors detected" in completed.stdout


def test_count_reads_raw_zip_and_mixed_top_level_files_without_materializing(
    tmp_path: Path,
):
    root = tmp_path / "species_tree"
    name = "single_copy_iqtree_dna"
    raw = root / name
    (raw / "nested").mkdir(parents=True)
    (raw / "BUSCO1.dna.nwk").write_text("(A,B);\n", encoding="utf-8")
    (raw / ".hidden.nwk").write_text("hidden\n", encoding="utf-8")
    (raw / "nested" / "BUSCO2.dna.nwk").write_text("(A,C);\n", encoding="utf-8")

    assert STORE.count_matching_files(root, name, "*.nwk") == 1
    STORE.pack_directory(root, name)
    archive_path = root / f"{name}.zip"
    archive_stat = archive_path.stat()

    completed = subprocess.run(
        [
            sys.executable,
            str(MODULE_PATH),
            "count",
            "--root",
            str(root),
            "--directory",
            name,
            "--pattern",
            "*.nwk",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    assert completed.stdout.strip() == "1"
    assert not raw.exists()
    assert archive_path.stat().st_ino == archive_stat.st_ino
    assert archive_path.stat().st_mtime_ns == archive_stat.st_mtime_ns

    raw.mkdir()
    (raw / "BUSCO1.dna.nwk").write_text("live override\n", encoding="utf-8")
    (raw / "BUSCO3.dna.nwk").write_text("(B,C);\n", encoding="utf-8")
    assert STORE.count_matching_files(root, name, "*.nwk") == 2


def test_cli_converts_all_managed_directories_both_ways(tmp_path: Path):
    root = tmp_path / "species_tree"
    for name in ("single_copy_iqtree_dna", "single_copy_trimal"):
        directory = root / name
        directory.mkdir(parents=True)
        (directory / "marker.txt").write_text(name, encoding="utf-8")

    packed = subprocess.run(
        [
            "bash",
            str(WRAPPER_PATH),
            "convert-storage",
            "--root",
            str(root),
            "--to",
            "zip",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    packed_records = json.loads(packed.stdout)
    assert {record["directory"] for record in packed_records} == set(
        STORE.MANAGED_DIRECTORIES
    )
    assert (root / "single_copy_iqtree_dna.zip").is_file()
    assert (root / "single_copy_trimal.zip").is_file()

    subprocess.run(
        [
            "bash",
            str(WRAPPER_PATH),
            "convert-storage",
            "--root",
            str(root),
            "--to",
            "raw",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert (root / "single_copy_iqtree_dna" / "marker.txt").is_file()
    assert (root / "single_copy_trimal" / "marker.txt").is_file()
    assert not (root / "single_copy_iqtree_dna.zip").exists()


def test_multi_directory_pack_continues_after_one_directory_error(tmp_path: Path):
    root = tmp_path / "species_tree"
    root.mkdir()
    with zipfile.ZipFile(root / "single_copy_cds_fasta.zip", "w") as archive:
        archive.writestr("single_copy_cds_fasta/unmanaged.txt", "unmarked")
    raw = root / "single_copy_iqtree_dna"
    raw.mkdir()
    (raw / "BUSCO1.dna.nwk").write_text("(A,B);\n", encoding="utf-8")

    completed = subprocess.run(
        [
            sys.executable,
            str(MODULE_PATH),
            "pack",
            "--root",
            str(root),
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 1
    assert "single_copy_cds_fasta" in completed.stderr
    assert (root / "single_copy_iqtree_dna.zip").is_file()
    assert not raw.exists()


def test_genome_evolution_archives_each_high_count_stage_after_last_use():
    text = CORE_PATH.read_text(encoding="utf-8")
    for name in STORE.MANAGED_DIRECTORIES:
        variable = "dir_" + name.replace("_cds_fasta", "_fasta")
        assert f'species_tree_materialize_directory "${{{variable}}}"' in text
        assert f'species_tree_archive_directory "${{{variable}}}"' in text
    assert "trap cleanup_genome_evolution_on_exit EXIT" in text
    assert "species_tree_archive_managed_directories" in text
    assert "species_tree_materialize_managed_directories_for_files_mode" in text

    for name in STORE.MANAGED_DIRECTORIES:
        variable = "dir_" + name.replace("_cds_fasta", "_fasta")
        archive_call = f'species_tree_archive_directory "${{{variable}}}"'
        archive_index = text.index(archive_call)
        remaining_references = [
            match.start()
            for match in re.finditer(rf"\$\{{{re.escape(variable)}\}}", text)
            if match.start() >= archive_index + len(archive_call)
        ]
        assert not remaining_references, f"{variable} is referenced after its archive point"


def test_progress_summary_reports_managed_species_tree_storage():
    text = PROGRESS_CORE_PATH.read_text(encoding="utf-8")
    assert 'dir_species_tree="${gg_workspace_output_dir}/species_tree"' in text
    assert 'species_tree_output_store.py" status' in text


def test_invalid_compression_level_is_rejected():
    with pytest.raises(STORE.SpeciesTreeArchiveError, match="from 0 through 9"):
        STORE._validate_options("adaptive", 10)

"""Exercise native input discovery without starting scientific processing."""

import os
import subprocess
from pathlib import Path

import pytest

CORE = Path(__file__).resolve().parents[1] / "core/gg_transcriptome_generation_core.sh"


def discover(root):
    return subprocess.run(
        ["bash", "-c", 'GG_CORE_SOURCE_ONLY=1 source "$1"; discover_native_fastq_files "$2"',
         "test", str(CORE), str(root)],
        capture_output=True, text=True, check=False,
    )


def test_valid_file_links_preserve_input_names_and_sort_order(tmp_path):
    source = tmp_path / "outside"
    source.mkdir()
    reads = tmp_path / "species"
    reads.mkdir()
    (source / "raw-one").write_bytes(b"@read\nA\n+\nI\n")
    (source / "raw-two").write_bytes(b"@read\nT\n+\nI\n")
    (reads / "RUN_2.fastq.gz").symlink_to(source / "raw-two")
    (reads / "RUN_1.fastq.gz").symlink_to("../outside/raw-one")
    (reads / "RUN_3.fq").write_bytes(b"@read\nG\n+\nI\n")
    result = discover(reads)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == [str(reads / name) for name in (
        "RUN_1.fastq.gz", "RUN_2.fastq.gz", "RUN_3.fq")]
    assert (reads / "RUN_1.fastq.gz").is_symlink()


def test_filters_hidden_nested_and_unsupported_entries(tmp_path):
    (tmp_path / "nested").mkdir()
    (tmp_path / "nested/ignored.fastq").write_bytes(b"nested")
    (tmp_path / "directory.fastq").mkdir()
    for name in (".hidden.fastq", "notes.txt", "upper.FASTQ"):
        (tmp_path / name).write_bytes(b"ignored")
    for name in ("a.fq", "b.fastq", "c.fq.gz", "d.fastq.gz"):
        (tmp_path / name).write_bytes(b"read")
    result = discover(tmp_path)
    assert result.returncode == 0, result.stderr
    assert [Path(path).name for path in result.stdout.splitlines()] == [
        "a.fq", "b.fastq", "c.fq.gz", "d.fastq.gz"]


@pytest.mark.parametrize("target_kind", ["missing", "directory", "loop"])
def test_invalid_fastq_link_fails_without_partial_inventory(tmp_path, target_kind):
    (tmp_path / "a.fastq").write_bytes(b"valid")
    target = tmp_path / "target"
    if target_kind == "directory":
        target.mkdir()
    elif target_kind == "loop":
        target.symlink_to("z.fastq.gz")
    (tmp_path / "z.fastq.gz").symlink_to(target)
    result = discover(tmp_path)
    assert result.returncode != 0
    assert result.stdout == ""
    assert "Native FASTQ input discovery failed:" in result.stderr


def test_species_directory_link_and_line_break_are_rejected(tmp_path):
    real = tmp_path / "real"
    real.mkdir()
    link = tmp_path / "link"
    link.symlink_to(real, target_is_directory=True)
    assert discover(link).returncode != 0
    (real / "bad\nname.fastq").write_bytes(b"read")
    assert "line break" in discover(real).stderr


def test_empty_inventory_is_empty_and_non_regular_file_is_ignored(tmp_path):
    os.mkfifo(tmp_path / "pipe.fastq")
    result = discover(tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout == ""

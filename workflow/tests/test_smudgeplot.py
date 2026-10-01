import argparse
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "support/run_smudgeplot.py"
spec = importlib.util.spec_from_file_location("run_smudgeplot", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def write_reads(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("@read\nACGT\n+\nIIII\n")
    return path


def analysis_args(tmp_path, **overrides):
    values = dict(reads=[str(write_reads(tmp_path / "reads.fastq"))],
                  output_dir=str(tmp_path / "result"), database_dir=None,
                  kmer_length=21, lower_count="4", aggregation_distance=2,
                  threads=2, memory_gb=2, title="Fixture")
    values.update(overrides)
    return argparse.Namespace(**values)


def fake_tools(monkeypatch, *, fail_at=None, cutoff="4", corrupt=None, mutate=None):
    calls = []
    monkeypatch.setattr(module.shutil, "which", lambda name: f"/bin/{name}")

    def run(command, *, cwd, env, check, stdin, stdout, stderr, text):
        assert check and text and stdin == subprocess.DEVNULL
        assert env["MPLBACKEND"] == "Agg"
        calls.append(command)
        task = command[1] if command[0] == "smudgeplot" else command[0]
        if task == fail_at:
            raise subprocess.CalledProcessError(1, command)
        if task == "FastK":
            database = Path(next(arg[2:] for arg in command if arg.startswith("-N")))
            database.with_suffix(".ktab").write_bytes(b"temporary table")
        elif task == "Histex":
            stdout.write("1\t20\n4\t10\n")
        elif task == "hetmers":
            (cwd / "kmerpairs.smu").write_text("12\t12\t10000\n")
        elif task == "all":
            for suffix in ("", "_log10"):
                (cwd / f"smudgeplot_smudgeplot{suffix}.pdf").write_bytes(
                    b"not a PDF" if corrupt == "pdf" else b"%PDF-1.7\nfixture\n")
            (cwd / "smudgeplot_smudgeplot_report.json").write_text(
                "invalid json" if corrupt == "json" else '{"haploid_coverage": 12, "smudges": []}')
            if mutate:
                mutate.write_text("changed during counting")
        return subprocess.CompletedProcess(command, 0, stdout=cutoff if task == "cutoff" else "")

    monkeypatch.setattr(module.subprocess, "run", run)
    return calls


@pytest.mark.parametrize("alias", ["same_path", "symlink", "hardlink"])
def test_duplicate_reads_fail_before_counting(tmp_path, alias):
    reads = write_reads(tmp_path / "reads.fastq")
    other = tmp_path / "other.fastq"
    if alias == "symlink":
        other.symlink_to(reads)
    elif alias == "hardlink":
        os.link(reads, other)
    else:
        other = reads
    with pytest.raises(ValueError, match="duplicate read inputs"):
        module.run_analysis(argparse.Namespace(reads=[str(reads), str(other)]))


def test_empty_read_input_fails_before_counting(tmp_path):
    reads = tmp_path / "reads.fastq"
    reads.touch()
    with pytest.raises(ValueError, match="nonempty file"):
        module.run_analysis(argparse.Namespace(reads=[str(reads)]))


@pytest.mark.parametrize("value", ["0", "-4", "bad"])
def test_resource_and_cutoff_values_must_be_positive(value):
    with pytest.raises(argparse.ArgumentTypeError):
        module.positive_int(value)


def test_discovery_follows_links_but_excludes_hidden_and_unrelated_files(tmp_path):
    source = write_reads(tmp_path / "source" / "first.fastq.gz")
    second = write_reads(tmp_path / "reads" / "nested" / "second.fq")
    (tmp_path / "reads" / "link.fastq.gz").symlink_to(source)
    write_reads(tmp_path / "reads" / ".hidden.fastq")
    write_reads(tmp_path / "reads" / ".cache" / "hidden.fastq")
    write_reads(tmp_path / "reads" / "notes.txt")
    (tmp_path / "alias").symlink_to(tmp_path / "reads", target_is_directory=True)
    assert module.discover_reads(tmp_path / "alias") == sorted([source.resolve(), second.resolve()])


def test_discovery_rejects_broken_fastq_link(tmp_path):
    (tmp_path / "broken.fastq").symlink_to(tmp_path / "absent.fastq")
    with pytest.raises(OSError):
        module.discover_reads(tmp_path)


@pytest.mark.parametrize("cycle", [True, False])
def test_discovery_rejects_directory_cycles_and_repeated_directory_aliases(tmp_path, cycle):
    root = tmp_path / "reads"
    write_reads(root / "nested" / "reads.fastq")
    (root / "alias").symlink_to(root if cycle else root / "nested", target_is_directory=True)
    with pytest.raises(ValueError, match="repeated or cyclic"):
        module.discover_reads(root)


def test_discovery_propagates_directory_errors(tmp_path, monkeypatch):
    def broken_walk(root, *, followlinks, onerror):
        onerror(PermissionError("cannot list reads"))
        yield
    monkeypatch.setattr(module.os, "walk", broken_walk)
    with pytest.raises(PermissionError, match="cannot list reads"):
        module.discover_reads(tmp_path)


def test_list_cli_preserves_filename_boundaries_and_needs_no_runtime(tmp_path):
    path = write_reads(tmp_path / "line\nbreak.fastq")
    result = subprocess.run([sys.executable, str(SCRIPT), "--list-reads", str(tmp_path)],
                            capture_output=True, check=True)
    assert result.stdout == os.fsencode(path.resolve()) + b"\0"


@pytest.mark.parametrize("field", ["reads", "output_dir", "database_dir"])
def test_unsafe_native_paths_fail_before_creating_outputs(tmp_path, field):
    args = analysis_args(tmp_path)
    unsafe = tmp_path / "with spaces"
    if field == "reads":
        args.reads = [str(write_reads(unsafe / "read.fastq"))]
    else:
        setattr(args, field, str(unsafe))
    with pytest.raises(ValueError, match="whitespace or shell metacharacters"):
        module.run_analysis(args)
    assert not Path(args.output_dir).exists()


@pytest.mark.parametrize("missing", ["FastK", "Histex", "Logex", "Symmex", "Fastrm", "smudgeplot"])
def test_missing_native_command_fails_before_counting(tmp_path, monkeypatch, missing):
    args = analysis_args(tmp_path)
    monkeypatch.setattr(module.shutil, "which", lambda name: None if name == missing else name)
    with pytest.raises(ValueError, match=f"required command is missing: {missing}"):
        module.run_analysis(args)
    assert not Path(args.output_dir).exists()


@pytest.mark.parametrize("lower", ["auto", "4"])
def test_successful_run_validates_outputs_and_cleans_temporary_database(tmp_path, monkeypatch, lower):
    args = analysis_args(tmp_path, lower_count=lower)
    calls = fake_tools(monkeypatch)
    metadata = module.run_analysis(args)
    assert metadata["lower_count_used"] == 4
    assert metadata["commands"] == calls
    assert not list(Path(args.output_dir).glob("fastk-*"))
    assert json.loads((Path(args.output_dir) / "run_metadata.json").read_text()) == metadata
    assert ("cutoff" in [command[1] for command in calls if command[0] == "smudgeplot"]) == (lower == "auto")


def test_explicit_database_is_retained(tmp_path, monkeypatch):
    args = analysis_args(tmp_path, database_dir=str(tmp_path / "database"))
    fake_tools(monkeypatch)
    module.run_analysis(args)
    assert (tmp_path / "database" / "reads.ktab").is_file()


@pytest.mark.parametrize("fail_at", ["FastK", "Histex", "hetmers", "all"])
def test_failure_cleans_temporary_database_and_does_not_claim_success(tmp_path, monkeypatch, fail_at):
    args = analysis_args(tmp_path)
    fake_tools(monkeypatch, fail_at=fail_at)
    with pytest.raises(subprocess.CalledProcessError):
        module.run_analysis(args)
    assert not list(Path(args.output_dir).glob("fastk-*"))
    assert not (Path(args.output_dir) / "run_metadata.json").exists()


@pytest.mark.parametrize("cutoff", ["", "bad", "0", "-1", "4\n5"])
def test_invalid_auto_cutoff_stops_before_pair_counting(tmp_path, monkeypatch, cutoff):
    args = analysis_args(tmp_path, lower_count="auto")
    calls = fake_tools(monkeypatch, cutoff=cutoff)
    with pytest.raises(ValueError, match="cutoff did not return a positive integer"):
        module.run_analysis(args)
    assert not any(command[:2] == ["smudgeplot", "hetmers"] for command in calls)


@pytest.mark.parametrize("corrupt", ["pdf", "json"])
def test_invalid_nonempty_outputs_are_rejected(tmp_path, monkeypatch, corrupt):
    args = analysis_args(tmp_path)
    fake_tools(monkeypatch, corrupt=corrupt)
    with pytest.raises(ValueError):
        module.run_analysis(args)
    assert not (Path(args.output_dir) / "run_metadata.json").exists()


def test_input_mutation_is_not_recorded_as_a_successful_run(tmp_path, monkeypatch):
    args = analysis_args(tmp_path)
    fake_tools(monkeypatch, mutate=Path(args.reads[0]))
    with pytest.raises(ValueError, match="changed during analysis"):
        module.run_analysis(args)
    assert not (Path(args.output_dir) / "run_metadata.json").exists()


@pytest.mark.parametrize("directory", ["output_dir", "database_dir"])
def test_nonempty_directories_preserve_prior_data(tmp_path, monkeypatch, directory):
    args = analysis_args(tmp_path, database_dir=str(tmp_path / "database"))
    fake_tools(monkeypatch)
    path = Path(getattr(args, directory))
    path.mkdir()
    (path / "precious").write_bytes(b"old result")
    with pytest.raises(ValueError, match="directory must be empty"):
        module.run_analysis(args)
    assert (path / "precious").read_bytes() == b"old result"
    assert not (Path(args.output_dir) / "commands.log").exists()


def test_concurrent_output_claim_does_not_overwrite_other_run(tmp_path, monkeypatch):
    args = analysis_args(tmp_path)
    calls = fake_tools(monkeypatch)
    original_open = Path.open

    def competing_open(path, mode="r", *positional, **keywords):
        if path.name == "commands.log":
            with original_open(path, "w") as handle:
                handle.write("other run owns this output\n")
        return original_open(path, mode, *positional, **keywords)

    monkeypatch.setattr(Path, "open", competing_open)
    with pytest.raises(FileExistsError):
        module.run_analysis(args)
    assert not calls
    assert (Path(args.output_dir) / "commands.log").read_text() == "other run owns this output\n"


def test_auto_inference_failure_preserves_error_and_does_not_lower_cutoff(tmp_path, monkeypatch):
    args = analysis_args(tmp_path, lower_count="auto")
    calls = fake_tools(monkeypatch, fail_at="all", cutoff="10")
    with pytest.raises(ValueError, match="automatic cutoff 10") as failure:
        module.run_analysis(args)
    assert isinstance(failure.value.__cause__, subprocess.CalledProcessError)
    assert "--lower-count" in str(failure.value)
    assert sum(command[:2] == ["smudgeplot", "hetmers"] for command in calls) == 1
    assert not (Path(args.output_dir) / "run_metadata.json").exists()


def test_same_size_input_mutation_with_restored_mtime_is_rejected(tmp_path, monkeypatch):
    args = analysis_args(tmp_path)
    reads = Path(args.reads[0])
    original = reads.stat()
    run = fake_tools(monkeypatch)
    fake_run = module.subprocess.run

    def mutate_then_restore_mtime(command, **keywords):
        result = fake_run(command, **keywords)
        if command[:2] == ["smudgeplot", "all"]:
            reads.write_bytes(reads.read_bytes().replace(b"ACGT", b"TGCA"))
            os.utime(reads, ns=(original.st_atime_ns, original.st_mtime_ns))
        return result

    monkeypatch.setattr(module.subprocess, "run", mutate_then_restore_mtime)
    with pytest.raises(ValueError, match="changed during analysis"):
        module.run_analysis(args)
    assert run
    assert not (Path(args.output_dir) / "run_metadata.json").exists()

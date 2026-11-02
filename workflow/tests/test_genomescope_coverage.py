import json
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from genomescope_coverage import integer_coverage, plan, run_fits


def archive(tmp_path, *, h=12.01, k=21):
    read = tmp_path / "reads.fastq"
    read.write_text("synthetic input identity\n")
    stat = read.stat()
    metadata = dict(kmer_length=k, lower_count_requested="auto", lower_count_used=10,
                    reads=[dict(path=str(read.resolve()), size_bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)])
    path = tmp_path / "smudgeplot.zip"
    with zipfile.ZipFile(path, "w") as handle:
        handle.writestr("Species.smudgeplot/run_metadata.json", json.dumps(metadata))
        handle.writestr("Species.smudgeplot/smudgeplot_smudgeplot_report.json", json.dumps(dict(haploid_coverage=h)))
    return path, read


def test_native_h_drives_initialization_and_sensitivity_without_inferring_ploidy(tmp_path):
    path, read = archive(tmp_path)
    result = plan(-1, "smudgeplot", path, "Species", 21, [read])
    assert result["initial_coverage"] == 12
    assert result["smudgeplot"]["haploid_coverage"] == 12.01
    assert len(result["smudgeplot"]["archive_sha256"]) == 64
    assert [fit["initial_coverage"] for fit in result["fits"]] == [12, 6, 24, -1]
    assert "ploidy" not in result


@pytest.mark.parametrize("value", [None, True, "12", 0, -1, float("nan"), float("inf"), 2147483648])
def test_unusable_h_is_rejected(tmp_path, value):
    path, read = archive(tmp_path, h=value)
    with pytest.raises(ValueError):
        plan(-1, "smudgeplot", path, "Species", 21, [read])


def test_rounding_and_redundant_sensitivity_starts(tmp_path):
    assert integer_coverage(6.5) == 7
    path, read = archive(tmp_path, h=1.01)
    assert [fit["initial_coverage"] for fit in plan(-1, "smudgeplot", path, "Species", 21, [read])["fits"]] == [1, 2, -1]


@pytest.mark.parametrize("mismatch", ["k", "reads", "changed", "empty"])
def test_unmatched_smudgeplot_cannot_supply_h(tmp_path, mismatch):
    path, read = archive(tmp_path)
    if mismatch == "changed":
        read.write_text("modified\n")
    reads = [read]
    if mismatch == "reads":
        other = tmp_path / "other.fastq"
        other.write_text("other")
        reads = [other]
    if mismatch == "empty":
        reads = []
    with pytest.raises(ValueError, match="same"):
        plan(-1, "smudgeplot", path, "Species", 31 if mismatch == "k" else 21, reads)


def test_manual_override_and_native_auto_do_not_need_smudgeplot():
    assert plan(66, "smudgeplot", None, "Species", 21, [])["selected_source"] == "manual"
    assert plan(-1, "auto", None, "Species", 21, [])["initial_coverage"] == -1
    with pytest.raises(ValueError, match="requires"):
        plan(-1, "smudgeplot", None, "Species", 21, [])
    with pytest.raises(ValueError, match="source"):
        plan(-1, "other", None, "Species", 21, [])


def test_failed_sensitivity_is_retained_and_does_not_replace_primary(tmp_path, monkeypatch):
    path, read = archive(tmp_path)
    selection = plan(-1, "smudgeplot", path, "Species", 21, [read])
    output = tmp_path / "out"
    def fake_run(command, **kwargs):
        directory = Path(command[command.index("-o") + 1])
        initial = int(command[command.index("-l") + 1])
        assert command[command.index("-p") + 1] == "4"
        if initial == 6:
            return subprocess.CompletedProcess(command, 1)
        for name in ("summary.txt", "transformed_linear_plot.png", "model.txt"):
            (directory / name).write_text(str(initial))
        (directory / "report.json").write_text(json.dumps(dict(model_summary={"parameters": {}}, kcov=12 if initial == -1 else initial, genome_haploid_length={"avg": 1000000})))
        return subprocess.CompletedProcess(command, 0)
    monkeypatch.setattr(subprocess, "run", fake_run)
    run_fits(selection, output, read, "genomescope.R", 21, 4, -1)
    results = json.loads((output / "coverage_selection.json").read_text())["results"]
    assert [result["status"] for result in results] == ["success", "failed", "success", "success"]
    assert (output / "summary.txt").read_text() == "12"
    # Native tools can exit zero without producing a fitted model.
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs: subprocess.CompletedProcess(args, 0))
    with pytest.raises(ValueError, match="primary fit failed"):
        run_fits(selection, tmp_path / "failed", read, "genomescope.R", 21, 2, -1)
    assert (tmp_path / "failed/coverage_selection.json").exists()


def test_reads_changed_during_counting_cannot_publish(tmp_path):
    path, read = archive(tmp_path)
    selection = plan(-1, "smudgeplot", path, "Species", 21, [read])
    read.write_text("changed during counting")
    with pytest.raises(ValueError, match="changed during"):
        run_fits(selection, tmp_path / "out", read, "missing_program", 21, 2, -1)


@pytest.mark.parametrize("failure", ["missing_program", "invalid_json"])
def test_execution_and_json_errors_are_retained_as_failed_fits(tmp_path, monkeypatch, failure):
    path, read = archive(tmp_path)
    selection = plan(12, "smudgeplot", None, "Species", 21, [])
    def fake_run(command, **kwargs):
        if failure == "missing_program":
            raise FileNotFoundError("program unavailable")
        directory = Path(command[command.index("-o") + 1])
        for name in ("summary.txt", "transformed_linear_plot.png", "model.txt", "report.json"):
            (directory / name).write_text("invalid JSON")
        return subprocess.CompletedProcess(command, 0)
    monkeypatch.setattr(subprocess, "run", fake_run)
    with pytest.raises(ValueError, match="primary fit failed"):
        run_fits(selection, tmp_path / "out", read, "program", 21, 2, -1)
    result = json.loads((tmp_path / "out/coverage_selection.json").read_text())["results"][0]
    assert result["status"] == "failed" and result["error"]


@pytest.mark.parametrize("location", ["primary", "twice_H"])
def test_stale_fit_outputs_are_rejected_before_any_fit_runs(tmp_path, monkeypatch, location):
    path, read = archive(tmp_path)
    selection = plan(-1, "smudgeplot", path, "Species", 21, [read])
    output = tmp_path / "out"
    directory = output if location == "primary" else output / "sensitivity" / location
    directory.mkdir(parents=True)
    stale = directory / "report.json"
    stale.write_text("old fitted output")
    def unexpected_run(*args, **kwargs):
        pytest.fail("No fits should run when a previous result would be overwritten")
    monkeypatch.setattr(subprocess, "run", unexpected_run)
    with pytest.raises(ValueError, match="already exists"):
        run_fits(selection, output, read, "program", 21, 2, -1)
    assert stale.read_text() == "old fitted output"
    assert not (output / "command.txt").exists()


@pytest.mark.parametrize("bad_length", [None, [], "unknown", {"avg": None}, {"avg": float("nan")}])
def test_malformed_sensitivity_report_retains_diagnostics_and_primary(tmp_path, monkeypatch, bad_length):
    path, read = archive(tmp_path)
    selection = plan(-1, "smudgeplot", path, "Species", 21, [read])
    output = tmp_path / "out"
    def fake_run(command, **kwargs):
        directory = Path(command[command.index("-o") + 1])
        initial = int(command[command.index("-l") + 1])
        for name in ("summary.txt", "transformed_linear_plot.png", "model.txt"):
            (directory / name).write_text(str(initial))
        report = dict(model_summary={"parameters": {}}, kcov=12,
                      genome_haploid_length=bad_length if initial == 6 else {"avg": 1000000})
        (directory / "report.json").write_text(json.dumps(report))
        return subprocess.CompletedProcess(command, 0)
    monkeypatch.setattr(subprocess, "run", fake_run)
    run_fits(selection, output, read, "program", 21, 2, -1)
    results = json.loads((output / "coverage_selection.json").read_text())["results"]
    assert [result["status"] for result in results] == ["success", "failed", "success", "success"]
    assert results[1]["error"]
    assert (output / "summary.txt").read_text() == "12"

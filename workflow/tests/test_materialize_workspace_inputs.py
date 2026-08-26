import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "support" / "materialize_workspace_inputs.py"
SPEC = importlib.util.spec_from_file_location("materialize_workspace_inputs", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def _workspace(tmp_path):
    workspace = tmp_path / "workspace"
    generated = workspace / "output" / "input_generation"
    inputs = workspace / "input"
    inputs.mkdir(parents=True)
    for name in ("species_cds", "species_gff"):
        source = generated / name
        source.mkdir(parents=True)
        (source / f"{name}.txt").write_text(f"{name}\n", encoding="utf-8")
        (inputs / name).symlink_to(Path("../output/input_generation") / name)
    return workspace


def _run(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *map(str, args)],
        text=True,
        capture_output=True,
        check=False,
    )


def test_apply_materializes_internal_symlinks_and_preserves_recovery_links(tmp_path):
    workspace = _workspace(tmp_path)
    backup = tmp_path / "recovery"
    completed = _run(
        "apply",
        "--workspace",
        workspace,
        "--backup-root",
        backup,
        "--name",
        "species_cds",
        "--name",
        "species_gff",
    )
    assert completed.returncode == 0, completed.stdout
    payload = json.loads(completed.stdout)
    assert payload["status"] == "ok"

    for name in ("species_cds", "species_gff"):
        destination = workspace / "input" / name
        assert destination.is_dir()
        assert not destination.is_symlink()
        assert (backup / f"{name}.original-symlink").is_symlink()
        assert (backup / f"{name}.receipt.json").is_file()
        source_file = workspace / "output" / "input_generation" / name / f"{name}.txt"
        destination_file = destination / f"{name}.txt"
        assert destination_file.read_bytes() == source_file.read_bytes()
        assert os.stat(destination_file).st_ino == os.stat(source_file).st_ino

    verified = _run(
        "verify",
        "--workspace",
        workspace,
        "--backup-root",
        backup,
        "--name",
        "species_cds",
        "--name",
        "species_gff",
    )
    assert verified.returncode == 0, verified.stdout


def test_plan_rejects_symlink_that_does_not_target_generated_input(tmp_path):
    workspace = _workspace(tmp_path)
    wrong = workspace / "wrong"
    wrong.mkdir()
    destination = workspace / "input" / "species_cds"
    destination.unlink()
    destination.symlink_to(wrong)

    completed = _run(
        "plan",
        "--workspace",
        workspace,
        "--backup-root",
        tmp_path / "recovery",
        "--name",
        "species_cds",
    )
    assert completed.returncode == 1
    assert "does not point to its generated source" in completed.stdout

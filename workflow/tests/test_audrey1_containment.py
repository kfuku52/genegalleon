import json
import shlex
import subprocess
from pathlib import Path

import pytest
from shell_static_helpers import WORKFLOW_DIR

from workflow.support.local_input_manifest_binds import source_files

UTIL = WORKFLOW_DIR / "support" / "gg_util.sh"


def run_site_command(tmp_path: Path, bind: str, setup: str = "") -> subprocess.CompletedProcess[str]:
    command = (
        f"source {shlex.quote(str(UTIL))}; "
        "hostname() { printf 'audrey1\\n'; }; "
        f"export GG_CONTAINER_PROJECT_ROOT_BIND={shlex.quote(bind)}; "
        + setup + "gg_site_container_shell_command singularity singularity_command; "
        "printf 'command=%s\\n' \"${singularity_command[*]}\"; "
        "printf 'bind=%s\\n' \"${SINGULARITY_BINDPATH:-}\""
    )
    return subprocess.run(
        ["bash", "-euo", "pipefail", "-c", command],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )


def test_audrey1_isolates_shm_and_keeps_project_paths_visible(tmp_path):
    root = tmp_path / "project"
    root.mkdir()
    result = run_site_command(tmp_path, f"{root}:{root}")
    assert result.returncode == 0, result.stderr
    assert "site profile = audrey1" in result.stdout
    assert "command=singularity exec --contain" in result.stdout
    assert f"bind={root}:{root}" in result.stdout


def test_audrey1_rejects_unreviewed_bind(tmp_path):
    root = tmp_path / "project"
    root.mkdir()
    result = run_site_command(tmp_path, f"{root}:/unrelated")
    assert result.returncode != 0
    assert "invalid GG_CONTAINER_PROJECT_ROOT_BIND" in result.stderr


def test_audrey1_existing_unbound_runs_keep_previous_runtime(tmp_path):
    result = run_site_command(tmp_path, "")
    assert result.returncode == 0, result.stderr
    assert "command=singularity exec\n" in result.stdout


def test_native_array_mounts_only_declared_sources_read_only(tmp_path):
    root = tmp_path / "project"
    root.mkdir()
    workspace = root / "work"
    workspace.mkdir()
    external = tmp_path / "source with spaces.fa"
    external.write_text(">seq\nATG\n")
    (external.parent / "unrelated.fa").write_text("private")
    manifest = workspace / "download.tsv"
    manifest.write_text("provider\tid\tcds_url\nlocal\tTest_species\t" + external.as_uri() + "\n")
    setup = (f"gg_workspace_dir={shlex.quote(str(workspace))}; "
             "export GG_INPUT_INPUT_GENERATION_MODE=array_prepare GG_INPUT_DOWNLOAD_MANIFEST=/workspace/download.tsv; ")
    result = run_site_command(tmp_path, f"{root}:{root}", setup)
    assert result.returncode == 0, result.stderr
    assert f"{external}:{external}:ro" in result.stdout
    assert "unrelated.fa" not in result.stdout
    external.unlink()
    result = run_site_command(tmp_path, f"{root}:{root}", setup)
    assert result.returncode != 0
    assert "No such file" in result.stderr


def test_non_array_project_run_does_not_mount_manifest_sources(tmp_path):
    root = tmp_path / "project"
    root.mkdir()
    result = run_site_command(tmp_path, f"{root}:{root}",
                              "export GG_INPUT_INPUT_GENERATION_MODE=single GG_INPUT_DOWNLOAD_MANIFEST=/missing; ")
    assert result.returncode == 0, result.stderr


def test_frozen_plan_sources_are_exact_files_and_deduplicated(tmp_path):
    source = tmp_path / "input.fa"
    source.write_text(">seq\nATG\n")
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({"tasks": [{"manifest_row": {"cds_url": source.as_uri(),
        "gff_url": "https://example.org/annotation.gff"}, "cds_path": "/workspace/input.fa",
        "input_sha256": {str(source): "sealed-digest"}}]}))
    assert source_files(plan, tmp_path, plan=True) == [source]


@pytest.mark.parametrize("name", ["unsafe,source.fa", "unsafe:source.fa", "unsafe\nsource.fa"])
def test_local_binds_reject_argument_delimiters(tmp_path, name):
    source = tmp_path / name
    source.touch()
    manifest = tmp_path / "plan.tsv"
    manifest.write_text("cds_url\n" + source.as_uri() + "\n")
    with pytest.raises(ValueError, match="safely"):
        source_files(manifest, tmp_path)


def test_local_binds_reject_directory_and_remote_file_authority(tmp_path):
    manifest = tmp_path / "plan.tsv"
    manifest.write_text("local_cds_path\n" + str(tmp_path) + "\n")
    with pytest.raises(ValueError, match="regular file"):
        source_files(manifest, tmp_path)
    manifest.write_text("cds_url\nfile://remote-host/data.fa\n")
    with pytest.raises(ValueError, match="this host"):
        source_files(manifest, tmp_path)

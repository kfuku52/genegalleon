from pathlib import Path
import shlex
import subprocess


GG_UTIL_PATH = Path(__file__).resolve().parents[1] / "script" / "gg_util.sh"


def run_bash(cmd: str, cwd: Path):
    return subprocess.run(
        ["bash", "-lc", cmd],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )


def test_workspace_pfam_le_dir_is_under_db_dedicated_folder(tmp_path):
    project_dir = tmp_path / "project"
    project_dir.mkdir()
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"workspace_pfam_le_dir {shlex.quote(str(project_dir))}"
    )

    completed = run_bash(command, cwd=tmp_path)
    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == str(project_dir / "db" / "pfam" / "Pfam_LE")


def test_ensure_pfam_le_db_migrates_legacy_workspace_layout(tmp_path):
    project_dir = tmp_path / "project"
    legacy_dir = project_dir / "db" / "Pfam_LE"
    new_dir = project_dir / "db" / "pfam" / "Pfam_LE"
    legacy_dir.mkdir(parents=True)
    (legacy_dir / "Pfam.pal").write_text("dummy\n")

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"ensure_pfam_le_db {shlex.quote(str(project_dir))}"
    )
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == str(new_dir)
    assert (new_dir / "Pfam.pal").is_file()
    assert not legacy_dir.exists()


def test_mv_out_accepts_pipe_input(tmp_path):
    outfile = tmp_path / "nested" / "out.txt"
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"printf 'hello\\n' | mv_out {shlex.quote(str(outfile))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert outfile.read_text() == "hello\n"


def test_cp_out_accepts_pipe_input(tmp_path):
    outfile = tmp_path / "nested" / "out.txt"
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"printf 'world\\n' | cp_out {shlex.quote(str(outfile))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert outfile.read_text() == "world\n"


def test_cp_out_single_argument_without_pipe_fails(tmp_path):
    outfile = tmp_path / "out.txt"
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"cp_out {shlex.quote(str(outfile))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode != 0
    assert "at least 2 arguments are required unless stdin is piped" in completed.stdout


def test_ensure_jaspar_file_latest_prefers_highest_local_release(tmp_path):
    project_dir = tmp_path / "project"
    jaspar_dir = project_dir / "db" / "jaspar"
    jaspar_dir.mkdir(parents=True)
    old_file = jaspar_dir / "JASPAR2022_CORE_plants_non-redundant_pfms_meme.txt"
    new_file = jaspar_dir / "JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt"
    old_file.write_text("old\n")
    new_file.write_text("new\n")

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"GG_JASPAR_SKIP_REMOTE_LOOKUP=1; "
        f"ensure_jaspar_file {shlex.quote(str(project_dir))} latest"
    )
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == str(new_file)
    marker = jaspar_dir / "latest_core_plants_non-redundant_pfms_meme.filename"
    assert marker.read_text().strip() == new_file.name


def test_ensure_jaspar_file_latest_uses_cached_marker_without_remote_lookup(tmp_path):
    project_dir = tmp_path / "project"
    jaspar_dir = project_dir / "db" / "jaspar"
    jaspar_dir.mkdir(parents=True)
    filename = "JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt"
    target_file = jaspar_dir / filename
    marker = jaspar_dir / "latest_core_plants_non-redundant_pfms_meme.filename"
    target_file.write_text("cached\n")
    marker.write_text(f"{filename}\n")

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"GG_JASPAR_SKIP_REMOTE_LOOKUP=1; "
        f"ensure_jaspar_file {shlex.quote(str(project_dir))} auto"
    )
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == str(target_file)

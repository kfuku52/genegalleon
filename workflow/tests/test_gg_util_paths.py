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

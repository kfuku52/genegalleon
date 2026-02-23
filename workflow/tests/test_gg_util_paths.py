from pathlib import Path
import shlex
import subprocess


GG_UTIL_PATH = Path(__file__).resolve().parents[1] / "support" / "gg_util.sh"


def run_bash(cmd: str, cwd: Path):
    return subprocess.run(
        ["bash", "-lc", cmd],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )


def test_workspace_pfam_le_dir_is_under_downloads_dedicated_folder(tmp_path):
    project_dir = tmp_path / "project"
    project_dir.mkdir()
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"workspace_pfam_le_dir {shlex.quote(str(project_dir))}"
    )

    completed = run_bash(command, cwd=tmp_path)
    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == str(project_dir / "downloads" / "pfam" / "Pfam_LE")


def test_ensure_pfam_le_db_uses_new_workspace_layout_without_migrating_legacy_dir(tmp_path):
    project_dir = tmp_path / "project"
    legacy_dir = project_dir / "downloads" / "Pfam_LE"
    new_dir = project_dir / "downloads" / "pfam" / "Pfam_LE"
    legacy_dir.mkdir(parents=True)
    new_dir.mkdir(parents=True)
    (legacy_dir / "Pfam.pal").write_text("dummy\n")
    (new_dir / "Pfam.pal").write_text("dummy\n")

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"ensure_pfam_le_db {shlex.quote(str(project_dir))}"
    )
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == str(new_dir)
    assert (new_dir / "Pfam.pal").is_file()
    assert legacy_dir.is_dir()


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


def test_cp_out_creates_destination_dir_when_target_has_trailing_slash(tmp_path):
    src = tmp_path / "src.txt"
    src.write_text("abc\n")
    dest_dir = tmp_path / "nested" / "dest_dir"
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"cp_out {shlex.quote(str(src))} {shlex.quote(str(dest_dir) + '/')}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    copied = dest_dir / src.name
    assert copied.read_text() == "abc\n"


def test_mv_out_creates_destination_dir_for_multi_source_move(tmp_path):
    src1 = tmp_path / "a.txt"
    src2 = tmp_path / "b.txt"
    src1.write_text("a\n")
    src2.write_text("b\n")
    dest_dir = tmp_path / "nested" / "mv_dest"
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"mv_out {shlex.quote(str(src1))} {shlex.quote(str(src2))} {shlex.quote(str(dest_dir))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert not src1.exists()
    assert not src2.exists()
    assert (dest_dir / "a.txt").read_text() == "a\n"
    assert (dest_dir / "b.txt").read_text() == "b\n"


def test_ensure_jaspar_file_latest_prefers_highest_local_release(tmp_path):
    project_dir = tmp_path / "project"
    jaspar_dir = project_dir / "downloads" / "jaspar"
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
    jaspar_dir = project_dir / "downloads" / "jaspar"
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


def test_recreate_dir_rejects_root_path(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "recreate_dir /"
    )
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode != 0
    assert "Refusing to recreate unsafe directory path: /" in completed.stdout


def test_is_species_set_identical_ignores_hidden_files(tmp_path):
    dir1 = tmp_path / "d1"
    dir2 = tmp_path / "d2"
    dir1.mkdir()
    dir2.mkdir()
    (dir1 / "Arabidopsis_thaliana.fa").write_text(">a\nATG\n")
    (dir2 / "Arabidopsis_thaliana.fa").write_text(">a\nATG\n")
    (dir1 / ".DS_Store").write_text("x\n")

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"if is_species_set_identical {shlex.quote(str(dir1))} {shlex.quote(str(dir2))}; then echo OK; else echo NG; fi"
    )
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert "OK" in completed.stdout
    assert "NG" not in completed.stdout


def test_is_species_set_identical_returns_nonzero_when_directory_is_missing(tmp_path):
    dir1 = tmp_path / "d1"
    dir2 = tmp_path / "d2"
    dir1.mkdir()

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"is_species_set_identical {shlex.quote(str(dir1))} {shlex.quote(str(dir2))}"
    )
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode != 0
    assert "Directory not found for species-set comparison" in completed.stdout


def test_check_if_species_files_unique_ignores_hidden_files(tmp_path):
    species_dir = tmp_path / "species"
    species_dir.mkdir()
    (species_dir / "Arabidopsis_thaliana.fa").write_text(">a\nATG\n")
    (species_dir / ".DS_Store").write_text("x\n")

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"check_if_species_files_unique {shlex.quote(str(species_dir))}"
    )
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert "Species files are unique in" in completed.stdout


def test_gg_find_fasta_files_excludes_hidden_files(tmp_path):
    fasta_dir = tmp_path / "fasta"
    fasta_dir.mkdir()
    visible = fasta_dir / "Arabidopsis_thaliana.fa"
    hidden = fasta_dir / ".Arabidopsis_thaliana.fa"
    visible.write_text(">a\nATG\n")
    hidden.write_text(">h\nATG\n")

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"gg_find_fasta_files {shlex.quote(str(fasta_dir))} 1"
    )
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    out_lines = [line for line in completed.stdout.splitlines() if line.strip()]
    assert str(visible) in out_lines
    assert str(hidden) not in out_lines

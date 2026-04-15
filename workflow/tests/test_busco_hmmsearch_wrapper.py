from pathlib import Path
import os
import shlex
import stat
import subprocess


REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_DIR = REPO_ROOT / "workflow"
SUPPORT_DIR = WORKFLOW_DIR / "support"
GG_BUSCO_PATH = SUPPORT_DIR / "gg_busco.sh"
HMMSEARCH_WRAPPER_PATH = SUPPORT_DIR / "gg_wrapper_bin" / "hmmsearch"


def _write_executable(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _run_bash(command: str, cwd: Path, env: dict[str, str] | None = None):
    run_env = os.environ.copy()
    if env:
        run_env.update(env)
    return subprocess.run(
        ["bash", "-c", command],
        cwd=cwd,
        capture_output=True,
        text=True,
        env=run_env,
        check=False,
    )


def test_hmmsearch_wrapper_creates_modified_fas_symlink_when_missing(tmp_path: Path):
    real_hmmsearch = tmp_path / "bin" / "hmmsearch-real"
    observed_args = tmp_path / "observed.args"
    _write_executable(
        real_hmmsearch,
        f"""#!/usr/bin/env bash
set -euo pipefail
printf '%s\\n' "$@" > {shlex.quote(str(observed_args))}
[[ -e "${{!#}}" ]]
""",
    )

    metaeuk_dir = tmp_path / "metaeuk_output" / "initial_results"
    metaeuk_dir.mkdir(parents=True)
    fallback_input = metaeuk_dir / "busco_infile_cds.fa.fas"
    expected_input = metaeuk_dir / "busco_infile_cds.fa.modified.fas"
    fallback_input.write_text(">seq1\nAAAA\n", encoding="utf-8")

    completed = subprocess.run(
        [
            "bash",
            str(HMMSEARCH_WRAPPER_PATH),
            "--domtblout",
            str(tmp_path / "domtblout.tsv"),
            "--cpu",
            "1",
            str(tmp_path / "marker.hmm"),
            str(expected_input),
        ],
        capture_output=True,
        text=True,
        env={
            **os.environ,
            "GG_REAL_HMMSEARCH": str(real_hmmsearch),
            "GG_BUSCO_METAEUK_MODIFIED_FAS_COMPAT": "1",
        },
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert expected_input.is_symlink()
    assert expected_input.resolve() == fallback_input.resolve()
    assert observed_args.read_text(encoding="utf-8").splitlines()[-1] == str(expected_input)
    assert "GeneGalleon BUSCO compatibility: created" in completed.stderr


def test_gg_run_busco_with_metaeuk_modified_fas_compat_wraps_busco_env(tmp_path: Path):
    bin_dir = tmp_path / "bin"
    real_hmmsearch = bin_dir / "hmmsearch"
    busco = bin_dir / "busco"
    observed_wrapper = tmp_path / "observed.wrapper"
    observed_real = tmp_path / "observed.real"
    observed_compat = tmp_path / "observed.compat"
    _write_executable(
        real_hmmsearch,
        """#!/usr/bin/env bash
set -euo pipefail
exit 0
""",
    )
    _write_executable(
        busco,
        f"""#!/usr/bin/env bash
set -euo pipefail
command -v hmmsearch > {shlex.quote(str(observed_wrapper))}
printf '%s\\n' "${{GG_REAL_HMMSEARCH:-}}" > {shlex.quote(str(observed_real))}
printf '%s\\n' "${{GG_BUSCO_METAEUK_MODIFIED_FAS_COMPAT:-}}" > {shlex.quote(str(observed_compat))}
""",
    )

    command = (
        f"source {shlex.quote(str(GG_BUSCO_PATH))}; "
        f"gg_support_dir={shlex.quote(str(SUPPORT_DIR))}; "
        "gg_run_busco_with_metaeuk_modified_fas_compat "
        "--in input.fasta --mode transcriptome --out busco_tmp"
    )

    completed = _run_bash(
        command,
        cwd=tmp_path,
        env={
            "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}",
        },
    )

    assert completed.returncode == 0, completed.stderr
    assert observed_wrapper.read_text(encoding="utf-8").strip() == str(HMMSEARCH_WRAPPER_PATH)
    assert observed_real.read_text(encoding="utf-8").strip() == str(real_hmmsearch)
    assert observed_compat.read_text(encoding="utf-8").strip() == "1"


def test_gg_busco_stderr_matches_known_metaeuk_modified_fas_bug(tmp_path: Path):
    stderr_log = tmp_path / "busco.stderr.log"
    stderr_log.write_text(
        "Error: Failed to open sequence file "
        "/tmp/run_eukaryota_odb12/metaeuk_output/initial_results/busco_infile_cds.fa.modified.fas for reading\n",
        encoding="utf-8",
    )

    command = (
        f"source {shlex.quote(str(GG_BUSCO_PATH))}; "
        f"gg_busco_stderr_matches_known_metaeuk_modified_fas_bug {shlex.quote(str(stderr_log))}"
    )

    completed = _run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr

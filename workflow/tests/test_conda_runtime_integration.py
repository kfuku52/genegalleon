import os
import subprocess
import textwrap
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
GG_UTIL = REPO_ROOT / "workflow" / "support" / "gg_util.sh"


def _install_fake_micromamba(tmp_path: Path) -> Path:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    micromamba = bin_dir / "micromamba"
    micromamba.write_text(
        textwrap.dedent(
            """\
            #!/usr/bin/env bash
            if [[ "${1:-}" == "shell" && "${2:-}" == "hook" && "${3:-}" == "--shell" && "${4:-}" == "bash" ]]; then
              cat <<'EOF'
            micromamba() {
              command micromamba "$@"
            }
            EOF
              exit 0
            fi
            if [[ "${1:-}" == "activate" && "${2:-}" == "base" ]]; then
              exit 0
            fi
            if [[ "${1:-}" == "deactivate" ]]; then
              exit 0
            fi
            exit 64
            """
        ),
        encoding="utf-8",
    )
    micromamba.chmod(0o755)
    return bin_dir


def test_conda_shell_reinitializes_in_nested_bash(tmp_path: Path):
    fake_bin = _install_fake_micromamba(tmp_path)
    env = {
        key: value
        for key, value in os.environ.items()
        if not (
            key.startswith("CONDA")
            or key.startswith("MAMBA")
            or key == "GG_CONDA_SHELL_INITIALIZED"
        )
    }
    env["PATH"] = os.pathsep.join((str(fake_bin), "/usr/bin", "/bin"))

    script = textwrap.dedent(
        """\
        set -euo pipefail
        source "$1"
        gg_activate_conda_env base
        declare -F conda >/dev/null
        if env | grep -q '^GG_CONDA_SHELL_INITIALIZED='; then
          echo "The conda initialization marker leaked into the child environment." >&2
          exit 1
        fi

        GG_CONDA_SHELL_INITIALIZED=1 bash --noprofile --norc -c '
          set -euo pipefail
          source "$1"
          if declare -F conda >/dev/null 2>&1; then
            echo "The parent conda function unexpectedly leaked into the child Bash process." >&2
            exit 1
          fi
          gg_activate_conda_env base
          declare -F conda >/dev/null
          if env | grep -q "^GG_CONDA_SHELL_INITIALIZED="; then
            echo "The inherited conda initialization marker remained exported." >&2
            exit 1
          fi
          gg_deactivate_conda_env
        ' nested-bash "$1"
        """
    )
    completed = subprocess.run(
        ["bash", "--noprofile", "--norc", "-c", script, "outer-bash", str(GG_UTIL)],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr

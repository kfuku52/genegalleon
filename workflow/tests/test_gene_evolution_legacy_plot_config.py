"""Legacy project settings must reach terminal stages under Bash nounset."""
import os
import subprocess
from pathlib import Path

import pytest

CORE = Path(__file__).resolve().parents[1] / "core/gg_gene_evolution_core.sh"
PARAMETERS = {
    "synteny_search_window": "20",
    "treevis_sequence_similarity": "1",
    "treevis_sequence_similarity_width_mm": "15",
    "treevis_cis_similarity": "1",
    "treevis_cis_similarity_width_mm": "15",
    "treevis_synteny_similarity": "1",
    "treevis_synteny_similarity_width_mm": "15",
}


@pytest.mark.parametrize("custom", [False, True])
def test_legacy_plot_configuration_initializes_before_analysis(tmp_path, custom):
    env = {key: value for key, value in os.environ.items() if key not in PARAMETERS}
    env.update(GG_CORE_SOURCE_ONLY="1", GG_TASK_CPUS="1", GG_MEM_TOTAL_GB="2",
               gg_workspace_dir=str(tmp_path / "workspace"),
               gg_support_dir=str(CORE.parent.parent / "support"))
    expected = {key: ("0" if "width" not in key and key != "synteny_search_window" else "7")
                for key in PARAMETERS} if custom else PARAMETERS
    if custom:
        env.update(expected)
    # Source the actual core rather than a duplicated configuration fragment.
    # Expanding these terminal-stage arguments used to abort legacy runs late.
    command = 'source "$1"; shift; for name in "$@"; do printf "%s=%s\\n" "$name" "${!name}"; done'
    result = subprocess.run(["bash", "-eu", "-c", command, "_", str(CORE), *PARAMETERS],
                            cwd=tmp_path, env=env, text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr
    observed = dict(line.split("=", 1) for line in result.stdout.splitlines()
                    if line.partition("=")[0] in PARAMETERS)
    assert observed == expected
    assert not (tmp_path / "workspace/output/query2family").exists()

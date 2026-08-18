import subprocess
from pathlib import Path


WORKFLOW_DIR = Path(__file__).resolve().parents[1]
CORE_PATH = WORKFLOW_DIR / "core" / "gg_genome_evolution_core.sh"
EXECUTION_RUNTIME_PATH = WORKFLOW_DIR / "support" / "gg_util" / "07_execution_runtime.sh"


def _run_gating_helpers(tmp_path: Path) -> subprocess.CompletedProcess[str]:
    script = r'''
set -eu
source "$1"
probe_dir=$2

printf '(A:1,B:1);\n' > "${probe_dir}/available.nwk"

any_available=1
disable_if_no_nonempty_input_file any_available \
  "${probe_dir}/missing.nwk" "${probe_dir}/available.nwk"
printf 'any_available=%s\n' "${any_available}"

all_missing=1
disable_if_no_nonempty_input_file all_missing \
  "${probe_dir}/missing-1.nwk" "${probe_dir}/missing-2.nwk"
printf 'all_missing=%s\n' "${all_missing}"

mkdir "${probe_dir}/empty-busco"
empty_busco=1
disable_if_no_matching_input_file empty_busco \
  "${probe_dir}/empty-busco" '*busco.full.tsv'
printf 'empty_busco=%s\n' "${empty_busco}"

mkdir "${probe_dir}/unrelated-busco"
printf 'not a BUSCO table\n' > "${probe_dir}/unrelated-busco/README.txt"
unrelated_busco=1
disable_if_no_matching_input_file unrelated_busco \
  "${probe_dir}/unrelated-busco" '*busco.full.tsv'
printf 'unrelated_busco=%s\n' "${unrelated_busco}"

mkdir "${probe_dir}/empty-table-busco"
: > "${probe_dir}/empty-table-busco/A.busco.full.tsv"
empty_table_busco=1
disable_if_no_matching_input_file empty_table_busco \
  "${probe_dir}/empty-table-busco" '*busco.full.tsv'
printf 'empty_table_busco=%s\n' "${empty_table_busco}"

mkdir "${probe_dir}/valid-busco"
printf 'species\tstatus\n' > "${probe_dir}/valid-busco/A.busco.full.tsv"
valid_busco=1
disable_if_no_matching_input_file valid_busco \
  "${probe_dir}/valid-busco" '*busco.full.tsv'
printf 'valid_busco=%s\n' "${valid_busco}"
'''
    return subprocess.run(
        ["bash", "-s", "--", str(EXECUTION_RUNTIME_PATH), str(tmp_path)],
        input=script,
        text=True,
        capture_output=True,
        check=False,
    )


def test_input_gating_distinguishes_any_file_and_matching_directory_content(
    tmp_path: Path,
):
    completed = _run_gating_helpers(tmp_path)

    assert completed.returncode == 0, completed.stderr + completed.stdout
    assert "any_available=1" in completed.stdout
    assert "all_missing=0" in completed.stdout
    assert "empty_busco=0" in completed.stdout
    assert "unrelated_busco=0" in completed.stdout
    assert "empty_table_busco=0" in completed.stdout
    assert "valid_busco=1" in completed.stdout


def test_species_tree_plot_stages_do_not_mutate_the_shared_user_flag():
    core = CORE_PATH.read_text(encoding="utf-8")
    start = core.index('task="Species tree plotting"')
    end = core.index('task="Time-constrained tree preparation"', start)
    plotting_stages = core[start:end]

    assert "run_species_tree_comparison_plot=${run_plot_species_trees}" in plotting_stages
    assert "run_species_tree_busco_plot=${run_plot_species_trees}" in plotting_stages
    assert 'disable_if_no_input_file "run_plot_species_trees"' not in plotting_stages
    assert 'disable_if_no_nonempty_input_file "run_plot_species_trees"' not in plotting_stages
    assert 'disable_if_no_matching_input_file "run_plot_species_trees"' not in plotting_stages

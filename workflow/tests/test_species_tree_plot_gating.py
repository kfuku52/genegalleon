import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

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


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="requires the GeneGalleon R runtime")
@pytest.mark.parametrize(("tips", "pages"), [(12, 1), (374, 2)])
def test_species_tree_comparison_pdf_paginates_without_an_extra_render(tmp_path, tips, pages):
    def balanced_tree(labels):
        if len(labels) == 1:
            return labels[0] + ":1"
        mid = len(labels) // 2
        return f"({balanced_tree(labels[:mid])},{balanced_tree(labels[mid:])})95:1"

    tree = tmp_path / "fixture.nwk"
    tree.write_text(balanced_tree([f"Taxon_{i}" for i in range(tips)]) + ";\n")
    script = WORKFLOW_DIR / "support" / "plot_species_trees.r"
    labels = ("iqtree_dna", "iqtree_pep", "astral_dna", "astral_pep")
    result = subprocess.run(
        ["Rscript", str(script), *(f"--{label}_nwk={tree}" for label in labels)],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    pdf = (tmp_path / "species_trees.pdf").read_bytes()
    assert pdf.startswith(b"%PDF-") and pdf.rstrip().endswith(b"%%EOF")
    assert len(re.findall(rb"/Type\s*/Page\b", pdf)) == pages
    assert not (tmp_path / "Rplots.pdf").exists()


def test_plot_only_rejects_changed_recorded_tree(tmp_path):
    core = CORE_PATH.read_text()
    start = core.index("verify_species_tree_plot_inputs() {")
    end = core.index("\ncase ", start)
    validator = core[start:end]
    tree = "output/species_tree"
    for name, step, relative in (
        ("concat_iqtree_pep_root", "species_tree_root_concat_iqtree_pep",
         f"{tree}/concatenated_iqtree_pep/concat.pep.rooted.nwk"),
        ("undated_summary", "species_tree_select_undated_summary",
         f"{tree}/species_tree_summary/undated_species_tree.nwk"),
    ):
        path = tmp_path / relative
        path.parent.mkdir(parents=True)
        body = b"(A:1,B:1);\n"
        path.write_bytes(body)
        provenance = tmp_path / "output/artifact_provenance/genome_evolution"
        provenance.mkdir(parents=True, exist_ok=True)
        (provenance / f"species_tree.{name}.json").write_text(json.dumps({
            "schema_version": 1, "step": step, "family_id": "all_buscos", "inputs": [],
            "outputs": [{"scope": "workspace", "path": relative, "artifact_type": "file",
                         "size_bytes": len(body), "sha256": hashlib.sha256(body).hexdigest()}],
        }))
    script = (
        'gg_support_dir=$1; gg_workspace_dir=$2; input_sequence_mode=protein\n'
        + validator + "\nverify_species_tree_plot_inputs"
    )
    cmd = ["bash", "-c", script, "--", str(WORKFLOW_DIR / "support"), str(tmp_path)]
    valid = subprocess.run(cmd, capture_output=True, text=True)
    assert valid.returncode == 0, valid.stderr
    assert "Verified 2 stored species-tree artifacts" in valid.stdout
    path.write_bytes(b"(A:1,C:1);\n")
    invalid = subprocess.run(cmd, capture_output=True, text=True)
    assert invalid.returncode == 3
    assert "recorded tree changed" in invalid.stderr

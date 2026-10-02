import subprocess
from pathlib import Path

import pytest

WORKFLOW_DIR = Path(__file__).resolve().parents[1]
CORE_PATH = WORKFLOW_DIR / "core" / "gg_gene_evolution_core.sh"
EXECUTION_RUNTIME_PATH = WORKFLOW_DIR / "support" / "gg_util" / "07_execution_runtime.sh"


@pytest.mark.parametrize("rows,accepted", [
    ("HOG0000004\t1\nHOG0000005\t20\n", True),
    ("HOG0000004\t1\nHOG0000004\t20\n", False),
    ("HOG0000004\t1\nHOG0000004\t1\n", False),
    ("\t1\nHOG0000004\t20\n", False),
])
def test_orthogroup_array_catalog_is_unambiguous_before_production(tmp_path, rows, accepted):
    catalog = tmp_path / "Orthogroups.GeneCount.selected.tsv"
    catalog.write_text("Orthogroup\tTotal\n" + rows, encoding="utf-8")
    text = CORE_PATH.read_text(encoding="utf-8")
    start = text.index("    if ! num_orthogroups=$(awk -F")
    end = text.index("    if [[ ${num_orthogroups} -le 0 ]]", start)
    script = 'file_orthogroup_genecount_selected="$1"\n' + text[start:end] + '\nprintf "families=%s\\n" "$num_orthogroups"\n'
    result = subprocess.run(["bash", "-c", script, "catalog-preflight", str(catalog)],
                            text=True, capture_output=True, check=False)
    if accepted:
        assert result.returncode == 0, result.stderr
        assert result.stdout == "families=2\n"
    else:
        assert result.returncode != 0
        assert "Refusing to start gg_gene_evolution with an ambiguous orthogroup catalog" in result.stderr
        assert "families=" not in result.stdout


def _orthogroup_extraction_preflight_source() -> str:
    text = CORE_PATH.read_text(encoding="utf-8")
    start = text.index('task="Orthogroup extraction with NWKIT"')
    end = text.index(
        'if [[ "${mode_gene_evolution}" == "query2family" && ${orthogroup_extraction_needs_update}',
        start,
    )
    return text[start:end]


def _run_orthogroup_extraction_preflight(
    tmp_path: Path,
    *,
    create_rooted_tree: bool,
) -> subprocess.CompletedProcess[str]:
    query = tmp_path / "query.txt"
    query.write_text("Species_one_gene1\n", encoding="utf-8")
    alignment = tmp_path / "trimmed.fa"
    alignment.write_text(">Species_one_gene1\nATG\n", encoding="utf-8")
    unrooted_tree = tmp_path / "unrooted.nwk"
    unrooted_tree.write_text("(Species_one_gene1);\n", encoding="utf-8")
    rooted_tree = tmp_path / "rooted.nwk"
    if create_rooted_tree:
        rooted_tree.write_text("(Species_one_gene1);\n", encoding="utf-8")
    missing_generax_tree = tmp_path / "generax.nhx"
    capture = tmp_path / "provenance-args.txt"

    script = "\n".join(
        [
            "set -eu",
            r'''
execution_runtime_path=$1
file_query_gene=$2
file_og_trimmed_aln_analysis=$3
file_og_iqtree_tree=$4
file_og_rooted_tree=$5
file_og_rooted_tree_analysis=$6
capture_path=$7
source "${execution_runtime_path}"

mode_gene_evolution=query2family
run_orthogroup_extraction=1
dir_output_active=/output/query2family
og_id=OG0001
gg_workspace_dir=/workspace
file_og_query_blast=/output/query2family/query_blast/OG0001.tsv
file_og_orthogroup_extraction_nwk=/output/query2family/orthogroup_extraction_nwk/OG0001.nwk
file_og_orthogroup_extraction_rooted_nwk=/output/query2family/orthogroup_extraction_rooted_nwk/OG0001.nwk
file_og_orthogroup_extraction_fasta=/output/query2family/orthogroup_extraction_fasta/OG0001.fa.gz
input_sequence_mode=cds
genetic_code=1
tree_rooting_method=mv

gg_artifact_prepare_stage() {
  local needs_update_variable=$1
  shift 2
  printf -v "${needs_update_variable}" '%s' 0
  printf '%s\n' "$@" > "${capture_path}"
}
''',
            _orthogroup_extraction_preflight_source(),
            "printf 'run=%s original=%s\\n' \"${run_orthogroup_extraction}\" \"${run_orthogroup_extraction_original}\"",
        ]
    )
    return subprocess.run(
        [
            "bash",
            "-s",
            "--",
            str(EXECUTION_RUNTIME_PATH),
            str(query),
            str(alignment),
            str(unrooted_tree),
            str(rooted_tree),
            str(missing_generax_tree),
            str(capture),
        ],
        input=script,
        text=True,
        capture_output=True,
        check=False,
    )


def test_fresh_run_uses_pre_generax_rooted_tree_for_extraction(tmp_path: Path):
    completed = _run_orthogroup_extraction_preflight(
        tmp_path,
        create_rooted_tree=True,
    )

    assert completed.returncode == 0, completed.stderr + completed.stdout
    assert "run=1 original=1" in completed.stdout
    captured = (tmp_path / "provenance-args.txt").read_text(encoding="utf-8")
    assert f"rooted_tree={tmp_path / 'rooted.nwk'}" in captured
    assert str(tmp_path / "generax.nhx") not in captured


def test_requested_extraction_fails_instead_of_silently_skipping_missing_root(
    tmp_path: Path,
):
    completed = _run_orthogroup_extraction_preflight(
        tmp_path,
        create_rooted_tree=False,
    )

    assert completed.returncode != 0
    assert "Disabled run_orthogroup_extraction" in completed.stdout
    assert (
        "Orthogroup extraction was requested but its required inputs are unavailable. "
        "Refusing to continue without extraction."
    ) in completed.stderr
    assert not (tmp_path / "provenance-args.txt").exists()

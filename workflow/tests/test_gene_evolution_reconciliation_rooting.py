import subprocess
from pathlib import Path

WORKFLOW_DIR = Path(__file__).resolve().parents[1]
CORE_PATH = WORKFLOW_DIR / "core" / "gg_gene_evolution_core.sh"


def _tree_rooting_block() -> str:
    text = CORE_PATH.read_text(encoding="utf-8")
    start = text.index('task="Gene tree rooting"')
    end = text.index('task="Orthogroup extraction with NWKIT"', start)
    return text[start:end]


def _run_reconciliation_rooting(
    tmp_path: Path,
    *,
    create_species_tree: bool,
) -> subprocess.CompletedProcess[str]:
    unrooted_tree = tmp_path / "unrooted.nwk"
    unrooted_tree.write_text("(A_g,B_g);\n", encoding="utf-8")
    species_tree = tmp_path / "species.nwk"
    if create_species_tree:
        species_tree.write_text("(A,B);\n", encoding="utf-8")
    species_map = tmp_path / "species-map.tsv"
    species_map.write_text("gene\tspecies\nA_g\tA\nB_g\tB\n", encoding="utf-8")
    rooted_tree = tmp_path / "rooted.nwk"
    rooting_log = tmp_path / "rooting.log"
    call_capture = tmp_path / "nwkit-calls.txt"
    provenance_capture = tmp_path / "provenance-args.txt"

    harness = r'''
set -euo pipefail
tmp_dir=$1
file_og_unrooted_tree_analysis=$2
species_tree_pruned=$3
species_label_map_tsv=$4
file_og_rooted_tree=$5
file_og_rooted_log=$6
nwkit_call_capture=$7
provenance_capture=$8
cd "${tmp_dir}"

run_tree_root=1
tree_rooting_method=reconciliation
species_label_parser=taxonomic
species_label_regex='^([^_]+)'
dir_output_active="${tmp_dir}/output"
og_id=OG0001
gg_workspace_dir="${tmp_dir}/workspace"
memory_notung=1
notung_jar=/unused/notung.jar
mkdir -p "${dir_output_active}/artifact_provenance"

disable_if_no_input_file() {
  local variable_name=$1
  shift
  local input_path
  for input_path in "$@"; do
    if [[ ! -s "${input_path}" ]]; then
      printf -v "${variable_name}" '%s' 0
    fi
  done
}

gg_artifact_prepare_stage() {
  local needs_update_variable=$1
  shift 2
  printf -v "${needs_update_variable}" '%s' 1
  printf '%s\n' "$@" > "${provenance_capture}"
}

gg_step_start() { :; }
gg_step_skip() { :; }
gg_artifact_record() { :; }
mv_out() { mv -- "$1" "$2"; }

nwkit() {
  if [[ "$1" == "root" ]]; then
    printf '%s\n' "$@" > "${nwkit_call_capture}"
    printf '(A_g,B_g);\n'
    return
  fi
  if [[ "$1" == "label" ]]; then
    local output_path=""
    local previous=""
    local argument
    for argument in "$@"; do
      if [[ "${previous}" == "--outfile" ]]; then
        output_path=${argument}
        break
      fi
      previous=${argument}
    done
    cat > "${output_path}"
    return
  fi
  return 1
}
'''
    script = "\n".join((harness, _tree_rooting_block()))
    return subprocess.run(
        [
            "bash",
            "-s",
            "--",
            str(tmp_path),
            str(unrooted_tree),
            str(species_tree),
            str(species_map),
            str(rooted_tree),
            str(rooting_log),
            str(call_capture),
            str(provenance_capture),
        ],
        input=script,
        text=True,
        capture_output=True,
        check=False,
    )


def test_reconciliation_rooting_passes_species_inputs_to_nwkit(tmp_path: Path):
    completed = _run_reconciliation_rooting(
        tmp_path,
        create_species_tree=True,
    )

    assert completed.returncode == 0, completed.stderr + completed.stdout
    root_call = (tmp_path / "nwkit-calls.txt").read_text(
        encoding="utf-8"
    ).splitlines()
    assert root_call == [
        "root",
        "--method",
        "reconciliation",
        "--infile",
        str(tmp_path / "unrooted.nwk"),
        "--species-tree",
        str(tmp_path / "species.nwk"),
        "--species-parser",
        "taxonomic",
        "--species-regex",
        "^([^_]+)",
        "--species-map-tsv",
        str(tmp_path / "species-map.tsv"),
    ]
    assert (tmp_path / "rooted.nwk").read_text(encoding="utf-8") == "(A_g,B_g);\n"
    assert "tree_rooting_method=reconciliation" in (
        tmp_path / "rooting.log"
    ).read_text(encoding="utf-8")
    provenance = (tmp_path / "provenance-args.txt").read_text(encoding="utf-8")
    assert f"species_tree={tmp_path / 'species.nwk'}" in provenance
    assert f"species_map={tmp_path / 'species-map.tsv'}" in provenance


def test_reconciliation_rooting_requires_species_tree(tmp_path: Path):
    completed = _run_reconciliation_rooting(
        tmp_path,
        create_species_tree=False,
    )

    assert completed.returncode != 0
    assert "tree_rooting_method=reconciliation requires species tree" in completed.stdout
    assert not (tmp_path / "nwkit-calls.txt").exists()

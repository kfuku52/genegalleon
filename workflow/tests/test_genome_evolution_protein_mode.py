from pathlib import Path
import os
import re
import shlex
import stat
import subprocess

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_DIR = REPO_ROOT / "workflow"
CORE_PATH = WORKFLOW_DIR / "core" / "gg_genome_evolution_core.sh"
ENTRYPOINT_PATH = WORKFLOW_DIR / "gg_genome_evolution_entrypoint.sh"
SYSTEM_BASH_MAJOR = int(
    subprocess.run(
        ["/bin/bash", "-c", "printf '%s' \"${BASH_VERSINFO[0]}\""],
        capture_output=True,
        text=True,
        check=False,
    ).stdout
    or "0"
)


def _write_executable(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _prepare_stub_binaries(tmp_path: Path) -> Path:
    bin_dir = tmp_path / "bin"
    capture_dir = tmp_path / "capture"
    capture_dir.mkdir(parents=True, exist_ok=True)

    _write_executable(
        bin_dir / "parallel",
        """#!/usr/bin/env bash
set -euo pipefail
func=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --jobs)
      shift 2
      ;;
    :::)
      shift
      break
      ;;
    *)
      if [[ -z "${func}" ]]; then
        func="$1"
      fi
      shift
      ;;
  esac
done
for item in "$@"; do
  bash -lc "${func} \"${1}\"" _ "${item}"
done
""",
    )

    _write_executable(
        bin_dir / "conda",
        """#!/usr/bin/env bash
set -euo pipefail
if [[ "${1:-}" == "shell.bash" && "${2:-}" == "hook" ]]; then
  cat <<'EOF'
conda() {
  return 0
}
EOF
  exit 0
fi
exit 0
""",
    )

    _write_executable(
        bin_dir / "ulimit",
        """#!/usr/bin/env bash
exit 0
""",
    )

    _write_executable(
        bin_dir / "orthofinder",
        f"""#!/usr/bin/env bash
set -euo pipefail
capture_dir={shlex.quote(str(capture_dir))}
input_dir=""
output_dir=""
run_name=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -f)
      input_dir="$2"
      shift 2
      ;;
    -o)
      output_dir="$2"
      shift 2
      ;;
    -n)
      run_name="$2"
      shift 2
      ;;
    *)
      shift
      ;;
  esac
done
results_dir="${{output_dir}}/Results_${{run_name}}"
mkdir -p "${{results_dir}}/Phylogenetic_Hierarchical_Orthogroups"
find "${{input_dir}}" -maxdepth 1 -type f ! -name '.*' | sort > "${{capture_dir}}/input_files.txt"
>"${{capture_dir}}/proteins.fasta"
while IFS= read -r fasta; do
  cat "${{fasta}}" >> "${{capture_dir}}/proteins.fasta"
done < "${{capture_dir}}/input_files.txt"
python - "${{capture_dir}}/input_files.txt" "${{results_dir}}/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" <<'PY'
import pathlib
import sys

input_files = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8").splitlines()
outfile = pathlib.Path(sys.argv[2])
species = []
for path in input_files:
    name = pathlib.Path(path).name
    if "_" in name:
        species.append(name.split("_", 2)[0] + "_" + name.split("_", 2)[1])
    else:
        species.append(pathlib.Path(name).stem)
with outfile.open("w", encoding="utf-8", newline="") as handle:
    cols = ["HOG", "OG", "Gene Tree Parent Clade", *species]
    handle.write("\\t".join(cols) + "\\n")
    genes = [f"{{sp}}_gene1" for sp in species]
    handle.write("\\t".join(["N0.HOG0000001", "OG0000001", "n0", *genes]) + "\\n")
PY
""",
    )

    return bin_dir


def _load_entrypoint_defaults() -> dict[str, str]:
    text = ENTRYPOINT_PATH.read_text(encoding="utf-8")
    in_block = False
    defaults: dict[str, str] = {}
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if line.startswith("### Start: Modify this block"):
            in_block = True
            continue
        if line.startswith("### End: Modify this block"):
            break
        if not in_block or not line or line.startswith("#"):
            continue
        match = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)=(.*?)(?:\s+#.*)?$", line)
        if not match:
            continue
        key = match.group(1)
        value = match.group(2).strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
            value = value[1:-1]
        defaults[key] = value
    return defaults


def _run_core(tmp_path: Path) -> subprocess.CompletedProcess[str]:
    workspace = tmp_path / "workspace"
    workspace.mkdir(parents=True, exist_ok=True)
    bin_dir = _prepare_stub_binaries(tmp_path)
    env = {
        key: value
        for key, value in os.environ.items()
        if not key.startswith("CONDA") and not key.startswith("_CE_")
    }
    env["PATH"] = f"{bin_dir}{os.pathsep}{env['PATH']}"
    env["gg_workspace_dir"] = str(workspace)
    env["GG_TASK_CPUS"] = "1"
    env["GG_MEM_PER_CPU_GB"] = "1"
    env["GG_MEM_TOTAL_GB"] = "1"
    env["GG_JOB_ID"] = "1"
    env["GG_ARRAY_TASK_ID"] = "1"
    env.update(_load_entrypoint_defaults())
    env.update(
        {
            "input_sequence_mode": "protein",
            "run_species_busco": "0",
            "run_species_get_busco_summary": "0",
            "run_individual_get_fasta": "0",
            "run_individual_mafft": "0",
            "run_individual_trimal": "0",
            "run_concat_alignment": "0",
            "run_concat_iqtree_protein": "0",
            "run_concat_iqtree_dna": "0",
            "run_individual_iqtree_pep": "0",
            "run_astral_pep": "0",
            "run_individual_iqtree_dna": "0",
            "run_astral_dna": "0",
            "run_plot_species_trees": "0",
            "run_constrained_tree": "0",
            "run_plot_constrained_tree": "0",
            "run_mcmctree1": "0",
            "run_mcmctree2": "0",
            "run_convert_tree_format": "0",
            "run_plot_mcmctreer": "0",
            "run_orthofinder": "1",
            "run_og_selection": "0",
            "run_orthogroup_method_comparison": "0",
            "run_genome_busco": "0",
            "run_genome_get_busco_summary": "0",
            "run_busco_getfasta": "0",
            "run_busco_mafft": "0",
            "run_busco_trimal": "0",
            "run_busco_iqtree_dna": "0",
            "run_busco_iqtree_pep": "0",
            "run_busco_notung_root_dna": "0",
            "run_busco_notung_root_pep": "0",
            "run_busco_root_dna": "0",
            "run_busco_root_pep": "0",
            "run_busco_grampa_dna": "0",
            "run_busco_grampa_pep": "0",
            "run_orthogroup_grampa": "0",
            "run_cafe": "0",
            "run_go_enrichment": "0",
            "change_direction_go": "increase",
            "target_branch_go": "",
        }
    )

    command = (
        "set -euo pipefail; "
        "mapfile() { "
        "  if [[ ${1:-} == '-t' ]]; then shift; fi; "
        "  local target_array=${1:?}; "
        "  local line=''; "
        "  eval \"${target_array}=()\"; "
        "  while IFS= read -r line; do "
        "    eval \"${target_array}+=(\\\"\\${line}\\\")\"; "
        "  done; "
        "}; "
        "enable -n ulimit; "
        f"source {shlex.quote(str(CORE_PATH))}"
    )
    return subprocess.run(
        ["bash", "-c", command],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )


@pytest.mark.skipif(SYSTEM_BASH_MAJOR < 4, reason="gg_genome_evolution_core.sh requires bash 4+ features such as local -n")
def test_genome_evolution_protein_mode_prefers_species_protein_inputs(tmp_path: Path):
    workspace = tmp_path / "workspace"
    species_protein_dir = workspace / "input" / "species_protein"
    species_cds_dir = workspace / "input" / "species_cds"
    species_code_dir = workspace / "input" / "species_genetic_code"
    species_protein_dir.mkdir(parents=True)
    species_cds_dir.mkdir(parents=True)
    species_code_dir.mkdir(parents=True)

    (species_protein_dir / "Tetrahymena_thermophila_pep.fa").write_text(
        ">Tetrahymena_thermophila_gene1\nMPEP\n",
        encoding="utf-8",
    )
    (species_cds_dir / "Tetrahymena_thermophila_cds.fa").write_text(
        ">Tetrahymena_thermophila_gene1\nATGTAA\n",
        encoding="utf-8",
    )
    (species_code_dir / "species_genetic_code.tsv").write_text(
        "species\tgenetic_code\nTetrahymena_thermophila\t6\n",
        encoding="utf-8",
    )

    completed = _run_core(tmp_path)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "species_genetic_code.tsv is ignored because species_protein inputs are provided" in completed.stdout
    proteins = (tmp_path / "capture" / "proteins.fasta").read_text(encoding="utf-8")
    assert ">Tetrahymena_thermophila_gene1" in proteins
    assert "MPEP" in proteins
    assert "MQ" not in proteins
    assert "M*" not in proteins


@pytest.mark.skipif(SYSTEM_BASH_MAJOR < 4, reason="gg_genome_evolution_core.sh requires bash 4+ features such as local -n")
def test_genome_evolution_protein_mode_translates_species_cds_with_species_specific_codes(tmp_path: Path):
    workspace = tmp_path / "workspace"
    species_cds_dir = workspace / "input" / "species_cds"
    species_code_dir = workspace / "input" / "species_genetic_code"
    species_cds_dir.mkdir(parents=True)
    species_code_dir.mkdir(parents=True)

    (species_cds_dir / "Tetrahymena_thermophila_cds.fa").write_text(
        ">Tetrahymena_thermophila_gene1\nATGTAA\n",
        encoding="utf-8",
    )
    (species_cds_dir / "Arabidopsis_thaliana_cds.fa").write_text(
        ">Arabidopsis_thaliana_gene1\nATGTAA\n",
        encoding="utf-8",
    )
    (species_code_dir / "species_genetic_code.tsv").write_text(
        "species\tgenetic_code\nTetrahymena_thermophila\t6\n",
        encoding="utf-8",
    )

    completed = _run_core(tmp_path)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    proteins = (tmp_path / "capture" / "proteins.fasta").read_text(encoding="utf-8")
    assert ">Tetrahymena_thermophila_gene1" in proteins
    assert ">Arabidopsis_thaliana_gene1" in proteins
    assert "MQ" in proteins
    assert "M*" in proteins

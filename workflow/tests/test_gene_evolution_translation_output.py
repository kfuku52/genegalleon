import re
import subprocess
from pathlib import Path

import pytest


CORE_PATH = Path(__file__).resolve().parents[1] / "core" / "gg_gene_evolution_core.sh"
BASH_HAS_MAPFILE = subprocess.run(
    ["bash", "-c", "type mapfile >/dev/null 2>&1"],
    check=False,
).returncode == 0
pytestmark = pytest.mark.skipif(
    not BASH_HAS_MAPFILE,
    reason="the workflow requires Bash 4 mapfile; macOS Bash 3 cannot exercise this helper",
)


def _shell_function_source(name: str) -> str:
    text = CORE_PATH.read_text(encoding="utf-8")
    match = re.search(
        rf"^{re.escape(name)}\(\) \{{\n.*?^\}}$",
        text,
        re.DOTALL | re.MULTILINE,
    )
    assert match is not None
    return match.group(0)


def _run_translation(cds_fasta: Path, protein_out: Path, *, fail_final_write: bool):
    script = "\n".join(
        [
            "set -u",
            _shell_function_source("translate_orthogroup_cds_to_protein_fasta"),
            r'''
cds_fasta=$1
protein_out=$2
fail_final_write=$3
og_id=OG0001
GG_TASK_CPUS=1
genetic_code=1

ensure_parent_dir() {
  mkdir -p -- "$(dirname -- "$1")"
}
mv_out() {
  ensure_parent_dir "$2"
  mv -- "$1" "$2"
}
gg_species_name_from_path_or_dot() {
  printf '%s\n' sp1
}
gg_lookup_species_genetic_code() {
  printf '%s\n' 1
}
gg_prepare_cds_fasta_stream() {
  cat
}
gg_count_fasta_records() {
  grep -c '^>' "$1"
}
seqkit() {
  local command=$1
  shift
  if [[ "${command}" == "grep" ]]; then
    printf '>sp1_gene1\nATG\n'
    return 0
  fi
  if [[ "${command}" == "translate" ]]; then
    cat
    return 0
  fi
  if [[ "${command}" != "seq" ]]; then
    return 2
  fi
  if [[ " $* " == *" --name "* ]]; then
    printf '>sp1_gene1\n'
    return 0
  fi
  if [[ " $* " == *" --remove-gaps "* ]]; then
    cat
    return 0
  fi
  local source_path=""
  local destination_path=""
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --threads)
        shift 2
        ;;
      --out-file)
        destination_path=$2
        shift 2
        ;;
      *)
        source_path=$1
        shift
        ;;
    esac
  done
  if [[ "${fail_final_write}" == "yes" ]]; then
    printf 'partial\n' > "${destination_path}"
    return 9
  fi
  cp -- "${source_path}" "${destination_path}"
}

translate_orthogroup_cds_to_protein_fasta "${cds_fasta}" "${protein_out}" ignored.tsv
''',
        ]
    )
    return subprocess.run(
        [
            "bash",
            "-s",
            "--",
            str(cds_fasta),
            str(protein_out),
            "yes" if fail_final_write else "no",
        ],
        input=script,
        text=True,
        capture_output=True,
        check=False,
        cwd=cds_fasta.parents[1],
    )


def test_translation_creates_missing_output_parent_and_atomically_publishes(tmp_path):
    cds_fasta = tmp_path / "input" / "OG0001_cds.fa"
    cds_fasta.parent.mkdir()
    cds_fasta.write_text(">sp1_gene1\nATG\n", encoding="utf-8")
    protein_out = tmp_path / "output" / "protein_fasta" / "OG0001_pep.fa"

    completed = _run_translation(cds_fasta, protein_out, fail_final_write=False)

    assert completed.returncode == 0, completed.stderr + completed.stdout
    assert protein_out.read_text(encoding="utf-8") == ">sp1_gene1\nATG\n"
    assert not list(protein_out.parent.glob("OG0001_pep.fa.partial.*"))


def test_translation_removes_partial_output_when_seqkit_write_fails(tmp_path):
    cds_fasta = tmp_path / "input" / "OG0001_cds.fa"
    cds_fasta.parent.mkdir()
    cds_fasta.write_text(">sp1_gene1\nATG\n", encoding="utf-8")
    protein_out = tmp_path / "output" / "protein_fasta" / "OG0001_pep.fa"

    completed = _run_translation(cds_fasta, protein_out, fail_final_write=True)

    assert completed.returncode != 0
    assert not protein_out.exists()
    assert not list(protein_out.parent.glob("OG0001_pep.fa.partial.*"))
    assert not (tmp_path / "OG0001.translated.pep.tmp.fasta").exists()
    assert "Failed to write translated protein FASTA" in completed.stderr

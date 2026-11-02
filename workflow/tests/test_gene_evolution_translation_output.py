import gzip
import os
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


@pytest.mark.parametrize("stage", ["jcvi", "genome_translation"])
def test_downstream_fasta_prefix_processors_match_literal_species(tmp_path, stage):
    species = "Plant_cf._alpha"
    near_match = "Plant_cfX_alpha_g2"
    payload = f">{species}_g1\nATG\n>{near_match}\nAAA\n>other_gene3\nCCC\n"
    if stage == "jcvi":
        source = (CORE_PATH.parent / "gg_genome_annotation_core.sh").read_text()
        processor = source.split('seqkit seq "${file_sp_cds}" |', 1)[1].split('> tmp.jcvi_input.cds.fasta', 1)[0]
        expected = ["g1", near_match, "other_gene3"]
    else:
        source = (CORE_PATH.parent / "gg_genome_evolution_core.sh").read_text()
        processor = source.split('seqkit translate --transl-table "${species_code}" --threads "${GG_TASK_CPUS}" |', 1)[1]
        processor = processor.split("sed -e '/^1 1$/d'", 1)[0].rstrip().removesuffix("|")
        expected = [f"{species}_g1", f"{species}_{near_match}", f"{species}_other_gene3"]
    result = subprocess.run(["bash", "-c", processor], input=payload, text=True, capture_output=True,
                            cwd=tmp_path, env=dict(os.environ, sp_ub=species), check=True)
    assert [line[1:] for line in result.stdout.splitlines() if line.startswith(">")] == expected
    assert [line for line in result.stdout.splitlines() if not line.startswith(">")] == ["ATG", "AAA", "CCC"]


@pytest.mark.parametrize("known_species", ["", "Standard_species", "Ciliate_species"])
def test_native_preparation_uses_each_species_code_and_preserves_internal_stops(tmp_path, known_species):
    cds = tmp_path / "mixed.fa"
    cds.write_text("".join(f">{species}_gene1\nATGTAATAGTGATGG\n" for species in
                           ([known_species] if known_species else ["Standard_species", "Ciliate_species"])))
    codes = tmp_path / "codes.tsv"
    codes.write_text("species\tgenetic_code\nStandard_species\t1\nCiliate_species\t6\n")
    protein = tmp_path / "mixed.pep.fa.gz"
    env = os.environ.copy()
    env.update(GG_CORE_SOURCE_ONLY="1", GG_TASK_CPUS="1", GG_TASK_MEM_GB="2",
               gg_workspace_dir=str(tmp_path / "workspace"), gg_support_dir=str(CORE_PATH.parent.parent / "support"))
    result = subprocess.run([
        "bash", "-c",
        'source "$1"; og_id=mixed; translate_orthogroup_cds_to_protein_fasta "$2" "$3" "$4" "$5"',
        "_", str(CORE_PATH), str(cds), str(protein), str(codes), known_species,
    ], cwd=tmp_path, env=env, text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr
    with gzip.open(protein, "rt") as handle:
        observed = handle.read()
    if known_species != "Ciliate_species":
        assert ">Standard_species_gene1\nM***W\n" in observed
    if known_species != "Standard_species":
        assert ">Ciliate_species_gene1\nMQQ*W\n" in observed
    assert not (tmp_path / "workspace" / "output" / "query2family").exists()


def test_declared_species_cannot_silently_drop_records(tmp_path):
    cds = tmp_path / "mixed.fa"
    cds.write_text(">Standard_species_gene1\nATGTAA\n>Ciliate_species_gene2\nATGTAA\n")
    table = tmp_path / "codes.tsv"
    table.write_text("species\tgenetic_code\nStandard_species\t1\nCiliate_species\t6\n")
    protein = tmp_path / "out.fa.gz"
    env = os.environ.copy()
    env.update(GG_CORE_SOURCE_ONLY="1", GG_TASK_CPUS="1", GG_MEM_TOTAL_GB="2",
               gg_workspace_dir=str(tmp_path / "workspace"), gg_support_dir=str(CORE_PATH.parent.parent / "support"))
    result = subprocess.run([
        "bash", "-c", 'source "$1"; og_id=check; translate_orthogroup_cds_to_protein_fasta "$2" "$3" "$4" Standard_species',
        "_", str(CORE_PATH), str(cds), str(protein), str(table),
    ], cwd=tmp_path, env=env, text=True, capture_output=True)
    assert result.returncode != 0
    assert "different number of sequences" in result.stdout
    assert not protein.exists()


@pytest.mark.parametrize("header,accepted", [("Standard_species.v1_gene1", True),
                                           ("Standard_speciesXv1_gene1", False)])
def test_declared_species_regex_metacharacters_are_literal(tmp_path, header, accepted):
    cds, table, protein = tmp_path / "input.fa", tmp_path / "codes.tsv", tmp_path / "protein.fa.gz"
    cds.write_text(">" + header + "\nATGAAATAA\n")
    table.write_text("species\tgenetic_code\nStandard_species.v1\t1\n")
    env = dict(os.environ, GG_CORE_SOURCE_ONLY="1", GG_TASK_CPUS="1", GG_TASK_MEM_GB="2",
               gg_workspace_dir=str(tmp_path / "workspace"), gg_support_dir=str(CORE_PATH.parent.parent / "support"))
    result = subprocess.run(["bash", "-c", 'source "$1"; og_id=literal; '
                             'translate_orthogroup_cds_to_protein_fasta "$2" "$3" "$4" Standard_species.v1',
                             "_", str(CORE_PATH), str(cds), str(protein), str(table)],
                            env=env, cwd=tmp_path, text=True, capture_output=True)
    assert (result.returncode == 0) is accepted, result.stdout + result.stderr
    assert protein.exists() is accepted


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
  if [[ "${destination_path}" == *.gz ]]; then
    gzip -c -- "${source_path}" > "${destination_path}"
  else
    cp -- "${source_path}" "${destination_path}"
  fi
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


def test_translation_keeps_gzip_suffix_on_staging_file(tmp_path):
    cds_fasta = tmp_path / "input" / "OG0001_cds.fa"
    cds_fasta.parent.mkdir()
    cds_fasta.write_text(">sp1_gene1\nATG\n", encoding="utf-8")
    protein_out = tmp_path / "output" / "protein_fasta" / "OG0001_pep.fa.gz"

    completed = _run_translation(cds_fasta, protein_out, fail_final_write=False)

    assert completed.returncode == 0, completed.stderr + completed.stdout
    with gzip.open(protein_out, "rt", encoding="utf-8") as handle:
        assert handle.read() == ">sp1_gene1\nATG\n"
    assert not list(protein_out.parent.glob("OG0001_pep.fa.gz.partial.*"))

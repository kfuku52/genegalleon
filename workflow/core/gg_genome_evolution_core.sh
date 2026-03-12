#!/usr/bin/env bash
set -euo pipefail

# Unified genome evolution pipeline (single-file form).
# Former inline stages are now written directly in this file.

gg_core_bootstrap="/script/support/gg_core_bootstrap.sh"
if [[ ! -s "${gg_core_bootstrap}" ]]; then
  gg_core_bootstrap="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../support/gg_core_bootstrap.sh"
fi
# shellcheck disable=SC1090
source "${gg_core_bootstrap}"
unset gg_core_bootstrap
gg_source_common_params_from_core "${BASH_SOURCE[0]:-$0}"

#run a busco analysis for all species in species_cds and extract shared complete BUSCO genes among all species.
#then make cds fasta files including all species' cds sequenes for each BUSCO genes and make alignment
#then concatenate all fasta files for further downstream analysis -> task = trimal in gg_species_tree.

### Start: Job-supplied configuration ###

# Configuration variables are provided by gg_genome_evolution_entrypoint.sh.
genetic_code="${genetic_code:-${GG_COMMON_GENETIC_CODE:-1}}"
input_sequence_mode="${input_sequence_mode:-cds}"
busco_lineage="${busco_lineage:-${GG_COMMON_BUSCO_LINEAGE:-auto}}"
species_tree_rooting="${species_tree_rooting:-taxonomy}"
annotation_species="${annotation_species:-${GG_COMMON_REFERENCE_SPECIES:-auto}}"
omark_db_path="${omark_db_path:-auto}"
# Backward-compatible defaults for launchers that predate the shared protein/OMArk flags.
run_cds_translation="${run_cds_translation:-1}"
run_species_omark="${run_species_omark:-0}"
run_species_get_omark_summary="${run_species_get_omark_summary:-1}"
mcmctree_divergence_time_constraints_str="${mcmctree_divergence_time_constraints_str:-}"
grampa_h1="${grampa_h1:-}"
target_branch_go="${target_branch_go:-}"
mcmctree_divergence_time_constraints=()
if [[ -n "${mcmctree_divergence_time_constraints_str:-}" ]]; then
  IFS='|' read -r -a mcmctree_divergence_time_constraints <<< "${mcmctree_divergence_time_constraints_str}"
fi

### End: Job-supplied configuration ###

### ----------------------------------------------------------------------- ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

gg_bootstrap_core_runtime "${BASH_SOURCE[0]:-$0}" "base" 1 1
delete_tmp_dir=${delete_tmp_dir:-1}
busco_lineage_resolved=""
omark_db_resolved=""

copy_busco_tables() {
  local busco_root_dir="$1"
  local lineage="$2"
  local file_full="$3"
  local file_short="$4"
  local run_dir="${busco_root_dir}/run_${lineage}"
  local full_src=""
  local short_src=""
  local -a short_candidates=()

  if [[ -s "${run_dir}/full_table.tsv" ]]; then
    full_src="${run_dir}/full_table.tsv"
  elif [[ -s "${run_dir}/full_table.tsv.gz" ]]; then
    gzip -cd "${run_dir}/full_table.tsv.gz" > "./tmp.busco.full_table.tsv"
    full_src="./tmp.busco.full_table.tsv"
  fi

  if [[ -s "${run_dir}/short_summary.txt" ]]; then
    short_src="${run_dir}/short_summary.txt"
  else
    mapfile -t short_candidates < <(find "${run_dir}" -maxdepth 1 -type f -name "short_summary*.txt" | sort)
    if [[ ${#short_candidates[@]} -gt 0 ]]; then
      short_src="${short_candidates[0]}"
    fi
  fi

  if [[ -z "${full_src}" || -z "${short_src}" ]]; then
    echo "BUSCO outputs were not found under ${run_dir}."
    echo "Expected full_table.tsv(.gz) and short_summary*.txt."
    rm -f -- "./tmp.busco.full_table.tsv"
    return 1
  fi

  cp_out "${full_src}" "${file_full}"
  cp_out "${short_src}" "${file_short}"
  rm -f -- "./tmp.busco.full_table.tsv"
}

resolve_busco_lineage_for_species_set() {
  if [[ -n "${busco_lineage_resolved}" ]]; then
    return 0
  fi
  if [[ $# -eq 0 ]]; then
    echo "No species names were provided for BUSCO lineage resolution." >&2
    return 1
  fi
  if ! busco_lineage_resolved=$(gg_resolve_busco_lineage "${gg_workspace_dir}" "${busco_lineage}" "$@"); then
    echo "Failed to resolve BUSCO lineage from request: ${busco_lineage}" >&2
    return 1
  fi
  echo "Resolved BUSCO lineage for species set (${#} species): ${busco_lineage_resolved}"
}

get_busco_summary_gene_count() {
  local summary_table=$1
  local line_count=0
  if [[ ! -s "${summary_table}" ]]; then
    echo 0
    return 0
  fi
  line_count=$(wc -l < "${summary_table}")
  if [[ ${line_count} -le 1 ]]; then
    echo 0
  else
    echo $((line_count - 1))
  fi
}

normalize_busco_table_naming() {
  local full_dir=$1
  local short_dir=$2
  # BUSCO outputs are produced directly in canonical naming:
  #   <species>.busco.full.tsv / <species>.busco.short.txt
  # Keep this function as a stable call point for summary steps.
  if [[ -n "${full_dir}" ]]; then
    ensure_dir "${full_dir}"
  fi
  if [[ -n "${short_dir}" ]]; then
    ensure_dir "${short_dir}"
  fi
}

busco_output_exists_for_species() {
  local search_dir=$1
  local species_name=$2
  local name_glob=$3
  local busco_file busco_base busco_species

  if [[ -z "${search_dir}" || ! -d "${search_dir}" ]]; then
    return 1
  fi

  while IFS= read -r busco_file; do
    busco_base=$(basename "${busco_file}")
    busco_species=$(gg_species_name_from_path_or_dot "${busco_base}")
    if [[ "${busco_species}" == "${species_name}" ]]; then
      return 0
    fi
  done < <(find "${search_dir}" -maxdepth 1 -type f ! -name '.*' -name "${name_glob}" 2> /dev/null | sort)

  return 1
}

remove_busco_outputs_for_species() {
  local search_dir=$1
  local species_name=$2
  local name_glob=$3
  local busco_file busco_base busco_species

  if [[ -z "${search_dir}" || ! -d "${search_dir}" ]]; then
    return 0
  fi

  while IFS= read -r busco_file; do
    busco_base=$(basename "${busco_file}")
    busco_species=$(gg_species_name_from_path_or_dot "${busco_base}")
    if [[ "${busco_species}" == "${species_name}" ]]; then
      rm -f -- "${busco_file}"
    fi
  done < <(find "${search_dir}" -maxdepth 1 -type f ! -name '.*' -name "${name_glob}" 2> /dev/null | sort)
}

trim_ascii_whitespace() {
  local s=$1
  s="${s#"${s%%[![:space:]]*}"}"
  s="${s%"${s##*[![:space:]]}"}"
  printf '%s' "${s}"
}

species_tree_rooting_looks_like_label_list() {
  local token=""
  if [[ $# -eq 0 ]]; then
    return 1
  fi
  for token in "$@"; do
    if [[ -z "${token}" ]]; then
      continue
    fi
    if [[ ! "${token}" =~ ^[[:alnum:].-]+_[[:alnum:]_.-]+$ ]]; then
      return 1
    fi
  done
  return 0
}

input_sequence_mode=$(printf '%s' "${input_sequence_mode}" | tr '[:upper:]' '[:lower:]')
if [[ "${input_sequence_mode}" != "cds" && "${input_sequence_mode}" != "protein" ]]; then
  echo "Invalid input_sequence_mode: ${input_sequence_mode}"
  echo 'input_sequence_mode must be either "cds" or "protein".'
  exit 1
fi

species_protein_ready=0
species_protein_source=""
species_tree_sequence_label="cds"
species_tree_busco_mode="transcriptome"
single_copy_fasta_glob="*.cds.fa.gz"
single_copy_aln_glob="*.cds.aln.fa.gz"
single_copy_trimal_glob="*.trimal.fa.gz"
single_copy_fasta_suffix=".cds.fa.gz"
single_copy_aln_suffix=".cds.aln.fa.gz"
single_copy_trimal_suffix=".trimal.fa.gz"
genome_busco_fasta_glob="*.busco.cds.fa.gz"
genome_busco_aln_glob="*.busco.cds.aln.fa.gz"
genome_busco_trimal_glob="*.busco.trimal.fa.gz"
genome_busco_fasta_suffix=".busco.cds.fa.gz"
genome_busco_aln_suffix=".busco.cds.aln.fa.gz"
genome_busco_trimal_suffix=".busco.trimal.fa.gz"
if [[ "${input_sequence_mode}" == "protein" ]]; then
  species_tree_sequence_label="pep"
  species_tree_busco_mode="proteins"
  single_copy_fasta_glob="*.pep.fa.gz"
  single_copy_aln_glob="*.pep.aln.fa.gz"
  single_copy_trimal_glob="*.pep.trimal.fa.gz"
  single_copy_fasta_suffix=".pep.fa.gz"
  single_copy_aln_suffix=".pep.aln.fa.gz"
  single_copy_trimal_suffix=".pep.trimal.fa.gz"
  genome_busco_fasta_glob="*.busco.pep.fa.gz"
  genome_busco_aln_glob="*.busco.pep.aln.fa.gz"
  genome_busco_trimal_glob="*.busco.pep.trimal.fa.gz"
  genome_busco_fasta_suffix=".busco.pep.fa.gz"
  genome_busco_aln_suffix=".busco.pep.aln.fa.gz"
  genome_busco_trimal_suffix=".busco.pep.trimal.fa.gz"
fi

species_genetic_code_table_path() {
  echo "${gg_workspace_input_dir}/species_genetic_code/species_genetic_code.tsv"
}

species_protein_input_dir_path() {
  echo "${gg_workspace_input_dir}/species_protein"
}

prepare_species_genetic_code_table() {
  local cds_dir=$1
  local default_code=$2
  local outfile=$3
  local input_table
  input_table=$(species_genetic_code_table_path)
  python - "${cds_dir}" "${default_code}" "${outfile}" "${input_table}" <<'PY'
import csv
import gzip
import os
import sys
from pathlib import Path


def species_from_filename(path: Path) -> str:
    name = path.name
    if "_" not in name:
        return path.stem
    return name.split("_", 2)[0] + "_" + name.split("_", 2)[1]


cds_dir = Path(sys.argv[1])
default_code = sys.argv[2].strip()
outfile = Path(sys.argv[3])
input_table = Path(sys.argv[4])
fasta_suffixes = (".fa", ".fa.gz", ".fas", ".fas.gz", ".fasta", ".fasta.gz", ".fna", ".fna.gz")

species = sorted(
    species_from_filename(path)
    for path in cds_dir.iterdir()
    if path.is_file() and not path.name.startswith(".") and path.name.endswith(fasta_suffixes)
)
if not species:
    sys.stderr.write(f"No species CDS files were found in: {cds_dir}\n")
    sys.exit(1)

overrides = {}
if input_table.exists():
    with input_table.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(
            (line for line in handle if line.strip() and not line.lstrip().startswith("#")),
            delimiter="\t",
        )
        if reader.fieldnames is None or "species" not in reader.fieldnames or "genetic_code" not in reader.fieldnames:
            sys.stderr.write(
                f"{input_table} must be a tab-delimited file with at least species and genetic_code columns.\n"
            )
            sys.exit(1)
        for row in reader:
            sp = (row.get("species") or "").strip().replace(" ", "_")
            code = (row.get("genetic_code") or "").strip()
            if not sp or not code:
                sys.stderr.write(f"Empty species or genetic_code entry was detected in {input_table}.\n")
                sys.exit(1)
            try:
                code_int = int(code)
            except ValueError:
                sys.stderr.write(f"Invalid genetic_code for {sp}: {code}\n")
                sys.exit(1)
            if code_int <= 0:
                sys.stderr.write(f"genetic_code must be a positive integer for {sp}: {code}\n")
                sys.exit(1)
            if sp in overrides:
                sys.stderr.write(f"Duplicate species entry in {input_table}: {sp}\n")
                sys.exit(1)
            overrides[sp] = str(code_int)

outfile.parent.mkdir(parents=True, exist_ok=True)
with outfile.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["species", "genetic_code", "source"])
    for sp in species:
        if sp in overrides:
            writer.writerow([sp, overrides[sp], "species_genetic_code.tsv"])
        else:
            writer.writerow([sp, default_code, "default"])

unknown = sorted(sp for sp in overrides if sp not in set(species))
if unknown:
    sys.stderr.write(
        "Warning: species_genetic_code.tsv entries were ignored because no matching species_cds file was found: "
        + ", ".join(unknown)
        + "\n"
    )
PY
}

lookup_species_genetic_code() {
  local species_name=$1
  local table_path=$2
  local code=""
  if [[ -s "${table_path}" ]]; then
    code=$(awk -F'\t' -v sp="${species_name}" 'NR>1 && $1==sp {print $2; exit}' "${table_path}")
  fi
  if [[ -z "${code}" ]]; then
    code="${genetic_code}"
  fi
  printf '%s' "${code}"
}

species_tree_rooting_method=""
species_tree_rooting_value=""

parse_species_tree_rooting() {
  local raw_config=$1
  local -n method_ref=$2
  local -n value_ref=$3
  local -a fields=()
  local -a normalized_fields=()
  local idx
  local method=""

  IFS=',' read -r -a fields <<< "${raw_config}"
  if [[ ${#fields[@]} -eq 0 ]]; then
    echo "species_tree_rooting is empty."
    echo 'Expected one of: "outgroup,GENUS_SPECIES[,GENUS_SPECIES...]", "midpoint", "mad", "mv", "taxonomy[,ncbi[,opentree,timetree...]]".'
    return 1
  fi

  for idx in "${!fields[@]}"; do
    normalized_fields[idx]=$(trim_ascii_whitespace "${fields[idx]}")
  done
  method=$(printf '%s' "${normalized_fields[0]}" | tr '[:upper:]' '[:lower:]')
  if [[ -z "${method}" ]]; then
    echo "species_tree_rooting is missing a rooting method: ${raw_config}"
    return 1
  fi

  case "${method}" in
    outgroup)
      if [[ ${#normalized_fields[@]} -lt 2 ]]; then
        echo "species_tree_rooting=${raw_config} is invalid."
        echo 'Outgroup rooting requires at least one species label after "outgroup,".'
        return 1
      fi
      value_ref=""
      for ((idx = 1; idx < ${#normalized_fields[@]}; idx++)); do
        if [[ -z "${normalized_fields[idx]}" ]]; then
          continue
        fi
        if [[ -n "${value_ref}" ]]; then
          value_ref+=","
        fi
        value_ref+="${normalized_fields[idx]}"
      done
      if [[ -z "${value_ref}" ]]; then
        echo "species_tree_rooting=${raw_config} is invalid."
        echo 'Outgroup rooting requires one or more non-empty species labels after "outgroup,".'
        return 1
      fi
      ;;
    midpoint | mad | mv)
      if [[ ${#normalized_fields[@]} -gt 1 ]]; then
        echo "species_tree_rooting=${raw_config} is invalid."
        echo "Method ${method} does not accept additional comma-separated fields."
        return 1
      fi
      value_ref=""
      ;;
    md)
      if [[ ${#normalized_fields[@]} -gt 1 ]]; then
        echo "species_tree_rooting=${raw_config} is invalid."
        echo 'Method md does not accept additional comma-separated fields.'
        return 1
      fi
      method="mv"
      value_ref=""
      ;;
    taxonomy)
      value_ref=""
      for ((idx = 1; idx < ${#normalized_fields[@]}; idx++)); do
        if [[ -z "${normalized_fields[idx]}" ]]; then
          continue
        fi
        if [[ -n "${value_ref}" ]]; then
          value_ref+=","
        fi
        value_ref+="${normalized_fields[idx]}"
      done
      ;;
    *)
      if species_tree_rooting_looks_like_label_list "${normalized_fields[@]}"; then
        method="outgroup"
        value_ref=""
        for ((idx = 0; idx < ${#normalized_fields[@]}; idx++)); do
          if [[ -z "${normalized_fields[idx]}" ]]; then
            continue
          fi
          if [[ -n "${value_ref}" ]]; then
            value_ref+=","
          fi
          value_ref+="${normalized_fields[idx]}"
        done
        echo "species_tree_rooting=${raw_config} uses legacy outgroup-label syntax; interpreting it as outgroup,${value_ref}."
      else
        echo "Invalid species_tree_rooting: ${raw_config}"
        echo 'species_tree_rooting must be one of "outgroup,GENUS_SPECIES[,GENUS_SPECIES...]", "midpoint", "mad", "mv", or "taxonomy[,ncbi[,opentree,timetree...]]".'
        return 1
      fi
      ;;
  esac

  method_ref="${method}"
  return 0
}

parse_mcmctree_constraint_record() {
  local raw_record=$1
  local __out_name=$2
  local -n __out_ref="${__out_name}"
  local idx

  IFS=',' read -r -a __out_ref <<< "${raw_record}"
  if [[ ${#__out_ref[@]} -ne 4 ]]; then
    echo "Invalid mcmctree divergence-time constraint: ${raw_record}"
    echo "Expected exactly 4 comma-separated fields: left_species,right_species,lower_bound,upper_bound"
    return 1
  fi

  for idx in 0 1 2 3; do
    __out_ref[${idx}]=$(trim_ascii_whitespace "${__out_ref[${idx}]}")
  done

  if [[ -z "${__out_ref[0]}" || -z "${__out_ref[1]}" || -z "${__out_ref[2]}" || -z "${__out_ref[3]}" ]]; then
    echo "Invalid mcmctree divergence-time constraint (empty field): ${raw_record}"
    echo "Expected non-empty fields: left_species,right_species,lower_bound,upper_bound"
    return 1
  fi

  return 0
}

clear_directory_contents_safe() {
  local target_dir=$1
  local target_real
  local entries=()

  if [[ -z "${target_dir}" || ! -d "${target_dir}" ]]; then
    echo "Refusing to clear missing directory: ${target_dir}" >&2
    return 1
  fi

  target_real=$(cd "${target_dir}" && pwd -P)
  if [[ -z "${target_real}" || "${target_real}" == "/" ]]; then
    echo "Refusing to clear unsafe directory: ${target_dir}" >&2
    return 1
  fi

  shopt -s dotglob nullglob
  entries=("${target_real}"/*)
  shopt -u dotglob nullglob
  if [[ ${#entries[@]} -gt 0 ]]; then
    rm -rf -- "${entries[@]}"
  fi
}
# Directories
dir_sp_cds="${gg_workspace_input_dir}/species_cds"
dir_sp_protein_input="$(species_protein_input_dir_path)"
file_species_genetic_code="$(species_genetic_code_table_path)"
file_species_genetic_code_resolved="${gg_workspace_downloads_dir}/tmp/species_genetic_code.resolved.tsv"
dir_og_rooted_tree="${gg_workspace_output_dir}/orthogroup/rooted_tree"
annotation_species_resolved=""
annotation_species_candidates=()
file_go_annotation=""
mapfile -t annotation_species_candidates < <(gg_species_names_from_annotation_dir "${gg_workspace_output_dir}/species_cds_annotation")
if [[ ${#annotation_species_candidates[@]} -eq 0 && "${input_sequence_mode}" == "protein" ]]; then
  mapfile -t annotation_species_candidates < <(gg_species_names_from_fasta_dir "${dir_sp_protein_input}")
fi
if [[ ${#annotation_species_candidates[@]} -eq 0 ]]; then
  mapfile -t annotation_species_candidates < <(gg_species_names_from_fasta_dir "${dir_sp_cds}")
fi
if annotation_species_resolved=$(gg_resolve_annotation_species "${annotation_species}" "${annotation_species_candidates[@]}"); then
  if [[ -n "${annotation_species_resolved}" ]]; then
    file_go_annotation="${gg_workspace_output_dir}/species_cds_annotation/${annotation_species_resolved}_annotation.tsv"
  fi
fi

# Species tree
dir_species_tree="${gg_workspace_output_dir}/species_tree"
if [[ "${input_sequence_mode}" == "protein" ]]; then
  dir_species_busco_full="${gg_workspace_output_dir}/species_protein_busco_full"
  dir_species_busco_short="${gg_workspace_output_dir}/species_protein_busco_short"
else
  dir_species_busco_full="${gg_workspace_output_dir}/species_cds_busco_full"
  dir_species_busco_short="${gg_workspace_output_dir}/species_cds_busco_short"
fi
dir_single_copy_fasta="${dir_species_tree}/single_copy_cds_fasta"
dir_single_copy_mafft="${dir_species_tree}/single_copy_mafft"
dir_single_copy_trimal="${dir_species_tree}/single_copy_trimal"
dir_single_copy_iqtree_dna="${dir_species_tree}/single_copy_iqtree_dna"
dir_single_copy_iqtree_pep="${dir_species_tree}/single_copy_iqtree_pep"
dir_astral_dna="${dir_species_tree}/single_copy_astral_dna"
dir_astral_pep="${dir_species_tree}/single_copy_astral_pep"
dir_species_tree_summary="${dir_species_tree}/species_tree_summary"
dir_concat_fasta="${dir_species_tree}/concatenated_alignment"
dir_concat_iqtree_dna="${dir_species_tree}/concatenated_iqtree_dna"
dir_concat_iqtree_pep="${dir_species_tree}/concatenated_iqtree_pep"
dir_mcmctree2="${dir_species_tree}/mcmctree_main"
dir_tmp="${dir_species_tree}/tmp"
dir_nwkit_download_dir="${gg_workspace_downloads_dir}/nwkit_downloads"

# Orthogroup
dir_sp_protein="${gg_workspace_downloads_dir}/tmp/species_protein"
dir_orthofinder="${gg_workspace_output_dir}/orthofinder"
dir_orthofinder_og="${dir_orthofinder}/Orthogroups"
dir_orthofinder_filtered="${dir_orthofinder}/Orthogroups_filtered"
dir_orthofinder_hog2og="${dir_orthofinder}/hog2og"

# Genome evolution
dir_genome_evolution="${gg_workspace_output_dir}/genome_evolution"
dir_species_omamer="${dir_genome_evolution}/omamer_search"
dir_species_omark="${dir_genome_evolution}/omark"
dir_busco_fasta="${dir_genome_evolution}/busco_cds_fasta"
dir_busco_mafft="${dir_genome_evolution}/busco_mafft"
dir_busco_trimal="${dir_genome_evolution}/busco_trimal"
dir_busco_iqtree_dna="${dir_genome_evolution}/busco_iqtree_dna"
dir_busco_iqtree_pep="${dir_genome_evolution}/busco_iqtree_pep"
dir_busco_notung_dna="${dir_genome_evolution}/busco_notung_dna"
dir_busco_notung_pep="${dir_genome_evolution}/busco_notung_pep"
dir_busco_rooted_txt_dna="${dir_genome_evolution}/busco_rooted_txt_dna"
dir_busco_rooted_txt_pep="${dir_genome_evolution}/busco_rooted_txt_pep"
dir_busco_rooted_nwk_dna="${dir_genome_evolution}/busco_rooted_nwk_dna"
dir_busco_rooted_nwk_pep="${dir_genome_evolution}/busco_rooted_nwk_pep"
dir_cafe="${dir_genome_evolution}/cafe"
dir_cafe_orthogroup_selection="${dir_cafe}/orthogroup_selection"
dir_cafe_output="${dir_cafe}/cafe_output"

# Output files
# Species tree
file_species_busco_summary_table="${dir_species_tree}/busco_summary_table/busco_summary.tsv"
file_astral_tree_dna_q1="${dir_astral_dna}/single_copy.astral.dna.q1.nwk"     # Quartet supports for the main topology; The lengths of terminal branches are not computed.
file_astral_tree_dna="${dir_astral_dna}/single_copy.astral.dna.optimized.nwk" # ASTRAL topology with branch lengths optimized by IQ-TREE.
file_astral_log_dna="${dir_astral_dna}/single_copy.astral.dna.log"
file_astral_tree_pep_q1="${dir_astral_pep}/single_copy.astral.pep.q1.nwk"     # Quartet supports for the main topology; The lengths of terminal branches are not computed.
file_astral_tree_pep="${dir_astral_pep}/single_copy.astral.pep.optimized.nwk" # ASTRAL topology with branch lengths optimized by IQ-TREE.
file_astral_log_pep="${dir_astral_pep}/single_copy.astral.pep.log"
file_concat_cds="${dir_concat_fasta}/concat.cds.trimal.fa.gz"
file_concat_pep="${dir_concat_fasta}/concat.pep.trimal.fa.gz"
file_concat_iqtree_dna="${dir_concat_iqtree_dna}/concat.dna.iqtree.nwk"
file_concat_iqtree_pep="${dir_concat_iqtree_pep}/concat.pep.iqtree.nwk"
file_concat_iqtree_dna_root="${dir_concat_iqtree_dna}/concat.dna.rooted.nwk"
file_concat_iqtree_pep_root="${dir_concat_iqtree_pep}/concat.pep.rooted.nwk"
file_constrained_tree="${dir_species_tree}/constrained_tree/constrained.nwk"
file_plot_constrained_tree="${dir_species_tree}/constrained_tree/constrained_tree_constraints.pdf"
file_iq2mc_prefix="${dir_species_tree}/mcmctree_parameter_estimation/iq2mc"
file_iq2mc_ctl="${file_iq2mc_prefix}.mcmctree.ctl"
file_iq2mc_hessian="${file_iq2mc_prefix}.mcmctree.hessian"
file_iq2mc_rooted_tree="${file_iq2mc_prefix}.rooted.nwk"
file_iq2mc_dummy_phy="${file_iq2mc_prefix}.dummy.phy"
file_mcmctree_raw_output="${dir_mcmctree2}/iq2mc.mcmctree.out"
file_mcmctree_figtree_tre="${dir_mcmctree2}/FigTree.tre"
file_mcmctree_dated_nwk="${dir_mcmctree2}/dated_species_tree.nwk"
file_plot_species_trees="${dir_species_tree_summary}/species_trees.pdf"
file_plot_mcmctree_pdf="${file_mcmctree_dated_nwk%.*}.pdf"
file_dated_species_tree="${dir_species_tree_summary}/dated_species_tree.nwk"
file_dated_species_tree_pdf="${dir_species_tree_summary}/dated_species_tree.pdf"
file_undated_species_tree="${dir_species_tree_summary}/undated_species_tree.nwk"

# Orthogroup
file_orthofinder_done_marker="${dir_orthofinder_hog2og}/README.txt"
file_orthogroup_selection="${dir_orthofinder_filtered}/Orthogroups.selected.tsv"
file_orthogroup_method_comparison="${dir_orthofinder}/orthogroup_method_comparison/orthogroup_method_comparison.pdf"

# Genome evolution
file_orthogroup_genecount_selected="${dir_orthofinder_filtered}/Orthogroups.GeneCount.selected.tsv"
file_genome_busco_summary_table="${dir_genome_evolution}/busco_summary_table/busco_summary.tsv"
file_species_omark_summary_table="${dir_genome_evolution}/omark_summary_table/omark_summary.tsv"
file_busco_grampa_dna="${dir_genome_evolution}/grampa_busco_dna/grampa_summary.tsv"
file_busco_grampa_pep="${dir_genome_evolution}/grampa_busco_pep/grampa_summary.tsv"
file_orthogroup_grampa="${dir_genome_evolution}/grampa_orthogroup/grampa_summary.tsv"
file_gene_id="${dir_orthofinder_filtered}/Orthogroups.selected.tsv"
file_cafe_summary_all_pdf="${dir_cafe}/summary_plot/summary_all.pdf"
file_cafe_summary_significant_pdf="${dir_cafe}/summary_plot/summary_significant.pdf"
file_go_enrichment_significant="${dir_cafe}/go_enrichment/enrichment_significant_${change_direction_go}_${target_branch_go}_significant_go.tsv"

# Runtime helpers
shared_species_busco_stage_done=0
shared_busco_summary_stage_done=0
shared_species_omark_stage_done=0
shared_omark_summary_stage_done=0

sync_genome_busco_summary_table_from_shared() {
  if [[ "${file_genome_busco_summary_table}" == "${file_species_busco_summary_table}" ]]; then
    return 0
  fi
  if [[ ! -s "${file_species_busco_summary_table}" ]]; then
    return 1
  fi
  ensure_parent_dir "${file_genome_busco_summary_table}"
  if [[ -s "${file_genome_busco_summary_table}" ]] && cmp -s "${file_species_busco_summary_table}" "${file_genome_busco_summary_table}"; then
    return 0
  fi
  cp_out "${file_species_busco_summary_table}" "${file_genome_busco_summary_table}"
}

run_shared_species_busco_stage() {
  local task="BUSCO analysis of species-wise input files"
  local source_species_input_dir=""
  local -a source_species_input_fasta=()
  local -a species_input_fasta=()
  local -a input_species_set=()
  local -a busco_output_files=()
  local seq_full seq_file sp_ub file_sp_busco_full file_sp_busco_short
  local input_species busco_file busco_base busco_species busco_species_found
  local dir_busco_db="" dir_busco_lineage=""
  local full_found=0 short_found=0
  local missing_busco_outputs=0

  if [[ ${shared_species_busco_stage_done} -eq 1 ]]; then
    return 0
  fi
  shared_species_busco_stage_done=1

  if [[ ${run_species_busco} -ne 1 ]]; then
    gg_step_skip "${task}"
    return 0
  fi

  source_species_input_dir=$(effective_species_input_source_dir_path)
  ensure_dir "${dir_species_busco_full}"
  ensure_dir "${dir_species_busco_short}"
  mapfile -t source_species_input_fasta < <(gg_find_fasta_files "${source_species_input_dir}" 1)
  echo "Number of ${input_sequence_mode} files for BUSCO: ${#source_species_input_fasta[@]}"
  if [[ ${#source_species_input_fasta[@]} -eq 0 ]]; then
    echo "No ${input_sequence_mode} file found. Exiting."
    exit 1
  fi
  mapfile -t input_species_set < <(gg_species_names_from_fasta_dir "${source_species_input_dir}")
  if ! resolve_busco_lineage_for_species_set "${input_species_set[@]}"; then
    exit 1
  fi
  mapfile -t busco_output_files < <(
    find "${dir_species_busco_full}" "${dir_species_busco_short}" -maxdepth 1 -type f \
      \( -name "*busco.full.tsv" -o -name "*busco.short.txt" \) \
      2> /dev/null | sort
  )
  for busco_file in "${busco_output_files[@]}"; do
    busco_base=$(basename "${busco_file}")
    busco_species=$(gg_species_name_from_path_or_dot "${busco_base}")
    busco_species_found=0
    for input_species in "${input_species_set[@]}"; do
      if [[ "${input_species}" == "${busco_species}" ]]; then
        busco_species_found=1
        break
      fi
    done
    if [[ ${busco_species_found} -eq 0 ]]; then
      echo "Removing stale BUSCO output for species not in current input: ${busco_file}"
      rm -f -- "${busco_file}"
    fi
  done
  for input_species in "${input_species_set[@]}"; do
    if busco_output_exists_for_species "${dir_species_busco_full}" "${input_species}" "*busco.full.tsv"; then
      full_found=1
    else
      full_found=0
    fi
    if busco_output_exists_for_species "${dir_species_busco_short}" "${input_species}" "*busco.short.txt"; then
      short_found=1
    else
      short_found=0
    fi
    if [[ ${full_found} -ne 1 || ${short_found} -ne 1 ]]; then
      missing_busco_outputs=1
      break
    fi
  done
  if [[ ${missing_busco_outputs} -ne 1 ]]; then
    for seq_full in "${source_species_input_fasta[@]}"; do
      echo "Skipped BUSCO: $(basename "${seq_full}")"
    done
    echo "$(date): End: ${task}"
    return 0
  fi

  prepare_species_tree_input_dir
  mapfile -t species_input_fasta < <(gg_find_fasta_files "${species_tree_input_dir}" 1)
  for seq_full in "${species_input_fasta[@]}"; do
    seq_file=$(basename "${seq_full}")
    sp_ub=$(gg_species_name_from_path_or_dot "${seq_file}")
    file_sp_busco_full="${dir_species_busco_full}/${sp_ub}.busco.full.tsv"
    file_sp_busco_short="${dir_species_busco_short}/${sp_ub}.busco.short.txt"
    if busco_output_exists_for_species "${dir_species_busco_full}" "${sp_ub}" "*busco.full.tsv"; then
      full_found=1
    else
      full_found=0
    fi
    if busco_output_exists_for_species "${dir_species_busco_short}" "${sp_ub}" "*busco.short.txt"; then
      short_found=1
    else
      short_found=0
    fi
    if [[ ${full_found} -ne 1 || ${short_found} -ne 1 ]]; then
      gg_step_start "${task}: ${seq_file}"
      remove_busco_outputs_for_species "${dir_species_busco_full}" "${sp_ub}" "*busco.full.tsv"
      remove_busco_outputs_for_species "${dir_species_busco_short}" "${sp_ub}" "*busco.short.txt"
      seqkit seq --threads "${GG_TASK_CPUS}" "${species_tree_input_dir}/${seq_file}" --out-file "tmp.busco_input.fasta"

      if ! dir_busco_db=$(ensure_busco_download_path "${gg_workspace_dir}" "${busco_lineage_resolved}"); then
        echo "Failed to prepare BUSCO dataset: ${busco_lineage_resolved}"
        exit 1
      fi
      dir_busco_lineage="${dir_busco_db}/lineages/${busco_lineage_resolved}"

      busco \
        --in "tmp.busco_input.fasta" \
        --mode "${species_tree_busco_mode}" \
        --out "busco_tmp" \
        --cpu "${GG_TASK_CPUS}" \
        --force \
        --evalue 1e-03 \
        --limit 20 \
        --lineage_dataset "${dir_busco_lineage}" \
        --download_path "${dir_busco_db}" \
        --offline

      if copy_busco_tables "./busco_tmp" "${busco_lineage_resolved}" "${file_sp_busco_full}" "${file_sp_busco_short}"; then
        rm -rf -- "./busco_tmp"
      else
        echo "Failed to locate normalized BUSCO outputs for ${sp_ub}. Exiting."
        exit 1
      fi
    else
      echo "Skipped BUSCO: ${seq_file}"
    fi
  done
  echo "$(date): End: ${task}"
}

run_shared_busco_summary_stage() {
  local task="Collecting IDs of common BUSCO genes"
  local num_busco_ids=0
  local source_species_input_dir=""

  if [[ ${shared_busco_summary_stage_done} -eq 1 ]]; then
    sync_genome_busco_summary_table_from_shared || true
    return 0
  fi
  shared_busco_summary_stage_done=1

  if [[ ${run_species_get_busco_summary} -ne 1 ]]; then
    gg_step_skip "${task}"
    return 0
  fi

  source_species_input_dir=$(effective_species_input_source_dir_path)
  normalize_busco_table_naming "${dir_species_busco_full}" "${dir_species_busco_short}"
  if ! is_species_set_identical "${source_species_input_dir}" "${dir_species_busco_full}"; then
    echo "Exiting due to species-set mismatch between ${source_species_input_dir} and ${dir_species_busco_full}"
    exit 1
  fi
  if [[ ! -s "${file_species_busco_summary_table}" ]]; then
    gg_step_start "${task}"
    ensure_parent_dir "${file_species_busco_summary_table}"

    python "${gg_support_dir}/collect_common_BUSCO_genes.py" \
      --busco_outdir "${dir_species_busco_full}" \
      --ncpu "${GG_TASK_CPUS}" \
      --outfile "tmp.busco_summary_table.tsv"
    mv_out "tmp.busco_summary_table.tsv" "${file_species_busco_summary_table}"

    num_busco_ids=$(get_busco_summary_gene_count "${file_species_busco_summary_table}")
    echo "Number of BUSCO genes: ${num_busco_ids}"
    echo "$(date): End: ${task}"
  else
    gg_step_skip "${task}"
  fi
  sync_genome_busco_summary_table_from_shared || true
}

run_shared_species_omark_stage() {
  local task="OMArk analysis of species-wise protein input files"
  local source_species_input_dir=""
  local -a source_species_input_fasta=()
  local -a species_input_fasta=()
  local -a input_species_set=()
  local -a existing_species_dirs=()
  local protein_full protein_file sp_ub omamer_out omark_outdir sum_file omamer_query
  local input_species existing_dir existing_species existing_found
  local omark_db_file=""
  local missing_omark_outputs=0

  if [[ ${shared_species_omark_stage_done} -eq 1 ]]; then
    return 0
  fi
  shared_species_omark_stage_done=1

  if [[ ${run_species_omark} -ne 1 ]]; then
    gg_step_skip "${task}"
    return 0
  fi

  ensure_dir "${dir_species_omamer}"
  ensure_dir "${dir_species_omark}"
  source_species_input_dir=$(effective_species_input_source_dir_path)
  mapfile -t source_species_input_fasta < <(gg_find_fasta_files "${source_species_input_dir}" 1)
  echo "Number of protein files for OMArk: ${#source_species_input_fasta[@]}"
  if [[ ${#source_species_input_fasta[@]} -eq 0 ]]; then
    echo "No protein file found for OMArk. Exiting."
    exit 1
  fi

  mapfile -t input_species_set < <(gg_species_names_from_fasta_dir "${source_species_input_dir}")

  mapfile -t existing_species_dirs < <(find "${dir_species_omark}" -maxdepth 1 -mindepth 1 -type d ! -name '.*' | sort)
  for existing_dir in "${existing_species_dirs[@]}"; do
    existing_species=$(basename "${existing_dir}")
    existing_found=0
    for input_species in "${input_species_set[@]}"; do
      if [[ "${input_species}" == "${existing_species}" ]]; then
        existing_found=1
        break
      fi
    done
    if [[ ${existing_found} -eq 0 ]]; then
      echo "Removing stale OMArk output for species not in current input: ${existing_dir}"
      rm -rf -- "${existing_dir}"
      rm -f -- "${dir_species_omamer}/${existing_species}.omamer"
    fi
  done
  for input_species in "${input_species_set[@]}"; do
    omamer_out="${dir_species_omamer}/${input_species}.omamer"
    sum_file="${dir_species_omark}/${input_species}/${input_species}.sum"
    if [[ ! -s "${omamer_out}" || ! -s "${sum_file}" ]]; then
      missing_omark_outputs=1
      break
    fi
  done
  if [[ ${missing_omark_outputs} -ne 1 ]]; then
    for protein_full in "${source_species_input_fasta[@]}"; do
      echo "Skipped OMArk: $(basename "${protein_full}")"
    done
    echo "$(date): End: ${task}"
    return 0
  fi

  if ! omark_db_file=$(ensure_omark_database "${gg_workspace_dir}" "${omark_db_path}"); then
    echo "Failed to prepare OMArk database: ${omark_db_path}"
    exit 1
  fi
  omark_db_resolved="${omark_db_file}"
  prepare_species_protein_tmp
  mapfile -t species_input_fasta < <(gg_find_fasta_files "${dir_sp_protein}" 1)

  for protein_full in "${species_input_fasta[@]}"; do
    protein_file=$(basename "${protein_full}")
    sp_ub=$(gg_species_name_from_path_or_dot "${protein_file}")
    omamer_out="${dir_species_omamer}/${sp_ub}.omamer"
    omark_outdir="${dir_species_omark}/${sp_ub}"
    sum_file="${omark_outdir}/${sp_ub}.sum"

    if [[ -s "${omamer_out}" && -s "${sum_file}" ]]; then
      echo "Skipped OMArk: ${protein_file}"
      continue
    fi

    gg_step_start "${task}: ${protein_file}"
    ensure_dir "${omark_outdir}"
    if [[ ! -s "${omamer_out}" ]]; then
      omamer_query="${omark_outdir}/${sp_ub}.query.fa"
      stage_species_protein_fasta "${protein_full}" "${omamer_query}"
      omamer search \
        --db "${omark_db_file}" \
        --query "${omamer_query}" \
        --out "${omamer_out}"
    fi
    if [[ ! -s "${omamer_out}" ]]; then
      echo "OMAmer search output was not created: ${omamer_out}"
      exit 1
    fi

    if [[ ! -s "${sum_file}" ]]; then
      omark \
        -f "${omamer_out}" \
        -d "${omark_db_file}" \
        -o "${omark_outdir}"
    fi
    if [[ ! -s "${sum_file}" ]]; then
      echo "OMArk summary output was not created: ${sum_file}"
      exit 1
    fi
  done
  echo "$(date): End: ${task}"
}

run_shared_omark_summary_stage() {
  local task="Summarizing OMArk species quality results"

  if [[ ${shared_omark_summary_stage_done} -eq 1 ]]; then
    return 0
  fi
  shared_omark_summary_stage_done=1

  if [[ ${run_species_get_omark_summary} -ne 1 ]]; then
    gg_step_skip "${task}"
    return 0
  fi

  if [[ ! -s "${file_species_omark_summary_table}" ]]; then
    gg_step_start "${task}"
    ensure_parent_dir "${file_species_omark_summary_table}"

    python "${gg_support_dir}/summarize_omark.py" \
      --omark_outdir "${dir_species_omark}" \
      --outfile "tmp.omark_summary.tsv"
    mv_out "tmp.omark_summary.tsv" "${file_species_omark_summary_table}"
    echo "$(date): End: ${task}"
  else
    gg_step_skip "${task}"
  fi
}

species_protein_input_has_files() {
  local protein_files=()
  mapfile -t protein_files < <(gg_find_fasta_files "${dir_sp_protein_input}" 1)
  [[ ${#protein_files[@]} -gt 0 ]]
}

species_cds_input_has_files() {
  local cds_files=()
  mapfile -t cds_files < <(gg_find_fasta_files "${dir_sp_cds}" 1)
  [[ ${#cds_files[@]} -gt 0 ]]
}

effective_species_input_source_dir_path() {
  if [[ "${input_sequence_mode}" == "protein" ]] && species_protein_input_has_files; then
    echo "${dir_sp_protein_input}"
  else
    echo "${dir_sp_cds}"
  fi
}

compute_shared_protein_input_signature() {
  local input_file
  local -a input_files=()
  local -a stat_lines=()
  local metadata_source="species_cds"
  if [[ -s "${file_species_genetic_code}" ]]; then
    input_files+=( "${file_species_genetic_code}" )
  fi
  if [[ "${input_sequence_mode}" == "protein" ]] && species_protein_input_has_files; then
    metadata_source="species_protein"
    while IFS= read -r input_file; do
      input_files+=( "${input_file}" )
    done < <(gg_find_fasta_files "${dir_sp_protein_input}" 1)
  else
    while IFS= read -r input_file; do
      input_files+=( "${input_file}" )
    done < <(gg_find_fasta_files "${dir_sp_cds}" 1)
  fi

  if [[ ${#input_files[@]} -gt 0 ]]; then
    if stat --version > /dev/null 2>&1; then
      mapfile -t stat_lines < <(stat -c '%n:%s:%Y' "${input_files[@]}")
    else
      mapfile -t stat_lines < <(stat -f '%N:%z:%m' "${input_files[@]}")
    fi
  fi

  {
    printf 'input_sequence_mode=%s\n' "${input_sequence_mode}"
    printf 'genetic_code=%s\n' "${genetic_code}"
    printf 'metadata_source=%s\n' "${metadata_source}"
    printf 'species_tree_busco_mode=%s\n' "${species_tree_busco_mode}"
    printf '%s\n' "${stat_lines[@]}"
  } | cksum | awk '{print $1}'
}

print_effective_genome_evolution_config_summary() {
  local species_protein_input_available=0
  local species_cds_input_available=0
  local species_genetic_code_table_present=0

  if species_protein_input_has_files; then
    species_protein_input_available=1
  fi
  if species_cds_input_has_files; then
    species_cds_input_available=1
  fi
  if [[ -s "${file_species_genetic_code}" ]]; then
    species_genetic_code_table_present=1
  fi

  gg_print_registered_config_summary \
    "gg_genome_evolution_entrypoint.sh" \
    "effective config summary (gg_genome_evolution_core.sh)" \
    genetic_code \
    busco_lineage \
    annotation_species \
    annotation_species_resolved \
    omark_db_path \
    omark_db_resolved \
    species_tree_rooting_method \
    species_tree_rooting_value \
    species_tree_busco_mode \
    species_tree_sequence_label \
    species_protein_input_available \
    species_cds_input_available \
    species_genetic_code_table_present
}

refresh_dir_for_shared_protein_input_signature() {
  local target_dir=$1
  local description=$2
  local signature=$3
  local stamp_file="${target_dir}/.shared_protein_input_signature"
  local previous_signature=""

  ensure_dir "${target_dir}"
  if [[ -s "${stamp_file}" ]]; then
    previous_signature=$(< "${stamp_file}")
  fi
  if [[ -n "${previous_signature}" && "${previous_signature}" != "${signature}" ]]; then
    echo "Shared protein input signature changed for ${description}. Clearing derived outputs in ${target_dir}"
    if ! clear_directory_contents_safe "${target_dir}"; then
      echo "Failed to clear ${description} directory after input signature change: ${target_dir}"
      exit 1
    fi
  fi
  ensure_dir "${target_dir}"
  printf '%s\n' "${signature}" > "${stamp_file}"
}

cleanup_species_protein_tmp() {
  local cleanup_target
  for cleanup_target in "${dir_sp_protein}" "${dir_sp_protein}_orthofinder" "${dir_sp_protein}_core" "${dir_sp_protein}_additional"; do
    if [[ -d "${cleanup_target}" ]]; then
      echo "Removing temporary species_protein directory: ${cleanup_target}"
      rm -rf -- "${cleanup_target}"
    fi
  done
  species_protein_ready=0
}

stage_species_protein_fasta() {
  local source_fasta=$1
  local staged_fasta=$2
  ensure_parent_dir "${staged_fasta}"
  seqkit seq --threads "${GG_TASK_CPUS}" "${source_fasta}" --out-file "${staged_fasta}"
  if [[ ! -s "${staged_fasta}" ]]; then
    echo "Temporary species_protein FASTA was not created: ${staged_fasta}"
    exit 1
  fi
}

prepare_species_protein_orthofinder_dir() {
  local source_dir=$1
  local target_dir=$2
  local source_fasta source_name source_stem staged_fasta
  if [[ -e "${target_dir}" ]]; then
    rm -rf -- "${target_dir}"
  fi
  ensure_dir "${target_dir}"
  while IFS= read -r source_fasta; do
    source_name=$(basename "${source_fasta}")
    source_stem=${source_name%.gz}
    source_stem=${source_stem%.*}
    staged_fasta="${target_dir}/${source_stem}.fa"
    stage_species_protein_fasta "${source_fasta}" "${staged_fasta}"
  done < <(gg_find_fasta_files "${source_dir}" 1)
}

prepare_species_protein_tmp() {
  local cds_path cds sp_ub translated_file species_code protein_path
  local -a cds_files=()
  local -a protein_files=()

  if [[ ${species_protein_ready} -eq 1 ]]; then
    return 0
  fi

  gg_step_start "Prepare temporary species_protein FASTA files"
  cleanup_species_protein_tmp
  ensure_dir "${dir_sp_protein}"

  if [[ "${input_sequence_mode}" == "protein" ]] && species_protein_input_has_files; then
    check_species_protein_dir "${dir_sp_protein_input}"
    check_if_species_files_unique "${dir_sp_protein_input}"
    mapfile -t protein_files < <(gg_find_fasta_files "${dir_sp_protein_input}" 1)
    for protein_path in "${protein_files[@]}"; do
      sp_ub=$(gg_species_name_from_path "$(basename "${protein_path}")")
      translated_file="${sp_ub}.fa.gz"
      echo "Copying protein FASTA: $(basename "${protein_path}") -> ${translated_file}"
      seqkit seq --threads "${GG_TASK_CPUS}" "${protein_path}" --out-file "${dir_sp_protein}/${translated_file}"
    done
    if [[ -s "${file_species_genetic_code}" ]]; then
      echo "species_genetic_code.tsv is ignored because species_protein inputs are provided: ${file_species_genetic_code}"
    fi
    species_protein_source="species_protein"
    species_protein_ready=1
    return 0
  fi
  if [[ "${input_sequence_mode}" != "protein" ]] && species_protein_input_has_files; then
    echo "Ignoring species_protein inputs in cds mode: ${dir_sp_protein_input}"
    echo "cds mode always generates temporary species_protein FASTA files from species_cds."
  fi

  if ! species_cds_input_has_files; then
    echo "protein mode requires either species_protein inputs or species_cds inputs for fallback translation."
    echo "Checked: ${dir_sp_protein_input} and ${dir_sp_cds}"
    exit 1
  fi
  if [[ ${run_cds_translation} -ne 1 ]]; then
    echo "run_cds_translation must be 1 when species proteins need to be generated from species_cds."
    exit 1
  fi

  check_species_cds "${gg_workspace_dir}"
  check_if_species_files_unique "${dir_sp_cds}"
  prepare_species_genetic_code_table "${dir_sp_cds}" "${genetic_code}" "${file_species_genetic_code_resolved}"
  mapfile -t cds_files < <(gg_find_fasta_files "${dir_sp_cds}" 1)
  for cds_path in "${cds_files[@]}"; do
    cds=$(basename "${cds_path}")
    if zgrep -q -e "^>.*[[:blank:]]" "${cds_path}"; then
      echo "Space (\" \") is detected. Please remove all annotation info after spaces in sequence names. Exiting: ${cds}"
      exit 1
    fi
    sp_ub=$(gg_species_name_from_path "${cds}")
    translated_file="${sp_ub}.fa.gz"
    species_code=$(lookup_species_genetic_code "${sp_ub}" "${file_species_genetic_code_resolved}")
    echo "Translation started: ${cds} (genetic_code=${species_code}) -> ${translated_file}"

    seqkit seq --remove-gaps --threads "${GG_TASK_CPUS}" "${cds_path}" |
      gg_prepare_cds_fasta_stream "${GG_TASK_CPUS}" "${species_code}" |
      seqkit translate --transl-table "${species_code}" --threads "${GG_TASK_CPUS}" |
      sed -e "s/^>${sp_ub}[-_\.]/>/" -e "s/^>/>${sp_ub}_/" |
      sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' |
      seqkit seq --threads "${GG_TASK_CPUS}" --out-file "${dir_sp_protein}/${translated_file}"
  done
  species_protein_source="species_cds"
  species_protein_ready=1
}

prepare_species_tree_input_dir() {
  if [[ "${input_sequence_mode}" == "protein" ]]; then
    prepare_species_protein_tmp
    species_tree_input_dir="${dir_sp_protein}"
  else
    species_tree_input_dir="${dir_sp_cds}"
  fi
}

trap cleanup_species_protein_tmp EXIT

# Runtime setup
if [[ "${input_sequence_mode}" == "protein" ]]; then
  if ! species_protein_input_has_files && ! species_cds_input_has_files; then
    echo "protein mode requires either species_protein inputs or species_cds inputs."
    echo "Checked: ${dir_sp_protein_input} and ${dir_sp_cds}"
    exit 1
  fi
  if [[ "${undated_species_tree}" == "iqtree_dna" || "${undated_species_tree}" == "astral_dna" ]]; then
    echo "protein mode does not support undated_species_tree=${undated_species_tree}."
    echo 'Use iqtree_pep or astral_pep instead.'
    exit 1
  fi
  if [[ ${run_concat_iqtree_dna} -eq 1 || ${run_individual_iqtree_dna} -eq 1 || ${run_astral_dna} -eq 1 ]]; then
    echo "Disabling DNA-only species-tree steps in protein mode: run_concat_iqtree_dna, run_individual_iqtree_dna, run_astral_dna"
    run_concat_iqtree_dna=0
    run_individual_iqtree_dna=0
    run_astral_dna=0
  fi
  if [[ ${run_constrained_tree} -eq 1 || ${run_plot_constrained_tree} -eq 1 || ${run_mcmctree1} -eq 1 || ${run_mcmctree2} -eq 1 || ${run_convert_tree_format} -eq 1 || ${run_plot_mcmctreer} -eq 1 ]]; then
    echo "Disabling CDS-only dating steps in protein mode: run_constrained_tree, run_plot_constrained_tree, run_mcmctree1, run_mcmctree2, run_convert_tree_format, run_plot_mcmctreer"
    run_constrained_tree=0
    run_plot_constrained_tree=0
    run_mcmctree1=0
    run_mcmctree2=0
    run_convert_tree_format=0
    run_plot_mcmctreer=0
  fi
  if [[ ${run_busco_iqtree_dna} -eq 1 || ${run_busco_notung_root_dna} -eq 1 || ${run_busco_root_dna} -eq 1 || ${run_busco_grampa_dna} -eq 1 ]]; then
    echo "Disabling DNA-only BUSCO genome-evolution steps in protein mode: run_busco_iqtree_dna, run_busco_notung_root_dna, run_busco_root_dna, run_busco_grampa_dna"
    run_busco_iqtree_dna=0
    run_busco_notung_root_dna=0
    run_busco_root_dna=0
    run_busco_grampa_dna=0
  fi
else
  check_species_cds "${gg_workspace_dir}"
  check_if_species_files_unique "${dir_sp_cds}"
fi
shared_protein_input_signature=$(compute_shared_protein_input_signature)
refresh_dir_for_shared_protein_input_signature "${dir_species_tree}" "species_tree" "${shared_protein_input_signature}"
refresh_dir_for_shared_protein_input_signature "${dir_orthofinder}" "orthofinder" "${shared_protein_input_signature}"
refresh_dir_for_shared_protein_input_signature "${dir_genome_evolution}" "genome_evolution" "${shared_protein_input_signature}"
memory_notung=${GG_MEM_PER_CPU_GB}

ensure_dir "${dir_species_tree_summary}"
ensure_dir "${dir_tmp}"
cd "${dir_tmp}"

enable_all_run_flags_for_debug_mode
orthogroup_annotation_method=$(echo "${orthogroup_annotation_method:-mmseqs2}" | tr '[:upper:]' '[:lower:]')
if [[ "${orthogroup_annotation_method}" != "blastp" && "${orthogroup_annotation_method}" != "mmseqs2" ]]; then
  echo "Invalid orthogroup_annotation_method: ${orthogroup_annotation_method}"
  echo 'orthogroup_annotation_method must be either "blastp" or "mmseqs2". Exiting.'
  exit 1
fi

if ! parse_species_tree_rooting "${species_tree_rooting}" species_tree_rooting_method species_tree_rooting_value; then
  exit 1
fi
echo "Resolved species_tree_rooting method: ${species_tree_rooting_method}"
if [[ -n "${species_tree_rooting_value}" ]]; then
  echo "Resolved species_tree_rooting value: ${species_tree_rooting_value}"
fi
print_effective_genome_evolution_config_summary

root_species_tree() {
  local infile=$1
  local outfile=$2
  local tree_description=$3
  local root_log="${dir_tmp}/tmp.nwkit.root.$$.log"
  local outgroup_label_list=()
  local root_method="${species_tree_rooting_method}"
  local root_value="${species_tree_rooting_value}"
  local missing_outgroup=0
  local root_exit_code=0
  local -a nwkit_root_args=()
  ensure_parent_dir "${outfile}"
  rm -f -- "${outfile}" "${root_log}"

  nwkit_root_args=(--method "${root_method}" --infile "${infile}" --outfile "${outfile}")
  if [[ "${root_method}" == "outgroup" ]]; then
    mapfile -t outgroup_label_list < <(printf '%s' "${root_value}" | tr ',' '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' -e '/^$/d')
    for outgroup_label in "${outgroup_label_list[@]}"; do
      if ! grep -q -F -- "${outgroup_label}" "${infile}"; then
        missing_outgroup=1
        break
      fi
    done
    if [[ ${missing_outgroup} -eq 1 ]]; then
      echo "Error: Outgroup labels (${root_value}) are not present in ${tree_description}."
      return 1
    fi
    nwkit_root_args+=(--outgroup "${root_value}")
  elif [[ "${root_method}" == "taxonomy" ]]; then
    ensure_dir "${dir_nwkit_download_dir}"
    nwkit_root_args+=(--download_dir "${dir_nwkit_download_dir}")
    if [[ -n "${root_value}" ]]; then
      nwkit_root_args+=(--taxonomy_source "${root_value}")
    fi
  fi

  if nwkit root "${nwkit_root_args[@]}" 2> "${root_log}"; then
    root_exit_code=0
  else
    root_exit_code=$?
  fi

  if [[ ${root_exit_code} -eq 0 && -s "${outfile}" ]]; then
    rm -f -- "${root_log}"
    return 0
  fi

  echo "Error: Failed to root ${tree_description}."
  [[ -s "${root_log}" ]] && cat "${root_log}"
  rm -f -- "${root_log}"
  return 1
}

normalize_iq2mc_constraint_tree() {
  local infile=$1
  local tmpfile

  if [[ ! -s "${infile}" ]]; then
    echo "Error: Constraint tree file is missing or empty: ${infile}"
    return 1
  fi

  tmpfile=$(mktemp "${dir_tmp}/tmp.iq2mc.constraint.XXXXXX")
  # IQ-TREE's dated-tree parser rejects spaces after commas inside B(...) labels.
  tr -d '\r' < "${infile}" | sed -E 's/,[[:space:]]+/,/g; s/[[:space:]]+$//' > "${tmpfile}"

  if [[ ! -s "${tmpfile}" ]]; then
    echo "Error: Failed to normalize IQ2MC constraint tree: ${infile}"
    rm -f -- "${tmpfile}"
    return 1
  fi

  mv_out "${tmpfile}" "${infile}"
  return 0
}

iq2mc_option_supported() {
  local candidate=$1
  local resolved_candidate
  local probe_prefix
  local probe_aln
  local probe_tree
  local probe_rc=1
  local probe_hessian=""
  resolved_candidate=$(command -v "${candidate}" 2> /dev/null || true)

  # Avoid pipefail+SIGPIPE false negatives from `... | grep -q`.
  if grep -q -- "--mcmc-bds" < <("${candidate}" -h 2>&1 || true); then
    return 0
  fi
  if command -v strings > /dev/null 2>&1 && [[ -n "${resolved_candidate}" ]]; then
    if grep -q -- "--mcmc-bds" < <(strings "${resolved_candidate}" 2> /dev/null || true); then
      return 0
    fi
  fi

  # Fallback: run a tiny dating job to detect support reliably.
  probe_prefix="${dir_tmp}/tmp.iq2mc_probe.$$"
  probe_aln="${probe_prefix}.fa"
  probe_tree="${probe_prefix}.nwk"
  probe_hessian="${probe_prefix}.mcmctree.hessian"
  cat > "${probe_aln}" << 'EOF'
>a
ACGT
>b
ACGT
>c
ACGT
EOF
  cat > "${probe_tree}" << 'EOF'
((a:0.1,b:0.1):0.1,c:0.1);
EOF
  if [[ -s "${probe_aln}" && -s "${probe_tree}" ]]; then
    if "${candidate}" \
      -s "${probe_aln}" \
      -te "${probe_tree}" \
      -m JC \
      -n 0 \
      -T 1 \
      --dating mcmctree \
      --mcmc-bds 1,1,0.5 \
      --prefix "${probe_prefix}" \
      --redo > /dev/null 2>&1; then
      probe_rc=0
    else
      probe_rc=$?
    fi
  fi
  if [[ ${probe_rc} -eq 0 && -s "${probe_hessian}" ]]; then
    rm -f -- "${probe_aln}" "${probe_tree}" "${probe_prefix}".* 2> /dev/null || true
    return 0
  fi
  rm -f -- "${probe_aln}" "${probe_tree}" "${probe_prefix}".* 2> /dev/null || true

  return 1
}

resolve_iq2mc_binary() {
  local candidate
  for candidate in iqtree3 iqtree; do
    if ! command -v "${candidate}" > /dev/null 2>&1; then
      continue
    fi
    if iq2mc_option_supported "${candidate}"; then
      echo "${candidate}"
      return 0
    fi
  done
  return 1
}

count_newick_tips() {
  local tree_file=$1
  python - "${tree_file}" << 'PY'
import sys
from Bio import Phylo
try:
    tree = Phylo.read(sys.argv[1], "newick")
    print(len(tree.get_terminals()))
except Exception:
    print(0)
PY
}

build_astral_input() {
  local tree_dir=$1
  local merged_file=$2
  local min_tips=$3
  local kept skipped ntips tree_file
  local -a tree_files

  kept=0
  skipped=0
  : > "${merged_file}"

  mapfile -t tree_files < <(find "${tree_dir}" -maxdepth 1 -type f -name "*.nwk" | sort)
  if [[ ${#tree_files[@]} -eq 0 ]]; then
    echo "No gene-tree files were found for ASTRAL in: ${tree_dir}"
    rm -f -- "${merged_file}"
    return 1
  fi

  for tree_file in "${tree_files[@]}"; do
    ntips=$(count_newick_tips "${tree_file}")
    if [[ "${ntips}" -ge "${min_tips}" ]]; then
      cat "${tree_file}" >> "${merged_file}"
      echo "" >> "${merged_file}"
      kept=$((kept + 1))
    else
      skipped=$((skipped + 1))
      echo "Skipped for ASTRAL (tips=${ntips} < ${min_tips}): $(basename "${tree_file}")"
    fi
  done

  echo "ASTRAL input tree filtering summary: kept=${kept}, skipped=${skipped}, min_tips=${min_tips}"
  if [[ ${kept} -eq 0 ]]; then
    rm -f -- "${merged_file}"
    return 1
  fi
  return 0
}

optimize_astral_tree_branch_lengths() {
  local astral_support_tree=$1
  local concat_alignment=$2
  local model=$3
  local optimized_outfile=$4
  local tag=$5
  local tmp_topology="tmp.astral.${tag}.topology.nwk"
  local iqtree_prefix="tmp.astral.${tag}.optimized"
  local tmp_concat_alignment="tmp.astral.${tag}.concat.fasta"
  local tmp_rooted_optimized_tree="tmp.astral.${tag}.optimized.rooted.nwk"
  local tmp_optimized_tree_with_support="tmp.astral.${tag}.optimized.with_support.nwk"

  if [[ ! -s "${astral_support_tree}" ]]; then
    echo "Warning: Cannot optimize ASTRAL tree branch lengths. Missing tree: ${astral_support_tree}"
    return 1
  fi
  if [[ ! -s "${concat_alignment}" ]]; then
    echo "Warning: Cannot optimize ASTRAL tree branch lengths. Missing alignment: ${concat_alignment}"
    return 1
  fi

  rm -f -- "${tmp_topology}" "${iqtree_prefix}".* "${tmp_concat_alignment}" "${tmp_rooted_optimized_tree}" "${tmp_optimized_tree_with_support}"
  if ! nwkit drop \
    --infile "${astral_support_tree}" \
    --outfile "${tmp_topology}" \
    --target intnode \
    --name yes \
    --support yes; then
    echo "Warning: Failed to sanitize ASTRAL topology for IQ-TREE branch-length optimization."
    rm -f -- "${tmp_topology}" "${iqtree_prefix}".* "${tmp_concat_alignment}" "${tmp_rooted_optimized_tree}" "${tmp_optimized_tree_with_support}"
    return 1
  fi
  seqkit seq --threads 1 "${concat_alignment}" --out-file "${tmp_concat_alignment}"
  if [[ ! -s "${tmp_concat_alignment}" ]]; then
    echo "Warning: Failed to prepare alignment for IQ-TREE branch-length optimization: ${concat_alignment}"
    rm -f -- "${tmp_topology}" "${iqtree_prefix}".* "${tmp_concat_alignment}" "${tmp_rooted_optimized_tree}" "${tmp_optimized_tree_with_support}"
    return 1
  fi

  iqtree \
    -s "${tmp_concat_alignment}" \
    -te "${tmp_topology}" \
    -m "${model}" \
    -n 0 \
    -T "${GG_TASK_CPUS}" \
    --prefix "${iqtree_prefix}" \
    --seed 12345 \
    --redo

  if [[ ! -s "${iqtree_prefix}.treefile" ]]; then
    echo "Warning: IQ-TREE optimization did not produce a tree for ${tag}."
    rm -f -- "${tmp_topology}" "${iqtree_prefix}".* "${tmp_concat_alignment}" "${tmp_rooted_optimized_tree}" "${tmp_optimized_tree_with_support}"
    return 1
  fi

  local root_exit_code=0
  if root_species_tree \
    "${iqtree_prefix}.treefile" \
    "${tmp_rooted_optimized_tree}" \
    "optimized ASTRAL ${tag} tree"; then
    root_exit_code=0
  else
    root_exit_code=$?
  fi
  if [[ ${root_exit_code} -ne 0 || ! -s "${tmp_rooted_optimized_tree}" ]]; then
    echo "Warning: Failed to save optimized ASTRAL ${tag} tree."
    rm -f -- "${tmp_topology}" "${iqtree_prefix}".* "${tmp_concat_alignment}" "${tmp_rooted_optimized_tree}" "${tmp_optimized_tree_with_support}"
    return 1
  fi

  local transfer_exit_code=0
  if nwkit transfer \
    --infile "${tmp_rooted_optimized_tree}" \
    --infile2 "${astral_support_tree}" \
    --target intnode \
    --support yes \
    --name no \
    --length no \
    --outfile "${tmp_optimized_tree_with_support}"; then
    transfer_exit_code=0
  else
    transfer_exit_code=$?
  fi
  local copy_exit_code=0
  if [[ ${transfer_exit_code} -eq 0 && -s "${tmp_optimized_tree_with_support}" ]]; then
    if cp_out "${tmp_optimized_tree_with_support}" "${optimized_outfile}"; then
      copy_exit_code=0
    else
      copy_exit_code=$?
    fi
  else
    echo "Warning: Failed to transfer ASTRAL support values to optimized ${tag} tree. Keeping branch-length-optimized tree without support transfer."
    if cp_out "${tmp_rooted_optimized_tree}" "${optimized_outfile}"; then
      copy_exit_code=0
    else
      copy_exit_code=$?
    fi
  fi
  rm -f -- "${tmp_topology}" "${iqtree_prefix}".* "${tmp_concat_alignment}" "${tmp_rooted_optimized_tree}" "${tmp_optimized_tree_with_support}"

  if [[ ${copy_exit_code} -ne 0 || ! -s "${optimized_outfile}" ]]; then
    echo "Warning: Failed to write optimized ASTRAL ${tag} tree."
    return 1
  fi
  return 0
}

iq2mc_binary=""
if [[ ${run_mcmctree1} -eq 1 || ${run_mcmctree2} -eq 1 ]]; then
  echo "IQ2MC dating will run. Checking tools and calibration species."
  if ! iq2mc_binary=$(resolve_iq2mc_binary); then
    echo "Error: IQ2MC-capable IQ-TREE was not found."
    echo "Please install an IQ-TREE build that supports --dating mcmctree (and --mcmc-bds)."
    exit 1
  fi
  if ! iq2mc_option_supported "${iq2mc_binary}"; then
    echo "Warning: Could not verify --mcmc-bds support in ${iq2mc_binary} by help/strings checks."
    echo "Warning: Proceeding with ${iq2mc_binary}; IQ2MC step 2 will fail if this binary is incompatible."
  fi
  if ! command -v mcmctree > /dev/null 2>&1; then
    echo "Error: mcmctree command was not found."
    echo "Please install the IQ2MC-compatible mcmctree binary."
    exit 1
  fi
  echo "Using IQ2MC binary: ${iq2mc_binary}"
  outgroup_label_list=()
  if [[ "${species_tree_rooting_method}" == "outgroup" ]]; then
    mapfile -t outgroup_label_list < <(printf '%s' "${species_tree_rooting_value}" | tr ',' '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' -e '/^$/d')
    for outgroup_label in "${outgroup_label_list[@]}"; do
      stop_if_species_not_found_in "${dir_sp_cds}" "${outgroup_label}"
    done
  fi
  if [[ ${timetree_constraint} -eq 0 ]]; then
    if [[ ${#mcmctree_divergence_time_constraints[@]} -eq 0 ]]; then
      echo "timetree_constraint=0 requires mcmctree_divergence_time_constraints_str with one or more records."
      echo "Example: Arabidopsis_thaliana,Oryza_sativa,130,-|Arabidopsis_thaliana,Amborella_trichopoda,150,200"
      exit 1
    fi
    mcmctree_params=()
    if ! parse_mcmctree_constraint_record "${mcmctree_divergence_time_constraints[0]}" mcmctree_params; then
      exit 1
    fi
    stop_if_species_not_found_in "${dir_sp_cds}" "${mcmctree_params[0]}"
    stop_if_species_not_found_in "${dir_sp_cds}" "${mcmctree_params[1]}"
  fi
  echo ""
fi

task="BUSCO analysis of species-wise input files"
run_shared_species_busco_stage

task="Collecting IDs of common BUSCO genes"
run_shared_busco_summary_stage

task="Generating fasta files for individual single-copy genes"
ensure_dir "${dir_single_copy_fasta}"
num_busco_ids=$(get_busco_summary_gene_count "${file_species_busco_summary_table}")
singlecopy_fasta_files=()
mapfile -t singlecopy_fasta_files < <(gg_find_file_basenames "${dir_single_copy_fasta}" "${single_copy_fasta_glob}")
num_singlecopy_fasta=${#singlecopy_fasta_files[@]}
if [[ ${num_busco_ids} -ne ${num_singlecopy_fasta} && ${run_individual_get_fasta} -eq 1 ]]; then
  prepare_species_tree_input_dir
  gg_step_start "${task}"

  generate_single_copy_fasta() {
    local busco_id
    busco_id=$(awk -v row="$1" 'NR==row {print $1; exit}' "${file_species_busco_summary_table}")
    local remove_nonsingle=$2
    local outfile1="${dir_single_copy_fasta}/${busco_id}${single_copy_fasta_suffix}"
    if [[ -s "${outfile1}" ]]; then
      return 0
    fi
    local genes=()
    IFS=$'\t' read -r -a genes <<< "$(sed -n "${1}P" "${file_species_busco_summary_table}" | cut -f 4-)"
    echo "busco_id: ${busco_id}"
    if [[ ${strictly_single_copy_only} -eq 1 ]]; then
      if [[ "${genes[@]}" == *" -"* || "${genes[@]}" == *"- "* ]]; then
        echo "Skipping. ${busco_id} has missing gene(s)."
        return 0
      fi
      if [[ "${genes[@]}" == *","* ]]; then
        echo "Skipping. ${busco_id} has duplicated gene(s)."
        return 0
      fi
    fi
    local genes1=()
    read -r -a genes1 <<< "$(printf "%s " "${genes[@]}" | sed -e "s/[[:space:]]-[[:space:]]/ /g" -e "s/[^[[:space:]]]*,[^[[:space:]]]*//g" -e "s/[[:space:]]+/ /g")"
    local num_gene=${#genes1[@]}
    if [[ ${num_gene} -eq 0 ]]; then
      echo "Skipping. ${busco_id} has no genes in the selected species."
      return 0
    fi
    if [[ ! -s "${outfile1}" ]]; then
      if [[ "${input_sequence_mode}" == "protein" ]]; then
        gg_seqkit_grep_by_patterns_from_infile_list 1 "species_tree_input_fasta_list.txt" "${genes1[@]}" |
          seqkit replace --pattern " .*" --replacement "" --ignore-case --threads 1 |
          sed -e "s/_/|/" -e "s/_.*//" -e "s/|/_/" |
          seqkit seq --threads 1 --out-file "${outfile1}"
      else
        local pattern_args=()
        for gene in "${genes1[@]}"; do
          pattern_args+=(--pattern "${gene}")
        done
        seqkit grep --threads 1 "${pattern_args[@]}" --infile-list "species_tree_input_fasta_list.txt" |
          seqkit replace --pattern X --replacement N --by-seq --ignore-case --threads 1 |
          seqkit replace --pattern " .*" --replacement "" --ignore-case --threads 1 |
          cdskit pad |
          sed -e "s/_/|/" -e "s/_.*//" -e "s/|/_/" |
          seqkit seq --threads 1 --out-file "${outfile1}"
      fi
      if [[ ! -s "${outfile1}" ]]; then
        echo "File is empty. Removing: ${outfile1}"
        rm -f -- "${outfile1}"
      fi
    fi
    local fasta_genes=()
    if [[ -s "${outfile1}" ]]; then
      mapfile -t fasta_genes < <(seqkit seq --name --threads 1 "${outfile1}")
    fi
    local num_seq=${#fasta_genes[@]}
    # this block needs to be disabeld for ${strictly_single_copy_only} -eq 1, because orthogroups won't be complete
    if [[ ${num_gene} -ne ${num_seq} && ${strictly_single_copy_only} -eq 1 ]]; then
      echo "${busco_id}: Error. Number of genes and sequences did not match."
      echo "Genes in the orthogroup or BLAST hit:"
      printf '%s\n' "${genes1[@]}"
      echo ""
      echo "Genes in the generated FASTA:"
      printf '%s\n' "${fasta_genes[@]}"
      echo ""
      echo "Check duplicated sequence names in the species input FASTA files. Exiting."
      rm -f -- "${outfile1}"
      exit 1
    fi
  }

  gg_find_fasta_files "${species_tree_input_dir}" 1 > species_tree_input_fasta_list.txt
  num_busco_ids=$(get_busco_summary_gene_count "${file_species_busco_summary_table}")
  for ((i = 2; i <= num_busco_ids + 1; i++)); do # starting from 2 because the line 1 is header.
    wait_until_jobn_le ${GG_TASK_CPUS}
    generate_single_copy_fasta ${i} ${strictly_single_copy_only} &
  done
  wait_for_background_jobs
  rm -f -- species_tree_input_fasta_list.txt
  rm -f -- tmp.*
else
  gg_step_skip "${task}"
fi

task="MAFFT alignment"
ensure_dir "${dir_single_copy_mafft}"
num_busco_ids=$(get_busco_summary_gene_count "${file_species_busco_summary_table}")
mafft_fasta_files=()
mapfile -t mafft_fasta_files < <(gg_find_file_basenames "${dir_single_copy_mafft}" "${single_copy_aln_glob}")
num_mafft_fasta=${#mafft_fasta_files[@]}
if [[ ${num_busco_ids} -ne ${num_mafft_fasta} && ${run_individual_mafft} -eq 1 ]]; then
  gg_step_start "${task}"

  run_mafft() {
    local infile=$1
    local infile_base=${infile%%.*}
    local outfile="${dir_single_copy_mafft}/${infile_base}${single_copy_aln_suffix}"
    if [[ -s "${outfile}" ]]; then
      return 0
    fi
    local infile_path="${dir_single_copy_fasta}/${infile}"
    local num_seq=$(gg_count_fasta_records "${infile_path}")
    if [[ ${num_seq} -lt 2 ]]; then
      echo "Skipped. At least 2 sequences are necessary for MAFFT: ${infile}"
      return 0
    fi
    echo "$(date): start mafft: ${infile_base}"
    if [[ "${input_sequence_mode}" == "protein" ]]; then
      seqkit seq --threads 1 "${infile_path}" --out-file "tmp.${infile_base}.input.pep.fasta"
      mafft \
        --auto \
        --amino \
        --thread 1 \
        "tmp.${infile_base}.input.pep.fasta" \
        > "tmp.${infile_base}.pep.aln.fasta"
      if [[ -s "tmp.${infile_base}.pep.aln.fasta" ]]; then
        seqkit seq --threads 1 "tmp.${infile_base}.pep.aln.fasta" --out-file "tmp.${infile_base}.pep.aln.out.fa.gz"
        mv_out "tmp.${infile_base}.pep.aln.out.fa.gz" "${outfile}"
      fi
    else
      seqkit seq --threads 1 "${infile_path}" --out-file "tmp.${infile_base}.input.cds.fasta"
      cdskit mask \
        --seqfile "tmp.${infile_base}.input.cds.fasta" \
        --outfile "tmp.${infile_base}.cds.fasta"
      seqkit translate \
        --allow-unknown-codon \
        --transl-table "${genetic_code}" \
        --threads 1 \
        "tmp.${infile_base}.cds.fasta" \
        > "tmp.${infile_base}.pep.fasta"
      mafft \
        --auto \
        --thread 1 \
        "tmp.${infile_base}.pep.fasta" \
        > "tmp.${infile_base}.pep.aln.fasta"
      cdskit backalign \
        --seqfile "tmp.${infile_base}.cds.fasta" \
        --aa_aln "tmp.${infile_base}.pep.aln.fasta" \
        --codontable "${genetic_code}" \
        --outfile "tmp.${infile_base}.cds.aln.fasta"
      if [[ -s "tmp.${infile_base}.cds.aln.fasta" ]]; then
        seqkit seq --threads 1 "tmp.${infile_base}.cds.aln.fasta" --out-file "tmp.${infile_base}.cds.aln.out.fa.gz"
        mv_out "tmp.${infile_base}.cds.aln.out.fa.gz" "${outfile}"
      fi
    fi
    rm -f -- "tmp.${infile_base}"*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_single_copy_fasta}" "${single_copy_fasta_glob}")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    run_mafft "${input_alignment_file}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

task="TrimAl"
ensure_dir "${dir_single_copy_trimal}"
trimal_fasta_files=()
mapfile -t trimal_fasta_files < <(gg_find_file_basenames "${dir_single_copy_trimal}" "${single_copy_trimal_glob}")
num_trimal_fasta=${#trimal_fasta_files[@]}
if [[ ${num_busco_ids} -ne ${num_trimal_fasta} && ${run_individual_trimal} -eq 1 ]]; then
  gg_step_start "${task}"

  run_trimal() {
    local infile=$1
    local infile_base=${infile%%.*}
    local outfile="${dir_single_copy_trimal}/${infile_base}${single_copy_trimal_suffix}"
    if [[ ! -s "${outfile}" ]]; then
      if [[ "${input_sequence_mode}" == "protein" ]]; then
        seqkit seq --threads 1 "${dir_single_copy_mafft}/${infile}" --out-file "tmp.${infile_base}.pep.aln.fasta"
        trimal \
          -in "tmp.${infile_base}.pep.aln.fasta" \
          -out "tmp.${infile_base}.trimal.fasta" \
          -automated1
      else
        seqkit seq --remove-gaps --threads 1 "${dir_single_copy_mafft}/${infile}" > "tmp.${infile_base}.degap.fasta"
        seqkit translate --transl-table "${genetic_code}" --threads 1 "${dir_single_copy_mafft}/${infile}" > "tmp.${infile_base}.pep.fasta"
        trimal \
          -in tmp.${infile_base}.pep.fasta \
          -backtrans tmp.${infile_base}.degap.fasta \
          -out tmp.${infile_base}.trimal.fasta \
          -ignorestopcodon \
          -automated1
      fi
      if [[ -s "tmp.${infile_base}.trimal.fasta" ]]; then
        seqkit seq --threads 1 "tmp.${infile_base}.trimal.fasta" --out-file "tmp.${infile_base}.trimal.out.fa.gz"
        mv_out "tmp.${infile_base}.trimal.out.fa.gz" "${outfile}"
      fi
      rm -f -- "tmp.${infile_base}."*
    fi
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_single_copy_mafft}" "${single_copy_aln_glob}")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    run_trimal "${input_alignment_file}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

task="Concatenating single-copy CDS fasta files"
concat_alignment_ready=0
if [[ "${input_sequence_mode}" == "protein" ]]; then
  if [[ -s "${file_concat_pep}" ]]; then
    concat_alignment_ready=1
  fi
else
  if [[ -s "${file_concat_cds}" && -s "${file_concat_pep}" ]]; then
    concat_alignment_ready=1
  fi
fi
if [[ ${concat_alignment_ready} -eq 0 && ${run_concat_alignment} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_parent_dir "${file_concat_pep}"
  ensure_dir "${dir_concat_fasta}"
  if [[ "${input_sequence_mode}" != "protein" ]]; then
    ensure_parent_dir "${file_concat_cds}"
  fi
  mapfile -t trimal_files < <(find "${dir_single_copy_trimal}" -maxdepth 1 -type f -name "${single_copy_trimal_glob}" | sort)
  if [[ ${#trimal_files[@]} -eq 0 ]]; then
    echo "No trimmed single-copy FASTA files were found in: ${dir_single_copy_trimal}"
    exit 1
  fi

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    if [[ ${strictly_single_copy_only} -eq 0 ]]; then
      concat_pep_tmp="tmp.concat.pep.fa.gz"
      seqkit concat --full --fill "-" --threads "${GG_TASK_CPUS}" "${trimal_files[@]}" |
        seqkit sort --threads "${GG_TASK_CPUS}" |
        seqkit seq --threads "${GG_TASK_CPUS}" --out-file "${concat_pep_tmp}"
    else
      concat_pep_tmp="tmp.concat.pep.fa.gz"
      seqkit concat --threads "${GG_TASK_CPUS}" "${trimal_files[@]}" |
        seqkit sort --threads "${GG_TASK_CPUS}" |
        seqkit seq --threads "${GG_TASK_CPUS}" --out-file "${concat_pep_tmp}"
    fi
    mv_out "${concat_pep_tmp}" "${file_concat_pep}"
  else
    if [[ ${strictly_single_copy_only} -eq 0 ]]; then
      concat_cds_tmp="tmp.concat.cds.fa.gz"
      seqkit concat --full --fill "-" --threads "${GG_TASK_CPUS}" "${trimal_files[@]}" |
        seqkit sort --threads "${GG_TASK_CPUS}" |
        seqkit seq --threads "${GG_TASK_CPUS}" --out-file "${concat_cds_tmp}"
    else
      concat_cds_tmp="tmp.concat.cds.fa.gz"
      seqkit concat --threads "${GG_TASK_CPUS}" "${trimal_files[@]}" |
        seqkit sort --threads "${GG_TASK_CPUS}" |
        seqkit seq --threads "${GG_TASK_CPUS}" --out-file "${concat_cds_tmp}"
    fi
    mv_out "${concat_cds_tmp}" "${file_concat_cds}"
    seqkit translate --transl-table "${genetic_code}" --threads "${GG_TASK_CPUS}" "${file_concat_cds}" |
      seqkit seq --threads "${GG_TASK_CPUS}" --out-file "tmp.concat.pep.fa.gz"
    mv_out "tmp.concat.pep.fa.gz" "${file_concat_pep}"
  fi
else
  gg_step_skip "${task}"
fi

task="IQ-TREE of the concatenated alignment with a Protein evolution model"
if [[ ! -s "${file_concat_iqtree_pep}" && ${run_concat_iqtree_protein} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_concat_iqtree_pep}"

  ntaxa=$(gg_count_fasta_records "${file_concat_pep}")
  bootstrap_params=(--ufboot 1000 --bnni)
  if [[ ${ntaxa} -lt 4 ]]; then
    bootstrap_params=()
  fi
  iqtree_mem_args=()
  if [[ -n "${GG_MEM_TOTAL_GB:-}" ]]; then
    iqtree_mem_args=(-mem "${GG_MEM_TOTAL_GB}G")
  fi

  cd "${dir_concat_iqtree_pep}"
  concat_pep_local="tmp.concat.pep.input.fasta"
  seqkit seq --threads 1 "${file_concat_pep}" --out-file "./${concat_pep_local}"
  if iqtree \
    -s "${concat_pep_local}" \
    -st AA \
    -m "${protein_model}" \
    -pre "${concat_pep_local}" \
    -nt "${GG_TASK_CPUS}" \
    "${iqtree_mem_args[@]}" \
    -seed 12345 \
    "${bootstrap_params[@]}"; then
    iqtree_exit_status=0
  else
    iqtree_exit_status=$?
  fi

  if [[ ${iqtree_exit_status} -ne 0 ]]; then
    echo "Error. IQ-TREE exited with non-zero status: ${iqtree_exit_status}"
  else
    echo "IQ-TREE successfully exited with zero status: ${iqtree_exit_status}"
    if [[ ${ntaxa} -lt 4 ]]; then
      out_nwk="${concat_pep_local}.treefile"
    else
      out_nwk="${concat_pep_local}.contree"
    fi
    cp_out "${out_nwk}" "${file_concat_iqtree_pep}"
    rm -f -- "./${concat_pep_local}"
  fi
  cd "${dir_tmp}"
else
  gg_step_skip "${task}"
fi

task="Rooting of IQ-TREE's protein tree"
if [[ ! -s "${file_concat_iqtree_pep_root}" && ${run_concat_iqtree_protein} -eq 1 ]]; then
  gg_step_start "${task}"

  root_species_tree \
    "${file_concat_iqtree_pep}" \
    "${file_concat_iqtree_pep_root}" \
    "concatenated protein tree"

  if [[ -s "${file_concat_iqtree_pep_root}" ]]; then
    Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="${file_concat_iqtree_pep_root}"
  fi

  if [[ -s "${file_concat_iqtree_pep_root}" && ${undated_species_tree} == 'iqtree_pep' ]]; then
    echo "Undated species tree is copied,"
    echo "from: ${file_concat_iqtree_pep_root}"
    echo "to: ${file_undated_species_tree}"
    cp_out "${file_concat_iqtree_pep_root}" "${file_undated_species_tree}"
    Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="${file_undated_species_tree}"
  fi
else
  gg_step_skip "${task}"
fi

task="IQ-TREE of the concatenated alignment with a DNA evolution model"
if [[ ! -s "${file_concat_iqtree_dna}" && ${run_concat_iqtree_dna} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_concat_iqtree_dna}"

  ntaxa=$(gg_count_fasta_records "${file_concat_cds}")
  bootstrap_params=(--ufboot 1000 --bnni)
  if [[ ${ntaxa} -lt 4 ]]; then
    bootstrap_params=()
  fi
  iqtree_mem_args=()
  if [[ -n "${GG_MEM_TOTAL_GB:-}" ]]; then
    iqtree_mem_args=(-mem "${GG_MEM_TOTAL_GB}G")
  fi

  cd "${dir_concat_iqtree_dna}"
  concat_cds_local="tmp.concat.cds.input.fasta"
  seqkit seq --threads 1 "${file_concat_cds}" --out-file "./${concat_cds_local}"
  if iqtree \
    -s "${concat_cds_local}" \
    -st DNA \
    -m "${nucleotide_model}" \
    -pre "${concat_cds_local}" \
    -nt "${GG_TASK_CPUS}" \
    "${iqtree_mem_args[@]}" \
    -seed 12345 \
    "${bootstrap_params[@]}"; then
    iqtree_exit_status=0
  else
    iqtree_exit_status=$?
  fi

  if [[ ${iqtree_exit_status} -ne 0 ]]; then
    echo "Error. IQ-TREE exited with non-zero status: ${iqtree_exit_status}"
  else
    echo "IQ-TREE successfully exited with zero status: ${iqtree_exit_status}"
    if [[ ${ntaxa} -lt 4 ]]; then
      out_nwk="${concat_cds_local}.treefile"
    else
      out_nwk="${concat_cds_local}.contree"
    fi
    cp_out "${out_nwk}" "${file_concat_iqtree_dna}"
    rm -f -- "./${concat_cds_local}"
  fi
  cd "${dir_tmp}"
else
  gg_step_skip "${task}"
fi

task="Rooting of IQ-TREE's DNA tree"
if [[ ! -s "${file_concat_iqtree_dna_root}" && ${run_concat_iqtree_dna} -eq 1 ]]; then
  gg_step_start "${task}"

  root_species_tree \
    "${file_concat_iqtree_dna}" \
    "${file_concat_iqtree_dna_root}" \
    "concatenated DNA tree"

  if [[ -s "${file_concat_iqtree_dna_root}" ]]; then
    Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="${file_concat_iqtree_dna_root}"
  fi

  if [[ -s "${file_concat_iqtree_dna_root}" && ${undated_species_tree} == 'iqtree_dna' ]]; then
    echo "Undated species tree is copied,"
    echo "from: ${file_concat_iqtree_dna_root}"
    echo "to: ${file_undated_species_tree}"
    cp_out "${file_concat_iqtree_dna_root}" "${file_undated_species_tree}"
    Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="${file_undated_species_tree}"
  fi
else
  gg_step_skip "${task}"
fi

task="IQ-TREE for individual single-copy protein trees"
ensure_dir "${dir_single_copy_iqtree_pep}"
iqtree_pep_tree_files=()
mapfile -t iqtree_pep_tree_files < <(gg_find_file_basenames "${dir_single_copy_iqtree_pep}" "*.nwk")
num_iqtree_pep=${#iqtree_pep_tree_files[@]}
if [[ ${num_busco_ids} -ne ${num_iqtree_pep} && ${run_individual_iqtree_pep} -eq 1 ]]; then
  gg_step_start "${task}"

  run_iqtree_pep() {
    local infile=$1
    local infile_base=${infile%%.*}
    local outfile="${dir_single_copy_iqtree_pep}/${infile_base}.pep.nwk"
    if [[ -s "${outfile}" ]]; then
      return 0
    fi
    local num_seq=$(gg_count_fasta_records "${dir_single_copy_trimal}/${infile}")
    if [[ ${num_seq} -lt 3 ]]; then
      echo "Skipped. At least 3 sequences are necessary for IQ-TREE: ${infile}"
      return 0
    fi
    if [[ "${input_sequence_mode}" == "protein" ]]; then
      seqkit seq --threads 1 "${dir_single_copy_trimal}/${infile}" > "tmp.${infile_base}.pep.fasta"
    else
      seqkit translate --transl-table "${genetic_code}" --threads 1 "${dir_single_copy_trimal}/${infile}" \
        > "tmp.${infile_base}.pep.fasta"
    fi
    iqtree \
      -s "tmp.${infile_base}.pep.fasta" \
      -m "${protein_model}" \
      -T 1 \
      --prefix "tmp.${infile_base}" \
      --seed 12345 \
      --redo
    mv_out tmp."${infile_base}".treefile "${outfile}"
    rm -f -- "tmp.${infile_base}."*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_single_copy_trimal}" "${single_copy_trimal_glob}")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    run_iqtree_pep "${input_alignment_file}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

task="ASTRAL of individual single-copy protein trees"
if [[ (! -s "${file_astral_tree_pep}" || ! -s "${file_astral_log_pep}") && ${run_astral_pep} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_astral_pep}"

  if compgen -G "tmp.astral.*" > /dev/null; then
    rm -f -- tmp.astral.*
  fi

  if ! build_astral_input "${dir_single_copy_iqtree_pep}" "tmp.astral.merged.iqtree.nwk" "${astral_min_tips}"; then
    echo "Skipped. No eligible protein gene trees for ASTRAL after filtering."
  else
    astral-hybrid \
      --input "tmp.astral.merged.iqtree.nwk" \
      --output "tmp.astral.out.tree" \
      --mode 3 \
      --support 2 \
      --thread "${GG_TASK_CPUS}" \
      2> "tmp.astral.log.txt"

    labels=("pp1" "pp2" "pp3" "f1" "f2" "f3" "q1" "q2" "q3") # https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md
    for i in "${!labels[@]}"; do
      tmp_label_tree="tmp.astral.pep.${labels[i]}.tmp.nwk"
      python "${gg_support_dir}/extract_astral_support_label.py" \
        --infile "tmp.astral.out.tree" \
        --outfile "${tmp_label_tree}" \
        --label_key "${labels[i]}"
      root_species_tree \
        "${tmp_label_tree}" \
        "single_copy.astral.pep.${labels[i]}.nwk" \
        "ASTRAL protein tree (${labels[i]})"
      if [[ -s "single_copy.astral.pep.${labels[i]}.nwk" ]]; then
        Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="single_copy.astral.pep.${labels[i]}.nwk"
      fi
      rm -f -- "${tmp_label_tree}"
    done

    if [[ -s "single_copy.astral.pep.q1.nwk" ]]; then
      shopt -s nullglob
      astral_pep_outputs=(single_copy.astral.pep.*)
      shopt -u nullglob
      if [[ ${#astral_pep_outputs[@]} -eq 0 ]]; then
        echo "ASTRAL protein outputs were expected but not found."
        exit 1
      fi
      mv_out "${astral_pep_outputs[@]}" "${dir_astral_pep}"
      mv_out "tmp.astral.log.txt" "${file_astral_log_pep}"
      echo "For more information on support values (e.g., f1, f2, pp1, q1, ...), please refer to: https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md" > "${dir_astral_pep}/README.txt"
    fi

    if [[ -s "${file_astral_tree_pep_q1}" ]]; then
      if optimize_astral_tree_branch_lengths "${file_astral_tree_pep_q1}" "${file_concat_pep}" "${protein_model}" "${file_astral_tree_pep}" "pep"; then
        Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="${file_astral_tree_pep}"
      else
        echo "Warning: Falling back to unoptimized ASTRAL protein tree."
        cp_out "${file_astral_tree_pep_q1}" "${file_astral_tree_pep}"
      fi
    fi
  fi
  rm -f -- tmp.astral.*

  if [[ -s "${file_astral_tree_pep}" && ${undated_species_tree} == 'astral_pep' ]]; then
    echo "Undated species tree is copied,"
    echo "from: ${file_astral_tree_pep}"
    echo "to: ${file_undated_species_tree}"
    nwkit drop --name yes --target intnode --infile "${file_astral_tree_pep}" | nwkit label --target intnode --force yes --outfile "tmp.undated_species_tree.nwk"
    if [[ -s "tmp.undated_species_tree.nwk" ]]; then
      mv_out "tmp.undated_species_tree.nwk" "${file_undated_species_tree}"
    fi
    if [[ ! -s "${file_undated_species_tree}" ]]; then
      echo "Warning: Failed to convert optimized ASTRAL protein tree with nwkit drop|label. Copying optimized tree instead."
      cp_out "${file_astral_tree_pep}" "${file_undated_species_tree}"
    fi
    if [[ -s "${file_undated_species_tree}" ]]; then
      Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="${file_undated_species_tree}"
    fi
  fi
else
  gg_step_skip "${task}"
fi

task="IQ-TREE for individual single-copy DNA trees"
ensure_dir "${dir_single_copy_iqtree_dna}"
iqtree_dna_tree_files=()
mapfile -t iqtree_dna_tree_files < <(gg_find_file_basenames "${dir_single_copy_iqtree_dna}" "*.nwk")
num_iqtree_dna=${#iqtree_dna_tree_files[@]}
if [[ ${num_busco_ids} -ne ${num_iqtree_dna} && ${run_individual_iqtree_dna} -eq 1 ]]; then
  gg_step_start "${task}"

  run_iqtree_dna() {
    local infile=$1
    local infile_base=${infile%%.*}
    local outfile="${dir_single_copy_iqtree_dna}/${infile_base}.dna.nwk"
    if [[ -s "${outfile}" ]]; then
      return 0
    fi
    local num_seq=$(gg_count_fasta_records "${dir_single_copy_trimal}/${infile}")
    if [[ ${num_seq} -lt 3 ]]; then
      echo "Skipped. At least 3 sequences are necessary for IQ-TREE: ${infile}"
      return 0
    fi
    cp_out "${dir_single_copy_trimal}/${infile}" "./tmp.${infile_base}.input.fasta"
    iqtree \
      -s "./tmp.${infile_base}.input.fasta" \
      -m "${nucleotide_model}" \
      -T 1 \
      --prefix "tmp.${infile_base}" \
      --seed 12345 \
      --redo
    mv_out tmp."${infile_base}".treefile "${outfile}"
    rm -f -- "tmp.${infile_base}."*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_single_copy_trimal}" "*.trimal.fa.gz")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    run_iqtree_dna "${input_alignment_file}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

task="ASTRAL of individual single-copy DNA trees"
if [[ (! -s "${file_astral_tree_dna}" || ! -s "${file_astral_log_dna}") && ${run_astral_dna} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_astral_dna}"

  if compgen -G "tmp.astral.*" > /dev/null; then
    rm -f -- tmp.astral.*
  fi

  if ! build_astral_input "${dir_single_copy_iqtree_dna}" "tmp.astral.merged.iqtree.nwk" "${astral_min_tips}"; then
    echo "Skipped. No eligible DNA gene trees for ASTRAL after filtering."
  else
    astral-hybrid \
      --input "tmp.astral.merged.iqtree.nwk" \
      --output "tmp.astral.out.tree" \
      --mode 3 \
      --support 2 \
      --thread "${GG_TASK_CPUS}" \
      2> "tmp.astral.log.txt"

    labels=("pp1" "pp2" "pp3" "f1" "f2" "f3" "q1" "q2" "q3") # https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md
    for i in "${!labels[@]}"; do
      tmp_label_tree="tmp.astral.dna.${labels[i]}.tmp.nwk"
      python "${gg_support_dir}/extract_astral_support_label.py" \
        --infile "tmp.astral.out.tree" \
        --outfile "${tmp_label_tree}" \
        --label_key "${labels[i]}"
      root_species_tree \
        "${tmp_label_tree}" \
        "single_copy.astral.dna.${labels[i]}.nwk" \
        "ASTRAL DNA tree (${labels[i]})"
      if [[ -s "single_copy.astral.dna.${labels[i]}.nwk" ]]; then
        Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="single_copy.astral.dna.${labels[i]}.nwk"
      fi
      rm -f -- "${tmp_label_tree}"
    done

    if [[ -s "single_copy.astral.dna.q1.nwk" ]]; then
      shopt -s nullglob
      astral_dna_outputs=(single_copy.astral.dna.*)
      shopt -u nullglob
      if [[ ${#astral_dna_outputs[@]} -eq 0 ]]; then
        echo "ASTRAL DNA outputs were expected but not found."
        exit 1
      fi
      mv_out "${astral_dna_outputs[@]}" "${dir_astral_dna}"
      mv_out "tmp.astral.log.txt" "${file_astral_log_dna}"
      echo "For more information on support values (e.g., f1, f2, pp1, q1, ...), please refer to: https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md" > "${dir_astral_dna}/README.txt"
    fi

    if [[ -s "${file_astral_tree_dna_q1}" ]]; then
      if optimize_astral_tree_branch_lengths "${file_astral_tree_dna_q1}" "${file_concat_cds}" "${nucleotide_model}" "${file_astral_tree_dna}" "dna"; then
        Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="${file_astral_tree_dna}"
      else
        echo "Warning: Falling back to unoptimized ASTRAL DNA tree."
        cp_out "${file_astral_tree_dna_q1}" "${file_astral_tree_dna}"
      fi
    fi
  fi
  rm -f -- tmp.astral.*

  if [[ -s "${file_astral_tree_dna}" && ${undated_species_tree} == 'astral_dna' ]]; then
    echo "Undated species tree is copied,"
    echo "from: ${file_astral_tree_dna}"
    echo "to: ${file_undated_species_tree}"
    nwkit drop --name yes --target intnode --infile "${file_astral_tree_dna}" | nwkit label --target intnode --force yes --outfile "tmp.undated_species_tree.nwk"
    if [[ -s "tmp.undated_species_tree.nwk" ]]; then
      mv_out "tmp.undated_species_tree.nwk" "${file_undated_species_tree}"
    fi
    if [[ ! -s "${file_undated_species_tree}" ]]; then
      echo "Warning: Failed to convert optimized ASTRAL DNA tree with nwkit drop|label. Copying optimized tree instead."
      cp_out "${file_astral_tree_dna}" "${file_undated_species_tree}"
    fi
    if [[ -s "${file_undated_species_tree}" ]]; then
      Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="${file_undated_species_tree}"
    fi
  fi
else
  gg_step_skip "${task}"
fi

task="Species tree plotting"
disable_if_no_input_file "run_plot_species_trees" "${file_concat_iqtree_dna_root}" "${file_concat_iqtree_pep_root}" "${file_astral_tree_dna}" "${file_astral_tree_pep}"
if [[ ! -s "${file_plot_species_trees}" && ${run_plot_species_trees} -eq 1 ]]; then
  gg_step_start "${task}"

  Rscript "${gg_support_dir}/plot_species_trees.r" \
    --iqtree_dna_nwk="${file_concat_iqtree_dna_root}" \
    --iqtree_pep_nwk="${file_concat_iqtree_pep_root}" \
    --iqtree_dna_log="${dir_concat_iqtree_dna}/tmp.concat.cds.input.fasta.log" \
    --iqtree_pep_log="${dir_concat_iqtree_pep}/tmp.concat.pep.input.fasta.log" \
    --astral_dna_nwk="${file_astral_tree_dna}" \
    --astral_pep_nwk="${file_astral_tree_pep}" \
    --astral_dna_log="${file_astral_log_dna}" \
    --astral_pep_log="${file_astral_log_pep}"

  if [[ -s "species_trees.pdf" ]]; then
    echo "Output file found for the task: ${task}"
    mv_out "species_trees.pdf" "${file_plot_species_trees}"
  fi
else
  gg_step_skip "${task}"
fi

task="Time-constrained tree preparation"
disable_if_no_input_file "run_constrained_tree" "${file_undated_species_tree}"
if [[ ! -s "${file_constrained_tree}" && ${run_constrained_tree} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_parent_dir "${file_constrained_tree}"
  ensure_dir "${dir_nwkit_download_dir}"
  if [[ ${timetree_constraint} -eq 1 ]]; then
    nwkit mcmctree \
      --download_dir "${dir_nwkit_download_dir}" \
      --infile "${file_undated_species_tree}" \
      --timetree "ci" \
      --min_clade_prop 0.2 \
      --outfile "tmp.constrained.tree.nwk"
    if [[ -s "tmp.constrained.tree.nwk" ]]; then
      mv_out "tmp.constrained.tree.nwk" "${file_constrained_tree}"
    fi
  else
    if [[ ${#mcmctree_divergence_time_constraints[@]} -eq 0 ]]; then
      echo "timetree_constraint=0 requires mcmctree_divergence_time_constraints_str with one or more records."
      echo "Example: Arabidopsis_thaliana,Oryza_sativa,130,-|Arabidopsis_thaliana,Amborella_trichopoda,150,200"
      exit 1
    fi
    tree_string=$(< "${file_undated_species_tree}")
    for mdtc in "${mcmctree_divergence_time_constraints[@]}"; do
      echo "applying ${mdtc}"
      mcmctree_params=()
      if ! parse_mcmctree_constraint_record "${mdtc}" mcmctree_params; then
        exit 1
      fi
      nwkit_args=(
        --download_dir "${dir_nwkit_download_dir}"
        --left_species "${mcmctree_params[0]}"
        --right_species "${mcmctree_params[1]}"
      )
      if [[ "${mcmctree_params[2]}" != "-" ]]; then
        nwkit_args+=(--lower_bound "${mcmctree_params[2]}")
      fi
      if [[ "${mcmctree_params[3]}" != "-" ]]; then
        nwkit_args+=(--upper_bound "${mcmctree_params[3]}")
      fi
      echo "nwkit mcmctree params: ${nwkit_args[*]}"
      tree_string=$(printf '%s\n' "${tree_string}" | nwkit mcmctree "${nwkit_args[@]}")
    done
    printf '%s\n' "${tree_string}" > "tmp.constrained.tree.nwk"
    if [[ -s "tmp.constrained.tree.nwk" ]]; then
      mv_out "tmp.constrained.tree.nwk" "${file_constrained_tree}"
    fi
  fi
  if [[ -s "${file_constrained_tree}" ]]; then
    if ! normalize_iq2mc_constraint_tree "${file_constrained_tree}"; then
      echo "Error: Failed to normalize the constrained tree for IQ2MC."
      exit 1
    fi
  fi
else
  gg_step_skip "${task}"
fi

task="Constrained range plotting"
disable_if_no_input_file "run_plot_constrained_tree" "${file_constrained_tree}"
if [[ ! -s "${file_plot_constrained_tree}" && ${run_plot_constrained_tree} -eq 1 ]]; then
  gg_step_start "${task}"
  Rscript "${gg_support_dir}/plot_constrained_tree.r" \
    --infile="${file_constrained_tree}" \
    --outfile="tmp.constrained_tree_plot.pdf"
  if [[ -s "tmp.constrained_tree_plot.pdf" ]]; then
    mv_out "tmp.constrained_tree_plot.pdf" "${file_plot_constrained_tree}"
  fi
  if [[ -s "${file_plot_constrained_tree}" ]]; then
    echo "Output file found for the task: ${task}"
  fi
else
  gg_step_skip "${task}"
fi

task="IQ2MC step 2 (IQ-TREE Hessian/control generation)"
disable_if_no_input_file "run_mcmctree1" "${file_constrained_tree}" "${file_concat_cds}"
if [[ (! -s "${file_iq2mc_ctl}" || ! -s "${file_iq2mc_hessian}" || ! -s "${file_iq2mc_rooted_tree}" || ! -s "${file_iq2mc_dummy_phy}") && ${run_mcmctree1} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "$(dirname "${file_iq2mc_ctl}")"

  if ! normalize_iq2mc_constraint_tree "${file_constrained_tree}"; then
    echo "Error: Failed to normalize the constrained tree for IQ2MC."
    exit 1
  fi
  if ! clear_directory_contents_safe "$(dirname "${file_iq2mc_ctl}")"; then
    echo "Error: Failed to clear IQ2MC working directory safely."
    exit 1
  fi
  cd "$(dirname "${file_iq2mc_ctl}")"
  seqkit seq --threads 1 "${file_concat_cds}" --out-file "./tmp.iq2mc.concat.cds.fasta"

  if ! "${iq2mc_binary}" \
    -s "./tmp.iq2mc.concat.cds.fasta" \
    -m "${nucleotide_model}" \
    -te "${file_constrained_tree}" \
    --dating mcmctree \
    --mcmc-bds "${mcmc_birth_death_sampling}" \
    --mcmc-clock "${mcmc_clock_model}" \
    --mcmc-iter "${mcmc_burnin},${mcmc_sampfreq},${mcmc_nsample}" \
    -T "${GG_TASK_CPUS}" \
    --prefix iq2mc; then
    echo "Error: IQ2MC step 2 failed. Deleting generated files."
    rm -f -- "${file_iq2mc_prefix}".*
  elif [[ ! -s "${file_iq2mc_ctl}" || ! -s "${file_iq2mc_hessian}" || ! -s "${file_iq2mc_rooted_tree}" || ! -s "${file_iq2mc_dummy_phy}" ]]; then
    echo "Error: IQ2MC step 2 did not generate all expected files."
    rm -f -- "${file_iq2mc_prefix}".*
  fi
  rm -f -- "./tmp.iq2mc.concat.cds.fasta"
  cd "${dir_tmp}"
else
  gg_step_skip "${task}"
fi

task="IQ2MC step 3 (MCMCtree dating run)"
disable_if_no_input_file "run_mcmctree2" "${file_iq2mc_ctl}" "${file_iq2mc_hessian}" "${file_iq2mc_rooted_tree}" "${file_iq2mc_dummy_phy}"
if [[ ! -s "${file_mcmctree_raw_output}" && ${run_mcmctree2} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_mcmctree2}"

  if ! clear_directory_contents_safe "${dir_mcmctree2}"; then
    echo "Error: Failed to clear MCMCtree working directory safely."
    exit 1
  fi
  cd "${dir_mcmctree2}"
  cp_out "${file_iq2mc_ctl}" ./
  cp_out "${file_iq2mc_hessian}" ./
  cp_out "${file_iq2mc_rooted_tree}" ./
  cp_out "${file_iq2mc_dummy_phy}" ./

  # Ensure MCMCtree emits CI-rich summaries (including 95% HPD annotations in FigTree output).
  ctl_basename="$(basename "${file_iq2mc_ctl}")"
  sed -i -e "s/FossilErrprint/FossilErr\\nprint/g" "${ctl_basename}"
  if ! grep -q "^print[[:space:]]*=" "${ctl_basename}"; then
    echo "print = 1" >> "${ctl_basename}"
  fi

  if ! mcmctree "${ctl_basename}"; then
    echo "Error: IQ2MC step 3 failed."
    rm -f -- "${file_mcmctree_raw_output}"
  fi
  cd "${dir_tmp}"
else
  gg_step_skip "${task}"
fi

if [[ ! -s "${file_mcmctree_figtree_tre}" && -s "${file_mcmctree_raw_output}" ]]; then
  echo "Generating ${file_mcmctree_figtree_tre} from ${file_mcmctree_raw_output}"
  awk '
  /Species tree for FigTree/ {print; in_figtree=1; next}
  in_figtree && /^\(\(/ {print; count++; if (count >= 3) exit}
  ' "${file_mcmctree_raw_output}" > "tmp.mcmctree2.txt"
  if [[ -s "tmp.mcmctree2.txt" ]]; then
    mv_out "tmp.mcmctree2.txt" "${file_mcmctree_figtree_tre}"
  fi
  if [[ ! -s "${file_mcmctree_figtree_tre}" ]]; then
    echo "Warning: Failed to extract FigTree content from ${file_mcmctree_raw_output}. Copying raw file instead."
    cp_out "${file_mcmctree_raw_output}" "${file_mcmctree_figtree_tre}"
  fi
fi

task="Convert tree format"
disable_if_no_input_file "run_convert_tree_format" "${file_mcmctree_figtree_tre}"
if [[ ! -s "${file_mcmctree_dated_nwk}" && ${run_convert_tree_format} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_parent_dir "${file_mcmctree_dated_nwk}"

  if grep -q -e "UTREE" "${file_mcmctree_figtree_tre}"; then
    grep -e "UTREE" "${file_mcmctree_figtree_tre}" |
      sed -e "s/.*UTREE 1 = //" -e "s/;.*/;/" \
        > "${dir_mcmctree2}/mcmctree_95CI.nwk"

    grep -e "UTREE" "${file_mcmctree_figtree_tre}" |
      sed -e "s/.*UTREE 1 = //" -e "s/;.*/;/" -e "s/[[:space:]]*\[&95%={[0-9.]*,[[:space:]][0-9.]*}\][[:space:]]*//g" -e "s/:[[:space:]]/:/g" \
        > "${dir_mcmctree2}/mcmctree_no95CI.nwk"
  else
    tree_line="$(awk '/^\(\(/ {line=$0} END {print line}' "${file_mcmctree_figtree_tre}")"
    if [[ -n "${tree_line}" ]]; then
      echo "${tree_line}" > "${dir_mcmctree2}/mcmctree_95CI.nwk"
      echo "${tree_line}" |
        sed -e "s/[[:space:]]*\[&95%={[0-9.]*,[[:space:]][0-9.]*}\][[:space:]]*//g" -e "s/:[[:space:]]/:/g" \
          > "${dir_mcmctree2}/mcmctree_no95CI.nwk"
    else
      echo "Error: Failed to detect a tree string in ${file_mcmctree_figtree_tre}"
      rm -f -- "${dir_mcmctree2}/mcmctree_95CI.nwk" "${dir_mcmctree2}/mcmctree_no95CI.nwk"
    fi
  fi

  if [[ -s "${dir_mcmctree2}/mcmctree_no95CI.nwk" ]]; then
    Rscript -e "library(ape); t=read.tree(\"${dir_mcmctree2}/mcmctree_no95CI.nwk\"); \
    t[['node.label']]=paste0('s',1:(length(t[['tip.label']])-1)); \
    write.tree(t, \"${file_mcmctree_dated_nwk}\")"
  else
    echo "Error: Missing mcmctree_no95CI.nwk. Skipping tree conversion."
  fi

  if [[ ! -s "${file_dated_species_tree}" ]]; then
    echo "Dated species tree for gg_pipeline is not placed yet: ${file_dated_species_tree}"
    if [[ -s "${file_mcmctree_dated_nwk}" ]]; then
      echo "Copying from: ${file_mcmctree_dated_nwk}"
      echo "Copying to: ${file_dated_species_tree}"
      cp_out "${file_mcmctree_dated_nwk}" "${file_dated_species_tree}"
      if [[ -s "${file_plot_mcmctree_pdf}" ]]; then
        echo "Copying from: ${file_plot_mcmctree_pdf}"
        echo "Copying to: ${file_dated_species_tree_pdf}"
        cp_out "${file_plot_mcmctree_pdf}" "${file_dated_species_tree_pdf}"
      fi
      echo "Please manually check whether the species tree is valid."
    fi
  else
    echo "Dated species tree for gg_pipeline is already placed: ${file_dated_species_tree}"
    echo "If necessary, please replace the file with: ${file_mcmctree_dated_nwk}"
  fi
else
  gg_step_skip "${task}"
fi

task="Dated species tree plotting"
disable_if_no_input_file "run_plot_mcmctreer" "${file_mcmctree_dated_nwk}"
if [[ ! -s "${file_plot_mcmctree_pdf}" && ${run_plot_mcmctreer} -eq 1 ]]; then
  gg_step_start "${task}"

  Rscript "${gg_support_dir}/plot_mcmctreer.r" \
    --infile="${file_mcmctree_dated_nwk}" \
    --outfile="tmp.plot_mcmctreer.pdf"
  if [[ -s "tmp.plot_mcmctreer.pdf" ]]; then
    mv_out "tmp.plot_mcmctreer.pdf" "${file_plot_mcmctree_pdf}"
  fi

  if [[ -s "${file_plot_mcmctree_pdf}" ]]; then
    echo "Output file found for the task: ${task}"
    echo "Output file: ${file_plot_mcmctree_pdf}"
    if [[ -s "${file_dated_species_tree}" ]]; then
      echo "Copying from: ${file_plot_mcmctree_pdf}"
      echo "Copying to: ${file_dated_species_tree_pdf}"
      cp_out "${file_plot_mcmctree_pdf}" "${file_dated_species_tree_pdf}"
    fi
  fi
else
  gg_step_skip "${task}"
fi

remove_empty_subdirs "${dir_species_tree}"
if [[ ${delete_tmp_dir} -eq 1 ]]; then
  if [[ -d "${dir_tmp}" ]]; then
    echo "Removing tmp directory: ${dir_tmp}"
    rm -rf -- "${dir_tmp}"
  fi
fi

# Orthogroup inference

orthofinder_output_directory_cleanup() {
  local target_dir=$1
  local _threads=${2:-1}
  if [[ -d "${target_dir}" ]]; then
    remove_empty_subdirs "${target_dir}"
  fi
}

ensure_dir "${dir_tmp}"
cd "${dir_tmp}"

task="OrthoFinder"
if [[ ! -s "${file_orthofinder_done_marker}" && ${run_orthofinder} -eq 1 ]]; then
  dir_sp_protein_orthofinder="${dir_sp_protein}_orthofinder"
  gg_step_start "${task}"
  prepare_species_protein_tmp
  prepare_species_protein_orthofinder_dir "${dir_sp_protein}" "${dir_sp_protein_orthofinder}"
  ensure_dir "${dir_orthofinder}"
  ensure_dir "${dir_orthofinder_hog2og}"

  orthofinder_algorithm_threads=$((${GG_TASK_CPUS} / 4))
  if [[ ${orthofinder_algorithm_threads} -eq 0 ]]; then
    orthofinder_algorithm_threads=1
  fi
  param_species_tree=()
  species_tree=""
  if [[ -s "${dir_species_tree_summary}/dated_species_tree.nwk" ]]; then
    species_tree="${dir_species_tree_summary}/dated_species_tree.nwk"
  elif [[ -s "${dir_species_tree_summary}/undated_species_tree.nwk" ]]; then
    species_tree="${dir_species_tree_summary}/undated_species_tree.nwk"
  fi
  if [[ -n "${species_tree}" && -s "${species_tree}" ]]; then
    echo "OrthoFinder will use the species tree: ${species_tree}"
    param_species_tree=(-s "${species_tree}")
  fi
  echo "OrthoFinder will use ${GG_TASK_CPUS} threads for diamond search."
  echo "OrthoFinder will use ${orthofinder_algorithm_threads} threads for the OrthoFinder algorithm."

  protein_files=()
  mapfile -t protein_files < <(find "${dir_sp_protein_orthofinder}" -maxdepth 1 -type f ! -name '.*' \( -name "*.fa" -o -name "*.fa.gz" -o -name "*.fas" -o -name "*.fas.gz" -o -name "*.fasta" -o -name "*.fasta.gz" -o -name "*.fna" -o -name "*.fna.gz" \) | sort)
  num_sp=${#protein_files[@]}
  if [[ ${num_sp} -eq 0 ]]; then
    echo "No protein FASTA files were found in: ${dir_sp_protein_orthofinder}. Exiting."
    exit 1
  fi
  if [[ ${#param_species_tree[@]} -gt 0 ]]; then
    species_ids=()
    for protein_file in "${protein_files[@]}"; do
      species_base=$(basename "${protein_file}")
      species_base=${species_base%.gz}
      species_base=${species_base%.*}
      species_ids+=("${species_base}")
    done
    mapfile -t missing_species < <(
      python - "${species_tree}" "${species_ids[@]}" << 'PY'
import re
import sys

tree_file = sys.argv[1]
species = sys.argv[2:]
with open(tree_file) as f:
    tree = f.read().strip()

leaves = set(re.findall(r'(?<=[(,])([^:(),]+?)(?=[:),])', tree))
missing = [sp for sp in species if sp not in leaves]
for sp in missing:
    print(sp)
PY
    )
    if [[ ${#missing_species[@]} -gt 0 ]]; then
      echo "Species tree is missing ${#missing_species[@]} species: ${missing_species[*]}"
      echo "Running OrthoFinder without species tree constraints."
      param_species_tree=()
    fi
  fi
  if [[ ${num_sp} -gt ${max_orthofinder_core_species} ]]; then
    echo "The number of species (${num_sp}) is greater than the maximum number of core species (${max_orthofinder_core_species}) for OrthoFinder."
    echo "OrthoFinder will be run for 2 rounds (using --assign): For details, see https://github.com/davidemms/OrthoFinder"

    species_cds_core=()
    mapfile -t species_cds_core < <(printf "%s\n" "${protein_files[@]##*/}" | shuf -n "${max_orthofinder_core_species}")
    echo "Core CDS files: ${species_cds_core[@]}"
    if [[ -e "${dir_sp_protein}_core" ]]; then
      rm -rf -- "${dir_sp_protein}_core"
    fi
    mkdir -p "${dir_sp_protein}_core"
    for sp_cds_core in "${species_cds_core[@]}"; do
      cp_out "${dir_sp_protein_orthofinder}"/"${sp_cds_core}" "${dir_sp_protein}"_core
    done

    species_cds_additional=()
    mapfile -t species_cds_additional < <(printf "%s\n" "${protein_files[@]##*/}" | grep -v -F -x -f <(printf "%s\n" "${species_cds_core[@]}"))
    echo "Additional CDS files: ${species_cds_additional[@]}"
    if [[ -e "${dir_sp_protein}_additional" ]]; then
      rm -rf -- "${dir_sp_protein}_additional"
    fi
    mkdir -p "${dir_sp_protein}_additional"
    for sp_cds_additional in "${species_cds_additional[@]}"; do
      cp_out "${dir_sp_protein_orthofinder}"/"${sp_cds_additional}" "${dir_sp_protein}"_additional
    done

    if [[ -e "${dir_orthofinder}"/core ]]; then
      rm -rf -- "${dir_orthofinder}/core"
    fi

    if [[ -e "${dir_orthofinder}/species_tree_core.nwk" ]]; then
      rm -f -- "${dir_orthofinder}/species_tree_core.nwk"
    fi
    core_species_names=()
    mapfile -t core_species_names < <(
      find "${dir_sp_protein}_core" -maxdepth 1 -type f ! -name '.*' | sort |
        awk '{name=$0; sub(/^.*\//, "", name); sub(/\.[^.]*$/, "", name); print name}'
    )
    core_species_regex=$(printf '%s|' "${core_species_names[@]}")
    core_species_regex=${core_species_regex%|}
    if [[ -z "${core_species_regex}" ]]; then
      echo "Failed to build core species regex. No core protein files were detected. Exiting."
      exit 1
    fi
    nwkit prune --invert_match yes --pattern "${core_species_regex}" --infile "${species_tree}" --outfile "${dir_orthofinder}/species_tree_core.nwk"

    if orthofinder \
      -t "${GG_TASK_CPUS}" \
      -a "${orthofinder_algorithm_threads}" \
      -M "msa" \
      -S "diamond" \
      -f "${dir_sp_protein}_core" \
      -n "core" \
      -o "${dir_orthofinder}/core" \
      -s "${dir_orthofinder}/species_tree_core.nwk"; then
      orthofinder_core_exit_code=0
    else
      orthofinder_core_exit_code=$?
    fi

    if [[ ${orthofinder_core_exit_code} -ne 0 ]]; then
      echo "OrthoFinder failed in the core-species run. Exiting."
      exit 1
    fi

    if orthofinder \
      -t "${GG_TASK_CPUS}" \
      -a "${orthofinder_algorithm_threads}" \
      -M "msa" \
      -S "diamond" \
      -n "all" \
      --core "${dir_orthofinder}/core/Results_core" \
      --assign "${dir_sp_protein}_additional" \
      "${param_species_tree[@]}"; then
      orthofinder_main_exit_code=0
    else
      orthofinder_main_exit_code=$?
    fi

    if [[ ${orthofinder_main_exit_code} -ne 0 ]]; then
      echo "OrthoFinder failed in the all-species run. Exiting."
      exit 1
    fi

    shopt -s nullglob
    orthofinder_all_outputs=("${dir_orthofinder}"/core/Results_all/*)
    orthofinder_core_outputs=("${dir_orthofinder}"/core/Results_core/*)
    shopt -u nullglob
    if [[ ${#orthofinder_all_outputs[@]} -eq 0 || ${#orthofinder_core_outputs[@]} -eq 0 ]]; then
      echo "OrthoFinder core/main output files were expected but not found after completion."
      exit 1
    fi
    mv_out "${orthofinder_all_outputs[@]}" "${dir_orthofinder}"
    mv_out "${orthofinder_core_outputs[@]}" "${dir_orthofinder}/core"
    shopt -s nullglob
    orthofinder_result_dirs=("${dir_orthofinder}/core/Results_"*)
    shopt -u nullglob
    if [[ ${#orthofinder_result_dirs[@]} -gt 0 ]]; then
      rm -rf -- "${orthofinder_result_dirs[@]}"
    fi
    orthofinder_output_directory_cleanup "${dir_orthofinder}/core" "${GG_TASK_CPUS}"
  else
    echo "The number of species (${num_sp}) is less than or equal to the maximum number of core species (${max_orthofinder_core_species}) for OrthoFinder."
    echo "OrthoFinder will be run for 1 round."

    if orthofinder \
      -t "${GG_TASK_CPUS}" \
      -a "${orthofinder_algorithm_threads}" \
      -M "msa" \
      -S "diamond" \
      -f "${dir_sp_protein_orthofinder}" \
      -n "main" \
      -o "${dir_orthofinder}/main" \
      "${param_species_tree[@]}"; then
      orthofinder_main_exit_code=0
    else
      orthofinder_main_exit_code=$?
    fi

    if [[ ${orthofinder_main_exit_code} -ne 0 ]]; then
      echo "OrthoFinder failed in the all-species run. Exiting."
      exit 1
    fi

    shopt -s nullglob
    orthofinder_main_outputs=("${dir_orthofinder}"/main/Results_main/*)
    shopt -u nullglob
    if [[ ${#orthofinder_main_outputs[@]} -eq 0 ]]; then
      echo "OrthoFinder main output files were expected but not found after completion."
      exit 1
    fi
    mv_out "${orthofinder_main_outputs[@]}" "${dir_orthofinder}"
    rm -rf -- "${dir_orthofinder}/main"
  fi

  echo "OrthoFinder finished successfully."
  orthofinder_output_directory_cleanup "${dir_orthofinder}" "${GG_TASK_CPUS}"

  hog_table="${dir_orthofinder}/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"
  if [[ ! -s "${hog_table}" ]]; then
    hog_candidates=()
    mapfile -t hog_candidates < <(find "${dir_orthofinder}/Phylogenetic_Hierarchical_Orthogroups" -maxdepth 1 -type f -name 'N*.tsv' | sort -V)
    if [[ ${#hog_candidates[@]} -eq 0 ]]; then
      echo "No HOG tables were found in ${dir_orthofinder}/Phylogenetic_Hierarchical_Orthogroups. Exiting."
      exit 1
    fi
    hog_table="${hog_candidates[0]}"
    echo "N0.tsv was not found. Falling back to HOG table: ${hog_table}"
  fi

  python "${gg_support_dir}/orthogroup_table_formatter.py" \
    --file_orthogroup_table "${hog_table}" \
    --dir_out "${dir_orthofinder_hog2og}" \
    --mode "hog2og"
else
  gg_step_skip "${task}"
fi

task="OMArk analysis of species-wise protein input files"
run_shared_species_omark_stage

task="Summarizing OMArk species quality results"
run_shared_omark_summary_stage

task="Selecting orthogroups based on gene and species numbers"
if [[ ! -s "${file_orthogroup_selection}" && ${run_og_selection} -eq 1 ]]; then
  gg_step_start "${task}"
  prepare_species_protein_tmp
  ensure_dir "${dir_orthofinder_filtered}"
  if [[ "${orthogroup_annotation_method}" == "blastp" ]]; then
    if ! uniprot_db_prefix=$(ensure_uniprot_sprot_blast_db "${gg_workspace_dir}"); then
      echo "Failed to prepare UniProt Swiss-Prot BLASTP DB. Exiting."
      exit 1
    fi
  else
    if ! uniprot_db_prefix=$(ensure_uniprot_sprot_mmseqs_db "${gg_workspace_dir}"); then
      echo "Failed to prepare UniProt Swiss-Prot MMseqs2 DB. Exiting."
      exit 1
    fi
  fi

  if [[ ${orthogroup_table} == "OG" ]]; then
    dir_orthogroup_selection_input="${dir_orthofinder_og}"
  elif [[ ${orthogroup_table} == "HOG" ]]; then
    dir_orthogroup_selection_input="${dir_orthofinder_hog2og}"
  else
    echo "Unsupported orthogroup_table: ${orthogroup_table}. Allowed values are OG or HOG."
    exit 1
  fi

  if [[ ! -s "${dir_orthogroup_selection_input}/Orthogroups.tsv" || ! -s "${dir_orthogroup_selection_input}/Orthogroups.GeneCount.tsv" ]]; then
    echo "Orthogroup source files were not found in ${dir_orthogroup_selection_input}. Exiting."
    exit 1
  fi
  cp_out "${dir_orthogroup_selection_input}/Orthogroups.tsv" "${dir_orthofinder_filtered}/Orthogroups.tsv"
  cp_out "${dir_orthogroup_selection_input}/Orthogroups.GeneCount.tsv" "${dir_orthofinder_filtered}/Orthogroups.GeneCount.tsv"

  if python "${gg_support_dir}/orthogroup_selection.py" \
    --dir_orthofinder_og "${dir_orthofinder_filtered}" \
    --dir_species_protein "${dir_sp_protein}" \
    --min_gene_num "${min_num_gene}" \
    --max_gene_num "${max_num_gene}" \
    --min_species_num "${min_num_species}" \
    --min_percent_species_coverage "${min_percent_species_coverage}" \
    --remove_unannotated 'yes' \
    --gene_size_quantiles '0.05,0.25,0.5,0.75,0.95' \
    --annotation_search_method "${orthogroup_annotation_method}" \
    --path_search_db "${uniprot_db_prefix}" \
    --evalue '1e-2' \
    --ncpu "${GG_TASK_CPUS}"; then
    exit_code=0
  else
    exit_code=$?
  fi

  if [[ ${exit_code} -eq 0 ]]; then
    echo "Orthogroup selection finished successfully."
  else
    echo "Orthogroup selection failed. Exiting."
    exit 1
  fi
else
  gg_step_skip "${task}"
fi

task="Orthogroup method comparison"
disable_if_no_input_file "run_orthogroup_method_comparison" "${file_orthofinder_done_marker}"
if [[ ! -s "${file_orthogroup_method_comparison}" && ${run_orthogroup_method_comparison} -eq 1 ]]; then
  gg_step_start "${task}"

  if python "${gg_support_dir}/orthogroup_method_comparison.py" \
    --orthofinder_og_genecount "${dir_orthofinder_og}/Orthogroups.GeneCount.tsv" \
    --orthofinder_hog_genecount "${dir_orthofinder_hog2og}/Orthogroups.GeneCount.tsv"; then
    exit_code=0
  else
    exit_code=$?
  fi
  if [[ ${exit_code} -eq 0 ]]; then
    echo "Orthogroup method comparison finished successfully."
    mv_out orthogroup_histogram.pdf "${file_orthogroup_method_comparison}"
  else
    echo "Orthogroup method comparison failed. Exiting."
    exit 1
  fi
else
  gg_step_skip "${task}"
fi

cleanup_species_protein_tmp
trap - EXIT

# Genome evolution
if [[ ${run_busco_notung_root_dna} -eq 1 || ${run_busco_notung_root_pep} -eq 1 || ${run_busco_root_dna} -eq 1 || ${run_busco_root_pep} -eq 1 || ${run_busco_grampa_dna} -eq 1 || ${run_busco_grampa_pep} -eq 1 ]]; then
  if [[ ! -s "${notung_jar}" ]]; then
    echo "Notung jar was not found: ${notung_jar}"
    echo "Disabling NOTUNG-dependent BUSCO tasks: run_busco_notung_root_dna, run_busco_notung_root_pep, run_busco_root_dna, run_busco_root_pep, run_busco_grampa_dna, run_busco_grampa_pep"
    run_busco_notung_root_dna=0
    run_busco_notung_root_pep=0
    run_busco_root_dna=0
    run_busco_root_pep=0
    run_busco_grampa_dna=0
    run_busco_grampa_pep=0
  fi
fi

if [[ ${run_orthogroup_grampa} -eq 1 ]]; then
  if [[ ! -s "${file_orthogroup_genecount_selected}" ]]; then
    echo "Disabling run_orthogroup_grampa because required file is missing: ${file_orthogroup_genecount_selected}"
    run_orthogroup_grampa=0
  elif [[ ! -d "${dir_og_rooted_tree}" ]]; then
    echo "Disabling run_orthogroup_grampa because required directory is missing: ${dir_og_rooted_tree}"
    run_orthogroup_grampa=0
  fi
fi

if [[ -z "${grampa_h1}" ]]; then
  if [[ ${run_busco_grampa_dna} -eq 1 || ${run_busco_grampa_pep} -eq 1 || ${run_orthogroup_grampa} -eq 1 ]]; then
    echo "Disabling GRAMPA tasks because grampa_h1 is empty. Set grampa_h1 in gg_genome_evolution_entrypoint.sh to enable them."
  fi
  run_busco_grampa_dna=0
  run_busco_grampa_pep=0
  run_orthogroup_grampa=0
fi

if [[ -z "${target_branch_go}" && ${run_go_enrichment} -eq 1 ]]; then
  echo "Disabling run_go_enrichment because target_branch_go is empty. Set target_branch_go in gg_genome_evolution_entrypoint.sh to enable it."
  run_go_enrichment=0
fi

ensure_dir "${dir_tmp}"
cd "${dir_tmp}"

task="Generating fasta files for individual single-copy genes"
sync_genome_busco_summary_table_from_shared || true
disable_if_no_input_file "run_busco_getfasta" "${file_genome_busco_summary_table}"
if [[ ${run_busco_getfasta} -eq 1 ]]; then
  prepare_species_tree_input_dir
  gg_step_start "${task}"
  ensure_dir "${dir_busco_fasta}"
  busco_rows=()
  mapfile -t busco_rows < <(tail -n +2 "${file_genome_busco_summary_table}")
  num_busco_ids=${#busco_rows[@]}

  generate_single_copy_fasta() {
    local busco_idx=$1
    local busco_row="${busco_rows[${busco_idx}]}"
    local busco_id
    local outfile2
    local cols=()
    local genes=()
    local genes2=()
    IFS=$'\t' read -r -a cols <<< "${busco_row}"
    busco_id="${cols[0]:-}"
    if [[ -z "${busco_id}" ]]; then
      return 0
    fi
    if [[ ${#cols[@]} -gt 3 ]]; then
      genes=("${cols[@]:3}")
    fi
    outfile2="${dir_busco_fasta}/${busco_id}${genome_busco_fasta_suffix}"
    if [[ -s "${outfile2}" ]]; then
      return 0
    fi
    echo "busco_id: ${busco_id}"
    if [[ ! -s "${outfile2}" ]]; then
      mapfile -t genes2 < <(gg_busco_gene_tokens "split_duplicated" "${genes[@]}")
      if [[ ${#genes2[@]} -eq 0 ]]; then
        echo "Skipping. ${busco_id} has no genes in the selected species."
        return 0
      fi
      if [[ "${input_sequence_mode}" == "protein" ]]; then
        gg_seqkit_grep_by_patterns_from_infile_list 1 "species_tree_input_fasta_list.txt" "${genes2[@]}" |
          seqkit replace --pattern " .*" --replacement "" --ignore-case --threads 1 |
          seqkit seq --threads 1 --out-file "${outfile2}"
      else
        gg_seqkit_grep_by_patterns_from_infile_list 1 "species_tree_input_fasta_list.txt" "${genes2[@]}" |
          gg_prepare_cds_fasta_stream 1 |
          seqkit seq --threads 1 --out-file "${outfile2}"
      fi
      if [[ ! -s "${outfile2}" ]]; then
        echo "File is empty. Removing: ${outfile2}"
        rm -f -- "${outfile2}"
      fi
    fi
  }

  gg_find_fasta_files "${species_tree_input_dir}" 1 > species_tree_input_fasta_list.txt
  for ((busco_idx = 0; busco_idx < num_busco_ids; busco_idx++)); do
    wait_until_jobn_le ${GG_TASK_CPUS}
    generate_single_copy_fasta "${busco_idx}" &
  done
  wait_for_background_jobs
  rm -f -- species_tree_input_fasta_list.txt
  find . -maxdepth 1 -type f -name 'tmp.*' -delete
else
  gg_step_skip "${task}"
fi

task="In-frame mafft alignment of duplicate-containing BUSCO genes"
if [[ ${run_busco_mafft} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_busco_mafft}"

  run_mafft() {
    infile=$1
    infile_base=${infile%%.*}
    outfile=${dir_busco_mafft}/${infile_base}${genome_busco_aln_suffix}
    if [[ -s "${outfile}" ]]; then
      return 0
    fi
    echo "$(date): start mafft: ${infile_base}"

    if [[ "${input_sequence_mode}" == "protein" ]]; then
      seqkit seq --threads 1 "${dir_busco_fasta}/${infile}" --out-file "tmp.${infile_base}.input.pep.fasta"
      num_seq=$(gg_count_fasta_records "tmp.${infile_base}.input.pep.fasta")
      if [[ ${num_seq} -lt 2 ]]; then
        echo "Skipped MAFFT because fewer than 2 sequences were found: ${infile}"
        seqkit seq --threads 1 "tmp.${infile_base}.input.pep.fasta" --out-file "tmp.${infile_base}.pep.out.fa.gz"
        mv_out "tmp.${infile_base}.pep.out.fa.gz" "${outfile}"
        rm -f -- "tmp.${infile_base}"*
        return 0
      fi
      mafft \
        --auto \
        --amino \
        --thread 1 \
        "tmp.${infile_base}.input.pep.fasta" \
        > "tmp.${infile_base}.pep.aln.fasta"
      if [[ -s "tmp.${infile_base}.pep.aln.fasta" ]]; then
        seqkit seq --threads 1 "tmp.${infile_base}.pep.aln.fasta" --out-file "tmp.${infile_base}.pep.aln.out.fa.gz"
        mv_out "tmp.${infile_base}.pep.aln.out.fa.gz" "${outfile}"
      fi
    else
      seqkit seq --threads 1 "${dir_busco_fasta}/${infile}" --out-file "tmp.${infile_base}.input.cds.fasta"
      cdskit mask \
        --seqfile "tmp.${infile_base}.input.cds.fasta" \
        --outfile "tmp.${infile_base}.cds.fasta"

      num_seq=$(gg_count_fasta_records "tmp.${infile_base}.cds.fasta")
      if [[ ${num_seq} -lt 2 ]]; then
        echo "Skipped MAFFT/backalign because fewer than 2 sequences were found: ${infile}"
        seqkit seq --threads 1 "tmp.${infile_base}.cds.fasta" --out-file "tmp.${infile_base}.cds.out.fa.gz"
        mv_out "tmp.${infile_base}.cds.out.fa.gz" "${outfile}"
        rm -f -- "tmp.${infile_base}"*
        return 0
      fi

      seqkit translate \
        --allow-unknown-codon \
        --transl-table "${genetic_code}" \
        --threads 1 \
        "tmp.${infile_base}.cds.fasta" \
        > "tmp.${infile_base}.pep.fasta"

      mafft \
        --auto \
        --amino \
        --thread 1 \
        "tmp.${infile_base}.pep.fasta" \
        > "tmp.${infile_base}.pep.aln.fasta"

      cdskit backalign \
        --seqfile "tmp.${infile_base}.cds.fasta" \
        --aa_aln "tmp.${infile_base}.pep.aln.fasta" \
        --codontable "${genetic_code}" \
        --outfile "tmp.${infile_base}.cds.aln.fasta"

      if [[ -s "tmp.${infile_base}.cds.aln.fasta" ]]; then
        seqkit seq --threads 1 "tmp.${infile_base}.cds.aln.fasta" --out-file "tmp.${infile_base}.cds.aln.out.fa.gz"
        mv_out "tmp.${infile_base}.cds.aln.out.fa.gz" "${outfile}"
      fi
    fi
    rm -f -- "tmp.${infile_base}"*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_busco_fasta}" "${genome_busco_fasta_glob}")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    run_mafft "${input_alignment_file}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

task="TrimAl of duplicate-containing BUSCO genes"
if [[ ${run_busco_trimal} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_busco_trimal}"

  run_trimal() {
    infile=$1
    infile_base=${infile%%.*}
    outfile="${dir_busco_trimal}/${infile_base}${genome_busco_trimal_suffix}"
    if [[ -s "${outfile}" ]]; then
      return 0
    fi
    if [[ "${input_sequence_mode}" == "protein" ]]; then
      seqkit seq --threads 1 "${dir_busco_mafft}/${infile}" --out-file "tmp.${infile_base}.pep.aln.fasta"
      trimal \
        -in "tmp.${infile_base}.pep.aln.fasta" \
        -out "tmp.${infile_base}.trimal.fasta" \
        -automated1
    else
      seqkit seq --remove-gaps --threads 1 "${dir_busco_mafft}/${infile}" > "tmp.${infile_base}.degap.fasta"
      seqkit translate --transl-table "${genetic_code}" --threads 1 "${dir_busco_mafft}/${infile}" > "tmp.${infile_base}.pep.fasta"

      trimal \
        -in "tmp.${infile_base}.pep.fasta" \
        -backtrans "tmp.${infile_base}.degap.fasta" \
        -out "tmp.${infile_base}.trimal.fasta" \
        -ignorestopcodon \
        -automated1
    fi

    if [[ -s "tmp.${infile_base}.trimal.fasta" ]]; then
      seqkit seq --threads 1 "tmp.${infile_base}.trimal.fasta" --out-file "tmp.${infile_base}.trimal.out.fa.gz"
      mv_out "tmp.${infile_base}.trimal.out.fa.gz" "${outfile}"
    fi
    rm -f -- "tmp.${infile_base}."*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_busco_mafft}" "${genome_busco_aln_glob}")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    run_trimal "${input_alignment_file}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

task="IQ-TREE for duplicate-containing BUSCO DNA trees"
if [[ ${run_busco_iqtree_dna} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_busco_iqtree_dna}"

  busco_iqtree_dna() {
    infile=$1
    indir=$2
    outdir=$3
    infile_base=${infile%%.*}
    outfile="${outdir}/${infile_base}.busco.nwk"
    if [[ -s "${outfile}" ]]; then
      return 0
    fi
    num_seq=$(gg_count_fasta_records "${dir_busco_trimal}/${infile}")
    if [[ ${num_seq} -lt 3 ]]; then
      echo "Skipped. At least 3 sequences are necessary for IQ-TREE: ${infile}"
      return 0
    fi

    seqkit seq --threads 1 "${indir}/${infile}" --out-file "./tmp.${infile_base}.input.fasta"

    iqtree \
      -s "./tmp.${infile_base}.input.fasta" \
      -m "${nucleotide_model}" \
      -T 1 \
      --prefix "tmp.${infile_base}" \
      --seed 12345 \
      --redo

    mv_out tmp."${infile_base}".treefile "${outfile}"
    rm -f -- "tmp.${infile_base}."*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_busco_trimal}" "*.busco.trimal.fa.gz")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    busco_iqtree_dna "${input_alignment_file}" "${dir_busco_trimal}" "${dir_busco_iqtree_dna}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

task="IQ-TREE for duplicate-containing BUSCO protein trees"
if [[ ${run_busco_iqtree_pep} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_busco_iqtree_pep}"

  busco_iqtree_pep() {
    infile=$1
    indir=$2
    outdir=$3
    infile_base=${infile%%.*}
    outfile="${outdir}/${infile_base}.busco.nwk"
    if [[ -s "${outfile}" ]]; then
      return 0
    fi
    num_seq=$(gg_count_fasta_records "${dir_busco_trimal}/${infile}")
    if [[ ${num_seq} -lt 3 ]]; then
      echo "Skipped. At least 3 sequences are necessary for IQ-TREE: ${infile}"
      return 0
    fi

    if [[ "${input_sequence_mode}" == "protein" ]]; then
      seqkit seq --threads 1 "${indir}/${infile}" > "./tmp.${infile_base}.input.fasta"
    else
      seqkit translate --transl-table "${genetic_code}" --threads 1 "${indir}/${infile}" > "./tmp.${infile_base}.input.fasta"
    fi

    iqtree \
      -s "./tmp.${infile_base}.input.fasta" \
      -m "${protein_model}" \
      -st AA \
      -T 1 \
      --prefix "tmp.${infile_base}" \
      --seed 12345 \
      --redo

    mv_out tmp."${infile_base}".treefile "${outfile}"
    rm -f -- "tmp.${infile_base}."*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_busco_trimal}" "${genome_busco_trimal_glob}")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    busco_iqtree_pep "${input_alignment_file}" "${dir_busco_trimal}" "${dir_busco_iqtree_pep}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

busco_notung() {
  infile=$1
  indir=$2
  outdir=$3
  busco_id=${infile%%.*}
  outfile="${outdir}/${busco_id}.busco.notung.root.zip"
  if [[ -s "${outfile}" ]]; then
    return 0
  fi
  if [[ -e "./${busco_id}.notung.root" ]]; then
    rm -rf -- "./${busco_id}.notung.root"
  fi
  java -jar -Xmx${memory_notung}g ${notung_jar} \
    -s "${file_dated_species_tree}" \
    -g "${indir}/${infile}" \
    --root \
    --infertransfers "false" \
    --treeoutput newick \
    --log \
    --treestats \
    --events \
    --parsable \
    --speciestag prefix \
    --allopt \
    --maxtrees 1000 \
    --nolosses \
    --outputdir "./${busco_id}.notung.root"
  if [[ -e "${busco_id}.notung.root/${busco_id}.busco.nwk.rooting.0" ]]; then
    zip -rq "${busco_id}.busco.notung.root.zip" "${busco_id}.notung.root"
    mv_out "${busco_id}.busco.notung.root.zip" "${outfile}"
    rm -rf -- "./${busco_id}.notung.root"
  fi
}

task="NOTUNG rooting of duplicate-containing BUSCO DNA trees"
disable_if_no_input_file "run_busco_notung_root_dna" "${file_dated_species_tree}"
if [[ ${run_busco_notung_root_dna} -eq 1 ]]; then
  gg_step_start "${task}"

  infiles=()
  mapfile -t infiles < <(gg_find_file_basenames "${dir_busco_iqtree_dna}")
  for infile in "${infiles[@]}"; do
    wait_until_jobn_le $((${GG_TASK_CPUS} / 2))
    busco_notung "${infile}" "${dir_busco_iqtree_dna}" "${dir_busco_notung_dna}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

task="NOTUNG rooting of duplicate-containing BUSCO protein trees"
disable_if_no_input_file "run_busco_notung_root_pep" "${file_dated_species_tree}"
if [[ ${run_busco_notung_root_pep} -eq 1 ]]; then
  gg_step_start "${task}"

  infiles=()
  mapfile -t infiles < <(gg_find_file_basenames "${dir_busco_iqtree_pep}")
  for infile in "${infiles[@]}"; do
    wait_until_jobn_le $((${GG_TASK_CPUS} / 2))
    busco_notung "${infile}" "${dir_busco_iqtree_pep}" "${dir_busco_notung_pep}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

busco_species_tree_assisted_gene_tree_rooting() {
  infile=$1
  indir=$2
  intreedir=$3
  outdir_txt=$4
  outdir_nwk=$5
  busco_id=${infile%%.*}
  intree="${intreedir}/${busco_id}.busco.nwk"
  outfile_txt="${outdir_txt}/${busco_id}.busco.root.txt"
  outfile_nwk="${outdir_nwk}/${busco_id}.busco.root.nwk"
  if [[ -s "${outfile_txt}" && -s "${outfile_nwk}" ]]; then
    return 0
  fi
  echo "Start NOTUNG root: ${busco_id}"
  if [[ -e "./${busco_id}.notung.root" ]]; then
    rm -rf -- "./${busco_id}.notung.root"
  fi
  cp_out "${indir}"/"${infile}" .
  unzip -q "${infile}"

  Rscript "${gg_support_dir}/species_tree_guided_gene_tree_rooting.r" \
    "--notung_root_zip=./${infile}" \
    "--in_tree=${intree}" \
    "--out_tree=${busco_id}.root.nwk" \
    "--ncpu=${GG_TASK_CPUS}" \
    2>&1 | tee "${busco_id}.root.txt"

  if [[ -s "${busco_id}.root.nwk" ]]; then
    mv_out "${busco_id}".root.txt "${outfile_txt}"
    mv_out "${busco_id}".root.nwk "${outfile_nwk}"
    rm -f -- "${infile}"
    rm -rf -- "${busco_id}.notung.root"
  fi
}

task="Species-tree-guided gene tree rooting of duplicate-containing BUSCO DNA trees"
disable_if_no_input_file "run_busco_root_dna" "${file_dated_species_tree}"
if [[ ${run_busco_root_dna} -eq 1 ]]; then
  gg_step_start "${task}"

  infiles=()
  mapfile -t infiles < <(gg_find_file_basenames "${dir_busco_notung_dna}")
  for infile in "${infiles[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    busco_species_tree_assisted_gene_tree_rooting "${infile}" "${dir_busco_notung_dna}" "${dir_busco_iqtree_dna}" "${dir_busco_rooted_txt_dna}" "${dir_busco_rooted_nwk_dna}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

task="Species-tree-guided gene tree rooting of duplicate-containing BUSCO protein trees"
disable_if_no_input_file "run_busco_root_pep" "${file_dated_species_tree}"
if [[ ${run_busco_root_pep} -eq 1 ]]; then
  gg_step_start "${task}"

  infiles=()
  mapfile -t infiles < <(gg_find_file_basenames "${dir_busco_notung_pep}")
  for infile in "${infiles[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    busco_species_tree_assisted_gene_tree_rooting "${infile}" "${dir_busco_notung_pep}" "${dir_busco_iqtree_pep}" "${dir_busco_rooted_txt_pep}" "${dir_busco_rooted_nwk_pep}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

busco_grampa() {
  indir=$1
  outdir=$2
  outfile=$3
  ensure_dir "${outdir}"
  local nwk_files=()
  shopt -s nullglob
  nwk_files=("${indir}"/*.nwk)
  shopt -u nullglob
  if [[ ${#nwk_files[@]} -eq 0 ]]; then
    echo "Skipping Grampa because no rooted gene trees were found in: ${indir}"
    return 0
  fi
  if [[ -e "./grampa_out" ]]; then
    rm -rf -- "./grampa_out"
  fi
  nwkit drop --infile "${file_dated_species_tree}" --target 'intnode' --name 'yes' |
    sed -e "s/_/-/g" \
      > "grampa_input_species_tree.nwk"

  > "grampa_input_gene_trees.nwk"
  > "busco_genetree_filenames.txt"
  for nwk_file in "${nwk_files[@]}"; do
    transformed_tree=$(
      sed -E "s/([(,])([^_)(,:]+)_([^_)(,:]+)_([^)(,:]+)([)(,:])/\1\4|||\2\-\3\5/g" "${nwk_file}" |
        sed -e "s/_/-/g" -e "s/|||/_/g" |
        sed -E 's/\)([^():;,]+):/\):/g; s/\)([^():;,]+);/\);/g' |
        tr -d '\r\n'
    )
    if [[ "${transformed_tree}" == *"("* && "${transformed_tree}" == *")"* && "${transformed_tree}" == *";"* ]]; then
      printf "%s\n" "${transformed_tree}" >> "grampa_input_gene_trees.nwk"
      printf "%s\n" "${nwk_file##*/}" >> "busco_genetree_filenames.txt"
    fi
  done

  valid_tree_count=$(awk 'NF>0{n++} END{print n+0}' "grampa_input_gene_trees.nwk")
  if [[ ${valid_tree_count} -eq 0 ]]; then
    echo "Skipping Grampa because no valid rooted gene trees were prepared from: ${indir}"
    rm -f -- "grampa_input_species_tree.nwk" "grampa_input_gene_trees.nwk" "busco_genetree_filenames.txt"
    return 0
  fi
  echo "Number of rooted gene trees passed to Grampa: ${valid_tree_count}"

  grampa_args=(
    -s "grampa_input_species_tree.nwk"
    -g "grampa_input_gene_trees.nwk"
    -o "./grampa_out"
    -p "${GG_TASK_CPUS}"
    -v -1
    --maps
  )
  if [[ -n "${grampa_h1}" ]]; then
    local grampa_h1_normalized=${grampa_h1//_/-}
    grampa_h1_normalized=${grampa_h1_normalized//[[:space:]]/-}
    grampa_args+=(-h1 "${grampa_h1_normalized}")
  fi

  grampa.py "${grampa_args[@]}"

  grampa_filtered_file=""
  if [[ -s "./grampa_out/grampa_trees_filtered.txt" ]]; then
    grampa_filtered_file="./grampa_out/grampa_trees_filtered.txt"
  elif [[ -s "./grampa_out/grampa-trees-filtered.txt" ]]; then
    grampa_filtered_file="./grampa_out/grampa-trees-filtered.txt"
  fi

  if [[ -n "${grampa_filtered_file}" ]]; then
    # Filtered gene trees are not analyzed by Grampa, but it does not disturb gene tree IDs.
    # For example, if the 136th gene tree is filtered out, it is replaced with a placeholder text in grampa_trees_filtered.txt.
    # And GT-136 does not appear in the Grampa outputs. Still, the 137th gene tree is labeled correctly as GT-137.
    sed -e "s/$/;/" "${grampa_filtered_file}" > "grampa_input_gene_trees_filtered.nwk"
  fi

  grampa_det_file=""
  grampa_out_file=""
  grampa_checknums_file=""
  for candidate in "./grampa_out/grampa_det.txt" "./grampa_out/grampa-detailed.txt"; do
    if [[ -s "${candidate}" ]]; then
      grampa_det_file="${candidate}"
      break
    fi
  done
  for candidate in "./grampa_out/grampa_out.txt" "./grampa_out/grampa-scores.txt"; do
    if [[ -s "${candidate}" ]]; then
      grampa_out_file="${candidate}"
      break
    fi
  done
  for candidate in "./grampa_out/grampa_checknums.txt" "./grampa_out/grampa-checknums.txt"; do
    if [[ -s "${candidate}" ]]; then
      grampa_checknums_file="${candidate}"
      break
    fi
  done

  if [[ -z "${grampa_det_file}" || -z "${grampa_out_file}" || -z "${grampa_checknums_file}" ]]; then
    echo "Grampa output files are missing. Ending Grampa."
    return 0
  fi

  python "${gg_support_dir}/parse_grampa.py" \
    --grampa_det "${grampa_det_file}" \
    --grampa_out "${grampa_out_file}" \
    --gene_trees "./grampa_input_gene_trees.nwk" \
    --species_tree "./grampa_input_species_tree.nwk" \
    --ncpu "${GG_TASK_CPUS}" \
    --sorted_gene_tree_file_names "./busco_genetree_filenames.txt"

  if [[ -s "${grampa_checknums_file}" && -s "${grampa_det_file}" && -s "${grampa_out_file}" && -s "grampa_summary.tsv" ]]; then
    if ! awk -F'\t' '/^The MUL-tree with the minimum parsimony score/ {print $NF; found=1} END{exit(found?0:1)}' "${grampa_out_file}" > "${outdir}/best_mul_tree.nwk"; then
      awk -F'\t' '
        /^MT-/ {
          score = $NF + 0
          if (!seen || score < best_score) {
            best_score = score
            best_tree = $(NF-1)
            seen = 1
          }
        }
        END {
          if (seen) {
            print best_tree
          }
        }
      ' "${grampa_out_file}" > "${outdir}/best_mul_tree.nwk"
    fi
    cp_out "${grampa_checknums_file}" "${outdir}/grampa_checknums.txt"
    cp_out "${grampa_det_file}" "${outdir}/grampa_det.txt"
    cp_out "${grampa_out_file}" "${outdir}/grampa_out.txt"
    if [[ -s "./grampa_input_gene_trees_filtered.nwk" ]]; then
      mv_out "./grampa_input_gene_trees_filtered.nwk" "${outdir}"
    fi
    mv_out "./grampa_input_species_tree.nwk" "${outdir}"
    mv_out "./grampa_input_gene_trees.nwk" "${outdir}"
    mv_out "./grampa_summary.tsv" "${outfile}"
    rm -rf -- "./grampa_out"
  else
    echo "Grampa output files are missing. Ending Grampa."
  fi
}

task="BUSCO-based Grampa analysis for the polyploidization history with BUSCO DNA trees"
disable_if_no_input_file "run_busco_grampa_dna" "${file_dated_species_tree}"
if [[ ! -s "${file_busco_grampa_dna}" && ${run_busco_grampa_dna} -eq 1 ]]; then
  gg_step_start "${task}"

  busco_grampa "${dir_busco_rooted_nwk_dna}" "$(dirname "${file_busco_grampa_dna}")" "${file_busco_grampa_dna}"
else
  gg_step_skip "${task}"
fi

task="BUSCO-based Grampa analysis for the polyploidization history with BUSCO protein trees"
disable_if_no_input_file "run_busco_grampa_pep" "${file_dated_species_tree}"
if [[ ! -s "${file_busco_grampa_pep}" && ${run_busco_grampa_pep} -eq 1 ]]; then
  gg_step_start "${task}"

  busco_grampa "${dir_busco_rooted_nwk_pep}" "$(dirname "${file_busco_grampa_pep}")" "${file_busco_grampa_pep}"
else
  gg_step_skip "${task}"
fi

task="Orthogroup-based Grampa analysis for polyploidization history"
if [[ ! -s "${file_orthogroup_grampa}" && ${run_orthogroup_grampa} -eq 1 ]]; then
  gg_step_start "${task}"

  og_ids=()
  mapfile -t og_ids < <(
    python -c 'import pandas,sys; d=pandas.read_csv(sys.argv[1],sep="\t",header=0); ids=d.loc[(d["Total"]>=int(sys.argv[2]))&(d["Total"]<=int(sys.argv[3])),"Orthogroup"].astype(str).tolist(); print("\n".join(ids))' \
      "${file_orthogroup_genecount_selected}" "${min_gene_orthogroup_grampa}" "${max_gene_orthogroup_grampa}"
  )
  file_names=()
  mapfile -t file_names < <(gg_find_file_basenames "${dir_og_rooted_tree}")
  echo "Number of files in ${dir_og_rooted_tree}: ${#file_names[@]}"
  echo "Number of selected orthogroups with ${min_gene_orthogroup_grampa}<=gene number<=${max_gene_orthogroup_grampa}: ${#og_ids[@]}"
  if [[ -e ./tmp.orthogroup_grampa_indir ]]; then
    rm -rf -- ./tmp.orthogroup_grampa_indir
  fi
  mkdir -p ./tmp.orthogroup_grampa_indir
  for file_name in "${file_names[@]}"; do
    for og_id in "${og_ids[@]}"; do
      if [[ "${file_name}" == "${og_id}"* ]]; then
        cp_out "${dir_og_rooted_tree}"/"${file_name}" ./tmp.orthogroup_grampa_indir
        mapfile -t og_ids < <(printf "%s\n" "${og_ids[@]}" | grep -v -Fx -- "${og_id}" || true)
        break
      fi
    done
  done

  busco_grampa "./tmp.orthogroup_grampa_indir" "$(dirname "${file_orthogroup_grampa}")" "${file_orthogroup_grampa}"
else
  gg_step_skip "${task}"
fi

task='CAFE analysis'
disable_if_no_input_file "run_cafe" "${file_orthogroup_genecount_selected}" "${file_dated_species_tree}"
if [[ (! -s "${file_cafe_summary_all_pdf}" || ! -s "${file_cafe_summary_significant_pdf}") && ${run_cafe} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_cafe}"

  if [[ ! -s "${dir_cafe_output}/Gamma_asr.tre" || ! -s "${dir_cafe_output}/Gamma_count.tab" || ! -s "${dir_cafe_output}/Gamma_change.tab" ]]; then
    if [[ -d "${dir_cafe_orthogroup_selection}" ]]; then
      rm -rf -- "${dir_cafe_orthogroup_selection}"
    fi
    python "${gg_support_dir}/cafe_orthogroup_selection.py" \
      --genecount "${file_orthogroup_genecount_selected}" \
      --dated_species_tree "${file_dated_species_tree}" \
      --output_dir "${dir_cafe_orthogroup_selection}" \
      --max_size_differential "${max_size_differential_cafe}"

    if [[ -d "${dir_cafe_output}" ]]; then
      rm -rf -- "${dir_cafe_output}"
    fi
    cafe5 \
      --infile "${dir_cafe_orthogroup_selection}/cafe_input.tsv" \
      --tree "${file_dated_species_tree}" \
      --n_gamma_cats "${n_gamma_cats_cafe}" \
      --pvalue 0.05 \
      --cores "${GG_TASK_CPUS}" \
      --output_prefix "${dir_cafe_output}"
  else
    echo "CAFE output files already exist. Skipping CAFE run."
  fi

  if [[ -s "${dir_cafe_output}/Gamma_asr.tre" && -s "${dir_cafe_output}/Gamma_count.tab" && -s "${dir_cafe_output}/Gamma_change.tab" ]]; then
    if [[ -d "${dir_cafe}/each_family_plot" ]]; then
      rm -rf -- "${dir_cafe}/each_family_plot"
    fi
    if ! Rscript "${gg_support_dir}/cafe_plot_each_family.r" \
      "${dir_cafe_output}/Gamma_asr.tre" \
      "${dir_cafe_output}/Gamma_count.tab" \
      "${dir_cafe_output}/Gamma_change.tab" \
      "${dir_cafe}/each_family_plot" \
      "${GG_TASK_CPUS}"; then
      echo "Error in Rscript cafe_plot_each_family.r. Exiting."
      exit 1
    fi

    if [[ -d "$(dirname "${file_cafe_summary_all_pdf}")" ]]; then
      rm -rf -- "$(dirname "${file_cafe_summary_all_pdf}")"
    fi
    Rscript "${gg_support_dir}/cafe_plot_summary.r" \
      "${dir_cafe_output}/Gamma_asr.tre" \
      "${dir_cafe_output}/Gamma_change.tab" \
      "$(dirname "${file_cafe_summary_all_pdf}")"

    Rscript "${gg_support_dir}/cafe_plot_branch_id.r" \
      "${dir_cafe_output}/Gamma_asr.tre" \
      "$(dirname "${file_cafe_summary_all_pdf}")"
  else
    echo "CAFE did not finish successfully. Exiting."
    exit 1
  fi

  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task="GO enrichment analysis"
disable_if_no_input_file "run_go_enrichment" "${dir_cafe_output}/Gamma_change.tab" "${dir_cafe_output}/Gamma_branch_probabilities.tab" "${file_gene_id}" "${file_go_annotation}"
if [[ ! -s "${file_go_enrichment_significant}" && ${run_go_enrichment} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "$(dirname "${file_go_enrichment_significant}")"
  if ! Rscript "${gg_support_dir}/cafe_go_enrichment.r" \
    "${dir_cafe_output}/Gamma_change.tab" \
    "${dir_cafe_output}/Gamma_branch_probabilities.tab" \
    "${file_gene_id}" \
    "${file_go_annotation}" \
    "$(dirname "${file_go_enrichment_significant}")" \
    "${target_branch_go}" \
    "${change_direction_go}" \
    "${go_category}"; then
    echo "Error in Rscript cafe_go_enrichment.r. Exiting."
    exit 1
  fi
  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

if [[ ${delete_tmp_dir} -eq 1 ]]; then
  if [[ -d "${dir_tmp}" && -n "${dir_tmp:-}" && "${dir_tmp}" != "/" ]]; then
    echo "Removing tmp directory: ${dir_tmp}"
    rm -rf -- "${dir_tmp}"
  elif [[ -n "${dir_tmp:-}" && "${dir_tmp}" == "/" ]]; then
    echo "Refusing to delete unsafe tmp directory: '${dir_tmp}'"
  fi
fi

echo "$(date): end"

#!/usr/bin/env bash
set -euo pipefail

gg_core_bootstrap="/script/support/gg_core_bootstrap.sh"
if [[ ! -s "${gg_core_bootstrap}" ]]; then
  gg_core_bootstrap="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../support/gg_core_bootstrap.sh"
fi
# shellcheck disable=SC1090
source "${gg_core_bootstrap}"
unset gg_core_bootstrap

### Start: Job-supplied configuration ###
# Configuration variables are provided by gg_hgt_entrypoint.sh.
### End: Job-supplied configuration ###

gg_bootstrap_core_runtime "${BASH_SOURCE[0]:-$0}" "base" 1 1

run_hgt_eval="${run_hgt_eval:-1}"
hgt_min_branch_score="${hgt_min_branch_score:-0.45}"
hgt_use_taxonomy_db="${hgt_use_taxonomy_db:-1}"
hgt_contamination_dir="${hgt_contamination_dir:-}"

file_orthogroup_db="${gg_workspace_output_dir}/orthogroup/gg_orthogroup.db"
dir_hgt="${gg_workspace_output_dir}/hgt"
default_hgt_contamination_dir="${gg_workspace_output_dir}/species_cds_contamination_removal_tsv"

enable_all_run_flags_for_debug_mode

if [[ "${run_hgt_eval}" != "0" && "${run_hgt_eval}" != "1" ]]; then
  echo "Invalid binary flag value: run_hgt_eval=${run_hgt_eval} (expected 0 or 1)"
  exit 1
fi
if [[ "${hgt_use_taxonomy_db}" != "0" && "${hgt_use_taxonomy_db}" != "1" ]]; then
  echo "Invalid binary flag value: hgt_use_taxonomy_db=${hgt_use_taxonomy_db} (expected 0 or 1)"
  exit 1
fi
if ! [[ "${hgt_min_branch_score}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "Invalid hgt_min_branch_score: ${hgt_min_branch_score}"
  exit 1
fi

mkdir -p "${dir_hgt}"

if [[ ${run_hgt_eval} -eq 1 ]]; then
  if [[ ! -s "${file_orthogroup_db}" ]]; then
    echo "Skipping HGT evaluation because the orthogroup database was not found: ${file_orthogroup_db}"
    exit 0
  fi

  hgt_taxonomy_dbfile=""
  if [[ ${hgt_use_taxonomy_db} -eq 1 ]]; then
    if ensure_ete_taxonomy_db "${gg_workspace_dir}"; then
      hgt_taxonomy_dbfile=$(workspace_taxonomy_dbfile "${gg_workspace_dir}")
    else
      echo "Warning: Failed to prepare the ETE taxonomy DB. Continuing with best-hit name heuristics only." >&2
    fi
  fi

  contamination_arg=""
  if [[ -n "${hgt_contamination_dir}" ]]; then
    if [[ -d "${hgt_contamination_dir}" ]]; then
      contamination_arg="${hgt_contamination_dir}"
    else
      echo "Warning: HGT contamination directory was provided but not found. Skipping contamination input: ${hgt_contamination_dir}" >&2
    fi
  elif [[ -d "${default_hgt_contamination_dir}" ]]; then
    contamination_arg="${default_hgt_contamination_dir}"
  fi

  python "${gg_support_dir}/score_hgt_candidates.py" \
    --dbpath "${file_orthogroup_db}" \
    --branch_out "${dir_hgt}/hgt_branch_candidates.tsv" \
    --gene_out "${dir_hgt}/hgt_gene_candidates.tsv" \
    --orthogroup_out "${dir_hgt}/hgt_orthogroup_summary.tsv" \
    --min_branch_score "${hgt_min_branch_score}" \
    --dir_contamination_tsv "${contamination_arg}" \
    --taxonomy_dbfile "${hgt_taxonomy_dbfile}"
fi

echo "$(date): Exiting Singularity environment"

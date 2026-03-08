#!/usr/bin/env bash
set -euo pipefail

gg_core_bootstrap="/script/support/gg_core_bootstrap.sh"
if [[ ! -s "${gg_core_bootstrap}" ]]; then
  gg_core_bootstrap="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../support/gg_core_bootstrap.sh"
fi
# shellcheck disable=SC1090
source "${gg_core_bootstrap}"
unset gg_core_bootstrap

gg_bootstrap_core_runtime "${BASH_SOURCE[0]:-$0}" "base" 1 1

mode_transcriptome_assembly="${mode_transcriptome_assembly:-auto}"
ncpu_progress_summary="${ncpu_progress_summary:-${GG_TASK_CPUS:-1}}"

gg_workspace_output_dir=$(workspace_output_root "${gg_workspace_dir}")
gg_workspace_input_dir=$(workspace_input_root "${gg_workspace_dir}")
dir_orthogroup="${gg_workspace_output_dir}/orthogroup"
dir_transcriptome_assembly="${gg_workspace_output_dir}/transcriptome_assembly"
file_orthogroup_genecount_selected="${gg_workspace_output_dir}/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv"

if [[ -d "${dir_orthogroup}" ]]; then
  echo ""
  echo "Checking directory: ${dir_orthogroup}"
  if [[ ! -s "${file_orthogroup_genecount_selected}" ]]; then
    echo "Skipping orthogroup summary because the gene-count table was not found: ${file_orthogroup_genecount_selected}"
  elif [[ ! -d "${dir_orthogroup}/amas_original" || ! -d "${dir_orthogroup}/amas_cleaned" ]]; then
    echo "Skipping orthogroup summary because required AMAS directories were not found under: ${dir_orthogroup}"
  else
    python "${gg_support_dir}/orthogroup_output_summary.py" \
      --dir_og "${dir_orthogroup}" \
      --genecount "${file_orthogroup_genecount_selected}" \
      --ncpu "${ncpu_progress_summary}" \
      --out orthogroup_summary.tsv
  fi
fi

if [[ -d "${dir_transcriptome_assembly}" ]]; then
  echo ""
  echo "Checking directory: ${dir_transcriptome_assembly}"
  python "${gg_support_dir}/transcriptome_assembly_output_summary.py" \
    --dir_transcriptome_assembly "${dir_transcriptome_assembly}" \
    --gg_workspace_input_dir "${gg_workspace_input_dir}" \
    --mode "${mode_transcriptome_assembly}" \
    --ncpu "${ncpu_progress_summary}" \
    --out transcriptome_assembly_summary.tsv
fi

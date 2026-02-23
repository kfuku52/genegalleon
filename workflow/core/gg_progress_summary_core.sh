#!/usr/bin/env bash
set -euo pipefail

dir_pg="/workspace"
dir_script="/script/support"
source "${dir_script}/gg_util.sh"
gg_source_home_bashrc
gg_prepare_cmd_runtime "${dir_pg}" "base" 1 1

mode_transcriptome_assembly="${mode_transcriptome_assembly:-auto}"
ncpu_progress_summary="${ncpu_progress_summary:-${NSLOTS:-1}}"

dir_pg_output=$(workspace_output_root "${dir_pg}")
dir_pg_input=$(workspace_input_root "${dir_pg}")
dir_orthogroup="${dir_pg_output}/orthogroup"
dir_transcriptome_assembly="${dir_pg_output}/transcriptome_assembly"
file_orthogroup_genecount_selected="${dir_pg_output}/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv"

if [[ -d "${dir_orthogroup}" ]]; then
  echo ""
  echo "Checking directory: ${dir_orthogroup}"
  if [[ ! -s "${file_orthogroup_genecount_selected}" ]]; then
    echo "Skipping orthogroup summary because the gene-count table was not found: ${file_orthogroup_genecount_selected}"
  elif [[ ! -d "${dir_orthogroup}/amas_original" || ! -d "${dir_orthogroup}/amas_cleaned" ]]; then
    echo "Skipping orthogroup summary because required AMAS directories were not found under: ${dir_orthogroup}"
  else
    python "${dir_script}/orthogroup_output_summary.py" \
      --dir_og "${dir_orthogroup}" \
      --genecount "${file_orthogroup_genecount_selected}" \
      --ncpu "${ncpu_progress_summary}" \
      --out orthogroup_summary.tsv
  fi
fi

if [[ -d "${dir_transcriptome_assembly}" ]]; then
  echo ""
  echo "Checking directory: ${dir_transcriptome_assembly}"
  python "${dir_script}/transcriptome_assembly_output_summary.py" \
    --dir_transcriptome_assembly "${dir_transcriptome_assembly}" \
    --dir_pg_input "${dir_pg_input}" \
    --mode "${mode_transcriptome_assembly}" \
    --ncpu "${ncpu_progress_summary}" \
    --out transcriptome_assembly_summary.tsv
fi

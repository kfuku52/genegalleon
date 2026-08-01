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
gene_family_output_storage=$(printf '%s' "${gene_family_output_storage:-${GG_COMMON_GENE_FAMILY_OUTPUT_STORAGE:-zip}}" | tr '[:upper:]' '[:lower:]')
gene_family_zip_compression=$(printf '%s' "${gene_family_zip_compression:-${GG_COMMON_GENE_FAMILY_ZIP_COMPRESSION:-adaptive}}" | tr '[:upper:]' '[:lower:]')
gene_family_zip_compression_level="${gene_family_zip_compression_level:-${GG_COMMON_GENE_FAMILY_ZIP_COMPRESSION_LEVEL:-6}}"
gene_family_zip_workers="${gene_family_zip_workers:-${GG_COMMON_GENE_FAMILY_ZIP_WORKERS:-1}}"
gene_family_tmp_retention_days="${gene_family_tmp_retention_days:-${GG_COMMON_GENE_FAMILY_TMP_RETENTION_DAYS:-7}}"
gene_family_tmp_max_dirs="${gene_family_tmp_max_dirs:-${GG_COMMON_GENE_FAMILY_TMP_MAX_DIRS:-100}}"
gene_family_tmp_max_bytes="${gene_family_tmp_max_bytes:-${GG_COMMON_GENE_FAMILY_TMP_MAX_BYTES:-107374182400}}"
gene_family_tmp_max_files="${gene_family_tmp_max_files:-${GG_COMMON_GENE_FAMILY_TMP_MAX_FILES:-100000}}"
case "${gene_family_output_storage}" in
  zip|files)
    ;;
  raw)
    gene_family_output_storage="files"
    ;;
  *)
    echo "Invalid gene_family_output_storage: ${gene_family_output_storage}; expected zip, files, or raw." >&2
    exit 1
    ;;
esac
case "${gene_family_zip_compression}" in
  adaptive|deflate|store)
    ;;
  *)
    echo "Invalid gene_family_zip_compression: ${gene_family_zip_compression}; expected adaptive, deflate, or store." >&2
    exit 1
    ;;
esac
if [[ ! "${gene_family_zip_compression_level}" =~ ^[0-9]+$ || ${gene_family_zip_compression_level} -gt 9 ]]; then
  echo "Invalid gene_family_zip_compression_level: ${gene_family_zip_compression_level}; expected an integer from 0 through 9." >&2
  exit 1
fi
if [[ ! "${gene_family_zip_workers}" =~ ^[0-9]+$ || ${gene_family_zip_workers} -lt 1 || ${gene_family_zip_workers} -gt 4 ]]; then
  echo "Invalid gene_family_zip_workers: ${gene_family_zip_workers}; expected an integer from 1 through 4." >&2
  exit 1
fi
gene_family_archive_write_args=(
  --compression "${gene_family_zip_compression}"
  --compression-level "${gene_family_zip_compression_level}"
  --workers "${gene_family_zip_workers}"
)
if [[ ! "${gene_family_tmp_retention_days}" =~ ^[0-9]+$ ]]; then
  echo "Invalid gene_family_tmp_retention_days: ${gene_family_tmp_retention_days}; expected a non-negative integer." >&2
  exit 1
fi
if [[ ! "${gene_family_tmp_max_dirs}" =~ ^[0-9]+$ ]]; then
  echo "Invalid gene_family_tmp_max_dirs: ${gene_family_tmp_max_dirs}; expected a non-negative integer." >&2
  exit 1
fi
if [[ ! "${gene_family_tmp_max_bytes}" =~ ^[0-9]+$ ]]; then
  echo "Invalid gene_family_tmp_max_bytes: ${gene_family_tmp_max_bytes}; expected a non-negative integer." >&2
  exit 1
fi
if [[ ! "${gene_family_tmp_max_files}" =~ ^[0-9]+$ ]]; then
  echo "Invalid gene_family_tmp_max_files: ${gene_family_tmp_max_files}; expected a non-negative integer." >&2
  exit 1
fi

gg_workspace_output_dir=$(workspace_output_root "${gg_workspace_dir}")
gg_workspace_input_dir=$(workspace_input_root "${gg_workspace_dir}")
dir_orthogroup="${gg_workspace_output_dir}/orthogroup"
dir_query2family="${gg_workspace_output_dir}/query2family"
dir_transcriptome_assembly="${gg_workspace_output_dir}/transcriptome_assembly"
dir_query_gene="${gg_workspace_input_dir}/query_gene"
file_orthogroup_genecount_selected="${gg_workspace_output_dir}/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv"

assert_gene_family_storage_ready() {
  local output_root=$1
  local conversion_marker="${output_root}/.gg_archives/storage-conversion.pending"
  if [[ -e "${conversion_marker}" || -L "${conversion_marker}" ]]; then
    echo "Gene-family storage conversion is in progress or needs to be resumed: ${conversion_marker}" >&2
    echo "Refusing to create a progress summary from a converting store." >&2
    exit 1
  fi
}

cleanup_gene_family_tmp_if_enabled() {
  local output_root=$1
  if [[
    "${gene_family_output_storage}" != "zip"
    || (
      ${gene_family_tmp_retention_days} -eq 0
      && ${gene_family_tmp_max_dirs} -eq 0
      && ${gene_family_tmp_max_bytes} -eq 0
      && ${gene_family_tmp_max_files} -eq 0
    )
  ]]; then
    return 0
  fi
  if ! python "${gg_support_dir}/gene_family_output_store.py" cleanup-tmp \
    --root "${output_root}" \
    --older-than-days "${gene_family_tmp_retention_days}" \
    --max-directories "${gene_family_tmp_max_dirs}" \
    --max-bytes "${gene_family_tmp_max_bytes}" \
    --max-files "${gene_family_tmp_max_files}" \
    --nonblocking
  then
    echo "Warning: Failed to clean stale gene-family temporary directories below ${output_root}." >&2
  fi
}

if [[ -d "${dir_orthogroup}" ]]; then
  assert_gene_family_storage_ready "${dir_orthogroup}"
  echo ""
  echo "Checking directory: ${dir_orthogroup}"
  if [[ ! -s "${file_orthogroup_genecount_selected}" ]]; then
    echo "Skipping orthogroup summary because the gene-count table was not found: ${file_orthogroup_genecount_selected}"
  else
    python "${gg_support_dir}/orthogroup_output_summary.py" \
      --dir_og "${dir_orthogroup}" \
      --genecount "${file_orthogroup_genecount_selected}" \
      --ncpu "${ncpu_progress_summary}" \
      --out orthogroup_summary.tsv
    if [[ "${gene_family_output_storage}" == "zip" ]]; then
      python "${gg_support_dir}/gene_family_output_store.py" archive-completed \
        --root "${dir_orthogroup}" \
        --mode orthogroup \
        --genecount "${file_orthogroup_genecount_selected}" \
        "${gene_family_archive_write_args[@]}" \
        --min-files 1 \
        --nonblocking
    fi
  fi
  cleanup_gene_family_tmp_if_enabled "${dir_orthogroup}"
fi

if [[ -d "${dir_query2family}" ]]; then
  assert_gene_family_storage_ready "${dir_query2family}"
  echo ""
  echo "Checking directory: ${dir_query2family}"
  if [[ ! -d "${dir_query_gene}" ]]; then
    echo "Skipping query2family summary because the query_gene input directory was not found: ${dir_query_gene}"
  elif [[ -z "$(find "${dir_query_gene}" -mindepth 1 -maxdepth 1 -type f ! -name '.*' -print -quit 2>/dev/null)" ]]; then
    echo "Skipping query2family summary because the query_gene input directory has no visible files: ${dir_query_gene}"
  else
    python "${gg_support_dir}/query2family_output_summary.py" \
      --dir_query2family "${dir_query2family}" \
      --dir_query_gene "${dir_query_gene}" \
      --ncpu "${ncpu_progress_summary}" \
      --out query2family_summary.tsv
    if [[ "${gene_family_output_storage}" == "zip" ]]; then
      python "${gg_support_dir}/gene_family_output_store.py" archive-completed \
        --root "${dir_query2family}" \
        --mode query2family \
        --query-dir "${dir_query_gene}" \
        "${gene_family_archive_write_args[@]}" \
        --min-files 1 \
        --nonblocking
    fi
  fi
  cleanup_gene_family_tmp_if_enabled "${dir_query2family}"
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

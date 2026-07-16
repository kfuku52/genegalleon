#!/usr/bin/env bash
set -euo pipefail

gg_core_bootstrap="/script/support/gg_core_bootstrap.sh"
if [[ ! -s "${gg_core_bootstrap}" ]]; then
  gg_core_bootstrap="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../support/gg_core_bootstrap.sh"
fi
# shellcheck disable=SC1090
source "${gg_core_bootstrap}"
unset gg_core_bootstrap
gg_source_common_params_from_core "${BASH_SOURCE[0]:-$0}"

### Start: Job-supplied configuration ###
# Configuration variables are provided by gg_transcriptome_generation_entrypoint.sh.
busco_lineage="${busco_lineage:-${GG_COMMON_BUSCO_LINEAGE:-auto}}"
contamination_removal_rank="${contamination_removal_rank:-domain}"
contamination_removal_target_taxon="${contamination_removal_target_taxon:-}"
### End: Job-supplied configuration ###

### ----------------------------------------------------------------------- ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

gg_bootstrap_core_runtime "${BASH_SOURCE[0]:-$0}" "base" 1 1
genetic_code="${genetic_code:-${GG_COMMON_GENETIC_CODE:-1}}"
# shellcheck disable=SC1090
source "${gg_support_dir}/gg_busco.sh"
delete_tmp_dir=${delete_tmp_dir:-1}
mode_transcriptome_assembly=$(echo "${mode_transcriptome_assembly:-sraid}" | tr '[:upper:]' '[:lower:]')
requested_assembly_method=$(printf '%s' "${assembly_method:-auto}" | tr '[:upper:]' '[:lower:]' | tr '_' '-')
amalgkit_ncbi_metadata_max_concurrency="${amalgkit_ncbi_metadata_max_concurrency:-20}"
amalgkit_ncbi_download_max_concurrency="${amalgkit_ncbi_download_max_concurrency:-20}"
amalgkit_aws_download_max_concurrency="${amalgkit_aws_download_max_concurrency:-20}"
amalgkit_gcp_download_max_concurrency="${amalgkit_gcp_download_max_concurrency:-20}"
amalgkit_quant_backend="${amalgkit_quant_backend:-auto}"
amalgkit_oarfish_seq_tech="${amalgkit_oarfish_seq_tech:-auto}"
amalgkit_oarfish_options="${amalgkit_oarfish_options:-}"
amalgkit_sra_strategy_query="${amalgkit_sra_strategy_query:-\"RNA-seq\"[Strategy] OR \"EST\"[Strategy] OR \"CLONE\"[Strategy]}"
busco_lineage_resolved=""
contamination_removal_rank_for_remove_contaminated_sequences="$(
  gg_normalize_contamination_removal_rank_for_remove_contaminated_sequences "${contamination_removal_rank}"
)"
contamination_removal_rank_for_amalgkit="$(
  gg_normalize_contamination_removal_rank_for_amalgkit "${contamination_removal_rank}"
)"
effective_assembly_method=""
detected_metadata_run_count=0
detected_short_read_run_count=0
detected_pacbio_run_count=0
detected_ont_cdna_run_count=0
detected_ont_direct_rna_run_count=0
detected_long_read_unknown_run_count=0
detected_has_long_reads=0
detected_has_short_reads=0
detected_has_pacbio=0
detected_has_ont=0
detected_has_ont_cdna=0
detected_has_ont_direct_rna=0
detected_has_long_read_unknown=0
detected_input_class="unknown"
detected_metadata_instrument_field="none"
classified_short_single_fastq_files=()
classified_short_left_fastq_files=()
classified_short_right_fastq_files=()
classified_long_fastq_files=()
classified_pacbio_fastq_files=()
classified_ont_fastq_files=()

gg_core_stage_library="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/stages/gg_transcriptome_generation_core_functions.sh"
# shellcheck disable=SC1090
source "${gg_core_stage_library}"
unset gg_core_stage_library

































# Setting modes
if [[ ${gg_debug_mode:-0} -eq 1 ]]; then
  enable_all_run_flags_for_debug_mode
  echo "gg debug mode: max_assembly_input_fastq_size is set to 30,000,000 bp."
  max_assembly_input_fastq_size="30,000,000"
fi

dir_transcriptome_assembly_output="${gg_workspace_output_dir}/transcriptome_assembly"
dir_input_fastq="${gg_workspace_input_dir}/species_rnaseq"
dir_input_sra_list="${gg_workspace_input_dir}/query_sra_id"
dir_input_amalgkit_metadata="${gg_workspace_input_dir}/amalgkit_metadata"
dir_generated_amalgkit_metadata="${dir_transcriptome_assembly_output}/amalgkit_metadata"
dir_amalgkit_quant="${dir_transcriptome_assembly_output}/amalgkit_quant"

if [[ "${mode_transcriptome_assembly}" != "auto" && "${mode_transcriptome_assembly}" != "sraid" && "${mode_transcriptome_assembly}" != "fastq" && "${mode_transcriptome_assembly}" != "metadata" ]]; then
  echo "Invalid mode_transcriptome_assembly: ${mode_transcriptome_assembly}"
  echo 'mode_transcriptome_assembly must be one of "auto", "sraid", "fastq", or "metadata". Exiting.'
  exit 1
fi

fastq_mode_dirs=()
if [[ -d "${dir_input_fastq}" ]]; then
  mapfile -t fastq_mode_dirs < <(find "${dir_input_fastq}" -mindepth 1 -maxdepth 1 -type d ! -name '.*' | sort)
fi
sra_mode_files=()
if [[ -d "${dir_input_sra_list}" ]]; then
  mapfile -t sra_mode_files < <(find "${dir_input_sra_list}" -mindepth 1 -maxdepth 1 -type f ! -name '.*' | sort)
fi
metadata_mode_files=()
if [[ -d "${dir_input_amalgkit_metadata}" ]]; then
  mapfile -t metadata_mode_files < <(find "${dir_input_amalgkit_metadata}" -mindepth 1 -maxdepth 1 -type f ! -name '.*' | sort)
fi

available_transcriptome_modes=()
if [[ ${#fastq_mode_dirs[@]} -gt 0 ]]; then
  available_transcriptome_modes+=("fastq")
fi
if [[ ${#sra_mode_files[@]} -gt 0 ]]; then
  available_transcriptome_modes+=("sraid")
fi
if [[ ${#metadata_mode_files[@]} -gt 0 ]]; then
  available_transcriptome_modes+=("metadata")
fi

selected_transcriptome_mode="${mode_transcriptome_assembly}"
if [[ "${mode_transcriptome_assembly}" == "auto" ]]; then
  if [[ ${#available_transcriptome_modes[@]} -eq 0 ]]; then
    echo "No valid transcriptome input mode is available. Skipping transcriptome assembly workflow."
    exit 0
  fi
  if [[ ${#available_transcriptome_modes[@]} -gt 1 ]]; then
    echo "mode_transcriptome_assembly=auto is ambiguous. Multiple input layouts are available: ${available_transcriptome_modes[*]}"
    echo 'Set mode_transcriptome_assembly explicitly to one of "sraid", "fastq", or "metadata". Exiting.'
    exit 1
  fi
  selected_transcriptome_mode="${available_transcriptome_modes[0]}"
  echo "Auto-detected transcriptome input mode: ${selected_transcriptome_mode}"
else
  case "${selected_transcriptome_mode}" in
    "fastq")
      if [[ ! -d "${dir_input_fastq}" ]]; then
        echo "Missing mode_transcriptome_assembly=fastq input directory: ${dir_input_fastq}"
        exit 1
      fi
      if [[ ${#fastq_mode_dirs[@]} -eq 0 ]]; then
        echo "Input directory is empty for mode_transcriptome_assembly=fastq: ${dir_input_fastq}"
        exit 1
      fi
      ;;
    "sraid")
      if [[ ! -d "${dir_input_sra_list}" ]]; then
        echo "Missing mode_transcriptome_assembly=sraid input directory: ${dir_input_sra_list}"
        exit 1
      fi
      if [[ ${#sra_mode_files[@]} -eq 0 ]]; then
        echo "Input directory is empty for mode_transcriptome_assembly=sraid: ${dir_input_sra_list}"
        exit 1
      fi
      ;;
    "metadata")
      if [[ ! -d "${dir_input_amalgkit_metadata}" ]]; then
        echo "Missing mode_transcriptome_assembly=metadata input directory: ${dir_input_amalgkit_metadata}"
        exit 1
      fi
      if [[ ${#metadata_mode_files[@]} -eq 0 ]]; then
        echo "Input directory is empty for mode_transcriptome_assembly=metadata: ${dir_input_amalgkit_metadata}"
        exit 1
      fi
      ;;
  esac
fi
if [[ ! "${GG_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] || [[ ${GG_ARRAY_TASK_ID} -lt 1 ]]; then
  echo "Invalid GG_ARRAY_TASK_ID value (must be a positive integer): ${GG_ARRAY_TASK_ID}"
  exit 1
fi

if [[ "${selected_transcriptome_mode}" == "fastq" ]]; then
  echo 'Mode: fastq input'
  dirs=( "${fastq_mode_dirs[@]}" )
  id=$((GG_ARRAY_TASK_ID - 1))
  if [[ ${id} -ge ${#dirs[@]} ]]; then
    echo "Input fastq directory not found, probably due to the out-of-range specification of array jobs. Exiting."
    exit 1
  fi
  dir_species_fastq="${dirs[${id}]}"
  sp_ub="$(basename "${dir_species_fastq}")"
  files_fastq=()
  mapfile -t files_fastq < <(find "${dir_species_fastq}" -maxdepth 1 -type f ! -name '.*' \( -name "*.fq" -o -name "*.fastq" -o -name "*.fq.gz" -o -name "*.fastq.gz" \) | sort)
  echo "Input fastq directory: ${dir_species_fastq}"
  echo "Species: ${sp_ub}"
  printf 'Input fastq files: %s\n' "${files_fastq[*]}"
  if [[ ! -d "${dir_species_fastq}" || ${#files_fastq[@]} -eq 0 ]]; then
    echo "No FASTQ files were found for species ${sp_ub} in: ${dir_species_fastq}. Exiting."
    exit 1
  fi
elif [[ "${selected_transcriptome_mode}" == "sraid" ]]; then
  echo 'Mode: sraid input'
  files=( "${sra_mode_files[@]}" )
  id=$((GG_ARRAY_TASK_ID - 1))
  if [[ ${id} -ge ${#files[@]} ]]; then
    echo "Input SRA list file not found, probably due to the out-of-range specification of array jobs. Exiting."
    exit 1
  fi
  file_input_sra_list="${files[${id}]}"
  if [[ ! -f "${file_input_sra_list}" ]]; then
    echo "Input SRA list file not found, probably due to the out-of-range specification of array jobs. Exiting."
    exit 1
  fi
  sra_ids=()
  mapfile -t sra_ids < "${file_input_sra_list}"
  sp_ub="$(gg_species_name_from_path_or_dot "${file_input_sra_list}")"
  echo "Input SRA list file: ${file_input_sra_list}"
  echo "Species: ${sp_ub}"
  echo "Number of input SRA IDs: ${#sra_ids[@]}"
  printf 'Input SRA IDs: %s\n' "${sra_ids[*]}"
  if [[ ${#sra_ids[@]} -eq 0 ]]; then
    echo "SRA IDs not found, probably due to the out-of-range specification of array jobs. Exiting."
    exit 1
  fi
elif [[ "${selected_transcriptome_mode}" == "metadata" ]]; then
  echo 'Mode: metadata input'
  files=( "${metadata_mode_files[@]}" )
  id=$((GG_ARRAY_TASK_ID - 1))
  if [[ ${id} -ge ${#files[@]} ]]; then
    echo "Input metadata file not found, probably due to the out-of-range specification of array jobs. Exiting."
    exit 1
  fi
  file_metadata="${files[${id}]}"
  sp_ub="$(gg_species_name_from_path_or_dot "${file_metadata}")"
  echo "Input metadata file: ${file_metadata}"
  echo "Species: ${sp_ub}"
  if [[ ! -f "${file_metadata}" ]]; then
    echo "Input metadata file not found, probably due to the out-of-range specification of array jobs. Exiting."
    exit 1
  fi
fi

dir_tmp="${dir_transcriptome_assembly_output}/tmp/${GG_ARRAY_TASK_ID}_${sp_ub}"
dir_amalgkit_getfastq_sp="${dir_transcriptome_assembly_output}/amalgkit_getfastq/${sp_ub}"
dir_amalgkit_download_dir="${gg_workspace_downloads_dir}"
dir_amalgkit_download_lock_dir="${dir_amalgkit_download_dir}/locks"
dir_mmseqs2_db="${gg_workspace_downloads_dir}/mmseqs2"
file_input_amalgkit_metadata="${dir_input_amalgkit_metadata}/${sp_ub}_metadata.tsv"
file_generated_amalgkit_metadata="${dir_generated_amalgkit_metadata}/${sp_ub}_metadata.tsv"
if [[ "${selected_transcriptome_mode}" == "metadata" ]]; then
  file_amalgkit_metadata="${file_input_amalgkit_metadata}"
else
  file_amalgkit_metadata="${file_generated_amalgkit_metadata}"
fi
file_amalgkit_read_technology="${dir_transcriptome_assembly_output}/amalgkit_read_technology/${sp_ub}_read_technology.tsv"
file_amalgkit_read_technology_summary_sh="${dir_tmp}/metadata/read_technology.summary.sh"
file_amalgkit_getfastq_legacy_safely_removed_flag=${dir_transcriptome_assembly_output}/amalgkit_getfastq/${sp_ub}_safely_removed.txt
file_isoform="${dir_transcriptome_assembly_output}/assembled_transcripts_with_isoforms/${sp_ub}_isoform.fa.gz"
file_corset_clusters="${dir_transcriptome_assembly_output}/corset_clusters/${sp_ub}_corset.clusters.tsv"
file_corset_counts="${dir_transcriptome_assembly_output}/corset_counts/${sp_ub}_corset.counts.tsv"
file_longestcds="${dir_transcriptome_assembly_output}/longest_cds/${sp_ub}_longestCDS.fa.gz"
file_longestcds_transcript="${dir_transcriptome_assembly_output}/longest_cds_transcript/${sp_ub}_longestCDS.transcript.fa.gz"
file_longestcds_fx2tab="${dir_transcriptome_assembly_output}/longest_cds_fx2tab/${sp_ub}_longestCDS.fx2tab_cds.tsv"
file_longestcds_mmseqs2taxonomy="${dir_transcriptome_assembly_output}/longest_cds_mmseqs2taxonomy/${sp_ub}_longestCDS.mmseqs2taxonomy.tsv"
file_longestcds_contamination_removal_fasta="${dir_transcriptome_assembly_output}/longest_cds_contamination_removal_fasta/${sp_ub}_longestCDS_contamination_removal.fa.gz"
file_longestcds_contamination_removal_tsv="${dir_transcriptome_assembly_output}/longest_cds_contamination_removal_tsv/${sp_ub}_longestCDS_contamination_removal.tsv"
file_assembly_stat="${dir_transcriptome_assembly_output}/assembly_stat/${sp_ub}_assembly_stat.tsv"
file_busco_full_cdna_isoforms="${dir_transcriptome_assembly_output}/busco_full_cdna_isoforms/${sp_ub}_busco.full.tsv"
file_busco_short_cdna_isoforms="${dir_transcriptome_assembly_output}/busco_short_cdna_isoforms/${sp_ub}_busco.short.txt"
file_busco_full_longest_cds="${dir_transcriptome_assembly_output}/busco_full_longest_cds/${sp_ub}_busco.full.tsv"
file_busco_short_longest_cds="${dir_transcriptome_assembly_output}/busco_short_longest_cds/${sp_ub}_busco.short.txt"
file_busco_full_longest_cds_filtered="${dir_transcriptome_assembly_output}/busco_full_longest_cds_contamination_removal/${sp_ub}_busco.full.tsv"
file_busco_short_longest_cds_filtered="${dir_transcriptome_assembly_output}/busco_short_longest_cds_contamination_removal/${sp_ub}_busco.short.txt"
file_amalgkit_merge_efflen="${dir_transcriptome_assembly_output}/amalgkit_merge/${sp_ub}/${sp_ub}_eff_length.tsv"
file_amalgkit_merge_count="${dir_transcriptome_assembly_output}/amalgkit_merge/${sp_ub}/${sp_ub}_est_counts.tsv"
file_amalgkit_merge_tpm="${dir_transcriptome_assembly_output}/amalgkit_merge/${sp_ub}/${sp_ub}_tpm.tsv"
file_amalgkit_merge_metadata="${dir_transcriptome_assembly_output}/amalgkit_merge/${sp_ub}/${sp_ub}_metadata.tsv"
file_multispecies_summary="${dir_transcriptome_assembly_output}/annotation_summary/assembly_stat_summary.pdf"

ensure_dir "${dir_tmp}"
ensure_dir "${dir_amalgkit_download_dir}"
ensure_dir "${dir_amalgkit_download_lock_dir}"
cd "${dir_tmp}"
amalgkit_ncbi_metadata_max_concurrency="$(normalize_amalgkit_download_limit_value "${amalgkit_ncbi_metadata_max_concurrency}" "amalgkit_ncbi_metadata_max_concurrency")"
amalgkit_ncbi_download_max_concurrency="$(normalize_amalgkit_download_limit_value "${amalgkit_ncbi_download_max_concurrency}" "amalgkit_ncbi_download_max_concurrency")"
amalgkit_aws_download_max_concurrency="$(normalize_amalgkit_download_limit_value "${amalgkit_aws_download_max_concurrency}" "amalgkit_aws_download_max_concurrency")"
amalgkit_gcp_download_max_concurrency="$(normalize_amalgkit_download_limit_value "${amalgkit_gcp_download_max_concurrency}" "amalgkit_gcp_download_max_concurrency")"


gg_core_execution_stage_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/stages/gg_transcriptome_generation"
# shellcheck source=/dev/null
source "${gg_core_execution_stage_dir}/01_metadata_and_reads.sh"
# shellcheck source=/dev/null
source "${gg_core_execution_stage_dir}/02_assembly.sh"
# shellcheck source=/dev/null
source "${gg_core_execution_stage_dir}/03_cds_and_quality.sh"
# shellcheck source=/dev/null
source "${gg_core_execution_stage_dir}/04_quant_and_summary.sh"
unset gg_core_execution_stage_dir

#!/usr/bin/env bash
set -euo pipefail

# Load shared defaults when available.
gg_core_self="${BASH_SOURCE[0]:-/script/core/gg_transcriptome_generation_core.sh}"
gg_common_params_file=""
for gg_common_params_candidate in \
  "$(cd "$(dirname "${gg_core_self}")" && pwd)/gg_common_params.sh" \
  "$(cd "$(dirname "${gg_core_self}")" && pwd)/../gg_common_params.sh" \
  "/script/gg_common_params.sh"
do
  if [[ -s "${gg_common_params_candidate}" ]]; then
    gg_common_params_file="${gg_common_params_candidate}"
    break
  fi
done
if [[ -n "${gg_common_params_file}" ]]; then
  # shellcheck disable=SC1090
  source "${gg_common_params_file}"
fi
unset gg_common_params_file gg_common_params_candidate gg_core_self

### Start: Job-supplied configuration ###
# Configuration variables are provided by gg_transcriptome_generation_entrypoint.sh.
busco_lineage="${busco_lineage:-${GG_COMMON_BUSCO_LINEAGE:-embryophyta_odb12}}"
contamination_removal_rank="${contamination_removal_rank:-${GG_COMMON_CONTAMINATION_REMOVAL_RANK:-phylum}}"
### End: Job-supplied configuration ###

### ----------------------------------------------------------------------- ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

dir_pg="/workspace"
dir_script="/script/support"
source "${dir_script}/gg_util.sh" # Load utility functions
gg_source_home_bashrc
gg_prepare_cmd_runtime "${dir_pg}" "base" 1 1
delete_tmp_dir=${delete_tmp_dir:-1}

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

# Setting modes
if [[ ${gg_debug_mode:-0} -eq 1 ]]; then
  enable_all_run_flags_for_debug_mode
  echo "gg debug mode: max_assembly_input_fastq_size is set to 30,000,000 bp."
  max_assembly_input_fastq_size="30,000,000"
fi

dir_transcriptome_assembly_input="${dir_pg_input}/transcriptome_assembly"
dir_transcriptome_assembly_output="${dir_pg_output}/transcriptome_assembly"
dir_input_fastq="${dir_pg_input}/species_rnaseq"
dir_input_sra_list="${dir_pg_input}/query_sra_id"
dir_amalgkit_metadata="${dir_transcriptome_assembly_input}/amalgkit_metadata"
dir_amalgkit_quant="${dir_transcriptome_assembly_output}/amalgkit_quant"

if [[ ${mode_fastq} -eq 1 && ! -d "${dir_input_fastq}" ]]; then
  echo "Missing mode_fastq input directory: ${dir_input_fastq}"
  mode_fastq=0
fi
if [[ ${mode_fastq} -eq 1 ]]; then
  fastq_mode_dirs=()
  mapfile -t fastq_mode_dirs < <(find "${dir_input_fastq}" -mindepth 1 -maxdepth 1 -type d ! -name '.*' | sort)
  if [[ ${#fastq_mode_dirs[@]} -eq 0 ]]; then
    echo "Input directory is empty for mode_fastq: ${dir_input_fastq}"
    mode_fastq=0
  fi
fi
if [[ ${mode_sraid} -eq 1 ]]; then
  if [[ ! -d "${dir_input_sra_list}" ]]; then
    echo "Missing mode_sraid input directory: ${dir_input_sra_list}"
    mode_sraid=0
  else
    sra_mode_files=()
    mapfile -t sra_mode_files < <(find "${dir_input_sra_list}" -mindepth 1 -maxdepth 1 -type f ! -name '.*' | sort)
    if [[ ${#sra_mode_files[@]} -eq 0 ]]; then
      echo "Input directory is empty for mode_sraid: ${dir_input_sra_list}"
      mode_sraid=0
    fi
  fi
fi
if [[ ${mode_metadata} -eq 1 ]]; then
  if [[ ! -d "${dir_amalgkit_metadata}" ]]; then
    echo "Missing mode_metadata input directory: ${dir_amalgkit_metadata}"
    mode_metadata=0
  else
    metadata_mode_files=()
    mapfile -t metadata_mode_files < <(find "${dir_amalgkit_metadata}" -mindepth 1 -maxdepth 1 -type f ! -name '.*' | sort)
    if [[ ${#metadata_mode_files[@]} -eq 0 ]]; then
      echo "Input directory is empty for mode_metadata: ${dir_amalgkit_metadata}"
      mode_metadata=0
    fi
  fi
fi
if [[ ${mode_fastq} -eq 0 && ${mode_sraid} -eq 0 && ${mode_metadata} -eq 0 ]]; then
  echo "No valid transcriptome input mode is available. Skipping transcriptome assembly workflow."
  exit 0
fi
if [[ ! "${SGE_TASK_ID}" =~ ^[0-9]+$ ]] || [[ ${SGE_TASK_ID} -lt 1 ]]; then
  echo "Invalid SGE_TASK_ID value (must be a positive integer): ${SGE_TASK_ID}"
  exit 1
fi

if [[ ${mode_fastq} -eq 1 ]]; then
  echo 'Mode: fastq input'
  if [[ ! -d "${dir_input_fastq}" ]]; then
    echo "Input directory does not exist: ${dir_input_fastq}"
    exit 1
  fi
  dirs=()
  mapfile -t dirs < <(find "${dir_input_fastq}" -mindepth 1 -maxdepth 1 -type d ! -name '.*' | sort)
  if [[ ${#dirs[@]} -eq 0 ]]; then
    echo "Input directory is empty: ${dir_input_fastq}"
    exit 1
  fi
  id=$((SGE_TASK_ID-1))
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
  echo "Input fastq files: ${files_fastq[@]}"
  if [[ ! -d "${dir_species_fastq}" || ${#files_fastq[@]} -eq 0 ]]; then
    echo "No FASTQ files were found for species ${sp_ub} in: ${dir_species_fastq}. Exiting."
    exit 1
  fi
elif [[ ${mode_sraid} -eq 1 ]]; then
  echo 'Mode: sraid input'
  if [[ ! -d "${dir_input_sra_list}" ]]; then
    echo "Input directory does not exist: ${dir_input_sra_list}"
    exit 1
  fi
  files=()
  mapfile -t files < <(find "${dir_input_sra_list}" -mindepth 1 -maxdepth 1 -type f ! -name '.*' | sort)
  if [[ ${#files[@]} -eq 0 ]]; then
    echo "Input directory is empty: ${dir_input_sra_list}"
    exit 1
  fi
  id=$((SGE_TASK_ID-1))
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
  echo "Input SRA IDs: ${sra_ids[@]}"
  if [[ ${#sra_ids[@]} -eq 0 ]]; then
    echo "SRA IDs not found, probably due to the out-of-range specification of array jobs. Exiting."
    exit 1
  fi
elif [[ ${mode_metadata} -eq 1 ]]; then
  echo 'Mode: metadata input'
  if [[ ! -d "${dir_amalgkit_metadata}" ]]; then
    echo "Input directory does not exist: ${dir_amalgkit_metadata}"
    exit 1
  fi
  files=()
  mapfile -t files < <(find "${dir_amalgkit_metadata}" -mindepth 1 -maxdepth 1 -type f ! -name '.*' | sort)
  if [[ ${#files[@]} -eq 0 ]]; then
    echo "Input directory is empty: ${dir_amalgkit_metadata}"
    exit 1
  fi
  id=$((SGE_TASK_ID-1))
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

dir_tmp="${dir_transcriptome_assembly_output}/tmp/${SGE_TASK_ID}_${sp_ub}"
dir_amalgkit_getfastq_sp="${dir_transcriptome_assembly_output}/amalgkit_getfastq/${sp_ub}"
file_amalgkit_metadata="${dir_amalgkit_metadata}/${sp_ub}_metadata.tsv"
file_amalgkit_getfastq_safely_removed_flag=${dir_transcriptome_assembly_output}/amalgkit_getfastq/${sp_ub}_safely_removed.txt
file_rRNA_contamination_report="${dir_transcriptome_assembly_output}/amalgkit_metadata_rRNA_mapping_rate/${sp_ub}_rRNA_mapping_rate.tsv"
file_isoform="${dir_transcriptome_assembly_output}/assembled_transcripts_with_isoforms/${sp_ub}_isoform.fa.gz"
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
cd "${dir_tmp}"

task="amalgkit metadata/integrate"
if [[ ! -s "${file_amalgkit_metadata}" && ${run_amalgkit_metadata_or_integrate} -eq 1 ]]; then
  gg_step_start "${task}"
  if [[ -e "./metadata" ]]; then
    rm -rf -- "./metadata"
  fi

  if [[ ${mode_sraid} -eq 1 ]]; then
    search_string=$(tr -s "\n" < "${file_input_sra_list}" | sed ':a;N;$!ba;s/\n/ OR /g' | sed -e "s/ OR $//" -e "s/^/(/" -e "s/$/)/")
    search_string="${search_string} AND \"Illumina\"[Platform] AND (\"RNA-seq\"[Strategy] OR \"EST\"[Strategy])"
    echo "Entrez search string: ${search_string}"

    amalgkit metadata \
    --out_dir "./" \
    --search_string "${search_string}"

    sp_space="${sp_ub//_/ }"
    {
        head -n 1 "./metadata/metadata.tsv";
        grep -F -- "${sp_space}" "./metadata/metadata.tsv" || true;
    } | sed -e "s/\t\t\tno\t/\tyes\tyes\tno\t/g" > "./metadata.tsv"
    if [[ $(wc -l < "./metadata.tsv") -le 1 ]]; then
      echo "No metadata rows matched species '${sp_space}' in ./metadata/metadata.tsv. Exiting."
      exit 1
    fi

  elif [[ ${mode_fastq} -eq 1 ]]; then
    amalgkit integrate \
    --out_dir "./" \
    --fastq_dir "${dir_species_fastq}" \
    --threads "${NSLOTS}" \
    --remove_tmp yes

    mv_out ./metadata_private_fastq.tsv ./metadata/metadata.tsv
    python -c "import sys,pandas as pd; d=pd.read_csv('./metadata/metadata.tsv',sep='\t',header=0); d.loc[:,'scientific_name']=sys.argv[1]; d.to_csv('./metadata.tsv',sep='\t',index=False)" "${sp_ub}"
  fi

  if [[ -s "./metadata.tsv" ]]; then
    mv_out "./metadata.tsv" "${file_amalgkit_metadata}"
  else
    echo "metadata.tsv not found. Exiting."
    exit 1
  fi

  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task="amalgkit getfastq"
amalgkit_fastq_files=()
if [[ -d "${dir_amalgkit_getfastq_sp}" ]]; then
  mapfile -t amalgkit_fastq_files < <(find "${dir_amalgkit_getfastq_sp}" -type f -name "*.amalgkit.fastq.gz" | sort)
fi
echo "Number of amalgkit getfastq fastq files: ${#amalgkit_fastq_files[@]}"
echo "is_fastq_requiring_downstream_analysis_done: $(is_fastq_requiring_downstream_analysis_done)"
if [[ ( ${#amalgkit_fastq_files[@]} -eq 0 && ${run_amalgkit_getfastq} -eq 1 ) && $(is_fastq_requiring_downstream_analysis_done) -eq 0 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_amalgkit_getfastq_sp}"

  if [[ -e "${file_amalgkit_getfastq_safely_removed_flag}" ]]; then
    rm -f -- "${file_amalgkit_getfastq_safely_removed_flag}"
  fi

	  if amalgkit getfastq \
	    --out_dir "${dir_tmp}" \
	    --metadata "${file_amalgkit_metadata}" \
	    --threads "${NSLOTS}" \
    --remove_sra yes \
    --remove_tmp yes \
    --read_name 'trinity' \
    --aws yes \
	    --ncbi yes \
	    --redo no; then
	    echo "amalgkit getfastq safely finished."
	    shopt -s nullglob
	    getfastq_outputs=( "${dir_tmp}"/getfastq/* )
	    shopt -u nullglob
	    if [[ ${#getfastq_outputs[@]} -eq 0 ]]; then
	      echo "amalgkit getfastq finished but no output files were found under: ${dir_tmp}/getfastq"
	      exit 1
	    fi
	    mv_out "${getfastq_outputs[@]}" "${dir_amalgkit_getfastq_sp}"
	    rm -rf -- "${dir_tmp}/getfastq"
  else
    status_amalgkit=$?
    echo "amalgkit getfastq did not safely finish. Exiting."
    echo "amalgkit getfastq exit code: ${status_amalgkit}"
    exit 1
  fi
  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task='rRNA contamination report'
if [[ ( ! -s "${file_rRNA_contamination_report}" ) && ${run_rRNA_contamination_report} -eq 1 ]]; then
  gg_step_start "${task}"
  if [[ ! -d "${dir_amalgkit_getfastq_sp}" ]]; then
    echo "amalgkit getfastq output directory not found: ${dir_amalgkit_getfastq_sp}. Exiting."
    exit 1
  fi
  if [[ -z "$(find "${dir_amalgkit_getfastq_sp}" -type f -name "*.amalgkit.fastq.gz" -print -quit 2>/dev/null)" ]]; then
    echo "No amalgkit getfastq FASTQ files were found in: ${dir_amalgkit_getfastq_sp}. Exiting."
    exit 1
  fi
  if ! silva_rrna_ref=$(ensure_silva_rrna_ref_db "${dir_pg}"); then
    echo "Failed to prepare SILVA rRNA reference. Exiting."
    exit 1
  fi

  recreate_dir "./fasta"
  recreate_dir "./index"
  recreate_dir "./quant"
  recreate_dir "./merge"
  if [[ -e "./getfastq" ]]; then
    rm -rf -- "./getfastq"
  fi
  ln -s "${dir_amalgkit_getfastq_sp}" "./getfastq"
  file_reference_fasta_link="./fasta/${sp_ub}_dummy_rRNA_for_kallisto_index.fa.gz"
  ln -s "${silva_rrna_ref}" "${file_reference_fasta_link}"
  set +e
  kallisto index --make-unique -i "./index/${sp_ub}.idx" "${file_reference_fasta_link}"
  status_kallisto=$?
  set -e
  if [[ ${status_kallisto} -ne 0 ]]; then
    echo "kallisto index failed with exit code ${status_kallisto}"
    exit 1
  fi

  set +e
  amalgkit quant \
  --out_dir "./" \
  --threads "${NSLOTS}" \
  --metadata "${file_amalgkit_metadata}" \
  --clean_fastq no \
  --fasta_dir "./fasta" \
  --build_index no
  exit_code_amalgkit_quant=$?

  amalgkit merge \
  --out_dir "./" \
  --metadata "${file_amalgkit_metadata}"
  exit_code_amalgkit_merge=$?
  set -e

  if [[ ${exit_code_amalgkit_quant} -ne 0 || ${exit_code_amalgkit_merge} -ne 0 ]]; then
    echo "amalgkit quant for rRNA contamination rate estimation failed with exit code ${exit_code_amalgkit_quant}"
    echo "If you got an error 'FileNotFoundError: Could not find index file.', the specified species name may not be the same as the one in the metadata.tsv file."
    echo "In that case, please correct the species name in the metadata.tsv file (or in your input files/directories) and run this pipeline again."
    echo "Exiting."
    exit 1
  else
    echo "amalgkit quant for rRNA contamination rate estimation finished successfully"
    mv_out ./merge/metadata.tsv "${file_rRNA_contamination_report}"
    rm -f -- "${file_reference_fasta_link}"
    rm -rf -- "./fasta"
    rm -f -- "./getfastq" # Do not put -r, otherwise the original getfastq files will be deleted.
  fi
else
  gg_step_skip "${task}"
fi

task='De novo transcriptome assembly'
if [[ ! -s "${file_isoform}" && ${run_assembly} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_parent_dir "${file_isoform}"
  if [[ ! -d "${dir_amalgkit_getfastq_sp}" ]]; then
    echo "amalgkit getfastq output directory not found: ${dir_amalgkit_getfastq_sp}. Exiting."
    exit 1
  fi
  if [[ -z "$(find "${dir_amalgkit_getfastq_sp}" -type f -name "*.amalgkit.fastq.gz" -print -quit 2>/dev/null)" ]]; then
    echo "No amalgkit getfastq FASTQ files were found in: ${dir_amalgkit_getfastq_sp}. Exiting."
    exit 1
  fi

  mapfile -t files_right < <(find "${dir_amalgkit_getfastq_sp}" -type f -name "*_2.amalgkit.fastq.gz" | sort)
  if [[ ${#files_right[@]} -eq 0 ]]; then
    echo "Paired-end samples were not detected. SE reads will be used for transcriptome assembly."
    lib_layout='single'
  else
    echo "Paired-end samples were detected. PE reads will be used for transcriptome assembly and SE reads, if any, will not be used."
    lib_layout='paired'
  fi

  selected_fastq_dir="./tmp_selected_fastq"
  if [[ -e "${selected_fastq_dir}" ]]; then
    rm -rf -- "${selected_fastq_dir}"
  fi
  mkdir -p "${selected_fastq_dir}"
  if [[ ${assembly_method} == "rnaSPAdes" && ${protocol_rna_seq} == "mixed" ]]; then
    if [[ ${lib_layout} == 'single' ]]; then
      mapfile -t candidate_fastq_files < <(find "${dir_amalgkit_getfastq_sp}" -type f -name "*.amalgkit.fastq.gz" | sort)
    elif [[ ${lib_layout} == 'paired' ]]; then
      mapfile -t candidate_fastq_files < <(find "${dir_amalgkit_getfastq_sp}" -type f -name "*_1.amalgkit.fastq.gz" | sort)
    fi
    n_pair=${#candidate_fastq_files[@]}
    if [[ ${n_pair} -ge 10 ]]; then
      echo "Selecting top 9 ${lib_layout}-end fastq files by size, as rnaSPAdes accepts fewer than 10 libraries with different protocols."
      echo "For details, please refer to the rnaSPAdes manual: https://ablab.github.io/spades/rna.html"
      echo "Selected fastq files:"
      if [[ ${lib_layout} == 'single' ]]; then
        find "${dir_amalgkit_getfastq_sp}" -type f -name "*.amalgkit.fastq.gz" -exec du -b {} + \
        | sort -nr | head -n 9 | cut -f2 | while read -r file; do
          cp_out "${file}" "${selected_fastq_dir}"; echo "${file}"
        done
      elif [[ ${lib_layout} == 'paired' ]]; then
        find "${dir_amalgkit_getfastq_sp}" -type f -name "*_1.amalgkit.fastq.gz" -exec du -b {} + \
        | sort -nr | head -n 9 | cut -f2 | while read -r file1; do
          file2="${file1/_1.amalgkit.fastq.gz/_2.amalgkit.fastq.gz}"
          cp_out "${file1}" "${selected_fastq_dir}"; echo "${file1}"
          cp_out "${file2}" "${selected_fastq_dir}"; echo "${file2}"
        done
      fi
    else
      echo "All ${lib_layout}-end fastq files will be used for transcriptome assembly."
      selected_fastq_dir=${dir_amalgkit_getfastq_sp}
    fi
  else
    echo "All ${lib_layout}-end fastq files will be used for transcriptome assembly."
    selected_fastq_dir=${dir_amalgkit_getfastq_sp}
  fi

  if [[ ${lib_layout} == 'single' ]]; then
    total_fastq_len=$(get_total_fastq_len "${selected_fastq_dir}" "*.amalgkit.fastq.gz")
  elif [[ ${lib_layout} == 'paired' ]]; then
    total_fastq_len1=$(get_total_fastq_len "${selected_fastq_dir}" "*_1.amalgkit.fastq.gz")
    total_fastq_len2=$(get_total_fastq_len "${selected_fastq_dir}" "*_2.amalgkit.fastq.gz")
    total_fastq_len=$((${total_fastq_len1}+${total_fastq_len2}))
  fi
  max_assembly_input_fastq_size="${max_assembly_input_fastq_size//,/}"
  if [[ ${total_fastq_len} -gt ${max_assembly_input_fastq_size} ]]; then
    echo "Total ${lib_layout} fastq length is ${total_fastq_len} bp, which is greater than ${max_assembly_input_fastq_size} bp."
    echo "Only ${max_assembly_input_fastq_size} bp will be used."
    assembly_input_fastq_dir="./tmp_assembly_input_fastq"
    if [[ -e "${assembly_input_fastq_dir}" ]]; then
      rm -rf -- "${assembly_input_fastq_dir}"
    fi
    mkdir -p "${assembly_input_fastq_dir}"
    proportion=$(awk -v max="${max_assembly_input_fastq_size}" -v total="${total_fastq_len}" 'BEGIN {printf "%.3f\n", max/total}')
    echo "Proportion of fastq reads to be used: ${proportion} (${max_assembly_input_fastq_size}/${total_fastq_len})"
    files=()
    if [[ ${lib_layout} == 'single' ]]; then
      mapfile -t files < <(find "${selected_fastq_dir}" -type f -name "*.amalgkit.fastq.gz" | sort)
    elif [[ ${lib_layout} == 'paired' ]]; then
      mapfile -t files1 < <(find "${selected_fastq_dir}" -type f -name "*_1.amalgkit.fastq.gz" | sort)
      mapfile -t files2 < <(find "${selected_fastq_dir}" -type f -name "*_2.amalgkit.fastq.gz" | sort)
      files=( "${files1[@]}" "${files2[@]}" )
    fi
    for file in "${files[@]}"; do
      seqkit sample --proportion "${proportion}" --rand-seed 11 --out-file "${assembly_input_fastq_dir}/$(basename "${file}")" "${file}"
    done
    echo "Total fastq length of the subsampled fastq files: $(get_total_fastq_len "${assembly_input_fastq_dir}" "*.amalgkit.fastq.gz") bp"
  else
    echo "Total fastq length is ${total_fastq_len} bp, which is less than ${max_assembly_input_fastq_size} bp. All fastq reads will be used."
    assembly_input_fastq_dir=${selected_fastq_dir}
  fi

  nslots_assembly=$((${NSLOTS}-${assembly_cpu_offset}))
  memory_assembly=$((${MEM_PER_HOST}-${assembly_ram_offset}))
  if [[ ${nslots_assembly} -lt 1 ]]; then
    echo "Adjusted nslots_assembly from ${nslots_assembly} to 1 because it must be >=1."
    nslots_assembly=1
  fi
  if [[ ${memory_assembly} -lt 1 ]]; then
    echo "Adjusted memory_assembly from ${memory_assembly}G to 1G because it must be >=1G."
    memory_assembly=1
  fi
  bflyHeapSpaceMax=${MEM_PER_SLOT} # GB
  echo "Number of offset CPUs for transcriptome assembly is ${assembly_cpu_offset}."
  echo "${NSLOTS} CPUs and ${MEM_PER_HOST}G RAM are allocated to this job."
  echo "${nslots_assembly} CPUs and ${memory_assembly}G RAM are used for transcriptome assembly."

  files_single=()
  files_left=()
  files_right=()
  mapfile -t files_single < <(find "${assembly_input_fastq_dir}" -type f -name "*.amalgkit.fastq.gz" | sort)
  mapfile -t files_left < <(find "${assembly_input_fastq_dir}" -type f -name "*_1.amalgkit.fastq.gz" | sort)
  mapfile -t files_right < <(find "${assembly_input_fastq_dir}" -type f -name "*_2.amalgkit.fastq.gz" | sort)

  if [[ ${lib_layout} == 'single' && ${#files_single[@]} -eq 0 ]]; then
    echo "No single-end fastq files were found for transcriptome assembly. Exiting."
    exit 1
  fi
  if [[ ${lib_layout} == 'paired' ]]; then
    if [[ ${#files_left[@]} -eq 0 || ${#files_right[@]} -eq 0 || ${#files_left[@]} -ne ${#files_right[@]} ]]; then
      echo "Paired-end input files were not detected correctly in ${assembly_input_fastq_dir}."
      echo "Detected left/right counts: ${#files_left[@]}/${#files_right[@]}. Exiting."
      exit 1
    fi
  fi

  if [[ ${assembly_method} == 'Trinity' ]]; then
    if [[ ${lib_layout} == 'single' ]]; then
      in_single="$(IFS=","; echo "${files_single[*]}")"
    elif [[ ${lib_layout} == 'paired' ]]; then
      in_left="$(IFS=","; echo "${files_left[*]}")"
      in_right="$(IFS=","; echo "${files_right[*]}")"
    fi
    if [[ ${lib_layout} == 'single' ]]; then
      Trinity \
      --seqType fq \
      --CPU "${nslots_assembly}" \
      --max_memory "${memory_assembly}G" \
      --min_contig_length 200 \
      --output trinity \
      --full_cleanup \
      --NO_SEQTK \
      --bflyHeapSpaceMax "${bflyHeapSpaceMax}G" \
      --bflyGCThreads 1 \
      --bflyCalculateCPU \
      --single "${in_single}"
    else
      Trinity \
      --seqType fq \
      --CPU "${nslots_assembly}" \
      --max_memory "${memory_assembly}G" \
      --min_contig_length 200 \
      --output trinity \
      --full_cleanup \
      --NO_SEQTK \
      --bflyHeapSpaceMax "${bflyHeapSpaceMax}G" \
      --bflyGCThreads 1 \
      --bflyCalculateCPU \
      --left "${in_left}" \
      --right "${in_right}"
    fi
    # For --NO_SEQTK, see https://github.com/trinityrnaseq/trinityrnaseq/issues/787
    if [[ -s "${dir_tmp}/trinity.Trinity.fasta" ]]; then
      seqkit seq --threads "${NSLOTS}" "${dir_tmp}/trinity.Trinity.fasta" --out-file "tmp.isoform.fa.gz"
      mv_out "tmp.isoform.fa.gz" "${file_isoform}"
    fi
  elif [[ ${assembly_method} == 'rnaSPAdes' ]]; then
    rnaspades_input_args=()
    if [[ ${protocol_rna_seq} == "same" ]]; then
      if [[ ${lib_layout} == 'single' ]]; then
        for i in "${!files_single[@]}"; do
          rnaspades_input_args+=( --s1 "${files_single[i]}" )
        done
      elif [[ ${lib_layout} == 'paired' ]]; then
        for i in "${!files_left[@]}"; do
          rnaspades_input_args+=( --pe1-1 "${files_left[i]}" --pe1-2 "${files_right[i]}" )
        done
      fi
    elif [[ ${protocol_rna_seq} == "mixed" ]]; then
      if [[ ${lib_layout} == 'single' ]]; then
        for i in "${!files_single[@]}"; do
          j=$((i + 1))
          rnaspades_input_args+=( "--s${j}" "${files_single[i]}" )
        done
      elif [[ ${lib_layout} == 'paired' ]]; then
        for i in "${!files_left[@]}"; do
          j=$((i + 1))
          rnaspades_input_args+=( "--pe${j}-1" "${files_left[i]}" "--pe${j}-2" "${files_right[i]}" )
        done
      fi
    else
      echo "Invalid value for 'protocol_rna_seq'. Please specify either 'same' or 'mixed'."
      echo "Exiting."
      exit 1
    fi
    if [[ -d "${dir_tmp}/rnaspades_output" ]]; then
      rm -rf -- "${dir_tmp}/rnaspades_output"
    fi
    rnaspades.py \
    --threads "${nslots_assembly}" \
    --memory "${memory_assembly}" \
    -o rnaspades_output \
    "${rnaspades_input_args[@]}"
    if [[ -s "${dir_tmp}/rnaspades_output/transcripts.fasta" ]]; then
      seqkit seq --threads "${NSLOTS}" "${dir_tmp}/rnaspades_output/transcripts.fasta" --out-file "tmp.isoform.fa.gz"
      mv_out "tmp.isoform.fa.gz" "${file_isoform}"
    fi
  else
    echo "Invalid value for 'assembly_method'. Please specify either 'Trinity' or 'rnaSPAdes'."
    echo "Exiting."
    exit 1
  fi

  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task='Longest CDS extraction'
disable_if_no_input_file "run_longestcds" "${file_isoform}"
if [[ ( ! -s "${file_longestcds}" || ! -s "${file_longestcds_transcript}" ) && ${run_longestcds} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_parent_dir "${file_longestcds}"
  ensure_parent_dir "${file_longestcds_transcript}"

  if [[ ${orf_aggregation_level} = "p" ]]; then
    aggregate_expression="\.${orf_aggregation_level}[0-9].*"
  else
    aggregate_expression="\-${orf_aggregation_level}[0-9].*"
  fi
  if [[ ${orf_aggregation_level} == "c" && ${assembly_method} == "rnaSPAdes" ]]; then
    echo "The aggregation level 'c' is not supported for rnaSPAdes assemblies. Please set 'orf_aggregation_level' to either 'i' or 'g'."
    echo "Exiting"
    exit 1
  fi

  echo "scientific_name: ${sp_ub}"
  if [[ ${assembly_method} == 'Trinity' ]]; then
    seqkit sort --by-length --reverse --threads "${NSLOTS}" "${file_isoform}" \
    | sed -e "s/_/-/g" -e "s/^>[[:space:]]*/>${sp_ub}_/" -e "s/TRINITY-//" -e "s/[[:space:]].*//" \
    > "${sp_ub}.tmp.renamed.fa"
  elif [[ ${assembly_method} == 'rnaSPAdes' ]]; then
    seqkit sort --by-length --reverse --threads "${NSLOTS}" "${file_isoform}" \
    | sed -e "s/^>.*_g\([0-9]\+\)_i\([0-9]\+\)$/>${sp_ub}_g\1-i\2/" \
    > "${sp_ub}.tmp.renamed.fa"
  else
    echo "Invalid value for 'assembly_method'. Please specify either 'Trinity' or 'rnaSPAdes'."
    echo "Exiting."
    exit 1
  fi

  longestcds_tmp="${sp_ub}.tmp.longestcds.raw.fa"
  if ! cdskit longestcds \
    --seqfile "${sp_ub}.tmp.renamed.fa" \
    --outfile "${longestcds_tmp}" \
    --codontable 1 \
    --annotate_seqname no \
    --threads "${NSLOTS}"; then
    echo "Error: cdskit longestcds failed."
    echo "Please install a cdskit version that supports the 'longestcds' subcommand."
    exit 1
  fi
  if [[ ! -s "${longestcds_tmp}" ]]; then
    echo "Error: No CDS output generated by cdskit longestcds."
    exit 1
  fi

  seqkit rmdup --by-seq "${longestcds_tmp}" \
  | cdskit pad \
  | seqkit sort --by-name --threads "${NSLOTS}" \
  | cdskit aggregate -x "${aggregate_expression}" \
  > "${sp_ub}.tmp.aggregated.fa"

  awk '/^>/ {sub(/^>/, "", $0); sub(/[[:space:]].*$/, "", $0); sub(/\.p[0-9].*$/, "", $0); print}' \
  "${sp_ub}.tmp.aggregated.fa" > "${sp_ub}.tmp.aggregated.transcript_id.txt"

  if [[ -s "${sp_ub}.tmp.aggregated.fa" && -s "${sp_ub}.tmp.aggregated.transcript_id.txt" ]]; then
    echo "Output file detected for the task: ${task}"
    seqkit grep \
    --threads "${NSLOTS}" \
    --pattern-file "${sp_ub}.tmp.aggregated.transcript_id.txt" \
    "${sp_ub}.tmp.renamed.fa" \
    --out-file "${sp_ub}.tmp.longestcds.transcript.fa.gz"
    mv_out "${sp_ub}.tmp.longestcds.transcript.fa.gz" "${file_longestcds_transcript}"
    sed -e "s/[[:space:]].*//" -e "s/${aggregate_expression}//" "${sp_ub}.tmp.aggregated.fa" \
    | seqkit seq --threads "${NSLOTS}" --out-file "${sp_ub}.tmp.longestcds.fa.gz"
    mv_out "${sp_ub}.tmp.longestcds.fa.gz" "${file_longestcds}"
  else
    echo "Output file not detected for the task: ${task}"
  fi
  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task="seqkit fx2tab for the longest CDS sequences"
disable_if_no_input_file "run_longestcds_fx2tab" "${file_longestcds}"
if [[ ! -s "${file_longestcds_fx2tab}" && ${run_longestcds_fx2tab} -eq 1 ]]; then
  gg_step_start "${task}"

  seqkit fx2tab \
  --threads "${NSLOTS}" \
  --length \
  --name \
  --gc \
  --gc-skew \
  --header-line \
  --only-id \
  "${file_longestcds}" \
  > "tmp.cds_length.tsv"

  if [[ -s "tmp.cds_length.tsv" ]]; then
    echo "Output file detected for the task: ${task}"
    mv_out "tmp.cds_length.tsv" "${file_longestcds_fx2tab}"
  fi
else
  gg_step_skip "${task}"
fi

task="MMseqs2 Taxonomy of the CDS sequences"
disable_if_no_input_file "run_longestcds_mmseqs2taxonomy" "${file_longestcds}"
if [[ ! -s "${file_longestcds_mmseqs2taxonomy}" && ${run_longestcds_mmseqs2taxonomy} -eq 1 ]]; then
  gg_step_start "${task}"

  if ! ensure_mmseqs_uniref90_db "${dir_pg_db}/mmseqs2" "${NSLOTS}"; then
    echo "Failed to prepare MMseqs2 UniRef90 DB. Exiting."
    exit 1
  fi

  if [[ ! -s "queryDB" ]]; then
    mmseqs createdb "${file_longestcds}" queryDB
  fi

  if [[ ! -e "tmp_mmseqs2" ]]; then
    mkdir -p "tmp_mmseqs2"
  fi

  mmseqs taxonomy "queryDB" "${dir_pg_db}/mmseqs2/UniRef90_DB" "output_prefix" "tmp" \
  --split-mode 2 \
  --split-memory-limit $((${MEM_PER_HOST}*3/4))G \
  --majority 0.5 \
  --lca-mode 3 \
  --vote-mode 1 \
  --tax-lineage 2 \
  --orf-filter 1 \
  --threads "${NSLOTS}"

  mmseqs createtsv "queryDB" "output_prefix" "result.tsv" --threads "${NSLOTS}"

  if [[ -s "result.tsv" ]]; then
    echo "Output file detected for the task: ${task}"
    mv_out "result.tsv" "${file_longestcds_mmseqs2taxonomy}"
    rm -f -- queryDB*
    rm -f -- output_prefix.*
    rm -rf -- "tmp_mmseqs2"
  fi
else
  gg_step_skip "${task}"
fi

task="Contaminated sequence removal from the CDS sequences"
disable_if_no_input_file "run_longestcds_contamination_removal" "${file_longestcds}" "${file_longestcds_fx2tab}" "${file_longestcds_mmseqs2taxonomy}"
if [[ ( ! -s "${file_longestcds_contamination_removal_fasta}" || ! -s "${file_longestcds_contamination_removal_tsv}" ) && ${run_longestcds_contamination_removal} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_parent_dir "${file_longestcds_contamination_removal_fasta}"

  if ! ensure_ete_taxonomy_db "${dir_pg}"; then
    echo "Failed to prepare ETE taxonomy DB. Exiting."
    exit 1
  fi

  python "${dir_script}/remove_contaminated_sequences.py" \
  --fasta_file "${file_longestcds}" \
  --mmseqs2taxonomy_tsv "${file_longestcds_mmseqs2taxonomy}" \
  --fx2tab_tsv "${file_longestcds_fx2tab}" \
  --species_name "${sp_ub}" \
  --rank "${contamination_removal_rank}" \
  --ncpu "${NSLOTS}" \
  --rename_seq "no" \
  --verbose "no"

  if [[ -s "clean_sequences.fa" && -s "lineage_compatibility.tsv" ]]; then
    echo "Output file detected for the task: ${task}"
    seqkit seq --threads "${NSLOTS}" "clean_sequences.fa" --out-file "tmp.longestcds.clean.fa.gz"
    mv_out "tmp.longestcds.clean.fa.gz" "${file_longestcds_contamination_removal_fasta}"
    rm -f -- "clean_sequences.fa"
    mv_out "lineage_compatibility.tsv" "${file_longestcds_contamination_removal_tsv}"
  fi
else
  gg_step_skip "${task}"
fi

task='BUSCO for cDNA isoforms (isoform.fasta)'
disable_if_no_input_file "run_busco1" "${file_isoform}"
if [[ ( ! -s "${file_busco_full_cdna_isoforms}" || ! -s "${file_busco_short_cdna_isoforms}" ) && ${run_busco1} -eq 1 ]]; then
  gg_step_start "${task}"

  seqkit seq --threads "${NSLOTS}" "${file_isoform}" --out-file "busco_infile_cdna.fa"

  if ! dir_busco_db=$(ensure_busco_download_path "${dir_pg}" "${busco_lineage}"); then
    echo "Failed to prepare BUSCO dataset: ${busco_lineage}"
    exit 1
  fi
  dir_busco_lineage="${dir_busco_db}/lineages/${busco_lineage}"

  busco \
  --in "busco_infile_cdna.fa" \
  --mode transcriptome \
  --out "busco_tmp" \
  --cpu "${NSLOTS}" \
  --force \
  --evalue 1e-03 \
  --limit 20 \
  --lineage_dataset "${dir_busco_lineage}" \
  --download_path "${dir_busco_db}" \
  --offline

  if ! copy_busco_tables "./busco_tmp" "${busco_lineage}" "${file_busco_full_cdna_isoforms}" "${file_busco_short_cdna_isoforms}"; then
    echo "Failed to locate normalized BUSCO outputs for cDNA isoforms. Exiting."
    exit 1
  fi
  rm -rf -- "./busco_tmp"
  rm -f -- "busco_infile_cdna.fa"

  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task='BUSCO for longest CDS'
disable_if_no_input_file "run_busco2" "${file_longestcds}"
if [[ ( ! -s "${file_busco_full_longest_cds}" || ! -s "${file_busco_short_longest_cds}" ) && ${run_busco2} -eq 1 ]]; then
  gg_step_start "${task}"

  seqkit seq --threads "${NSLOTS}" "${file_longestcds}" --out-file "busco_infile_cds.fa"

  if ! dir_busco_db=$(ensure_busco_download_path "${dir_pg}" "${busco_lineage}"); then
    echo "Failed to prepare BUSCO dataset: ${busco_lineage}"
    exit 1
  fi
  dir_busco_lineage="${dir_busco_db}/lineages/${busco_lineage}"

  busco \
  --in "busco_infile_cds.fa" \
  --mode transcriptome \
  --out "busco_tmp" \
  --cpu "${NSLOTS}" \
  --force \
  --evalue 1e-03 \
  --limit 20 \
  --lineage_dataset "${dir_busco_lineage}" \
  --download_path "${dir_busco_db}" \
  --offline

  if ! copy_busco_tables "./busco_tmp" "${busco_lineage}" "${file_busco_full_longest_cds}" "${file_busco_short_longest_cds}"; then
    echo "Failed to locate normalized BUSCO outputs for longest CDS. Exiting."
    exit 1
  fi
  rm -rf -- "./busco_tmp"
  rm -f -- "busco_infile_cds.fa"

  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task='BUSCO for contamination-removed longest CDS'
disable_if_no_input_file "run_busco3" "${file_longestcds_contamination_removal_fasta}"
if [[ ( ! -s "${file_busco_full_longest_cds_filtered}" || ! -s "${file_busco_short_longest_cds_filtered}" ) && ${run_busco3} -eq 1 ]]; then
  gg_step_start "${task}"

  seqkit seq --threads "${NSLOTS}" "${file_longestcds_contamination_removal_fasta}" --out-file "busco_infile_cds.fa"

  if ! dir_busco_db=$(ensure_busco_download_path "${dir_pg}" "${busco_lineage}"); then
    echo "Failed to prepare BUSCO dataset: ${busco_lineage}"
    exit 1
  fi
  dir_busco_lineage="${dir_busco_db}/lineages/${busco_lineage}"

  busco \
  --in "busco_infile_cds.fa" \
  --mode transcriptome \
  --out "busco_tmp" \
  --cpu "${NSLOTS}" \
  --force \
  --evalue 1e-03 \
  --limit 20 \
  --lineage_dataset "${dir_busco_lineage}" \
  --download_path "${dir_busco_db}" \
  --offline

  if ! copy_busco_tables "./busco_tmp" "${busco_lineage}" "${file_busco_full_longest_cds_filtered}" "${file_busco_short_longest_cds_filtered}"; then
    echo "Failed to locate normalized BUSCO outputs for contamination-removed longest CDS. Exiting."
    exit 1
  fi
  rm -rf -- "./busco_tmp"
  rm -f -- "busco_infile_cds.fa"

  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task='Assembly statistics'
disable_if_no_input_file "run_assembly_stat" "${file_isoform}"
if [[ "${run_longestcds}" -eq 1 ]]; then disable_if_no_input_file "run_assembly_stat" "${file_longestcds}"; fi
if [[ "${run_longestcds_contamination_removal}" -eq 1 ]]; then disable_if_no_input_file "run_assembly_stat" "${file_longestcds_contamination_removal_fasta}"; fi
if [[ ! -s "${file_assembly_stat}" && ${run_assembly_stat} -eq 1 ]]; then
  gg_step_start "${task}"

  input_files=( "${file_isoform}" )
  if [[ -s "${file_longestcds}" ]]; then
    input_files+=( "${file_longestcds}" )
  fi
  if [[ -s "${file_longestcds_contamination_removal_fasta}" ]]; then
    input_files+=( "${file_longestcds_contamination_removal_fasta}" )
  fi

  seqkit stats \
  --all \
  --tabular \
  --threads "${NSLOTS}" \
  --out-file assembly_stat.tsv \
  "${input_files[@]}"

  if [[ -s assembly_stat.tsv ]]; then
    echo "The task ${task} was successful."
    mv_out assembly_stat.tsv "${file_assembly_stat}"
  fi
  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task='amalgkit quant'
disable_if_no_input_file "run_amalgkit_quant" "${file_amalgkit_metadata}"
if [[ ( ! -s "${file_amalgkit_merge_efflen}" || ! -s "${file_amalgkit_merge_count}" || ! -s "${file_amalgkit_merge_tpm}" ) && ${run_amalgkit_quant} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_amalgkit_quant}/${sp_ub}"
  if [[ ! -d "${dir_amalgkit_getfastq_sp}" ]]; then
    echo "amalgkit getfastq output directory not found: ${dir_amalgkit_getfastq_sp}. Exiting."
    exit 1
  fi
  if [[ -z "$(find "${dir_amalgkit_getfastq_sp}" -type f -name "*.amalgkit.fastq.gz" -print -quit 2>/dev/null)" ]]; then
    echo "No amalgkit getfastq FASTQ files were found in: ${dir_amalgkit_getfastq_sp}. Exiting."
    exit 1
  fi

  if [[ -e "./getfastq" ]]; then
    rm -rf -- "./getfastq"
  fi
  ln -s "${dir_amalgkit_getfastq_sp}" "./getfastq"
  recreate_dir "./fasta"
  recreate_dir "./index"
  recreate_dir "./quant"

  file_reference_fasta_link="./fasta/${sp_ub}_for_kallisto_index.fasta"
  if [[ ${kallisto_reference} == 'species_cds' ]]; then
    kallisto_ref_candidates=()
    mapfile -t kallisto_ref_candidates < <(find "${dir_pg_input}/species_cds" -maxdepth 1 -type f -name "${sp_ub}_*" | sort)
    if [[ ${#kallisto_ref_candidates[@]} -eq 0 ]]; then
      echo "No species_cds reference file matched ${sp_ub}_* in ${dir_pg_input}/species_cds. Exiting."
      exit 1
    fi
    if [[ ${#kallisto_ref_candidates[@]} -gt 1 ]]; then
      echo "Multiple species_cds reference files matched ${sp_ub}_*. Using the first one: ${kallisto_ref_candidates[0]}"
    fi
    file_kallisto_reference_fasta="${kallisto_ref_candidates[0]}"
  elif [[ ${kallisto_reference} == 'longest_transcript' ]]; then
    file_kallisto_reference_fasta=${file_longestcds_transcript}
  elif [[ ${kallisto_reference} == 'longest_cds' ]]; then
    file_kallisto_reference_fasta=${file_longestcds}
  elif [[ ${kallisto_reference} == 'contamination_removed_longest_cds' ]]; then
    file_kallisto_reference_fasta=${file_longestcds_contamination_removal_fasta}
  else
    echo "Please check the input parameter. kallisto_reference must not be: ${kallisto_reference}"
    exit 1
  fi

  echo "kallisto reference = ${kallisto_reference}: ${file_kallisto_reference_fasta}"
  if [[ -e "${file_kallisto_reference_fasta}" ]]; then
    ln -s "${file_kallisto_reference_fasta}" "${file_reference_fasta_link}"
  else
    echo "kallisto reference fasta file was not found in: ${file_kallisto_reference_fasta}"
    exit 1
  fi

  if amalgkit quant \
    --out_dir "./" \
    --threads "${NSLOTS}" \
    --metadata "${file_amalgkit_metadata}" \
    --clean_fastq no \
    --fasta_dir "./fasta" \
    --build_index yes; then
    exit_code_amalgkit_quant=0
  else
    exit_code_amalgkit_quant=$?
  fi

  if [[ ${exit_code_amalgkit_quant} -ne 0 ]]; then
    echo "amalgkit quant failed with exit code ${exit_code_amalgkit_quant}"
    exit 1
	  else
	    echo "amalgkit quant finished successfully"
	    if [[ ! -e "${dir_amalgkit_quant}/${sp_ub}" ]]; then
	      mkdir -p "${dir_amalgkit_quant}/${sp_ub}"
	    fi
	    shopt -s nullglob
	    quant_outputs=( ./quant/* )
	    shopt -u nullglob
	    if [[ ${#quant_outputs[@]} -eq 0 ]]; then
	      echo "amalgkit quant finished but no files were found in ./quant."
	      exit 1
	    fi
	    mv_out "${quant_outputs[@]}" "${dir_amalgkit_quant}/${sp_ub}"
	    rm -rf -- "./quant"
	    rm -f -- "./getfastq" # Do not put -r, otherwise the original getfastq files will be deleted.
	  fi

  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task='amalgkit merge'
disable_if_no_input_file "run_amalgkit_merge" "${file_amalgkit_metadata}"
if [[ ( ! -s "${file_amalgkit_merge_efflen}" || ! -s "${file_amalgkit_merge_count}" || ! -s "${file_amalgkit_merge_tpm}" || ! -s "${file_amalgkit_merge_metadata}" ) && ${run_amalgkit_merge} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ -e ./quant ]]; then
    rm -rf -- ./quant
  fi
  ln -s "${dir_amalgkit_quant}/${sp_ub}" "./quant"

  amalgkit merge \
  --out_dir "./" \
  --metadata "${file_amalgkit_metadata}"

  #sp_metadata=$(python -c "import pandas; d=pandas.read_csv('${file_amalgkit_metadata}',sep='\t',header=0); print(d.at[0,'scientific_name'].replace(' ','_'))")
  if [[ -s "./merge/${sp_ub}/${sp_ub}_eff_length.tsv" ]]; then
    echo "Copying amalgkit merge outputs from: ./merge/${sp_ub}"
    mv_out "./merge/${sp_ub}" "$(dirname "$(dirname "${file_amalgkit_merge_tpm}")")"
    mv_out "./merge/metadata.tsv" "${file_amalgkit_merge_metadata}"
    rm -rf -- "./merge"
    rm -f -- "./quant" # Do not put -r, otherwise the original quant files will be deleted.
  else
    echo "amalgkit merge outputs were not found in: ./merge/${sp_ub}"
  fi
  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task='Multispecies summary'
if is_output_older_than_inputs "^file_" "${file_multispecies_summary}"; then
  summary_flag=0
else
  summary_flag=$?
fi
if [[ ${run_multispecies_summary} -eq 1 && ${summary_flag} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "$(dirname "${file_multispecies_summary}")"
  cd "$(dirname "${file_multispecies_summary}")"

  python "${dir_script}/collect_common_BUSCO_genes.py" \
  --busco_outdir "$(dirname "${file_busco_full_longest_cds}")" \
  --ncpu "${NSLOTS}" \
  --outfile busco_table.tsv

  for dir_busco in \
    "$(dirname "${file_busco_full_cdna_isoforms}")" \
    "$(dirname "${file_busco_full_longest_cds}")" \
    "$(dirname "${file_busco_full_longest_cds_filtered}")"
  do
    if [[ ! -d "${dir_busco}" || -z "$(find "${dir_busco}" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
      echo "Skipping. No BUSCO output was found in: ${dir_busco}"
      continue
    fi
    Rscript "${dir_script}/annotation_summary.r" \
    --dir_species_cds_busco="${dir_busco}" \
    --tree_annotation_dir="${dir_script}/tree_annotation" \
    --min_og_species='auto'
    mv_out "annotation_summary.tsv" "$(basename "${dir_busco}").tsv"
    mv_out "busco_cds.svg" "$(basename "${dir_busco}").svg"
    mv_out "busco_cds.pdf" "$(basename "${dir_busco}").pdf"
  done

  Rscript "${dir_script}/multispecies_transcriptome_summary.r" \
  --dir_assembly_stat="$(dirname "${file_assembly_stat}")" \
  --dir_amalgkit_metadata="${dir_amalgkit_metadata}" \
  --dir_amalgkit_merge="$(dirname "$(dirname "${file_amalgkit_merge_tpm}")")" \
  --dir_busco_isoform="$(dirname "${file_busco_full_cdna_isoforms}")" \
  --dir_busco_longest_cds="$(dirname "${file_busco_full_longest_cds}")" \
  --dir_myscript="${dir_script}"

  if [[ -e "Rplots.pdf" ]]; then
    rm -f -- "Rplots.pdf"
  fi
  cd "${dir_tmp}"

  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

if [[ ${remove_amalgkit_fastq_after_completion} -eq 1 && $(is_fastq_requiring_downstream_analysis_done) -eq 1 ]]; then
  echo "remove_amalgkit_fastq_after_completion=1: All necessary output files were detected. amalgkit getfastq outputs will be removed."
  if [[ -e "${dir_amalgkit_getfastq_sp}" ]]; then
    rm -rf -- "${dir_amalgkit_getfastq_sp}"
    ensure_parent_dir "${file_amalgkit_getfastq_safely_removed_flag}"
    echo "Fastq files for this species have been safely removed." > "tmp.amalgkit_getfastq_removed.flag.txt"
    mv_out "tmp.amalgkit_getfastq_removed.flag.txt" "${file_amalgkit_getfastq_safely_removed_flag}"
  fi
else
  echo "fastp fastq files will not be removed."
fi

remove_empty_subdirs "${dir_transcriptome_assembly_output}"
if [[ ${delete_tmp_dir} -eq 1 ]]; then
  echo "delete_tmp_dir=1: Normal completion. Deleting ${dir_tmp}"
  if [[ -n "${dir_tmp:-}" && "${dir_tmp}" != "/" ]]; then
    rm -rf -- "${dir_tmp}"
  else
    echo "Refusing to delete unsafe tmp directory: '${dir_tmp:-}'"
  fi
else
  echo "Tmp directory will not be deleted: ${dir_tmp}"
fi

echo "$(date): Exiting Singularity environment"

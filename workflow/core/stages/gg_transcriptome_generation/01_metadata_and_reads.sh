# shellcheck shell=bash
# Sourced by gg_transcriptome_generation_core.sh.

task="amalgkit metadata/integrate"
if [[ ! -s "${file_amalgkit_metadata}" && ${run_amalgkit_metadata_or_integrate} -eq 1 ]]; then
  gg_step_start "${task}"
  if [[ -e "./metadata" ]]; then
    rm -rf -- "./metadata"
  fi

  if [[ "${selected_transcriptome_mode}" == "sraid" ]]; then
    search_string=$(build_entrez_or_search_string_from_file "${file_input_sra_list}")
    if [[ -n "${amalgkit_sra_strategy_query}" ]]; then
      search_string="${search_string} AND (${amalgkit_sra_strategy_query})"
    fi
    echo "Entrez strategy filter: ${amalgkit_sra_strategy_query:-disabled}"
    if ! run_amalgkit_metadata_query "${search_string}"; then
      echo "amalgkit metadata did not safely finish. Exiting."
      exit 1
    fi

    sp_space="${sp_ub//_/ }"
    extract_sraid_metadata_rows_for_species "./metadata/metadata.tsv" "${sp_space}" "./metadata.tsv"
    if [[ -n "${amalgkit_sra_strategy_query}" ]]; then
      missing_requested_sra_ids=()
      mapfile -t missing_requested_sra_ids < <(extract_requested_accessions_missing_from_metadata "./metadata.tsv" "${file_input_sra_list}")
      if [[ ${#missing_requested_sra_ids[@]} -gt 0 ]]; then
        missing_accessions_csv=$(csv_join_from_array "${missing_requested_sra_ids[@]}")
        echo "Strategy-filtered metadata lookup is missing ${#missing_requested_sra_ids[@]} requested accession(s): ${missing_accessions_csv}"
        echo "Retrying the missing accessions without the Entrez strategy filter and keeping transcriptomic rows such as lib_source=TRANSCRIPTOMIC even when lib_strategy=OTHER."
        printf '%s\n' "${missing_requested_sra_ids[@]}" > "./metadata_missing_accessions.txt"
        relaxed_search_string=$(build_entrez_or_search_string_from_file "./metadata_missing_accessions.txt")
        rm -rf -- "./metadata"
        if ! run_amalgkit_metadata_query "${relaxed_search_string}"; then
          echo "Relaxed accession-driven amalgkit metadata retry did not safely finish. Exiting."
          exit 1
        fi
        extract_transcriptomic_rows_for_requested_accessions "./metadata/metadata.tsv" "./metadata_missing_accessions.txt" "./metadata.relaxed.tsv"
        if metadata_table_has_data_rows "./metadata.relaxed.tsv"; then
          merge_metadata_tables_by_run "./metadata.tsv" "./metadata.relaxed.tsv" "./metadata.merged.tsv"
          mv_out "./metadata.merged.tsv" "./metadata.tsv"
          echo "Relaxed accession-driven metadata fallback retained $(($(wc -l < "./metadata.relaxed.tsv") - 1)) transcriptomic run(s)."
        else
          echo "Relaxed accession-driven metadata fallback did not recover transcriptomic metadata for the missing accessions."
        fi
      fi
    fi
    if ! metadata_table_has_data_rows "./metadata.tsv"; then
      echo "No metadata rows matched species '${sp_space}' in ./metadata/metadata.tsv. Exiting."
      echo "If the accession exists in SRA but was excluded by the Entrez strategy clause, relax amalgkit_sra_strategy_query or set it empty to disable strategy filtering."
      echo 'If the run is small-RNA/miRNA rather than transcriptome assembly input, prefer a different accession or provide explicit metadata/FASTQ inputs instead of mode_transcriptome_assembly="sraid".'
      exit 1
    fi

  elif [[ "${selected_transcriptome_mode}" == "fastq" ]]; then
    amalgkit integrate \
      --out_dir "./" \
      --download_dir "${dir_amalgkit_download_dir}" \
      --fastq_dir "${dir_species_fastq}" \
      --threads "${GG_TASK_CPUS}" \
      --remove_tmp yes

    repair_private_fastq_metadata_scientific_names \
      "./metadata_private_fastq.tsv" \
      "${sp_ub}" \
      "${gg_support_dir}"
    mv_out "./metadata_private_fastq.tsv" "./metadata.tsv"
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

ensure_parent_dir "${file_amalgkit_read_technology}"
ensure_parent_dir "${file_amalgkit_read_technology_summary_sh}"
if [[ -e "${file_amalgkit_read_technology_summary_sh}" ]]; then
  rm -f -- "${file_amalgkit_read_technology_summary_sh}"
fi
detect_transcriptome_read_technology_from_metadata "${file_amalgkit_metadata}" "${file_amalgkit_read_technology}" "${file_amalgkit_read_technology_summary_sh}"
configure_transcriptome_runtime_from_detected_metadata

task="amalgkit getfastq"
amalgkit_fastq_files=()
if [[ -d "${dir_amalgkit_getfastq_sp}" ]]; then
  mapfile -t amalgkit_fastq_files < <(find "${dir_amalgkit_getfastq_sp}" -type f -name "*.amalgkit.fastq.gz" | sort)
fi
echo "Number of amalgkit getfastq fastq files: ${#amalgkit_fastq_files[@]}"
echo "is_fastq_requiring_downstream_analysis_done: $(is_fastq_requiring_downstream_analysis_done)"
if [[ (${#amalgkit_fastq_files[@]} -eq 0 && ${run_amalgkit_getfastq} -eq 1) && $(is_fastq_requiring_downstream_analysis_done) -eq 0 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_amalgkit_getfastq_sp}"

  clear_getfastq_safely_removed_markers

  if [[ "${amalgkit_contam_filter}" == "yes" ]]; then
    if ! awk -F '\t' 'NR==1 {found=0; for (i=1; i<=NF; i++) if ($i=="taxid") found=1; exit(found ? 0 : 1)}' "${file_amalgkit_metadata}"; then
      echo "amalgkit_contam_filter=yes requires a taxid column in metadata: ${file_amalgkit_metadata}. Exiting."
      exit 1
    fi
  fi
  if ! run_amalgkit_getfastq_or_fallback; then
    exit 1
  fi
  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

# shellcheck shell=bash
# Sourced by gg_transcriptome_generation_core.sh.

task='amalgkit quant'
disable_if_no_input_file "run_amalgkit_quant" "${file_amalgkit_metadata}"
if [[ (! -s "${file_amalgkit_merge_efflen}" || ! -s "${file_amalgkit_merge_count}" || ! -s "${file_amalgkit_merge_tpm}") && ${run_amalgkit_quant} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_amalgkit_quant}/${sp_ub}"
  if [[ ! -d "${dir_amalgkit_getfastq_sp}" ]]; then
    echo "amalgkit getfastq output directory not found: ${dir_amalgkit_getfastq_sp}. Exiting."
    exit 1
  fi
  if [[ -z "$(find "${dir_amalgkit_getfastq_sp}" -type f -name "*.amalgkit.fastq.gz" -print -quit 2> /dev/null)" ]]; then
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
  if [[ ${kallisto_reference} == 'species_cds' ]]; then
    kallisto_ref_candidates=()
    mapfile -t kallisto_ref_candidates < <(gg_find_species_files_by_label "${gg_workspace_input_dir}/species_cds" "${sp_ub}")
    if [[ ${#kallisto_ref_candidates[@]} -eq 0 ]]; then
      echo "No species_cds reference file matched species label ${sp_ub} in ${gg_workspace_input_dir}/species_cds. Exiting."
      exit 1
    fi
    if [[ ${#kallisto_ref_candidates[@]} -gt 1 ]]; then
      echo "Multiple species_cds reference files matched species label ${sp_ub}. Exiting."
      printf '  %s\n' "${kallisto_ref_candidates[@]}"
      exit 1
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
    while IFS= read -r alias_name; do
      [[ -n "${alias_name}" ]] || continue
      echo "Staged exact quant reference from metadata scientific_name: ./fasta/${alias_name}"
    done < <(
      stage_quant_reference_fasta_aliases \
        "${file_amalgkit_metadata}" \
        "${file_kallisto_reference_fasta}" \
        "./fasta" \
        "${sp_ub}" \
        "${gg_support_dir}"
    )
  else
    echo "kallisto reference fasta file was not found in: ${file_kallisto_reference_fasta}"
    exit 1
  fi

  amalgkit_quant_cmd=(
    amalgkit quant
    --out_dir "./"
    --threads "${GG_TASK_CPUS}"
    --metadata "${file_amalgkit_metadata}"
    --clean_fastq no
    --fasta_dir "./fasta"
    --build_index yes
    --quant_backend "${amalgkit_quant_backend}"
    --oarfish_seq_tech "${amalgkit_oarfish_seq_tech}"
  )
  if [[ -n "${amalgkit_oarfish_options}" ]]; then
    amalgkit_quant_cmd+=(--oarfish_options "${amalgkit_oarfish_options}")
  fi

  if "${amalgkit_quant_cmd[@]}"; then
    exit_code_amalgkit_quant=0
  else
    exit_code_amalgkit_quant=$?
  fi

  if [[ ${exit_code_amalgkit_quant} -ne 0 ]]; then
    echo "amalgkit quant failed with exit code ${exit_code_amalgkit_quant}"
    exit 1
  else
    echo "amalgkit quant finished successfully"
    shopt -s nullglob
    quant_outputs=(./quant/*)
    shopt -u nullglob
    if [[ ${#quant_outputs[@]} -eq 0 ]]; then
      echo "amalgkit quant finished but no files were found in ./quant."
      exit 1
    fi
    mv_out_replace_dir "./quant" "${dir_amalgkit_quant}/${sp_ub}"
    rm -rf -- "./quant"
    rm -f -- "./getfastq" # Do not put -r, otherwise the original getfastq files will be deleted.
  fi

  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task='amalgkit merge'
disable_if_no_input_file "run_amalgkit_merge" "${file_amalgkit_metadata}"
if [[ (! -s "${file_amalgkit_merge_efflen}" || ! -s "${file_amalgkit_merge_count}" || ! -s "${file_amalgkit_merge_tpm}" || ! -s "${file_amalgkit_merge_metadata}") && ${run_amalgkit_merge} -eq 1 ]]; then
  gg_step_start "${task}"
  merge_output_prefix=""
  merge_output_dir=""
  merge_metadata_file="./metadata.merge_species_level.tsv"

  if [[ -e ./quant ]]; then
    rm -rf -- ./quant
  fi
  ln -s "${dir_amalgkit_quant}/${sp_ub}" "./quant"
  stage_amalgkit_merge_metadata_for_species \
    "${file_amalgkit_metadata}" \
    "${merge_metadata_file}" \
    "${sp_ub}" \
    "${gg_support_dir}"

  amalgkit merge \
    --out_dir "./" \
    --metadata "${merge_metadata_file}"

  if merge_output_prefix=$(resolve_amalgkit_merge_output_prefix \
    "${merge_metadata_file}" \
    "./merge" \
    "${sp_ub}"); then
    merge_output_dir="./merge/${merge_output_prefix}"
    echo "Copying amalgkit merge outputs from: ${merge_output_dir}"
    mv_out "${merge_output_dir}/${merge_output_prefix}_eff_length.tsv" "${file_amalgkit_merge_efflen}"
    mv_out "${merge_output_dir}/${merge_output_prefix}_est_counts.tsv" "${file_amalgkit_merge_count}"
    mv_out "${merge_output_dir}/${merge_output_prefix}_tpm.tsv" "${file_amalgkit_merge_tpm}"
    mv_out "./merge/metadata.tsv" "${file_amalgkit_merge_metadata}"
    rm -rf -- "./merge"
    rm -f -- "./quant" # Do not put -r, otherwise the original quant files will be deleted.
    rm -f -- "${merge_metadata_file}"
  else
    echo "amalgkit merge outputs were not found for the metadata scientific_name prefix under: ./merge"
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
  cd "$(dirname "${file_multispecies_summary}")" || exit 1

  python "${gg_support_dir}/collect_common_BUSCO_genes.py" \
    --busco_outdir "$(dirname "${file_busco_full_longest_cds}")" \
    --ncpu "${GG_TASK_CPUS}" \
    --outfile busco_table.tsv

  for dir_busco in \
    "$(dirname "${file_busco_full_cdna_isoforms}")" \
    "$(dirname "${file_busco_full_longest_cds}")" \
    "$(dirname "${file_busco_full_longest_cds_filtered}")"; do
    if [[ ! -d "${dir_busco}" || -z "$(find "${dir_busco}" -mindepth 1 -print -quit 2> /dev/null)" ]]; then
      echo "Skipping. No BUSCO output was found in: ${dir_busco}"
      continue
    fi
    Rscript "${gg_support_dir}/annotation_summary.r" \
      --dir_species_cds_busco="${dir_busco}" \
      --min_og_species='auto'
    mv_out "annotation_summary.tsv" "$(basename "${dir_busco}").tsv"
    mv_out "busco_cds.svg" "$(basename "${dir_busco}").svg"
    mv_out "busco_cds.pdf" "$(basename "${dir_busco}").pdf"
  done

  Rscript "${gg_support_dir}/multispecies_transcriptome_summary.r" \
    --dir_assembly_stat="$(dirname "${file_assembly_stat}")" \
    --dir_amalgkit_metadata="${dir_input_amalgkit_metadata}" \
    --dir_amalgkit_merge="$(dirname "$(dirname "${file_amalgkit_merge_tpm}")")" \
    --dir_busco_isoform="$(dirname "${file_busco_full_cdna_isoforms}")" \
    --dir_busco_longest_cds="$(dirname "${file_busco_full_longest_cds}")" \
    --dir_myscript="${gg_support_dir}"

  if [[ -e "Rplots.pdf" ]]; then
    rm -f -- "Rplots.pdf"
  fi
  cd "${dir_tmp}" || exit 1

  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

if [[ ${remove_amalgkit_fastq_after_completion} -eq 1 && $(is_fastq_requiring_downstream_analysis_done) -eq 1 ]]; then
  echo "remove_amalgkit_fastq_after_completion=1: All necessary output files were detected. amalgkit getfastq outputs will be removed."
  if [[ -e "${dir_amalgkit_getfastq_sp}" ]]; then
    safe_delete_getfastq_fastq_files
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

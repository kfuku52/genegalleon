# shellcheck shell=bash
# Sourced by gg_genome_evolution_core.sh.

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
      if gg_orthogroup_file_matches_id "${file_name}" "${og_id}"; then
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

if [[ ${run_cafe} -eq 1 ]]; then
  if [[ -s "${file_orthogroup_copy_number}" ]]; then
    disable_if_no_input_file "run_cafe" "${file_orthogroup_copy_number}" "${file_dated_species_tree}"
  else
    disable_if_no_input_file "run_cafe" "${file_orthogroup_genecount_selected}" "${file_dated_species_tree}"
  fi
fi
if [[ ${run_orthogroup_copy_number_trait_pgls} -eq 1 ]]; then
  if [[ -s "${file_orthogroup_copy_number}" ]]; then
    disable_if_no_input_file "run_orthogroup_copy_number_trait_pgls" "${file_orthogroup_copy_number}" "${file_dated_species_tree}" "${file_trait}"
  else
    disable_if_no_input_file "run_orthogroup_copy_number_trait_pgls" "${file_orthogroup_genecount_selected}" "${file_dated_species_tree}" "${file_trait}"
  fi
fi

task="Orthogroup copy-number matrix preparation"
if [[ ! -s "${file_orthogroup_copy_number}" && ( ${run_cafe} -eq 1 || ${run_orthogroup_copy_number_trait_pgls} -eq 1 ) ]]; then
  gg_step_start "${task}"
  if [[ -d "${dir_orthogroup_copy_number}" ]]; then
    rm -rf -- "${dir_orthogroup_copy_number}"
  fi
  python "${gg_support_dir}/prepare_orthogroup_copy_number.py" \
    --genecount "${file_orthogroup_genecount_selected}" \
    --dated_species_tree "${file_dated_species_tree}" \
    --output_dir "${dir_orthogroup_copy_number}" \
    --max_size_differential "${orthogroup_copy_number_max_size_differential}"
  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task='CAFE analysis'
if [[ (! -s "${file_cafe_summary_all_pdf}" || ! -s "${file_cafe_summary_significant_pdf}") && ${run_cafe} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_cafe}"

  if [[ ! -s "${dir_cafe_output}/Gamma_asr.tre" || ! -s "${dir_cafe_output}/Gamma_count.tab" || ! -s "${dir_cafe_output}/Gamma_change.tab" ]]; then
    if [[ -d "${dir_cafe_output}" ]]; then
      rm -rf -- "${dir_cafe_output}"
    fi
    cafe5 \
      --infile "${file_orthogroup_copy_number}" \
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

task="Orthogroup copy-number trait PGLS"
disable_if_no_input_file "run_orthogroup_copy_number_trait_pgls" "${file_orthogroup_copy_number}" "${file_dated_species_tree}" "${file_trait}"
if [[ (! -s "${file_orthogroup_copy_number_matrix}" || ! -s "${file_orthogroup_copy_number_trait_pgls}" || ! -s "${file_orthogroup_copy_number_trait_pgls_summary_pdf}") && ${run_orthogroup_copy_number_trait_pgls} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_orthogroup_copy_number_trait_pgls}"
  if ! Rscript "${gg_support_dir}/orthogroup_copy_number_trait_pgls.r" \
    --file_orthogroup_copy_number="${file_orthogroup_copy_number}" \
    --file_sptree="${file_dated_species_tree}" \
    --file_trait="${file_trait}" \
    --outdir="${dir_orthogroup_copy_number_trait_pgls}" \
    --trait="${orthogroup_copy_number_trait}" \
    --min_species="${orthogroup_copy_number_trait_min_species}" \
    --family_ids="${orthogroup_copy_number_trait_family_ids}" \
    --family_file="${orthogroup_copy_number_trait_family_file}" \
    --max_families="${orthogroup_copy_number_trait_max_families}" \
    --p_adjust_method="${orthogroup_copy_number_trait_p_adjust_method}" \
    --alpha="${orthogroup_copy_number_trait_alpha}" \
    --plot_top_n="${orthogroup_copy_number_trait_plot_top_n}" \
    --verbose=0; then
    echo "Error in Rscript orthogroup_copy_number_trait_pgls.r. Exiting."
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

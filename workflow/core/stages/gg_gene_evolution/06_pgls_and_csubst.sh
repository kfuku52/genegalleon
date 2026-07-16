# shellcheck shell=bash
# Sourced by gg_gene_evolution_core.sh.

task="Gene tree E-PGLS analysis"
gg_step_skip "${task}"

task="Species tree PGLS analysis"
disable_if_no_input_file "run_pgls_species_tree" "${file_sp_trait}" "${species_tree_pruned}" "${file_og_expression}"
if [[ ! -s "${file_og_species_pgls}" && ${run_pgls_species_tree} -eq 1 ]]; then
  gg_step_start "${task}"
  pgls_merge_replicates="yes"
  if [[ ${pgls_use_phenocov} -eq 1 ]]; then
    pgls_merge_replicates="no"
  fi

  Rscript "${gg_support_dir}/species_tree_pgls.r" \
    --file_sptree="${species_tree_pruned}" \
    --file_exp="${file_og_expression}" \
    --file_trait="${file_sp_trait}" \
    --replicate_sep="_" \
    --exp_value_type="${exp_value_type}" \
    --merge_replicates="${pgls_merge_replicates}" \
    2>&1 | tee pgls.log

  mv_out species_tree_PGLS.tsv "${file_og_species_pgls}"
  mv_out species_tree_PGLS.barplot.pdf "${file_og_species_pgls_plot}"
else
  gg_step_skip "${task}"
fi

task="IQ-TREE ancestral codon sequence reconstruction for CSUBST"
disable_if_no_input_file "run_iqtree_anc" "${file_og_trimmed_aln_analysis}" "${file_og_rooted_tree_analysis}"
if [[ ! -s "${file_og_iqtree_anc}" && ${run_iqtree_anc} -eq 1 ]]; then
  gg_step_start "${task}"

  shopt -s nullglob
  csubst_cleanup_paths=(csubst.* csubst_search csubst_scan)
  shopt -u nullglob
  if [[ ${#csubst_cleanup_paths[@]} -gt 0 ]]; then
    rm -rf -- "${csubst_cleanup_paths[@]}"
  fi
  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file "tmp.csubst.fasta"

  nwkit drop --target intnode --support yes --name no --infile "${file_og_rooted_tree_analysis}" |
    nwkit intersection --seqin tmp.csubst.fasta --seqout csubst.fasta |
    nwkit sanitize --remove_singleton yes --resolve_polytomy no \
      > csubst.nwk

  rm -f -- tmp.csubst.fasta

  build_iqtree_mem_args
  iqtree \
    -s csubst.fasta \
    -te csubst.nwk \
    -m "${codon_model}" \
    -T AUTO \
    --threads-max "${GG_TASK_CPUS}" \
    --seqtype "CODON${genetic_code}" \
    --prefix csubst \
    --ancestral \
    --rate \
    "${IQTREE_MEM_ARGS[@]}" \
    --seed 12345 \
    --redo

  if [[ -s csubst.rate && -s csubst.state && -s csubst.treefile ]]; then
    if [[ -e "${og_id}.iqtree.anc" ]]; then
      rm -rf -- "${og_id}.iqtree.anc"
    fi
    mkdir -p "${og_id}.iqtree.anc"
    shopt -s nullglob
    csubst_outputs=(csubst.*)
    shopt -u nullglob
    if [[ ${#csubst_outputs[@]} -eq 0 ]]; then
      echo "Expected csubst output files were not found."
      exit 1
    fi
    mv_out "${csubst_outputs[@]}" "${og_id}.iqtree.anc"
    zip -rq "${og_id}.iqtree.anc.zip" "${og_id}.iqtree.anc"
    mv_out "${og_id}.iqtree.anc.zip" "${file_og_iqtree_anc}"
    rm -rf -- "${og_id}.iqtree.anc"
  fi
else
  gg_step_skip "${task}"
fi

task="CSUBST"
disable_if_no_input_file "run_csubst" "${file_og_iqtree_anc}"
if [[ (! -s "${file_og_csubst_b}" || ! -s "${file_og_csubst_cb_stats}") && ${run_csubst} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ -s "${file_sp_trait}" ]]; then
    echo "CSUBST foreground specification file: ${file_sp_trait}"
    first_trait_header=""
    IFS= read -r first_trait_header < "${file_sp_trait}" || true
    if [[ "${first_trait_header}" == *" "* ]]; then
      echo "Column names should not contain spaces: ${file_sp_trait}"
      echo "Exiting."
      exit 1
    fi
    write_species_trait_foreground_regex_table "${file_sp_trait}" "foreground.tsv"
    foreground_params=(--foreground foreground.tsv --fg_format 2)
  else
    echo 'Foreground specification file was not found. CSUBST will run without it.'
    foreground_params=()
  fi
  shopt -s nullglob
  csubst_cleanup_paths=(csubst.*)
  shopt -u nullglob
  if [[ ${#csubst_cleanup_paths[@]} -gt 0 ]]; then
    rm -rf -- "${csubst_cleanup_paths[@]}"
  fi
  if [[ -z "${og_id:-}" ]]; then
    echo "og_id is empty. Refusing to run cleanup pattern."
    exit 1
  fi
  shopt -s nullglob
  og_cleanup_paths=("${og_id}".*)
  shopt -u nullglob
  if [[ ${#og_cleanup_paths[@]} -gt 0 ]]; then
    rm -rf -- "${og_cleanup_paths[@]}"
  fi
  unzip -q "${file_og_iqtree_anc}"
  csubst_input_base="./${og_id}.iqtree.anc/csubst"
  csubst_search_dir="csubst_search"

  csubst search \
    --genetic_code "${genetic_code}" \
    --alignment_file "${csubst_input_base}.fasta" \
    --rooted_tree_file "${csubst_input_base}.nwk" \
    --iqtree_treefile "${csubst_input_base}.treefile" \
    --iqtree_state "${csubst_input_base}.state" \
    --iqtree_rate "${csubst_input_base}.rate" \
    --iqtree_iqtree "${csubst_input_base}.iqtree" \
    --iqtree_log "${csubst_input_base}.log" \
    --iqtree_model "${codon_model}" \
    --max_arity "${csubst_max_arity}" \
    --exhaustive_until "${csubst_exhaustive_until}" \
    --cutoff_stat "${csubst_cutoff_stat}" \
    --max_combination "${csubst_max_combination}" \
    --fg_exclude_wg "${csubst_fg_exclude_wg}" \
    --fg_stem_only "${csubst_fg_stem_only}" \
    --nonsyn_recode "${csubst_nonsyn_recode}" \
    --expectation_method "codon_model" \
    --threads "${GG_TASK_CPUS}" \
    --float_type 32 \
    "${foreground_params[@]}"

  csubst_b_src=""
  csubst_cb_stats_src=""
  for candidate in "${csubst_search_dir}/csubst_b.tsv" "csubst_b.tsv"; do
    if [[ -s "${candidate}" ]]; then
      csubst_b_src="${candidate}"
      break
    fi
  done
  for candidate in "${csubst_search_dir}/csubst_cb_stats.tsv" "csubst_cb_stats.tsv"; do
    if [[ -s "${candidate}" ]]; then
      csubst_cb_stats_src="${candidate}"
      break
    fi
  done

  if [[ -n "${csubst_b_src}" && -n "${csubst_cb_stats_src}" ]]; then
    echo "CSUBST was successful."
    mv_out "${csubst_b_src}" "${file_og_csubst_b}"
    mv_out "${csubst_cb_stats_src}" "${file_og_csubst_cb_stats}"
    if [[ ${csubst_max_arity} -gt 2 ]]; then
      for ((i = 2; i <= csubst_max_arity; i++)); do
        csubst_cb_src=""
        for candidate in "${csubst_search_dir}/csubst_cb_${i}.tsv" "csubst_cb_${i}.tsv"; do
          if [[ -e "${candidate}" ]]; then
            csubst_cb_src="${candidate}"
            break
          fi
        done
        if [[ -n "${csubst_cb_src}" ]]; then
          my_csubst_file=file_og_csubst_cb_${i}
          mv_out "${csubst_cb_src}" "${!my_csubst_file}"
        fi
      done
    fi
  else
    echo "CSUBST failed."
  fi
else
  gg_step_skip "${task}"
fi

task="CSUBST scan"
disable_if_no_input_file "run_csubst_scan" "${file_og_iqtree_anc}" "${file_sp_trait}"
if [[ (! -s "${file_og_csubst_scan}" || ! -s "${file_og_csubst_scan_units}") && ${run_csubst_scan} -eq 1 ]]; then
  gg_step_start "${task}"

  echo "CSUBST scan foreground specification file: ${file_sp_trait}"
  first_trait_header=""
  IFS= read -r first_trait_header < "${file_sp_trait}" || true
  if [[ "${first_trait_header}" == *" "* ]]; then
    echo "Column names should not contain spaces: ${file_sp_trait}"
    echo "Exiting."
    exit 1
  fi
  write_species_trait_foreground_regex_table "${file_sp_trait}" "foreground.tsv"

  if ! csubst scan -h >/dev/null 2>&1; then
    echo "csubst scan is unavailable in this runtime. Rebuild the GeneGalleon container with a csubst version that provides the scan subcommand."
    exit 1
  fi

  shopt -s nullglob
  csubst_cleanup_paths=(csubst.* csubst_scan)
  shopt -u nullglob
  if [[ ${#csubst_cleanup_paths[@]} -gt 0 ]]; then
    rm -rf -- "${csubst_cleanup_paths[@]}"
  fi
  if [[ -e "${og_id}.iqtree.anc" ]]; then
    rm -rf -- "${og_id}.iqtree.anc"
  fi
  unzip -q "${file_og_iqtree_anc}"
  csubst_input_base="./${og_id}.iqtree.anc/csubst"
  csubst_scan_dir="csubst_scan"

  csubst scan \
    --genetic_code "${genetic_code}" \
    --alignment_file "${csubst_input_base}.fasta" \
    --rooted_tree_file "${csubst_input_base}.nwk" \
    --iqtree_treefile "${csubst_input_base}.treefile" \
    --iqtree_state "${csubst_input_base}.state" \
    --iqtree_rate "${csubst_input_base}.rate" \
    --iqtree_iqtree "${csubst_input_base}.iqtree" \
    --iqtree_log "${csubst_input_base}.log" \
    --iqtree_model "${codon_model}" \
    --foreground foreground.tsv \
    --fg_format 2 \
    --scan_match "${csubst_scan_match}" \
    --scan_min_event_pp "${csubst_scan_min_event_pp}" \
    --scan_min_support "${csubst_scan_min_support}" \
    --scan_rate_event_mode "${csubst_scan_rate_event_mode}" \
    --scan_rate_length "${csubst_scan_rate_length}" \
    --scan_rate_exposure "${csubst_scan_rate_exposure}" \
    --scan_other_scope "${csubst_scan_other_scope}" \
    --scan_pvalue_calibration "${csubst_scan_pvalue_calibration}" \
    --scan_n_permutations "${csubst_scan_n_permutations}" \
    --scan_site_plot "${csubst_scan_site_plot}" \
    --tree_site_plot_format "${csubst_scan_tree_site_plot_format}" \
    --tree_site_plot_max_sites "${csubst_scan_tree_site_plot_max_sites}" \
    --nonsyn_recode "${csubst_nonsyn_recode}" \
    --threads "${GG_TASK_CPUS}" \
    --float_type 32 \
    --outdir "${csubst_scan_dir}" \
    --output_prefix csubst

  if [[ -s "${csubst_scan_dir}/csubst_scan.tsv" && -s "${csubst_scan_dir}/csubst_scan_units.tsv" ]]; then
    echo "CSUBST scan was successful."
    mv_out "${csubst_scan_dir}/csubst_scan.tsv" "${file_og_csubst_scan}"
    mv_out "${csubst_scan_dir}/csubst_scan_units.tsv" "${file_og_csubst_scan_units}"
    shopt -s nullglob
    csubst_scan_foreground_branch_files=("${csubst_scan_dir}"/csubst_foreground_branch*.txt)
    shopt -u nullglob
    for csubst_scan_foreground_branch_file in "${csubst_scan_foreground_branch_files[@]}"; do
      csubst_scan_foreground_branch_name=$(basename "${csubst_scan_foreground_branch_file}")
      csubst_scan_foreground_branch_suffix="${csubst_scan_foreground_branch_name#csubst_foreground_branch}"
      mv_out \
        "${csubst_scan_foreground_branch_file}" \
        "${dir_output_active}/csubst_scan_foreground_branch/${og_id}_csubst_foreground_branch${csubst_scan_foreground_branch_suffix}"
    done
    if [[ -e "${csubst_scan_dir}/csubst_scan.tree_site.${csubst_scan_tree_site_plot_format}" ]]; then
      mv_out "${csubst_scan_dir}/csubst_scan.tree_site.${csubst_scan_tree_site_plot_format}" "${file_og_csubst_scan_plot}"
    fi
    if [[ -e "${csubst_scan_dir}/csubst.log" ]]; then
      mv_out "${csubst_scan_dir}/csubst.log" "${file_og_csubst_scan_log}"
    fi
  else
    echo "CSUBST scan failed."
  fi
  rm -rf -- "${og_id}.iqtree.anc" "${csubst_scan_dir}"
else
  gg_step_skip "${task}"
fi

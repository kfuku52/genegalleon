# shellcheck shell=bash
# Sourced by gg_gene_evolution_core.sh.

task="Expression matrix preparation"
disable_if_no_input_file "run_generate_expression_matrix" "${file_og_trimmed_aln_analysis}"
if [[ ! -s "${file_og_expression}" && ${run_generate_expression_matrix} -eq 1 ]]; then
  gg_step_start "${task}"
  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file "${og_id}.trait_matrix_input.fasta"

  python "${gg_support_dir}/get_trait_matrix.py" \
    --dir_trait "${dir_sp_expression}" \
    --seqfile "${og_id}.trait_matrix_input.fasta" \
    --ncpu "${GG_TASK_CPUS}" \
    --outfile expression_matrix.tsv
  rm -f -- "${og_id}.trait_matrix_input.fasta"
  if [[ -s expression_matrix.tsv ]]; then
    mv_out expression_matrix.tsv "${file_og_expression}"
  fi
else
  gg_step_skip "${task}"
fi

task="Promoter fasta generation"
disable_if_no_input_file "run_extract_promoter_fasta" "${file_og_gff_info}"
if [[ ! -s "${file_og_promoter_fasta}" && ${run_extract_promoter_fasta} -eq 1 ]]; then
  gg_step_start "${task}"

  python "${gg_support_dir}/get_promoter_fasta.py" \
    --dir_genome "${dir_sp_genome}" \
    --geneinfo_tsv "${file_og_gff_info}" \
    --seqkit_exe "seqkit" \
    --outfile "${og_id}.promoter.tmp.fasta" \
    --promoter_bp "${promoter_bp}" \
    --ncpu "${GG_TASK_CPUS}"
  if [[ -s "${og_id}.promoter.tmp.fasta" ]]; then
    seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.promoter.tmp.fasta" --out-file "${og_id}.promoter.out.fa.gz"
    mv_out "${og_id}.promoter.out.fa.gz" "${file_og_promoter_fasta}"
    rm -f -- "${og_id}.promoter.tmp.fasta"
  fi
else
  gg_step_skip "${task}"
fi

task="fimo"
jaspar_path=""
if [[ ${run_fimo} -eq 1 ]]; then
  if ! jaspar_path=$(ensure_jaspar_file "${gg_workspace_dir}" "${jaspar_file}") || [[ -z "${jaspar_path}" ]]; then
    echo "Failed to prepare JASPAR motif file (${jaspar_file}). Exiting."
    exit 1
  fi
fi
disable_if_no_input_file "run_fimo" "${file_og_promoter_fasta}" "${jaspar_path}"
if [[ ! -s "${file_og_fimo}" && ${run_fimo} -eq 1 ]]; then
  gg_step_start "${task}"
  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_promoter_fasta}" --out-file "${og_id}.fimo.input.fasta"

  fimo \
    --oc "fimo_out" \
    "${jaspar_path}" \
    "${og_id}.fimo.input.fasta"
  rm -f -- "${og_id}.fimo.input.fasta"

  fimo_result_table=""
  if [[ -s "./fimo_out/fimo.tsv" ]]; then
    fimo_result_table="./fimo_out/fimo.tsv"
  elif [[ -s "./fimo_out/fimo.txt" ]]; then
    fimo_result_table="./fimo_out/fimo.txt"
  fi

  if [[ -n "${fimo_result_table}" ]]; then
    mv_out "${fimo_result_table}" "${file_og_fimo}"
    rm -rf -- "./fimo_out"
  else
    echo "FIMO result table was not detected (expected fimo.tsv or fimo.txt). Keeping fimo_out for inspection."
  fi
else
  gg_step_skip "${task}"
fi

task="Tree pruning"
disable_if_no_input_file "run_tree_pruning" "${file_og_expression}" "${file_og_untrimmed_aln_analysis}" "${file_og_trimmed_aln_analysis}" "${file_og_unrooted_tree_analysis}" "${file_og_rooted_tree_analysis}"
if [[ (! -s "${file_og_untrimmed_aln_pruned}" || ! -s "${file_og_trimmed_aln_pruned}" || ! -s "${file_og_unrooted_tree_pruned}" || ! -s "${file_og_rooted_tree_pruned}") ]]; then
  is_all_outputs_exist=0
else
  if [[ ${run_tree_dating} -eq 1 && ! -s "${file_og_dated_tree_pruned}" ]]; then
    is_all_outputs_exist=0
  else
    is_all_outputs_exist=1
  fi
fi
if [[ ${is_all_outputs_exist} -eq 0 && ${run_tree_pruning} -eq 1 ]]; then
  gg_step_start "${task}"

  cut -f 1 "${file_og_expression}" | tail -n +2 > target_genes.txt

  if [[ -s "${file_og_untrimmed_aln_analysis}" ]]; then
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_untrimmed_aln_analysis}" --out-file "${og_id}.untrimmed.input.fasta"
    awk '
    NR==FNR {
      sub(/\r$/, "", $1)
      if ($1 != "") {
        keep[$1] = 1
      }
      next
    }
    /^>/ {
      id = substr($0, 2)
      sub(/[[:space:]].*$/, "", id)
      keep_seq = (id in keep)
    }
    keep_seq {
      print
    }
    ' target_genes.txt "${og_id}.untrimmed.input.fasta" > "${og_id}.untrimmed.pruned.tmp.fasta"
    rm -f -- "${og_id}.untrimmed.input.fasta"
    seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.untrimmed.pruned.tmp.fasta" --out-file "${og_id}.untrimmed.pruned.out.fa.gz"
    mv_out "${og_id}.untrimmed.pruned.out.fa.gz" "${file_og_untrimmed_aln_pruned}"
    rm -f -- "${og_id}.untrimmed.pruned.tmp.fasta"
  fi

  if [[ -s "${file_og_trimmed_aln_analysis}" ]]; then
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file "${og_id}.trimmed.input.fasta"
    awk '
    NR==FNR {
      sub(/\r$/, "", $1)
      if ($1 != "") {
        keep[$1] = 1
      }
      next
    }
    /^>/ {
      id = substr($0, 2)
      sub(/[[:space:]].*$/, "", id)
      keep_seq = (id in keep)
    }
    keep_seq {
      print
    }
    ' target_genes.txt "${og_id}.trimmed.input.fasta" > "${og_id}.trimmed.pruned.tmp.fasta"
    rm -f -- "${og_id}.trimmed.input.fasta"
    seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.trimmed.pruned.tmp.fasta" --out-file "${og_id}.trimmed.pruned.out.fa.gz"
    mv_out "${og_id}.trimmed.pruned.out.fa.gz" "${file_og_trimmed_aln_pruned}"
    rm -f -- "${og_id}.trimmed.pruned.tmp.fasta"
  fi

  mapfile -t prune_genes < <(sed -e '/^[[:space:]]*$/d' target_genes.txt)
  prune_pattern=""
  if [[ ${#prune_genes[@]} -gt 0 ]]; then
    prune_pattern=$(
      printf '%s\n' "${prune_genes[@]}" |
        sed -e 's/[][(){}.^$+*?|\\-]/\\&/g' |
        paste -sd'|' -
    )
  fi

  if [[ -s "${file_og_unrooted_tree_analysis}" ]]; then
    if [[ ${#prune_genes[@]} -eq 0 ]]; then
      cat "${file_og_unrooted_tree_analysis}"
    else
      nwkit prune \
        --infile "${file_og_unrooted_tree_analysis}" \
        --pattern "^(${prune_pattern})$" \
        --invert_match yes
    fi |
      nwkit drop --target root --length yes |
      nwkit sanitize --remove_singleton yes --resolve_polytomy no \
        > "${og_id}.unrooted.pruned.tmp.nwk"
    mv_out "${og_id}.unrooted.pruned.tmp.nwk" "${file_og_unrooted_tree_pruned}"
  fi

  if [[ -s "${file_og_rooted_tree_analysis}" ]]; then
    if [[ ${#prune_genes[@]} -eq 0 ]]; then
      cat "${file_og_rooted_tree_analysis}"
    else
      nwkit prune \
        --infile "${file_og_rooted_tree_analysis}" \
        --pattern "^(${prune_pattern})$" \
        --invert_match yes
    fi |
      nwkit drop --target root --length yes |
      nwkit sanitize --remove_singleton yes --resolve_polytomy no \
        > "${og_id}.rooted.pruned.tmp.nwk"
    mv_out "${og_id}.rooted.pruned.tmp.nwk" "${file_og_rooted_tree_pruned}"
  fi

  if [[ -s "${file_og_dated_tree_analysis}" ]]; then
    if [[ ${#prune_genes[@]} -eq 0 ]]; then
      cat "${file_og_dated_tree_analysis}"
    else
      nwkit prune \
        --infile "${file_og_dated_tree_analysis}" \
        --pattern "^(${prune_pattern})$" \
        --invert_match yes
    fi |
      nwkit drop --target root --length yes |
      nwkit sanitize --remove_singleton yes --resolve_polytomy no \
        > "${og_id}.dated.pruned.tmp.nwk"
    mv_out "${og_id}.dated.pruned.tmp.nwk" "${file_og_dated_tree_pruned}"
  fi
else
  gg_step_skip "${task}"
fi
if [[ ${run_tree_pruning} -eq 1 ]]; then
  num_gene_before_pruning=$(gg_count_fasta_records "${file_og_trimmed_aln_analysis}")
  num_gene_after_pruning=$(gg_count_fasta_records "${file_og_trimmed_aln_pruned}")
  echo "Number of genes before pruning: ${num_gene_before_pruning}"
  echo "Number of genes after pruning: ${num_gene_after_pruning}"
  if [[ ${num_gene_after_pruning} -lt 3 ]]; then
    echo 'This is not sufficient for tree-based analysis (<3). Exiting.'
    exit 0
  fi
  set_analysis_file untrimmed_aln "${file_og_untrimmed_aln_pruned}"
  set_analysis_file trimmed_aln "${file_og_trimmed_aln_pruned}"
  set_analysis_file unrooted_tree "${file_og_unrooted_tree_pruned}"
  set_analysis_file rooted_tree "${file_og_rooted_tree_pruned}"
  set_analysis_file dated_tree "${file_og_dated_tree_pruned}"
  if [[ -s "${file_og_rooted_tree_analysis}" ]]; then
    assert_strictly_bifurcating_tree "${file_og_rooted_tree_analysis}" "Rooted analysis tree"
  fi
  if [[ -s "${file_og_dated_tree_analysis}" ]]; then
    assert_strictly_bifurcating_tree "${file_og_dated_tree_analysis}" "Dated analysis tree"
  fi
fi
if [[ -s "${file_og_expression}" && ${run_l1ou} -eq 1 ]]; then
  # This block should be run after tree pruning.
  num_gene_trait=$(($(wc -l < "${file_og_expression}") - 1)) # -1 for header
  num_gene_tree=$(gg_count_fasta_records "${file_og_trimmed_aln_analysis}")
  if [[ ${num_gene_trait} -eq ${num_gene_tree} ]]; then
    echo "num_gene_trait (${num_gene_trait}) and num_gene_tree (${num_gene_tree}) matched."
  else
    echo "num_gene_trait (${num_gene_trait}) and num_gene_tree (${num_gene_tree}) did not match."
    if [[ ${run_tree_pruning} -ne 1 && ${run_l1ou} -eq 1 ]]; then
      echo "Set run_tree_pruning=1 to run phylogenetic comparative analysis. Exiting."
      exit 1
    fi
  fi
fi

task="Checking whether downstream output files are correctly unpruned/pruned."
if [[ ${check_pruned} -eq 1 ]]; then
  gg_step_start "${task}"
  file_to_remove=(
    "${file_og_mapdnds_dn}"
    "${file_og_mapdnds_ds}"
    "${file_og_hyphy_dnds}"
    "${file_og_gff_info}"
    "${file_og_scm_intron_summary}"
    "${file_og_scm_intron_plot}"
    "${file_og_l1ou_fit_rdata}"
    "${file_og_l1ou_fit_conv_rdata}"
    "${file_og_l1ou_fit_tree}"
    "${file_og_l1ou_fit_regime}"
    "${file_og_l1ou_fit_leaf}"
    "${file_og_l1ou_fit_plot}"
    "${file_og_iqtree_anc}"
    "${file_og_csubst_b}"
    "${file_og_csubst_cb_stats}"
    "${file_og_csubst_scan}"
    "${file_og_csubst_scan_units}"
    "${file_og_csubst_scan_foreground_branch}"
    "${file_og_csubst_scan_plot}"
    "${file_og_csubst_scan_log}"
    "${file_og_gene_pgls}"
    "${file_og_gene_pgls_plot}"
    "${file_og_species_pgls}"
    "${file_og_species_pgls_plot}"
  )
  for ((i = 2; i <= csubst_max_arity; i++)); do
    varname="file_og_csubst_cb_${i}"
    if [[ -n "${!varname:-}" ]]; then
      file_to_remove+=("${!varname}")
    fi
  done
  remove_flag=0
  if [[ ${run_tree_pruning} -eq 0 ]]; then
    if [[ -s "${file_og_untrimmed_aln_pruned}" ]]; then
      remove_flag=1
    fi
  fi
  if [[ ${run_tree_pruning} -eq 1 ]]; then
    if [[ -s "${file_og_stat_tree}" ]]; then
      num_gene_pruned=$(gg_count_fasta_records "${file_og_trimmed_aln_pruned}")
      num_stat_branch_row=$(wc -l < "${file_og_stat_branch}")
      if [[ $((${num_gene_pruned} * 2)) -lt ${num_stat_branch_row} ]]; then
        remove_flag=1
      fi
    fi
  fi
  if [[ ${remove_flag} -eq 1 ]]; then
    echo "Downstream output files are inconsistent with the run_tree_pruning setting."
    for file in "${file_to_remove[@]}"; do
      echo "Not found: ${file}"
      if [[ -e "${file}" ]]; then
        echo "Deleting: ${file}"
        rm -f -- "${file}"
      fi
    done
    echo "Completed the deletion of inconsistent files."
  else
    echo "Downstream output files are consistent with the run_tree_pruning setting."
  fi
else
  gg_step_skip "${task}"
fi

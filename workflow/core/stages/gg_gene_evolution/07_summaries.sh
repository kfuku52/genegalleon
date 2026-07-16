# shellcheck shell=bash
# Sourced by gg_gene_evolution_core.sh.

task="summary statistics"
if is_output_older_than_inputs "^file_og_" "${file_og_tree_plot}"; then
  summary_flag=0
else
  summary_flag=$?
fi
task="Synteny neighborhood grouping"
if [[ ${treevis_synteny} -eq 1 ]] && { [[ ${run_summary} -eq 1 ]] || [[ ${run_tree_plot} -eq 1 ]]; }; then
  synteny_source_dir="${dir_sp_cds}"
  synteny_sequence_mode="cds"
  if [[ "${input_sequence_mode}" == "protein" ]] && species_protein_input_has_files; then
    synteny_source_dir="${dir_sp_protein_input}"
    synteny_sequence_mode="protein"
  fi
  synteny_needs_update=0
  if [[ ! -s "${file_og_synteny}" || "${file_og_primary_fasta}" -nt "${file_og_synteny}" ]]; then
    synteny_needs_update=1
  fi
  if [[ ${synteny_needs_update} -eq 1 ]]; then
    gg_step_start "${task}"
    if [[ ! -d "${synteny_source_dir}" ]]; then
      echo "Sequence directory not found. Skipping synteny panel input generation: ${synteny_source_dir}"
    elif [[ ! -d "${dir_sp_gff}" ]]; then
      echo "species_gff directory not found. Skipping synteny panel input generation: ${dir_sp_gff}"
    elif [[ ! -s "${file_og_primary_fasta}" ]]; then
      echo "Focal sequence fasta file not found. Skipping synteny panel input generation: ${file_og_primary_fasta}"
    else
      synteny_evalue="${query_blast_evalue}"
      if [[ "$(printf '%s' "${query_blast_evalue}" | tr '[:upper:]' '[:lower:]')" == "auto" ]]; then
        synteny_evalue_query_fasta="${og_id}.synteny.evalue.aa.tmp.fasta"
        if ! prepare_synteny_evalue_fasta "${synteny_evalue_query_fasta}"; then
          rm -f -- "${synteny_evalue_query_fasta}"
          echo "Failed to prepare amino-acid FASTA for synteny query_blast_evalue=auto: ${file_og_primary_fasta}"
          exit 1
        fi
        if ! resolve_query_blast_evalue "${query_blast_evalue}" "${synteny_evalue_query_fasta}" "${query_blast_auto_evalue_maxlen_cutoffs}"; then
          rm -f -- "${synteny_evalue_query_fasta}"
          exit 1
        fi
        synteny_evalue="${effective_query_blast_evalue}"
        echo "synteny auto E-value: query_count=${query_blast_query_num_seqs} min_aa_len=${query_blast_query_min_aa_len} avg_aa_len=${query_blast_query_avg_aa_len} max_aa_len=${query_blast_query_max_aa_len}"
        echo "synteny auto E-value: cutoffs=${query_blast_auto_evalue_maxlen_cutoffs} effective_synteny_evalue=${synteny_evalue}"
        rm -f -- "${synteny_evalue_query_fasta}"
      fi
      python "${gg_support_dir}/synteny_neighbors.py" \
        --focal_cds_fasta "${file_og_primary_fasta}" \
        --dir_sp_cds "${synteny_source_dir}" \
        --dir_sp_gff "${dir_sp_gff}" \
        --cache_dir "${gg_workspace_output_dir}/species_gff_info" \
        --lock_dir "${file_og_parameters_dir}/synteny_locks" \
        --gff2genestat_script "${gg_support_dir}/gff2genestat.py" \
        --input_sequence_mode "${synteny_sequence_mode}" \
        --window "${treevis_synteny_window}" \
        --evalue "${synteny_evalue}" \
        --genetic_code "${genetic_code}" \
        --threads "${GG_TASK_CPUS}" \
        --outfile "${file_og_synteny}"
      if [[ ! -s "${file_og_synteny}" ]]; then
        echo "No synteny links were generated: ${file_og_synteny}"
      fi
    fi
  else
    gg_step_skip "${task}"
  fi
else
  gg_step_skip "${task}"
fi
task="summary statistics"
disable_if_no_input_file "run_summary" "${file_og_rooted_tree_analysis}"
if [[ (${summary_flag} -eq 1 || ! -s "${file_og_stat_branch}" || ! -s "${file_og_stat_tree}") && ${run_summary} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ -s "${file_og_notung_reconcil}" ]]; then
    unzip -qf "${file_og_notung_reconcil}"
  fi
  notung_root_log_for_summary="PLACEHOLDER"
  if [[ -d "./${og_id}.notung.root" ]]; then
    notung_log_candidates=()
    mapfile -t notung_log_candidates < <(find "./${og_id}.notung.root" -maxdepth 1 -type f -name "*.ntglog" | sort)
    if [[ ${#notung_log_candidates[@]} -gt 0 ]]; then
      notung_root_log_for_summary="${notung_log_candidates[0]}"
    fi
  fi
  notung_reconcil_stats_for_summary="PLACEHOLDER"
  if [[ -d "./${og_id}.notung_reconcile" ]]; then
    reconcil_stats_candidates=()
    mapfile -t reconcil_stats_candidates < <(find "./${og_id}.notung_reconcile" -maxdepth 1 -type f -name "*.reconciled*.parsable.txt" | sort)
    if [[ ${#reconcil_stats_candidates[@]} -gt 0 ]]; then
      notung_reconcil_stats_for_summary="${reconcil_stats_candidates[0]}"
    fi
  fi
  if [[ ${run_tree_pruning} -eq 1 ]]; then
    generax2orthogroup_statistics="PLACEHOLDER" # generax nhx should be pruned to get used here.
  else
    generax2orthogroup_statistics=${file_og_generax_nhx}
  fi
  summary_unaligned_fasta="PLACEHOLDER"
  summary_trimmed_fasta="PLACEHOLDER"
  summary_promoter_fasta="PLACEHOLDER"
  summary_tmp_files=()
  if [[ -s "${file_og_primary_fasta}" ]]; then
    summary_unaligned_fasta="${og_id}.summary.unaligned.fasta"
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_primary_fasta}" --out-file "${summary_unaligned_fasta}"
    summary_tmp_files+=("${summary_unaligned_fasta}")
  fi
  if [[ -s "${file_og_trimmed_aln_analysis}" ]]; then
    summary_trimmed_fasta="${og_id}.summary.trimmed.fasta"
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file "${summary_trimmed_fasta}"
    summary_tmp_files+=("${summary_trimmed_fasta}")
  fi
  if [[ -s "${file_og_promoter_fasta}" ]]; then
    summary_promoter_fasta="${og_id}.summary.promoter.fasta"
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_promoter_fasta}" --out-file "${summary_promoter_fasta}"
    summary_tmp_files+=("${summary_promoter_fasta}")
  fi

  python "${gg_support_dir}/orthogroup_statistics.py" \
    --species_tree "${species_tree_pruned}" \
    --unaligned_aln "${summary_unaligned_fasta}" \
    --trimal_aln "${summary_trimmed_fasta}" \
    --unrooted_tree "${file_og_unrooted_tree_analysis}" \
    --rooted_tree "${file_og_rooted_tree_analysis}" \
    --rooting_log "${file_og_rooted_log}" \
    --notung_root_log "${notung_root_log_for_summary}" \
    --notung_reconcil_stats "${notung_reconcil_stats_for_summary}" \
    --dated_tree "${file_og_dated_tree_analysis}" \
    --dated_log "${file_og_dated_tree_log}" \
    --generax_nhx "${generax2orthogroup_statistics}" \
    --hyphy_dnds_json "${file_og_hyphy_dnds}" \
    --hyphy_relax_json "${file_og_hyphy_relax}" \
    --hyphy_relax_reversed_json "${file_og_hyphy_relax_reversed}" \
    --l1ou_tree "${file_og_l1ou_fit_tree}" \
    --l1ou_regime "${file_og_l1ou_fit_regime}" \
    --l1ou_leaf "${file_og_l1ou_fit_leaf}" \
    --expression "${file_og_expression}" \
    --mapdnds_tree_dn "${file_og_mapdnds_dn}" \
    --mapdnds_tree_ds "${file_og_mapdnds_ds}" \
    --codeml_tsv "${file_og_codeml_two_ratio}" \
    --character_gff "${file_og_gff_info}" \
    --fimo "${file_og_fimo}" \
    --promoter_fasta "${summary_promoter_fasta}" \
    --scm_intron "${file_og_scm_intron_summary}" \
    --csubst_b "${file_og_csubst_b}" \
    --gene_pgls_stats "${file_og_gene_pgls}" \
    --species_pgls_stats "${file_og_species_pgls}" \
    --rpsblast "${file_og_rpsblast}" \
    --uniprot "${file_og_uniprot_annotation}" \
    --cdskit_localize "${file_og_cdskit_localize}" \
    --synteny "${file_og_synteny}" \
    --ncpu "${GG_TASK_CPUS}" \
    --clade_ortholog_prefix "${treevis_clade_ortholog_prefix}"
  if [[ ${#summary_tmp_files[@]} -gt 0 ]]; then
    rm -f -- "${summary_tmp_files[@]}"
  fi

  #--csubst_cb_stats ${file_og_csubst_cb_stats} \ # Does not support --arity 3 or larger

  cp_out orthogroup.branch.tsv "${file_og_stat_branch}"
  cp_out orthogroup.tree.tsv "${file_og_stat_tree}"

else
  gg_step_skip "${task}"
fi

task="Query marker annotation"
if [[ "${mode_gene_evolution}" == "query2family" && ${treevis_query_marker} -eq 1 && -s "${file_og_stat_branch}" ]]; then
  query_marker_needs_update=0
  if ! awk -F $'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "query_marker") found = 1; exit(found ? 0 : 1) }' "${file_og_stat_branch}"; then
    query_marker_needs_update=1
  elif [[ -s "${file_query_gene}" && "${file_query_gene}" -nt "${file_og_stat_branch}" ]]; then
    query_marker_needs_update=1
  elif [[ -s "${file_og_query_aa_fasta}" && "${file_og_query_aa_fasta}" -nt "${file_og_stat_branch}" ]]; then
    query_marker_needs_update=1
  elif [[ -s "${file_og_query_blast}" && "${file_og_query_blast}" -nt "${file_og_stat_branch}" ]]; then
    query_marker_needs_update=1
  fi
  if [[ ${query_marker_needs_update} -eq 1 ]]; then
    gg_step_start "${task}"
    python "${gg_support_dir}/annotate_stat_branch_query_markers.py" \
      --stat_branch "${file_og_stat_branch}" \
      --query_gene "${file_query_gene}" \
      --query_aa_fasta "${file_og_query_aa_fasta}" \
      --query_blast "${file_og_query_blast}" \
      --min_query_blast_coverage "${query_blast_coverage}" \
      --outfile "orthogroup.branch.query_marker.tsv"
    cp_out "orthogroup.branch.query_marker.tsv" "${file_og_stat_branch}"
    summary_flag=1
  else
    gg_step_skip "${task}"
  fi
else
  gg_step_skip "${task}"
fi

if [[ ${treevis_synteny} -eq 1 && -s "${file_og_synteny}" ]]; then
  if [[ ! -s "${file_og_tree_plot}" || "${file_og_synteny}" -nt "${file_og_tree_plot}" ]]; then
    summary_flag=1
  fi
fi

task="stat_branch2tree_plot"
disable_if_no_input_file "run_tree_plot" "${file_og_stat_branch}" "${file_og_stat_tree}"
if [[ ${run_tree_plot} -eq 1 ]]; then
  if ! Rscript -e "if (!requireNamespace('ggimage', quietly=TRUE)) quit(status=1)" > /dev/null 2>&1; then
    echo "ggimage package is unavailable. Disabling run_tree_plot."
    run_tree_plot=0
  fi
fi
if ([[ ${summary_flag} -eq 1 || ! -s "${file_og_tree_plot}" ]]) && [[ ${run_tree_plot} -eq 1 ]]; then
  gg_step_start "${task}"

  num_tip_treeplot=$(
    awk -F $'\t' '
        NR==1 {
          col=0
          for (i=1; i<=NF; i++) {
            if ($i=="so_event") {
              col=i
              break
            }
          }
          next
        }
        (col>0 && $col=="L") {n++}
        END {print n+0}
      ' "${file_og_stat_branch}"
  )
  panel11_trimmed_aln="${file_og_trimmed_aln_analysis}"
  panel11_trimmed_n=$(gg_count_fasta_records "${panel11_trimmed_aln}")
  if [[ ${panel11_trimmed_n} -lt ${num_tip_treeplot} ]]; then
    for candidate in \
      "${file_og_clipkit}" \
      "${file_og_orthogroup_extraction_fasta}" \
      "${file_og_maxalign}" \
      "${file_og_mafft}" \
      "${file_og_primary_fasta}"; do
      candidate_n=$(gg_count_fasta_records "${candidate}")
      if [[ ${candidate_n} -ge ${num_tip_treeplot} ]]; then
        panel11_trimmed_aln="${candidate}"
        panel11_trimmed_n=${candidate_n}
        break
      fi
    done
  fi
  panel11_untrimmed_aln="${file_og_untrimmed_aln_analysis}"
  panel11_untrimmed_n=$(gg_count_fasta_records "${panel11_untrimmed_aln}")
  if [[ ${panel11_untrimmed_n} -lt ${num_tip_treeplot} ]]; then
    for candidate in \
      "${file_og_orthogroup_extraction_fasta}" \
      "${file_og_mafft}" \
      "${file_og_primary_fasta}"; do
      candidate_n=$(gg_count_fasta_records "${candidate}")
      if [[ ${candidate_n} -ge ${num_tip_treeplot} ]]; then
        panel11_untrimmed_aln="${candidate}"
        panel11_untrimmed_n=${candidate_n}
        break
      fi
    done
  fi
  echo "Tree plot alignment inputs: tips=${num_tip_treeplot}, trimmed=${panel11_trimmed_n} (${panel11_trimmed_aln}), untrimmed=${panel11_untrimmed_n} (${panel11_untrimmed_aln})"

  if [[ ${treevis_clade_ortholog} -eq 1 ]]; then
    ortholog_prefix=${treevis_clade_ortholog_prefix}
  else
    ortholog_prefix=""
  fi
  cb_path=${file_og_csubst_cb_2/cb_2/cb_ARITY}

  tree_plot_panel_args=(
    "--panel1=tree,${treevis_branch_length},${treevis_support_value},${treevis_branch_color},L"
    "--panel2=heatmap,${treevis_heatmap_transform},abs,_,expression_"
    "--panel3=pointplot,no,rel,_,expression_"
    "--panel4=cluster_membership,${treevis_max_intergenic_dist}"
    "--panel5=synteny,${file_og_synteny},${treevis_synteny_window}"
    "--panel6=tiplabel"
  )
  panel_index=7
  if [[ "${mode_gene_evolution}" == "query2family" && ${treevis_query_marker} -eq 1 ]]; then
    tree_plot_panel_args+=("--panel${panel_index}=categorical,query_marker,Query,-")
    panel_index=$((panel_index + 1))
  fi
  tree_plot_panel_args+=(
    "--panel${panel_index}=signal_peptide"
  )
  panel_index=$((panel_index + 1))
  tree_plot_panel_args+=(
    "--panel${panel_index}=transmembrane_domain"
  )
  panel_index=$((panel_index + 1))
  tree_plot_panel_args+=(
    "--panel${panel_index}=intron_number"
  )
  panel_index=$((panel_index + 1))
  tree_plot_panel_args+=(
    "--panel${panel_index}=domain,${file_og_rpsblast}"
  )
  panel_index=$((panel_index + 1))
  tree_plot_panel_args+=(
    "--panel${panel_index}=alignment,${panel11_trimmed_aln},${panel11_untrimmed_aln}"
  )
  panel_index=$((panel_index + 1))
  tree_plot_panel_args+=(
    "--panel${panel_index}=fimo,${promoter_bp},${fimo_qvalue}"
  )
  panel_index=$((panel_index + 1))
  tree_plot_panel_args+=(
    "--panel${panel_index}=meme,${file_og_meme}"
  )
  panel_index=$((panel_index + 1))
  tree_plot_panel_args+=(
    "--panel${panel_index}=ortholog,${ortholog_prefix},${file_og_dated_tree}"
  )

  TREEVIS_SPECIES_PARSER="${species_label_parser}" \
  Rscript "${gg_support_dir}/stat_branch2tree_plot.r" \
    --stat_branch="${file_og_stat_branch}" \
    --max_delta_intron_present="${treevis_retrotransposition_delta_intron}" \
    --width="7.2" \
    --rel_widths="" \
    "${tree_plot_panel_args[@]}" \
    --show_branch_id="yes" \
    --event_method="${treevis_event_method}" \
    --species_color_table="PLACEHOLDER" \
    --pie_chart_value_transformation="${treevis_pie_chart_value_transformation}" \
    --long_branch_display="${treevis_long_branch_display}" \
    --long_branch_ref_quantile="${treevis_long_branch_ref_quantile}" \
    --long_branch_detect_ratio="${treevis_long_branch_detect_ratio}" \
    --long_branch_cap_ratio="${treevis_long_branch_cap_ratio}" \
    --long_branch_tail_shrink="${treevis_long_branch_tail_shrink}" \
    --long_branch_max_fraction="${treevis_long_branch_max_fraction}" \
    --protein_convergence="100,100,yes,3-${csubst_max_arity},${cb_path},${csubst_cutoff_stat}"

  if [[ -e "df_fimo.tsv" ]]; then
    mv_out "df_fimo.tsv" "${file_og_fimo_collapsed}"
  fi
  mv_out stat_branch2tree_plot.pdf "${file_og_tree_plot}"
else
  gg_step_skip "${task}"
fi

# Copy parameter files and codes to ${file_og_parameters_dir} for record
mkdir -p "${file_og_parameters_dir}"
file_params=(
  "${file_sp_trait}"
  "${species_tree}"
  "${species_tree_pruned}"
)
for file_from in "${file_params[@]}"; do
  file_to="${file_og_parameters_dir}/$(basename "${file_from}")"
  if [[ ! -e "${file_from}" ]]; then
    continue
  fi
  lock_file="${file_to}.lock"
  (
    if ! gg_shared_lock_acquire "${lock_file}" "parameter artifact copy (${file_to})"; then
      exit 1
    fi
    gg_shared_lock_start_heartbeat "${lock_file}"
    heartbeat_pid=${GG_SHARED_LOCK_HEARTBEAT_PID:-}
    cleanup_parameter_copy_lock() {
      gg_shared_lock_stop_heartbeat "${heartbeat_pid}"
      gg_shared_lock_release "${lock_file}"
    }
    trap cleanup_parameter_copy_lock EXIT
    filesize_from=$(stat -c%s "${file_from}")
    filesize_to=0
    if [[ -s "${file_to}" ]]; then
      filesize_to=$(stat -c%s "${file_to}")
    fi
    if [[ ${filesize_from} -ne ${filesize_to} ]]; then
      echo "Storing important files for record: ${file_to}"
      cp_out "${file_from}" "${file_to}"
    fi
  ) || exit 1
done

cd "${gg_workspace_dir}" || exit 1
remove_empty_subdirs "${dir_output_active}"

if [[ -s "${file_og_stat_branch}" && -s "${file_og_stat_tree}" && -s "${file_og_tree_plot}" && ${gg_debug_mode:-0} -eq 0 ]]; then
  echo "Output files detected."
  echo "${file_og_stat_branch}"
  echo "${file_og_stat_tree}"
  echo "${file_og_tree_plot}"
  echo "$(date): Exiting Singularity environment"
  exit 8
elif [[ -s "${file_og_stat_branch}" && -s "${file_og_stat_tree}" && -s "${file_og_tree_plot}" && ${gg_debug_mode:-0} -eq 1 ]]; then
  echo "Output files detected & debug mode."
else
  echo "Output files not found."
fi

###################
echo "$(date): Exiting Singularity environment"

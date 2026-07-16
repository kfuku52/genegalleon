# shellcheck shell=bash
# Sourced by gg_gene_evolution_core.sh.

task="Parameter estimation for mapdNdS"
disable_if_no_input_file "run_mapdnds_parameter_estimation" "${file_og_rooted_tree_analysis}" "${file_og_trimmed_aln_analysis}"
if [[ ! -s "${file_og_mapdnds_parameter}" && ${run_mapdnds_parameter_estimation} -eq 1 ]]; then
  gg_step_start "${task}"

  nwkit drop --target intnode --support yes --name yes \
    --infile "${file_og_rooted_tree_analysis}" |
    nwkit sanitize --remove_singleton yes --resolve_polytomy no \
      > mapdnds_input.nwk

  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file ./mapdnds_input.fasta

  # F3X4+G4 shouldn not be changed otherwise iqtree2mapnh.py has to be updated.
  build_iqtree_mem_args
  iqtree \
    -s mapdnds_input.fasta \
    -m "GY+F3X4+G4" \
    -te mapdnds_input.nwk \
    -T AUTO \
    --threads-max "${GG_TASK_CPUS}" \
    --seqtype "CODON${genetic_code}" \
    --prefix "${og_id}.iqtree2mapdNdS" \
    "${IQTREE_MEM_ARGS[@]}" \
    --ancestral \
    --seed 12345 \
    --redo

  python "${gg_support_dir}/iqtree2mapnh.py" \
    --iqtree "${og_id}.iqtree2mapdNdS.iqtree" \
    --log "${og_id}.iqtree2mapdNdS.log" \
    --state "${og_id}.iqtree2mapdNdS.state" \
    --alignment mapdnds_input.fasta \
    --treefile "${og_id}.iqtree2mapdNdS.treefile" \
    --rooted_tree mapdnds_input.nwk \
    --genetic_code "${genetic_code}"

  if [[ -s "iqtree2mapnh.params" && -s "iqtree2mapnh.nwk" ]]; then
    echo "iqtree2mapnh was successfully completed."
    if [[ -e "${og_id}.mapdnds_parameter" ]]; then
      rm -rf -- "${og_id}.mapdnds_parameter"
    fi
    mkdir -p "${og_id}.mapdnds_parameter"
    mv_out "iqtree2mapnh.params" ./"${og_id}".mapdnds_parameter
    mv_out "iqtree2mapnh.nwk" ./"${og_id}".mapdnds_parameter
    zip -r "${og_id}.mapdnds_parameter.zip" "${og_id}.mapdnds_parameter"
    mv_out "${og_id}.mapdnds_parameter.zip" "${file_og_mapdnds_parameter}"
  else
    echo "iqtree2mapnh.params was not generated."
  fi
else
  gg_step_skip "${task}"
fi

task="mapdNdS main run"
disable_if_no_input_file "run_mapdnds" "${file_og_mapdnds_parameter}" "${file_og_trimmed_aln_analysis}"
if [[ (! -s "${file_og_mapdnds_dn}" || ! -s "${file_og_mapdnds_ds}") && ${run_mapdnds} -eq 1 ]]; then
  gg_step_start "${task}"

  unzip -o "${file_og_mapdnds_parameter}"
  cd "${dir_tmp}/${og_id}.mapdnds_parameter" || exit 1
  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file ./mapdnds_input.fasta
  normalize_mapnh_params_for_mapnh_v1 "iqtree2mapnh.params" "${genetic_code}"

  mapnh_exit_code=0
  if ! mapnh \
    SEQ=mapdnds_input.fasta \
    TREE=iqtree2mapnh.nwk \
    OUT=${og_id} \
    param=iqtree2mapnh.params \
    2>&1 | tee mapnh.log.txt; then
    mapnh_exit_code=${PIPESTATUS[0]:-1}
  fi

  if [[ -s "${og_id}.dN.dnd" && -s "${og_id}.dS.dnd" ]]; then
    echo "mapnh successfully generated dN and dS trees."
    mv_out "${og_id}.dN.dnd" "${file_og_mapdnds_dn}"
    mv_out "${og_id}.dS.dnd" "${file_og_mapdnds_ds}"
  else
    echo "mapnh did not generate dN/dS trees (exit code: ${mapnh_exit_code})."
    echo "mapnh output and HyPhy output are managed separately; no cross-substitution is applied."
  fi
  cd "${dir_tmp}" || exit 1
else
  gg_step_skip "${task}"
fi

task="CodeML two-ratio model"
disable_if_no_input_file "run_codeml_two_ratio" "${file_og_rooted_tree_analysis}" "${file_og_trimmed_aln_analysis}" "${file_sp_trait}"
if [[ ! -s "${file_og_codeml_two_ratio}" && "${run_codeml_two_ratio}" -eq 1 ]]; then
  gg_step_start "${task}"

  binarize_species_trait "${file_sp_trait}" species_trait_binary.tsv
  write_species_trait_foreground_regex_table species_trait_binary.tsv foreground.tsv
  IFS=$'\t' read -r -a colname_array < foreground.tsv

  if [[ ${#colname_array[@]} -le 1 ]]; then
    echo "No trait columns were detected in ${file_sp_trait}. Skipping ${task}."
  else
    for ((i = 1; i < ${#colname_array[@]}; i++)); do
      trait="${colname_array[$i]}"
      echo "Processing trait: ${trait}"
      awk -F'\t' -v trait_col="$((i + 1))" 'NR>1 && $trait_col == 1 { print $1 }' foreground.tsv > "foreground_${trait}.txt"

      target_spnode=$(paste -sd'|' "foreground_${trait}.txt") # Bar-separated list of target species nodes
      echo "Regular expression for CodeML foreground node search: ${target_spnode}"

      nwkit drop \
        --infile "${file_og_rooted_tree_analysis}" \
        --target "intnode" \
        --support "yes" \
        --name "yes" |
        nwkit mark \
          --pattern "${target_spnode}" \
          --insert_txt "#1" \
          --insert_sep "" \
          --target "mrca" \
          --target_only_clade "yes" \
          --outfile "codeml_input_${trait}.nwk"
      echo "CodeML input tree: $(< "codeml_input_${trait}.nwk")"

      bash "${gg_support_dir}/shorten_fasta_newick_names.sh" \
        "${file_og_trimmed_aln_analysis}" "codeml_input2_${trait}.fasta" "codeml_input_${trait}.nwk" "codeml_input2_${trait}.nwk" 90

      if grep -q "#1:0;$" "codeml_input_${trait}.nwk"; then
        exit_code1=0
      else
        exit_code1=1
      fi
      if grep -q "#1" "codeml_input_${trait}.nwk"; then
        exit_code2=0
      else
        exit_code2=1
      fi
      flag_unanalyzable=0
      if [[ ${exit_code1} -eq 0 || ${exit_code2} -eq 1 ]]; then
        flag_unanalyzable=1
      fi
      if [[ ${flag_unanalyzable} -eq 1 ]]; then
        echo "Target species tree node (${target_spnode}) is the gene tree root node. Generating an empty output file."
        codeml_out_treelength="NA"
        codeml_out_treelength_dn="NA"
        codeml_out_treelength_ds="NA"
        codeml_out_kappa="NA"
        codeml_out_background_omega="NA"
        codeml_out_foreground_omega="NA"
        codeml_out_time="NA"
      else
        python -c 'from pathlib import Path; import sys; template = Path(sys.argv[1]).read_text(encoding="utf-8"); rendered = template.replace("__SEQFILE__", sys.argv[2]).replace("__TREEFILE__", sys.argv[3]).replace("__ICODE__", sys.argv[4]); Path(sys.argv[5]).write_text(rendered, encoding="utf-8")' \
          "${gg_support_dir}/codeml/codeml_two_ratio.ctl.template" \
          "codeml_input2_${trait}.fasta" \
          "codeml_input2_${trait}.nwk" \
          "$((genetic_code - 1))" \
          "my_codeml_${trait}.ctl"

        codeml "my_codeml_${trait}.ctl"
        codeml_out_treelength=$(awk '/^tree length =/ {sub(/^tree length =[[:space:]]*/, "", $0); print; exit}' mlc)
        codeml_out_treelength_dn=$(awk '/^tree length for dN:/ {sub(/^tree length for dN:[[:space:]]*/, "", $0); print; exit}' mlc)
        codeml_out_treelength_ds=$(awk '/^tree length for dS:/ {sub(/^tree length for dS:[[:space:]]*/, "", $0); print; exit}' mlc)
        codeml_out_kappa=$(awk '/^kappa \(ts\/tv\) =/ {sub(/^kappa \(ts\/tv\) =[[:space:]]*/, "", $0); print; exit}' mlc)
        read -r -a codeml_out_omegas <<< "$(awk '/^w \(dN\/dS\) for branches:/ {sub(/^w \(dN\/dS\) for branches:[[:space:]]*/, "", $0); print; exit}' mlc)"
        codeml_out_background_omega=${codeml_out_omegas[0]:-}
        codeml_out_foreground_omega=${codeml_out_omegas[1]:-}
        codeml_out_time=$(awk '/^Time used:/ {sub(/^Time used:[[:space:]]*/, "", $0); print; exit}' mlc)
      fi
      if [[ -n "${codeml_out_background_omega}" && -n "${codeml_out_foreground_omega}" ]]; then
        echo "The task '${task}' has completed successfully for trait '${trait}'."
        printf 'tree_length_%s\ttree_length_dn_%s\ttree_length_ds_%s\tkappa_%s\tbackground_omega_%s\tforeground_omega_%s\tcodeml_time_%s\n' \
          "${trait}" "${trait}" "${trait}" "${trait}" "${trait}" "${trait}" "${trait}" > "file_og_codeml_two_ratio_${trait}.tsv"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
          "${codeml_out_treelength}" "${codeml_out_treelength_dn}" "${codeml_out_treelength_ds}" "${codeml_out_kappa}" "${codeml_out_background_omega}" "${codeml_out_foreground_omega}" "${codeml_out_time}" >> "file_og_codeml_two_ratio_${trait}.tsv"
      else
        echo "The task '${task}' failed for trait '${trait}'."
      fi
    done

    # Combine all codeml output files
    missing_files=()
    codeml_output_files=()
    for ((i = 1; i < ${#colname_array[@]}; i++)); do
      trait="${colname_array[$i]}"
      codeml_output_file="file_og_codeml_two_ratio_${trait}.tsv"
      if [[ -s "${codeml_output_file}" ]]; then
        codeml_output_files+=("${codeml_output_file}")
      else
        missing_files+=("${codeml_output_file}")
      fi
    done
    if [[ ${#missing_files[@]} -gt 0 ]]; then
      echo "The following codeml output files are missing:"
      for f in "${missing_files[@]}"; do
        echo "${f}"
      done
      echo "The task has failed: ${task}"
    else
      echo "All codeml two-ratio model output files are generated. Combining them into a single tsv file: ${file_og_codeml_two_ratio}"
      header_files=()
      data_files=()
      for f in "${codeml_output_files[@]}"; do
        base=$(basename "$f")
        head -n1 "$f" > "header_${base}"
        tail -n1 "$f" > "data_${base}"
        header_files+=("header_${base}")
        data_files+=("data_${base}")
      done
      paste -d$'\t' "${header_files[@]}" > combined_header.tsv
      paste -d$'\t' "${data_files[@]}" > combined_data.tsv
      cat combined_header.tsv combined_data.tsv > "${og_id}.codeml.two_ratio.tmp.tsv"
      mv_out "${og_id}.codeml.two_ratio.tmp.tsv" "${file_og_codeml_two_ratio}"
      echo "The task has completed successfully: ${task}"
    fi
  fi
else
  gg_step_skip "${task}"
fi

task="HyPhy dN-dS estimation"
disable_if_no_input_file "run_hyphy_dnds" "${file_og_rooted_tree_analysis}" "${file_og_trimmed_aln_analysis}"
if [[ ! -s "${file_og_hyphy_dnds}" && ${run_hyphy_dnds} -eq 1 ]]; then
  gg_step_start "${task}"

  nwkit drop --target intnode --support yes --name yes \
    --infile "${file_og_rooted_tree_analysis}" |
    nwkit sanitize --remove_singleton yes --resolve_polytomy no \
      > "hyphy_input.nwk"

  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file "hyphy_input.fasta"

  hyphy_genetic_code=$(get_hyphy_genetic_code "${genetic_code}")

  hyphy "${gg_support_dir}/hyphy/FitMG94.bf" \
    --alignment "hyphy_input.fasta" \
    --tree "hyphy_input.nwk" \
    --code "${hyphy_genetic_code}" \
    --frequencies "CF3x4" \
    --type "local" \
    --lrt "No" \
    --rooted "Yes" \
    --CPU "${GG_TASK_CPUS}"

  # --lrt "Yes" took too long time for some genes. 20 sec vs 10 min in a small tree.

  hyphy_dnds_json=""
  if [[ -s "hyphy_input.fasta.FITTER.json" ]]; then
    hyphy_dnds_json="hyphy_input.fasta.FITTER.json"
  else
    hyphy_dnds_candidates=()
    mapfile -t hyphy_dnds_candidates < <(find "." -maxdepth 1 -type f -name "hyphy_input.fasta*.json" | sort)
    if [[ ${#hyphy_dnds_candidates[@]} -gt 0 ]]; then
      hyphy_dnds_json="${hyphy_dnds_candidates[0]}"
    fi
  fi
  if [[ -z "${hyphy_dnds_json}" ]]; then
    echo "HyPhy FitMG94 output JSON was not detected. Exiting."
    exit 1
  fi
  mv_out "${hyphy_dnds_json}" "${file_og_hyphy_dnds}"
else
  gg_step_skip "${task}"
fi

# HyPhy RELAX


task="HyPhy RELAX"
disable_if_no_input_file "run_hyphy_relax" "${file_og_rooted_tree_analysis}" "${file_og_trimmed_aln_analysis}" "${file_sp_trait}"
if [[ ! -s "${file_og_hyphy_relax}" && ${run_hyphy_relax} -eq 1 ]]; then
  gg_step_start "${task}"
  run_hyphy_relax_for_all_traits 1 "${file_og_hyphy_relax}"
  if [[ -s "${file_og_hyphy_relax}" ]]; then
    echo "The task has completed successfully: ${task}"
  else
    echo "The task has failed: ${task}"
  fi
else
  gg_step_skip "${task}"
fi

task="HyPhy RELAX with reversed foreground/background"
disable_if_no_input_file "run_hyphy_relax_reversed" "${file_og_rooted_tree_analysis}" "${file_og_trimmed_aln_analysis}" "${file_sp_trait}"
if [[ ! -s "${file_og_hyphy_relax_reversed}" && ${run_hyphy_relax_reversed} -eq 1 ]]; then
  gg_step_start "${task}"
  run_hyphy_relax_for_all_traits 0 "${file_og_hyphy_relax_reversed}"
  if [[ -s "${file_og_hyphy_relax_reversed}" ]]; then
    echo "The task has completed successfully: ${task}"
  else
    echo "The task has failed: ${task}"
  fi
else
  gg_step_skip "${task}"
fi

task="Stochastic character mapping of intron evolution"
disable_if_no_input_file "run_scm_intron" "${file_og_gff_info}" "${file_og_dated_tree_analysis}"
if [[ ! -s "${file_og_scm_intron_summary}" && ${run_scm_intron} -eq 1 ]]; then
  gg_step_start "${task}"

  Rscript "${gg_support_dir}/scm_intron_evolution.r" \
    --tree_file="${file_og_dated_tree_analysis}" \
    --trait_file="${file_og_gff_info}" \
    --intron_gain_rate="${intron_gain_rate}" \
    --retrotransposition_rate="${retrotransposition_rate}" \
    --nrep=1000 \
    --nslots="${GG_TASK_CPUS}"

  cp_out intron_evolution_summary.tsv "${file_og_scm_intron_summary}"
  if [[ -e intron_evolution_plot.pdf ]]; then
    cp_out intron_evolution_plot.pdf "${file_og_scm_intron_plot}"
  fi
else
  gg_step_skip "${task}"
fi

task="l1ou"
disable_if_no_input_file "run_l1ou" "${file_og_trimmed_aln_analysis}" "${file_og_expression}" "${file_og_dated_tree_analysis}"
if [[ (! -s "${file_og_l1ou_fit_rdata}" || ! -s "${file_og_l1ou_fit_tree}" || ! -s "${file_og_l1ou_fit_regime}" || ! -s "${file_og_l1ou_fit_leaf}") && ${run_l1ou} -eq 1 ]]; then
  gg_step_start "${task}"

  num_gene=$(gg_count_fasta_records "${file_og_trimmed_aln_analysis}")
  if [[ ${num_gene} -ge ${large_tree_num_gene} ]]; then
    max_nshift=${large_tree_max_nshift}
  else
    max_nshift=0
  fi

  fit_ind_file=''
  if [[ ${l1ou_use_fit_file} -eq 1 && -s "${file_og_l1ou_fit_rdata}" ]]; then
    fit_ind_file=${file_og_l1ou_fit_rdata}
  fi

  l1ou_cmd=(
    Rscript "${gg_support_dir}/detect_OU_shift_kfl1ou.r"
    --max_nshift="${max_nshift}"
    --tree_file="${file_og_dated_tree_analysis}"
    --trait_file="${file_og_expression}"
    --nslots="${GG_TASK_CPUS}"
    --criterion="${l1ou_criterion}"
    --nbootstrap="${l1ou_nbootstrap}"
    --fit_ind_file="${fit_ind_file}"
    --fit_conv_file=''
    --alpha_upper="${l1ou_alpha_upper}"
    --detect_convergence="${l1ou_convergence}"
    --replicate_sep="_"
  )
  "${l1ou_cmd[@]}"

  mv_out fit_ind.RData "${file_og_l1ou_fit_rdata}"
  mv_out l1ou_tree.tsv "${file_og_l1ou_fit_tree}"
  mv_out l1ou_regime.tsv "${file_og_l1ou_fit_regime}"
  mv_out l1ou_leaf.tsv "${file_og_l1ou_fit_leaf}"
  mv_out l1ou_plot.pdf "${file_og_l1ou_fit_plot}"
  if [[ ${l1ou_nbootstrap} -gt 0 ]]; then
    cp_out l1ou_bootstrap.tsv "${dir_l1ou_bootstrap}"/"${l1ou_bootstrap}"
  fi
  if [[ ${l1ou_convergence} -eq 1 ]]; then
    cp_out fit_conv.RData "${file_og_l1ou_fit_conv_rdata}"
  fi

else
  gg_step_skip "${task}"
fi

# shellcheck shell=bash
# Sourced by gg_genome_evolution_core.sh.

task="Materialize selected undated species tree summary"
if [[ ! -s "${file_undated_species_tree}" && ${species_tree_requested_for_orthofinder} -eq 1 ]]; then
  selected_undated_species_tree_source=""
  selected_undated_species_tree_is_astral=0
  case "${undated_species_tree}" in
    iqtree_pep)
      selected_undated_species_tree_source="${file_concat_iqtree_pep_root}"
      ;;
    iqtree_dna)
      selected_undated_species_tree_source="${file_concat_iqtree_dna_root}"
      ;;
    astral_pep)
      selected_undated_species_tree_source="${file_astral_tree_pep}"
      selected_undated_species_tree_is_astral=1
      ;;
    astral_dna)
      selected_undated_species_tree_source="${file_astral_tree_dna}"
      selected_undated_species_tree_is_astral=1
      ;;
  esac

  if [[ -s "${selected_undated_species_tree_source}" ]]; then
    gg_step_start "${task}"
    echo "Rebuilding missing undated species tree summary from cached source:"
    echo "from: ${selected_undated_species_tree_source}"
    echo "to: ${file_undated_species_tree}"
    if [[ ${selected_undated_species_tree_is_astral} -eq 1 ]]; then
      rm -f -- "tmp.undated_species_tree.nwk"
      if ! nwkit drop --name yes --target intnode --infile "${selected_undated_species_tree_source}" |
        nwkit label --target intnode --force yes --outfile "tmp.undated_species_tree.nwk"; then
        echo "Warning: Failed to convert cached ASTRAL tree with nwkit drop|label. Copying cached tree instead."
      fi
      if [[ -s "tmp.undated_species_tree.nwk" ]]; then
        mv_out "tmp.undated_species_tree.nwk" "${file_undated_species_tree}"
      else
        cp_out "${selected_undated_species_tree_source}" "${file_undated_species_tree}"
      fi
    else
      cp_out "${selected_undated_species_tree_source}" "${file_undated_species_tree}"
    fi
    if [[ -s "${file_undated_species_tree}" ]]; then
      Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="${file_undated_species_tree}"
    fi
  else
    echo "Selected undated species tree source is not available: ${selected_undated_species_tree_source}"
  fi
else
  gg_step_skip "${task}"
fi

task="Species tree plotting"
if [[ ${run_plot_species_trees} -eq 1 ]]; then
  if [[ "${input_sequence_mode}" == "protein" ]]; then
    disable_if_no_input_file "run_plot_species_trees" "${file_concat_iqtree_pep_root}" "${file_astral_tree_pep}"
  else
    disable_if_no_input_file "run_plot_species_trees" "${file_concat_iqtree_dna_root}" "${file_concat_iqtree_pep_root}" "${file_astral_tree_dna}" "${file_astral_tree_pep}"
  fi
fi
if [[ ! -s "${file_plot_species_trees}" && ${run_plot_species_trees} -eq 1 ]]; then
  gg_step_start "${task}"

  Rscript "${gg_support_dir}/plot_species_trees.r" \
    --iqtree_dna_nwk="${file_concat_iqtree_dna_root}" \
    --iqtree_pep_nwk="${file_concat_iqtree_pep_root}" \
    --iqtree_dna_log="${dir_concat_iqtree_dna}/tmp.concat.cds.input.fasta.log" \
    --iqtree_pep_log="${dir_concat_iqtree_pep}/tmp.concat.pep.input.fasta.log" \
    --astral_dna_nwk="${file_astral_tree_dna}" \
    --astral_pep_nwk="${file_astral_tree_pep}" \
    --astral_dna_log="${file_astral_log_dna}" \
    --astral_pep_log="${file_astral_log_pep}"

  if [[ -s "species_trees.pdf" ]]; then
    echo "Output file found for the task: ${task}"
    mv_out "species_trees.pdf" "${file_plot_species_trees}"
  fi
else
  gg_step_skip "${task}"
fi

task="BUSCO stacked bar species tree plotting"
disable_if_no_input_file "run_plot_species_trees" "${file_undated_species_tree}" "${dir_species_busco_full}"
if [[ (! -s "${file_species_tree_busco_cds_pdf}" || ! -s "${file_species_tree_busco_cds_svg}" || ! -s "${file_species_tree_busco_summary}") && ${run_plot_species_trees} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_species_tree_summary}"
  cd "${dir_species_tree_summary}" || exit 1

  Rscript "${gg_support_dir}/annotation_summary.r" \
    --dir_species_tree="${dir_species_tree}" \
    --dir_species_cds_busco="${dir_species_busco_full}" \
    --dir_species_genome_busco="" \
    --dir_species_annotation="" \
    --dir_species_cds_fx2tab="" \
    --dir_species_genome_fx2tab="" \
    --file_species_trait="" \
    --file_orthogroup_gene_count="" \
    --min_og_species='auto'

  if [[ -e "Rplots.pdf" ]]; then
    rm -f -- "Rplots.pdf"
  fi
  cd "${dir_tmp}" || exit 1

  if [[ -s "${file_species_tree_busco_cds_pdf}" && -s "${file_species_tree_busco_cds_svg}" ]]; then
    echo "Output file found for the task: ${task}"
  fi
else
  gg_step_skip "${task}"
fi

task="Time-constrained tree preparation"
disable_if_no_input_file "run_constrained_tree" "${file_undated_species_tree}"
if [[ ! -s "${file_constrained_tree}" && ${run_constrained_tree} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_parent_dir "${file_constrained_tree}"
  ensure_dir "${dir_nwkit_download_dir}"
  if [[ ${timetree_constraint} -eq 1 ]]; then
    nwkit_args=(
      --download_dir "${dir_nwkit_download_dir}"
      --infile "${file_undated_species_tree}"
      --timetree "ci"
      --min_clade_prop 0.2
      --outfile "tmp.constrained.tree.nwk"
      --species-parser "${species_label_parser}"
    )
    if [[ -n "${species_label_regex}" ]]; then
      nwkit_args+=(--species-regex "${species_label_regex}")
    fi
    if [[ -n "${species_label_map_tsv}" ]]; then
      nwkit_args+=(--species-map-tsv "${species_label_map_tsv}")
    fi
    nwkit mcmctree "${nwkit_args[@]}"
    if [[ -s "tmp.constrained.tree.nwk" ]]; then
      mv_out "tmp.constrained.tree.nwk" "${file_constrained_tree}"
    fi
  else
    if [[ ${#mcmctree_divergence_time_constraints[@]} -eq 0 ]]; then
      echo "timetree_constraint=0 requires mcmctree_divergence_time_constraints_str with one or more records."
      echo "Example: Arabidopsis_thaliana,Oryza_sativa,130,-|Arabidopsis_thaliana,Amborella_trichopoda,150,200"
      exit 1
    fi
    tree_string=$(< "${file_undated_species_tree}")
    for mdtc in "${mcmctree_divergence_time_constraints[@]}"; do
      echo "applying ${mdtc}"
      mcmctree_params=()
      if ! parse_mcmctree_constraint_record "${mdtc}" mcmctree_params; then
        exit 1
      fi
      nwkit_args=(
        --download_dir "${dir_nwkit_download_dir}"
        --left_species "${mcmctree_params[0]}"
        --right_species "${mcmctree_params[1]}"
        --species-parser "${species_label_parser}"
      )
      if [[ -n "${species_label_regex}" ]]; then
        nwkit_args+=(--species-regex "${species_label_regex}")
      fi
      if [[ -n "${species_label_map_tsv}" ]]; then
        nwkit_args+=(--species-map-tsv "${species_label_map_tsv}")
      fi
      if [[ "${mcmctree_params[2]}" != "-" ]]; then
        nwkit_args+=(--lower_bound "${mcmctree_params[2]}")
      fi
      if [[ "${mcmctree_params[3]}" != "-" ]]; then
        nwkit_args+=(--upper_bound "${mcmctree_params[3]}")
      fi
      echo "nwkit mcmctree params: ${nwkit_args[*]}"
      tree_string=$(printf '%s\n' "${tree_string}" | nwkit mcmctree "${nwkit_args[@]}")
    done
    printf '%s\n' "${tree_string}" > "tmp.constrained.tree.nwk"
    if [[ -s "tmp.constrained.tree.nwk" ]]; then
      mv_out "tmp.constrained.tree.nwk" "${file_constrained_tree}"
    fi
  fi
  if [[ -s "${file_constrained_tree}" ]]; then
    if ! normalize_iq2mc_constraint_tree "${file_constrained_tree}"; then
      echo "Error: Failed to normalize the constrained tree for IQ2MC."
      exit 1
    fi
  fi
else
  gg_step_skip "${task}"
fi

task="Constrained range plotting"
disable_if_no_input_file "run_plot_constrained_tree" "${file_constrained_tree}"
if [[ ! -s "${file_plot_constrained_tree}" && ${run_plot_constrained_tree} -eq 1 ]]; then
  gg_step_start "${task}"
  Rscript "${gg_support_dir}/plot_constrained_tree.r" \
    --infile="${file_constrained_tree}" \
    --outfile="tmp.constrained_tree_plot.pdf"
  if [[ -s "tmp.constrained_tree_plot.pdf" ]]; then
    mv_out "tmp.constrained_tree_plot.pdf" "${file_plot_constrained_tree}"
  fi
  if [[ -s "${file_plot_constrained_tree}" ]]; then
    echo "Output file found for the task: ${task}"
  fi
else
  gg_step_skip "${task}"
fi

task="IQ2MC step 2 (IQ-TREE Hessian/control generation)"
disable_if_no_input_file "run_mcmctree1" "${file_constrained_tree}" "${file_concat_cds}"
if [[ (! -s "${file_iq2mc_ctl}" || ! -s "${file_iq2mc_hessian}" || ! -s "${file_iq2mc_rooted_tree}" || ! -s "${file_iq2mc_dummy_phy}") && ${run_mcmctree1} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "$(dirname "${file_iq2mc_ctl}")"

  if ! normalize_iq2mc_constraint_tree "${file_constrained_tree}"; then
    echo "Error: Failed to normalize the constrained tree for IQ2MC."
    exit 1
  fi
  if ! clear_directory_contents_safe "$(dirname "${file_iq2mc_ctl}")"; then
    echo "Error: Failed to clear IQ2MC working directory safely."
    exit 1
  fi
  mcmctree_time_scale_factor=$(resolve_mcmctree_time_scale_factor)
  if [[ "${mcmctree_time_scale_factor}" != "1" ]]; then
    echo "MCMCTree calibration ages will be divided by ${mcmctree_time_scale_factor} for internal IQ2MC/MCMCTree execution."
    echo "Public GeneGalleon tree outputs remain in the original time unit."
  fi

  iq2mc_work_dir=$(mktemp -d "${dir_tmp}/tmp.iq2mc.work.XXXXXX")
  iq2mc_scaled_constraint_tree="${iq2mc_work_dir}/iq2mc.scaled.constraint.nwk"
  scale_mcmctree_calibrations_file "${file_constrained_tree}" "${iq2mc_scaled_constraint_tree}" "${mcmctree_time_scale_factor}" "down"

  cd "${iq2mc_work_dir}" || exit 1
  seqkit seq --threads 1 "${file_concat_cds}" --out-file "./tmp.iq2mc.concat.cds.fasta"

  if ! "${iq2mc_binary}" \
    -s "./tmp.iq2mc.concat.cds.fasta" \
    -m "${nucleotide_model}" \
    -te "${iq2mc_scaled_constraint_tree}" \
    --dating mcmctree \
    --mcmc-bds "${mcmc_birth_death_sampling}" \
    --mcmc-clock "${mcmc_clock_model}" \
    --mcmc-iter "${mcmc_burnin},${mcmc_sampfreq},${mcmc_nsample}" \
    -T "${GG_TASK_CPUS}" \
    --prefix iq2mc; then
    echo "Error: IQ2MC step 2 failed. Deleting generated files."
    rm -f -- "${file_iq2mc_prefix}".*
  elif [[ ! -s "${iq2mc_work_dir}/iq2mc.mcmctree.ctl" || ! -s "${iq2mc_work_dir}/iq2mc.mcmctree.hessian" || ! -s "${iq2mc_work_dir}/iq2mc.rooted.nwk" || ! -s "${iq2mc_work_dir}/iq2mc.dummy.phy" ]]; then
    echo "Error: IQ2MC step 2 did not generate all expected files."
  else
    scale_mcmctree_ctl_rootage_file "${iq2mc_work_dir}/iq2mc.mcmctree.ctl" "${file_iq2mc_ctl}" "${mcmctree_time_scale_factor}" "up"
    cp_out "${iq2mc_work_dir}/iq2mc.mcmctree.hessian" "${file_iq2mc_hessian}"
    cp_out "${iq2mc_work_dir}/iq2mc.dummy.phy" "${file_iq2mc_dummy_phy}"
    scale_mcmctree_calibrations_file "${iq2mc_work_dir}/iq2mc.rooted.nwk" "${file_iq2mc_rooted_tree}" "${mcmctree_time_scale_factor}" "up"
  fi
  if [[ ! -s "${file_iq2mc_ctl}" || ! -s "${file_iq2mc_hessian}" || ! -s "${file_iq2mc_rooted_tree}" || ! -s "${file_iq2mc_dummy_phy}" ]]; then
    rm -f -- "${file_iq2mc_prefix}".*
  fi
  rm -f -- "./tmp.iq2mc.concat.cds.fasta"
  cd "${dir_tmp}" || exit 1
  if [[ "${GG_KEEP_MCMCTREE_RAW_DEBUG:-0}" == "1" ]]; then
    echo "Keeping internal scaled IQ2MC debug directory until tmp cleanup: ${iq2mc_work_dir}"
    echo "Set delete_tmp_dir=0 to preserve it after the workflow."
  else
    rm -rf -- "${iq2mc_work_dir}"
  fi
else
  gg_step_skip "${task}"
fi

task="IQ2MC step 3 (MCMCtree dating run)"
disable_if_no_input_file "run_mcmctree2" "${file_iq2mc_ctl}" "${file_iq2mc_hessian}" "${file_iq2mc_rooted_tree}" "${file_iq2mc_dummy_phy}"
if [[ ! -s "${file_mcmctree_figtree_tre}" && ${run_mcmctree2} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_mcmctree2}"

  if ! clear_directory_contents_safe "${dir_mcmctree2}"; then
    echo "Error: Failed to clear MCMCtree working directory safely."
    exit 1
  fi
  mcmctree_time_scale_factor=$(resolve_mcmctree_time_scale_factor)
  cd "${dir_mcmctree2}" || exit 1
  cp_out "${file_iq2mc_ctl}" ./
  cp_out "${file_iq2mc_hessian}" ./
  cp_out "${file_iq2mc_rooted_tree}" ./
  cp_out "${file_iq2mc_dummy_phy}" ./

  mcmctree_work_dir=$(mktemp -d "${dir_tmp}/tmp.mcmctree.work.XXXXXX")
  scale_mcmctree_ctl_rootage_file "${file_iq2mc_ctl}" "${mcmctree_work_dir}/$(basename "${file_iq2mc_ctl}")" "${mcmctree_time_scale_factor}" "down"
  cp_out "${file_iq2mc_hessian}" "${mcmctree_work_dir}/"
  cp_out "${file_iq2mc_dummy_phy}" "${mcmctree_work_dir}/"
  scale_mcmctree_calibrations_file "${file_iq2mc_rooted_tree}" "${mcmctree_work_dir}/$(basename "${file_iq2mc_rooted_tree}")" "${mcmctree_time_scale_factor}" "down"
  cd "${mcmctree_work_dir}" || exit 1

  # Ensure MCMCtree emits CI-rich summaries (including 95% HPD annotations in FigTree output).
  ctl_basename="$(basename "${file_iq2mc_ctl}")"
  sed -i -e "s/FossilErrprint/FossilErr\\nprint/g" "${ctl_basename}"
  if ! grep -q "^print[[:space:]]*=" "${ctl_basename}"; then
    echo "print = 1" >> "${ctl_basename}"
  fi
  normalize_mcmctree_ctl_for_installed_paml "${ctl_basename}"

  if ! mcmctree "${ctl_basename}"; then
    echo "Error: IQ2MC step 3 failed."
    rm -f -- "${file_mcmctree_raw_output}"
  elif [[ ! -s "${mcmctree_work_dir}/$(basename "${file_mcmctree_raw_output}")" ]]; then
    echo "Error: IQ2MC step 3 did not generate $(basename "${file_mcmctree_raw_output}")."
    rm -f -- "${file_mcmctree_raw_output}" "${file_mcmctree_figtree_tre}"
  elif extract_scaled_mcmctree_figtree "${mcmctree_work_dir}/$(basename "${file_mcmctree_raw_output}")" "${file_mcmctree_figtree_tre}" "${mcmctree_time_scale_factor}"; then
    {
      echo "GeneGalleon ran MCMCTree in an internal scaled time unit."
      echo "Raw scaled MCMCTree output is not retained by default."
      echo "The FigTree block below is converted back to the original public time unit."
      cat "${file_mcmctree_figtree_tre}"
    } > "${file_mcmctree_raw_output}"
  else
    echo "Error: Failed to extract an original-unit FigTree tree block from MCMCTree output."
    rm -f -- "${file_mcmctree_raw_output}" "${file_mcmctree_figtree_tre}"
  fi
  cd "${dir_tmp}" || exit 1
  if [[ "${GG_KEEP_MCMCTREE_RAW_DEBUG:-0}" == "1" ]]; then
    echo "Keeping internal scaled MCMCTree debug directory until tmp cleanup: ${mcmctree_work_dir}"
    echo "Set delete_tmp_dir=0 to preserve it after the workflow."
  else
    rm -rf -- "${mcmctree_work_dir}"
  fi
else
  gg_step_skip "${task}"
fi

if [[ ! -s "${file_mcmctree_figtree_tre}" && -s "${file_mcmctree_raw_output}" ]]; then
  echo "Generating ${file_mcmctree_figtree_tre} from ${file_mcmctree_raw_output}"
  awk '
  /Species tree for FigTree/ {print; in_figtree=1; next}
  in_figtree && /^\(\(/ {print; count++; if (count >= 3) exit}
  ' "${file_mcmctree_raw_output}" > "tmp.mcmctree2.txt"
  if [[ -s "tmp.mcmctree2.txt" ]]; then
    mv_out "tmp.mcmctree2.txt" "${file_mcmctree_figtree_tre}"
  fi
  if [[ ! -s "${file_mcmctree_figtree_tre}" ]]; then
    echo "Warning: Failed to extract FigTree content from ${file_mcmctree_raw_output}. Copying raw file instead."
    cp_out "${file_mcmctree_raw_output}" "${file_mcmctree_figtree_tre}"
  fi
fi

task="Convert tree format"
disable_if_no_input_file "run_convert_tree_format" "${file_mcmctree_figtree_tre}"
if [[ ! -s "${file_mcmctree_dated_nwk}" && ${run_convert_tree_format} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_parent_dir "${file_mcmctree_dated_nwk}"

  if grep -q -e "UTREE" "${file_mcmctree_figtree_tre}"; then
    grep -e "UTREE" "${file_mcmctree_figtree_tre}" |
      sed -e "s/.*UTREE 1 = //" -e "s/;.*/;/" \
        > "${dir_mcmctree2}/mcmctree_95CI.nwk"

    grep -e "UTREE" "${file_mcmctree_figtree_tre}" |
      sed -e "s/.*UTREE 1 = //" -e "s/;.*/;/" -e "s/[[:space:]]*\[&95%={[0-9.]*,[[:space:]][0-9.]*}\][[:space:]]*//g" -e "s/:[[:space:]]/:/g" \
        > "${dir_mcmctree2}/mcmctree_no95CI.nwk"
  else
    tree_line="$(awk '/^\(\(/ {line=$0} END {print line}' "${file_mcmctree_figtree_tre}")"
    if [[ -n "${tree_line}" ]]; then
      echo "${tree_line}" > "${dir_mcmctree2}/mcmctree_95CI.nwk"
      echo "${tree_line}" |
        sed -e "s/[[:space:]]*\[&95%={[0-9.]*,[[:space:]][0-9.]*}\][[:space:]]*//g" -e "s/:[[:space:]]/:/g" \
          > "${dir_mcmctree2}/mcmctree_no95CI.nwk"
    else
      echo "Error: Failed to detect a tree string in ${file_mcmctree_figtree_tre}"
      rm -f -- "${dir_mcmctree2}/mcmctree_95CI.nwk" "${dir_mcmctree2}/mcmctree_no95CI.nwk"
    fi
  fi

  if [[ -s "${dir_mcmctree2}/mcmctree_no95CI.nwk" ]]; then
    Rscript -e "library(ape); t=read.tree(\"${dir_mcmctree2}/mcmctree_no95CI.nwk\"); \
    t[['node.label']]=paste0('s',1:(length(t[['tip.label']])-1)); \
    write.tree(t, \"${file_mcmctree_dated_nwk}\")"
  else
    echo "Error: Missing mcmctree_no95CI.nwk. Skipping tree conversion."
  fi

  if [[ ! -s "${file_dated_species_tree}" ]]; then
    echo "Dated species tree for gg_pipeline is not placed yet: ${file_dated_species_tree}"
    if [[ -s "${file_mcmctree_dated_nwk}" ]]; then
      echo "Copying from: ${file_mcmctree_dated_nwk}"
      echo "Copying to: ${file_dated_species_tree}"
      cp_out "${file_mcmctree_dated_nwk}" "${file_dated_species_tree}"
      if [[ -s "${file_plot_mcmctree_pdf}" ]]; then
        echo "Copying from: ${file_plot_mcmctree_pdf}"
        echo "Copying to: ${file_dated_species_tree_pdf}"
        cp_out "${file_plot_mcmctree_pdf}" "${file_dated_species_tree_pdf}"
      fi
      echo "Please manually check whether the species tree is valid."
    fi
  else
    echo "Dated species tree for gg_pipeline is already placed: ${file_dated_species_tree}"
    echo "If necessary, please replace the file with: ${file_mcmctree_dated_nwk}"
  fi
else
  gg_step_skip "${task}"
fi

task="Dated species tree plotting"
disable_if_no_input_file "run_plot_mcmctreer" "${file_mcmctree_dated_nwk}"
if [[ ! -s "${file_plot_mcmctree_pdf}" && ${run_plot_mcmctreer} -eq 1 ]]; then
  gg_step_start "${task}"

  Rscript "${gg_support_dir}/plot_mcmctreer.r" \
    --infile="${file_mcmctree_dated_nwk}" \
    --outfile="tmp.plot_mcmctreer.pdf"
  if [[ -s "tmp.plot_mcmctreer.pdf" ]]; then
    mv_out "tmp.plot_mcmctreer.pdf" "${file_plot_mcmctree_pdf}"
  fi

  if [[ -s "${file_plot_mcmctree_pdf}" ]]; then
    echo "Output file found for the task: ${task}"
    echo "Output file: ${file_plot_mcmctree_pdf}"
    if [[ -s "${file_dated_species_tree}" ]]; then
      echo "Copying from: ${file_plot_mcmctree_pdf}"
      echo "Copying to: ${file_dated_species_tree_pdf}"
      cp_out "${file_plot_mcmctree_pdf}" "${file_dated_species_tree_pdf}"
    fi
  fi
else
  gg_step_skip "${task}"
fi

remove_empty_subdirs "${dir_species_tree}"
if [[ ${delete_tmp_dir} -eq 1 ]]; then
  if [[ -d "${dir_tmp}" ]]; then
    echo "Removing tmp directory: ${dir_tmp}"
    rm -rf -- "${dir_tmp}"
  fi
fi

# Orthogroup inference





ensure_dir "${dir_tmp}"
cd "${dir_tmp}" || exit 1

# shellcheck shell=bash
# Sourced by gg_genome_evolution_core.sh.

task="Concatenating single-copy CDS fasta files"
concat_alignment_ready=0
if [[ "${input_sequence_mode}" == "protein" ]]; then
  if [[ -s "${file_concat_pep}" ]]; then
    concat_alignment_ready=1
  fi
else
  if [[ -s "${file_concat_cds}" && -s "${file_concat_pep}" ]]; then
    concat_alignment_ready=1
  fi
fi
if [[ ${concat_alignment_ready} -eq 0 && ${run_concat_alignment} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_parent_dir "${file_concat_pep}"
  ensure_dir "${dir_concat_fasta}"
  if [[ "${input_sequence_mode}" != "protein" ]]; then
    ensure_parent_dir "${file_concat_cds}"
  fi
  mapfile -t trimal_files < <(find "${dir_single_copy_trimal}" -maxdepth 1 -type f -name "${single_copy_trimal_glob}" | sort)
  if [[ ${#trimal_files[@]} -eq 0 ]]; then
    echo "No trimmed single-copy FASTA files were found in: ${dir_single_copy_trimal}"
    exit 1
  fi

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    if [[ ${strictly_single_copy_only} -eq 0 ]]; then
      concat_pep_tmp="tmp.concat.pep.fa.gz"
      seqkit concat --full --fill "-" --threads "${GG_TASK_CPUS}" "${trimal_files[@]}" |
        seqkit sort --threads "${GG_TASK_CPUS}" |
        seqkit seq --threads "${GG_TASK_CPUS}" --out-file "${concat_pep_tmp}"
    else
      concat_pep_tmp="tmp.concat.pep.fa.gz"
      seqkit concat --threads "${GG_TASK_CPUS}" "${trimal_files[@]}" |
        seqkit sort --threads "${GG_TASK_CPUS}" |
        seqkit seq --threads "${GG_TASK_CPUS}" --out-file "${concat_pep_tmp}"
    fi
    mv_out "${concat_pep_tmp}" "${file_concat_pep}"
  else
    if [[ ${strictly_single_copy_only} -eq 0 ]]; then
      concat_cds_tmp="tmp.concat.cds.fa.gz"
      seqkit concat --full --fill "-" --threads "${GG_TASK_CPUS}" "${trimal_files[@]}" |
        seqkit sort --threads "${GG_TASK_CPUS}" |
        seqkit seq --threads "${GG_TASK_CPUS}" --out-file "${concat_cds_tmp}"
    else
      concat_cds_tmp="tmp.concat.cds.fa.gz"
      seqkit concat --threads "${GG_TASK_CPUS}" "${trimal_files[@]}" |
        seqkit sort --threads "${GG_TASK_CPUS}" |
        seqkit seq --threads "${GG_TASK_CPUS}" --out-file "${concat_cds_tmp}"
    fi
    mv_out "${concat_cds_tmp}" "${file_concat_cds}"
    seqkit translate --transl-table "${genetic_code}" --threads "${GG_TASK_CPUS}" "${file_concat_cds}" |
      seqkit seq --threads "${GG_TASK_CPUS}" --out-file "tmp.concat.pep.fa.gz"
    mv_out "tmp.concat.pep.fa.gz" "${file_concat_pep}"
  fi
else
  gg_step_skip "${task}"
fi

task="IQ-TREE of the concatenated alignment with a Protein evolution model"
if [[ ! -s "${file_concat_iqtree_pep}" && ${run_concat_iqtree_protein} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_concat_iqtree_pep}"

  ntaxa=$(gg_count_fasta_records "${file_concat_pep}")
  bootstrap_params=(--ufboot 1000 --bnni)
  if [[ ${ntaxa} -lt 4 ]]; then
    bootstrap_params=()
  fi
  iqtree_mem_args=()
  if [[ -n "${GG_MEM_TOOL_GB:-}" ]]; then
    iqtree_mem_args=(-mem "${GG_MEM_TOOL_GB}G")
  fi

  cd "${dir_concat_iqtree_pep}" || exit 1
  concat_pep_local="tmp.concat.pep.input.fasta"
  seqkit seq --threads 1 "${file_concat_pep}" --out-file "./${concat_pep_local}"
  if iqtree \
    -s "${concat_pep_local}" \
    -st AA \
    -m "${protein_model}" \
    -pre "${concat_pep_local}" \
    -nt "${GG_TASK_CPUS}" \
    "${iqtree_mem_args[@]}" \
    -seed 12345 \
    "${bootstrap_params[@]}"; then
    iqtree_exit_status=0
  else
    iqtree_exit_status=$?
  fi

  if [[ ${iqtree_exit_status} -ne 0 ]]; then
    echo "Error. IQ-TREE exited with non-zero status: ${iqtree_exit_status}"
  else
    echo "IQ-TREE successfully exited with zero status: ${iqtree_exit_status}"
    if [[ ${ntaxa} -lt 4 ]]; then
      out_nwk="${concat_pep_local}.treefile"
    else
      out_nwk="${concat_pep_local}.contree"
    fi
    cp_out "${out_nwk}" "${file_concat_iqtree_pep}"
    rm -f -- "./${concat_pep_local}"
  fi
  cd "${dir_tmp}" || exit 1
else
  gg_step_skip "${task}"
fi

task="Rooting of IQ-TREE's protein tree"
if [[ ! -s "${file_concat_iqtree_pep_root}" && ${run_concat_iqtree_protein} -eq 1 ]]; then
  gg_step_start "${task}"

  root_species_tree \
    "${file_concat_iqtree_pep}" \
    "${file_concat_iqtree_pep_root}" \
    "concatenated protein tree"
  restore_rooted_tree_internal_support \
    "${file_concat_iqtree_pep_root}" \
    "${file_concat_iqtree_pep}" \
    "concatenated protein tree"

  if [[ -s "${file_concat_iqtree_pep_root}" ]]; then
    Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="${file_concat_iqtree_pep_root}"
  fi

  if [[ -s "${file_concat_iqtree_pep_root}" && ${undated_species_tree} == 'iqtree_pep' ]]; then
    echo "Undated species tree is copied,"
    echo "from: ${file_concat_iqtree_pep_root}"
    echo "to: ${file_undated_species_tree}"
    cp_out "${file_concat_iqtree_pep_root}" "${file_undated_species_tree}"
    Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="${file_undated_species_tree}"
  fi
else
  gg_step_skip "${task}"
fi

task="IQ-TREE of the concatenated alignment with a DNA evolution model"
if [[ ! -s "${file_concat_iqtree_dna}" && ${run_concat_iqtree_dna} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_concat_iqtree_dna}"

  ntaxa=$(gg_count_fasta_records "${file_concat_cds}")
  bootstrap_params=(--ufboot 1000 --bnni)
  if [[ ${ntaxa} -lt 4 ]]; then
    bootstrap_params=()
  fi
  iqtree_mem_args=()
  if [[ -n "${GG_MEM_TOOL_GB:-}" ]]; then
    iqtree_mem_args=(-mem "${GG_MEM_TOOL_GB}G")
  fi

  cd "${dir_concat_iqtree_dna}" || exit 1
  concat_cds_local="tmp.concat.cds.input.fasta"
  seqkit seq --threads 1 "${file_concat_cds}" --out-file "./${concat_cds_local}"
  if iqtree \
    -s "${concat_cds_local}" \
    -st DNA \
    -m "${nucleotide_model}" \
    -pre "${concat_cds_local}" \
    -nt "${GG_TASK_CPUS}" \
    "${iqtree_mem_args[@]}" \
    -seed 12345 \
    "${bootstrap_params[@]}"; then
    iqtree_exit_status=0
  else
    iqtree_exit_status=$?
  fi

  if [[ ${iqtree_exit_status} -ne 0 ]]; then
    echo "Error. IQ-TREE exited with non-zero status: ${iqtree_exit_status}"
  else
    echo "IQ-TREE successfully exited with zero status: ${iqtree_exit_status}"
    if [[ ${ntaxa} -lt 4 ]]; then
      out_nwk="${concat_cds_local}.treefile"
    else
      out_nwk="${concat_cds_local}.contree"
    fi
    cp_out "${out_nwk}" "${file_concat_iqtree_dna}"
    rm -f -- "./${concat_cds_local}"
  fi
  cd "${dir_tmp}" || exit 1
else
  gg_step_skip "${task}"
fi

task="Rooting of IQ-TREE's DNA tree"
if [[ ! -s "${file_concat_iqtree_dna_root}" && ${run_concat_iqtree_dna} -eq 1 ]]; then
  gg_step_start "${task}"

  root_species_tree \
    "${file_concat_iqtree_dna}" \
    "${file_concat_iqtree_dna_root}" \
    "concatenated DNA tree"
  restore_rooted_tree_internal_support \
    "${file_concat_iqtree_dna_root}" \
    "${file_concat_iqtree_dna}" \
    "concatenated DNA tree"

  if [[ -s "${file_concat_iqtree_dna_root}" ]]; then
    Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="${file_concat_iqtree_dna_root}"
  fi

  if [[ -s "${file_concat_iqtree_dna_root}" && ${undated_species_tree} == 'iqtree_dna' ]]; then
    echo "Undated species tree is copied,"
    echo "from: ${file_concat_iqtree_dna_root}"
    echo "to: ${file_undated_species_tree}"
    cp_out "${file_concat_iqtree_dna_root}" "${file_undated_species_tree}"
    Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="${file_undated_species_tree}"
  fi
else
  gg_step_skip "${task}"
fi

task="IQ-TREE for individual single-copy protein trees"
ensure_dir "${dir_single_copy_iqtree_pep}"
iqtree_pep_tree_files=()
mapfile -t iqtree_pep_tree_files < <(gg_find_file_basenames "${dir_single_copy_iqtree_pep}" "*.nwk")
num_iqtree_pep=${#iqtree_pep_tree_files[@]}
if [[ ${num_busco_ids} -ne ${num_iqtree_pep} && ${run_individual_iqtree_pep} -eq 1 ]]; then
  gg_step_start "${task}"

  run_iqtree_pep() {
    local infile=$1
    local infile_base=${infile%%.*}
    local outfile="${dir_single_copy_iqtree_pep}/${infile_base}.pep.nwk"
    if [[ -s "${outfile}" ]]; then
      return 0
    fi
    local num_seq
    num_seq=$(gg_count_fasta_records "${dir_single_copy_trimal}/${infile}")
    if [[ ${num_seq} -lt 3 ]]; then
      echo "Skipped. At least 3 sequences are necessary for IQ-TREE: ${infile}"
      return 0
    fi
    if [[ "${input_sequence_mode}" == "protein" ]]; then
      seqkit seq --threads 1 "${dir_single_copy_trimal}/${infile}" > "tmp.${infile_base}.pep.fasta"
    else
      seqkit translate --transl-table "${genetic_code}" --threads 1 "${dir_single_copy_trimal}/${infile}" \
        > "tmp.${infile_base}.pep.fasta"
    fi
    iqtree \
      -s "tmp.${infile_base}.pep.fasta" \
      -m "${protein_model}" \
      -T 1 \
      "${iqtree_parallel_mem_args[@]}" \
      --prefix "tmp.${infile_base}" \
      --seed 12345 \
      --redo
    mv_out tmp."${infile_base}".treefile "${outfile}"
    rm -f -- "tmp.${infile_base}."*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_single_copy_trimal}" "${single_copy_trimal_glob}")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    run_iqtree_pep "${input_alignment_file}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

task="ASTRAL of individual single-copy protein trees"
if [[ (! -s "${file_astral_tree_pep}" || ! -s "${file_astral_log_pep}") && ${run_astral_pep} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_astral_pep}"

  if compgen -G "tmp.astral.*" > /dev/null; then
    rm -f -- tmp.astral.*
  fi

  if ! build_astral_input "${dir_single_copy_iqtree_pep}" "tmp.astral.merged.iqtree.nwk" "${astral_min_tips}"; then
    echo "Skipped. No eligible protein gene trees for ASTRAL after filtering."
  else
    astral-hybrid \
      --input "tmp.astral.merged.iqtree.nwk" \
      --output "tmp.astral.out.tree" \
      --mode 3 \
      --support 2 \
      --thread "${GG_TASK_CPUS}" \
      2> "tmp.astral.log.txt"

    labels=("pp1" "pp2" "pp3" "f1" "f2" "f3" "q1" "q2" "q3") # https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md
    for i in "${!labels[@]}"; do
      tmp_label_tree="tmp.astral.pep.${labels[i]}.tmp.nwk"
      python "${gg_support_dir}/extract_astral_support_label.py" \
        --infile "tmp.astral.out.tree" \
        --outfile "${tmp_label_tree}" \
        --label_key "${labels[i]}"
      save_astral_label_tree \
        "${tmp_label_tree}" \
        "single_copy.astral.pep.${labels[i]}.nwk" \
        "ASTRAL protein tree (${labels[i]})"
      if [[ -s "single_copy.astral.pep.${labels[i]}.nwk" ]] && ! astral_label_tree_rooting_is_deferred; then
        Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="single_copy.astral.pep.${labels[i]}.nwk"
      fi
      rm -f -- "${tmp_label_tree}"
    done

    if [[ -s "single_copy.astral.pep.q1.nwk" ]]; then
      shopt -s nullglob
      astral_pep_outputs=(single_copy.astral.pep.*)
      shopt -u nullglob
      if [[ ${#astral_pep_outputs[@]} -eq 0 ]]; then
        echo "ASTRAL protein outputs were expected but not found."
        exit 1
      fi
      mv_out "${astral_pep_outputs[@]}" "${dir_astral_pep}"
      mv_out "tmp.astral.log.txt" "${file_astral_log_pep}"
      echo "For more information on support values (e.g., f1, f2, pp1, q1, ...), please refer to: https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md" > "${dir_astral_pep}/README.txt"
    fi

    if [[ -s "${file_astral_tree_pep_q1}" ]]; then
      if optimize_astral_tree_branch_lengths "${file_astral_tree_pep_q1}" "${file_concat_pep}" "${protein_model}" "${file_astral_tree_pep}" "pep"; then
        transfer_deferred_astral_label_tree_roots "${file_astral_tree_pep}" "${dir_astral_pep}" "single_copy.astral.pep" "ASTRAL protein support trees"
        Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="${file_astral_tree_pep}"
      else
        echo "Warning: Falling back to unoptimized ASTRAL protein tree."
        cp_out "${file_astral_tree_pep_q1}" "${file_astral_tree_pep}"
      fi
    fi
  fi
  rm -f -- tmp.astral.*

  if [[ -s "${file_astral_tree_pep}" && ${undated_species_tree} == 'astral_pep' ]]; then
    echo "Undated species tree is copied,"
    echo "from: ${file_astral_tree_pep}"
    echo "to: ${file_undated_species_tree}"
    nwkit drop --name yes --target intnode --infile "${file_astral_tree_pep}" | nwkit label --target intnode --force yes --outfile "tmp.undated_species_tree.nwk"
    if [[ -s "tmp.undated_species_tree.nwk" ]]; then
      mv_out "tmp.undated_species_tree.nwk" "${file_undated_species_tree}"
    fi
    if [[ ! -s "${file_undated_species_tree}" ]]; then
      echo "Warning: Failed to convert optimized ASTRAL protein tree with nwkit drop|label. Copying optimized tree instead."
      cp_out "${file_astral_tree_pep}" "${file_undated_species_tree}"
    fi
    if [[ -s "${file_undated_species_tree}" ]]; then
      Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="${file_undated_species_tree}"
    fi
  fi
else
  gg_step_skip "${task}"
fi

task="IQ-TREE for individual single-copy DNA trees"
ensure_dir "${dir_single_copy_iqtree_dna}"
iqtree_dna_tree_files=()
mapfile -t iqtree_dna_tree_files < <(gg_find_file_basenames "${dir_single_copy_iqtree_dna}" "*.nwk")
num_iqtree_dna=${#iqtree_dna_tree_files[@]}
if [[ ${num_busco_ids} -ne ${num_iqtree_dna} && ${run_individual_iqtree_dna} -eq 1 ]]; then
  gg_step_start "${task}"

  run_iqtree_dna() {
    local infile=$1
    local infile_base=${infile%%.*}
    local outfile="${dir_single_copy_iqtree_dna}/${infile_base}.dna.nwk"
    if [[ -s "${outfile}" ]]; then
      return 0
    fi
    local num_seq
    num_seq=$(gg_count_fasta_records "${dir_single_copy_trimal}/${infile}")
    if [[ ${num_seq} -lt 3 ]]; then
      echo "Skipped. At least 3 sequences are necessary for IQ-TREE: ${infile}"
      return 0
    fi
    cp_out "${dir_single_copy_trimal}/${infile}" "./tmp.${infile_base}.input.fasta"
    iqtree \
      -s "./tmp.${infile_base}.input.fasta" \
      -m "${nucleotide_model}" \
      -T 1 \
      "${iqtree_parallel_mem_args[@]}" \
      --prefix "tmp.${infile_base}" \
      --seed 12345 \
      --redo
    mv_out tmp."${infile_base}".treefile "${outfile}"
    rm -f -- "tmp.${infile_base}."*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_single_copy_trimal}" "*.trimal.fa.gz")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    run_iqtree_dna "${input_alignment_file}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

task="ASTRAL of individual single-copy DNA trees"
if [[ (! -s "${file_astral_tree_dna}" || ! -s "${file_astral_log_dna}") && ${run_astral_dna} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_astral_dna}"

  if compgen -G "tmp.astral.*" > /dev/null; then
    rm -f -- tmp.astral.*
  fi

  if ! build_astral_input "${dir_single_copy_iqtree_dna}" "tmp.astral.merged.iqtree.nwk" "${astral_min_tips}"; then
    echo "Skipped. No eligible DNA gene trees for ASTRAL after filtering."
  else
    astral-hybrid \
      --input "tmp.astral.merged.iqtree.nwk" \
      --output "tmp.astral.out.tree" \
      --mode 3 \
      --support 2 \
      --thread "${GG_TASK_CPUS}" \
      2> "tmp.astral.log.txt"

    labels=("pp1" "pp2" "pp3" "f1" "f2" "f3" "q1" "q2" "q3") # https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md
    for i in "${!labels[@]}"; do
      tmp_label_tree="tmp.astral.dna.${labels[i]}.tmp.nwk"
      python "${gg_support_dir}/extract_astral_support_label.py" \
        --infile "tmp.astral.out.tree" \
        --outfile "${tmp_label_tree}" \
        --label_key "${labels[i]}"
      save_astral_label_tree \
        "${tmp_label_tree}" \
        "single_copy.astral.dna.${labels[i]}.nwk" \
        "ASTRAL DNA tree (${labels[i]})"
      if [[ -s "single_copy.astral.dna.${labels[i]}.nwk" ]] && ! astral_label_tree_rooting_is_deferred; then
        Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="single_copy.astral.dna.${labels[i]}.nwk"
      fi
      rm -f -- "${tmp_label_tree}"
    done

    if [[ -s "single_copy.astral.dna.q1.nwk" ]]; then
      shopt -s nullglob
      astral_dna_outputs=(single_copy.astral.dna.*)
      shopt -u nullglob
      if [[ ${#astral_dna_outputs[@]} -eq 0 ]]; then
        echo "ASTRAL DNA outputs were expected but not found."
        exit 1
      fi
      mv_out "${astral_dna_outputs[@]}" "${dir_astral_dna}"
      mv_out "tmp.astral.log.txt" "${file_astral_log_dna}"
      echo "For more information on support values (e.g., f1, f2, pp1, q1, ...), please refer to: https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md" > "${dir_astral_dna}/README.txt"
    fi

    if [[ -s "${file_astral_tree_dna_q1}" ]]; then
      if optimize_astral_tree_branch_lengths "${file_astral_tree_dna_q1}" "${file_concat_cds}" "${nucleotide_model}" "${file_astral_tree_dna}" "dna"; then
        transfer_deferred_astral_label_tree_roots "${file_astral_tree_dna}" "${dir_astral_dna}" "single_copy.astral.dna" "ASTRAL DNA support trees"
        Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="${file_astral_tree_dna}"
      else
        echo "Warning: Falling back to unoptimized ASTRAL DNA tree."
        cp_out "${file_astral_tree_dna_q1}" "${file_astral_tree_dna}"
      fi
    fi
  fi
  rm -f -- tmp.astral.*

  if [[ -s "${file_astral_tree_dna}" && ${undated_species_tree} == 'astral_dna' ]]; then
    echo "Undated species tree is copied,"
    echo "from: ${file_astral_tree_dna}"
    echo "to: ${file_undated_species_tree}"
    nwkit drop --name yes --target intnode --infile "${file_astral_tree_dna}" | nwkit label --target intnode --force yes --outfile "tmp.undated_species_tree.nwk"
    if [[ -s "tmp.undated_species_tree.nwk" ]]; then
      mv_out "tmp.undated_species_tree.nwk" "${file_undated_species_tree}"
    fi
    if [[ ! -s "${file_undated_species_tree}" ]]; then
      echo "Warning: Failed to convert optimized ASTRAL DNA tree with nwkit drop|label. Copying optimized tree instead."
      cp_out "${file_astral_tree_dna}" "${file_undated_species_tree}"
    fi
    if [[ -s "${file_undated_species_tree}" ]]; then
      Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="${file_undated_species_tree}"
    fi
  fi
else
  gg_step_skip "${task}"
fi

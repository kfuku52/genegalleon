# shellcheck shell=bash
# Sourced by gg_genome_evolution_core.sh.

task="Generating FASTA files for duplicate-aware BUSCO genes"
sync_genome_busco_summary_table_from_shared || true
disable_if_no_input_file "run_busco_dupaware_extract_fasta" "${file_genome_busco_summary_table}"
if [[ ${run_busco_dupaware_extract_fasta} -eq 1 ]]; then
  prepare_species_tree_input_dir
  gg_step_start "${task}"
  ensure_dir "${dir_busco_fasta}"
  busco_rows=()
  mapfile -t busco_rows < <(tail -n +2 "${file_genome_busco_summary_table}")
  num_busco_ids=${#busco_rows[@]}

  generate_genome_dupaware_busco_fasta() {
    local busco_idx=$1
    local busco_row="${busco_rows[${busco_idx}]}"
    local busco_id
    local outfile2
    local cols=()
    local genes=()
    local genes2=()
    IFS=$'\t' read -r -a cols <<< "${busco_row}"
    busco_id="${cols[0]:-}"
    if [[ -z "${busco_id}" ]]; then
      return 0
    fi
    if [[ ${#cols[@]} -gt 3 ]]; then
      genes=("${cols[@]:3}")
    fi
    outfile2="${dir_busco_fasta}/${busco_id}${genome_busco_fasta_suffix}"
    if [[ -s "${outfile2}" ]]; then
      return 0
    fi
    echo "busco_id: ${busco_id}"
    if [[ ! -s "${outfile2}" ]]; then
      mapfile -t genes2 < <(gg_busco_gene_tokens "split_duplicated" "${genes[@]}")
      if [[ ${#genes2[@]} -eq 0 ]]; then
        echo "Skipping. ${busco_id} has no genes in the selected species."
        return 0
      fi
      if [[ "${input_sequence_mode}" == "protein" ]]; then
        gg_seqkit_grep_by_patterns_from_infile_list 1 "species_tree_input_fasta_list.txt" "${genes2[@]}" |
          seqkit replace --pattern " .*" --replacement "" --ignore-case --threads 1 |
          seqkit seq --threads 1 --out-file "${outfile2}"
      else
        gg_seqkit_grep_by_patterns_from_infile_list 1 "species_tree_input_fasta_list.txt" "${genes2[@]}" |
          gg_prepare_cds_fasta_stream 1 |
          seqkit seq --threads 1 --out-file "${outfile2}"
      fi
      if [[ ! -s "${outfile2}" ]]; then
        echo "File is empty. Removing: ${outfile2}"
        rm -f -- "${outfile2}"
      fi
    fi
  }

  gg_find_fasta_files "${species_tree_input_dir}" 1 > species_tree_input_fasta_list.txt
  for ((busco_idx = 0; busco_idx < num_busco_ids; busco_idx++)); do
    wait_until_jobn_le ${GG_TASK_CPUS}
    generate_genome_dupaware_busco_fasta "${busco_idx}" &
  done
  wait_for_background_jobs
  rm -f -- species_tree_input_fasta_list.txt
else
  gg_step_skip "${task}"
fi

task="In-frame mafft alignment of duplicate-containing BUSCO genes"
if [[ ${run_busco_dupaware_mafft} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_busco_mafft}"

  run_mafft() {
    infile=$1
    infile_base=${infile%%.*}
    outfile=${dir_busco_mafft}/${infile_base}${genome_busco_aln_suffix}
    if [[ -s "${outfile}" ]]; then
      return 0
    fi
    echo "$(date): start mafft: ${infile_base}"

    if [[ "${input_sequence_mode}" == "protein" ]]; then
      seqkit seq --threads 1 "${dir_busco_fasta}/${infile}" --out-file "tmp.${infile_base}.input.pep.fasta"
      num_seq=$(gg_count_fasta_records "tmp.${infile_base}.input.pep.fasta")
      if [[ ${num_seq} -lt 2 ]]; then
        echo "Skipped MAFFT because fewer than 2 sequences were found: ${infile}"
        seqkit seq --threads 1 "tmp.${infile_base}.input.pep.fasta" --out-file "tmp.${infile_base}.pep.out.fa.gz"
        mv_out "tmp.${infile_base}.pep.out.fa.gz" "${outfile}"
        rm -f -- "tmp.${infile_base}"*
        return 0
      fi
      mafft \
        --auto \
        --amino \
        --thread 1 \
        "tmp.${infile_base}.input.pep.fasta" \
        > "tmp.${infile_base}.pep.aln.fasta"
      if [[ -s "tmp.${infile_base}.pep.aln.fasta" ]]; then
        seqkit seq --threads 1 "tmp.${infile_base}.pep.aln.fasta" --out-file "tmp.${infile_base}.pep.aln.out.fa.gz"
        mv_out "tmp.${infile_base}.pep.aln.out.fa.gz" "${outfile}"
      fi
    else
      seqkit seq --threads 1 "${dir_busco_fasta}/${infile}" --out-file "tmp.${infile_base}.input.cds.fasta"
      cdskit mask \
        --seqfile "tmp.${infile_base}.input.cds.fasta" \
        --outfile "tmp.${infile_base}.cds.fasta"

      num_seq=$(gg_count_fasta_records "tmp.${infile_base}.cds.fasta")
      if [[ ${num_seq} -lt 2 ]]; then
        echo "Skipped MAFFT/backalign because fewer than 2 sequences were found: ${infile}"
        seqkit seq --threads 1 "tmp.${infile_base}.cds.fasta" --out-file "tmp.${infile_base}.cds.out.fa.gz"
        mv_out "tmp.${infile_base}.cds.out.fa.gz" "${outfile}"
        rm -f -- "tmp.${infile_base}"*
        return 0
      fi

      seqkit translate \
        --allow-unknown-codon \
        --transl-table "${genetic_code}" \
        --threads 1 \
        "tmp.${infile_base}.cds.fasta" \
        > "tmp.${infile_base}.pep.fasta"

      mafft \
        --auto \
        --amino \
        --thread 1 \
        "tmp.${infile_base}.pep.fasta" \
        > "tmp.${infile_base}.pep.aln.fasta"

      cdskit backalign \
        --seqfile "tmp.${infile_base}.cds.fasta" \
        --aa_aln "tmp.${infile_base}.pep.aln.fasta" \
        --codontable "${genetic_code}" \
        --outfile "tmp.${infile_base}.cds.aln.fasta"

      if [[ -s "tmp.${infile_base}.cds.aln.fasta" ]]; then
        seqkit seq --threads 1 "tmp.${infile_base}.cds.aln.fasta" --out-file "tmp.${infile_base}.cds.aln.out.fa.gz"
        mv_out "tmp.${infile_base}.cds.aln.out.fa.gz" "${outfile}"
      fi
    fi
    rm -f -- "tmp.${infile_base}"*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_busco_fasta}" "${genome_busco_fasta_glob}")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    run_mafft "${input_alignment_file}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

task="TrimAl of duplicate-containing BUSCO genes"
if [[ ${run_busco_dupaware_trimal} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_busco_trimal}"

  run_trimal() {
    infile=$1
    infile_base=${infile%%.*}
    outfile="${dir_busco_trimal}/${infile_base}${genome_busco_trimal_suffix}"
    if [[ -s "${outfile}" ]]; then
      return 0
    fi
    if [[ "${input_sequence_mode}" == "protein" ]]; then
      seqkit seq --threads 1 "${dir_busco_mafft}/${infile}" --out-file "tmp.${infile_base}.pep.aln.fasta"
      trimal \
        -in "tmp.${infile_base}.pep.aln.fasta" \
        -out "tmp.${infile_base}.trimal.fasta" \
        -automated1
    else
      seqkit seq --remove-gaps --threads 1 "${dir_busco_mafft}/${infile}" > "tmp.${infile_base}.degap.fasta"
      seqkit translate --transl-table "${genetic_code}" --threads 1 "${dir_busco_mafft}/${infile}" > "tmp.${infile_base}.pep.fasta"

      trimal \
        -in "tmp.${infile_base}.pep.fasta" \
        -backtrans "tmp.${infile_base}.degap.fasta" \
        -out "tmp.${infile_base}.trimal.fasta" \
        -ignorestopcodon \
        -automated1
    fi

    if [[ -s "tmp.${infile_base}.trimal.fasta" ]]; then
      seqkit seq --threads 1 "tmp.${infile_base}.trimal.fasta" --out-file "tmp.${infile_base}.trimal.out.fa.gz"
      mv_out "tmp.${infile_base}.trimal.out.fa.gz" "${outfile}"
    fi
    rm -f -- "tmp.${infile_base}."*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_busco_mafft}" "${genome_busco_aln_glob}")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    run_trimal "${input_alignment_file}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

task="IQ-TREE for duplicate-containing BUSCO DNA trees"
if [[ ${run_busco_dupaware_iqtree_dna} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_busco_iqtree_dna}"

  busco_iqtree_dna() {
    infile=$1
    indir=$2
    outdir=$3
    infile_base=${infile%%.*}
    outfile="${outdir}/${infile_base}.busco.nwk"
    if [[ -s "${outfile}" ]]; then
      return 0
    fi
    num_seq=$(gg_count_fasta_records "${dir_busco_trimal}/${infile}")
    if [[ ${num_seq} -lt 3 ]]; then
      echo "Skipped. At least 3 sequences are necessary for IQ-TREE: ${infile}"
      return 0
    fi

    seqkit seq --threads 1 "${indir}/${infile}" --out-file "./tmp.${infile_base}.input.fasta"

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
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_busco_trimal}" "*.busco.trimal.fa.gz")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    busco_iqtree_dna "${input_alignment_file}" "${dir_busco_trimal}" "${dir_busco_iqtree_dna}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

task="IQ-TREE for duplicate-containing BUSCO protein trees"
if [[ ${run_busco_dupaware_iqtree_pep} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_busco_iqtree_pep}"

  busco_iqtree_pep() {
    infile=$1
    indir=$2
    outdir=$3
    infile_base=${infile%%.*}
    outfile="${outdir}/${infile_base}.busco.nwk"
    if [[ -s "${outfile}" ]]; then
      return 0
    fi
    num_seq=$(gg_count_fasta_records "${dir_busco_trimal}/${infile}")
    if [[ ${num_seq} -lt 3 ]]; then
      echo "Skipped. At least 3 sequences are necessary for IQ-TREE: ${infile}"
      return 0
    fi

    if [[ "${input_sequence_mode}" == "protein" ]]; then
      seqkit seq --threads 1 "${indir}/${infile}" > "./tmp.${infile_base}.input.fasta"
    else
      seqkit translate --transl-table "${genetic_code}" --threads 1 "${indir}/${infile}" > "./tmp.${infile_base}.input.fasta"
    fi

    iqtree \
      -s "./tmp.${infile_base}.input.fasta" \
      -m "${protein_model}" \
      -st AA \
      -T 1 \
      "${iqtree_parallel_mem_args[@]}" \
      --prefix "tmp.${infile_base}" \
      --seed 12345 \
      --redo

    mv_out tmp."${infile_base}".treefile "${outfile}"
    rm -f -- "tmp.${infile_base}."*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_busco_trimal}" "${genome_busco_trimal_glob}")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    busco_iqtree_pep "${input_alignment_file}" "${dir_busco_trimal}" "${dir_busco_iqtree_pep}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi


task="NOTUNG rooting of duplicate-containing BUSCO DNA trees"
disable_if_no_input_file "run_busco_dupaware_notung_root_dna" "${file_dated_species_tree}"
if [[ ${run_busco_dupaware_notung_root_dna} -eq 1 ]]; then
  gg_step_start "${task}"

  infiles=()
  mapfile -t infiles < <(gg_find_file_basenames "${dir_busco_iqtree_dna}")
  for infile in "${infiles[@]}"; do
    wait_until_jobn_le $((${GG_TASK_CPUS} / 2))
    busco_notung "${infile}" "${dir_busco_iqtree_dna}" "${dir_busco_notung_dna}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

task="NOTUNG rooting of duplicate-containing BUSCO protein trees"
disable_if_no_input_file "run_busco_dupaware_notung_root_pep" "${file_dated_species_tree}"
if [[ ${run_busco_dupaware_notung_root_pep} -eq 1 ]]; then
  gg_step_start "${task}"

  infiles=()
  mapfile -t infiles < <(gg_find_file_basenames "${dir_busco_iqtree_pep}")
  for infile in "${infiles[@]}"; do
    wait_until_jobn_le $((${GG_TASK_CPUS} / 2))
    busco_notung "${infile}" "${dir_busco_iqtree_pep}" "${dir_busco_notung_pep}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi


task="Species-tree-guided gene tree rooting of duplicate-containing BUSCO DNA trees"
disable_if_no_input_file "run_busco_dupaware_root_dna" "${file_dated_species_tree}"
if [[ ${run_busco_dupaware_root_dna} -eq 1 ]]; then
  gg_step_start "${task}"

  infiles=()
  mapfile -t infiles < <(gg_find_file_basenames "${dir_busco_notung_dna}")
  for infile in "${infiles[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    busco_species_tree_assisted_gene_tree_rooting "${infile}" "${dir_busco_notung_dna}" "${dir_busco_iqtree_dna}" "${dir_busco_rooted_txt_dna}" "${dir_busco_rooted_nwk_dna}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

task="Species-tree-guided gene tree rooting of duplicate-containing BUSCO protein trees"
disable_if_no_input_file "run_busco_dupaware_root_pep" "${file_dated_species_tree}"
if [[ ${run_busco_dupaware_root_pep} -eq 1 ]]; then
  gg_step_start "${task}"

  infiles=()
  mapfile -t infiles < <(gg_find_file_basenames "${dir_busco_notung_pep}")
  for infile in "${infiles[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    busco_species_tree_assisted_gene_tree_rooting "${infile}" "${dir_busco_notung_pep}" "${dir_busco_iqtree_pep}" "${dir_busco_rooted_txt_pep}" "${dir_busco_rooted_nwk_pep}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi


task="BUSCO-based Grampa analysis for the polyploidization history with BUSCO DNA trees"
disable_if_no_input_file "run_busco_dupaware_grampa_dna" "${file_dated_species_tree}"
if [[ ! -s "${file_busco_grampa_dna}" && ${run_busco_dupaware_grampa_dna} -eq 1 ]]; then
  gg_step_start "${task}"

  busco_grampa "${dir_busco_rooted_nwk_dna}" "$(dirname "${file_busco_grampa_dna}")" "${file_busco_grampa_dna}"
else
  gg_step_skip "${task}"
fi

task="BUSCO-based Grampa analysis for the polyploidization history with BUSCO protein trees"
disable_if_no_input_file "run_busco_dupaware_grampa_pep" "${file_dated_species_tree}"
if [[ ! -s "${file_busco_grampa_pep}" && ${run_busco_dupaware_grampa_pep} -eq 1 ]]; then
  gg_step_start "${task}"

  busco_grampa "${dir_busco_rooted_nwk_pep}" "$(dirname "${file_busco_grampa_pep}")" "${file_busco_grampa_pep}"
else
  gg_step_skip "${task}"
fi

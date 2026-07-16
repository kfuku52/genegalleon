# shellcheck shell=bash
# Sourced by gg_genome_evolution_core.sh.

task="BUSCO analysis of species-wise input files"
run_shared_species_busco_stage

task="Collecting IDs of common BUSCO genes"
run_shared_busco_summary_stage

task="Generating FASTA files for duplicate-aware BUSCO genes"
ensure_dir "${dir_single_copy_fasta}"
num_busco_ids=$(get_busco_summary_gene_count "${file_species_busco_summary_table}")
singlecopy_fasta_files=()
mapfile -t singlecopy_fasta_files < <(gg_find_file_basenames "${dir_single_copy_fasta}" "${single_copy_fasta_glob}")
num_singlecopy_fasta=${#singlecopy_fasta_files[@]}
if [[ ${num_busco_ids} -ne ${num_singlecopy_fasta} && ${run_extract_species_tree_fasta} -eq 1 ]]; then
  prepare_species_tree_input_dir
  gg_step_start "${task}"

  generate_dupaware_busco_fasta() {
    local busco_id
    busco_id=$(awk -v row="$1" 'NR==row {print $1; exit}' "${file_species_busco_summary_table}")
    local remove_nonsingle=$2
    local outfile1="${dir_single_copy_fasta}/${busco_id}${single_copy_fasta_suffix}"
    if [[ -s "${outfile1}" ]]; then
      return 0
    fi
    local genes=()
    IFS=$'\t' read -r -a genes <<< "$(sed -n "${1}P" "${file_species_busco_summary_table}" | cut -f 4-)"
    echo "busco_id: ${busco_id}"
    if [[ ${strictly_single_copy_only} -eq 1 ]]; then
      local gene_token
      for gene_token in "${genes[@]}"; do
        if [[ "${gene_token}" == "-" ]]; then
          echo "Skipping. ${busco_id} has missing gene(s)."
          return 0
        fi
        if [[ "${gene_token}" == *","* ]]; then
          echo "Skipping. ${busco_id} has duplicated gene(s)."
          return 0
        fi
      done
    fi
    local genes1=()
    read -r -a genes1 <<< "$(printf "%s " "${genes[@]}" | sed -e "s/[[:space:]]-[[:space:]]/ /g" -e "s/[^[[:space:]]]*,[^[[:space:]]]*//g" -e "s/[[:space:]]+/ /g")"
    local num_gene=${#genes1[@]}
    if [[ ${num_gene} -eq 0 ]]; then
      echo "Skipping. ${busco_id} has no genes in the selected species."
      return 0
    fi
    if [[ ! -s "${outfile1}" ]]; then
      if [[ "${input_sequence_mode}" == "protein" ]]; then
        gg_seqkit_grep_by_patterns_from_infile_list 1 "species_tree_input_fasta_list.txt" "${genes1[@]}" |
          seqkit replace --pattern " .*" --replacement "" --ignore-case --threads 1 |
          gg_fasta_relabel_headers_to_species |
          seqkit seq --threads 1 --out-file "${outfile1}"
      else
        local pattern_args=()
        for gene in "${genes1[@]}"; do
          pattern_args+=(--pattern "${gene}")
        done
        seqkit grep --threads 1 "${pattern_args[@]}" --infile-list "species_tree_input_fasta_list.txt" |
          seqkit replace --pattern X --replacement N --by-seq --ignore-case --threads 1 |
          seqkit replace --pattern " .*" --replacement "" --ignore-case --threads 1 |
          cdskit pad |
          gg_fasta_relabel_headers_to_species |
          seqkit seq --threads 1 --out-file "${outfile1}"
      fi
      if [[ ! -s "${outfile1}" ]]; then
        echo "File is empty. Removing: ${outfile1}"
        rm -f -- "${outfile1}"
      fi
    fi
    local fasta_genes=()
    if [[ -s "${outfile1}" ]]; then
      mapfile -t fasta_genes < <(seqkit seq --name --threads 1 "${outfile1}")
    fi
    local num_seq=${#fasta_genes[@]}
    # this block needs to be disabeld for ${strictly_single_copy_only} -eq 1, because orthogroups won't be complete
    if [[ ${num_gene} -ne ${num_seq} && ${strictly_single_copy_only} -eq 1 ]]; then
      echo "${busco_id}: Error. Number of genes and sequences did not match."
      echo "Genes in the orthogroup or BLAST hit:"
      printf '%s\n' "${genes1[@]}"
      echo ""
      echo "Genes in the generated FASTA:"
      printf '%s\n' "${fasta_genes[@]}"
      echo ""
      echo "Check duplicated sequence names in the species input FASTA files. Exiting."
      rm -f -- "${outfile1}"
      exit 1
    fi
  }

  gg_find_fasta_files "${species_tree_input_dir}" 1 > species_tree_input_fasta_list.txt
  num_busco_ids=$(get_busco_summary_gene_count "${file_species_busco_summary_table}")
  for ((i = 2; i <= num_busco_ids + 1; i++)); do # starting from 2 because the line 1 is header.
    wait_until_jobn_le ${GG_TASK_CPUS}
    generate_dupaware_busco_fasta ${i} ${strictly_single_copy_only} &
  done
  wait_for_background_jobs
  rm -f -- species_tree_input_fasta_list.txt
else
  gg_step_skip "${task}"
fi

task="MAFFT alignment"
ensure_dir "${dir_single_copy_mafft}"
num_busco_ids=$(get_busco_summary_gene_count "${file_species_busco_summary_table}")
mafft_fasta_files=()
mapfile -t mafft_fasta_files < <(gg_find_file_basenames "${dir_single_copy_mafft}" "${single_copy_aln_glob}")
num_mafft_fasta=${#mafft_fasta_files[@]}
if [[ ${num_busco_ids} -ne ${num_mafft_fasta} && ${run_individual_mafft} -eq 1 ]]; then
  gg_step_start "${task}"

  run_mafft() {
    local infile=$1
    local infile_base=${infile%%.*}
    local outfile="${dir_single_copy_mafft}/${infile_base}${single_copy_aln_suffix}"
    if [[ -s "${outfile}" ]]; then
      return 0
    fi
    local infile_path="${dir_single_copy_fasta}/${infile}"
    local num_seq
    num_seq=$(gg_count_fasta_records "${infile_path}")
    if [[ ${num_seq} -lt 2 ]]; then
      echo "Skipped. At least 2 sequences are necessary for MAFFT: ${infile}"
      return 0
    fi
    echo "$(date): start mafft: ${infile_base}"
    if [[ "${input_sequence_mode}" == "protein" ]]; then
      seqkit seq --threads 1 "${infile_path}" --out-file "tmp.${infile_base}.input.pep.fasta"
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
      seqkit seq --threads 1 "${infile_path}" --out-file "tmp.${infile_base}.input.cds.fasta"
      cdskit mask \
        --seqfile "tmp.${infile_base}.input.cds.fasta" \
        --outfile "tmp.${infile_base}.cds.fasta"
      seqkit translate \
        --allow-unknown-codon \
        --transl-table "${genetic_code}" \
        --threads 1 \
        "tmp.${infile_base}.cds.fasta" \
        > "tmp.${infile_base}.pep.fasta"
      mafft \
        --auto \
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
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_single_copy_fasta}" "${single_copy_fasta_glob}")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    run_mafft "${input_alignment_file}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

task="TrimAl"
ensure_dir "${dir_single_copy_trimal}"
trimal_fasta_files=()
mapfile -t trimal_fasta_files < <(gg_find_file_basenames "${dir_single_copy_trimal}" "${single_copy_trimal_glob}")
num_trimal_fasta=${#trimal_fasta_files[@]}
if [[ ${num_busco_ids} -ne ${num_trimal_fasta} && ${run_individual_trimal} -eq 1 ]]; then
  gg_step_start "${task}"

  run_trimal() {
    local infile=$1
    local infile_base=${infile%%.*}
    local outfile="${dir_single_copy_trimal}/${infile_base}${single_copy_trimal_suffix}"
    if [[ ! -s "${outfile}" ]]; then
      if [[ "${input_sequence_mode}" == "protein" ]]; then
        seqkit seq --threads 1 "${dir_single_copy_mafft}/${infile}" --out-file "tmp.${infile_base}.pep.aln.fasta"
        trimal \
          -in "tmp.${infile_base}.pep.aln.fasta" \
          -out "tmp.${infile_base}.trimal.fasta" \
          -automated1
      else
        seqkit seq --remove-gaps --threads 1 "${dir_single_copy_mafft}/${infile}" > "tmp.${infile_base}.degap.fasta"
        seqkit translate --transl-table "${genetic_code}" --threads 1 "${dir_single_copy_mafft}/${infile}" > "tmp.${infile_base}.pep.fasta"
        trimal \
          -in tmp.${infile_base}.pep.fasta \
          -backtrans tmp.${infile_base}.degap.fasta \
          -out tmp.${infile_base}.trimal.fasta \
          -ignorestopcodon \
          -automated1
      fi
      if [[ -s "tmp.${infile_base}.trimal.fasta" ]]; then
        seqkit seq --threads 1 "tmp.${infile_base}.trimal.fasta" --out-file "tmp.${infile_base}.trimal.out.fa.gz"
        mv_out "tmp.${infile_base}.trimal.out.fa.gz" "${outfile}"
      fi
      rm -f -- "tmp.${infile_base}."*
    fi
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_single_copy_mafft}" "${single_copy_aln_glob}")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    run_trimal "${input_alignment_file}" &
  done
  wait_for_background_jobs
else
  gg_step_skip "${task}"
fi

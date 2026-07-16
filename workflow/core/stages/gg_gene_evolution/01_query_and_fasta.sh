# shellcheck shell=bash
# Sourced by gg_gene_evolution_core.sh.

task="Species tree availability check"
if [[ ! -s "${species_tree_pruned}" ]]; then
  echo "$(date): Warning: ${task}: species tree file was not found."
  echo "Missing: ${species_tree_pruned}"
else
  gg_step_skip "${task}"
fi

task="Query fasta generation"
if [[ ! -s "${file_og_query_aa_fasta}" && ${run_extract_query_fasta} -eq 1 ]]; then
  gg_step_start "${task}"
  if [[ "$(head -c 1 "${file_query_gene}")" == ">" ]]; then
    seqtype=$(seqkit stats --tabular "${file_query_gene}" | awk 'NR>1 {print $3}')
    if [[ ${seqtype} == "DNA" ]]; then
      echo "DNA sequences were detected. The file will be treated as in-frame CDS sequences, translated into amino acids, and used as a ${query_blast_method} query: ${file_query_gene}"
      seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${GG_TASK_CPUS}" "${file_query_gene}" > "${og_id}.query.aa.tmp.fasta"
      seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.query.aa.tmp.fasta" --out-file "${og_id}.query.aa.out.fa.gz"
      mv_out "${og_id}.query.aa.out.fa.gz" "${file_og_query_aa_fasta}"
      rm -f -- "${og_id}.query.aa.tmp.fasta"
    elif [[ ${seqtype} == "Protein" ]]; then
      echo "Amino acid sequences were detected. The file will be used as a ${query_blast_method} query: ${file_query_gene}"
      seqkit seq --threads "${GG_TASK_CPUS}" "${file_query_gene}" --out-file "${og_id}.query.aa.out.fa.gz"
      mv_out "${og_id}.query.aa.out.fa.gz" "${file_og_query_aa_fasta}"
    else
      echo "Unsupported sequence type '${seqtype}' in '${file_query_gene}'. Only \"DNA\" or \"Protein\" are allowed. Exiting."
      exit 1
    fi
  else
    echo "Gene IDs were detected. Extracting in-frame CDS sequences from species_cds: ${file_query_gene}"
    cp_out "${file_query_gene}" "${dir_output_active}/query_gene/$(basename "${file_query_gene}")"
    mapfile -t genes < <(sed -e '/^[[:space:]]*$/d' "${file_query_gene}")
    mapfile -t cds_files < <(gg_find_fasta_files "${dir_sp_cds}" 1)
    if [[ -e pattern.txt ]]; then
      rm -f -- pattern.txt
    fi
    touch pattern.txt
    for gene in "${genes[@]}"; do
      echo "${gene}" >> pattern.txt
      if [[ "${gene}" == *"−"* ]]; then
        echo "Query sequence name contains minus sign. Searching the sequence name with hyphen as well: ${gene}"
        echo "${gene//−/-}" >> pattern.txt # Replace minus signs ("−") with hyphens ("-") and add to pattern.txt
      fi
    done
    if [[ -e "${og_id}.query.cds.fasta" ]]; then
      rm -f -- "${og_id}.query.cds.fasta"
    fi
    if [[ -e "${og_id}.query.cds.2.fasta" ]]; then
      rm -f -- "${og_id}.query.cds.2.fasta"
    fi
    touch "${og_id}.query.cds.fasta"
    query_hits_tmp_dir="./tmp.query_hits"
    if [[ -d "${query_hits_tmp_dir}" ]]; then
      rm -rf -- "${query_hits_tmp_dir}"
    fi
    mkdir -p "${query_hits_tmp_dir}"
    for file_cds in "${cds_files[@]}"; do
      wait_until_jobn_le ${GG_TASK_CPUS}
      (
        sp_ub=$(gg_species_name_from_path "${file_cds}")
        query_hits_tmp_file="${query_hits_tmp_dir}/$(basename "${file_cds}").hits.fasta"
        seqkit grep --threads "${GG_TASK_CPUS}" --ignore-case --pattern-file <(awk -v sp="${sp_ub}" '{print $0; print sp "_" $0}' pattern.txt) "${file_cds}" |
          sed -e "s/^>${sp_ub}_/>/" -e "s/^>${sp_ub}-/>/" -e "s/^>${sp_ub}[[:space:]]/>/" -e "s/^>${sp_ub}\./>/" -e "s/^>/>${sp_ub}_/" \
            > "${query_hits_tmp_file}"
      ) &
    done
    wait_for_background_jobs
    shopt -s nullglob
    query_hits_tmp_files=("${query_hits_tmp_dir}"/*.hits.fasta)
    shopt -u nullglob
    for query_hits_tmp_file in "${query_hits_tmp_files[@]}"; do
      cat "${query_hits_tmp_file}" >> "${og_id}.query.cds.fasta"
    done
    rm -rf -- "${query_hits_tmp_dir}"
    gg_prepare_cds_fasta_stream "${GG_TASK_CPUS}" "${genetic_code}" < "${og_id}.query.cds.fasta" |
      sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
        > "${og_id}.query.cds.2.fasta"
    num_query=${#genes[@]}
    num_result=$(grep -c -e "^>" "${og_id}.query.cds.2.fasta")
    echo "Number of gene names in query: ${num_query}"
    echo "Number of gene names in extracted fasta: ${num_result}"
    if [[ ${num_query} -ne ${num_result} ]]; then
      echo "Some gene names were not found in species_cds."
      for gene_name in "${genes[@]}"; do
        if ! awk -v gene="${gene_name}" '
	                  /^>/ {
	                    header=$0
	                    sub(/^>/, "", header)
	                    sub(/[[:space:]].*$/, "", header)
	                    if (header == gene) {
	                      found=1
	                      exit
	                    }
	                  }
	                  END { exit(found ? 0 : 1) }
	                ' "${og_id}.query.cds.2.fasta"; then
          echo "Query gene not found in species_cds: ${gene_name}"
        fi
      done
      echo "Exiting."
      exit 1
    fi
    if [[ -s "${og_id}.query.cds.2.fasta" ]]; then
      echo "Translating in-frame CDS sequences to amino acid sequences: ${og_id}.query.cds.2.fasta"
      seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${GG_TASK_CPUS}" "${og_id}.query.cds.2.fasta" > "${og_id}.query.aa.tmp.fasta"
      seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.query.aa.tmp.fasta" --out-file "${og_id}.query.aa.out.fa.gz"
      mv_out "${og_id}.query.aa.out.fa.gz" "${file_og_query_aa_fasta}"
      rm -f -- "${og_id}.query.aa.tmp.fasta"
      rm -f -- "${og_id}.query.cds.2.fasta"
    fi
  fi
else
  gg_step_skip "${task}"
fi

task="In-frame query BLAST (${query_blast_method})"
if [[ ! -s "${file_og_query_blast}" && ${run_query_blast} -eq 1 && "${mode_gene_evolution}" == "query2family" ]]; then
  gg_step_start "${task}"

  if [[ ${query_blast_method} == "tblastn" ]]; then
    if ! type makeblastdb > /dev/null 2>&1; then
      echo "makeblastdb was not found but query_blast_method=tblastn. Exiting."
      exit 1
    fi
    if ! type tblastn > /dev/null 2>&1; then
      echo "tblastn was not found but query_blast_method=tblastn. Exiting."
      exit 1
    fi
  elif [[ ${query_blast_method} == "diamond" ]]; then
    if ! type diamond > /dev/null 2>&1; then
      echo "diamond was not found but query_blast_method=diamond. Exiting."
      exit 1
    fi
    echo "DIAMOND mode selected: species CDS will be translated to proteins because diamond makedb/blastp use protein reference databases."
  fi

  export BLASTDB_LMDB_MAP_SIZE=100000000
  check_species_cds "${gg_workspace_dir}"
  check_if_species_files_unique "${dir_sp_cds}"

  if [[ -e "${og_id}".blastQuery.fasta ]]; then
    rm -f -- "${og_id}.blastQuery.fasta"
  fi
  touch "${og_id}.blastQuery.fasta"

  db_files=()
  ensure_dir "${dir_sp_blastdb}"
  cds_files=()
  mapfile -t cds_files < <(gg_find_fasta_files "${dir_sp_cds}" 1)
  cds_spp=()
  for cds_file in "${cds_files[@]}"; do
    cds_spp+=("$(gg_species_name_from_path "${cds_file}")")
  done
  mapfile -t cds_spp < <(printf "%s\n" "${cds_spp[@]}" | sort -u)
  filter_translated_fasta_for_diamond() {
    awk '
      BEGIN {
        seen = 0
        dropped = 0
        header = ""
        seq = ""
      }
      /^>/ {
        if (seen) {
          if (seq != "") {
            print header
            print seq
          } else {
            dropped++
          }
        }
        header = $0
        seq = ""
        seen = 1
        next
      }
      {
        gsub(/\*/, "", $0)
        if ($0 != "") {
          seq = seq $0
        }
      }
      END {
        if (seen) {
          if (seq != "") {
            print header
            print seq
          } else {
            dropped++
          }
        }
        if (dropped > 0) {
          printf("Dropped %d translated protein records with empty sequence after stop-codon removal.\n", dropped) > "/dev/stderr"
        }
      }
    '
  }
  for sp in "${cds_spp[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    echo "sp: ${sp}"
    sp_cds_candidates=()
    for cds_file in "${cds_files[@]}"; do
      if [[ "$(gg_species_name_from_path "${cds_file}")" == "${sp}" ]]; then
        sp_cds_candidates+=("${cds_file}")
      fi
    done
    if [[ ${#sp_cds_candidates[@]} -eq 0 ]]; then
      echo "No CDS file was found for species: ${sp}. Skipping."
      continue
    fi
    mapfile -t sp_cds_candidates < <(printf "%s\n" "${sp_cds_candidates[@]}" | sort)
    sp_cds=${sp_cds_candidates[0]}
    sp_cds_blastdb="${dir_sp_blastdb}/$(basename "${sp_cds}")"
    db_files+=("${sp_cds_blastdb}")
    if [[ ${query_blast_method} == "tblastn" ]]; then
      echo "makeblastdb input CDS file: ${sp_cds}"
      echo "makeblastdb output database file: ${sp_cds_blastdb}"
      if [[ ! -e "${sp_cds_blastdb}".nhr || ! -e "${sp_cds_blastdb}".nin || ! -e "${sp_cds_blastdb}".nsq || ! -e "${sp_cds_blastdb}".ndb ]]; then
        db_lock_file="${sp_cds_blastdb}.tblastn.build.lock"
        (
          if ! gg_shared_lock_acquire "${db_lock_file}" "TBLASTN database build (${sp})"; then
            exit 1
          fi
          gg_shared_lock_start_heartbeat "${db_lock_file}"
          heartbeat_pid=${GG_SHARED_LOCK_HEARTBEAT_PID:-}
          cleanup_tblastn_db_lock() {
            gg_shared_lock_stop_heartbeat "${heartbeat_pid}"
            gg_shared_lock_release "${db_lock_file}"
          }
          trap cleanup_tblastn_db_lock EXIT
          if [[ ! -e "${sp_cds_blastdb}".nhr || ! -e "${sp_cds_blastdb}".nin || ! -e "${sp_cds_blastdb}".nsq || ! -e "${sp_cds_blastdb}".ndb ]]; then
            if zgrep -q -e "^>.*[[:blank:]]" "${sp_cds}"; then
              echo "Space is detected. Please remove all annotation info after spaces in sequence names. Exiting: ${sp_cds}"
              exit 1
            fi
            if zgrep -q -e "^>.*[|]" "${sp_cds}"; then
              echo "Bar (|) is detected. Bars in sequence names will be replaced with underlines (_): ${sp_cds}"
            fi
            echo "Generating BLAST database: ${sp_cds}"
            echo "Generating BLAST database: ${sp_cds}" >&2
            if [[ ${sp_cds} == *.gz ]]; then
              seqkit seq --threads "${GG_TASK_CPUS}" "${sp_cds}" | makeblastdb -dbtype nucl -title "${sp_cds}" -out "${sp_cds_blastdb}"
            else
              makeblastdb -dbtype nucl -in "${sp_cds}" -out "${sp_cds_blastdb}"
            fi
          fi
        ) || exit 1
      fi
    elif [[ ${query_blast_method} == "diamond" ]]; then
      sp_cds_diamond_fasta="${sp_cds_blastdb}.diamond.fasta"
      echo "diamond input CDS file: ${sp_cds}"
      echo "diamond translated protein file: ${sp_cds_diamond_fasta}"
      echo "diamond database file: ${sp_cds_blastdb}.dmnd"
      if [[ ! -e "${sp_cds_blastdb}".dmnd ]]; then
        db_lock_file="${sp_cds_blastdb}.diamond.build.lock"
        (
          if ! gg_shared_lock_acquire "${db_lock_file}" "DIAMOND database build (${sp})"; then
            exit 1
          fi
          gg_shared_lock_start_heartbeat "${db_lock_file}"
          heartbeat_pid=${GG_SHARED_LOCK_HEARTBEAT_PID:-}
          cleanup_diamond_db_lock() {
            gg_shared_lock_stop_heartbeat "${heartbeat_pid}"
            gg_shared_lock_release "${db_lock_file}"
          }
          trap cleanup_diamond_db_lock EXIT
          if [[ ! -e "${sp_cds_blastdb}".dmnd ]]; then
            if zgrep -q -e "^>.*[[:blank:]]" "${sp_cds}"; then
              echo "Space is detected. Please remove all annotation info after spaces in sequence names. Exiting: ${sp_cds}"
              exit 1
            fi
            if zgrep -q -e "^>.*[|]" "${sp_cds}"; then
              echo "Bar (|) is detected. Bars in sequence names will be replaced with underlines (_): ${sp_cds}"
            fi
            echo "Generating DIAMOND database: ${sp_cds}"
            echo "Generating DIAMOND database: ${sp_cds}" >&2
            if [[ ${sp_cds} == *.gz ]]; then
              seqkit seq --remove-gaps --threads "${GG_TASK_CPUS}" "${sp_cds}" |
                seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${GG_TASK_CPUS}" |
                filter_translated_fasta_for_diamond \
                  > "${sp_cds_diamond_fasta}"
            else
              seqkit seq --remove-gaps --threads "${GG_TASK_CPUS}" "${sp_cds}" |
                seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${GG_TASK_CPUS}" |
                filter_translated_fasta_for_diamond \
                  > "${sp_cds_diamond_fasta}"
            fi
            if [[ "$(head -c 1 "${sp_cds_diamond_fasta}")" != '>' ]]; then
              sed -e "1d" "${sp_cds_diamond_fasta}" > "${sp_cds_diamond_fasta}.tmp"
              mv_out "${sp_cds_diamond_fasta}.tmp" "${sp_cds_diamond_fasta}"
            fi
            if [[ ! -s "${sp_cds_diamond_fasta}" ]]; then
              echo "Translated FASTA for DIAMOND is empty: ${sp_cds_diamond_fasta}. Exiting."
              exit 1
            fi
            if ! diamond makedb --in "${sp_cds_diamond_fasta}" --db "${sp_cds_blastdb}"; then
              echo "diamond makedb failed for ${sp_cds}. Exiting."
              exit 1
            fi
            rm -f -- "${sp_cds_diamond_fasta}"
          fi
        ) || exit 1
      fi
    fi
  done
  wait_for_background_jobs
  echo "db_files: ${db_files[*]}"
  query_aa_local="${og_id}.query.aa.tmp.for_blast.fasta"
  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_query_aa_fasta}" --out-file "${query_aa_local}"
  resolve_query_blast_evalue "${query_blast_evalue}" "${query_aa_local}" "${query_blast_auto_evalue_maxlen_cutoffs}"
  if [[ "$(printf '%s' "${query_blast_evalue}" | tr '[:upper:]' '[:lower:]')" == "auto" ]]; then
    echo "query BLAST auto E-value: query_count=${query_blast_query_num_seqs} min_aa_len=${query_blast_query_min_aa_len} avg_aa_len=${query_blast_query_avg_aa_len} max_aa_len=${query_blast_query_max_aa_len}"
    echo "query BLAST auto E-value: cutoffs=${query_blast_auto_evalue_maxlen_cutoffs} effective_query_blast_evalue=${effective_query_blast_evalue}"
  else
    echo "query BLAST E-value: effective_query_blast_evalue=${effective_query_blast_evalue}"
  fi

  outfmt="qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore frames qlen slen"
  if [[ ${query_blast_method} == "tblastn" ]]; then
    db_files_str=$(printf " %s" "${db_files[@]}")
    db_files_str="${db_files_str# }"
    echo "Running tblastn."
    if ! tblastn \
      -query "${query_aa_local}" \
      -db "${db_files_str}" \
      -out blast_out.tsv \
      -db_gencode "${genetic_code}" \
      -evalue "${effective_query_blast_evalue}" \
      -max_target_seqs 50000 \
      -outfmt "6 ${outfmt}" \
      -num_threads "${GG_TASK_CPUS}"; then
      echo "tblastn failed. Exiting."
      exit 1
    fi
  elif [[ ${query_blast_method} == "diamond" ]]; then
    echo "Running diamond blastp."
    rm -f -- blast_out.tsv
    touch blast_out.tsv
    for db_file in "${db_files[@]}"; do
      if [[ ! -e "${db_file}".dmnd ]]; then
        echo "DIAMOND database file is missing: ${db_file}.dmnd. Exiting."
        exit 1
      fi
      tmp_diamond_out="$(basename "${db_file}").diamond.out.tsv"
      if ! diamond blastp \
        --query "${query_aa_local}" \
        --db "${db_file}" \
        --out "${tmp_diamond_out}" \
        --evalue "${effective_query_blast_evalue}" \
        --max-target-seqs 50000 \
        --threads "${GG_TASK_CPUS}" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen; then
        echo "diamond blastp failed for database: ${db_file}. Exiting."
        exit 1
      fi
      if [[ -s "${tmp_diamond_out}" ]]; then
        awk -F '\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,"0/1",$13,$14}' "${tmp_diamond_out}" >> blast_out.tsv
      fi
      rm -f -- "${tmp_diamond_out}"
    done
  fi
  rm -f -- "${query_aa_local}"

  python "${gg_support_dir}/annotate_blast_coverage.py" \
    --in blast_out.tsv \
    --ncpu "${GG_TASK_CPUS}" \
    --outfmt-columns "${outfmt}" \
    --frame-filter "0/1" \
    --out blast_out_inframe.tmp3.tsv

  if [[ -s blast_out_inframe.tmp3.tsv ]]; then
    mv_out blast_out_inframe.tmp3.tsv "${file_og_query_blast}"
  else
    echo "No query BLAST hits were detected after in-frame filtering. Exiting."
    exit 1
  fi
else
  gg_step_skip "${task}"
fi

task="Fasta generation"
if [[ ! -s "${file_og_primary_fasta}" && ${run_extract_primary_fasta} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ "${mode_gene_evolution}" == "orthogroup" ]]; then
    genes=()
    read -r -a genes <<< "$(awk -v og="${og_id}" '$1==og {$1=""; sub(/^[[:space:]]+/, "", $0); gsub(",", "", $0); gsub(/\t/, " ", $0); sub(/[[:space:]]*$/, "", $0); gsub(/\047|"/, "", $0); print; exit}' "${file_og}")"
  elif [[ "${mode_gene_evolution}" == "query2family" ]]; then
    python "${gg_support_dir}/extract_gene_id_from_blast_table.py" \
      --infile "${file_og_query_blast}" \
      --outfile gene_id_list.txt \
      --min_query_blast_coverage "${query_blast_coverage}" \
      --max_num_gene_blast_hit_retrieval "${max_num_gene_blast_hit_retrieval}"
    mapfile -t genes < gene_id_list.txt
  fi
  if [[ -e pattern.txt ]]; then
    rm -f -- pattern.txt
  fi
  touch pattern.txt
  for gene in "${genes[@]}"; do
    echo "${gene}" >> pattern.txt
  done

  num_gene=${#genes[@]}
  if [[ "${input_sequence_mode}" == "protein" ]] && species_protein_input_has_files; then
    protein_files=()
    check_species_protein_dir "${dir_sp_protein_input}"
    check_if_species_files_unique "${dir_sp_protein_input}"
    mapfile -t protein_files < <(gg_find_fasta_files "${dir_sp_protein_input}" 1)
    echo "Number of protein files in ${dir_sp_protein_input}: ${#protein_files[@]}"
    if [[ -s "${file_species_genetic_code}" ]]; then
      echo "species_genetic_code.tsv is ignored because species_protein inputs are provided: ${file_species_genetic_code}"
    fi
    rm -f -- "${og_id}.pep.fasta"
    touch "${og_id}.pep.fasta"
    for protein_path in "${protein_files[@]}"; do
      seqkit grep --threads "${GG_TASK_CPUS}" --pattern-file pattern.txt "${protein_path}" \
        >> "${og_id}.pep.fasta"
    done
    seqkit replace --pattern " .*" --replacement "" --ignore-case --threads "${GG_TASK_CPUS}" "${og_id}.pep.fasta" |
      seqkit replace --pattern "\+" --replacement "_" --ignore-case --threads "${GG_TASK_CPUS}" |
      sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
        > "${og_id}.pep.2.fasta"
    fasta_genes=()
    mapfile -t fasta_genes < <(awk '/^>/ {sub(/^>/, "", $0); print}' "${og_id}.pep.2.fasta")
    num_seq=${#fasta_genes[@]}
    echo "Number of genes in the orthogroup or BLAST hit: ${num_gene}"
    echo "Number of sequences in the protein fasta: ${num_seq}"
    if [[ ${num_gene} -eq ${num_seq} ]]; then
      echo "Number of genes and sequences matched. Protein fasta generation completed!"
      seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.pep.2.fasta" --out-file "${og_id}.pep.out.fa.gz"
      mv_out "${og_id}.pep.out.fa.gz" "${file_og_pep_fasta}"
    else
      echo "Number of genes and sequences did not match."
      echo "Genes in the orthogroup or BLAST hit:"
      printf '%s\n' "${genes[@]}"
      echo ""
      echo "Genes in the generated protein FASTA:"
      printf '%s\n' "${fasta_genes[@]}" | sort | tr '\n' ' '
      echo ""
      echo "There may be duplicated or missing sequences."
      echo "Exiting."
      exit 1
    fi
  else
    if [[ ! -s "${file_og_cds_fasta}" ]]; then
      cds_files=()
      mapfile -t cds_files < <(gg_find_fasta_files "${dir_sp_cds}" 1)
      echo "Number of CDS files in ${dir_sp_cds}: ${#cds_files[@]}"
      if [[ ${#cds_files[@]} -eq 0 ]]; then
        echo "No species_cds FASTA files were found for focal fasta generation: ${dir_sp_cds}"
        exit 1
      fi
      if [[ -e "${og_id}.cds.fasta" ]]; then
        rm -f -- "${og_id}.cds.fasta"
      fi
      touch "${og_id}.cds.fasta"
      for file_cds in "${cds_files[@]}"; do
        seqkit grep --threads "${GG_TASK_CPUS}" --pattern-file pattern.txt "${file_cds}" \
          >> "${og_id}.cds.fasta"
      done

      seqkit replace --pattern "X" --replacement "N" --by-seq --ignore-case --threads "${GG_TASK_CPUS}" "${og_id}.cds.fasta" |
        seqkit replace --pattern " .*" --replacement "" --ignore-case --threads "${GG_TASK_CPUS}" |
        seqkit replace --pattern "\+" --replacement "_" --ignore-case --threads "${GG_TASK_CPUS}" |
        cdskit pad --codontable "${genetic_code}" |
        sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
          > "${og_id}.cds.2.fasta"

      fasta_genes=()
      mapfile -t fasta_genes < <(awk '/^>/ {sub(/^>/, "", $0); print}' "${og_id}.cds.2.fasta")
      num_seq=${#fasta_genes[@]}
      echo "Number of genes in the orthogroup or BLAST hit: ${num_gene}"
      echo "Number of sequences in the fasta: ${num_seq}"
      if [[ ${num_gene} -eq ${num_seq} ]]; then
        echo "Number of genes and sequences matched. Fasta generation completed!"
        seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.cds.2.fasta" --out-file "${og_id}.cds.out.fa.gz"
        mv_out "${og_id}.cds.out.fa.gz" "${file_og_cds_fasta}"
      else
        echo "Number of genes and sequences did not match."
        echo "Genes in the orthogroup or BLAST hit:"
        printf '%s\n' "${genes[@]}"
        echo ""
        echo "Genes in the generated FASTA:"
        printf '%s\n' "${fasta_genes[@]}" | sort | tr '\n' ' '
        echo ""
        echo "There may be duplicated or missing sequences."
        echo "If you have recently replaced species_cds files, please make sure to remove species_cds_blastdb before rerunning."
        echo "Exiting."
        exit 1
      fi
    fi
    if [[ "${input_sequence_mode}" == "protein" ]]; then
      gg_prepare_species_genetic_code_table "${dir_sp_cds}" "${genetic_code}" "${file_species_genetic_code_resolved}" "${file_species_genetic_code}"
      translate_orthogroup_cds_to_protein_fasta "${file_og_cds_fasta}" "${file_og_pep_fasta}" "${file_species_genetic_code_resolved}"
    fi
  fi
else
  gg_step_skip "${task}"
fi

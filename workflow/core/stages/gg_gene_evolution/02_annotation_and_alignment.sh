# shellcheck shell=bash
# Sourced by gg_gene_evolution_core.sh.

task="Protein RPS-BLAST"
disable_if_no_input_file "run_rps_blast" "${file_og_primary_fasta}"
if [[ ! -s "${file_og_rpsblast}" && ${run_rps_blast} -eq 1 ]]; then
  gg_step_start "${task}"
  if ! dir_rpsblastdb=$(ensure_pfam_le_db "${gg_workspace_dir}"); then
    echo "Failed to prepare Pfam_LE DB. Exiting."
    exit 1
  fi
  if [[ -s "${dir_rpsblastdb}/Pfam.pal" ]]; then
    db_rpsblast="${dir_rpsblastdb}/Pfam"
  else
    rps_db_candidates=()
    mapfile -t rps_db_candidates < <(find "${dir_rpsblastdb}" -maxdepth 1 -type f -name "*.loo" | sort)
    if [[ ${#rps_db_candidates[@]} -eq 0 ]]; then
      echo "No RPS-BLAST DB index (*.loo) was found in: ${dir_rpsblastdb}. Exiting."
      exit 1
    fi
    if [[ ${#rps_db_candidates[@]} -gt 1 ]]; then
      echo "Multiple RPS-BLAST DB indices were found. Using the first sorted entry: ${rps_db_candidates[0]}"
    fi
    db_rpsblast="${rps_db_candidates[0]%.loo}"
  fi
  echo "db_rpsblast: ${db_rpsblast}"

  if [[ -e "${og_id}.rpsblast.tmp.tsv" ]]; then
    rm -f -- "${og_id}.rpsblast.tmp.tsv"
  fi

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    seqkit seq --remove-gaps --threads "${GG_TASK_CPUS}" "${file_og_pep_fasta}" > ungapped_translated_cds.fas
  else
    seqkit seq --remove-gaps --threads "${GG_TASK_CPUS}" "${file_og_cds_fasta}" |
      seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${GG_TASK_CPUS}" \
        > ungapped_translated_cds.fas
  fi

  if [[ "$(head -c 1 ungapped_translated_cds.fas)" != '>' ]]; then
    sed -e "1d" ungapped_translated_cds.fas > ungapped_translated_cds2.fas
    mv_out ungapped_translated_cds2.fas ungapped_translated_cds.fas
  fi

  outfmt="qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle"

  if ! rpsblast \
    -query ungapped_translated_cds.fas \
    -db "${db_rpsblast}" \
    -out "${og_id}.rpsblast.tmp.tsv" \
    -evalue 0.01 \
    -outfmt "6 ${outfmt}" \
    -num_threads "${GG_TASK_CPUS}"; then
    echo "RPS-BLAST failed. Exiting."
    exit 1
  fi

  genes=()
  mapfile -t genes < <(awk '/^>/ {sub(/^>/, "", $0); sub(/^[[:space:]]*/, "", $0); sub(/[[:space:]].*$/, "", $0); print}' ungapped_translated_cds.fas)
  for gene in "${genes[@]}"; do
    if ! awk -F '\t' -v gene="${gene}" '$1 == gene {found=1; exit} END {exit(found ? 0 : 1)}' "${og_id}.rpsblast.tmp.tsv"; then
      echo "${gene}: no hit in RPS-BLAST. Appending qlen to output tsv."
      qlen=$(seqkit fx2tab --length ungapped_translated_cds.fas | awk -F '\t' -v gene="${gene}" '$1 == gene {print $NF}')
      printf '%s\t\t\t\t\t\t\t\t\t\t\t\t%s\t\t\n' "${gene}" "${qlen}" >> "${og_id}.rpsblast.tmp.tsv"
    else
      echo "${gene}: RPS-BLAST hit found."
    fi
  done
  {
    printf '%s\n' "${outfmt}" | tr ' ' '\t'
    cat "${og_id}.rpsblast.tmp.tsv"
  } > "${og_id}.rpsblast.tsv"

  cp_out "${og_id}.rpsblast.tsv" "${file_og_rpsblast}"
else
  gg_step_skip "${task}"
fi

task="Gene trait extraction from gff files"
disable_if_no_input_file "run_collect_gff_info" "${file_og_primary_fasta}"
if [[ ! -s "${file_og_gff_info}" && ${run_collect_gff_info} -eq 1 ]]; then
  gg_step_start "${task}"
  if [[ -e gff2genestat.tsv ]]; then
    rm -f -- gff2genestat.tsv
  fi
  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_primary_fasta}" --out-file "${og_id}.gff2genestat_input.fasta"

  python "${gg_support_dir}/gff2genestat.py" \
    --dir_gff "${dir_sp_gff}" \
    --feature "CDS" \
    --multiple_hits "longest" \
    --seqfile "${og_id}.gff2genestat_input.fasta" \
    --ncpu "${GG_TASK_CPUS}" \
    --outfile gff2genestat.tsv
  rm -f -- "${og_id}.gff2genestat_input.fasta"

  if [[ -s gff2genestat.tsv ]]; then
    mv_out gff2genestat.tsv "${file_og_gff_info}"
  fi
else
  gg_step_skip "${task}"
fi

task="UniProt annotation (${uniprot_annotation_method})"
disable_if_no_input_file "run_uniprot_annotation" "${file_og_primary_fasta}"
if [[ ! -s "${file_og_uniprot_annotation}" && ${run_uniprot_annotation} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    seqkit seq --remove-gaps --only-id --threads "${GG_TASK_CPUS}" "${file_og_pep_fasta}" > uniprot.query.pep.fas
  else
    seqkit seq --remove-gaps --only-id --threads "${GG_TASK_CPUS}" "${file_og_cds_fasta}" |
      seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${GG_TASK_CPUS}" \
        > uniprot.query.pep.fas
  fi

  if [[ "${uniprot_annotation_method}" == "blastp" ]]; then
    if ! uniprot_db_prefix=$(ensure_uniprot_sprot_blast_db "${gg_workspace_dir}"); then
      echo "Failed to prepare UniProt Swiss-Prot BLASTP DB. Exiting."
      exit 1
    fi
    if ! validate_uniprot_sprot_db_prefix "${uniprot_db_prefix}" "blastp"; then
      echo "Invalid UniProt Swiss-Prot BLASTP DB prefix. Exiting."
      exit 1
    fi

    blastp \
      -query uniprot.query.pep.fas \
      -num_threads "${GG_TASK_CPUS}" \
      -db "${uniprot_db_prefix}" \
      -out uniprot.search.tsv \
      -outfmt "6 qseqid sseqid pident length evalue bitscore qlen" \
      -max_target_seqs 1 \
      -evalue 1e-2
  else
    if ! uniprot_db_prefix=$(ensure_uniprot_sprot_mmseqs_db "${gg_workspace_dir}"); then
      echo "Failed to prepare UniProt Swiss-Prot MMseqs2 DB. Exiting."
      exit 1
    fi
    if ! validate_uniprot_sprot_db_prefix "${uniprot_db_prefix}" "mmseqs2"; then
      echo "Invalid UniProt Swiss-Prot MMseqs2 DB prefix. Exiting."
      exit 1
    fi

    mmseqs createdb "uniprot.query.pep.fas" "uniprot.queryDB"
    mmseqs search "uniprot.queryDB" "${uniprot_db_prefix}.mmseqs" "uniprot.resultDB" "tmp_mmseqs2_uniprot" \
      --threads "${GG_TASK_CPUS}" \
      --split-memory-limit "$(gg_memory_fraction_gb "${GG_MEM_TOOL_GB}" 3 4)G" \
      --max-seqs 1 \
      -e 1e-2 \
      -s 7.5
    mmseqs convertalis "uniprot.queryDB" "${uniprot_db_prefix}.mmseqs" "uniprot.resultDB" "uniprot.search.tsv" \
      --threads "${GG_TASK_CPUS}" \
      --format-output "query,target,pident,alnlen,evalue,bits,qlen"
    rm -f -- uniprot.queryDB* uniprot.resultDB*
    rm -rf -- "tmp_mmseqs2_uniprot"
  fi

  uniprot_meta_tsv=""
  if ! uniprot_meta_tsv=$(ensure_uniprot_sprot_metadata_tsv "${gg_workspace_dir}" "${uniprot_db_prefix}" 2>/dev/null); then
    echo "Warning: UniProt Swiss-Prot metadata TSV is unavailable for prefix: ${uniprot_db_prefix}" >&2
    uniprot_meta_tsv=""
  fi

  python "${gg_support_dir}/reformat_uniprot_diamond.py" \
    --diamond_tsv uniprot.search.tsv \
    --query_fasta uniprot.query.pep.fas \
    --uniprot_fasta "${uniprot_db_prefix}.pep" \
    --uniprot_meta_tsv "${uniprot_meta_tsv}" \
    --outfile uniprot.annotation.tsv

  cp_out uniprot.annotation.tsv "${file_og_uniprot_annotation}"
else
  gg_step_skip "${task}"
fi

task="cdskit localize"
disable_if_no_input_file "run_cdskit_localize" "${file_og_primary_fasta}"
if [[ ! -s "${file_og_cdskit_localize}" && ${run_cdskit_localize} -eq 1 ]]; then
  gg_step_start "${task}"

  cdskit_localize_species_dir="${dir_sp_cds}"
  if [[ "${input_sequence_mode}" == "protein" ]] && species_protein_input_has_files; then
    cdskit_localize_species_dir="${dir_sp_protein_input}"
  fi
  cdskit_localize_species_names=()
  mapfile -t cdskit_localize_species_names < <(gg_species_names_from_fasta_dir "${cdskit_localize_species_dir}")
  cdskit_localize_group_resolved=$(
    gg_resolve_cdskit_localize_organism_group \
      "${cdskit_localize_organism_group}" \
      "${gg_workspace_dir}" \
      "${busco_lineage}" \
      "${cdskit_localize_species_names[@]}"
  )

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    seqkit seq --remove-gaps --only-id --threads "${GG_TASK_CPUS}" "${file_og_pep_fasta}" > "cdskit_localize.input.pep.fasta"
    cdskit_localize_input="cdskit_localize.input.pep.fasta"
    cdskit_localize_seqtype="protein"
  else
    gg_prepare_cdskit_localize_cds_input \
      "${file_og_cds_fasta}" \
      "cdskit_localize.input.cds.fasta" \
      "${GG_TASK_CPUS}" \
      "${genetic_code}"
    cdskit_localize_input="cdskit_localize.input.cds.fasta"
    cdskit_localize_seqtype="dna"
  fi

  gg_run_cdskit_localize \
    "${cdskit_localize_input}" \
    "${cdskit_localize_seqtype}" \
    "cdskit_localize.tsv" \
    "${cdskit_localize_model}" \
    "${cdskit_localize_group_resolved}" \
    "${cdskit_localize_include_features}" \
    "${cdskit_localize_no_model_download}" \
    "${GG_TASK_CPUS}" \
    "${genetic_code}"
  if [[ -s "cdskit_localize.tsv" ]]; then
    mv_out "cdskit_localize.tsv" "${file_og_cdskit_localize}"
  fi
  rm -f -- "cdskit_localize.input.cds.fasta" "cdskit_localize.input.pep.fasta"
else
  gg_step_skip "${task}"
fi

task="In-frame mafft alignment"
disable_if_no_input_file "run_mafft" "${file_og_primary_fasta}"
if [[ ! -s "${file_og_mafft}" && ${run_mafft} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_pep_fasta}" --out-file tmp.pep.input.fasta
    mafft \
      --auto \
      --amino \
      --thread "${GG_TASK_CPUS}" \
      --quiet \
      tmp.pep.input.fasta \
      > "${og_id}.cds.aln.fasta"
  else
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_cds_fasta}" --out-file tmp.cds.input.fasta
    cdskit mask --seqfile tmp.cds.input.fasta --codontable "${genetic_code}" --outfile tmp.cds.fasta

    seqkit translate \
      --allow-unknown-codon \
      --transl-table "${genetic_code}" \
      --threads "${GG_TASK_CPUS}" \
      tmp.cds.fasta \
      > tmp.pep.fasta

    mafft \
      --auto \
      --amino \
      --thread "${GG_TASK_CPUS}" \
      --quiet \
      tmp.pep.fasta \
      > tmp.pep.aln.fasta

    cdskit backalign \
      --seqfile tmp.cds.fasta \
      --aa_aln tmp.pep.aln.fasta \
      --codontable "${genetic_code}" \
      --outfile "${og_id}.cds.aln.fasta"
  fi

  seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.cds.aln.fasta" --out-file "${og_id}.cds.aln.out.fa.gz"
  mv_out "${og_id}.cds.aln.out.fa.gz" "${file_og_mafft}"
  rm -f -- tmp.cds.input.fasta tmp.cds.fasta tmp.pep.fasta tmp.pep.aln.fasta tmp.pep.input.fasta
else
  gg_step_skip "${task}"
fi

task="AMAS for original alignment"
disable_if_no_input_file "run_amas_original" "${file_og_untrimmed_aln_analysis}"
if [[ ! -s "${file_og_amas_original}" && ${run_amas_original} -eq 1 ]]; then
  gg_step_start "${task}"
  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_untrimmed_aln_analysis}" --out-file "${og_id}.amas.original.input.fasta"

  AMAS.py summary \
    --in-format fasta \
    --data-type "${amas_data_type}" \
    --in-files "${og_id}.amas.original.input.fasta"

  mv_out summary.txt "${file_og_amas_original}"
  rm -f -- "${og_id}.amas.original.input.fasta"
else
  gg_step_skip "${task}"
fi

task="MaxAlign"
disable_if_no_input_file "run_maxalign" "${file_og_untrimmed_aln_analysis}"
if [[ ! -s "${file_og_maxalign}" && ${run_maxalign} -eq 1 ]]; then
  gg_step_start "${task}"

  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_untrimmed_aln_analysis}" --out-file "${og_id}.cds.aln.fasta"
  maxalign_keep_regex=""
  if [[ "${mode_gene_evolution}" == "query2family" && ${retain_query_in_maxalign} -eq 0 ]]; then
    echo "Query sequence(s) is NOT necessarily retained in MaxAlign."
  elif [[ "${mode_gene_evolution}" == "query2family" && ${retain_query_in_maxalign} -eq 1 ]]; then
    echo "Query sequence(s) is retained in MaxAlign."
    maxalign_keep_regex=$(
      python - "${file_query_gene}" << 'PY'
import re
import sys

infile = sys.argv[1]
gene_ids = []
with open(infile, 'r', encoding='utf-8', errors='replace') as handle:
    first_char = handle.read(1)
    handle.seek(0)
    if first_char == '>':
        for line in handle:
            if not line.startswith('>'):
                continue
            gene = line[1:].strip().split()[0]
            if gene:
                gene_ids.append(gene)
    else:
        for line in handle:
            gene = line.strip()
            if gene:
                gene_ids.append(gene)

normalized_ids = []
seen = set()
for gene in gene_ids:
    for candidate in (gene, gene.replace('−', '-')):
        if candidate and candidate not in seen:
            seen.add(candidate)
            normalized_ids.append(candidate)

patterns = [f"(?i:^{re.escape(gene)}$)" for gene in normalized_ids]
print(','.join(patterns))
PY
    )
    if [[ -z "${maxalign_keep_regex}" ]]; then
      echo "Warning: No query IDs were parsed for MaxAlign --keep. Running without keep constraints."
    fi
  else
    maxalign_keep_regex=""
  fi

  maxalign_cmd=(
    cdskit maxalign
    --seqfile "${og_id}.cds.aln.fasta"
    --outfile "${og_id}.maxalign.output.fasta"
  )
  if [[ -n "${maxalign_keep_regex}" ]]; then
    maxalign_cmd+=(--keep "${maxalign_keep_regex}")
  fi
  "${maxalign_cmd[@]}"

  echo "Number of sequences before MaxAlign: $(gg_count_fasta_records "${og_id}.cds.aln.fasta")"
  echo "Number of sequences after MaxAlign: $(gg_count_fasta_records "${og_id}.maxalign.output.fasta")"

  seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.maxalign.output.fasta" --out-file "${og_id}.maxalign.out.fa.gz"
  mv_out "${og_id}.maxalign.out.fa.gz" "${file_og_maxalign}"
  rm -f -- "${og_id}.maxalign.output.fasta"
else
  gg_step_skip "${task}"
fi
if [[ ${run_maxalign} -eq 1 ]]; then
  switch_alignment_analysis_source "${file_og_maxalign}"
fi

task="TrimAl"
disable_if_no_input_file "run_trimal" "${file_og_untrimmed_aln_analysis}"
if [[ ! -s "${file_og_trimal}" && ${run_trimal} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_untrimmed_aln_analysis}" --out-file untrimmed.pep.fasta
    trimal \
      -in untrimmed.pep.fasta \
      -out "${og_id}.cds.trimal.tmp2.fasta" \
      -automated1
  else
    seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${GG_TASK_CPUS}" "${file_og_untrimmed_aln_analysis}" |
      sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
        > untrimmed.pep.fasta

    seqkit seq --remove-gaps --threads "${GG_TASK_CPUS}" \
      "${file_og_untrimmed_aln_analysis}" \
      > untrimmed.cds.degap.fasta

    trimal \
      -in untrimmed.pep.fasta \
      -backtrans untrimmed.cds.degap.fasta \
      -out "${og_id}.cds.trimal.tmp1.fasta" \
      -ignorestopcodon \
      -automated1

    cdskit rmseq \
      --seqfile "${og_id}.cds.trimal.tmp1.fasta" \
      --problematic_percent 100 |
      cdskit hammer \
        --seqfile "-" \
        --codontable "${genetic_code}" \
        --nail 4 \
        --outfile "${og_id}.cds.trimal.tmp2.fasta"
  fi

  if [[ -s "${og_id}.cds.trimal.tmp2.fasta" ]]; then
    echo "Copying. Output file detected for the task: ${task}"
    seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.cds.trimal.tmp2.fasta" --out-file "${og_id}.cds.trimal.out.fa.gz"
    mv_out "${og_id}.cds.trimal.out.fa.gz" "${file_og_trimal}"
  fi
else
  gg_step_skip "${task}"
fi
if [[ ${run_trimal} -eq 1 ]]; then
  set_analysis_file trimmed_aln "${file_og_trimal}"
fi

task="ClipKIT"
disable_if_no_input_file "run_clipkit" "${file_og_untrimmed_aln_analysis}"
if [[ ! -s "${file_og_clipkit}" && ${run_clipkit} -eq 1 ]]; then
  gg_step_start "${task}"

  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_untrimmed_aln_analysis}" --out-file "${og_id}.cds.clipkit.input.fasta"

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    clipkit \
      "${og_id}.cds.clipkit.input.fasta" \
      --mode smart-gap \
      --sequence_type aa \
      --input_file_format "fasta" \
      --output_file_format "fasta" \
      --output "${og_id}.cds.clipkit.hammer.fasta" \
      --log
  else
    clipkit \
      "${og_id}.cds.clipkit.input.fasta" \
      --mode smart-gap \
      --sequence_type nt \
      --codon \
      --input_file_format "fasta" \
      --output_file_format "fasta" \
      --output "${og_id}.cds.clipkit.tmp.fasta" \
      --log

    cdskit hammer \
      --codontable "${genetic_code}" \
      --nail 4 \
      --seqfile "${og_id}.cds.clipkit.tmp.fasta" |
      cdskit rmseq \
        --problematic_percent 100 \
        --outfile "${og_id}.cds.clipkit.hammer.fasta"
  fi

  if [[ -s "${og_id}.cds.clipkit.hammer.fasta" ]]; then
    echo "Copying. Output file detected for the task: ${task}"
    seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.cds.clipkit.hammer.fasta" --out-file "${og_id}.cds.clipkit.out.fa.gz"
    mv_out "${og_id}.cds.clipkit.out.fa.gz" "${file_og_clipkit}"
    if [[ -s "${og_id}.cds.clipkit.tmp.fasta.log" ]]; then
      cp_out "${og_id}.cds.clipkit.tmp.fasta.log" "${file_og_clipkit_log}"
    elif [[ -s "${og_id}.cds.clipkit.hammer.fasta.log" ]]; then
      cp_out "${og_id}.cds.clipkit.hammer.fasta.log" "${file_og_clipkit_log}"
    fi
  fi
  rm -f -- "${og_id}.cds.clipkit.input.fasta"
else
  gg_step_skip "${task}"
fi
if [[ ${run_clipkit} -eq 1 ]]; then
  set_analysis_file trimmed_aln "${file_og_clipkit}"
fi

task="AMAS for cleaned alignment"
disable_if_no_input_file "run_amas_cleaned" "${file_og_trimmed_aln_analysis}"
if [[ ! -s "${file_og_amas_cleaned}" && ${run_amas_cleaned} -eq 1 ]]; then
  gg_step_start "${task}"
  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file "${og_id}.amas.cleaned.input.fasta"

  AMAS.py summary \
    --in-format fasta \
    --data-type "${amas_data_type}" \
    --in-files "${og_id}.amas.cleaned.input.fasta"

  mv_out summary.txt "${file_og_amas_cleaned}"
  rm -f -- "${og_id}.amas.cleaned.input.fasta"
else
  gg_step_skip "${task}"
fi
if [[ ${run_maxalign} -eq 1 ]]; then
  # This code block should be placed immediately after "AMAS for cleaned alignment".
  # orthogroup_summary.tsv will not include necessary info otherwise.
  num_gene_before_maxalign=$(gg_count_fasta_records "${file_og_mafft}")
  num_gene_after_maxalign=$(gg_count_fasta_records "${file_og_maxalign}")
  echo "Number of genes before MaxAlign: ${num_gene_before_maxalign}"
  echo "Number of genes after MaxAlign: ${num_gene_after_maxalign}"
  if [[ ${num_gene_after_maxalign} -lt 3 ]]; then
    echo "This is not sufficient for tree-based analysis (<3). Exiting."
    exit 1
  fi
fi

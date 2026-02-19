#!/usr/bin/env bash

### Start: Modify this block to tailor your analysis ###

# CDS analysis
run_get_gff_info=0 # Collect gene information from workspace/input/species_gff
run_busco_cds=0 # Gene set completeness analysis. https://gitlab.com/ezlab/busco
run_uniprot_annotation=0 # DIAMOND-based CDS annotation against UniProt Swiss-Prot.
run_cds_fx2tab=0 # Compute sequence stats of the CDS sequences. https://bioinf.shenwei.me/seqkit/usage/#fx2tab-tab2fx
run_cds_mmseqs2taxonomy=0 # Taxonomic assignment of CDS sequences. 16CPUS * 8GRAM should work. https://github.com/timkahlke/BASTA
run_cds_contamination_removal=0 # Removal of taxonomically non-compatible contaminated seqeunces using the result of MMseqs2 Taxonomy.
run_annotation=0 # Summarize per-gene annotation info including BUSCO, UniProt/DIAMOND, Orthofinder, gene location (requires workspace/input/species_gff), and gene expression (requires workspace/input/species_expression).
run_wgd_ksd=0 # NOT WORKING YET! Obtain dS distribution in paralog pairs to infer WGD events. https://github.com/arzwa/wgd
# Genome analysis
run_busco_genome=0 # Gene set completeness analysis. https://gitlab.com/ezlab/busco
run_subphaser=0 # Infer polyploid subgenome structures. https://github.com/zhangrengang/SubPhaser
run_genome_fx2tab=0 # Compute sequence stats of the reference genome. https://bioinf.shenwei.me/seqkit/usage/#fx2tab-tab2fx
run_scaffold_histogram=0 # Generate a scaffold size histogram.
run_genome_mmseqs2taxonomy=0 # Taxonomic assignment of genome assembly. Estimated memory consumption: 120G. https://github.com/timkahlke/BASTA
run_genome_contamination_removal=0 # Removal of taxonomically non-compatible contaminated seqeunces using the result of MMseqs2 Taxonomy.
run_jcvi_dotplot=0 # JCVI's MCscan pipeline to generate self-self dotplot for synteny https://github.com/tanghaibao/jcvi
# DNA-seq analysis
run_genomescope=0 # Estimate genome size and heterozygosity using DNA-seq k-mer frequency https://github.com/tbenavi1/genomescope2.0
# Summary
run_multispecies_summary=1 # Generate plots and tables that summarize annotation info across species.

busco_lineage="embryophyta_odb12" # See here for available datasets: https://busco-data.ezlab.org/v5/data/lineages/
genetic_code=1
contamination_removal_rank="phylum" # {phylum, class, order, family, genus, species} The rank of taxonomic assignment to be used for contamination removal.

### End: Modify this block to tailor your analysis ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

dir_pg="/workspace"
dir_myscript="/script/script"
source "${dir_myscript}/gg_util.sh" # Load utility functions
gg_source_home_bashrc
gg_prepare_cmd_runtime "${dir_pg}" "base" 0 1

enable_all_run_flags_for_debug_mode

resolve_species_file() {
  local search_dir=$1
  local species_prefix=$2
  local label=$3
  local required=${4:-0}
  local matches=()
  if [[ ! -d "${search_dir}" ]]; then
    if [[ ${required} -eq 1 ]]; then
      echo "${label} directory not found: ${search_dir}. Exiting." >&2
      exit 1
    fi
    echo ""
    return
  fi
  mapfile -t matches < <(find "${search_dir}" -maxdepth 1 -type f -name "${species_prefix}*" | sort)
  if [[ ${#matches[@]} -eq 0 ]]; then
    if [[ ${required} -eq 1 ]]; then
      echo "No ${label} file matched prefix '${species_prefix}' in: ${search_dir}. Exiting." >&2
      exit 1
    fi
    echo ""
    return
  fi
  if [[ ${#matches[@]} -gt 1 ]]; then
    echo "Multiple ${label} files matched '${species_prefix}' in ${search_dir}. Using: ${matches[0]}" >&2
  fi
  echo "${matches[0]}"
}

dir_og="${dir_pg_output}/orthogroup"
dir_sp_cds="${dir_pg_input}/species_cds"
dir_sp_expression="${dir_pg_input}/species_expression"
dir_sp_gff="${dir_pg_input}/species_gff"
dir_sp_genome="${dir_pg_input}/species_genome"
dir_sp_dnaseq="${dir_pg_input}/species_dnaseq"
dir_sp_subphaser_cfg="${dir_pg_input}/species_genome_subphaser_cfg" # Directory for SubPhaser's species-wise CFGFILEs.
dir_mmseqs2_db="${dir_pg_output}/db_mmseqs2"

dir_sp_gff_info="${dir_pg_output}/species_cds_gff_info"
dir_sp_cds_busco_full="${dir_pg_output}/species_cds_busco_full"
dir_sp_cds_busco_short="${dir_pg_output}/species_cds_busco_short"
dir_sp_genome_busco_full="${dir_pg_output}/species_genome_busco_full"
dir_sp_genome_busco_short="${dir_pg_output}/species_genome_busco_short"
dir_sp_uniprot_annotation="${dir_pg_output}/species_cds_uniprot_annotation"
dir_sp_annotation="${dir_pg_output}/species_cds_annotation"
dir_sp_wgd_ksd="${dir_pg_output}/species_cds_wgd_ksd"
dir_sp_subphaser="${dir_pg_output}/species_genome_subphaser"
dir_sp_cds_fx2tab="${dir_pg_output}/species_cds_fx2tab"
dir_sp_genome_fx2tab="${dir_pg_output}/species_genome_fx2tab"
dir_sp_scaffold_histogram="${dir_pg_output}/species_genome_scaffold_histogram"
dir_sp_cds_mmseqs2taxonomy="${dir_pg_output}/species_cds_mmseqs2taxonomy"
dir_sp_cds_contamination_removal_fasta="${dir_pg_output}/species_cds_contamination_removal_fasta"
dir_sp_cds_contamination_removal_tsv="${dir_pg_output}/species_cds_contamination_removal_tsv"
dir_sp_genome_mmseqs2taxonomy="${dir_pg_output}/species_genome_mmseqs2taxonomy"
dir_sp_genome_contamination_removal_fasta="${dir_pg_output}/species_genome_contamination_removal_fasta"
dir_sp_genome_contamination_removal_tsv="${dir_pg_output}/species_genome_contamination_removal_tsv"
dir_sp_genomescope="${dir_pg_output}/species_dnaseq_genomescope"
dir_sp_jcvi_dotplot="${dir_pg_output}/species_jcvi_dotplot"
dir_multispecies_summary="${dir_pg_output}/annotation_summary"
dir_tmp="${dir_pg_output}/tmp"

infiles=()
mapfile -t infiles < <(find "${dir_sp_cds}" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fa.gz" -o -name "*.fas" -o -name "*.fas.gz" -o -name "*.fasta" -o -name "*.fasta.gz" -o -name "*.fna" -o -name "*.fna.gz" \) | sort)
if [[ ${#infiles[@]} -eq 0 ]]; then
  echo "No input fasta files were detected in: ${dir_sp_cds}. Exiting."
  exit 1
fi
task_index=$((SGE_TASK_ID-1))
if [[ ${task_index} -lt 0 || ${task_index} -ge ${#infiles[@]} ]]; then
  echo "SGE_TASK_ID=${SGE_TASK_ID} is out of range for ${#infiles[@]} input fasta files. Exiting."
  exit 1
fi
file_sp_cds="${infiles[${task_index}]}"
sp_ub=$(gg_species_name_from_path "${file_sp_cds}")
dir_sp_tmp="${dir_tmp}/${SGE_TASK_ID}_${sp_ub}"
echo "${#infiles[@]} input fasta files were detected in: ${dir_sp_cds}"
echo "Processing ${SGE_TASK_ID}th file: ${file_sp_cds}"
echo "Scientific name: ${sp_ub}"
echo "Working directory: ${dir_sp_tmp}"

file_orthogroup="${dir_pg_output}/orthofinder/Orthogroups/Orthogroups.tsv"
file_sp_expression="$(resolve_species_file "${dir_sp_expression}" "${sp_ub}" "expression")"
file_sp_gff="$(resolve_species_file "${dir_sp_gff}" "${sp_ub}" "gff")"
file_sp_genome="$(resolve_species_file "${dir_sp_genome}" "${sp_ub}" "genome")"
file_sp_subphaser_cfg="$(resolve_species_file "${dir_sp_subphaser_cfg}" "${sp_ub}" "subphaser cfg")"
if [[ -s ${file_sp_cds} ]]; then echo "CDS file found: ${file_sp_cds}"; else echo "CDS file not found."; fi
if [[ -s ${file_sp_expression} ]]; then echo "Expression file found: ${file_sp_expression}"; else echo "Expression file not found."; fi
if [[ -s ${file_sp_gff} ]]; then echo "GFF file found: ${file_sp_gff}"; else echo "GFF file not found."; fi
if [[ -s ${file_sp_genome} ]]; then echo "Genome fasta file found: ${file_sp_genome}"; else echo "Genome fasta file not found."; fi
if [[ -s ${file_sp_subphaser_cfg} ]]; then echo "SubPhaser's CFGFILE found: ${file_sp_subphaser_cfg}"; else echo "SubPhaser's CFGFILE not found."; fi

file_sp_gff_info="${dir_sp_gff_info}/${sp_ub}_gff_info.tsv"
file_sp_cds_busco_full="${dir_sp_cds_busco_full}/${sp_ub}_busco.full.tsv"
file_sp_cds_busco_short="${dir_sp_cds_busco_short}/${sp_ub}_busco.short.txt"
file_sp_genome_busco_full="${dir_sp_genome_busco_full}/${sp_ub}_busco.full.tsv"
file_sp_genome_busco_short="${dir_sp_genome_busco_short}/${sp_ub}_busco.short.txt"
file_sp_uniprot_annotation="${dir_sp_uniprot_annotation}/${sp_ub}_uniprot.tsv"
file_sp_annotation="${dir_sp_annotation}/${sp_ub}_annotation.tsv"
file_sp_wgd_ksd="${dir_sp_wgd_ksd}/${sp_ub}_wgd_ksd.tsv"
file_sp_subphaser="${dir_sp_subphaser}/${sp_ub}_subphaser.zip"
file_sp_cds_fx2tab="${dir_sp_cds_fx2tab}/${sp_ub}_fx2tab_cds.tsv"
file_sp_genome_fx2tab="${dir_sp_genome_fx2tab}/${sp_ub}_fx2tab_genome.tsv"
file_sp_scaffold_histogram="${dir_sp_scaffold_histogram}/${sp_ub}_scaffold_histogram.pdf"
file_sp_cds_mmseqs2taxonomy="${dir_sp_cds_mmseqs2taxonomy}/${sp_ub}_mmseqs2taxonomy.tsv"
file_sp_cds_contamination_removal_fasta="${dir_sp_cds_contamination_removal_fasta}/${sp_ub}_contamination_removal.fa.gz"
file_sp_cds_contamination_removal_tsv="${dir_sp_cds_contamination_removal_tsv}/${sp_ub}_contamination_removal.tsv"
file_sp_genome_mmseqs2taxonomy="${dir_sp_genome_mmseqs2taxonomy}/${sp_ub}_mmseqs2taxonomy.tsv"
file_sp_genome_contamination_removal_fasta="${dir_sp_genome_contamination_removal_fasta}/${sp_ub}_contamination_removal.fa.gz"
file_sp_genome_contamination_removal_tsv="${dir_sp_genome_contamination_removal_tsv}/${sp_ub}_contamination_removal.tsv"
file_sp_genomescope="${dir_sp_genomescope}/${sp_ub}_genomescope.zip"
file_sp_jcvi_dotplot="${dir_sp_jcvi_dotplot}/${sp_ub}_jcvi_dotplot.zip"
file_multispecies_summary="${dir_multispecies_summary}/annotation_summary.tsv"

ensure_dir "${dir_sp_tmp}"
cd "${dir_sp_tmp}"
ensure_dir "${dir_tmp}"
species_cds_validation_signature=$(
  {
    if stat --version >/dev/null 2>&1; then
      stat -c '%n:%s:%Y' "${infiles[@]}"
    else
      stat -f '%N:%z:%m' "${infiles[@]}"
    fi
  } | cksum | awk '{print $1}'
)
species_cds_validation_stamp="${dir_tmp}/species_cds_validation.${species_cds_validation_signature}.ok"
species_cds_validation_lock="${species_cds_validation_stamp}.lock"
if command -v flock >/dev/null 2>&1; then
  exec 9> "${species_cds_validation_lock}"
  flock 9
  if [[ ! -s ${species_cds_validation_stamp} ]]; then
    check_species_cds ${dir_pg}
    check_if_species_files_unique ${dir_sp_cds}
    touch "${species_cds_validation_stamp}"
  fi
  flock -u 9
  exec 9>&-
else
  species_cds_validation_lock_dir="${species_cds_validation_lock}.d"
  while ! mkdir "${species_cds_validation_lock_dir}" 2>/dev/null; do
    sleep 1
  done
  if [[ ! -s ${species_cds_validation_stamp} ]]; then
    check_species_cds ${dir_pg}
    check_if_species_files_unique ${dir_sp_cds}
    touch "${species_cds_validation_stamp}"
  fi
  rmdir "${species_cds_validation_lock_dir}" 2>/dev/null || true
fi

task="Gene trait extraction from gff files"
disable_if_no_input_file "run_get_gff_info" ${file_sp_gff}
if [[ ! -s ${file_sp_gff_info} && ${run_get_gff_info} -eq 1 ]]; then
  gg_step_start "${task}"
  if [[ -e gff2genestat.tsv ]]; then
    rm gff2genestat.tsv
  fi
  mkdir input_gff
  cp_out ${file_sp_gff} ./input_gff/

  python ${dir_myscript}/gff2genestat.py \
  --dir_gff ./input_gff \
  --feature 'CDS' \
  --multiple_hits 'longest' \
  --seqfile ${file_sp_cds} \
  --ncpu ${NSLOTS} \
  --outfile gff2genestat.tsv

  if [[ -s gff2genestat.tsv ]]; then
    mv_out gff2genestat.tsv ${file_sp_gff_info}
  fi
else
  gg_step_skip "${task}"
fi

task='BUSCO of species_cds'
disable_if_no_input_file "run_busco_cds" ${file_sp_cds}
if [[ ( ! -s ${file_sp_cds_busco_full} || ! -s ${file_sp_cds_busco_short} ) && ${run_busco_cds} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ -e "./busco_tmp" ]]; then
    rm -r "./busco_tmp"
  fi

  if [[ "${file_sp_cds}" == *gz ]]; then
    busco_infile=$(basename ${file_sp_cds} .gz)
    seqkit seq --threads ${NSLOTS} "${file_sp_cds}" --out-file "${busco_infile}"
  else
    busco_infile=${file_sp_cds}
  fi

  dir_busco_db=$(ensure_busco_download_path "${dir_pg}" "${busco_lineage}")
  if [[ $? -ne 0 ]]; then
    echo "Failed to prepare BUSCO dataset: ${busco_lineage}"
    exit 1
  fi
  dir_busco_lineage="${dir_busco_db}/lineages/${busco_lineage}"

  busco \
  --in "${busco_infile}" \
  --mode transcriptome \
  --out "busco_tmp" \
  --cpu ${NSLOTS} \
  --force \
  --evalue 1e-03 \
  --limit 20 \
  --lineage_dataset ${dir_busco_lineage} \
  --download_path ${dir_busco_db} \
  --offline
  busco_exit_status=$?

  if [[ ${busco_exit_status} -eq 0 ]]; then
    cp_out ./busco_tmp/run_${busco_lineage}/full_table.tsv ${file_sp_cds_busco_full}
    cp_out ./busco_tmp/run_${busco_lineage}/short_summary.txt ${file_sp_cds_busco_short}
    rm -r './busco_tmp'
  fi
  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task='BUSCO of species_genome'
disable_if_no_input_file "run_busco_genome" ${file_sp_genome}
if [[ ( ! -s ${file_sp_genome_busco_full} || ! -s ${file_sp_genome_busco_short} ) && ${run_busco_genome} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ -e "./busco_tmp" ]]; then
    rm -r "./busco_tmp"
  fi
  seqkit seq --threads ${NSLOTS} ${file_sp_genome} > "busco_genome_input.fa"
  dir_busco_db=$(ensure_busco_download_path "${dir_pg}" "${busco_lineage}")
  if [[ $? -ne 0 ]]; then
    echo "Failed to prepare BUSCO dataset: ${busco_lineage}"
    exit 1
  fi
  dir_busco_lineage="${dir_busco_db}/lineages/${busco_lineage}"

  busco \
  --in "busco_genome_input.fa" \
  --mode "genome" \
  --out "busco_tmp" \
  --cpu ${NSLOTS} \
  --force \
  --evalue 1e-03 \
  --limit 20 \
  --lineage_dataset ${dir_busco_lineage} \
  --download_path ${dir_busco_db} \
  --offline
  busco_exit_status=$?

  if [[ ${busco_exit_status} -eq 0 ]]; then
    cp_out ./busco_tmp/run_${busco_lineage}/full_table.tsv ${file_sp_genome_busco_full}
    cp_out ./busco_tmp/run_${busco_lineage}/short_summary.txt ${file_sp_genome_busco_short}
    rm -r './busco_tmp'
    rm "busco_genome_input.fa"
  fi
  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task="UniProt annotation by DIAMOND"
disable_if_no_input_file "run_uniprot_annotation" ${file_sp_cds}
if [[ ! -s ${file_sp_uniprot_annotation} ]] && [[ ${run_uniprot_annotation} -eq 1 ]]; then
  gg_step_start "${task}"
  uniprot_db_prefix=$(ensure_uniprot_sprot_db "${dir_pg}")
  if [[ $? -ne 0 ]]; then
    echo "Failed to prepare UniProt Swiss-Prot DB. Exiting."
    exit 1
  fi

  seqkit seq --remove-gaps --only-id --threads ${NSLOTS} ${file_sp_cds} \
  | cdskit pad --codontable ${genetic_code} \
  | cdskit mask --codontable ${genetic_code} --stopcodon yes --ambiguouscodon yes --maskchar 'N' \
  | seqkit translate --transl-table ${genetic_code} --threads ${NSLOTS} \
  > uniprot.query.pep.fas

  diamond blastp \
  --query uniprot.query.pep.fas \
  --threads ${NSLOTS} \
  --db ${uniprot_db_prefix} \
  --very-sensitive \
  --out uniprot.diamond.tsv \
  --outfmt 6 qseqid sseqid pident length evalue bitscore qlen \
  --max-target-seqs 1 \
  --evalue 1e-2

  python ${dir_myscript}/reformat_uniprot_diamond.py \
  --diamond_tsv uniprot.diamond.tsv \
  --query_fasta uniprot.query.pep.fas \
  --uniprot_fasta ${uniprot_db_prefix}.pep \
  --outfile uniprot.annotation.tsv

  cp_out uniprot.annotation.tsv ${file_sp_uniprot_annotation}
else
  gg_step_skip "${task}"
fi

task="seqkit fx2tab for the CDS sequences"
disable_if_no_input_file "run_cds_fx2tab" ${file_sp_cds}
if [[ ! -s ${file_sp_cds_fx2tab} && ${run_cds_fx2tab} -eq 1 ]]; then
  gg_step_start "${task}"

  seqkit fx2tab \
  --threads ${NSLOTS} \
  --length \
  --name \
  --gc \
  --gc-skew \
  --header-line \
  --only-id \
  ${file_sp_cds} \
  > "tmp.cds_length.tsv"

  if [[ -s "tmp.cds_length.tsv" ]]; then
    echo "Output file detected for the task: ${task}"
    mv_out "tmp.cds_length.tsv" ${file_sp_cds_fx2tab}
  fi
else
  gg_step_skip "${task}"
fi

task="MMseqs2 Taxonomy of the CDS sequences"
disable_if_no_input_file "run_cds_mmseqs2taxonomy" ${file_sp_cds}
if [[ ! -s ${file_sp_cds_mmseqs2taxonomy} && ${run_cds_mmseqs2taxonomy} -eq 1  && ${gg_debug_mode:-0} -eq 0 ]]; then
  gg_step_start "${task}"

  ensure_mmseqs_uniref90_db "${dir_mmseqs2_db}" "${NSLOTS}"
  if [[ $? -ne 0 ]]; then
    echo "Failed to prepare MMseqs2 UniRef90 DB. Exiting."
    exit 1
  fi

  if [[ ! -s "queryDB" ]]; then
    mmseqs createdb ${file_sp_cds} queryDB
  fi

  if [[ ! -e "tmp_mmseqs2" ]]; then
    mkdir -p "tmp_mmseqs2"
  fi

  mmseqs taxonomy "queryDB" "${dir_mmseqs2_db}/UniRef90_DB" "output_prefix" "tmp" \
  --split-mode 2 \
  --split-memory-limit $((${MEM_PER_HOST}*3/4))G \
  --majority 0.5 \
  --lca-mode 3 \
  --vote-mode 1 \
  --tax-lineage 2 \
  --orf-filter 1 \
  --threads ${NSLOTS}

  mmseqs createtsv "queryDB" "output_prefix" "result.tsv" --threads ${NSLOTS}

  if [[ -s "result.tsv" ]]; then
    echo "Output file detected for the task: ${task}"
    mv_out "result.tsv" ${file_sp_cds_mmseqs2taxonomy}
    rm queryDB*
    rm output_prefix.*
    rm -r "tmp_mmseqs2"
  fi
else
  gg_step_skip "${task}"
fi

task="Contaminated sequence removal from the CDS sequences"
disable_if_no_input_file "run_cds_contamination_removal" ${file_sp_cds} ${file_sp_cds_fx2tab} ${file_sp_cds_mmseqs2taxonomy}
if [[ ( ! -s ${file_sp_cds_contamination_removal_fasta} || ! -s ${file_sp_cds_contamination_removal_tsv} ) && ${run_cds_contamination_removal} -eq 1 ]]; then
  gg_step_start "${task}"

  if ! ensure_ete_taxonomy_db "${dir_pg}"; then
    echo "Failed to prepare ETE taxonomy DB. Exiting."
    exit 1
  fi

  python ${dir_myscript}/remove_contaminated_sequences.py \
  --fasta_file ${file_sp_cds} \
  --mmseqs2taxonomy_tsv ${file_sp_cds_mmseqs2taxonomy} \
  --fx2tab_tsv ${file_sp_cds_fx2tab} \
  --species_name ${sp_ub} \
  --rank "${contamination_removal_rank}" \
  --ncpu ${NSLOTS} \
  --rename_seq "no" \
  --verbose "no"

  if [[ -s "clean_sequences.fa" && -s "lineage_compatibility.tsv" ]]; then
    echo "Output file detected for the task: ${task}"
    ensure_parent_dir "${file_sp_cds_contamination_removal_fasta}"
    seqkit seq --threads ${NSLOTS} "clean_sequences.fa" --out-file "tmp.cds.clean.fa.gz"
    mv_out "tmp.cds.clean.fa.gz" ${file_sp_cds_contamination_removal_fasta}
    rm "clean_sequences.fa"
    mv_out "lineage_compatibility.tsv" ${file_sp_cds_contamination_removal_tsv}
  fi
else
  gg_step_skip "${task}"
fi

task="Merge annotations"
disable_if_no_input_file "run_annotation" ${file_sp_cds}
if [[ ! -s ${file_sp_annotation} ]] && [[ ${run_annotation} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ ! -n ${file_sp_expression} ]]; then
    file_sp_expression="placeholder_path"
  fi
  if [[ ! -n ${file_sp_gff} ]]; then
    file_sp_gff="placeholder_path"
  fi

  ensure_parent_dir "${file_sp_annotation}"
  python ${dir_myscript}/merge_cds_annotation.py \
  --scientific_name "${sp_ub}" \
  --cds_fasta ${file_sp_cds} \
  --uniprot_tsv ${file_sp_uniprot_annotation} \
  --busco_tsv ${file_sp_cds_busco_full} \
  --expression_tsv ${file_sp_expression} \
  --gff_info ${file_sp_gff_info} \
  --orthogroup_tsv ${file_orthogroup} \
  --mmseqs2taxonomy_tsv ${file_sp_cds_mmseqs2taxonomy} \
  --fx2tab ${file_sp_cds_fx2tab} \
  --ncpu ${NSLOTS} \
  --out_tsv "tmp.annotation.tsv"
  if [[ -s "tmp.annotation.tsv" ]]; then
    mv_out "tmp.annotation.tsv" ${file_sp_annotation}
  fi
else
  gg_step_skip "${task}"
fi

task="wgd ksd"
disable_if_no_input_file "run_wgd_ksd" ${file_sp_cds}
if [[ ! -s ${file_sp_wgd_ksd} && ${run_wgd_ksd} -eq 1 && ${gg_debug_mode:-0} -eq 0 ]]; then
  gg_step_start "${task}"

  seqkit seq --remove-gaps --upper-case --threads ${NSLOTS} ${file_sp_cds} \
  | gg_prepare_cds_fasta_stream "${NSLOTS}" "${genetic_code}" \
  | cdskit mask --codontable ${genetic_code} --stopcodon no --ambiguouscodon yes --maskchar 'N' \
  > "tmp.${sp_ub}.nuc.fasta"

  wgd dmd --ignorestop  --nostrictcds "tmp.${sp_ub}.nuc.fasta"
  wgd --verbosity debug ksd \
  --aligner 'mafft' \
  --wm 'alc' \
  --n_threads ${NSLOTS} \
  --max_pairwise 1000 \
  "./wgd_dmd/tmp.${sp_ub}.nuc.fasta.mcl" \
  "tmp.${sp_ub}.nuc.fasta"

else
  gg_step_skip "${task}"
fi

task="SubPhaser"
disable_if_no_input_file "run_subphaser" ${file_sp_genome} ${file_sp_subphaser_cfg}
if [[ ! -s ${file_sp_subphaser} && ${run_subphaser} -eq 1 ]]; then
  gg_step_start "${task}"

  subphaser \
  -genomes ${file_sp_genome} \
  -sg_cfgs ${file_sp_subphaser_cfg} \
  -no_label \
  -prefix "${sp_ub}." \
  -outdir "subphaser" \
  -tmpdir "tmp" \
  -ncpu ${NSLOTS} \
  -max_memory ${MEM_PER_HOST}
  subphaser_exit_status=$?
  echo "Subphaser exit status = ${subphaser_exit_status}"

  if [[ ${subphaser_exit_status} -eq 0 ]]; then
    echo "Zipping and copying SubPhaser's output files."
    zip -rq "${sp_ub}.subphaser.zip" "${sp_ub}.subphaser"
    cp_out "${sp_ub}.subphaser.zip" ${file_sp_subphaser}
  else
    echo "SubPhaser's exit status is not 0. Skipped zipping and copying SubPhaser's output files."
  fi
else
  gg_step_skip "${task}"
fi

task="seqkit fx2tab for the reference genome"
disable_if_no_input_file "run_genome_fx2tab" ${file_sp_genome}
if [[ ! -s ${file_sp_genome_fx2tab} && ${run_genome_fx2tab} -eq 1 ]]; then
  gg_step_start "${task}"

  seqkit fx2tab \
  --threads ${NSLOTS} \
  --length \
  --name \
  --gc \
  --gc-skew \
  --header-line \
  --only-id \
  ${file_sp_genome} \
  > "tmp.scaffold_length.tsv"

  if [[ -s "tmp.scaffold_length.tsv" ]]; then
    echo "Output file detected for the task: ${task}"
    mv_out "tmp.scaffold_length.tsv" ${file_sp_genome_fx2tab}
  fi
else
  gg_step_skip "${task}"
fi

task="scaffold size histogram"
disable_if_no_input_file "run_scaffold_histogram" ${file_sp_genome_fx2tab}
if [[ ! -s ${file_sp_scaffold_histogram} && ${run_scaffold_histogram} -eq 1 ]]; then
  gg_step_start "${task}"

  python ${dir_myscript}/scaffold_size_histogram.py \
  --min_scaffold_size 1000000 \
  --fx2tab ${file_sp_genome_fx2tab}

  if [[ -s "scaffold_size_histogram.pdf" ]]; then
    echo "Output file detected for the task: ${task}"
    mv_out "scaffold_size_histogram.pdf" ${file_sp_scaffold_histogram}
  fi
else
  gg_step_skip "${task}"
fi

task="MMseqs2 Taxonomy of the reference genome"
disable_if_no_input_file "run_genome_mmseqs2taxonomy" ${file_sp_genome}
if [[ ! -s ${file_sp_genome_mmseqs2taxonomy} && ${run_genome_mmseqs2taxonomy} -eq 1 ]]; then
  gg_step_start "${task}"

  ensure_mmseqs_uniref90_db "${dir_mmseqs2_db}" "${NSLOTS}"
  if [[ $? -ne 0 ]]; then
    echo "Failed to prepare MMseqs2 UniRef90 DB. Exiting."
    exit 1
  fi

  if [[ ! -s "queryDB" ]]; then
    mmseqs createdb ${file_sp_genome} queryDB
  fi

  if [[ ! -e "tmp_mmseqs2" ]]; then
    mkdir -p "tmp_mmseqs2"
  fi

  mmseqs taxonomy "queryDB" "${dir_mmseqs2_db}/UniRef90_DB" "output_prefix" "tmp" \
  --split-mode 2 \
  --split-memory-limit $((${MEM_PER_HOST}*3/4))G \
  --majority 0.5 \
  --lca-mode 3 \
  --vote-mode 1 \
  --tax-lineage 2 \
  --orf-filter 1 \
  --threads ${NSLOTS}

  mmseqs createtsv "queryDB" "output_prefix" "result.tsv" --threads ${NSLOTS}

  if [[ -s "result.tsv" ]]; then
    echo "Output file detected for the task: ${task}"
    mv_out "result.tsv" ${file_sp_genome_mmseqs2taxonomy}
    rm queryDB*
    rm output_prefix.*
    rm -r "tmp_mmseqs2"
  fi
else
  gg_step_skip "${task}"
fi

task="Contaminated sequence removal from the reference genome"
disable_if_no_input_file "run_genome_contamination_removal" ${file_sp_genome} ${file_sp_genome_mmseqs2taxonomy} ${file_sp_genome_fx2tab}
if [[ ( ! -s ${file_sp_genome_contamination_removal_fasta} || ! -s ${file_sp_genome_contamination_removal_tsv} ) && ${run_genome_contamination_removal} -eq 1 ]]; then
  gg_step_start "${task}"

  if ! ensure_ete_taxonomy_db "${dir_pg}"; then
    echo "Failed to prepare ETE taxonomy DB. Exiting."
    exit 1
  fi

  python ${dir_myscript}/remove_contaminated_sequences.py \
  --fasta_file ${file_sp_genome} \
  --mmseqs2taxonomy_tsv ${file_sp_genome_mmseqs2taxonomy} \
  --fx2tab_tsv ${file_sp_genome_fx2tab} \
  --species_name ${sp_ub} \
  --rank "phylum" \
  --ncpu ${NSLOTS} \
  --rename_seq "yes" \
  --rename_prefix "scaffold" \
  --verbose "no"

  if [[ -s "clean_sequences.fa" && -s "lineage_compatibility.tsv" ]]; then
    echo "Output file detected for the task: ${task}"
    ensure_parent_dir "${file_sp_genome_contamination_removal_fasta}"
    seqkit seq --threads ${NSLOTS} "clean_sequences.fa" --out-file "tmp.genome.clean.fa.gz"
    mv_out "tmp.genome.clean.fa.gz" ${file_sp_genome_contamination_removal_fasta}
    rm "clean_sequences.fa"
    mv_out "lineage_compatibility.tsv" ${file_sp_genome_contamination_removal_tsv}
  fi
else
  gg_step_skip "${task}"
fi

task="GenomeScope"
if [[ ! -e "${dir_sp_dnaseq}/${sp_ub}" ]]; then echo "dir_sp_dnaseq/sp not found. Skipping ${task}"; run_genomescope=0; fi
if [[ ! -s ${file_sp_genomescope} && ${run_genomescope} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ -e ${sp_ub}.genomescope ]]; then
    rm -r ${sp_ub}.genomescope
  fi
  mkdir ${sp_ub}.genomescope

  find ${dir_sp_dnaseq}/${sp_ub} | grep -e ".fq$" -e ".fastq$" -e ".fq.gz$" -e ".fastq.gz$" > input_fastq_files.txt
  klength=21
  kmer_lower=1
  kmer_upper=1000 # Could be automatically calculated with the file from species_genome
  mkdir ./tmp
  kmc -k${klength} -t${NSLOTS} -m${MEM_PER_HOST} -ci${kmer_lower} -cs${kmer_upper} @input_fastq_files.txt tmp.reads ./tmp
  kmc_tools transform tmp.reads histogram tmp.reads.histo -cx${kmer_upper}
  rm -r ./tmp
  /usr/local/bin/genomescope2.0/genomescope.R -i "tmp.reads.histo" -o "${sp_ub}.genomescope" -k ${klength}

  if [[ -s "${sp_ub}.genomescope/transformed_linear_plot.png" ]]; then
    echo "GenomeScope output file was detected. Start compressing."
    mv_out "tmp.reads.histo" "${sp_ub}.genomescope/kmc.histo.tsv"
    zip -rq "${sp_ub}.genomescope.zip" "${sp_ub}.genomescope"
    mv_out "${sp_ub}.genomescope.zip" "${file_sp_genomescope}"
  fi
else
  gg_step_skip "${task}"
fi

task="JCVI synteny dotplot"
disable_if_no_input_file "run_jcvi_dotplot" ${file_sp_cds} ${file_sp_gff} ${file_sp_genome_fx2tab}
if [[ ! -s ${file_sp_jcvi_dotplot} && ${run_jcvi_dotplot} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ -e ${sp_ub}.jcvi_dotplot ]]; then
    rm -r ${sp_ub}.jcvi_dotplot
  fi
  mkdir ${sp_ub}.jcvi_dotplot
  cd ${sp_ub}.jcvi_dotplot

  minimum_scaffold_size=1000000
  jcvi_cscore="0.7"
  gff_feature="mRNA"
  sp=${sp_ub//_/ }
  gff_attribute=$(python ${dir_myscript}/gff_attribute_identifier.py --gff_feature ${gff_feature} --fasta_file ${file_sp_cds} --gff_file ${file_sp_gff})
  selected_scaffolds=$(python -c "import sys,pandas;d=pandas.read_csv(sys.argv[1],sep='\t',header=0);print(','.join(d.loc[(d.loc[:,'length']>=int(sys.argv[2])),'#id']))" ${file_sp_genome_fx2tab} ${minimum_scaffold_size})
  echo "Detected GFF attribute that matches to the FASTA sequence names: ${gff_attribute}"
  echo "Selected sequences for the dot-plot visualization: ${selected_scaffolds}"

  seqkit seq ${file_sp_cds} | sed -e "s/^>${sp_ub}_/>/" > tmp.jcvi_input.cds.fasta
  python -m jcvi.formats.gff bed --primary_only --type=${gff_feature} --key=${gff_attribute} ${file_sp_gff} -o tmp.species1.bed
  python -m jcvi.formats.fasta format tmp.jcvi_input.cds.fasta species1.cds
  touch species1.bed
  selected_scaffold_list=()
  mapfile -t selected_scaffold_list < <(printf '%s' "${selected_scaffolds}" | tr ',' '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' -e '/^$/d')
  for scaffold in "${selected_scaffold_list[@]}"; do
      grep -e "^${scaffold}[[:space:]]" tmp.species1.bed >> species1.bed
  done
  cp_out species1.bed species2.bed
  cp_out species1.cds species2.cds

  python -m jcvi.compara.catalog ortholog species1 species2 --no_strip_names --no_dotplot --cscore=${jcvi_cscore}
  python -m jcvi.graphics.dotplot species1.species2.anchors --notex --figsize 7x7 --genomenames="${sp}_${sp}" --style='white' --nmax=100000
  python -m jcvi.compara.synteny depth --histogram species1.species2.anchors
  python -m jcvi.compara.synteny screen --minspan=30 --minsize=0 --simple species1.species2.anchors species1.species2.anchors.new

  echo """# y, xstart, xend, rotation, color, label, va,  bed
0.6,     0.15,    0.95,       0,      , species1, top, species1.bed
0.4,     0.15,    0.95,       0,      , species2, top, species2.bed
# edges
e, 0, 1, species1.species2.anchors.simple""" > layout

  echo """${selected_scaffolds}
  ${selected_scaffolds}""" > seqids
  python -m jcvi.graphics.karyotype --notex --figsize 7x7 seqids layout

  rm species2.*
  rm species1.bed species1.cds
  rm seqids layout
  rm species1.species2.last species1.species2.last.filtered
  rm tmp.scaffold_length.tsv
  rm tmp.species1.bed
  rm tmp.jcvi_input.cds.fasta
  mv_out karyotype.pdf ${sp_ub}.karyotype.pdf
  shopt -s nullglob
  for file in *; do
    mv_out "${file}" "${file/species1.species2/${sp_ub}}"
  done
  shopt -u nullglob

  cd ${dir_sp_tmp}

  if [[ -s "${sp_ub}.jcvi_dotplot/${sp_ub}.pdf" ]]; then
    echo "JCVI dotplot output file was detected. Start compressing."
    zip -rq "${sp_ub}.jcvi_dotplot.zip" "${sp_ub}.jcvi_dotplot"
    mv_out "${sp_ub}.jcvi_dotplot.zip" "${file_sp_jcvi_dotplot}"
  fi
else
  gg_step_skip "${task}"
fi

# Repeat-annotation steps were removed together with the dedicated
# `repeat` conda environment.

task="Multispecies annotation summary"
is_output_older_than_inputs "^file_sp_" ${file_multispecies_summary}; summary_flag=$?
if [[ ${run_multispecies_summary} -eq 1 && ${summary_flag} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_multispecies_summary}"
  cd ${dir_multispecies_summary}

  Rscript ${dir_myscript}/annotation_summary.r \
  --dir_species_tree="${dir_pg_output}/species_tree" \
  --dir_species_cds_busco="${dir_sp_cds_busco_full}" \
  --dir_species_genome_busco="${dir_sp_genome_busco_full}" \
  --dir_species_annotation="${dir_sp_annotation}" \
  --dir_species_cds_fx2tab="${dir_sp_cds_fx2tab}" \
  --dir_species_genome_fx2tab="${dir_sp_genome_fx2tab}" \
  --file_species_trait="${dir_pg_input}/species_trait/species_trait.tsv" \
  --file_orthogroup_gene_count="${dir_pg_output}/orthofinder/Orthogroups/Orthogroups.GeneCount.tsv" \
  --tree_annotation_dir="${dir_myscript}/tree_annotation" \
  --min_og_species='auto'

  if [[ -e "Rplots.pdf" ]]; then
    rm "Rplots.pdf"
  fi

  cd ${dir_sp_tmp}
else
  gg_step_skip "${task}"
fi

remove_empty_subdirs ${dir_pg_output}
if [[ ${delete_tmp_dir} -eq 1 ]]; then
    echo "Deleting tmp directory: ${dir_sp_tmp}"
    rm -rf ${dir_sp_tmp}
fi

echo "$(date): Exiting Singularity environment"

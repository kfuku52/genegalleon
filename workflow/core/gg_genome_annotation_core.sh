#!/usr/bin/env bash
set -euo pipefail

# Load shared defaults when available.
gg_core_self="${BASH_SOURCE[0]:-/script/core/gg_genome_annotation_core.sh}"
gg_common_params_file=""
for gg_common_params_candidate in \
  "$(cd "$(dirname "${gg_core_self}")" && pwd)/gg_common_params.sh" \
  "$(cd "$(dirname "${gg_core_self}")" && pwd)/../gg_common_params.sh" \
  "/script/gg_common_params.sh"
do
  if [[ -s "${gg_common_params_candidate}" ]]; then
    gg_common_params_file="${gg_common_params_candidate}"
    break
  fi
done
if [[ -n "${gg_common_params_file}" ]]; then
  # shellcheck disable=SC1090
  source "${gg_common_params_file}"
fi
unset gg_common_params_file gg_common_params_candidate gg_core_self

### Start: Job-supplied configuration ###
# Configuration variables are provided by gg_genome_annotation_entrypoint.sh.
busco_lineage="${busco_lineage:-${GG_COMMON_BUSCO_LINEAGE:-embryophyta_odb12}}"
genetic_code="${genetic_code:-${GG_COMMON_GENETIC_CODE:-1}}"
contamination_removal_rank="${contamination_removal_rank:-${GG_COMMON_CONTAMINATION_REMOVAL_RANK:-phylum}}"
### End: Job-supplied configuration ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

dir_pg="/workspace"
dir_script="/script/support"
source "${dir_script}/gg_util.sh" # Load utility functions
gg_source_home_bashrc
gg_prepare_cmd_runtime "${dir_pg}" "base" 0 1
delete_tmp_dir=${delete_tmp_dir:-1}

copy_busco_tables() {
  local busco_root_dir="$1"
  local lineage="$2"
  local file_full="$3"
  local file_short="$4"
  local run_dir="${busco_root_dir}/run_${lineage}"
  local full_src=""
  local short_src=""
  local -a short_candidates=()

  if [[ -s "${run_dir}/full_table.tsv" ]]; then
    full_src="${run_dir}/full_table.tsv"
  elif [[ -s "${run_dir}/full_table.tsv.gz" ]]; then
    gzip -cd "${run_dir}/full_table.tsv.gz" > "./tmp.busco.full_table.tsv"
    full_src="./tmp.busco.full_table.tsv"
  fi

  if [[ -s "${run_dir}/short_summary.txt" ]]; then
    short_src="${run_dir}/short_summary.txt"
  else
    mapfile -t short_candidates < <(find "${run_dir}" -maxdepth 1 -type f -name "short_summary*.txt" | sort)
    if [[ ${#short_candidates[@]} -gt 0 ]]; then
      short_src="${short_candidates[0]}"
    fi
  fi

  if [[ -z "${full_src}" || -z "${short_src}" ]]; then
    echo "BUSCO outputs were not found under ${run_dir}."
    echo "Expected full_table.tsv(.gz) and short_summary*.txt."
    rm -f -- "./tmp.busco.full_table.tsv"
    return 1
  fi

  cp_out "${full_src}" "${file_full}"
  cp_out "${short_src}" "${file_short}"
  rm -f -- "./tmp.busco.full_table.tsv"
}

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

dir_sp_cds="${dir_pg_input}/species_cds"
dir_sp_dnaseq="${dir_pg_input}/species_dnaseq"
dir_mmseqs2_db="${dir_pg_db}/mmseqs2"
dir_tmp="${dir_pg_output}/tmp"

if [[ ! -d "${dir_sp_cds}" ]]; then
  echo "Input directory not found: ${dir_sp_cds}. Exiting."
  exit 1
fi
infiles=()
mapfile -t infiles < <(find "${dir_sp_cds}" -maxdepth 1 -type f ! -name '.*' \( -name "*.fa" -o -name "*.fa.gz" -o -name "*.fas" -o -name "*.fas.gz" -o -name "*.fasta" -o -name "*.fasta.gz" -o -name "*.fna" -o -name "*.fna.gz" \) | sort)
if [[ ${#infiles[@]} -eq 0 ]]; then
  echo "No input fasta files were detected in: ${dir_sp_cds}. Exiting."
  exit 1
fi
if [[ ! "${SGE_TASK_ID}" =~ ^[0-9]+$ ]] || [[ ${SGE_TASK_ID} -lt 1 ]]; then
  echo "Invalid SGE_TASK_ID value (must be a positive integer): ${SGE_TASK_ID}"
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
file_sp_expression="$(resolve_species_file "${dir_pg_input}/species_expression" "${sp_ub}" "expression")"
file_sp_gff="$(resolve_species_file "${dir_pg_input}/species_gff" "${sp_ub}" "gff")"
file_sp_genome="$(resolve_species_file "${dir_pg_input}/species_genome" "${sp_ub}" "genome")"
file_sp_subphaser_cfg="$(resolve_species_file "${dir_pg_input}/species_genome_subphaser_cfg" "${sp_ub}" "subphaser cfg")"
if [[ -s "${file_sp_cds}" ]]; then echo "CDS file found: ${file_sp_cds}"; else echo "CDS file not found."; fi
if [[ -s "${file_sp_expression}" ]]; then echo "Expression file found: ${file_sp_expression}"; else echo "Expression file not found."; fi
if [[ -s "${file_sp_gff}" ]]; then echo "GFF file found: ${file_sp_gff}"; else echo "GFF file not found."; fi
if [[ -s "${file_sp_genome}" ]]; then echo "Genome fasta file found: ${file_sp_genome}"; else echo "Genome fasta file not found."; fi
if [[ -s "${file_sp_subphaser_cfg}" ]]; then echo "SubPhaser's CFGFILE found: ${file_sp_subphaser_cfg}"; else echo "SubPhaser's CFGFILE not found."; fi

file_sp_gff_info="${dir_pg_output}/species_gff_info/${sp_ub}_gff_info.tsv"
file_sp_cds_busco_full="${dir_pg_output}/species_cds_busco_full/${sp_ub}_busco.full.tsv"
file_sp_cds_busco_short="${dir_pg_output}/species_cds_busco_short/${sp_ub}_busco.short.txt"
file_sp_genome_busco_full="${dir_pg_output}/species_genome_busco_full/${sp_ub}_busco.full.tsv"
file_sp_genome_busco_short="${dir_pg_output}/species_genome_busco_short/${sp_ub}_busco.short.txt"
file_sp_uniprot_annotation="${dir_pg_output}/species_cds_uniprot_annotation/${sp_ub}_uniprot.tsv"
file_sp_annotation="${dir_pg_output}/species_cds_annotation/${sp_ub}_annotation.tsv"
file_sp_wgd_ksd="${dir_pg_output}/species_cds_wgd_ksd/${sp_ub}_wgd_ksd.tsv"
file_sp_subphaser="${dir_pg_output}/species_genome_subphaser/${sp_ub}_subphaser.zip"
file_sp_cds_fx2tab="${dir_pg_output}/species_cds_fx2tab/${sp_ub}_fx2tab_cds.tsv"
file_sp_genome_fx2tab="${dir_pg_output}/species_genome_fx2tab/${sp_ub}_fx2tab_genome.tsv"
file_sp_scaffold_histogram="${dir_pg_output}/species_genome_scaffold_histogram/${sp_ub}_scaffold_histogram.pdf"
file_sp_cds_mmseqs2taxonomy="${dir_pg_output}/species_cds_mmseqs2taxonomy/${sp_ub}_mmseqs2taxonomy.tsv"
file_sp_cds_contamination_removal_fasta="${dir_pg_output}/species_cds_contamination_removal_fasta/${sp_ub}_contamination_removal.fa.gz"
file_sp_cds_contamination_removal_tsv="${dir_pg_output}/species_cds_contamination_removal_tsv/${sp_ub}_contamination_removal.tsv"
file_sp_genome_mmseqs2taxonomy="${dir_pg_output}/species_genome_mmseqs2taxonomy/${sp_ub}_mmseqs2taxonomy.tsv"
file_sp_genome_contamination_removal_fasta="${dir_pg_output}/species_genome_contamination_removal_fasta/${sp_ub}_contamination_removal.fa.gz"
file_sp_genome_contamination_removal_tsv="${dir_pg_output}/species_genome_contamination_removal_tsv/${sp_ub}_contamination_removal.tsv"
file_sp_genomescope="${dir_pg_output}/species_dnaseq_genomescope/${sp_ub}_genomescope.zip"
file_sp_jcvi_dotplot="${dir_pg_output}/species_jcvi_dotplot/${sp_ub}_jcvi_dotplot.zip"
file_multispecies_summary="${dir_pg_output}/annotation_summary/annotation_summary.tsv"

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
  if [[ ! -s "${species_cds_validation_stamp}" ]]; then
    check_species_cds "${dir_pg}"
    check_if_species_files_unique "${dir_sp_cds}"
    touch "${species_cds_validation_stamp}"
  fi
  flock -u 9
  exec 9>&-
else
  species_cds_validation_lock_dir="${species_cds_validation_lock}.d"
  while ! mkdir "${species_cds_validation_lock_dir}" 2>/dev/null; do
    sleep 1
  done
  if [[ ! -s "${species_cds_validation_stamp}" ]]; then
    check_species_cds "${dir_pg}"
    check_if_species_files_unique "${dir_sp_cds}"
    touch "${species_cds_validation_stamp}"
  fi
  rmdir "${species_cds_validation_lock_dir}" 2>/dev/null || true
fi

task="Gene trait extraction from gff files"
disable_if_no_input_file "run_get_gff_info" "${file_sp_gff}"
if [[ ! -s "${file_sp_gff_info}" && ${run_get_gff_info} -eq 1 ]]; then
  gg_step_start "${task}"
  if [[ -e gff2genestat.tsv ]]; then
    rm -f -- gff2genestat.tsv
  fi
  if [[ -e input_gff ]]; then
    rm -rf -- input_gff
  fi
  mkdir -p input_gff
  cp_out "${file_sp_gff}" ./input_gff/

  python "${dir_script}/gff2genestat.py" \
  --dir_gff ./input_gff \
  --feature 'CDS' \
  --multiple_hits 'longest' \
  --seqfile "${file_sp_cds}" \
  --ncpu "${NSLOTS}" \
  --outfile gff2genestat.tsv

  if [[ -s gff2genestat.tsv ]]; then
    mv_out gff2genestat.tsv "${file_sp_gff_info}"
  fi
else
  gg_step_skip "${task}"
fi

task='BUSCO of species_cds'
disable_if_no_input_file "run_busco_cds" "${file_sp_cds}"
if [[ ( ! -s "${file_sp_cds_busco_full}" || ! -s "${file_sp_cds_busco_short}" ) && ${run_busco_cds} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ -e "./busco_tmp" ]]; then
    rm -rf -- "./busco_tmp"
  fi

  if [[ "${file_sp_cds}" == *gz ]]; then
    busco_infile=$(basename "${file_sp_cds}" .gz)
    seqkit seq --threads "${NSLOTS}" "${file_sp_cds}" --out-file "${busco_infile}"
  else
    busco_infile=${file_sp_cds}
  fi

  if ! dir_busco_db=$(ensure_busco_download_path "${dir_pg}" "${busco_lineage}"); then
    echo "Failed to prepare BUSCO dataset: ${busco_lineage}"
    exit 1
  fi
  dir_busco_lineage="${dir_busco_db}/lineages/${busco_lineage}"

  if busco \
    --in "${busco_infile}" \
    --mode transcriptome \
    --out "busco_tmp" \
    --cpu "${NSLOTS}" \
    --force \
    --evalue 1e-03 \
    --limit 20 \
    --lineage_dataset "${dir_busco_lineage}" \
    --download_path "${dir_busco_db}" \
    --offline; then
    if copy_busco_tables "./busco_tmp" "${busco_lineage}" "${file_sp_cds_busco_full}" "${file_sp_cds_busco_short}"; then
      rm -rf -- "./busco_tmp"
    else
      echo "Failed to locate normalized BUSCO outputs for CDS BUSCO. Exiting."
      exit 1
    fi
  else
    busco_exit_status=$?
    echo "BUSCO of species_cds failed with exit code ${busco_exit_status}. Continuing without BUSCO outputs."
  fi
  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task='BUSCO of species_genome'
disable_if_no_input_file "run_busco_genome" "${file_sp_genome}"
if [[ ( ! -s "${file_sp_genome_busco_full}" || ! -s "${file_sp_genome_busco_short}" ) && ${run_busco_genome} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ -e "./busco_tmp" ]]; then
    rm -rf -- "./busco_tmp"
  fi
  seqkit seq --threads "${NSLOTS}" "${file_sp_genome}" > "busco_genome_input.fa"
  if ! dir_busco_db=$(ensure_busco_download_path "${dir_pg}" "${busco_lineage}"); then
    echo "Failed to prepare BUSCO dataset: ${busco_lineage}"
    exit 1
  fi
  dir_busco_lineage="${dir_busco_db}/lineages/${busco_lineage}"

  if busco \
    --in "busco_genome_input.fa" \
    --mode "genome" \
    --out "busco_tmp" \
    --cpu "${NSLOTS}" \
    --force \
    --evalue 1e-03 \
    --limit 20 \
    --lineage_dataset "${dir_busco_lineage}" \
    --download_path "${dir_busco_db}" \
    --offline; then
    if copy_busco_tables "./busco_tmp" "${busco_lineage}" "${file_sp_genome_busco_full}" "${file_sp_genome_busco_short}"; then
      rm -rf -- "./busco_tmp"
      rm -f -- "busco_genome_input.fa"
    else
      echo "Failed to locate normalized BUSCO outputs for genome BUSCO. Exiting."
      exit 1
    fi
  else
    busco_exit_status=$?
    echo "BUSCO of species_genome failed with exit code ${busco_exit_status}. Continuing without BUSCO outputs."
  fi
  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task="UniProt annotation by DIAMOND"
disable_if_no_input_file "run_uniprot_annotation" "${file_sp_cds}"
if [[ ! -s "${file_sp_uniprot_annotation}" ]] && [[ ${run_uniprot_annotation} -eq 1 ]]; then
  gg_step_start "${task}"
  if ! uniprot_db_prefix=$(ensure_uniprot_sprot_db "${dir_pg}"); then
    echo "Failed to prepare UniProt Swiss-Prot DB. Exiting."
    exit 1
  fi

  seqkit seq --remove-gaps --only-id --threads "${NSLOTS}" "${file_sp_cds}" \
  | cdskit pad --codontable "${genetic_code}" \
  | cdskit mask --codontable "${genetic_code}" --stopcodon yes --ambiguouscodon yes --maskchar 'N' \
  | seqkit translate --transl-table "${genetic_code}" --threads "${NSLOTS}" \
  > uniprot.query.pep.fas

  diamond blastp \
  --query uniprot.query.pep.fas \
  --threads "${NSLOTS}" \
  --db "${uniprot_db_prefix}" \
  --very-sensitive \
  --out uniprot.diamond.tsv \
  --outfmt 6 qseqid sseqid pident length evalue bitscore qlen \
  --max-target-seqs 1 \
  --evalue 1e-2

  python "${dir_script}/reformat_uniprot_diamond.py" \
  --diamond_tsv uniprot.diamond.tsv \
  --query_fasta uniprot.query.pep.fas \
  --uniprot_fasta "${uniprot_db_prefix}.pep" \
  --outfile uniprot.annotation.tsv

  cp_out uniprot.annotation.tsv "${file_sp_uniprot_annotation}"
else
  gg_step_skip "${task}"
fi

task="seqkit fx2tab for the CDS sequences"
disable_if_no_input_file "run_cds_fx2tab" "${file_sp_cds}"
if [[ ! -s "${file_sp_cds_fx2tab}" && ${run_cds_fx2tab} -eq 1 ]]; then
  gg_step_start "${task}"

  seqkit fx2tab \
  --threads "${NSLOTS}" \
  --length \
  --name \
  --gc \
  --gc-skew \
  --header-line \
  --only-id \
  "${file_sp_cds}" \
  > "tmp.cds_length.tsv"

  if [[ -s "tmp.cds_length.tsv" ]]; then
    echo "Output file detected for the task: ${task}"
    mv_out "tmp.cds_length.tsv" "${file_sp_cds_fx2tab}"
  fi
else
  gg_step_skip "${task}"
fi

task="MMseqs2 Taxonomy of the CDS sequences"
disable_if_no_input_file "run_cds_mmseqs2taxonomy" "${file_sp_cds}"
if [[ ! -s "${file_sp_cds_mmseqs2taxonomy}" && ${run_cds_mmseqs2taxonomy} -eq 1  && ${gg_debug_mode:-0} -eq 0 ]]; then
  gg_step_start "${task}"

  if ! ensure_mmseqs_uniref90_db "${dir_mmseqs2_db}" "${NSLOTS}"; then
    echo "Failed to prepare MMseqs2 UniRef90 DB. Exiting."
    exit 1
  fi

  if [[ ! -s "queryDB" ]]; then
    mmseqs createdb "${file_sp_cds}" queryDB
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
  --threads "${NSLOTS}"

  mmseqs createtsv "queryDB" "output_prefix" "result.tsv" --threads "${NSLOTS}"

  if [[ -s "result.tsv" ]]; then
    echo "Output file detected for the task: ${task}"
    mv_out "result.tsv" "${file_sp_cds_mmseqs2taxonomy}"
    rm -f -- queryDB*
    rm -f -- output_prefix.*
    rm -rf -- "tmp_mmseqs2"
  fi
else
  gg_step_skip "${task}"
fi

task="Contaminated sequence removal from the CDS sequences"
disable_if_no_input_file "run_cds_contamination_removal" "${file_sp_cds}" "${file_sp_cds_fx2tab}" "${file_sp_cds_mmseqs2taxonomy}"
if [[ ( ! -s "${file_sp_cds_contamination_removal_fasta}" || ! -s "${file_sp_cds_contamination_removal_tsv}" ) && ${run_cds_contamination_removal} -eq 1 ]]; then
  gg_step_start "${task}"

  if ! ensure_ete_taxonomy_db "${dir_pg}"; then
    echo "Failed to prepare ETE taxonomy DB. Exiting."
    exit 1
  fi

  python "${dir_script}/remove_contaminated_sequences.py" \
  --fasta_file "${file_sp_cds}" \
  --mmseqs2taxonomy_tsv "${file_sp_cds_mmseqs2taxonomy}" \
  --fx2tab_tsv "${file_sp_cds_fx2tab}" \
  --species_name "${sp_ub}" \
  --rank "${contamination_removal_rank}" \
  --ncpu "${NSLOTS}" \
  --rename_seq "no" \
  --verbose "no"

  if [[ -s "clean_sequences.fa" && -s "lineage_compatibility.tsv" ]]; then
    echo "Output file detected for the task: ${task}"
    ensure_parent_dir "${file_sp_cds_contamination_removal_fasta}"
    seqkit seq --threads "${NSLOTS}" "clean_sequences.fa" --out-file "tmp.cds.clean.fa.gz"
    mv_out "tmp.cds.clean.fa.gz" "${file_sp_cds_contamination_removal_fasta}"
    rm -f -- "clean_sequences.fa"
    mv_out "lineage_compatibility.tsv" "${file_sp_cds_contamination_removal_tsv}"
  fi
else
  gg_step_skip "${task}"
fi

task="Merge annotations"
disable_if_no_input_file "run_annotation" "${file_sp_cds}"
if [[ ! -s "${file_sp_annotation}" ]] && [[ ${run_annotation} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ -z "${file_sp_expression}" ]]; then
    file_sp_expression="placeholder_path"
  fi
  if [[ -z "${file_sp_gff}" ]]; then
    file_sp_gff="placeholder_path"
  fi

  ensure_parent_dir "${file_sp_annotation}"
  python "${dir_script}/merge_cds_annotation.py" \
  --scientific_name "${sp_ub}" \
  --cds_fasta "${file_sp_cds}" \
  --uniprot_tsv "${file_sp_uniprot_annotation}" \
  --busco_tsv "${file_sp_cds_busco_full}" \
  --expression_tsv "${file_sp_expression}" \
  --gff_info "${file_sp_gff_info}" \
  --orthogroup_tsv "${file_orthogroup}" \
  --mmseqs2taxonomy_tsv "${file_sp_cds_mmseqs2taxonomy}" \
  --fx2tab "${file_sp_cds_fx2tab}" \
  --ncpu "${NSLOTS}" \
  --out_tsv "tmp.annotation.tsv"
  if [[ -s "tmp.annotation.tsv" ]]; then
    mv_out "tmp.annotation.tsv" "${file_sp_annotation}"
  fi
else
  gg_step_skip "${task}"
fi

task="wgd ksd"
disable_if_no_input_file "run_wgd_ksd" "${file_sp_cds}"
if [[ ! -s "${file_sp_wgd_ksd}" && ${run_wgd_ksd} -eq 1 && ${gg_debug_mode:-0} -eq 0 ]]; then
  gg_step_start "${task}"

  seqkit seq --remove-gaps --upper-case --threads "${NSLOTS}" "${file_sp_cds}" \
  | gg_prepare_cds_fasta_stream "${NSLOTS}" "${genetic_code}" \
  | cdskit mask --codontable "${genetic_code}" --stopcodon no --ambiguouscodon yes --maskchar 'N' \
  > "tmp.${sp_ub}.nuc.fasta"

  wgd dmd --ignorestop  --nostrictcds "tmp.${sp_ub}.nuc.fasta"
  wgd --verbosity debug ksd \
  --aligner 'mafft' \
  --wm 'alc' \
  --n_threads "${NSLOTS}" \
  --max_pairwise 1000 \
  "./wgd_dmd/tmp.${sp_ub}.nuc.fasta.mcl" \
  "tmp.${sp_ub}.nuc.fasta"

else
  gg_step_skip "${task}"
fi

task="SubPhaser"
disable_if_no_input_file "run_subphaser" "${file_sp_genome}" "${file_sp_subphaser_cfg}"
if [[ ! -s "${file_sp_subphaser}" && ${run_subphaser} -eq 1 ]]; then
  gg_step_start "${task}"

  if subphaser \
    -genomes "${file_sp_genome}" \
    -sg_cfgs "${file_sp_subphaser_cfg}" \
    -no_label \
    -prefix "${sp_ub}." \
    -outdir "subphaser" \
    -tmpdir "tmp" \
    -ncpu "${NSLOTS}" \
    -max_memory "${MEM_PER_HOST}"; then
    subphaser_exit_status=0
  else
    subphaser_exit_status=$?
  fi
  echo "Subphaser exit status = ${subphaser_exit_status}"

  if [[ ${subphaser_exit_status} -eq 0 ]]; then
    echo "Zipping and copying SubPhaser's output files."
    zip -rq "${sp_ub}.subphaser.zip" "${sp_ub}.subphaser"
    cp_out "${sp_ub}.subphaser.zip" "${file_sp_subphaser}"
  else
    echo "SubPhaser's exit status is not 0. Skipped zipping and copying SubPhaser's output files."
  fi
else
  gg_step_skip "${task}"
fi

task="seqkit fx2tab for the reference genome"
disable_if_no_input_file "run_genome_fx2tab" "${file_sp_genome}"
if [[ ! -s "${file_sp_genome_fx2tab}" && ${run_genome_fx2tab} -eq 1 ]]; then
  gg_step_start "${task}"

  seqkit fx2tab \
  --threads "${NSLOTS}" \
  --length \
  --name \
  --gc \
  --gc-skew \
  --header-line \
  --only-id \
  "${file_sp_genome}" \
  > "tmp.scaffold_length.tsv"

  if [[ -s "tmp.scaffold_length.tsv" ]]; then
    echo "Output file detected for the task: ${task}"
    mv_out "tmp.scaffold_length.tsv" "${file_sp_genome_fx2tab}"
  fi
else
  gg_step_skip "${task}"
fi

task="scaffold size histogram"
disable_if_no_input_file "run_scaffold_histogram" "${file_sp_genome_fx2tab}"
if [[ ! -s "${file_sp_scaffold_histogram}" && ${run_scaffold_histogram} -eq 1 ]]; then
  gg_step_start "${task}"

  python "${dir_script}/scaffold_size_histogram.py" \
  --min_scaffold_size 1000000 \
  --fx2tab "${file_sp_genome_fx2tab}"

  if [[ -s "scaffold_size_histogram.pdf" ]]; then
    echo "Output file detected for the task: ${task}"
    mv_out "scaffold_size_histogram.pdf" "${file_sp_scaffold_histogram}"
  fi
else
  gg_step_skip "${task}"
fi

task="MMseqs2 Taxonomy of the reference genome"
disable_if_no_input_file "run_genome_mmseqs2taxonomy" "${file_sp_genome}"
if [[ ! -s "${file_sp_genome_mmseqs2taxonomy}" && ${run_genome_mmseqs2taxonomy} -eq 1 ]]; then
  gg_step_start "${task}"

  if ! ensure_mmseqs_uniref90_db "${dir_mmseqs2_db}" "${NSLOTS}"; then
    echo "Failed to prepare MMseqs2 UniRef90 DB. Exiting."
    exit 1
  fi

  if [[ ! -s "queryDB" ]]; then
    mmseqs createdb "${file_sp_genome}" queryDB
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
  --threads "${NSLOTS}"

  mmseqs createtsv "queryDB" "output_prefix" "result.tsv" --threads "${NSLOTS}"

  if [[ -s "result.tsv" ]]; then
    echo "Output file detected for the task: ${task}"
    mv_out "result.tsv" "${file_sp_genome_mmseqs2taxonomy}"
    rm -f -- queryDB*
    rm -f -- output_prefix.*
    rm -rf -- "tmp_mmseqs2"
  fi
else
  gg_step_skip "${task}"
fi

task="Contaminated sequence removal from the reference genome"
disable_if_no_input_file "run_genome_contamination_removal" "${file_sp_genome}" "${file_sp_genome_mmseqs2taxonomy}" "${file_sp_genome_fx2tab}"
if [[ ( ! -s "${file_sp_genome_contamination_removal_fasta}" || ! -s "${file_sp_genome_contamination_removal_tsv}" ) && ${run_genome_contamination_removal} -eq 1 ]]; then
  gg_step_start "${task}"

  if ! ensure_ete_taxonomy_db "${dir_pg}"; then
    echo "Failed to prepare ETE taxonomy DB. Exiting."
    exit 1
  fi

  python "${dir_script}/remove_contaminated_sequences.py" \
  --fasta_file "${file_sp_genome}" \
  --mmseqs2taxonomy_tsv "${file_sp_genome_mmseqs2taxonomy}" \
  --fx2tab_tsv "${file_sp_genome_fx2tab}" \
  --species_name "${sp_ub}" \
  --rank "phylum" \
  --ncpu "${NSLOTS}" \
  --rename_seq "yes" \
  --rename_prefix "scaffold" \
  --verbose "no"

  if [[ -s "clean_sequences.fa" && -s "lineage_compatibility.tsv" ]]; then
    echo "Output file detected for the task: ${task}"
    ensure_parent_dir "${file_sp_genome_contamination_removal_fasta}"
    seqkit seq --threads "${NSLOTS}" "clean_sequences.fa" --out-file "tmp.genome.clean.fa.gz"
    mv_out "tmp.genome.clean.fa.gz" "${file_sp_genome_contamination_removal_fasta}"
    rm -f -- "clean_sequences.fa"
    mv_out "lineage_compatibility.tsv" "${file_sp_genome_contamination_removal_tsv}"
  fi
else
  gg_step_skip "${task}"
fi

task="GenomeScope"
if [[ ! -e "${dir_sp_dnaseq}/${sp_ub}" ]]; then echo "dir_sp_dnaseq/sp not found. Skipping ${task}"; run_genomescope=0; fi
if [[ ! -s "${file_sp_genomescope}" && ${run_genomescope} -eq 1 ]]; then
  gg_step_start "${task}"

	if [[ -e "${sp_ub}.genomescope" ]]; then
	  rm -rf -- "${sp_ub}.genomescope"
	fi
		mkdir -p "${sp_ub}.genomescope"

	fastq_files=()
	mapfile -t fastq_files < <(find "${dir_sp_dnaseq}/${sp_ub}" -type f ! -name '.*' \( -name "*.fq" -o -name "*.fastq" -o -name "*.fq.gz" -o -name "*.fastq.gz" \) | sort)
	if [[ ${#fastq_files[@]} -eq 0 ]]; then
	  echo "No FASTQ files were found in: ${dir_sp_dnaseq}/${sp_ub}. Skipping ${task}."
	else
	  printf '%s\n' "${fastq_files[@]}" > input_fastq_files.txt
	  klength=21
	  kmer_lower=1
	  kmer_upper=1000 # Could be automatically calculated with the file from species_genome
	  if [[ -e ./tmp ]]; then
	    rm -rf -- ./tmp
	  fi
	  mkdir -p ./tmp
	  kmc -k${klength} -t${NSLOTS} -m${MEM_PER_HOST} -ci${kmer_lower} -cs${kmer_upper} @input_fastq_files.txt tmp.reads ./tmp
	  kmc_tools transform tmp.reads histogram tmp.reads.histo -cx${kmer_upper}
	  rm -rf -- ./tmp
	  /usr/local/bin/genomescope2.0/genomescope.R -i "tmp.reads.histo" -o "${sp_ub}.genomescope" -k ${klength}

	  if [[ -s "${sp_ub}.genomescope/transformed_linear_plot.png" ]]; then
	    echo "GenomeScope output file was detected. Start compressing."
	    mv_out "tmp.reads.histo" "${sp_ub}.genomescope/kmc.histo.tsv"
	    zip -rq "${sp_ub}.genomescope.zip" "${sp_ub}.genomescope"
	    mv_out "${sp_ub}.genomescope.zip" "${file_sp_genomescope}"
	  fi
	fi
else
	gg_step_skip "${task}"
fi

task="JCVI synteny dotplot"
disable_if_no_input_file "run_jcvi_dotplot" "${file_sp_cds}" "${file_sp_gff}" "${file_sp_genome_fx2tab}"
if [[ ! -s "${file_sp_jcvi_dotplot}" && ${run_jcvi_dotplot} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ -e "${sp_ub}.jcvi_dotplot" ]]; then
    rm -rf -- "${sp_ub}.jcvi_dotplot"
  fi
  mkdir -p "${sp_ub}.jcvi_dotplot"
  cd "${sp_ub}.jcvi_dotplot"

  minimum_scaffold_size=1000000
  jcvi_cscore="0.7"
  gff_feature="mRNA"
  sp=${sp_ub//_/ }
  gff_attribute=$(python "${dir_script}/gff_attribute_identifier.py" --gff_feature "${gff_feature}" --fasta_file "${file_sp_cds}" --gff_file "${file_sp_gff}")
  selected_scaffolds=$(python -c "import sys,pandas as pd; d=pd.read_csv(sys.argv[1],sep='\t',header=0); d['length']=pd.to_numeric(d['length'],errors='coerce'); ids=d.loc[d['length']>=int(sys.argv[2]),'#id'].dropna().astype(str).tolist(); print(','.join(ids))" "${file_sp_genome_fx2tab}" "${minimum_scaffold_size}")
  if [[ -z "${selected_scaffolds//[[:space:],]/}" ]]; then
    echo "No scaffolds >= ${minimum_scaffold_size} bp were found. Falling back to top 20 longest scaffolds."
    selected_scaffolds=$(python -c "import sys,pandas as pd; d=pd.read_csv(sys.argv[1],sep='\t',header=0); d['length']=pd.to_numeric(d['length'],errors='coerce'); d=d.dropna(subset=['#id','length']).sort_values('length',ascending=False); ids=d.loc[:,'#id'].astype(str).head(int(sys.argv[2])).tolist(); print(','.join(ids))" "${file_sp_genome_fx2tab}" 20)
  fi
  if [[ -z "${selected_scaffolds//[[:space:],]/}" ]]; then
    echo "No scaffolds were selected from ${file_sp_genome_fx2tab}. Exiting."
    exit 1
  fi
  echo "Detected GFF attribute that matches to the FASTA sequence names: ${gff_attribute}"
  echo "Selected sequences for the dot-plot visualization: ${selected_scaffolds}"

  seqkit seq "${file_sp_cds}" | sed -e "s/^>${sp_ub}_/>/" > tmp.jcvi_input.cds.fasta
  python -m jcvi.formats.gff bed --primary_only --type="${gff_feature}" --key="${gff_attribute}" "${file_sp_gff}" -o tmp.species1.bed
  python -m jcvi.formats.fasta format tmp.jcvi_input.cds.fasta species1.cds
  touch species1.bed
  selected_scaffold_list=()
  mapfile -t selected_scaffold_list < <(printf '%s' "${selected_scaffolds}" | tr ',' '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' -e '/^$/d')
  for scaffold in "${selected_scaffold_list[@]}"; do
      awk -v scaffold="${scaffold}" '$1 == scaffold' tmp.species1.bed >> species1.bed
  done
  cp_out species1.bed species2.bed
  cp_out species1.cds species2.cds

  python -m jcvi.compara.catalog ortholog species1 species2 --no_strip_names --no_dotplot --cscore="${jcvi_cscore}"
  python -m jcvi.graphics.dotplot species1.species2.anchors --notex --figsize 7x7 --genomenames="${sp}_${sp}" --style='white' --nmax=100000
  python -m jcvi.compara.synteny depth --histogram species1.species2.anchors
  python -m jcvi.compara.synteny screen --minspan=30 --minsize=0 --simple species1.species2.anchors species1.species2.anchors.new

  printf '%s\n' \
    '# y, xstart, xend, rotation, color, label, va,  bed' \
    '0.6,     0.15,    0.95,       0,      , species1, top, species1.bed' \
    '0.4,     0.15,    0.95,       0,      , species2, top, species2.bed' \
    '# edges' \
    'e, 0, 1, species1.species2.anchors.simple' \
    > layout

  printf '%s\n%s\n' "${selected_scaffolds}" "${selected_scaffolds}" > seqids
  python -m jcvi.graphics.karyotype --notex --figsize 7x7 seqids layout

  rm -f -- species2.*
  rm -f -- species1.bed species1.cds
  rm -f -- seqids layout
  rm -f -- species1.species2.last species1.species2.last.filtered
  rm -f -- tmp.scaffold_length.tsv
  rm -f -- tmp.species1.bed
  rm -f -- tmp.jcvi_input.cds.fasta
  mv_out karyotype.pdf "${sp_ub}".karyotype.pdf
  shopt -s nullglob
  for file in *; do
    renamed_file="${file/species1.species2/${sp_ub}}"
    if [[ "${renamed_file}" != "${file}" ]]; then
      mv_out "${file}" "${renamed_file}"
    fi
  done
  shopt -u nullglob

  cd "${dir_sp_tmp}"

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
if is_output_older_than_inputs "^file_sp_" "${file_multispecies_summary}"; then
  summary_flag=0
else
  summary_flag=$?
fi
if [[ ${run_multispecies_summary} -eq 1 && ${summary_flag} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "$(dirname "${file_multispecies_summary}")"
  cd "$(dirname "${file_multispecies_summary}")"

  Rscript "${dir_script}/annotation_summary.r" \
  --dir_species_tree="${dir_pg_output}/species_tree" \
  --dir_species_cds_busco="${dir_pg_output}/species_cds_busco_full" \
  --dir_species_genome_busco="${dir_pg_output}/species_genome_busco_full" \
  --dir_species_annotation="${dir_pg_output}/species_cds_annotation" \
  --dir_species_cds_fx2tab="${dir_pg_output}/species_cds_fx2tab" \
  --dir_species_genome_fx2tab="${dir_pg_output}/species_genome_fx2tab" \
  --file_species_trait="${dir_pg_input}/species_trait/species_trait.tsv" \
  --file_orthogroup_gene_count="${dir_pg_output}/orthofinder/Orthogroups/Orthogroups.GeneCount.tsv" \
  --tree_annotation_dir="${dir_script}/tree_annotation" \
  --min_og_species='auto'

  if [[ -e "Rplots.pdf" ]]; then
    rm -f -- "Rplots.pdf"
  fi

  cd "${dir_sp_tmp}"
else
  gg_step_skip "${task}"
fi

remove_empty_subdirs "${dir_pg_output}"
if [[ ${delete_tmp_dir} -eq 1 ]]; then
    echo "Deleting tmp directory: ${dir_sp_tmp}"
    if [[ -n "${dir_sp_tmp:-}" && "${dir_sp_tmp}" != "/" ]]; then
      rm -rf -- "${dir_sp_tmp}"
    else
      echo "Refusing to delete unsafe tmp directory: '${dir_sp_tmp:-}'"
    fi
fi

echo "$(date): Exiting Singularity environment"

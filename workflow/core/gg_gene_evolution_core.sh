#!/usr/bin/env bash
set -euo pipefail

gg_core_bootstrap="/script/support/gg_core_bootstrap.sh"
if [[ ! -s "${gg_core_bootstrap}" ]]; then
  gg_core_bootstrap="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../support/gg_core_bootstrap.sh"
fi
# shellcheck disable=SC1090
source "${gg_core_bootstrap}"
unset gg_core_bootstrap
gg_source_common_params_from_core "${BASH_SOURCE[0]:-$0}"

### Start: Job-supplied configuration ###
# Configuration variables are provided by gg_gene_evolution_entrypoint.sh.
busco_lineage="${busco_lineage:-${GG_COMMON_BUSCO_LINEAGE:-auto}}"
genetic_code="${genetic_code:-${GG_COMMON_GENETIC_CODE:-1}}"
annotation_species="${annotation_species:-${GG_COMMON_REFERENCE_SPECIES:-auto}}"
input_sequence_mode="${input_sequence_mode:-${GG_COMMON_INPUT_SEQUENCE_MODE:-cds}}"
species_label_parser="${species_label_parser:-${GG_COMMON_SPECIES_LABEL_PARSER:-taxonomic}}"
species_label_regex="${species_label_regex:-${GG_COMMON_SPECIES_LABEL_REGEX:-}}"
species_label_map_tsv="${species_label_map_tsv:-${GG_COMMON_SPECIES_LABEL_MAP_TSV:-}}"
cdskit_localize_model="${cdskit_localize_model:-targeting5-perox-deeploc21-et-v1}"
cdskit_localize_organism_group="${cdskit_localize_organism_group:-auto}"
cdskit_localize_include_features="${cdskit_localize_include_features:-0}"
cdskit_localize_no_model_download="${cdskit_localize_no_model_download:-0}"
run_extract_query_fasta="${run_extract_query_fasta:-1}"
run_extract_primary_fasta="${run_extract_primary_fasta:-1}"
run_generate_expression_matrix="${run_generate_expression_matrix:-0}"
run_collect_gff_info="${run_collect_gff_info:-0}"
run_extract_promoter_fasta="${run_extract_promoter_fasta:-0}"
treevis_query_marker="${treevis_query_marker:-1}"
query_blast_evalue="${query_blast_evalue:-auto}"
query_blast_auto_evalue_maxlen_cutoffs="${query_blast_auto_evalue_maxlen_cutoffs:-40:1000,80:100,150:10,300:1,inf:0.01}"

# Substitution model in CSUBST and mapdNdS
if [[ ${genetic_code} -eq 1 ]]; then
  codon_model="ECMK07+F+R4" # Standard genetic code
else
  codon_model="GY+F+R4" # Non-standard genetic code
fi
### End: Job-supplied configuration ###

### ----------------------------------------------------------------------- ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

gg_bootstrap_core_runtime "${BASH_SOURCE[0]:-$0}" "base" 1 1
delete_preexisting_tmp_dir=${delete_preexisting_tmp_dir:-1}
delete_tmp_dir=${delete_tmp_dir:-1}

gg_core_stage_library="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/stages/gg_gene_evolution_core_functions.sh"
# shellcheck disable=SC1090
source "${gg_core_stage_library}"
unset gg_core_stage_library












# Setting modes
if [[ ${gg_debug_mode:-0} -eq 1 ]]; then
  enable_all_run_flags_for_debug_mode "gg debug mode: All run_* variables are forced to set 1, except for too-time-consuming tasks."
  #run_orthogroup_extraction=0; echo "gg debug mode: run_orthogroup_extraction=${run_orthogroup_extraction}"
  run_codeml_two_ratio=0
  echo "gg debug mode: run_codeml_two_ratio=${run_codeml_two_ratio}"
  run_hyphy_relax=0
  echo "gg debug mode: run_hyphy_relax=${run_hyphy_relax}"
  run_hyphy_relax_reversed=0
  echo "gg debug mode: run_hyphy_relax_reversed=${run_hyphy_relax_reversed}"
  csubst_cutoff_stat="OCNany2spe,0|omegaCany2spe,1"
  echo "gg debug mode: csubst_cutoff_stat=${csubst_cutoff_stat}"
fi
query_blast_method=$(echo "${query_blast_method}" | tr '[:upper:]' '[:lower:]')
mode_gene_evolution=$(echo "${mode_gene_evolution:-query2family}" | tr '[:upper:]' '[:lower:]')
gene_evolution_profile=$(echo "${gene_evolution_profile:-default}" | tr '[:upper:]' '[:lower:]')
input_sequence_mode=$(gg_normalize_input_sequence_mode "${input_sequence_mode}")
csubst_nonsyn_recode=$(echo "${csubst_nonsyn_recode:-${GG_COMMON_CSUBST_NONSYN_RECODE:-no}}" | tr '[:upper:]' '[:lower:]')
csubst_scan_match=$(echo "${csubst_scan_match:-any2spe}" | tr '[:upper:]' '[:lower:]')
csubst_scan_min_event_pp="${csubst_scan_min_event_pp:-0.5}"
csubst_scan_min_support="${csubst_scan_min_support:-2}"
csubst_scan_rate_event_mode=$(echo "${csubst_scan_rate_event_mode:-posterior_sum}" | tr '[:upper:]' '[:lower:]')
csubst_scan_rate_length=$(echo "${csubst_scan_rate_length:-n_rescaled}" | tr '[:upper:]' '[:lower:]')
csubst_scan_rate_exposure=$(echo "${csubst_scan_rate_exposure:-q_weighted}" | tr '[:upper:]' '[:lower:]')
csubst_scan_other_scope=$(echo "${csubst_scan_other_scope:-all}" | tr '[:upper:]' '[:lower:]')
csubst_scan_pvalue_calibration=$(echo "${csubst_scan_pvalue_calibration:-full_scan}" | tr '[:upper:]' '[:lower:]')
csubst_scan_n_permutations="${csubst_scan_n_permutations:-1000}"
csubst_scan_site_plot=$(echo "${csubst_scan_site_plot:-yes}" | tr '[:upper:]' '[:lower:]')
csubst_scan_tree_site_plot_format=$(echo "${csubst_scan_tree_site_plot_format:-pdf}" | tr '[:upper:]' '[:lower:]')
csubst_scan_tree_site_plot_max_sites="${csubst_scan_tree_site_plot_max_sites:-30}"
uniprot_annotation_method=$(echo "${uniprot_annotation_method:-mmseqs2}" | tr '[:upper:]' '[:lower:]')
tree_rooting_method=$(echo "${tree_rooting_method}" | tr '[:upper:]' '[:lower:]')
apply_gene_evolution_profile
if [[ "${mode_gene_evolution}" != "orthogroup" && "${mode_gene_evolution}" != "query2family" ]]; then
  echo "Invalid mode_gene_evolution: ${mode_gene_evolution}"
  echo 'mode_gene_evolution must be either "orthogroup" or "query2family". Exiting.'
  exit 1
fi
if [[ "${tree_rooting_method}" != "notung" && "${tree_rooting_method}" != "midpoint" && "${tree_rooting_method}" != "mad" && "${tree_rooting_method}" != "md" ]]; then
  echo "Invalid tree_rooting_method: ${tree_rooting_method}"
  echo "tree_rooting_method must be one of notung, midpoint, mad, md. Exiting."
  exit 1
fi
if [[ "${uniprot_annotation_method}" != "blastp" && "${uniprot_annotation_method}" != "mmseqs2" ]]; then
  echo "Invalid uniprot_annotation_method: ${uniprot_annotation_method}"
  echo 'uniprot_annotation_method must be either "blastp" or "mmseqs2". Exiting.'
  exit 1
fi
case "${csubst_nonsyn_recode}" in
  no|3di20|dayhoff6|sr6|kgb6|sr4|dayhoff9|dayhoff12|dayhoff15|dayhoff18|srchisq6|kgbauto6)
    ;;
  *)
    echo "Invalid csubst_nonsyn_recode: ${csubst_nonsyn_recode}"
    echo 'csubst_nonsyn_recode must be one of no, 3di20, dayhoff6, sr6, kgb6, sr4, dayhoff9, dayhoff12, dayhoff15, dayhoff18, srchisq6, kgbauto6. Exiting.'
    exit 1
    ;;
esac
case "${csubst_scan_rate_event_mode}" in
  posterior_sum|called)
    ;;
  *)
    echo "Invalid csubst_scan_rate_event_mode: ${csubst_scan_rate_event_mode}"
    echo 'csubst_scan_rate_event_mode must be either "posterior_sum" or "called". Exiting.'
    exit 1
    ;;
esac
case "${csubst_scan_rate_length}" in
  raw|sn_rescaled|n_rescaled)
    ;;
  *)
    echo "Invalid csubst_scan_rate_length: ${csubst_scan_rate_length}"
    echo 'csubst_scan_rate_length must be one of raw, sn_rescaled, n_rescaled. Exiting.'
    exit 1
    ;;
esac
case "${csubst_scan_rate_exposure}" in
  q_weighted|state_aware|raw_branch_length)
    ;;
  *)
    echo "Invalid csubst_scan_rate_exposure: ${csubst_scan_rate_exposure}"
    echo 'csubst_scan_rate_exposure must be one of q_weighted, state_aware, raw_branch_length. Exiting.'
    exit 1
    ;;
esac
case "${csubst_scan_other_scope}" in
  all|sister)
    ;;
  *)
    echo "Invalid csubst_scan_other_scope: ${csubst_scan_other_scope}"
    echo 'csubst_scan_other_scope must be either "all" or "sister". Exiting.'
    exit 1
    ;;
esac
case "${csubst_scan_pvalue_calibration}" in
  none|candidate_fixed|full_scan)
    ;;
  *)
    echo "Invalid csubst_scan_pvalue_calibration: ${csubst_scan_pvalue_calibration}"
    echo 'csubst_scan_pvalue_calibration must be one of none, candidate_fixed, full_scan. Exiting.'
    exit 1
    ;;
esac
case "${csubst_scan_site_plot}" in
  yes|no)
    ;;
  *)
    echo "Invalid csubst_scan_site_plot: ${csubst_scan_site_plot}"
    echo 'csubst_scan_site_plot must be either "yes" or "no". Exiting.'
    exit 1
    ;;
esac
case "${csubst_scan_tree_site_plot_format}" in
  pdf|png|svg)
    ;;
  *)
    echo "Invalid csubst_scan_tree_site_plot_format: ${csubst_scan_tree_site_plot_format}"
    echo 'csubst_scan_tree_site_plot_format must be one of pdf, png, svg. Exiting.'
    exit 1
    ;;
esac
apply_gene_evolution_input_sequence_mode
if [[ "${mode_gene_evolution}" == "query2family" && ${run_query_blast} -eq 1 ]]; then
  if [[ "${query_blast_method}" != "tblastn" && "${query_blast_method}" != "diamond" ]]; then
    echo "Invalid query_blast_method: ${query_blast_method}"
    echo 'query_blast_method must be either "tblastn" or "diamond". Exiting.'
    exit 1
  fi
fi
case "${mode_gene_evolution}" in
  "orthogroup")
    dir_output_active="${gg_workspace_output_dir}/orthogroup"
    file_orthogroup_genecount_selected="${gg_workspace_output_dir}/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv"
    if [[ ! -s "${file_orthogroup_genecount_selected}" ]]; then
      echo "Orthogroup gene-count table not found: ${file_orthogroup_genecount_selected}"
      exit 1
    fi
    if [[ ! "${GG_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] || [[ ${GG_ARRAY_TASK_ID} -lt 1 ]]; then
      echo "Invalid GG_ARRAY_TASK_ID value (must be a positive integer): ${GG_ARRAY_TASK_ID}"
      exit 1
    fi
    num_orthogroups=$(awk 'END { print (NR > 0 ? NR - 1 : 0) }' "${file_orthogroup_genecount_selected}")
    if [[ ${num_orthogroups} -le 0 ]]; then
      echo "No orthogroup rows were found in: ${file_orthogroup_genecount_selected}"
      exit 1
    fi
    if [[ ${GG_ARRAY_TASK_ID} -gt ${num_orthogroups} ]]; then
      echo "GG_ARRAY_TASK_ID=${GG_ARRAY_TASK_ID} is out of range for ${num_orthogroups} orthogroups in: ${file_orthogroup_genecount_selected}"
      exit 1
    fi
    og_id=$(awk -F'\t' -v row="${GG_ARRAY_TASK_ID}" 'NR == (row + 1) { print $1; exit }' "${file_orthogroup_genecount_selected}")
    if [[ -z "${og_id}" ]]; then
      echo "Failed to read OrthoGroup ID at row ${GG_ARRAY_TASK_ID} from: ${file_orthogroup_genecount_selected}"
      exit 1
    fi
    echo "OrthoGroup ID: ${og_id}"
    run_extract_query_fasta=0
    run_query_blast=0
    run_orthogroup_extraction=0
    ;;
  "query2family")
    dir_output_active="${gg_workspace_output_dir}/query2family"
    dir_genelist="${gg_workspace_input_dir}/query_gene"
    if [[ ! -d "${dir_genelist}" ]]; then
      echo "Input directory does not exist: ${dir_genelist}"
      exit 1
    fi
    files=()
    mapfile -t files < <(find "${dir_genelist}" -mindepth 1 -maxdepth 1 -type f ! -name '.*' | sort)
    if [[ ${#files[@]} -eq 0 ]]; then
      echo "Input directory is empty: ${dir_genelist}"
      exit 1
    fi
    if [[ ! "${GG_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] || [[ ${GG_ARRAY_TASK_ID} -lt 1 ]]; then
      echo "Invalid GG_ARRAY_TASK_ID value (must be a positive integer): ${GG_ARRAY_TASK_ID}"
      exit 1
    fi
    idx=$((GG_ARRAY_TASK_ID - 1))
    if [[ ${idx} -ge ${#files[@]} ]]; then
      echo "Input genelist file not found, probably due to out-of-range GG_ARRAY_TASK_ID=${GG_ARRAY_TASK_ID}. Number of files: ${#files[@]}"
      exit 1
    fi
    file_query_gene="${files[${idx}]}"
    og_id="$(basename "${file_query_gene}")"
    echo "mode_gene_evolution=query2family: ${#files[@]} input genelist file(s) were detected in ${dir_genelist}/"
    echo "mode_gene_evolution=query2family: Input genelist file = ${file_query_gene}"
    echo "output file prefix: ${og_id}"
    if [[ ! -f "${file_query_gene}" ]]; then
      echo "Input genelist file not found, probably due to out-of-range GG_ARRAY_TASK_ID=${GG_ARRAY_TASK_ID}: ${file_query_gene}"
      exit 1
    fi
    ;;
  *)
    echo "Invalid mode_gene_evolution: ${mode_gene_evolution}. Exiting."
    exit 1
    ;;
esac
if [[ -z "$og_id" ]]; then
  echo "og_id is empty. Exiting."
  exit 1
fi

dir_sp_genome="${gg_workspace_input_dir}/species_genome"
dir_sp_gff="${gg_workspace_input_dir}/species_gff"
dir_sp_expression="${gg_workspace_input_dir}/species_expression"
dir_sp_cds="${gg_workspace_input_dir}/species_cds"
dir_sp_protein_input="$(gg_species_protein_input_dir_path "${gg_workspace_input_dir}")"
dir_sp_blastdb="${gg_workspace_output_dir}/species_cds_blastdb"
file_species_genetic_code="$(gg_species_genetic_code_table_path "${gg_workspace_input_dir}")"
file_species_genetic_code_resolved="${dir_output_active}/parameters/${og_id}_species_genetic_code.resolved.tsv"
annotation_species_resolved=""
treevis_clade_ortholog_prefix=""
annotation_species_candidates=()
if [[ "${input_sequence_mode}" == "protein" ]]; then
  mapfile -t annotation_species_candidates < <(gg_species_names_from_fasta_dir "${dir_sp_protein_input}")
fi
if [[ ${#annotation_species_candidates[@]} -eq 0 ]]; then
  mapfile -t annotation_species_candidates < <(gg_species_names_from_fasta_dir "${dir_sp_cds}")
fi
if annotation_species_resolved=$(gg_resolve_annotation_species "${annotation_species}" "${annotation_species_candidates[@]}"); then
  if [[ -n "${annotation_species_resolved}" ]]; then
    treevis_clade_ortholog_prefix="${annotation_species_resolved}_"
  fi
fi
file_sp_trait="${gg_workspace_input_dir}/species_trait/species_trait.tsv"
file_og="${gg_workspace_output_dir}/orthofinder/Orthogroups_filtered/Orthogroups.selected.tsv"
file_og_parameters_dir="${dir_output_active}/parameters"
dir_species_tree="${gg_workspace_output_dir}/species_tree"
dir_species_tree_summary="${dir_species_tree}/species_tree_summary"
if [[ -s "${dir_species_tree_summary}/dated_species_tree.nwk" ]]; then
  species_tree_basename="dated_species_tree"
  species_tree="${dir_species_tree_summary}/${species_tree_basename}.nwk"
elif [[ -s "${dir_species_tree_summary}/undated_species_tree.nwk" ]]; then
  species_tree_basename="undated_species_tree"
  species_tree="${dir_species_tree_summary}/${species_tree_basename}.nwk"
else
  species_tree_basename="undated_species_tree"
  species_tree="${dir_species_tree_summary}/${species_tree_basename}.nwk"
fi
species_tree_generax="${file_og_parameters_dir}/${species_tree_basename}.generax.nwk" # generated later
species_tree_pruned="${file_og_parameters_dir}/${species_tree_basename}.pruned.nwk"
ensure_dir "${file_og_parameters_dir}"
notung_jar="/usr/local/bin/Notung.jar"
dir_rpsblastdb="/usr/local/db/Pfam_LE"

# Directory PATHs
# Directories for temporary files
dir_tmp="${dir_output_active}/tmp/${GG_ARRAY_TASK_ID}_${og_id}" #_${RANDOM}

# File PATHs
# Alignment and gene tree preparation and others
file_og_query_aa_fasta="${dir_output_active}/query_aa_fasta/${og_id}_query.aa.fa.gz"
file_og_query_blast="${dir_output_active}/query_blast/${og_id}_query_blast.tsv"
file_og_cds_fasta="${dir_output_active}/cds_fasta/${og_id}_cds.fa.gz"
file_og_pep_fasta="${dir_output_active}/protein_fasta/${og_id}_pep.fa.gz"
file_og_rpsblast="${dir_output_active}/rpsblast/${og_id}_rpsblast.tsv"
file_og_uniprot_annotation="${dir_output_active}/uniprot_annotation/${og_id}_uniprot.tsv"
file_og_cdskit_localize="${dir_output_active}/cdskit_localize/${og_id}_cdskit_localize.tsv"
file_og_mafft="${dir_output_active}/mafft/${og_id}_cds.aln.fa.gz"
file_og_maxalign="${dir_output_active}/maxalign/${og_id}_cds.maxalign.fa.gz"
file_og_trimal="${dir_output_active}/trimal/${og_id}_cds.trimal.fa.gz"
file_og_clipkit="${dir_output_active}/clipkit/${og_id}_cds.clipkit.fa.gz"
file_og_clipkit_log="${dir_output_active}/clipkit_log/${og_id}_cds.clipkit.log"
file_og_iqtree_tree="${dir_output_active}/iqtree_tree/${og_id}_iqtree.nwk"
file_og_iqtree_generax_ufboot="${dir_output_active}/generax_ufboot_tree/${og_id}_generax.ufboot.nwk"
file_og_orthogroup_extraction_nwk="${dir_output_active}/orthogroup_extraction_nwk/${og_id}_orthogroup_extraction.nwk"
file_og_orthogroup_extraction_rooted_nwk="${dir_output_active}/orthogroup_extraction_rooted_nwk/${og_id}_orthogroup_extraction.rooted.nwk"
file_og_orthogroup_extraction_fasta="${dir_output_active}/orthogroup_extraction_fasta/${og_id}_orthogroup_extraction.fa.gz"
file_og_generax_nhx="${dir_output_active}/generax_tree/${og_id}_generax.nhx"
file_og_generax_nwk="${dir_output_active}/generax_nwk/${og_id}_generax.nwk"
file_og_generax_xml="${dir_output_active}/generax_xml/${og_id}_generax.xml"
file_og_rooted_tree="${dir_output_active}/rooted_tree/${og_id}_root.nwk"
file_og_rooted_log="${dir_output_active}/rooted_tree_log/${og_id}_root.txt"
file_og_notung_reconcil="${dir_output_active}/notung_reconcile/${og_id}_notung_reconcile.zip"
file_og_dated_tree="${dir_output_active}/dated_tree/${og_id}_dated.nwk"
file_og_dated_tree_log="${dir_output_active}/dated_tree_log/${og_id}_dated.log.txt"
file_og_mapdnds_parameter="${dir_output_active}/mapdnds_parameter/${og_id}_parameter.zip"
file_og_mapdnds_dn="${dir_output_active}/mapdnds_dn_tree/${og_id}_mapdNdS.dN.nwk"
file_og_mapdnds_ds="${dir_output_active}/mapdnds_ds_tree/${og_id}_mapdNdS.dS.nwk"
file_og_codeml_two_ratio="${dir_output_active}/codeml_two_ratio/${og_id}_codeml.two_ratio.tsv"
file_og_hyphy_dnds="${dir_output_active}/hyphy_dnds/${og_id}_hyphy.dnds.json"
file_og_hyphy_relax="${dir_output_active}/hyphy_relax/${og_id}_hyphy.relax.json"
file_og_hyphy_relax_reversed="${dir_output_active}/hyphy_relax_reversed/${og_id}_hyphy.relax.reversed.json"
file_og_expression="${dir_output_active}/character_expression/${og_id}_expression.tsv"
file_og_gff_info="${dir_output_active}/character_gff_info/${og_id}_gff.tsv"
file_og_scm_intron_summary="${dir_output_active}/scm_intron_summary/${og_id}_scm.intron.tsv"
file_og_scm_intron_plot="${dir_output_active}/scm_intron_plot/${og_id}_scm.intron.pdf"
# Cis-regulatory motif
file_og_promoter_fasta="${dir_output_active}/promoter_fasta/${og_id}_promoter.fa.gz"
file_og_meme="${dir_output_active}/meme/${og_id}_meme.xml"
file_og_fimo="${dir_output_active}/fimo/${og_id}_fimo.tsv"
file_og_fimo_collapsed="${dir_output_active}/fimo_collapsed/${og_id}_fimo.collapsed.tsv"
file_og_l1ou_fit_rdata="${dir_output_active}/l1ou_fit_rdata/${og_id}_l1ou.RData"
file_og_l1ou_fit_conv_rdata="${dir_output_active}/l1ou_fit_conv_rdata/${og_id}_l1ou.conv.RData"
file_og_l1ou_fit_tree="${dir_output_active}/l1ou_fit_tree/${og_id}_l1ou.tree.tsv"
file_og_l1ou_fit_regime="${dir_output_active}/l1ou_fit_regime/${og_id}_l1ou.regime.tsv"
file_og_l1ou_fit_leaf="${dir_output_active}/l1ou_fit_leaf/${og_id}_l1ou.leaf.tsv"
file_og_l1ou_fit_plot="${dir_output_active}/l1ou_fit_plot/${og_id}_l1ou.pdf"
# Protein convergence analysis
file_og_iqtree_anc="${dir_output_active}/iqtree_anc/${og_id}_iqtree.anc.zip"
file_og_csubst_b="${dir_output_active}/csubst_b/${og_id}_csubst_b.tsv"
file_og_csubst_cb_2="${dir_output_active}/csubst_cb_2/${og_id}_csubst_cb_2.tsv"
file_og_csubst_cb_stats="${dir_output_active}/csubst_cb_stats/${og_id}_csubst_cb_stats.tsv"
file_og_csubst_scan="${dir_output_active}/csubst_scan/${og_id}_csubst_scan.tsv"
file_og_csubst_scan_units="${dir_output_active}/csubst_scan_units/${og_id}_csubst_scan_units.tsv"
file_og_csubst_scan_foreground_branch="${dir_output_active}/csubst_scan_foreground_branch/${og_id}_csubst_foreground_branch.txt"
file_og_csubst_scan_plot="${dir_output_active}/csubst_scan_plot/${og_id}_csubst_scan.tree_site.${csubst_scan_tree_site_plot_format}"
file_og_csubst_scan_log="${dir_output_active}/csubst_scan_log/${og_id}_csubst_scan.log"
if [[ ${csubst_max_arity} -gt 2 ]]; then
  for ((i = 3; i <= csubst_max_arity; i++)); do
    declare file_og_csubst_cb_${i}="${dir_output_active}/csubst_cb_${i}/${og_id}.csubst_cb_${i}.tsv"
  done
fi
# PGLS output
file_og_gene_pgls="${dir_output_active}/pgls_gene_tree/${og_id}_gene_tree_PGLS.tsv"
file_og_gene_pgls_plot="${dir_output_active}/pgls_gene_tree_plot/${og_id}_gene_PGLS.barplot.pdf"
file_og_species_pgls="${dir_output_active}/pgls_species_tree/${og_id}_species_PGLS.tsv"
file_og_species_pgls_plot="${dir_output_active}/pgls_species_tree_plot/${og_id}_species_PGLS.barplot.pdf"
# Summary
file_og_stat_branch="${dir_output_active}/stat_branch/${og_id}_stat.branch.tsv"
file_og_stat_tree="${dir_output_active}/stat_tree/${og_id}_stat.tree.tsv"
file_og_amas_original="${dir_output_active}/amas_original/${og_id}_amas.original.tsv"
file_og_amas_cleaned="${dir_output_active}/amas_cleaned/${og_id}_amas.cleaned.tsv"
file_og_tree_plot="${dir_output_active}/tree_plot/${og_id}_tree_plot.pdf"
file_og_synteny="${dir_output_active}/synteny/${og_id}_synteny.tsv"
# Pruned datasets
file_og_untrimmed_aln_pruned="${dir_output_active}/pruned_untrimmed_alignment/${og_id}_cds.untrimmed.pruned.fa.gz"
file_og_trimmed_aln_pruned="${dir_output_active}/pruned_trimmed_alignment/${og_id}_cds.trimmed.pruned.fa.gz"
file_og_unrooted_tree_pruned="${dir_output_active}/pruned_unrooted_tree/${og_id}_unrooted.pruned.nwk"
file_og_rooted_tree_pruned="${dir_output_active}/pruned_rooted_tree/${og_id}_rooted.pruned.nwk"
file_og_dated_tree_pruned="${dir_output_active}/pruned_dated_tree/${og_id}_dated.pruned.nwk"
file_og_primary_fasta="${file_og_cds_fasta}"
amas_data_type="dna"
if [[ "${input_sequence_mode}" == "protein" ]]; then
  file_og_primary_fasta="${file_og_pep_fasta}"
  amas_data_type="aa"
fi

# Define intermediate files for downstream analysis.
# These variables are updated via helper functions to keep routing logic explicit.






prepare_species_tree_pruned || true
set_default_analysis_files

memory_notung=$(gg_memory_fraction_gb "${GG_MEM_TOOL_GB}" 1 4)

echo "Checking parameter conflicts..."
if [[ ${run_trimal} -eq 1 && ${run_clipkit} -eq 1 ]]; then
  echo 'Both of ${run_trimal} and ${run_clipkit} are set to 1.'
  echo '${run_trimal} is deactivated. ${run_clipkit} is still active.'
  run_trimal=0
fi
if [[ ! -d "${dir_sp_gff}" ]] || [[ -z "$(find "${dir_sp_gff}" -mindepth 1 -maxdepth 1 -print -quit 2> /dev/null)" ]]; then
  if [[ ${run_collect_gff_info} -eq 1 ]]; then
    echo "\${run_collect_gff_info} is deactivated. Empty input: ${dir_sp_gff}"
    run_collect_gff_info=0
  fi
  if [[ ${run_scm_intron} -eq 1 ]]; then
    echo "\${run_scm_intron} is deactivated. Empty input: ${dir_sp_gff}"
    run_scm_intron=0
  fi
fi
if [[ -d "${dir_sp_expression}" ]] && [[ -n "$(find "${dir_sp_expression}" -mindepth 1 -maxdepth 1 -print -quit 2> /dev/null)" ]]; then
  echo "\${dir_sp_expression} is not empty. Continued: ${dir_sp_expression}"
else
  echo "\${dir_sp_expression} is empty: ${dir_sp_expression}"
  echo '${run_generate_expression_matrix}, ${run_tree_pruning}, and other options are deactivated.'
  run_tree_pruning=0
  run_generate_expression_matrix=0
  run_l1ou=0
fi
if [[ -d "${dir_sp_genome}" ]] && [[ -n "$(find "${dir_sp_genome}" -mindepth 1 -maxdepth 1 -print -quit 2> /dev/null)" ]]; then
  echo "\${dir_sp_genome} is not empty. Continued: ${dir_sp_genome}"
else
  echo "\${dir_sp_genome} is empty: ${dir_sp_genome}"
  echo '${run_extract_promoter_fasta} and ${run_fimo} are deactivated.'
  run_extract_promoter_fasta=0
  run_fimo=0
fi
echo "Checking preexisting tmp directory."
if [[ -e "${dir_tmp}" && ${delete_preexisting_tmp_dir} -eq 1 ]]; then
  echo "$(date): Deleting preexisting ${dir_tmp}"
  shopt -s nullglob
  stale_tmp_paths=("${dir_output_active}/tmp/${GG_ARRAY_TASK_ID}_"*)
  shopt -u nullglob
  if [[ ${#stale_tmp_paths[@]} -gt 0 ]]; then
    rm -rf -- "${stale_tmp_paths[@]}"
  fi
fi
if [[ ! -e "${dir_tmp}" ]]; then
  echo "Making ${dir_tmp}"
  mkdir -p "${dir_tmp}"
fi
cd "${dir_tmp}"
echo "Working at: $(pwd)"

trap cleanup_tmp_dir_on_normal_exit EXIT


gg_core_execution_stage_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/stages/gg_gene_evolution"
# shellcheck source=/dev/null
source "${gg_core_execution_stage_dir}/01_query_and_fasta.sh"
# shellcheck source=/dev/null
source "${gg_core_execution_stage_dir}/02_annotation_and_alignment.sh"
# shellcheck source=/dev/null
source "${gg_core_execution_stage_dir}/03_gene_tree.sh"
# shellcheck source=/dev/null
source "${gg_core_execution_stage_dir}/04_traits_motifs_and_pruning.sh"
# shellcheck source=/dev/null
source "${gg_core_execution_stage_dir}/05_selection_models.sh"
# shellcheck source=/dev/null
source "${gg_core_execution_stage_dir}/06_pgls_and_csubst.sh"
# shellcheck source=/dev/null
source "${gg_core_execution_stage_dir}/07_summaries.sh"
unset gg_core_execution_stage_dir

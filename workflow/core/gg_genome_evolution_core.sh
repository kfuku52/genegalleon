#!/usr/bin/env bash
set -euo pipefail

# Unified genome evolution pipeline (single-file form).
# Former inline stages are now written directly in this file.

gg_core_bootstrap="/script/support/gg_core_bootstrap.sh"
if [[ ! -s "${gg_core_bootstrap}" ]]; then
  gg_core_bootstrap="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../support/gg_core_bootstrap.sh"
fi
# shellcheck disable=SC1090
source "${gg_core_bootstrap}"
unset gg_core_bootstrap
gg_source_common_params_from_core "${BASH_SOURCE[0]:-$0}"

#run a busco analysis for all species in species_cds and extract shared complete BUSCO genes among all species.
#then make cds fasta files including all species' cds sequenes for each BUSCO genes and make alignment
#then concatenate all fasta files for further downstream analysis -> task = trimal in gg_species_tree.

### Start: Job-supplied configuration ###

# Configuration variables are provided by gg_genome_evolution_entrypoint.sh.
genetic_code="${genetic_code:-${GG_COMMON_GENETIC_CODE:-1}}"
input_sequence_mode="${input_sequence_mode:-${GG_COMMON_INPUT_SEQUENCE_MODE:-cds}}"
busco_lineage="${busco_lineage:-${GG_COMMON_BUSCO_LINEAGE:-auto}}"
species_tree_rooting="${species_tree_rooting:-taxonomy}"
annotation_species="${annotation_species:-${GG_COMMON_REFERENCE_SPECIES:-auto}}"
species_label_parser="${species_label_parser:-${GG_COMMON_SPECIES_LABEL_PARSER:-taxonomic}}"
species_label_regex="${species_label_regex:-${GG_COMMON_SPECIES_LABEL_REGEX:-}}"
species_label_map_tsv="${species_label_map_tsv:-${GG_COMMON_SPECIES_LABEL_MAP_TSV:-}}"
omark_db_path="${omark_db_path:-auto}"
run_cds_translation="${run_cds_translation:-1}"
run_species_omark="${run_species_omark:-0}"
run_build_species_busco_summary="${run_build_species_busco_summary:-1}"
run_build_species_omark_summary="${run_build_species_omark_summary:-1}"
run_extract_species_tree_fasta="${run_extract_species_tree_fasta:-1}"
run_single_copy_ortholog_decay_plot="${run_single_copy_ortholog_decay_plot:-1}"
orthogroup_decay_replicates="${orthogroup_decay_replicates:-1000}"
orthogroup_decay_species_counts="${orthogroup_decay_species_counts:-auto}"
orthogroup_decay_seed="${orthogroup_decay_seed:-1}"
orthofinder_core_filters="${orthofinder_core_filters:-busco_complete_pct:ge:80,num_seq:le:100000}"
orthofinder_core_rank="${orthofinder_core_rank:-num_seq:asc,busco_complete_pct:desc}"
orthofinder_core_method="${orthofinder_core_method:-max-pd}"
run_busco_dupaware_extract_fasta="${run_busco_dupaware_extract_fasta:-0}"
run_busco_dupaware_mafft="${run_busco_dupaware_mafft:-0}"
run_busco_dupaware_trimal="${run_busco_dupaware_trimal:-0}"
run_busco_dupaware_iqtree_dna="${run_busco_dupaware_iqtree_dna:-0}"
run_busco_dupaware_iqtree_pep="${run_busco_dupaware_iqtree_pep:-0}"
run_busco_dupaware_notung_root_dna="${run_busco_dupaware_notung_root_dna:-0}"
run_busco_dupaware_notung_root_pep="${run_busco_dupaware_notung_root_pep:-0}"
run_busco_dupaware_root_dna="${run_busco_dupaware_root_dna:-0}"
run_busco_dupaware_root_pep="${run_busco_dupaware_root_pep:-0}"
run_busco_dupaware_grampa_dna="${run_busco_dupaware_grampa_dna:-0}"
run_busco_dupaware_grampa_pep="${run_busco_dupaware_grampa_pep:-0}"
mcmctree_divergence_time_constraints_str="${mcmctree_divergence_time_constraints_str:-}"
grampa_h1="${grampa_h1:-}"
target_branch_go="${target_branch_go:-}"
orthogroup_copy_number_max_size_differential="${orthogroup_copy_number_max_size_differential:-9999999}"
run_orthogroup_copy_number_trait_pgls="${run_orthogroup_copy_number_trait_pgls:-0}"
orthogroup_copy_number_trait="${orthogroup_copy_number_trait:-all}"
orthogroup_copy_number_trait_min_species="${orthogroup_copy_number_trait_min_species:-4}"
orthogroup_copy_number_trait_family_ids="${orthogroup_copy_number_trait_family_ids:-}"
orthogroup_copy_number_trait_family_file="${orthogroup_copy_number_trait_family_file:-}"
orthogroup_copy_number_trait_max_families="${orthogroup_copy_number_trait_max_families:-all}"
orthogroup_copy_number_trait_p_adjust_method="${orthogroup_copy_number_trait_p_adjust_method:-BH}"
orthogroup_copy_number_trait_alpha="${orthogroup_copy_number_trait_alpha:-0.05}"
orthogroup_copy_number_trait_plot_top_n="${orthogroup_copy_number_trait_plot_top_n:-50}"
file_trait="${file_trait:-auto}"
mcmctree_divergence_time_constraints=()
if [[ -n "${mcmctree_divergence_time_constraints_str:-}" ]]; then
  IFS='|' read -r -a mcmctree_divergence_time_constraints <<< "${mcmctree_divergence_time_constraints_str}"
fi

### End: Job-supplied configuration ###

### ----------------------------------------------------------------------- ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

gg_bootstrap_core_runtime "${BASH_SOURCE[0]:-$0}" "base" 1 1
# shellcheck disable=SC1090
source "${gg_support_dir}/gg_busco.sh"
delete_tmp_dir=${delete_tmp_dir:-1}
busco_lineage_resolved=""
omark_db_resolved=""

gg_core_stage_library="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/stages/gg_genome_evolution_core_functions.sh"
# shellcheck disable=SC1090
source "${gg_core_stage_library}"
unset gg_core_stage_library



input_sequence_mode=$(gg_normalize_input_sequence_mode "${input_sequence_mode}")

species_protein_ready=0
species_protein_source=""
species_tree_sequence_label="cds"
species_tree_busco_mode="transcriptome"
single_copy_fasta_glob="*.cds.fa.gz"
single_copy_aln_glob="*.cds.aln.fa.gz"
single_copy_trimal_glob="*.trimal.fa.gz"
single_copy_fasta_suffix=".cds.fa.gz"
single_copy_aln_suffix=".cds.aln.fa.gz"
single_copy_trimal_suffix=".trimal.fa.gz"
genome_busco_fasta_glob="*.busco.cds.fa.gz"
genome_busco_aln_glob="*.busco.cds.aln.fa.gz"
genome_busco_trimal_glob="*.busco.trimal.fa.gz"
genome_busco_fasta_suffix=".busco.cds.fa.gz"
genome_busco_aln_suffix=".busco.cds.aln.fa.gz"
genome_busco_trimal_suffix=".busco.trimal.fa.gz"
if [[ "${input_sequence_mode}" == "protein" ]]; then
  species_tree_sequence_label="pep"
  species_tree_busco_mode="proteins"
  single_copy_fasta_glob="*.pep.fa.gz"
  single_copy_aln_glob="*.pep.aln.fa.gz"
  single_copy_trimal_glob="*.pep.trimal.fa.gz"
  single_copy_fasta_suffix=".pep.fa.gz"
  single_copy_aln_suffix=".pep.aln.fa.gz"
  single_copy_trimal_suffix=".pep.trimal.fa.gz"
  genome_busco_fasta_glob="*.busco.pep.fa.gz"
  genome_busco_aln_glob="*.busco.pep.aln.fa.gz"
  genome_busco_trimal_glob="*.busco.pep.trimal.fa.gz"
  genome_busco_fasta_suffix=".busco.pep.fa.gz"
  genome_busco_aln_suffix=".busco.pep.aln.fa.gz"
  genome_busco_trimal_suffix=".busco.pep.trimal.fa.gz"
fi





species_tree_rooting_method=""
species_tree_rooting_value=""



# Directories
dir_sp_cds="${gg_workspace_input_dir}/species_cds"
dir_sp_protein_input="$(species_protein_input_dir_path)"
file_species_genetic_code="$(species_genetic_code_table_path)"
file_species_genetic_code_resolved="${gg_workspace_downloads_dir}/tmp/species_genetic_code.resolved.tsv"
dir_og_rooted_tree="${gg_workspace_output_dir}/orthogroup/rooted_tree"
annotation_species_resolved=""
annotation_species_candidates=()
file_go_annotation=""
mapfile -t annotation_species_candidates < <(gg_species_names_from_annotation_dir "${gg_workspace_output_dir}/species_cds_annotation")
if [[ ${#annotation_species_candidates[@]} -eq 0 && "${input_sequence_mode}" == "protein" ]]; then
  mapfile -t annotation_species_candidates < <(gg_species_names_from_fasta_dir "${dir_sp_protein_input}")
fi
if [[ ${#annotation_species_candidates[@]} -eq 0 ]]; then
  mapfile -t annotation_species_candidates < <(gg_species_names_from_fasta_dir "${dir_sp_cds}")
fi
if annotation_species_resolved=$(gg_resolve_annotation_species "${annotation_species}" "${annotation_species_candidates[@]}"); then
  if [[ -n "${annotation_species_resolved}" ]]; then
    file_go_annotation="${gg_workspace_output_dir}/species_cds_annotation/${annotation_species_resolved}_annotation.tsv"
  fi
fi

# Species tree
dir_species_tree="${gg_workspace_output_dir}/species_tree"
if [[ "${input_sequence_mode}" == "protein" ]]; then
  dir_species_busco_full="${gg_workspace_output_dir}/species_protein_busco_full"
  dir_species_busco_short="${gg_workspace_output_dir}/species_protein_busco_short"
else
  dir_species_busco_full="${gg_workspace_output_dir}/species_cds_busco_full"
  dir_species_busco_short="${gg_workspace_output_dir}/species_cds_busco_short"
fi
dir_single_copy_fasta="${dir_species_tree}/single_copy_cds_fasta"
dir_single_copy_mafft="${dir_species_tree}/single_copy_mafft"
dir_single_copy_trimal="${dir_species_tree}/single_copy_trimal"
dir_single_copy_iqtree_dna="${dir_species_tree}/single_copy_iqtree_dna"
dir_single_copy_iqtree_pep="${dir_species_tree}/single_copy_iqtree_pep"
dir_astral_dna="${dir_species_tree}/single_copy_astral_dna"
dir_astral_pep="${dir_species_tree}/single_copy_astral_pep"
dir_species_tree_summary="${dir_species_tree}/species_tree_summary"
dir_concat_fasta="${dir_species_tree}/concatenated_alignment"
dir_concat_iqtree_dna="${dir_species_tree}/concatenated_iqtree_dna"
dir_concat_iqtree_pep="${dir_species_tree}/concatenated_iqtree_pep"
dir_mcmctree2="${dir_species_tree}/mcmctree_main"
dir_tmp="${dir_species_tree}/tmp"
dir_nwkit_download_dir="${gg_workspace_downloads_dir}/nwkit_downloads"

# Orthogroup
dir_sp_protein="${gg_workspace_downloads_dir}/tmp/species_protein"
dir_orthofinder="${gg_workspace_output_dir}/orthofinder"
dir_orthofinder_og="${dir_orthofinder}/Orthogroups"
dir_orthofinder_filtered="${dir_orthofinder}/Orthogroups_filtered"
dir_orthofinder_hog2og="${dir_orthofinder}/hog2og"
dir_orthogroup_decay="${dir_orthofinder}/single_copy_ortholog_decay"

# Genome evolution
dir_genome_evolution="${gg_workspace_output_dir}/genome_evolution"
dir_species_omamer="${dir_genome_evolution}/omamer_search"
dir_species_omark="${dir_genome_evolution}/omark"
dir_busco_fasta="${dir_genome_evolution}/busco_cds_fasta"
dir_busco_mafft="${dir_genome_evolution}/busco_mafft"
dir_busco_trimal="${dir_genome_evolution}/busco_trimal"
dir_busco_iqtree_dna="${dir_genome_evolution}/busco_iqtree_dna"
dir_busco_iqtree_pep="${dir_genome_evolution}/busco_iqtree_pep"
dir_busco_notung_dna="${dir_genome_evolution}/busco_notung_dna"
dir_busco_notung_pep="${dir_genome_evolution}/busco_notung_pep"
dir_busco_rooted_txt_dna="${dir_genome_evolution}/busco_rooted_txt_dna"
dir_busco_rooted_txt_pep="${dir_genome_evolution}/busco_rooted_txt_pep"
dir_busco_rooted_nwk_dna="${dir_genome_evolution}/busco_rooted_nwk_dna"
dir_busco_rooted_nwk_pep="${dir_genome_evolution}/busco_rooted_nwk_pep"
dir_cafe="${dir_genome_evolution}/cafe"
dir_orthogroup_copy_number="${dir_genome_evolution}/orthogroup_copy_number"
dir_cafe_output="${dir_cafe}/cafe_output"
dir_orthogroup_copy_number_trait_pgls="${dir_orthogroup_copy_number}/trait_pgls"
if [[ "${file_trait}" == "auto" ]]; then
  file_trait="${gg_workspace_input_dir}/species_trait/species_trait.tsv"
fi

# Output files
# Species tree
file_species_busco_summary_table="${dir_species_tree}/busco_summary_table/busco_summary.tsv"
file_astral_tree_dna_q1="${dir_astral_dna}/single_copy.astral.dna.q1.nwk"     # Quartet supports for the main topology; The lengths of terminal branches are not computed.
file_astral_tree_dna="${dir_astral_dna}/single_copy.astral.dna.optimized.nwk" # ASTRAL topology with branch lengths optimized by IQ-TREE.
file_astral_log_dna="${dir_astral_dna}/single_copy.astral.dna.log"
file_astral_tree_pep_q1="${dir_astral_pep}/single_copy.astral.pep.q1.nwk"     # Quartet supports for the main topology; The lengths of terminal branches are not computed.
file_astral_tree_pep="${dir_astral_pep}/single_copy.astral.pep.optimized.nwk" # ASTRAL topology with branch lengths optimized by IQ-TREE.
file_astral_log_pep="${dir_astral_pep}/single_copy.astral.pep.log"
file_concat_cds="${dir_concat_fasta}/concat.cds.trimal.fa.gz"
file_concat_pep="${dir_concat_fasta}/concat.pep.trimal.fa.gz"
file_concat_iqtree_dna="${dir_concat_iqtree_dna}/concat.dna.iqtree.nwk"
file_concat_iqtree_pep="${dir_concat_iqtree_pep}/concat.pep.iqtree.nwk"
file_concat_iqtree_dna_root="${dir_concat_iqtree_dna}/concat.dna.rooted.nwk"
file_concat_iqtree_pep_root="${dir_concat_iqtree_pep}/concat.pep.rooted.nwk"
file_constrained_tree="${dir_species_tree}/constrained_tree/constrained.nwk"
file_plot_constrained_tree="${dir_species_tree}/constrained_tree/constrained_tree_constraints.pdf"
file_iq2mc_prefix="${dir_species_tree}/mcmctree_parameter_estimation/iq2mc"
file_iq2mc_ctl="${file_iq2mc_prefix}.mcmctree.ctl"
file_iq2mc_hessian="${file_iq2mc_prefix}.mcmctree.hessian"
file_iq2mc_rooted_tree="${file_iq2mc_prefix}.rooted.nwk"
file_iq2mc_dummy_phy="${file_iq2mc_prefix}.dummy.phy"
file_mcmctree_raw_output="${dir_mcmctree2}/iq2mc.mcmctree.out"
file_mcmctree_figtree_tre="${dir_mcmctree2}/FigTree.tre"
file_mcmctree_dated_nwk="${dir_mcmctree2}/dated_species_tree.nwk"
file_plot_species_trees="${dir_species_tree_summary}/species_trees.pdf"
file_species_tree_busco_cds_pdf="${dir_species_tree_summary}/busco_cds.pdf"
file_species_tree_busco_cds_svg="${dir_species_tree_summary}/busco_cds.svg"
file_species_tree_busco_summary="${dir_species_tree_summary}/annotation_summary.tsv"
file_plot_mcmctree_pdf="${file_mcmctree_dated_nwk%.*}.pdf"
file_dated_species_tree="${dir_species_tree_summary}/dated_species_tree.nwk"
file_dated_species_tree_pdf="${dir_species_tree_summary}/dated_species_tree.pdf"
file_undated_species_tree="${dir_species_tree_summary}/undated_species_tree.nwk"

# Orthogroup
file_orthofinder_done_marker="${dir_orthofinder_hog2og}/README.txt"
file_orthofinder_core_candidates="${dir_orthofinder}/orthofinder_core_species.candidates.tsv"
file_orthofinder_core_selected="${dir_orthofinder}/orthofinder_core_species.selected.tsv"
file_orthofinder_core_selected_list="${dir_orthofinder}/orthofinder_core_species.selected_files.txt"
file_orthofinder_core_species_tree="${dir_orthofinder}/species_tree_core.nwk"
file_orthogroup_selection="${dir_orthofinder_filtered}/Orthogroups.selected.tsv"
file_orthogroup_method_comparison="${dir_orthofinder}/orthogroup_method_comparison/orthogroup_method_comparison.pdf"
file_single_copy_ortholog_decay_plot="${dir_orthogroup_decay}/single_copy_ortholog_decay_plot.pdf"
file_single_copy_ortholog_decay_summary="${dir_orthogroup_decay}/single_copy_ortholog_decay_summary.tsv"

# Genome evolution
file_orthogroup_genecount_selected="${dir_orthofinder_filtered}/Orthogroups.GeneCount.selected.tsv"
file_genome_busco_summary_table="${dir_genome_evolution}/busco_summary_table/busco_summary.tsv"
file_species_omark_summary_table="${dir_genome_evolution}/omark_summary_table/omark_summary.tsv"
file_busco_grampa_dna="${dir_genome_evolution}/grampa_busco_dna/grampa_summary.tsv"
file_busco_grampa_pep="${dir_genome_evolution}/grampa_busco_pep/grampa_summary.tsv"
file_orthogroup_grampa="${dir_genome_evolution}/grampa_orthogroup/grampa_summary.tsv"
file_gene_id="${dir_orthofinder_filtered}/Orthogroups.selected.tsv"
file_cafe_summary_all_pdf="${dir_cafe}/summary_plot/summary_all.pdf"
file_cafe_summary_significant_pdf="${dir_cafe}/summary_plot/summary_significant.pdf"
file_orthogroup_copy_number="${dir_orthogroup_copy_number}/orthogroup_copy_number.tsv"
file_orthogroup_copy_number_matrix="${dir_orthogroup_copy_number_trait_pgls}/orthogroup_copy_number_matrix.tsv"
file_orthogroup_copy_number_trait_pgls="${dir_orthogroup_copy_number_trait_pgls}/orthogroup_copy_number_trait_pgls.tsv"
file_orthogroup_copy_number_trait_pgls_significant="${dir_orthogroup_copy_number_trait_pgls}/orthogroup_copy_number_trait_pgls.significant.tsv"
file_orthogroup_copy_number_trait_pgls_summary_pdf="${dir_orthogroup_copy_number_trait_pgls}/orthogroup_copy_number_trait_pgls.summary.pdf"
file_go_enrichment_significant="${dir_cafe}/go_enrichment/enrichment_significant_${change_direction_go}_${target_branch_go}_significant_go.tsv"

# Runtime helpers
shared_species_busco_stage_done=0
shared_busco_summary_stage_done=0
shared_species_omark_stage_done=0
shared_omark_summary_stage_done=0



















trap cleanup_species_protein_tmp EXIT

# Runtime setup
if [[ "${input_sequence_mode}" == "protein" ]]; then
  if ! species_protein_input_has_files && ! species_cds_input_has_files; then
    echo "protein mode requires either species_protein inputs or species_cds inputs."
    echo "Checked: ${dir_sp_protein_input} and ${dir_sp_cds}"
    exit 1
  fi
  if [[ "${undated_species_tree}" == "iqtree_dna" || "${undated_species_tree}" == "astral_dna" ]]; then
    echo "protein mode does not support undated_species_tree=${undated_species_tree}."
    echo 'Use iqtree_pep or astral_pep instead.'
    exit 1
  fi
  if [[ ${run_concat_iqtree_dna} -eq 1 || ${run_individual_iqtree_dna} -eq 1 || ${run_astral_dna} -eq 1 ]]; then
    echo "Disabling DNA-only species-tree steps in protein mode: run_concat_iqtree_dna, run_individual_iqtree_dna, run_astral_dna"
    run_concat_iqtree_dna=0
    run_individual_iqtree_dna=0
    run_astral_dna=0
  fi
  if [[ ${run_constrained_tree} -eq 1 || ${run_plot_constrained_tree} -eq 1 || ${run_mcmctree1} -eq 1 || ${run_mcmctree2} -eq 1 || ${run_convert_tree_format} -eq 1 || ${run_plot_mcmctreer} -eq 1 ]]; then
    echo "Disabling CDS-only dating steps in protein mode: run_constrained_tree, run_plot_constrained_tree, run_mcmctree1, run_mcmctree2, run_convert_tree_format, run_plot_mcmctreer"
    run_constrained_tree=0
    run_plot_constrained_tree=0
    run_mcmctree1=0
    run_mcmctree2=0
    run_convert_tree_format=0
    run_plot_mcmctreer=0
  fi
  if [[ ${run_busco_dupaware_iqtree_dna} -eq 1 || ${run_busco_dupaware_notung_root_dna} -eq 1 || ${run_busco_dupaware_root_dna} -eq 1 || ${run_busco_dupaware_grampa_dna} -eq 1 ]]; then
    echo "Disabling DNA-only duplicate-aware BUSCO genome-evolution steps in protein mode: run_busco_dupaware_iqtree_dna, run_busco_dupaware_notung_root_dna, run_busco_dupaware_root_dna, run_busco_dupaware_grampa_dna"
    run_busco_dupaware_iqtree_dna=0
    run_busco_dupaware_notung_root_dna=0
    run_busco_dupaware_root_dna=0
    run_busco_dupaware_grampa_dna=0
  fi
else
  check_species_cds "${gg_workspace_dir}"
  check_if_species_files_unique "${dir_sp_cds}"
fi
shared_protein_input_signature=$(compute_shared_protein_input_signature)
species_tree_requested_for_orthofinder=0
if species_tree_summary_generation_requested; then
  species_tree_requested_for_orthofinder=1
fi
refresh_species_tree_for_shared_protein_input_signature "${shared_protein_input_signature}"
refresh_dir_for_shared_protein_input_signature "${dir_orthofinder}" "orthofinder" "${shared_protein_input_signature}"
refresh_dir_for_shared_protein_input_signature "${dir_genome_evolution}" "genome_evolution" "${shared_protein_input_signature}"
memory_notung=$(gg_memory_fraction_gb "${GG_MEM_TOOL_GB}" 1 "${GG_TASK_CPUS}")
memory_iqtree_parallel=$(gg_memory_fraction_gb "${GG_MEM_TOOL_GB}" 1 "${GG_TASK_CPUS}")
iqtree_full_mem_args=(-mem "${GG_MEM_TOOL_GB}G")
iqtree_parallel_mem_args=(-mem "${memory_iqtree_parallel}G")

ensure_dir "${dir_species_tree_summary}"
ensure_dir "${dir_tmp}"
cd "${dir_tmp}"

enable_all_run_flags_for_debug_mode
orthogroup_annotation_method=$(echo "${orthogroup_annotation_method:-mmseqs2}" | tr '[:upper:]' '[:lower:]')
if [[ "${orthogroup_annotation_method}" != "blastp" && "${orthogroup_annotation_method}" != "mmseqs2" ]]; then
  echo "Invalid orthogroup_annotation_method: ${orthogroup_annotation_method}"
  echo 'orthogroup_annotation_method must be either "blastp" or "mmseqs2". Exiting.'
  exit 1
fi

if ! parse_species_tree_rooting "${species_tree_rooting}" species_tree_rooting_method species_tree_rooting_value; then
  exit 1
fi
echo "Resolved species_tree_rooting method: ${species_tree_rooting_method}"
if [[ -n "${species_tree_rooting_value}" ]]; then
  echo "Resolved species_tree_rooting value: ${species_tree_rooting_value}"
fi
print_effective_genome_evolution_config_summary







mcmctree_time_scale_factor_cache=""












iq2mc_binary=""
if [[ ${run_mcmctree1} -eq 1 || ${run_mcmctree2} -eq 1 ]]; then
  echo "IQ2MC dating will run. Checking tools and calibration species."
  if ! iq2mc_binary=$(resolve_iq2mc_binary); then
    echo "Error: IQ2MC-capable IQ-TREE was not found."
    echo "Please install an IQ-TREE build that supports --dating mcmctree (and --mcmc-bds)."
    exit 1
  fi
  if ! iq2mc_option_supported "${iq2mc_binary}"; then
    echo "Warning: Could not verify --mcmc-bds support in ${iq2mc_binary} by help/strings checks."
    echo "Warning: Proceeding with ${iq2mc_binary}; IQ2MC step 2 will fail if this binary is incompatible."
  fi
  if ! command -v mcmctree > /dev/null 2>&1; then
    echo "Error: mcmctree command was not found."
    echo "Please install the IQ2MC-compatible mcmctree binary."
    exit 1
  fi
  echo "Using IQ2MC binary: ${iq2mc_binary}"
  outgroup_label_list=()
  if [[ "${species_tree_rooting_method}" == "outgroup" ]]; then
    mapfile -t outgroup_label_list < <(printf '%s' "${species_tree_rooting_value}" | tr ',' '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' -e '/^$/d')
    for outgroup_label in "${outgroup_label_list[@]}"; do
      stop_if_species_not_found_in "${dir_sp_cds}" "${outgroup_label}"
    done
  fi
  if [[ ${timetree_constraint} -eq 0 ]]; then
    if [[ ${#mcmctree_divergence_time_constraints[@]} -eq 0 ]]; then
      echo "timetree_constraint=0 requires mcmctree_divergence_time_constraints_str with one or more records."
      echo "Example: Arabidopsis_thaliana,Oryza_sativa,130,-|Arabidopsis_thaliana,Amborella_trichopoda,150,200"
      exit 1
    fi
    mcmctree_params=()
    if ! parse_mcmctree_constraint_record "${mcmctree_divergence_time_constraints[0]}" mcmctree_params; then
      exit 1
    fi
    stop_if_species_not_found_in "${dir_sp_cds}" "${mcmctree_params[0]}"
    stop_if_species_not_found_in "${dir_sp_cds}" "${mcmctree_params[1]}"
  fi
  echo ""
fi


gg_core_execution_stage_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/stages/gg_genome_evolution"
# shellcheck source=/dev/null
source "${gg_core_execution_stage_dir}/01_busco_inputs.sh"
# shellcheck source=/dev/null
source "${gg_core_execution_stage_dir}/02_concatenated_species_tree.sh"
# shellcheck source=/dev/null
source "${gg_core_execution_stage_dir}/03_species_tree_summary_and_dating.sh"
# shellcheck source=/dev/null
source "${gg_core_execution_stage_dir}/04_orthofinder_and_selection.sh"
# shellcheck source=/dev/null
source "${gg_core_execution_stage_dir}/05_duplicate_aware_trees.sh"
# shellcheck source=/dev/null
source "${gg_core_execution_stage_dir}/06_family_evolution.sh"
unset gg_core_execution_stage_dir

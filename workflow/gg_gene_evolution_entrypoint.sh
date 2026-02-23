#!/usr/bin/env bash
set -euo pipefail

# SLURM in NIG supercomputer
#SBATCH -J gg_gene_evolution
#SBATCH -c 4 # Number of CPUs
#SBATCH --mem-per-cpu=8G # RAM per CPU in GB
#SBATCH -t 2976:00:00 # maximum time in d-hh:mm:ss format. NIG supercomputer epyc/medium MaxTime=2976:00:00
#SBATCH --output=gg_gene_evolution_%A_%a.out
#SBATCH --error=gg_gene_evolution_%A_%a.err
#SBATCH -p medium # partition name, cluster environment specific
#SBATCH --chdir=.
#SBATCH -a 1 # Array job, 1-N
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>
##SBATCH -N 1-1 # Number of nodes
##SBATCH -n 1 # Number of tasks

## UGE
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -l s_vmem=8G
#$ -l mem_req=8G
#$ -l epyc
#$ -l d_rt=124:00:00:00
#$ -l s_rt=124:00:00:00
#$ -t 1

## PBS in BIAS at NIBB
##PBS -S /bin/bash
##PBS -l ncpus=4
##PBS -l mem=32G
##PBS -J 1
##PBS -q small
##PBS -V

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# N = Number of rows (excluding the header) in workspace/output/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv for mode_orthogroup=1
# N = Number of files in workspace/input/query_gene for mode_query2family=1

echo "$(date): Starting"

# Change these directories for your custom-made analysis
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dir_pg="${script_dir}/../workspace" # pg input and output directory
dir_script="${script_dir}" # directory where core/gg_*_core.sh and support/ locate
gg_image="${script_dir}/../genegalleon.sif" # path to the singularity image

# Load shared defaults when available.
if [[ -s "${script_dir}/gg_common_params.sh" ]]; then
  # shellcheck disable=SC1091
  source "${script_dir}/gg_common_params.sh"
fi

### Start: Modify this block to tailor your analysis ###

# Mode
mode_orthogroup="${mode_orthogroup:-0}" # Analyze OrthoFinder orthogroups
mode_query2family="${mode_query2family:-1}" # Analyze all homologs of input genelist in ${dir_pg}/input/query_gene

# Workflow
run_get_query_fasta="${run_get_query_fasta:-1}" # Activated if mode_query2family=1. Generate amino acid fasta file for query BLAST.
run_query_blast="${run_query_blast:-1}" # Activated if mode_query2family=1.
run_get_fasta="${run_get_fasta:-1}" # Generate in-frame CDS fasta file.
run_rps_blast="${run_rps_blast:-1}" # RPS-BLAST protein domain search.
run_uniprot_annotation="${run_uniprot_annotation:-0}" # DIAMOND-based annotation against UniProt Swiss-Prot.
run_mafft="${run_mafft:-1}" # In-frame nucleotide alignment using MAFFT.
run_amas_original="${run_amas_original:-1}" # Alignment statistics before MaxAlign and TrimAl using AMAS.
run_maxalign="${run_maxalign:-0}" # Remove anomalous sequences by cdskit maxalign.
run_trimal="${run_trimal:-0}" # Remove less-alignable codon sites.
run_clipkit="${run_clipkit:-1}" # Remove less-alignable codon sites.
run_amas_cleaned="${run_amas_cleaned:-1}" # Alignment statistics after MaxAlign and TrimAl using AMAS.
run_iqtree="${run_iqtree:-1}" # Maximum-likelihood phylogenetic reconstruction.
run_tree_root="${run_tree_root:-1}" # Root gene tree using tree_rooting_method.
tree_rooting_method="${tree_rooting_method:-mad}" # notung|midpoint|mad|md; md is mapped to nwkit method "mv".
run_orthogroup_extraction="${run_orthogroup_extraction:-0}" # Activated if mode_query2family=1.
run_generax="${run_generax:-0}" # GeneRax off by default for local/smoke environments without MPI setup.
run_notung_reconcil="${run_notung_reconcil:-0}" # Run NOTUNG for RADTE.
run_tree_dating="${run_tree_dating:-0}" # Species-tree-guided divergence time estimation with RADTE.
run_get_expression_matrix="${run_get_expression_matrix:-0}" # Generate trait matrix of gene expression level.
run_get_gff_info="${run_get_gff_info:-0}" # Generate trait matrix on introns and gene positions.
run_get_promoter_fasta="${run_get_promoter_fasta:-0}" # Generate a promoter fasta file from reference genomes.
run_fimo="${run_fimo:-0}" # Identify promoter motif sequences.
run_tree_pruning="${run_tree_pruning:-0}" # If 1, discard genes without expression data.
check_pruned="${check_pruned:-0}" # Delete downstream outputs if inconsistent to run_tree_pruning.
run_mapdnds_parameter_estimation="${run_mapdnds_parameter_estimation:-0}" # Parameter estimation for mapdNdS.
run_mapdnds="${run_mapdnds:-0}" # Stochastic substitution mapping to estimate dN/dS by mapdNdS.
run_codeml_two_ratio="${run_codeml_two_ratio:-0}" # Run codeml two-ratio model for dN/dS.
run_hyphy_dnds="${run_hyphy_dnds:-0}" # ML dN/dS estimation by HyPhy FitMG94.bf.
run_hyphy_relax="${run_hyphy_relax:-0}" # Run HyPhy RELAX.
run_hyphy_relax_reversed="${run_hyphy_relax_reversed:-0}" # Run HyPhy RELAX with reversed foreground/background.
run_scm_intron="${run_scm_intron:-0}" # Stochastic character mapping of intron traits.
run_phylogeneticem="${run_phylogeneticem:-0}" # OU modeling of gene expression by PhylogeneticEM.
run_l1ou="${run_l1ou:-0}" # OU modeling of gene expression by l1ou.
run_pgls_species_tree="${run_pgls_species_tree:-0}" # PGLS with species tree.
run_iqtree_anc="${run_iqtree_anc:-0}" # Ancestral state reconstruction required for CSUBST.
run_csubst="${run_csubst:-0}" # Protein convergence analysis with CSUBST.
run_summary="${run_summary:-1}" # Generate summary tables.
run_tree_plot="${run_tree_plot:-1}" # Tree visualization pdf.

# Sequence selection in mode_query2family=1
query_blast_method="${query_blast_method:-diamond}" # diamond|tblastn
query_blast_evalue="${query_blast_evalue:-0.01}" # BLAST E-value threshold.
query_blast_coverage="${query_blast_coverage:-0.25}" # BLAST coverage threshold.
max_num_gene_blast_hit_retrieval="${max_num_gene_blast_hit_retrieval:-5000}" # Maximum number of genes to retrieve.
retain_query_in_maxalign="${retain_query_in_maxalign:-1}" # BOOL.

# Phylogeny reconstruction and reconciliation
iqtree_fast_mode_gt="${iqtree_fast_mode_gt:-2000}" # INTEGER.
generax_model="${generax_model:-GTR+G4}" # GeneRax substitution model.
generax_rec_model="${generax_rec_model:-UndatedDL}" # "UndatedDTL" or "UndatedDL"
radte_max_age="${radte_max_age:-1000}" # Upper limit of estimated divergence time in MY.

# species_expression data (value in input files)
exp_value_type="${exp_value_type:-log2p1}"
pgls_use_phenocov="${pgls_use_phenocov:-0}" # BOOL.

# Promoter cis-element analysis
promoter_bp="${promoter_bp:-2000}" # Promoter length in bp
fimo_qvalue="${fimo_qvalue:-0.05}" # False discovery rate threshold for FIMO motif search
jaspar_file="${jaspar_file:-latest}" # "latest"/"auto" or explicit JASPAR filename in ${dir_jaspardb}

# Ornstein-Uhlenbeck modeling of gene expression evolution
clade_collapse_similarity_method="${clade_collapse_similarity_method:-pearson}"
clade_collapse_similarity_threshold="${clade_collapse_similarity_threshold:-0.99}"
require_internal_node_labels="${require_internal_node_labels:-1}" # BOOL. Require non-empty unique internal node labels in input trees for OU scripts.
l1ou_criterion="${l1ou_criterion:-AICc}" # "pBIC", "mBIC", "BIC", or "AICc"
l1ou_nbootstrap="${l1ou_nbootstrap:-0}" # INTEGER.
l1ou_use_fit_file="${l1ou_use_fit_file:-1}" # BOOL.
l1ou_alpha_upper="${l1ou_alpha_upper:-PhylogeneticEM}"
l1ou_convergence="${l1ou_convergence:-1}" # BOOL.
large_tree_num_gene="${large_tree_num_gene:-1000}" # INTEGER.
large_tree_max_nshift="${large_tree_max_nshift:-10}" # INTEGER.
phylogeneticem_use_fit_file="${phylogeneticem_use_fit_file:-1}" # BOOL.

# CSUBST options
csubst_max_arity="${csubst_max_arity:-10}"
csubst_exhaustive_until="${csubst_exhaustive_until:-1}"
csubst_cutoff_stat="${csubst_cutoff_stat:-OCNany2spe,2.0|omegaCany2spe,5.0}"
csubst_max_combination="${csubst_max_combination:-10000}"
csubst_fg_exclude_wg="${csubst_fg_exclude_wg:-no}"
csubst_fg_stem_only="${csubst_fg_stem_only:-yes}"

# Intron and chromosomal character evolution
intron_gain_rate="${intron_gain_rate:-0.0001}" # FLOAT.
retrotransposition_rate="${retrotransposition_rate:-0.001}" # FLOAT.

# Tree visualization
treevis_event_method="${treevis_event_method:-species_overlap}" # "auto", "generax", or "species_overlap"
treevis_clade_ortholog="${treevis_clade_ortholog:-1}" # BOOL.
treevis_support_value="${treevis_support_value:-support_unrooted}" # "support_unrooted", "dup_conf_score", "no"
treevis_branch_length="${treevis_branch_length:-bl_rooted}" # "bl_dated", "bl_rooted", "mapdnds_omega"
treevis_branch_color="${treevis_branch_color:-l1ou_regime}" # "species", "no", or *_regime
treevis_retrotransposition_delta_intron="${treevis_retrotransposition_delta_intron:--0.5}"
treevis_heatmap_transform="${treevis_heatmap_transform:-no}" # "no", "log2", "log10p1", "log2p1"
treevis_pie_chart_value_transformation="${treevis_pie_chart_value_transformation:-identity}" # identity|delog2|delog2p1|delog10|delog10p1
treevis_max_intergenic_dist="${treevis_max_intergenic_dist:-100000}" # Maximum distance between genes in bp.
treevis_synteny="${treevis_synteny:-1}" # BOOL.
treevis_synteny_window="${treevis_synteny_window:-5}" # INTEGER.
treevis_long_branch_display="${treevis_long_branch_display:-auto}" # "auto" or "no"
treevis_long_branch_ref_quantile="${treevis_long_branch_ref_quantile:-0.95}"
treevis_long_branch_detect_ratio="${treevis_long_branch_detect_ratio:-5}"
treevis_long_branch_cap_ratio="${treevis_long_branch_cap_ratio:-2.5}"
treevis_long_branch_tail_shrink="${treevis_long_branch_tail_shrink:-0.02}"
treevis_long_branch_max_fraction="${treevis_long_branch_max_fraction:-0.1}"

### End: Modify this block to tailor your analysis ###

# Misc
exit_if_running=0 # Exit without main analysis if the same SGE_TASK_ID is already running.
delete_tmp_dir="${delete_tmp_dir:-1}" # After this run, delete tmp directory created for each job. Set 0 when debugging.
delete_preexisting_tmp_dir="${delete_preexisting_tmp_dir:-1}" # Before starting this job, delete tmp directory created by previous run.

# Forward config variables (including external overrides) into container environment.
forward_config_vars_to_container_env() {
  local job_script=$1
  local var_name
  local var_names
  var_names=$(
    sed -n '/^### Start: Modify this block to tailor your analysis ###$/,/^### End: Modify this block to tailor your analysis ###$/p' "${job_script}" \
      | sed -n -E 's/^[[:space:]]*([A-Za-z_][A-Za-z0-9_]*)=.*/\1/p' \
      | sort -u
  )
  while IFS= read -r var_name; do
    if [[ -z "${var_name}" ]]; then
      continue
    fi
    if [[ -n "${!var_name+x}" ]]; then
      export "${var_name}"
      export "SINGULARITYENV_${var_name}=${!var_name}"
      export "APPTAINERENV_${var_name}=${!var_name}"
    fi
  done <<< "${var_names}"

  # gg_debug_mode is read by gg_gene_evolution_core.sh but defined outside the config block.
  if [[ -n "${gg_debug_mode+x}" ]]; then
    export gg_debug_mode
    export "SINGULARITYENV_gg_debug_mode=${gg_debug_mode}"
    export "APPTAINERENV_gg_debug_mode=${gg_debug_mode}"
  fi

  # Cleanup controls are defined outside the config block but consumed by core script.
  if [[ -n "${delete_tmp_dir+x}" ]]; then
    export delete_tmp_dir
    export "SINGULARITYENV_delete_tmp_dir=${delete_tmp_dir}"
    export "APPTAINERENV_delete_tmp_dir=${delete_tmp_dir}"
  fi
  if [[ -n "${delete_preexisting_tmp_dir+x}" ]]; then
    export delete_preexisting_tmp_dir
    export "SINGULARITYENV_delete_preexisting_tmp_dir=${delete_preexisting_tmp_dir}"
    export "APPTAINERENV_delete_preexisting_tmp_dir=${delete_preexisting_tmp_dir}"
  fi
}

forward_config_vars_to_container_env "${BASH_SOURCE[0]}"
unset -f forward_config_vars_to_container_env

source "${dir_script}/support/gg_util.sh" # loading utility functions
gg_scheduler_runtime_prelude
unset_singularity_envs
exit_if_running_qstat
if ! set_singularity_command; then
  exit 1
fi
variable_SGEnizer
set_singularityenv
gg_print_scheduler_runtime_summary

mkdir -p "${dir_pg}"
cd "${dir_pg}"
set +e
${singularity_command} "${gg_image}" < "${dir_script}/core/gg_gene_evolution_core.sh"
cmd_exit_code=$?
set -e
if [[ ${cmd_exit_code} -eq 8 ]]; then
  echo "Output files were detected. No more run of gg_gene_evolution_core.sh is necessary."
elif [[ ${cmd_exit_code} -ne 0 ]]; then
  echo "gg_gene_evolution_core.sh failed with code ${cmd_exit_code}."
  exit "${cmd_exit_code}"
fi
if ! gg_trigger_versions_dump "$(basename "${BASH_SOURCE[0]}")"; then
  echo "Warning: gg_versions trigger failed."
fi
echo "$(date): Ending"

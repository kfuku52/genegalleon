#!/usr/bin/env bash
set -euo pipefail

# SLURM in NIG supercomputer
#SBATCH -J gg_genome_evolution
#SBATCH -c 4 # Number of CPUs
#SBATCH --mem-per-cpu=8G # RAM per CPU in GB
#SBATCH -t 2976:00:00 # maximum time in d-hh:mm:ss format. NIG supercomputer epyc/medium MaxTime=2976:00:00
#SBATCH --output=gg_genome_evolution_%j.out
#SBATCH --error=gg_genome_evolution_%j.err
#SBATCH -p medium # partition name, cluster environment specific
#SBATCH --chdir=.
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>
##SBATCH -N 1-1 # Number of nodes
##SBATCH -n 1 # Number of tasks

#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -l s_vmem=8G
#$ -l mem_req=8G
#$ -l epyc
#$ -l d_rt=124:00:00:00
#$ -l s_rt=124:00:00:00

echo "$(date): Starting"

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

# Species-tree workflow flags
run_species_busco="${run_species_busco:-1}"
run_species_get_busco_summary="${run_species_get_busco_summary:-1}"
run_individual_get_fasta="${run_individual_get_fasta:-1}"
run_individual_mafft="${run_individual_mafft:-1}"
run_individual_trimal="${run_individual_trimal:-1}"
run_concat_alignment="${run_concat_alignment:-1}"
run_concat_iqtree_protein="${run_concat_iqtree_protein:-1}"
run_concat_iqtree_dna="${run_concat_iqtree_dna:-1}"
run_individual_iqtree_pep="${run_individual_iqtree_pep:-1}"
run_astral_pep="${run_astral_pep:-1}"
run_individual_iqtree_dna="${run_individual_iqtree_dna:-1}"
run_astral_dna="${run_astral_dna:-1}"
run_plot_species_trees="${run_plot_species_trees:-1}" # Plotting 4 species trees
run_constrained_tree="${run_constrained_tree:-1}" # Introduce divergence time constraints for IQ2MC input
run_plot_constrained_tree="${run_plot_constrained_tree:-1}" # Plot node-wise constrained ranges in constrained.nwk
run_mcmctree1="${run_mcmctree1:-1}" # IQ2MC step 2
run_mcmctree2="${run_mcmctree2:-1}" # IQ2MC step 3
run_convert_tree_format="${run_convert_tree_format:-1}"
run_plot_mcmctreer="${run_plot_mcmctreer:-1}" # Plot dated species tree

# Orthogroup workflow flags
run_cds_translation="${run_cds_translation:-1}" # Internal helper for temporary protein FASTA generation.
run_orthofinder="${run_orthofinder:-1}" # OrthoFinder run
run_og_selection="${run_og_selection:-1}" # Selecting orthogroups for downstream analyses
run_orthogroup_method_comparison="${run_orthogroup_method_comparison:-1}" # Method comparison plot

# Genome-evolution workflow flags
run_genome_busco="${run_genome_busco:-1}"
run_genome_get_busco_summary="${run_genome_get_busco_summary:-1}"
run_busco_getfasta="${run_busco_getfasta:-1}"
run_busco_mafft="${run_busco_mafft:-1}"
run_busco_trimal="${run_busco_trimal:-1}"
run_busco_iqtree_dna="${run_busco_iqtree_dna:-1}"
run_busco_iqtree_pep="${run_busco_iqtree_pep:-1}"
run_busco_notung_root_dna="${run_busco_notung_root_dna:-1}"
run_busco_notung_root_pep="${run_busco_notung_root_pep:-1}"
run_busco_root_dna="${run_busco_root_dna:-1}"
run_busco_root_pep="${run_busco_root_pep:-1}"
run_busco_grampa_dna="${run_busco_grampa_dna:-1}"
run_busco_grampa_pep="${run_busco_grampa_pep:-1}"
run_orthogroup_grampa="${run_orthogroup_grampa:-1}" # Requires gg_gene_evolution rooted trees
run_cafe="${run_cafe:-0}"
run_go_enrichment="${run_go_enrichment:-0}"

# Shared parameters
strictly_single_copy_only="${strictly_single_copy_only:-0}"
bootstrap_params="${bootstrap_params:--bb 1000 -bnni}"
nucleotide_model="${nucleotide_model:-GTR+R4}"
protein_model="${protein_model:-LG+R4}"
notung_jar="${notung_jar:-/usr/local/bin/Notung.jar}"

# Species-tree parameters
undated_species_tree="${undated_species_tree:-astral_pep}" # {iqtree_dna,iqtree_pep,astral_dna,astral_pep}
astral_min_tips="${astral_min_tips:-4}"
timetree_constraint="${timetree_constraint:-1}"
mcmc_burnin="${mcmc_burnin:-20000}"
mcmc_sampfreq="${mcmc_sampfreq:-100}"
mcmc_nsample="${mcmc_nsample:-20000}"
mcmc_birth_death_sampling="${mcmc_birth_death_sampling:-1,1,0.5}" # birth,death,sampling_fraction
mcmc_clock_model="${mcmc_clock_model:-IND}" # {EQUAL, IND, CORR}

# Orthogroup parameters
orthogroup_table="${orthogroup_table:-HOG}" # "OG" or "HOG"
min_num_gene="${min_num_gene:-4}"
min_num_species="${min_num_species:-2}"
max_orthofinder_core_species="${max_orthofinder_core_species:-50}"
min_percent_species_coverage="${min_percent_species_coverage:-50}"
max_num_gene="${max_num_gene:-1000}"

# Genome-evolution parameters
min_gene_orthogroup_grampa="${min_gene_orthogroup_grampa:-5}"
max_gene_orthogroup_grampa="${max_gene_orthogroup_grampa:-50}"
max_size_differential_cafe="${max_size_differential_cafe:-9999999}"
n_gamma_cats_cafe="${n_gamma_cats_cafe:-4}"
change_direction_go="${change_direction_go:-increase}" # "increase" or "decrease"
go_category="${go_category:-BP,MF,CC}" # BP, MF, CC
delete_tmp_dir="${delete_tmp_dir:-1}" # After normal completion, delete tmp directories. Set 0 when debugging.

### End: Modify this block to tailor your analysis ###

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

  # gg_debug_mode is read by gg_*_core.sh scripts but defined outside the config block.
  if [[ -n "${gg_debug_mode+x}" ]]; then
    export gg_debug_mode
    export "SINGULARITYENV_gg_debug_mode=${gg_debug_mode}"
    export "APPTAINERENV_gg_debug_mode=${gg_debug_mode}"
  fi
}

forward_config_vars_to_container_env "${BASH_SOURCE[0]}"
unset -f forward_config_vars_to_container_env

source "${dir_script}/support/gg_util.sh" # loading utility functions
gg_scheduler_runtime_prelude
unset_singularity_envs
if ! set_singularity_command; then
  exit 1
fi
variable_SGEnizer
set_singularityenv
gg_print_scheduler_runtime_summary

mkdir -p "${dir_pg}"
cd "${dir_pg}"
${singularity_command} "${gg_image}" < "${dir_script}/core/gg_genome_evolution_core.sh"
if ! gg_trigger_versions_dump "$(basename "${BASH_SOURCE[0]}")"; then
  echo "Warning: gg_versions trigger failed."
fi

echo "$(date): Ending"

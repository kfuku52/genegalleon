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
run_species_busco=1
run_species_get_busco_summary=1
run_individual_get_fasta=1
run_individual_mafft=1
run_individual_trimal=1
run_concat_alignment=1
run_concat_iqtree_protein=1
run_concat_iqtree_dna=1
run_individual_iqtree_pep=1
run_astral_pep=1
run_individual_iqtree_dna=1
run_astral_dna=1
run_plot_species_trees=1 # Plotting 4 species trees
run_constrained_tree=1 # Introduce divergence time constraints for IQ2MC input
run_plot_constrained_tree=1 # Plot node-wise constrained ranges in constrained.nwk
run_mcmctree1=1 # IQ2MC step 2
run_mcmctree2=1 # IQ2MC step 3
run_convert_tree_format=1
run_plot_mcmctreer=1 # Plot dated species tree

# Orthogroup workflow flags
run_cds_translation=1 # Internal helper for temporary protein FASTA generation.
run_orthofinder=1 # OrthoFinder run
run_og_selection=1 # Selecting orthogroups for downstream analyses
run_orthogroup_method_comparison=1 # Method comparison plot

# Genome-evolution workflow flags
run_genome_busco=1
run_genome_get_busco_summary=1
run_busco_getfasta=1
run_busco_mafft=1
run_busco_trimal=1
run_busco_iqtree_dna=1
run_busco_iqtree_pep=1
run_busco_notung_root_dna=1
run_busco_notung_root_pep=1
run_busco_root_dna=1
run_busco_root_pep=1
run_busco_grampa_dna=1
run_busco_grampa_pep=1
run_orthogroup_grampa=1 # Requires gg_gene_evolution rooted trees
run_cafe=0
run_go_enrichment=0

# Shared parameters
strictly_single_copy_only=0
bootstrap_params="-bb 1000 -bnni"
nucleotide_model="GTR+R4"
protein_model="LG+R4"
notung_jar="/usr/local/bin/Notung.jar"

# Species-tree parameters
undated_species_tree="astral_pep" # {iqtree_dna,iqtree_pep,astral_dna,astral_pep}
astral_min_tips=4
timetree_constraint=1
mcmc_burnin=20000
mcmc_sampfreq=100
mcmc_nsample=20000
mcmc_birth_death_sampling="1,1,0.5" # birth,death,sampling_fraction
mcmc_clock_model="IND" # {EQUAL, IND, CORR}

# Orthogroup parameters
orthogroup_table="HOG" # "OG" or "HOG"
orthogroup_annotation_method="mmseqs2" # blastp|mmseqs2 for representative-gene UniProt annotation in orthogroup selection.
min_num_gene=4
min_num_species=2
max_orthofinder_core_species=50
min_percent_species_coverage=50
max_num_gene=1000

# Genome-evolution parameters
min_gene_orthogroup_grampa=5
max_gene_orthogroup_grampa=50
max_size_differential_cafe=9999999
n_gamma_cats_cafe=4
change_direction_go="increase" # "increase" or "decrease"
go_category="BP,MF,CC" # BP, MF, CC
delete_tmp_dir=1 # After normal completion, delete tmp directories. Set 0 when debugging.

### End: Modify this block to tailor your analysis ###

source "${dir_script}/support/gg_util.sh" # loading utility functions
# Forward config variables (including external overrides) into container environment.
forward_config_vars_to_container_env "${BASH_SOURCE[0]}"
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

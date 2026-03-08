#!/usr/bin/env bash

# SLURM in NIG supercomputer
#SBATCH -J gg_genome_evolution
#SBATCH -c 4
#SBATCH --mem-per-cpu=8G
#SBATCH -t 2976:00:00
#SBATCH --output=gg_genome_evolution_%j.out
#SBATCH --error=gg_genome_evolution_%j.err
#SBATCH -p epyc
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

set -euo pipefail

echo "$(date): Starting"

# Resolve workflow paths for local and scheduler-spooled execution.
gg_bootstrap_submit_dir="${SLURM_SUBMIT_DIR:-${PBS_O_WORKDIR:-${PWD:-}}}"
gg_bootstrap_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for bootstrap_path in \
  "${gg_bootstrap_submit_dir}/support/gg_entrypoint_bootstrap.sh" \
  "${gg_bootstrap_submit_dir}/workflow/support/gg_entrypoint_bootstrap.sh" \
  "${gg_bootstrap_script_dir}/support/gg_entrypoint_bootstrap.sh"
do
  if [[ -s "${bootstrap_path}" ]]; then
    # shellcheck disable=SC1090
    source "${bootstrap_path}"
    break
  fi
done
unset gg_bootstrap_submit_dir gg_bootstrap_script_dir bootstrap_path
if ! declare -F gg_entrypoint_initialize >/dev/null 2>&1; then
  echo "Failed to locate gg_entrypoint_bootstrap.sh from BASH_SOURCE[0]=${BASH_SOURCE[0]}" >&2
  exit 1
fi
if ! gg_entrypoint_initialize "${BASH_SOURCE[0]}" 1; then
  exit 1
fi
gg_entrypoint_name="gg_genome_evolution_entrypoint.sh"

### Start: Modify this block to tailor your analysis ###

# Species-tree workflow flags
run_species_busco=1 # Run BUSCO on species inputs for species-tree marker discovery.
run_species_get_busco_summary=1 # Summarize BUSCO results across species-tree inputs.
run_individual_get_fasta=1 # Extract per-ortholog FASTA files for species-tree reconstruction.
run_individual_mafft=1 # Align individual species-tree ortholog sets with MAFFT.
run_individual_trimal=1 # Trim individual species-tree ortholog alignments.
run_concat_alignment=1 # Concatenate selected species-tree ortholog alignments.
run_concat_iqtree_protein=1 # Build a concatenated protein species tree with IQ-TREE.
run_concat_iqtree_dna=1 # Build a concatenated DNA species tree with IQ-TREE.
run_individual_iqtree_pep=1 # Build per-ortholog protein trees for ASTRAL input.
run_astral_pep=1 # Infer a protein species tree with ASTRAL.
run_individual_iqtree_dna=1 # Build per-ortholog DNA trees for ASTRAL input.
run_astral_dna=1 # Infer a DNA species tree with ASTRAL.
run_plot_species_trees=1 # Plotting 4 species trees
run_constrained_tree=1 # Introduce divergence time constraints for IQ2MC input
run_plot_constrained_tree=1 # Plot node-wise constrained ranges in constrained.nwk
run_mcmctree1=1 # IQ2MC step 2
run_mcmctree2=1 # IQ2MC step 3
run_convert_tree_format=1 # Convert dated species trees into downstream exchange formats.
run_plot_mcmctreer=1 # Plot dated species tree

# Orthogroup workflow flags
run_cds_translation=1 # Internal helper for temporary protein FASTA generation.
run_orthofinder=1 # OrthoFinder run
run_og_selection=1 # Selecting orthogroups for downstream analyses
run_orthogroup_method_comparison=1 # Method comparison plot

# Genome-evolution workflow flags
run_genome_busco=1 # Run BUSCO on genome-derived CDS sets for genome-evolution analyses.
run_genome_get_busco_summary=1 # Summarize BUSCO results across genome-derived CDS sets.
run_busco_getfasta=1 # Extract BUSCO ortholog FASTA files from genome annotations.
run_busco_mafft=1 # Align BUSCO ortholog FASTA files with MAFFT.
run_busco_trimal=1 # Trim BUSCO ortholog alignments.
run_busco_iqtree_dna=1 # Build BUSCO DNA gene trees with IQ-TREE.
run_busco_iqtree_pep=1 # Build BUSCO protein gene trees with IQ-TREE.
run_busco_notung_root_dna=1 # Root BUSCO DNA trees with NOTUNG.
run_busco_notung_root_pep=1 # Root BUSCO protein trees with NOTUNG.
run_busco_root_dna=1 # Run standard rooting on BUSCO DNA trees.
run_busco_root_pep=1 # Run standard rooting on BUSCO protein trees.
run_busco_grampa_dna=1 # Run GRAMPA on rooted BUSCO DNA trees.
run_busco_grampa_pep=1 # Run GRAMPA on rooted BUSCO protein trees.
run_orthogroup_grampa=1 # Requires gg_gene_evolution rooted trees
run_cafe=0 # Run CAFE family-size evolution analysis.
run_go_enrichment=0 # Run GO enrichment on selected orthogroups or shifts.

# Shared parameters
strictly_single_copy_only=0 # Restrict marker selection to strictly single-copy orthologs only.
bootstrap_params="-bb 1000 -bnni" # Extra IQ-TREE bootstrap parameters.
nucleotide_model="GTR+R4" # IQ-TREE nucleotide substitution model.
protein_model="LG+R4" # IQ-TREE protein substitution model.
notung_jar="/usr/local/bin/Notung.jar" # Path to the Notung JAR used for rooting.

# Species-tree parameters
undated_species_tree="astral_pep" # {iqtree_dna,iqtree_pep,astral_dna,astral_pep}
species_tree_rooting="taxonomy" # taxonomy[,ncbi[,opentree,timetree...]] | outgroup,GENUS_SPECIES[,GENUS_SPECIES...] | midpoint | mad | mv
astral_min_tips=4 # Minimum tip count required for per-gene trees used by ASTRAL.
timetree_constraint=1 # Use TimeTree confidence intervals for species-tree dating when set to 1.
mcmctree_divergence_time_constraints_str="" # Used only when timetree_constraint=0. Example: "Arabidopsis_thaliana,Oryza_sativa,130,-|Arabidopsis_thaliana,Amborella_trichopoda,150,200"
mcmc_burnin=20000 # Burn-in iterations for MCMCTree.
mcmc_sampfreq=100 # Sampling frequency for MCMCTree.
mcmc_nsample=20000 # Number of posterior samples retained by MCMCTree.
mcmc_birth_death_sampling="1,1,0.5" # birth,death,sampling_fraction
mcmc_clock_model="IND" # {EQUAL, IND, CORR}

# Orthogroup parameters
orthogroup_table="HOG" # "OG" or "HOG"
orthogroup_annotation_method="mmseqs2" # blastp|mmseqs2 for representative-gene UniProt annotation in orthogroup selection.
min_num_gene=4 # Minimum total gene count required for an orthogroup.
min_num_species=2 # Minimum number of species required for an orthogroup.
max_orthofinder_core_species=50 # Maximum number of species retained in the core OrthoFinder set.
min_percent_species_coverage=50 # Minimum percent species coverage required for orthogroup selection.
max_num_gene=1000 # Maximum total gene count allowed for an orthogroup.

# Genome-evolution parameters
min_gene_orthogroup_grampa=5 # Minimum gene count required for GRAMPA-ready orthogroups.
max_gene_orthogroup_grampa=50 # Maximum gene count allowed for GRAMPA-ready orthogroups.
grampa_h1="" # Optional GRAMPA H1 hypothesis. Leave empty to skip GRAMPA steps. Example: "2" or "x,y,z".
max_size_differential_cafe=9999999 # Maximum family-size difference modeled by CAFE.
n_gamma_cats_cafe=4 # Number of gamma categories used by CAFE.
target_branch_go="" # Optional GO-enrichment target branch. Leave empty to skip GO enrichment. Example: "<1>" or "Arabidopsis_thaliana".
change_direction_go="increase" # "increase" or "decrease"
go_category="BP,MF,CC" # BP, MF, CC
delete_tmp_dir=1 # After normal completion, delete tmp directories. Set 0 when debugging.

### End: Modify this block to tailor your analysis ###

source "${gg_support_dir}/gg_util.sh" # loading utility functions
# Forward config variables (including external overrides) into container environment.
forward_config_vars_to_container_env "${gg_entrypoint_name}"
if ! gg_entrypoint_prepare_container_runtime 0; then
  exit 1
fi
gg_entrypoint_activate_container_runtime

gg_entrypoint_enter_workspace
${singularity_command} "${gg_container_image_path}" < "${gg_core_dir}/gg_genome_evolution_core.sh"
gg_require_versions_dump "${gg_entrypoint_name}"

echo "$(date): Ending"

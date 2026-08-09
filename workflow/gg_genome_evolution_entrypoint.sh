#!/usr/bin/env bash

# Scheduler header notes:
# - Keep sections ordered as SLURM -> UGE -> PBS across entrypoints.
# - Update job name, CPU count, memory, walltime, log paths, and array size together.
# - UGE resources default to SHIROKANE AGE; SLURM/PBS site-specific lines remain examples.

# SLURM
# Common parameters: job name, cores per task, memory per core, walltime, log files, and working directory.
#SBATCH -J gg_genome_evolution
#SBATCH -c 4
#SBATCH --mem-per-cpu=8G
#SBATCH -t 14-00:00:00
#SBATCH --output=gg_genome_evolution_entrypoint.sh_%j.out
#SBATCH --error=gg_genome_evolution_entrypoint.sh_%j.err
#SBATCH --chdir=.
#SBATCH --ignore-pbs
# Site-specific partition example.
#SBATCH -p epyc,rome,medium
# Optional notifications and single-node examples.
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>

## UGE
# SHIROKANE AGE defaults: shell, working directory, slot count, memory per slot, and ljob.
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -l s_vmem=8G
#$ -l ljob

## PBS
# Common parameters: shell, CPU count, total memory, and exported environment.
#PBS -S /bin/bash
#PBS -l ncpus=4
#PBS -l mem=32G
# Site-specific queue example.
#PBS -q small
#PBS -V

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# Fixed to 1

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

# Input-preparation workflow flags
run_cds_translation=1 # Internal helper for temporary protein FASTA generation when protein mode falls back to species_cds.

# Species-tree workflow flags
run_species_busco=1 # Run the shared multi-species BUSCO stage used for species-tree marker discovery and BUSCO-based genome-evolution analyses.
run_build_species_busco_summary=1 # Build the shared BUSCO summary table used by species-tree and BUSCO-based genome-evolution stages.
run_extract_species_tree_fasta=1 # Extract per-ortholog FASTA files for species-tree reconstruction.
run_individual_mafft=1 # Align individual species-tree ortholog sets with MAFFT.
run_individual_trimal=1 # Trim individual species-tree ortholog alignments.
run_concat_alignment=1 # Concatenate selected species-tree ortholog alignments.
run_concat_iqtree_protein=1 # Build a concatenated protein species tree with IQ-TREE.
run_concat_iqtree_dna=1 # Build a concatenated DNA species tree with IQ-TREE.
run_individual_iqtree_pep=1 # Build per-ortholog protein trees for ASTRAL input.
run_individual_iqtree_dna=1 # Build per-ortholog DNA trees for ASTRAL input.
run_astral_pep=1 # Infer a protein species tree with ASTRAL.
run_astral_dna=1 # Infer a DNA species tree with ASTRAL.
run_plot_species_trees=1 # Plot the four candidate species trees: concatenated protein, concatenated DNA, ASTRAL protein, and ASTRAL DNA.
run_constrained_tree=1 # Add divergence-time calibration constraints to the selected undated species tree for IQ2MC/MCMCTree.
run_plot_constrained_tree=1 # Plot node-wise divergence-time constraint ranges in constrained.nwk.
run_mcmctree1=1 # Run IQ2MC parameter-estimation/MCMCTree preparation step.
run_mcmctree2=1 # Run the main MCMCTree dating step and export dated species-tree outputs.
run_convert_tree_format=1 # Convert dated species trees into downstream exchange formats.
run_plot_mcmctreer=1 # Plot the MCMCTree dated species tree as a publication-style summary figure.

# Orthogroup and species-protein QC workflow flags
run_orthofinder=1 # Run OrthoFinder to infer orthogroups from the configured species CDS/protein inputs.
run_species_omark=0 # Run OMArk proteome quality assessment for each species using shared protein inputs.
run_build_species_omark_summary=1 # Build the shared OMArk summary table for species-wise proteome quality assessment.
run_og_selection=1 # Select orthogroups that pass gene-count, alignment, and annotation filters for downstream analyses.
run_orthogroup_method_comparison=1 # Plot comparison among orthogroup/species-tree inference methods.
run_single_copy_ortholog_decay_plot=1 # Line plot with SD of orthogroup retention when randomly subsampling species.

# Genome-evolution workflow flags
run_busco_dupaware_extract_fasta=0 # Extract duplicate-aware BUSCO ortholog FASTA files from genome annotations.
run_busco_dupaware_mafft=0 # Align duplicate-aware BUSCO ortholog FASTA files with MAFFT.
run_busco_dupaware_trimal=0 # Trim duplicate-aware BUSCO ortholog alignments.
run_busco_dupaware_iqtree_dna=0 # Build duplicate-aware BUSCO DNA gene trees with IQ-TREE.
run_busco_dupaware_iqtree_pep=0 # Build duplicate-aware BUSCO protein gene trees with IQ-TREE.
run_busco_dupaware_notung_root_dna=0 # Root duplicate-aware BUSCO DNA trees with NOTUNG.
run_busco_dupaware_notung_root_pep=0 # Root duplicate-aware BUSCO protein trees with NOTUNG.
run_busco_dupaware_root_dna=0 # Run standard rooting on duplicate-aware BUSCO DNA trees.
run_busco_dupaware_root_pep=0 # Run standard rooting on duplicate-aware BUSCO protein trees.
run_busco_dupaware_grampa_dna=0 # Run GRAMPA on rooted duplicate-aware BUSCO DNA trees.
run_busco_dupaware_grampa_pep=0 # Run GRAMPA on rooted duplicate-aware BUSCO protein trees.
run_orthogroup_grampa=1 # Run GRAMPA on orthogroup gene trees; requires rooted trees produced by gg_gene_evolution.
run_cafe=0 # Run CAFE family-size evolution analysis on selected orthogroup gene-count tables and the dated species tree.
run_orthogroup_copy_number_trait_pgls=0 # Test associations between orthogroup copy numbers and species traits with species-tree PGLS.
run_go_enrichment=0 # Run GO enrichment for branches or orthogroups selected by family-size change tests.

# Shared parameters
input_sequence_mode="${input_sequence_mode:-${GG_COMMON_INPUT_SEQUENCE_MODE:-cds}}" # {cds,protein}; protein mode uses species_protein inputs or per-species CDS->protein translation with optional species_genetic_code/species_genetic_code.tsv overrides.
strictly_single_copy_only=0 # Restrict marker selection to strictly single-copy orthologs only.
bootstrap_params="-bb 1000 -bnni" # Extra IQ-TREE bootstrap parameters.
nucleotide_model="GTR+R4" # IQ-TREE nucleotide substitution model.
protein_model="LG+R4" # IQ-TREE protein substitution model.
notung_jar="/usr/local/bin/Notung.jar" # Path to the Notung JAR used for rooting.

# Species-tree parameters
species_tree_output_storage="${species_tree_output_storage:-${GG_COMMON_SPECIES_TREE_OUTPUT_STORAGE:-zip}}" # zip|files|raw; ZIP mode archives high-file-count single-copy stage directories between uses.
species_tree_zip_compression="${species_tree_zip_compression:-${GG_COMMON_SPECIES_TREE_ZIP_COMPRESSION:-adaptive}}" # adaptive|deflate|store.
species_tree_zip_compression_level="${species_tree_zip_compression_level:-${GG_COMMON_SPECIES_TREE_ZIP_COMPRESSION_LEVEL:-6}}" # Deflate level 0-9.
undated_species_tree="astral_pep" # {iqtree_dna,iqtree_pep,astral_dna,astral_pep}; species-tree source copied to undated_species_tree.nwk for dating and downstream summaries.
species_tree_rooting="taxonomy" # taxonomy[,ncbi[,opentree,timetree...]] | outgroup,GENUS_SPECIES[,GENUS_SPECIES...] | midpoint | mad | mv; selects how species trees are rooted before dating, using taxonomy providers, explicit outgroups, or topology/branch-length rooting methods.
species_busco_parallel_jobs="auto" # auto uses up to four concurrent species BUSCO jobs within GG_TASK_CPUS; set a positive integer to cap memory use.
species_busco_memory_gb_per_job=4 # Minimum tool-memory budget per concurrent BUSCO species job; parallelism is capped by GG_MEM_TOOL_GB / this value.
astral_min_tips=4 # Minimum tip count required for per-gene trees used by ASTRAL.
timetree_constraint=1 # Use TimeTree confidence intervals for species-tree dating when set to 1.
mcmctree_divergence_time_constraints_str="" # Used only when timetree_constraint=0. Example: "Arabidopsis_thaliana,Oryza_sativa,130,-|Arabidopsis_thaliana,Amborella_trichopoda,150,200"
mcmc_burnin=20000 # Burn-in iterations for MCMCTree.
mcmc_sampfreq=100 # Sampling frequency for MCMCTree.
mcmc_nsample=20000 # Number of posterior samples retained by MCMCTree.
mcmc_birth_death_sampling="1,1,0.5" # MCMCTree birth-death-sampling prior as birth,death,sampling_fraction.
mcmc_clock_model="IND" # {EQUAL, IND, CORR}; MCMCTree clock model, where EQUAL is strict-clock, IND uses independent rates, and CORR uses autocorrelated rates.

# Orthogroup and species-protein QC parameters
omark_db_path="auto" # Path to the OMArk OMAmer LUCA.h5 database, or "auto" to download it under workspace/downloads/omark/.
orthogroup_table="HOG" # "OG" or "HOG"; choose whether orthogroup selection starts from OrthoFinder orthogroups or hierarchical orthogroups.
orthogroup_annotation_method="mmseqs2" # blastp|mmseqs2 for representative-gene UniProt annotation in orthogroup selection.
min_num_gene=4 # Minimum total gene count required for an orthogroup.
min_num_species=2 # Minimum number of species required for an orthogroup.
max_orthofinder_core_species=50 # Maximum number of species retained in the core OrthoFinder set.
orthofinder_core_filters="busco_complete_pct:ge:80,num_seq:le:100000" # Comma-separated nwkit sample filters for two-round OrthoFinder core selection.
orthofinder_core_rank="num_seq:asc,busco_complete_pct:desc" # Comma-separated nwkit sample rank keys used after filtering.
orthofinder_core_method="max-pd" # nwkit sample method for tree-aware OrthoFinder core selection.
orthofinder_algorithm_threads="auto" # auto uses GG_TASK_CPUS for OrthoFinder analysis phases; set a positive integer to cap them.
min_percent_species_coverage=50 # Minimum percent species coverage required for orthogroup selection.
max_num_gene=1000 # Maximum total gene count allowed for an orthogroup.
orthogroup_decay_replicates=1000 # Random species subsets per species-count value for single-copy ortholog decay plotting.
orthogroup_decay_species_counts="auto" # Species-count values for single-copy ortholog decay plotting: auto, comma list, or ranges such as 1-50:5.
orthogroup_decay_seed=1 # Random seed for single-copy ortholog decay plotting.

# Genome-evolution parameters
min_gene_orthogroup_grampa=5 # Minimum gene count required for GRAMPA-ready orthogroups.
max_gene_orthogroup_grampa=50 # Maximum gene count allowed for GRAMPA-ready orthogroups.
grampa_h1="" # Optional GRAMPA H1 hypothesis. Leave empty to skip GRAMPA steps. Example: "2" or "x,y,z".
orthogroup_copy_number_max_size_differential=9999999 # Maximum family-size difference retained in the shared orthogroup copy-number matrix.
n_gamma_cats_cafe=4 # Number of gamma categories used by CAFE.
orthogroup_copy_number_trait="all" # Trait column name(s) in species_trait.tsv to test against orthogroup copy numbers, or "all".
orthogroup_copy_number_trait_min_species=4 # Minimum number of tree-matched species required for each orthogroup copy-number trait PGLS fit.
orthogroup_copy_number_trait_family_ids="" # Optional comma/space-separated orthogroup IDs to test; empty means use max_families.
orthogroup_copy_number_trait_family_file="" # Optional file listing orthogroup IDs to test.
orthogroup_copy_number_trait_max_families="all" # Maximum orthogroups tested: all|auto|0 for unlimited, or a non-negative integer.
orthogroup_copy_number_trait_p_adjust_method="BH" # P-value adjustment method passed to p.adjust for orthogroup copy-number trait PGLS.
orthogroup_copy_number_trait_alpha=0.05 # Adjusted P-value cutoff used for orthogroup_copy_number_trait_pgls.significant.tsv and summary plot guide line.
orthogroup_copy_number_trait_plot_top_n=50 # Number of strongest orthogroup copy-number trait associations shown in the summary plot.
file_trait="auto" # Species trait table path for orthogroup copy-number trait PGLS, or auto for workspace/input/species_trait/species_trait.tsv.
target_branch_go="" # Optional GO-enrichment target branch. Leave empty to skip GO enrichment. Example: "<1>" or "Arabidopsis_thaliana".
change_direction_go="increase" # "increase" or "decrease"; family-size direction tested for GO enrichment on target_branch_go.
go_category="BP,MF,CC" # GO aspects included in enrichment: BP biological process, MF molecular function, CC cellular component.

### End: Modify this block to tailor your analysis ###

# Misc
delete_tmp_dir=1 # After normal completion, delete tmp directories. Set 0 when debugging.

source "${gg_support_dir}/gg_util.sh" # loading utility functions
# Forward config variables (including external overrides) into container environment.
gg_apply_registered_env_overrides "${gg_entrypoint_name}"
forward_config_vars_to_container_env "${gg_entrypoint_name}"
if ! gg_entrypoint_prepare_container_runtime 0; then
  exit 1
fi
gg_entrypoint_activate_container_runtime

gg_entrypoint_enter_workspace
gg_runtime_core_script="$(gg_prepare_entrypoint_runtime_snapshot "${gg_entrypoint_name}" "${gg_core_dir}/gg_genome_evolution_core.sh")"
gg_run_container_shell_script "${gg_container_image_path}" "${gg_runtime_core_script}"
gg_require_versions_dump "${gg_entrypoint_name}"

echo "$(date): Ending"

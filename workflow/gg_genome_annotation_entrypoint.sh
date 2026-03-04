#!/usr/bin/env bash
set -euo pipefail

# SLURM in NIG supercomputer
#SBATCH -J gg_genome_annotation
#SBATCH -c 4 # Number of CPUs
#SBATCH --mem-per-cpu=8G # RAM per CPU in GB
#SBATCH -t 2976:00:00 # maximum time in d-hh:mm:ss format. NIG supercomputer epyc/medium MaxTime=2976:00:00
#SBATCH --output=gg_genome_annotation_%A_%a.out
#SBATCH --error=gg_genome_annotation_%A_%a.err
#SBATCH -p medium # partition name, cluster environment specific
#SBATCH --chdir=.
#SBATCH -a 1 # Array job, 1-N
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
#$ -t 1

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# N = Number of fasta files in workspace/input/species_cds

echo "$(date): Starting"

# Change these directories for your custom-made analysis
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dir_pg="${script_dir}/../workspace" # pg input and output directory
dir_script="${script_dir}" # directory where gg_util.sh and core/gg_*_core.sh locate
gg_image="${script_dir}/../genegalleon.sif" # path to the singularity image

# Load shared defaults when available.
if [[ -s "${script_dir}/gg_common_params.sh" ]]; then
  # shellcheck disable=SC1091
  source "${script_dir}/gg_common_params.sh"
fi

### Start: Modify this block to tailor your analysis ###

# CDS analysis
run_get_gff_info=0 # Collect gene information from workspace/input/species_gff
run_busco_cds=0 # Gene set completeness analysis
run_uniprot_annotation=0 # CDS annotation against UniProt Swiss-Prot.
uniprot_annotation_method="mmseqs2" # blastp|mmseqs2 for UniProt Swiss-Prot annotation search engine.
run_cds_fx2tab=0 # Sequence stats of CDS sequences
run_cds_mmseqs2taxonomy=0 # Taxonomic assignment of CDS sequences
run_cds_contamination_removal=0 # Removal of contaminated sequences
run_annotation=0 # Per-gene annotation summary
run_wgd_ksd=0 # WGD inference by dS distribution

# Genome analysis
run_busco_genome=0 # Gene set completeness analysis
run_subphaser=0 # Subgenome structure inference
run_genome_fx2tab=0 # Sequence stats of reference genome
run_scaffold_histogram=0 # Scaffold histogram
run_genome_mmseqs2taxonomy=0 # Taxonomic assignment of genome assembly
run_genome_contamination_removal=0 # Removal of contaminated sequences
run_jcvi_dotplot=0 # Self-self synteny dotplot

# DNA-seq analysis
run_genomescope=0 # GenomeScope

# Summary
run_multispecies_summary=1 # Multi-species summary plots and tables

### End: Modify this block to tailor your analysis ###

# Misc
exit_if_running=0 # Exit without main analysis if the same SGE_TASK_ID is already running.
delete_tmp_dir=1 # After this run, delete tmp directory created for each job. Set 0 when debugging.

source "${dir_script}/support/gg_util.sh" # loading utility functions
# Forward config variables (including external overrides) into container environment.
forward_config_vars_to_container_env "${BASH_SOURCE[0]}" "delete_tmp_dir"
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
${singularity_command} "${gg_image}" < "${dir_script}/core/gg_genome_annotation_core.sh"
if ! gg_trigger_versions_dump "$(basename "${BASH_SOURCE[0]}")"; then
  echo "Warning: gg_versions trigger failed."
fi

echo "$(date): Ending"

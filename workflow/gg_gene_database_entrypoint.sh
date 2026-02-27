#!/usr/bin/env bash
set -euo pipefail

# SLURM in NIG supercomputer
#SBATCH -J gg_gene_database
#SBATCH -c 1 # Number of CPUs
#SBATCH --mem-per-cpu=64G # RAM per CPU in GB
#SBATCH -t 2976:00:00 # maximum time in d-hh:mm:ss format. NIG supercomputer epyc/medium MaxTime=2976:00:00
#SBATCH --output=gg_gene_database_%j.out
#SBATCH --error=gg_gene_database_%j.err
#SBATCH -p medium # partition name, cluster environment specific
#SBATCH --chdir=.
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>
##SBATCH -N 1-1 # Number of nodes
##SBATCH -n 1 # Number of tasks

#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 1
#$ -l s_vmem=64G
#$ -l mem_req=64G
#$ -l epyc
#$ -l d_rt=124:00:00:00
#$ -l s_rt=124:00:00:00

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# Fixed to 1

echo "$(date): Starting"

# Change these directories for your custom-made analysis
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dir_pg="${script_dir}/../workspace" # pg input and output directory
dir_script="${script_dir}" # directory where gg_util.sh and core/gg_*_core.sh locate
gg_image="${script_dir}/../genegalleon.sif" # path to the singularity image

### Start: Modify this block to tailor your analysis ###

run_database_prep="${run_database_prep:-1}"

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
${singularity_command} "${gg_image}" < "${dir_script}/core/gg_gene_database_core.sh"
if ! gg_trigger_versions_dump "$(basename "${BASH_SOURCE[0]}")"; then
  echo "Warning: gg_versions trigger failed."
fi

echo "$(date): Ending"

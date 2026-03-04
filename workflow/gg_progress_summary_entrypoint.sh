#!/usr/bin/env bash
set -euo pipefail

# SLURM in NIG supercomputer
#SBATCH -J gg_progress_summary
#SBATCH -c 4 # Number of CPUs
#SBATCH --mem-per-cpu=8G # RAM per CPU in GB
#SBATCH -t 2976:00:00 # maximum time in d-hh:mm:ss format. NIG supercomputer epyc/medium MaxTime=2976:00:00
#SBATCH --output=gg_progress_summary_%j.out
#SBATCH --error=gg_progress_summary_%j.err
#SBATCH -p medium # partition name, cluster environment specific
#SBATCH --chdir=.
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>
##SBATCH -N 1-1 # Number of nodes
##SBATCH -n 1 # Number of tasks

## UGE in NIG supercomputer
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -l s_vmem=8G
#$ -l mem_req=8G
#$ -l epyc
#$ -l d_rt=124:00:00:00
#$ -l s_rt=124:00:00:00

## PBS in BIAS at NIBB
##PBS -S /bin/bash
##PBS -l ncpus=4
##PBS -l mem=32G
##PBS -q small
##PBS -V

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# Fixed to 1

echo "$(date): Starting"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dir_script="${script_dir}"
dir_pg="${script_dir}/../workspace"
gg_image="${script_dir}/../genegalleon.sif"

### Start: Modify this block to tailor your analysis ###
mode_transcriptome_assembly="auto" # {"auto", "sraid", "fastq", "metadata"}
ncpu_progress_summary=""
### End: Modify this block to tailor your analysis ###

source "${script_dir}/support/gg_util.sh"
gg_scheduler_runtime_prelude
unset_singularity_envs
if ! set_singularity_command; then
  exit 1
fi
variable_SGEnizer
: "${ncpu_progress_summary:=${NSLOTS:-1}}"
set_singularityenv
gg_print_scheduler_runtime_summary

forward_config_vars_to_container_env "${BASH_SOURCE[0]}"

mkdir -p "${dir_pg}"
cd "${dir_pg}"
${singularity_command} "${gg_image}" < "${dir_script}/core/gg_progress_summary_core.sh"
if ! gg_trigger_versions_dump "$(basename "${BASH_SOURCE[0]}")"; then
  echo "Warning: gg_versions trigger failed."
fi

echo "$(date): Ending"

#!/usr/bin/env bash
set -eo pipefail

# SLURM in NIG supercomputer
#SBATCH -J gg_cdsAnnotation
#SBATCH -c 4 # Number of CPUs
#SBATCH --mem-per-cpu=8G # RAM per CPU in GB
#SBATCH -t 2976:00:00 # maximum time in d-hh:mm:ss format. NIG supercomputer epyc/medium MaxTime=2976:00:00
#SBATCH --output=gg_cdsAnnotation_%A_%a.out
#SBATCH --error=gg_cdsAnnotation_%A_%a.err
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
#$ -l d_rt=30:00:00:00
#$ -l s_rt=30:00:00:00
#$ -t 1

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# N = Number of fasta files in workspace/input/species_cds

echo "$(date): Starting"
ulimit -s unlimited 2>/dev/null || true

# Change these directories for your custom-made analysis
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dir_pg="${script_dir}/../workspace" # pg input and output directory
dir_script="${script_dir}" # directory where gg_util.sh and gg_*_cmd.sh locate
gg_image="${script_dir}/../genegalleon.sif" # path to the singularity image

# Misc
exit_if_running=0 # Exit without main analysis if the same SGE_TASK_ID is already running.
delete_tmp_dir=1 # After this run, delete tmp directory created for each job. Set 0 when debugging.

source "${dir_script}/script/gg_util.sh" # loading utility functions
unset_singularity_envs
set_singularity_command
variable_SGEnizer
set_singularityenv

cd "${dir_pg}"
${singularity_command} "${gg_image}" < "${dir_script}/gg_cdsAnnotation_cmd.sh"

echo "$(date): Ending"

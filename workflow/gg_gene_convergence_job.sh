#!/usr/bin/env bash
set -eo pipefail

# SLURM in NIG supercomputer
#SBATCH -J gg_gene_convergence
#SBATCH -c 8 # Number of CPUs
#SBATCH --mem-per-cpu=32G # RAM per CPU in GB
#SBATCH -t 60-00:00:00 # maximum time in d-hh:mm:ss format
#SBATCH --output=gg_gene_convergence_%A_%a.out
#SBATCH --error=gg_gene_convergence_%A_%a.err
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
#$ -pe def_slot 8
#$ -l s_vmem=32G
#$ -l mem_req=32G
#$ -l medium-ubuntu
#$ -l d_rt=30:00:00:00
#$ -l s_rt=30:00:00:00
#$ -t 1

## PBS in BIAS at NIBB
##PBS -s /bin/bash
##PBS -l ncpus=1
##PBS -l mem=4G
##PBS -J 1
##PBS -q small
##PBS -V

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# Fixed to 1

if [[ -n "${PBS_O_WORKDIR:-}" ]]; then # Instead of -cwd
  cd "${PBS_O_WORKDIR}"
  export PATH="${PATH}:/bio/package/singularity/singularity_3.0/bin"
fi

echo "$(date): Starting"
ulimit -s unlimited 2>/dev/null || true

# Change these directories for your custom-made analysis
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dir_pg="${script_dir}/../workspace" # pg input and output directory
dir_script="${script_dir}" # directory where gg_*_cmd.sh and script/ locate
gg_image="${script_dir}/../genegalleon.sif" # path to the singularity image

# Misc
exit_if_running=0
delete_tmp_dir=1

source "${dir_script}/script/gg_util.sh" # loading utility functions
unset_singularity_envs
set_singularity_command
variable_SGEnizer
set_singularityenv

mkdir -p "${dir_pg}"
cd "${dir_pg}"
${singularity_command} "${gg_image}" < "${dir_script}/gg_gene_convergence_cmd.sh"
if ! gg_trigger_versions_dump "$(basename "${BASH_SOURCE[0]}")"; then
  echo "Warning: gg_versions trigger failed."
fi

echo "$(date): Ending"

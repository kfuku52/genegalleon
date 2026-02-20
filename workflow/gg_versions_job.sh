#!/usr/bin/env bash
set -eo pipefail

# SLURM in NIG supercomputer
#SBATCH -J gg_versions
#SBATCH -c 1 # Number of CPUs
#SBATCH --mem-per-cpu=8G # RAM per CPU in GB
#SBATCH -t 2976:00:00 # maximum time in d-hh:mm:ss format. NIG supercomputer epyc/medium MaxTime=2976:00:00
#SBATCH --output=gg_versions_%A_%a.out
#SBATCH --error=gg_versions_%A_%a.err
#SBATCH -p medium # partition name, cluster environment specific
#SBATCH --chdir=.
#SBATCH -a 1 # Array job, 1-N
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>
##SBATCH -N 1-1 # Number of nodes
##SBATCH -n 1 # Number of tasks

#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 1
#$ -l s_vmem=4G
#$ -l mem_req=4G
##$ -l epyc
#$ -l d_rt=30:00:00:00
#$ -l s_rt=30:00:00:00
#$ -t 1

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# Fixed to 1

echo "$(date): Starting"
#ulimit -s unlimited

# Change these directories for your custom-made analysis
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dir_pg="${script_dir}/../workspace" # pg input and output directory
dir_script="${script_dir}" # directory where gg_util.sh and gg_*_cmd.sh locate
gg_image="${script_dir}/../genegalleon.sif" # path to the singularity image
file_version="${script_dir}/../VERSION"
gg_version="unknown"
if [[ -s "${file_version}" ]]; then
  gg_version="$(head -n 1 "${file_version}" | tr -d '\r')"
fi
echo "### genegalleon version ###"
echo "${gg_version}"
echo ''
echo ''

source "${dir_script}/script/gg_util.sh" # loading utility functions
unset_singularity_envs
set_singularity_command
variable_SGEnizer
set_singularityenv
export SINGULARITYENV_GG_VERSION="${gg_version}"

if ! command -v singularity >/dev/null 2>&1; then
	echo "singularity command not found on host. Exiting."
	exit 1
fi

cd "${dir_pg}"
if ! gg_trigger_versions_dump "$(basename "${BASH_SOURCE[0]}")"; then
  echo "Warning: gg_versions trigger failed."
fi

echo "$(date): Ending"

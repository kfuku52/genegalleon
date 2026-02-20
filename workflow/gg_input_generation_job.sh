#!/usr/bin/env bash
set -eo pipefail

# SLURM in NIG supercomputer
#SBATCH -J gg_input_generation
#SBATCH -c 2 # Number of CPUs
#SBATCH --mem-per-cpu=8G # RAM per CPU in GB
#SBATCH -t 2976:00:00 # maximum time in d-hh:mm:ss format. NIG supercomputer epyc/medium MaxTime=2976:00:00
#SBATCH --output=gg_input_generation_%A_%a.out
#SBATCH --error=gg_input_generation_%A_%a.err
#SBATCH -p medium # partition name, cluster environment specific
#SBATCH --chdir=.
#SBATCH -a 1 # Array job, 1-N
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>
##SBATCH -N 1-1 # Number of nodes
##SBATCH -n 1 # Number of tasks

#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 2
#$ -l s_vmem=8G
#$ -l mem_req=8G
##$ -l epyc
#$ -l d_rt=30:00:00:00
#$ -l s_rt=30:00:00:00
#$ -t 1

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# Fixed to 1

echo "$(date): Starting"
ulimit -s unlimited 2>/dev/null || true

# Change these directories for your custom-made analysis
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dir_pg="${script_dir}/../workspace" # pg input and output directory
dir_script="${script_dir}" # directory where gg_util.sh and gg_*_cmd.sh locate
gg_image="${script_dir}/../genegalleon.sif" # path to the singularity image

source "${dir_script}/script/gg_util.sh" # loading utility functions
unset_singularity_envs
set_singularity_command
variable_SGEnizer
set_singularityenv

if [[ -n "${GG_INPUT_DATASET_HOST_PATH}" ]]; then
  if [[ ! -d "${GG_INPUT_DATASET_HOST_PATH}" ]]; then
    echo "GG_INPUT_DATASET_HOST_PATH was set but directory was not found: ${GG_INPUT_DATASET_HOST_PATH}"
    exit 1
  fi
  dataset_host_path=$(cd "$(dirname "${GG_INPUT_DATASET_HOST_PATH}")" && pwd -P)/$(basename "${GG_INPUT_DATASET_HOST_PATH}")
  if [[ -n "${SINGULARITY_BIND}" ]]; then
    SINGULARITY_BIND="${dataset_host_path}:/external/gfe_dataset,${SINGULARITY_BIND}"
  else
    SINGULARITY_BIND="${dataset_host_path}:/external/gfe_dataset"
  fi
  if [[ -n "${SINGULARITY_BINDPATH}" ]]; then
    SINGULARITY_BINDPATH="${dataset_host_path}:/external/gfe_dataset,${SINGULARITY_BINDPATH}"
  else
    SINGULARITY_BINDPATH="${dataset_host_path}:/external/gfe_dataset"
  fi
  export SINGULARITY_BIND
  export SINGULARITY_BINDPATH
  export SINGULARITYENV_GG_DATASET_ROOT="/external/gfe_dataset"
  echo "Binding external dataset root: ${dataset_host_path} -> /external/gfe_dataset"
fi

if ! command -v singularity >/dev/null 2>&1; then
  echo "singularity command not found on host. Exiting."
  exit 1
fi

cd "${dir_pg}"
${singularity_command} "${gg_image}" < "${dir_script}/gg_input_generation_cmd.sh"
if ! gg_trigger_versions_dump "$(basename "${BASH_SOURCE[0]}")"; then
  echo "Warning: gg_versions trigger failed."
fi

echo "$(date): Ending"

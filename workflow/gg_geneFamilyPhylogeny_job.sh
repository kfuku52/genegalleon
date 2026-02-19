#!/usr/bin/env bash
set -eo pipefail

# SLURM in NIG supercomputer
#SBATCH -J gg_geneFamilyPhylogeny
#SBATCH -c 4 # Number of CPUs
#SBATCH --mem-per-cpu=8G # RAM per CPU in GB
#SBATCH -t 2976:00:00 # maximum time in d-hh:mm:ss format. NIG supercomputer epyc/medium MaxTime=2976:00:00
#SBATCH --output=gg_geneFamilyPhylogeny_%A_%a.out
#SBATCH --error=gg_geneFamilyPhylogeny_%A_%a.err
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
#$ -pe def_slot 4
#$ -l s_vmem=4G
#$ -l mem_req=4G
#$ -l epyc
#$ -l d_rt=30:00:00:00
#$ -l s_rt=30:00:00:00
#$ -t 1

## PBS in BIAS at NIBB
##PBS -s /bin/bash
##PBS -l ncpus=4
##PBS -l mem=4G
##PBS -J 1-3
##PBS -q small
##PBS -V

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# N = Number of rows (excluding the header) in workspace/output/orthofinder/Orthogroups/Orthogroups.GeneCount.selected.tsv for mode_orthogroup=1
# N = Number of files in workspace/input/query2family_input for mode_query2family=1

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
exit_if_running=0 # Exit without main analysis if the same SGE_TASK_ID is already running.
delete_tmp_dir=1 # After this run, delete tmp directory created for each job. Set 0 when debugging.
delete_preexisting_tmp_dir=1 # Before starting this job, delete tmp directory created by previous run.
export SINGULARITYENV_run_rps_blast="${SINGULARITYENV_run_rps_blast:-${run_rps_blast:-1}}" # Set to 0 explicitly to skip Pfam domain annotation.
export run_rps_blast="${run_rps_blast:-1}" # The Docker-backed singularity shim forwards this variable directly.
export SINGULARITYENV_run_generax="${SINGULARITYENV_run_generax:-0}" # GeneRax requires MPI launcher setup not available in local smoke runs.
export run_generax="${run_generax:-0}" # Forwarded directly by the Docker-backed singularity shim.
export SINGULARITYENV_run_midpoint_root="${SINGULARITYENV_run_midpoint_root:-1}" # Provide rooted trees without GeneRax in local validation.
export run_midpoint_root="${run_midpoint_root:-1}" # Forwarded directly by the Docker-backed singularity shim.

source "${dir_script}/script/gg_util.sh" # loading utility functions
unset_singularity_envs
exit_if_running_qstat
set_singularity_command
variable_SGEnizer
set_singularityenv

cd "${dir_pg}"
set +e
${singularity_command} "${gg_image}" < "${dir_script}/gg_geneFamilyPhylogeny_cmd.sh"
first_run_exit_code=$?
set -e
if [[ ${first_run_exit_code} -eq 8 ]]; then
  echo "Output files were detected. No more run of gg_geneFamilyPhylogeny_cmd.sh is necessary."
else
  if [[ ${first_run_exit_code} -ne 0 ]]; then
    echo "First run exited with code ${first_run_exit_code}. Retrying once."
  fi
  # Jobs are misteriously killed after GeneRax, so run twice to complete all steps if output files were not detected.
  set +e
  ${singularity_command} "${gg_image}" < "${dir_script}/gg_geneFamilyPhylogeny_cmd.sh"
  second_run_exit_code=$?
  set -e
  if [[ ${second_run_exit_code} -eq 8 ]]; then
    echo "Second run reported output files already detected (exit code 8)."
  elif [[ ${second_run_exit_code} -ne 0 ]]; then
    echo "Second run failed with code ${second_run_exit_code}."
    exit "${second_run_exit_code}"
  fi
fi
echo "$(date): Ending"

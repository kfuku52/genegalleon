#!/usr/bin/env bash
set -eo pipefail

# SLURM in NIG supercomputer
#SBATCH -J gg_test
#SBATCH -c 4 # Number of CPUs
#SBATCH --mem-per-cpu=8G # RAM per CPU in GB
#SBATCH -t 2976:00:00 # maximum time in d-hh:mm:ss format. NIG supercomputer epyc/medium MaxTime=2976:00:00
#SBATCH --output=gg_test_%A_%a.out
#SBATCH --error=gg_test_%A_%a.err
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
#$ -pe def_slot 2
#$ -l s_vmem=4G
#$ -l mem_req=4G
#$ -l epyc
#$ -l d_rt=30:00:00:00
#$ -l s_rt=30:00:00:00
#$ -t 1-3

## PBS in BIAS at NIBB
##PBS -s /bin/bash
##PBS -l ncpus=4
##PBS -l mem=4G
##PBS -J 1-3
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

exit_if_running=0 # Exit without main analysis if the same SGE_TASK_ID is already running.
delete_tmp_dir=0 # After this run, delete tmp directory created for each job. Set 0 when debugging.
delete_preexisting_tmp_dir=1 # Before starting this job, delete tmp directory created by previous run.

export NSLOTS=4
export SINGULARITYENV_pg_debug_mode="${SINGULARITYENV_pg_debug_mode:-0}"
export SINGULARITYENV_delete_tmp_dir=0
export SINGULARITYENV_run_rps_blast="${SINGULARITYENV_run_rps_blast:-0}"
export run_rps_blast="${run_rps_blast:-0}"
export SINGULARITYENV_run_generax="${SINGULARITYENV_run_generax:-0}"
export run_generax="${run_generax:-0}"
export SINGULARITYENV_run_midpoint_root="${SINGULARITYENV_run_midpoint_root:-1}"
export run_midpoint_root="${run_midpoint_root:-1}"

source "${dir_script}/script/gg_util.sh" # loading utility functions
unset_singularity_envs
exit_if_running_qstat
set_singularity_command
variable_SGEnizer
set_singularityenv

script_files=(
#gg_test_cmd.sh
#gg_versions_cmd.sh
#gg_orthofinder_cmd.sh
gg_geneFamilyPhylogeny_cmd.sh
gg_geneFamilyPhylogeny_cmd.sh
#gg_genomeEvolution_cmd.sh
#gg_cdsAnnotation_cmd.sh
#gg_speciesTree_cmd.sh
#gg_createSpeciesPromoterFasta_cmd.sh
#gg_orthogroupDatabasePrep_cmd.sh
#gg_transcriptomeAssembly_cmd.sh
)

cd "${dir_pg}"
for script_file in "${script_files[@]}"; do
  echo "Starting: ${script_file}"
  prefix="${script_file}"
  command="${singularity_command} \"${gg_image}\" < \"${dir_script}/${script_file}\""
  echo "full command: ${command}"
  set +e
  ${singularity_command} "${gg_image}" < "${dir_script}/${script_file}" \
    > >(tee "${prefix}.out") \
    2> >(tee "${prefix}.err" >&2)
  script_rc=$?
  set -e
  if [[ ${script_rc} -eq 8 ]]; then
    echo "Skipped: ${script_file} (output files were already detected, exit code 8)"
  elif [[ ${script_rc} -ne 0 ]]; then
    echo "Failed: ${script_file} (exit code ${script_rc})"
    exit "${script_rc}"
  fi
  echo "Ending: ${script_file}"
  echo ""
  echo ""
  echo ""
done

echo "$(date): Ending"

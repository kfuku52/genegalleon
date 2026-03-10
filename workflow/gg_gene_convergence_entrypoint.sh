#!/usr/bin/env bash

# Scheduler header notes:
# - Keep sections ordered as SLURM -> UGE -> PBS across entrypoints.
# - Update job name, CPU count, memory, walltime, log paths, and array size together.
# - Site-specific partition/queue/resource lines stay commented examples by default.

# SLURM
# Common parameters: job name, cores per task, memory per core, walltime, log files, and working directory.
#SBATCH -J gg_gene_convergence
#SBATCH -c 8
#SBATCH --mem-per-cpu=32G
#SBATCH -t 60-00:00:00
#SBATCH --output=gg_gene_convergence_%j.out
#SBATCH --error=gg_gene_convergence_%j.err
#SBATCH --chdir=.
# Site-specific partition example.
#SBATCH -p epyc
# Optional notifications and single-node examples.
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>

## UGE
# Common parameters: shell, working directory, slot count, memory per slot, and runtime limits.
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 8
#$ -l s_vmem=32G
#$ -l mem_req=32G
# Site-specific resource example.
#$ -l medium-ubuntu
#$ -l d_rt=60:00:00:00
#$ -l s_rt=60:00:00:00

## PBS
# Common parameters: shell, CPU count, total memory, and exported environment.
#PBS -S /bin/bash
#PBS -l ncpus=8
#PBS -l mem=256G
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
if ! gg_entrypoint_initialize "${BASH_SOURCE[0]}" 0; then
  exit 1
fi
gg_entrypoint_name="gg_gene_convergence_entrypoint.sh"

### Start: Modify this block to tailor your analysis ###

arity_range="2-10" # Branch-combination arity range to scan.
trait="all" # Trait column name(s) to analyze, or "all".
skip_lower_order="yes" # Skip lower-order combinations once higher-order hits are retained.
min_fg_stem_ratio="0.5" # Minimum fraction of foreground branches that must be stem branches.
min_OCNany2spe="1.8" # Minimum OCNany2spe cutoff for candidate combinations.
min_omegaCany2spe="3.0" # Minimum omegaCany2spe cutoff for candidate combinations.
min_OCNCoD=0 # Minimum OCNCoD cutoff for candidate combinations.
max_per_K=100 # Maximum number of combinations retained per arity K.
file_trait="auto" # Trait table path, or "auto" to use the workspace default.
dir_orthogroup="auto" # Orthogroup result directory, or "auto" to detect it from the workspace.
dir_orthofinder="auto" # OrthoFinder result directory, or "auto" to detect it from the workspace.
dir_out="auto" # Output directory for convergence results, or "auto" to use the default workspace path.

### End: Modify this block to tailor your analysis ###

# Misc
exit_if_running=0

source "${gg_support_dir}/gg_util.sh" # loading utility functions
# Forward config variables (including external overrides) into container environment.
forward_config_vars_to_container_env "${gg_entrypoint_name}"
gg_export_var_to_container_env_if_set "PYMOL_HEADLESS"
gg_export_var_to_container_env_if_set "QT_QPA_PLATFORM"
if ! gg_entrypoint_prepare_container_runtime 1; then
  exit 1
fi
gg_entrypoint_activate_container_runtime

gg_entrypoint_enter_workspace
gg_run_container_shell_script "${gg_container_image_path}" "${gg_core_dir}/gg_gene_convergence_core.sh"
gg_require_versions_dump "${gg_entrypoint_name}"

echo "$(date): Ending"

#!/usr/bin/env bash

# Scheduler header notes:
# - Keep sections ordered as SLURM -> UGE -> PBS across entrypoints.
# - Update job name, CPU count, memory, walltime, log paths, and array size together.
# - UGE resources default to SHIROKANE AGE; SLURM/PBS site-specific lines remain examples.

# SLURM
# Common parameters: job name, cores per task, memory per core, walltime, log files, and working directory.
#SBATCH -J gg_progress_summary
#SBATCH -c 4
#SBATCH --mem-per-cpu=8G
#SBATCH -t 3:00:00
#SBATCH --output=gg_progress_summary_entrypoint.sh_%j.out
#SBATCH --error=gg_progress_summary_entrypoint.sh_%j.err
#SBATCH --chdir=.
#SBATCH --ignore-pbs
# Site-specific partition example.
#SBATCH -p epyc,rome,medium
# Optional notifications and single-node examples.
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>

## UGE
# SHIROKANE AGE defaults: shell, working directory, slot count, memory per slot, and ljob.
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -l s_vmem=8G
#$ -l ljob

## PBS
# Common parameters: shell, CPU count, total memory, and exported environment.
#PBS -S /bin/bash
#PBS -l ncpus=4
#PBS -l mem=32G
# Site-specific queue example.
#PBS -q small
#PBS -V

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# Fixed to 1

set -euo pipefail

echo "$(date): Starting"

# Resolve workflow paths for local and scheduler-spooled execution.
gg_bootstrap_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
gg_bootstrap_checked_bases=""
for gg_bootstrap_base in \
  "${SLURM_SUBMIT_DIR:-}" \
  "${PBS_O_WORKDIR:-}" \
  "${PWD:-}" \
  "${gg_bootstrap_script_dir}"
do
  [[ -n "${gg_bootstrap_base}" ]] || continue
  case ":${gg_bootstrap_checked_bases}:" in
    *":${gg_bootstrap_base}:"*) continue ;;
  esac
  gg_bootstrap_checked_bases="${gg_bootstrap_checked_bases:+${gg_bootstrap_checked_bases}:}${gg_bootstrap_base}"
  for bootstrap_path in \
    "${gg_bootstrap_base}/support/gg_entrypoint_bootstrap.sh" \
    "${gg_bootstrap_base}/workflow/support/gg_entrypoint_bootstrap.sh"
  do
    if [[ -s "${bootstrap_path}" ]]; then
      # shellcheck disable=SC1090
      source "${bootstrap_path}"
      break
    fi
  done
  if declare -F gg_entrypoint_initialize >/dev/null 2>&1; then
    break
  fi
done
unset gg_bootstrap_base gg_bootstrap_checked_bases gg_bootstrap_script_dir bootstrap_path
if ! declare -F gg_entrypoint_initialize >/dev/null 2>&1; then
  echo "Failed to locate gg_entrypoint_bootstrap.sh from BASH_SOURCE[0]=${BASH_SOURCE[0]}" >&2
  exit 1
fi
if ! gg_entrypoint_initialize "${BASH_SOURCE[0]}" 0; then
  exit 1
fi
gg_entrypoint_name="gg_progress_summary_entrypoint.sh"

### Start: Modify this block to tailor your analysis ###

# Mode
mode_transcriptome_assembly="auto" # {"auto", "sraid", "fastq", "metadata"}; input mode used to interpret transcriptome-generation progress directories, with auto detecting the available workspace/input layout.
gene_family_output_storage="${gene_family_output_storage:-${GG_COMMON_GENE_FAMILY_OUTPUT_STORAGE:-zip}}" # zip|files|raw; raw aliases files, while zip flushes completed live gene-family artifacts after archive-aware summaries.
gene_family_zip_compression="${gene_family_zip_compression:-${GG_COMMON_GENE_FAMILY_ZIP_COMPRESSION:-adaptive}}" # adaptive|deflate|store.
gene_family_zip_compression_level="${gene_family_zip_compression_level:-${GG_COMMON_GENE_FAMILY_ZIP_COMPRESSION_LEVEL:-6}}" # Deflate level 0-9.
gene_family_zip_workers="${gene_family_zip_workers:-${GG_COMMON_GENE_FAMILY_ZIP_WORKERS:-1}}" # Concurrent ZIP writers, capped at 4.
gene_family_final_zip_max_bytes="${gene_family_final_zip_max_bytes:-${GG_COMMON_GENE_FAMILY_FINAL_ZIP_MAX_BYTES:-0}}" # 0 permits one final ZIP of any size; a positive byte limit retains named part ZIPs for larger subdirectories.
gene_family_tmp_retention_days="${gene_family_tmp_retention_days:-${GG_COMMON_GENE_FAMILY_TMP_RETENTION_DAYS:-7}}" # Failed task directories below query2family/orthogroup tmp roots are removed after this many days when their family lock is idle; 0 disables the age limit.
gene_family_tmp_max_dirs="${gene_family_tmp_max_dirs:-${GG_COMMON_GENE_FAMILY_TMP_MAX_DIRS:-100}}" # Maximum failed task directories retained per query2family/orthogroup root; oldest excess directories are removed and 0 disables the count limit.
gene_family_tmp_max_bytes="${gene_family_tmp_max_bytes:-${GG_COMMON_GENE_FAMILY_TMP_MAX_BYTES:-107374182400}}" # Maximum bytes retained across inactive failed task directories; 0 disables the byte limit.
gene_family_tmp_max_files="${gene_family_tmp_max_files:-${GG_COMMON_GENE_FAMILY_TMP_MAX_FILES:-100000}}" # Maximum files retained across inactive failed task directories; 0 disables the file-count limit.

# Runtime parameters
ncpu_progress_summary="" # Number of CPU threads used by summary scripts; empty falls back to GG_TASK_CPUS.

### End: Modify this block to tailor your analysis ###

source "${gg_support_dir}/gg_util.sh"
if ! gg_entrypoint_prepare_container_runtime 0; then
  exit 1
fi
: "${ncpu_progress_summary:=${GG_TASK_CPUS:-1}}"
gg_entrypoint_activate_container_runtime

gg_apply_registered_env_overrides "${gg_entrypoint_name}"
forward_config_vars_to_container_env "${gg_entrypoint_name}"

gg_entrypoint_enter_workspace
gg_runtime_core_script="$(gg_prepare_entrypoint_runtime_snapshot "${gg_entrypoint_name}" "${gg_core_dir}/gg_progress_summary_core.sh")"
gg_run_container_shell_script "${gg_container_image_path}" "${gg_runtime_core_script}"
gg_require_versions_dump "${gg_entrypoint_name}"

echo "$(date): Ending"

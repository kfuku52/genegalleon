#!/usr/bin/env bash
set -eo pipefail

# SLURM in NIG supercomputer
#SBATCH -J gg_build
#SBATCH -c 2 # Number of CPUs
#SBATCH --mem-per-cpu=32G # RAM per CPU in GB
#SBATCH -t 2976:00:00 # maximum time in d-hh:mm:ss format. NIG supercomputer epyc/medium MaxTime=2976:00:00
#SBATCH --output=gg_build_%A_%a.out
#SBATCH --error=gg_build_%A_%a.err
#SBATCH -p epyc # partition name, cluster environment specific
#SBATCH --chdir=.
#SBATCH -a 1 # Array job, 1-N
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>
##SBATCH -N 1-1 # Number of nodes
##SBATCH -n 1 # Number of tasks

#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 2
#$ -l s_vmem=32G
#$ -l mem_req=32G
#$ -l epyc
#$ -l d_rt=31:00:00:00
#$ -l s_rt=31:00:00:00
#$ -t 1

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# Fixed to 1

echo "$(date): Starting"
ulimit -s unlimited 2>/dev/null || true

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dir_script="${script_dir}"
dir_pg="${script_dir}/../workspace"
cd "${script_dir}"
source "${dir_script}/script/gg_util.sh"

SECONDS=0
export TMPDIR=$(pwd)
my_date=$(date '+%Y%m%d')
host_arch="$(uname -m)"
case "${host_arch}" in
  aarch64|arm64)
    arch_tag="arm"
    ;;
  x86_64|amd64)
    arch_tag="amd"
    ;;
  *)
    arch_tag="$(echo "${host_arch}" | tr '[:upper:]' '[:lower:]' | tr -cd '[:alnum:]_-')"
    ;;
esac
build_dir="${BUILD_DIR:-genegalleon}"
if [[ ! -d "${build_dir}" && -d pg ]]; then
  # Backward compatibility for legacy local build directories.
  build_dir="pg"
fi

sif_name="genegalleon_${my_date}_${arch_tag}.sif"
zip_pid=""
if [[ ! -d "${build_dir}" ]]; then
  echo "Build directory not found: ${build_dir}"
  exit 1
fi

zip -ryq genegalleon_${my_date}.zip "${build_dir}" &
zip_pid=$!
singularity build "${sif_name}" "${build_dir}"
echo "Generated ${sif_name}"
echo "$SECONDS sec elapsed in singularity build"
if [[ -n "${zip_pid}" ]]; then
  wait "${zip_pid}"
fi
gg_image="${script_dir}/${sif_name}"
if ! gg_trigger_versions_dump "$(basename "${BASH_SOURCE[0]}")"; then
  echo "Warning: gg_versions trigger failed."
fi

echo "$(date): Ending"

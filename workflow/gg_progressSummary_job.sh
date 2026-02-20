#!/usr/bin/env bash
set -eo pipefail

# SLURM in NIG supercomputer
#SBATCH -J gg_progressSummary
#SBATCH -c 4 # Number of CPUs
#SBATCH --mem-per-cpu=8G # RAM per CPU in GB
#SBATCH -t 2976:00:00 # maximum time in d-hh:mm:ss format. NIG supercomputer epyc/medium MaxTime=2976:00:00
#SBATCH --output=gg_progressSummary_%A_%a.out
#SBATCH --error=gg_progressSummary_%A_%a.err
#SBATCH -p medium # partition name, cluster environment specific
#SBATCH --chdir=.
#SBATCH -a 1 # Array job, 1-N
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>
##SBATCH -N 1-1 # Number of nodes
##SBATCH -n 1 # Number of tasks

## UGE in NIG supercomputer
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 1
##$ -l s_vmem=4G
#$ -l mem_req=4G
#$ -l epyc
##$ -l d_rt=30:00:00:00
##$ -l s_rt=30:00:00:00
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

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dir_script="${script_dir}"
dir_pg="${script_dir}/../workspace"
gg_image="${script_dir}/../genegalleon.sif"
source "${script_dir}/script/gg_util.sh"
dir_pg_input=$(workspace_input_root "${dir_pg}")
dir_pg_output=$(workspace_output_root "${dir_pg}")
dir_orthogroup="${dir_pg_output}/orthogroup"
dir_transcriptome_assembly="${dir_pg_output}/transcriptome_assembly"
file_genecount="${dir_pg_output}/orthofinder/Orthogroups/Orthogroups.GeneCount.selected.tsv"
mode_transcriptome_assembly="auto" # {"auto", "sraid", "fastq", "metadata"} — Set to the same mode as in gg_transcriptomeAssembly_cmd.sh, or use "auto".
ncpu_progress_summary="${NSLOTS:-1}"

if [[ -d "${dir_orthogroup}" ]]; then
	echo ""
	echo "Checking directory: ${dir_orthogroup}"
	if [[ ! -s "${file_genecount}" ]]; then
		echo "Skipping orthogroup summary because the gene-count table was not found: ${file_genecount}"
	elif [[ ! -d "${dir_orthogroup}/amas.original" || ! -d "${dir_orthogroup}/amas.cleaned" ]]; then
		echo "Skipping orthogroup summary because required AMAS directories were not found under: ${dir_orthogroup}"
	else
		python "${script_dir}/script/orthogroup_output_summary.py" \
		--dir_og "${dir_orthogroup}" \
		--genecount "${file_genecount}" \
		--ncpu "${ncpu_progress_summary}" \
		--out orthogroup_summary.tsv
	fi
fi

if [[ -d "${dir_transcriptome_assembly}" ]]; then
	echo ""
	echo "Checking directory: ${dir_transcriptome_assembly}"
	python "${script_dir}/script/transcriptome_assembly_output_summary.py" \
	--dir_transcriptome_assembly "${dir_transcriptome_assembly}" \
	--mode "${mode_transcriptome_assembly}" \
	--ncpu "${ncpu_progress_summary}" \
	--out transcriptome_assembly_summary.tsv
fi
if ! gg_trigger_versions_dump "$(basename "${BASH_SOURCE[0]}")"; then
  echo "Warning: gg_versions trigger failed."
fi

echo "$(date): Ending"

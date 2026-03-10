#!/usr/bin/env bash

# SLURM in NIG supercomputer
#SBATCH -J gg_genome_annotation
#SBATCH -c 4
#SBATCH --mem-per-cpu=8G
#SBATCH -t 2976:00:00
#SBATCH --output=gg_genome_annotation_%A_%a.out
#SBATCH --error=gg_genome_annotation_%A_%a.err
#SBATCH -p epyc
#SBATCH --chdir=.
#SBATCH -a 1
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
#$ -l d_rt=124:00:00:00
#$ -l s_rt=124:00:00:00
#$ -t 1

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# N = Number of fasta files in workspace/input/species_cds

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
if ! gg_entrypoint_initialize "${BASH_SOURCE[0]}" 1; then
  exit 1
fi
gg_entrypoint_name="gg_genome_annotation_entrypoint.sh"

### Start: Modify this block to tailor your analysis ###

# CDS analysis
run_get_gff_info=0 # Collect gene information from workspace/input/species_gff
run_busco_cds=0 # Gene set completeness analysis
run_uniprot_annotation=0 # CDS annotation against UniProt Swiss-Prot.
uniprot_annotation_method="mmseqs2" # blastp|mmseqs2 for UniProt Swiss-Prot annotation search engine.
run_cds_fx2tab=0 # Sequence stats of CDS sequences
run_cds_mmseqs2taxonomy=0 # Taxonomic assignment of CDS sequences
run_cds_contamination_removal=0 # Removal of contaminated sequences
run_annotation=0 # Per-gene annotation summary
run_wgd_ksd=0 # WGD inference by dS distribution

# Genome analysis
run_busco_genome=0 # Gene set completeness analysis
run_subphaser=0 # Subgenome structure inference
run_genome_fx2tab=0 # Sequence stats of reference genome
run_scaffold_histogram=0 # Scaffold histogram
run_genome_mmseqs2taxonomy=0 # Taxonomic assignment of genome assembly
run_genome_contamination_removal=0 # Removal of contaminated sequences
run_jcvi_dotplot=0 # Self-self synteny dotplot

# DNA-seq analysis
run_genomescope=0 # GenomeScope

# Summary
run_multispecies_summary=1 # Multi-species summary plots and tables

### End: Modify this block to tailor your analysis ###

# Misc
contamination_removal_rank="domain" # Taxonomic rank for contamination removal. Canonical value is domain; GeneGalleon normalizes tool-specific synonyms automatically.
exit_if_running=0 # Exit without main analysis if the same GG_ARRAY_TASK_ID is already running.
delete_tmp_dir=1 # After this run, delete tmp directory created for each job. Set 0 when debugging.

source "${gg_support_dir}/gg_util.sh" # loading utility functions
# Forward config variables (including external overrides) into container environment.
forward_config_vars_to_container_env "${gg_entrypoint_name}" "delete_tmp_dir"
if ! gg_entrypoint_prepare_container_runtime 1; then
  exit 1
fi
gg_entrypoint_activate_container_runtime

gg_entrypoint_enter_workspace
gg_run_container_shell_script "${gg_container_image_path}" "${gg_core_dir}/gg_genome_annotation_core.sh"
gg_require_versions_dump "${gg_entrypoint_name}"

echo "$(date): Ending"

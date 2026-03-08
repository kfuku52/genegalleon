#!/usr/bin/env bash

# SLURM in NIG supercomputer
#SBATCH -J gg_transcriptome_generation
#SBATCH -c 4
#SBATCH --mem-per-cpu=32G
#SBATCH -t 2976:00:00
#SBATCH --output=gg_transcriptome_generation_%A_%a.out
#SBATCH --error=gg_transcriptome_generation_%A_%a.err
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
#$ -l s_vmem=32G
#$ -l mem_req=32G
#$ -l epyc
#$ -l d_rt=124:00:00:00
#$ -l s_rt=124:00:00:00
#$ -t 1

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# N = Number of directories in workspace/input/species_rnaseq for mode_fastq=1
# N = Number of files in workspace/input/query_sra_id for mode_sraid=1

set -euo pipefail

echo "$(date): Starting"

# In many cases, 4 cores and mem_req=16G (=64G) worked well in per-SRA Trinity assembly.
# Increase per-core RAM to 32G if Trinity fails.

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
gg_entrypoint_name="gg_transcriptome_generation_entrypoint.sh"

### Start: Modify this block to tailor your analysis ###

mode_sraid=1 # Need input at workspace/input/query_sra_id.
mode_fastq=0 # Need input at workspace/input/species_rnaseq.
mode_metadata=0 # Need input at workspace/input/amalgkit_metadata.

run_amalgkit_metadata_or_integrate=1 # Metadata retrieval.
run_amalgkit_getfastq=1 # fastq generation from NCBI SRA.
run_assembly=1 # Transcriptome assembly with Trinity or rnaSPAdes.
run_longestcds=1 # Longest CDS extraction.
run_longestcds_fx2tab=1 # Sequence stats for longest CDS.
run_longestcds_mmseqs2taxonomy=0 # MMseqs2 taxonomy.
run_longestcds_contamination_removal=0 # Contamination removal.
run_busco_isoforms=1 # BUSCO for transcriptome assembly with isoforms.
run_busco_longest_cds=1 # BUSCO for longest CDS.
run_busco_contamination_removed_longest_cds=0 # BUSCO for contamination-removed longest CDS.
run_assembly_stat=1 # seqkit stat for assembly and CDS files.
run_amalgkit_quant=1 # Expression quantification.
run_amalgkit_merge=1 # Expression merge.
run_multispecies_summary=1 # Multi-species summary.

amalgkit_rrna_filter="yes" # read-level rRNA removal in amalgkit getfastq. Observed to finish for ~4.2 Gbp total / ~14 million reads within a 350G job before contamination filtering; exact peak RAM and elapsed time for this step alone were not logged.
amalgkit_contam_filter="no" # read-level contamination removal in amalgkit getfastq. Rank follows contamination_removal_rank below. Setting yes may require >350G RAM for large private FASTQ inputs (for example, ~4.2 Gbp total, ~14 million reads).
amalgkit_filter_order="fastp_first" # {fastp_first,rrna_first}
remove_amalgkit_fastq_after_completion=1 # Delete per-species amalgkit FASTQ files after downstream completion.
max_assembly_input_fastq_size="30,000,000,000" # Maximum total FASTQ length in bp used for transcriptome assembly.
assembly_method="rnaSPAdes" # {Trinity,rnaSPAdes}
protocol_rna_seq="mixed" # {same,mixed}
kallisto_reference="longest_cds" # {species_cds,longest_transcript,longest_cds,contamination_removed_longest_cds}
orf_aggregation_level="i" # {c,g,i,p}
assembly_cpu_offset=0 # Number of CPU cores reserved from NSLOTS before launching the assembler.
assembly_ram_offset=4 # Amount of RAM in GB reserved from MEM_PER_HOST before launching the assembler.
contamination_removal_rank="domain" # Taxonomic rank for contamination removal. Canonical value is domain; GeneGalleon normalizes tool-specific synonyms automatically.

### End: Modify this block to tailor your analysis ###

delete_tmp_dir=1 # After this run, delete tmp directory created for each job. Set 0 when debugging.

source "${gg_support_dir}/gg_util.sh" # loading utility functions
# Forward config variables (including external overrides) into container environment.
forward_config_vars_to_container_env "${gg_entrypoint_name}" "delete_tmp_dir"
if ! gg_entrypoint_prepare_container_runtime 0; then
  exit 1
fi
gg_entrypoint_activate_container_runtime

gg_entrypoint_enter_workspace
${singularity_command} "${gg_container_image_path}" < "${gg_core_dir}/gg_transcriptome_generation_core.sh"
gg_require_versions_dump "${gg_entrypoint_name}"

echo "$(date): Ending"

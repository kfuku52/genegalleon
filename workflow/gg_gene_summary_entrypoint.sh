#!/usr/bin/env bash

# Scheduler header notes:
# - Keep sections ordered as SLURM -> UGE -> PBS across entrypoints.
# - Update job name, CPU count, memory, walltime, log paths, and working directory together.
# - UGE resources default to SHIROKANE AGE; SLURM/PBS site-specific lines remain examples.

# SLURM
# Common parameters: job name, cores per task, memory per core, walltime, log files, and working directory.
#SBATCH -J gg_gene_summary
#SBATCH -c 4
#SBATCH --mem-per-cpu=16G
#SBATCH -t 1-00:00:00
#SBATCH --output=gg_gene_summary_entrypoint.sh_%j.out
#SBATCH --error=gg_gene_summary_entrypoint.sh_%j.err
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
#$ -l s_vmem=16G
#$ -l ljob

## PBS
# Common parameters: shell, CPU count, total memory, and exported environment.
#PBS -S /bin/bash
#PBS -l ncpus=4
#PBS -l mem=64G
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
gg_entrypoint_name="gg_gene_summary_entrypoint.sh"

### Start: Modify this block to tailor your analysis ###

# Gene-family source
gene_family_source="${gene_family_source:-query2family}" # query2family|orthogroup; selects which gene-family output source to summarize.

# Workflow flags
run_family_completion_summary="${run_family_completion_summary:-1}" # Generate per-family completion/alignment summary tables.
run_presence_absence_summary="${run_presence_absence_summary:-1}" # Generate species x gene-family presence/absence and copy-number matrices.
run_gene_family_database_build="${run_gene_family_database_build:-0}" # Prepare or refresh gg_orthogroup.db for the selected source.
run_csubst_scan_aa_change_summary="${run_csubst_scan_aa_change_summary:-0}" # Generate CSUBST scan AA-change summary tables and plots from gg_orthogroup.db.
run_hgt_candidate_summary="${run_hgt_candidate_summary:-0}" # Summarize GeneRax-first HGT candidate evidence from gg_orthogroup.db.
run_hgt_summary_plots="${run_hgt_summary_plots:-0}" # Generate HGT overview, taxonomy-flow, and per-family HGT tree plots.
run_csubst_site_convergence_summary="${run_csubst_site_convergence_summary:-0}" # Run site-level convergence screening for the selected source.

# Output and shared summary parameters
summary_output_dir="${summary_output_dir:-auto}" # Output directory for source-specific gene summaries.

# Gene-family species-tree presence/absence parameters
presence_absence_species_tree="${presence_absence_species_tree:-auto}" # Species tree path, or auto to use query2family/species_tree outputs.
presence_absence_species_tree_ci="${presence_absence_species_tree_ci:-auto}" # Dated species-tree CI path, or auto to use mcmctree_95CI.nwk when available.
presence_absence_species_tree_support="${presence_absence_species_tree_support:-auto}" # Species tree with numeric branch-support labels, or auto to use species_tree outputs.
presence_absence_busco_table="${presence_absence_busco_table:-auto}" # BUSCO full-table directory or species summary table for the right-side stacked bars, or auto to use workspace outputs.
presence_absence_include_incomplete="${presence_absence_include_incomplete:-0}" # Include query files without stat_branch as NA columns.
presence_absence_heatmap_value="${presence_absence_heatmap_value:-presence}" # presence|copy_number for the PDF/SVG heatmap fill.
presence_absence_plot_width="${presence_absence_plot_width:-7.2}" # Width in inches for the presence/absence summary plot.
presence_absence_max_families="${presence_absence_max_families:-auto}" # Maximum plotted families; auto means all query2family queries and first 100 orthogroups. Use 0/all for no cap.
presence_absence_family_ids="${presence_absence_family_ids:-}" # Optional comma/space-separated family IDs to plot, overriding the default subset.
presence_absence_family_file="${presence_absence_family_file:-}" # Optional file listing family IDs to plot, one ID per line or first column.

# HGT parameters
hgt_summary_use_taxonomy_db="${hgt_summary_use_taxonomy_db:-1}" # Resolve UniProt best-hit taxonomic distances with the local ETE taxonomy DB when available.
hgt_summary_contamination_dir="${hgt_summary_contamination_dir:-}" # Optional directory containing species_cds_contamination_removal_tsv files; empty auto-detects the workspace default.
hgt_summary_taxonomy_flow_rank="${hgt_summary_taxonomy_flow_rank:-phylum}" # Taxonomic rank used to collapse recipient/best-hit lineages in the flow plot.
hgt_summary_taxonomy_flow_max_categories="${hgt_summary_taxonomy_flow_max_categories:-12}" # Maximum recipient and best-hit categories retained before collapsing to Other.
hgt_summary_tree_plot_width="${hgt_summary_tree_plot_width:-24}" # Width in inches for HGT-specific tree plot PDFs.
hgt_summary_promoter_bp="${hgt_summary_promoter_bp:-2000}" # Promoter length used when re-rendering FIMO panels in HGT tree plots.
hgt_summary_fimo_qvalue="${hgt_summary_fimo_qvalue:-0.05}" # FIMO q-value threshold used when re-rendering HGT tree plots.
hgt_summary_output_dir="${hgt_summary_output_dir:-auto}" # HGT output directory; auto is source-specific.

# CSUBST site convergence parameters
csubst_site_arity_range="${csubst_site_arity_range:-2-10}" # Branch-combination arity range to scan.
csubst_site_trait="${csubst_site_trait:-all}" # Trait column name(s) to analyze, or "all".
csubst_site_skip_lower_order="${csubst_site_skip_lower_order:-yes}" # Skip lower-order combinations once higher-order hits are retained.
csubst_site_min_fg_stem_ratio="${csubst_site_min_fg_stem_ratio:-0.5}" # Minimum fraction of foreground branches that must be stem branches.
csubst_site_min_ocn_any2spe="${csubst_site_min_ocn_any2spe:-1.8}" # Minimum OCNany2spe cutoff for candidate combinations.
csubst_site_min_omega_c_any2spe="${csubst_site_min_omega_c_any2spe:-3.0}" # Minimum omegaCany2spe cutoff for candidate combinations.
csubst_site_min_ocn_cod="${csubst_site_min_ocn_cod:-0}" # Minimum OCNCoD cutoff for candidate combinations.
csubst_site_max_candidates_per_arity="${csubst_site_max_candidates_per_arity:-100}" # Maximum number of combinations retained per arity K.
csubst_site_nonsyn_recode="${csubst_site_nonsyn_recode:-}" # Optional CSUBST recoding override; empty uses GG_COMMON_CSUBST_NONSYN_RECODE.
csubst_site_trait_file="${csubst_site_trait_file:-auto}" # Trait table path, or auto to use the workspace default.
csubst_site_orthofinder_dir="${csubst_site_orthofinder_dir:-auto}" # OrthoFinder result directory, or auto to detect it from the workspace.
csubst_site_output_dir="${csubst_site_output_dir:-auto}" # Output directory for convergence results; auto is source-specific.

### End: Modify this block to tailor your analysis ###

source "${gg_support_dir}/gg_util.sh"
gg_apply_registered_env_overrides "${gg_entrypoint_name}"
forward_config_vars_to_container_env "${gg_entrypoint_name}"
gg_export_var_to_container_env_if_set "PYMOL_HEADLESS"
gg_export_var_to_container_env_if_set "QT_QPA_PLATFORM"
if ! gg_entrypoint_prepare_container_runtime 0; then
  exit 1
fi
gg_entrypoint_activate_container_runtime

gg_entrypoint_enter_workspace
gg_runtime_core_script="$(gg_prepare_entrypoint_runtime_snapshot "${gg_entrypoint_name}" "${gg_core_dir}/gg_gene_summary_core.sh")"
gg_run_container_shell_script "${gg_container_image_path}" "${gg_runtime_core_script}"
gg_require_versions_dump "${gg_entrypoint_name}"

echo "$(date): Ending"

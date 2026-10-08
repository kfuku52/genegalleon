#!/usr/bin/env bash
set -euo pipefail

gg_core_bootstrap="/script/support/gg_core_bootstrap.sh"
if [[ ! -s "${gg_core_bootstrap}" ]]; then
  gg_core_bootstrap="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../support/gg_core_bootstrap.sh"
fi
# shellcheck disable=SC1090
source "${gg_core_bootstrap}"
unset gg_core_bootstrap

### Start: Job-supplied configuration ###
# Configuration variables are provided by gg_gene_summary_entrypoint.sh.
### End: Job-supplied configuration ###

gg_bootstrap_core_runtime "${BASH_SOURCE[0]:-$0}" "base" 1 1

hgt_gene_family_mode="${hgt_gene_family_mode:-orthogroup}"

run_hgt_eval="${run_hgt_eval:-1}"
run_hgt_plot="${run_hgt_plot:-1}"
run_hgt_focus="${run_hgt_focus:-1}"
hgt_use_taxonomy_db="${hgt_use_taxonomy_db:-1}"
hgt_contamination_dir="${hgt_contamination_dir:-}"
hgt_taxonomy_flow_rank="${hgt_taxonomy_flow_rank:-phylum}"
hgt_taxonomy_flow_max_categories="${hgt_taxonomy_flow_max_categories:-12}"
hgt_species_tree="${hgt_species_tree:-auto}"
hgt_species_trait="${hgt_species_trait:-auto}"
hgt_focus_event_tsv="${hgt_focus_event_tsv:-auto}"
hgt_focus_event_gene_tsv="${hgt_focus_event_gene_tsv:-auto}"
hgt_focus_filter_audit_tsv="${hgt_focus_filter_audit_tsv:-}"
hgt_focus_context_annotations_tsv="${hgt_focus_context_annotations_tsv:-}"
hgt_focus_require_shared_pfam="${hgt_focus_require_shared_pfam:-1}"
hgt_focus_allow_both_no_pfam="${hgt_focus_allow_both_no_pfam:-0}"
hgt_focus_min_shared_pfam_coverage="${hgt_focus_min_shared_pfam_coverage:-0.5}"
hgt_focus_require_length_ratio="${hgt_focus_require_length_ratio:-1}"
hgt_focus_min_length_ratio="${hgt_focus_min_length_ratio:-0.5}"
hgt_focus_direction_filter="${hgt_focus_direction_filter:-any}"
hgt_focus_species_taxonomy="${hgt_focus_species_taxonomy:-auto}"
case "${hgt_focus_direction_filter}" in
  any|non_arthropoda_to_insecta) ;;
  *) echo "Invalid focused HGT direction filter: ${hgt_focus_direction_filter}" >&2; exit 1 ;;
esac
hgt_transfer_tree_max_edges="${hgt_transfer_tree_max_edges:-200}"
hgt_transfer_arrow_alpha="${hgt_transfer_arrow_alpha:-0.55}"
hgt_tree_width_mm="${hgt_tree_width_mm:-60}"
hgt_promoter_bp="${hgt_promoter_bp:-2000}"
hgt_fimo_qvalue="${hgt_fimo_qvalue:-0.05}"

dir_orthogroup="${hgt_gene_family_dir:-${gg_workspace_output_dir}/orthogroup}"
file_orthogroup_db="${hgt_db_path:-${dir_orthogroup}/gg_orthogroup.db}"
dir_hgt="${hgt_output_dir:-${gg_workspace_output_dir}/hgt}"
default_hgt_contamination_dir="${gg_workspace_output_dir}/species_cds_contamination_removal_tsv"
file_hgt_branch="${dir_hgt}/hgt_branch_candidates.tsv"
file_hgt_gene="${dir_hgt}/hgt_gene_candidates.tsv"
file_hgt_orthogroup="${dir_hgt}/hgt_orthogroup_summary.tsv"
file_hgt_events="${dir_hgt}/hgt_transfer_events.tsv"
file_hgt_event_genes="${dir_hgt}/hgt_transfer_event_genes.tsv"
file_hgt_readme="${dir_hgt}/README.md"
dir_hgt_plot="${dir_hgt}/plots"
dir_hgt_tree_plot="${dir_hgt}/tree_plot"
dir_hgt_tree_input="${dir_hgt}/tree_plot_input"
dir_hgt_tmp="${dir_hgt}/tmp"
dir_hgt_provenance="${dir_hgt}/artifact_provenance"
dir_hgt_trait_focus="${dir_hgt}/trait_focus"
file_hgt_overview_pdf="${dir_hgt_plot}/hgt_branch_overview.pdf"
file_hgt_taxonomy_flow_pdf="${dir_hgt_plot}/hgt_taxonomy_flow.pdf"
file_hgt_transfer_tree_pdf="${dir_hgt_plot}/hgt_transfer_tree.pdf"
file_hgt_transfer_edges="${dir_hgt_plot}/hgt_transfer_edges.tsv"

enable_all_run_flags_for_debug_mode

if [[ "${run_hgt_eval}" != "0" && "${run_hgt_eval}" != "1" ]]; then
  echo "Invalid binary flag value: run_hgt_eval=${run_hgt_eval} (expected 0 or 1)"
  exit 1
fi
if [[ "${run_hgt_plot}" != "0" && "${run_hgt_plot}" != "1" ]]; then
  echo "Invalid binary flag value: run_hgt_plot=${run_hgt_plot} (expected 0 or 1)"
  exit 1
fi
if [[ "${run_hgt_focus}" != "0" && "${run_hgt_focus}" != "1" ]]; then
  echo "Invalid binary flag value: run_hgt_focus=${run_hgt_focus} (expected 0 or 1)"
  exit 1
fi
if [[ "${hgt_use_taxonomy_db}" != "0" && "${hgt_use_taxonomy_db}" != "1" ]]; then
  echo "Invalid binary flag value: hgt_use_taxonomy_db=${hgt_use_taxonomy_db} (expected 0 or 1)"
  exit 1
fi
for hgt_focus_flag in hgt_focus_require_shared_pfam hgt_focus_allow_both_no_pfam hgt_focus_require_length_ratio; do
  if [[ "${!hgt_focus_flag}" != "0" && "${!hgt_focus_flag}" != "1" ]]; then
    echo "Invalid binary flag value: ${hgt_focus_flag}=${!hgt_focus_flag} (expected 0 or 1)" >&2
    exit 1
  fi
done
for hgt_focus_fraction in hgt_focus_min_shared_pfam_coverage hgt_focus_min_length_ratio; do
  if ! python -c 'import math, sys; c = float(sys.argv[1]); sys.exit(not (math.isfinite(c) and 0 <= c <= 1))' "${!hgt_focus_fraction}"; then
    echo "Invalid ${hgt_focus_fraction}: ${!hgt_focus_fraction} (expected a finite fraction from 0 to 1)" >&2
    exit 1
  fi
done
if ! [[ "${hgt_taxonomy_flow_max_categories}" =~ ^[0-9]+$ ]]; then
  echo "Invalid hgt_taxonomy_flow_max_categories: ${hgt_taxonomy_flow_max_categories}"
  exit 1
fi
if ! [[ "${hgt_transfer_tree_max_edges}" =~ ^[0-9]+$ ]]; then
  echo "Invalid hgt_transfer_tree_max_edges: ${hgt_transfer_tree_max_edges}"
  exit 1
fi
if ! python -c 'import math, sys; a = float(sys.argv[1]); sys.exit(not (math.isfinite(a) and 0 <= a <= 1))' "${hgt_transfer_arrow_alpha}"; then
  echo "Invalid hgt_transfer_arrow_alpha: ${hgt_transfer_arrow_alpha} (expected a finite value from 0 to 1)" >&2
  exit 1
fi
if ! [[ "${hgt_tree_width_mm}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "Invalid hgt_tree_width_mm: ${hgt_tree_width_mm}"
  exit 1
fi
if ! [[ "${hgt_promoter_bp}" =~ ^[0-9]+$ ]]; then
  echo "Invalid hgt_promoter_bp: ${hgt_promoter_bp}"
  exit 1
fi
if ! [[ "${hgt_fimo_qvalue}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "Invalid hgt_fimo_qvalue: ${hgt_fimo_qvalue}"
  exit 1
fi
mkdir -p "${dir_hgt}"

resolve_hgt_species_tree() {
  local candidate
  if [[ "${hgt_species_tree}" == "none" ]]; then
    return 1
  fi
  if [[ "${hgt_species_tree}" != "auto" ]]; then
    if [[ -s "${hgt_species_tree}" ]]; then
      printf '%s\n' "${hgt_species_tree}"
      return 0
    fi
    return 1
  fi
  local candidates=(
    "${dir_orthogroup}/parameters/dated_species_tree.pruned.nwk"
    "${dir_orthogroup}/parameters/undated_species_tree.pruned.nwk"
    "${dir_orthogroup}/parameters/dated_species_tree.nwk"
    "${dir_orthogroup}/parameters/undated_species_tree.nwk"
    "${dir_orthogroup}/species_tree.nwk"
    "${gg_workspace_output_dir}/species_tree/species_tree_summary/dated_species_tree.nwk"
    "${gg_workspace_output_dir}/species_tree/species_tree_summary/undated_species_tree.nwk"
    "${gg_workspace_output_dir}/species_tree/mcmctree_main/dated_species_tree.nwk"
  )
  for candidate in "${candidates[@]}"; do
    if [[ -s "${candidate}" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  return 1
}

hgt_species_tree_path=""
if ! hgt_species_tree_path=$(resolve_hgt_species_tree); then
  if [[ "${hgt_species_tree}" != "auto" && "${hgt_species_tree}" != "none" ]]; then
    echo "Warning: HGT transfer-tree species tree was not found: ${hgt_species_tree}" >&2
  elif [[ "${hgt_species_tree}" == "auto" ]]; then
    echo "Warning: No species tree was found for the HGT transfer-tree plot; the edge TSV will still be written." >&2
  fi
fi

hgt_majority_ortholog_prefix() {
  local orthogroup_id=$1
  local gene_tsv=$2
  if [[ ! -s "${gene_tsv}" ]]; then
    return 0
  fi
  awk -F $'\t' -v og="${orthogroup_id}" '
      NR==1 {
        for (i=1; i<=NF; i++) {
          if ($i=="orthogroup") og_col=i
          if ($i=="gene_taxon") tax_col=i
        }
        next
      }
      (og_col>0 && tax_col>0 && $og_col==og && $tax_col!="") {
        counts[$tax_col]++
      }
      END {
        max_count=0
        best=""
        for (taxon in counts) {
          if (counts[taxon] > max_count || (counts[taxon] == max_count && taxon < best)) {
            max_count=counts[taxon]
            best=taxon
          }
        }
        gsub(/ /, "_", best)
        if (best != "") {
          print best "_"
        }
      }
    ' "${gene_tsv}"
}

hgt_count_tree_tips() {
  local stat_branch=$1
  awk -F $'\t' '
      NR==1 {
        for (i=1; i<=NF; i++) {
          if ($i=="so_event") so_col=i
          if ($i=="num_leaf") leaf_col=i
        }
        next
      }
      ((so_col>0 && $so_col=="L") || (leaf_col>0 && $leaf_col=="1")) {n++}
      END {print n+0}
    ' "${stat_branch}"
}

hgt_select_alignment_input() {
  local num_tip=$1
  shift
  local candidate=""
  local candidate_n=0
  for candidate in "$@"; do
    candidate_n=$(gg_count_fasta_records "${candidate}")
    if [[ ${candidate_n} -ge ${num_tip} ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  printf '%s\n' "$1"
}

hgt_materialization_root="${dir_hgt_tmp}/materialized"
hgt_materialization_global_lock="${dir_hgt_tmp}/materialized.cleanup.lock"
hgt_materialization_run_dir=""
hgt_materialization_lock_held=0
hgt_materialization_locking_available=0
if command -v flock >/dev/null 2>&1; then
  hgt_materialization_locking_available=1
fi

hgt_cleanup_materialization_run() {
  local target="${hgt_materialization_run_dir:-}"
  if [[ -n "${target}" && ${hgt_materialization_locking_available} -eq 1 ]]; then
    mkdir -p -- "${dir_hgt_tmp}"
    exec 198> "${hgt_materialization_global_lock}"
    flock -x 198
  fi
  if [[ ${hgt_materialization_lock_held:-0} -eq 1 ]]; then
    flock -u 197 2>/dev/null || true
    exec 197>&-
    hgt_materialization_lock_held=0
  fi
  if [[ -z "${target}" ]]; then
    return 0
  fi
  case "${target}" in
    "${hgt_materialization_root}/"*)
      rm -rf -- "${target}"
      rmdir -- "${hgt_materialization_root}" 2>/dev/null || true
      ;;
    *)
      echo "Refusing to remove unexpected HGT materialization run: ${target}" >&2
      return 1
      ;;
  esac
  hgt_materialization_run_dir=""
  if [[ ${hgt_materialization_locking_available} -eq 1 ]]; then
    flock -u 198 2>/dev/null || true
    exec 198>&-
  fi
}

hgt_cleanup_stale_materialization_runs() {
  local candidate=""
  if [[ ! -d "${hgt_materialization_root}" || -L "${hgt_materialization_root}" ]]; then
    return 0
  fi
  for candidate in "${hgt_materialization_root}"/*; do
    [[ -d "${candidate}" && ! -L "${candidate}" ]] || continue
    exec 196> "${candidate}/.run.lock"
    if flock -n -x 196; then
      case "${candidate}" in
        "${hgt_materialization_root}/"*)
          rm -rf -- "${candidate}"
          ;;
      esac
      flock -u 196 2>/dev/null || true
    fi
    exec 196>&-
  done
  rmdir -- "${hgt_materialization_root}" 2>/dev/null || true
}

hgt_prepare_materialization_run() {
  mkdir -p -- "${dir_hgt_tmp}"
  if [[ ${hgt_materialization_locking_available} -eq 0 ]]; then
    mkdir -p "${hgt_materialization_root}"
    hgt_materialization_run_dir="${hgt_materialization_root}/${GG_JOB_ID:-local}_${GG_ARRAY_TASK_ID:-1}_$$"
    mkdir -p "${hgt_materialization_run_dir}"
    return 0
  fi
  exec 198> "${hgt_materialization_global_lock}"
  flock -x 198
  hgt_cleanup_stale_materialization_runs
  mkdir -p "${hgt_materialization_root}"
  hgt_materialization_run_dir="${hgt_materialization_root}/${GG_JOB_ID:-local}_${GG_ARRAY_TASK_ID:-1}_$$"
  mkdir -p "${hgt_materialization_run_dir}"
  exec 197> "${hgt_materialization_run_dir}/.run.lock"
  flock -x 197
  hgt_materialization_lock_held=1
  flock -u 198
  exec 198>&-
}

hgt_cleanup_family_materialization() {
  local target="${1:-}"
  [[ -n "${target}" ]] || return 0
  case "${target}" in
    "${hgt_materialization_run_dir}/"*)
      rm -rf -- "${target}"
      ;;
    *)
      echo "Refusing to remove unexpected HGT family materialization: ${target}" >&2
      return 1
      ;;
  esac
}

trap hgt_cleanup_materialization_run EXIT

contamination_arg=""
if [[ -n "${hgt_contamination_dir}" ]]; then
  if [[ -d "${hgt_contamination_dir}" ]]; then
    contamination_arg="${hgt_contamination_dir}"
  else
    echo "HGT contamination directory was provided but not found: ${hgt_contamination_dir}" >&2
    exit 1
  fi
elif [[ -d "${default_hgt_contamination_dir}" ]]; then
  contamination_arg="${default_hgt_contamination_dir}"
fi

hgt_taxonomy_db_candidate=""
if [[ ${hgt_use_taxonomy_db} -eq 1 ]]; then
  hgt_taxonomy_db_candidate=$(workspace_taxonomy_dbfile "${gg_workspace_dir}")
fi

hgt_eval_provenance_args=()
gg_artifact_contract_init \
  hgt_eval_provenance_args \
  "hgt_evaluation" \
  "all_gene_families" \
  "${dir_hgt_provenance}/hgt_evaluation.json"
hgt_eval_provenance_args+=(
  --input "gene_family_database=${file_orthogroup_db}"
  --input "hgt_candidate_scorer=${gg_support_dir}/score_hgt_candidates.py"
  --input "scaffold_taxonomy_helper=${gg_support_dir}/scaffold_taxonomy.py"
  --input "species_tree_reader=${gg_support_dir}/hgt_species_tree.py"
  --output "branch_candidates=${file_hgt_branch}"
  --output "gene_candidates=${file_hgt_gene}"
  --output "gene_family_summary=${file_hgt_orthogroup}"
  --parameter "use_taxonomy_db=${hgt_use_taxonomy_db}"
  --parameter "schema_version=4"
)
gg_artifact_add_input_if_present hgt_eval_provenance_args "scaffold_taxonomy" "${gg_workspace_output_dir}/species_scaffold_taxonomy"
gg_artifact_add_input_if_present hgt_eval_provenance_args "species_tree" "${hgt_species_tree_path}"
gg_artifact_add_input_if_present hgt_eval_provenance_args "contamination_tables" "${contamination_arg}"
gg_artifact_add_input_if_present hgt_eval_provenance_args "taxonomy_database" "${hgt_taxonomy_db_candidate}"
gg_artifact_prepare_stage hgt_eval_needs_update run_hgt_eval "${hgt_eval_provenance_args[@]}" || exit $?

if [[ ${run_hgt_eval} -eq 1 && ${hgt_eval_needs_update} -eq 1 ]]; then
  if [[ ! -s "${file_orthogroup_db}" ]]; then
    echo "Skipping HGT evaluation because the orthogroup database was not found: ${file_orthogroup_db}"
  else
    hgt_taxonomy_dbfile=""
    if [[ ${hgt_use_taxonomy_db} -eq 1 ]]; then
      if ensure_ete_taxonomy_db "${gg_workspace_dir}"; then
        hgt_taxonomy_dbfile=$(workspace_taxonomy_dbfile "${gg_workspace_dir}")
        gg_artifact_add_input_if_present hgt_eval_provenance_args "taxonomy_database" "${hgt_taxonomy_dbfile}"
      else
        echo "Warning: Failed to prepare the ETE taxonomy DB. Continuing with best-hit name heuristics only." >&2
      fi
    fi

    python "${gg_support_dir}/score_hgt_candidates.py" \
      --dbpath "${file_orthogroup_db}" \
      --branch_out "${file_hgt_branch}" \
      --gene_out "${file_hgt_gene}" \
      --orthogroup_out "${file_hgt_orthogroup}" \
      --dir_contamination_tsv "${contamination_arg}" \
      --dir_scaffold_taxonomy "${gg_workspace_output_dir}/species_scaffold_taxonomy" \
      --species_tree "${hgt_species_tree_path}" \
      --taxonomy_dbfile "${hgt_taxonomy_dbfile}"
    gg_artifact_record "${hgt_eval_provenance_args[@]}"
  fi
fi

hgt_context_provenance_args=()
gg_artifact_contract_init \
  hgt_context_provenance_args \
  "hgt_transfer_context" \
  "all_gene_families" \
  "${dir_hgt_provenance}/hgt_transfer_context.json"
hgt_context_provenance_args+=(
  --input "branch_candidates=${file_hgt_branch}"
  --input "gene_candidates=${file_hgt_gene}"
  --input "context_summarizer=${gg_support_dir}/summarize_hgt_transfer_context.py"
  --input "scaffold_taxonomy_helper=${gg_support_dir}/scaffold_taxonomy.py"
  --input "species_tree_reader=${gg_support_dir}/hgt_species_tree.py"
  --input-gene-family-store "reconciliations=${dir_orthogroup}"
  --output "transfer_events=${file_hgt_events}"
  --output "transfer_event_genes=${file_hgt_event_genes}"
  --parameter "schema_version=1"
)
gg_artifact_add_input_if_present hgt_context_provenance_args "species_tree" "${hgt_species_tree_path}"
if [[ ${run_hgt_eval} -eq 1 && -s "${file_hgt_branch}" && -s "${file_hgt_gene}" ]]; then
  gg_artifact_prepare_stage hgt_context_needs_update run_hgt_eval "${hgt_context_provenance_args[@]}" || exit $?
  if [[ ${hgt_context_needs_update} -eq 1 ]]; then
    python "${gg_support_dir}/summarize_hgt_transfer_context.py" \
      --branch_tsv "${file_hgt_branch}" \
      --gene_tsv "${file_hgt_gene}" \
      --dir_gene_family "${dir_orthogroup}" \
      --species_tree "${hgt_species_tree_path}" \
      --event_out "${file_hgt_events}" \
      --event_gene_out "${file_hgt_event_genes}"
    gg_artifact_record "${hgt_context_provenance_args[@]}"
  fi
fi

hgt_output_readme_run=1
hgt_output_readme_provenance_args=()
gg_artifact_contract_init \
  hgt_output_readme_provenance_args \
  "hgt_output_readme" \
  "all_candidates" \
  "${dir_hgt_provenance}/hgt_output_readme.json"
hgt_output_readme_provenance_args+=(
  --input "branch_candidates=${file_hgt_branch}"
  --input "gene_candidates=${file_hgt_gene}"
  --input "orthogroup_summary=${file_hgt_orthogroup}"
  --input "readme_generator=${gg_support_dir}/write_hgt_output_readme.py"
  --input "table_schema=${gg_support_dir}/score_hgt_candidates.py"
  --input "scaffold_schema=${gg_support_dir}/scaffold_taxonomy.py"
  --input "transfer_context_schema=${gg_support_dir}/summarize_hgt_transfer_context.py"
  --output "readme=${file_hgt_readme}"
  --parameter "schema_version=2"
)
gg_artifact_add_input_if_present hgt_output_readme_provenance_args "transfer_events" "${file_hgt_events}"
gg_artifact_add_input_if_present hgt_output_readme_provenance_args "transfer_event_genes" "${file_hgt_event_genes}"
gg_artifact_prepare_stage \
  hgt_output_readme_needs_update \
  hgt_output_readme_run \
  "${hgt_output_readme_provenance_args[@]}" || exit $?
if [[ ${hgt_output_readme_run} -eq 1 && ${hgt_output_readme_needs_update} -eq 1 ]]; then
  python "${gg_support_dir}/write_hgt_output_readme.py" \
    --output "${file_hgt_readme}" \
    --branch_tsv "${file_hgt_branch}" \
    --gene_tsv "${file_hgt_gene}" \
    --orthogroup_tsv "${file_hgt_orthogroup}" \
    --event_tsv "${file_hgt_events}" \
    --event_gene_tsv "${file_hgt_event_genes}"
  gg_artifact_record "${hgt_output_readme_provenance_args[@]}"
fi

hgt_summary_plot_provenance_args=()
hgt_species_trait_path=""
if [[ "${hgt_species_trait}" == "auto" ]]; then
  if [[ -s "${gg_workspace_input_dir}/species_trait/species_trait.tsv" ]]; then
    hgt_species_trait_path="${gg_workspace_input_dir}/species_trait/species_trait.tsv"
  fi
elif [[ "${hgt_species_trait}" != "none" ]]; then
  if [[ ! -s "${hgt_species_trait}" ]]; then
    echo "Species trait file not found: ${hgt_species_trait}" >&2
    exit 1
  fi
  hgt_species_trait_path="${hgt_species_trait}"
fi
gg_artifact_contract_init \
  hgt_summary_plot_provenance_args \
  "hgt_summary_plot" \
  "all_candidates" \
  "${dir_hgt_provenance}/hgt_summary_plot.json"
hgt_summary_plot_provenance_args+=(
  --input "plotter=${gg_support_dir}/plot_hgt_summary.py"
  --input "species_tree_reader=${gg_support_dir}/hgt_species_tree.py"
  --input "taxonomy_resolver=${gg_support_dir}/score_hgt_candidates.py"
  --input "branch_candidates=${file_hgt_branch}"
  --input "gene_candidates=${file_hgt_gene}"
  --output "branch_overview=${file_hgt_overview_pdf}"
  --output "taxonomy_flow=${file_hgt_taxonomy_flow_pdf}"
  --output "transfer_tree=${file_hgt_transfer_tree_pdf}"
  --output "transfer_edges=${file_hgt_transfer_edges}"
  --parameter "use_taxonomy_db=${hgt_use_taxonomy_db}"
  --parameter "taxonomy_flow_rank=${hgt_taxonomy_flow_rank}"
  --parameter "taxonomy_flow_max_categories=${hgt_taxonomy_flow_max_categories}"
  --parameter "transfer_tree_max_edges=${hgt_transfer_tree_max_edges}"
  --parameter "transfer_arrow_alpha=${hgt_transfer_arrow_alpha}"
)
gg_artifact_add_input_if_present hgt_summary_plot_provenance_args "taxonomy_database" "${hgt_taxonomy_db_candidate}"
gg_artifact_add_input_if_present hgt_summary_plot_provenance_args "species_tree" "${hgt_species_tree_path}"
gg_artifact_add_input_if_present hgt_summary_plot_provenance_args "species_trait" "${hgt_species_trait_path}"
gg_artifact_add_input_if_present hgt_summary_plot_provenance_args "transfer_events" "${file_hgt_events}"
if [[ -n "${hgt_species_trait_path}" ]]; then
  gg_artifact_add_input_if_present hgt_summary_plot_provenance_args "species_trait_schema" "${hgt_species_trait_path}.schema.json"
  gg_artifact_add_input_if_present hgt_summary_plot_provenance_args "species_trait_metadata" "${hgt_species_trait_path}.metadata.json"
  hgt_summary_plot_provenance_args+=(--input "trait_contract=${gg_support_dir}/species_trait_contract.py" --input "trait_schema=${gg_support_dir}/species_trait_schema.py")
fi
gg_artifact_prepare_stage hgt_summary_plot_needs_update run_hgt_plot "${hgt_summary_plot_provenance_args[@]}" || exit $?

if [[ ${run_hgt_plot} -eq 1 && ${hgt_summary_plot_needs_update} -eq 1 ]]; then
  if [[ ! -s "${file_hgt_branch}" || ! -s "${file_hgt_gene}" ]]; then
    echo "Skipping HGT summary plotting because candidate tables were not found: ${file_hgt_branch}, ${file_hgt_gene}"
  else
    mkdir -p "${dir_hgt_plot}"
    hgt_taxonomy_dbfile=""
    if [[ ${hgt_use_taxonomy_db} -eq 1 ]] && ensure_ete_taxonomy_db "${gg_workspace_dir}"; then
      hgt_taxonomy_dbfile=$(workspace_taxonomy_dbfile "${gg_workspace_dir}")
      gg_artifact_add_input_if_present hgt_summary_plot_provenance_args "taxonomy_database" "${hgt_taxonomy_dbfile}"
    fi
    hgt_transfer_plot_args=()
    if [[ -s "${file_hgt_events}" ]]; then
      hgt_transfer_plot_args+=(--transfer_event_tsv "${file_hgt_events}")
    fi
    python "${gg_support_dir}/plot_hgt_summary.py" \
      --branch_tsv "${file_hgt_branch}" \
      --gene_tsv "${file_hgt_gene}" \
      --overview_pdf "${file_hgt_overview_pdf}" \
      --taxonomy_flow_pdf "${file_hgt_taxonomy_flow_pdf}" \
      --taxonomy_dbfile "${hgt_taxonomy_dbfile}" \
      --flow_rank "${hgt_taxonomy_flow_rank}" \
      --flow_max_categories "${hgt_taxonomy_flow_max_categories}" \
      --transfer_tree_pdf "${file_hgt_transfer_tree_pdf}" \
      --transfer_edges_tsv "${file_hgt_transfer_edges}" \
      --species_tree "${hgt_species_tree_path}" \
      --species_trait "${hgt_species_trait_path}" \
      --transfer_tree_max_edges "${hgt_transfer_tree_max_edges}" \
      --transfer_arrow_alpha "${hgt_transfer_arrow_alpha}" \
      "${hgt_transfer_plot_args[@]}"
    gg_artifact_record "${hgt_summary_plot_provenance_args[@]}"
  fi
fi

# Apply query-Pfam once to the shared input cohort before category-1 trait selection.
if [[ ${run_hgt_focus} -eq 1 && -n "${hgt_species_trait_path}" ]]; then
  hgt_focus_query_taxonomy_path=""
  if [[ -f "${hgt_taxonomy_db_candidate}" ]]; then
    hgt_focus_query_taxonomy_path="${hgt_taxonomy_db_candidate}"
  fi
  hgt_focus_events="${file_hgt_events}"
  hgt_focus_links="${file_hgt_event_genes}"
  [[ "${hgt_focus_event_tsv}" == "auto" ]] || hgt_focus_events="${hgt_focus_event_tsv}"
  [[ "${hgt_focus_event_gene_tsv}" == "auto" ]] || hgt_focus_links="${hgt_focus_event_gene_tsv}"
  if [[ ! -s "${hgt_focus_events}" || ! -s "${hgt_focus_links}" || -z "${hgt_species_tree_path}" ]]; then
    if [[ "${hgt_focus_event_tsv}" != "auto" || "${hgt_focus_event_gene_tsv}" != "auto" ]]; then
      echo "Trait-focused HGT results require the supplied event tables and analysis species tree." >&2
      exit 1
    fi
    echo "Skipping trait-focused HGT results: event context tables or analysis species tree unavailable."
  else
    hgt_focus_provenance_args=()
    hgt_focus_taxonomy_path="${hgt_focus_species_taxonomy}"
    if [[ "${hgt_focus_taxonomy_path}" == "auto" ]]; then
      hgt_focus_taxonomy_path="${gg_workspace_output_dir}/species_taxonomy/species_taxonomy.tsv"
    fi
    if [[ "${hgt_focus_direction_filter}" != "any" && ! -s "${hgt_focus_taxonomy_path}" ]]; then
      echo "Focused species-branch direction filtering requires existing species taxonomy: ${hgt_focus_taxonomy_path}" >&2
      exit 1
    fi
    gg_artifact_contract_init hgt_focus_provenance_args "hgt_trait_focus" "all_category1_targets" \
      "${dir_hgt_provenance}/hgt_trait_focus.json"
    hgt_focus_provenance_args+=(
      --input "events=${hgt_focus_events}"
      --input "event_genes=${hgt_focus_links}"
      --input "species_tree=${hgt_species_tree_path}"
      --input "species_trait=${hgt_species_trait_path}"
      --input "focus_helper=${gg_support_dir}/focus_hgt_traits.py"
      --input "direction_helper=${gg_support_dir}/focus_hgt_direction.py"
      --input "origin_review_helper=${gg_support_dir}/focus_hgt_origin.py"
      --input "pfam_filter=${gg_support_dir}/focus_hgt_pfam.py"
      --input "pfam_background_profile=${gg_support_dir}/focus_hgt_gene_trees.py"
      --input "family_store=${gg_support_dir}/gene_family_output_store.py"
      --input "plotter=${gg_support_dir}/plot_hgt_summary.py"
      --input "species_tree_reader=${gg_support_dir}/hgt_species_tree.py"
      --input "trait_contract=${gg_support_dir}/species_trait_contract.py"
      --input "trait_schema=${gg_support_dir}/species_trait_schema.py"
      --output "result_bundle=${dir_hgt_trait_focus}"
      --parameter "schema_version=1"
      --parameter "plots=${run_hgt_plot}"
      --parameter "transfer_arrow_alpha=${hgt_transfer_arrow_alpha}"
      --parameter "require_shared_pfam=${hgt_focus_require_shared_pfam}"
      --parameter "allow_both_no_pfam=${hgt_focus_allow_both_no_pfam}"
      --parameter "min_shared_pfam_coverage=${hgt_focus_min_shared_pfam_coverage}"
      --parameter "require_length_ratio=${hgt_focus_require_length_ratio}"
      --parameter "min_length_ratio=${hgt_focus_min_length_ratio}"
      --parameter "direction_filter=${hgt_focus_direction_filter}"
    )
    if [[ "${hgt_focus_direction_filter}" != "any" ]]; then
      hgt_focus_provenance_args+=(--input "direction_species_taxonomy=${hgt_focus_taxonomy_path}")
    else
      gg_artifact_add_input_if_present hgt_focus_provenance_args "origin_species_taxonomy" "${hgt_focus_taxonomy_path}"
    fi
    gg_artifact_add_input_if_present hgt_focus_provenance_args "context_mmseqs2_classifications" "${gg_workspace_output_dir}/species_cds_mmseqs2taxonomy"
    gg_artifact_add_input_if_present hgt_focus_provenance_args "context_gene_host_labels" "${gg_workspace_output_dir}/species_scaffold_taxonomy"
    if [[ ${hgt_focus_require_shared_pfam} -eq 1 || ${hgt_focus_require_length_ratio} -eq 1 ]]; then
      hgt_focus_provenance_args+=(
        --input-gene-family-subdir "pfam_query_hits=${dir_orthogroup}::rpsblast"
      )
    fi
    if [[ ${run_hgt_plot} -eq 1 ]]; then
      hgt_focus_provenance_args+=(
        --input "gene_tree_focus_helper=${gg_support_dir}/focus_hgt_gene_trees.py"
        --input "gene_tree_config=${gg_support_dir}/gene_tree_plot_config.py"
        --input "gene_context_annotations=${gg_support_dir}/focus_hgt_context_annotations.py"
        --input "gene_context_taxonomy_schema=${gg_support_dir}/scaffold_taxonomy.py"
        --input "gene_context=${gg_support_dir}/focus_hgt_context.py"
        --input "focused_figures=${gg_support_dir}/focus_hgt_figures.py"
        --input "gene_tree_plotter=${gg_support_dir}/stat_branch2tree_plot.r"
        --input "gene_tree_renderer=${gg_support_dir}/treevis"
        --input-gene-family-subdir "gene_tree_stats=${dir_orthogroup}::stat_branch"
      )
      for hgt_focus_subdir in artifact_provenance synteny rpsblast clipkit orthogroup_extraction_fasta maxalign mafft cds_fasta protein_fasta dated_tree fimo meme promoter_fasta; do
        hgt_focus_provenance_args+=(--input-gene-family-subdir "focus_${hgt_focus_subdir}=${dir_orthogroup}::${hgt_focus_subdir}")
      done
      gg_artifact_add_input_if_present hgt_focus_provenance_args "gff_coordinates" "${gg_workspace_output_dir}/species_gff_info"
      gg_artifact_add_input_if_present hgt_focus_provenance_args "context_query_lineage_database" "${hgt_taxonomy_db_candidate}"
      if [[ -n "${hgt_focus_context_annotations_tsv}" ]]; then
        hgt_focus_provenance_args+=(--input "context_annotations=${hgt_focus_context_annotations_tsv}")
      fi
      gg_artifact_add_input_if_present hgt_focus_provenance_args "filtering_audit" "${hgt_focus_filter_audit_tsv}"
    fi
    gg_artifact_add_input_if_present hgt_focus_provenance_args "trait_schema_input" "${hgt_species_trait_path}.schema.json"
    gg_artifact_add_input_if_present hgt_focus_provenance_args "trait_metadata" "${hgt_species_trait_path}.metadata.json"
    gg_artifact_prepare_stage hgt_focus_needs_update run_hgt_focus "${hgt_focus_provenance_args[@]}" || exit $?
    if [[ ${hgt_focus_needs_update} -eq 1 ]]; then
      hgt_focus_existing_taxonomy="${hgt_focus_taxonomy_path}"
      [[ -f "${hgt_focus_existing_taxonomy}" ]] || hgt_focus_existing_taxonomy=""
      python "${gg_support_dir}/focus_hgt_traits.py" \
        --event_tsv "${hgt_focus_events}" --event_gene_tsv "${hgt_focus_links}" \
        --species_tree "${hgt_species_tree_path}" --species_trait "${hgt_species_trait_path}" \
        --output_dir "${dir_hgt_trait_focus}" --plots "${run_hgt_plot}" \
        --gene_family_root "${dir_orthogroup}" \
        --gff_info_root "${gg_workspace_output_dir}/species_gff_info" \
        --filter_audit_tsv "${hgt_focus_filter_audit_tsv}" \
        --context_annotations_tsv "${hgt_focus_context_annotations_tsv}" \
        --mmseqs2_taxonomy_dir "${gg_workspace_output_dir}/species_cds_mmseqs2taxonomy" \
        --scaffold_taxonomy_dir "${gg_workspace_output_dir}/species_scaffold_taxonomy" \
        --taxonomy_dbfile "${hgt_focus_query_taxonomy_path}" \
        --require_shared_pfam "${hgt_focus_require_shared_pfam}" \
        --allow_both_no_pfam "${hgt_focus_allow_both_no_pfam}" \
        --min_shared_pfam_coverage "${hgt_focus_min_shared_pfam_coverage}" \
        --require_length_ratio "${hgt_focus_require_length_ratio}" \
        --min_length_ratio "${hgt_focus_min_length_ratio}" \
        --direction_filter "${hgt_focus_direction_filter}" \
        --species_taxonomy "${hgt_focus_existing_taxonomy}" \
        --transfer_arrow_alpha "${hgt_transfer_arrow_alpha}"
      gg_artifact_record "${hgt_focus_provenance_args[@]}"
    fi
  fi
fi

if [[ -s "${file_hgt_branch}" && -s "${file_hgt_gene}" ]]; then
  hgt_orthogroups=()
  while IFS= read -r og_id; do
    [[ -z "${og_id}" ]] && continue
    hgt_orthogroups+=("${og_id}")
  done < <(
    awk -F $'\t' '
        NR==1 {
          for (i=1; i<=NF; i++) {
            if ($i=="orthogroup") col=i
          }
          next
        }
        (col>0 && $col!="") {seen[$col]=1}
        END {
          for (k in seen) print k
        }
      ' "${file_hgt_branch}" | sort
  )

  hgt_ggimage_available=-1
  for og_id in "${hgt_orthogroups[@]}"; do
    [[ -z "${og_id}" ]] && continue
    file_hgt_stat_branch="${dir_hgt_tree_input}/${og_id}_hgt_stat.branch.tsv"
    file_hgt_tree_plot="${dir_hgt_tree_plot}/${og_id}_hgt_tree_plot.pdf"
    file_hgt_tree_manifest="${dir_hgt_provenance}/${og_id}.hgt_tree_plot.json"
    hgt_tree_plot_run="${run_hgt_plot}"
    if [[ ${hgt_tree_plot_run} -ne 1 \
      && ! -e "${file_hgt_stat_branch}" \
      && ! -e "${file_hgt_tree_plot}" \
      && ! -e "${file_hgt_tree_manifest}" ]]; then
      continue
    fi

    if [[ -z "${hgt_materialization_run_dir}" ]]; then
      mkdir -p "${dir_hgt_tree_plot}" "${dir_hgt_tree_input}" "${dir_hgt_tmp}"
      hgt_prepare_materialization_run
    fi
    dir_og_input_root="${dir_orthogroup}"
    dir_hgt_materialized="${hgt_materialization_run_dir}/${og_id}"
    hgt_current_materialized=""
    hgt_uses_gene_family_store=0
    if [[ -d "${dir_orthogroup}/.gg_store" || -d "${dir_orthogroup}/.gg_archives" ]]; then
      hgt_uses_gene_family_store=1
      hgt_current_materialized="${dir_hgt_materialized}"
      materialize_args=(
        materialize-family
        --root "${dir_orthogroup}"
        --mode "${hgt_gene_family_mode}"
        --family-id "${og_id}"
        --destination-root "${dir_hgt_materialized}"
        --subdirs "stat_branch,synteny,rpsblast,clipkit,orthogroup_extraction_fasta,maxalign,mafft,protein_fasta,cds_fasta,dated_tree,fimo,meme"
      )
      if [[ "${hgt_gene_family_mode}" == "query2family" ]]; then
        materialize_args+=(--query-dir "${gg_workspace_input_dir}/query_gene")
      fi
      if ! python "${gg_support_dir}/gene_family_output_store.py" "${materialize_args[@]}"; then
        echo "Skipping HGT tree plot for ${og_id}: failed to materialize ZIP-backed inputs." >&2
        hgt_cleanup_family_materialization "${hgt_current_materialized}"
        hgt_current_materialized=""
        continue
      fi
      dir_og_input_root="${dir_hgt_materialized}"
    fi
    file_og_stat_branch="${dir_og_input_root}/stat_branch/${og_id}_stat.branch.tsv"
    file_og_synteny="${dir_og_input_root}/synteny/${og_id}_synteny.tsv"
    file_og_rpsblast="${dir_og_input_root}/rpsblast/${og_id}_rpsblast.tsv"
    file_og_clipkit="${dir_og_input_root}/clipkit/${og_id}_cds.clipkit.fa.gz"
    file_og_orthogroup_extraction_fasta="${dir_og_input_root}/orthogroup_extraction_fasta/${og_id}_orthogroup_extraction.fa.gz"
    file_og_maxalign="${dir_og_input_root}/maxalign/${og_id}_cds.maxalign.fa.gz"
    file_og_mafft="${dir_og_input_root}/mafft/${og_id}_cds.aln.fa.gz"
    file_og_pep_fasta="${dir_og_input_root}/protein_fasta/${og_id}_pep.fa.gz"
    file_og_cds_fasta="${dir_og_input_root}/cds_fasta/${og_id}_cds.fa.gz"
    file_og_dated_tree="${dir_og_input_root}/dated_tree/${og_id}_dated.nwk"
    file_og_fimo="${dir_og_input_root}/fimo/${og_id}_fimo.tsv"
    file_og_meme="${dir_og_input_root}/meme/${og_id}_meme.xml"
    if [[ ! -s "${file_og_stat_branch}" ]]; then
      echo "Skipping HGT tree plot for ${og_id}: stat_branch not found (${file_og_stat_branch})"
      hgt_cleanup_family_materialization "${hgt_current_materialized}"
      hgt_current_materialized=""
      continue
    fi

    hgt_tree_plot_provenance_args=()
    gg_artifact_contract_init \
      hgt_tree_plot_provenance_args \
      "hgt_tree_plot" \
      "${og_id}" \
      "${file_hgt_tree_manifest}"
    hgt_tree_plot_provenance_args+=(
      --input "branch_candidates=${file_hgt_branch}"
      --input "gene_candidates=${file_hgt_gene}"
      --output "annotated_stat_branch=${file_hgt_stat_branch}"
      --output "tree_plot=${file_hgt_tree_plot}"
      --parameter "column_layout=physical-mm-v2-compact-legends"
      --parameter "tree_width_mm=${hgt_tree_width_mm}"
      --parameter "promoter_bp=${hgt_promoter_bp}"
      --parameter "fimo_qvalue=${hgt_fimo_qvalue}"
    )
    hgt_family_input_specs=(
      "stat_branch|stat_branch|${og_id}_stat.branch.tsv|${file_og_stat_branch}"
      "synteny|synteny|${og_id}_synteny.tsv|${file_og_synteny}"
      "rpsblast|rpsblast|${og_id}_rpsblast.tsv|${file_og_rpsblast}"
      "clipkit_alignment|clipkit|${og_id}_cds.clipkit.fa.gz|${file_og_clipkit}"
      "extracted_alignment|orthogroup_extraction_fasta|${og_id}_orthogroup_extraction.fa.gz|${file_og_orthogroup_extraction_fasta}"
      "maxalign_alignment|maxalign|${og_id}_cds.maxalign.fa.gz|${file_og_maxalign}"
      "mafft_alignment|mafft|${og_id}_cds.aln.fa.gz|${file_og_mafft}"
      "protein_fasta|protein_fasta|${og_id}_pep.fa.gz|${file_og_pep_fasta}"
      "cds_fasta|cds_fasta|${og_id}_cds.fa.gz|${file_og_cds_fasta}"
      "dated_tree|dated_tree|${og_id}_dated.nwk|${file_og_dated_tree}"
      "fimo|fimo|${og_id}_fimo.tsv|${file_og_fimo}"
      "meme|meme|${og_id}_meme.xml|${file_og_meme}"
    )
    for hgt_family_input_spec in "${hgt_family_input_specs[@]}"; do
      IFS='|' read -r hgt_input_label hgt_input_subdir hgt_input_name hgt_input_path <<< "${hgt_family_input_spec}"
      [[ -e "${hgt_input_path}" ]] || continue
      if [[ ${hgt_uses_gene_family_store} -eq 1 ]]; then
        hgt_tree_plot_provenance_args+=(
          --input-gene-family-artifact
          "${hgt_input_label}=${dir_orthogroup}::${hgt_input_subdir}::${hgt_input_name}"
        )
      else
        hgt_tree_plot_provenance_args+=(--input "${hgt_input_label}=${hgt_input_path}")
      fi
    done
    gg_artifact_prepare_stage \
      hgt_tree_plot_needs_update \
      hgt_tree_plot_run \
      "${hgt_tree_plot_provenance_args[@]}" || exit $?
    if [[ ${hgt_tree_plot_run} -ne 1 || ${hgt_tree_plot_needs_update} -ne 1 ]]; then
      hgt_cleanup_family_materialization "${hgt_current_materialized}"
      hgt_current_materialized=""
      continue
    fi

    if [[ ${hgt_ggimage_available} -lt 0 ]]; then
      if Rscript -e "if (!requireNamespace('ggimage', quietly=TRUE)) quit(status=1)" > /dev/null 2>&1; then
        hgt_ggimage_available=1
      else
        hgt_ggimage_available=0
        echo "ggimage package is unavailable. Skipping HGT tree plots." >&2
      fi
    fi
    if [[ ${hgt_ggimage_available} -eq 0 ]]; then
      hgt_cleanup_family_materialization "${hgt_current_materialized}"
      hgt_current_materialized=""
      continue
    fi

    ortholog_prefix=$(hgt_majority_ortholog_prefix "${og_id}" "${file_hgt_gene}")
    if [[ -z "${ortholog_prefix}" ]]; then
      ortholog_prefix="HGT_UNRESOLVED_"
    fi

    rm -f -- "${file_hgt_stat_branch}" "${file_hgt_tree_plot}"
    python "${gg_support_dir}/annotate_hgt_tree_plot.py" \
      --stat_branch "${file_og_stat_branch}" \
      --branch_tsv "${file_hgt_branch}" \
      --gene_tsv "${file_hgt_gene}" \
      --orthogroup "${og_id}" \
      --outfile "${file_hgt_stat_branch}"

    num_tip_treeplot=$(hgt_count_tree_tips "${file_hgt_stat_branch}")
    panel_trimmed_aln=$(hgt_select_alignment_input "${num_tip_treeplot}" \
      "${file_og_clipkit}" \
      "${file_og_orthogroup_extraction_fasta}" \
      "${file_og_maxalign}" \
      "${file_og_mafft}" \
      "${file_og_pep_fasta}" \
      "${file_og_cds_fasta}")
    panel_untrimmed_aln=$(hgt_select_alignment_input "${num_tip_treeplot}" \
      "${file_og_orthogroup_extraction_fasta}" \
      "${file_og_mafft}" \
      "${file_og_pep_fasta}" \
      "${file_og_cds_fasta}")

    (
      cd "${dir_hgt_tmp}"
      rm -f -- stat_branch2tree_plot.pdf
      Rscript "${gg_support_dir}/stat_branch2tree_plot.r" \
        --stat_branch="${file_hgt_stat_branch}" \
        --max_delta_intron_present="-0.5" \
        --panel_widths_mm="tree:${hgt_tree_width_mm}" \
        --panel1="tree,bl_rooted,support_unrooted,species,L" \
        --panel2="heatmap,no,abs,_,expression_,Expression" \
        --panel3="pointplot,no,rel,_,expression_" \
        --panel4="heatmap,no,colrel,_,hgt_,HGT evidence (column max=1)" \
        --panel5="cluster_membership,100000" \
        --panel6="synteny,${file_og_synteny},5" \
        --panel7="tiplabel" \
        --panel8="categorical,besthit_lca_rank_display,Hit LCA,-" \
        --panel9="signal_peptide" \
        --panel10="transmembrane_domain" \
        --panel11="intron_number" \
        --panel12="domain,${file_og_rpsblast}" \
        --panel13="alignment,${panel_trimmed_aln},${panel_untrimmed_aln}" \
        --panel14="fimo,${hgt_promoter_bp},${hgt_fimo_qvalue}" \
        --panel15="meme,${file_og_meme}" \
        --panel16="ortholog,${ortholog_prefix},${file_og_dated_tree}" \
        --show_branch_id="yes" \
        --event_method="generax" \
        --species_color_table="PLACEHOLDER" \
        --pie_chart_value_transformation="identity" \
        --long_branch_display="auto" \
        --long_branch_ref_quantile="0.95" \
        --long_branch_detect_ratio="5" \
        --long_branch_cap_ratio="2.5" \
        --long_branch_tail_shrink="0.02" \
        --long_branch_max_fraction="0.1"
      if [[ -s "stat_branch2tree_plot.pdf" ]]; then
        mv_out "stat_branch2tree_plot.pdf" "${file_hgt_tree_plot}"
      else
        echo "Warning: HGT tree plot was not generated for ${og_id}."
      fi
    )
    if [[ -s "${file_hgt_stat_branch}" && -s "${file_hgt_tree_plot}" ]]; then
      gg_artifact_record "${hgt_tree_plot_provenance_args[@]}"
    else
      echo "Warning: HGT tree plot outputs are incomplete for ${og_id}; provenance was not recorded." >&2
    fi
    hgt_cleanup_family_materialization "${hgt_current_materialized}"
    hgt_current_materialized=""
  done
fi

hgt_cleanup_materialization_run
trap - EXIT
echo "$(date): Exiting Singularity environment"

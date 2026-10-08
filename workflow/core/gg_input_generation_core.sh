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
# Configuration variables are provided by gg_input_generation_entrypoint.sh.
### End: Job-supplied configuration ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

gg_bootstrap_core_runtime "${BASH_SOURCE[0]:-$0}" "base" 0 1
# shellcheck disable=SC1090
source "${gg_support_dir}/gg_busco.sh"

config_file="${config_file:-gg_input_generation_entrypoint.sh}"
input_generation_mode="${input_generation_mode:-single}"
run_gene_model_refinement="${run_gene_model_refinement:-0}"
gene_model_refinement_dir="${gene_model_refinement_dir:-}"
gene_model_refinement_policy="${gene_model_refinement_policy:-conserved}"
gene_model_refinement_isoform_adoption="${gene_model_refinement_isoform_adoption:-rna_required}"
gene_model_refinement_mode="${gene_model_refinement_mode:-conservative}"
gene_model_refinement_inputs="${gene_model_refinement_inputs:-}"
gene_model_refinement_edges="${gene_model_refinement_edges:-}"
gene_model_refinement_rescue_dir="${gene_model_refinement_rescue_dir:-}"
gene_model_refinement_rna="${gene_model_refinement_rna:-}"
gene_model_refinement_min_margin="${gene_model_refinement_min_margin:-0.10}"
gene_model_refinement_min_support="${gene_model_refinement_min_support:-2}"
gene_model_refinement_candidate_limit="${gene_model_refinement_candidate_limit:-32}"
gene_model_refinement_padding="${gene_model_refinement_padding:-2000}"
run_gene_model_rescue="${run_gene_model_rescue:-0}"
run_gene_model_rescue_swissprot="${run_gene_model_rescue_swissprot:-1}"
gene_model_rescue_swissprot_dir="${gene_model_rescue_swissprot_dir:-}"
gene_model_rescue_tree="${gene_model_rescue_tree:-auto}"
gene_model_rescue_guide_markers="${gene_model_rescue_guide_markers:-200}"
gene_model_rescue_guide_k="${gene_model_rescue_guide_k:-5}"
gene_model_rescue_guide_sketch_size="${gene_model_rescue_guide_sketch_size:-256}"
gene_model_rescue_guide_occupancy="${gene_model_rescue_guide_occupancy:-0.8}"
gene_model_rescue_guide_minimum_shared="${gene_model_rescue_guide_minimum_shared:-50}"
gene_model_rescue_guide_dir="${gene_model_rescue_guide_dir:-}"
gene_model_rescue_guide_cache="${gene_model_rescue_guide_cache:-}"
gene_model_rescue_dir="${gene_model_rescue_dir:-}"
gene_model_rescue_prediction_cache="${gene_model_rescue_prediction_cache:-}"
gene_model_rescue_common_references="${gene_model_rescue_common_references:-5}"
gene_model_rescue_nearest_references="${gene_model_rescue_nearest_references:-3}"
gene_model_rescue_minimum_busco="${gene_model_rescue_minimum_busco:-90}"
gene_model_rescue_minimum_coverage="${gene_model_rescue_minimum_coverage:-0.95}"
gene_model_rescue_minimum_identity="${gene_model_rescue_minimum_identity:-0.5}"
gene_model_rescue_max_interval="${gene_model_rescue_max_interval:-200000}"
gene_model_rescue_max_intron="${gene_model_rescue_max_intron:-20000}"
gene_model_species_profiles="${gene_model_species_profiles:-}"
gene_model_genome_index_cache="${gene_model_genome_index_cache:-}"
[[ -z "${gene_model_genome_index_cache}" ]] || export GG_GENOME_INDEX_CACHE="${gene_model_genome_index_cache}"
gene_model_rescue_genome_fallback="${gene_model_rescue_genome_fallback:-1}"
gene_model_rescue_max_genome_queries="${gene_model_rescue_max_genome_queries:-20000}"
gene_model_rescue_unanchored_min_species="${gene_model_rescue_unanchored_min_species:-2}"
gene_model_rescue_terminal_max_extension="${gene_model_rescue_terminal_max_extension:-300}"
gene_model_rescue_terminal_max_unaligned_c_overhang="${gene_model_rescue_terminal_max_unaligned_c_overhang:-2}"
gene_model_rescue_gemoma_jar="${gene_model_rescue_gemoma_jar:-}"
gene_model_rescue_gemoma_java="${gene_model_rescue_gemoma_java:-java}"
require_cds="${require_cds:-0}"
require_gff="${require_gff:-0}"
require_genome="${require_genome:-0}"
run_species_busco="${run_species_busco:-1}"
busco_timeout_seconds="${busco_timeout_seconds:-0}"
species_busco_parallel_jobs="${species_busco_parallel_jobs:-auto}"
species_busco_memory_gb_per_job="${species_busco_memory_gb_per_job:-4}"
run_cds_fx2tab="${run_cds_fx2tab:-1}"
run_multispecies_summary="${run_multispecies_summary:-1}"
run_generate_species_trait="${run_generate_species_trait:-0}"
trait_profile="${trait_profile:-none}"
busco_lineage="${busco_lineage:-${GG_COMMON_BUSCO_LINEAGE:-auto}}"
busco_lineage_resolved=""

species_cds_dir="${species_cds_dir:-}"
species_cds_fx2tab_dir="${species_cds_fx2tab_dir:-}"
species_busco_full_dir="${species_busco_full_dir:-}"
species_busco_short_dir="${species_busco_short_dir:-}"
species_gff_dir="${species_gff_dir:-}"
species_genome_dir="${species_genome_dir:-}"
species_summary_output="${species_summary_output:-}"
resolved_manifest_output="${resolved_manifest_output:-}"
species_trait_output="${species_trait_output:-}"
task_plan_output="${task_plan_output:-}"
resume_from_task_plan="${resume_from_task_plan:-}"
resume_from_task_plan_sha256="${resume_from_task_plan_sha256:-}"
resume_from_input_generation_root="${resume_from_input_generation_root:-}"
resume_fallback_task_plan="${resume_fallback_task_plan:-}"
resume_fallback_task_plan_sha256="${resume_fallback_task_plan_sha256:-}"
resume_fallback_input_generation_root="${resume_fallback_input_generation_root:-}"
trait_plan="${trait_plan:-}"
trait_database_sources="${trait_database_sources:-}"
trait_download_dir="${trait_download_dir:-}"
trait_download_timeout="${trait_download_timeout:-120}"
trait_species_source="${trait_species_source:-download_manifest}"
trait_databases="${trait_databases:-auto}"
gbif_api="${gbif_api:-}"
gbif_page_size="${gbif_page_size:-}"
gbif_max_occurrences_per_species="${gbif_max_occurrences_per_species:-}"
gbif_grid_degrees="${gbif_grid_degrees:-}"
gbif_min_match_confidence="${gbif_min_match_confidence:-}"
gbif_max_coordinate_uncertainty_m="${gbif_max_coordinate_uncertainty_m:-}"
gbif_min_distance_from_known_centroid_m="${gbif_min_distance_from_known_centroid_m:-}"
gbif_year_min="${gbif_year_min:-}"
gbif_year_max="${gbif_year_max:-}"
gbif_countries="${gbif_countries:-}"
gbif_include_basis_of_record="${gbif_include_basis_of_record:-}"
gbif_exclude_basis_of_record="${gbif_exclude_basis_of_record:-}"
gbif_include_establishment_means="${gbif_include_establishment_means:-}"
gbif_missing_date="${gbif_missing_date:-}"
gbif_missing_uncertainty="${gbif_missing_uncertainty:-}"
gbif_missing_centroid_distance="${gbif_missing_centroid_distance:-}"
gbif_use_cache="${gbif_use_cache:-}"
gbif_require_complete="${gbif_require_complete:-}"
gbif_occurrence_file="${gbif_occurrence_file:-}"
gbif_taxon_map="${gbif_taxon_map:-}"
gbif_download_metadata="${gbif_download_metadata:-}"
gene_grouping_mode="${gene_grouping_mode:-rescue_overlap}"
gff_repair_mode="${gff_repair_mode:-safe}"
format_contract_version=24

run_species_taxonomy="${run_species_taxonomy:-1}"
taxonomy_species_tree="${taxonomy_species_tree:-auto}"
taxonomy_ranks="${taxonomy_ranks:-all}"
taxonomy_plot_clades="${taxonomy_plot_clades:-0}"
taxonomy_taxid_map="${taxonomy_taxid_map:-}"
taxonomy_taxid_override="${taxonomy_taxid_override:-}"
# Resolve explicit paths before downstream stages change the working directory.
case "${taxonomy_species_tree}" in auto|/*) ;; *) taxonomy_species_tree="${PWD}/${taxonomy_species_tree}" ;; esac
case "${taxonomy_taxid_map}" in ""|/*) ;; *) taxonomy_taxid_map="${PWD}/${taxonomy_taxid_map}" ;; esac
case "${run_species_taxonomy}" in
  0|1) ;;
  *) echo "run_species_taxonomy must be 0 or 1" >&2; exit 1 ;;
esac

enable_all_run_flags_for_debug_mode

case "${trait_profile}" in
  ""|none)
    trait_profile="none"
    ;;
  gift_starter)
    run_generate_species_trait=1
    if [[ -z "${trait_databases}" || "${trait_databases}" == "auto" ]]; then
      trait_databases="gift"
    fi
    ;;
  gbif_distribution)
    run_generate_species_trait=1
    if [[ -z "${trait_databases}" || "${trait_databases}" == "auto" ]]; then
      trait_databases="gbif"
    fi
    ;;
  *)
    echo "Invalid trait_profile: ${trait_profile} (allowed: none|gift_starter|gbif_distribution)"
    exit 1
    ;;
esac

download_manifest_explicit=0
if [[ -n "${download_manifest}" ]]; then
  download_manifest_explicit=1
fi

case "${input_generation_mode}" in
  single|array_prepare|array_worker|array_finalize|rescue_prepare|rescue_synteny|rescue_models|rescue_finalize|refinement_prepare|refinement_catalog|refinement_correspondence|refinement_predict|refinement_finalize) ;;
  *)
    echo "Invalid input_generation_mode: ${input_generation_mode} (allowed: single|array_prepare|array_worker|array_finalize|rescue_prepare|rescue_synteny|rescue_models|rescue_finalize|refinement_prepare|refinement_catalog|refinement_correspondence|refinement_predict|refinement_finalize)"
    exit 1
    ;;
esac
for requirement in require_cds require_gff require_genome; do
  case "${!requirement}" in
    0|1) ;;
    *) echo "${requirement} must be 0 or 1" >&2; exit 1 ;;
  esac
done
if [[ ! "${busco_timeout_seconds}" =~ ^(0|[1-9][0-9]*)$ ]]; then
  echo "busco_timeout_seconds must be 0 or a positive integer." >&2
  exit 1
fi

case "${provider}" in
  refseq|genbank)
    echo "Provider '${provider}' is treated as alias of 'ncbi'."
    provider="ncbi"
    ;;
  all|ensembl|ensemblplants|ensemblmetazoa|ensemblprotists|phycocosm|phytozome|ncbi|coge|cngb|gwh|flybase|wormbase|vectorbase|fernbase|veupathdb|dictybase|insectbase|direct|local) ;;
  *)
    echo "Invalid provider: ${provider} (allowed: all|ensembl|ensemblplants|ensemblmetazoa|ensemblprotists|phycocosm|phytozome|ncbi|coge|cngb|gwh|flybase|wormbase|vectorbase|fernbase|veupathdb|dictybase|insectbase|direct|local)"
    exit 1
    ;;
esac

for binary_flag_name in \
  run_format_inputs \
  run_validate_inputs \
  run_cds_fx2tab \
  run_species_busco \
  run_gene_model_refinement \
  run_gene_model_rescue \
  run_gene_model_rescue_swissprot \
  gene_model_rescue_genome_fallback \
  run_multispecies_summary \
  run_generate_species_trait \
  strict \
  overwrite \
  download_only \
  dry_run
do
  binary_flag_value="${!binary_flag_name}"
  if [[ "${binary_flag_value}" != "0" && "${binary_flag_value}" != "1" ]]; then
    echo "Invalid binary flag value: ${binary_flag_name}=${binary_flag_value} (expected 0 or 1)"
    exit 1
  fi
done

if ! [[ "${download_timeout}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "Invalid download_timeout: ${download_timeout}"
  exit 1
fi
if ! [[ "${trait_download_timeout}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "Invalid trait_download_timeout: ${trait_download_timeout}"
  exit 1
fi
if [[ "${trait_species_source}" != "download_manifest" && "${trait_species_source}" != "species_cds" ]]; then
  echo "Invalid trait_species_source: ${trait_species_source} (allowed: download_manifest|species_cds)"
  exit 1
fi
case "${gene_grouping_mode}" in
  strict|rescue_overlap) ;;
  *)
    echo "Invalid gene_grouping_mode: ${gene_grouping_mode} (allowed: strict|rescue_overlap)"
    exit 1
    ;;
esac
case "${gff_repair_mode}" in
  off|safe|strict) ;;
  *)
    echo "Invalid gff_repair_mode: ${gff_repair_mode} (allowed: off|safe|strict)"
    exit 1
    ;;
esac
if [[ "${input_generation_mode}" != "single" && ${dry_run} -eq 1 ]]; then
  echo "dry_run=1 is not supported with input_generation_mode=${input_generation_mode}"
  exit 1
fi
if [[ "${input_generation_mode}" != "single" && ${download_only} -eq 1 ]]; then
  echo "download_only=1 is only supported in input_generation_mode=single"
  exit 1
fi
if [[ "${input_generation_mode}" == array_* && ${run_format_inputs} -ne 1 ]]; then
  echo "run_format_inputs must be 1 when input_generation_mode=${input_generation_mode}"
  exit 1
fi

input_generation_root="${gg_workspace_output_dir}/input_generation"
gene_model_refinement_dir="${gene_model_refinement_dir:-${input_generation_root}/gene_model_refinement}"
case "${gene_model_refinement_dir}" in /*) ;; *) gene_model_refinement_dir="${PWD}/${gene_model_refinement_dir}" ;; esac
gene_model_rescue_dir="${gene_model_rescue_dir:-${input_generation_root}/gene_model_rescue}"
case "${gene_model_rescue_dir}" in /*) ;; *) gene_model_rescue_dir="${PWD}/${gene_model_rescue_dir}" ;; esac
case "${gene_model_rescue_swissprot_dir}" in ""|/*) ;; *) gene_model_rescue_swissprot_dir="${PWD}/${gene_model_rescue_swissprot_dir}" ;; esac
gene_model_rescue_guide_dir="${gene_model_rescue_guide_dir:-${gene_model_rescue_dir}/guide_tree}"
gene_model_rescue_guide_cache="${gene_model_rescue_guide_cache:-$(dirname "${gene_model_rescue_dir}")/busco_guide_sketch_cache}"
case "${gene_model_rescue_guide_dir}" in /*) ;; *) gene_model_rescue_guide_dir="${PWD}/${gene_model_rescue_guide_dir}" ;; esac
case "${gene_model_rescue_guide_cache}" in /*) ;; *) gene_model_rescue_guide_cache="${PWD}/${gene_model_rescue_guide_cache}" ;; esac
case "${gene_model_rescue_tree}" in auto|/*) ;; *) gene_model_rescue_tree="${PWD}/${gene_model_rescue_tree}" ;; esac
case "${gene_model_rescue_gemoma_jar}" in ""|/*) ;; *) gene_model_rescue_gemoma_jar="${PWD}/${gene_model_rescue_gemoma_jar}" ;; esac
input_generation_tmp_root="${input_generation_root}/tmp"
# Shared across arrays and preserved when task scratch directories are cleaned.
export download_limit_dir="${download_limit_dir:-${gg_workspace_dir}/.gg_cache/input_download_limits}"
input_generation_provenance_dir="${input_generation_root}/artifact_provenance"
download_tmp_root=$(gg_task_tmp_path "${input_generation_tmp_root}") || exit 1
ensure_dir "${download_tmp_root}"
dir_species_summary_shards="${input_generation_tmp_root}/species_summary_shards"
dir_task_stats_shards="${input_generation_tmp_root}/task_stats_shards"
dir_task_meta_shards="${input_generation_tmp_root}/task_meta_shards"
file_busco_lineage_resolved="${input_generation_tmp_root}/busco_lineage.resolved.txt"

if [[ -z "${download_dir}" ]]; then
  download_dir="${input_generation_tmp_root}/input_download_cache"
fi
if [[ -z "${summary_output}" ]]; then
  summary_output="${input_generation_root}/gg_input_generation_runs.tsv"
fi

if [[ -z "${species_cds_dir}" ]]; then
  species_cds_dir="${input_generation_root}/species_cds"
fi
if [[ -z "${species_cds_fx2tab_dir}" ]]; then
  species_cds_fx2tab_dir="${input_generation_root}/species_cds_fx2tab"
fi
if [[ -z "${species_busco_full_dir}" ]]; then
  species_busco_full_dir="${input_generation_root}/species_cds_busco_full"
fi
if [[ -z "${species_busco_short_dir}" ]]; then
  species_busco_short_dir="${input_generation_root}/species_cds_busco_short"
fi
if [[ -z "${species_gff_dir}" ]]; then
  species_gff_dir="${input_generation_root}/species_gff"
fi
if [[ -z "${species_genome_dir}" ]]; then
  species_genome_dir="${input_generation_root}/species_genome"
fi
if [[ -z "${species_summary_output}" ]]; then
  species_summary_output="${input_generation_root}/gg_input_generation_species.tsv"
fi
if [[ -z "${resolved_manifest_output}" ]]; then
  resolved_manifest_output="${input_generation_root}/download_plan.resolved.tsv"
fi
if [[ -z "${species_trait_output}" ]]; then
  species_trait_output="${gg_workspace_input_dir}/species_trait/species_trait.tsv"
fi
if [[ -z "${task_plan_output}" ]]; then
  task_plan_output="${input_generation_tmp_root}/task_plan.json"
fi
if [[ -z "${trait_plan}" ]]; then
  trait_plan="${gg_workspace_input_dir}/input_generation/trait_plan.tsv"
fi
if [[ -z "${trait_database_sources}" ]]; then
  trait_database_sources="${gg_workspace_input_dir}/input_generation/trait_database_sources.tsv"
fi
if [[ -z "${trait_download_dir}" ]]; then
  trait_download_dir="${gg_workspace_downloads_dir}/trait_datasets"
fi

discover_input_generation_manifests() {
  local manifest_dir="${gg_workspace_input_dir}/input_generation"
  if [[ ! -d "${manifest_dir}" ]]; then
    return 0
  fi
  python - "${manifest_dir}" <<'PY'
import csv
import sys
from pathlib import Path


manifest_dir = Path(sys.argv[1])
excluded_markers = ("audit", "resolved", "summary", "template", "snapshot", "options")


def delimited_manifest_shape(path):
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open("rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader, [])
        has_data = any(any(str(value).strip() for value in row) for row in reader)
    return header, has_data


def xlsx_manifest_shape(path):
    from openpyxl import load_workbook

    workbook = load_workbook(path, read_only=True, data_only=True)
    try:
        sheet = workbook.active
        rows = sheet.iter_rows(values_only=True)
        header = list(next(rows, ()) or ())
        has_data = any(any(value is not None and str(value).strip() for value in row) for row in rows)
        return header, has_data
    finally:
        workbook.close()


for path in sorted(manifest_dir.iterdir()):
    if not path.is_file() or path.stat().st_size == 0:
        continue
    suffix = path.suffix.lower()
    if suffix not in (".csv", ".tsv", ".xlsx"):
        continue
    stem = path.stem.lower()
    if not stem.startswith("download_plan"):
        continue
    if any(marker in stem for marker in excluded_markers):
        continue
    try:
        if suffix == ".xlsx":
            header, has_data = xlsx_manifest_shape(path)
        else:
            header, has_data = delimited_manifest_shape(path)
    except Exception:
        print(path.resolve())
        continue
    normalized_header = {str(value or "").strip().lower() for value in header}
    if has_data and {"provider", "id"}.issubset(normalized_header):
        print(path.resolve())
PY
}

if [[ -z "${download_manifest}" && "${input_generation_mode}" != array_worker && "${input_generation_mode}" != array_finalize && "${input_generation_mode}" != rescue_* && "${input_generation_mode}" != refinement_* ]]; then
  default_download_manifests=()
  while IFS= read -r discovered_manifest; do
    [[ -n "${discovered_manifest}" ]] || continue
    default_download_manifests+=( "${discovered_manifest}" )
  done < <(discover_input_generation_manifests)
  unset discovered_manifest
  if [[ ${#default_download_manifests[@]} -gt 1 ]]; then
    echo "Multiple input-generation download manifests were discovered. Set download_manifest explicitly:"
    printf '  %s\n' "${default_download_manifests[@]}"
    exit 1
  fi
  if [[ ${#default_download_manifests[@]} -eq 1 ]]; then
    download_manifest="${default_download_manifests[0]}"
    echo "Auto-selected download_manifest: ${download_manifest}"
  fi
  unset default_download_manifests
fi

file_multispecies_summary="${input_generation_root}/annotation_summary/annotation_summary.tsv"

num_species_cds=""
num_species_cds_fx2tab=""
num_species_gff=""
num_species_genome=""
num_species_busco_full=""
num_species_busco_short=""
num_species_trait=""
num_trait_columns=""
cds_sequences_before=""
cds_sequences_after=""
cds_first_sequence_name=""
stage_format_status="not_run"
stage_validate_status="not_run"
stage_cds_fx2tab_status="not_run"
stage_species_busco_status="not_run"
stage_multispecies_summary_status="not_run"
stage_trait_status="not_run"
write_run_summary_on_exit=1
cleanup_input_generation_tmp=0

run_started_epoch=$(date +%s)
run_started_iso=$(date -u '+%Y-%m-%dT%H:%M:%SZ')

sanitize_tsv_value() {
  local value="$1"
  value=$(printf '%s' "${value}" | tr '\t\r\n' '   ')
  printf '%s' "${value}"
}

manifest_data_row_count() {
  local manifest_path="$1"
  if [[ ! -s "${manifest_path}" ]]; then
    echo 0
    return 0
  fi
  if [[ "${manifest_path##*.}" == "xlsx" ]]; then
    python - "${manifest_path}" <<'PY'
import sys
from pathlib import Path

path = Path(sys.argv[1])
try:
    from openpyxl import load_workbook
except Exception:
    print(1)
    raise SystemExit(0)

try:
    workbook = load_workbook(path, read_only=True, data_only=True)
    sheet = workbook.active
    row_iter = sheet.iter_rows(values_only=True)
    next(row_iter, None)
    count = 0
    for row in row_iter:
        if row is None:
            continue
        nonempty = False
        for value in row:
            if value is None:
                continue
            if str(value).strip() != "":
                nonempty = True
                break
        if nonempty:
            count += 1
    print(count)
except Exception:
    print(1)
PY
    return 0
  fi
  awk '
    NR == 1 {next}
    /^[[:space:]]*$/ {next}
    /^[[:space:]]*#/ {next}
    {n++}
    END {print n+0}
  ' "${manifest_path}"
}

ensure_selected_download_manifest_has_rows() {
  local manifest_rows=""
  if [[ -z "${download_manifest}" ]]; then
    return 0
  fi
  manifest_rows="$(manifest_data_row_count "${download_manifest}")"
  if [[ ${manifest_rows} -gt 0 ]]; then
    return 0
  fi
  echo "Download manifest has no data rows: ${download_manifest}"
  if [[ ${download_manifest_explicit} -eq 1 ]]; then
    echo "Explicitly provided download_manifest is empty."
    return 1
  fi
  echo "Ignoring empty auto-selected download_manifest."
  download_manifest=""
  return 0
}

read_stats_json_field() {
  local stats_file="$1"
  local field_name="$2"
  python - "${stats_file}" "${field_name}" <<'PY'
import json
import sys

path = sys.argv[1]
field = sys.argv[2]
with open(path, "rt", encoding="utf-8") as handle:
    data = json.load(handle)
value = data.get(field, "")
if value is None:
    value = ""
print(value)
PY
}

read_stats_json_fields() {
  local stats_file="$1"
  shift
  local fields_file
  local field_name
  local field_value
  local -a field_names=("$@")
  fields_file=$(mktemp "${download_tmp_root}/json-fields.XXXXXX") || return $?
  if ! python - "${stats_file}" "${field_names[@]}" > "${fields_file}" <<'PYFIELDS'
import json
import sys
with open(sys.argv[1], encoding="utf-8") as handle:
    data = json.load(handle)
values = ["" if data.get(field) is None else str(data.get(field, "")).rstrip("\n")
          for field in (item.split("=", 1)[-1] for item in sys.argv[2:])]
if any("\0" in value for value in values):
    raise ValueError("NUL is unsupported in shell metadata fields")
sys.stdout.buffer.write(b"".join(value.encode("utf-8") + b"\0" for value in values))
PYFIELDS
  then
    rm -f -- "${fields_file}"
    return 1
  fi
  for field_name in "${field_names[@]}"; do
    IFS= read -r -d '' field_value || { rm -f -- "${fields_file}"; return 1; }
    printf -v "${field_name%%=*}" '%s' "${field_value}"
  done < "${fields_file}"
  rm -f -- "${fields_file}"
}

input_generation_effective_input_dir_path() {
  if [[ -n "${input_dir}" ]]; then
    printf '%s\n' "${input_dir}"
    return 0
  fi
  if [[ -z "${download_manifest}" ]]; then
    printf '%s\n' ""
    return 0
  fi
  if [[ "${provider}" == "all" ]]; then
    printf '%s\n' "${download_dir}"
    return 0
  fi
  case "${provider}" in
    ensembl)
      printf '%s\n' "${download_dir}/Ensembl/original_files"
      ;;
    ensemblplants)
      printf '%s\n' "${download_dir}/20230216_EnsemblPlants/original_files"
      ;;
    ensemblmetazoa)
      printf '%s\n' "${download_dir}/EnsemblMetazoa/original_files"
      ;;
    ensemblprotists)
      printf '%s\n' "${download_dir}/EnsemblProtists/original_files"
      ;;
    phycocosm)
      printf '%s\n' "${download_dir}/PhycoCosm/species_wise_original"
      ;;
    phytozome)
      printf '%s\n' "${download_dir}/Phytozome/species_wise_original"
      ;;
    ncbi|refseq|genbank)
      printf '%s\n' "${download_dir}/NCBI_Genome/species_wise_original"
      ;;
    coge)
      printf '%s\n' "${download_dir}/CoGe/species_wise_original"
      ;;
    cngb)
      printf '%s\n' "${download_dir}/CNGB/species_wise_original"
      ;;
    gwh)
      printf '%s\n' "${download_dir}/GWH/species_wise_original"
      ;;
    flybase)
      printf '%s\n' "${download_dir}/FlyBase/species_wise_original"
      ;;
    wormbase)
      printf '%s\n' "${download_dir}/WormBase/species_wise_original"
      ;;
    vectorbase)
      printf '%s\n' "${download_dir}/VectorBase/species_wise_original"
      ;;
    fernbase)
      printf '%s\n' "${download_dir}/FernBase/species_wise_original"
      ;;
    veupathdb)
      printf '%s\n' "${download_dir}/VEuPathDB/species_wise_original"
      ;;
    dictybase)
      printf '%s\n' "${download_dir}/dictyBase/species_wise_original"
      ;;
    insectbase)
      printf '%s\n' "${download_dir}/InsectBase/species_wise_original"
      ;;
    direct)
      printf '%s\n' "${download_dir}/Direct/species_wise_original"
      ;;
    local)
      printf '%s\n' "${download_dir}/Local/species_wise_original"
      ;;
    *)
      printf '%s\n' ""
      ;;
  esac
}

count_nonhidden_matching_files() {
  local search_dir=$1
  local pattern=$2
  if [[ ! -d "${search_dir}" ]]; then
    echo 0
    return 0
  fi
  find "${search_dir}" -maxdepth 1 -type f ! -name '.*' -name "${pattern}" | wc -l | awk '{print $1}'
}

list_nonhidden_matching_files() {
  local search_dir=$1
  shift
  if [[ ! -d "${search_dir}" ]]; then
    return 0
  fi
  find "${search_dir}" -maxdepth 1 -type f ! -name '.*' \( "$@" \) | sort
}

task_plan_task_count() {
  local plan_file=$1
  if [[ ! -s "${plan_file}" ]]; then
    echo 0
    return 0
  fi
  python - "${plan_file}" <<'PY'
import json
import sys

with open(sys.argv[1], "rt", encoding="utf-8") as handle:
    payload = json.load(handle)
print(int(payload.get("task_count", len(payload.get("tasks") or []))))
PY
}

task_plan_species_names() {
  local plan_file=$1
  python - "${plan_file}" <<'PY'
import json
import sys

with open(sys.argv[1], "rt", encoding="utf-8") as handle:
    payload = json.load(handle)
for task in payload.get("tasks") or []:
    species = str(task.get("species_prefix") or "").strip()
    if species:
        print(species)
PY
}

resolve_busco_lineage_from_task_plan() {
  local plan_file=$1
  local -a species=()
  while IFS= read -r species_name; do
    [[ -n "${species_name}" ]] || continue
    species+=( "${species_name}" )
  done < <(task_plan_species_names "${plan_file}")
  if [[ ${#species[@]} -eq 0 ]]; then
    echo "No species names were discovered in task plan: ${plan_file}" >&2
    return 1
  fi
  if ! busco_lineage_resolved=$(gg_resolve_busco_lineage "${gg_workspace_dir}" "${busco_lineage}" "${species[@]}"); then
    echo "Failed to resolve BUSCO lineage from request: ${busco_lineage}" >&2
    return 1
  fi
  printf '%s\n' "${busco_lineage_resolved}" > "${file_busco_lineage_resolved}"
  echo "Resolved BUSCO lineage for task plan (${#species[@]} species): ${busco_lineage_resolved}"
}

ensure_shared_busco_lineage_ready() {
  local source_plan=${1:-}
  if [[ -n "${busco_lineage_resolved}" ]]; then
    return 0
  fi
  if [[ -s "${file_busco_lineage_resolved}" ]]; then
    busco_lineage_resolved=$(tr -d '\r\n' < "${file_busco_lineage_resolved}")
    if [[ -n "${busco_lineage_resolved}" ]]; then
      return 0
    fi
  fi
  if [[ -n "${source_plan}" ]]; then
    resolve_busco_lineage_from_task_plan "${source_plan}"
    return 0
  fi
  return 1
}

write_gg_input_generation_summary_on_exit() {
  local exit_code=${1:-$?}
  local run_ended_iso
  local run_duration_sec
  local header
  local expected_header_line
  local existing_header_line
  local legacy_output
  local legacy_stamp
  local row

  if [[ ${write_run_summary_on_exit} -ne 1 ]]; then
    return 0
  fi

  run_ended_iso=$(date -u '+%Y-%m-%dT%H:%M:%SZ')
  run_duration_sec=$(( $(date +%s) - run_started_epoch ))

  if [[ ${exit_code} -ne 0 ]]; then
    if [[ "${stage_format_status}" == "running" ]]; then stage_format_status="failed"; fi
    if [[ "${stage_validate_status}" == "running" ]]; then stage_validate_status="failed"; fi
    if [[ "${stage_cds_fx2tab_status}" == "running" ]]; then stage_cds_fx2tab_status="failed"; fi
    if [[ "${stage_species_busco_status}" == "running" ]]; then stage_species_busco_status="failed"; fi
    if [[ "${stage_multispecies_summary_status}" == "running" ]]; then stage_multispecies_summary_status="failed"; fi
    if [[ "${stage_trait_status}" == "running" ]]; then stage_trait_status="failed"; fi
  fi

  ensure_parent_dir "${summary_output}"
  header="started_utc\tended_utc\tduration_sec\texit_code\tprovider\tinput_generation_mode\trun_format_inputs\trun_validate_inputs\trun_cds_fx2tab\trun_species_busco\trun_multispecies_summary\trun_generate_species_trait\tbusco_lineage\tbusco_lineage_resolved\tstrict\toverwrite\tdownload_only\tdry_run\tdownload_timeout\tinput_dir\tdownload_manifest\tdownload_dir\ttask_plan_output\tspecies_cds_dir\tspecies_cds_fx2tab_dir\tspecies_busco_full_dir\tspecies_busco_short_dir\tspecies_gff_dir\tspecies_genome_dir\tspecies_trait_output\tfile_multispecies_summary\tnum_species_cds\tnum_species_cds_fx2tab\tnum_species_gff\tnum_species_genome\tnum_species_busco_full\tnum_species_busco_short\tnum_species_trait\tnum_trait_columns\tcds_sequences_before\tcds_sequences_after\tcds_first_sequence_name\tstage_format_status\tstage_validate_status\tstage_cds_fx2tab_status\tstage_species_busco_status\tstage_multispecies_summary_status\tstage_trait_status\tconfig_file"
  expected_header_line=$(printf '%b' "${header}")
  if [[ -s "${summary_output}" ]]; then
    existing_header_line=$(head -n 1 "${summary_output}" || true)
    if [[ "${existing_header_line}" != "${expected_header_line}" ]]; then
      legacy_stamp=$(date -u '+%Y%m%dT%H%M%SZ')
      legacy_output="${summary_output%.tsv}.legacy.${legacy_stamp}.tsv"
      if [[ "${legacy_output}" == "${summary_output}" ]]; then
        legacy_output="${summary_output}.legacy.${legacy_stamp}"
      fi
      mv -- "${summary_output}" "${legacy_output}"
      echo "Summary header format changed. Archived previous summary: ${legacy_output}"
    fi
  fi
  if [[ ! -s "${summary_output}" ]]; then
    printf '%b\n' "${header}" > "${summary_output}"
  fi
  row="$(sanitize_tsv_value "${run_started_iso}")"
  row="${row}\t$(sanitize_tsv_value "${run_ended_iso}")"
  row="${row}\t$(sanitize_tsv_value "${run_duration_sec}")"
  row="${row}\t$(sanitize_tsv_value "${exit_code}")"
  row="${row}\t$(sanitize_tsv_value "${provider}")"
  row="${row}\t$(sanitize_tsv_value "${input_generation_mode}")"
  row="${row}\t$(sanitize_tsv_value "${run_format_inputs}")"
  row="${row}\t$(sanitize_tsv_value "${run_validate_inputs}")"
  row="${row}\t$(sanitize_tsv_value "${run_cds_fx2tab}")"
  row="${row}\t$(sanitize_tsv_value "${run_species_busco}")"
  row="${row}\t$(sanitize_tsv_value "${run_multispecies_summary}")"
  row="${row}\t$(sanitize_tsv_value "${run_generate_species_trait}")"
  row="${row}\t$(sanitize_tsv_value "${busco_lineage}")"
  row="${row}\t$(sanitize_tsv_value "${busco_lineage_resolved}")"
  row="${row}\t$(sanitize_tsv_value "${strict}")"
  row="${row}\t$(sanitize_tsv_value "${overwrite}")"
  row="${row}\t$(sanitize_tsv_value "${download_only}")"
  row="${row}\t$(sanitize_tsv_value "${dry_run}")"
  row="${row}\t$(sanitize_tsv_value "${download_timeout}")"
  row="${row}\t$(sanitize_tsv_value "${input_dir}")"
  row="${row}\t$(sanitize_tsv_value "${download_manifest}")"
  row="${row}\t$(sanitize_tsv_value "${download_dir}")"
  row="${row}\t$(sanitize_tsv_value "${task_plan_output}")"
  row="${row}\t$(sanitize_tsv_value "${species_cds_dir}")"
  row="${row}\t$(sanitize_tsv_value "${species_cds_fx2tab_dir}")"
  row="${row}\t$(sanitize_tsv_value "${species_busco_full_dir}")"
  row="${row}\t$(sanitize_tsv_value "${species_busco_short_dir}")"
  row="${row}\t$(sanitize_tsv_value "${species_gff_dir}")"
  row="${row}\t$(sanitize_tsv_value "${species_genome_dir}")"
  row="${row}\t$(sanitize_tsv_value "${species_trait_output}")"
  row="${row}\t$(sanitize_tsv_value "${file_multispecies_summary}")"
  row="${row}\t$(sanitize_tsv_value "${num_species_cds}")"
  row="${row}\t$(sanitize_tsv_value "${num_species_cds_fx2tab}")"
  row="${row}\t$(sanitize_tsv_value "${num_species_gff}")"
  row="${row}\t$(sanitize_tsv_value "${num_species_genome}")"
  row="${row}\t$(sanitize_tsv_value "${num_species_busco_full}")"
  row="${row}\t$(sanitize_tsv_value "${num_species_busco_short}")"
  row="${row}\t$(sanitize_tsv_value "${num_species_trait}")"
  row="${row}\t$(sanitize_tsv_value "${num_trait_columns}")"
  row="${row}\t$(sanitize_tsv_value "${cds_sequences_before}")"
  row="${row}\t$(sanitize_tsv_value "${cds_sequences_after}")"
  row="${row}\t$(sanitize_tsv_value "${cds_first_sequence_name}")"
  row="${row}\t$(sanitize_tsv_value "${stage_format_status}")"
  row="${row}\t$(sanitize_tsv_value "${stage_validate_status}")"
  row="${row}\t$(sanitize_tsv_value "${stage_cds_fx2tab_status}")"
  row="${row}\t$(sanitize_tsv_value "${stage_species_busco_status}")"
  row="${row}\t$(sanitize_tsv_value "${stage_multispecies_summary_status}")"
  row="${row}\t$(sanitize_tsv_value "${stage_trait_status}")"
  row="${row}\t$(sanitize_tsv_value "${config_file}")"
  printf '%b\n' "${row}" >> "${summary_output}"
}

array_lock_paths=()
array_lock_tokens=()
array_lock_modes=()

input_generation_lock() {
  local path=$1 mode=$2 token
  local -a contention_args=(--nonblocking)
  # Shared readers briefly serialize while registering their ownership. Retry
  # that gate contention instead of rejecting a compatible parallel worker.
  if [[ ${mode} == shared ]]; then
    contention_args=(--timeout 30)
  fi
  token=$(python "${gg_support_dir}/shared_namespace_lock.py" "acquire-${mode}" "${path}" --owner-pid "$$" "${contention_args[@]}") || {
    echo "Input-generation location is already in use: ${path}" >&2
    return 1
  }
  array_lock_paths+=("${path}")
  array_lock_tokens+=("${token}")
  array_lock_modes+=("${mode}")
}

input_generation_on_exit() {
  local status=$? index interrupted=0
  trap - EXIT
  (( status < 128 )) || interrupted=1
  write_gg_input_generation_summary_on_exit "${status}" || status=1
  if [[ ${interrupted} -eq 1 ]]; then
    echo "Interrupted input generation: retaining namespace locks until the job and its children are reconciled." >&2
    exit "${status}"
  fi
  for ((index=${#array_lock_paths[@]}-1; index>=0; index--)); do
    python "${gg_support_dir}/shared_namespace_lock.py" "release-${array_lock_modes[index]}" \
      "${array_lock_paths[index]}" --token "${array_lock_tokens[index]}" || status=1
  done
  exit "${status}"
}
trap input_generation_on_exit EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

prepare_input_generation_tmp_dirs() {
  ensure_dir "${input_generation_tmp_root}"
  ensure_dir "${dir_species_summary_shards}"
  ensure_dir "${dir_task_stats_shards}"
  ensure_dir "${dir_task_meta_shards}"
}

handle_obsolete_input_generation_output() {
  local path=$1
  local description=$2
  if [[ ! -e "${path}" ]]; then
    return 0
  fi
  case "${artifact_stale_policy:-stop}" in
    stop)
      echo "Stale artifact detected: ${description}" >&2
      echo "Reason: output belongs to a species that is absent from the current input set." >&2
      echo "Path: ${path}" >&2
      echo "No artifact files were modified. Use artifact_stale_policy=rebuild to remove obsolete outputs or artifact_stale_policy=reuse to keep them." >&2
      return 3
      ;;
    reuse)
      echo "Warning: keeping obsolete output because artifact_stale_policy=reuse: ${path}" >&2
      return 0
      ;;
    rebuild)
      echo "Removing obsolete output because artifact_stale_policy=rebuild: ${path}"
      rm -f -- "${path}"
      return 0
      ;;
    *)
      echo "Invalid artifact_stale_policy=${artifact_stale_policy:-}; expected stop, reuse, or rebuild." >&2
      return 2
      ;;
  esac
}

remove_formatted_species_outputs_for_rebuild() {
  local species_label=$1
  local managed_dir=""
  local managed_path=""
  for managed_dir in "${species_cds_dir}" "${species_gff_dir}" "${species_genome_dir}"; do
    while IFS= read -r managed_path; do
      [[ -n "${managed_path}" ]] || continue
      case "${managed_path}" in
        "${managed_dir}/"*) rm -f -- "${managed_path}" ;;
        *)
          echo "Refusing to remove an input-generation output outside ${managed_dir}: ${managed_path}" >&2
          return 1
          ;;
      esac
    done < <(gg_find_species_files_by_label "${managed_dir}" "${species_label}")
  done
}

clean_input_generation_shards() {
  if [[ -d "${dir_species_summary_shards}" ]]; then rm -rf -- "${dir_species_summary_shards}"; fi
  if [[ -d "${dir_task_stats_shards}" ]]; then rm -rf -- "${dir_task_stats_shards}"; fi
  if [[ -d "${dir_task_meta_shards}" ]]; then rm -rf -- "${dir_task_meta_shards}"; fi
  ensure_dir "${dir_species_summary_shards}"
  ensure_dir "${dir_task_stats_shards}"
  ensure_dir "${dir_task_meta_shards}"
}

validate_required_formatted_outputs() {
  local summary_file="$1"
  local expected_count="${2:-0}"
  if [[ ${require_cds} -ne 1 && ${require_gff} -ne 1 && ${require_genome} -ne 1 ]]; then
    return 0
  fi
  local cmd=(python "${gg_support_dir}/validate_required_species_outputs.py" --species-summary "${summary_file}")
  [[ ${require_cds} -ne 1 ]] || cmd+=(--require-cds)
  [[ ${require_gff} -ne 1 ]] || cmd+=(--require-gff)
  [[ ${require_genome} -ne 1 ]] || cmd+=(--require-genome)
  [[ ${expected_count} -le 0 ]] || cmd+=(--expected-task-count "${expected_count}")
  "${cmd[@]}"
}

run_format_stage_single() {
  local task="Format species inputs"
  local format_stats_file=""
  local cmd=()
  local cmd_status=0
  local existing_cds=()
  local existing_gff=()
  local existing_genome=()
  local format_needs_update=0
  local format_force_overwrite=${overwrite}
  local format_work_dir=""
  local format_cds_dir="${species_cds_dir}"
  local format_gff_dir="${species_gff_dir}"
  local format_genome_dir="${species_genome_dir}"
  local format_summary_output="${species_summary_output}"
  local format_resolved_output="${resolved_manifest_output}"
  local format_provenance_manifest="${input_generation_provenance_dir}/format.single.json"
  local -a format_provenance_args=()

  if [[ ${run_format_inputs} -eq 1 ]] && ! ensure_selected_download_manifest_has_rows; then
    stage_format_status="failed"
    exit 1
  fi

  if [[ -z "${download_manifest}" && -z "${input_dir}" ]]; then
    while IFS= read -r path; do
      [[ -n "${path}" ]] || continue
      existing_cds+=( "${path}" )
    done < <(list_nonhidden_matching_files "${species_cds_dir}" -name "*.fa" -o -name "*.fa.gz" -o -name "*.fas" -o -name "*.fas.gz" -o -name "*.fasta" -o -name "*.fasta.gz" -o -name "*.fna" -o -name "*.fna.gz")
    while IFS= read -r path; do
      [[ -n "${path}" ]] || continue
      existing_gff+=( "${path}" )
    done < <(list_nonhidden_matching_files "${species_gff_dir}" -name "*.gff" -o -name "*.gff.gz" -o -name "*.gff3" -o -name "*.gff3.gz" -o -name "*.gtf" -o -name "*.gtf.gz")
    while IFS= read -r path; do
      [[ -n "${path}" ]] || continue
      existing_genome+=( "${path}" )
    done < <(list_nonhidden_matching_files "${species_genome_dir}" -name "*.fa" -o -name "*.fa.gz" -o -name "*.fas" -o -name "*.fas.gz" -o -name "*.fasta" -o -name "*.fasta.gz" -o -name "*.fna" -o -name "*.fna.gz")
    if [[ ${#existing_cds[@]} -gt 0 || ${#existing_gff[@]} -gt 0 || ${#existing_genome[@]} -gt 0 ]]; then
      echo "Skipping ${task} because no input source was provided and existing species inputs were detected."
      stage_format_status="skipped"
      run_format_inputs=0
      return 0
    fi
    if [[ ${run_format_inputs} -ne 1 ]]; then
      gg_step_skip "${task}"
      stage_format_status="skipped"
      return 0
    fi
    echo "No input source was specified for formatting."
    echo "Set one of input_dir / download_manifest."
    stage_format_status="failed"
    exit 1
  fi

  if [[ ${download_only} -eq 0 && ${dry_run} -eq 0 ]]; then
    gg_artifact_contract_init format_provenance_args "input_generation_format" "all_species" "${format_provenance_manifest}"
    gg_artifact_add_input_if_present format_provenance_args "input_directory" "${input_dir}"
    gg_artifact_add_input_if_present format_provenance_args "download_manifest" "${download_manifest}"
    format_provenance_args+=(
      --output "species_cds=${species_cds_dir}"
      --output "species_gff=${species_gff_dir}"
      --output "species_genome=${species_genome_dir}"
      --output "species_summary=${species_summary_output}"
      --parameter "provider=${provider}"
      --parameter "gene_grouping_mode=${gene_grouping_mode}"
      --parameter "gff_repair_mode=${gff_repair_mode}"
      --parameter "genetic_code=${GG_COMMON_GENETIC_CODE:-1}"
      --parameter "format_contract_version=${format_contract_version}"
      --parameter "strict=${strict}"
    )
    if [[ -n "${download_manifest}" ]]; then
      format_provenance_args+=(--optional-output "resolved_manifest=${resolved_manifest_output}")
    fi
    if [[ ${overwrite} -eq 1 ]]; then
      format_needs_update=1
    else
      gg_artifact_prepare_stage format_needs_update run_format_inputs "${format_provenance_args[@]}" || return $?
    fi
    if [[ ${format_needs_update} -ne 1 || ${run_format_inputs} -ne 1 ]]; then
      gg_step_skip "${task}"
      stage_format_status="skipped"
      return 0
    fi
    if [[ -s "${format_provenance_manifest}" ]]; then
      format_force_overwrite=1
    fi
  elif [[ ${run_format_inputs} -ne 1 ]]; then
    gg_step_skip "${task}"
    stage_format_status="skipped"
    return 0
  fi

  gg_step_start "${task}"
  stage_format_status="running"

  if [[ ${download_only} -eq 0 && ${dry_run} -eq 0 ]]; then
    ensure_dir "${download_tmp_root}"
    format_work_dir=$(mktemp -d "${download_tmp_root}/formatted-inputs.XXXXXX") || return $?
    format_cds_dir="${format_work_dir}/species_cds"
    format_gff_dir="${format_work_dir}/species_gff"
    format_genome_dir="${format_work_dir}/species_genome"
    format_summary_output="${format_work_dir}/species_summary.tsv"
    format_resolved_output="${format_work_dir}/resolved_manifest.tsv"
    echo "Staging formatted-input outputs before publication."
  fi

  if ! ensure_ete_taxonomy_db "${gg_workspace_dir}"; then
    echo "Warning: Failed to prepare ETE taxonomy DB for species_summary taxonomy metadata. Continuing without taxid/genetic code annotation." >&2
  fi

  format_stats_file="${download_tmp_root}/gg_input_generation_stats.json"
  ensure_parent_dir "${format_stats_file}"
  rm -f -- "${format_stats_file}"
  cmd=(python "${gg_support_dir}/format_species_inputs.py")
  cmd+=(--provider "${provider}")
  cmd+=(--species-cds-dir "${format_cds_dir}")
  cmd+=(--species-gff-dir "${format_gff_dir}")
  cmd+=(--species-genome-dir "${format_genome_dir}")
  cmd+=(--species-summary-output "${format_summary_output}")
  cmd+=(--stats-output "${format_stats_file}")
  cmd+=(--gene-grouping-mode "${gene_grouping_mode}")
  cmd+=(--gff-repair-mode "${gff_repair_mode}")
  cmd+=(--genetic-code "${GG_COMMON_GENETIC_CODE:-1}")

  if [[ -n "${download_manifest}" ]]; then
    cmd+=(--download-manifest "${download_manifest}")
    cmd+=(--download-dir "${download_dir}")
    cmd+=(--resolved-manifest-output "${format_resolved_output}")
  fi
  if [[ -n "${input_dir}" ]]; then
    cmd+=(--input-dir "${input_dir}")
  fi
  if [[ ${overwrite} -eq 1 ]]; then
    cmd+=(--overwrite)
  elif [[ ${format_force_overwrite} -eq 1 ]]; then
    cmd+=(--overwrite-formatted)
  fi
  if [[ ${strict} -eq 1 ]]; then
    cmd+=(--strict)
  fi
  if [[ ${download_only} -eq 1 ]]; then
    cmd+=(--download-only)
  fi
  if [[ ${dry_run} -eq 1 ]]; then
    cmd+=(--dry-run)
  fi
  cmd+=(--download-timeout "${download_timeout}")
  cmd+=(--jobs "${GG_TASK_CPUS:-1}")
  if [[ -n "${auth_bearer_token_env}" ]]; then
    cmd+=(--auth-bearer-token-env "${auth_bearer_token_env}")
  fi
  if [[ -n "${http_header}" ]]; then
    cmd+=(--http-header "${http_header}")
  fi

  echo "Running: ${cmd[*]}"
  if (
    trap 'if [[ -n "${format_work_dir}" ]]; then rm -rf -- "${format_work_dir}"; fi' EXIT
    "${cmd[@]}" || exit $?
    if [[ -n "${format_work_dir}" ]]; then
      python "${gg_support_dir}/relocate_formatted_inputs.py" \
        --mapping "${format_cds_dir}" "${species_cds_dir}" \
        --mapping "${format_gff_dir}" "${species_gff_dir}" \
        --mapping "${format_genome_dir}" "${species_genome_dir}" \
        --summary "${format_summary_output}" || exit $?
      local format_publish_args=(
        "${format_cds_dir}" "${species_cds_dir}"
        "${format_gff_dir}" "${species_gff_dir}"
        "${format_genome_dir}" "${species_genome_dir}"
        "${format_summary_output}" "${species_summary_output}"
      )
      if [[ -f "${format_resolved_output}" ]]; then
        format_publish_args+=("${format_resolved_output}" "${resolved_manifest_output}")
      fi
      mv_out_bundle "${format_publish_args[@]}" || exit $?
    fi
  ); then
    cmd_status=0
  else
    cmd_status=$?
  fi
  if [[ ${cmd_status} -ne 0 ]]; then
    stage_format_status="failed"
    echo "Failed: ${task} (exit=${cmd_status})"
    exit "${cmd_status}"
  fi

  if [[ -s "${format_stats_file}" ]]; then
    read_stats_json_fields "${format_stats_file}" num_species_cds=num_species_cds_files num_species_gff=num_species_gff_files num_species_genome=num_species_genome_files cds_sequences_before cds_sequences_after cds_first_sequence_name || exit $?
    rm -f -- "${format_stats_file}"
  fi
  if [[ ${download_only} -eq 0 && ${dry_run} -eq 0 ]]; then
    gg_artifact_record "${format_provenance_args[@]}"
  fi
  stage_format_status="ok"
}

run_validate_stage() {
  local task="Validate formatted species inputs"
  local mapping_stats_file=""
  local cmd=()
  local cds_files=()
  local gff_files=()
  local genome_files=()
  local set_status=0
  local validation_reuse_args=()

  if [[ ${run_validate_inputs} -ne 1 || ${run_format_inputs} -ne 1 || ${download_only} -ne 0 || ${dry_run} -ne 0 ]]; then
    gg_step_skip "${task}"
    stage_validate_status="skipped"
    return 0
  fi

  gg_step_start "${task}"
  stage_validate_status="running"
  if [[ "${input_generation_mode}" == array_finalize ]]; then
    validation_reuse_args=(--reuse-validation-root "${input_generation_root}"
      --reuse-task-plan "${task_plan_output}" --format-contract-version "${format_contract_version}")
  fi
  while IFS= read -r path; do
    [[ -n "${path}" ]] || continue
    cds_files+=( "${path}" )
  done < <(list_nonhidden_matching_files "${species_cds_dir}" -name "*.fa" -o -name "*.fa.gz" -o -name "*.fas" -o -name "*.fas.gz" -o -name "*.fasta" -o -name "*.fasta.gz" -o -name "*.fna" -o -name "*.fna.gz")
  while IFS= read -r path; do
    [[ -n "${path}" ]] || continue
    gff_files+=( "${path}" )
  done < <(list_nonhidden_matching_files "${species_gff_dir}" -name "*.gff" -o -name "*.gff.gz" -o -name "*.gff3" -o -name "*.gff3.gz" -o -name "*.gtf" -o -name "*.gtf.gz")
  while IFS= read -r path; do
    [[ -n "${path}" ]] || continue
    genome_files+=( "${path}" )
  done < <(list_nonhidden_matching_files "${species_genome_dir}" -name "*.fa" -o -name "*.fa.gz" -o -name "*.fas" -o -name "*.fas.gz" -o -name "*.fasta" -o -name "*.fasta.gz" -o -name "*.fna" -o -name "*.fna.gz")
  num_species_cds=${#cds_files[@]}
  num_species_gff=${#gff_files[@]}
  num_species_genome=${#genome_files[@]}
  echo "Detected ${#cds_files[@]} species CDS files, ${#gff_files[@]} species GFF files, and ${#genome_files[@]} species genome files."

  if [[ ${#cds_files[@]} -eq 0 ]]; then
    echo "No species CDS files were detected in: ${species_cds_dir}"
    if [[ ${strict} -eq 1 ]]; then
      stage_validate_status="failed"
      exit 1
    fi
  else
    check_species_cds_dir "${species_cds_dir}"
  fi

  if [[ ${#gff_files[@]} -eq 0 ]]; then
    echo "No species GFF files were detected in: ${species_gff_dir}"
    if [[ ${strict} -eq 1 ]]; then
      stage_validate_status="failed"
      exit 1
    fi
  fi

  if [[ ${#genome_files[@]} -eq 0 ]]; then
    echo "No species genome files were detected in: ${species_genome_dir}"
    if [[ ${strict} -eq 1 ]]; then
      stage_validate_status="failed"
      exit 1
    fi
  fi

  if [[ ${#cds_files[@]} -gt 0 && ${#gff_files[@]} -gt 0 ]]; then
    if is_species_set_identical "${species_cds_dir}" "${species_gff_dir}"; then
      set_status=0
    else
      set_status=$?
    fi
    if [[ ${set_status} -ne 0 && ${strict} -eq 1 ]]; then
      echo "Species set mismatch between species_cds and species_gff. Exiting due to strict mode."
      stage_validate_status="failed"
      exit 1
    fi
    if [[ ${set_status} -eq 0 ]]; then
      mapping_stats_file="${download_tmp_root}/gg_input_generation_mapping_stats.json"
      longest_cds_stats_file="${download_tmp_root}/gg_input_generation_longest_cds_stats.json"
      ensure_parent_dir "${mapping_stats_file}"
      rm -f -- "${mapping_stats_file}"
      rm -f -- "${longest_cds_stats_file}"
      if [[ "${input_generation_mode}" == array_finalize ]]; then
        cmd=(python "${gg_support_dir}/input_validation_reuse.py")
        cmd+=(--species-cds-dir "${species_cds_dir}" --species-gff-dir "${species_gff_dir}")
        cmd+=(--species-genome-dir "${species_genome_dir}" --species-summary "${species_summary_output}")
        cmd+=(--nthreads "${GG_TASK_CPUS:-1}" --mapping-stats-output "${mapping_stats_file}")
        cmd+=(--ownership-stats-output "${longest_cds_stats_file}" "${validation_reuse_args[@]}")
        [[ ${strict} -ne 1 ]] || cmd+=(--strict)
        echo "Running: ${cmd[*]}"
        if ! "${cmd[@]}"; then
          stage_validate_status="failed"
          echo "Failed: Validate formatted species inputs"
          exit 1
        fi
      else
        cmd=(python "${gg_support_dir}/validate_cds_gff_mapping.py")
        cmd+=(--species-cds-dir "${species_cds_dir}")
        cmd+=(--species-gff-dir "${species_gff_dir}")
        cmd+=(--species-genome-dir "${species_genome_dir}")
        cmd+=(--species-summary "${species_summary_output}")
        cmd+=(--nthreads "${GG_TASK_CPUS:-1}")
        cmd+=(--stats-output "${mapping_stats_file}")
        cmd+=("${validation_reuse_args[@]}")
        if [[ ${strict} -eq 1 ]]; then
          cmd+=(--strict)
        fi
        echo "Running: ${cmd[*]}"
        if "${cmd[@]}"; then
          :
        else
          stage_validate_status="failed"
          echo "Failed: Validate CDS-to-GFF mapping compatibility"
          exit 1
        fi
        if [[ ! -s "${species_summary_output}" ]]; then
          stage_validate_status="failed"
          echo "Species summary not found for longest CDS validation: ${species_summary_output}"
          exit 1
        fi
        cmd=(python "${gg_support_dir}/validate_longest_cds_selection.py")
        cmd+=(--species-cds-dir "${species_cds_dir}")
        cmd+=(--species-summary "${species_summary_output}")
        cmd+=(--nthreads "${GG_TASK_CPUS:-1}")
        cmd+=(--stats-output "${longest_cds_stats_file}")
        cmd+=("${validation_reuse_args[@]}")
        echo "Running: ${cmd[*]}"
        if "${cmd[@]}"; then
          :
        else
          stage_validate_status="failed"
          echo "Failed: Validate longest CDS representative selection"
          exit 1
        fi
      fi
      rm -f -- "${mapping_stats_file}"
      # Retain source-ownership coverage alongside the selection result.
    fi
  fi
  stage_validate_status="ok"
}

run_validate_stage_one_worker() {
  local task="Validate formatted species input"
  local task_summary_file="${dir_species_summary_shards}/${GG_ARRAY_TASK_ID}.tsv"
  local mapping_stats_file="${dir_task_stats_shards}/${GG_ARRAY_TASK_ID}.mapping.json"
  local longest_stats_file="${dir_task_stats_shards}/${GG_ARRAY_TASK_ID}.longest.json"
  local cmd=()

  if [[ ${run_validate_inputs} -ne 1 || ${run_format_inputs} -ne 1 || ${download_only} -ne 0 || ${dry_run} -ne 0 ]]; then
    gg_step_skip "${task}"
    stage_validate_status="skipped"
    return 0
  fi
  if [[ ! -s "${task_summary_file}" ]]; then
    echo "Species summary shard not found for array-worker validation: ${task_summary_file}"
    stage_validate_status="failed"
    exit 1
  fi

  if [[ ${overwrite} -ne 1 ]] && python "${gg_support_dir}/input_generation_stage_resume.py" check \
    --task-plan "${task_plan_output}" --root "${input_generation_root}" \
    --format-contract-version "${format_contract_version}" \
    --task-index "${GG_ARRAY_TASK_ID}" --stage validate
  then
    echo "Reused verified CDS/GFF validation: task ${GG_ARRAY_TASK_ID}"
    stage_validate_status="ok"
    return 0
  fi

  gg_step_start "${task}"
  stage_validate_status="running"
  rm -f -- "${mapping_stats_file}" "${longest_stats_file}"
  # CDS-only tasks are supported when GFF is optional. Validate supplied GFFs
  # and preserve their QC, but do not require annotation solely to run BUSCO.
  if [[ -n "${gff_output_path:-}" ]]; then
    if [[ ! -s "${gff_output_path}" ]]; then
      stage_validate_status="failed"
      echo "Formatted GFF is missing for task ${GG_ARRAY_TASK_ID}: ${gff_output_path}" >&2
      exit 1
    fi
    cmd=(python "${gg_support_dir}/validate_cds_gff_mapping.py")
    cmd+=(--species-cds-dir "${species_cds_dir}")
    cmd+=(--species-gff-dir "${species_gff_dir}")
    cmd+=(--species-genome-dir "${species_genome_dir}")
    cmd+=(--species-summary "${task_summary_file}")
    cmd+=(--nthreads 1)
    cmd+=(--stats-output "${mapping_stats_file}")
    if [[ ${strict} -eq 1 ]]; then
      cmd+=(--strict)
    fi
    echo "Running: ${cmd[*]}"
    if ! "${cmd[@]}"; then
      stage_validate_status="failed"
      echo "Failed: ${task} (CDS-to-GFF mapping)"
      exit 1
    fi
  elif [[ ${require_gff} -eq 1 ]]; then
    stage_validate_status="failed"
    echo "Required formatted GFF is missing for task ${GG_ARRAY_TASK_ID}" >&2
    exit 1
  else
    echo "No GFF supplied for task ${GG_ARRAY_TASK_ID}; validating CDS selection only."
  fi

  cmd=(python "${gg_support_dir}/validate_longest_cds_selection.py")
  cmd+=(--species-cds-dir "${species_cds_dir}")
  cmd+=(--species-summary "${task_summary_file}")
  cmd+=(--summary-outputs-only)
  cmd+=(--nthreads 1)
  cmd+=(--stats-output "${longest_stats_file}")
  echo "Running: ${cmd[*]}"
  if ! "${cmd[@]}"; then
    stage_validate_status="failed"
    echo "Failed: ${task} (longest CDS selection)"
    exit 1
  fi
  # Keep both mapping and source-ownership QC, bound to the validation
  # checkpoint and worker completion receipt.
  python "${gg_support_dir}/input_generation_stage_resume.py" record \
    --task-plan "${task_plan_output}" --root "${input_generation_root}" \
    --format-contract-version "${format_contract_version}" \
    --task-index "${GG_ARRAY_TASK_ID}" --stage validate
  stage_validate_status="ok"
}

run_cds_fx2tab_for_one_file() {
  local seq_full=$1
  local species_name=$2
  local seq_file=""
  local file_sp_cds_fx2tab=""
  local tmp_fx2tab_tsv=""
  local fx2tab_needs_update=0
  local -a fx2tab_provenance_args=()

  seq_file=$(basename "${seq_full}")
  file_sp_cds_fx2tab="${species_cds_fx2tab_dir}/${species_name}_fx2tab_cds.tsv"

  gg_artifact_contract_init fx2tab_provenance_args "input_generation_cds_fx2tab" "${species_name}" "${input_generation_provenance_dir}/fx2tab.${species_name}.json"
  fx2tab_provenance_args+=(
    --input "species_cds=${seq_full}"
    --output "fx2tab=${file_sp_cds_fx2tab}"
    --parameter "length=yes"
    --parameter "name=yes"
    --parameter "gc=yes"
    --parameter "gc_skew=yes"
    --parameter "only_id=yes"
  )
  if [[ ${overwrite} -eq 1 ]]; then
    fx2tab_needs_update=1
  else
    gg_artifact_prepare_stage fx2tab_needs_update run_cds_fx2tab "${fx2tab_provenance_args[@]}" || return $?
  fi
  if [[ ${fx2tab_needs_update} -ne 1 || ${run_cds_fx2tab} -ne 1 ]]; then
    echo "Skipped fx2tab: ${seq_file}"
    return 0
  fi

  ensure_dir "${species_cds_fx2tab_dir}"
  rm -f -- "${file_sp_cds_fx2tab}"
  tmp_fx2tab_tsv=$(mktemp "${download_tmp_root}/fx2tab.${species_name}.XXXXXX.tsv")
  if seqkit fx2tab \
    --threads "${GG_TASK_CPUS:-1}" \
    --length \
    --name \
    --gc \
    --gc-skew \
    --header-line \
    --only-id \
    "${seq_full}" \
    > "${tmp_fx2tab_tsv}"; then
    :
  else
    echo "Failed: seqkit fx2tab for ${seq_full}"
    rm -f -- "${tmp_fx2tab_tsv}"
    exit 1
  fi

  if [[ ! -s "${tmp_fx2tab_tsv}" ]]; then
    echo "seqkit fx2tab produced no output for: ${seq_full}"
    rm -f -- "${tmp_fx2tab_tsv}"
    exit 1
  fi

  mv -- "${tmp_fx2tab_tsv}" "${file_sp_cds_fx2tab}"
  gg_artifact_record "${fx2tab_provenance_args[@]}"
}

run_cds_fx2tab_stage_all() {
  local task="seqkit fx2tab for species CDS files"
  local source_species_input_fasta=()
  local input_species_set=()
  local fx2tab_output_files=()
  local fx2tab_file fx2tab_base fx2tab_species fx2tab_species_found input_species
  local seq_full seq_file sp_ub

  if [[ ${download_only} -ne 0 || ${dry_run} -ne 0 ]]; then
    gg_step_skip "${task}"
    stage_cds_fx2tab_status="skipped"
    return 0
  fi

  gg_step_start "${task}"
  stage_cds_fx2tab_status="running"
  ensure_dir "${species_cds_fx2tab_dir}"
  while IFS= read -r path; do
    [[ -n "${path}" ]] || continue
    source_species_input_fasta+=( "${path}" )
  done < <(gg_find_fasta_files "${species_cds_dir}" 1)
  echo "Number of CDS files for fx2tab: ${#source_species_input_fasta[@]}"
  if [[ ${#source_species_input_fasta[@]} -eq 0 ]]; then
    if [[ ${run_cds_fx2tab} -eq 1 ]]; then
      echo "No CDS file found. Exiting."
      stage_cds_fx2tab_status="failed"
      exit 1
    fi
    gg_step_skip "${task}"
    stage_cds_fx2tab_status="skipped"
    return 0
  fi
  while IFS= read -r species_name; do
    [[ -n "${species_name}" ]] || continue
    input_species_set+=( "${species_name}" )
  done < <(gg_species_names_from_fasta_dir "${species_cds_dir}")
  while IFS= read -r path; do
    [[ -n "${path}" ]] || continue
    fx2tab_output_files+=( "${path}" )
  done < <(find "${species_cds_fx2tab_dir}" -maxdepth 1 -type f -name "*_fx2tab_cds.tsv" 2> /dev/null | sort)
  if [[ ${#fx2tab_output_files[@]} -gt 0 ]]; then
    for fx2tab_file in "${fx2tab_output_files[@]}"; do
      fx2tab_base=$(basename "${fx2tab_file}")
      fx2tab_species=$(gg_species_name_from_path_or_dot "${fx2tab_base}")
      fx2tab_species_found=0
      for input_species in "${input_species_set[@]}"; do
        if [[ "${input_species}" == "${fx2tab_species}" ]]; then
          fx2tab_species_found=1
          break
        fi
      done
      if [[ ${fx2tab_species_found} -eq 0 ]]; then
        handle_obsolete_input_generation_output "${fx2tab_file}" "fx2tab output for removed species" || return $?
      fi
    done
  fi

  for seq_full in "${source_species_input_fasta[@]}"; do
    seq_file=$(basename "${seq_full}")
    sp_ub=$(gg_species_name_from_path_or_dot "${seq_file}")
    gg_step_start "${task}: ${seq_file}"
    run_cds_fx2tab_for_one_file "${seq_full}" "${sp_ub}"
  done
  num_species_cds_fx2tab=$(count_nonhidden_matching_files "${species_cds_fx2tab_dir}" "*_fx2tab_cds.tsv")
  if [[ ${run_cds_fx2tab} -eq 1 ]]; then
    stage_cds_fx2tab_status="ok"
  else
    stage_cds_fx2tab_status="skipped"
  fi
}

run_cds_fx2tab_stage_one_worker() {
  local task="seqkit fx2tab for species CDS files"
  local task_meta_file="${dir_task_meta_shards}/${GG_ARRAY_TASK_ID}.json"
  local species_prefix=""
  local cds_output_path=""

  if [[ ${download_only} -ne 0 || ${dry_run} -ne 0 ]]; then
    gg_step_skip "${task}"
    stage_cds_fx2tab_status="skipped"
    return 0
  fi

  gg_step_start "${task}"
  stage_cds_fx2tab_status="running"
  if [[ ! -s "${task_meta_file}" ]]; then
    echo "Missing task metadata shard for fx2tab worker: ${task_meta_file}"
    stage_cds_fx2tab_status="failed"
    exit 1
  fi
  read_stats_json_fields "${task_meta_file}" species_prefix cds_output_path || exit $?
  if [[ -z "${species_prefix}" || -z "${cds_output_path}" ]]; then
    echo "Task metadata shard is missing fx2tab fields: ${task_meta_file}"
    stage_cds_fx2tab_status="failed"
    exit 1
  fi

  run_cds_fx2tab_for_one_file "${cds_output_path}" "${species_prefix}"
  if [[ ${run_cds_fx2tab} -eq 1 ]]; then
    stage_cds_fx2tab_status="ok"
  else
    stage_cds_fx2tab_status="skipped"
  fi
}

run_species_busco_for_one_file() {
  local seq_full=$1
  local species_name=$2
  local busco_threads=${3:-${GG_TASK_CPUS}}
  local seq_file=""
  local file_sp_busco_full=""
  local file_sp_busco_short=""
  local file_sp_busco_proteins=""
  local dir_busco_db=""
  local dir_busco_lineage=""
  local busco_work_root=""
  local busco_input_fasta=""
  local busco_output_dir=""
  local busco_needs_update=0
  local -a busco_provenance_args=()

  seq_file=$(basename "${seq_full}")
  file_sp_busco_full="${species_busco_full_dir}/${species_name}.busco.full.tsv"
  file_sp_busco_short="${species_busco_short_dir}/${species_name}.busco.short.txt"
  file_sp_busco_proteins="${species_busco_full_dir}/single_copy/${species_name}.json.gz"

  if [[ -z "${busco_lineage_resolved}" ]]; then
    if [[ "${input_generation_mode}" == "array_worker" ]]; then
      ensure_shared_busco_lineage_ready "${task_plan_output}" || return 1
    else
      echo "BUSCO lineage must be resolved before running the serial species stage." >&2
      return 1
    fi
  fi

  gg_artifact_contract_init busco_provenance_args "input_generation_species_busco" "${species_name}" "${input_generation_provenance_dir}/busco.${species_name}.json"
  busco_provenance_args+=(
    --input "species_cds=${seq_full}"
    --output "busco_full=${file_sp_busco_full}"
    --output "busco_short=${file_sp_busco_short}"
    --output "busco_single_copy=${file_sp_busco_proteins}"
    --parameter "busco_lineage_request=${busco_lineage}"
    --parameter "busco_lineage_resolved=${busco_lineage_resolved}"
    --parameter "busco_mode=transcriptome"
    --parameter "evalue=1e-03"
    --parameter "limit=20"
  )
  if [[ ${overwrite} -eq 1 ]]; then
    busco_needs_update=1
  else
    gg_artifact_prepare_stage busco_needs_update run_species_busco "${busco_provenance_args[@]}" || return $?
  fi
  if [[ ${busco_needs_update} -ne 1 || ${run_species_busco} -ne 1 ]]; then
    echo "Skipped BUSCO: ${seq_file}"
    return 0
  fi
  remove_busco_outputs_for_species "${species_busco_full_dir}" "${species_name}" "*busco.full.tsv"
  remove_busco_outputs_for_species "${species_busco_short_dir}" "${species_name}" "*busco.short.txt"
  busco_work_root=$(mktemp -d "${download_tmp_root}/busco.${species_name}.XXXXXX")
  busco_input_fasta="${busco_work_root}/input.fasta"
  busco_output_dir="${busco_work_root}/busco_tmp"
  seqkit seq --threads "${busco_threads}" "${seq_full}" --out-file "${busco_input_fasta}"

  if ! dir_busco_db=$(ensure_busco_download_path "${gg_workspace_dir}" "${busco_lineage_resolved}"); then
    echo "Failed to prepare BUSCO dataset: ${busco_lineage_resolved}"
    rm -rf -- "${busco_work_root}"
    exit 1
  fi
  dir_busco_lineage="${dir_busco_db}/lineages/${busco_lineage_resolved}"

  (
    cd "${busco_work_root}"
    GG_BUSCO_TIMEOUT_SECONDS="${busco_timeout_seconds}" gg_run_busco_with_metaeuk_modified_fas_compat \
      --in "input.fasta" \
      --mode "transcriptome" \
      --out "busco_tmp" \
      --cpu "${busco_threads}" \
      --force \
      --evalue 1e-03 \
      --limit 20 \
      --lineage_dataset "${dir_busco_lineage}" \
      --download_path "${dir_busco_db}" \
      --offline 2>&1 | awk -v max_lines=10000 '
        NR <= max_lines { print; fflush(); next }
        { tail_lines[NR % 50] = $0 }
        END {
          if (NR > max_lines) {
            print "BUSCO output exceeded 10000 lines; intermediate lines suppressed (" NR " total)."
            first = NR > 50 ? NR - 49 : 1
            for (i = first; i <= NR; i++) print tail_lines[i % 50]
          }
        }'
  )

  if copy_busco_tables "${busco_output_dir}" "${busco_lineage_resolved}" "${file_sp_busco_full}" "${file_sp_busco_short}"; then
    python "${gg_support_dir}/busco_guide_tree.py" preserve \
      --run-dir "${busco_output_dir}/run_${busco_lineage_resolved}" \
      --full "${file_sp_busco_full}" --short "${file_sp_busco_short}" \
      --input "${seq_full}" --species "${species_name}" --output "${file_sp_busco_proteins}" || return $?
    rm -rf -- "${busco_work_root}"
    gg_artifact_record "${busco_provenance_args[@]}"
  else
    echo "Failed to locate normalized BUSCO outputs for ${species_name}. Exiting."
    rm -rf -- "${busco_work_root}"
    exit 1
  fi
}

run_species_busco_stage_all() {
  local task="BUSCO analysis of species CDS files"
  local source_species_input_fasta=()
  local input_species_set=()
  local busco_output_files=()
  local busco_file busco_base busco_species busco_species_found input_species
  local seq_full seq_file sp_ub
  local busco_stage_relevant=${run_species_busco}

  gg_step_start "${task}"
  stage_species_busco_status="running"
  while IFS= read -r path; do
    [[ -n "${path}" ]] || continue
    source_species_input_fasta+=( "${path}" )
  done < <(gg_find_fasta_files "${species_cds_dir}" 1)
  echo "Number of CDS files for BUSCO: ${#source_species_input_fasta[@]}"
  if [[ ${#source_species_input_fasta[@]} -eq 0 ]]; then
    if [[ ${run_species_busco} -eq 1 ]]; then
      echo "No CDS file found. Exiting."
      stage_species_busco_status="failed"
      exit 1
    fi
  fi
  while IFS= read -r species_name; do
    [[ -n "${species_name}" ]] || continue
    input_species_set+=( "${species_name}" )
  done < <(gg_species_names_from_fasta_dir "${species_cds_dir}")
  while IFS= read -r path; do
    [[ -n "${path}" ]] || continue
    busco_output_files+=( "${path}" )
  done < <(
    find "${species_busco_full_dir}" "${species_busco_short_dir}" -maxdepth 1 -type f \
      \( -name "*busco.full.tsv" -o -name "*busco.short.txt" \) \
      2> /dev/null | sort
  )
  if [[ ${#busco_output_files[@]} -gt 0 ]] || compgen -G "${input_generation_provenance_dir}/busco.*.json" > /dev/null; then
    busco_stage_relevant=1
  fi
  if [[ ${busco_stage_relevant} -ne 1 ]]; then
    gg_step_skip "${task}"
    stage_species_busco_status="skipped"
    return 0
  fi
  normalize_busco_table_naming "${species_busco_full_dir}" "${species_busco_short_dir}"
  if [[ ${#busco_output_files[@]} -gt 0 ]]; then
    for busco_file in "${busco_output_files[@]}"; do
      busco_base=$(basename "${busco_file}")
      busco_species=$(gg_species_name_from_path_or_dot "${busco_base}")
      busco_species_found=0
      for input_species in "${input_species_set[@]}"; do
        if [[ "${input_species}" == "${busco_species}" ]]; then
          busco_species_found=1
          break
        fi
      done
      if [[ ${busco_species_found} -eq 0 ]]; then
        handle_obsolete_input_generation_output "${busco_file}" "BUSCO output for removed species" || return $?
      fi
    done
  fi
  if [[ ${#source_species_input_fasta[@]} -eq 0 ]]; then
    echo "No current CDS inputs are available for tracked BUSCO outputs." >&2
    stage_species_busco_status="failed"
    return 2
  fi

  if ! ensure_shared_busco_lineage_ready; then
    if ! busco_lineage_resolved=$(gg_resolve_busco_lineage "${gg_workspace_dir}" "${busco_lineage}" "${input_species_set[@]}"); then
      echo "Failed to resolve BUSCO lineage from request: ${busco_lineage}" >&2
      stage_species_busco_status="failed"
      exit 1
    fi
    echo "Resolved BUSCO lineage for species set (${#input_species_set[@]} species): ${busco_lineage_resolved}"
  fi
  local busco_jobs=1
  local busco_memory_job_cap=1
  if [[ "${species_busco_parallel_jobs}" == "auto" ]]; then
    busco_jobs=${GG_TASK_CPUS}
    [[ ${busco_jobs} -gt 4 ]] && busco_jobs=4
    [[ ${busco_jobs} -gt ${#source_species_input_fasta[@]} ]] && busco_jobs=${#source_species_input_fasta[@]}
  elif [[ "${species_busco_parallel_jobs}" =~ ^[0-9]+$ && ${species_busco_parallel_jobs} -ge 1 ]]; then
    busco_jobs=${species_busco_parallel_jobs}
    [[ ${busco_jobs} -gt ${GG_TASK_CPUS} ]] && busco_jobs=${GG_TASK_CPUS}
    [[ ${busco_jobs} -gt ${#source_species_input_fasta[@]} ]] && busco_jobs=${#source_species_input_fasta[@]}
  else
    echo "Invalid species_busco_parallel_jobs=${species_busco_parallel_jobs}; expected auto or a positive integer." >&2
    stage_species_busco_status="failed"
    return 2
  fi
  if ! busco_memory_job_cap=$(gg_memory_parallel_job_cap "${GG_MEM_TOOL_GB}" "${species_busco_memory_gb_per_job}"); then
    stage_species_busco_status="failed"
    return 2
  fi
  if [[ ${busco_jobs} -gt ${busco_memory_job_cap} ]]; then
    echo "Capping BUSCO species parallelism at ${busco_memory_job_cap} job(s) for ${GG_MEM_TOOL_GB}G tool memory (${species_busco_memory_gb_per_job}G/job)."
    busco_jobs=${busco_memory_job_cap}
  fi
  [[ ${busco_jobs} -lt 1 ]] && busco_jobs=1
  local busco_threads_per_job=$((GG_TASK_CPUS / busco_jobs))
  [[ ${busco_threads_per_job} -lt 1 ]] && busco_threads_per_job=1
  echo "BUSCO species parallelism: jobs=${busco_jobs}, threads_per_job=${busco_threads_per_job}"
  for seq_full in "${source_species_input_fasta[@]}"; do
    seq_file=$(basename "${seq_full}")
    sp_ub=$(gg_species_name_from_path_or_dot "${seq_file}")
    gg_step_start "${task}: ${seq_file}"
    if [[ ${busco_jobs} -eq 1 ]]; then
      run_species_busco_for_one_file "${seq_full}" "${sp_ub}" "${busco_threads_per_job}"
    else
      wait_until_jobn_le "${busco_jobs}"
      run_species_busco_for_one_file "${seq_full}" "${sp_ub}" "${busco_threads_per_job}" &
      gg_background_register "$!"
    fi
  done
  if [[ ${busco_jobs} -gt 1 ]]; then
    wait_for_background_jobs
  fi
  num_species_busco_full=$(count_nonhidden_matching_files "${species_busco_full_dir}" "*busco.full.tsv")
  num_species_busco_short=$(count_nonhidden_matching_files "${species_busco_short_dir}" "*busco.short.txt")
  if [[ ${run_species_busco} -eq 1 ]]; then
    stage_species_busco_status="ok"
  else
    stage_species_busco_status="skipped"
  fi
}

run_species_busco_stage_one_worker() {
  local task="BUSCO analysis of species CDS files"
  local task_meta_file="${dir_task_meta_shards}/${GG_ARRAY_TASK_ID}.json"
  local species_prefix=""
  local cds_output_path=""

  gg_step_start "${task}"
  stage_species_busco_status="running"
  if [[ ! -s "${task_meta_file}" ]]; then
    echo "Missing task metadata shard for BUSCO worker: ${task_meta_file}"
    stage_species_busco_status="failed"
    exit 1
  fi
  read_stats_json_fields "${task_meta_file}" species_prefix cds_output_path || exit $?
  if [[ -z "${species_prefix}" || -z "${cds_output_path}" ]]; then
    echo "Task metadata shard is missing BUSCO fields: ${task_meta_file}"
    stage_species_busco_status="failed"
    exit 1
  fi
  normalize_busco_table_naming "${species_busco_full_dir}" "${species_busco_short_dir}"
  run_species_busco_for_one_file "${cds_output_path}" "${species_prefix}"
  if [[ ${run_species_busco} -eq 1 ]]; then
    stage_species_busco_status="ok"
  else
    stage_species_busco_status="skipped"
  fi
}

run_multispecies_summary_stage() {
  local task="Generate multispecies BUSCO summary"
  local cmd=()
  local cmd_status=0
  local summary_needs_update=0
  local -a summary_provenance_args=()

  normalize_busco_table_naming "${species_busco_full_dir}" "${species_busco_short_dir}"
  num_species_busco_full=$(count_nonhidden_matching_files "${species_busco_full_dir}" "*busco.full.tsv")
  num_species_busco_short=$(count_nonhidden_matching_files "${species_busco_short_dir}" "*busco.short.txt")
  if [[ ${run_cds_fx2tab} -eq 1 && -d "${species_cds_fx2tab_dir}" ]]; then
    num_species_cds_fx2tab=$(count_nonhidden_matching_files "${species_cds_fx2tab_dir}" "*_fx2tab_cds.tsv")
  fi

  gg_artifact_contract_init summary_provenance_args "input_generation_multispecies_summary" "all_species" "${input_generation_provenance_dir}/multispecies_summary.json"
  summary_provenance_args+=(
    --input "species_busco_full=${species_busco_full_dir}"
    --input "adapter=${gg_support_dir}/annotation_summary.r"
    --input "busco_plot_metadata=${gg_support_dir}/busco_plot_metadata.r"
    --output "summary=${file_multispecies_summary}"
    --parameter "min_og_species=auto"
    --parameter "include_fx2tab=${run_cds_fx2tab}"
  )
  if [[ ${run_cds_fx2tab} -eq 1 ]]; then
    gg_artifact_add_input_if_present summary_provenance_args "species_cds_fx2tab" "${species_cds_fx2tab_dir}"
  fi
  gg_artifact_add_input_if_present summary_provenance_args "species_trait" "${species_trait_output}"
  if [[ ${overwrite} -eq 1 ]]; then
    summary_needs_update=1
  else
    gg_artifact_prepare_stage summary_needs_update run_multispecies_summary "${summary_provenance_args[@]}" || return $?
  fi
  if [[ ${summary_needs_update} -ne 1 || ${run_multispecies_summary} -ne 1 ]]; then
    gg_step_skip "${task}"
    stage_multispecies_summary_status="skipped"
    return 0
  fi
  if [[ "${num_species_busco_full}" == "0" ]]; then
    echo "No species BUSCO full tables were found. Skipping multispecies summary generation."
    gg_step_skip "${task}"
    stage_multispecies_summary_status="skipped"
    return 0
  fi
  if ! is_species_set_identical "${species_cds_dir}" "${species_busco_full_dir}"; then
    echo "Exiting due to species-set mismatch between ${species_cds_dir} and ${species_busco_full_dir}"
    stage_multispecies_summary_status="failed"
    exit 1
  fi
  if [[ ${run_cds_fx2tab} -eq 1 && "${num_species_cds_fx2tab:-0}" != "0" ]] \
    && ! is_species_set_identical "${species_cds_dir}" "${species_cds_fx2tab_dir}"; then
    echo "Exiting due to species-set mismatch between ${species_cds_dir} and ${species_cds_fx2tab_dir}"
    stage_multispecies_summary_status="failed"
    exit 1
  fi

  gg_step_start "${task}"
  stage_multispecies_summary_status="running"
  ensure_dir "$(dirname "${file_multispecies_summary}")"
  cd "$(dirname "${file_multispecies_summary}")"

  cmd=(Rscript "${gg_support_dir}/annotation_summary.r")
  cmd+=(--dir_species_cds_busco="${species_busco_full_dir}")
  if [[ ${run_cds_fx2tab} -eq 1 ]]; then
    cmd+=(--dir_species_cds_fx2tab="${species_cds_fx2tab_dir}")
  fi
  if [[ -s "${species_trait_output}" ]]; then
    cmd+=(--file_species_trait="${species_trait_output}")
  fi
  cmd+=(--min_og_species=auto)
  echo "Running: ${cmd[*]}"
  if "${cmd[@]}"; then
    cmd_status=0
  else
    cmd_status=$?
  fi
  if [[ ${cmd_status} -ne 0 ]]; then
    stage_multispecies_summary_status="failed"
    echo "Failed: ${task} (exit=${cmd_status})"
    exit "${cmd_status}"
  fi
  if [[ -e "Rplots.pdf" ]]; then
    rm -f -- "Rplots.pdf"
  fi
  gg_artifact_record "${summary_provenance_args[@]}"
  cd "${gg_workspace_dir}"
  stage_multispecies_summary_status="ok"
}

run_trait_stage() {
  local task="Generate species trait table"
  local trait_manifest_path=""
  local trait_manifest_default=""
  local trait_stats_file=""
  local cmd=()
  local cmd_status=0
  local trait_needs_update=0
  local -a trait_provenance_args=()

  gg_step_start "${task}"
  stage_trait_status="running"

  trait_manifest_path="${download_manifest}"
  if [[ "${input_generation_mode}" == array_finalize ]]; then
    trait_manifest_path="${task_plan_output}.manifest.tsv"
    python "${gg_support_dir}/input_generation_array_state.py" export-manifest --task-plan "${task_plan_output}" --outfile "${trait_manifest_path}"
  elif [[ -z "${trait_manifest_path}" ]]; then
    trait_manifest_default="${gg_workspace_input_dir}/input_generation/download_plan.xlsx"
    if [[ -s "${trait_manifest_default}" ]]; then
      trait_manifest_path="${trait_manifest_default}"
    fi
  fi

  local -a gbif_cli_args=()
  [[ -z "${gbif_api}" ]] || gbif_cli_args+=(--gbif-api "${gbif_api}")
  [[ -z "${gbif_page_size}" ]] || gbif_cli_args+=(--gbif-page-size "${gbif_page_size}")
  [[ -z "${gbif_max_occurrences_per_species}" ]] || gbif_cli_args+=(--gbif-max-occurrences-per-species "${gbif_max_occurrences_per_species}")
  [[ -z "${gbif_grid_degrees}" ]] || gbif_cli_args+=(--gbif-grid-degrees "${gbif_grid_degrees}")
  [[ -z "${gbif_min_match_confidence}" ]] || gbif_cli_args+=(--gbif-min-match-confidence "${gbif_min_match_confidence}")
  [[ -z "${gbif_max_coordinate_uncertainty_m}" ]] || gbif_cli_args+=(--gbif-max-coordinate-uncertainty-m "${gbif_max_coordinate_uncertainty_m}")
  [[ -z "${gbif_min_distance_from_known_centroid_m}" ]] || gbif_cli_args+=(--gbif-min-distance-from-known-centroid-m "${gbif_min_distance_from_known_centroid_m}")
  [[ -z "${gbif_year_min}" ]] || gbif_cli_args+=(--gbif-year-min "${gbif_year_min}")
  [[ -z "${gbif_year_max}" ]] || gbif_cli_args+=(--gbif-year-max "${gbif_year_max}")
  [[ -z "${gbif_countries}" ]] || gbif_cli_args+=(--gbif-countries "${gbif_countries}")
  [[ -z "${gbif_include_basis_of_record}" ]] || gbif_cli_args+=(--gbif-include-basis-of-record "${gbif_include_basis_of_record}")
  [[ -z "${gbif_exclude_basis_of_record}" ]] || gbif_cli_args+=(--gbif-exclude-basis-of-record "${gbif_exclude_basis_of_record}")
  [[ -z "${gbif_include_establishment_means}" ]] || gbif_cli_args+=(--gbif-include-establishment-means "${gbif_include_establishment_means}")
  [[ -z "${gbif_missing_date}" ]] || gbif_cli_args+=(--gbif-missing-date "${gbif_missing_date}")
  [[ -z "${gbif_missing_uncertainty}" ]] || gbif_cli_args+=(--gbif-missing-uncertainty "${gbif_missing_uncertainty}")
  [[ -z "${gbif_missing_centroid_distance}" ]] || gbif_cli_args+=(--gbif-missing-centroid-distance "${gbif_missing_centroid_distance}")
  [[ -z "${gbif_use_cache}" ]] || gbif_cli_args+=(--gbif-use-cache "${gbif_use_cache}")
  [[ -z "${gbif_require_complete}" ]] || gbif_cli_args+=(--gbif-require-complete "${gbif_require_complete}")
  [[ -z "${gbif_occurrence_file}" ]] || gbif_cli_args+=(--gbif-occurrence-file "${gbif_occurrence_file}")
  [[ -z "${gbif_taxon_map}" ]] || gbif_cli_args+=(--gbif-taxon-map "${gbif_taxon_map}")
  [[ -z "${gbif_download_metadata}" ]] || gbif_cli_args+=(--gbif-download-metadata "${gbif_download_metadata}")
  local gbif_source_identity
  local -a gbif_identity_cmd=(python "${gg_support_dir}/generate_species_trait.py" --database-sources "${trait_database_sources}" --trait-plan "${trait_plan}" --databases "${trait_databases}" --output "${species_trait_output}")
  if (( ${#gbif_cli_args[@]} > 0 )); then
    gbif_identity_cmd+=("${gbif_cli_args[@]}")
  fi
  gbif_identity_cmd+=(--print-gbif-input-identity)
  gbif_source_identity=$("${gbif_identity_cmd[@]}") || return $?
  gg_artifact_contract_init trait_provenance_args "input_generation_species_trait" "all_species" "${input_generation_provenance_dir}/species_trait.json"
  gg_artifact_add_input_if_present trait_provenance_args "download_manifest" "${trait_manifest_path}"
  if [[ "${trait_species_source}" == "species_cds" ]]; then
    gg_artifact_add_input_if_present trait_provenance_args "species_cds" "${species_cds_dir}"
  fi
  gg_artifact_add_input_if_present trait_provenance_args "trait_plan" "${trait_plan}"
  gg_artifact_add_input_if_present trait_provenance_args "database_sources" "${trait_database_sources}"
  trait_provenance_args+=(
    --input "trait_adapter=${gg_support_dir}/generate_species_trait.py"
    --input "gift_retrieval=${gg_support_dir}/gift_retrieval.py"
    --input "public_plant_traits=${gg_support_dir}/public_plant_traits.py"
    --input "trait_schema_adapter=${gg_support_dir}/species_trait_schema.py"
    --input "gift_reviewed_mappings=${gg_support_dir}/gift_species_mappings.tsv"
  )
  local gift_mapping_inputs=""
  local gift_mapping_path=""
  local gift_mapping_index=0
  gift_mapping_inputs=$(python "${gg_support_dir}/generate_species_trait.py" --print-gift-mapping-inputs --database-sources "${trait_database_sources}") || return $?
  while IFS= read -r gift_mapping_path; do
    [[ -n "${gift_mapping_path}" ]] || continue
    trait_provenance_args+=(--input "gift_custom_mapping_${gift_mapping_index}=${gift_mapping_path}")
    gift_mapping_index=$((gift_mapping_index + 1))
  done <<< "${gift_mapping_inputs}"
  trait_provenance_args+=(
    --output "species_trait=${species_trait_output}"
    --output "species_trait_schema=${species_trait_output}.schema.json"
    --output "trait_metadata=${species_trait_output}.metadata.json"
    --output "gbif_quality=${species_trait_output}.gbif-quality.tsv"
    --output "gbif_observations=${species_trait_output}.gbif-observations.tsv"
    --input "trait_generator=${gg_support_dir}/generate_species_trait.py"
    --input "gbif_adapter=${gg_support_dir}/gbif_observations.py"
    --input "trait_contract=${gg_support_dir}/species_trait_contract.py"
    --parameter "gbif_source_identity=${gbif_source_identity}"
    --parameter "trait_profile=${trait_profile}"
    --parameter "trait_species_source=${trait_species_source}"
    --parameter "trait_databases=${trait_databases}"
    --parameter "gbif_api=${gbif_api}"
    --parameter "gbif_page_size=${gbif_page_size}"
    --parameter "gbif_max_occurrences_per_species=${gbif_max_occurrences_per_species}"
    --parameter "gbif_grid_degrees=${gbif_grid_degrees}"
    --parameter "gbif_min_match_confidence=${gbif_min_match_confidence}"
    --parameter "gbif_max_coordinate_uncertainty_m=${gbif_max_coordinate_uncertainty_m}"
    --parameter "gbif_min_distance_from_known_centroid_m=${gbif_min_distance_from_known_centroid_m}"
    --parameter "gbif_year_min=${gbif_year_min}"
    --parameter "gbif_year_max=${gbif_year_max}"
    --parameter "gbif_countries=${gbif_countries}"
    --parameter "gbif_include_basis_of_record=${gbif_include_basis_of_record}"
    --parameter "gbif_exclude_basis_of_record=${gbif_exclude_basis_of_record}"
    --parameter "gbif_include_establishment_means=${gbif_include_establishment_means}"
    --parameter "gbif_missing_date=${gbif_missing_date}"
    --parameter "gbif_missing_uncertainty=${gbif_missing_uncertainty}"
    --parameter "gbif_missing_centroid_distance=${gbif_missing_centroid_distance}"
    --parameter "gbif_use_cache=${gbif_use_cache}"
    --parameter "gbif_require_complete=${gbif_require_complete}"
    --parameter "gbif_occurrence_file=${gbif_occurrence_file}"
    --parameter "gbif_taxon_map=${gbif_taxon_map}"
    --parameter "gbif_download_metadata=${gbif_download_metadata}"
    --parameter "strict=${strict}"
  )
  if [[ ${dry_run} -eq 0 ]]; then
    if [[ ${overwrite} -eq 1 ]]; then
      trait_needs_update=1
    else
      gg_artifact_prepare_stage trait_needs_update run_generate_species_trait "${trait_provenance_args[@]}" || return $?
    fi
    if [[ ${trait_needs_update} -ne 1 || ${run_generate_species_trait} -ne 1 ]]; then
      gg_step_skip "${task}"
      stage_trait_status="skipped"
      return 0
    fi
  fi

  trait_stats_file="${download_tmp_root}/gg_input_generation_trait_stats.json"
  ensure_parent_dir "${trait_stats_file}"
  rm -f -- "${trait_stats_file}"
  cmd=(python "${gg_support_dir}/generate_species_trait.py")
  cmd+=(--species-source "${trait_species_source}")
  cmd+=(--species-cds-dir "${species_cds_dir}")
  cmd+=(--trait-plan "${trait_plan}")
  cmd+=(--database-sources "${trait_database_sources}")
  cmd+=(--databases "${trait_databases}")
  cmd+=(--downloads-dir "${trait_download_dir}")
  cmd+=(--output "${species_trait_output}")
  cmd+=(--download-timeout "${trait_download_timeout}")
  cmd+=(--stats-output "${trait_stats_file}")
  if (( ${#gbif_cli_args[@]} > 0 )); then
    cmd+=("${gbif_cli_args[@]}")
  fi
  if [[ -n "${trait_manifest_path}" ]]; then
    cmd+=(--download-manifest "${trait_manifest_path}")
  fi
  if [[ ${strict} -eq 1 ]]; then
    cmd+=(--strict)
  fi
  if [[ ${dry_run} -eq 1 ]]; then
    cmd+=(--dry-run)
  fi

  echo "Running: ${cmd[*]}"
  if "${cmd[@]}"; then
    cmd_status=0
  else
    cmd_status=$?
  fi
  if [[ ${cmd_status} -ne 0 ]]; then
    stage_trait_status="failed"
    echo "Failed: ${task} (exit=${cmd_status})"
    exit "${cmd_status}"
  fi

  if [[ -s "${trait_stats_file}" ]]; then
    read_stats_json_fields "${trait_stats_file}" num_species_trait=num_species_with_any_trait num_trait_columns || exit $?
    rm -f -- "${trait_stats_file}"
  fi
  if [[ ${dry_run} -eq 0 ]]; then
    gg_artifact_record "${trait_provenance_args[@]}"
  fi
  stage_trait_status="ok"
}

run_array_prepare_mode() {
  local task="Prepare input-generation array task plan"
  local effective_input_dir=""
  local cmd=()
  local cmd_status=0
  local expected_tasks=0
  local phase_started
  local taxonomy_dataset_status=ok

  write_run_summary_on_exit=0
  rm -f -- "${task_plan_output}.prepared.json"
  prepare_input_generation_tmp_dirs

  gg_step_start "${task}"
  stage_format_status="running"

  if ! ensure_selected_download_manifest_has_rows; then
    stage_format_status="failed"
    exit 1
  fi


  effective_input_dir=$(input_generation_effective_input_dir_path)
  if [[ -z "${effective_input_dir}" ]]; then
    echo "Failed to resolve effective input_dir for array task planning."
    stage_format_status="failed"
    exit 1
  fi

  cmd=(python "${gg_support_dir}/plan_input_generation_tasks.py")
  cmd+=(--provider "${provider}")
  cmd+=(--input-dir "${effective_input_dir}")
  cmd+=(--outfile "${task_plan_output}")
  if [[ -n "${download_manifest}" ]]; then
    cmd+=(--download-manifest "${download_manifest}" --download-dir "${download_dir}" --stage-downloads)
  fi
  cmd+=(--gene-grouping-mode "${gene_grouping_mode}")
  cmd+=(--gff-repair-mode "${gff_repair_mode}")
  cmd+=(--genetic-code "${GG_COMMON_GENETIC_CODE:-1}")
  if [[ ${strict} -eq 1 ]]; then
    cmd+=(--strict)
  fi
  [[ ${require_gff} -ne 1 ]] || cmd+=(--require-gff)
  [[ ${require_genome} -ne 1 ]] || cmd+=(--require-genome)
  echo "Running: ${cmd[*]}"
  if "${cmd[@]}"; then
    cmd_status=0
  else
    cmd_status=$?
  fi
  if [[ ${cmd_status} -ne 0 ]]; then
    stage_format_status="failed"
    echo "Failed: task-plan generation (exit=${cmd_status})"
    exit "${cmd_status}"
  fi

  python "${gg_support_dir}/input_generation_array_state.py" claim-workspace --task-plan "${task_plan_output}" --workspace "${input_generation_root}" "${array_output_args[@]}" --prepare
  if [[ -n "${download_manifest}" ]]; then
    cmd=(python "${gg_support_dir}/stage_input_generation_downloads.py" --task-plan "${task_plan_output}"
      --jobs "${GG_TASK_CPUS:-1}" --download-timeout "${download_timeout}")
    [[ ${require_gff} -ne 1 ]] || cmd+=(--require-gff)
    [[ ${require_genome} -ne 1 ]] || cmd+=(--require-genome)
    [[ -z "${auth_bearer_token_env}" ]] || cmd+=(--auth-bearer-token-env "${auth_bearer_token_env}")
    [[ -z "${http_header}" ]] || cmd+=(--http-header "${http_header}")
    "${cmd[@]}"
  fi
  expected_tasks=$(task_plan_task_count "${task_plan_output}")
  num_species_cds="${expected_tasks}"
  num_species_gff=""
  num_species_genome=""

  if [[ ${run_species_busco} -eq 1 ]]; then
    phase_started=$(python "${gg_support_dir}/performance_metrics.py" start)
    ensure_parent_dir "${file_busco_lineage_resolved}"
    resolve_busco_lineage_from_task_plan "${task_plan_output}"
    ensure_busco_download_path "${gg_workspace_dir}" "${busco_lineage_resolved}" >/dev/null || exit 1
    python "${gg_support_dir}/performance_metrics.py" elapsed --phase busco_dataset_prepare --started "${phase_started}"
  fi
  phase_started=$(python "${gg_support_dir}/performance_metrics.py" start)
  if ! ensure_ete_taxonomy_db "${gg_workspace_dir}"; then
    taxonomy_dataset_status=failed
    echo "Warning: Failed to prepare ETE taxonomy DB before array workers." >&2
  fi
  python "${gg_support_dir}/performance_metrics.py" elapsed --phase taxonomy_dataset_prepare --started "${phase_started}" --status "${taxonomy_dataset_status}"
  local resume_prefix donor_plan donor_sha donor_root
  for resume_prefix in resume_from resume_fallback; do
    donor_plan=${resume_prefix}_task_plan
    donor_sha=${resume_prefix}_task_plan_sha256
    donor_root=${resume_prefix}_input_generation_root
    if [[ -n "${!donor_plan}${!donor_sha}${!donor_root}" ]]; then
      [[ -n "${resume_from_task_plan}" && -n "${!donor_plan}" && -n "${!donor_sha}" && -n "${!donor_root}" ]] || {
        echo "Stage resume requires each donor plan, SHA-256 and output root." >&2
        exit 1
      }
      python "${gg_support_dir}/input_generation_stage_resume.py" check-source \
        --task-plan "${task_plan_output}" --root "${input_generation_root}" \
        --format-contract-version "${format_contract_version}" \
        --source-plan "${!donor_plan}" --source-plan-sha256 "${!donor_sha}" --source-root "${!donor_root}"
    fi
  done
  local prepared_cmd=(python "${gg_support_dir}/input_generation_array_state.py" prepared --task-plan "${task_plan_output}")
  if [[ -n "${download_manifest}" ]]; then
    local staged_file
    for staged_file in "${task_plan_output}.tasks/"*.json "${task_plan_output}.tasks/"*.resolved.tsv; do
      prepared_cmd+=(--file "${staged_file}")
    done
  fi
  [[ ${run_species_busco} -ne 1 ]] || prepared_cmd+=(--file "${file_busco_lineage_resolved}")
  "${prepared_cmd[@]}"
  stage_format_status="ok"
}

run_array_worker_mode() {
  local task="Format a single input-generation species task"
  local task_stats_file=""
  local task_meta_file=""
  local task_summary_file=""
  local cmd=()
  local cmd_status=0
  local describe_cmd=()
  local array_task_lock_token
  local format_needs_update=0
  local format_force_overwrite=${overwrite}
  local species_prefix=""
  local cds_input_path=""
  local gff_input_path=""
  local gbff_input_path=""
  local genome_input_path=""
  local cds_output_path=""
  local gff_output_path=""
  local genome_output_path=""
  local format_provenance_manifest=""
  local -a format_provenance_args=()

  write_run_summary_on_exit=0
  prepare_input_generation_tmp_dirs

  if [[ ! -s "${task_plan_output}" ]]; then
    echo "Task plan not found for array worker: ${task_plan_output}"
    exit 1
  fi
  if [[ ! "${GG_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]]; then
    echo "Invalid GG_ARRAY_TASK_ID value (must be a positive integer): ${GG_ARRAY_TASK_ID}"
    exit 1
  fi

  GG_ARRAY_TASK_ID=$(python "${gg_support_dir}/input_generation_array_state.py" index --task-plan "${task_plan_output}" --task-index "${GG_ARRAY_TASK_ID}")

  # Namespace ownership prevents duplicate/requeued tasks from sharing outputs.
  ensure_dir "${task_plan_output}.locks"
  input_generation_lock "${task_plan_output}.locks/${GG_ARRAY_TASK_ID}.lock" exclusive
  array_task_lock_token=${array_lock_tokens[${#array_lock_tokens[@]}-1]}
  python "${gg_support_dir}/input_generation_array_state.py" invalidate --task-plan "${task_plan_output}" --task-index "${GG_ARRAY_TASK_ID}"
  gg_step_start "${task}"
  stage_format_status="running"
  task_stats_file="${dir_task_stats_shards}/${GG_ARRAY_TASK_ID}.json"
  task_meta_file="${dir_task_meta_shards}/${GG_ARRAY_TASK_ID}.json"
  task_summary_file="${dir_species_summary_shards}/${GG_ARRAY_TASK_ID}.tsv"

  describe_cmd=(python "${gg_support_dir}/input_generation_stage_resume.py" start-worker)
  describe_cmd+=(--task-plan "${task_plan_output}" --task-index "${GG_ARRAY_TASK_ID}")
  describe_cmd+=(--root "${input_generation_root}" --format-contract-version "${format_contract_version}")
  describe_cmd+=(--target-lock-token "${array_task_lock_token}" --download-timeout "${download_timeout}")
  [[ ${overwrite} -ne 1 ]] || describe_cmd+=(--overwrite)
  [[ -z "${http_header}" ]] || describe_cmd+=(--http-header "${http_header}")
  [[ -z "${auth_bearer_token_env}" ]] || describe_cmd+=(--auth-bearer-token-env "${auth_bearer_token_env}")
  if [[ -n "${resume_from_task_plan}" && ${overwrite} -ne 1 ]]; then
    describe_cmd+=(--source-plan "${resume_from_task_plan}" --source-plan-sha256 "${resume_from_task_plan_sha256}"
      --source-root "${resume_from_input_generation_root}")
    if [[ -n "${resume_fallback_task_plan}" ]]; then
      describe_cmd+=(--fallback-source-plan "${resume_fallback_task_plan}"
        --fallback-source-plan-sha256 "${resume_fallback_task_plan_sha256}"
        --fallback-source-root "${resume_fallback_input_generation_root}")
    fi
  fi
  if ! "${describe_cmd[@]}"; then
    stage_format_status="failed"
    echo "Failed to describe/import input-generation array task ${GG_ARRAY_TASK_ID}."
    exit 1
  fi
  read_stats_json_fields "${task_meta_file}" species_prefix cds_input_path=cds_path \
    gff_input_path=gff_path gbff_input_path=gbff_path genome_input_path=genome_path \
    cds_output_path gff_output_path genome_output_path || exit $?
  if [[ -z "${species_prefix}" || -z "${cds_output_path}" || ( ${require_gff} -eq 1 && -z "${gff_output_path}" ) ]]; then
    stage_format_status="failed"
    echo "Array task description is missing required species or output paths: ${task_meta_file}"
    exit 1
  fi
  if [[ ${require_genome} -eq 1 && ( -z "${genome_input_path}" || ! -s "${genome_input_path}" ) && ( -z "${gbff_input_path}" || ! -s "${gbff_input_path}" ) ]]; then
    stage_format_status="failed"
    echo "Required genome input is missing for ${species_prefix}" >&2
    exit 1
  fi
  if [[ ${require_gff} -eq 1 && ( -z "${gff_input_path}" || ! -s "${gff_input_path}" ) && ( -z "${gbff_input_path}" || ! -s "${gbff_input_path}" ) ]]; then
    stage_format_status="failed"
    echo "Required GFF input is missing for ${species_prefix}" >&2
    exit 1
  fi
  if [[ ${overwrite} -ne 1 ]] && python "${gg_support_dir}/input_generation_stage_resume.py" check \
    --task-plan "${task_plan_output}" --root "${input_generation_root}" \
    --format-contract-version "${format_contract_version}" \
    --task-index "${GG_ARRAY_TASK_ID}" --stage format
  then
    echo "Reused verified input formatting: ${species_prefix}"
    stage_format_status="ok"
  else
    rm -f -- "${task_stats_file}" "${task_summary_file}"
    format_provenance_manifest="${input_generation_provenance_dir}/format.${species_prefix}.json"
    gg_artifact_contract_init format_provenance_args "input_generation_format" "${species_prefix}" "${format_provenance_manifest}"
    gg_artifact_add_input_if_present format_provenance_args "cds_input" "${cds_input_path}"
    gg_artifact_add_input_if_present format_provenance_args "gff_input" "${gff_input_path}"
    gg_artifact_add_input_if_present format_provenance_args "gbff_input" "${gbff_input_path}"
    gg_artifact_add_input_if_present format_provenance_args "genome_input" "${genome_input_path}"
    format_provenance_args+=(--output "formatted_cds=${cds_output_path}")
    if [[ -n "${gff_output_path}" ]]; then
      if [[ ${require_gff} -eq 1 || -n "${gff_input_path}" ]]; then
        format_provenance_args+=(--output "formatted_gff=${gff_output_path}")
      else
        format_provenance_args+=(--optional-output "formatted_gff=${gff_output_path}")
      fi
    fi
    format_provenance_args+=(
      --parameter "provider=${provider}"
      --parameter "gene_grouping_mode=${gene_grouping_mode}"
      --parameter "gff_repair_mode=${gff_repair_mode}"
      --parameter "genetic_code=${GG_COMMON_GENETIC_CODE:-1}"
      --parameter "format_contract_version=${format_contract_version}"
      --parameter "strict=${strict}"
    )
    if [[ -n "${genome_output_path}" ]]; then
      format_provenance_args+=(--optional-output "formatted_genome=${genome_output_path}")
    fi
    if [[ ${overwrite} -eq 1 ]]; then
      format_needs_update=1
    else
      gg_artifact_prepare_stage format_needs_update run_format_inputs "${format_provenance_args[@]}" || exit $?
    fi
    if [[ ${format_needs_update} -eq 1 && -s "${format_provenance_manifest}" ]]; then
      format_force_overwrite=1
      remove_formatted_species_outputs_for_rebuild "${species_prefix}"
    fi

    if ! ensure_ete_taxonomy_db "${gg_workspace_dir}"; then
      echo "Warning: Failed to prepare ETE taxonomy DB for species_summary taxonomy metadata. Continuing without taxid/genetic code annotation." >&2
    fi

    cmd=(python "${gg_support_dir}/run_input_generation_task.py")
    cmd+=(--task-plan "${task_plan_output}")
    cmd+=(--task-index "${GG_ARRAY_TASK_ID}")
    cmd+=(--species-cds-dir "${species_cds_dir}")
    cmd+=(--species-gff-dir "${species_gff_dir}")
    cmd+=(--species-genome-dir "${species_genome_dir}")
    cmd+=(--species-summary-output "${task_summary_file}")
    cmd+=(--stats-output "${task_stats_file}")
    cmd+=(--task-meta-output "${task_meta_file}")
    if [[ ${format_force_overwrite} -eq 1 ]]; then
      cmd+=(--overwrite)
    fi
    if [[ "${artifact_stale_policy:-stop}" == "reuse" ]]; then
      cmd+=(--reuse-existing)
    fi
    echo "Running: ${cmd[*]}"
    if "${cmd[@]}"; then
      cmd_status=0
    else
      cmd_status=$?
    fi
    if [[ ${cmd_status} -ne 0 ]]; then
      stage_format_status="failed"
      echo "Failed: ${task} (exit=${cmd_status})"
      exit "${cmd_status}"
    fi
    if [[ ${require_genome} -eq 1 && ( -z "${genome_output_path}" || ! -s "${genome_output_path}" ) ]]; then
      stage_format_status="failed"
      echo "Required formatted genome is missing for ${species_prefix}" >&2
      exit 1
    fi
    validate_required_formatted_outputs "${task_summary_file}" 1 || {
      stage_format_status="failed"
      exit 1
    }
    if [[ ${format_needs_update} -eq 1 ]]; then
      gg_artifact_record "${format_provenance_args[@]}"
    fi
    python "${gg_support_dir}/input_generation_stage_resume.py" record \
      --task-plan "${task_plan_output}" --root "${input_generation_root}" \
    --format-contract-version "${format_contract_version}" \
      --task-index "${GG_ARRAY_TASK_ID}" --stage format
    stage_format_status="ok"
  fi

  if [[ -s "${task_stats_file}" ]]; then
    read_stats_json_fields "${task_stats_file}" num_species_cds=num_species_cds_files num_species_gff=num_species_gff_files num_species_genome=num_species_genome_files cds_sequences_before cds_sequences_after cds_first_sequence_name || exit $?
  fi

  run_validate_stage_one_worker
  run_cds_fx2tab_stage_one_worker
  run_species_busco_stage_one_worker
  stage_multispecies_summary_status="skipped"
  stage_trait_status="skipped"
  local receipt_cmd=(python "${gg_support_dir}/input_generation_array_state.py" complete
    --task-plan "${task_plan_output}" --task-index "${GG_ARRAY_TASK_ID}"
    --file "${task_plan_output}.settings.json" --file "${task_stats_file}" --file "${task_summary_file}" --file "${task_meta_file}"
    --file "${cds_output_path}")
  [[ -z "${gff_output_path}" || ! -s "${gff_output_path}" ]] || receipt_cmd+=(--file "${gff_output_path}")
  if [[ ${stage_validate_status} == "ok" ]]; then
    [[ -s "${dir_task_stats_shards}/${GG_ARRAY_TASK_ID}.longest.json" ]] || {
      echo "CDS selection/source ownership QC is missing for task ${GG_ARRAY_TASK_ID}" >&2
      exit 1
    }
    receipt_cmd+=(--file "${dir_task_stats_shards}/${GG_ARRAY_TASK_ID}.longest.json")
  fi
  if [[ ${run_validate_inputs} -eq 1 && -n "${gff_output_path}" ]]; then
    [[ -s "${dir_task_stats_shards}/${GG_ARRAY_TASK_ID}.mapping.json" ]] || {
      echo "CDS/GFF mapping QC is missing for task ${GG_ARRAY_TASK_ID}" >&2
      exit 1
    }
    receipt_cmd+=(--file "${dir_task_stats_shards}/${GG_ARRAY_TASK_ID}.mapping.json")
  fi
  local raw_input_path
  for raw_input_path in "${cds_input_path}" "${gff_input_path}" "${gbff_input_path}" "${genome_input_path}"; do
    [[ -z "${raw_input_path}" ]] || receipt_cmd+=(--file "${raw_input_path}")
  done
  if [[ -s "${task_plan_output}.tasks/${GG_ARRAY_TASK_ID}.resolved.tsv" ]]; then
    receipt_cmd+=(--file "${task_plan_output}.tasks/${GG_ARRAY_TASK_ID}.resolved.tsv" --file "${task_plan_output}.tasks/${GG_ARRAY_TASK_ID}.json")
  fi
  [[ -z "${genome_output_path}" ]] || receipt_cmd+=(--file "${genome_output_path}")
  if [[ ${run_cds_fx2tab} -eq 1 ]]; then
    receipt_cmd+=(--file "${species_cds_fx2tab_dir}/${species_prefix}_fx2tab_cds.tsv")
  fi
  if [[ ${run_species_busco} -eq 1 ]]; then
    receipt_cmd+=(--file "${species_busco_full_dir}/single_copy/${species_prefix}.json.gz"
      --file "${species_busco_full_dir}/${species_prefix}.busco.full.tsv"
      --file "${species_busco_short_dir}/${species_prefix}.busco.short.txt")
  fi
  "${receipt_cmd[@]}"
}

run_array_finalize_mode() {
  local canonical_summary_output="${species_summary_output}"
  local species_summary_output=""
  local staged_resolved_manifest=""
  local mapping_qc_output="${input_generation_root}/species_mapping_qc.tsv"
  local staged_mapping_qc=""
  ensure_parent_dir "${canonical_summary_output}"
  ensure_parent_dir "${resolved_manifest_output}"
  species_summary_output=$(mktemp "${canonical_summary_output}.array.XXXXXX")
  staged_resolved_manifest=$(mktemp "${resolved_manifest_output}.array.XXXXXX")
  staged_mapping_qc=$(mktemp "${mapping_qc_output}.array.XXXXXX")
  local task="Finalize input-generation array outputs"
  local merge_stats_file="${input_generation_tmp_root}/merged_task_stats.json"
  local expected_tasks=0
  local merged_rows=0
  local task_stats_files=0
  local cmd=()
  local cmd_status=0

  prepare_input_generation_tmp_dirs
  if [[ ! -s "${task_plan_output}" ]]; then
    echo "Task plan not found for array finalize: ${task_plan_output}"
    exit 1
  fi
  expected_tasks=$(task_plan_task_count "${task_plan_output}")
  if [[ ${expected_tasks} -le 0 ]]; then
    echo "Task plan does not contain any tasks: ${task_plan_output}"
    exit 1
  fi

  gg_step_start "${task}"
  stage_format_status="running"
  cmd=(python "${gg_support_dir}/merge_input_generation_shards.py")
  cmd+=(--species-summary-shard-dir "${dir_species_summary_shards}")
  cmd+=(--species-summary-output "${species_summary_output}")
  cmd+=(--task-stats-dir "${dir_task_stats_shards}")
  cmd+=(--aggregate-stats-output "${merge_stats_file}")
  cmd+=(--mapping-qc-output "${staged_mapping_qc}")
  cmd+=(--expected-task-count "${expected_tasks}" --task-plan "${task_plan_output}")
  cmd+=(--resolved-manifest-output "${staged_resolved_manifest}")
  echo "Running: ${cmd[*]}"
  if "${cmd[@]}"; then
    cmd_status=0
  else
    cmd_status=$?
  fi
  if [[ ${cmd_status} -ne 0 ]]; then
    stage_format_status="failed"
    echo "Failed: ${task} (exit=${cmd_status})"
    exit "${cmd_status}"
  fi

  read_stats_json_fields "${merge_stats_file}" num_species_cds=num_species_cds_files num_species_gff=num_species_gff_files num_species_genome=num_species_genome_files cds_sequences_before cds_sequences_after cds_first_sequence_name merged_rows=merged_species_summary_rows task_stats_files || exit $?
  if [[ "${task_stats_files}" != "${expected_tasks}" ]]; then
    echo "Merged task stats count does not match expected tasks: ${task_stats_files} != ${expected_tasks}"
    stage_format_status="failed"
    exit 1
  fi
  if [[ "${merged_rows}" == "0" ]]; then
    echo "Merged species summary is empty after array finalize."
    stage_format_status="failed"
    exit 1
  fi
  validate_required_formatted_outputs "${species_summary_output}" "${expected_tasks}" || {
    stage_format_status="failed"
    exit 1
  }
  stage_format_status="ok"

  run_validate_stage
  if [[ ${run_cds_fx2tab} -eq 1 ]]; then
    num_species_cds_fx2tab=$(count_nonhidden_matching_files "${species_cds_fx2tab_dir}" "*_fx2tab_cds.tsv")
    if [[ "${num_species_cds_fx2tab}" != "${expected_tasks}" ]]; then
      echo "fx2tab output count mismatch after array workers: cds=${num_species_cds_fx2tab}, expected=${expected_tasks}"
      stage_cds_fx2tab_status="failed"
      exit 1
    fi
    if ! is_species_set_identical "${species_cds_dir}" "${species_cds_fx2tab_dir}"; then
      echo "Exiting due to species-set mismatch between ${species_cds_dir} and ${species_cds_fx2tab_dir}"
      stage_cds_fx2tab_status="failed"
      exit 1
    fi
    stage_cds_fx2tab_status="ok"
  else
    stage_cds_fx2tab_status="skipped"
  fi
  if [[ ${run_species_busco} -eq 1 ]]; then
    num_species_busco_full=$(count_nonhidden_matching_files "${species_busco_full_dir}" "*busco.full.tsv")
    num_species_busco_short=$(count_nonhidden_matching_files "${species_busco_short_dir}" "*busco.short.txt")
    if [[ "${num_species_busco_full}" != "${expected_tasks}" || "${num_species_busco_short}" != "${expected_tasks}" ]]; then
      echo "BUSCO output count mismatch after array workers: full=${num_species_busco_full}, short=${num_species_busco_short}, expected=${expected_tasks}"
      stage_species_busco_status="failed"
      exit 1
    fi
    stage_species_busco_status="ok"
  else
    stage_species_busco_status="skipped"
  fi
  # Optional shared stages must also succeed before publishing canonical tables.
  run_trait_stage
  run_multispecies_summary_stage
  mv -- "${species_summary_output}" "${canonical_summary_output}"
  species_summary_output="${canonical_summary_output}"
  mv -- "${staged_mapping_qc}" "${mapping_qc_output}"
  if [[ -s "${staged_resolved_manifest}" ]]; then
    mv -- "${staged_resolved_manifest}" "${resolved_manifest_output}"
  else
    rm -f -- "${staged_resolved_manifest}"
  fi
  # Preserve the immutable plan and completion receipts for audit/retry.
  cleanup_input_generation_tmp=0
}

# Serialize shared stages against active species workers in this output workspace.
prepare_gene_model_rescue() {
  local rescue_tree="${gene_model_rescue_tree}"
  local -a rescue_args=()
  if [[ "${rescue_tree}" == auto ]]; then
    python "${gg_support_dir}/busco_guide_tree.py" build \
      --cds-dir "${species_cds_dir}" --full-dir "${species_busco_full_dir}" --short-dir "${species_busco_short_dir}" \
      --output "${gene_model_rescue_guide_dir}" --cache "${gene_model_rescue_guide_cache}" \
      --markers "${gene_model_rescue_guide_markers}" --k "${gene_model_rescue_guide_k}" \
      --sketch-size "${gene_model_rescue_guide_sketch_size}" --occupancy "${gene_model_rescue_guide_occupancy}" \
      --minimum-shared "${gene_model_rescue_guide_minimum_shared}" --nearest "${gene_model_rescue_nearest_references}" \
      --cpus "${GG_TASK_CPUS}" || return $?
    rescue_tree="${gene_model_rescue_guide_dir}/guide_tree.nwk"
    rescue_args+=(--guide-tree-receipt "${gene_model_rescue_guide_dir}/receipt.json")
  fi
  [[ -s "${rescue_tree}" ]] || {
    echo "Gene-model rescue requires an initial species tree: ${rescue_tree}" >&2
    return 1
  }
  rescue_args=(plan "${rescue_args[@]}" --cds-dir "${species_cds_dir}" --gff-dir "${species_gff_dir}"
    --genome-dir "${species_genome_dir}" --busco-dir "${species_busco_short_dir}"
    --tree "${rescue_tree}" --output "${gene_model_rescue_dir}"
    --common-references "${gene_model_rescue_common_references}" --nearest-references "${gene_model_rescue_nearest_references}"
    --minimum-busco "${gene_model_rescue_minimum_busco}" --minimum-coverage "${gene_model_rescue_minimum_coverage}"
    --minimum-identity "${gene_model_rescue_minimum_identity}" --max-interval "${gene_model_rescue_max_interval}"
    --max-intron "${gene_model_rescue_max_intron}" --genome-fallback "${gene_model_rescue_genome_fallback}"
    --max-genome-queries "${gene_model_rescue_max_genome_queries}"
    --unanchored-min-species "${gene_model_rescue_unanchored_min_species}"
    --terminal-max-extension "${gene_model_rescue_terminal_max_extension}"
    --terminal-max-unaligned-c-overhang "${gene_model_rescue_terminal_max_unaligned_c_overhang}")
  [[ -z "${gene_model_rescue_prediction_cache}" ]] || rescue_args+=(--prediction-cache "${gene_model_rescue_prediction_cache}")
  [[ -z "${gene_model_rescue_gemoma_jar}" ]] || rescue_args+=(--gemoma-jar "${gene_model_rescue_gemoma_jar}" --gemoma-java "${gene_model_rescue_gemoma_java}")
  [[ ! -s "${gg_workspace_input_dir}/species_genetic_code/species_genetic_code.tsv" ]] || \
    rescue_args+=(--genetic-codes "${gg_workspace_input_dir}/species_genetic_code/species_genetic_code.tsv")
  [[ -z "${gene_model_species_profiles}" ]] || rescue_args+=(--species-profiles "${gene_model_species_profiles}")
  python "${gg_support_dir}/rescue_gene_models.py" "${rescue_args[@]}"
}

gene_model_rescue_busco_species() {
  # First-pass BUSCO and provenance remain frozen. Re-evaluate changed species
  # into a separate namespace using the original shared lineage.
  local before_busco_full="${species_busco_full_dir}"
  local before_busco_short="${species_busco_short_dir}"
  local species_busco_full_dir="${gene_model_rescue_dir}/qc/species_cds_busco_full"
  local species_busco_short_dir="${gene_model_rescue_dir}/qc/species_cds_busco_short"
  local input_generation_provenance_dir="${gene_model_rescue_dir}/qc/artifact_provenance"
  local run_species_busco=1
  local busco_lineage_resolved
  local rescue_species=$1 rescue_cds=$2 changed=$3
  busco_lineage_resolved=$(python - "${gene_model_rescue_dir}/plan.json" <<'PY'
import json
import sys
plan = json.load(open(sys.argv[1]))
print(next(iter(plan['request']['sources'].values()))['quality']['lineage'])
PY
)
  ensure_dir "${species_busco_full_dir}"
  ensure_dir "${species_busco_short_dir}"
  ensure_dir "${input_generation_provenance_dir}"
  if [[ "${changed}" == 0 ]]; then
    cp -- "${before_busco_full}/${rescue_species}.busco.full.tsv" "${species_busco_full_dir}/${rescue_species}.busco.full.tsv"
    cp -- "${before_busco_short}/${rescue_species}.busco.short.txt" "${species_busco_short_dir}/${rescue_species}.busco.short.txt"
  else
    run_species_busco_for_one_file "${rescue_cds}" "${rescue_species}" "${GG_TASK_CPUS}"
  fi
}

annotate_gene_model_rescue_swissprot() {
  [[ ${run_gene_model_rescue_swissprot} -eq 1 ]] || return 0
  local rescue_root=$1 uniprot_prefix uniprot_meta
  local evidence_dir="${gene_model_rescue_swissprot_dir:-${rescue_root%/}.swissprot}"
  uniprot_prefix=$(ensure_uniprot_sprot_mmseqs_db "${gg_workspace_dir}") || return $?
  uniprot_meta=$(ensure_uniprot_sprot_metadata_tsv "${gg_workspace_dir}" "${uniprot_prefix}") || return $?
  python "${gg_support_dir}/rescue_swissprot_evidence.py" --rescue-output "${rescue_root}" \
    --output "${evidence_dir}" --db-prefix "${uniprot_prefix}" --metadata "${uniprot_meta}" \
    --cache "${gg_workspace_dir}/downloads/rescue_swissprot_cache" --cpus "${GG_TASK_CPUS}" \
    --memory-gb "$(gg_memory_fraction_gb "${GG_MEM_TOOL_GB}" 3 4)"
}

finish_gene_model_rescue() {
  python "${gg_support_dir}/rescue_gene_models.py" finalize --output "${gene_model_rescue_dir}" --cpus "${GG_TASK_CPUS}"
  local rescue_index rescue_species rescue_cds changed qc_complete rescue_qc_rows
  rescue_qc_rows=$(python "${gg_support_dir}/rescue_gene_models.py" qc-inputs --output "${gene_model_rescue_dir}" --cpus "${GG_TASK_CPUS}") || return $?
  while IFS=$'\t' read -r rescue_index rescue_species rescue_cds changed qc_complete; do
    if [[ "${qc_complete}" != 1 ]]; then
      gene_model_rescue_busco_species "${rescue_species}" "${rescue_cds}" "${changed}"
      python "${gg_support_dir}/rescue_gene_models.py" worker-complete --output "${gene_model_rescue_dir}" --task-index "${rescue_index}" --cpus "${GG_TASK_CPUS}"
    fi
  done <<< "${rescue_qc_rows}"
  python "${gg_support_dir}/rescue_gene_models.py" qc --output "${gene_model_rescue_dir}" --busco-dir "${gene_model_rescue_dir}/qc/species_cds_busco_short" --cpus "${GG_TASK_CPUS}"
  annotate_gene_model_rescue_swissprot "${gene_model_rescue_dir}"
  echo "Augmented CDS/GFF inputs: ${gene_model_rescue_dir}/augmented/inputs.tsv"
}

prepare_gene_model_refinement() {
  local -a refinement_args=(plan --output "${gene_model_refinement_dir}"
    --policy "${gene_model_refinement_policy}" --mode "${gene_model_refinement_mode}"
    --isoform-adoption "${gene_model_refinement_isoform_adoption}"
    --min-margin "${gene_model_refinement_min_margin}" --min-support "${gene_model_refinement_min_support}"
    --candidate-limit "${gene_model_refinement_candidate_limit}" --padding "${gene_model_refinement_padding}"
    --minimum-coverage "${gene_model_rescue_minimum_coverage}" --minimum-identity "${gene_model_rescue_minimum_identity}"
    --max-intron "${gene_model_rescue_max_intron}" --max-interval "${gene_model_rescue_max_interval}")
  if [[ -n "${gene_model_refinement_inputs}" ]]; then
    refinement_args+=(--inputs "${gene_model_refinement_inputs}")
  else
    local anchor_dir="${gene_model_refinement_rescue_dir:-${gene_model_rescue_dir}}"
    if [[ ! -s "${anchor_dir}/plan.json" ]]; then
      [[ -z "${gene_model_refinement_rescue_dir}" ]] || { echo "Frozen rescue plan missing: ${anchor_dir}" >&2; return 1; }
      prepare_gene_model_rescue
    fi
    refinement_args+=(--rescue-output "${anchor_dir}")
  fi
  [[ -z "${gene_model_refinement_edges}" ]] || refinement_args+=(--edges "${gene_model_refinement_edges}")
  [[ -z "${gene_model_refinement_rna}" ]] || refinement_args+=(--rna "${gene_model_refinement_rna}")
  [[ -z "${gene_model_species_profiles}" ]] || refinement_args+=(--species-profiles "${gene_model_species_profiles}")
  python "${gg_support_dir}/gene_model_refinement.py" "${refinement_args[@]}"
}

finish_gene_model_refinement() {
  python "${gg_support_dir}/gene_model_refinement.py" finalize --output "${gene_model_refinement_dir}" --cpus "${GG_TASK_CPUS}"
  python "${gg_support_dir}/gene_model_refinement.py" qc --output "${gene_model_refinement_dir}"
  local refinement_review_dir="${gene_model_refinement_dir%/}.review"
  local anchor_dir="${gene_model_refinement_rescue_dir:-${gene_model_rescue_dir}}"
  local -a swissprot_plot_args=()
  local -a busco_phase_args=()
  if [[ -z "${gene_model_refinement_inputs}" && -s "${anchor_dir}/augmented/receipt.json" ]]; then
    busco_phase_args+=(--three-stage)
  fi
  if [[ ${run_gene_model_rescue_swissprot} -eq 1 && -z "${gene_model_refinement_inputs}" && -s "${anchor_dir}/augmented/receipt.json" ]]; then
    annotate_gene_model_rescue_swissprot "${anchor_dir}"
    swissprot_plot_args+=(--rescue-swissprot-dir "${gene_model_rescue_swissprot_dir:-${anchor_dir%/}.swissprot}")
  fi
  python "${gg_support_dir}/plot_gene_model_refinement.py" --output "${gene_model_refinement_dir}" \
    --report "${refinement_review_dir}" --cds-dir "${species_cds_dir}"
  if [[ ${run_species_busco} -eq 1 ]]; then
    local refinement_busco_db refinement_busco_jobs=1 refinement_busco_memory_cap
    ensure_shared_busco_lineage_ready "${task_plan_output}"
    refinement_busco_db=$(ensure_busco_download_path "${gg_workspace_dir}" "${busco_lineage_resolved}")
    if [[ "${species_busco_parallel_jobs}" == auto ]]; then
      refinement_busco_jobs=${GG_TASK_CPUS}
      [[ ${refinement_busco_jobs} -le 4 ]] || refinement_busco_jobs=4
    elif [[ "${species_busco_parallel_jobs}" =~ ^[1-9][0-9]*$ ]]; then
      refinement_busco_jobs=${species_busco_parallel_jobs}
      [[ ${refinement_busco_jobs} -le ${GG_TASK_CPUS} ]] || refinement_busco_jobs=${GG_TASK_CPUS}
    else
      echo "Invalid species_busco_parallel_jobs: ${species_busco_parallel_jobs}" >&2
      return 2
    fi
    refinement_busco_memory_cap=$(gg_memory_parallel_job_cap "${GG_MEM_TOOL_GB}" "${species_busco_memory_gb_per_job}")
    [[ ${refinement_busco_jobs} -le ${refinement_busco_memory_cap} ]] || refinement_busco_jobs=${refinement_busco_memory_cap}
    python "${gg_support_dir}/gene_model_refinement_busco.py" --output "${gene_model_refinement_dir}" \
      --report "${refinement_review_dir}/busco" --cds-dir "${species_cds_dir}" \
      "${swissprot_plot_args[@]}" \
      "${busco_phase_args[@]}" \
      --lineage "${refinement_busco_db}/lineages/${busco_lineage_resolved}" --download-path "${refinement_busco_db}" \
      --jobs "${refinement_busco_jobs}" --cpus "$((GG_TASK_CPUS / refinement_busco_jobs))"
  fi
  echo "Selected CDS/protein/GFF inputs: ${gene_model_refinement_dir}/effective/inputs.tsv"
  echo "Refinement review and paired BUSCO comparison: ${refinement_review_dir}"
}

ensure_dir "${input_generation_root}"
export GG_PERFORMANCE_DIR="${input_generation_root}/tmp/performance/$$"
array_lock_mode=exclusive
if [[ "${input_generation_mode}" == array_worker || "${input_generation_mode}" == rescue_synteny || "${input_generation_mode}" == rescue_models || "${input_generation_mode}" == refinement_catalog || "${input_generation_mode}" == refinement_predict ]]; then
  array_lock_mode=shared
fi
input_generation_lock "${input_generation_root}/.array-phase.lock" "${array_lock_mode}"
if [[ ${run_gene_model_rescue} -eq 1 || ${run_gene_model_refinement} -eq 1 || "${input_generation_mode}" == rescue_* || "${input_generation_mode}" == refinement_* ]]; then
  ensure_dir "${gene_model_rescue_dir}"
  input_generation_lock "${gene_model_rescue_dir}/.array-phase.lock" "${array_lock_mode}"
fi

if [[ ${run_gene_model_refinement} -eq 1 || "${input_generation_mode}" == refinement_* ]]; then
  ensure_dir "${gene_model_refinement_dir}"
  input_generation_lock "${gene_model_refinement_dir}/.array-phase.lock" "${array_lock_mode}"
fi

# Custom output directories can be shared across workspace paths. Lock their
# canonical locations as well, and bind array outputs to one plan at prepare.
array_output_args=(--file "${species_cds_dir}" --file "${species_gff_dir}" --file "${species_genome_dir}")
[[ ${run_cds_fx2tab} -ne 1 ]] || array_output_args+=(--file "${species_cds_fx2tab_dir}")
[[ ${run_species_busco} -ne 1 ]] || array_output_args+=(--file "${species_busco_full_dir}" --file "${species_busco_short_dir}")
array_output_locks=$(python "${gg_support_dir}/input_generation_array_state.py" output-lock-paths --task-plan "${task_plan_output}" "${array_output_args[@]}")
while IFS= read -r array_output_lock; do
  input_generation_lock "${array_output_lock}" "${array_lock_mode}"
done <<< "${array_output_locks}"

if [[ "${input_generation_mode}" == array_* ]]; then
  array_settings_cmd=(python "${gg_support_dir}/input_generation_array_state.py" configure --task-plan "${task_plan_output}")
  for array_setting in provider download_limit_dir gene_grouping_mode gff_repair_mode strict busco_lineage busco_timeout_seconds \
    run_validate_inputs run_cds_fx2tab run_species_busco run_generate_species_trait run_multispecies_summary \
    run_species_taxonomy taxonomy_species_tree taxonomy_ranks taxonomy_plot_clades taxonomy_taxid_map \
    species_cds_dir species_gff_dir species_genome_dir species_cds_fx2tab_dir species_busco_full_dir species_busco_short_dir \
    species_summary_output resolved_manifest_output species_trait_output file_multispecies_summary \
    trait_profile trait_species_source trait_databases trait_plan trait_database_sources trait_download_dir trait_download_timeout \
    gbif_api gbif_page_size gbif_max_occurrences_per_species gbif_grid_degrees gbif_min_match_confidence \
    gbif_max_coordinate_uncertainty_m gbif_min_distance_from_known_centroid_m \
    gbif_year_min gbif_year_max gbif_countries gbif_include_basis_of_record gbif_exclude_basis_of_record gbif_include_establishment_means gbif_missing_date gbif_missing_uncertainty gbif_missing_centroid_distance gbif_use_cache gbif_require_complete gbif_occurrence_file gbif_taxon_map gbif_download_metadata; do
    array_settings_cmd+=(--setting "${array_setting}=${!array_setting}")
  done
  array_settings_cmd+=(--setting "genetic_code=${GG_COMMON_GENETIC_CODE:-1}")
  # Do not add empty resume fields to old immutable settings documents.
  for resume_prefix in resume_from resume_fallback; do
    donor_plan=${resume_prefix}_task_plan
    donor_sha=${resume_prefix}_task_plan_sha256
    donor_root=${resume_prefix}_input_generation_root
    if [[ -n "${!donor_plan}${!donor_sha}${!donor_root}" ]]; then
      for array_setting in "${donor_plan}" "${donor_sha}" "${donor_root}"; do
        array_settings_cmd+=(--setting "${array_setting}=${!array_setting}")
      done
    fi
  done
  [[ ${require_cds} -ne 1 ]] || array_settings_cmd+=(--setting "require_cds=1")
  [[ ${require_gff} -ne 1 ]] || array_settings_cmd+=(--setting "require_gff=1")
  [[ ${require_genome} -ne 1 ]] || array_settings_cmd+=(--setting "require_genome=1")
  if [[ ${run_species_taxonomy} -eq 1 ]]; then
    [[ -z "${taxonomy_taxid_map}" ]] || array_settings_cmd+=(--file "${taxonomy_taxid_map}")
    [[ "${taxonomy_species_tree}" == auto ]] || array_settings_cmd+=(--file "${taxonomy_species_tree}")
  fi
  if [[ ${run_generate_species_trait} -eq 1 ]]; then
    array_settings_cmd+=(--file "${trait_plan}" --file "${trait_database_sources}")
    gbif_input_files=$(python "${gg_support_dir}/generate_species_trait.py" \
      --database-sources "${trait_database_sources}" --trait-plan "${trait_plan}" --databases "${trait_databases}" \
      --gbif-occurrence-file "${gbif_occurrence_file}" --gbif-taxon-map "${gbif_taxon_map}" \
      --gbif-download-metadata "${gbif_download_metadata}" --print-gbif-input-files) || exit $?
    while IFS= read -r gbif_input_file; do
      [[ -z "${gbif_input_file}" ]] || array_settings_cmd+=(--file "${gbif_input_file}")
    done <<< "${gbif_input_files}"
  fi
  [[ "${input_generation_mode}" != array_prepare ]] || array_settings_cmd+=(--prepare)
  "${array_settings_cmd[@]}"
  if [[ "${input_generation_mode}" != array_prepare ]]; then
    prepare_check_cmd=(python "${gg_support_dir}/input_generation_array_state.py" check-prepared --task-plan "${task_plan_output}")
    [[ "${input_generation_mode}" != array_worker ]] || prepare_check_cmd+=(--task-index "${GG_ARRAY_TASK_ID}")
    "${prepare_check_cmd[@]}" || {
      echo "Array prepare did not complete for this plan and settings. Run array_prepare first."; exit 1;
    }
    python "${gg_support_dir}/input_generation_array_state.py" claim-workspace --task-plan "${task_plan_output}" --workspace "${input_generation_root}" "${array_output_args[@]}"
  fi
fi

case "${input_generation_mode}" in
  single)
    run_format_stage_single
    if [[ ${download_only} -eq 0 && ${dry_run} -eq 0 ]]; then
      validate_required_formatted_outputs "${species_summary_output}"
    fi
    run_validate_stage
    run_cds_fx2tab_stage_all
    run_species_busco_stage_all
    run_trait_stage
    run_multispecies_summary_stage
    cleanup_input_generation_tmp=1
    ;;
  array_prepare)
    run_array_prepare_mode
    ;;
  array_worker)
    run_array_worker_mode
    ;;
  array_finalize)
    run_array_finalize_mode
    ;;
  refinement_prepare)
    prepare_gene_model_refinement
    ;;
  refinement_catalog|refinement_predict|refinement_correspondence)
    write_run_summary_on_exit=0
    refinement_command="${input_generation_mode#refinement_}"
    refinement_args=(--output "${gene_model_refinement_dir}" --cpus "${GG_TASK_CPUS}")
    if [[ "${refinement_command}" != correspondence ]]; then
      [[ "${GG_ARRAY_TASK_ID}" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid refinement worker index" >&2; exit 1; }
      refinement_args+=(--task-index "${GG_ARRAY_TASK_ID}")
    fi
    [[ -z "${gene_model_rescue_comparison_cache}" ]] || refinement_args+=(--comparison-cache "${gene_model_rescue_comparison_cache}")
    python "${gg_support_dir}/gene_model_refinement.py" "${refinement_command}" "${refinement_args[@]}"
    if [[ "${refinement_command}" == correspondence ]]; then
      python "${gg_support_dir}/gene_model_refinement.py" select --output "${gene_model_refinement_dir}"
    fi
    ;;
  refinement_finalize)
    finish_gene_model_refinement
    ;;
  rescue_prepare)
    prepare_gene_model_rescue
    ;;
  rescue_synteny|rescue_models)
    write_run_summary_on_exit=0
    if [[ "${input_generation_mode}" == rescue_models ]]; then
      [[ "${GG_ARRAY_TASK_ID}" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid rescue worker index" >&2; exit 1; }
      ensure_dir "${gene_model_rescue_dir}/core-worker-locks"
      input_generation_lock "${gene_model_rescue_dir}/core-worker-locks/${GG_ARRAY_TASK_ID}.lock" exclusive
    fi
    rescue_subcommand=synteny
    [[ "${input_generation_mode}" != rescue_models ]] || rescue_subcommand=rescue
    rescue_execution_args=()
    if [[ "${rescue_subcommand}" == synteny && -n "${gene_model_rescue_comparison_cache:-}" ]]; then
      rescue_execution_args+=(--comparison-cache "${gene_model_rescue_comparison_cache}")
    fi
    if [[ "${rescue_subcommand}" == rescue && "${gene_model_rescue_interval_workers:-0}" != 0 ]]; then
      rescue_execution_args+=(--interval-workers "${gene_model_rescue_interval_workers}")
    fi
    python "${gg_support_dir}/rescue_gene_models.py" "${rescue_subcommand}" --output "${gene_model_rescue_dir}" \
      --task-index "${GG_ARRAY_TASK_ID}" --cpus "${GG_TASK_CPUS}" "${rescue_execution_args[@]}"
    if [[ "${input_generation_mode}" == rescue_models ]]; then
      rescue_qc_rows=$(python "${gg_support_dir}/rescue_gene_models.py" qc-inputs --output "${gene_model_rescue_dir}" --task-index "${GG_ARRAY_TASK_ID}") || exit $?
      while IFS=$'\t' read -r rescue_index rescue_species rescue_cds changed qc_complete; do
        gene_model_rescue_busco_species "${rescue_species}" "${rescue_cds}" "${changed}"
      done <<< "${rescue_qc_rows}"
      python "${gg_support_dir}/rescue_gene_models.py" worker-complete --output "${gene_model_rescue_dir}" --task-index "${GG_ARRAY_TASK_ID}"
    fi
    ;;
  rescue_finalize)
    finish_gene_model_rescue
    ;;
esac

# Species taxonomy uses the current input set and preserves completed output on failure.
if [[ ${run_species_taxonomy} -eq 1 && ( "${input_generation_mode}" == single || "${input_generation_mode}" == array_finalize ) && ${dry_run} -ne 1 && ${download_only} -ne 1 ]]; then
  ensure_ete_taxonomy_db "${gg_workspace_dir}" || exit 1
  taxonomy_summary_args=()
  [[ ! -s "${species_summary_output}" ]] || taxonomy_summary_args+=(--species-summary "${species_summary_output}")
  python "${gg_support_dir}/species_taxonomy.py" \
    --workspace "${gg_workspace_dir}" \
    --species-tree "${taxonomy_species_tree}" \
    --ranks "${taxonomy_ranks}" \
    --plot-clades "${taxonomy_plot_clades}" \
    --taxid-map "${taxonomy_taxid_map}" \
    --taxid-override "${taxonomy_taxid_override}" \
    --species-dir "${species_cds_dir}" \
    "${taxonomy_summary_args[@]}" || exit $?
fi

if [[ ${run_gene_model_rescue} -eq 1 && ( "${input_generation_mode}" == single || "${input_generation_mode}" == array_finalize ) && ${dry_run} -ne 1 && ${download_only} -ne 1 ]]; then
  prepare_gene_model_rescue
  if [[ "${input_generation_mode}" == single ]]; then
    rescue_execution_args=()
    [[ -z "${gene_model_rescue_comparison_cache:-}" ]] || rescue_execution_args+=(--comparison-cache "${gene_model_rescue_comparison_cache}")
    [[ "${gene_model_rescue_interval_workers:-0}" == 0 ]] || rescue_execution_args+=(--interval-workers "${gene_model_rescue_interval_workers}")
    python "${gg_support_dir}/rescue_gene_models.py" run --output "${gene_model_rescue_dir}" --cpus "${GG_TASK_CPUS}" "${rescue_execution_args[@]}"
    finish_gene_model_rescue
  fi
fi

if [[ ${run_gene_model_refinement} -eq 1 && ( "${input_generation_mode}" == single || "${input_generation_mode}" == array_finalize ) && ${dry_run} -ne 1 && ${download_only} -ne 1 ]]; then
  # A combined array run freezes refinement only after rescue finalization;
  # the augmented inputs do not exist at the initial formatting finalizer.
  if [[ "${input_generation_mode}" == single || ${run_gene_model_rescue} -ne 1 ]]; then
    prepare_gene_model_refinement
  fi
  if [[ "${input_generation_mode}" == single ]]; then
    refinement_execution_args=()
    [[ -z "${gene_model_rescue_comparison_cache}" ]] || refinement_execution_args+=(--comparison-cache "${gene_model_rescue_comparison_cache}")
    python "${gg_support_dir}/gene_model_refinement.py" run --output "${gene_model_refinement_dir}" --cpus "${GG_TASK_CPUS}" "${refinement_execution_args[@]}"
    finish_gene_model_refinement
  fi
fi

if [[ ${cleanup_input_generation_tmp} -eq 1 && -d "${download_tmp_root}" ]]; then
  echo "Removing temporary input_generation directory: ${download_tmp_root}"
  rm -rf -- "${download_tmp_root}"
fi

echo "$(date): Exiting Singularity environment"

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
run_species_busco="${run_species_busco:-1}"
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
gbif_max_distance_from_centroid_m="${gbif_max_distance_from_centroid_m:-}"
gene_grouping_mode="${gene_grouping_mode:-rescue_overlap}"
gff_repair_mode="${gff_repair_mode:-safe}"

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
  single|array_prepare|array_worker|array_finalize) ;;
  *)
    echo "Invalid input_generation_mode: ${input_generation_mode} (allowed: single|array_prepare|array_worker|array_finalize)"
    exit 1
    ;;
esac

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
if [[ "${input_generation_mode}" != "single" && ${run_format_inputs} -ne 1 ]]; then
  echo "run_format_inputs must be 1 when input_generation_mode=${input_generation_mode}"
  exit 1
fi

input_generation_root="${gg_workspace_output_dir}/input_generation"
input_generation_tmp_root="${input_generation_root}/tmp"
input_generation_provenance_dir="${input_generation_root}/artifact_provenance"
download_tmp_root="${input_generation_tmp_root}"
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

if [[ -z "${download_manifest}" ]]; then
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
  local exit_code=$?
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

trap write_gg_input_generation_summary_on_exit EXIT

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

  if [[ ${format_force_overwrite} -eq 1 && -s "${format_provenance_manifest}" ]]; then
    echo "Clearing managed formatted-input outputs before provenance rebuild."
    rm -rf -- "${species_cds_dir}" "${species_gff_dir}" "${species_genome_dir}"
    rm -f -- "${species_summary_output}" "${resolved_manifest_output}"
  fi

  if ! ensure_ete_taxonomy_db "${gg_workspace_dir}"; then
    echo "Warning: Failed to prepare ETE taxonomy DB for species_summary taxonomy metadata. Continuing without taxid/genetic code annotation." >&2
  fi

  format_stats_file="${download_tmp_root}/gg_input_generation_stats.json"
  ensure_parent_dir "${format_stats_file}"
  rm -f -- "${format_stats_file}"
  cmd=(python "${gg_support_dir}/format_species_inputs.py")
  cmd+=(--provider "${provider}")
  cmd+=(--species-cds-dir "${species_cds_dir}")
  cmd+=(--species-gff-dir "${species_gff_dir}")
  cmd+=(--species-genome-dir "${species_genome_dir}")
  cmd+=(--species-summary-output "${species_summary_output}")
  cmd+=(--stats-output "${format_stats_file}")
  cmd+=(--gene-grouping-mode "${gene_grouping_mode}")
  cmd+=(--gff-repair-mode "${gff_repair_mode}")

  if [[ -n "${download_manifest}" ]]; then
    cmd+=(--download-manifest "${download_manifest}")
    cmd+=(--download-dir "${download_dir}")
    cmd+=(--resolved-manifest-output "${resolved_manifest_output}")
  fi
  if [[ -n "${input_dir}" ]]; then
    cmd+=(--input-dir "${input_dir}")
  fi
  if [[ ${format_force_overwrite} -eq 1 ]]; then
    cmd+=(--overwrite)
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

  if [[ -s "${format_stats_file}" ]]; then
    num_species_cds="$(read_stats_json_field "${format_stats_file}" "num_species_cds_files")"
    num_species_gff="$(read_stats_json_field "${format_stats_file}" "num_species_gff_files")"
    num_species_genome="$(read_stats_json_field "${format_stats_file}" "num_species_genome_files")"
    cds_sequences_before="$(read_stats_json_field "${format_stats_file}" "cds_sequences_before")"
    cds_sequences_after="$(read_stats_json_field "${format_stats_file}" "cds_sequences_after")"
    cds_first_sequence_name="$(read_stats_json_field "${format_stats_file}" "cds_first_sequence_name")"
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

  if [[ ${run_validate_inputs} -ne 1 || ${run_format_inputs} -ne 1 || ${download_only} -ne 0 || ${dry_run} -ne 0 ]]; then
    gg_step_skip "${task}"
    stage_validate_status="skipped"
    return 0
  fi

  gg_step_start "${task}"
  stage_validate_status="running"
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
      cmd=(python "${gg_support_dir}/validate_cds_gff_mapping.py")
      cmd+=(--species-cds-dir "${species_cds_dir}")
      cmd+=(--species-gff-dir "${species_gff_dir}")
      cmd+=(--species-summary "${species_summary_output}")
      cmd+=(--nthreads "${GG_TASK_CPUS:-1}")
      cmd+=(--stats-output "${mapping_stats_file}")
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
      echo "Running: ${cmd[*]}"
      if "${cmd[@]}"; then
        :
      else
        stage_validate_status="failed"
        echo "Failed: Validate longest CDS representative selection"
        exit 1
      fi
      rm -f -- "${mapping_stats_file}"
      rm -f -- "${longest_cds_stats_file}"
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

  gg_step_start "${task}"
  stage_validate_status="running"
  rm -f -- "${mapping_stats_file}" "${longest_stats_file}"
  cmd=(python "${gg_support_dir}/validate_cds_gff_mapping.py")
  cmd+=(--species-cds-dir "${species_cds_dir}")
  cmd+=(--species-gff-dir "${species_gff_dir}")
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
  rm -f -- "${mapping_stats_file}" "${longest_stats_file}"
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
  tmp_fx2tab_tsv=$(mktemp "${input_generation_tmp_root}/fx2tab.${species_name}.XXXXXX.tsv")
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
  species_prefix=$(read_stats_json_field "${task_meta_file}" "species_prefix")
  cds_output_path=$(read_stats_json_field "${task_meta_file}" "cds_output_path")
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
  busco_work_root=$(mktemp -d "${input_generation_tmp_root}/busco.${species_name}.XXXXXX")
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
    gg_run_busco_with_metaeuk_modified_fas_compat \
      --in "input.fasta" \
      --mode "transcriptome" \
      --out "busco_tmp" \
      --cpu "${busco_threads}" \
      --force \
      --evalue 1e-03 \
      --limit 20 \
      --lineage_dataset "${dir_busco_lineage}" \
      --download_path "${dir_busco_db}" \
      --offline
  )

  if copy_busco_tables "${busco_output_dir}" "${busco_lineage_resolved}" "${file_sp_busco_full}" "${file_sp_busco_short}"; then
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
  species_prefix=$(read_stats_json_field "${task_meta_file}" "species_prefix")
  cds_output_path=$(read_stats_json_field "${task_meta_file}" "cds_output_path")
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
  if [[ -z "${trait_manifest_path}" ]]; then
    trait_manifest_default="${gg_workspace_input_dir}/input_generation/download_plan.xlsx"
    if [[ -s "${trait_manifest_default}" ]]; then
      trait_manifest_path="${trait_manifest_default}"
    fi
  fi

  gg_artifact_contract_init trait_provenance_args "input_generation_species_trait" "all_species" "${input_generation_provenance_dir}/species_trait.json"
  gg_artifact_add_input_if_present trait_provenance_args "download_manifest" "${trait_manifest_path}"
  if [[ "${trait_species_source}" == "species_cds" ]]; then
    gg_artifact_add_input_if_present trait_provenance_args "species_cds" "${species_cds_dir}"
  fi
  gg_artifact_add_input_if_present trait_provenance_args "trait_plan" "${trait_plan}"
  gg_artifact_add_input_if_present trait_provenance_args "database_sources" "${trait_database_sources}"
  trait_provenance_args+=(
    --output "species_trait=${species_trait_output}"
    --parameter "trait_profile=${trait_profile}"
    --parameter "trait_species_source=${trait_species_source}"
    --parameter "trait_databases=${trait_databases}"
    --parameter "gbif_api=${gbif_api}"
    --parameter "gbif_page_size=${gbif_page_size}"
    --parameter "gbif_max_occurrences_per_species=${gbif_max_occurrences_per_species}"
    --parameter "gbif_grid_degrees=${gbif_grid_degrees}"
    --parameter "gbif_min_match_confidence=${gbif_min_match_confidence}"
    --parameter "gbif_max_coordinate_uncertainty_m=${gbif_max_coordinate_uncertainty_m}"
    --parameter "gbif_max_distance_from_centroid_m=${gbif_max_distance_from_centroid_m}"
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
  if [[ -n "${gbif_api}" ]]; then
    cmd+=(--gbif-api "${gbif_api}")
  fi
  if [[ -n "${gbif_page_size}" ]]; then
    cmd+=(--gbif-page-size "${gbif_page_size}")
  fi
  if [[ -n "${gbif_max_occurrences_per_species}" ]]; then
    cmd+=(--gbif-max-occurrences-per-species "${gbif_max_occurrences_per_species}")
  fi
  if [[ -n "${gbif_grid_degrees}" ]]; then
    cmd+=(--gbif-grid-degrees "${gbif_grid_degrees}")
  fi
  if [[ -n "${gbif_min_match_confidence}" ]]; then
    cmd+=(--gbif-min-match-confidence "${gbif_min_match_confidence}")
  fi
  if [[ -n "${gbif_max_coordinate_uncertainty_m}" ]]; then
    cmd+=(--gbif-max-coordinate-uncertainty-m "${gbif_max_coordinate_uncertainty_m}")
  fi
  if [[ -n "${gbif_max_distance_from_centroid_m}" ]]; then
    cmd+=(--gbif-max-distance-from-centroid-m "${gbif_max_distance_from_centroid_m}")
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
    num_species_trait="$(read_stats_json_field "${trait_stats_file}" "num_species_with_any_trait")"
    num_trait_columns="$(read_stats_json_field "${trait_stats_file}" "num_trait_columns")"
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

  write_run_summary_on_exit=0
  prepare_input_generation_tmp_dirs
  clean_input_generation_shards

  gg_step_start "${task}"
  stage_format_status="running"

  if ! ensure_selected_download_manifest_has_rows; then
    stage_format_status="failed"
    exit 1
  fi

  if [[ -n "${download_manifest}" ]]; then
    cmd=(python "${gg_support_dir}/format_species_inputs.py")
    cmd+=(--provider "${provider}")
    cmd+=(--download-manifest "${download_manifest}")
    cmd+=(--download-dir "${download_dir}")
    cmd+=(--resolved-manifest-output "${resolved_manifest_output}")
    cmd+=(--download-only)
    cmd+=(--download-timeout "${download_timeout}")
    cmd+=(--jobs "${GG_TASK_CPUS:-1}")
    cmd+=(--gene-grouping-mode "${gene_grouping_mode}")
    cmd+=(--gff-repair-mode "${gff_repair_mode}")
    if [[ ${overwrite} -eq 1 ]]; then
      cmd+=(--overwrite)
    fi
    if [[ ${strict} -eq 1 ]]; then
      cmd+=(--strict)
    fi
    if [[ -n "${auth_bearer_token_env}" ]]; then
      cmd+=(--auth-bearer-token-env "${auth_bearer_token_env}")
    fi
    if [[ -n "${http_header}" ]]; then
      cmd+=(--http-header "${http_header}")
    fi
    echo "Running: ${cmd[*]}"
    if "${cmd[@]}"; then
      cmd_status=0
    else
      cmd_status=$?
    fi
    if [[ ${cmd_status} -ne 0 ]]; then
      stage_format_status="failed"
      echo "Failed: download-manifest prepare stage (exit=${cmd_status})"
      exit "${cmd_status}"
    fi
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
  cmd+=(--gene-grouping-mode "${gene_grouping_mode}")
  cmd+=(--gff-repair-mode "${gff_repair_mode}")
  if [[ ${strict} -eq 1 ]]; then
    cmd+=(--strict)
  fi
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

  expected_tasks=$(task_plan_task_count "${task_plan_output}")
  num_species_cds="${expected_tasks}"
  num_species_gff="${expected_tasks}"
  num_species_genome=""

  if [[ ${run_species_busco} -eq 1 ]]; then
    ensure_parent_dir "${file_busco_lineage_resolved}"
    resolve_busco_lineage_from_task_plan "${task_plan_output}"
  fi
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
  if [[ ! "${GG_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] || [[ ${GG_ARRAY_TASK_ID} -lt 1 ]]; then
    echo "Invalid GG_ARRAY_TASK_ID value (must be a positive integer): ${GG_ARRAY_TASK_ID}"
    exit 1
  fi

  gg_step_start "${task}"
  stage_format_status="running"
  task_stats_file="${dir_task_stats_shards}/${GG_ARRAY_TASK_ID}.json"
  task_meta_file="${dir_task_meta_shards}/${GG_ARRAY_TASK_ID}.json"
  task_summary_file="${dir_species_summary_shards}/${GG_ARRAY_TASK_ID}.tsv"
  rm -f -- "${task_stats_file}" "${task_meta_file}" "${task_summary_file}"

  describe_cmd=(python "${gg_support_dir}/run_input_generation_task.py")
  describe_cmd+=(--task-plan "${task_plan_output}")
  describe_cmd+=(--task-index "${GG_ARRAY_TASK_ID}")
  describe_cmd+=(--species-cds-dir "${species_cds_dir}")
  describe_cmd+=(--species-gff-dir "${species_gff_dir}")
  describe_cmd+=(--species-genome-dir "${species_genome_dir}")
  describe_cmd+=(--task-meta-output "${task_meta_file}")
  describe_cmd+=(--describe-only)
  if ! "${describe_cmd[@]}"; then
    stage_format_status="failed"
    echo "Failed to describe input-generation array task ${GG_ARRAY_TASK_ID}."
    exit 1
  fi
  species_prefix=$(read_stats_json_field "${task_meta_file}" "species_prefix")
  cds_input_path=$(read_stats_json_field "${task_meta_file}" "cds_path")
  gff_input_path=$(read_stats_json_field "${task_meta_file}" "gff_path")
  gbff_input_path=$(read_stats_json_field "${task_meta_file}" "gbff_path")
  genome_input_path=$(read_stats_json_field "${task_meta_file}" "genome_path")
  cds_output_path=$(read_stats_json_field "${task_meta_file}" "cds_output_path")
  gff_output_path=$(read_stats_json_field "${task_meta_file}" "gff_output_path")
  genome_output_path=$(read_stats_json_field "${task_meta_file}" "genome_output_path")
  if [[ -z "${species_prefix}" || -z "${cds_output_path}" || -z "${gff_output_path}" ]]; then
    stage_format_status="failed"
    echo "Array task description is missing required species or output paths: ${task_meta_file}"
    exit 1
  fi
  format_provenance_manifest="${input_generation_provenance_dir}/format.${species_prefix}.json"
  gg_artifact_contract_init format_provenance_args "input_generation_format" "${species_prefix}" "${format_provenance_manifest}"
  gg_artifact_add_input_if_present format_provenance_args "cds_input" "${cds_input_path}"
  gg_artifact_add_input_if_present format_provenance_args "gff_input" "${gff_input_path}"
  gg_artifact_add_input_if_present format_provenance_args "gbff_input" "${gbff_input_path}"
  gg_artifact_add_input_if_present format_provenance_args "genome_input" "${genome_input_path}"
  format_provenance_args+=(
    --output "formatted_cds=${cds_output_path}"
    --output "formatted_gff=${gff_output_path}"
    --parameter "provider=${provider}"
    --parameter "gene_grouping_mode=${gene_grouping_mode}"
    --parameter "gff_repair_mode=${gff_repair_mode}"
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
  if [[ ${format_needs_update} -eq 1 ]]; then
    gg_artifact_record "${format_provenance_args[@]}"
  fi

  if [[ -s "${task_stats_file}" ]]; then
    num_species_cds="$(read_stats_json_field "${task_stats_file}" "num_species_cds_files")"
    num_species_gff="$(read_stats_json_field "${task_stats_file}" "num_species_gff_files")"
    num_species_genome="$(read_stats_json_field "${task_stats_file}" "num_species_genome_files")"
    cds_sequences_before="$(read_stats_json_field "${task_stats_file}" "cds_sequences_before")"
    cds_sequences_after="$(read_stats_json_field "${task_stats_file}" "cds_sequences_after")"
    cds_first_sequence_name="$(read_stats_json_field "${task_stats_file}" "cds_first_sequence_name")"
  fi
  stage_format_status="ok"

  run_validate_stage_one_worker
  run_cds_fx2tab_stage_one_worker
  run_species_busco_stage_one_worker
  stage_multispecies_summary_status="skipped"
  stage_trait_status="skipped"
}

run_array_finalize_mode() {
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
  cmd+=(--expected-task-count "${expected_tasks}")
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

  num_species_cds="$(read_stats_json_field "${merge_stats_file}" "num_species_cds_files")"
  num_species_gff="$(read_stats_json_field "${merge_stats_file}" "num_species_gff_files")"
  num_species_genome="$(read_stats_json_field "${merge_stats_file}" "num_species_genome_files")"
  cds_sequences_before="$(read_stats_json_field "${merge_stats_file}" "cds_sequences_before")"
  cds_sequences_after="$(read_stats_json_field "${merge_stats_file}" "cds_sequences_after")"
  cds_first_sequence_name="$(read_stats_json_field "${merge_stats_file}" "cds_first_sequence_name")"
  merged_rows="$(read_stats_json_field "${merge_stats_file}" "merged_species_summary_rows")"
  task_stats_files="$(read_stats_json_field "${merge_stats_file}" "task_stats_files")"
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
  run_trait_stage
  run_multispecies_summary_stage
  cleanup_input_generation_tmp=1
}

case "${input_generation_mode}" in
  single)
    run_format_stage_single
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
esac

if [[ ${cleanup_input_generation_tmp} -eq 1 && -d "${download_tmp_root}" ]]; then
  echo "Removing temporary input_generation directory: ${download_tmp_root}"
  rm -rf -- "${download_tmp_root}"
fi

echo "$(date): Exiting Singularity environment"

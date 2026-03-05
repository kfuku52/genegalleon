#!/usr/bin/env bash
set -euo pipefail

### Start: Job-supplied configuration ###
# Configuration variables are provided by gg_input_generation_entrypoint.sh.
### End: Job-supplied configuration ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

dir_pg="/workspace"
dir_script="/script/support"
source "${dir_script}/gg_util.sh" # Load utility functions
gg_source_home_bashrc
gg_prepare_cmd_runtime "${dir_pg}" "base" 0 1

config_file="${config_file:-gg_input_generation_entrypoint.sh}"

species_cds_dir="${species_cds_dir:-}"
species_gff_dir="${species_gff_dir:-}"
species_genome_dir="${species_genome_dir:-}"
species_summary_output="${species_summary_output:-}"
resolved_manifest_output="${resolved_manifest_output:-}"

apply_env_override() {
  local var_name=$1
  local env_name=$2
  local env_value="${!env_name:-}"
  if [[ -n "${env_value}" ]]; then
    printf -v "${var_name}" '%s' "${env_value}"
  fi
}

apply_env_override provider GG_INPUT_PROVIDER
apply_env_override strict GG_INPUT_STRICT
apply_env_override overwrite GG_INPUT_OVERWRITE
apply_env_override download_only GG_INPUT_DOWNLOAD_ONLY
apply_env_override dry_run GG_INPUT_DRY_RUN
apply_env_override download_timeout GG_INPUT_DOWNLOAD_TIMEOUT
apply_env_override input_dir GG_INPUT_INPUT_DIR
apply_env_override download_manifest GG_INPUT_DOWNLOAD_MANIFEST
apply_env_override download_dir GG_INPUT_DOWNLOAD_DIR
apply_env_override summary_output GG_INPUT_SUMMARY_OUTPUT
apply_env_override auth_bearer_token_env GG_INPUT_AUTH_BEARER_TOKEN_ENV
apply_env_override http_header GG_INPUT_HTTP_HEADER
apply_env_override run_format_inputs GG_INPUT_RUN_FORMAT_INPUTS
apply_env_override run_validate_inputs GG_INPUT_RUN_VALIDATE_INPUTS
apply_env_override species_cds_dir GG_INPUT_SPECIES_CDS_DIR
apply_env_override species_gff_dir GG_INPUT_SPECIES_GFF_DIR
apply_env_override species_genome_dir GG_INPUT_SPECIES_GENOME_DIR
apply_env_override species_summary_output GG_INPUT_SPECIES_SUMMARY_OUTPUT
apply_env_override resolved_manifest_output GG_INPUT_RESOLVED_MANIFEST_OUTPUT

enable_all_run_flags_for_debug_mode

download_manifest_explicit=0
if [[ -n "${download_manifest}" ]]; then
  download_manifest_explicit=1
fi

case "${provider}" in
  all|ensembl|ensemblplants|phycocosm|phytozome|ncbi|refseq|genbank|coge|cngb|flybase|wormbase|vectorbase|local) ;;
  *)
    echo "Invalid provider: ${provider} (allowed: all|ensembl|ensemblplants|phycocosm|phytozome|ncbi|refseq|genbank|coge|cngb|flybase|wormbase|vectorbase|local)"
    exit 1
    ;;
esac

for binary_flag_name in \
  run_format_inputs \
  run_validate_inputs \
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

if [[ -z "${download_dir}" ]]; then
  download_dir="${dir_pg_output}/input_generation/tmp/input_download_cache"
fi
download_tmp_root="${dir_pg_output}/input_generation/tmp"
if [[ -z "${summary_output}" ]]; then
  summary_output="${dir_pg_output}/input_generation/gg_input_generation_runs.tsv"
fi

if [[ -z "${species_cds_dir}" ]]; then
  species_cds_dir="${dir_pg_output}/input_generation/species_cds"
fi
if [[ -z "${species_gff_dir}" ]]; then
  species_gff_dir="${dir_pg_output}/input_generation/species_gff"
fi
if [[ -z "${species_genome_dir}" ]]; then
  species_genome_dir="${dir_pg_output}/input_generation/species_genome"
fi
if [[ -z "${species_summary_output}" ]]; then
  species_summary_output="${dir_pg_output}/input_generation/gg_input_generation_species.tsv"
fi
if [[ -z "${resolved_manifest_output}" ]]; then
  resolved_manifest_output="${dir_pg_output}/input_generation/download_plan.resolved.tsv"
fi
num_species_cds=""
num_species_gff=""
num_species_genome=""
cds_sequences_before=""
cds_sequences_after=""
cds_first_sequence_name=""
stage_format_status="not_run"
stage_validate_status="not_run"

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
    next(row_iter, None)  # header
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

  run_ended_iso=$(date -u '+%Y-%m-%dT%H:%M:%SZ')
  run_duration_sec=$(( $(date +%s) - run_started_epoch ))

  if [[ ${exit_code} -ne 0 ]]; then
    if [[ "${stage_format_status}" == "running" ]]; then stage_format_status="failed"; fi
    if [[ "${stage_validate_status}" == "running" ]]; then stage_validate_status="failed"; fi
  fi

  ensure_parent_dir "${summary_output}"
  header="started_utc\tended_utc\tduration_sec\texit_code\tprovider\trun_format_inputs\trun_validate_inputs\tstrict\toverwrite\tdownload_only\tdry_run\tdownload_timeout\tinput_dir\tdownload_manifest\tdownload_dir\tspecies_cds_dir\tspecies_gff_dir\tspecies_genome_dir\tnum_species_cds\tnum_species_gff\tnum_species_genome\tcds_sequences_before\tcds_sequences_after\tcds_first_sequence_name\tstage_format_status\tstage_validate_status\tconfig_file"
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
  row="${row}\t$(sanitize_tsv_value "${run_format_inputs}")"
  row="${row}\t$(sanitize_tsv_value "${run_validate_inputs}")"
  row="${row}\t$(sanitize_tsv_value "${strict}")"
  row="${row}\t$(sanitize_tsv_value "${overwrite}")"
  row="${row}\t$(sanitize_tsv_value "${download_only}")"
  row="${row}\t$(sanitize_tsv_value "${dry_run}")"
  row="${row}\t$(sanitize_tsv_value "${download_timeout}")"
  row="${row}\t$(sanitize_tsv_value "${input_dir}")"
  row="${row}\t$(sanitize_tsv_value "${download_manifest}")"
  row="${row}\t$(sanitize_tsv_value "${download_dir}")"
  row="${row}\t$(sanitize_tsv_value "${species_cds_dir}")"
  row="${row}\t$(sanitize_tsv_value "${species_gff_dir}")"
  row="${row}\t$(sanitize_tsv_value "${species_genome_dir}")"
  row="${row}\t$(sanitize_tsv_value "${num_species_cds}")"
  row="${row}\t$(sanitize_tsv_value "${num_species_gff}")"
  row="${row}\t$(sanitize_tsv_value "${num_species_genome}")"
  row="${row}\t$(sanitize_tsv_value "${cds_sequences_before}")"
  row="${row}\t$(sanitize_tsv_value "${cds_sequences_after}")"
  row="${row}\t$(sanitize_tsv_value "${cds_first_sequence_name}")"
  row="${row}\t$(sanitize_tsv_value "${stage_format_status}")"
  row="${row}\t$(sanitize_tsv_value "${stage_validate_status}")"
  row="${row}\t$(sanitize_tsv_value "${config_file}")"
  printf '%b\n' "${row}" >> "${summary_output}"
}

trap write_gg_input_generation_summary_on_exit EXIT

if [[ ${run_format_inputs} -eq 1 ]]; then
  task="Format species inputs"
  gg_step_start "${task}"
  stage_format_status="running"

  if [[ -n "${download_manifest}" ]]; then
    manifest_rows="$(manifest_data_row_count "${download_manifest}")"
    if [[ ${manifest_rows} -eq 0 ]]; then
      echo "Download manifest has no data rows: ${download_manifest}"
      if [[ ${download_manifest_explicit} -eq 1 ]]; then
        echo "Explicitly provided download_manifest is empty."
        stage_format_status="failed"
        exit 1
      fi
      echo "Ignoring empty auto-selected download_manifest."
      download_manifest=""
    fi
  fi

  if [[ -z "${download_manifest}" && -z "${input_dir}" ]]; then
    mapfile -t existing_cds < <(find "${species_cds_dir}" -maxdepth 1 -type f ! -name '.*' \( -name "*.fa" -o -name "*.fa.gz" -o -name "*.fas" -o -name "*.fas.gz" -o -name "*.fasta" -o -name "*.fasta.gz" -o -name "*.fna" -o -name "*.fna.gz" \) | sort)
    mapfile -t existing_gff < <(find "${species_gff_dir}" -maxdepth 1 -type f ! -name '.*' \( -name "*.gff" -o -name "*.gff.gz" -o -name "*.gff3" -o -name "*.gff3.gz" -o -name "*.gtf" -o -name "*.gtf.gz" \) | sort)
    mapfile -t existing_genome < <(find "${species_genome_dir}" -maxdepth 1 -type f ! -name '.*' \( -name "*.fa" -o -name "*.fa.gz" -o -name "*.fas" -o -name "*.fas.gz" -o -name "*.fasta" -o -name "*.fasta.gz" -o -name "*.fna" -o -name "*.fna.gz" \) | sort)
    if [[ ${#existing_cds[@]} -gt 0 || ${#existing_gff[@]} -gt 0 || ${#existing_genome[@]} -gt 0 ]]; then
      echo "Skipping ${task} because no input source was provided and existing species inputs were detected."
      stage_format_status="skipped"
      run_format_inputs=0
    else
      echo "No input source was specified for formatting."
      echo "Set one of input_dir / download_manifest."
      stage_format_status="failed"
      exit 1
    fi
  fi

  if [[ "${stage_format_status}" == "running" ]]; then
    format_stats_file="${download_tmp_root}/gg_input_generation_stats.json"
    ensure_parent_dir "${format_stats_file}"
    rm -f -- "${format_stats_file}"
    cmd=(python "${dir_script}/format_species_inputs.py")
    cmd+=(--provider "${provider}")
    cmd+=(--species-cds-dir "${species_cds_dir}")
    cmd+=(--species-gff-dir "${species_gff_dir}")
    cmd+=(--species-genome-dir "${species_genome_dir}")
    cmd+=(--species-summary-output "${species_summary_output}")
    cmd+=(--stats-output "${format_stats_file}")

    if [[ -n "${download_manifest}" ]]; then
      cmd+=(--download-manifest "${download_manifest}")
      cmd+=(--download-dir "${download_dir}")
      cmd+=(--resolved-manifest-output "${resolved_manifest_output}")
    fi
    if [[ -n "${input_dir}" ]]; then
      cmd+=(--input-dir "${input_dir}")
    fi
    if [[ ${overwrite} -eq 1 ]]; then
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
    cmd+=(--jobs "${NSLOTS:-1}")
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
    stage_format_status="ok"
  fi
else
  gg_step_skip "Format species inputs"
  stage_format_status="skipped"
fi

if [[ ${run_validate_inputs} -eq 1 && ${run_format_inputs} -eq 1 && ${download_only} -eq 0 && ${dry_run} -eq 0 ]]; then
  task="Validate formatted species inputs"
  gg_step_start "${task}"
  stage_validate_status="running"
  mapfile -t cds_files < <(find "${species_cds_dir}" -maxdepth 1 -type f ! -name '.*' \( -name "*.fa" -o -name "*.fa.gz" -o -name "*.fas" -o -name "*.fas.gz" -o -name "*.fasta" -o -name "*.fasta.gz" -o -name "*.fna" -o -name "*.fna.gz" \) | sort)
  mapfile -t gff_files < <(find "${species_gff_dir}" -maxdepth 1 -type f ! -name '.*' \( -name "*.gff" -o -name "*.gff.gz" -o -name "*.gff3" -o -name "*.gff3.gz" -o -name "*.gtf" -o -name "*.gtf.gz" \) | sort)
  mapfile -t genome_files < <(find "${species_genome_dir}" -maxdepth 1 -type f ! -name '.*' \( -name "*.fa" -o -name "*.fa.gz" -o -name "*.fas" -o -name "*.fas.gz" -o -name "*.fasta" -o -name "*.fasta.gz" -o -name "*.fna" -o -name "*.fna.gz" \) | sort)
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
    check_species_cds "${dir_pg}"
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
  fi
  stage_validate_status="ok"
else
  gg_step_skip "Validate formatted species inputs"
  stage_validate_status="skipped"
fi

if [[ -d "${download_tmp_root}" ]]; then
  echo "Removing temporary input_generation directory: ${download_tmp_root}"
  rm -rf -- "${download_tmp_root}"
fi

echo "$(date): Exiting Singularity environment"

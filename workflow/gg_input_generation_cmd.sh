#!/usr/bin/env bash

### Start: Modify this block to tailor your analysis ###

run_build_manifest=1
run_format_inputs=1
run_validate_inputs=1

provider="all" # all|ensemblplants|phycocosm|phytozome
strict=0
overwrite=0
download_only=0
dry_run=0
download_timeout=120

# Optional download request settings.
auth_bearer_token_env="" # e.g., GFE_DOWNLOAD_BEARER_TOKEN
http_header="" # e.g., "User-Agent: genegalleon-inputprep"

# Optional runtime config file. If empty, auto-load from:
#   /workspace/input/manifest/inputprep_profile.sh
# You can also override by exporting GG_INPUTPREP_CONFIG.
config_file=""

# Local raw dataset root for manifest build and direct formatting.
# Example: /workspace/input/raw/gfe_dataset
# If empty and GG_DATASET_ROOT is exported, that value is used.
dataset_root=""

# Optional provider-specific direct input directory. Leave empty in most cases.
# Example: /workspace/input/raw/gfe_dataset/PhycoCosm/species_wise_original
input_dir=""

# Optional download manifest to fetch inputs before formatting.
# Example: /workspace/input/manifest/download_manifest.tsv
download_manifest=""

# Download cache for --download-manifest mode.
download_dir=""

# Output manifest path for run_build_manifest.
manifest_output=""

# Run summary TSV output path.
summary_output=""

### End: Modify this block to tailor your analysis ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

dir_pg="/workspace"
dir_myscript="/script/script"
source "${dir_myscript}/gg_util.sh" # Load utility functions
gg_source_home_bashrc
gg_prepare_cmd_runtime "${dir_pg}" "base" 0 1

if [[ -n "${GG_INPUTPREP_CONFIG}" ]]; then
  config_file="${GG_INPUTPREP_CONFIG}"
fi
if [[ -z "${config_file}" ]]; then
  config_file="${dir_pg_input}/manifest/inputprep_profile.sh"
fi
if [[ -s "${config_file}" ]]; then
  echo "Loading config file: ${config_file}"
  # shellcheck disable=SC1090
  source "${config_file}"
fi

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
apply_env_override dataset_root GG_INPUT_DATASET_ROOT
apply_env_override input_dir GG_INPUT_INPUT_DIR
apply_env_override download_manifest GG_INPUT_DOWNLOAD_MANIFEST
apply_env_override download_dir GG_INPUT_DOWNLOAD_DIR
apply_env_override manifest_output GG_INPUT_MANIFEST_OUTPUT
apply_env_override summary_output GG_INPUT_SUMMARY_OUTPUT
apply_env_override auth_bearer_token_env GG_INPUT_AUTH_BEARER_TOKEN_ENV
apply_env_override http_header GG_INPUT_HTTP_HEADER
apply_env_override run_build_manifest GG_INPUT_RUN_BUILD_MANIFEST
apply_env_override run_format_inputs GG_INPUT_RUN_FORMAT_INPUTS
apply_env_override run_validate_inputs GG_INPUT_RUN_VALIDATE_INPUTS

case "${provider}" in
  all|ensemblplants|phycocosm|phytozome) ;;
  *)
    echo "Invalid provider: ${provider} (allowed: all|ensemblplants|phycocosm|phytozome)"
    exit 1
    ;;
esac

for binary_flag_name in \
  run_build_manifest \
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

if [[ -z "${dataset_root}" ]]; then
  if [[ -n "${GG_DATASET_ROOT}" ]]; then
    dataset_root="${GG_DATASET_ROOT}"
  fi
fi
if [[ -z "${download_dir}" ]]; then
  download_dir="${dir_pg_output}/input_download_cache"
fi
if [[ -z "${manifest_output}" ]]; then
  manifest_output="${dir_pg_input}/manifest/download_manifest.tsv"
fi
if [[ -z "${summary_output}" ]]; then
  summary_output="${dir_pg_output}/inputprep/inputprep_runs.tsv"
fi

species_cds_dir="${dir_pg_input}/species_cds"
species_gff_dir="${dir_pg_input}/species_gff"
num_species_cds=""
num_species_gff=""
stage_manifest_status="not_run"
stage_format_status="not_run"
stage_validate_status="not_run"

run_started_epoch=$(date +%s)
run_started_iso=$(date -u '+%Y-%m-%dT%H:%M:%SZ')

sanitize_tsv_value() {
  local value="$1"
  value=$(printf '%s' "${value}" | tr '\t\r\n' '   ')
  printf '%s' "${value}"
}

write_inputprep_summary_on_exit() {
  local exit_code=$?
  local run_ended_iso
  local run_duration_sec
  local header
  local row

  run_ended_iso=$(date -u '+%Y-%m-%dT%H:%M:%SZ')
  run_duration_sec=$(( $(date +%s) - run_started_epoch ))

  if [[ ${exit_code} -ne 0 ]]; then
    if [[ "${stage_manifest_status}" == "running" ]]; then stage_manifest_status="failed"; fi
    if [[ "${stage_format_status}" == "running" ]]; then stage_format_status="failed"; fi
    if [[ "${stage_validate_status}" == "running" ]]; then stage_validate_status="failed"; fi
  fi

  ensure_parent_dir "${summary_output}"
  header="started_utc\tended_utc\tduration_sec\texit_code\tprovider\trun_build_manifest\trun_format_inputs\trun_validate_inputs\tstrict\toverwrite\tdownload_only\tdry_run\tdownload_timeout\tdataset_root\tinput_dir\tdownload_manifest\tdownload_dir\tmanifest_output\tspecies_cds_dir\tspecies_gff_dir\tnum_species_cds\tnum_species_gff\tstage_manifest_status\tstage_format_status\tstage_validate_status\tconfig_file"
  if [[ ! -s "${summary_output}" ]]; then
    printf '%b\n' "${header}" > "${summary_output}"
  fi
  row="$(sanitize_tsv_value "${run_started_iso}")"
  row="${row}\t$(sanitize_tsv_value "${run_ended_iso}")"
  row="${row}\t$(sanitize_tsv_value "${run_duration_sec}")"
  row="${row}\t$(sanitize_tsv_value "${exit_code}")"
  row="${row}\t$(sanitize_tsv_value "${provider}")"
  row="${row}\t$(sanitize_tsv_value "${run_build_manifest}")"
  row="${row}\t$(sanitize_tsv_value "${run_format_inputs}")"
  row="${row}\t$(sanitize_tsv_value "${run_validate_inputs}")"
  row="${row}\t$(sanitize_tsv_value "${strict}")"
  row="${row}\t$(sanitize_tsv_value "${overwrite}")"
  row="${row}\t$(sanitize_tsv_value "${download_only}")"
  row="${row}\t$(sanitize_tsv_value "${dry_run}")"
  row="${row}\t$(sanitize_tsv_value "${download_timeout}")"
  row="${row}\t$(sanitize_tsv_value "${dataset_root}")"
  row="${row}\t$(sanitize_tsv_value "${input_dir}")"
  row="${row}\t$(sanitize_tsv_value "${download_manifest}")"
  row="${row}\t$(sanitize_tsv_value "${download_dir}")"
  row="${row}\t$(sanitize_tsv_value "${manifest_output}")"
  row="${row}\t$(sanitize_tsv_value "${species_cds_dir}")"
  row="${row}\t$(sanitize_tsv_value "${species_gff_dir}")"
  row="${row}\t$(sanitize_tsv_value "${num_species_cds}")"
  row="${row}\t$(sanitize_tsv_value "${num_species_gff}")"
  row="${row}\t$(sanitize_tsv_value "${stage_manifest_status}")"
  row="${row}\t$(sanitize_tsv_value "${stage_format_status}")"
  row="${row}\t$(sanitize_tsv_value "${stage_validate_status}")"
  row="${row}\t$(sanitize_tsv_value "${config_file}")"
  printf '%b\n' "${row}" >> "${summary_output}"
}

trap write_inputprep_summary_on_exit EXIT

if [[ ${run_build_manifest} -eq 1 ]]; then
  task="Build input download manifest"
  gg_step_start "${task}"
  stage_manifest_status="running"
  if [[ -z "${dataset_root}" ]]; then
    dataset_root="${dir_pg_input}/raw/gfe_dataset"
  fi
  if [[ ! -d "${dataset_root}" ]]; then
    echo "Skipping ${task}. dataset_root was not found: ${dataset_root}"
    stage_manifest_status="skipped"
  else
    ensure_parent_dir "${manifest_output}"
    cmd=(python "${dir_myscript}/build_download_manifest.py")
    cmd+=(--provider "${provider}")
    cmd+=(--dataset-root "${dataset_root}")
    cmd+=(--output "${manifest_output}")
    if [[ ${strict} -eq 1 ]]; then
      cmd+=(--strict)
    fi
    echo "Running: ${cmd[*]}"
    "${cmd[@]}"; cmd_status=$?
    if [[ ${cmd_status} -ne 0 ]]; then
      stage_manifest_status="failed"
      echo "Failed: ${task} (exit=${cmd_status})"
      exit "${cmd_status}"
    fi
    stage_manifest_status="ok"
  fi
else
  gg_step_skip "Build input download manifest"
  stage_manifest_status="skipped"
fi

if [[ -z "${download_manifest}" && -s "${manifest_output}" && ${run_build_manifest} -eq 1 ]]; then
  download_manifest="${manifest_output}"
  echo "Using manifest generated in this run: ${download_manifest}"
fi

if [[ ${run_format_inputs} -eq 1 ]]; then
  task="Format species inputs"
  gg_step_start "${task}"
  stage_format_status="running"
  if [[ -z "${download_manifest}" && -z "${input_dir}" && -z "${dataset_root}" ]]; then
    echo "No input source was specified for formatting."
    echo "Set one of dataset_root / input_dir / download_manifest (or GG_DATASET_ROOT)."
    stage_format_status="failed"
    exit 1
  fi

  cmd=(python "${dir_myscript}/format_species_inputs.py")
  cmd+=(--provider "${provider}")
  cmd+=(--species-cds-dir "${species_cds_dir}")
  cmd+=(--species-gff-dir "${species_gff_dir}")

  if [[ -n "${download_manifest}" ]]; then
    cmd+=(--download-manifest "${download_manifest}")
    cmd+=(--download-dir "${download_dir}")
  fi
  if [[ -n "${input_dir}" ]]; then
    cmd+=(--input-dir "${input_dir}")
  elif [[ -n "${dataset_root}" && -z "${download_manifest}" ]]; then
    cmd+=(--dataset-root "${dataset_root}")
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
  if [[ -n "${auth_bearer_token_env}" ]]; then
    cmd+=(--auth-bearer-token-env "${auth_bearer_token_env}")
  fi
  if [[ -n "${http_header}" ]]; then
    cmd+=(--http-header "${http_header}")
  fi

  echo "Running: ${cmd[*]}"
  "${cmd[@]}"; cmd_status=$?
  if [[ ${cmd_status} -ne 0 ]]; then
    stage_format_status="failed"
    echo "Failed: ${task} (exit=${cmd_status})"
    exit "${cmd_status}"
  fi
  stage_format_status="ok"
else
  gg_step_skip "Format species inputs"
  stage_format_status="skipped"
fi

if [[ ${run_validate_inputs} -eq 1 && ${run_format_inputs} -eq 1 && ${download_only} -eq 0 && ${dry_run} -eq 0 ]]; then
  task="Validate formatted species inputs"
  gg_step_start "${task}"
  stage_validate_status="running"
  mapfile -t cds_files < <(find "${species_cds_dir}" -maxdepth 1 -type f | grep $(get_fasta_extensions_for_grep))
  mapfile -t gff_files < <(find "${species_gff_dir}" -maxdepth 1 -type f | grep -E '\.gff3?(\.gz)?$|\.gtf(\.gz)?$')
  num_species_cds=${#cds_files[@]}
  num_species_gff=${#gff_files[@]}
  echo "Detected ${#cds_files[@]} species CDS files and ${#gff_files[@]} species GFF files."

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

  if [[ ${#cds_files[@]} -gt 0 && ${#gff_files[@]} -gt 0 ]]; then
    is_species_set_identical "${species_cds_dir}" "${species_gff_dir}"
    set_status=$?
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

echo "$(date): Exiting Singularity environment"

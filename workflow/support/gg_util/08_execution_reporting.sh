# shellcheck shell=bash
# Sourced by workflow/support/gg_util.sh.

gg_trigger_versions_dump() {
  local trigger_name=${1:-unknown_job}
  local versions_script
  local version_file
  local gg_version
  local dir_output
  local versions_dir
  local container_key_seed
  local container_key_hash
  local inspect_snapshot
  local image_file_hash
  local versions_script_hash=""
  local log_file
  local tmp_log_file
  local lock_file
  local failed_log_file
  local had_flock=0
  local singularity_bin
  local container_runtime_bin
  local versions_exit_code=0
  local block_exit_code=0
  local had_errexit=0

  if [[ -z "${gg_workflow_dir:-}" || -z "${gg_support_dir:-}" || -z "${gg_workspace_dir:-}" || -z "${gg_container_image_path:-}" ]]; then
    echo "gg_trigger_versions_dump: gg_workflow_dir/gg_support_dir/gg_workspace_dir/gg_container_image_path are required." >&2
    return 1
  fi

  version_file="${gg_workflow_dir}/../VERSION"
  gg_version="${SINGULARITYENV_GG_VERSION:-}"
  if [[ -z "${gg_version}" ]]; then
    gg_version="$(gg_read_repo_version "${version_file}")"
  fi
  export SINGULARITYENV_GG_VERSION="${gg_version}"
  export APPTAINERENV_GG_VERSION="${gg_version}"

  versions_script="${gg_support_dir}/gg_versions.sh"
  if [[ ! -s "${versions_script}" ]]; then
    echo "gg_trigger_versions_dump: versions script not found: ${versions_script}" >&2
    return 1
  fi

  if ! gg_container_shell_command_is_set; then
    set_singularity_command
  fi
  container_runtime_bin="$(gg_container_shell_command_runtime_bin || true)"
  if [[ -z "${container_runtime_bin}" ]]; then
    if ! container_runtime_bin=$(gg_detect_container_runtime_binary); then
      echo "gg_trigger_versions_dump: container runtime command not found. Skipping version dump." >&2
      return 0
    fi
  fi
  if ! command -v "${container_runtime_bin}" >/dev/null 2>&1; then
    echo "gg_trigger_versions_dump: container runtime command not found (${container_runtime_bin}). Skipping version dump." >&2
    return 0
  fi

  if ! gg_container_bind_destination_exists "${gg_workspace_dir}:/workspace" \
    || ! gg_container_bind_destination_exists "${gg_workflow_dir}:/script"; then
    gg_normalize_scheduler_env
    set_singularityenv
  fi

  dir_output=$(workspace_output_root "${gg_workspace_dir}")
  versions_dir="${dir_output}/versions"
  ensure_dir "${versions_dir}"

  container_key_seed="gg_container_image_path=${gg_container_image_path};runtime=${container_runtime_bin};gg_version=${gg_version}"
  if [[ -s "${versions_script}" ]]; then
    if command -v sha256sum >/dev/null 2>&1; then
      versions_script_hash=$(sha256sum "${versions_script}" | awk '{print $1}')
      container_key_seed="${container_key_seed};versions_script_sha256=${versions_script_hash}"
    elif command -v shasum >/dev/null 2>&1; then
      versions_script_hash=$(shasum -a 256 "${versions_script}" | awk '{print $1}')
      container_key_seed="${container_key_seed};versions_script_sha256=${versions_script_hash}"
    else
      versions_script_hash=$(cksum "${versions_script}" | awk '{print $1 "-" $2}')
      container_key_seed="${container_key_seed};versions_script_cksum=${versions_script_hash}"
    fi
  fi
  if [[ -s "${gg_container_image_path}" ]]; then
    if command -v sha256sum >/dev/null 2>&1; then
      image_file_hash=$(sha256sum "${gg_container_image_path}" | awk '{print $1}')
      container_key_seed="${container_key_seed};image_sha256=${image_file_hash}"
    elif command -v shasum >/dev/null 2>&1; then
      image_file_hash=$(shasum -a 256 "${gg_container_image_path}" | awk '{print $1}')
      container_key_seed="${container_key_seed};image_sha256=${image_file_hash}"
    else
      image_file_hash=$(cksum "${gg_container_image_path}" | awk '{print $1 "-" $2}')
      container_key_seed="${container_key_seed};image_cksum=${image_file_hash}"
    fi
  fi
  inspect_snapshot=""
  if inspect_snapshot=$("${container_runtime_bin}" inspect "${gg_container_image_path}" 2>/dev/null); then
    if [[ -n "${inspect_snapshot}" ]]; then
      container_key_seed="${container_key_seed};inspect=${inspect_snapshot}"
    fi
  fi
  if command -v sha256sum >/dev/null 2>&1; then
    container_key_hash=$(printf '%s' "${container_key_seed}" | sha256sum | awk '{print $1}')
  elif command -v shasum >/dev/null 2>&1; then
    container_key_hash=$(printf '%s' "${container_key_seed}" | shasum -a 256 | awk '{print $1}')
  else
    container_key_hash=$(printf '%s' "${container_key_seed}" | cksum | awk '{print $1}')
  fi
  if [[ -z "${container_key_hash}" ]]; then
    container_key_hash=$(echo "${gg_container_image_path}" | tr '[:space:]/:' '_' | tr -cd '[:alnum:]_.-')
    if [[ -z "${container_key_hash}" ]]; then
      container_key_hash="unknown_container"
    fi
  fi
  log_file="${versions_dir}/container.${container_key_hash}.versions.log"
  tmp_log_file="${log_file}.tmp.$$"
  lock_file="${log_file}.lock"

  if command -v flock >/dev/null 2>&1; then
    exec 9>"${lock_file}"
    flock 9
    had_flock=1
  fi
  if [[ -s "${log_file}" ]]; then
    echo "gg_trigger_versions_dump: skipped existing ${log_file}"
    if [[ ${had_flock} -eq 1 ]]; then
      flock -u 9
      exec 9>&-
    fi
    return 0
  fi

  if [[ $- == *e* ]]; then
    had_errexit=1
    set +e
  fi
  {
    gg_print_version_summary "genegalleon versions" "${container_runtime_bin}" "${inspect_snapshot}"
    echo "$(date): Triggered gg_versions by ${trigger_name}"
    gg_run_container_shell_script "${gg_container_image_path}" "${versions_script}" || {
      cmd_rc=$?
      if [[ ${versions_exit_code} -eq 0 ]]; then
        versions_exit_code=${cmd_rc}
      fi
    }
    echo ""
    echo "### gg_container ###"
    singularity_bin="$(command -v "${container_runtime_bin}")"
    if [[ "${singularity_bin}" == *"/gg_wrapper_bin/"* ]]; then
      echo "Skipping container inspect/version under Docker-backed singularity shim: ${singularity_bin}"
    else
      "${container_runtime_bin}" inspect "${gg_container_image_path}" || {
        cmd_rc=$?
        if [[ ${versions_exit_code} -eq 0 ]]; then
          versions_exit_code=${cmd_rc}
        fi
      }
      echo ""
      echo "### ${container_runtime_bin} version ###"
      "${container_runtime_bin}" version || {
        cmd_rc=$?
        if [[ ${versions_exit_code} -eq 0 ]]; then
          versions_exit_code=${cmd_rc}
        fi
      }
    fi
    echo ""
    echo "### Host OS info ###"
    if [[ -f /etc/os-release ]]; then
      cat /etc/os-release
    else
      echo "/etc/os-release was not found. Falling back to uname output."
      uname -a
    fi
    echo ""
    echo "$(date): gg_versions trigger completed"
  } > "${tmp_log_file}" 2>&1
  # The status belongs to the redirected command group above, not to a
  # standalone echo/printf command.
  # shellcheck disable=SC2320
  block_exit_code=$?
  if [[ ${block_exit_code} -ne 0 && ${versions_exit_code} -eq 0 ]]; then
    versions_exit_code=${block_exit_code}
  fi
  if [[ ${had_errexit} -eq 1 ]]; then
    set -e
  fi

  if [[ ${versions_exit_code} -ne 0 ]]; then
    failed_log_file="${versions_dir}/container.${container_key_hash}.versions.failed.$(date '+%Y%m%d_%H%M%S').log"
    if [[ -s "${tmp_log_file}" ]]; then
      mv_out "${tmp_log_file}" "${failed_log_file}"
    else
      : > "${failed_log_file}"
    fi
    if [[ ${had_flock} -eq 1 ]]; then
      flock -u 9
      exec 9>&-
    fi
    echo "gg_trigger_versions_dump: failed (exit=${versions_exit_code}). Log: ${failed_log_file}" >&2
    return "${versions_exit_code}"
  fi

  if [[ -s "${tmp_log_file}" ]]; then
    mv_out "${tmp_log_file}" "${log_file}"
  fi
  if [[ ${had_flock} -eq 1 ]]; then
    flock -u 9
    exec 9>&-
  fi
  echo "gg_trigger_versions_dump: wrote ${log_file}"
  return 0
}

gg_require_versions_dump() {
  local trigger_name=${1:-unknown_job}
  local versions_exit_code=0
  local had_errexit=0

  if [[ $- == *e* ]]; then
    had_errexit=1
    set +e
  fi
  gg_trigger_versions_dump "${trigger_name}"
  versions_exit_code=$?
  if [[ ${had_errexit} -eq 1 ]]; then
    set -e
  fi
  if [[ ${versions_exit_code} -ne 0 ]]; then
    echo "gg_require_versions_dump: gg_versions trigger failed for ${trigger_name}." >&2
    return "${versions_exit_code}"
  fi
  return 0
}

enable_all_run_flags_for_debug_mode() {
  local message=${1:-"gg debug mode: All run_* variables are forced to set 1."}
  local debug_mode=${gg_debug_mode:-0}
  local flags=()
  local flag

  if [[ "${debug_mode}" -ne 1 ]]; then
    return 0
  fi

  echo "${message}"
  while IFS= read -r flag; do
    flags+=( "${flag}" )
  done < <(set | sed -n -E 's/^(run_[A-Za-z0-9_]+)=.*/\1/p')
  for flag in "${flags[@]}"; do
    printf -v "${flag}" '%s' 1
    export "${flag?}"
  done
  return 0
}

ensure_gg_scheduler_defaults() {
  local echo_header="ensure_gg_scheduler_defaults: "
  if [[ -z "${GG_TASK_CPUS:-}" ]]; then
    echo "${echo_header}GG_TASK_CPUS is undefined or empty. GG_TASK_CPUS=1"
    GG_TASK_CPUS=1
  fi
  if [[ -z "${GG_JOB_ID:-}" ]]; then
    echo "${echo_header}GG_JOB_ID is undefined or empty. GG_JOB_ID=1"
    GG_JOB_ID=1
  fi
  if [[ -z "${GG_ARRAY_TASK_ID:-}" ]]; then
    echo "${echo_header}GG_ARRAY_TASK_ID is undefined or empty. GG_ARRAY_TASK_ID=1"
    GG_ARRAY_TASK_ID=1
  fi
  if [[ -z "${GG_MEM_PER_CPU_GB:-}" ]]; then
    echo "${echo_header}GG_MEM_PER_CPU_GB is undefined or empty. GG_MEM_PER_CPU_GB=3"
    GG_MEM_PER_CPU_GB=3
  fi
  if [[ -z "${GG_MEM_TOTAL_GB:-}" ]]; then
    GG_MEM_TOTAL_GB=$((GG_MEM_PER_CPU_GB * GG_TASK_CPUS))
    echo "${echo_header}GG_MEM_TOTAL_GB is undefined or empty. GG_MEM_TOTAL_GB=${GG_MEM_TOTAL_GB}"
  fi
  gg_normalize_memory_budget
  gg_sync_legacy_scheduler_aliases
}

ensure_scheduler_defaults() {
  ensure_gg_scheduler_defaults "$@"
}

print_gg_container_starting_message() {
  local gg_workspace_input_dir_resolved
  local gg_workspace_output_dir_resolved
  local gg_workspace_downloads_dir_resolved
  local gg_workspace_layout_local
  ensure_gg_scheduler_defaults
  gg_workspace_layout_local=$(gg_resolve_workspace_layout "${gg_workspace_dir}")
  gg_workspace_input_dir_resolved=$(workspace_input_root "${gg_workspace_dir}")
  gg_workspace_output_dir_resolved=$(workspace_output_root "${gg_workspace_dir}")
  gg_workspace_downloads_dir_resolved=$(workspace_downloads_root "${gg_workspace_dir}")
  echo "$(date): Starting gg Singularity/Apptainer environment"
  echo "pwd: $(pwd)"
  echo "gg_workspace_dir: ${gg_workspace_dir}"
  echo "gg_workspace_layout: ${gg_workspace_layout_local}"
  echo "gg_workspace_input_dir: ${gg_workspace_input_dir_resolved}"
  echo "gg_workspace_output_dir: ${gg_workspace_output_dir_resolved}"
  echo "gg_workspace_downloads_dir: ${gg_workspace_downloads_dir_resolved}"
  echo "python: $(command -v python)"
  echo "GG_TASK_CPUS: ${GG_TASK_CPUS}"
  echo "GG_MEM_PER_CPU_GB: ${GG_MEM_PER_CPU_GB}"
  echo "GG_MEM_TOTAL_GB: ${GG_MEM_TOTAL_GB}"
  echo "GG_MEM_TOOL_RESERVE_GB: ${GG_MEM_TOOL_RESERVE_GB}"
  echo "GG_MEM_TOOL_GB: ${GG_MEM_TOOL_GB}"
  echo "GG_JOB_ID: ${GG_JOB_ID}"
  echo "GG_ARRAY_TASK_ID: ${GG_ARRAY_TASK_ID}"
  echo "ulimit -Hn: $(ulimit -Hn)"
  echo "ulimit -Sn: $(ulimit -Sn)"
}

recreate_dir() {
  local dir=$1
  if [[ -z "${dir}" || "${dir}" == "/" ]]; then
    echo "Refusing to recreate unsafe directory path: ${dir}"
    return 1
  fi
  if [[ -e "${dir}" ]]; then
    echo "Deleting: ${dir}"
    rm -rf -- "${dir}"
  fi
  echo "Creating: ${dir}"
  mkdir -p "${dir}"
}

gg_task_token_contains_sge_id() {
  local task_token=$1
  local target_task_id=$2
  local token_parts=()
  IFS=',' read -r -a token_parts <<< "${task_token}"
  local part
  for part in "${token_parts[@]}"; do
    if [[ "${part}" == "${target_task_id}" ]]; then
      return 0
    fi
    if [[ "${part}" =~ ^([0-9]+)-([0-9]+)(:([0-9]+))?$ ]]; then
      local start_id=${BASH_REMATCH[1]}
      local end_id=${BASH_REMATCH[2]}
      local step_id=${BASH_REMATCH[4]:-1}
      if (( target_task_id >= start_id && target_task_id <= end_id )); then
        if (( (target_task_id - start_id) % step_id == 0 )); then
          return 0
        fi
      fi
    fi
  done
  return 1
}

gg_is_sge_task_present_in_qstat() {
  local job_id=$1
  local target_task_id=$2
  if ! command -v qstat >/dev/null 2>&1; then
    return 2
  fi
  if [[ -z "${job_id}" ]]; then
    return 2
  fi
  local qstat_task_tokens=()
  local qstat_token
  while IFS= read -r qstat_token; do
    qstat_task_tokens+=( "${qstat_token}" )
  done < <(qstat 2>/dev/null | awk -v jid="${job_id}" 'NR>2 && $1==jid {print $NF}')
  if [[ ${#qstat_task_tokens[@]} -eq 0 ]]; then
    return 1
  fi
  local token
  for token in "${qstat_task_tokens[@]}"; do
    if gg_task_token_contains_sge_id "${token}" "${target_task_id}"; then
      return 0
    fi
  done
  return 1
}

gg_get_slurm_array_master_job_id() {
  if [[ -n "${SLURM_ARRAY_JOB_ID:-}" ]]; then
    echo "${SLURM_ARRAY_JOB_ID}"
    return 0
  fi
  local id_candidate="${SLURM_JOB_ID:-${GG_JOB_ID:-}}"
  if [[ -z "${id_candidate}" ]]; then
    return 0
  fi
  if [[ "${id_candidate}" =~ ^([0-9]+)_([0-9]+)$ ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    echo "${id_candidate}"
  fi
}

gg_task_token_contains_slurm_id() {
  local job_token=$1
  local array_job_id=$2
  local target_task_id=$3
  local suffix
  local expr
  local expr_parts=()
  local part

  if [[ "${job_token}" == "${array_job_id}_${target_task_id}" ]]; then
    return 0
  fi
  if [[ "${job_token}" != "${array_job_id}_"* ]]; then
    return 1
  fi

  suffix="${job_token#${array_job_id}_}"
  if [[ "${suffix}" == \[*\] ]]; then
    expr="${suffix#[}"
    expr="${expr%]}"
  else
    expr="${suffix}"
  fi

  IFS=',' read -r -a expr_parts <<< "${expr}"
  for part in "${expr_parts[@]}"; do
    part="${part%%%*}" # remove SLURM throttling suffix, e.g. 1-100%10
    if [[ "${part}" == "${target_task_id}" ]]; then
      return 0
    fi
    if [[ "${part}" =~ ^([0-9]+)-([0-9]+)(:([0-9]+))?$ ]]; then
      local start_id=${BASH_REMATCH[1]}
      local end_id=${BASH_REMATCH[2]}
      local step_id=${BASH_REMATCH[4]:-1}
      if (( target_task_id >= start_id && target_task_id <= end_id )); then
        if (( (target_task_id - start_id) % step_id == 0 )); then
          return 0
        fi
      fi
    fi
  done
  return 1
}

gg_is_sge_task_present_in_squeue() {
  local array_job_id=$1
  local target_task_id=$2
  if ! command -v squeue >/dev/null 2>&1; then
    return 2
  fi
  if [[ -z "${array_job_id}" ]]; then
    return 2
  fi

  # Fast path: query the specific array task directly.
  if [[ -n "$(squeue -h -j "${array_job_id}_${target_task_id}" -o "%i" 2>/dev/null)" ]]; then
    return 0
  fi

  local slurm_job_tokens=()
  local slurm_token
  while IFS= read -r slurm_token; do
    slurm_job_tokens+=( "${slurm_token}" )
  done < <(squeue -h -j "${array_job_id}" -o "%i" 2>/dev/null)
  if [[ ${#slurm_job_tokens[@]} -eq 0 ]]; then
    return 1
  fi
  local token
  for token in "${slurm_job_tokens[@]}"; do
    if gg_task_token_contains_slurm_id "${token}" "${array_job_id}" "${target_task_id}"; then
      return 0
    fi
  done
  return 1
}

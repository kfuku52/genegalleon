#!/usr/bin/env bash
set -euo pipefail

gg_resolve_core_support_dir() {
  local core_script_path=${1:-${BASH_SOURCE[1]:-}}
  local core_dir=""
  local candidate=""

  if [[ -n "${core_script_path}" ]]; then
    core_dir="$(cd "$(dirname "${core_script_path}")" && pwd 2>/dev/null || true)"
  fi

  for candidate in \
    "${gg_support_dir:-}" \
    "/script/support" \
    "${core_dir}/../support"
  do
    [[ -n "${candidate}" ]] || continue
    if [[ -s "${candidate}/gg_util.sh" ]]; then
      (cd "${candidate}" && pwd -P)
      return 0
    fi
  done

  return 1
}

gg_source_common_params_from_core() {
  local core_script_path=${1:-${BASH_SOURCE[1]:-}}
  local core_dir=""
  local candidate=""

  if [[ -n "${core_script_path}" ]]; then
    core_dir="$(cd "$(dirname "${core_script_path}")" && pwd 2>/dev/null || true)"
  fi

  for candidate in \
    "${core_dir}/gg_common_params.sh" \
    "${core_dir}/../gg_common_params.sh" \
    "/script/gg_common_params.sh"
  do
    if [[ -s "${candidate}" ]]; then
      # shellcheck disable=SC1090
      source "${candidate}"
      GG_CORE_COMMON_PARAMS_LOADED=1
      return 0
    fi
  done

  return 0
}

gg_configure_python_pycacheprefix_from_core() {
  [[ -z "${PYTHONPYCACHEPREFIX:-}" ]] || return 0
  local storage_support="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
  # shellcheck disable=SC1090
  source "${storage_support}/gg_tmp_storage.sh"
  local selected runtime_tmp
  selected=$(gg_resolve_tmp_root) || return 1
  runtime_tmp=$(gg_private_runtime_tmp "${selected}") || return 1
  export PYTHONPYCACHEPREFIX="${runtime_tmp}/pycache"
}

gg_bootstrap_core_runtime() {
  local core_script_path=${1:-${BASH_SOURCE[1]:-}}
  local conda_env=${2:-}
  local set_unlimited_stack=${3:-1}
  local print_start_message=${4:-1}
  local resolved_support_dir=""

  # Apptainer/Singularity normally bind the host home directory. Disable the
  # Python user site again inside the core process so host ~/.local packages
  # cannot shadow the versions installed in the GeneGalleon image.
  export PYTHONNOUSERSITE=1
  export GG_RESOURCE_OWNER_PID=${BASHPID:-$$}


  if [[ "${GG_CORE_COMMON_PARAMS_LOADED:-0}" -ne 1 ]]; then
    gg_source_common_params_from_core "${core_script_path}"
  fi

  resolved_support_dir="$(gg_resolve_core_support_dir "${core_script_path}")" || {
    echo "Failed to resolve gg_util.sh from core_script_path=${core_script_path:-unset}" >&2
    return 1
  }

  gg_support_dir="${resolved_support_dir}"
  gg_workflow_dir="$(cd "${gg_support_dir}/.." && pwd -P)"
  gg_core_dir="${gg_workflow_dir}/core"
  : "${gg_workspace_dir:=/workspace}"
  export GG_CONTENT_DIGEST_CACHE="${GG_CONTENT_DIGEST_CACHE:-${gg_workspace_dir}/.gg_cache/content_digests.sqlite3}"
  gg_configure_python_pycacheprefix_from_core

  # shellcheck disable=SC1090
  source "${gg_support_dir}/gg_util.sh"
  if [[ "${GG_COMMON_TMP_ROOT:-auto}" == auto ]]; then
    GG_COMMON_TMP_ROOT=$(gg_resolve_tmp_root) || return 1
    export GG_COMMON_TMP_ROOT
  fi
  if [[ "${GG_COMMON_TMP_ROOT}" == workspace && -z "${GG_TMP_TASK_ROOT:-}" ]]; then
    local runtime_tmp
    runtime_tmp=$(gg_private_runtime_tmp workspace) || return 1
    export TMPDIR="${runtime_tmp}/tmp" TMP="${runtime_tmp}/tmp" TEMP="${runtime_tmp}/tmp"
  fi
  gg_prepare_cmd_runtime "${gg_workspace_dir}" "${conda_env}" "${set_unlimited_stack}" "${print_start_message}"
}

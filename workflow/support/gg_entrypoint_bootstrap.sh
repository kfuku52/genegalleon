#!/usr/bin/env bash
set -euo pipefail

gg_search_workflow_dir_from_base() {
  local base_dir=${1:-}
  local search_dir=""
  local parent_dir=""

  [[ -n "${base_dir}" ]] || return 1
  search_dir="$(cd "$(dirname "${base_dir}")" 2>/dev/null && pwd -P || true)"
  if [[ -d "${base_dir}" ]]; then
    search_dir="$(cd "${base_dir}" 2>/dev/null && pwd -P || true)"
  fi

  while [[ -n "${search_dir}" ]]; do
    if [[ -f "${search_dir}/support/gg_util.sh" ]]; then
      (cd "${search_dir}" && pwd -P)
      return 0
    fi
    if [[ -f "${search_dir}/workflow/support/gg_util.sh" ]]; then
      (cd "${search_dir}/workflow" && pwd -P)
      return 0
    fi
    parent_dir="$(dirname "${search_dir}")"
    if [[ "${parent_dir}" == "${search_dir}" ]]; then
      break
    fi
    search_dir="${parent_dir}"
  done

  return 1
}

gg_find_workflow_dir() {
  local entrypoint_path=${1:-}
  local bootstrap_path=${2:-${BASH_SOURCE[0]:-}}
  local base_dir=""

  for base_dir in "${entrypoint_path}" "${bootstrap_path}"; do
    if gg_search_workflow_dir_from_base "${base_dir}"; then
      return 0
    fi
  done

  for base_dir in "${SLURM_SUBMIT_DIR:-}" "${PBS_O_WORKDIR:-}" "${PWD:-}"; do
    [[ -n "${base_dir}" ]] || continue
    if gg_search_workflow_dir_from_base "${base_dir}"; then
      return 0
    fi
  done

  return 1
}

gg_set_workflow_dir() {
  local entrypoint_path=${1:-${BASH_SOURCE[1]:-}}
  local bootstrap_path=${2:-${BASH_SOURCE[0]:-}}
  local resolved_workflow_dir=""

  resolved_workflow_dir="$(gg_find_workflow_dir "${entrypoint_path}" "${bootstrap_path}")" || {
    echo "Failed to resolve workflow directory from entrypoint_path=${entrypoint_path:-unset}" >&2
    return 1
  }

  gg_workflow_dir="${resolved_workflow_dir}"
}

gg_source_common_params_if_available() {
  local resolved_workflow_dir=${1:-${gg_workflow_dir:-}}

  if [[ -n "${resolved_workflow_dir}" && -s "${resolved_workflow_dir}/gg_common_params.sh" ]]; then
    # shellcheck disable=SC1090
    source "${resolved_workflow_dir}/gg_common_params.sh"
  fi
}

gg_configure_python_pycacheprefix() {
  local default_pycache_prefix=""

  if [[ -n "${PYTHONPYCACHEPREFIX:-}" ]]; then
    return 0
  fi

  default_pycache_prefix="${TMPDIR:-/tmp}/genegalleon_pycache"
  mkdir -p -- "${default_pycache_prefix}" 2>/dev/null || true
  export PYTHONPYCACHEPREFIX="${default_pycache_prefix}"
}

gg_entrypoint_initialize() {
  local entrypoint_path=${1:-${BASH_SOURCE[1]:-}}
  local load_common_params=${2:-0}

  if ! gg_set_workflow_dir "${entrypoint_path}"; then
    return 1
  fi
  # shellcheck disable=SC1090
  source "${gg_workflow_dir}/gg_path_defaults.sh"
  if [[ "${load_common_params}" -eq 1 ]]; then
    gg_source_common_params_if_available "${gg_workflow_dir}"
  fi
  gg_configure_python_pycacheprefix
}

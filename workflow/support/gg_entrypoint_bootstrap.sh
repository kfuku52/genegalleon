#!/usr/bin/env bash
set -euo pipefail

gg_find_workflow_dir() {
  local entrypoint_path=${1:-}
  local search_dir=""
  local parent_dir=""

  if [[ -n "${entrypoint_path}" ]]; then
    search_dir="$(cd "$(dirname "${entrypoint_path}")" && pwd 2>/dev/null || true)"
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
  fi

  for base_dir in "${SLURM_SUBMIT_DIR:-}" "${PBS_O_WORKDIR:-}" "${PWD:-}"; do
    [[ -n "${base_dir}" ]] || continue
    search_dir="${base_dir}"
    while true; do
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
  done

  return 1
}

gg_set_workflow_dir() {
  local entrypoint_path=${1:-${BASH_SOURCE[1]:-}}
  local resolved_workflow_dir=""

  resolved_workflow_dir="$(gg_find_workflow_dir "${entrypoint_path}")" || {
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

gg_entrypoint_print_config_summary_if_available() {
  local entrypoint_name=${1:-entrypoint}
  shift || true

  if declare -F gg_print_entrypoint_config_summary >/dev/null 2>&1; then
    gg_print_entrypoint_config_summary "${entrypoint_name}" "$@"
  else
    echo "Config summary helper is unavailable in gg_util.sh; skipping entrypoint config summary for ${entrypoint_name}."
  fi
}

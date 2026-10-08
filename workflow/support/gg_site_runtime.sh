#!/usr/bin/env bash
set -euo pipefail

gg_detect_site_profile() {
  local hostname_text

  if [[ -n "${GG_SITE_PROFILE:-}" ]]; then
    echo "${GG_SITE_PROFILE}"
    return 0
  fi

  if [[ "${SGE_ROOT:-}" == "/home/geadmin/N1GE" ]]; then
    echo "shirokane"
    return 0
  fi

  hostname_text="$(hostname 2>/dev/null || true)"
  if [[ "${hostname_text}" =~ ^a[0-9]+$ ]]; then
    echo "nig"
    return 0
  fi
  case "${hostname_text}" in
    audrey1)
      echo "audrey1"
      return 0
      ;;
    at*|m*|igt*|it*)
      echo "nig"
      return 0
      ;;
    *.nhr.fau.de)
      echo "nhr-fau"
      return 0
      ;;
    *)
      echo "default"
      return 0
      ;;
  esac
}

gg_initialize_environment_modules() {
  local module_init

  if type module >/dev/null 2>&1; then
    return 0
  fi

  for module_init in \
    /etc/profile.d/modules.sh \
    /usr/share/Modules/init/bash \
    /usr/share/lmod/lmod/init/bash
  do
    if [[ -s "${module_init}" ]]; then
      # shellcheck disable=SC1090
      source "${module_init}"
      if type module >/dev/null 2>&1; then
        return 0
      fi
    fi
  done

  echo "Environment Modules could not be initialized." >&2
  return 1
}

gg_shirokane_load_apptainer_module() {
  local modulefiles_dir="${GG_SHIROKANE_MODULEFILES_DIR:-/usr/local/package/modulefiles}"

  if command -v apptainer >/dev/null 2>&1; then
    return 0
  fi
  if ! gg_initialize_environment_modules; then
    echo "SHIROKANE requires the Apptainer module on an AGE compute node." >&2
    echo "Submit this command with qsub or enter a compute node with qlogin." >&2
    return 1
  fi
  if [[ ! -d "${modulefiles_dir}" ]]; then
    echo "SHIROKANE module directory is unavailable on this host: ${modulefiles_dir}" >&2
    echo "Apptainer is provided on AGE compute nodes, not on the login node." >&2
    return 1
  fi

  if ! module use "${modulefiles_dir}"; then
    echo "Failed to add the SHIROKANE module directory: ${modulefiles_dir}" >&2
    return 1
  fi
  if ! module load apptainer; then
    echo "Failed to load the SHIROKANE Apptainer module." >&2
    return 1
  fi
  if ! command -v apptainer >/dev/null 2>&1; then
    echo "The SHIROKANE Apptainer module loaded without providing the apptainer command." >&2
    return 1
  fi
}

gg_site_prepend_path_once() {
  local path_entry=${1:-}

  if [[ -z "${path_entry}" ]]; then
    echo "Cannot prepend an empty PATH entry." >&2
    return 1
  fi

  case ":${PATH:-}:" in
    *":${path_entry}:"*)
      return 0
      ;;
  esac

  export PATH="${path_entry}${PATH:+:${PATH}}"
}

# Optional path overrides are used by runtime-discovery tests and by sites
# whose package tree is mounted somewhere other than /opt/pkg.
# shellcheck disable=SC2120
gg_prepend_container_runtime_path() {
  local package_root=${1:-/opt/pkg}
  local legacy_runtime_dir=${2:-/bio/package/singularity/singularity_3.0/bin}
  local runtime_name
  local runtime_bin

  if command -v apptainer >/dev/null 2>&1 \
    || command -v singularity >/dev/null 2>&1; then
    return 0
  fi

  for runtime_name in apptainer singularity; do
    for runtime_bin in \
      "${package_root}/${runtime_name}/bin/${runtime_name}" \
      "${package_root}/${runtime_name}"/*/bin/"${runtime_name}"
    do
      if [[ ! -x "${runtime_bin}" ]]; then
        continue
      fi
      gg_site_prepend_path_once "${runtime_bin%/*}"
      return 0
    done
  done

  for runtime_name in apptainer singularity; do
    runtime_bin="${legacy_runtime_dir}/${runtime_name}"
    if [[ -x "${runtime_bin}" ]]; then
      gg_site_prepend_path_once "${legacy_runtime_dir}"
      return 0
    fi
  done

  return 0
}

# Backward-compatible NIG entry point. Runtime discovery itself is generic and
# can be reused by validation commands on any host with a versioned package tree.
# shellcheck disable=SC2120
gg_nig_prepend_container_runtime_path() {
  gg_prepend_container_runtime_path "${1:-/opt/pkg}" "${2:-/bio/package/singularity/singularity_3.0/bin}"
}

gg_site_scheduler_prelude() {
  case "$(gg_detect_site_profile)" in
    shirokane)
      gg_shirokane_load_apptainer_module || return 1
      ;;
    nig)
      if [[ -n "${PBS_O_WORKDIR:-}" ]]; then
        cd "${PBS_O_WORKDIR}" || return 1
      fi
      gg_nig_prepend_container_runtime_path || return 1
      ;;
    *)
      if [[ -n "${PBS_O_WORKDIR:-}" ]]; then
        cd "${PBS_O_WORKDIR}" || return 1
      fi
      ;;
  esac
}

gg_set_command_array() {
  local out_var=${1:-}
  shift || true

  if [[ -z "${out_var}" ]]; then
    echo "gg_set_command_array: output variable name is required." >&2
    return 1
  fi
  if [[ ! "${out_var}" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]]; then
    echo "gg_set_command_array: invalid output variable name: ${out_var}" >&2
    return 1
  fi

  local arg
  local quoted_arg=""
  eval "${out_var}=()"
  for arg in "$@"; do
    printf -v quoted_arg '%q' "${arg}"
    eval "${out_var}+=( ${quoted_arg} )"
  done
}

gg_bind_native_array_sources() {
  local bind_args_var=$1
  local mode=${GG_INPUT_INPUT_GENERATION_MODE:-}
  local source_path="" kind="" source_list="" source_file="" quoted_bind=""
  local mount_entry="" destination="" filtered_mounts="" matched=""
  local -a source_files=() existing_mounts=() source_list_args=()
  gg_set_command_array "${bind_args_var}" || return 1
  case "${mode}" in
    array_prepare)
      source_path=${GG_INPUT_DOWNLOAD_MANIFEST:-}
      kind=--manifest
      ;;
    array_worker|array_finalize)
      source_path=${GG_INPUT_TASK_PLAN_OUTPUT:-}
      kind=--plan
      ;;
    *) return 0 ;;
  esac
  # Automatic in-workspace discovery keeps its existing behavior. Explicit
  # manifests/plans can also reference source files outside the project bind.
  [[ -n "${source_path}" ]] || return 0
  if [[ "${source_path}" == /workspace/* ]]; then
    source_path="${gg_workspace_dir}/${source_path#/workspace/}"
  fi
  gg_sync_container_bind_envs || return 1
  source_list_args=( "${kind}" "${source_path}" --workspace "${gg_workspace_dir}"
    --coalesce-project-sources "${GG_CONTAINER_PROJECT_ROOT_BIND%%:*}" )
  IFS=',' read -r -a existing_mounts <<< "${GG_CONTAINER_BIND_MOUNTS:-}"
  for mount_entry in "${existing_mounts[@]}"; do
    destination=${mount_entry#*:}
    destination=${destination%%:*}
    [[ -z "${destination}" ]] || source_list_args+=( --existing-bind-destination "${destination}" )
  done
  source_list=$(python "${BASH_SOURCE[0]%/*}/local_input_manifest_binds.py" \
    "${source_list_args[@]}") || return 1
  while IFS= read -r source_file; do
    [[ -n "${source_file}" ]] || continue
    source_files+=( "${source_file}" )
    printf -v quoted_bind '%q' "${source_file}:${source_file}:ro"
    eval "${bind_args_var}+=( --bind ${quoted_bind} )"
  done <<< "${source_list}"
  [[ ${#source_files[@]} -gt 0 ]] || return 0
  # Avoid Linux's per-environment-variable length limit for large manifests.
  # Existing non-source binds keep their precedence; an exact source bind is
  # removed so only the reviewed read-only command argument can provide it.
  gg_sync_container_bind_envs || return 1
  IFS=',' read -r -a existing_mounts <<< "${GG_CONTAINER_BIND_MOUNTS:-}"
  for mount_entry in "${existing_mounts[@]}"; do
    destination=${mount_entry#*:}
    destination=${destination%%:*}
    matched=0
    for source_file in "${source_files[@]}"; do
      if [[ "${destination}" == "${source_file}" ]]; then
        matched=1
        break
      fi
    done
    [[ "${matched}" == 0 ]] || continue
    filtered_mounts="${filtered_mounts:+${filtered_mounts},}${mount_entry}"
  done
  unset SINGULARITY_BIND SINGULARITY_BINDPATH APPTAINER_BIND APPTAINER_BINDPATH
  GG_CONTAINER_BIND_MOUNTS="${filtered_mounts}"
  export GG_CONTAINER_BIND_MOUNTS
  gg_sync_container_bind_envs || return 1
}

gg_bind_requested_project_root() {
  local project_root=${GG_CONTAINER_PROJECT_ROOT_BIND%%:*}
  if [[ "${project_root}" != /* || "${GG_CONTAINER_PROJECT_ROOT_BIND}" != "${project_root}:${project_root}" || \
        "${project_root}" == *','* || "${project_root}" == *'/../'* || "${project_root}" == *'/./'* || \
        ! -d "${project_root}" || -L "${project_root}" ]]; then
    echo "set_singularity_command: invalid GG_CONTAINER_PROJECT_ROOT_BIND" >&2
    return 1
  fi
  gg_add_container_bind_mount "${GG_CONTAINER_PROJECT_ROOT_BIND}"
}

gg_site_container_shell_command() {
  local runtime_bin=$1
  local out_var=${2:-}
  local echo_header="set_singularity_command: "
  local site_profile
  local -a source_bind_args=()

  site_profile="$(gg_detect_site_profile)"
  case "${site_profile}" in
    audrey1)
      echo "${echo_header}site profile = audrey1"
      if [[ -n "${GG_CONTAINER_PROJECT_ROOT_BIND:-}" ]]; then
        gg_bind_requested_project_root || return 1
        gg_bind_native_array_sources source_bind_args || return 1
        gg_set_command_array "${out_var}" "${runtime_bin}" exec --contain "${source_bind_args[@]}" || return 1
      else
        gg_set_command_array "${out_var}" "${runtime_bin}" exec || return 1
      fi
      ;;
    shirokane)
      echo "${echo_header}site profile = shirokane"
      gg_set_command_array "${out_var}" "${runtime_bin}" exec || return 1
      ;;
    nig)
      echo "${echo_header}site profile = nig"
      if [[ -e /var/spool/uge ]]; then
        gg_add_container_bind_mount "/var/spool/uge:/var/spool/uge"
      fi
      if [[ -e /var/spool/age ]]; then
        gg_add_container_bind_mount "/var/spool/age:/var/spool/age"
      fi
      if [[ -e /opt/pkg ]]; then
        gg_add_container_bind_mount "/opt/pkg:/opt/pkg"
      fi
      if [[ -e /home/geadmin/UGER/uger/spool ]]; then
        gg_add_container_bind_mount "/home/geadmin/UGER/uger/spool:/home/geadmin/UGER/uger/spool"
      fi
      if [[ -n "${GG_CONTAINER_PROJECT_ROOT_BIND:-}" ]]; then
        gg_bind_requested_project_root || return 1
        gg_bind_native_array_sources source_bind_args || return 1
      fi
      gg_set_command_array "${out_var}" "${runtime_bin}" exec "${source_bind_args[@]}" || return 1
      ;;
    nhr-fau)
      echo "${echo_header}site profile = nhr-fau"
      gg_set_command_array "${out_var}" "${runtime_bin}" exec --contain || return 1
      ;;
    *)
      echo "${echo_header}site profile = default"
      gg_set_command_array "${out_var}" "${runtime_bin}" exec || return 1
      ;;
  esac
}

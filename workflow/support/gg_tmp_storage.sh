#!/usr/bin/env bash

# Resolve on the execution host, before any container mounts are constructed.
gg_tmp_is_nig() {
  if [[ -n "${GG_SITE_PROFILE:-}" ]]; then
    [[ "${GG_SITE_PROFILE}" == nig ]]
    return
  fi
  [[ "${SGE_ROOT:-}" != /home/geadmin/N1GE ]] || return 1
  local host
  host=$(hostname 2>/dev/null) || return 1
  [[ "${host}" =~ ^(a[0-9]+|m[0-9]+|at[0-9]+|igt[0-9]+|it[0-9]+)(\.ib\.cluster|\.ddbj\.nig\.ac\.jp)?$ ]] \
    && [[ -d /lustre10/home ]]
}

gg_tmp_validate_nig_data1() {
  local root=${1:-/data1}
  [[ -d "${root}" && -w "${root}" && -x "${root}" ]] || {
    echo "NIG auto scratch requires an existing writable /data1; select workspace explicitly to use it." >&2
    return 1
  }
  local scratch_device system_device
  scratch_device=$(stat -c %d "${root}") || return 1
  system_device=$(stat -c %d /) || return 1
  [[ "${scratch_device}" != "${system_device}" ]] || {
    echo "NIG auto scratch refuses /data1 on the operating-system filesystem." >&2
    return 1
  }
}

gg_resolve_tmp_root() {
  local requested="${GG_COMMON_TMP_ROOT:-auto}"
  case "${requested}" in
    auto)
      if gg_tmp_is_nig; then
        gg_tmp_validate_nig_data1 /data1 || return 1
        requested=/data1
      else
        requested=workspace
      fi
      ;;
    env)
      requested="${TMPDIR:-}"
      [[ -n "${requested}" ]] || {
        echo "GG_COMMON_TMP_ROOT=env requires TMPDIR on the execution node." >&2
        return 1
      }
      ;;
  esac
  if [[ "${requested}" == workspace ]]; then
    printf '%s\n' workspace
    return 0
  fi
  if [[ "${requested}" != /* || "${requested}" == *[:,]* || "${requested}" == *$'\n'* ]]; then
    echo "Scratch root must be an absolute path without colons, commas or newlines: ${requested}" >&2
    return 1
  fi
  if [[ ! -d "${requested}" || ! -w "${requested}" || ! -x "${requested}" ]]; then
    echo "Scratch root must be an existing writable directory: ${requested}" >&2
    return 1
  fi
  local resolved
  resolved=$(cd -P -- "${requested}" && printf '%s.' "$PWD") || return 1
  resolved=${resolved%.}
  if [[ "${resolved}" == *[:,]* || "${resolved}" == *$'\n'* ]]; then
    echo "Resolved scratch root contains a container bind delimiter: ${resolved}" >&2
    return 1
  fi
  printf '%s\n' "${resolved}"
}

# Keep auxiliary tempfile and bytecode writes off the system /tmp as well.
gg_private_runtime_tmp() {
  local root=$1 path
  if [[ "${root}" == workspace ]]; then
    root="${gg_workspace_dir:?Workspace is required for temporary storage}"
    (umask 077; mkdir -p -- "${root}") || return 1
    root=$(cd -P -- "${root}" && pwd -P) || return 1
  fi
  path="${root%/}/.genegalleon-runtime-$(id -u)"
  if [[ -L "${path}" || ( -e "${path}" && ( ! -d "${path}" || ! -O "${path}" ) ) ]]; then
    echo "Refusing unsafe GeneGalleon runtime temporary directory: ${path}" >&2
    return 1
  fi
  (umask 077; mkdir -p -- "${path}") || return 1
  [[ ! -L "${path}" && -d "${path}" && -O "${path}" ]] || return 1
  chmod 700 -- "${path}" || return 1
  local child
  for child in tmp pycache; do
    if [[ -L "${path}/${child}" || ( -e "${path}/${child}" && ( ! -d "${path}/${child}" || ! -O "${path}/${child}" ) ) ]]; then
      echo "Refusing unsafe runtime temporary child: ${path}/${child}" >&2
      return 1
    fi
    (umask 077; mkdir -p -- "${path}/${child}") || return 1
    chmod 700 -- "${path}/${child}" || return 1
  done
  printf '%s\n' "${path}"
}

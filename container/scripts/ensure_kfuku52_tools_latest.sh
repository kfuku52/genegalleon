#!/usr/bin/env bash
set -euo pipefail

env_name=${1:-base}
repo_owner=${KFU52_REPO_OWNER:-kfuku52}
tools_csv=${KFU52_TOOLS:-amalgkit,cdskit,csubst,nwkit}
repo_ref=${KFU52_REPO_REF:-master}
amalgkit_branch_candidates=${KFU52_AMALGKIT_BRANCH_CANDIDATES:-master,kfdevel}
amalgkit_auto_select_ref=${KFU52_AMALGKIT_AUTO_SELECT_REF:-1}
amalgkit_repo_ref_override=${KFU52_AMALGKIT_REPO_REF:-}
max_attempts=${KFU52_MAX_ATTEMPTS:-5}
retry_sleep_sec=${KFU52_RETRY_SLEEP_SEC:-15}

log() {
  echo "[ensure_kfuku52_tools_latest] $*"
}

print_tool_version() {
  local tool="$1"
  micromamba run -n "${env_name}" python - "${tool}" <<'PY'
import importlib
import sys

tool_name = sys.argv[1]
try:
    module = importlib.import_module(tool_name)
    print(getattr(module, "__version__", "unknown"))
except Exception:
    print("not-installed")
PY
}

install_latest_from_github() {
  local tool="$1"
  local tool_ref="$2"
  local spec
  spec="git+https://github.com/${repo_owner}/${tool}.git@${tool_ref}"
  log "Installing ${tool} from GitHub ${tool_ref}: ${spec}"

  local attempt
  for (( attempt=1; attempt<=max_attempts; attempt++ )); do
    if micromamba run -n "${env_name}" \
      pip install --no-cache-dir --upgrade --force-reinstall --no-deps --no-build-isolation "${spec}"; then
      return 0
    fi
    if [[ "${attempt}" -lt "${max_attempts}" ]]; then
      log "Retrying ${tool} install in ${retry_sleep_sec}s (${attempt}/${max_attempts} failed)"
      sleep "${retry_sleep_sec}"
    fi
  done

  return 1
}

branch_exists_on_remote() {
  local remote_url="$1"
  local branch="$2"
  git ls-remote --heads --exit-code "${remote_url}" "${branch}" >/dev/null 2>&1
}

get_remote_branch_commit_epoch() {
  local remote_url="$1"
  local branch="$2"
  local tmp_dir
  local commit_epoch=""

  tmp_dir=$(mktemp -d)
  if git clone --quiet --depth 1 --single-branch --branch "${branch}" "${remote_url}" "${tmp_dir}/repo" >/dev/null 2>&1; then
    commit_epoch=$(git -C "${tmp_dir}/repo" log -1 --format=%ct 2>/dev/null || true)
  fi
  rm -rf -- "${tmp_dir}"

  echo "${commit_epoch}"
}

resolve_amalgkit_ref_by_latest_commit() {
  local remote_url
  remote_url="https://github.com/${repo_owner}/amalgkit.git"

  local best_ref=""
  local best_epoch=-1
  local ref
  local epoch
  local refs=()
  mapfile -t refs < <(printf '%s' "${amalgkit_branch_candidates}" | tr ',' '\n')
  for ref in "${refs[@]}"; do
    ref=$(printf '%s' "${ref}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
    if [[ -z "${ref}" ]]; then
      continue
    fi
    if ! branch_exists_on_remote "${remote_url}" "${ref}"; then
      log "Skipping missing amalgkit branch on remote: ${ref}" >&2
      continue
    fi
    epoch=$(get_remote_branch_commit_epoch "${remote_url}" "${ref}")
    if [[ ! "${epoch}" =~ ^[0-9]+$ ]]; then
      log "Could not resolve commit timestamp for amalgkit branch: ${ref}" >&2
      continue
    fi
    log "amalgkit branch ${ref} latest commit epoch: ${epoch}" >&2
    if (( epoch > best_epoch )); then
      best_epoch=${epoch}
      best_ref=${ref}
    fi
  done

  echo "${best_ref}"
}

resolve_tool_repo_ref() {
  local tool="$1"

  if [[ "${tool}" != "amalgkit" ]]; then
    echo "${repo_ref}"
    return 0
  fi

  if [[ -n "${amalgkit_repo_ref_override}" ]]; then
    log "Using amalgkit ref override from KFU52_AMALGKIT_REPO_REF: ${amalgkit_repo_ref_override}" >&2
    echo "${amalgkit_repo_ref_override}"
    return 0
  fi

  if [[ "${amalgkit_auto_select_ref}" != "1" ]]; then
    log "amalgkit auto-ref selection disabled; using KFU52_REPO_REF=${repo_ref}" >&2
    echo "${repo_ref}"
    return 0
  fi

  local selected_ref
  selected_ref=$(resolve_amalgkit_ref_by_latest_commit)
  if [[ -n "${selected_ref}" ]]; then
    log "Selected amalgkit branch by latest commit: ${selected_ref}" >&2
    echo "${selected_ref}"
    return 0
  fi

  log "Falling back to KFU52_REPO_REF for amalgkit: ${repo_ref}" >&2
  echo "${repo_ref}"
}

ensure_tool_available() {
  local tool="$1"
  if ! micromamba run -n "${env_name}" bash -lc "command -v ${tool} >/dev/null 2>&1"; then
    log "ERROR: command '${tool}' was not found in env '${env_name}' after GitHub install."
    exit 1
  fi
}

main() {
  local tools_normalized
  tools_normalized=$(echo "${tools_csv}" | tr ',' ' ')

  local tool
  for tool in ${tools_normalized}; do
    local before_version
    local after_version
    local tool_ref
    before_version=$(print_tool_version "${tool}")
    log "Current ${tool} version: ${before_version}"

    tool_ref=$(resolve_tool_repo_ref "${tool}")
    install_latest_from_github "${tool}" "${tool_ref}"
    ensure_tool_available "${tool}"

    after_version=$(print_tool_version "${tool}")
    log "${tool} version after GitHub upgrade: ${after_version}"
  done
}

main "$@"

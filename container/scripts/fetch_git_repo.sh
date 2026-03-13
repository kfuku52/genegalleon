#!/usr/bin/env bash
set -euo pipefail

repo_url="${1:-}"
repo_ref="${2:-}"
dest_dir="${3:-}"
fallback_repo_url="${4:-}"

if [[ -z "${repo_url}" || -z "${dest_dir}" ]]; then
  echo "usage: $0 <repo-url> <repo-ref> <dest-dir> [fallback-repo-url]" >&2
  exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
download_url_script="${script_dir}/download_url.sh"
max_attempts=${FETCH_GIT_MAX_ATTEMPTS:-5}
retry_sleep_sec=${FETCH_GIT_RETRY_SLEEP_SEC:-15}

log() {
  echo "[fetch_git_repo] $*" >&2
}

github_archive_url() {
  local url="$1"
  local ref="$2"
  if [[ "${url}" =~ ^https://github\.com/([^/]+)/([^/]+)(\.git)?$ ]]; then
    local owner="${BASH_REMATCH[1]}"
    local repo_name="${BASH_REMATCH[2]%.git}"
    printf 'https://codeload.github.com/%s/%s/tar.gz/%s\n' "${owner}" "${repo_name}" "${ref}"
    return 0
  fi
  return 1
}

gitlab_archive_url() {
  local url="$1"
  local ref="$2"
  if [[ "${url}" =~ ^https://gitlab\.com/(.+)(\.git)?$ ]]; then
    local project_path="${BASH_REMATCH[1]%.git}"
    local repo_name="${project_path##*/}"
    printf 'https://gitlab.com/%s/-/archive/%s/%s-%s.tar.gz\n' "${project_path}" "${ref}" "${repo_name}" "${ref}"
    return 0
  fi
  return 1
}

download_archive() {
  local url="$1"
  local ref="$2"
  local dest="$3"
  local archive_url=""
  local archive_path

  if archive_url="$(github_archive_url "${url}" "${ref}")"; then
    :
  elif archive_url="$(gitlab_archive_url "${url}" "${ref}")"; then
    :
  else
    return 1
  fi

  archive_path="$(mktemp)"
  rm -rf -- "${dest}"
  mkdir -p "${dest}"
  if ! bash "${download_url_script}" "${archive_url}" "${archive_path}"; then
    rm -f -- "${archive_path}"
    rm -rf -- "${dest}"
    return 1
  fi
  if ! tar -xzf "${archive_path}" -C "${dest}" --strip-components=1; then
    rm -f -- "${archive_path}"
    rm -rf -- "${dest}"
    return 1
  fi
  rm -f -- "${archive_path}"
  return 0
}

git_clone_ref() {
  local url="$1"
  local ref="$2"
  local dest="$3"

  rm -rf -- "${dest}"
  if [[ -z "${ref}" ]]; then
    git clone --depth 1 --recursive "${url}" "${dest}" >/dev/null
    return 0
  fi

  mkdir -p "${dest}"
  git -C "${dest}" init -q
  git -C "${dest}" remote add origin "${url}"
  git -C "${dest}" fetch --depth 1 origin "${ref}" >/dev/null
  git -C "${dest}" checkout -q --detach FETCH_HEAD
  if [[ -f "${dest}/.gitmodules" ]]; then
    git -C "${dest}" submodule update --init --recursive >/dev/null
  fi
}

fetch_from_source() {
  local url="$1"
  local ref="$2"
  local dest="$3"
  local attempt

  if [[ -n "${ref}" ]]; then
    if download_archive "${url}" "${ref}" "${dest}"; then
      log "Fetched ${url}@${ref} via archive"
      return 0
    fi
    log "Archive fetch failed for ${url}@${ref}; falling back to git"
  fi

  for (( attempt=1; attempt<=max_attempts; attempt++ )); do
    if git_clone_ref "${url}" "${ref}" "${dest}"; then
      if [[ -n "${ref}" ]]; then
        log "Fetched ${url}@${ref} via git"
      else
        log "Fetched ${url} default branch via git"
      fi
      return 0
    fi
    rm -rf -- "${dest}"
    if [[ "${attempt}" -lt "${max_attempts}" ]]; then
      if [[ -n "${ref}" ]]; then
        log "Retrying git fetch for ${url}@${ref} in ${retry_sleep_sec}s (${attempt}/${max_attempts} failed)"
      else
        log "Retrying git fetch for ${url} in ${retry_sleep_sec}s (${attempt}/${max_attempts} failed)"
      fi
      sleep "${retry_sleep_sec}"
    fi
  done

  return 1
}

main() {
  local sources=("${repo_url}")
  if [[ -n "${fallback_repo_url}" && "${fallback_repo_url}" != "${repo_url}" ]]; then
    sources+=("${fallback_repo_url}")
  fi

  local source_url
  for source_url in "${sources[@]}"; do
    if fetch_from_source "${source_url}" "${repo_ref}" "${dest_dir}"; then
      return 0
    fi
  done

  if [[ -n "${repo_ref}" ]]; then
    log "ERROR: failed to fetch ${repo_url}@${repo_ref}"
  else
    log "ERROR: failed to fetch ${repo_url}"
  fi
  exit 1
}

main "$@"

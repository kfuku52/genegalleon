#!/usr/bin/env bash
set -euo pipefail

repo_root=${1:-.}
dockerfile="${repo_root}/container/Dockerfile"
dropped_envs_file="${repo_root}/container/spec/dropped_envs.txt"

if [[ ! -f "${dockerfile}" ]]; then
  echo "Dockerfile not found: ${dockerfile}"
  exit 1
fi

dropped_envs=()
if [[ -f "${dropped_envs_file}" ]]; then
  while IFS= read -r line; do
    env_name=$(echo "${line%%#*}" | xargs)
    [[ -n "${env_name}" ]] || continue
    dropped_envs+=("${env_name}")
  done < "${dropped_envs_file}"
fi

is_dropped_env() {
  local env_name=$1
  local dropped
  for dropped in "${dropped_envs[@]-}"; do
    if [[ "${env_name}" == "${dropped}" ]]; then
      return 0
    fi
  done
  return 1
}

search_lines() {
  local pattern=$1
  shift
  local matches
  local status

  if matches="$(grep -nE -- "${pattern}" "$@")"; then
    printf '%s\n' "${matches}"
  else
    status=$?
    if [[ ${status} -ne 1 ]]; then
      return "${status}"
    fi
  fi
}

core_scripts=("${repo_root}/workflow/core"/gg_*_core.sh)
activate_matches="$(
  search_lines "conda activate [A-Za-z0-9_.-]+" "${core_scripts[@]}"
)"
bootstrap_matches="$(
  search_lines \
    'gg_bootstrap_core_runtime[[:space:]]+[^[:space:]]+[[:space:]]+"[A-Za-z0-9_.-]+"' \
    "${core_scripts[@]}"
)"
prepare_matches="$(
  search_lines \
    'gg_prepare_cmd_runtime[[:space:]]+[^[:space:]]+[[:space:]]+"[A-Za-z0-9_.-]+"' \
    "${core_scripts[@]}"
)"
pipeline_env_output="$(
  {
    printf '%s\n' "${activate_matches}" \
      | sed -E -n 's/.*conda activate ([A-Za-z0-9_.-]+).*/\1/p'
    printf '%s\n' "${bootstrap_matches}" \
      | sed -E -n 's/.*gg_bootstrap_core_runtime[[:space:]]+[^[:space:]]+[[:space:]]+"([A-Za-z0-9_.-]+)".*/\1/p'
    printf '%s\n' "${prepare_matches}" \
      | sed -E -n 's/.*gg_prepare_cmd_runtime[[:space:]]+[^[:space:]]+[[:space:]]+"([A-Za-z0-9_.-]+)".*/\1/p'
  } | sort -u
)"

pipeline_envs=()
while IFS= read -r env_name; do
  [[ -n "${env_name}" ]] || continue
  is_dropped_env "${env_name}" && continue
  pipeline_envs+=("${env_name}")
done <<< "${pipeline_env_output}"

docker_env_matches="$(
  search_lines "/opt/pg/scripts/install_env.sh [A-Za-z0-9_.-]+" "${dockerfile}"
)"
docker_env_output="$(
  printf '%s\n' "${docker_env_matches}" \
    | sed -E -n 's/.*install_env.sh ([A-Za-z0-9_.-]+).*/\1/p' \
    | sort -u
)"

docker_envs=()
while IFS= read -r env_name; do
  [[ -n "${env_name}" ]] || continue
  is_dropped_env "${env_name}" && continue
  docker_envs+=("${env_name}")
done <<< "${docker_env_output}"

contains_env() {
  local target=$1
  shift
  local candidate
  for candidate in "$@"; do
    if [[ "${candidate}" == "${target}" ]]; then
      return 0
    fi
  done
  return 1
}

is_pipeline_env_covered() {
  local env_name=$1
  contains_env "${env_name}" "${docker_envs[@]-}"
}

missing=0
for env_name in "${pipeline_envs[@]-}"; do
  if ! is_pipeline_env_covered "${env_name}"; then
    echo "Missing conda env in Dockerfile: ${env_name}"
    missing=1
  fi
done

unused=0
for env_name in "${docker_envs[@]-}"; do
  if contains_env "${env_name}" "${pipeline_envs[@]-}"; then
    continue
  fi
  if ! contains_env "${env_name}" "${pipeline_envs[@]-}"; then
    echo "Unused conda env in Dockerfile: ${env_name}"
    unused=1
  fi
done

if [[ ${missing} -ne 0 || ${unused} -ne 0 ]]; then
  exit 1
fi

echo "Conda environment coverage matches between pipeline scripts and container/Dockerfile."

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

pipeline_envs=()
while IFS= read -r env_name; do
  [[ -n "${env_name}" ]] || continue
  is_dropped_env "${env_name}" && continue
  pipeline_envs+=("${env_name}")
done < <(
  rg -n "conda activate [A-Za-z0-9_.-]+" "${repo_root}/workflow"/gg_*_cmd.sh \
    | sed -E 's/.*conda activate ([A-Za-z0-9_.-]+).*/\1/' \
    | sort -u
)

docker_envs=()
while IFS= read -r env_name; do
  [[ -n "${env_name}" ]] || continue
  is_dropped_env "${env_name}" && continue
  docker_envs+=("${env_name}")
done < <(
  rg -n "/opt/pg/scripts/install_env.sh [A-Za-z0-9_.-]+" "${dockerfile}" \
    | sed -E 's/.*install_env.sh ([A-Za-z0-9_.-]+).*/\1/' \
    | sort -u
)

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

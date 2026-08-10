#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
runtime="${GG_TEST_RUNTIME:-auto}"
sif_path="${GENEGALLEON_SIF:-${gg_container_image_path:-${repo_root}/genegalleon.sif}}"
docker_image="${GG_CONTAINER_DOCKER_IMAGE:-local/genegalleon:dev}"

usage() {
  cat <<'EOF'
Usage:
  bash workflow/tests/run_in_runtime.sh <command> [args...]

Selects an authoritative GeneGalleon validation runtime:
  GG_TEST_RUNTIME=auto    Prefer an available SIF runtime, otherwise Docker.
  GG_TEST_RUNTIME=sif     Require Apptainer/Singularity and a GeneGalleon SIF.
  GG_TEST_RUNTIME=docker  Require Docker and GG_CONTAINER_DOCKER_IMAGE.

Docker defaults to local/genegalleon:dev. SIF defaults to ./genegalleon.sif.
EOF
}

if [[ $# -eq 0 ]]; then
  usage >&2
  exit 2
fi

case "${runtime}" in
  auto|sif|docker)
    ;;
  *)
    echo "GG_TEST_RUNTIME must be one of: auto, sif, docker" >&2
    exit 2
    ;;
esac

has_sif_runtime() {
  [[ -f "${sif_path}" ]] \
    && { command -v apptainer >/dev/null 2>&1 || command -v singularity >/dev/null 2>&1; }
}

has_docker_runtime() {
  command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1
}

if [[ "${runtime}" == "auto" ]]; then
  if has_sif_runtime; then
    runtime="sif"
  elif has_docker_runtime; then
    runtime="docker"
  else
    cat >&2 <<EOF
No GeneGalleon validation runtime is available.

Provide ${sif_path} with Apptainer/Singularity, or start Docker and build:
  BUILD_SIF=0 IMAGE_SOURCE=local IMAGE=local/genegalleon TAG=dev bash ./gg_container_build_entrypoint.sh
EOF
    exit 1
  fi
fi

if [[ "${runtime}" == "sif" ]]; then
  exec bash "${script_dir}/run_in_sif.sh" "$@"
fi

if ! has_docker_runtime; then
  echo "Docker is unavailable; cannot use GG_TEST_RUNTIME=docker." >&2
  exit 1
fi
if ! docker image inspect "${docker_image}" >/dev/null 2>&1; then
  cat >&2 <<EOF
GeneGalleon Docker image was not found: ${docker_image}

Build the current repository with:
  BUILD_SIF=0 IMAGE_SOURCE=local IMAGE=local/genegalleon TAG=dev bash ./gg_container_build_entrypoint.sh
EOF
  exit 1
fi

bind_args=(--volume "${repo_root}:${repo_root}")
if [[ -n "${GENEGALLEON_DOCKER_EXTRA_BINDS:-}" ]]; then
  while IFS= read -r bind_path; do
    [[ -z "${bind_path}" ]] && continue
    if [[ "${bind_path}" != /* || "${bind_path}" == *:* || ! -e "${bind_path}" ]]; then
      echo "Invalid GENEGALLEON_DOCKER_EXTRA_BINDS path: ${bind_path}" >&2
      exit 1
    fi
    bind_args+=(--volume "${bind_path}:${bind_path}:ro")
  done <<< "${GENEGALLEON_DOCKER_EXTRA_BINDS}"
fi

exec docker run --rm -i \
  "${bind_args[@]}" \
  --workdir "${repo_root}" \
  "${docker_image}" \
  "$@"

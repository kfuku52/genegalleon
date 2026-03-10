#!/usr/bin/env bash
set -euo pipefail

# Convert a built Docker image to Apptainer/Singularity SIF.
# Run this on a host that has apptainer/singularity installed.
#
# Example:
#   IMAGE=ghcr.io/your-org/genegalleon TAG=20260211 ./container/apptainer_from_docker.sh
#   SOURCE=docker-daemon IMAGE=local/genegalleon TAG=dev ./container/apptainer_from_docker.sh

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

IMAGE=${IMAGE:-ghcr.io/example/genegalleon}
TAG=${TAG:-dev}
ENGINE=${ENGINE:-auto} # auto | apptainer | singularity
SOURCE=${SOURCE:-docker} # docker | docker-daemon
OUT=${OUT:-${repo_root}/genegalleon.sif}

resolve_container_engine() {
  if [[ "${ENGINE}" == "apptainer" || "${ENGINE}" == "singularity" ]]; then
    if command -v "${ENGINE}" >/dev/null 2>&1; then
      echo "${ENGINE}"
      return 0
    fi
    echo "Container engine not found: ${ENGINE}" >&2
    return 1
  fi
  if command -v apptainer >/dev/null 2>&1; then
    echo apptainer
    return 0
  fi
  if command -v singularity >/dev/null 2>&1; then
    echo singularity
    return 0
  fi
  echo "Neither apptainer nor singularity was found on PATH." >&2
  return 1
}

resolved_engine="$(resolve_container_engine)"
case "${SOURCE}" in
  docker)
    source_uri="docker://${IMAGE}:${TAG}"
    ;;
  docker-daemon)
    source_uri="docker-daemon:${IMAGE}:${TAG}"
    ;;
  *)
    echo "SOURCE must be one of: docker, docker-daemon" >&2
    exit 1
    ;;
esac

"${resolved_engine}" build "${OUT}" "${source_uri}"
echo "Generated ${OUT}"

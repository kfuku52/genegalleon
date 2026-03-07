#!/usr/bin/env bash
set -euo pipefail

# Convert a built Docker image to Apptainer/Singularity SIF.
# Run this on a host that has apptainer/singularity installed.
#
# Example:
#   IMAGE=ghcr.io/your-org/genegalleon TAG=20260211 OUT=genegalleon_20260211_amd.sif ./container/apptainer_from_docker.sh

IMAGE=${IMAGE:-ghcr.io/example/genegalleon}
TAG=${TAG:-dev}
ENGINE=${ENGINE:-auto} # auto | apptainer | singularity

build_date="$(date '+%Y%m%d')"
host_arch="$(uname -m)"
case "${host_arch}" in
  aarch64|arm64)
    arch_tag="arm"
    ;;
  x86_64|amd64)
    arch_tag="amd"
    ;;
  *)
    arch_tag="$(echo "${host_arch}" | tr '[:upper:]' '[:lower:]' | tr -cd '[:alnum:]_-')"
    ;;
esac
OUT=${OUT:-genegalleon_${build_date}_${arch_tag}.sif}

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
"${resolved_engine}" build "${OUT}" "docker://${IMAGE}:${TAG}"
echo "Generated ${OUT}"

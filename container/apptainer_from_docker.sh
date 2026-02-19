#!/usr/bin/env bash
set -euo pipefail

# Convert a built Docker image to Apptainer/Singularity SIF.
# Run this on a host that has apptainer/singularity installed.
#
# Example:
#   IMAGE=ghcr.io/your-org/genegalleon TAG=20260211 OUT=genegalleon_20260211_amd.sif ./container/apptainer_from_docker.sh

IMAGE=${IMAGE:-ghcr.io/example/genegalleon}
TAG=${TAG:-dev}
ENGINE=${ENGINE:-apptainer} # apptainer | singularity

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

if ! command -v "${ENGINE}" >/dev/null 2>&1; then
  echo "Container engine not found: ${ENGINE}"
  exit 1
fi

"${ENGINE}" build "${OUT}" "docker://${IMAGE}:${TAG}"
echo "Generated ${OUT}"

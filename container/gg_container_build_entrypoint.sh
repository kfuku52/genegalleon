#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
buildx_script="${script_dir}/buildx.sh"
apptainer_script="${script_dir}/apptainer_from_docker.sh"

if [[ ! -f "${buildx_script}" ]]; then
  echo "Missing build script: ${buildx_script}"
  exit 1
fi
if [[ ! -f "${apptainer_script}" ]]; then
  echo "Missing conversion script: ${apptainer_script}"
  exit 1
fi

host_arch="$(uname -m)"
case "${host_arch}" in
  aarch64|arm64)
    default_platform="linux/arm64"
    ;;
  x86_64|amd64)
    default_platform="linux/amd64"
    ;;
  *)
    echo "Unsupported host architecture for default PLATFORMS: ${host_arch}"
    echo "Set PLATFORMS explicitly (example: linux/amd64 or linux/arm64)."
    exit 1
    ;;
esac

default_build_sif=1
if [[ "$(uname -s)" == "Darwin" ]]; then
  # Apptainer is typically unavailable on macOS hosts.
  default_build_sif=0
fi

IMAGE=${IMAGE:-local/genegalleon}
TAG=${TAG:-dev}
PLATFORMS=${PLATFORMS:-${default_platform}}
MODE=${MODE:-load} # push | load
BUILD_SIF=${BUILD_SIF:-${default_build_sif}} # 1 | 0
ENGINE=${ENGINE:-apptainer} # apptainer | singularity
OUT=${OUT:-${repo_root}/genegalleon.sif}

if [[ "${MODE}" != "push" && "${MODE}" != "load" ]]; then
  echo "MODE must be one of: push, load"
  exit 1
fi
if [[ "${MODE}" == "load" && "${PLATFORMS}" == *","* ]]; then
  echo "MODE=load requires a single platform. Set PLATFORMS=linux/amd64 or linux/arm64."
  exit 1
fi
if [[ "${BUILD_SIF}" != "0" && "${BUILD_SIF}" != "1" ]]; then
  echo "BUILD_SIF must be 0 or 1."
  exit 1
fi

echo "[gg_container_build] repo_root=${repo_root}"
echo "[gg_container_build] image=${IMAGE}:${TAG}"
echo "[gg_container_build] platforms=${PLATFORMS} mode=${MODE}"

echo "[gg_container_build] step 1/2: docker buildx build"
(
  cd "${repo_root}"
  IMAGE="${IMAGE}" TAG="${TAG}" PLATFORMS="${PLATFORMS}" MODE="${MODE}" \
    bash "${buildx_script}"
)

if [[ "${BUILD_SIF}" == "1" ]]; then
  mkdir -p "$(dirname "${OUT}")"
  echo "[gg_container_build] step 2/2: convert docker image to SIF"
  (
    cd "${repo_root}"
    IMAGE="${IMAGE}" TAG="${TAG}" ENGINE="${ENGINE}" OUT="${OUT}" \
      bash "${apptainer_script}"
  )
  echo "[gg_container_build] done: ${OUT}"
else
  echo "[gg_container_build] BUILD_SIF=0, skipped SIF conversion."
fi

#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
buildx_script="${script_dir}/buildx.sh"
apptainer_script="${script_dir}/apptainer_from_docker.sh"
apptainer_local_build_script="${script_dir}/apptainer_local_build.sh"

if [[ ! -f "${buildx_script}" ]]; then
  echo "Missing build script: ${buildx_script}"
  exit 1
fi
if [[ ! -f "${apptainer_script}" ]]; then
  echo "Missing conversion script: ${apptainer_script}"
  exit 1
fi
if [[ ! -f "${apptainer_local_build_script}" ]]; then
  echo "Missing native build script: ${apptainer_local_build_script}"
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

image_was_default=0
tag_was_default=0
if [[ -z "${IMAGE+x}" ]]; then
  image_was_default=1
fi
if [[ -z "${TAG+x}" ]]; then
  tag_was_default=1
fi

IMAGE_SOURCE=${IMAGE_SOURCE:-auto} # auto | local | public
LOCAL_IMAGE_DEFAULT=${LOCAL_IMAGE_DEFAULT:-local/genegalleon}
LOCAL_TAG_DEFAULT=${LOCAL_TAG_DEFAULT:-dev}
PUBLIC_IMAGE_DEFAULT=${PUBLIC_IMAGE_DEFAULT:-ghcr.io/kfuku52/genegalleon}
PUBLIC_TAG_DEFAULT=${PUBLIC_TAG_DEFAULT:-latest}

case "${IMAGE_SOURCE}" in
  auto|local)
    default_image="${LOCAL_IMAGE_DEFAULT}"
    default_tag="${LOCAL_TAG_DEFAULT}"
    ;;
  public)
    default_image="${PUBLIC_IMAGE_DEFAULT}"
    default_tag="${PUBLIC_TAG_DEFAULT}"
    ;;
  *)
    echo "IMAGE_SOURCE must be one of: auto, local, public"
    exit 1
    ;;
esac

IMAGE=${IMAGE:-${default_image}}
TAG=${TAG:-${default_tag}}
PLATFORMS=${PLATFORMS:-${default_platform}}
MODE=${MODE:-load} # push | load
BUILD_SIF=${BUILD_SIF:-${default_build_sif}} # 1 | 0
ENGINE=${ENGINE:-auto} # auto | apptainer | singularity
OUT=${OUT:-${repo_root}/genegalleon.sif}
FALLBACK_REMOTE_IMAGE=${FALLBACK_REMOTE_IMAGE:-${PUBLIC_IMAGE_DEFAULT}}
FALLBACK_REMOTE_TAG=${FALLBACK_REMOTE_TAG:-${PUBLIC_TAG_DEFAULT}}

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

resolve_container_engine() {
  if [[ "${ENGINE}" == "apptainer" || "${ENGINE}" == "singularity" ]]; then
    if command -v "${ENGINE}" >/dev/null 2>&1; then
      echo "${ENGINE}"
      return 0
    fi
    echo "Requested container engine not found on PATH: ${ENGINE}" >&2
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

docker_buildx_available() {
  command -v docker >/dev/null 2>&1 && docker buildx ls >/dev/null 2>&1
}

image_looks_remote() {
  local image_ref="$1"
  local first_component
  if [[ "${image_ref}" != */* ]]; then
    return 0
  fi
  first_component="${image_ref%%/*}"
  [[ "${first_component}" == "localhost" || "${first_component}" == *.* || "${first_component}" == *:* ]]
}

echo "[gg_container_build] repo_root=${repo_root}"
echo "[gg_container_build] image=${IMAGE}:${TAG}"
echo "[gg_container_build] image_source=${IMAGE_SOURCE}"
echo "[gg_container_build] platforms=${PLATFORMS} mode=${MODE}"

did_buildx=0
did_native_local_build=0
resolved_image_source="${IMAGE_SOURCE}"
sif_source="docker"
case "${IMAGE_SOURCE}" in
  local)
    if docker_buildx_available; then
      echo "[gg_container_build] step 1/2: docker buildx build"
      (
        cd "${repo_root}"
        IMAGE="${IMAGE}" TAG="${TAG}" PLATFORMS="${PLATFORMS}" MODE="${MODE}" \
          bash "${buildx_script}"
      )
      did_buildx=1
      if [[ "${MODE}" == "load" ]]; then
        sif_source="docker-daemon"
      fi
    elif [[ "${BUILD_SIF}" == "1" ]]; then
      resolved_engine="$(resolve_container_engine)"
      echo "[gg_container_build] docker buildx not found; using ${resolved_engine} for native local build."
      echo "[gg_container_build] engine=${resolved_engine}"
      echo "[gg_container_build] step 1/1: native local build from repository"
      (
        cd "${repo_root}"
        IMAGE="${IMAGE}" TAG="${TAG}" PLATFORMS="${PLATFORMS}" ENGINE="${resolved_engine}" OUT="${OUT}" \
          bash "${apptainer_local_build_script}"
      )
      did_native_local_build=1
    else
      echo "[gg_container_build] IMAGE_SOURCE=local without docker buildx requires BUILD_SIF=1."
      echo "[gg_container_build] Native local builds produce a SIF directly."
      exit 1
    fi
    ;;
  public)
    if [[ "${BUILD_SIF}" != "1" ]]; then
      echo "[gg_container_build] IMAGE_SOURCE=public requires BUILD_SIF=1."
      exit 1
    fi
    if ! image_looks_remote "${IMAGE}"; then
      echo "[gg_container_build] IMAGE_SOURCE=public requires a registry image reference."
      echo "[gg_container_build] Example: IMAGE=ghcr.io/kfuku52/genegalleon TAG=latest"
      exit 1
    fi
    echo "[gg_container_build] skipping local Docker build and using published image ${IMAGE}:${TAG}."
    ;;
  auto)
    if [[ "${BUILD_SIF}" == "1" ]] && image_looks_remote "${IMAGE}"; then
      resolved_image_source="public"
      echo "[gg_container_build] IMAGE points at a registry image; skipping local Docker build."
      echo "[gg_container_build] using Apptainer/Singularity to pull docker://${IMAGE}:${TAG} directly."
    elif docker_buildx_available; then
      echo "[gg_container_build] step 1/2: docker buildx build"
      (
        cd "${repo_root}"
        IMAGE="${IMAGE}" TAG="${TAG}" PLATFORMS="${PLATFORMS}" MODE="${MODE}" \
          bash "${buildx_script}"
      )
      did_buildx=1
      resolved_image_source="local"
      if [[ "${MODE}" == "load" ]]; then
        sif_source="docker-daemon"
      fi
    elif [[ "${BUILD_SIF}" == "1" ]]; then
      if image_looks_remote "${IMAGE}"; then
        resolved_image_source="public"
        echo "[gg_container_build] docker buildx not found; skipping local Docker build."
        echo "[gg_container_build] using Apptainer/Singularity to pull docker://${IMAGE}:${TAG} directly."
      elif [[ "${image_was_default}" == "1" ]]; then
        resolved_image_source="public"
        echo "[gg_container_build] docker buildx not found; default local image cannot be built on this host."
        IMAGE="${FALLBACK_REMOTE_IMAGE}"
        if [[ "${tag_was_default}" == "1" || "${TAG}" == "${LOCAL_TAG_DEFAULT}" ]]; then
          TAG="${FALLBACK_REMOTE_TAG}"
        fi
        echo "[gg_container_build] falling back to published image ${IMAGE}:${TAG}."
      else
        echo "[gg_container_build] docker buildx not found, so local Dockerfile builds are unavailable."
        echo "[gg_container_build] Either install docker buildx, or use a registry image such as:"
        echo "[gg_container_build]   IMAGE_SOURCE=public IMAGE=${FALLBACK_REMOTE_IMAGE} TAG=${FALLBACK_REMOTE_TAG} bash ./gg_container_build_entrypoint.sh"
        exit 1
      fi
    else
      echo "[gg_container_build] docker buildx not found and BUILD_SIF=0."
      echo "[gg_container_build] Nothing can be built on this host without Docker."
      exit 1
    fi
    ;;
esac

echo "[gg_container_build] resolved_image_source=${resolved_image_source}"

if [[ "${BUILD_SIF}" == "1" ]]; then
  if [[ "${did_native_local_build}" == "1" ]]; then
    echo "[gg_container_build] done: ${OUT}"
    exit 0
  fi
  resolved_engine="$(resolve_container_engine)"
  echo "[gg_container_build] engine=${resolved_engine}"
  echo "[gg_container_build] sif_source=${sif_source}"
  mkdir -p "$(dirname "${OUT}")"
  if [[ "${did_buildx}" == "1" ]]; then
    echo "[gg_container_build] step 2/2: convert docker image to SIF"
  else
    echo "[gg_container_build] step 1/1: build SIF from registry image"
  fi
  (
    cd "${repo_root}"
    IMAGE="${IMAGE}" TAG="${TAG}" ENGINE="${resolved_engine}" OUT="${OUT}" SOURCE="${sif_source}" \
      bash "${apptainer_script}"
  )
  echo "[gg_container_build] done: ${OUT}"
else
  echo "[gg_container_build] BUILD_SIF=0, skipped SIF conversion."
fi

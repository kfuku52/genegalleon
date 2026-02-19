#!/usr/bin/env bash
set -euo pipefail

# Example:
#   IMAGE=ghcr.io/your-org/genegalleon TAG=20260211 MODE=push ./container/buildx.sh
#   IMAGE=local/genegalleon TAG=dev PLATFORMS=linux/arm64 MODE=load ./container/buildx.sh
#   KFU52_AMALGKIT_BRANCH_CANDIDATES=master,kfdevel MODE=push ./container/buildx.sh

IMAGE=${IMAGE:-ghcr.io/example/genegalleon}
TAG=${TAG:-dev}
PLATFORMS=${PLATFORMS:-linux/amd64,linux/arm64}
MODE=${MODE:-push} # push | load
NOTUNG_DOWNLOAD_PAGE=${NOTUNG_DOWNLOAD_PAGE:-https://amberjack.compbio.cs.cmu.edu/Notung/download29.html}
KFU52_REPO_REF=${KFU52_REPO_REF:-master}
KFU52_AMALGKIT_AUTO_SELECT_REF=${KFU52_AMALGKIT_AUTO_SELECT_REF:-1}
KFU52_AMALGKIT_BRANCH_CANDIDATES=${KFU52_AMALGKIT_BRANCH_CANDIDATES:-master,kfdevel}
KFU52_AMALGKIT_REPO_REF=${KFU52_AMALGKIT_REPO_REF:-}
KFTOOLS_REPO_URL=${KFTOOLS_REPO_URL:-https://github.com/kfuku52/kftools.git}
RKFTOOLS_REPO_URL=${RKFTOOLS_REPO_URL:-https://github.com/kfuku52/rkftools.git}
KFTOOLS_REPO_REF=${KFTOOLS_REPO_REF:-${KFU52_REPO_REF}}
RKFTOOLS_REPO_REF=${RKFTOOLS_REPO_REF:-${KFU52_REPO_REF}}
CACHE_DIR=${CACHE_DIR:-.buildx-cache}
USE_LOCAL_CACHE=${USE_LOCAL_CACHE:-1}
CACHE_FROM=${CACHE_FROM:-}
CACHE_TO=${CACHE_TO:-}
FORCE_LOAD_IN_CI=${FORCE_LOAD_IN_CI:-0}

if [[ "${MODE}" != "push" && "${MODE}" != "load" ]]; then
  echo "MODE must be one of: push, load"
  exit 1
fi

if [[ "${MODE}" == "load" ]]; then
  echo "[buildx] MODE=load exports a local image tarball and is typically slower than MODE=push."
  if [[ -n "${CI:-}" && "${FORCE_LOAD_IN_CI}" != "1" ]]; then
    echo "[buildx] CI detected; switching MODE=push for faster non-interactive builds."
    echo "[buildx] Set FORCE_LOAD_IN_CI=1 to keep MODE=load in CI."
    MODE="push"
  fi
fi

if [[ "${MODE}" == "load" && "${PLATFORMS}" == *","* ]]; then
  echo "MODE=load supports a single platform only. Set PLATFORMS=linux/amd64 or linux/arm64."
  exit 1
fi

if ! docker buildx ls >/dev/null 2>&1; then
  echo "docker buildx is required."
  exit 1
fi

# Ensure conda env coverage stays consistent with pipeline script usage.
container/scripts/check_env_coverage.sh .

output_flag="--push"
if [[ "${MODE}" == "load" ]]; then
  output_flag="--load"
fi

build_date=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
vcs_ref=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")

cache_args=()
cache_dir_new=""
buildx_driver=$(docker buildx inspect --bootstrap 2>/dev/null | awk -F': *' '/^Driver:/ {print $2; exit}' | xargs)
if [[ -n "${CACHE_FROM}" || -n "${CACHE_TO}" ]]; then
  if [[ -n "${CACHE_FROM}" ]]; then
    cache_args+=(--cache-from "${CACHE_FROM}")
  fi
  if [[ -n "${CACHE_TO}" ]]; then
    cache_args+=(--cache-to "${CACHE_TO}")
  fi
elif [[ "${USE_LOCAL_CACHE}" == "1" ]]; then
  if [[ "${buildx_driver}" == "docker" ]]; then
    echo "[buildx] Driver '${buildx_driver}' does not support local cache export; continuing without --cache-to/--cache-from."
    echo "[buildx] For external cache, use a docker-container builder or set CACHE_FROM/CACHE_TO explicitly."
  else
    cache_dir_new="${CACHE_DIR}.new"
    mkdir -p "${CACHE_DIR}"
    rm -rf "${cache_dir_new}"
    cache_args+=(--cache-from "type=local,src=${CACHE_DIR}")
    cache_args+=(--cache-to "type=local,dest=${cache_dir_new},mode=max")
  fi
fi

if [[ ${#cache_args[@]} -gt 0 ]]; then
  docker buildx build \
    --platform "${PLATFORMS}" \
    --file container/Dockerfile \
    --build-arg BUILD_DATE="${build_date}" \
    --build-arg VCS_REF="${vcs_ref}" \
    --build-arg NOTUNG_DOWNLOAD_PAGE="${NOTUNG_DOWNLOAD_PAGE}" \
    --build-arg KFU52_REPO_REF="${KFU52_REPO_REF}" \
    --build-arg KFU52_AMALGKIT_AUTO_SELECT_REF="${KFU52_AMALGKIT_AUTO_SELECT_REF}" \
    --build-arg KFU52_AMALGKIT_BRANCH_CANDIDATES="${KFU52_AMALGKIT_BRANCH_CANDIDATES}" \
    --build-arg KFU52_AMALGKIT_REPO_REF="${KFU52_AMALGKIT_REPO_REF}" \
    --build-arg KFTOOLS_REPO_URL="${KFTOOLS_REPO_URL}" \
    --build-arg RKFTOOLS_REPO_URL="${RKFTOOLS_REPO_URL}" \
    --build-arg KFTOOLS_REPO_REF="${KFTOOLS_REPO_REF}" \
    --build-arg RKFTOOLS_REPO_REF="${RKFTOOLS_REPO_REF}" \
    --tag "${IMAGE}:${TAG}" \
    "${cache_args[@]}" \
    ${output_flag} \
    .
else
  docker buildx build \
    --platform "${PLATFORMS}" \
    --file container/Dockerfile \
    --build-arg BUILD_DATE="${build_date}" \
    --build-arg VCS_REF="${vcs_ref}" \
    --build-arg NOTUNG_DOWNLOAD_PAGE="${NOTUNG_DOWNLOAD_PAGE}" \
    --build-arg KFU52_REPO_REF="${KFU52_REPO_REF}" \
    --build-arg KFU52_AMALGKIT_AUTO_SELECT_REF="${KFU52_AMALGKIT_AUTO_SELECT_REF}" \
    --build-arg KFU52_AMALGKIT_BRANCH_CANDIDATES="${KFU52_AMALGKIT_BRANCH_CANDIDATES}" \
    --build-arg KFU52_AMALGKIT_REPO_REF="${KFU52_AMALGKIT_REPO_REF}" \
    --build-arg KFTOOLS_REPO_URL="${KFTOOLS_REPO_URL}" \
    --build-arg RKFTOOLS_REPO_URL="${RKFTOOLS_REPO_URL}" \
    --build-arg KFTOOLS_REPO_REF="${KFTOOLS_REPO_REF}" \
    --build-arg RKFTOOLS_REPO_REF="${RKFTOOLS_REPO_REF}" \
    --tag "${IMAGE}:${TAG}" \
    ${output_flag} \
    .
fi

if [[ -n "${cache_dir_new}" && -d "${cache_dir_new}" ]]; then
  rm -rf "${CACHE_DIR}"
  mv "${cache_dir_new}" "${CACHE_DIR}"
fi

echo "Built ${IMAGE}:${TAG} for ${PLATFORMS}"

#!/usr/bin/env bash
set -euo pipefail

# Example:
#   IMAGE=ghcr.io/your-org/genegalleon TAG=20260211 MODE=push ./container/buildx.sh
#   IMAGE=local/genegalleon TAG=dev PLATFORMS=linux/arm64 MODE=load ./container/buildx.sh
#   KFU52_AMALGKIT_BRANCH_CANDIDATES=master,kfdevel,devel MODE=push ./container/buildx.sh

IMAGE=${IMAGE:-ghcr.io/example/genegalleon}
TAG=${TAG:-dev}
PLATFORMS=${PLATFORMS:-linux/amd64,linux/arm64}
MODE=${MODE:-push} # push | load
NOTUNG_DOWNLOAD_PAGE=${NOTUNG_DOWNLOAD_PAGE:-https://amberjack.compbio.cs.cmu.edu/Notung/Notung-2.9.1.5.zip}
NOTUNG_DOWNLOAD_HOST_IP=${NOTUNG_DOWNLOAD_HOST_IP:-128.2.205.60}
KFU52_REPO_REF=${KFU52_REPO_REF:-master}
KFU52_AMALGKIT_AUTO_SELECT_REF=${KFU52_AMALGKIT_AUTO_SELECT_REF:-1}
KFU52_AMALGKIT_BRANCH_CANDIDATES=${KFU52_AMALGKIT_BRANCH_CANDIDATES:-master,kfdevel,devel}
KFU52_AMALGKIT_REPO_REF=${KFU52_AMALGKIT_REPO_REF:-}
KFU52_AMALGKIT_REPO_SHA=${KFU52_AMALGKIT_REPO_SHA:-1ce96e3148cace8f96a67fdb46df99400a9e3f19}
KFU52_CDSKIT_REPO_SHA=${KFU52_CDSKIT_REPO_SHA:-9cd69c01c9736dec2a1eea63e63d994fbc4c812c}
KFU52_CSUBST_REPO_SHA=${KFU52_CSUBST_REPO_SHA:-613771ce6d19d583a7adee0bcc63abed312cdad4}
KFU52_NWKIT_REPO_SHA=${KFU52_NWKIT_REPO_SHA:-106a8cfd97cad993772aa4c28e7d9a87d667ff86}
BUSCO_REPO_URL=${BUSCO_REPO_URL:-https://gitlab.com/ezlab/busco.git}
BUSCO_MIRROR_REPO_URL=${BUSCO_MIRROR_REPO_URL:-}
BUSCO_REPO_REF=${BUSCO_REPO_REF:-6.0.0}
BUSCO_REPO_SHA=${BUSCO_REPO_SHA:-6278721a1916f6da310e03ec9674099028c927a4}
PAML_REPO_URL=${PAML_REPO_URL:-https://github.com/iqtree/paml.git}
PAML_REPO_REF=${PAML_REPO_REF:-master}
PAML_REPO_SHA=${PAML_REPO_SHA:-8daeead6b55523f375d9ac56dcfac38373ef8a2e}
KFL1OU_REPO_URL=${KFL1OU_REPO_URL:-https://github.com/kfuku52/kfl1ou.git}
KFL1OU_REPO_REF=${KFL1OU_REPO_REF:-}
KFL1OU_REPO_SHA=${KFL1OU_REPO_SHA:-1bf3028f204a6d58e697f58461c82ecfc7c29802}
KFTOOLS_REPO_URL=${KFTOOLS_REPO_URL:-https://github.com/kfuku52/kftools.git}
RKFTOOLS_REPO_URL=${RKFTOOLS_REPO_URL:-https://github.com/kfuku52/rkftools.git}
RADTE_REPO_URL=${RADTE_REPO_URL:-https://github.com/kfuku52/RADTE.git}
KFTOOLS_REPO_REF=${KFTOOLS_REPO_REF:-${KFU52_REPO_REF}}
RKFTOOLS_REPO_REF=${RKFTOOLS_REPO_REF:-${KFU52_REPO_REF}}
RADTE_REPO_REF=${RADTE_REPO_REF:-}
KFTOOLS_REPO_SHA=${KFTOOLS_REPO_SHA:-4918fed6146b9cef1df66b5ce33de70b74454547}
RKFTOOLS_REPO_SHA=${RKFTOOLS_REPO_SHA:-cf16e570300ec32909d8cd458119712d40bcf06f}
RADTE_REPO_SHA=${RADTE_REPO_SHA:-873c4acb22d3decedf417bb95e3d292abccbb386}
TESTNH_TARBALL_SHA256=${TESTNH_TARBALL_SHA256:-598337183d2cec9c61cd364fab255a270062844b0ba5172913f7cf97512c43e2}
CAFE5_TARBALL_SHA256=${CAFE5_TARBALL_SHA256:-71871bdc74c2ffc7c1c0f4500f4742f2ff46a15cfaba78dc179d21bb1ba67ba8}
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
gg_version="$(sed -n '1p' VERSION 2>/dev/null | tr -d '\r' || true)"
if [[ -z "${gg_version}" ]]; then
  gg_version="unknown"
fi

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
    rm -rf -- "${cache_dir_new}"
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
    --build-arg GG_VERSION="${gg_version}" \
    --build-arg NOTUNG_DOWNLOAD_PAGE="${NOTUNG_DOWNLOAD_PAGE}" \
    --build-arg NOTUNG_DOWNLOAD_HOST_IP="${NOTUNG_DOWNLOAD_HOST_IP}" \
    --build-arg KFU52_REPO_REF="${KFU52_REPO_REF}" \
    --build-arg KFU52_AMALGKIT_AUTO_SELECT_REF="${KFU52_AMALGKIT_AUTO_SELECT_REF}" \
    --build-arg KFU52_AMALGKIT_BRANCH_CANDIDATES="${KFU52_AMALGKIT_BRANCH_CANDIDATES}" \
    --build-arg KFU52_AMALGKIT_REPO_REF="${KFU52_AMALGKIT_REPO_REF}" \
    --build-arg KFU52_AMALGKIT_REPO_SHA="${KFU52_AMALGKIT_REPO_SHA}" \
    --build-arg KFU52_CDSKIT_REPO_SHA="${KFU52_CDSKIT_REPO_SHA}" \
    --build-arg KFU52_CSUBST_REPO_SHA="${KFU52_CSUBST_REPO_SHA}" \
    --build-arg KFU52_NWKIT_REPO_SHA="${KFU52_NWKIT_REPO_SHA}" \
    --build-arg BUSCO_REPO_URL="${BUSCO_REPO_URL}" \
    --build-arg BUSCO_MIRROR_REPO_URL="${BUSCO_MIRROR_REPO_URL}" \
    --build-arg BUSCO_REPO_REF="${BUSCO_REPO_REF}" \
    --build-arg BUSCO_REPO_SHA="${BUSCO_REPO_SHA}" \
    --build-arg PAML_REPO_URL="${PAML_REPO_URL}" \
    --build-arg PAML_REPO_REF="${PAML_REPO_REF}" \
    --build-arg PAML_REPO_SHA="${PAML_REPO_SHA}" \
    --build-arg KFL1OU_REPO_URL="${KFL1OU_REPO_URL}" \
    --build-arg KFL1OU_REPO_REF="${KFL1OU_REPO_REF}" \
    --build-arg KFL1OU_REPO_SHA="${KFL1OU_REPO_SHA}" \
    --build-arg KFTOOLS_REPO_URL="${KFTOOLS_REPO_URL}" \
    --build-arg RKFTOOLS_REPO_URL="${RKFTOOLS_REPO_URL}" \
    --build-arg RADTE_REPO_URL="${RADTE_REPO_URL}" \
    --build-arg KFTOOLS_REPO_REF="${KFTOOLS_REPO_REF}" \
    --build-arg RKFTOOLS_REPO_REF="${RKFTOOLS_REPO_REF}" \
    --build-arg RADTE_REPO_REF="${RADTE_REPO_REF}" \
    --build-arg KFTOOLS_REPO_SHA="${KFTOOLS_REPO_SHA}" \
    --build-arg RKFTOOLS_REPO_SHA="${RKFTOOLS_REPO_SHA}" \
    --build-arg RADTE_REPO_SHA="${RADTE_REPO_SHA}" \
    --build-arg TESTNH_TARBALL_SHA256="${TESTNH_TARBALL_SHA256}" \
    --build-arg CAFE5_TARBALL_SHA256="${CAFE5_TARBALL_SHA256}" \
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
    --build-arg GG_VERSION="${gg_version}" \
    --build-arg NOTUNG_DOWNLOAD_PAGE="${NOTUNG_DOWNLOAD_PAGE}" \
    --build-arg NOTUNG_DOWNLOAD_HOST_IP="${NOTUNG_DOWNLOAD_HOST_IP}" \
    --build-arg KFU52_REPO_REF="${KFU52_REPO_REF}" \
    --build-arg KFU52_AMALGKIT_AUTO_SELECT_REF="${KFU52_AMALGKIT_AUTO_SELECT_REF}" \
    --build-arg KFU52_AMALGKIT_BRANCH_CANDIDATES="${KFU52_AMALGKIT_BRANCH_CANDIDATES}" \
    --build-arg KFU52_AMALGKIT_REPO_REF="${KFU52_AMALGKIT_REPO_REF}" \
    --build-arg KFU52_AMALGKIT_REPO_SHA="${KFU52_AMALGKIT_REPO_SHA}" \
    --build-arg KFU52_CDSKIT_REPO_SHA="${KFU52_CDSKIT_REPO_SHA}" \
    --build-arg KFU52_CSUBST_REPO_SHA="${KFU52_CSUBST_REPO_SHA}" \
    --build-arg KFU52_NWKIT_REPO_SHA="${KFU52_NWKIT_REPO_SHA}" \
    --build-arg BUSCO_REPO_URL="${BUSCO_REPO_URL}" \
    --build-arg BUSCO_MIRROR_REPO_URL="${BUSCO_MIRROR_REPO_URL}" \
    --build-arg BUSCO_REPO_REF="${BUSCO_REPO_REF}" \
    --build-arg BUSCO_REPO_SHA="${BUSCO_REPO_SHA}" \
    --build-arg PAML_REPO_URL="${PAML_REPO_URL}" \
    --build-arg PAML_REPO_REF="${PAML_REPO_REF}" \
    --build-arg PAML_REPO_SHA="${PAML_REPO_SHA}" \
    --build-arg KFL1OU_REPO_URL="${KFL1OU_REPO_URL}" \
    --build-arg KFL1OU_REPO_REF="${KFL1OU_REPO_REF}" \
    --build-arg KFL1OU_REPO_SHA="${KFL1OU_REPO_SHA}" \
    --build-arg KFTOOLS_REPO_URL="${KFTOOLS_REPO_URL}" \
    --build-arg RKFTOOLS_REPO_URL="${RKFTOOLS_REPO_URL}" \
    --build-arg RADTE_REPO_URL="${RADTE_REPO_URL}" \
    --build-arg KFTOOLS_REPO_REF="${KFTOOLS_REPO_REF}" \
    --build-arg RKFTOOLS_REPO_REF="${RKFTOOLS_REPO_REF}" \
    --build-arg RADTE_REPO_REF="${RADTE_REPO_REF}" \
    --build-arg KFTOOLS_REPO_SHA="${KFTOOLS_REPO_SHA}" \
    --build-arg RKFTOOLS_REPO_SHA="${RKFTOOLS_REPO_SHA}" \
    --build-arg RADTE_REPO_SHA="${RADTE_REPO_SHA}" \
    --build-arg TESTNH_TARBALL_SHA256="${TESTNH_TARBALL_SHA256}" \
    --build-arg CAFE5_TARBALL_SHA256="${CAFE5_TARBALL_SHA256}" \
    --tag "${IMAGE}:${TAG}" \
    ${output_flag} \
    .
fi

if [[ -n "${cache_dir_new}" && -d "${cache_dir_new}" ]]; then
  rm -rf -- "${CACHE_DIR}"
  mv "${cache_dir_new}" "${CACHE_DIR}"
fi

echo "Built ${IMAGE}:${TAG} for ${PLATFORMS}"

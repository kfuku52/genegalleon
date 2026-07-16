#!/usr/bin/env bash
set -euo pipefail

# Example:
#   IMAGE=ghcr.io/your-org/genegalleon TAG=20260211 MODE=push ./container/buildx.sh
#   IMAGE=local/genegalleon TAG=dev PLATFORMS=linux/arm64 MODE=load ./container/buildx.sh
#   KFU52_AMALGKIT_REPO_REF=kfdevel MODE=push ./container/buildx.sh

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

IMAGE=${IMAGE:-ghcr.io/example/genegalleon}
TAG=${TAG:-dev}
PLATFORMS=${PLATFORMS:-linux/amd64,linux/arm64}
MODE=${MODE:-push} # push | load
NOTUNG_DOWNLOAD_PAGE=${NOTUNG_DOWNLOAD_PAGE:-https://amberjack.compbio.cs.cmu.edu/Notung/Notung-2.9.1.5.zip}
NOTUNG_DOWNLOAD_HOST_IP=${NOTUNG_DOWNLOAD_HOST_IP:-128.2.205.60}
KFU52_REPO_REF=${KFU52_REPO_REF:-master}
KFU52_AMALGKIT_AUTO_SELECT_REF=${KFU52_AMALGKIT_AUTO_SELECT_REF:-0}
KFU52_AMALGKIT_BRANCH_CANDIDATES=${KFU52_AMALGKIT_BRANCH_CANDIDATES:-master,kfdevel,devel}
KFU52_AMALGKIT_REPO_REF=${KFU52_AMALGKIT_REPO_REF-master}
KFU52_AMALGKIT_REPO_SHA=${KFU52_AMALGKIT_REPO_SHA:-}
KFU52_CDSKIT_REPO_SHA=${KFU52_CDSKIT_REPO_SHA:-}
KFU52_CSUBST_REPO_REF=${KFU52_CSUBST_REPO_REF:-master}
KFU52_CSUBST_REPO_SHA=${KFU52_CSUBST_REPO_SHA:-}
KFU52_NWKIT_REPO_SHA=${KFU52_NWKIT_REPO_SHA:-}
BUSCO_REPO_URL=${BUSCO_REPO_URL:-https://gitlab.com/ezlab/busco.git}
BUSCO_MIRROR_REPO_URL=${BUSCO_MIRROR_REPO_URL:-}
BUSCO_REPO_REF=${BUSCO_REPO_REF:-6.0.0}
BUSCO_REPO_SHA=${BUSCO_REPO_SHA:-6278721a1916f6da310e03ec9674099028c927a4}
PAML_REPO_URL=${PAML_REPO_URL:-https://github.com/iqtree/paml.git}
PAML_REPO_REF=${PAML_REPO_REF:-master}
PAML_REPO_SHA=${PAML_REPO_SHA:-8daeead6b55523f375d9ac56dcfac38373ef8a2e}
KFL1OU_REPO_URL=${KFL1OU_REPO_URL:-https://github.com/kfuku52/kfl1ou.git}
KFL1OU_REPO_REF=${KFL1OU_REPO_REF:-${KFU52_REPO_REF}}
KFL1OU_REPO_SHA=${KFL1OU_REPO_SHA:-}
KFTOOLS_REPO_URL=${KFTOOLS_REPO_URL:-https://github.com/kfuku52/kftools.git}
RKFTOOLS_REPO_URL=${RKFTOOLS_REPO_URL:-https://github.com/kfuku52/rkftools.git}
RADTE_REPO_URL=${RADTE_REPO_URL:-https://github.com/kfuku52/RADTE.git}
KFTOOLS_REPO_REF=${KFTOOLS_REPO_REF:-${KFU52_REPO_REF}}
RKFTOOLS_REPO_REF=${RKFTOOLS_REPO_REF:-${KFU52_REPO_REF}}
RADTE_REPO_REF=${RADTE_REPO_REF:-${KFU52_REPO_REF}}
KFTOOLS_REPO_SHA=${KFTOOLS_REPO_SHA:-}
RKFTOOLS_REPO_SHA=${RKFTOOLS_REPO_SHA:-}
RADTE_REPO_SHA=${RADTE_REPO_SHA:-}
TESTNH_TARBALL_SHA256=${TESTNH_TARBALL_SHA256:-598337183d2cec9c61cd364fab255a270062844b0ba5172913f7cf97512c43e2}
CAFE5_TARBALL_SHA256=${CAFE5_TARBALL_SHA256:-71871bdc74c2ffc7c1c0f4500f4742f2ff46a15cfaba78dc179d21bb1ba67ba8}
CACHE_DIR=${CACHE_DIR:-.buildx-cache}
# BuildKit already keeps an internal cache. Exporting this image's full cache
# to a local directory costs several gigabytes and is therefore opt-in.
USE_LOCAL_CACHE=${USE_LOCAL_CACHE:-0}
CACHE_FROM=${CACHE_FROM:-}
CACHE_TO=${CACHE_TO:-}
FORCE_LOAD_IN_CI=${FORCE_LOAD_IN_CI:-0}
SKIP_UNCHANGED_LOAD=${SKIP_UNCHANGED_LOAD:-1}
PRINT_BUILD_INPUT_HASH=${PRINT_BUILD_INPUT_HASH:-0}

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

resolve_source_sha() {
  local sha_var="$1"
  local repo_url="$2"
  local repo_ref="$3"
  local source_name="$4"
  local resolved_sha="${!sha_var}"

  if [[ -z "${resolved_sha}" ]]; then
    resolved_sha="$(bash container/scripts/resolve_git_branch_sha.sh "${repo_url}" "${repo_ref}")"
    printf -v "${sha_var}" '%s' "${resolved_sha}"
  fi
  echo "[buildx] Resolved ${source_name} ${repo_ref}: ${resolved_sha}"
}

resolve_source_sha KFU52_AMALGKIT_REPO_SHA https://github.com/kfuku52/amalgkit.git "${KFU52_AMALGKIT_REPO_REF}" amalgkit
resolve_source_sha KFU52_CDSKIT_REPO_SHA https://github.com/kfuku52/cdskit.git "${KFU52_REPO_REF}" cdskit
resolve_source_sha KFU52_CSUBST_REPO_SHA https://github.com/kfuku52/csubst.git "${KFU52_CSUBST_REPO_REF}" csubst
resolve_source_sha KFU52_NWKIT_REPO_SHA https://github.com/kfuku52/nwkit.git "${KFU52_REPO_REF}" nwkit
resolve_source_sha KFL1OU_REPO_SHA "${KFL1OU_REPO_URL}" "${KFL1OU_REPO_REF}" kfl1ou
resolve_source_sha KFTOOLS_REPO_SHA "${KFTOOLS_REPO_URL}" "${KFTOOLS_REPO_REF}" kftools
resolve_source_sha RKFTOOLS_REPO_SHA "${RKFTOOLS_REPO_URL}" "${RKFTOOLS_REPO_REF}" rkftools
resolve_source_sha RADTE_REPO_SHA "${RADTE_REPO_URL}" "${RADTE_REPO_REF}" RADTE

output_flag="--push"
if [[ "${MODE}" == "load" ]]; then
  output_flag="--load"
fi

vcs_ref=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")
gg_version="$(sed -n '1p' VERSION 2>/dev/null | tr -d '\r' || true)"
if [[ -z "${gg_version}" ]]; then
  gg_version="unknown"
fi

sha256_stream() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum | awk '{print $1}'
  else
    shasum -a 256 | awk '{print $1}'
  fi
}

sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

context_digest="$({
  printf '%s\n' .dockerignore container/Dockerfile container/pip-compatibility.requirements.txt
  find container/env container/spec container/testdata container/scripts -type f -print
  find workflow/support/treevis -type f -print
} | LC_ALL=C sort | while IFS= read -r path; do
  printf '%s\t%s\n' "${path}" "$(sha256_file "${path}")"
done | sha256_stream)"

build_input_hash="$({
  printf '%s\n' \
    "context=${context_digest}" \
    "platforms=${PLATFORMS}" \
    "vcs_ref=${vcs_ref}" \
    "gg_version=${gg_version}" \
    "notung_page=${NOTUNG_DOWNLOAD_PAGE}" \
    "notung_host=${NOTUNG_DOWNLOAD_HOST_IP}" \
    "kfu52_ref=${KFU52_REPO_REF}" \
    "amalgkit_auto=${KFU52_AMALGKIT_AUTO_SELECT_REF}" \
    "amalgkit_candidates=${KFU52_AMALGKIT_BRANCH_CANDIDATES}" \
    "amalgkit_ref=${KFU52_AMALGKIT_REPO_REF}" \
    "amalgkit_sha=${KFU52_AMALGKIT_REPO_SHA}" \
    "cdskit_sha=${KFU52_CDSKIT_REPO_SHA}" \
    "csubst_ref=${KFU52_CSUBST_REPO_REF}" \
    "csubst_sha=${KFU52_CSUBST_REPO_SHA}" \
    "nwkit_sha=${KFU52_NWKIT_REPO_SHA}" \
    "busco_url=${BUSCO_REPO_URL}" \
    "busco_mirror=${BUSCO_MIRROR_REPO_URL}" \
    "busco_ref=${BUSCO_REPO_REF}" \
    "busco_sha=${BUSCO_REPO_SHA}" \
    "paml_url=${PAML_REPO_URL}" \
    "paml_ref=${PAML_REPO_REF}" \
    "paml_sha=${PAML_REPO_SHA}" \
    "kfl1ou_url=${KFL1OU_REPO_URL}" \
    "kfl1ou_ref=${KFL1OU_REPO_REF}" \
    "kfl1ou_sha=${KFL1OU_REPO_SHA}" \
    "kftools_url=${KFTOOLS_REPO_URL}" \
    "kftools_ref=${KFTOOLS_REPO_REF}" \
    "kftools_sha=${KFTOOLS_REPO_SHA}" \
    "rkftools_url=${RKFTOOLS_REPO_URL}" \
    "rkftools_ref=${RKFTOOLS_REPO_REF}" \
    "rkftools_sha=${RKFTOOLS_REPO_SHA}" \
    "radte_url=${RADTE_REPO_URL}" \
    "radte_ref=${RADTE_REPO_REF}" \
    "radte_sha=${RADTE_REPO_SHA}" \
    "testnh_sha=${TESTNH_TARBALL_SHA256}" \
    "cafe5_sha=${CAFE5_TARBALL_SHA256}"
} | sha256_stream)"

if [[ "${PRINT_BUILD_INPUT_HASH}" == "1" ]]; then
  printf '%s\n' "${build_input_hash}"
  exit 0
fi

if [[ "${MODE}" == "load" && "${SKIP_UNCHANGED_LOAD}" == "1" ]]; then
  target_arch="${PLATFORMS#linux/}"
  image_state="$(docker image inspect \
    --format '{{.Architecture}} {{index .Config.Labels "io.genegalleon.build-input"}}' \
    "${IMAGE}:${TAG}" 2>/dev/null || true)"
  image_arch="${image_state%% *}"
  image_hash="${image_state#* }"
  if [[ "${image_state}" == *" "* && "${image_arch}" == "${target_arch}" && "${image_hash}" == "${build_input_hash}" ]]; then
    echo "[buildx] ${IMAGE}:${TAG} already matches build input ${build_input_hash}; skipping unchanged MODE=load build."
    echo "[buildx] Set SKIP_UNCHANGED_LOAD=0 to force BuildKit execution."
    exit 0
  fi
fi

build_date=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

cache_args=()
cache_dir_new=""
buildx_driver=""
if buildx_inspect="$(docker buildx inspect --bootstrap 2>/dev/null)"; then
  buildx_driver="$(printf '%s\n' "${buildx_inspect}" | awk -F': *' '/^Driver:/ {print $2; exit}' | xargs)"
else
  echo "[buildx] Could not inspect the builder; continuing without local cache configuration."
fi
if [[ -n "${CACHE_FROM}" || -n "${CACHE_TO}" ]]; then
  if [[ -n "${CACHE_FROM}" ]]; then
    cache_args+=(--cache-from "${CACHE_FROM}")
  fi
  if [[ -n "${CACHE_TO}" ]]; then
    cache_args+=(--cache-to "${CACHE_TO}")
  fi
elif [[ "${USE_LOCAL_CACHE}" == "1" ]]; then
  if [[ -z "${buildx_driver}" ]]; then
    echo "[buildx] Builder driver is unknown; continuing without --cache-to/--cache-from."
  elif [[ "${buildx_driver}" == "docker" ]]; then
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

run_build() {
  docker buildx build \
    --platform "${PLATFORMS}" \
    --file container/Dockerfile \
    --build-arg BUILD_DATE="${build_date}" \
    --build-arg VCS_REF="${vcs_ref}" \
    --build-arg GG_VERSION="${gg_version}" \
    --build-arg BUILD_INPUT_HASH="${build_input_hash}" \
    --build-arg NOTUNG_DOWNLOAD_PAGE="${NOTUNG_DOWNLOAD_PAGE}" \
    --build-arg NOTUNG_DOWNLOAD_HOST_IP="${NOTUNG_DOWNLOAD_HOST_IP}" \
    --build-arg KFU52_REPO_REF="${KFU52_REPO_REF}" \
    --build-arg KFU52_AMALGKIT_AUTO_SELECT_REF="${KFU52_AMALGKIT_AUTO_SELECT_REF}" \
    --build-arg KFU52_AMALGKIT_BRANCH_CANDIDATES="${KFU52_AMALGKIT_BRANCH_CANDIDATES}" \
    --build-arg KFU52_AMALGKIT_REPO_REF="${KFU52_AMALGKIT_REPO_REF}" \
    --build-arg KFU52_AMALGKIT_REPO_SHA="${KFU52_AMALGKIT_REPO_SHA}" \
    --build-arg KFU52_CDSKIT_REPO_SHA="${KFU52_CDSKIT_REPO_SHA}" \
    --build-arg KFU52_CSUBST_REPO_REF="${KFU52_CSUBST_REPO_REF}" \
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
    "$@" \
    ${output_flag} \
    .
}

if [[ ${#cache_args[@]} -gt 0 ]]; then
  run_build "${cache_args[@]}"
else
  run_build
fi

if [[ -n "${cache_dir_new}" && -d "${cache_dir_new}" ]]; then
  rm -rf -- "${CACHE_DIR}"
  mv "${cache_dir_new}" "${CACHE_DIR}"
fi

echo "Built ${IMAGE}:${TAG} for ${PLATFORMS}"

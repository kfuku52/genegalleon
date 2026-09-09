#!/usr/bin/env bash
# Build unmodified official IQ-TREE 3 with the local NWKIT consumer.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  echo "Usage: $0 NWKIT_SOURCE" >&2
  echo "Optional: BASE_IMAGE, IMAGE, GG_BUILD_JOBS, IQTREE_REPO_SHA" >&2
  exit 2
fi
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/source_branches.env"
nwkit_source="$(cd "$1" && pwd)"
if [[ ! -f "${nwkit_source}/nwkit/radte_iqtree.py" ]]; then
  echo "NWKIT_SOURCE must contain the IQ-TREE 3 dating adapter." >&2
  exit 1
fi
build_directory=$(mktemp -d)
trap 'rm -rf "${build_directory}"' EXIT
git clone --depth 1 --branch "${GG_SOURCE_IQTREE_REPO_REF}" \
  https://github.com/iqtree/iqtree3.git "${build_directory}/iqtree"
if [[ -n "${IQTREE_REPO_SHA:-}" ]]; then
  git -C "${build_directory}/iqtree" fetch --depth 1 origin "${IQTREE_REPO_SHA}"
  git -C "${build_directory}/iqtree" checkout --detach FETCH_HEAD
fi
git -C "${build_directory}/iqtree" submodule update --init --recursive --depth 1
revision=$(git -C "${build_directory}/iqtree" rev-parse HEAD)
echo "Building official IQ-TREE 3 source ${revision}"
docker build \
  --build-context "iqtree_source=${build_directory}/iqtree" \
  --build-context "nwkit_source=${nwkit_source}" \
  --build-arg "BASE_IMAGE=${BASE_IMAGE:-local/genegalleon:dev}" \
  --build-arg "GG_BUILD_JOBS=${GG_BUILD_JOBS:-2}" \
  --build-arg "IQTREE_SOURCE_REVISION=${revision}" \
  -f "${script_dir}/Dockerfile.iqtree3" \
  -t "${IMAGE:-local/genegalleon:iqtree3-dev}" \
  "${script_dir}"

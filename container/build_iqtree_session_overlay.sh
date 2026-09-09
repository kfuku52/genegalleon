#!/usr/bin/env bash
# Build the locally developed IQ-TREE/NWKIT protocol into a GeneGalleon runtime.
set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 IQTREE_SOURCE NWKIT_SOURCE" >&2
  echo "Optional: BASE_IMAGE, IMAGE, GG_BUILD_JOBS" >&2
  exit 2
fi
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
iqtree_source="$(cd "$1" && pwd)"
nwkit_source="$(cd "$2" && pwd)"
if [[ ! -f "${iqtree_source}/main/likelihoodsession.h" || ! -f "${nwkit_source}/nwkit/radte_iqtree_session.py" ]]; then
  echo "Both source directories must include the persistent likelihood protocol." >&2
  exit 1
fi
revision="$(git -C "${iqtree_source}" rev-parse HEAD)"
docker build \
  --build-context "iqtree_source=${iqtree_source}" \
  --build-context "nwkit_source=${nwkit_source}" \
  --build-arg "BASE_IMAGE=${BASE_IMAGE:-local/genegalleon:dev}" \
  --build-arg "GG_BUILD_JOBS=${GG_BUILD_JOBS:-2}" \
  --build-arg "IQTREE_SOURCE_REVISION=${revision}" \
  -f "${script_dir}/Dockerfile.iqtree-session" \
  -t "${IMAGE:-local/genegalleon:iqtree-session-dev}" \
  "${script_dir}"

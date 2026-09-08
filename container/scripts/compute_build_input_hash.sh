#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
# shellcheck source=../source_branches.env
source "${repo_root}/container/source_branches.env"

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

usage() {
  cat <<'EOF'
Usage:
  bash container/scripts/compute_build_input_hash.sh [--runtime] SECURITY_REFRESH_EPOCH

The default hash includes release metadata (GeneGalleon revision and version).
--runtime omits that metadata so an identical runtime can be reused while the
mounted GeneGalleon workflow changes.
EOF
}

hash_mode="full"
if [[ "${1:-}" == "--runtime" ]]; then
  hash_mode="runtime"
  shift
fi
if [[ $# -ne 1 || ! "$1" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]]; then
  usage >&2
  exit 2
fi
security_refresh_epoch="$1"
build_target="${GG_BUILD_TARGET:-${BUILD_TARGET:-runtime}}"
if [[ "${build_target}" != "runtime" && "${build_target}" != "development" ]]; then
  echo "BUILD_TARGET must be runtime or development." >&2
  exit 2
fi

# The image architecture is validated independently from Docker/SIF metadata.
# Excluding the requested platform list lets each image from one multi-platform
# Buildx invocation carry the same content-input identity.
vcs_revision="${GG_BUILD_VCS_REF:-${vcs_ref:-unknown}}"
version="${GG_BUILD_VERSION:-${gg_version:-unknown}}"
notung_download_page="${NOTUNG_DOWNLOAD_PAGE:-https://amberjack.compbio.cs.cmu.edu/Notung/Notung-2.9.1.5.zip}"
notung_download_host_ip="${NOTUNG_DOWNLOAD_HOST_IP:-128.2.205.60}"
notung_zip_sha256="${NOTUNG_ZIP_SHA256:-81cbff670ab4d2416c01eba503f81c454aa5a724b0982373dd17510113882ae6}"
kfu52_repo_ref="${KFU52_REPO_REF:-${GG_SOURCE_NWKIT_REPO_REF}}"
amalgkit_auto="${KFU52_AMALGKIT_AUTO_SELECT_REF:-0}"
amalgkit_candidates="${KFU52_AMALGKIT_BRANCH_CANDIDATES:-master,kfdevel,devel}"
amalgkit_ref="${KFU52_AMALGKIT_REPO_REF:-${GG_SOURCE_AMALGKIT_REPO_REF}}"
amalgkit_sha="${KFU52_AMALGKIT_REPO_SHA:-}"
cdskit_ref="${KFU52_CDSKIT_REPO_REF:-${GG_SOURCE_CDSKIT_REPO_REF}}"
cdskit_sha="${KFU52_CDSKIT_REPO_SHA:-}"
csubst_ref="${KFU52_CSUBST_REPO_REF:-${GG_SOURCE_CSUBST_REPO_REF}}"
csubst_sha="${KFU52_CSUBST_REPO_SHA:-}"
nwkit_ref="${KFU52_NWKIT_REPO_REF:-${GG_SOURCE_NWKIT_REPO_REF}}"
nwkit_sha="${KFU52_NWKIT_REPO_SHA:-}"
busco_url="${BUSCO_REPO_URL:-https://gitlab.com/ezlab/busco.git}"
busco_mirror="${BUSCO_MIRROR_REPO_URL:-}"
busco_ref="${BUSCO_REPO_REF:-${GG_SOURCE_BUSCO_REPO_REF}}"
busco_sha="${BUSCO_REPO_SHA:-}"
paml_url="${PAML_REPO_URL:-https://github.com/iqtree/paml.git}"
paml_ref="${PAML_REPO_REF:-${GG_SOURCE_PAML_REPO_REF}}"
paml_sha="${PAML_REPO_SHA:-}"
kfl1ou_url="${KFL1OU_REPO_URL:-https://github.com/kfuku52/kfl1ou.git}"
kfl1ou_ref="${KFL1OU_REPO_REF:-${GG_SOURCE_KFL1OU_REPO_REF}}"
kfl1ou_sha="${KFL1OU_REPO_SHA:-}"
kffractbias_url="${KFFRACTBIAS_REPO_URL:-https://github.com/kfuku52/kfFractBias.git}"
kffractbias_ref="${KFFRACTBIAS_REPO_REF:-${GG_SOURCE_KFFRACTBIAS_REPO_REF}}"
kffractbias_sha="${KFFRACTBIAS_REPO_SHA:-}"
kftools_url="${KFTOOLS_REPO_URL:-https://github.com/kfuku52/kftools.git}"
rkftools_url="${RKFTOOLS_REPO_URL:-https://github.com/kfuku52/rkftools.git}"
radte_url="${RADTE_REPO_URL:-https://github.com/kfuku52/RADTE.git}"
kftools_ref="${KFTOOLS_REPO_REF:-${GG_SOURCE_KFTOOLS_REPO_REF}}"
rkftools_ref="${RKFTOOLS_REPO_REF:-${GG_SOURCE_RKFTOOLS_REPO_REF}}"
radte_ref="${RADTE_REPO_REF:-${GG_SOURCE_RADTE_REPO_REF}}"
kftools_sha="${KFTOOLS_REPO_SHA:-}"
rkftools_sha="${RKFTOOLS_REPO_SHA:-}"
radte_sha="${RADTE_REPO_SHA:-}"
testnh_sha="${TESTNH_TARBALL_SHA256:-598337183d2cec9c61cd364fab255a270062844b0ba5172913f7cf97512c43e2}"
cafe5_sha="${CAFE5_TARBALL_SHA256:-71871bdc74c2ffc7c1c0f4500f4742f2ff46a15cfaba78dc179d21bb1ba67ba8}"

context_digest="$(
  cd "${repo_root}"
  python3 "${script_dir}/list_build_inputs.py" | while IFS= read -r path; do
    printf '%s\t%s\n' "${path}" "$(sha256_file "${path}")"
  done | sha256_stream
)"

{
  printf '%s\n' \
    "hash_mode=${hash_mode}" \
    "build_target=${build_target}" \
    "context=${context_digest}"
  if [[ "${hash_mode}" == "full" ]]; then
    printf '%s\n' \
      "vcs_ref=${vcs_revision}" \
      "gg_version=${version}"
  fi
  printf '%s\n' \
    "security_refresh_epoch=${security_refresh_epoch}" \
    "notung_page=${notung_download_page}" \
    "notung_host=${notung_download_host_ip}" \
    "notung_sha256=${notung_zip_sha256}" \
    "kfu52_ref=${kfu52_repo_ref}" \
    "amalgkit_auto=${amalgkit_auto}" \
    "amalgkit_candidates=${amalgkit_candidates}" \
    "amalgkit_ref=${amalgkit_ref}" \
    "amalgkit_sha=${amalgkit_sha}" \
    "cdskit_ref=${cdskit_ref}" \
    "cdskit_sha=${cdskit_sha}" \
    "csubst_ref=${csubst_ref}" \
    "csubst_sha=${csubst_sha}" \
    "nwkit_ref=${nwkit_ref}" \
    "nwkit_sha=${nwkit_sha}" \
    "busco_url=${busco_url}" \
    "busco_mirror=${busco_mirror}" \
    "busco_ref=${busco_ref}" \
    "busco_sha=${busco_sha}" \
    "paml_url=${paml_url}" \
    "paml_ref=${paml_ref}" \
    "paml_sha=${paml_sha}" \
    "kfl1ou_url=${kfl1ou_url}" \
    "kfl1ou_ref=${kfl1ou_ref}" \
    "kfl1ou_sha=${kfl1ou_sha}" \
    "kffractbias_url=${kffractbias_url}" \
    "kffractbias_ref=${kffractbias_ref}" \
    "kffractbias_sha=${kffractbias_sha}" \
    "kftools_url=${kftools_url}" \
    "kftools_ref=${kftools_ref}" \
    "kftools_sha=${kftools_sha}" \
    "rkftools_url=${rkftools_url}" \
    "rkftools_ref=${rkftools_ref}" \
    "rkftools_sha=${rkftools_sha}" \
    "radte_url=${radte_url}" \
    "radte_ref=${radte_ref}" \
    "radte_sha=${radte_sha}" \
    "testnh_sha=${testnh_sha}" \
    "cafe5_sha=${cafe5_sha}"
} | sha256_stream

#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
template_path="${script_dir}/apptainer_local_build.def.template"
coverage_script="${script_dir}/scripts/check_env_coverage.sh"

IMAGE=${IMAGE:-local/genegalleon}
TAG=${TAG:-dev}
PLATFORMS=${PLATFORMS:-}
ENGINE=${ENGINE:-auto} # auto | apptainer | singularity
OUT=${OUT:-${repo_root}/genegalleon.sif}
NATIVE_BUILD_FAKEROOT=${NATIVE_BUILD_FAKEROOT:-auto} # auto | always | never
NATIVE_BUILD_KEEP_WORKDIR=${NATIVE_BUILD_KEEP_WORKDIR:-0} # 0 | 1

NOTUNG_DOWNLOAD_PAGE=${NOTUNG_DOWNLOAD_PAGE:-https://amberjack.compbio.cs.cmu.edu/Notung/Notung-2.9.1.5.zip}
NOTUNG_DOWNLOAD_HOST_IP=${NOTUNG_DOWNLOAD_HOST_IP:-128.2.205.60}
KFU52_REPO_REF=${KFU52_REPO_REF:-master}
KFU52_AMALGKIT_AUTO_SELECT_REF=${KFU52_AMALGKIT_AUTO_SELECT_REF:-1}
KFU52_AMALGKIT_BRANCH_CANDIDATES=${KFU52_AMALGKIT_BRANCH_CANDIDATES:-master,kfdevel,devel}
KFU52_AMALGKIT_REPO_REF=${KFU52_AMALGKIT_REPO_REF:-}
KFU52_AMALGKIT_REPO_SHA=${KFU52_AMALGKIT_REPO_SHA:-}
KFU52_CDSKIT_REPO_SHA=${KFU52_CDSKIT_REPO_SHA:-}
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

if [[ ! -f "${template_path}" ]]; then
  echo "Definition template not found: ${template_path}"
  exit 1
fi
if [[ ! -f "${coverage_script}" ]]; then
  echo "Coverage preflight script not found: ${coverage_script}"
  exit 1
fi

resolve_container_engine() {
  if [[ "${ENGINE}" == "apptainer" || "${ENGINE}" == "singularity" ]]; then
    if command -v "${ENGINE}" >/dev/null 2>&1; then
      echo "${ENGINE}"
      return 0
    fi
    echo "Requested container engine not found: ${ENGINE}" >&2
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

escape_sed_replacement() {
  printf '%s' "$1" | sed -e 's/[&|]/\\&/g'
}

render_definition() {
  sed \
    -e "s|@@STAGING_ROOT@@|$(escape_sed_replacement "${staging_root}")|g" \
    -e "s|@@TARGET_ARCH@@|$(escape_sed_replacement "${target_arch}")|g" \
    -e "s|@@BUILD_DATE@@|$(escape_sed_replacement "${build_date}")|g" \
    -e "s|@@VCS_REF@@|$(escape_sed_replacement "${vcs_ref}")|g" \
    -e "s|@@GG_VERSION@@|$(escape_sed_replacement "${gg_version}")|g" \
    -e "s|@@LOCAL_IMAGE_REF@@|$(escape_sed_replacement "${IMAGE}")|g" \
    -e "s|@@LOCAL_IMAGE_TAG@@|$(escape_sed_replacement "${TAG}")|g" \
    -e "s|@@NOTUNG_DOWNLOAD_PAGE@@|$(escape_sed_replacement "${NOTUNG_DOWNLOAD_PAGE}")|g" \
    -e "s|@@NOTUNG_DOWNLOAD_HOST_IP@@|$(escape_sed_replacement "${NOTUNG_DOWNLOAD_HOST_IP}")|g" \
    -e "s|@@KFU52_REPO_REF@@|$(escape_sed_replacement "${KFU52_REPO_REF}")|g" \
    -e "s|@@KFU52_AMALGKIT_AUTO_SELECT_REF@@|$(escape_sed_replacement "${KFU52_AMALGKIT_AUTO_SELECT_REF}")|g" \
    -e "s|@@KFU52_AMALGKIT_BRANCH_CANDIDATES@@|$(escape_sed_replacement "${KFU52_AMALGKIT_BRANCH_CANDIDATES}")|g" \
    -e "s|@@KFU52_AMALGKIT_REPO_REF@@|$(escape_sed_replacement "${KFU52_AMALGKIT_REPO_REF}")|g" \
    -e "s|@@KFU52_AMALGKIT_REPO_SHA@@|$(escape_sed_replacement "${KFU52_AMALGKIT_REPO_SHA}")|g" \
    -e "s|@@KFU52_CDSKIT_REPO_SHA@@|$(escape_sed_replacement "${KFU52_CDSKIT_REPO_SHA}")|g" \
    -e "s|@@KFU52_CSUBST_REPO_SHA@@|$(escape_sed_replacement "${KFU52_CSUBST_REPO_SHA}")|g" \
    -e "s|@@KFU52_NWKIT_REPO_SHA@@|$(escape_sed_replacement "${KFU52_NWKIT_REPO_SHA}")|g" \
    -e "s|@@BUSCO_REPO_URL@@|$(escape_sed_replacement "${BUSCO_REPO_URL}")|g" \
    -e "s|@@BUSCO_MIRROR_REPO_URL@@|$(escape_sed_replacement "${BUSCO_MIRROR_REPO_URL}")|g" \
    -e "s|@@BUSCO_REPO_REF@@|$(escape_sed_replacement "${BUSCO_REPO_REF}")|g" \
    -e "s|@@BUSCO_REPO_SHA@@|$(escape_sed_replacement "${BUSCO_REPO_SHA}")|g" \
    -e "s|@@PAML_REPO_URL@@|$(escape_sed_replacement "${PAML_REPO_URL}")|g" \
    -e "s|@@PAML_REPO_REF@@|$(escape_sed_replacement "${PAML_REPO_REF}")|g" \
    -e "s|@@PAML_REPO_SHA@@|$(escape_sed_replacement "${PAML_REPO_SHA}")|g" \
    -e "s|@@KFL1OU_REPO_URL@@|$(escape_sed_replacement "${KFL1OU_REPO_URL}")|g" \
    -e "s|@@KFL1OU_REPO_REF@@|$(escape_sed_replacement "${KFL1OU_REPO_REF}")|g" \
    -e "s|@@KFL1OU_REPO_SHA@@|$(escape_sed_replacement "${KFL1OU_REPO_SHA}")|g" \
    -e "s|@@KFTOOLS_REPO_URL@@|$(escape_sed_replacement "${KFTOOLS_REPO_URL}")|g" \
    -e "s|@@RKFTOOLS_REPO_URL@@|$(escape_sed_replacement "${RKFTOOLS_REPO_URL}")|g" \
    -e "s|@@RADTE_REPO_URL@@|$(escape_sed_replacement "${RADTE_REPO_URL}")|g" \
    -e "s|@@KFTOOLS_REPO_REF@@|$(escape_sed_replacement "${KFTOOLS_REPO_REF}")|g" \
    -e "s|@@RKFTOOLS_REPO_REF@@|$(escape_sed_replacement "${RKFTOOLS_REPO_REF}")|g" \
    -e "s|@@RADTE_REPO_REF@@|$(escape_sed_replacement "${RADTE_REPO_REF}")|g" \
    -e "s|@@KFTOOLS_REPO_SHA@@|$(escape_sed_replacement "${KFTOOLS_REPO_SHA}")|g" \
    -e "s|@@RKFTOOLS_REPO_SHA@@|$(escape_sed_replacement "${RKFTOOLS_REPO_SHA}")|g" \
    -e "s|@@RADTE_REPO_SHA@@|$(escape_sed_replacement "${RADTE_REPO_SHA}")|g" \
    -e "s|@@TESTNH_TARBALL_SHA256@@|$(escape_sed_replacement "${TESTNH_TARBALL_SHA256}")|g" \
    -e "s|@@CAFE5_TARBALL_SHA256@@|$(escape_sed_replacement "${CAFE5_TARBALL_SHA256}")|g" \
    "${template_path}" > "${definition_path}"
}

case "${NATIVE_BUILD_FAKEROOT}" in
  auto|always|never)
    ;;
  *)
    echo "NATIVE_BUILD_FAKEROOT must be one of: auto, always, never"
    exit 1
    ;;
esac
if [[ "${NATIVE_BUILD_KEEP_WORKDIR}" != "0" && "${NATIVE_BUILD_KEEP_WORKDIR}" != "1" ]]; then
  echo "NATIVE_BUILD_KEEP_WORKDIR must be 0 or 1."
  exit 1
fi

host_arch="$(uname -m)"
case "${host_arch}" in
  aarch64|arm64)
    target_arch="arm64"
    platform="linux/arm64"
    ;;
  x86_64|amd64)
    target_arch="amd64"
    platform="linux/amd64"
    ;;
  *)
    echo "Unsupported host architecture for native local build: ${host_arch}"
    exit 1
    ;;
esac

if [[ -n "${PLATFORMS}" ]]; then
  if [[ "${PLATFORMS}" == *","* ]]; then
    echo "Native local build supports a single platform only. Set PLATFORMS=${platform}."
    exit 1
  fi
  if [[ "${PLATFORMS}" != "${platform}" ]]; then
    echo "Native local build only supports the host architecture."
    echo "Set PLATFORMS=${platform}."
    exit 1
  fi
fi

resolved_engine="$(resolve_container_engine)"
build_date="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
vcs_ref="$(git -C "${repo_root}" rev-parse --short HEAD 2>/dev/null || echo "unknown")"
gg_version="$(sed -n '1p' "${repo_root}/VERSION" 2>/dev/null | tr -d '\r' || true)"
if [[ -z "${gg_version}" ]]; then
  gg_version="unknown"
fi

work_dir="$(mktemp -d)"
staging_root="${work_dir}/context"
definition_path="${work_dir}/genegalleon_local.def"

cleanup() {
  if [[ "${NATIVE_BUILD_KEEP_WORKDIR}" == "1" ]]; then
    echo "[apptainer_local_build] kept workdir: ${work_dir}"
    return
  fi
  rm -rf -- "${work_dir}"
}
trap cleanup EXIT

mkdir -p "${staging_root}"
cp -R "${repo_root}/container/env" "${staging_root}/"
cp -R "${repo_root}/container/spec" "${staging_root}/"
cp -R "${repo_root}/container/testdata" "${staging_root}/"
cp -R "${repo_root}/container/scripts" "${staging_root}/"

(
  cd "${repo_root}"
  bash "${coverage_script}" .
)

render_definition

build_args=()
if [[ "${NATIVE_BUILD_FAKEROOT}" == "always" ]]; then
  build_args+=(--fakeroot)
fi

mkdir -p "$(dirname "${OUT}")"
echo "[apptainer_local_build] engine=${resolved_engine}"
echo "[apptainer_local_build] platform=${platform}"
echo "[apptainer_local_build] definition=${definition_path}"

if [[ ${#build_args[@]} -gt 0 ]]; then
  if ! "${resolved_engine}" build "${build_args[@]}" "${OUT}" "${definition_path}"; then
    echo "[apptainer_local_build] Native local build failed."
    echo "[apptainer_local_build] Retry with NATIVE_BUILD_FAKEROOT=always if your site supports fakeroot,"
    echo "[apptainer_local_build] or use IMAGE_SOURCE=public to pull the published image."
    exit 1
  fi
else
  if ! "${resolved_engine}" build "${OUT}" "${definition_path}"; then
    echo "[apptainer_local_build] Native local build failed."
    echo "[apptainer_local_build] Retry with NATIVE_BUILD_FAKEROOT=always if your site supports fakeroot,"
    echo "[apptainer_local_build] or use IMAGE_SOURCE=public to pull the published image."
    exit 1
  fi
fi

echo "Generated ${OUT}"

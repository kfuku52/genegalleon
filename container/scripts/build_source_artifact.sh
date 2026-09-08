#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ $# -lt 5 || $# -gt 6 ]]; then
  echo "Usage: $0 python|r|paml|radte SOURCE REPO_URL REVISION OUTPUT_DIR [MIRROR_URL]" >&2
  exit 2
fi
kind="$1"
source_name="$2"
repo_url="$3"
revision="$4"
output_dir="$5"
mirror_url="${6:-}"
jobs="${GG_BUILD_JOBS:-2}"
if [[ ! "${jobs}" =~ ^[1-9][0-9]*$ || ! "${source_name}" =~ ^[A-Za-z][A-Za-z0-9]*$ ]]; then
  echo "Invalid build job count or source name." >&2
  exit 2
fi
case "${kind}" in
  python|r|paml|radte) ;;
  *) echo "Unknown source artifact kind: ${kind}" >&2; exit 2 ;;
esac
if [[ ! "${revision}" =~ ^[0-9a-f]{40}$ ]]; then
  # Standard wrappers already resolve all branches once for every platform.
  # A direct stage build can still follow its configured moving branch.
  revision="$(bash "${script_dir}/resolve_git_branch_sha.sh" "${repo_url}" "${revision}")"
fi
work_dir="$(mktemp -d)"
trap 'rm -rf -- "${work_dir}"' EXIT
bash "${script_dir}/fetch_git_repo.sh" "${repo_url}" "${revision}" "${work_dir}/source" "${mirror_url}"
mkdir -p "${output_dir}"
echo "[source artifact] Building ${source_name} at ${revision} (${kind})"
case "${kind}" in
  python)
    mkdir -p "${output_dir}/wheels"
    micromamba run -n base python -m pip wheel \
      --no-deps --no-build-isolation --wheel-dir "${output_dir}/wheels" "${work_dir}/source"
    shopt -s nullglob
    wheels=("${output_dir}/wheels/"*.whl)
    if [[ ${#wheels[@]} -ne 1 ]]; then
      echo "Expected one wheel for ${source_name}, found ${#wheels[@]}." >&2
      exit 1
    fi
    ;;
  r)
    library="${output_dir}/rootfs/opt/conda/lib/R/library"
    mkdir -p "${library}"
    MAKEFLAGS="-j${jobs}" CMAKE_BUILD_PARALLEL_LEVEL="${jobs}" \
      micromamba run -n base R CMD INSTALL --library="${library}" "${work_dir}/source"
    ;;
  paml)
    make -C "${work_dir}/source/src" -j"${jobs}" mcmctree \
      || make -C "${work_dir}/source/src" -j"${jobs}"
    binary="${work_dir}/source/bin/mcmctree"
    if [[ ! -x "${binary}" ]]; then binary="${work_dir}/source/src/mcmctree"; fi
    install -D -m 0755 "${binary}" "${output_dir}/rootfs/usr/local/bin/mcmctree"
    mkdir -p "${output_dir}/rootfs/opt/conda/bin"
    ln -s /usr/local/bin/mcmctree "${output_dir}/rootfs/opt/conda/bin/mcmctree"
    ;;
  radte)
    if [[ ! -s "${work_dir}/source/radte.r" ]]; then
      echo "RADTE script was not found in ${repo_url}" >&2
      exit 1
    fi
    install -D -m 0755 "${work_dir}/source/radte.r" "${output_dir}/rootfs/usr/local/bin/radte.r"
    ;;
esac
printf '%s\t%s\n' "${source_name}" "${revision}" > "${output_dir}/source.tsv"

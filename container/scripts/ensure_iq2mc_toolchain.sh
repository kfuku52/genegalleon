#!/usr/bin/env bash
set -euo pipefail

env_name=${1:-base}
if [[ "${env_name}" == "base" ]]; then
  env_bin="/opt/conda/bin"
else
  env_bin="/opt/conda/envs/${env_name}/bin"
fi
jobs=${GG_BUILD_JOBS:-$(nproc)}
if [[ "${jobs}" -lt 1 ]]; then
  jobs=1
fi

log() {
  echo "[ensure_iq2mc_toolchain] $*"
}

check_iq2mc_support() {
  local bin=$1
  local resolved
  resolved=$(command -v "${bin}" 2>/dev/null || true)
  if [[ -z "${resolved}" && -x "${bin}" ]]; then
    resolved="${bin}"
  fi

  if "${bin}" -h 2>&1 | grep -q -- "--mcmc-bds"; then
    return 0
  fi
  if command -v strings >/dev/null 2>&1 && [[ -n "${resolved}" ]]; then
    if strings "${resolved}" 2>/dev/null | grep -q -- "--mcmc-bds"; then
      return 0
    fi
  fi
  return 1
}

install_into_path() {
  local src_bin=$1
  local dest_name=$2
  install -m 0755 "${src_bin}" "/usr/local/bin/${dest_name}"
  if [[ -d "${env_bin}" ]]; then
    install -m 0755 "${src_bin}" "${env_bin}/${dest_name}"
  fi
}

ensure_iqtree_commands() {
  # Keep iqtree3 as canonical command while exposing iqtree alias.
  local primary=""
  if [[ -x "${env_bin}/iqtree3" ]]; then
    primary="${env_bin}/iqtree3"
  elif [[ -x "${env_bin}/iqtree" ]]; then
    primary="${env_bin}/iqtree"
  elif command -v iqtree3 >/dev/null 2>&1; then
    primary=$(command -v iqtree3)
  elif command -v iqtree >/dev/null 2>&1; then
    primary=$(command -v iqtree)
  fi

  if [[ -z "${primary}" ]]; then
    return 1
  fi

  if [[ "${primary}" != "/usr/local/bin/iqtree3" ]]; then
    ln -sf "${primary}" /usr/local/bin/iqtree3
  fi
  ln -sf /usr/local/bin/iqtree3 /usr/local/bin/iqtree
  rm -f /usr/local/bin/iqtree2

  if [[ -d "${env_bin}" ]]; then
    if [[ "${primary}" != "${env_bin}/iqtree3" ]]; then
      ln -sf "${primary}" "${env_bin}/iqtree3"
    fi
    ln -sf "${env_bin}/iqtree3" "${env_bin}/iqtree"
    rm -f "${env_bin}/iqtree2"
  fi

  return 0
}

clone_branch() {
  local repo_url=$1
  local dest_dir=$2
  shift 2
  local branch
  for branch in "$@"; do
    if git clone --depth 1 --recursive --branch "${branch}" "${repo_url}" "${dest_dir}"; then
      log "Cloned ${repo_url} branch '${branch}'"
      return 0
    fi
    rm -rf "${dest_dir}"
  done
  return 1
}

build_mcmctree_iq2mc() {
  local tmpdir paml_src mcmctree_bin
  tmpdir=$(mktemp -d)
  paml_src="${tmpdir}/paml"

  if ! clone_branch \
    "https://github.com/iqtree/paml.git" \
    "${paml_src}" \
    "master"; then
    log "ERROR: Failed to clone IQ2MC-compatible paml source branch."
    rm -rf "${tmpdir}"
    exit 1
  fi

  make -C "${paml_src}/src" -j"${jobs}" mcmctree >/dev/null 2>&1 || make -C "${paml_src}/src" -j"${jobs}" >/dev/null

  mcmctree_bin=""
  for candidate in \
    "${paml_src}/bin/mcmctree" \
    "${paml_src}/src/mcmctree"; do
    if [[ -x "${candidate}" ]]; then
      mcmctree_bin="${candidate}"
      break
    fi
  done
  if [[ -z "${mcmctree_bin}" ]]; then
    log "ERROR: mcmctree build finished but binary was not found."
    rm -rf "${tmpdir}"
    exit 1
  fi

  install_into_path "${mcmctree_bin}" "mcmctree"
  rm -rf "${tmpdir}"
}

main() {
  if ! command -v git >/dev/null 2>&1; then
    log "ERROR: git is required to build IQ2MC toolchain."
    exit 1
  fi
  if ! command -v make >/dev/null 2>&1; then
    log "ERROR: make is required to build mcmctree."
    exit 1
  fi

  ensure_iqtree_commands || true
  if ! command -v iqtree3 >/dev/null 2>&1; then
    log "ERROR: iqtree3 command not found in the conda environment."
    exit 1
  fi
  if check_iq2mc_support "$(command -v iqtree3)"; then
    log "Conda IQ-TREE exposes IQ2MC flag (--mcmc-bds): $(command -v iqtree3)"
  else
    log "Conda IQ-TREE does not expose --mcmc-bds; keeping conda-provided binary without source reinstall."
  fi

  log "Installing IQ2MC-compatible mcmctree binary"
  build_mcmctree_iq2mc

  ensure_iqtree_commands || true
  if ! command -v iqtree3 >/dev/null 2>&1; then
    log "ERROR: iqtree3 command not found after command alias setup."
    exit 1
  fi
  if ! command -v mcmctree >/dev/null 2>&1; then
    log "ERROR: mcmctree command not found after IQ2MC installation."
    exit 1
  fi

  log "IQ2MC toolchain is ready: iqtree3=$(command -v iqtree3), iqtree=$(command -v iqtree), mcmctree=$(command -v mcmctree)"
}

main "$@"

#!/usr/bin/env bash
set -euo pipefail

tool_env_name=${GG_TOOL_ENV_NAME:-base}
r_env_name=${GG_R_ENV_NAME:-base}

env_prefix() {
  local env_name=$1
  if [[ "${env_name}" == "base" ]]; then
    echo "/opt/conda"
  else
    echo "/opt/conda/envs/${env_name}"
  fi
}

tool_env_bin="$(env_prefix "${tool_env_name}")/bin"

log() {
  echo "[install_nonconda_fallbacks] $*"
}

download_github_tag_tarball() {
  local owner=$1
  local repo=$2
  local tag=$3
  local dest=$4
  local url="https://codeload.github.com/${owner}/${repo}/tar.gz/refs/tags/${tag}"

  mkdir -p "${dest}"
  curl -fL --retry 5 --retry-all-errors --retry-delay 2 "${url}" \
    | tar -xzf - -C "${dest}" --strip-components=1
}

build_jobs() {
  local jobs=${GG_BUILD_JOBS:-$(nproc)}
  if [[ "${jobs}" -lt 1 ]]; then
    jobs=1
  fi
  echo "${jobs}"
}

resolve_binary() {
  local candidate
  for candidate in "$@"; do
    if [[ "${candidate}" == /* ]]; then
      if [[ -x "${candidate}" ]]; then
        echo "${candidate}"
        return 0
      fi
    elif command -v "${candidate}" >/dev/null 2>&1; then
      command -v "${candidate}"
      return 0
    elif [[ -x "${tool_env_bin}/${candidate}" ]]; then
      echo "${tool_env_bin}/${candidate}"
      return 0
    fi
  done
  return 1
}

install_mapnh() {
  if command -v mapnh >/dev/null 2>&1; then
    log "mapnh already available: $(command -v mapnh)"
    return
  fi

  local workdir
  workdir=$(mktemp -d)
  log "Building mapnh from BioPP/testnh v2.3.2 in ${workdir}"
  download_github_tag_tarball "BioPP" "testnh" "v2.3.2" "${workdir}/testnh"

  local cmake_file="${workdir}/testnh/CMakeLists.txt"
  sed -i -E 's/find_package[[:space:]]*\([[:space:]]*bpp-phyl[[:space:]]+[0-9.]+[[:space:]]+REQUIRED[[:space:]]*\)/find_package(bpp-phyl REQUIRED)/' "${cmake_file}"

  local jobs
  jobs=$(build_jobs)
  cmake -S "${workdir}/testnh" -B "${workdir}/testnh/build" -DCMAKE_BUILD_TYPE=Release
  cmake --build "${workdir}/testnh/build" --target mapnh -j"${jobs}"

  install -m 0755 "${workdir}/testnh/build/TestNH/mapnh" /usr/local/bin/mapnh
  rm -rf "${workdir}"

  if ! command -v mapnh >/dev/null 2>&1; then
    log "ERROR: mapnh was not installed correctly."
    exit 1
  fi
  log "Installed mapnh to /usr/local/bin/mapnh"
}

install_maxalign_wrapper() {
  local wrapper_tmp
  wrapper_tmp=$(mktemp)

  cat > "${wrapper_tmp}" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: maxalign.pl -f=<output_prefix> [options] <input_fasta>
  -f, --prefix             Output prefix (writes <prefix>heuristic.fsa)
  -i, --max-iterations     MaxAlign iterations (default: 999)
  -k, --keep-sequence      Keep sequence id (can be used multiple times)
  -v                       Accepted for compatibility (ignored)
  -V, --version            Show backend version
USAGE
}

find_backend() {
  local candidate
  for candidate in \
    maxalign-rs \
    __TOOL_ENV_BIN__/maxalign-rs \
    maxalign \
    __TOOL_ENV_BIN__/maxalign
  do
    if [[ "${candidate}" == /* ]]; then
      if [[ -x "${candidate}" ]]; then
        echo "${candidate}"
        return 0
      fi
    elif command -v "${candidate}" >/dev/null 2>&1; then
      command -v "${candidate}"
      return 0
    fi
  done
  return 1
}

iterations=999
prefix=""
input_file=""
declare -a keep_ids=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    -V|--version)
      backend=$(find_backend) || { echo "No maxalign backend found (maxalign-rs/maxalign)." >&2; exit 1; }
      exec "${backend}" --version
      ;;
    -i=*|--max-iterations=*)
      iterations="${1#*=}"
      shift
      ;;
    -i|--max-iterations)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
      iterations="$2"
      shift 2
      ;;
    -f=*|--prefix=*)
      prefix="${1#*=}"
      shift
      ;;
    -f|--prefix)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
      prefix="$2"
      shift 2
      ;;
    -k=*|--keep-sequence=*)
      keep_ids+=("${1#*=}")
      shift
      ;;
    -k|--keep-sequence)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
      keep_ids+=("$2")
      shift 2
      ;;
    -v|-v=*|--verbose|--verbose=*|--verbosity|--verbosity=*)
      if [[ "$1" == "-v" || "$1" == "--verbose" || "$1" == "--verbosity" ]]; then
        if [[ $# -ge 2 && "$2" != -* ]]; then
          shift 2
        else
          shift
        fi
      else
        shift
      fi
      ;;
    --)
      shift
      break
      ;;
    -*)
      echo "Unsupported option: $1" >&2
      usage >&2
      exit 2
      ;;
    *)
      if [[ -z "${input_file}" ]]; then
        input_file="$1"
      else
        echo "Unexpected argument: $1" >&2
        exit 2
      fi
      shift
      ;;
  esac
done

if [[ -z "${input_file}" && $# -gt 0 ]]; then
  input_file="$1"
  shift
fi
if [[ -n "${input_file}" && $# -gt 0 ]]; then
  echo "Unexpected extra arguments: $*" >&2
  exit 2
fi
if [[ -z "${input_file}" ]]; then
  usage >&2
  exit 2
fi
if [[ -z "${prefix}" ]]; then
  prefix="maxalign."
fi

if [[ -f "${input_file}" ]]; then
  while IFS= read -r header; do
    header="${header#>}"
    [[ -n "${header}" ]] || continue
    found=0
    for existing in "${keep_ids[@]}"; do
      if [[ "${existing}" == "${header}" ]]; then
        found=1
        break
      fi
    done
    if [[ ${found} -eq 0 ]]; then
      keep_ids+=("${header}")
    fi
  done < <(grep '^>+' "${input_file}" || true)
fi

backend=$(find_backend) || { echo "No maxalign backend found (maxalign-rs/maxalign)." >&2; exit 1; }
output_file="${prefix}heuristic.fsa"
declare -a cmd_args
cmd_args=(--max-iterations "${iterations}")
for keep_id in "${keep_ids[@]}"; do
  cmd_args+=(--keep-sequence "${keep_id}")
done

run_backend() {
  "${backend}" "${cmd_args[@]}" "$@"
}

if ! run_backend "${input_file}" "${output_file}"; then
  if ! run_backend --input "${input_file}" --output "${output_file}"; then
    if ! run_backend --infile "${input_file}" --outfile "${output_file}"; then
      echo "Failed to run maxalign backend with supported argument styles." >&2
      exit 1
    fi
  fi
fi

if [[ ! -s "${output_file}" ]]; then
  echo "Expected output was not created: ${output_file}" >&2
  exit 1
fi
EOF
  sed -i "s|__TOOL_ENV_BIN__|${tool_env_bin}|g" "${wrapper_tmp}"

  install -m 0755 "${wrapper_tmp}" /usr/local/bin/maxalign.pl
  if [[ -d "${tool_env_bin}" ]]; then
    install -m 0755 "${wrapper_tmp}" "${tool_env_bin}/maxalign.pl"
  fi
  rm -f "${wrapper_tmp}"
  log "Installed maxalign.pl compatibility wrapper"
}

install_astral_hybrid_wrapper() {
  local existing_wrapper
  if existing_wrapper=$(resolve_binary astral-hybrid); then
    log "astral-hybrid already available: ${existing_wrapper}"
    return
  fi

  local backend
  if ! backend=$(resolve_binary astral); then
    log "ERROR: astral backend was not found."
    exit 1
  fi

  local wrapper_tmp
  wrapper_tmp=$(mktemp)
  cat > "${wrapper_tmp}" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

input=""
output=""
branch_annotate="2"
declare -a pass_args=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input=*)
      input="${1#*=}"
      shift
      ;;
    --input)
      [[ $# -ge 2 ]] || { echo "Missing value for --input" >&2; exit 2; }
      input="$2"
      shift 2
      ;;
    --output=*)
      output="${1#*=}"
      shift
      ;;
    --output)
      [[ $# -ge 2 ]] || { echo "Missing value for --output" >&2; exit 2; }
      output="$2"
      shift 2
      ;;
    --support=*)
      branch_annotate="${1#*=}"
      shift
      ;;
    --support)
      [[ $# -ge 2 ]] || { echo "Missing value for --support" >&2; exit 2; }
      branch_annotate="$2"
      shift 2
      ;;
    --mode|--thread)
      if [[ $# -ge 2 ]]; then
        shift 2
      else
        shift
      fi
      ;;
    --mode=*|--thread=*)
      shift
      ;;
    --help|-h)
      cat <<'USAGE'
Usage: astral-hybrid --input FILE --output FILE [--support N] [--mode N] [--thread N]
Compatibility wrapper around ASTRAL v5 (`astral`).
USAGE
      exit 0
      ;;
    *)
      pass_args+=("$1")
      shift
      ;;
  esac
done

if [[ -z "${input}" || -z "${output}" ]]; then
  echo "astral-hybrid wrapper requires --input and --output." >&2
  exit 2
fi

exec __ASTRAL_BACKEND__ --input "${input}" --output "${output}" --branch-annotate "${branch_annotate}" "${pass_args[@]}"
EOF
  sed -i "s|__ASTRAL_BACKEND__|${backend}|g" "${wrapper_tmp}"

  install -m 0755 "${wrapper_tmp}" /usr/local/bin/astral-hybrid
  if [[ -d "${tool_env_bin}" ]]; then
    install -m 0755 "${wrapper_tmp}" "${tool_env_bin}/astral-hybrid"
  fi
  rm -f "${wrapper_tmp}"
  log "Installed astral-hybrid compatibility wrapper"
}

install_cafe5() {
  if command -v cafe5 >/dev/null 2>&1; then
    log "cafe5 already available: $(command -v cafe5)"
    return
  fi

  local workdir
  local jobs
  workdir=$(mktemp -d)
  jobs=$(build_jobs)
  log "Building CAFE5 v5.1 in ${workdir}"

  curl -fL --retry 5 --retry-all-errors --retry-delay 2 \
    "https://github.com/hahnlab/CAFE5/releases/download/v5.1/CAFE5-5.1.0.tar.gz" \
    | tar -xzf - -C "${workdir}"

  local src_dir="${workdir}/CAFE5"
  if [[ ! -d "${src_dir}" ]]; then
    log "ERROR: CAFE5 source directory was not found after extraction."
    exit 1
  fi

  (
    cd "${src_dir}"
    ./configure
    make -j"${jobs}"
  )

  if [[ ! -x "${src_dir}/bin/cafe5" ]]; then
    log "ERROR: CAFE5 binary was not built: ${src_dir}/bin/cafe5"
    exit 1
  fi

  install -m 0755 "${src_dir}/bin/cafe5" /usr/local/bin/cafe5
  if [[ -d "${tool_env_bin}" ]]; then
    install -m 0755 "${src_dir}/bin/cafe5" "${tool_env_bin}/cafe5"
  fi
  rm -rf "${workdir}"

  if ! command -v cafe5 >/dev/null 2>&1; then
    log "ERROR: cafe5 was not installed correctly."
    exit 1
  fi
  log "Installed cafe5 to /usr/local/bin/cafe5"
}

install_r_cran_packages() {
  local env_name=$1
  shift
  local packages=("$@")
  if [[ ${#packages[@]} -eq 0 ]]; then
    return
  fi

  local pkg_expr=""
  local pkg
  local jobs
  jobs=$(build_jobs)
  for pkg in "${packages[@]}"; do
    pkg_expr+="\"${pkg}\","
  done
  pkg_expr="${pkg_expr%,}"

  log "Ensuring CRAN packages in '${env_name}' with ${jobs} job(s): ${packages[*]}"
  MAKEFLAGS="-j${jobs}" CMAKE_BUILD_PARALLEL_LEVEL="${jobs}" \
    micromamba run -n "${env_name}" Rscript -e "options(repos=c(CRAN='https://cloud.r-project.org')); options(Ncpus=${jobs}L); dep_levels <- c('Depends','Imports','LinkingTo'); pkgs <- c(${pkg_expr}); missing <- pkgs[!vapply(pkgs, requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))]; if (length(missing) > 0) install.packages(missing, dependencies=dep_levels)"
  micromamba run -n "${env_name}" Rscript -e "pkgs <- c(${pkg_expr}); missing <- pkgs[!vapply(pkgs, requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))]; if (length(missing) > 0) stop(sprintf('Missing packages in env %s: %s', '${env_name}', paste(missing, collapse=', ')))"
}

install_r_l1ou() {
  local jobs
  jobs=$(build_jobs)
  log "Ensuring GitHub package in '${r_env_name}' with ${jobs} job(s): l1ou"
  MAKEFLAGS="-j${jobs}" CMAKE_BUILD_PARALLEL_LEVEL="${jobs}" \
    micromamba run -n "${r_env_name}" Rscript -e "options(repos=c(CRAN='https://cloud.r-project.org')); options(Ncpus=${jobs}L); dep_levels <- c('Depends','Imports','LinkingTo'); if (!requireNamespace('l1ou', quietly=TRUE)) { if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', dependencies=dep_levels); remotes::install_github('khabbazian/l1ou', dependencies=NA, upgrade='never'); }; if (!requireNamespace('l1ou', quietly=TRUE)) stop('Missing package in env ${r_env_name}: l1ou')"
}

verify_plotting_packages_in_r() {
  log "Verifying plotting packages are available in '${r_env_name}'"
  micromamba run -n "${r_env_name}" Rscript -e "pkgs <- c('Rphylopars','ape','aplot','cowplot','ggimage','ggmsa','ggplot2','ggrepel','ggtree','magick','phangorn','svglite','viridis','xml2'); missing <- pkgs[!vapply(pkgs, requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))]; if (length(missing) > 0) stop(sprintf('Missing packages in env ${r_env_name}: %s', paste(missing, collapse=', ')))"
}

main() {
  if ! command -v provean >/dev/null 2>&1; then
    log "ERROR: provean not found. Install it via apt before running this script."
    exit 1
  fi

  install_cafe5
  install_mapnh
  install_maxalign_wrapper
  install_astral_hybrid_wrapper
  install_r_cran_packages "${r_env_name}" PhylogeneticEM Rphylopars
  install_r_l1ou
  verify_plotting_packages_in_r
}

main "$@"

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

run_mapnh_smoke_test() {
  local mapnh_bin=$1
  local smoke_dir
  smoke_dir=$(mktemp -d)

  cat > "${smoke_dir}/tiny.fasta" <<'EOF'
>t1
ACGTTGCAACGTTGCAACGTTGCAACGTTGCA
>t2
ACGTTGCAACGTTGCAACGTTGCAACGTCGCA
>t3
ACGTTGCAACGTTGCAACGTTGCAACGTCGTA
>t4
ACGTTGCAACGTTGCAACGTTGCAACGCCGTA
EOF
  cat > "${smoke_dir}/tiny.nwk" <<'EOF'
((t1:0.1,t2:0.1):0.2,(t3:0.1,t4:0.1):0.2);
EOF
  cat > "${smoke_dir}/tiny.params" <<'EOF'
alphabet=DNA
input.sequence.file=$(SEQ)
input.sequence.format=Fasta
input.tree.file=$(TREE)
input.tree.format=Newick
model=JC69
rate_distribution=Constant()
map.type=Total
output.counts=PerType(prefix=$(OUT).)
output.tree_with_id.file=$(OUT).with_id.nwk
EOF

  if ! "${mapnh_bin}" \
    SEQ="${smoke_dir}/tiny.fasta" \
    TREE="${smoke_dir}/tiny.nwk" \
    OUT="${smoke_dir}/tiny_out" \
    param="${smoke_dir}/tiny.params" \
    > "${smoke_dir}/tiny.log" 2>&1; then
    log "mapnh smoke test failed for ${mapnh_bin}"
    sed -n '1,120p' "${smoke_dir}/tiny.log" >&2 || true
    rm -rf "${smoke_dir}"
    return 1
  fi

  rm -rf "${smoke_dir}"
}

install_mapnh() {
  if command -v mapnh >/dev/null 2>&1; then
    local existing_mapnh
    existing_mapnh=$(command -v mapnh)
    if run_mapnh_smoke_test "${existing_mapnh}"; then
      log "mapnh already available and passed smoke test: ${existing_mapnh}"
      return
    fi
    log "Existing mapnh is broken; rebuilding: ${existing_mapnh}"
  fi

  local workdir
  workdir=$(mktemp -d)
  log "Building mapnh from BioPP/testnh v2.3.2 in ${workdir}"
  download_github_tag_tarball "BioPP" "testnh" "v2.3.2" "${workdir}/testnh"

  local cmake_file="${workdir}/testnh/CMakeLists.txt"
  sed -i -E 's/find_package[[:space:]]*\([[:space:]]*bpp-phyl[[:space:]]+[0-9.]+[[:space:]]+REQUIRED[[:space:]]*\)/find_package(bpp-phyl REQUIRED)/' "${cmake_file}"
  # Bio++ >=12 with testnh v2.3.2 triggers a homogeneous-mode null dereference in MapNH.cpp.
  # Guard mapnh code paths so they use `model` when nonhomogeneous mode is disabled.
  local mapnh_cpp="${workdir}/testnh/TestNH/MapNH.cpp"
  sed -i 's/modelSet->getSubstitutionModel(0)/model ? model : modelSet->getSubstitutionModel(0)/g' "${mapnh_cpp}"

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
  if ! run_mapnh_smoke_test "/usr/local/bin/mapnh"; then
    log "ERROR: mapnh smoke test failed after installation."
    exit 1
  fi
  log "Installed mapnh to /usr/local/bin/mapnh"
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
  micromamba run -n "${r_env_name}" Rscript -e "pkgs <- c('Rphylopars','ape','aplot','cowplot','ggmsa','ggplot2','ggrepel','ggtree','phangorn','svglite','viridis','xml2'); missing <- pkgs[!vapply(pkgs, requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))]; if (length(missing) > 0) stop(sprintf('Missing packages in env ${r_env_name}: %s', paste(missing, collapse=', ')))"
}

main() {
  install_cafe5
  install_mapnh
  install_astral_hybrid_wrapper
  install_r_cran_packages "${r_env_name}" PhylogeneticEM Rphylopars
  install_r_l1ou
  verify_plotting_packages_in_r
}

main "$@"

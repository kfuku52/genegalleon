#!/usr/bin/env bash
set -euo pipefail

env_name=${1:-base}
fallback_spec="generax=2.1.3=*_2"
channel_args=(-c conda-forge -c bioconda -c defaults)

log() {
  echo "[ensure_generax_stable] $*"
}

print_generax_pkg() {
  micromamba list -n "${env_name}" | awk '$1=="generax"{print}'
}

run_smoke_test() {
  local workdir log_file rc
  workdir=$(mktemp -d)
  log_file="${workdir}/generax_smoke.log"

  cp /opt/pg/testdata/generax_smoke/* "${workdir}/"

  pushd "${workdir}" >/dev/null
  set +e
  micromamba run -n "${env_name}" \
    env OMPI_MCA_plm=isolated \
    mpiexec --allow-run-as-root -oversubscribe -np 1 \
    generax \
      --families families.txt \
      --species-tree species_tree.newick \
      --rec-model UndatedDL \
      --per-family-rates \
      --prefix smoke_out \
      --max-spr-radius 1 \
      > "${log_file}" 2>&1
  rc=$?
  set -e
  popd >/dev/null

  if [[ ${rc} -ne 0 ]]; then
    log "Smoke test failed with exit code ${rc}."
    tail -n 80 "${log_file}" || true
    rm -rf -- "${workdir}"
    return 1
  fi

  if [[ ! -s "${workdir}/smoke_out/results/Phy003AEDB_CUCME/geneTree.newick" ]]; then
    log "Smoke test did not produce expected output gene tree."
    tail -n 80 "${log_file}" || true
    rm -rf -- "${workdir}"
    return 1
  fi

  rm -rf -- "${workdir}"
  return 0
}

install_fallback() {
  log "Installing fallback GeneRax build: ${fallback_spec}"
  micromamba install -y -n "${env_name}" "${channel_args[@]}" "${fallback_spec}"
}

main() {
  if ! micromamba run -n "${env_name}" bash -lc "command -v generax >/dev/null 2>&1"; then
    log "ERROR: 'generax' was not found in env '${env_name}'."
    exit 1
  fi

  log "Installed GeneRax package: $(print_generax_pkg)"
  log "Running GeneRax runtime smoke test"

  if run_smoke_test; then
    log "Smoke test passed with current GeneRax build."
    return
  fi

  log "Current GeneRax build is unstable. Applying fallback."
  install_fallback
  log "GeneRax package after fallback: $(print_generax_pkg)"
  log "Re-running smoke test after fallback"

  if ! run_smoke_test; then
    log "ERROR: Smoke test still failed after fallback install."
    exit 1
  fi

  log "Fallback GeneRax build passed smoke test."
}

main "$@"

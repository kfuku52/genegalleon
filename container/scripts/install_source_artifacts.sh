#!/usr/bin/env bash
set -euo pipefail

artifact_root="${1:?Usage: install_source_artifacts.sh ARTIFACT_ROOT}"
sources=(amalgkit cdskit csubst nwkit BUSCO paml iqtree kfl1ou kfFractBias kftools rkftools)
wheels=()
shopt -s nullglob
for source_name in "${sources[@]}"; do
  artifact="${artifact_root}/${source_name}"
  IFS=$'\t' read -r recorded_source revision < "${artifact}/source.tsv"
  if [[ "${recorded_source}" != "${source_name}" || ! "${revision}" =~ ^[0-9a-f]{40}$ ]]; then
    echo "Invalid source artifact identity: ${source_name}" >&2
    exit 1
  fi
  wheels+=("${artifact}/wheels/"*.whl)
done
if [[ ${#wheels[@]} -ne 7 ]]; then
  echo "Expected seven upstream Python wheels, found ${#wheels[@]}." >&2
  exit 1
fi
# Artifacts are read-only BuildKit mounts, not COPY layers. No wheel archives or
# source checkouts remain in the runtime, and dependencies are never re-solved.
micromamba run -n base python -m pip install \
  --no-index --no-deps --force-reinstall "${wheels[@]}"
mkdir -p /opt/pg/logs
printf 'source\trevision\n' > /opt/pg/logs/source_revisions.tsv
for source_name in "${sources[@]}"; do
  artifact="${artifact_root}/${source_name}"
  if [[ -d "${artifact}/rootfs" ]]; then
    cp -a "${artifact}/rootfs/." /
  fi
  cat "${artifact}/source.tsv" >> /opt/pg/logs/source_revisions.tsv
done

# Both IQ-TREE entrypoints and the external worker come from one source build.
# The separate IQ2MC-compatible PAML artifact remains required by species dating.
micromamba run -n base python -m nwkit.iqtree_library check --interface library \
  > /opt/pg/logs/iqtree3_library_worker.json
micromamba run -n base python -m pip check
micromamba run -n base Rscript -e \
  'stopifnot(requireNamespace("kfl1ou", quietly = TRUE), requireNamespace("rkftools", quietly = TRUE))'

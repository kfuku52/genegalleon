#!/usr/bin/env bash
set -euo pipefail

artifact_root="${1:?Usage: install_source_artifacts.sh ARTIFACT_ROOT}"
sources=(amalgkit cdskit csubst nwkit BUSCO paml kfl1ou kfFractBias kftools rkftools RADTE)
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

# Preserve the existing IQ-TREE entrypoints while using the independently
# compiled IQ2MC-compatible mcmctree artifact.
if [[ -x /opt/conda/bin/iqtree3 ]]; then
  ln -sf /opt/conda/bin/iqtree3 /usr/local/bin/iqtree3
  ln -sf /opt/conda/bin/iqtree3 /opt/conda/bin/iqtree
elif [[ -x /opt/conda/bin/iqtree ]]; then
  ln -sf /opt/conda/bin/iqtree /usr/local/bin/iqtree3
else
  echo "The conda environment does not provide IQ-TREE." >&2
  exit 1
fi
ln -sf /usr/local/bin/iqtree3 /usr/local/bin/iqtree
if [[ ! -e /opt/conda/bin/iqtree3 ]]; then
  ln -sf /usr/local/bin/iqtree3 /opt/conda/bin/iqtree3
fi
rm -f /usr/local/bin/iqtree2 /opt/conda/bin/iqtree2
micromamba run -n base python -m pip check
micromamba run -n base Rscript -e \
  'stopifnot(requireNamespace("kfl1ou", quietly = TRUE), requireNamespace("rkftools", quietly = TRUE))'

#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
jobs="${GG_BUILD_JOBS:-2}"
if [[ ! "${jobs}" =~ ^[1-9][0-9]*$ ]]; then
  echo "GG_BUILD_JOBS must be a positive integer." >&2
  exit 2
fi

usage() {
  echo "Usage: $0 library REPO_URL REVISION ARTIFACT_DIR | worker ARTIFACT_DIR NWKIT_ARTIFACT_DIR" >&2
  exit 2
}

case "${1:-}" in
  library)
    [[ $# -eq 4 ]] || usage
    repo_url="$2"
    revision="$3"
    mkdir -p "$4"
    artifact="$(cd "$4" && pwd)"
    if [[ ! "${revision}" =~ ^[0-9a-f]{40}$ ]]; then
      revision="$(bash "${script_dir}/resolve_git_branch_sha.sh" "${repo_url}" "${revision}")"
    fi
    # Keep matching headers and CMake build metadata until the worker stage.
    # Git archives do not include IQ-TREE's recursive submodule contents.
    mkdir -p "${artifact}/build/source"
    git -C "${artifact}/build/source" init -q
    git -C "${artifact}/build/source" remote add origin "${repo_url}"
    git -C "${artifact}/build/source" fetch --depth 1 origin "${revision}"
    git -C "${artifact}/build/source" checkout -q --detach FETCH_HEAD
    git -C "${artifact}/build/source" submodule update --init --recursive --depth 1
    cmake -S "${artifact}/build/source" -B "${artifact}/build/cli" -DCMAKE_BUILD_TYPE=Release
    cmake --build "${artifact}/build/cli" --parallel "${jobs}"
    cmake -S "${artifact}/build/source" -B "${artifact}/build/library" \
      -DBUILD_LIB=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    cmake --build "${artifact}/build/library" --parallel "${jobs}"
    install -D -m 0755 "${artifact}/build/cli/iqtree3" "${artifact}/rootfs/opt/pg/iqtree3/iqtree3"
    for prefix in usr/local/bin opt/conda/bin; do
      mkdir -p "${artifact}/rootfs/${prefix}"
      ln -s /opt/pg/iqtree3/iqtree3 "${artifact}/rootfs/${prefix}/iqtree3"
      ln -s /opt/pg/iqtree3/iqtree3 "${artifact}/rootfs/${prefix}/iqtree"
    done
    # Ship corresponding upstream source, including submodules, and build
    # materials with the separate GPL program. None enters NWKIT's wheel.
    share="${artifact}/rootfs/usr/local/share/iqtree3"
    mkdir -p "${share}"
    cp "${artifact}/build/source/LICENSE" "${share}/LICENSE"
    tar -czf "${share}/source.tar.gz" --exclude=.git -C "${artifact}/build" source
    cp "${script_dir}/build_iqtree_artifact.sh" "${script_dir}/resolve_git_branch_sha.sh" "${share}/"
    printf 'iqtree\t%s\n' "${revision}" > "${artifact}/source.tsv"
    cp "${artifact}/source.tsv" "${share}/source.tsv"
    ;;
  worker)
    [[ $# -eq 3 ]] || usage
    artifact="$(cd "$2" && pwd)"
    nwkit_artifact="$(cd "$3" && pwd)"
    adapter="${artifact}/build/adapter"
    mkdir -p "${adapter}"
    # Extract exactly the adapter from the same wheel installed at runtime.
    micromamba run -n base python - "${nwkit_artifact}" "${adapter}" <<'PY'
from pathlib import Path
import sys
import zipfile

wheels = list((Path(sys.argv[1]) / "wheels").glob("*.whl"))
if len(wheels) != 1:
    raise SystemExit("Expected exactly one NWKIT wheel for the IQ-TREE adapter")
with zipfile.ZipFile(wheels[0]) as wheel:
    for name in ("nwkit/iqtree_library.py", "nwkit/data_iqtree/worker.cpp"):
        target = Path(sys.argv[2]) / name.removeprefix("nwkit/")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(wheel.read(name))
    licenses = [name for name in wheel.namelist() if name.endswith(".dist-info/licenses/LICENSE")]
    if len(licenses) != 1:
        raise SystemExit("Expected the NWKIT license in its wheel")
    (Path(sys.argv[2]) / "LICENSE.nwkit").write_bytes(wheel.read(licenses[0]))
PY
    micromamba run -n base python "${adapter}/iqtree_library.py" build \
      --build-dir "${artifact}/build/library" --prefix "${artifact}/rootfs/usr/local"
    share="${artifact}/rootfs/usr/local/share/iqtree3"
    cp -a "${adapter}" "${share}/adapter"
    cp "${nwkit_artifact}/source.tsv" "${share}/nwkit-source.tsv"
    mkdir -p "${artifact}/rootfs/opt/pg/logs"
    micromamba run -n base python "${adapter}/iqtree_library.py" check \
      --worker "${artifact}/rootfs/usr/local/bin/nwkit-iqtree-worker" \
      > "${artifact}/rootfs/opt/pg/logs/iqtree3_library_build.json"
    ;;
  *) usage ;;
esac

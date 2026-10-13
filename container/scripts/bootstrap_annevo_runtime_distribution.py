#!/usr/bin/env python3
"""Build fresh official-YAML ANNEVO in an immutable public GG Linux/amd64 image."""
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

OFFICIAL_REPOSITORY = "https://github.com/xjtu-omics/ANNEVO.git"


def run(argv, **kwargs):
    return subprocess.run(argv, check=True, text=True, **kwargs)


def immutable_amd64_reference(image):
    output = run(["docker", "buildx", "imagetools", "inspect", image, "--format", "{{json .Manifest}}"],
                 capture_output=True).stdout
    descriptor = json.loads(output)
    manifests = descriptor.get("manifests")
    if manifests:
        matches = [item for item in manifests if item.get("platform") == {"architecture": "amd64", "os": "linux"}]
        if len(matches) != 1:
            raise ValueError("Public GG base does not expose exactly one Linux/amd64 image")
        digest = matches[0]["digest"]
    else:
        digest = descriptor["digest"]
    if not re.fullmatch(r"sha256:[0-9a-f]{64}", digest):
        raise ValueError("GG base resolution must produce an immutable SHA-256 digest")
    repository = image.split("@", 1)[0]
    if ":" in repository.rsplit("/", 1)[-1]:
        repository = repository.rsplit(":", 1)[0]
    return repository + "@" + digest, digest


def resolve_source(source_ref):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._/-]*", source_ref):
        raise ValueError("Invalid official moving branch")
    data = run(["git", "ls-remote", "--exit-code", OFFICIAL_REPOSITORY, "refs/heads/" + source_ref],
               capture_output=True).stdout.splitlines()
    if len(data) != 1 or not re.fullmatch(r"[0-9a-f]{40}", data[0].split()[0]):
        raise ValueError("Official source branch did not resolve uniquely")
    return data[0].split()[0]


INNER = r'''
set -euo pipefail
mkdir -p /work/home
export HOME=/work/home MAMBA_ROOT_PREFIX=/work/mamba CONDA_PKGS_DIRS=/work/mamba/pkgs
export CUDA_VISIBLE_DEVICES=-1 PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
git clone --no-checkout https://github.com/xjtu-omics/ANNEVO.git /work/source
git -C /work/source checkout --detach "${ANNEVO_RESOLVED_COMMIT}"
test "$(git -C /work/source rev-parse HEAD)" = "${ANNEVO_RESOLVED_COMMIT}"
micromamba create --yes --prefix /work/environment --file /work/source/ANNEVO.yml
micromamba create --yes --prefix /work/pack-tools --override-channels -c conda-forge conda conda-pack
export PATH=/work/pack-tools/bin:${PATH} CONDA_EXE=/work/pack-tools/bin/conda
/opt/conda/bin/python - <<'PY'
import json,os,platform
from pathlib import Path
Path('/work/base-abi.json').write_text(json.dumps({'system':platform.system(),'machine':platform.machine(),
 'libc':{'name':platform.libc_ver()[0],'version':platform.libc_ver()[1]},
 'image_digest':os.environ['GG_BASE_DIGEST'],'image_reference':os.environ['GG_BASE_REFERENCE']},indent=2)+'\n')
PY
/opt/conda/bin/python /repo/container/scripts/build_annevo_runtime_bundle.py \
 --environment /work/environment --source /work/source --source-ref "${ANNEVO_SOURCE_REF}" \
 --base-abi-json /work/base-abi.json --profile linux-amd64-official-yaml \
 --conda-pack /work/pack-tools/bin/conda-pack --output /work/distribution
'''


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", choices=("docker", "apptainer"), default="docker")
    parser.add_argument("--base-image", default="ghcr.io/kfuku52/genegalleon:latest")
    parser.add_argument("--base-sha256", help="Required exact SIF SHA-256 with --engine apptainer")
    parser.add_argument("--apptainer", default="apptainer")
    parser.add_argument("--source-ref", default="main")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--work-parent", type=Path)
    args = parser.parse_args(argv)
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error("Output must be empty; never overwrite distribution bytes")
    if args.engine == "docker":
        image, digest = immutable_amd64_reference(args.base_image)
        run(["docker", "pull", "--platform", "linux/amd64", image])
    else:
        path = Path(args.base_image).resolve()
        if not args.base_sha256 or not re.fullmatch(r"[0-9a-f]{64}", args.base_sha256):
            parser.error("An exact --base-sha256 is required for the qualified GG SIF")
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
                digest.update(block)
        if digest.hexdigest() != args.base_sha256:
            parser.error("Qualified GG SIF SHA-256 changed")
        image, digest = str(path), "sha256:" + args.base_sha256
    commit = resolve_source(args.source_ref)
    repo = Path(__file__).resolve().parents[2]
    output.mkdir(parents=True, exist_ok=True)
    if args.work_parent:
        args.work_parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="annevo-fresh-", dir=args.work_parent))
    try:
        if args.engine == "docker":
            command = ["docker", "run", "--rm", "--platform", "linux/amd64", "--user", f"{os.getuid()}:{os.getgid()}",
                "--mount", f"type=bind,src={repo},dst=/repo,readonly", "--mount", f"type=bind,src={work},dst=/work",
                "--env", "GG_BASE_REFERENCE=" + image, "--env", "GG_BASE_DIGEST=" + digest,
                "--env", "ANNEVO_RESOLVED_COMMIT=" + commit, "--env", "ANNEVO_SOURCE_REF=" + args.source_ref,
                image, "bash", "-c", INNER]
        else:
            command = [args.apptainer, "exec", "--cleanenv", "--containall", "--no-home", "--no-mount", "hostfs,cwd", "--pwd", "/work",
                       "--bind", f"{repo}:/repo:ro,{work}:/work:rw", "--env", "GG_BASE_REFERENCE=" + image,
                       "--env", "GG_BASE_DIGEST=" + digest, "--env", "ANNEVO_RESOLVED_COMMIT=" + commit,
                       "--env", "ANNEVO_SOURCE_REF=" + args.source_ref, "--env", "CUDA_VISIBLE_DEVICES=-1",
                       image, "bash", "-c", INNER]
        run(command)
        distribution = work / "distribution"
        manifests = list(distribution.glob("*.manifest.json"))
        if len(manifests) != 1:
            raise ValueError("Fresh build did not produce one complete manifest")
        for path in sorted(distribution.iterdir()):
            shutil.move(str(path), output / path.name)
    except Exception:
        print("Failed fresh build preserved at " + str(work), file=sys.stderr)
        raise
    else:
        shutil.rmtree(work)
    manifest = next(output.glob("*.manifest.json"))
    print(json.dumps({"base_image": image, "official_source_commit": commit, "manifest": str(manifest),
                      "manifest_sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
                      "GPU_inference_performed": False}))


if __name__ == "__main__":
    main()

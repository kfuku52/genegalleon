#!/usr/bin/env python3
"""Publish immutable OCI runtime bytes; verify anonymously before moving alias."""
import argparse
import copy
import hashlib
import json
import re
import subprocess
from functools import partial
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

MANIFEST_TYPE = "application/vnd.genegalleon.annevo.manifest.v1+json"
PART_TYPE = "application/vnd.genegalleon.annevo.part.v1"
ARTIFACT_TYPE = "application/vnd.genegalleon.annevo.runtime.v1"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def oras(argv, *, check=True, cwd=None, executable="oras", registry_config=None):
    command = [executable, *argv]
    if registry_config:
        command += ["--registry-config", str(registry_config)]
    return subprocess.run(command, capture_output=True, text=True, check=check, cwd=cwd)


def anonymous_reader(repository):
    if not re.fullmatch(r"ghcr\.io/[a-z0-9_.-]+/[a-z0-9_.-]+", repository):
        raise ValueError("Publisher only supports bounded public GHCR repository names")
    namespace = repository.split("/", 1)[1]
    url = "https://ghcr.io/token?" + urlencode({"service": "ghcr.io", "scope": "repository:" + namespace + ":pull"})
    with urlopen(url, timeout=60) as response:
        token = json.load(response)["token"]

    def fetch(kind, digest):
        if not re.fullmatch(r"sha256:[0-9a-f]{64}", digest):
            raise ValueError("Anonymous verification needs a full digest")
        headers = {"Authorization": "Bearer " + token,
                   "Accept": "application/vnd.oci.image.manifest.v1+json"}
        return urlopen(Request("https://ghcr.io/v2/" + namespace + "/" + kind + "/" + digest, headers=headers), timeout=120)

    return fetch


def verify_anonymous(repository, digest, manifest_path, expected, fetch=None):
    fetch = fetch or anonymous_reader(repository)
    with fetch("manifests", digest) as response:
        raw = response.read()
    if "sha256:" + hashlib.sha256(raw).hexdigest() != digest:
        raise ValueError("OCI manifest digest mismatch")
    oci = json.loads(raw)
    if oci.get("artifactType") != ARTIFACT_TYPE:
        raise ValueError("OCI artifact type is not the GeneGalleon ANNEVO runtime")
    layers = oci.get("layers", [])
    metadata = [entry for entry in layers if entry.get("mediaType") == MANIFEST_TYPE]
    parts = [entry for entry in layers if entry.get("mediaType") == PART_TYPE]
    if len(metadata) != 1 or len(parts) != len(expected["parts"]) or len(layers) != len(parts) + 1:
        raise ValueError("OCI artifact does not contain the exact GG metadata/part layers")
    local = manifest_path.read_bytes()
    if metadata[0]["digest"] != "sha256:" + hashlib.sha256(local).hexdigest() or metadata[0]["size"] != len(local):
        raise ValueError("Immutable artifact metadata differs from this exact build")
    with fetch("blobs", metadata[0]["digest"]) as response:
        if response.read() != local:
            raise ValueError("Anonymous metadata bytes mismatch")
    descriptors = {entry["digest"]: entry for entry in parts}
    whole = hashlib.sha256()
    whole_size = 0
    for expected_part in expected["parts"]:
        descriptor = descriptors.get("sha256:" + expected_part["sha256"])
        if not descriptor or descriptor["size"] != expected_part["size"]:
            raise ValueError("OCI part descriptor does not match bundle SHA/size")
        part_sha = hashlib.sha256()
        size = 0
        with fetch("blobs", descriptor["digest"]) as response:
            for block in iter(lambda response=response: response.read(8 * 1024 * 1024), b""):
                size += len(block)
                part_sha.update(block)
                whole.update(block)
        if size != expected_part["size"] or part_sha.hexdigest() != expected_part["sha256"]:
            raise ValueError("Anonymous part bytes mismatch")
        whole_size += size
    if whole_size != expected["bundle"]["size"] or whole.hexdigest() != expected["bundle"]["sha256"]:
        raise ValueError("Anonymous joined archive differs from declared complete bundle")
    # Model downloads are not embedded in the artifact, but public first-use must be real.
    model_checks = {}
    for lineage, model in expected["models"].items():
        url = model["official_url"]
        prefix = "https://raw.githubusercontent.com/xjtu-omics/ANNEVO/" + expected["source"]["commit"] + "/saved_model/"
        if not url.startswith(prefix):
            raise ValueError("Model locator is not bound to the official source commit")
        size = 0
        h = hashlib.sha256()
        with urlopen(url, timeout=120) as response:
            for block in iter(lambda response=response: response.read(8 * 1024 * 1024), b""):
                size += len(block)
                h.update(block)
        if size != model["size"] or h.hexdigest() != model["sha256"]:
            raise ValueError("Anonymous official model SHA/size mismatch")
        model_checks[lineage] = {"size": size, "sha256": h.hexdigest()}
    return {"OCI_digest": digest, "bundle_sha256": whole.hexdigest(), "bundle_size": whole_size,
            "anonymous_parts_verified": len(parts), "anonymous_official_models_verified": model_checks}


def verify_gpu_qualification(expected, proof_path, prepared_path):
    """Bind explicit promotion to the actual cold/warm reserved-GPU qualification."""
    if proof_path is None or prepared_path is None:
        raise ValueError("Promotion requires the GPU qualification proof and prepared receipt")

    def read(path):
        path = Path(path)
        if path.is_symlink() or not path.is_file() or path.stat().st_size > 2 * 1024 * 1024:
            raise ValueError("Invalid GPU qualification evidence file")
        return json.loads(path.read_text())

    proof, selected = read(proof_path), read(prepared_path)
    if (proof.get("qualification_completed") is not True
            or type(proof.get("integration_returncode")) is not int
            or proof["integration_returncode"] != 0
            or any(proof.get(key) is not True for key in (
                "cold_first_install_verified", "warm_no_download_reuse_verified", "cold_warm_gff_bytes_equal"))):
        raise ValueError("GPU qualification did not complete a successful cold/warm comparison")
    observed = selected.get("manifest")
    if not isinstance(observed, dict):
        raise ValueError("Prepared receipt lacks the qualified manifest")
    canonical = (json.dumps(observed, sort_keys=True, separators=(",", ":")) + "\n").encode()
    if proof.get("publisher_manifest_sha256") != hashlib.sha256(canonical).hexdigest():
        raise ValueError("GPU qualification manifest differs from the prepared receipt")
    # The OCI resolver adds public part URLs; all scientific/runtime identities stay exact.
    normalized = copy.deepcopy(observed)
    for part in normalized.get("parts", []):
        part.pop("url", None)
    if normalized != expected:
        raise ValueError("GPU-qualified source/model/parts metadata differs from this candidate")
    model = expected["models"].get(selected.get("model_profile"))
    if (selected.get("runtime_digest") != expected["bundle"]["sha256"]
            or selected.get("source_commit") != expected["source"]["commit"]
            or not model or selected.get("model_sha256") != model["sha256"]
            or selected.get("model_size") != model["size"]):
        raise ValueError("GPU-qualified runtime/source/model selection differs from this candidate")
    probe = proof.get("gpu_probe", {})
    cuda = expected["cuda_profile"]
    if (probe.get("verified") is not True or probe.get("torch_version") != cuda["torch_version"]
            or probe.get("compiled_cuda") != cuda["compiled_cuda"]
            or not re.fullmatch(r"[0-9a-f]{64}", proof.get("annotation", {}).get("sha256", ""))):
        raise ValueError("GPU CUDA profile or successful annotation identity is missing")
    return {"GPU_qualification_proof_sha256": sha256(proof_path),
            "prepared_receipt_sha256": sha256(prepared_path),
            "qualified_bundle_sha256": selected["runtime_digest"],
            "cold_warm_gff_sha256": proof["annotation"]["sha256"]}


def publish(manifest_path, repository, moving_tag, *, oras_executable="oras", registry_config=None,
            promote=False, promote_existing=False, gpu_qualification_proof=None, prepared_receipt=None):
    manifest_path = Path(manifest_path).resolve()
    expected = json.loads(manifest_path.read_text())
    if expected.get("kind") != "genegalleon-annevo-runtime" or expected.get("platform") != "linux/amd64":
        raise ValueError("Publisher requires the Linux/amd64 GG ANNEVO runtime contract")
    if not re.fullmatch(r"annevo-runtime-linux-amd64-[A-Za-z0-9._-]+", moving_tag):
        raise ValueError("Moving artifact tag cannot collide with GeneGalleon image tags")
    qualification = None
    if promote:
        if not promote_existing:
            raise ValueError("Promotion requires an already published candidate and explicit promote-existing mode")
        qualification = verify_gpu_qualification(expected, gpu_qualification_proof, prepared_receipt)
    anonymous_reader(repository)  # verify existing public package without personal credentials
    files = [manifest_path.name + ":" + MANIFEST_TYPE]
    for part in expected["parts"]:
        name = part["name"]
        if Path(name).name != name:
            raise ValueError("Part name must be a plain basename")
        path = manifest_path.parent / name
        if path.stat().st_size != part["size"] or sha256(path) != part["sha256"]:
            raise ValueError("Local part no longer matches this manifest")
        files.append(name + ":" + PART_TYPE)
    immutable_tag = "annevo-runtime-linux-amd64-sha256-" + expected["bundle"]["sha256"]
    immutable = repository + ":" + immutable_tag
    run_oras = partial(oras, executable=oras_executable, registry_config=registry_config)
    existing = run_oras(["resolve", immutable], check=False)
    if existing.returncode == 0:
        digest = existing.stdout.strip()
    else:
        if promote_existing:
            raise ValueError("Qualified immutable artifact must already exist before promotion")
        # Do not turn auth/network failures into an overwrite decision.
        missing = existing.stderr.lower()
        if not any(message in missing for message in ["not found", "manifest_unknown", "404"]):
            raise ValueError("Immutable artifact lookup failed; publication is stopped: " + existing.stderr)
        run_oras(["push", immutable, "--artifact-type", ARTIFACT_TYPE,
              "--annotation", "org.opencontainers.image.source=https://github.com/kfuku52/genegalleon",
              "--annotation", "org.opencontainers.image.licenses=LicenseRef-ANNEVO-NC", *files], cwd=manifest_path.parent)
        digest = run_oras(["resolve", immutable]).stdout.strip()
    proof = verify_anonymous(repository, digest, manifest_path, expected)
    if promote:
        run_oras(["tag", repository + "@" + digest, moving_tag])
        observed = run_oras(["resolve", repository + ":" + moving_tag]).stdout.strip()
        if observed != digest:
            raise ValueError("Moving artifact tag verification failed")
    return {"immutable_ref": "oci://" + immutable, "immutable_digest": digest,
            "moving_ref": "oci://" + repository + ":" + moving_tag if promote else None,
            "moving_alias_updated": promote, "proof": proof, "GPU_qualification": qualification,
            "image_latest_and_sha_tags_changed": False}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--repository", default="ghcr.io/kfuku52/genegalleon")
    parser.add_argument("--moving-tag", default="annevo-runtime-linux-amd64-latest")
    parser.add_argument("--oras", default="oras", help="Private official ORAS executable")
    parser.add_argument("--registry-config", type=Path, help="Private transient credential file; never stored in publication evidence")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--immutable-only", action="store_true", help="Publish candidate bytes without changing the qualified default alias")
    mode.add_argument("--promote-existing", action="store_true", help="After exact candidate GPU qualification, reverify and promote an already published artifact")
    parser.add_argument("--gpu-qualification-proof", type=Path, help="Actual successful reserved-GPU cold/warm qualification exit receipt")
    parser.add_argument("--prepared-receipt", type=Path, help="Cold prepare receipt binding the qualified bundle/source/model")
    args = parser.parse_args(argv)
    print(json.dumps(publish(args.manifest, args.repository, args.moving_tag,
                             oras_executable=args.oras, registry_config=args.registry_config,
                             promote=args.promote_existing, promote_existing=args.promote_existing,
                             gpu_qualification_proof=args.gpu_qualification_proof, prepared_receipt=args.prepared_receipt)))


if __name__ == "__main__":
    main()

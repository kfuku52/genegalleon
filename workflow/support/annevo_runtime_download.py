"""Resolve and fetch optional ANNEVO runtime assets through GG's transport."""
from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import shutil
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlencode, urljoin, urlsplit
from urllib.request import Request

from format_species_download.locking import download_url_to_file
from format_species_network import limited_urlopen

DEFAULT_RUNTIME = "oci://ghcr.io/kfuku52/genegalleon:annevo-runtime-linux-amd64-latest"
MANIFEST_MEDIA = "application/vnd.genegalleon.annevo.manifest.v1+json"
SHA = re.compile(r"[0-9a-f]{64}\Z")


from runtime_asset_io import atomic_json as atomic_json
from runtime_asset_io import canonical_bytes as canonical_bytes
from runtime_asset_io import file_sha as file_sha


def platform_name():
    arch = {"x86_64": "amd64", "aarch64": "arm64"}.get(platform.machine())
    if platform.system() != "Linux" or arch is None:
        raise RuntimeError("ANNEVO runtimes require a qualified Linux CPU architecture")
    return "linux/" + arch


def _read(url, headers=None, maximum=2 * 1024 * 1024):
    with limited_urlopen(Request(url, headers=headers or {}), timeout=120) as response:
        raw = response.read(maximum + 1)
    if len(raw) > maximum:
        raise RuntimeError("ANNEVO metadata exceeds its bounded size")
    return raw


def registry_headers(reference):
    match = re.fullmatch(r"oci://(ghcr\.io)/([a-z0-9/_-]+):([A-Za-z0-9_.-]+)", reference)
    if not match:
        raise ValueError("ANNEVO OCI locator must name a GHCR repository and tag")
    host, repository, tag = match.groups()
    token_url = "https://" + host + "/token?" + urlencode(
        {"service": host, "scope": "repository:" + repository + ":pull"}
    )
    token = json.loads(_read(token_url)).get("token")
    if not isinstance(token, str) or not token:
        raise RuntimeError("Public ANNEVO runtime registry did not issue a pull token")
    return "https://" + host + "/v2/" + repository, tag, {"Authorization": "Bearer " + token}


def resolve_manifest(locator):
    """Return trusted-locator metadata and ephemeral download headers, never tokens on disk."""
    if locator.startswith("oci://"):
        base, tag, headers = registry_headers(locator)
        index_headers = {**headers, "Accept": "application/vnd.oci.image.manifest.v1+json"}
        try:
            index = json.loads(_read(base + "/manifests/" + tag, index_headers))
        except HTTPError as exc:
            raise RuntimeError("ANNEVO runtime artifact is unavailable; use a published manifest or prepopulate the cache") from exc
        descriptors = index.get("layers", [])
        metadata = [d for d in descriptors if d.get("mediaType") == MANIFEST_MEDIA]
        if len(metadata) != 1:
            raise RuntimeError("ANNEVO registry artifact lacks its unique runtime manifest")
        descriptor = metadata[0]
        digest = descriptor.get("digest", "").removeprefix("sha256:")
        if not SHA.fullmatch(digest):
            raise RuntimeError("Invalid ANNEVO registry metadata digest")
        raw = _read(base + "/blobs/sha256:" + digest, headers)
        if hashlib.sha256(raw).hexdigest() != digest or len(raw) != descriptor["size"]:
            raise RuntimeError("ANNEVO registry metadata checksum or size changed")
        manifest = json.loads(raw)
        layers = {d.get("digest", ""): d for d in descriptors}
        for part in manifest["parts"]:
            layer = layers.get("sha256:" + part["sha256"])
            if not layer or layer["size"] != part["size"]:
                raise RuntimeError("ANNEVO runtime part is missing from its registry artifact")
            part["url"] = base + "/blobs/" + layer["digest"]
        return manifest, headers
    parsed = urlsplit(locator)
    if parsed.scheme:
        if parsed.scheme not in {"https", "http", "file"} or parsed.username or parsed.password:
            raise ValueError("Invalid ANNEVO manifest locator")
        raw = _read(locator)
        base = locator
    else:
        path = Path(locator).absolute()
        if path.is_symlink() or not path.is_file() or path.stat().st_size > 2 * 1024 * 1024:
            raise ValueError("Invalid local ANNEVO manifest")
        raw = path.read_bytes()
        base = path.as_uri()
    manifest = json.loads(raw)
    for part in manifest["parts"]:
        part["url"] = urljoin(base, part.get("url", part["name"]))
    return manifest, {}


def validate_manifest(manifest, model_profile):
    if manifest.get("schema_version") != 1 or manifest.get("kind") != "genegalleon-annevo-runtime":
        raise ValueError("Unsupported ANNEVO runtime manifest")
    if manifest.get("format") != "tar.gz" or manifest.get("platform") != platform_name():
        raise ValueError("ANNEVO runtime platform or archive format does not match this container")
    if not re.fullmatch(r"[A-Za-z0-9_]+", model_profile):
        raise ValueError("Invalid ANNEVO model profile")
    source = manifest.get("source", {})
    if source.get("repository", "").removesuffix(".git") != "https://github.com/xjtu-omics/ANNEVO":
        raise ValueError("ANNEVO bundle must record its official source")
    if not re.fullmatch(r"[0-9a-f]{40}", source.get("commit", "")):
        raise ValueError("ANNEVO source revision is missing")
    if not SHA.fullmatch(source.get("native_tree_sha256", "")):
        raise ValueError("ANNEVO source tree identity is missing")
    abi = manifest.get("base_abi", {})
    libc = platform.libc_ver()
    required = abi.get("libc", {})
    def version(value):
        return tuple(int(item) for item in value.split("."))
    if (abi.get("system") != "Linux" or abi.get("machine") != platform.machine()
            or required.get("name") != libc[0]
            or version(required.get("version", "0")) > version(libc[1] or "0")):
        raise ValueError("ANNEVO runtime requires an incompatible base ABI")
    bundle, parts = manifest.get("bundle", {}), manifest.get("parts", [])
    if not SHA.fullmatch(bundle.get("sha256", "")) or not isinstance(bundle.get("size"), int) or bundle["size"] <= 0:
        raise ValueError("Invalid ANNEVO bundle identity")
    unpacked = manifest.get("payload_uncompressed_size")
    if not isinstance(unpacked, int) or not 0 < unpacked <= 16 * 1024 ** 3:
        raise ValueError("ANNEVO payload must record its bounded uncompressed size")
    if not parts or len(parts) > 32 or sum(p.get("size", 0) for p in parts) != bundle["size"]:
        raise ValueError("ANNEVO bundle part sizes do not match its total")
    names = set()
    for part in parts:
        name = part.get("name", "")
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", name) or name in {".", ".."} or name in names:
            raise ValueError("Unsafe or duplicate ANNEVO bundle part name")
        names.add(name)
        parsed = urlsplit(part.get("url", ""))
        if parsed.scheme not in {"https", "http", "file"} or parsed.username or parsed.password:
            raise ValueError("Invalid ANNEVO bundle part URL")
        if not SHA.fullmatch(part.get("sha256", "")) or not isinstance(part.get("size"), int) or part["size"] <= 0:
            raise ValueError("Invalid ANNEVO bundle part identity")
    models = manifest.get("models", {})
    model = models.get(model_profile)
    if not model or not SHA.fullmatch(model.get("sha256", "")) or not isinstance(model.get("size"), int) or model["size"] <= 0:
        raise ValueError("ANNEVO model profile is unavailable or lacks a checksum")
    if model.get("path") != "source/saved_model/ANNEVO_" + model_profile + ".pt":
        raise ValueError("Invalid ANNEVO model cache name")
    expected_url = "https://raw.githubusercontent.com/xjtu-omics/ANNEVO/" + source["commit"] + "/saved_model/ANNEVO_" + model_profile + ".pt"
    if model.get("official_url") != expected_url:
        raise ValueError("ANNEVO model is not bound to its official source revision")
    for key in ("dependency_manifest", "payload_manifest", "license"):
        item = manifest.get(key, {})
        path = Path(item.get("path", ""))
        if not path.parts or path.parts[0] not in {"source", "provenance"} or path.is_absolute() or ".." in path.parts or not SHA.fullmatch(item.get("sha256", "")):
            raise ValueError("Invalid ANNEVO " + key)
    if manifest.get("relocation", {}).get("rename_after_unpack") is not False:
        raise ValueError("ANNEVO environment must retain its final unpack prefix")
    relocation = manifest["relocation"]
    if any(relocation.get(key) != value for key, value in {
        "environment": "environment", "python": "environment/bin/python",
        "unpack": "environment/bin/conda-unpack",
    }.items()):
        raise ValueError("ANNEVO relocation entrypoints changed")
    cuda = manifest.get("cuda_profile", {})
    if set(cuda) != {"torch_version", "compiled_cuda"} or not all(isinstance(v, str) and v for v in cuda.values()):
        raise ValueError("ANNEVO requires a recorded CUDA-enabled Torch profile")
    return manifest


def fetch_verified(url, destination, expected_sha, expected_size, *, headers=None, no_download=False):
    destination = Path(destination)
    if destination.exists():
        if destination.is_symlink() or not destination.is_file():
            raise RuntimeError("Unsafe ANNEVO cached asset: " + str(destination))
        if destination.stat().st_size == expected_size and file_sha(destination) == expected_sha:
            return False
        raise RuntimeError("ANNEVO cached asset checksum changed: " + str(destination))
    if no_download:
        raise RuntimeError("ANNEVO asset is not cached and downloads are disabled")
    destination.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    parsed = urlsplit(url)
    if parsed.scheme == "file":
        from urllib.parse import unquote
        source = Path(unquote(parsed.path))
        if parsed.netloc or source.is_symlink() or not source.is_file():
            raise ValueError("Unsafe local ANNEVO asset")
        temporary = destination.with_name(destination.name + ".part")
        shutil.copyfile(source, temporary)
        if temporary.stat().st_size != expected_size or file_sha(temporary) != expected_sha:
            raise RuntimeError("Local ANNEVO asset checksum mismatch")
        os.replace(temporary, destination)
    else:
        if parsed.scheme not in {"https", "http"} or parsed.username or parsed.password:
            raise ValueError("Invalid ANNEVO asset URL")
        warnings = []
        try:
            download_url_to_file(url, destination, headers or {}, 1800, False, False, 86400, warnings, "ANNEVO optional runtime",
                                 expected_sha256=expected_sha, expected_size=expected_size)
        finally:
            import sys
            for warning in warnings:
                print(warning, file=sys.stderr)
    if destination.stat().st_size != expected_size or file_sha(destination) != expected_sha:
        raise RuntimeError("Downloaded ANNEVO asset checksum mismatch")
    destination.chmod(0o600)
    return True

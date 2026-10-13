"""Install optional ANNEVO generations without modifying the GG base environment."""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import stat
import subprocess
import tarfile
import tempfile
import uuid
from pathlib import Path

from annevo_runtime_download import (
    DEFAULT_RUNTIME,
    atomic_json,
    canonical_bytes,
    fetch_verified,
    file_sha,
    registry_headers,
    resolve_manifest,
    validate_manifest,
)
from shared_namespace_lock import namespace_lock


def private_root(path):
    path = Path(path).absolute()
    for ancestor in (path, *path.parents):
        if ancestor.is_symlink():
            raise RuntimeError("Symlinked ANNEVO cache path: " + str(ancestor))
    created = not path.exists()
    path.mkdir(mode=0o700, parents=True, exist_ok=True)
    if created:
        descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        try:
            os.fchmod(descriptor, 0o700)
        finally:
            os.close(descriptor)
    info = path.stat()
    if info.st_uid != os.getuid() or info.st_mode & 0o022:
        raise RuntimeError("ANNEVO cache must be owned and not writable by other users")
    return path


def signature(path):
    value = path.lstat()
    return [value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns,
            value.st_ctime_ns, stat.S_IMODE(value.st_mode)]


def inspect_tree(root, *, seal=False):
    files = {}
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root).as_posix()
        info = path.lstat()
        if path.is_symlink():
            if not path.resolve().is_relative_to(root):
                raise RuntimeError("ANNEVO runtime symlink escapes its generation: " + relative)
            files[relative] = {"type": "symlink", "target": os.readlink(path)}
        elif path.is_dir():
            if seal and stat.S_IMODE(info.st_mode) != 0o500:
                path.chmod(0o500)
            files[relative] = {"type": "directory"}
        elif path.is_file():
            digest = file_sha(path)
            desired_mode = 0o500 if info.st_mode & 0o111 else 0o400
            if seal and stat.S_IMODE(info.st_mode) != desired_mode:
                path.chmod(desired_mode)
            files[relative] = {"type": "file", "sha256": digest}
        else:
            raise RuntimeError("Unsupported ANNEVO runtime member: " + relative)
        files[relative]["signature"] = signature(path)
    if seal:
        root.chmod(0o500)
    return files


def verify_tree(root, recorded):
    current = {p.relative_to(root).as_posix(): p for p in root.rglob("*")}
    if set(current) != set(recorded):
        raise RuntimeError("ANNEVO runtime file inventory changed")
    for relative, expected in recorded.items():
        path = current[relative]
        if signature(path) != expected["signature"]:
            raise RuntimeError("ANNEVO runtime member changed: " + relative)
        if expected["type"] == "symlink":
            if not path.is_symlink() or os.readlink(path) != expected["target"] or not path.resolve().is_relative_to(root):
                raise RuntimeError("ANNEVO runtime link changed: " + relative)
        elif path.is_symlink():
            raise RuntimeError("ANNEVO runtime member became a symlink: " + relative)


def verify_ready(receipt):
    manifest = receipt["manifest"]
    profile = receipt.get("model_profile") or next(iter(manifest["models"]))
    validate_manifest(manifest, profile)
    if (receipt["runtime_digest"] != manifest["bundle"]["sha256"]
            or receipt["source_commit"] != manifest["source"]["commit"]):
        raise RuntimeError("ANNEVO installed receipt source identity changed")
    cache = private_root(receipt["cache_root"])
    key = hashlib.sha256(canonical_bytes(receipt["manifest"])).hexdigest()
    root = Path(receipt["generation"])
    expected_parent = cache / "runtimes" / key[:32]
    if (receipt["cache_key"] != key or root.parent != expected_parent
            or not root.name.startswith("g-")
            or Path(receipt["inventory_path"]) != expected_parent / (root.name + ".inventory.json")
            or Path(receipt["python"]) != root / "environment/bin/python"
            or Path(receipt["native_entrypoint"]) != root / "source/annotation.py"):
        raise RuntimeError("ANNEVO installed receipt paths or identity changed")
    if root.is_symlink() or root.stat().st_uid != os.getuid() or root.stat().st_mode & 0o222:
        raise RuntimeError("ANNEVO runtime generation is unsafe or writable")
    inventory_path = Path(receipt["inventory_path"])
    if inventory_path.is_symlink() or file_sha(inventory_path) != receipt["inventory_sha256"]:
        raise RuntimeError("ANNEVO installed inventory checksum changed")
    verify_tree(root, json.loads(inventory_path.read_text()))
    return receipt


def _safe_extract(archive, destination, maximum=16 * 1024 ** 3):
    total = 0
    seen = set()
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar:
            path = Path(member.name)
            if (not path.parts or path.is_absolute() or ".." in path.parts
                    or path.parts[0] not in {"environment", "source", "provenance"}
                    or member.size < 0):
                raise RuntimeError(f"Unsafe ANNEVO archive member: name={member.name!r}, normalized={path.as_posix()!r}")
            if path.as_posix() in seen:
                raise RuntimeError(f"Duplicate ANNEVO archive member: name={member.name!r}, normalized={path.as_posix()!r}")
            seen.add(path.as_posix())
            total += member.size
            if total > maximum or len(seen) > 200000 or not (member.isfile() or member.isdir() or member.issym() or member.islnk()):
                raise RuntimeError("ANNEVO archive exceeds its supported payload or member types")
            tar.extract(member, destination, filter="data")


def _check_payload(root, manifest):
    for field in ("payload_manifest", "dependency_manifest", "license"):
        path = root / manifest[field]["path"]
        if path.is_symlink() or not path.is_file() or file_sha(path) != manifest[field]["sha256"]:
            raise RuntimeError("ANNEVO bundle metadata changed: " + field)
    payload_path = root / manifest["payload_manifest"]["path"]
    document = json.loads(payload_path.read_text())
    if document.get("schema_version") != 1:
        raise RuntimeError("Unsupported ANNEVO payload inventory")
    entries = document.get("entries", [])
    declared = {r["path"]: r for r in entries}
    actual = {p.relative_to(root).as_posix(): p for p in root.rglob("*")}
    if len(declared) != len(entries) or set(actual) != set(declared) | {manifest["payload_manifest"]["path"]}:
        raise RuntimeError("ANNEVO bundle payload inventory is incomplete")
    for name, row in declared.items():
        path = actual[name]
        if not path.resolve().is_relative_to(root):
            raise RuntimeError("ANNEVO payload link escapes its environment")
        if row["type"] == "file":
            if path.is_symlink() or not path.is_file() or path.stat().st_size != row["size"] or file_sha(path) != row["sha256"]:
                raise RuntimeError("ANNEVO payload checksum mismatch: " + name)
        elif row["type"] == "symlink":
            if not path.is_symlink() or os.readlink(path) != row["target"]:
                raise RuntimeError("ANNEVO payload symlink mismatch: " + name)
        elif row["type"] == "directory":
            if path.is_symlink() or not path.is_dir():
                raise RuntimeError("ANNEVO payload directory mismatch: " + name)
        elif row["type"] == "hardlink":
            target = root / row["target"]
            if path.is_symlink() or not target.resolve().is_relative_to(root) or not os.path.samefile(path, target):
                raise RuntimeError("ANNEVO payload hardlink mismatch: " + name)
        else:
            raise RuntimeError("Unsupported ANNEVO payload inventory type")
    native = {name.removeprefix("source/"): {"sha256": row["sha256"], "size": row["size"]}
              for name, row in declared.items() if name.startswith("source/") and row["type"] == "file"}
    digest = hashlib.sha256((json.dumps(native, sort_keys=True, indent=2) + "\n").encode()).hexdigest()
    if digest != manifest["source"]["native_tree_sha256"]:
        raise RuntimeError("ANNEVO native source tree identity changed")


def runtime_environment(receipt):
    env = os.environ.copy()
    prefix = Path(receipt["python"]).parent.parent
    env["PATH"] = str(prefix / "bin") + os.pathsep + env.get("PATH", "")
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["PYTHONNOUSERSITE"] = "1"
    env.pop("PYTHONHOME", None)
    env.pop("PYTHONPATH", None)
    return env


def _cpu_check(receipt, manifest):
    code = "import json,torch,numpy,h5py,pandas,numba,torchmetrics; print(json.dumps({'torch_version':torch.__version__,'compiled_cuda':torch.version.cuda}))"
    completed = subprocess.run([receipt["python"], "-E", "-s", "-B", "-c", code],
                               env=runtime_environment(receipt), capture_output=True, text=True, timeout=180)
    if completed.returncode:
        raise RuntimeError("ANNEVO isolated dependency import failed: " + completed.stderr[-4000:])
    result = json.loads(completed.stdout.strip().splitlines()[-1])
    if result != manifest["cuda_profile"]:
        raise RuntimeError("ANNEVO installed CUDA dependency profile changed")
    completed = subprocess.run([receipt["python"], "-E", "-s", "-B", receipt["native_entrypoint"], "--help"],
                               env=runtime_environment(receipt), capture_output=True, text=True, timeout=180)
    if completed.returncode:
        raise RuntimeError("ANNEVO native CLI import failed: " + completed.stderr[-4000:])
    return result


def _install(root, key, manifest, headers, no_download, locator):
    object_root = root / "runtimes" / key[:32]
    private_root(object_root)
    ready_path = object_root / "READY.json"
    with namespace_lock(object_root / ".install.lock", exclusive=True, timeout=86400):
        if ready_path.is_symlink():
            raise RuntimeError("ANNEVO ready receipt is redirected")
        if ready_path.exists():
            ready = json.loads(ready_path.read_text())
            if ready.get("cache_key") != key or ready.get("manifest") != manifest:
                raise RuntimeError("ANNEVO abbreviated cache namespace collision")
            return verify_ready(ready), True
        if no_download:
            raise RuntimeError("ANNEVO environment is not installed and downloads are disabled")
        if not headers and any(p["url"].startswith("https://ghcr.io/") for p in manifest["parts"]):
            if not locator.startswith("oci://"):
                raise RuntimeError("GHCR asset authentication requires its explicit OCI locator")
            _, _, headers = registry_headers(locator)
        parts = []
        for part in manifest["parts"]:
            target = root / "archives" / part["sha256"]
            fetch_verified(part["url"], target, part["sha256"], part["size"], headers=headers)
            parts.append(target)
        needed = manifest["bundle"]["size"] + manifest["payload_uncompressed_size"] + 64 * 1024 ** 2
        if shutil.disk_usage(object_root).free < needed:
            raise RuntimeError("Insufficient free space for ANNEVO archive assembly and environment installation")
        generation = object_root / ("g-" + uuid.uuid4().hex[:16])
        generation.mkdir(mode=0o700)
        with tempfile.TemporaryDirectory(prefix="assemble-", dir=object_root) as assembly:
            archive = Path(assembly) / "runtime.tar.gz"
            digest = hashlib.sha256()
            with archive.open("xb") as output:
                for part in parts:
                    with part.open("rb") as source:
                        for chunk in iter(lambda: source.read(8 * 1024 * 1024), b""):
                            digest.update(chunk)
                            output.write(chunk)
            if digest.hexdigest() != manifest["bundle"]["sha256"] or archive.stat().st_size != manifest["bundle"]["size"]:
                raise RuntimeError("ANNEVO assembled bundle checksum mismatch")
            _safe_extract(archive, generation, manifest["payload_uncompressed_size"])
        _check_payload(generation, manifest)
        relocation = manifest["relocation"]
        python = generation / relocation["python"]
        unpack = generation / relocation["unpack"]
        if python.is_symlink() and not python.resolve().is_relative_to(generation):
            raise RuntimeError("ANNEVO Python points outside its environment")
        receipt = {"schema_version": 1, "kind": "genegalleon-annevo-installed-runtime",
                   "cache_root": str(root), "cache_key": key, "generation": str(generation),
                   "runtime_digest": manifest["bundle"]["sha256"], "manifest": manifest,
                   "source_commit": manifest["source"]["commit"], "python": str(python),
                   "source_dir": str(generation / "source"),
                   "native_entrypoint": str(generation / "source/annotation.py")}
        completed = subprocess.run([str(python), "-E", "-s", "-B", str(unpack)],
                                   env=runtime_environment(receipt), capture_output=True, text=True, timeout=600)
        if completed.returncode:
            raise RuntimeError("ANNEVO environment relocation failed: " + completed.stderr[-4000:])
        receipt["CPU_dependency_check"] = _cpu_check(receipt, manifest)
        inventory = inspect_tree(generation, seal=True)
        inventory_path = object_root / (generation.name + ".inventory.json")
        atomic_json(inventory_path, inventory)
        receipt["inventory_path"] = str(inventory_path)
        receipt["inventory_sha256"] = file_sha(inventory_path)
        atomic_json(ready_path, receipt)
        return receipt, False


def prepare(cache_root, model_profile, locator="", *, receipt_path=None, no_download=False, refresh=False):
    root = private_root(cache_root)
    locator = locator or DEFAULT_RUNTIME
    alias = root / "locators" / (hashlib.sha256(locator.encode()).hexdigest() + ".json")
    headers = {}
    with namespace_lock(root / ".resolve.lock", exclusive=True, timeout=86400):
        if alias.is_symlink():
            raise RuntimeError("ANNEVO cached locator is redirected")
        if alias.exists() and not refresh:
            manifest = json.loads(alias.read_text())["manifest"]
        else:
            if no_download:
                raise RuntimeError("ANNEVO runtime manifest is not cached and downloads are disabled")
            manifest, headers = resolve_manifest(locator)
            validate_manifest(manifest, model_profile)
            atomic_json(alias, {"locator": locator, "manifest": manifest})
    validate_manifest(manifest, model_profile)
    key = hashlib.sha256(canonical_bytes(manifest)).hexdigest()
    installed, reused = _install(root, key, manifest, headers, no_download, locator)
    model = manifest["models"][model_profile]
    model_path = root / "models" / model["sha256"] / Path(model["path"]).name
    with namespace_lock(root / (".model-" + model["sha256"] + ".lock"), exclusive=True, timeout=86400):
        fetch_verified(model["official_url"], model_path, model["sha256"], model["size"], no_download=no_download)
    selected = {**installed, "model_profile": model_profile, "model_path": str(model_path),
                "model_sha256": model["sha256"], "model_size": model["size"]}
    destination = Path(receipt_path) if receipt_path else root / "selected" / (model_profile + ".json")
    if destination.is_symlink():
        raise RuntimeError("ANNEVO selected receipt is redirected")
    raw = canonical_bytes(selected)
    if not destination.exists() or destination.read_bytes() != raw:
        atomic_json(destination, selected)
    return {**selected, "reused": reused, "receipt_path": str(destination),
            "ready_receipt_sha256": file_sha(destination)}

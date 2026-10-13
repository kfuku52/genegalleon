#!/usr/bin/env python3
"""Build a relocatable, source-bound ANNEVO environment artifact without weights.

Build inside the target GeneGalleon Linux/architecture/base-ABI container. The
input environment is read-only to this producer. conda-pack's strict checks are
required; a publisher supplies real public locators only after publication.
"""
import argparse
import copy
import gzip
import hashlib
import io
import json
import os
import posixpath
import re
import shutil
import stat
import subprocess
import tarfile
import tempfile
from pathlib import Path, PurePosixPath
from urllib.parse import urlsplit

OFFICIAL_REPOSITORY = "https://github.com/xjtu-omics/ANNEVO"
SCHEMA_VERSION = 1
DEFAULT_PART_BYTES = 1024 ** 3
ARCHITECTURES = {"x86_64": "amd64", "aarch64": "arm64"}
SHARED_OPENMP_FILES = {"lib/libiomp5.so", "lib/libomptarget.so"}


class BundleError(ValueError):
    pass


def sha256(path):
    path = Path(path)
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda stream=stream: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    after = path.stat()
    if (before.st_ino, before.st_size, before.st_mtime_ns) != (after.st_ino, after.st_size, after.st_mtime_ns):
        raise BundleError(f"Input changed during hashing: {path}")
    return digest.hexdigest()


def json_bytes(value):
    return (json.dumps(value, sort_keys=True, indent=2) + "\n").encode()


def checked_run(argv, **kwargs):
    result = subprocess.run(argv, check=False, capture_output=True, text=True, **kwargs)
    if result.returncode:
        raise BundleError(f"Command failed ({result.returncode}): {argv!r}\n{result.stderr}")
    return result.stdout.strip()


def safe_name(name):
    name = name.removeprefix("./").rstrip("/")
    path = PurePosixPath(name)
    if not name or not path.parts or path.is_absolute() or ".." in path.parts or "\x00" in name or "\\" in name:
        raise BundleError(f"Unsafe archive/source path: {name!r}")
    return str(path)


def inspect_source(source, source_ref):
    source = Path(source).resolve()
    origin = checked_run(["git", "-C", str(source), "remote", "get-url", "origin"])
    if origin.removesuffix(".git") != OFFICIAL_REPOSITORY:
        raise BundleError("Source origin must be the public official ANNEVO repository")
    if checked_run(["git", "-C", str(source), "status", "--porcelain"]):
        raise BundleError("Official source checkout must have no tracked or untracked changes")
    commit = checked_run(["git", "-C", str(source), "rev-parse", "HEAD"])
    if not re.fullmatch(r"[0-9a-f]{40}|[0-9a-f]{64}", commit):
        raise BundleError("Source commit is not a complete Git revision")
    tracked = checked_run(["git", "-C", str(source), "ls-files", "-z"]).split("\x00")
    tracked = [safe_name(name) for name in tracked if name]
    entrypoints = {"annotation.py", "prediction.py", "decoding.py"}
    notices = {name for name in tracked if "/" not in name and name.upper().startswith(("LICENSE", "COPYING", "NOTICE"))}
    if "LICENSE" not in notices or not entrypoints <= set(tracked):
        raise BundleError("Official source is missing its native entrypoints or LICENSE")
    selected = sorted(entrypoints | notices | {name for name in tracked if name.startswith("src/") and name.endswith(".py")})
    paths = {}
    for name in selected:
        path = source / name
        if path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(source):
            raise BundleError(f"Native source must be a regular contained file: {name}")
        paths[name] = {"sha256": sha256(path), "size": path.stat().st_size}
    models = {}
    for name in sorted(tracked):
        match = re.fullmatch(r"saved_model/ANNEVO_([A-Za-z0-9_]+)\.pt", name)
        if not match:
            continue
        path = source / name
        if path.is_symlink() or not path.is_file():
            raise BundleError(f"Official model is not a regular file: {name}")
        blob_size = int(checked_run(["git", "-C", str(source), "cat-file", "-s", f"HEAD:{name}"]))
        blob = checked_run(["git", "-C", str(source), "show", f"HEAD:{name}"], encoding="utf-8", errors="replace") if blob_size < 1024 else ""
        if blob.startswith("version https://git-lfs.github.com/spec/v1"):
            raise BundleError(f"Git LFS model contents must be materialized before packaging metadata: {name}")
        lineage = match.group(1)
        models[lineage] = {"official_url": f"https://raw.githubusercontent.com/xjtu-omics/ANNEVO/{commit}/{name}",
                           "path": "source/" + name, "sha256": sha256(path), "size": path.stat().st_size}
    if not models:
        raise BundleError("No materialized official model metadata was found")
    tree_digest = hashlib.sha256(json_bytes(paths)).hexdigest()
    return {"repository": OFFICIAL_REPOSITORY, "ref": source_ref, "commit": commit,
            "native_tree_sha256": tree_digest, "files": paths}, models


ENVIRONMENT_PROBE = r'''
import importlib.metadata as metadata,json,platform,sys,torch
assert not torch.cuda.is_available(), 'builder metadata probe must be CPU-hidden'
print(json.dumps({'system':platform.system(),'machine':platform.machine(),
 'libc':{'name':platform.libc_ver()[0],'version':platform.libc_ver()[1]},
 'python':platform.python_version(),'python_prefix':sys.prefix,
 'torch_version':torch.__version__,'compiled_cuda':torch.version.cuda,
 'packages':sorted([{'name':d.metadata['Name'],'version':d.version} for d in metadata.distributions()],key=lambda d:d['name'].lower())}))
'''


def inspect_environment(prefix, source, base_abi, staging):
    prefix = Path(prefix).resolve()
    if not (prefix / "conda-meta").is_dir() or not (prefix / "bin/python").is_file():
        raise BundleError("Input must be a complete Conda environment")
    env = {key: value for key, value in os.environ.items() if key not in {"PYTHONPATH", "PYTHONHOME", "LD_PRELOAD"}}
    env.update({"CUDA_VISIBLE_DEVICES": "-1", "PYTHONDONTWRITEBYTECODE": "1", "PYTHONNOUSERSITE": "1",
                "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                "NUMBA_CACHE_DIR": str(staging / "numba-cache"), "TORCH_HOME": str(staging / "torch-cache"),
                "XDG_CACHE_HOME": str(staging / "xdg-cache")})
    python = [str(prefix / "bin/python"), "-E", "-s", "-B"]
    observed = json.loads(checked_run(python + ["-c", ENVIRONMENT_PROBE], env=env, cwd=staging, timeout=120))
    expected = {key: base_abi.get(key) for key in ("system", "machine", "libc")}
    if expected != {key: observed[key] for key in expected} or observed["system"] != "Linux":
        raise BundleError(f"Environment probe does not match target GG Linux/base ABI: {expected!r}, {observed!r}")
    if observed["machine"] not in ARCHITECTURES or not observed["compiled_cuda"]:
        raise BundleError("A supported Linux architecture and CUDA-built torch are required; no CPU fallback")
    if Path(observed["python_prefix"]).resolve() != prefix:
        raise BundleError("Metadata probe escaped its dedicated Conda prefix")
    if not re.fullmatch(r"sha256:[0-9a-f]{64}", base_abi.get("image_digest", "")):
        raise BundleError("base_abi.image_digest must identify the qualified GG base image")
    for name in ("annotation.py", "prediction.py", "decoding.py"):
        checked_run(python + [str(Path(source).resolve() / name), "--help"], env=env, cwd=staging, timeout=120)
    packages = []
    for path in sorted((prefix / "conda-meta").glob("*.json")):
        record = json.loads(path.read_text())
        package = {key: record.get(key) for key in ("name", "version", "build", "channel", "url", "sha256", "md5")}
        for key in ("url", "channel"):
            value = package[key]
            if isinstance(value, str) and (urlsplit(value).username or urlsplit(value).password or "/t/" in value):
                raise BundleError("Credential-bearing package locators cannot enter a public bundle")
        packages.append(package)
    observed["conda_packages"] = packages
    return observed


def preserved_shared_files(prefix):
    """Preserve the effective prefix for the two demonstrated OpenMP collisions."""
    prefix = Path(prefix).resolve()
    owners = {}
    for path in sorted((prefix / "conda-meta").glob("*.json")):
        metadata = json.loads(path.read_text())
        for name in metadata.get("files", []):
            owners.setdefault(name, []).append((path, metadata))
    result = {}
    for name, providers in owners.items():
        if len(providers) < 2:
            continue
        if (name not in SHARED_OPENMP_FILES or len(providers) != 2
                or {item[1].get("name") for item in providers} != {"intel-openmp", "llvm-openmp"}):
            raise BundleError("Unrecognized managed-file overlap: " + name)
        path = prefix / name
        details = path.lstat()
        mode = stat.S_IMODE(details.st_mode)
        if not stat.S_ISREG(details.st_mode) or mode not in {0o644, 0o664, 0o755, 0o775} or details.st_size > 32 * 1024 ** 2:
            raise BundleError("Unrecognized shared-library physical type/mode/size: " + name)
        data = path.read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        proof = []
        for metadata_path, metadata in providers:
            entries = [entry for entry in metadata.get("paths_data", {}).get("paths", []) if entry.get("_path") == name]
            if len(entries) != 1:
                raise BundleError("Shared-library provider lacks exact path metadata: " + name)
            entry = entries[0]
            expected_type = "hardlink" if metadata["name"] == "intel-openmp" else "softlink"
            if (entry.get("path_type") != expected_type or entry.get("file_mode") is not None
                    or entry.get("prefix_placeholder") or not re.fullmatch(r"[0-9a-f]{64}", entry.get("sha256_in_prefix", ""))):
                raise BundleError("Unrecognized shared-library provider type or prefix placeholder: " + name)
            proof.append({"name": metadata["name"], "version": metadata["version"], "build": metadata["build"],
                          "metadata_sha256": sha256(metadata_path), "path_record": entry})
        if digest not in {provider["path_record"]["sha256_in_prefix"] for provider in proof}:
            raise BundleError("Physical shared-library SHA does not match a managed provider: " + name)
        if (str(prefix).encode() in data or b"_h_env_placehold" in data
                or b"/opt/anaconda1anaconda2anaconda3" in data):
            raise BundleError("Retained shared-library contains a prefix placeholder: " + name)
        if sha256(path) != digest or stat.S_IMODE(path.lstat().st_mode) != mode:
            raise BundleError("Physical shared-library changed during qualification: " + name)
        result[name] = {"sha256": digest, "size": len(data), "mode": mode, "providers": proof,
                        "matches_providers": [p["name"] for p in proof if p["path_record"]["sha256_in_prefix"] == digest]}
    return result


def pack_environment(prefix, output, conda_pack, exclude_files=()):
    argv = [conda_pack, "--prefix", str(Path(prefix).resolve()), "--output", str(output),
            "--format", "tar.gz", "--n-threads", "1"]
    for name in sorted(exclude_files):
        if name not in SHARED_OPENMP_FILES:
            raise BundleError("Only the demonstrated managed OpenMP collisions can be excluded")
        argv += ["--exclude", name]
    version = checked_run([conda_pack, "--version"], timeout=60)
    # Conda-pack creates a temporary archive before atomic publication. Keep it
    # beside the output on the qualified storage, away from contained /tmp RAM.
    checked_run(argv, env={**os.environ, "TMPDIR": str(Path(output).resolve().parent)}, timeout=3600)
    return {"argv": argv, "version": version, "strict_missing_files": True,
            "strict_editable_packages": True, "fixed_destination_prefix": False}


def environment_members(archive):
    seen = set()
    for original in archive:
        if original.name.rstrip("/") in {"", "."}:
            continue
        name = safe_name(original.name)
        if name in seen:
            raise BundleError(f"Duplicate environment archive path: {name}")
        seen.add(name)
        member = copy.copy(original)
        # TarInfo has already applied the input PAX overrides. Reusing them
        # would override our rewritten name/linkname and normalized metadata
        # when the output is read. The writer generates new long-path/link
        # headers from the validated attributes instead.
        member.pax_headers = {}
        member.name = "environment/" + name
        if member.issym():
            target = posixpath.normpath(posixpath.join(posixpath.dirname(name), member.linkname))
            if member.linkname.startswith("/") or target == ".." or target.startswith("../"):
                raise BundleError(f"Environment symlink escapes its prefix: {name}")
        elif member.islnk():
            member.linkname = "environment/" + safe_name(member.linkname)
        elif not (member.isfile() or member.isdir()):
            raise BundleError(f"Special file in environment archive: {name}")
        yield member, original
    if "bin/conda-unpack" not in seen or "bin/python" not in seen:
        raise BundleError("Relocatable Conda archive must contain bin/python and bin/conda-unpack")


def add_bytes(archive, name, data, mode=0o644):
    member = tarfile.TarInfo(name)
    member.size = len(data)
    member.mode = mode
    archive.addfile(member, io.BytesIO(data))
    return {"path": name, "type": "file", "sha256": hashlib.sha256(data).hexdigest(), "size": len(data), "mode": member.mode}


def write_payload(env_archive, source, source_info, dependencies, pack_info, output, prefix=None, shared_files=None):
    inventory = []
    emitted = {}
    with output.open("wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        with tarfile.open(fileobj=gz, mode="w|") as target, tarfile.open(env_archive, "r:gz") as environment:
            def ensure_parents(name):
                for parent in reversed(PurePosixPath(name).parents):
                    path = str(parent)
                    if path == ".":
                        continue
                    if path in emitted:
                        if emitted[path] != "directory":
                            raise BundleError("Archive parent is not a physical directory: " + path)
                        continue
                    directory = tarfile.TarInfo(path)
                    directory.type = tarfile.DIRTYPE
                    directory.mode = 0o755
                    target.addfile(directory)
                    inventory.append({"path": path, "type": "directory", "mode": 0o755})
                    emitted[path] = "directory"

            def add_record(name, data, mode=0o644):
                ensure_parents(name)
                record = add_bytes(target, name, data, mode)
                emitted[name] = "file"
                return record

            for name, proof in sorted((shared_files or {}).items()):
                path = Path(prefix) / name
                if (path.is_symlink() or not path.is_file() or stat.S_IMODE(path.stat().st_mode) != proof["mode"]
                        or path.stat().st_size != proof["size"]):
                    raise BundleError("Retained shared-library input identity changed: " + name)
                data = path.read_bytes()
                if hashlib.sha256(data).hexdigest() != proof["sha256"]:
                    raise BundleError("Retained shared-library input SHA changed: " + name)
                inventory.append(add_record("environment/" + name, data, 0o755 if proof["mode"] & 0o111 else 0o644))

            for member, original in environment_members(environment):
                ensure_parents(member.name)
                if member.isdir() and emitted.get(member.name) == "directory":
                    continue
                if member.name in emitted:
                    raise BundleError("Duplicate emitted payload path")
                if member.islnk() and emitted.get(member.linkname) not in {"file", "hardlink"}:
                    raise BundleError("Hardlink target must precede its link in the packed environment")
                member.uid = member.gid = member.mtime = 0
                member.uname = member.gname = ""
                member.mode = 0o755 if (member.isdir() or member.mode & 0o111) else 0o644
                record = {"path": member.name, "mode": member.mode}
                if member.isfile():
                    stream = environment.extractfile(original)
                    if stream is None:
                        raise BundleError("Missing packed member bytes")
                    # Spool only one member; CUDA shared libraries can exceed memory limits.
                    with tempfile.TemporaryFile(dir=Path(output).resolve().parent) as spool:
                        digest = hashlib.sha256()
                        for block in iter(lambda stream=stream: stream.read(8 * 1024 * 1024), b""):
                            digest.update(block)
                            spool.write(block)
                        if spool.tell() != member.size:
                            raise BundleError("Packed member size mismatch")
                        spool.seek(0)
                        target.addfile(member, spool)
                    record.update(type="file", sha256=digest.hexdigest(), size=member.size)
                else:
                    target.addfile(member)
                    record.update(type="directory" if member.isdir() else ("symlink" if member.issym() else "hardlink"))
                    if member.issym() or member.islnk():
                        record["target"] = member.linkname
                inventory.append(record)
                emitted[member.name] = record["type"]
            for name, expected in source_info["files"].items():
                path = Path(source) / name
                if sha256(path) != expected["sha256"] or path.stat().st_size != expected["size"]:
                    raise BundleError("Official source changed during packing")
                inventory.append(add_record("source/" + name, path.read_bytes()))
            dependency_record = add_record("provenance/dependencies.json", json_bytes(dependencies))
            inventory.append(dependency_record)
            inventory.append(add_record("provenance/source.json", json_bytes(source_info)))
            inventory.append(add_record("provenance/conda-pack.json", json_bytes(pack_info)))
            payload_record = add_record("provenance/payload-manifest.json", json_bytes({"schema_version": 1, "entries": inventory}))
    uncompressed_size = sum(record.get("size", 0) for record in inventory) + payload_record["size"]
    return dependency_record, payload_record, uncompressed_size


def verify_serialized_payload(archive, payload_record, maximum):
    """Read actual tar headers/bytes back before any distribution parts exist."""
    records = {}
    total = 0
    document = None
    with tarfile.open(archive, "r:gz") as packed:
        for member in packed:
            name = safe_name(member.name)
            if (name != member.name or PurePosixPath(name).parts[0] not in {"environment", "source", "provenance"}
                    or name in records or member.size < 0):
                raise BundleError(f"Unsafe or duplicate serialized payload member: {member.name!r}")
            total += member.size
            if total > maximum or len(records) >= 200000:
                raise BundleError("Serialized payload exceeds its declared bounds")
            if (member.uid, member.gid, member.mtime, member.uname, member.gname) != (0, 0, 0, "", ""):
                raise BundleError("Serialized payload metadata was not normalized: " + name)
            row = {"path": name, "mode": member.mode}
            if member.isfile():
                stream = packed.extractfile(member)
                digest = hashlib.sha256()
                size = 0
                inventory_bytes = bytearray() if name == payload_record["path"] else None
                if inventory_bytes is not None and member.size > 64 * 1024 ** 2:
                    raise BundleError("Payload inventory exceeds its supported size")
                for block in iter(lambda stream=stream: stream.read(8 * 1024 * 1024), b""):
                    digest.update(block)
                    size += len(block)
                    if inventory_bytes is not None:
                        inventory_bytes.extend(block)
                if size != member.size:
                    raise BundleError("Serialized payload file size changed: " + name)
                row.update(type="file", sha256=digest.hexdigest(), size=size)
                if inventory_bytes is not None:
                    if row != payload_record:
                        raise BundleError("Serialized payload inventory identity changed")
                    document = json.loads(inventory_bytes)
            elif member.isdir():
                row.update(type="directory")
            elif member.issym() or member.islnk():
                if member.issym():
                    target = posixpath.normpath(posixpath.join(posixpath.dirname(name), member.linkname))
                    parts = PurePosixPath(target).parts
                    if member.linkname.startswith("/") or not parts or parts[0] != "environment":
                        raise BundleError("Serialized payload symlink escapes its environment: " + name)
                    row.update(type="symlink", target=member.linkname)
                else:
                    target = safe_name(member.linkname)
                    if target != member.linkname or records.get(target, {}).get("type") not in {"file", "hardlink"}:
                        raise BundleError("Serialized payload hardlink target changed: " + name)
                    row.update(type="hardlink", target=member.linkname)
            else:
                raise BundleError("Unsupported serialized payload member: " + name)
            records[name] = row
    if total != maximum or not isinstance(document, dict) or document.get("schema_version") != 1:
        raise BundleError("Serialized payload inventory or total size changed")
    entries = document.get("entries", [])
    declared = {row["path"]: row for row in entries}
    actual = {name: row for name, row in records.items() if name != payload_record["path"]}
    if len(declared) != len(entries) or declared != actual:
        raise BundleError("Serialized payload differs from its complete inventory")


def split_archive(archive, destination, part_bytes):
    if not isinstance(part_bytes, int) or not 1 <= part_bytes < 2 * 1024 ** 3:
        raise BundleError("Each public part must be positive and strictly smaller than 2 GiB")
    parts = []
    with archive.open("rb") as stream:
        index = 0
        while True:
            block = stream.read(min(part_bytes, 8 * 1024 * 1024))
            if not block:
                break
            name = archive.name + f".part{index:04d}"
            path = destination / name
            digest = hashlib.sha256()
            size = 0
            with path.open("xb") as output:
                while block:
                    output.write(block)
                    digest.update(block)
                    size += len(block)
                    if size == part_bytes:
                        break
                    block = stream.read(min(part_bytes - size, 8 * 1024 * 1024))
            parts.append({"name": name, "sha256": digest.hexdigest(), "size": size})
            index += 1
    return parts


def build_bundle(prefix, source, output, base_abi, profile, source_ref="main", conda_pack="conda-pack", part_bytes=DEFAULT_PART_BYTES):
    output = Path(output).resolve()
    if output.is_relative_to(Path(prefix).resolve()) or output.is_relative_to(Path(source).resolve()):
        raise BundleError("Output must remain outside the read-only source and environment")
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise BundleError("Output directory must be empty; existing published bundle bytes are never overwritten")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}", profile):
        raise BundleError("Profile must be a bounded identifier")
    if not isinstance(source_ref, str) or not source_ref.strip() or len(source_ref) > 256:
        raise BundleError("Source ref must be a bounded nonempty build metadata string")
    source_info, models = inspect_source(source, source_ref)
    with tempfile.TemporaryDirectory(prefix=".annevo-build-", dir=output.parent) as directory:
        staging = Path(directory)
        dependencies = inspect_environment(prefix, source, base_abi, staging)
        shared_files = preserved_shared_files(prefix)
        env_archive = staging / "environment.tar.gz"
        pack_info = (pack_environment(prefix, env_archive, conda_pack, exclude_files=shared_files)
                     if shared_files else pack_environment(prefix, env_archive, conda_pack))
        pack_info["effective_prefix_shared_files"] = shared_files
        payload = staging / "payload.tar.gz"
        dependency_record, payload_record, uncompressed_size = write_payload(
            env_archive, source, source_info, dependencies, pack_info, payload, prefix=prefix, shared_files=shared_files)
        verify_serialized_payload(payload, payload_record, uncompressed_size)
        digest = sha256(payload)
        name = f"annevo-runtime-{digest}.tar.gz"
        payload.rename(staging / name)
        archive = staging / name
        # Partial construction stays outside the publication directory until all checks finish.
        parts = split_archive(archive, staging, part_bytes)
        license_record = source_info["files"]["LICENSE"]
        manifest = {"schema_version": SCHEMA_VERSION, "kind": "genegalleon-annevo-runtime", "format": "tar.gz",
                    "profile": profile, "platform": "linux/" + ARCHITECTURES[dependencies["machine"]],
                    "base_abi": base_abi, "source": {key: value for key, value in source_info.items() if key != "files"},
                    "models": models, "cuda_profile": {"torch_version": dependencies["torch_version"], "compiled_cuda": dependencies["compiled_cuda"]},
                    "bundle": {"name": name, "sha256": digest, "size": archive.stat().st_size}, "parts": parts,
                    "payload_uncompressed_size": uncompressed_size,
                    "dependency_manifest": dependency_record, "payload_manifest": payload_record,
                    "license": {"path": "source/LICENSE", **license_record},
                    "relocation": {"environment": "environment", "python": "environment/bin/python",
                                   "unpack": "environment/bin/conda-unpack", "rename_after_unpack": False},
                    "qualification": {"CPU_hidden_metadata_native_help": "passed", "serialized_payload_readback": "passed",
                                      "GPU_inference": "not_performed"}}
        for part in parts:
            shutil.move(str(staging / part["name"]), output / part["name"])
        manifest_path = output / (name + ".manifest.json")
        manifest_path.write_bytes(json_bytes(manifest))
        (output / (manifest_path.name + ".sha256")).write_text(sha256(manifest_path) + "  " + manifest_path.name + "\n")
    return manifest_path, manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--environment", required=True, type=Path)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--source-ref", default="main", help="Moving official branch resolved by the caller; one-off overrides are build metadata")
    parser.add_argument("--base-abi-json", required=True, type=Path)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--conda-pack", default="conda-pack", help="Builder-side executable; never install packaging tools into the input runtime")
    parser.add_argument("--part-bytes", type=int, default=DEFAULT_PART_BYTES)
    args = parser.parse_args(argv)
    try:
        path, manifest = build_bundle(args.environment, args.source, args.output, json.loads(args.base_abi_json.read_text()),
                                      args.profile, args.source_ref, args.conda_pack, args.part_bytes)
    except (BundleError, OSError, json.JSONDecodeError) as error:
        parser.exit(1, str(error) + "\n")
    print(json.dumps({"manifest": str(path), "manifest_sha256": sha256(path), "bundle": manifest["bundle"], "parts": len(manifest["parts"])}))


if __name__ == "__main__":
    main()

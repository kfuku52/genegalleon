"""Fresh image identities and conservative version-inventory cache contexts."""

import argparse
import contextlib
import errno
import fcntl
import hashlib
import json
import os
import platform
import shlex
import shutil
import stat
import struct
import sys
from pathlib import Path

from input_generation_array_state import _stat_identity, atomic_json, digest
from metadata_store import DATABASE, cache_policy, pack, read_bytes
from metadata_store import read_snapshot as metadata_snapshot
from performance_metrics import measure
from shared_namespace_lock import namespace_lock


def image_identity(path):
    path = Path(path)
    if sys.platform == "linux":
        with path.open("rb") as handle:
            before = os.fstat(handle.fileno())
            if not stat.S_ISREG(before.st_mode):
                raise ValueError("Version identity requires a regular image")
            # Linux UAPI: _IOWR('f', 134, struct fsverity_digest), header size 4.
            # Only measure already-enabled verity; never enable it or mutate SIFs.
            buffer = bytearray(struct.pack("=HH", 0, 64) + bytes(64))
            try:
                fcntl.ioctl(handle.fileno(), 0xC0046686, buffer, True)
            except OSError as exc:
                if exc.errno not in (errno.ENODATA, errno.ENOTTY, errno.EOPNOTSUPP):
                    raise
            else:
                algorithm, size = struct.unpack("=HH", buffer[:4])
                if (algorithm, size) not in ((1, 32), (2, 64)):
                    raise ValueError("Unsupported kernel verity digest")
                if _stat_identity(before) != _stat_identity(path.stat()):
                    raise ValueError("Image replaced while measuring verity")
                return "fsverity:" + str(algorithm) + ":" + buffer[4:4 + size].hex()
    return "sha256:" + digest(path)


def _environment_context(environ):
    diagnostic = {"PWD", "OLDPWD", "SHLVL", "_", "GG_PERFORMANCE_DIR", "GG_JOB_ID", "GG_ARRAY_TASK_ID",
                  "JOB_ID", "GG_RESOURCE_OWNER_PID", "GG_ARRAY_JOB_ID", "GG_TASK_CPUS", "GG_MEM_TOTAL_GB",
                  "GG_OBSERVATION_ATTEMPT_DIR", "GG_RESOURCE_EVENT_FILE", "GG_RESOURCE_TIMES_FILE",
                  "GG_RESOURCE_PROFILE"}
    context = {}
    for key, value in environ.items():
        suffix = key.removeprefix("SINGULARITYENV_").removeprefix("APPTAINERENV_")
        if suffix in diagnostic or suffix.startswith(("SLURM_", "SGE_")):
            continue
        context[key] = value
    return context


def environment_identity():
    # Preserve the public legacy identity; values may contain credentials.
    return hashlib.sha256(json.dumps(_environment_context(os.environ), sort_keys=True).encode()).hexdigest()


# Only this reviewed inventory script treats managed scratch as disposable
# storage rather than a software input. Unknown/custom scripts keep literal
# contexts. Changes to the script require reviewing this contract again.
MANAGED_SCRATCH_INVENTORY_SHA256 = "cf83779ddc3e0b509517196b88908905ed925970efea31be538a24e8854b328d"
BIND_CONTEXT_FIELDS = {"GG_CONTAINER_BIND_MOUNTS", "SINGULARITY_BIND", "SINGULARITY_BINDPATH",
                       "APPTAINER_BIND", "APPTAINER_BINDPATH"}
TEMP_CONTEXT_FIELDS = {"TMP", "TMPDIR", "TEMP", "APPTAINER_TMPDIR", "SINGULARITY_TMPDIR",
                       "PYTHONPYCACHEPREFIX"}


def _suffix(key):
    return key.removeprefix("SINGULARITYENV_").removeprefix("APPTAINERENV_")


def _scratch_binds(value, scratch, mount):
    entries = []
    matched = False
    for entry in value.split(","):
        fields = entry.split(":")
        # Apptainer's inherited bind list also uses source-only entries, whose
        # destination is the same path. Keep those literal; only the explicit
        # managed source:destination pair below may lose its scratch name.
        if len(fields) not in (1, 2, 3) or not fields[0]:
            return None
        if fields[:2] == [scratch, mount]:
            entries.append([{"managed_task_scratch": True}, mount, *fields[2:]])
            matched = True
        elif scratch in entry or mount in entry:
            # A scratch overlay of software, or an extra mount within scratch,
            # needs its literal identity, even when another bind is recognised.
            return None
        else:
            entries.append(fields)
    return entries, matched


def _managed_scratch_context(command, context, environ):
    scratch = environ.get("GG_TMP_HOST_ROOT", "")
    mount = environ.get("GG_TMP_MOUNT", "")
    if (not Path(scratch).is_absolute() or os.path.normpath(scratch) != scratch
            or mount != "/gg_tmp" or environ.get("GG_COMMON_TMP_ROOT") != scratch
            or not environ.get("GG_TMP_WORKSPACE_ID")):
        return None
    try:
        argv = shlex.split(command)
    except ValueError:
        return None
    if len(argv) < 2 or argv[1] != "exec":
        return None
    normal = dict(context)
    matched = False
    for key, value in list(normal.items()):
        suffix = _suffix(key)
        if suffix in {"GG_TMP_HOST_ROOT", "GG_COMMON_TMP_ROOT", "GG_TMP_MOUNT"}:
            expected = mount if suffix == "GG_TMP_MOUNT" else scratch
            if value != expected:
                return None
            del normal[key]
        elif suffix in TEMP_CONTEXT_FIELDS and any(
                value == root or value.startswith(root + "/") for root in (scratch, mount)):
            if ".." in Path(value).parts:
                return None
            del normal[key]
        elif suffix in BIND_CONTEXT_FIELDS and value:
            parsed = _scratch_binds(value, scratch, mount)
            if parsed is None:
                return None
            normal[key], found = parsed
            matched |= found
    normal_argv = []
    index = 0
    while index < len(argv):
        item = argv[index]
        if item in {"--bind", "-B"} and index + 1 < len(argv):
            parsed = _scratch_binds(argv[index + 1], scratch, mount)
            if parsed is None:
                return None
            bindings, found = parsed
            normal_argv.extend([item, {"binds": bindings}])
            matched |= found
            index += 2
        elif item.startswith("--bind="):
            parsed = _scratch_binds(item.split("=", 1)[1], scratch, mount)
            if parsed is None:
                return None
            bindings, found = parsed
            normal_argv.append({"--bind": bindings})
            matched |= found
            index += 1
        else:
            normal_argv.append(item)
            index += 1
    # Retain strict identities for custom PATH/library/config references into
    # scratch. Only the reviewed storage roles above may be normalised.
    remaining = [v for v in normal.values() if isinstance(v, str)]
    remaining.extend(v for v in normal_argv if isinstance(v, str))
    if not matched or any(scratch in value or mount in value for value in remaining):
        return None
    return normal, normal_argv


def version_context_identity(command, versions_script_sha256, *, environ=None):
    environ = os.environ if environ is None else environ
    context = _environment_context(environ)
    managed = (_managed_scratch_context(command, context, environ)
               if versions_script_sha256 == MANAGED_SCRATCH_INVENTORY_SHA256 else None)
    environment, invocation = managed if managed is not None else (context, command)
    uname = platform.uname()
    release = Path("/etc/os-release")
    host = {"node": uname.node, "system": uname.system, "release": uname.release,
            "machine": uname.machine, "os_release_sha256": digest(release) if release.is_file() else None}
    # Host/runtime/software selection stays explicit; no environment or command
    # values (which may contain credentials) are written to the inventory.
    value = {"schema_version": 1, "managed_scratch": managed is not None,
             "environment": environment, "command": invocation, "host": host}
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def verified_cache_bytes(path, key):
    try:
        if path.is_symlink() or path.with_suffix(".json").is_symlink():
            return None
        with metadata_snapshot():
            metadata = read_bytes(path.with_suffix(".json"))
            if metadata is None:
                return None
            saved = json.loads(metadata)
            raw = read_bytes(path)
            if raw is not None and saved == {"schema_version": 1, "key": key,
                                             "sha256": hashlib.sha256(raw).hexdigest()}:
                return raw
    except (OSError, ValueError, KeyError):
        pass
    return None


def cache_ready(path, key):
    return verified_cache_bytes(path, key) is not None


def copy_verified_cache(path, key, destination):
    raw = verified_cache_bytes(path, key)
    if raw is None:
        return False
    # Write the exact bytes verified above, not a second, potentially replaced
    # source read after checking the digest.
    destination.write_bytes(raw)
    return True


def private_cache(root):
    if root.is_symlink():
        raise ValueError("Symlinked version cache directory")
    root.parent.mkdir(parents=True, exist_ok=True)
    with namespace_lock(root.with_name(root.name + '.initialize.lock'), exclusive=True):
        try:
            root.mkdir(mode=0o700)
            created = True
        except FileExistsError:
            created = False
        # NFS servers can override mkdir's requested mode. Restrict only a new,
        # owned directory through its no-follow descriptor; unsafe existing
        # caches remain an error. Serialize creation before other workers read it.
        descriptor = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        try:
            info = os.fstat(descriptor)
            if info.st_uid != os.getuid():
                raise ValueError("Version cache must be owned by this user with mode 700")
            if created:
                os.fchmod(descriptor, 0o700)
            current = os.fstat(descriptor)
            path_info = root.lstat()
            if ((current.st_dev, current.st_ino) != (path_info.st_dev, path_info.st_ino)
                    or not stat.S_ISDIR(path_info.st_mode) or current.st_mode & 0o077):
                raise ValueError("Version cache must be owned by this user with mode 700")
        finally:
            os.close(descriptor)


SQLITE_VERSION_CACHE_LOCK_V1 = True


def cache_lock_path(cache_file):
    """Preserve legacy publisher coordination for SQLite or flat caches."""
    cache_file = Path(cache_file)
    database = cache_file.parent / DATABASE
    if database.is_symlink():
        raise ValueError("Symlinked version cache database")
    return cache_file.parent / ".versions-cache.lock" if database.is_file() else Path(str(cache_file) + ".lock")


def collection_lock_path(cache_file):
    """Bound metadata inodes while separating slow, unrelated collections."""
    cache_file = Path(cache_file)
    stripe = hashlib.sha256(cache_file.name.encode()).digest()[0] % 64
    return cache_file.parent / f".versions-collection.{stripe:02x}.lock"


def publish_cache(cache_file, source, key, *, lock_publication=False, lock_timeout=300):
    guard = (namespace_lock(cache_lock_path(cache_file), exclusive=True, timeout=lock_timeout)
             if lock_publication else contextlib.nullcontext())
    with guard:
        if cache_file.is_symlink():
            raise ValueError("Symlinked version cache file")
        temporary = cache_file.with_suffix(f".tmp.{os.getpid()}")
        try:
            shutil.copyfile(source, temporary)
            os.replace(temporary, cache_file)
            atomic_json(cache_file.with_suffix(".json"), {"schema_version": 1, "key": key,
                        "sha256": digest(cache_file)})
            if (cache_file.parent / DATABASE).is_file():
                # Publication, including SQLite packing, is the only cache
                # write phase. Other contexts may collect outside this guard.
                limit = cache_policy(cache_file.parent)
                pack(cache_file.parent, scratch="/tmp", writers_stopped=True, delete_originals=True)
                if limit is not None:
                    pack(cache_file.parent, scratch="/tmp", writers_stopped=True, cache_limit=limit)
        finally:
            temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("image", "environment", "context", "cache-ready", "publish-cache", "private-cache", "copy-cache", "cache-lock", "collection-lock"))
    parser.add_argument("--image", type=Path)
    parser.add_argument("--cache-file", type=Path)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--destination", type=Path)
    parser.add_argument("--key")
    parser.add_argument("--command")
    parser.add_argument("--versions-script-sha256")
    parser.add_argument("--if-ready", action="store_true")
    parser.add_argument("--lock-publication", action="store_true")
    parser.add_argument("--lock-timeout", type=int, default=300)
    args = parser.parse_args()
    if args.action == "image":
        if args.image is None:
            parser.error("--image is required")
        with measure("version_image_identity"):
            print(image_identity(args.image))
    elif args.action == "environment":
        print(environment_identity())
    elif args.action == "context":
        if args.command is None or args.versions_script_sha256 is None:
            parser.error("--command and --versions-script-sha256 are required")
        print(version_context_identity(args.command, args.versions_script_sha256))
    elif args.action == "private-cache":
        private_cache(args.cache_file)
    elif args.action == "cache-lock":
        print(cache_lock_path(args.cache_file))
    elif args.action == "collection-lock":
        print(collection_lock_path(args.cache_file))
    elif args.action == "cache-ready":
        raise SystemExit(0 if cache_ready(args.cache_file, args.key) else 1)
    elif args.action == "copy-cache":
        if not copy_verified_cache(args.cache_file, args.key, args.destination):
            if args.if_ready:
                raise SystemExit(3)
            raise ValueError("Version inventory is not verified for this context")
    else:
        publish_cache(args.cache_file, args.source, args.key,
                      lock_publication=args.lock_publication, lock_timeout=args.lock_timeout)


if __name__ == "__main__":
    main()

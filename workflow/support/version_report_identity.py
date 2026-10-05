"""Fresh image identities and conservative version-inventory cache contexts."""

import argparse
import errno
import fcntl
import hashlib
import json
import os
import shutil
import stat
import struct
import sys
from pathlib import Path

from input_generation_array_state import _stat_identity, atomic_json, digest
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


def environment_identity():
    diagnostic = {"PWD", "OLDPWD", "SHLVL", "_", "GG_PERFORMANCE_DIR", "GG_JOB_ID", "GG_ARRAY_TASK_ID",
                  "JOB_ID", "GG_RESOURCE_OWNER_PID"}
    context = {}
    for key, value in os.environ.items():
        suffix = key.removeprefix("SINGULARITYENV_").removeprefix("APPTAINERENV_")
        if suffix in diagnostic or suffix.startswith(("SLURM_", "SGE_")):
            continue
        context[key] = value
    # Store only the digest; environment values may contain credentials.
    return hashlib.sha256(json.dumps(context, sort_keys=True).encode()).hexdigest()


def cache_ready(path, key):
    try:
        if path.is_symlink() or path.with_suffix(".json").is_symlink():
            return False
        saved = json.loads(path.with_suffix(".json").read_text())
        return saved == {"schema_version": 1, "key": key, "sha256": digest(path)}
    except (OSError, ValueError, KeyError):
        return False


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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("image", "environment", "cache-ready", "publish-cache", "private-cache"))
    parser.add_argument("--image", type=Path)
    parser.add_argument("--cache-file", type=Path)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--key")
    args = parser.parse_args()
    if args.action == "image":
        if args.image is None:
            parser.error("--image is required")
        with measure("version_image_identity"):
            print(image_identity(args.image))
    elif args.action == "environment":
        print(environment_identity())
    elif args.action == "private-cache":
        private_cache(args.cache_file)
    elif args.action == "cache-ready":
        raise SystemExit(0 if cache_ready(args.cache_file, args.key) else 1)
    else:
        if args.cache_file.is_symlink():
            raise ValueError("Symlinked version cache file")
        temporary = args.cache_file.with_suffix(f".tmp.{os.getpid()}")
        try:
            shutil.copyfile(args.source, temporary)
            os.replace(temporary, args.cache_file)
            atomic_json(args.cache_file.with_suffix(".json"), {"schema_version": 1, "key": args.key,
                        "sha256": digest(args.cache_file)})
        finally:
            temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    main()

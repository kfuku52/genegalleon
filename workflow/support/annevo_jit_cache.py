"""Private, runtime/source/CPU-bound cache for ANNEVO's native Numba kernels."""
from __future__ import annotations

import hashlib
import json
import os
import platform
import stat
from pathlib import Path

from annevo_model_selection import protect_outputs
from annevo_runtime_download import atomic_json, canonical_bytes
from shared_namespace_lock import namespace_lock


def cpu_identity():
    value = {'machine': platform.machine(), 'system': platform.system(),
             'NUMBA_CPU_NAME': os.environ.get('NUMBA_CPU_NAME', ''),
             'NUMBA_CPU_FEATURES': os.environ.get('NUMBA_CPU_FEATURES', '')}
    path = Path('/proc/cpuinfo')
    if path.is_file():
        first = path.read_text().split('\n\n', 1)[0]
        value['features'] = {key.strip(): item.strip() for key, item in
                             (line.split(':', 1) for line in first.splitlines() if ':' in line)
                             if key.strip() in {'vendor_id', 'cpu family', 'model', 'flags', 'Features'}}
    return value


def private_directory(path):
    path = Path(path).absolute()
    if path != path.resolve():
        raise ValueError('ANNEVO JIT cache path is redirected')
    path.mkdir(mode=0o700, parents=True, exist_ok=True)
    info = path.lstat()
    if not stat.S_ISDIR(info.st_mode) or info.st_uid != os.getuid() or info.st_mode & 0o077:
        raise ValueError('ANNEVO JIT cache directory must be private (0700) and owned by this user: '
                         +str(path)+' (uid='+str(info.st_uid)+', mode='+oct(stat.S_IMODE(info.st_mode))
                         +'); select a filesystem that preserves private POSIX permissions')
    return path


def validate_members(path):
    for member in path.rglob('*'):
        info = member.lstat()
        if (not (stat.S_ISREG(info.st_mode) or stat.S_ISDIR(info.st_mode))
                or info.st_uid != os.getuid() or info.st_mode & 0o022):
            raise ValueError('ANNEVO JIT cache member is redirected, unsafe or owned by another user')


def prepare(root, runtime, protected):
    root = Path(root).absolute()
    generation = Path(runtime['generation']).resolve()
    if (root.resolve().is_relative_to(generation)
            or any(Path(p).resolve().is_relative_to(root.resolve()) for p in protected)):
        raise ValueError('ANNEVO JIT cache must be separate from inputs and the installed runtime')
    root = private_directory(root)
    identity = dict(schema_version=1, kind='genegalleon-annevo-numba-cache',
                    runtime_digest=runtime['runtime_digest'], source_commit=runtime['source_commit'],
                    generation=str(generation), cpu=cpu_identity())
    key = hashlib.sha256(canonical_bytes(identity)).hexdigest()
    directory = root / key
    metadata = directory / 'cache_identity.json'
    protect_outputs([directory, metadata], protected)
    with namespace_lock(root / (key + '.lock'), exclusive=True, timeout=86400):
        if directory.exists():
            directory = private_directory(directory)
            validate_members(directory)
            if not metadata.is_file() or json.loads(metadata.read_text()) != identity:
                raise ValueError('ANNEVO JIT cache identity changed')
        else:
            directory = private_directory(directory)
            atomic_json(metadata, identity)
            metadata.chmod(0o400)
    # Numba's own cache index additionally checks compiler/CPU compatibility.
    # Native atomic cache writes can be shared by this allocation's workers.
    return directory, key

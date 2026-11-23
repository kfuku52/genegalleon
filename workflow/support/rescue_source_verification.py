"""Pin the content generation of rescue source files inside the frozen plan."""
import errno
import fcntl
import os
import platform
import stat
import struct
from functools import lru_cache
from pathlib import Path

if __package__:
    from .input_generation_array_state import FreshDigestBatch, digest
else:
    from input_generation_array_state import FreshDigestBatch, digest


@lru_cache(maxsize=1)
def _lustre_devices():
    # Only use the Lustre UAPI on mounts and architectures with this ioctl ABI.
    if platform.system() != 'Linux' or platform.machine() not in {'x86_64', 'aarch64'}:
        return {}
    try:
        lines = Path('/proc/self/mountinfo').read_text().splitlines()
    except OSError:
        return {}
    devices = {}
    for line in lines:
        mount, filesystem = line.split(' - ', 1)
        kind, source, *_ = filesystem.split()
        if kind == 'lustre':
            major, minor = mount.split()[2].split(':')
            devices[os.makedev(int(major), int(minor))] = source
    return devices


def content_generation(path):
    """Observe Lustre's flushed data version; unsupported files use fresh hashes."""
    devices = _lustre_devices()
    if not devices or os.stat(path).st_dev not in devices:
        return None
    with os.fdopen(os.open(path, os.O_RDONLY | os.O_NONBLOCK), 'rb') as handle:
        info = os.fstat(handle.fileno())
        if not stat.S_ISREG(info.st_mode) or info.st_dev not in devices:
            return None
        # lustre_user.h: LL_IOC_DATA_VERSION, __u64 version, __u32 layout,
        # __u32 flags. LL_DV_RD_FLUSH includes other clients' dirty writes.
        try:
            # LL_IOC_PATH2FID encodes sizeof(long); its output is a 16-byte lu_fid.
            # The FID identifies the same object across client devices.
            fid = fcntl.ioctl(handle, 0x800866AD, bytes(16))
            raw = fcntl.ioctl(handle, 0x801066DA, struct.pack('=QII', 0, 0, 1))
        except OSError as exc:
            if exc.errno in {errno.ENOTTY, errno.EINVAL, errno.EOPNOTSUPP, errno.ENOSYS,
                             errno.EPERM, errno.EACCES, errno.ENODATA, errno.EOVERFLOW}:
                return None
            raise
        version, layout, _ = struct.unpack('=QII', raw)
        if version and any(fid):
            return [str(Path(path).resolve(strict=True)), devices[info.st_dev], fid.hex(),
                    info.st_size, info.st_mtime_ns, info.st_ctime_ns, version, layout]
    return None


def source_paths(request):
    return {source[role] for source in request['sources'].values() for role in ('fasta', 'gff', 'genome')}


def _digest_at_generation(path, generation):
    # Keep the Lustre object open while digest() reads. Replacing and restoring a
    # pathname between the two version checks must never bind another file's bytes.
    with os.fdopen(os.open(path, os.O_RDONLY | os.O_NONBLOCK), 'rb') as handle:
        pinned = '/proc/self/fd/' + str(handle.fileno())
        if content_generation(pinned) != generation:
            raise OSError('Rescue source changed before hashing: ' + path)
        sha = digest(pinned)
        if content_generation(pinned) != generation:
            raise OSError('Rescue source changed during hashing: ' + path)
        return sha


def freeze_sources(paths, expected=None, *, other_paths=()):
    """Hash one plan boundary and bind versioned source bytes on both sides."""
    paths = set(map(str, paths))
    versions = {}
    for path in sorted(paths):
        generation = content_generation(path)
        if generation is not None:
            versions[path] = generation
    # On a derived plan the worker will hash unversioned sources itself.
    plain = set(map(str, other_paths))
    if expected is None:
        plain.update(paths)
    batch = FreshDigestBatch()
    hashes = batch.read(sorted(plain - versions.keys()))
    records, objects = {}, {}
    for path, version in versions.items():
        identity = tuple(version[1:])  # Same Lustre FID/version can have multiple names.
        if identity not in objects:
            objects[identity] = _digest_at_generation(path, version)
        sha = objects[identity]
        if expected is not None and sha != expected[path]:
            raise ValueError('Frozen rescue input changed: ' + path)
        hashes[path] = sha
        records[path] = [sha, version]
    for path, version in versions.items():
        if content_generation(path) != version:
            raise OSError('Rescue source changed during planning: ' + path)
    batch.check()
    return hashes, records


def refresh_source_versions(plan):
    request = plan['request']
    _, records = freeze_sources(source_paths(request), request['files'])
    plan.pop('source_versions', None)
    if records:
        plan['source_versions'] = records


class RescueSourceBatch(FreshDigestBatch):
    """Reuse only plan-bound bytes. Fence subsequent use until publication."""

    def __init__(self, plan):
        super().__init__(scoped_reads=True)
        self.records = plan.get('source_versions', {})
        self.expected = plan['request']['files']
        self.live = {}

    def _hash(self, path):
        record = self.records.get(path)
        if not isinstance(record, list) or len(record) != 2 or record[0] != self.expected.get(path):
            return super()._hash(path)
        if path not in self.live:
            current = content_generation(path)
            if current is None:
                return super()._hash(path)
            sha = record[0] if current == record[1] else _digest_at_generation(path, current)
            self.live[path] = current, sha
        return self.live[path][1]

    def _check_paths(self, paths):
        for path in paths:
            super()._check_paths((path,))
            if path in self.live and content_generation(path) != self.live[path][0]:
                raise OSError('Frozen rescue source changed during use: ' + path)

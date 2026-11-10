"""Validate wrapped gzip FASTQ and bind the exact compressed bytes."""

import gzip
import hashlib
import os
import stat
import zlib


class FastqChanged(OSError):
    """The input changed while deriving its proof."""


class _DigestingReader:
    def __init__(self, handle):
        self.handle = handle
        self.digest = hashlib.sha256()

    def read(self, size=-1):
        data = self.handle.read(size)
        self.digest.update(data)
        return data


def _complete_fastq(handle):
    has_records = False
    readline = handle.readline
    while True:
        header = readline()
        if not header:
            return has_records
        if not header.startswith(b"@") or len(header.rstrip(b"\r\n")) < 2:
            return False
        sequence_length = 0
        while True:
            line = readline()
            if not line:
                return False
            if line.startswith(b"+"):
                break
            sequence_length += len(line.rstrip(b"\r\n"))
        if not sequence_length:
            return False
        quality_length = 0
        while quality_length < sequence_length:
            line = readline()
            if not line:
                return False
            quality_length += len(line.rstrip(b"\r\n"))
        if quality_length != sequence_length:
            return False
        has_records = True


def _identity(info):
    return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns


def _parse_gzip(file):
    try:
        with gzip.GzipFile(fileobj=file, mode="rb") as handle:
            return _complete_fastq(handle)
    except (EOFError, OSError, ValueError, zlib.error):
        return False


def proof_fastq(path):
    """Return structure validity, compressed size and SHA-256 for a stable file."""
    with os.fdopen(os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK), "rb") as raw:
        before = os.fstat(raw.fileno())
        if not stat.S_ISREG(before.st_mode) or not before.st_size:
            return False, before.st_size, None
        reader = _DigestingReader(raw)
        if not _parse_gzip(reader):
            return False, before.st_size, None
        while reader.read(1024 * 1024):
            pass
        try:
            current = path.stat(follow_symlinks=False)
        except FileNotFoundError as error:
            raise FastqChanged(path) from error
        if _identity(before) != _identity(os.fstat(raw.fileno())) or _identity(before) != _identity(current):
            raise FastqChanged(path)
        return True, before.st_size, reader.digest.hexdigest()

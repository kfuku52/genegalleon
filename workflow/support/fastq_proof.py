"""Validate getfastq reads and quantification statistics while hashing them."""

import csv
import gzip
import hashlib
import io
import math
import multiprocessing
import os
import stat
import zlib
from concurrent.futures import ProcessPoolExecutor


class FastqChanged(OSError):
    """The input changed while deriving its proof."""


class _DigestingReader(io.BufferedIOBase):
    def __init__(self, handle):
        self.handle = handle
        self.digest = hashlib.sha256()

    def read(self, size=-1):
        data = self.handle.read(size)
        self.digest.update(data)
        return data

    def readable(self):
        return True

    read1 = read


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
        with gzip.GzipFile(fileobj=file, mode="rb") as stream, io.BufferedReader(stream) as handle:
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


def proof_fastqs(paths, workers):
    """Prove independent files in manifest order using the task's CPU allocation."""
    workers = min(len(paths), workers)
    if workers <= 1:
        return [proof_fastq(path) for path in paths]
    # The workflow calls this from a fresh, single-threaded Python heredoc.
    # Fork keeps the worker entrypoint importable when __main__ is stdin.
    with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("fork")) as pool:
        return list(pool.map(proof_fastq, paths))


def _positive_number(value):
    try:
        number = float(value)
        return math.isfinite(number) and number > 0
    except (TypeError, ValueError):
        return False


def quant_stats_digest(path, run):
    """Return the size and SHA-256 of the same bytes used to validate the run."""
    with os.fdopen(os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK), "rb") as raw:
        before = os.fstat(raw.fileno())
        if not stat.S_ISREG(before.st_mode):
            raise ValueError(f"Expected regular getfastq statistics for run: {run}")
        reader = _DigestingReader(raw)
        with io.TextIOWrapper(reader, encoding="utf-8", newline="") as handle:
            rows = (row for row in csv.DictReader(handle, delimiter="\t")
                    if str(row.get("run") or "").strip() == run)
            row = next(rows, None)
            if row is None or next(rows, None) is not None:
                raise ValueError(f"Missing or ambiguous getfastq statistics for run: {run}")
        if not _positive_number(row.get("num_written")) or not any(
            _positive_number(row.get(column)) for column in ("bp_fastp_in", "bp_rrna_in", "bp_contam_in")
        ):
            raise ValueError(f"Found unusable getfastq statistics for run: {run}")
        if _identity(before) != _identity(os.fstat(raw.fileno())) or _identity(before) != _identity(path.lstat()):
            raise ValueError(f"getfastq statistics changed while validating: {path}")
        return before.st_size, reader.digest.hexdigest()

"""Explicit exact-byte transport identities; structural validation stays default."""
import hashlib
import os
import re
from pathlib import Path

from .local import validate_gzip_with_cache


def byte_integrity_contract(expected_sha256=None, expected_size=None):
    if expected_sha256 is None and expected_size is None:
        return None
    if (not isinstance(expected_sha256, str) or re.fullmatch(r"[0-9a-f]{64}", expected_sha256) is None
            or type(expected_size) is not int or expected_size <= 0):
        raise ValueError("Exact-byte downloads require both a full SHA-256 and a positive integer size")
    return {"sha256": expected_sha256, "size": expected_size}


def validate_download(path, *, expected_sha256=None, expected_size=None, **structural_options):
    contract = byte_integrity_contract(expected_sha256, expected_size)
    if contract is None:
        return validate_gzip_with_cache(path, **structural_options)
    path = Path(path)
    try:
        with path.open("rb") as source:
            before = os.fstat(source.fileno())
            if before.st_size != contract["size"]:
                return ValueError("downloaded bytes differ from the expected size")
            digest = hashlib.sha256()
            for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
                digest.update(block)
            after = os.fstat(source.fileno())
        current = path.stat()
        def signature(value):
            return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns
        if signature(before) != signature(after) or signature(after) != signature(current):
            return ValueError("downloaded bytes changed during integrity verification")
        if digest.hexdigest() != contract["sha256"]:
            return ValueError("downloaded bytes differ from the expected SHA-256")
    except OSError as error:
        return error
    return None

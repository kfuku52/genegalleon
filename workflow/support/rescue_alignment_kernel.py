"""Load the matching container's SIMD implementation of protein-only Gotoh QC."""
import ctypes
import hashlib
from functools import lru_cache
from pathlib import Path

from cdskit.protein_similarity import PairScore

_LIBRARY = Path('/usr/local/lib/libgg_rescue_protein_alignment.so')
_SOURCE = Path(__file__).with_name('rescue_protein_alignment.cpp')
MIN_FAST_LENGTH = 24
MAX_FAST_LENGTH = 16_384


@lru_cache(maxsize=1)
def compiled_kernel():
    """A checkout can also run against an image predating the native library."""
    if not _LIBRARY.is_file():
        return None
    lib = ctypes.CDLL(str(_LIBRARY))
    library_source = lib.gg_rescue_protein_alignment_source
    library_source.argtypes = ()
    library_source.restype = ctypes.c_char_p
    with _SOURCE.open('rb') as stream:
        source = hashlib.file_digest(stream, 'sha256').hexdigest().encode('ascii')
    if library_source() != source:
        return None
    kernel = lib.gg_rescue_protein_alignment_v1
    kernel.argtypes = (ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.POINTER(ctypes.c_int))
    kernel.restype = ctypes.c_int
    return kernel


def score_counts(n, m, matches, aligned, loss_a, loss_b):
    """Keep operation order and the empty-junction convention of CDSKit."""
    identity = matches / max(1, aligned)
    coverage_a, coverage_b = aligned / n, aligned / m
    score = min(coverage_a, coverage_b) * (0.85 * identity + 0.15)
    return PairScore(score, identity, coverage_a, coverage_b, 1.0, aligned, bool(loss_a), bool(loss_b))


def buffer():
    return (ctypes.c_int * 4)()

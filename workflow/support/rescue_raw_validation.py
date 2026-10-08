"""Bounded, context-local reuse of candidate-independent genomic validation.

The owner must supply one immutable indexed genome and a fixed validator. Only
DNA-derived overrides are retained. Candidate evidence, query metrics, terminal
completion and placement decisions always belong to the current candidate.
"""

from __future__ import annotations

import copy
import math
import struct
import sys
from collections import OrderedDict

DNA_FIELDS = ("cds", "sequence", "assembly_ambiguous_bases", "problems")


def _metric(value):
    if type(value) is int:
        return int, value
    if type(value) is float and math.isfinite(value):
        return float, struct.pack("!d", value)
    return None


def raw_validation_key(model):
    """Abstain on malformed input so the original validator raises its error."""
    if type(model) is not dict:
        return None
    seqid, strand, cds = model.get("seqid"), model.get("strand"), model.get("cds")
    if (type(seqid) is not str or type(strand) is not str or strand not in ("+", "-")
            or type(cds) not in (list, tuple) or not cds
            or type(model.get("frameshift")) is not bool):
        return None
    blocks = []
    for block in cds:
        if (type(block) not in (list, tuple) or len(block) != 3
                or any(type(value) is not int for value in block)):
            return None
        blocks.append((type(block), *block))
    coverage, identity = _metric(model.get("coverage")), _metric(model.get("identity"))
    if coverage is None or identity is None:
        return None
    return seqid, strand, type(cds), tuple(blocks), model["frameshift"], coverage, identity


def _retained_size(value):
    """Conservative accounting, including repeated references and LRU overhead."""
    size = sys.getsizeof(value)
    if isinstance(value, dict):
        return size + sum(_retained_size(k) + _retained_size(v) for k, v in value.items())
    if isinstance(value, (tuple, list)):
        return size + sum(_retained_size(v) for v in value)
    return size


class RawValidationMemo:
    """One bounded LRU per immutable indexed-genome/validator/code context."""

    def __init__(self, validator, genome, code, params, *, max_entries=8192,
                 max_bytes=32 * 1024 * 1024):
        if (type(max_entries) is not int or max_entries < 0
                or type(max_bytes) is not int or max_bytes < 0):
            raise ValueError("Raw validation memo limits must be nonnegative integers")
        self._validator, self._genome, self._code = validator, genome, code
        self._params = params
        self._scope_params = copy.deepcopy(params)
        self._cache = OrderedDict()
        self.max_entries, self.max_bytes = max_entries, max_bytes
        self._bytes = self._peak_bytes = self._peak_entries = 0
        self._counts = dict(calls=0, hits=0, misses=0, bypasses=0, evictions=0,
                            oversized=0, exceptions=0, scope_invalidations=0)

    def validate(self, model):
        self._counts["calls"] += 1
        if self._params != self._scope_params:
            self._cache.clear()
            self._bytes = 0
            self._scope_params = copy.deepcopy(self._params)
            self._counts["scope_invalidations"] += 1
        key = raw_validation_key(model) if self.max_entries and self.max_bytes else None
        if key is not None and key in self._cache:
            overlay, _ = self._cache[key]
            self._cache.move_to_end(key)
            self._counts["hits"] += 1
            return {**model, **copy.deepcopy(overlay)}
        self._counts["bypasses" if key is None else "misses"] += 1
        try:
            result = self._validator(model, self._genome, self._code, self._params)
        except Exception:  # noqa: BLE001 -- count and preserve the original exception exactly
            self._counts["exceptions"] += 1
            raise
        if key is None:
            return result
        # Future validator changes must not silently discard new derived fields.
        if (type(result) is not dict or set(result) != set(model) | set(DNA_FIELDS)
                or any(result[k] is not model[k] for k in model if k not in DNA_FIELDS)):
            self._counts["bypasses"] += 1
            return result
        overlay = {field: copy.deepcopy(result[field]) for field in DNA_FIELDS}
        weight = _retained_size(key) + _retained_size(overlay) + 256
        if weight > self.max_bytes:
            self._counts["oversized"] += 1
            return result
        while self._cache and (len(self._cache) >= self.max_entries
                               or self._bytes + weight > self.max_bytes):
            _, (_, removed) = self._cache.popitem(last=False)
            self._bytes -= removed
            self._counts["evictions"] += 1
        self._cache[key] = overlay, weight
        self._bytes += weight
        self._peak_bytes = max(self._peak_bytes, self._bytes)
        self._peak_entries = max(self._peak_entries, len(self._cache))
        return result

    def diagnostics(self):
        return {"schema": 1, **self._counts, "max_entries": self.max_entries,
                "max_bytes": self.max_bytes, "entries": len(self._cache),
                "accounted_bytes": self._bytes, "peak_accounted_bytes": self._peak_bytes,
                "peak_entries": self._peak_entries,
                "scope": "One immutable indexed genome, genetic code and validator; exact raw inputs",
                "shared_fields": list(DNA_FIELDS),
                "uncached": ["terminal_completion", "query", "evidence", "ownership", "placement"],
                "memory_accounting": "Conservative Python object estimate; excludes genome and caller records"}

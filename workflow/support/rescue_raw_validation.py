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

DNA_FIELDS = ("cds", "sequence", "genetic_code", "assembly_ambiguous_bases", "problems")
_ATOMIC_TYPES = frozenset((str, int, float, bool, type(None)))


def copy_validation_data(value):
    """Deep-copy validation overlays without dispatching every immutable leaf.

    Exact dicts/lists are common in CDS and completion metadata. Keep their
    mutable contents independent, preserve aliases/cycles, and delegate other
    types (including subclasses) to the ordinary deepcopy protocol. This is a
    copy operation, not a replacement for input or finite-number validation.
    """
    return _copy_validation_data(value, {})


def _copy_validation_data(value, memo):
    kind = type(value)
    if kind in _ATOMIC_TYPES:
        return value
    token = id(value)
    if token in memo:
        return memo[token]
    if kind is dict:
        result = {}
        memo[token] = result
        for key, item in value.items():
            result[key if type(key) in _ATOMIC_TYPES else _copy_validation_data(key, memo)] = (
                item if type(item) in _ATOMIC_TYPES else _copy_validation_data(item, memo))
    elif kind is list:
        result = []
        memo[token] = result
        result.extend(item if type(item) in _ATOMIC_TYPES else _copy_validation_data(item, memo) for item in value)
    else:
        return copy.deepcopy(value, memo)
    # Match deepcopy's lifetime guarantee when a custom copy protocol creates
    # temporary inputs or mutates the source graph during this operation.
    memo.setdefault(id(memo), []).append(value)
    return result


def _metric(value):
    if type(value) is int:
        return int, value
    if type(value) is float and math.isfinite(value):
        return float, struct.pack("!d", value)
    return None


def raw_validation_key(model, *, genomic_only=False):
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
    if genomic_only:
        if any(not 0 <= model[field] <= 1 for field in ("coverage", "identity")):
            return None
        return seqid, strand, type(cds), tuple(blocks), model["frameshift"]
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
                 max_bytes=32 * 1024 * 1024, candidate_homology=None):
        if (type(max_entries) is not int or max_entries < 0
                or type(max_bytes) is not int or max_bytes < 0):
            raise ValueError("Raw validation memo limits must be nonnegative integers")
        self._validator, self._genome, self._code = validator, genome, code
        # Opt in only for validators whose DNA overrides depend on metrics
        # solely through this homology policy. Arbitrary validators keep the
        # original exact-metric key and result contract.
        self._candidate_homology = candidate_homology
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
        key = (raw_validation_key(model, genomic_only=self._candidate_homology is not None)
               if self.max_entries and self.max_bytes else None)
        if key is not None and key in self._cache:
            overlay, _ = self._cache[key]
            self._cache.move_to_end(key)
            self._counts["hits"] += 1
            result = {**model, **copy_validation_data(overlay)}
            if self._candidate_homology is not None:
                result["problems"] = sorted(set(result["problems"]) | set(
                    self._candidate_homology(model["coverage"], model["identity"], self._params)))
            return result
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
        overlay = {field: copy_validation_data(result[field]) for field in DNA_FIELDS}
        if self._candidate_homology is not None:
            candidate_flags = set(self._candidate_homology(model["coverage"], model["identity"], self._params))
            overlay["problems"] = [p for p in overlay["problems"] if p not in candidate_flags]
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
                "scope": ("One immutable indexed genome, genetic code and validator; exact genomic inputs"
                          if self._candidate_homology is not None else
                          "One immutable indexed genome, genetic code and validator; exact raw inputs"),
                "homology_evaluated_per_candidate": self._candidate_homology is not None,
                "shared_fields": list(DNA_FIELDS),
                "uncached": ["terminal_completion", "query", "evidence", "ownership", "placement"],
                "memory_accounting": "Conservative Python object estimate; excludes genome and caller records"}

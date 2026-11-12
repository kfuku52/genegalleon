"""Reciprocal, unique protein matches from an existing tabular alignment."""
import math
from pathlib import Path


def reciprocal_matches(path, lengths_a, lengths_b, *, min_identity, coverage_a, coverage_b):
    """Retain two distinct partners per protein and abstain on near ties."""
    best_a, best_b = {}, {}

    def add(index, key, partner, score, metrics):
        hits = index.get(key, ())
        if any(other == partner and value >= score for value, other, _ in hits):
            return
        ranked = [(value, other, old) for value, other, old in hits if other != partner]
        ranked.append((score, partner, metrics))
        index[key] = sorted(ranked, key=lambda hit: (-hit[0], hit[1]))[:2]

    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 12:
                raise ValueError('Incomplete tabular protein alignment: ' + str(path))
            a, b = fields[:2]
            if a not in lengths_a or b not in lengths_b:
                raise ValueError('Alignment references an unknown protein: ' + str(path))
            identity, score = float(fields[2]) / 100, float(fields[11])
            a0, a1 = sorted(map(int, fields[6:8]))
            b0, b1 = sorted(map(int, fields[8:10]))
            if (not 0 <= identity <= 1 or not math.isfinite(score) or score <= 0
                    or not 0 < a0 <= a1 <= lengths_a[a] or not 0 < b0 <= b1 <= lengths_b[b]):
                raise ValueError('Invalid tabular protein alignment: ' + str(path))
            cov_a, cov_b = (a1 - a0 + 1) / lengths_a[a], (b1 - b0 + 1) / lengths_b[b]
            if identity < min_identity or cov_a < coverage_a or cov_b < coverage_b:
                continue
            metrics = {'identity': identity, 'coverage_a': cov_a, 'coverage_b': cov_b, 'bitscore': score}
            add(best_a, a, b, score, metrics)
            add(best_b, b, a, score, metrics)

    def unique(hits):
        return len(hits) == 1 or hits[1][0] < .95 * hits[0][0]

    for a, hits in sorted(best_a.items()):
        _, b, metrics = hits[0]
        if unique(hits) and unique(best_b[b]) and best_b[b][0][1] == a:
            yield a, b, metrics

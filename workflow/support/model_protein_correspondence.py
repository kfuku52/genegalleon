"""Reciprocal, unique protein matches from an existing tabular alignment."""
import math
from pathlib import Path


def _outranks(hit, other):
    return hit[0] > other[0] or (hit[0] == other[0] and hit[1] < other[1])


def _offer(index, key, hit):
    ranked = index.get(key)
    if ranked is None:
        index[key] = [hit]
        return
    first = ranked[0]
    if hit[1] == first[1]:
        if hit[0] > first[0]:
            ranked[0] = hit
        return
    if _outranks(hit, first):
        ranked[:] = [hit, first]
    elif len(ranked) == 1:
        ranked.append(hit)
    elif _outranks(hit, ranked[1]):
        ranked[1] = hit


def reciprocal_matches(path, lengths_a, lengths_b, *, min_identity, coverage_a, coverage_b):
    """Rank all reported hits, then require unique reciprocal winners to pass QC."""
    best_a, best_b = {}, {}
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 12:
                raise ValueError('Incomplete tabular protein alignment: ' + str(path))
            a, b = fields[0], fields[1]
            la, lb = lengths_a.get(a), lengths_b.get(b)
            if la is None or lb is None:
                raise ValueError('Alignment references an unknown protein: ' + str(path))
            identity, score = float(fields[2]) / 100, float(fields[11])
            a0, a1, b0, b1 = map(int, fields[6:10])
            if (not 0 <= identity <= 1 or not math.isfinite(score) or score <= 0
                    or not (0 < a0 <= la and 0 < a1 <= la and 0 < b0 <= lb and 0 < b1 <= lb)):
                raise ValueError('Invalid tabular protein alignment: ' + str(path))
            metrics = identity, (abs(a1 - a0) + 1) / la, (abs(b1 - b0) + 1) / lb
            _offer(best_a, a, (score, b, metrics))
            _offer(best_b, b, (score, a, metrics))

    def unique(hits):
        return len(hits) == 1 or hits[1][0] < .95 * hits[0][0]

    for a, hits in sorted(best_a.items()):
        score, b, (identity, cov_a, cov_b) = hits[0]
        if (identity >= min_identity and cov_a >= coverage_a and cov_b >= coverage_b
                and unique(hits) and unique(best_b[b]) and best_b[b][0][1] == a):
            yield a, b, {'identity': identity, 'coverage_a': cov_a, 'coverage_b': cov_b, 'bitscore': score}

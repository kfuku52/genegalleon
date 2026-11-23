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


def reciprocal_matches(path, lengths_a, lengths_b, *, min_identity, coverage_a, coverage_b, terminal_hits=None, target_is_a=True, anchors=None):
    """Rank all reported hits, then require unique reciprocal winners to pass QC."""
    best_a, best_b = {}, {}
    terminal = {} if terminal_hits is not None else None
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
            if terminal is not None:
                target, donor, t0, t1, lt, d0, d1, ld = ((a, b, a0, a1, la, b0, b1, lb) if target_is_a else
                                                          (b, a, b0, b1, lb, a0, a1, la))
                minimum = max(15, .10 * lt)
                # A syntenic hit can nominate a bounded projection at lower identity.
                # The projected coding sequence must still pass the configured QC.
                anchored = anchors is not None and donor in anchors.get(target, ())
                identity_floor = min(min_identity, max(.45, .90 * min_identity)) if anchored else min_identity
                if (t1 >= t0 and d1 >= d0 and ld - lt >= minimum and identity >= identity_floor
                        and (t1 - t0 + 1) >= .90 * lt and (d1 - d0 + 1) >= .65 * ld):
                    left, right = d0 - t0, (ld - d1) - (lt - t1)
                    side = None
                    if left >= minimum and abs(right) <= 10:
                        side = 'N'
                    elif right >= minimum and abs(left) <= 10:
                        side = 'C'
                    elif (d1 - d0) - (t1 - t0) >= max(15, .08 * lt) and max(abs(left), abs(right)) < minimum:
                        side = 'internal'
                    if side:
                        prev = terminal.get(target)
                        if prev is None or score > prev[0]:
                            terminal[target] = score, donor, side
            _offer(best_a, a, (score, b, metrics))
            _offer(best_b, b, (score, a, metrics))

    def unique(hits):
        return len(hits) == 1 or hits[1][0] < .95 * hits[0][0]

    if terminal is not None:
        rankings = best_a if target_is_a else best_b
        for target, (score, donor, side) in terminal.items():
            rank = rankings[target]
            if unique(rank) and (score, donor) == rank[0][:2]:
                terminal_hits.append((target, donor, side))
    for a, hits in sorted(best_a.items()):
        score, b, (identity, cov_a, cov_b) = hits[0]
        if (identity >= min_identity and cov_a >= coverage_a and cov_b >= coverage_b
                and unique(hits) and unique(best_b[b]) and best_b[b][0][1] == a):
            yield a, b, {'identity': identity, 'coverage_a': cov_a, 'coverage_b': cov_b, 'bitscore': score}

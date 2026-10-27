"""Deterministic, bounded order/orientation search for two or more JCVI tracks.

This changes display coordinates only. The legacy pairwise length solver is
deliberately separate: its objective, exactness claims and defaults are public.
"""

import itertools
import math
import random
from collections import defaultdict
from pathlib import Path

import numpy as np

EXACT_STATES = 20_000
MODES = {"none", "order", "order_orientation"}
TRACK_FIELDS = {"species", "bed", "seqids", "groups", "fixed_order", "orientation", "fixed_orientation", "colors", "labels"}
OPTIONS = {"mode", "seed", "restarts", "max_rounds"}


def validate_options(config):
    if config.get("mode", "order_orientation") not in MODES:
        raise ValueError("Layout mode must be none, order or order_orientation")
    for key, default, low, high in (("seed", 1, 0, 2**32 - 1), ("restarts", 4, 1, 16), ("max_rounds", 8, 1, 100)):
        value = config.get(key, default)
        if type(value) is not int or not low <= value <= high:
            raise ValueError(f"Layout {key} must be an integer in [{low}, {high}]")


def _identifiers(values, name):
    if (not isinstance(values, list) or not values or any(not isinstance(s, str) or not s.strip()
            or any(c in s for c in ",\n\r\0") or s.endswith("-") for s in values) or len(set(values)) != len(values)):
        raise ValueError(f"{name} must be a nonempty list of unique JCVI identifiers")
    return list(values)


def prepare_tracks(tracks, genomes):
    """Validate once and match JCVI's start/literal-accession rank ordering."""
    if (len(tracks) < 2 or len(tracks) != len(genomes)
            or any(not isinstance(t, dict) for t in tracks)):
        raise ValueError("Layout requires at least two matching tracks and genomes")
    species = _identifiers([t.get("species") for t in tracks], "Track species")
    selected, ranks, counts, units, orientations, locks = [], [], [], [], [], []
    for track, genes in zip(tracks, genomes, strict=True):
        if set(track) - TRACK_FIELDS:
            raise ValueError("Unknown track options: " + ", ".join(sorted(set(track) - TRACK_FIELDS)))
        if type(track.get("fixed_order", False)) is not bool:
            raise ValueError("fixed_order must be a boolean")
        grouped, lookup = defaultdict(list), {}
        for gene in genes:
            if gene.gene_id in lookup:
                raise ValueError("Layout BED gene IDs must be unique within each track")
            lookup[gene.gene_id] = gene.seqid
            grouped[gene.seqid].append(gene)
        order = _identifiers(track.get("seqids", list(grouped)), "Track seqids")
        if set(order) - set(grouped):
            raise ValueError("Selected layout chromosomes must have annotated genes")
        rank = {}
        for sid, members in grouped.items():
            for i, gene in enumerate(sorted(members, key=lambda g: (g.start, g.gene_id))):
                rank[gene.gene_id] = (sid, i)
        groups = track.get("groups", [])
        if not isinstance(groups, list):
            raise ValueError("groups must be a list of chromosome lists")
        membership = {}
        for number, group in enumerate(groups):
            group = _identifiers(group, "Chromosome group")
            if set(group) - set(order) or set(group) & set(membership):
                raise ValueError("Chromosome groups must be disjoint subsets of selected seqids")
            positions = sorted(order.index(s) for s in group)
            if positions != list(range(positions[0], positions[0] + len(group))):
                raise ValueError("Chromosome groups must be contiguous in the input order")
            membership.update((sid, number) for sid in group)
        blocks = []
        for sid in order:
            if blocks and sid in membership and membership.get(blocks[-1][-1]) == membership[sid]:
                blocks[-1].append(sid)
            else:
                blocks.append([sid])
        directions = dict.fromkeys(order, False)
        fixed = {}
        for key in ("orientation", "fixed_orientation"):
            values = track.get(key, {})
            if not isinstance(values, dict) or set(values) - set(order) or any(v not in ("+", "-") for v in values.values()):
                raise ValueError(f"{key} must map selected seqids to '+' or '-'")
            for sid, value in values.items():
                if key == "fixed_orientation":
                    if sid in track.get("orientation", {}) and directions[sid] != (value == "-"):
                        raise ValueError("Initial orientation conflicts with fixed_orientation")
                    fixed[sid] = value == "-"
                directions[sid] = value == "-"
        selected.append(order)
        ranks.append(rank)
        counts.append({s: len(grouped[s]) for s in order})
        units.append(blocks)
        orientations.append(directions)
        locks.append(fixed)
    return species, selected, ranks, counts, units, orientations, locks


def read_blocks(edges, ranks, selected):
    result, excluded = [], []
    seen = set()
    for edge in edges:
        if set(edge) != {"tracks", "blocks"}:
            raise ValueError("Each layout edge requires tracks and blocks")
        pair = edge["tracks"]
        if (not isinstance(pair, list) or len(pair) != 2 or any(type(i) is not int or not 0 <= i < len(ranks) for i in pair)
                or pair[0] == pair[1] or tuple(sorted(pair)) in seen):
            raise ValueError("Layout edges must connect distinct tracks, once per track pair")
        seen.add(tuple(sorted(pair)))
        a, b = pair
        rows, omitted = [], 0
        for line in Path(edge["blocks"]).read_text(encoding="utf-8").splitlines():
            if not line.strip() or (line.startswith("#") and "*" not in line.split()[0]) or line.startswith("StartGeneA"):
                continue
            fields = line.split()
            fields[0] = fields[0].split("*", 1)[-1]
            if len(fields) != 6 or fields[5] not in ("+", "-"):
                raise ValueError("Invalid JCVI simple block")
            try:
                if int(fields[4]) <= 0:
                    raise ValueError("Nonpositive block score")
                p, q = [ranks[a][g] for g in fields[:2]], [ranks[b][g] for g in fields[2:4]]
            except (KeyError, ValueError) as error:
                raise ValueError("Invalid layout block score or gene identifier") from error
            if p[0][0] != p[1][0] or q[0][0] != q[1][0]:
                raise ValueError("A layout ribbon cannot span multiple chromosomes")
            if p[0][0] not in selected[a] or q[0][0] not in selected[b]:
                omitted += 1
                continue
            if fields[5] == "-":
                q.reverse()  # The same interpretation as native JCVI SimpleFile.
            rows.append((p[0][0], p[0][1], p[1][1], q[0][0], q[0][1], q[1][1]))
        result.append((a, b, rows))
        excluded.append(omitted)
    if not result or not any(rows for _, _, rows in result):
        raise ValueError("Layout has no displayed syntenic blocks")
    connected = {i for a, b, rows in result if rows for i in (a, b)}
    if connected != set(range(len(ranks))):
        raise ValueError("Every layout track must have displayed syntenic connections")
    return result, excluded


def crossing_weight(x, y, weights):
    """Weighted midpoint inversions, excluding ties, in O(B log B) space/time.

    Fenwick storage is linear in block count; no all-block-pairs matrix is
    allocated for fragmented genomes. These are not Bézier area intersections.
    """
    order = np.lexsort((y, x))
    yy = np.searchsorted(np.unique(y), y) + 1
    tree = np.zeros(int(max(yy, default=0)) + 1)
    total, crossings, start = 0., 0., 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and x[order[stop]] == x[order[start]]:
            stop += 1
        for index in order[start:stop]:
            rank, prefix = int(yy[index]), 0.
            while rank:
                prefix += tree[rank]
                rank -= rank & -rank
            crossings += weights[index] * max(0., total - prefix)
        for index in order[start:stop]:
            rank = int(yy[index])
            while rank < len(tree):
                tree[rank] += weights[index]
                rank += rank & -rank
            total += weights[index]
        start = stop
    return float(crossings)


def optimize_layout(tracks, genomes, edges, *, mode="order_orientation", scale_mode="shared", seed=1,
                    restarts=4, max_rounds=8, xstart=.12, xend=.92, figsize=(20., 8.), ys=None, gaps=None):
    """Return one shared state per track, keeping track order and all blocks."""
    validate_options(dict(mode=mode, seed=seed, restarts=restarts, max_rounds=max_rounds))
    if scale_mode not in {"shared", "independent"} or not 0 <= xstart < xend <= 1:
        raise ValueError("Invalid layout scale or horizontal extent")
    species, initial, ranks, counts, units, directions, locks = prepare_tracks(tracks, genomes)
    blocks, excluded = read_blocks(edges, ranks, initial)
    ys = list(np.linspace(.8, .2, len(tracks))) if ys is None else list(ys)
    if (len(ys) != len(tracks) or any(not math.isfinite(y) or not 0 < y < 1 for y in ys)
            or any(a <= b for a, b in zip(ys, ys[1:], strict=False))
            or len(figsize) != 2 or any(not math.isfinite(v) or v <= 0 for v in figsize)):
        raise ValueError("Layout tracks need descending y positions and a positive figure size")
    gaps = ([min(.01, .01 * 16 / len(order) + .001) if len(order) > 16 else .01 for order in initial]
            if gaps is None else list(gaps))
    if len(gaps) != len(tracks) or any(not math.isfinite(g) or g < 0 for g in gaps):
        raise ValueError("Layout gaps must be finite nonnegative values, one per track")
    ratios = [(xend - xstart - gap * (len(order) - 1)) / sum(c.values())
              for order, c, gap in zip(initial, counts, gaps, strict=True)]
    if min(ratios) <= 0:
        raise ValueError("Layout chromosome gaps leave no space for genes")
    if scale_mode == "shared":
        ratios = [min(ratios)] * len(tracks)
    widths = [{s: n * ratio for s, n in c.items()} for c, ratio in zip(counts, ratios, strict=True)]
    cached = []
    for a, b, rows in blocks:
        sa, pa, qa, sb, pb, qb = zip(*rows, strict=True) if rows else ([], [], [], [], [], [])
        pa, qa, pb, qb = [np.asarray(x, dtype=float) for x in (pa, qa, pb, qb)]
        weights = (abs(pa - qa) * ratios[a] + abs(pb - qb) * ratios[b]) / 2
        cached.append((a, b, list(sa), pa, qa, list(sb), pb, qb, weights))
    total_weight = sum(float(w.sum()) for *_, w in cached)
    pair_weight = sum(float((w.sum()**2 - (w * w).sum()) / 2) for *_, w in cached)

    def track_offsets(orders):
        offsets = []
        for order, ww, gap in zip(orders, widths, gaps, strict=True):
            start, row = xstart, {}
            for sid in order:
                row[sid] = start
                start += ww[sid] + gap
            offsets.append(row)
        return offsets

    def measure(orders, orientation):
        offsets = track_offsets(orders)
        length, crossing, twist = 0., 0., 0.
        for a, b, sa, pa, qa, sb, pb, qb, w in cached:
            points = []
            for side, seqs, first, last in ((a, sa, pa, qa), (b, sb, pb, qb)):
                reverse = np.array([orientation[side][s] for s in seqs], dtype=bool)
                maximum = np.array([counts[side][s] - 1 for s in seqs])
                start = np.array([offsets[side][s] for s in seqs])
                points.append((start + ratios[side] * np.where(reverse, maximum - first, first),
                               start + ratios[side] * np.where(reverse, maximum - last, last)))
            (p, q), (r, s) = points
            x, y = (p + q) / 2, (r + s) / 2
            length += float(np.sum(w * np.hypot((x - y) * figsize[0], (ys[a] - ys[b]) * figsize[1])))
            twist += float(w[(q - p) * (s - r) < 0].sum())
            crossing += crossing_weight(x, y, w)
        parts = {"length": length / (total_weight * math.hypot(*figsize)) if total_weight else 0.,
                 "crossing": crossing / pair_weight if pair_weight else 0.,
                 "twist": twist / total_weight if total_weight else 0.}
        return {"total": sum(parts.values()), **parts}

    def arranged_units(order, side):
        membership = {s: j for j, group in enumerate(units[side]) for s in group}
        result = []
        for sid in order:
            if result and membership[result[-1][-1]] == membership[sid]:
                result[-1].append(sid)
            else:
                result.append([sid])
        return result

    def barycentric_order(state, side):
        """A global order proposal, always judged by the full objective.

        Adjacent swaps alone would need many rounds for fragmented assemblies.
        Both neighbors contribute to the same proposal for a middle track.
        """
        starts = track_offsets(state[0])
        sums = {s: [0., 0.] for s in initial[side]}
        for a, b, sa, pa, qa, sb, pb, qb, weights in cached:
            if side not in (a, b):
                continue
            moving, local, other, remote, fixed = ((sa, (pa + qa) / 2, sb, (pb + qb) / 2, b)
                if side == a else (sb, (pb + qb) / 2, sa, (pa + qa) / 2, a))
            for sid, rank, partner, partner_rank, weight in zip(moving, local, other, remote, weights, strict=True):
                if state[1][side][sid]:
                    rank = counts[side][sid] - 1 - rank
                if state[1][fixed][partner]:
                    partner_rank = counts[fixed][partner] - 1 - partner_rank
                desired = starts[fixed][partner] + partner_rank * ratios[fixed] - rank * ratios[side]
                sums[sid][0] += desired * weight
                sums[sid][1] += weight
        def center(sid):
            value, weight = sums[sid]
            return value / weight if weight else starts[side][sid]
        groups = [sorted(group, key=center) for group in arranged_units(state[0][side], side)]
        def group_center(group):
            weight = sum(sums[s][1] for s in group)
            return sum(sums[s][0] for s in group) / weight if weight else starts[side][group[0]]
        return [s for group in sorted(groups, key=group_center) for s in group]

    ordinal = [{s: i for i, s in enumerate(row)} for row in initial]

    def key(state):
        orders, orientation = state
        return (sum(orientation[i][s] != directions[i][s] for i, row in enumerate(initial) for s in row),
                tuple(tuple(ordinal[i][s] for s in row) for i, row in enumerate(orders)),
                tuple(orientation[i][s] for i, row in enumerate(initial) for s in row))

    def copy_state(state):
        return [list(row) for row in state[0]], [dict(row) for row in state[1]]

    baseline = ([list(row) for row in initial], [dict(row) for row in directions])
    before = measure(*baseline)
    best, best_score = copy_state(baseline), before
    evaluations = 1

    def consider(state):
        nonlocal best, best_score, evaluations
        score = measure(*state)
        evaluations += 1
        if score["total"] < best_score["total"] - 1e-12 or (abs(score["total"] - best_score["total"]) <= 1e-12 and key(state) < key(best)):
            best, best_score = copy_state(state), score
        return score

    estimate = 1
    for i, track in enumerate(tracks):
        if mode != "none" and not track.get("fixed_order", False):
            estimate *= math.factorial(len(units[i]))
            for group in units[i]:
                estimate *= math.factorial(len(group))
        if mode == "order_orientation":
            estimate *= 2 ** (len(initial[i]) - len(locks[i]))
        if estimate > EXACT_STATES:
            break
    exact = estimate <= EXACT_STATES
    if mode != "none" and exact:
        def states(i):
            group_orders = [units[i]] if tracks[i].get("fixed_order", False) else itertools.permutations(units[i])
            for arranged in group_orders:
                members = [[tuple(group)] if tracks[i].get("fixed_order", False) else itertools.permutations(group) for group in arranged]
                for parts in itertools.product(*members):
                    order = [s for part in parts for s in part]
                    free = [s for s in initial[i] if s not in locks[i]] if mode == "order_orientation" else []
                    for flips in itertools.product((False, True), repeat=len(free)):
                        orientation = dict(directions[i])
                        orientation.update(zip(free, flips, strict=True))
                        yield order, orientation
        # Materialize only the bounded per-track choices, not their joint product.
        for choice in itertools.product(*(list(states(i)) for i in range(len(tracks)))):
            consider(([row[0] for row in choice], [row[1] for row in choice]))
    elif mode != "none":
        rng = random.Random(seed)
        for restart in range(restarts):
            state = copy_state(baseline)
            if restart:
                for i, track in enumerate(tracks):
                    if not track.get("fixed_order", False):
                        groups = [list(group) for group in units[i]]
                        rng.shuffle(groups)
                        for group in groups:
                            rng.shuffle(group)
                        state[0][i] = [s for group in groups for s in group]
                    if mode == "order_orientation":
                        for sid in initial[i]:
                            if sid not in locks[i]:
                                state[1][i][sid] = bool(rng.getrandbits(1))
            current = consider(state)
            for _ in range(max_rounds):
                changed = False
                for i, track in enumerate(tracks):
                    candidates = []
                    if not track.get("fixed_order", False):
                        candidates.append(barycentric_order(state, i))
                        arranged = arranged_units(state[0][i], i)
                        for j in range(len(arranged) - 1):
                            order = arranged[:j] + [arranged[j + 1], arranged[j]] + arranged[j + 2:]
                            candidates.append([s for group in order for s in group])
                        for j, group in enumerate(arranged):
                            for k in range(len(group) - 1):
                                swapped = group[:k] + [group[k + 1], group[k]] + group[k + 2:]
                                order = arranged[:j] + [swapped] + arranged[j + 1:]
                                candidates.append([s for part in order for s in part])
                    # Best local move against the complete graph, not separate
                    # pairwise decisions for a middle track.
                    winner, score = None, current
                    for order in candidates:
                        trial = copy_state(state)
                        trial[0][i] = order
                        value = consider(trial)
                        if value["total"] < score["total"] - 1e-12:
                            winner, score = trial, value
                    if winner is not None:
                        state, current, changed = winner, score, True
                    if mode == "order_orientation":
                        # Accept supported flips sequentially against the current
                        # complete state, so ten independent improvements do not
                        # require ten outer rounds.
                        for sid in initial[i]:
                            if sid in locks[i]:
                                continue
                            trial = copy_state(state)
                            trial[1][i][sid] = not trial[1][i][sid]
                            score = consider(trial)
                            if score["total"] < current["total"] - 1e-12:
                                state, current, changed = trial, score, True
                if not changed:
                    break
    report = {"schema_version": 1, "method": "normalized_length_crossing_twist", "mode": mode,
              "solver": "preserve_input" if mode == "none" else ("exhaustive" if exact else "bounded_multistart_local_search"),
              "globally_optimal": mode != "none" and exact, "seed": seed, "restarts": restarts, "max_rounds": max_rounds,
              "evaluations": evaluations, "exact_state_limit": EXACT_STATES,
              "objective_before": before, "objective_after": best_score,
              "objective_definition": "sum of normalized width-weighted straight midpoint length, midpoint crossing-pair weight and oppositely directed endpoint-span weight",
              "block_weight": "mean endpoint width in native JCVI gene-rank coordinates; zero-width blocks retained with zero weight",
              "crossing_definition": "strict midpoint inversions per edge, excluding ties; not Bezier area overlap",
              "species_order": species, "input_order": dict(zip(species, initial, strict=True)),
              "display_order": dict(zip(species, best[0], strict=True)),
              "input_orientation": {s: {q: "-" if v else "+" for q, v in row.items()} for s, row in zip(species, directions, strict=True)},
              "display_orientation": {s: {q: "-" if v else "+" for q, v in row.items()} for s, row in zip(species, best[1], strict=True)},
              "orientation_changed": best[1] != directions, "coordinate_system": "gene_rank",
              "scale_mode": scale_mode, "track_ratios": ratios, "gaps": gaps, "xstart": xstart, "xend": xend,
              "figsize_inches": list(figsize), "ys": ys,
              "displayed_block_counts": [len(rows) for _, _, rows in blocks], "excluded_block_counts": excluded,
              "constraints": [{k: t[k] for k in ("groups", "fixed_order", "fixed_orientation") if k in t} for t in tracks]}
    return best[0], best[1], report

"""Width-weighted ribbon-length minimization in JCVI display space."""

import math
import random
from collections import Counter
from itertools import permutations
from pathlib import Path

import numpy as np

FIGSIZE = (20, 8)
XSTART, XEND = 0.12, 0.92
TRACK_DISTANCE = 0.4
EXACT_LIMIT = 20
JOINT_EXACT_WORK_LIMIT = 2_000_000
JOINT_RESTARTS = 6
JOINT_MAX_ROUNDS = 8


def track_geometry(seqids, genes):
    counts = Counter(g.seqid for g in genes)
    gap = min(0.01, 0.01 * 16 / len(seqids) + 0.001) if len(seqids) > 16 else 0.01
    span = XEND - XSTART - gap * (len(seqids) - 1)
    if span <= 0 or not seqids or any(counts[sid] == 0 for sid in seqids):
        raise ValueError("Invalid ribbon track geometry")
    ratio = span / sum(counts[sid] for sid in seqids)
    ranks = {}
    for sid in seqids:
        # JCVI Bed.order_in_chr sorts by start, then literal accession (not end
        # or natural accession order). Shared-start transcripts must agree too.
        ordered = sorted((g for g in genes if g.seqid == sid), key=lambda g: (g.start, g.gene_id))
        ranks.update((g.gene_id, rank) for rank, g in enumerate(ordered))
    widths = {sid: counts[sid] * ratio for sid in seqids}
    return widths, ranks, ratio, gap


def pair_track_geometry(selected, genomes, scale_mode="shared"):
    """Left-aligned tracks, optionally using the same width per gene."""
    if scale_mode not in {"shared", "independent"}:
        raise ValueError("karyotype-scale must be shared or independent")
    if len(selected) != 2 or len(genomes) != 2:
        raise ValueError("Pairwise geometry requires exactly two tracks")
    geometry = [track_geometry(seqids, genes) for seqids, genes in zip(selected, genomes, strict=True)]
    if scale_mode == "shared":
        ratio = min(track[2] for track in geometry)
        scaled = []
        for seqids, genes, (_, ranks, _, gap) in zip(selected, genomes, geometry, strict=True):
            counts = Counter(g.seqid for g in genes)
            scaled.append(({sid: counts[sid] * ratio for sid in seqids}, ranks, ratio, gap))
        geometry = scaled
    return geometry


def offsets(order, widths, gap):
    result = {}
    position = XSTART
    for sid in order:
        result[sid] = position
        position += widths[sid] + gap
    return result


def minimize_order(order, widths, gap, cost, exact_limit=EXACT_LIMIT):
    """Costs depend only on the preceding subset's width, not its internal order."""
    n = len(order)

    def objective(permutation):
        positions = offsets(permutation, widths, gap)
        return sum(float(cost(sid, positions[sid])) for sid in permutation)

    before = objective(order)
    if n <= exact_limit:
        states = 1 << n
        prefix = np.zeros(states)
        for mask in range(1, states):
            bit = mask & -mask
            prefix[mask] = prefix[mask ^ bit] + widths[order[bit.bit_length() - 1]] + gap
        costs = [np.asarray(cost(sid, XSTART + prefix)) for sid in order]
        distances = np.full(states, np.inf)
        distances[0] = 0
        previous = np.full(states, -1, dtype=np.int8)
        tie_codes = [0] * states
        for mask in range(1, states):
            bits = mask
            while bits:
                bit = bits & -bits
                i = bit.bit_length() - 1
                parent = mask ^ bit
                candidate = distances[parent] + costs[i][parent]
                code = tie_codes[parent] * (n + 1) + i + 1
                if candidate < distances[mask] or (candidate == distances[mask] and code < tie_codes[mask]):
                    distances[mask], previous[mask], tie_codes[mask] = candidate, i, code
                bits ^= bit
        mask = states - 1
        result = []
        while mask:
            i = int(previous[mask])
            result.append(order[i])
            mask ^= 1 << i
        result.reverse()
        method, optimal = "exact_subset_dynamic_programming", True
    else:
        # Explicitly reported local optimization for very fragmented assemblies.
        # No factorial/exponential allocation on hundreds of scaffolds.
        result = list(order)
        current = before
        for _ in range(n):
            improved = False
            for i in range(n - 1):
                candidate = list(result)
                candidate[i:i + 2] = reversed(candidate[i:i + 2])
                value = objective(candidate)
                if value < current:
                    result, current, improved = candidate, value, True
            if not improved:
                break
        method, optimal = "bounded_adjacent_swap", False
    return result, {"solver": method, "globally_optimal": optimal, "exact_limit": exact_limit,
                    "objective_before": before, "objective_after": objective(result)}


def order_by_ribbon_length(selected, genomes, simple, moving, scale_mode="shared"):
    if moving is None:
        return order_both_by_ribbon_length(selected, genomes, simple, scale_mode)
    fixed = 1 - moving
    geometry = pair_track_geometry(selected, genomes, scale_mode)
    widths, ranks, ratio, gap = geometry[moving]
    fixed_widths, fixed_ranks, fixed_ratio, fixed_gap = geometry[fixed]
    fixed_offsets = offsets(selected[fixed], fixed_widths, fixed_gap)
    lookup = [{g.gene_id: g.seqid for g in genes} for genes in genomes]
    connections = {sid: [] for sid in selected[moving]}
    for line in Path(simple).read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 6 or any(gene not in lookup[i // 2] for i, gene in enumerate(fields[:4])):
            raise ValueError(f"Invalid JCVI simple block: {line}")
        ends = [fields[:2], fields[2:4]]
        seqids = [lookup[i][ends[i][0]] for i in (0, 1)]
        if any(lookup[i][ends[i][1]] != seqids[i] for i in (0, 1)):
            raise ValueError("A ribbon spans multiple chromosomes")
        if seqids[moving] not in connections or seqids[fixed] not in fixed_offsets:
            continue
        a, b = (ranks[gene] for gene in ends[moving])
        c, d = (fixed_ranks[gene] for gene in ends[fixed])
        weight = (abs(a - b) * ratio + abs(c - d) * fixed_ratio) / 2
        connections[seqids[moving]].append((ratio * (a + b) / 2,
                                          fixed_offsets[seqids[fixed]] + fixed_ratio * (c + d) / 2, weight))

    def cost(sid, start):
        result = np.zeros_like(start, dtype=float)
        for local_center, fixed_center, weight in connections[sid]:
            result += weight * np.hypot((start + local_center - fixed_center) * FIGSIZE[0],
                                       TRACK_DISTANCE * FIGSIZE[1])
        return result

    ordered, solver = minimize_order(selected[moving], widths, gap, cost)
    result = [list(seqids) for seqids in selected]
    result[moving] = ordered
    metadata = {"method": "width_weighted_ribbon_length", "fixed_side": ("target", "query")[fixed],
                "weight": "mean_endpoint_width", "length": "straight_centerline_inches",
                "figsize_inches": list(FIGSIZE), "xstart": XSTART, "xend": XEND,
                "track_distance": TRACK_DISTANCE, "gaps": [g[3] for g in geometry],
                "scale_mode": scale_mode, "track_ratios": [g[2] for g in geometry],
                "track_spans": [sum(g[0].values()) + (len(g[0]) - 1) * g[3] for g in geometry],
                "displayed_block_count": sum(map(len, connections.values())),
                "orientation_changed": False, "input_order": dict(zip(("target", "query"), selected, strict=True)),
                "display_order": dict(zip(("target", "query"), result, strict=True)), **solver}
    if not math.isfinite(metadata["objective_after"]):
        raise ValueError("Non-finite ribbon objective")
    return result, metadata


def order_both_by_ribbon_length(selected, genomes, simple, scale_mode="shared"):
    """Optimize both tracks, recording the limits of a bounded joint search."""
    original = [list(order) for order in selected]
    indices = [{sid: i for i, sid in enumerate(order)} for order in original]

    def key(orders):
        return tuple(tuple(indices[side][sid] for sid in order) for side, order in enumerate(orders))

    best, metadata = order_by_ribbon_length(original, genomes, simple, 0, scale_mode)
    before = metadata["objective_before"]
    best_value = metadata["objective_after"]
    solves = 1

    def consider(orders, result_metadata):
        nonlocal best, best_value, metadata
        value = result_metadata["objective_after"]
        if value < best_value or (value == best_value and key(orders) < key(best)):
            best, best_value, metadata = orders, value, result_metadata

    # Exhaustively enumerate the smaller fixed track only for bounded joint
    # workloads. The other track is solved exactly for each permutation.
    candidates = []
    for fixed in (0, 1):
        n = len(original[1 - fixed])
        work = math.factorial(len(original[fixed])) * n * (1 << n) if n <= EXACT_LIMIT else math.inf
        candidates.append((work, fixed))
    work, fixed = min(candidates)
    if work <= JOINT_EXACT_WORK_LIMIT:
        for order in permutations(original[fixed]):
            start = [list(track) for track in original]
            start[fixed] = list(order)
            candidate, info = order_by_ribbon_length(start, genomes, simple, 1 - fixed, scale_mode)
            solves += 1
            consider(candidate, info)
        solver, optimal, starts, rounds = "exact_joint_permutation_subset_dp", True, 1, None
    else:
        rng = random.Random(1)
        starts = JOINT_RESTARTS
        rounds = JOINT_MAX_ROUNDS
        for restart in range(starts):
            current = [list(track) for track in original]
            if restart >= 2:
                for track in current:
                    rng.shuffle(track)
            last_value = math.inf
            for _ in range(rounds):
                prior = [list(track) for track in current]
                for moving in (restart % 2, 1 - restart % 2):
                    current, info = order_by_ribbon_length(current, genomes, simple, moving, scale_mode)
                    solves += 1
                    consider(current, info)
                value = info["objective_after"]
                if current == prior or value >= last_value:
                    break
                last_value = value
        solver, optimal = "multistart_alternating_track_optimization", False
    metadata.update(fixed_side=None, moving_sides=["target", "query"], solver=solver,
                    globally_optimal=optimal, objective_before=before, objective_after=best_value,
                    input_order=dict(zip(("target", "query"), original, strict=True)),
                    display_order=dict(zip(("target", "query"), best, strict=True)),
                    joint_exact_work_limit=JOINT_EXACT_WORK_LIMIT, restart_count=starts,
                    max_rounds=rounds, track_solves=solves, seed=1)
    return best, metadata

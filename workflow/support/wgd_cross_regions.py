"""Assess whether cross-synteny support depends on a single genomic region."""

from collections import Counter, defaultdict

REGION_SUMMARY = ("species_event_id", "branch_id", "mode", "matched_windows", "matched_families",
                  "regions", "largest_region_windows", "min_windows_after_one_region", "contigs_by_species",
                  "supports_before_holdout", "num_holdouts_losing_support", "holds_after_every_one_region")
REGION_DETAIL = ("species_event_id", "branch_id", "mode", "region", "matched_windows", "matched_families",
                 "num_original_blocks", "copy_contigs", "windows_after_exclusion", "supports_after_exclusion")
NESTED_FIELDS = ("younger_event_id", "younger_branch_id", "older_event_id", "older_branch_id", "species",
                 "younger_window_pairs", "older_window_pairs", "shared_window_pairs",
                 "younger_gene_pairs", "older_gene_pairs", "shared_gene_pairs",
                 "younger_families", "older_families", "shared_families")


def components(hits, reverse_windows=None):
    parents = list(range(len(hits)))

    def root(i):
        while parents[i] != i:
            parents[i] = parents[parents[i]]
            i = parents[i]
        return i

    seen = {}
    for i, (_, a, b, window, _, _, blocks, positions) in enumerate(hits):
        keys = (("scored_window", a, b, *window), *(("block", *v) for v in blocks),
                *(("copy_window", *v) for v in positions))
        if reverse_windows is not None:
            keys += tuple(("reversed_window", a, b, *v) for v in reverse_windows[hits[i]])
        for key in keys:
            parents[root(i)] = root(seen.setdefault(key, i))
    groups = defaultdict(list)
    for i, hit in enumerate(hits):
        groups[root(i)].append(hit)
    return sorted(groups.values(), key=lambda g: (-len({h[1:4] for h in g}), g[0]))


def support_after_holdouts(hits, daughters, floor, *, reverse_windows=None):
    groups = components(hits, reverse_windows)
    stats = []
    windows_by_pair, family_regions = Counter(), defaultdict(Counter)
    reverse_by_pair = Counter()
    for group in groups:
        windows = {h[1:4] for h in group}
        window_counts = Counter((a, b) for a, b, _ in windows)
        reverse = Counter((a, b) for a, b, _ in {
            (h[1], h[2], window) for h in group for window in reverse_windows[h]
        }) if reverse_windows is not None else Counter()
        reverse_by_pair.update(reverse)
        families = defaultdict(set)
        for _, a, b, _, family, *_ in group:
            families[a, b].add(family)
        windows_by_pair.update(window_counts)
        for pair, ids in families.items():
            family_regions[pair].update(ids)
        stats.append((group, window_counts, families, reverse))

    def pair_passes(windows, families, reverse):
        return min(windows, families) >= floor and (reverse_windows is None or reverse >= floor)

    initially_supported = {pair for pair, n in windows_by_pair.items()
                           if pair_passes(n, len(family_regions[pair]), reverse_by_pair[pair])}
    initial_taxa = Counter(sp for pair in initially_supported for sp in pair)

    def has_support(taxa):
        return all(any(taxa[sp] > 0 for sp in daughter) for daughter in daughters)

    total = windows_by_pair.total()
    result = []
    for i, (group, wc, families, rc) in enumerate(stats, 1):
        taxa = initial_taxa.copy()
        for pair in wc.keys() & initially_supported:
            lost_families = sum(family_regions[pair][f] == 1 for f in families[pair])
            if not pair_passes(windows_by_pair[pair] - wc[pair], len(family_regions[pair]) - lost_families,
                               reverse_by_pair[pair] - rc[pair]):
                taxa.subtract(pair)
        result.append({"region": i, "matched_windows": wc.total(), "matched_families": len({h[4] for h in group}),
                       "num_original_blocks": len({v for h in group for v in h[6]}),
                       "copy_contigs": ",".join(f"{sp}:{c}" for sp, c in sorted({v[:2] for h in group for v in h[7]})),
                       "windows_after_exclusion": total - wc.total(), "supports_after_exclusion": int(has_support(taxa))})
    return has_support(initial_taxa), result


def summarize_regions(candidates, traced_hits, daughters, floor, *, unique_phase_hits=None, reverse_windows=None):
    by_branch = defaultdict(list)
    for hit in sorted(traced_hits):
        by_branch[hit[0]].append(hit)
    unique_by_branch = defaultdict(list)
    if unique_phase_hits is not None:
        for hit in sorted(unique_phase_hits):
            unique_by_branch[hit[0]].append(hit)
    shared = defaultdict(list)
    if reverse_windows is not None:
        for hit in sorted(reverse_windows):
            shared[hit[0]].append(hit)
    summary, details = [], []
    for row in candidates:
        branch = row["species_event_id"]
        if branch not in daughters:
            continue
        hits = by_branch[branch]
        modes = [("paired", [h for h in hits if len(h[5]) == 4]), ("all", hits)]
        if unique_phase_hits is not None:
            modes.append(("unique_phase", unique_by_branch[branch]))
        if reverse_windows is not None:
            modes.append(("both_orientations", shared[branch]))
        for mode, selected in modes:
            supported, regions = support_after_holdouts(selected, daughters[branch], floor,
                                                         reverse_windows=reverse_windows if mode == "both_orientations" else None)
            contigs = defaultdict(set)
            for hit in selected:
                for sp, contig, _ in hit[7]:
                    contigs[sp].add(contig)
            keys = {"species_event_id": branch, "branch_id": row["branch_id"], "mode": mode}
            lost = sum(not r["supports_after_exclusion"] for r in regions) if supported else 0
            summary.append({**keys, "matched_windows": sum(r["matched_windows"] for r in regions),
                            "matched_families": len({h[4] for h in selected}), "regions": len(regions),
                            "largest_region_windows": regions[0]["matched_windows"] if regions else 0,
                            "min_windows_after_one_region": min((r["windows_after_exclusion"] for r in regions), default=0),
                            "contigs_by_species": ",".join(f"{sp}:{len(contigs[sp])}" for sp in sorted(contigs)),
                            "supports_before_holdout": int(supported), "num_holdouts_losing_support": lost,
                            "holds_after_every_one_region": int(supported and bool(regions) and lost == 0)})
            details.extend({**keys, **r} for r in regions)
    return summary, details


def nested_copy_overlap(tree, candidates, traced_hits):
    """Compare exact paired paths only, in species sampled at both nested branches."""
    from nwkit.clade_index import CladeIndex

    index = CladeIndex(tree)
    supported = {r["species_event_id"]: r["branch_id"] for r in candidates
                 if r["tree_support"] == "tree_score_supported"}
    records = defaultdict(dict)
    for event, a, b, _, family, genes, _, positions in traced_hits:
        if len(genes) != 4 or event not in supported:
            continue
        for sp, i in ((a, 0), (b, 1)):
            windows, pairs, families = records[event].setdefault(sp, (set(), set(), set()))
            windows.add(tuple(sorted((positions[i][1:], positions[i + 2][1:]))))
            pairs.add(tuple(sorted((genes[i], genes[i + 2]))))
            families.add(family)
    for node in tree.descendants("preorder"):
        young = index.clade_id_for_node(node)
        if young not in records:
            continue
        for parent in node.ancestors():
            old = index.clade_id_for_node(parent)
            if old not in records:
                continue
            for sp in sorted(records[young].keys() & records[old].keys()):
                counts = [v for x, y in zip(records[young][sp], records[old][sp], strict=True)
                          for v in (len(x), len(y), len(x & y))]
                yield dict(zip(NESTED_FIELDS, (young, supported[young], old, supported[old], sp, *counts), strict=True))

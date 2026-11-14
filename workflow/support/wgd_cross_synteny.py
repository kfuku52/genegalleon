"""Branch evidence from sparse, precomputed self and cross anchor blocks."""

import csv
import gzip
import itertools
from collections import defaultdict
from pathlib import Path

try:
    from wgd_evidence import read_table
except ImportError:
    from .wgd_evidence import read_table


ANCHOR_FLOOR = 5
WINDOW_GENES = 100
PAIR_FIELDS = ("species_a", "species_b", "bed_a", "bed_b", "anchors")
EMPTY = frozenset()
PAIR_OUTPUT = ("species_event_id", "species_a", "species_b", "status", "observable_windows",
               "paired_windows", "partial_only_windows", "matched_paired_windows",
               "matched_partial_only_windows", "matched_families")
WINDOW_OUTPUT = ("species_event_id", "species_a", "species_b", "contig", "window",
                 "paired", "partial_only", "matched_paired_families", "matched_partial_only_families")
COVERAGE_OUTPUT = ("species_event_id", "species_a", "species_b", "num_cross_anchor_pairs",
                   "num_archive_family_anchor_pairs", "num_observable_windows", "num_archive_family_windows",
                   "archive_family_window_floor_met")
NONREUSE_OUTPUT = ("species_event_id", "species_a", "species_b", "num_observable_windows",
                   "num_nonreused_anchor_windows", "num_archive_family_windows", "num_nonreused_archive_family_windows",
                   "nonreused_anchor_window_floor_met", "nonreused_archive_window_floor_met")
ORIENTATION_METRICS = ("status", "observable_windows", "matched_paired_windows",
                       "matched_partial_only_windows", "matched_families")
ORIENTATION_OUTPUT = ("species_event_id", "branch_id", "species_a", "species_b",
                      *(f"{side}_{field}" for side in ("a", "b") for field in ORIENTATION_METRICS),
                      "a_supports_pair", "b_supports_pair", "status")


def read_pairs(path, resolve, species, hash_input):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        header = next(csv.reader(stream, delimiter="\t"), [])
    if sorted(header) != sorted(PAIR_FIELDS):
        raise ValueError(f"Expected cross-synteny manifest columns: {', '.join(PAIR_FIELDS)}")
    beds, seen, result = {}, set(), []
    for row in read_table(path):
        a, b = row["species_a"], row["species_b"]
        if a not in species or b not in species:
            raise ValueError("Cross-synteny species must match exact species-tree tips")
        key = tuple(sorted((a, b)))
        if key in seen:
            raise ValueError(f"Duplicate cross-synteny pair: {key}")
        seen.add(key)
        for sp, field in ((a, "bed_a"), (b, "bed_b"), (None, "anchors")):
            if not row[field]:
                raise ValueError(f"Missing {field} for cross-synteny pair {key}")
            source = resolve(row[field])
            row[field] = str(source)
            if sp is not None and beds.setdefault(sp, source) != source:
                raise ValueError(f"Use one coordinate BED per cross-synteny species: {sp}")
            hash_input(source)
        result.append(row)
    if not result:
        raise ValueError("Cross-synteny manifest has no comparisons")
    return result


def coordinates(path, species):
    entries = {}
    with Path(path).open(encoding="utf-8", newline="") as stream:
        for row in csv.reader(stream, delimiter="\t", strict=True):
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 4:
                raise ValueError(f"Malformed cross-synteny BED row in {path}")
            chrom, start, end, gene = row[:4]
            if (not chrom or not gene.startswith(species + "_") or gene == species + "_"
                    or gene in entries or not 0 <= int(start) < int(end)):
                raise ValueError(f"Invalid or duplicate cross-synteny BED locus in {path}: {gene}")
            entries[gene] = chrom, int(start), int(end)
    if not entries:
        raise ValueError(f"Empty cross-synteny BED: {path}")
    ranks = defaultdict(int)
    result = {}
    for gene in sorted(entries, key=lambda name: (*entries[name], name)):
        chrom = entries[gene][0]
        result[gene] = chrom, ranks[chrom]
        ranks[chrom] += 1
    return result


def anchor_links(path, a, b, symmetric=False, *, reversed_input=False):
    links, neighbors = defaultdict(set), defaultdict(set)
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        block, contigs = 0, None
        for line in stream:
            line = line.strip()
            if line == "###":
                block, contigs = block + 1, None
            elif line and not line.startswith("#"):
                fields = line.split("\t")
                if len(fields) < 2 or not block:
                    raise ValueError(f"Anchor outside its block or BED in {path}: {line}")
                x, y = fields[:2]
                if reversed_input:
                    x, y = y, x
                if x not in a or y not in b:
                    raise ValueError(f"Anchor outside its block or BED in {path}: {line}")
                pair = a[x][0], b[y][0]
                if contigs is not None and contigs != pair:
                    raise ValueError(f"Anchor block spans multiple contig pairs in {path}")
                contigs = pair
                if symmetric and x == y:
                    continue
                links[x, y].add(block)
                neighbors[x].add(y)
                if symmetric:
                    links[y, x].add(block)
                    neighbors[y].add(x)
    return links, neighbors


def window_for(positions, gene):
    chrom, rank = positions[gene]
    return chrom, rank // WINDOW_GENES + 1


def qualifying(groups):
    for key, tuples in groups.items():
        if len(tuples) >= ANCHOR_FLOOR and all(len(set(slot)) >= ANCHOR_FLOOR for slot in zip(*tuples, strict=True)):
            yield key, tuples


def cross_window_groups(positions, cross):
    groups = defaultdict(set)
    for (x, y), ids in cross.items():
        window = window_for(positions, x)
        for block in ids:
            groups[window, block].add((x, y))
    return groups


def observable_windows(positions, cross):
    return {window for (window, _), _ in qualifying(cross_window_groups(positions, cross))}


def has_nonreused_anchors(tuples):
    # Most blocks already contain five gene-disjoint pairs. Fall back to an
    # exact bipartite match when an inexpensive greedy witness is insufficient.
    pairs = sorted(tuples)
    if len(pairs) < ANCHOR_FLOOR:
        return False
    used_a, used_b = set(), set()
    for a, b in pairs:
        if a not in used_a and b not in used_b:
            used_a.add(a)
            used_b.add(b)
            if len(used_a) >= ANCHOR_FLOOR:
                return True
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import maximum_bipartite_matching

    row, col = {}, {}
    indices = [(row.setdefault(a, len(row)), col.setdefault(b, len(col))) for a, b in pairs]
    ii, jj = zip(*indices, strict=True)
    graph = csr_matrix(([1] * len(indices), (ii, jj)), shape=(len(row), len(col)))
    return bool((maximum_bipartite_matching(graph, perm_type="column") >= 0).sum() >= ANCHOR_FLOOR)


def correspondence(positions, left, right, cross, provenance=None):
    aa, an = left
    bb, bn = right
    ab, cross_neighbors = cross
    groups = defaultdict(set)
    # Both paths must use the same two self blocks and two direct cross blocks.
    for x, ys in cross_neighbors.items():
        for xx in an.get(x, ()):
            if positions[x] >= positions[xx] or xx not in cross_neighbors:
                continue
            for y in ys:
                for yy in bn.get(y, EMPTY) & cross_neighbors[xx]:
                    for blocks in itertools.product(aa[x, xx], bb[y, yy], ab[x, y], ab[xx, yy]):
                        groups[window_for(positions, x), blocks].add((x, y, xx, yy))
    paired = defaultdict(set)
    for (window, blocks), genes in qualifying(groups):
        paired[window].update(genes)
        if provenance is not None:
            keys = tuple(zip((0, 1, 2, 2), blocks, strict=True))
            for gene_tuple in genes:
                provenance[window, gene_tuple].add(keys)

    partial = defaultdict(set)
    for reverse, (self_links, self_neighbors), (_, other_neighbors) in ((False, left, right), (True, right, left)):
        oriented = defaultdict(dict)
        for (x, y), ids in ab.items():
            a, b = (y, x) if reverse else (x, y)
            oriented[a][b] = ids
        groups = defaultdict(set)
        for x, others in oriented.items():
            for xx in self_neighbors.get(x, ()):
                other_path = oriented.get(xx, {})
                for y, cross_blocks in others.items():
                    window = window_for(positions, y if reverse else x)
                    if window in paired or not other_neighbors.get(y, EMPTY).isdisjoint(other_path):
                        continue
                    for blocks in itertools.product(self_links[x, xx], cross_blocks):
                        groups[window, blocks].add((x, y, xx))
        for (window, blocks), genes in qualifying(groups):
            partial[window].update(genes)
            if provenance is not None:
                keys = ((1 if reverse else 0, blocks[0]), (2, blocks[1]))
                for gene_tuple in genes:
                    provenance[window, gene_tuple].add(keys)
    return paired, partial


def cross_synteny_inputs(pairs):
    """Parse shared coordinates and self links once, without retaining cross links."""
    beds, self_blocks = {}, {}
    seen_genes = set()
    for row in pairs:
        for side in ("a", "b"):
            sp = row[f"species_{side}"]
            if sp not in beds:
                beds[sp] = coordinates(row[f"bed_{side}"], sp)
                if not seen_genes.isdisjoint(beds[sp]):
                    raise ValueError(f"Cross-synteny BED gene IDs are shared between species: {sp}")
                seen_genes.update(beds[sp])
        a, b = row["species_a"], row["species_b"]
        if a == b:
            self_blocks[a] = anchor_links(row["anchors"], beds[a], beds[a], symmetric=True)
    return beds, self_blocks


class CrossSynteny:
    """Collect exact copy-level agreement with the same duplication nodes used for scoring."""

    def __init__(self, tree, pairs, min_windows, age_species=(), trace_branches=(), audit_coverage=False,
                 *, reverse_cross=False, inputs=None):
        from nwkit.clade_index import CladeIndex

        self.min_windows = min_windows
        self.entries = defaultdict(list)
        self.matches = defaultdict(lambda: defaultdict(list))
        self.age_species = set(age_species)
        self.age_pairs = set()
        self.traced_hits = set()
        self.coverage = [] if audit_coverage else None
        self.coverage_genes = set() if audit_coverage else None
        self.archive_families = {}
        trace_branches = None if trace_branches is None else set(trace_branches)
        index = CladeIndex(tree)
        leaves = {leaf.name: leaf for leaf in tree.leaves()}
        beds, self_blocks = cross_synteny_inputs(pairs) if inputs is None else inputs
        self.daughters = {}
        for row in pairs:
            a, b = row["species_a"], row["species_b"]
            if a == b:
                continue
            if reverse_cross:
                a, b = b, a
            node = tree.common_ancestor([leaves[a], leaves[b]])
            cross = anchor_links(row["anchors"], beds[a], beds[b], reversed_input=reverse_cross)
            if node.is_root:
                continue
            branch = index.clade_id_for_node(node)
            self.daughters[branch] = tuple(set(child.leaf_names()) for child in node.children)
            observable = observable_windows(beds[a], cross[0])
            if self.coverage is not None:
                self.coverage.append((branch, a, b, beds[a], cross[0], len(observable)))
                self.coverage_genes.update(gene for pair in cross[0] for gene in pair)
            provenance = defaultdict(set) if trace_branches is None or branch in trace_branches else None
            if a in self_blocks and b in self_blocks:
                paired, partial = correspondence(beds[a], self_blocks[a], self_blocks[b], cross, provenance)
                status = "assessed" if len(observable) >= min_windows else "insufficient_cross_anchors"
            else:
                paired, partial = {}, {}
                status = "missing_self_comparison"
            entry = (a, b, status, observable, paired.keys(), partial.keys(), defaultdict(set), defaultdict(set))
            self.entries[branch].append(entry)
            age_pair = (a, b) if status == "assessed" and self.age_species.intersection((a, b)) else None
            for regions, observed, age in ((paired, entry[6], age_pair), (partial, entry[7], None)):
                for window, groups in regions.items():
                    for genes in groups:
                        regions_for_age = None if age is None else tuple(
                            (sp, tuple(sorted((x, y))), window_for(beds[sp], x), window_for(beds[sp], y))
                            for sp, x, y in zip(age, genes[:2], genes[2:], strict=True) if sp in self.age_species)
                        trace = (a, b, provenance[window, genes], beds) if provenance is not None and status == "assessed" else None
                        self.matches[branch][genes[0]].append((genes, window, observed, regions_for_age, trace))

    def observe_family(self, family, leaves):
        for gene in self.coverage_genes.intersection(leaves):
            self.archive_families[gene] = family

    def _covered_windows(self, audit_nonreuse):
        for event, a, b, positions, links, observable in self.coverage:
            both = {pair: blocks for pair, blocks in links.items() if pair[0] in self.archive_families
                    and self.archive_families[pair[0]] == self.archive_families.get(pair[1])}
            grouped = cross_window_groups(positions, both)
            eligible = list(qualifying(grouped)) if audit_nonreuse else qualifying(grouped)
            n = len({window for (window, _), _ in eligible})
            coverage = dict(zip(COVERAGE_OUTPUT, (event, a, b, len(links), len(both), observable, n,
                                                  int(n >= self.min_windows)), strict=True))
            nonreuse = None
            if audit_nonreuse:
                raw_count, archive_count = (
                    len({window for (window, _), pairs in rows if has_nonreused_anchors(pairs)})
                    for rows in (qualifying(cross_window_groups(positions, links)), eligible))
                nonreuse = dict(zip(NONREUSE_OUTPUT, (event, a, b, observable, raw_count, n, archive_count,
                                                     int(raw_count >= self.min_windows),
                                                     int(archive_count >= self.min_windows)), strict=True))
            yield coverage, nonreuse

    def coverage_rows(self):
        for coverage, _ in self._covered_windows(False):
            yield coverage

    def nonreused_coverage_rows(self):
        coverage, nonreuse = [], []
        for old, strict in self._covered_windows(True):
            coverage.append(old)
            nonreuse.append(strict)
        return coverage, nonreuse

    def observe(self, branch, family, node):
        targets = self.matches.get(branch)
        if not targets:
            return
        left, right = (set(child.leaf_names()) for child in node.children)
        for source, other in ((left, right), (right, left)):
            for gene in source.intersection(targets):
                for genes, window, observed, regions_for_age, trace in targets[gene]:
                    if genes[1] in source and all(g in other for g in genes[2:]):
                        observed[window].add(family)
                        for sp, pair, first, second in regions_for_age or ():
                            self.age_pairs.add((branch, sp, family, pair, *sorted((first, second))))
                        if trace is not None:
                            a, b, paths, beds = trace
                            scopes = ((a, a), (b, b), (a, b))
                            blocks = tuple(sorted((*scopes[side], block) for side, block in set().union(*paths)))
                            positions = tuple((sp, *window_for(beds[sp], gene)) for gene in genes
                                              for sp in (a if gene in beds[a] else b,))
                            self.traced_hits.add((branch, a, b, window, family, genes, blocks, positions))

    def summarize(self, *, include_windows=True):
        branches, pairs, windows = {}, [], []
        for branch, entries in self.entries.items():
            observed_species, supported_species = set(), set()
            for a, b, status, observable, paired, partial, full_hits, partial_hits in entries:
                assert (paired | partial) <= observable
                assert full_hits.keys() <= paired and partial_hits.keys() <= partial
                families = set().union(*full_hits.values(), *partial_hits.values())
                if status == "assessed":
                    observed_species.update((a, b))
                    if len(full_hits) + len(partial_hits) >= self.min_windows and len(families) >= self.min_windows:
                        supported_species.update((a, b))
                key = {"species_event_id": branch, "species_a": a, "species_b": b}
                pairs.append({**key, "status": status, "observable_windows": len(observable),
                              "paired_windows": len(paired), "partial_only_windows": len(partial),
                              "matched_paired_windows": len(full_hits), "matched_partial_only_windows": len(partial_hits),
                              "matched_families": len(families)})
                if include_windows:
                    for chrom, window in sorted(paired | partial):
                        k = chrom, window
                        windows.append({**key, "contig": chrom, "window": window, "paired": int(k in paired),
                                        "partial_only": int(k in partial),
                                        "matched_paired_families": ",".join(sorted(full_hits.get(k, ()))),
                                        "matched_partial_only_families": ",".join(sorted(partial_hits.get(k, ())))})
            daughters = self.daughters[branch]
            branches[branch] = {
                "observable": all(group & observed_species for group in daughters),
                "supported": all(group & supported_species for group in daughters),
                "observed_species": sorted(observed_species),
                "supported_species": sorted(supported_species),
            }
        return branches, pairs, windows


def orientation_rows(candidates, forward, reverse, floor):
    ids = {row["species_event_id"]: row["branch_id"] for row in candidates}
    reversed_by_pair = {(row["species_event_id"], row["species_b"], row["species_a"]): row for row in reverse}

    def passes(row):
        return (row["matched_paired_windows"] + row["matched_partial_only_windows"] >= floor
                and row["matched_families"] >= floor)

    for row in forward:
        event, a, b = (row[field] for field in ("species_event_id", "species_a", "species_b"))
        first, second = row, reversed_by_pair[event, a, b]
        assessed = tuple(item["status"] == "assessed" for item in (first, second))
        support = tuple(passes(item) if ok else "NA" for item, ok in zip((first, second), assessed, strict=True))
        status = ("not_assessed" if not any(assessed) else "orientation_dependent_observability" if not all(assessed)
                  else "orientation_dependent_support" if support[0] != support[1] else
                  "supports_in_both" if all(support) else "below_threshold_in_both")
        values = (event, ids[event], a, b, *(item[field] for item in (first, second) for field in ORIENTATION_METRICS),
                  *(int(value) if isinstance(value, bool) else value for value in support), status)
        yield dict(zip(ORIENTATION_OUTPUT, values, strict=True))

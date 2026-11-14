"""Check tree-matched paired anchors against split support from a separate archive."""

from collections import Counter, defaultdict
from contextlib import ExitStack
from pathlib import PurePosixPath
from zipfile import ZipFile

SUPPORT_FIELDS = ("species_event_id", "branch_id", "species_a", "species_b", "contig", "window",
                  "family_id", "genes", "status", "duplication_support", "first_copy_support", "second_copy_support")
COMPARISON_STATES = ("archive_missing", "different_leaves", "split_not_supported_or_missing", "all_splits_present")
SUPPORT_LEVELS = (50, 70, 95)
DUPLICATION_SUPPORT_FIELDS = ("species_event_id", "branch_id", "num_eligible_families", "num_supporting_families",
                              *(f"num_{state}" for state in COMPARISON_STATES),
                              *(f"num_all_splits_ge_{level}" for level in SUPPORT_LEVELS))


def archive_index(archive, suffix):
    result = {}
    for member in archive.infolist():
        name = PurePosixPath(member.filename).name
        if member.is_dir() or not name.endswith(suffix):
            continue
        family = name.removesuffix(suffix)
        if not family or family in result:
            raise ValueError(f"Duplicate or empty family in tree archive: {name}")
        result[family] = member
    if not result:
        raise ValueError(f"Tree archive contains no *{suffix} members")
    return result


def masks(tree, bits):
    result = {}
    for node in tree.traverse("postorder"):
        if node.is_leaf:
            result[node] = bits[node.name]
        else:
            result[node] = 0
            for child in node.children:
                result[node] |= result[child]
    return result


def split_id(mask, universe):
    return min(mask, universe ^ mask)


def support_splits(tree, bits):
    from math import isfinite

    by_node = masks(tree, bits)
    universe = by_node[tree]
    support = {}
    for node, mask in by_node.items():
        if node.is_root or node.is_leaf or "support" not in node.props:
            continue
        value = float(node.props["support"])
        if not isfinite(value) or not 0 <= value <= 100:
            raise ValueError("UFBoot support must be between 0 and 100")
        key = split_id(mask, universe)
        # A degree-two root can represent the same unrooted split twice.
        support[key] = value if key not in support or support[key] == value else None
    return support


def comparable_splits(archive, member, bits):
    from nwkit.util import read_tree

    if member is None:
        return "archive_missing", {}
    other = read_tree(archive.read(member).decode("utf-8"), 0, True, quiet=True)
    names = list(other.leaf_names())
    if len(names) != len(bits) or set(names) != bits.keys():
        return "different_leaves", {}
    return "comparable", support_splits(other, bits)


class DuplicationTreeSupport:
    """Stream one support tree per voting family; retain only branch counters."""

    def __init__(self, path):
        self.path = path
        self.counts = defaultdict(Counter)
        self.paired = defaultdict(dict)

    def __enter__(self):
        with ExitStack() as stack:
            self.archive = stack.enter_context(ZipFile(self.path))
            self.index = archive_index(self.archive, "_generax.ufboot.nwk")
            self.stack = stack.pop_all()
        return self

    def __exit__(self, *exc):
        return self.stack.__exit__(*exc)

    def observe(self, family, tree, duplications, paired_hits=()):
        leaves = {leaf.name: leaf for leaf in tree.leaves()}
        bits = {name: 1 << i for i, name in enumerate(sorted(leaves))}
        state, support = comparable_splits(self.archive, self.index.get(family), bits)
        by_node = masks(tree, bits) if support or paired_hits else {}
        best = dict.fromkeys(branch for branch, _ in duplications)
        if support:
            universe = by_node[tree]
            for branch, node in duplications:
                values = [support.get(split_id(by_node[n], universe)) for n in (node, *node.children)]
                if None not in values:
                    weakest = min(values)
                    best[branch] = weakest if best[branch] is None else max(best[branch], weakest)
        for branch, weakest in best.items():
            key = state if state != "comparable" else (
                "all_splits_present" if weakest is not None else "split_not_supported_or_missing")
            self.counts[branch][key] += 1
            if weakest is not None:
                for level in SUPPORT_LEVELS:
                    self.counts[branch][level] += weakest >= level
        for hit in paired_hits:
            nodes = paired_copy_nodes(tree, leaves, by_node, bits, hit[5])
            values = tuple(support.get(split_id(by_node[node], by_node[tree])) for node in nodes)
            status = state if state != "comparable" else (
                "all_splits_present" if None not in values else "split_not_supported_or_missing")
            self.paired[family][hit] = status, *values

    def rows(self, candidates):
        for row in candidates:
            branch = row["species_event_id"]
            counter = self.counts[branch]
            yield {**{key: row[key] for key in DUPLICATION_SUPPORT_FIELDS[:4]},
                   **{f"num_{state}": counter[state] for state in COMPARISON_STATES},
                   **{f"num_all_splits_ge_{level}": counter[level] for level in SUPPORT_LEVELS}}

    def paired_rows(self, candidates):
        branch_ids = {row["species_event_id"]: row["branch_id"] for row in candidates}
        for family, hits in sorted(self.paired.items()):
            for (event, a, b, (contig, window), _, genes, *_), values in sorted(hits.items()):
                yield dict(zip(SUPPORT_FIELDS, (event, branch_ids[event], a, b, contig, window, family,
                                               ",".join(genes), *values), strict=True))


def paired_copy_nodes(tree, leaves, by_node, bits, genes):
    x, y, xx, yy = (leaves[gene] for gene in genes)
    duplication = tree.common_ancestor([x, xx])
    first = bits[x.name] | bits[y.name]
    second = bits[xx.name] | bits[yy.name]
    copies = [next((c for c in duplication.children if by_node[c] & needed == needed), None)
              for needed in (first, second)]
    if None in copies or copies[0] is copies[1]:
        raise ValueError("Traced paired anchors do not split into the recorded duplication copies")
    return duplication, *copies

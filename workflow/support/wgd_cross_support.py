"""Check tree-matched paired anchors against split support from a separate archive."""

from collections import defaultdict
from pathlib import PurePosixPath
from zipfile import ZipFile

SUPPORT_FIELDS = ("species_event_id", "branch_id", "species_a", "species_b", "contig", "window",
                  "family_id", "genes", "status", "duplication_support", "first_copy_support", "second_copy_support")


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


def paired_tree_support(tree_archive, support_archive, candidates, traced_hits):
    from nwkit.util import read_tree

    branch_ids = {row["species_event_id"]: row["branch_id"] for row in candidates}
    by_family = defaultdict(list)
    for hit in sorted(traced_hits):
        if len(hit[5]) == 4:
            by_family[hit[4]].append(hit)
    with ZipFile(tree_archive) as primary, ZipFile(support_archive) as alternative:
        primary_index = archive_index(primary, "_root.nwk")
        support_index = archive_index(alternative, "_generax.ufboot.nwk")
        for family, hits in sorted(by_family.items()):
            tree = read_tree(primary.read(primary_index[family]).decode("utf-8"), "auto", True, quiet=True)
            leaves = {leaf.name: leaf for leaf in tree.leaves()}
            bits = {name: 1 << i for i, name in enumerate(sorted(leaves))}
            by_node = masks(tree, bits)
            universe = by_node[tree]
            support = {}
            state = "archive_missing"
            if family in support_index:
                other = read_tree(alternative.read(support_index[family]).decode("utf-8"), 0, True, quiet=True)
                names = list(other.leaf_names())
                state = "different_leaves"
                if len(names) == len(leaves) and set(names) == leaves.keys():
                    support = support_splits(other, bits)
                    state = "comparable"
            for event, a, b, (contig, window), _, genes, *_ in hits:
                nodes = paired_copy_nodes(tree, leaves, by_node, bits, genes)
                vals = tuple(support.get(split_id(by_node[node], universe)) for node in nodes)
                status = state if state != "comparable" else (
                    "all_splits_present" if None not in vals else "split_not_supported_or_missing")
                yield dict(zip(SUPPORT_FIELDS, (event, branch_ids[event], a, b, contig, window, family,
                                               ",".join(genes), status, *vals), strict=True))

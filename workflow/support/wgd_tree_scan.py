"""Branch-matched duplication scores from an archived set of rooted gene trees."""

import re
from collections import Counter
from pathlib import PurePosixPath
from zipfile import ZipFile

try:
    from wgd_taxon_holdout import DuplicationTaxonHoldout
except ImportError:
    from .wgd_taxon_holdout import DuplicationTaxonHoldout


def scan_duplications(species_tree, archive_path, minimum_score, observe=None, family_members=None,
                      observe_family=None, observe_duplications=None, *, holdout_output=None):
    from nwkit.clade_index import CladeIndex
    from nwkit.rooting_state import require_rooted
    from nwkit.util import assign_branch_ids, read_tree

    index = CladeIndex(species_tree)
    branch_ids = assign_branch_ids(species_tree)
    taxa = sorted(species_tree.leaf_names())
    bits = {name: 1 << i for i, name in enumerate(taxa)}
    masks = {}
    for node in species_tree.traverse("postorder"):
        masks[node] = bits[node.name] if node.is_leaf else sum(masks[child] for child in node.children)
    branches = {}
    for node in species_tree.descendants("preorder"):
        sisters = masks[node.up] ^ masks[node]
        if not node.is_leaf:
            branches[masks[node]] = (node, sisters, tuple(masks[c] for c in node.children),
                                    index.clade_id_for_node(node))
    holdout = DuplicationTaxonHoldout(taxa, branches) if holdout_output is not None else None
    prefix = re.compile(r"^(" + "|".join(map(re.escape, sorted(taxa, key=lambda name: (-len(name), name)))) + r")_(?=.)")
    eligible, hits = Counter(), Counter()
    families = set()
    with ZipFile(archive_path) as archive:
        for member in archive.infolist():
            name = PurePosixPath(member.filename).name
            if member.is_dir() or not name.endswith("_root.nwk"):
                continue
            family = name.removesuffix("_root.nwk")
            if not family or family in families:
                raise ValueError(f"Duplicate or empty gene-tree family: {name}")
            families.add(family)
            tree = read_tree(archive.read(member).decode("utf-8"), "auto", True, rooted="auto", quiet=True)
            require_rooted(tree, f"WGD duplication score requires rooted gene trees: {name}")
            leaves = set()
            expected = family_members(family) if family_members is not None else None
            if family_members is not None and expected is None:
                raise ValueError(f"Gene-tree family absent from membership table: {name}")
            observed, supported = {}, set()
            duplication_nodes = []
            for gene in tree.traverse("postorder"):
                if gene.is_leaf:
                    if gene.name in leaves:
                        raise ValueError(f"Duplicate gene or unmatched species in {name}: {gene.name}")
                    if expected is not None:
                        if gene.name not in expected:
                            raise ValueError(f"Gene absent from its family membership in {name}: {gene.name}")
                        taxon = expected[gene.name]
                    else:
                        match = prefix.match(gene.name or "")
                        if match is None:
                            raise ValueError(f"Duplicate gene or unmatched species in {name}: {gene.name}")
                        taxon = match[1]
                    leaves.add(gene.name)
                    observed[gene] = bits[taxon]
                    continue
                union = 0
                for child in gene.children:
                    union |= observed[child]
                observed[gene] = union
                if union not in branches or len(gene.children) != 2:
                    continue
                left, right = (observed[c] for c in gene.children)
                if left & right and all(left & daughter and right & daughter for daughter in branches[union][2]):
                    supported.add(union)
                    if observe is not None or observe_duplications is not None:
                        duplication_nodes.append((union, gene))
            full = observed[tree]
            if observe_family is not None:
                observe_family(family, leaves)
            for mask, (_, sisters, daughters, _) in branches.items():
                if full & sisters and all(full & daughter for daughter in daughters):
                    eligible[mask] += 1
                    if mask in supported:
                        hits[mask] += 1
            if holdout is not None:
                holdout.observe(tree, observed, supported)
            if observe is not None or observe_duplications is not None:
                family_nodes = []
                for mask, gene in duplication_nodes:
                    _, sisters, _, branch = branches[mask]
                    if full & sisters:
                        if observe is not None:
                            observe(branch, family, gene)
                        if observe_duplications is not None:
                            family_nodes.append((branch, gene))
                if family_nodes:
                    observe_duplications(family, tree, family_nodes)
    if not families:
        raise ValueError("Gene-tree archive contains no *_root.nwk families")
    rows = []
    for mask, (node, _, _, branch) in branches.items():
        n, support = eligible[mask], hits[mask]
        score = support / n if n else 0.0
        rows.append({
            "branch_id": branch_ids[node], "species_event_id": branch,
            "descendant_taxa": index.csv_for_node(node), "num_eligible_families": n,
            "num_supporting_families": support, "duplication_score": score,
            "minimum_duplication_score": minimum_score,
            "tree_support": "tree_score_supported" if n and score >= minimum_score else "below_threshold",
        })
    rows.sort(key=lambda r: (-r["duplication_score"], -r["num_supporting_families"], r["branch_id"]))
    if holdout is not None:
        holdout_output.extend(holdout.rows(rows))
    return rows, len(families)

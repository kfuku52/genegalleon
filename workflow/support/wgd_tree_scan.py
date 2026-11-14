"""Branch-matched duplication scores from an archived set of rooted gene trees."""

import re
from collections import Counter
from pathlib import PurePosixPath
from zipfile import ZipFile

try:
    from wgd_taxon_holdout import DuplicationTaxonHoldout
except ImportError:
    from .wgd_taxon_holdout import DuplicationTaxonHoldout

COVERAGE_FIELDS = ("species_event_id", "branch_id", "eligible_complete", "eligible_partial", "support_complete",
                   "support_partial_projected", "complete_case_score", "partial_case_score", "projected_score")
ROOT_FIELDS = ("species_event_id", "branch_id", "eligible_families", "rooted_supporting_families",
               "extra_families_under_alternate_root", "extra_families_with_other_rooted_branch",
               "any_orientation_score")


def alternate_root_votes(tree, observed, eligible_daughters):
    """Check the two possible child groups for each direction at a three-way junction."""
    above = {tree: 0}
    votes = set()

    def consider(left, right):
        combined = left | right
        daughters = eligible_daughters.get(combined)
        if daughters and left & right and all(left & d and right & d for d in daughters):
            votes.add(combined)

    for node in tree.traverse("preorder"):
        children = node.children
        if not children:
            continue
        outside = above.pop(node)
        groups = tuple(observed[c] for c in children)
        if len(groups) == 2:
            a, b = groups
            if outside:
                consider(a, outside)
                consider(b, outside)
            above[children[0]], above[children[1]] = outside | b, outside | a
        else:
            if not outside and len(groups) == 3:
                a, b, c = groups
                consider(a, b)
                consider(a, c)
                consider(b, c)
            prefix = outside
            for child, group in zip(children, groups, strict=True):
                above[child] = prefix
                prefix |= group
            suffix = 0
            for child, group in zip(reversed(children), reversed(groups), strict=True):
                above[child] |= suffix
                suffix |= group
    return votes


def scan_duplications(species_tree, archive_path, minimum_score, observe=None, family_members=None,
                      observe_family=None, observe_duplications=None, *, holdout_output=None, coverage_output=None,
                      root_output=None, observe_topology=None):
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
    complete, partial_hits = Counter(), Counter()
    alternate, alternate_elsewhere = Counter(), Counter()
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
            if observe_topology is not None:
                observe_topology(family, tree, expected)
            projected, root_eligible = {}, {}
            for mask, (_, sisters, daughters, _) in branches.items():
                if full & sisters and all(full & daughter for daughter in daughters):
                    eligible[mask] += 1
                    if mask in supported:
                        hits[mask] += 1
                    if coverage_output is not None:
                        if full & mask == mask:
                            complete[mask] += 1
                        else:
                            projected[full & mask] = mask, daughters
                    if root_output is not None:
                        root_eligible[mask] = daughters
            if root_eligible:
                new_votes = alternate_root_votes(tree, observed, root_eligible) - supported
                alternate.update(new_votes)
                if supported & root_eligible.keys():
                    alternate_elsewhere.update(new_votes)
            if projected:
                # Each eligible split still has observed daughters and a sister,
                # so projected clades cannot coincide with another eligible split.
                seen = set()
                for gene, union in observed.items():
                    if union not in projected or len(gene.children) != 2 or union in seen:
                        continue
                    mask, daughters = projected[union]
                    left, right = (observed[c] for c in gene.children)
                    if left & right and all(left & d and right & d for d in daughters):
                        partial_hits[mask] += 1
                        seen.add(union)
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
        if coverage_output is not None:
            full, partial = complete[mask], n - complete[mask]
            projected = partial_hits[mask]
            values = (branch, branch_ids[node], full, partial, support, projected,
                      support / full if full else None, projected / partial if partial else None,
                      (support + projected) / n if n else None)
            coverage_output.append(dict(zip(COVERAGE_FIELDS, values, strict=True)))
        if root_output is not None:
            values = (branch, branch_ids[node], n, support, alternate[mask], alternate_elsewhere[mask],
                      (support + alternate[mask]) / n if n else None)
            root_output.append(dict(zip(ROOT_FIELDS, values, strict=True)))
    rows.sort(key=lambda r: (-r["duplication_score"], -r["num_supporting_families"], r["branch_id"]))
    if holdout is not None:
        holdout_output.extend(holdout.rows(rows))
    return rows, len(families)

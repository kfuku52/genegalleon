"""Check the projected audit against explicitly pruned gene trees."""

from collections import defaultdict

from workflow.support.wgd_taxon_holdout import DuplicationTaxonHoldout


def test_taxon_holdout_matches_explicit_tree_pruning():
    from nwkit.util import read_tree

    taxa = tuple("ABCDEF")
    bits = {name: 1 << i for i, name in enumerate(taxa)}

    def taxon_masks(tree):
        observed = {}
        for node in tree.traverse("postorder"):
            if node.is_leaf:
                observed[node] = bits[node.name.split("_", 1)[0]]
            else:
                mask = 0
                for child in node.children:
                    mask |= observed[child]
                observed[node] = mask
        return observed

    species = read_tree("(((A,B),(C,D)),(E,F));", "auto", True, quiet=True)
    sm = taxon_masks(species)
    branches = {sm[node]: (node, sm[node.up] ^ sm[node], tuple(sm[c] for c in node.children), "".join(node.leaf_names()))
                for node in species.descendants() if not node.is_leaf}
    sources = (
        "(((A_x,C_x),(B_y,D_y)),(E_o,F_o));",  # No overlap between copies.
        "(((A_x,B_x,C_x,E_x),(A_y,D_y)),F_o);",  # Removing E exposes ABCD.
        "(((A_x,B_x,C_x),(A_y,B_y,D_y),E_x),F_o);",  # Removing E makes a multifurcation binary.
        "(((A_x,B_x,C_x,D_x),(A_y,C_y)),(E_o,F_o));",
        "(((A_x,B_x,C_x),(A_y,B_y,D_y)),E_o);",  # Removing E destroys family eligibility.
        "(((A_x,C_x),(A_y,D_y)),F_o);",  # Removing A destroys this vote and eligibility.
        "(((A_x,B_x,C_x),(A_y,D_y)),E_o);",  # Removing A loses overlap but not eligibility.
        "(((A_x,C_x),(A_y,C_y)),(E_o,F_o));",  # Removing B or D tests changes in species coverage.
        "(((A_x,C_x),(B_y,D_y)),(A_o,E_o));",
        "(((A_x,A_y),B_o),(C_o,F_o));",  # Duplication restricted to A.
        "((((A_x,B_x),(A_y,B_y)),((A_z,B_z),(A_w,B_w))),(C_o,E_o));",  # Only one vote per family.
    )

    def score(tree, branch, drop=0):
        source, sister, daughters, _ = branch
        retained = [leaf.name for leaf in tree.leaves() if bits[leaf.name.split("_", 1)[0]] != drop]
        if len(retained) < 3:
            return False, False
        tree = tree.copy()
        tree.prune(retained)
        observed = taxon_masks(tree)
        observed_groups = tuple(group & observed[tree] for group in (sister, *daughters))
        if not all(observed_groups):
            return False, False
        for node, node_mask in observed.items():
            if node_mask != sm[source] & ~drop or len(node.children) != 2:
                continue
            left, right = (observed[child] for child in node.children)
            if left & right and all(left & daughter and right & daughter for daughter in daughters):
                return True, True
        return True, False

    baseline = defaultdict(lambda: [0, 0])
    expected = defaultdict(lambda: [0, 0, 0, 0])
    audit = DuplicationTaxonHoldout(taxa, branches)
    for source in sources:
        tree = read_tree(source, "auto", True, quiet=True)
        votes = set()
        for mask, branch in branches.items():
            eligible, supported = score(tree, branch)
            baseline[mask][0] += eligible
            baseline[mask][1] += supported
            if supported:
                votes.add(mask)
            for name, bit in bits.items():
                included, vote = score(tree, branch, bit)
                for i, value in enumerate((included, vote, vote and not supported, included and supported and not vote)):
                    expected[branch[3], name][i] += value
        audit.observe(tree, taxon_masks(tree), votes)
    cases = []
    for mask, (_, _, _, event) in branches.items():
        eligible, support = baseline[mask]
        cases.append({"species_event_id": event, "branch_id": mask, "num_eligible_families": eligible,
                      "num_supporting_families": support, "duplication_score": support / eligible if eligible else 0})
    rows = {(row["species_event_id"], row["excluded_taxon"]): row for row in audit.rows(cases)}
    for key, row in rows.items():
        assert [row[field] for field in ("eligible_families", "supporting_families", "gained_votes",
                                         "lost_votes_among_eligible")] == expected[key]
    assert rows["ABCD", "E"]["gained_votes"] == 2
    assert rows["ABCD", "A"]["lost_votes_among_eligible"] > 0
    assert rows["AB", "A"]["branch_identifiable"] == 0
    assert rows["AB", "A"]["duplication_score"] == "NA"
    assert rows["AB", "F"]["branch_identifiable"] == 1


def test_reused_coverage_does_not_share_votes_between_gene_trees():
    from nwkit.util import read_tree

    audit = DuplicationTaxonHoldout(("A", "B", "C", "D"), {3: (None, 12, (1, 2), "AB")})
    bits = dict(zip("ABCD", (1, 2, 4, 8), strict=True))

    def observe(newick, votes):
        tree = read_tree(newick, "auto", True, quiet=True)
        counts = {}
        for node in tree.traverse("postorder"):
            counts[node] = bits[node.name[0]] if node.is_leaf else counts[node.children[0]] | counts[node.children[1]]
        audit.observe(tree, counts, votes)

    observe("(((A_x,B_x),(A_y,B_y)),(C_x,D_x));", {3})
    assert audit._eligibility.cache_info().hits == 0
    observe("((A_x,B_x),(C_x,D_x));", set())
    assert audit._eligibility.cache_info().hits == 1
    rows = list(audit.rows([dict(species_event_id="AB", branch_id=1, num_eligible_families=2,
                                num_supporting_families=1, duplication_score=.5)]))
    assert [(r["excluded_taxon"], r["eligible_families"], r["supporting_families"]) for r in rows] == [
        ("A", 0, 0), ("B", 0, 0), ("C", 2, 1), ("D", 2, 1),
    ]

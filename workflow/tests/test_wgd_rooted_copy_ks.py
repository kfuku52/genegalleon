from collections import defaultdict

import numpy as np
import pytest

from workflow.support.wgd_rooted_copy_ks import EXPECTED, RootedCopyKs, designs, rooted_copy_summary, rooted_fit


@pytest.mark.parametrize("index", range(15))
def test_every_exact_five_tip_topology_has_its_own_optimum(index):
    models = designs()
    key = list(models)[index]
    expected, alternative, gap, edge, status = rooted_fit(models[key] @ np.arange(1, 8))
    assert status == ("expected_best" if key == EXPECTED else "other_best")
    if key == EXPECTED:
        assert expected == pytest.approx(0, abs=1e-12)
        assert edge == pytest.approx(6)
        assert gap > 0
    else:
        assert alternative == pytest.approx(0, abs=1e-12)
        assert gap < 0


def test_polytomy_is_not_resolved_from_roundoff():
    distances = designs()[EXPECTED] @ [0.2, 0.5, 0.4, 0.7, 0.1, 0, 0]
    assert rooted_fit(distances)[-1] == "tied"


def test_fifth_copy_is_an_adjacent_single_species_copy_with_same_code():
    from nwkit.clade_index import CladeIndex
    from nwkit.util import read_tree

    tree = read_tree("(((A:1,B:1):1,(C:0.5,D:0.7):1):1,E:3);", "auto", True)
    event = CladeIndex(tree).clade_id_for_node(tree.common_ancestor(["A", "B"]))
    genes = ("A_a", "B_a", "A_b", "B_b")
    hit = (event, "A", "B", ("chr", 1), "f", genes, (), ())
    gene_tree = read_tree("((((A_a,B_a),(A_b,B_b)),(C_a,D_a)),E_a);", "auto", True)
    members = {g: g[0] for g in gene_tree.leaf_names()}
    checker = RootedCopyKs(tree, dict(A=1, B=1, C=1, D=1, E=1))
    checker.observe("f", gene_tree, members, {hit})
    assert checker.outgroups[event, "f", genes] == ("C", "C_a")
    fields = (event, 9, "A", "B", "chr", 1, "f", *genes)
    targets = defaultdict(set)
    assert list(checker.selection([(fields, ())], targets)) == [(fields, ("C", "C_a"))]
    assert targets == {"A": {"A_a", "A_b"}, "B": {"B_a", "B_b"}, "C": {"C_a"}}
    checker = RootedCopyKs(tree, dict(A=1, B=1, C=1, D=2, E=1))
    checker.observe("f", gene_tree, {**members, "C_not_in_archive": "C"}, {hit})
    assert checker.outgroups == {}


def test_rooted_summary_separates_quartets_from_families():
    base = dict(species_event_id="key", branch_id=1, species_a="A", species_b="B")
    row, = rooted_copy_summary([{**base, "family_id": family, "status": status}
                                 for family, status in [("f", "expected_best"), ("f", "other_best"),
                                                        ("f", "tied"), ("g", "cds_unavailable")]])
    assert (row["sampled_quartets"], row["complete_quartets"], row["complete_families"], row["mixed_families"]) == (4, 3, 1, 1)

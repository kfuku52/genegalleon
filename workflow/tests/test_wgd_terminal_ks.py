from collections import defaultdict

import pytest

from workflow.support.wgd_terminal_ks import TerminalKs, sample_anchors


def test_sample_anchors_is_block_balanced_deterministic_and_family_unique():
    blocks = ((("a", "b"), {1, 2}), (("b", "a"), {1, 2}), (("c", "d"), {3}),
              (("e", "f"), {4}), (("g", "h"), {4}), (("i", "j"), {5}))
    positions = dict.fromkeys("abcdefghij", ("ctg", 1))
    membership = dict(zip("abcdefghi", ("x", "x", "x", "x", "u", "u", "v", "v", "w"), strict=True))
    num, chosen = sample_anchors("Sp", positions, dict(blocks), membership, 3, 71)
    assert num == 4
    assert len(set(a.family for a in chosen)) == 3
    assert len(set(a.block for a in chosen)) == 2
    assert chosen == sample_anchors("Sp", positions, dict(reversed(blocks)), membership, 3, 71)[1]


def test_terminal_ks_joins_gene_topology_to_complete_d_s_quartets():
    from nwkit.util import read_tree

    species = read_tree("(((A:1,B:1):1,C:2):1,D:3);", "auto", True)
    complete = [(f"a{k}", f"b{k}") for k in range(9)]
    sources = {sp: {"genetic_code": 2 if sp == "D" else 1} for sp in "ABCD"}
    members = [{"family": str(i), "A": f"{a},{b}", "B": "s", "C": "o", "D": "d"}
               for i, (a, b) in enumerate(complete)]
    links = {pair: {i % 6 + 1} for i, pair in enumerate(complete)}
    positions = {"A": {gene: ("ctg" + gene[0], i) for i, pair in enumerate(complete) for gene in pair}}
    terminal = TerminalKs(species, members, "family", sources, (positions, {"A": (links, {})}), 9, 41)
    results = {}
    for i, (a, b) in enumerate(complete):
        if i == 8:  # This quartet has an exterior copy only on an incompatible code.
            gene = read_tree(f"((({a},{b}),s),d);", "auto", True)
        else:
            gene = read_tree(f"((({a},{b}),s),o);", "auto", True)
        terminal.observe(str(i), gene, {a: "A", b: "A", "s": "B", "o": "C", "d": "D"})
        if i < 8:
            pairs = terminal.records[next(x for x in terminal.sampled["A"] if x.family == str(i))][2]
            for j, (pair, value) in enumerate(zip(pairs, (.3, .4, .8, .4, .8, .8), strict=True)):
                results[pair] = ("saturated", "NA") if i == 7 and j == 2 else ("ok", value)
    pairs, summaries = terminal.summarize(results)
    by_taxon = {r["species"]: r for r in summaries}
    observed = defaultdict(list)
    for row in pairs:
        observed[row["quartet_status"]].append(row)
    assert len(observed["terminal_minimum"]) == 7
    assert observed["dS_incomplete"][0]["outgroup_a_status"] == "saturated"
    assert observed["not_requested"][0]["tree_status"] == "no_same_code_rooted_outgroup"
    a = by_taxon["A"]
    assert a["topology_compatible_families"] == 8 and a["six_pair_families"] == 7
    assert a["original_blocks"] == 6
    assert a["median_signed_margin"] == pytest.approx(.1)
    assert a["median_copy_ds"] == a["copy_lower"] == a["copy_upper"] == .3
    assert by_taxon["B"]["status"] == "unavailable_self_or_cds"

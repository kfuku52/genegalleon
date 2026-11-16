import pytest

from workflow.support.wgd_ks_outgroups import (
    CommonOutgroups,
    common_outgroup_requests,
    common_outgroup_rows,
    select_common_outgroups,
    select_sister_panels,
    sibling_panel_requests,
    sibling_panel_rows,
    topology_outgroup_rows,
)


def test_common_outgroups_select_before_distances_and_keep_the_same_families():
    from nwkit.clade_index import CladeIndex
    from nwkit.ksrate_model import KsTrio
    from nwkit.util import read_tree

    tree = read_tree("((((F,Y),X),O),(P,Q));", "auto", True, quiet=True)
    index = CladeIndex(tree)
    nodes = {index.clade_id_for_node(node): node for node in tree.traverse()}
    event = index.clade_id_for_node(tree.common_ancestor(["F", "Y"]))
    older = index.clade_id_for_node(nodes[event].up)
    trios = [KsTrio("F", "Y", "X", event, 2), KsTrio("F", "X", "O", older, 1)]
    taxa = "FYXOPQ"
    archive = {str(i): {name: f"{name}_{i}" for name in taxa if (name, i) != ("P", 0)} for i in range(20)}
    codes = {sp: 1 if sp != "Q" else 2 for sp in taxa}
    row, = select_common_outgroups([(event, "F")], trios, nodes, index, archive, codes, 7, 31)
    assert row.outgroups == ("O", "P")
    assert row.eligible == 19 and len(row.families) == 7 and "0" not in row.families
    requests = list(common_outgroup_requests([row], archive))
    assert len(dict(requests)) == 8 * len(row.families)
    assert all("Q" not in details[:2] for _, details in requests)
    none, = select_common_outgroups([(event, "F")], trios[1:], nodes, index, archive, codes, 7, 31)
    assert list(common_outgroup_rows([none], {}, {event: 2}))[0]["status"] == "missing_trio"


def test_common_outgroup_contrasts_keep_negative_values_and_exclude_incomplete_families():
    families = tuple(str(i) for i in range(30))
    row = CommonOutgroups("event", "F", "Y", "X", ("O", "P"), len(families), families)
    sampled = {}
    for outside, older_distance in (("O", 5.0), ("P", 8.0)):
        for key, value in zip(row.pair_keys(outside), (2.0, 3.0, 6.0, older_distance, 4.0), strict=True):
            sampled[key] = dict.fromkeys(families, value)
    sampled["P", "X"].pop("0")
    a, b = common_outgroup_rows([row], sampled, {"event": 2})
    assert a["num_measured_all_outgroups"] == b["num_measured_all_outgroups"] == 29
    assert a["median_signed_gap"] == pytest.approx(2.0) and a["status"] == "ordered"
    assert b["median_signed_gap"] == pytest.approx(-1.0) and b["status"] == "reversed"
    assert all(r["gap_lower"] == r["median_signed_gap"] == r["gap_upper"] for r in (a, b))
    sampled["P", "X"] = {"1": 7.0}
    for diagnostic in common_outgroup_rows([row], sampled, {"event": 2}):
        assert diagnostic["num_measured_all_outgroups"] == 1
        assert diagnostic["gap_lower"] is diagnostic["gap_upper"] is None
        assert diagnostic["status"] == "interval_unavailable"
    sampled["P", "X"] = {}
    assert {r["status"] for r in common_outgroup_rows([row], sampled, {"event": 2})} == {"no_complete_families"}


def test_sister_panel_keeps_every_comparison_on_the_same_measured_families():
    from nwkit.clade_index import CladeIndex
    from nwkit.util import read_tree

    tree = read_tree("(((F,(Y,Z)),(X,W)),(O,P));", "auto", True, quiet=True)
    event_node = tree.common_ancestor(["F", "Y", "Z"])
    event = CladeIndex(tree).clade_id_for_node(event_node)
    primary = CommonOutgroups(event, "F", "Y", "X", ("O", "P"), 30, ())
    archive = {str(i): {name: f"{name}_{i}" for name in "FYZXWOP" if (name, i) != ("W", 0)} for i in range(30)}
    codes = dict.fromkeys("FYZXWOP", 1)
    panel, = select_sister_panels([primary], {event: event_node}, archive, codes, 40, 9)
    assert set(panel.choices) == {("Y", "X"), ("Y", "W"), ("Z", "X"), ("Z", "W")}
    assert len(panel.families) == panel.eligible == 29 and "0" not in panel.families
    measured = {}
    for _, (a, b, family, _) in sibling_panel_requests([panel], archive):
        measured.setdefault((a, b), {})[family] = 1.0
    measured["P", "Z"].pop("1")
    rows = list(sibling_panel_rows([panel], measured, {event: 6}))
    assert len(rows) == 8
    assert {row["num_measured_all_outgroups"] for row in rows} == {28}
    assert {row["num_sister_combinations"] for row in rows} == {4}
    assert {row["status"] for row in rows} == {"order_unresolved"}


def test_rooted_filter_uses_only_original_complete_families_and_all_outgroups(tmp_path):
    from zipfile import ZipFile

    files = tmp_path / "rooted.zip"
    with ZipFile(files, "w") as archive:
        for family, newick in {
            "good": "((((F,Y),X),O),P);",
            "reverse_sisters": "((((F,X),Y),O),P);",
            "unresolved_sisters": "(((F,Y,X),O),P);",
            "wrong_one_outgroup": "((((F,P),Y),X),O);",
            "unresolved": "[&R]((F,Y,X,O),P);",
            "not_measured": "((((F,Y),X),O),P);",
        }.items():
            archive.writestr(f"{family}_root.nwk", newick)
    families = ("good", "reverse_sisters", "unresolved_sisters", "wrong_one_outgroup", "unresolved", "not_measured")
    archived = {family: {name: name for name in "FYXOP"} for family in families}
    selected = CommonOutgroups("event", "F", "Y", "X", ("O", "P"), len(families), families)
    sampled = {}
    for out in selected.outgroups:
        for key in selected.pair_keys(out):
            sampled[key] = dict.fromkeys(families, 1.)
    sampled[tuple(sorted(("X", "P")))].pop("not_measured")
    ordered = []
    rows = list(topology_outgroup_rows(iter([selected]), sampled, {"event": 2}, archived, files, sister_order=ordered))
    assert {r["num_measured_before_topology_filter"] for r in rows + ordered} == {5}
    assert {r["num_measured_all_outgroups"] for r in rows} == {3}
    assert {r["num_measured_all_outgroups"] for r in ordered} == {1}
    assert {r["median_signed_gap"] for r in rows} == {0.}

import pytest

from workflow.support.wgd_synteny_plan import closest_untried, main, plan_pairs


def species_tree(newick):
    from nwkit.util import read_tree

    return read_tree(newick, "auto", True, rooted="auto", quiet=True)


def test_nearest_untried_matches_complete_pair_order():
    first = [(0, "a"), (0, "d"), (1, "e")]
    second = [(0, "b"), (1, "c"), (1, "f"), (2, "g")]
    expected = sorted((d + e, tuple(sorted((x, y)))) for d, x in first for e, y in second)
    used = set()
    for item in expected:
        assert closest_untried(first, second, used) == closest_untried(second, first, used) == item
        used.add(item[1])
    assert closest_untried(first, second, used) is None


def test_plan_prioritizes_nearby_divisions_and_retries_only_missing_coverage():
    tree = species_tree("(((A:1,B:2):1,(C:1,D:3):1):1,E:4);")
    self_pairs = {(name, name) for name in "ABCDE"}
    first = ("A", "C")
    pairs, audit = plan_pairs(tree, self_pairs)
    assert [(r["target_species"], r["query_species"]) for r in pairs] == [("A", "B"), ("C", "D"), first]
    known = self_pairs | {first, ("A", "B"), ("C", "D"), ("A", "E")}
    assert not plan_pairs(tree, known)[0]
    retried, audit = plan_pairs(tree, known, {("C", "D")})
    assert [(r["target_species"], r["query_species"]) for r in retried] == [("B", "C")]
    assert {r["descendant_taxa"]: r["status"] for r in audit} == {
        "A,B": "unassessed_pairs_exhausted", "C,D": "covered_by_assessed_pairs", "A,B,C,D": "comparisons_planned",
    }
    missing = known - {("A", "A"), ("B", "B")}
    assert [r["status"] for r in plan_pairs(tree, missing, set())[1]].count("missing_self_synteny") == 2


def test_multifurcating_plan_covers_every_available_daughter_without_all_pairs():
    tree = species_tree("((A:1,B:2,C:3,D:4,F:5):1,E:6);")
    pairs, audit = plan_pairs(tree, {(name, name) for name in "ABCDE"})
    assert [(r["target_species"], r["query_species"]) for r in pairs] == [("A", "B"), ("C", "D")]
    row, = audit
    assert row["num_daughters_after_planned_pairs"] == 4
    assert row["num_daughters"] == 5
    assert row["status"] == "missing_self_synteny"


def test_sparse_replication_planning_preserves_a_focal_and_bounds_extra_pairs():
    tree = species_tree("(((A:1,B:2):1,(C:1,D:3):1):1,E:4);")
    seeds = {(name, name) for name in "ABCDE"}
    assessed = {("A", "C"), ("A", "B"), ("C", "D")}
    with pytest.raises(ValueError, match="assessed pair evidence"):
        plan_pairs(tree, seeds, replication_budget=1)
    known = seeds | assessed
    rows, audit = plan_pairs(tree, known, assessed, replication_budget=1)
    assert [(row["target_species"], row["query_species"]) for row in rows] == [("B", "C")]
    assert audit[-1]["status"] == "replication_comparison_planned"
    assert not plan_pairs(tree, known | {("B", "C")}, assessed | {("B", "C")}, replication_budget=1)[0]
    rows, audit = plan_pairs(tree, known | {("B", "C")}, assessed, replication_budget=1)
    assert [(row["target_species"], row["query_species"]) for row in rows] == [("A", "D")]


def test_planner_cli_writes_pairwise_schema_and_protects_inputs(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("((A:1,B:1):1,C:2);\n")
    manifest = tmp_path / "cross.tsv"
    manifest.write_text("species_a\tspecies_b\tbed_a\tbed_b\tanchors\n"
                        "A\tA\tA.bed\tA.bed\tA.anchors\nB\tB\tB.bed\tB.bed\tB.anchors\n")
    for name in ("A.bed", "B.bed", "A.anchors", "B.anchors"):
        (tmp_path / name).touch()
    pairs, audit = tmp_path / "next.tsv", tmp_path / "audit.tsv"
    args = ["--workspace", str(tmp_path), "--species-tree", str(tree), "--cross-synteny", str(manifest),
            "--audit-out", str(audit), "--pairs-out"]
    with pytest.raises(ValueError, match="must not replace"):
        main([*args, str(manifest)])
    main([*args, str(pairs)])
    header, row = pairs.read_text().splitlines()
    assert header == "analysis_id\ttarget_species\tquery_species"
    assert row.startswith("wgd-") and row.endswith("\tA\tB")
    assert "comparisons_planned" in audit.read_text()
    with pytest.raises(ValueError, match="Refusing to replace"):
        main([*args, str(pairs)])


def test_planner_accepts_archive_limited_pair_evidence(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("((A:1,B:1):1,C:2);\n")
    manifest = tmp_path / "cross.tsv"
    manifest.write_text("species_a\tspecies_b\tbed_a\tbed_b\tanchors\n"
                        "A\tA\tA.bed\tA.bed\tA.anchors\nB\tB\tB.bed\tB.bed\tB.anchors\n"
                        "A\tB\tA.bed\tB.bed\tAB.anchors\n")
    for name in ("A.bed", "B.bed", "A.anchors", "B.anchors", "AB.anchors"):
        (tmp_path / name).touch()
    evidence = tmp_path / "evidence.tsv"
    evidence.write_text("species_a\tspecies_b\tstatus\nA\tB\tinsufficient_archive_anchors\n")
    pairs, audit = tmp_path / "next.tsv", tmp_path / "audit.tsv"
    main(["--workspace", str(tmp_path), "--species-tree", str(tree), "--cross-synteny", str(manifest),
          "--pair-evidence", str(evidence), "--audit-out", str(audit), "--pairs-out", str(pairs)])
    assert len(pairs.read_text().splitlines()) == 1
    assert "unassessed_pairs_exhausted" in audit.read_text()

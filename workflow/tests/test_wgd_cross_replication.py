from workflow.support.wgd_cross_replication import replicated_copy_rows


def test_replication_counts_exact_focal_copies_not_just_family_or_windows():
    candidates = [{"species_event_id": "branch", "branch_id": 6}]
    pairs = [{"species_event_id": "branch", "species_a": a, "species_b": b, "status": status}
             for a, b, status in (("A", "B", "assessed"), ("C", "A", "assessed"), ("A", "D", "assessed"),
                                  ("A", "E", "insufficient_archive_anchors"))]

    def hit(a, b, window, family, genes):
        return "branch", a, b, ("chr", window), family, genes, (), ()

    h1 = hit("A", "B", 1, "f1", ("a", "b", "a2", "b2"))
    h2 = hit("C", "A", 3, "f1", ("c", "a2", "c2", "a"))
    hits = [h1, h2, hit("A", "B", 2, "f1", h1[5]),
            hit("A", "B", 1, "f2", ("a3", "b3", "a4", "b4")),
            hit("C", "A", 1, "f2", ("c3", "a4", "c4", "a5")),
            hit("A", "D", 1, "f1", ("a", "d", "a2")),
            hit("A", "E", 1, "f1", ("a", "e", "a2", "e2"))]
    rows = list(replicated_copy_rows(candidates, pairs, hits, [h1, h2]))
    assert len(rows) == 6
    assert {row["focal"] for row in rows} == {"A"}
    assert {row["comparator_b"] for row in rows} == {"C", "D"}
    baseline, strict = (row for row in rows if row["comparator_a"] == "B" and row["comparator_b"] == "C")
    assert (baseline["copy_pairs_a"], baseline["copy_pairs_b"], baseline["shared_copy_pairs"], baseline["shared_families"]) == (2, 2, 1, 1)
    assert strict["copy_pairs_a"] == strict["copy_pairs_b"] == strict["shared_copy_pairs"] == 1
    assert all(row["copy_pairs_b"] == row["shared_copy_pairs"] == 0 for row in rows if row["comparator_b"] == "D")

import pytest

from workflow.support.wgd_evidence import combine_node, positional_feature, summarize_events


def positions():
    return {"a": {"species": "A", "seqid": "chr1", "rank": 1, "locus_id": "l1", "start": 0, "end": 80},
            "b": {"species": "A", "seqid": "chr1", "rank": 2, "locus_id": "l2", "start": 100, "end": 180},
            "c": {"species": "A", "seqid": "chr2", "rank": 1, "locus_id": "l3", "start": 0, "end": 80}}


def test_missing_synteny_never_means_ssd():
    assert combine_node("event", [("a", "c")], positions(), {}, {})[0] == "unresolved"
    assert combine_node("event", [("a", "missing")], positions(), {}, {})[0] == "unresolved"


def test_terminal_tandem_is_positive_ssd_evidence_but_conflicts_stay_unresolved():
    assert combine_node("event", [("a", "b")], positions(), {}, {})[:2] == ("SSD-supported", "terminal_tandem_adjacency")
    assert combine_node("event", [("a", "b")], positions(), {}, {}, terminal_cherry=False)[0] == "unresolved"
    assert combine_node("event", [("a", "b")], positions(), {("a", "b"): [{"species_event_id": "older"}]}, {})[0] == "unresolved"


def test_annotation_isoforms_are_not_duplicate_origins():
    data = positions()
    data["b"]["locus_id"] = "l1"
    assert positional_feature("a", "b", data) == ("same_locus", 0)
    assert combine_node("event", [("a", "b")], data, {}, {})[1] == "same_locus_annotation_ambiguity"


def test_wgd_requires_same_branch_not_just_collinearity():
    anchors = {("a", "c"): [{"species_event_id": "older", "placement_status": "interval_supported",
                               "ks_status": "ok", "ks": "0.5", "block_id": "block1", "gene_a": "a", "gene_b": "c", "species": "A"}]}
    events = {"older": {"event_support": "WGD-supported"}}
    assert combine_node("younger", [("a", "c")], positions(), anchors, events)[0] == "unresolved"
    assert combine_node("older", [("a", "c")], positions(), anchors, events)[0] == "WGD-supported"


def test_local_segmental_blocks_uncalibrated_counts_and_missing_age_do_not_support_wgd():
    candidate = {"species_event_id": "event", "descendant_taxa": "A,B", "count_support": "count_supported_conditional",
                 "nuisance_bound_reached": "False", "background_nuisance_bound_reached": "False",
                 "branch_burst_nuisance_bound_reached": "False"}
    anchors = [{"species_event_id": "event", "species": name, "block_id": str(i),
                "gene_a": f"{name}{i}a", "gene_b": f"{name}{i}b", "placement_status": "interval_supported",
                "ks_status": "ok", "ks": "0.5"}
               for name in ("A", "B") for i in range(3)]
    genomic = {row[arm]: {"species": row["species"], "seqid": arm, "rank": i + 1, "locus_id": row[arm],
                         "start": i * 100, "end": i * 100 + 80}
               for i, row in enumerate(anchors) for arm in ("gene_a", "gene_b")}
    summaries = {name: {"num_annotated_loci": 100, "gene_span_coverage": 1} for name in ("A", "B")}
    assert summarize_events([candidate], anchors, summaries, positions=genomic)[0]["event_support"] == "unresolved"
    summaries = {name: {"num_annotated_loci": 20} for name in ("A", "B")}
    assert summarize_events([candidate], anchors, summaries, positions=genomic)[0]["event_support"] == "WGD-supported"
    assert summarize_events([{**candidate, "count_support": "not_calibrated"}], anchors, summaries, positions=genomic)[0]["event_support"] == "unresolved"
    assert summarize_events([{**candidate, "nuisance_bound_reached": "True"}], anchors, summaries, positions=genomic)[0]["event_support"] == "unresolved"
    for row in anchors:
        row["placement_status"] = "boundary_overlap"
    assert summarize_events([candidate], anchors, summaries, positions=genomic)[0]["event_support"] == "unresolved"


@pytest.mark.parametrize("coordinates", [(0, 150, 100, 200), (0, 200, 50, 100), (0, 100, 0, 100)])
def test_overlapping_or_nested_loci_are_not_terminal_tandem_ssd(coordinates):
    data = positions()
    data["a"].update(start=coordinates[0], end=coordinates[1])
    data["b"].update(start=coordinates[2], end=coordinates[3])
    assert positional_feature("a", "b", data)[0] == "overlapping_loci"
    assert combine_node("event", [("a", "b")], data, {}, {})[0] == "unresolved"


def test_internal_tandem_collinearity_conflict_is_not_wgd():
    row = {"species_event_id": "event", "placement_status": "interval_supported", "ks_status": "ok",
           "ks": "0.5", "gene_a": "a", "gene_b": "b", "species": "A"}
    result = combine_node("event", [("a", "b")], positions(), {("a", "b"): [row]},
                          {"event": {"event_support": "WGD-supported"}}, terminal_cherry=False)
    assert result[:2] == ("unresolved", "tandem_and_collinearity_conflict")


def test_missing_coordinates_cannot_establish_tandem_adjacency():
    data = positions()
    data["a"].pop("start")
    assert combine_node("event", [("a", "b")], data, {}, {})[0] == "unresolved"


@pytest.mark.parametrize("status,ks", [("saturated", "0.5"), ("ok", "NA"), ("ok", "-1")])
def test_ks_status_and_finite_value_are_both_required_for_node_support(status, ks):
    row = {"species_event_id": "event", "placement_status": "interval_supported", "ks_status": status,
           "ks": ks, "gene_a": "a", "gene_b": "c", "species": "A"}
    assert combine_node("event", [("a", "c")], positions(), {("a", "c"): [row]},
                        {"event": {"event_support": "WGD-supported"}})[0] == "unresolved"

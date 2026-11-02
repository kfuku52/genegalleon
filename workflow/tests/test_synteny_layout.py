import itertools
import json
import math
import random
from types import SimpleNamespace

import numpy as np
import pytest

from workflow.support import synteny_layout as layout
from workflow.support.synteny_multitrack import load_panels


def genes(prefix, seqids, count=8):
    return [SimpleNamespace(gene_id=f"{prefix}{sid}_{i}", seqid=sid, start=i * 10, end=i * 10 + 5)
            for sid in seqids for i in range(count)]


def block_file(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def test_signed_two_track_optimum_and_original_identifiers_are_preserved(tmp_path):
    genomes = [genes("a", ["chr1"]), genes("b", ["chr1"])]
    tracks = [{"species": "A", "fixed_orientation": {"chr1": "+"}}, {"species": "B"}]
    simple = block_file(tmp_path, "ab.simple", "achr1_0 achr1_7 bchr1_0 bchr1_7 8 -\n")
    original = json.dumps([[vars(g) for g in row] for row in genomes])
    orders, directions, report = layout.optimize_layout(tracks, genomes, [{"tracks": [0, 1], "blocks": simple}])
    assert orders == [["chr1"], ["chr1"]]
    assert directions == [{"chr1": False}, {"chr1": True}]
    assert report["globally_optimal"] and report["orientation_changed"]
    assert report["objective_after"]["twist"] == 0
    assert report["objective_before"]["twist"] == 1
    assert original == json.dumps([[vars(g) for g in row] for row in genomes])
    assert (tmp_path / "ab.simple").read_text().endswith("8 -\n")


def test_middle_track_is_one_joint_decision_using_both_edges(tmp_path):
    tracks = [{"species": "A", "fixed_orientation": {"chr1": "+"}}, {"species": "D"},
              {"species": "T", "fixed_orientation": {"chr1": "+"}}]
    genomes = [genes(p, ["chr1"]) for p in "adt"]
    first = block_file(tmp_path, "ad", "achr1_0 achr1_7 dchr1_0 dchr1_7 8 -\n")
    second = block_file(tmp_path, "dt", "dchr1_0 dchr1_1 tchr1_0 tchr1_1 2 +\n")
    edges = [{"tracks": [0, 1], "blocks": first}, {"tracks": [1, 2], "blocks": second}]
    _, directions, report = layout.optimize_layout(tracks, genomes, edges)
    assert directions == [{"chr1": False}, {"chr1": True}, {"chr1": False}]
    assert report["displayed_block_counts"] == [1, 1]
    assert report["species_order"] == ["A", "D", "T"]
    assert report["objective_after"]["total"] < report["objective_before"]["total"]
    # Ignoring the upper edge would choose the opposite middle-track direction.
    _, only_lower, _ = layout.optimize_layout(tracks[1:], genomes[1:], [{"tracks": [0, 1], "blocks": second}])
    assert only_lower[0] == {"chr1": False}


@pytest.mark.parametrize("seed", range(5))
def test_crossing_weight_matches_independent_quadratic_count_with_ties(seed):
    rng = random.Random(seed)
    x = np.array([rng.randrange(4) for _ in range(30)], dtype=float)
    y = np.array([rng.randrange(5) for _ in range(30)], dtype=float)
    w = np.array([rng.random() for _ in range(30)])
    expected = sum(w[i] * w[j] for i in range(len(x)) for j in range(i)
                   if (x[i] - x[j]) * (y[i] - y[j]) < 0)
    assert layout.crossing_weight(x, y, w) == pytest.approx(expected)


def test_joint_exact_solution_matches_independent_permutation_flip_enumeration(tmp_path):
    tracks = [{"species": "A", "seqids": ["1", "2"]}, {"species": "B", "seqids": ["1", "2"]}]
    genomes = [genes(p, ["1", "2"], 4) for p in "ab"]
    raw = [("1", 0, 3, "2", 0, 3, -1), ("2", 0, 2, "1", 1, 3, 1), ("1", 1, 2, "1", 0, 1, 1)]
    source = block_file(tmp_path, "ab", "".join(f"a{a}_{p} a{a}_{q} b{b}_{r} b{b}_{s} 4 {'+' if sign == 1 else '-'}\n"
                                                   for a, p, q, b, r, s, sign in raw))
    _, _, report = layout.optimize_layout(tracks, genomes, [{"tracks": [0, 1], "blocks": source}])
    def independent(order_a, order_b, flip_a, flip_b):
        ratio, gap = .79 / 8, .01
        positions, weights, lengths, twisted = [], [], [], []
        for a, p, q, b, r, s, sign in raw:
            if sign == -1:
                r, s = s, r
            weight = (abs(p - q) + abs(r - s)) * ratio / 2
            if flip_a[int(a) - 1]:
                p, q = 3 - p, 3 - q
            if flip_b[int(b) - 1]:
                r, s = 3 - r, 3 - s
            x = .12 + order_a.index(a) * (4 * ratio + gap) + (p + q) * ratio / 2
            y = .12 + order_b.index(b) * (4 * ratio + gap) + (r + s) * ratio / 2
            weights.append(weight)
            lengths.append(weight * math.hypot((x - y) * 20, .6 * 8))
            twisted.append(weight if (q - p) * (s - r) < 0 else 0)
            positions.append((x, y))
        crossing = sum(weights[i] * weights[j] for i in range(3) for j in range(i)
                       if (positions[i][0] - positions[j][0]) * (positions[i][1] - positions[j][1]) < 0)
        norm = sum(weights[i] * weights[j] for i in range(3) for j in range(i))
        return sum(lengths) / (sum(weights) * math.hypot(20, 8)) + sum(twisted) / sum(weights) + crossing / norm
    candidates = (independent(a, b, fa, fb) for a in itertools.permutations(["1", "2"])
                  for b in itertools.permutations(["1", "2"]) for fa in itertools.product((False, True), repeat=2)
                  for fb in itertools.product((False, True), repeat=2))
    assert report["objective_after"]["total"] == pytest.approx(min(candidates), abs=1e-12)
    assert report["globally_optimal"] is True


def test_bounded_search_preserves_groups_locks_and_seed(tmp_path, monkeypatch):
    monkeypatch.setattr(layout, "EXACT_STATES", 0)
    tracks = [{"species": "A", "seqids": ["a", "b", "c", "d"], "groups": [["a", "b"], ["c", "d"]],
               "fixed_orientation": {"a": "+"}}, {"species": "B", "fixed_order": True}]
    genomes = [genes("a", ["a", "b", "c", "d"]), genes("b", ["a", "b", "c", "d"])]
    source = block_file(tmp_path, "ab", "".join(f"a{x}_0 a{x}_7 b{y}_0 b{y}_7 8 -\n"
                                                   for x, y in zip("abcd", "dcba", strict=True)))
    edges = [{"tracks": [0, 1], "blocks": source}]
    one = layout.optimize_layout(tracks, genomes, edges, restarts=3, max_rounds=5)
    two = layout.optimize_layout(tracks, genomes, edges, restarts=3, max_rounds=5)
    assert one == two
    order, directions, report = one
    assert abs(order[0].index("a") - order[0].index("b")) == 1
    assert abs(order[0].index("c") - order[0].index("d")) == 1
    assert order[1] == list("abcd")
    assert directions[0]["a"] is False
    assert report["globally_optimal"] is False
    assert report["objective_after"]["total"] <= report["objective_before"]["total"] + 1e-12


@pytest.mark.parametrize("mode", ["none", "order", "order_orientation"])
def test_uninformative_zero_width_connections_keep_input_order_and_direction(tmp_path, mode):
    tracks = [{"species": "A", "seqids": list("dcba")}, {"species": "B", "seqids": list("badc")}]
    genomes = [genes(p, list("abcd"), 1) for p in "ab"]
    source = block_file(tmp_path, "ab", "".join(f"a{x}_0 a{x}_0 b{x}_0 b{x}_0 1 -\n" for x in "abcd"))
    order, direction, report = layout.optimize_layout(tracks, genomes, [{"tracks": [0, 1], "blocks": source}], mode=mode)
    assert order == [list("dcba"), list("badc")]
    assert not any(value for row in direction for value in row.values())
    assert report["objective_after"]["total"] == 0
    assert report["displayed_block_counts"] == [4]


@pytest.mark.parametrize("change", ["noncontiguous", "overlap", "bad-lock", "unknown", "bad-score", "cross-sequence"])
def test_invalid_layout_constraints_and_sources_fail(tmp_path, change):
    tracks = [{"species": "A", "seqids": list("abc")}, {"species": "B", "seqids": list("abc")}]
    genomes = [genes(p, list("abc")) for p in "ab"]
    text = "aa_0 aa_7 ba_0 ba_7 8 +\n"
    if change == "noncontiguous":
        tracks[0]["groups"] = [["a", "c"]]
    elif change == "overlap":
        tracks[0]["groups"] = [["a", "b"], ["b", "c"]]
    elif change == "bad-lock":
        tracks[0]["fixed_orientation"] = {"unknown": "+"}
    elif change == "unknown":
        tracks[0]["fixed_oder"] = True
    elif change == "bad-score":
        text = text.replace(" 8 +", " 0 +")
    else:
        text = text.replace("aa_7", "ab_7")
    with pytest.raises(ValueError):
        layout.optimize_layout(tracks, genomes, [{"tracks": [0, 1], "blocks": block_file(tmp_path, "ab", text)}])


def test_panel_schema_rejects_unknown_comparisons_species_and_duplicate_json(tmp_path):
    pair = {"analysis_id": "ad", "target": {"species": "A"}, "query": {"species": "D"}}
    valid = {"panels": [{"panel_id": "adt", "tracks": [{"species": "A"}, {"species": "D"}], "comparisons": ["ad"]}]}
    path = tmp_path / "panels.json"
    path.write_text(json.dumps(valid))
    assert load_panels(path, [pair]) == valid["panels"]
    for key, value in (("comparisons", ["missing"]), ("tracks", [{"species": "A"}, {"species": "T"}]), ("seed", False)):
        bad = json.loads(json.dumps(valid))
        bad["panels"][0][key] = value
        path.write_text(json.dumps(bad))
        with pytest.raises(ValueError):
            load_panels(path, [pair])
    path.write_text('{"panels":[],"panels":[]}')
    with pytest.raises(ValueError, match="Duplicate"):
        load_panels(path, [pair])

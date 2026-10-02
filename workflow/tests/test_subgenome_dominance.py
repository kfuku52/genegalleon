import csv
import itertools
import json

import numpy as np
import pytest

from workflow.support import subgenome_dominance as sd


def fixture_inputs(root, effect=1, scope="local"):
    root.mkdir(parents=True)
    mapping = ["gene_id\tgroup_id\tsubgenome\tassignment_scope\tassignment_basis\tevidence"]
    pairs = ["pair_id\tblock_id\tgene_id"]
    expr = ["gene_id\tleaf1\tleaf1tech\tleaf2"]
    retention = ["group_id\tblock_id\tlocus_id\tsubgenome\tcallable\tretained"]
    for i in range(12):
        for sg in ("A", "B"):
            gene = f"g{i}{sg}"
            mapping.append(f"{gene}\tgroup1\t{sg}\t{scope}\tsynteny\tvalidated_block")
            pairs.append(f"p{i}\tb{i}\t{gene}")
            value = 2 ** effect if sg == "A" else 1
            expr.append(f"{gene}\t{value}\t{value}\t{value}")
            retained = int(sg == "A") if effect else 1
            retention.append(f"group1\tb{i}\tl{i}\t{sg}\t1\t{retained}")
    for name, rows in [("mapping", mapping), ("pairs", pairs), ("expr", expr), ("retention", retention)]:
        (root / f"{name}.tsv").write_text("\n".join(rows) + "\n")
    (root / "samples.tsv").write_text("column\ttissue\tbiological_id\nleaf1\tleaf\tplant1\nleaf1tech\tleaf\tplant1\nleaf2\tleaf\tplant2\n")
    (root / "manifest.tsv").write_text("analysis_id\tspecies\tmapping_file\tretention_file\thomoeolog_file\texpression_file\tsamples_file\nfixture\tTest_species\tmapping.tsv\tretention.tsv\tpairs.tsv\texpr.tsv\tsamples.tsv\n")
    return sd.make_plan(root.parent, root / "manifest.tsv", 200, 7)


@pytest.mark.parametrize("effect", [0, 1, -1])
def test_known_expression_bias_and_null(tmp_path, effect):
    plan = fixture_inputs(tmp_path / "input", effect)
    result = sd.analyse(plan["analyses"][0], tmp_path / "results", 200, 7)
    expression = next(r for r in result["statistics"] if r["metric"] == "expression_log2_ratio")
    assert expression["effect"] == pytest.approx(effect)
    assert expression["n_loci"] == expression["n_blocks"] == 12
    assert expression["p_value"] == (1 if effect == 0 else pytest.approx(2 / 4096))
    assert result["genomewide_identity"] == "unresolved"
    assert result["genomewide_dominance"] == "not_tested"
    rows = list(csv.DictReader((tmp_path / "results/expression_pairs.tsv").open(), delimiter="\t"))
    assert all(r["biological_samples"] == "2" for r in rows)


def test_uncallable_retention_is_excluded_not_loss(tmp_path):
    plan = fixture_inputs(tmp_path / "input")
    path = tmp_path / "input/retention.tsv"
    path.write_text(path.read_text().replace("group1\tb0\tl0\tB\t1\t0", "group1\tb0\tl0\tB\t0\t0"))
    result = sd.analyse(plan["analyses"][0], tmp_path / "results", 200, 7)
    retention = result["statistics"][0]
    assert retention["n_opportunities"] == 12
    assert retention["n_loci"] == retention["retained_a"] == 11
    sd.verify_inputs(sd.make_plan(tmp_path, tmp_path / "input/manifest.tsv"))
    with pytest.raises(ValueError, match="Input changed"):
        sd.verify_inputs(plan)


def test_zeros_are_audited_without_pseudocount(tmp_path):
    plan = fixture_inputs(tmp_path / "input")
    path = tmp_path / "input/expr.tsv"
    path.write_text(path.read_text().replace("g0B\t1\t1\t1", "g0B\t1\t1\t0"))
    result = sd.analyse(plan["analyses"][0], tmp_path / "results", 200, 7)
    expression = result["statistics"][1]
    assert expression["n_opportunities"] == 12 and expression["n_loci"] == 11
    coverage = list(csv.DictReader((tmp_path / "results/expression_coverage.tsv").open(), delimiter="\t"))
    assert coverage[0]["positive_samples"] == "1"
    assert coverage[0]["positive_a_samples"] == "2" and coverage[0]["positive_b_samples"] == "1"
    detection = next(r for r in result["statistics"] if r["metric"] == "expression_detection_difference")
    assert detection["n_loci"] == 12
    assert detection["effect"] == pytest.approx(0.5 / 12)


@pytest.mark.parametrize("kind", ["duplicate", "expression_assignment", "reserved_group", "missing_expected_side", "negative", "unmatched_expression"])
def test_invalid_scientific_inputs_fail(tmp_path, kind):
    plan = fixture_inputs(tmp_path / "input")
    if kind in {"duplicate", "expression_assignment", "reserved_group"}:
        path = tmp_path / "input/mapping.tsv"
        content = path.read_text()
        if kind == "duplicate":
            path.write_text(content + content.splitlines()[1] + "\n")
        elif kind == "reserved_group":
            path.write_text(content.replace("group1", "ALL_GROUPS"))
        else:
            path.write_text(content.replace("\tsynteny\t", "\texpression\t"))
    elif kind == "missing_expected_side":
        path = tmp_path / "input/retention.tsv"
        path.write_text("\n".join(path.read_text().splitlines()[:-1]) + "\n")
    else:
        path = tmp_path / "input/expr.tsv"
        path.write_text(path.read_text().replace("g0B\t1\t1\t1", "g0B\t-1\t1\t1") if kind == "negative" else path.read_text().replace("\ng", "\nunmatched_g"))
    with pytest.raises(ValueError):
        sd.analyse(plan["analyses"][0], tmp_path / "results", 200, 7)


def test_block_dependence_and_mosaic_null():
    values = [1, -1] * 6
    result = sd.inference(values, [f"b{i}" for i in range(12)], 1000, np.random.default_rng(1))
    assert result["effect"] == 0 and result["p_value"] == 1
    result = sd.inference([1] * 100, ["same_block"] * 100, 1000, np.random.default_rng(1))
    assert result["n_blocks"] == 1 and result["p_value"] is None
    assert result["status"] == "insufficient_blocks"


def test_three_subgenomes_and_missing_assays(tmp_path):
    path = tmp_path / "mapping.tsv"
    path.write_text("gene_id\tgroup_id\tsubgenome\tassignment_scope\tassignment_basis\tevidence\n"
                    + "".join(f"g{sg}\tG\t{sg}\tlocal\tcurated\tindependent_phasing\n" for sg in "ABC"))
    analysis = {"mapping_file": str(path), "retention_file": "", "expression_file": "", "analysis_id": "x", "species": "test"}
    result = sd.analyse(analysis, tmp_path / "output", 200, 1)
    assert "not_estimable" in result["retention_status"]
    assert result["mapped_genes"] == 3 and not result["statistics"]
    assert json.loads((tmp_path / "output/summary.json").read_text())["genomewide_identity"] == "unresolved"


def test_independent_global_identity_enables_pooled_contrasts(tmp_path):
    plan = fixture_inputs(tmp_path / "input", scope="global")
    result = sd.analyse(plan["analyses"][0], tmp_path / "output", 200, 7)
    pooled = [r for r in result["statistics"] if r["group_id"] == "ALL_GROUPS"]
    assert len(pooled) == 3
    assert all(r["n_loci"] == 12 for r in pooled)
    assert [r["effect"] for r in pooled] == [1, 1, 0]
    assert result["genomewide_dominance"] == "contrasts_available"


def test_expression_units_are_explicit(tmp_path):
    plan = fixture_inputs(tmp_path / "input")
    manifest = tmp_path / "input/manifest.tsv"
    rows = table = list(csv.DictReader(manifest.open(), delimiter="\t"))
    sd.write_table(manifest, [{**table[0], "expression_unit": "counts"}], [*rows[0], "expression_unit"])
    with pytest.raises(ValueError, match="expression_unit"):
        sd.make_plan(tmp_path, manifest)
    sd.write_table(manifest, [{**table[0], "expression_unit": "FPKM"}], [*rows[0], "expression_unit"])
    plan = sd.make_plan(tmp_path, manifest)
    result = sd.analyse(plan["analyses"][0], tmp_path / "out", 200, 1)
    assert result["expression_unit"] == "FPKM" and result["statistics"][1]["effect"] == 1


def test_expression_stream_unchanged_by_retention_and_input_order(tmp_path):
    plan = fixture_inputs(tmp_path / "input")
    path = tmp_path / "input/expr.tsv"
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    for i, row in enumerate(rows):
        if row["gene_id"].endswith("A"):
            for column in ("leaf1", "leaf1tech", "leaf2"):
                row[column] = str(i + 2.5)
    sd.write_table(path, rows, rows[0].keys())
    analysis = plan["analyses"][0]
    first = sd.analyse(analysis, tmp_path / "first", 200, 7)
    # Other metric removal and file row ordering must leave this comparison invariant.
    sd.write_table(path, list(reversed(rows)), rows[0].keys())
    samples = tmp_path / "input/samples.tsv"
    lines = samples.read_text().splitlines()
    samples.write_text("\n".join([lines[0], *reversed(lines[1:])]) + "\n")
    second = sd.analyse({**analysis, "retention_file": ""}, tmp_path / "second", 200, 7)
    def expr(result):
        return [r for r in result["statistics"] if r["metric"].startswith("expression")]
    assert expr(first) == expr(second)
    assert sd.contrast_seed(7, "x", "metric", "G", "A", "B", "leaf", "bootstrap") != sd.contrast_seed(
        7, "x", "metric", "G", "A", "B", "leaf", "permutation")


@pytest.mark.parametrize("weights", [[0, 0, 0], [3, -3, 6, 0], [2, -1, 4, 2, -3, 1]])
def test_integer_dp_matches_independent_exhaustive_null(weights):
    observed = abs(sum(weights))
    expected = sum(abs(sum(s * w for s, w in zip(signs, weights, strict=True))) >= observed
                   for signs in itertools.product((-1, 1), repeat=len(weights))) / 2 ** len(weights)
    result = sd.sign_flip(np.array(weights, dtype=float), np.random.default_rng(1), 100, 100)
    assert result["test_method"] == "exact_integer_dp"
    assert result["p_value"] == expected
    assert result["null_draws"] == 2 ** len(weights)


def test_more_than_sixteen_integer_blocks_are_exact():
    result = sd.inference([1] * 30, [f"b{i}" for i in range(30)], 100, np.random.default_rng(1))
    assert result["test_method"] == "exact_integer_dp"
    assert result["p_value"] == 2 / 2 ** 30


def test_noninteger_meet_in_middle_matches_exhaustive_null():
    weights = np.array([.13, -.72, .49, .57, -.11, .38, .61, -.24, .16, .43, -.92, .55, .23, -.38, .69, .32, -.17])
    signs = ((np.arange(2 ** len(weights))[:, None] >> np.arange(len(weights))) & 1) * 2 - 1
    expected = np.mean(abs((signs * weights).sum(axis=1)) >= abs(weights.sum()) - 1e-12)
    result = sd.sign_flip(weights, np.random.default_rng(1), 100, 512)
    assert result["test_method"] == "exact_meet_in_middle"
    assert result["p_value"] == expected


def test_monte_carlo_stream_is_independent_of_bootstrap_and_block_order():
    weights = np.linspace(-.48, 1.43, 38)
    blocks = [f"b{i}" for i in range(len(weights))]
    def run(values, block_ids, replicates):
        return sd.inference(values, block_ids, replicates, np.random.default_rng(3),
                            null_rng=np.random.default_rng(4), permutation_replicates=1000, exact_max_states=1)
    first = run(weights, blocks, 200)
    second = run(weights, blocks, 800)
    reversed_rows = run(weights[::-1], blocks[::-1], 200)
    assert first["test_method"] == "monte_carlo" and first["null_draws"] == 1000
    assert first["p_value"] == second["p_value"]
    assert first == reversed_rows
    assert 0 <= first["mc_p_ci_low"] < first["mc_p_ci_high"] <= 1


def test_plot_config_and_font_are_hashed_in_plan(tmp_path):
    plan = fixture_inputs(tmp_path / "input")
    config = tmp_path / "plot.json"
    config.write_text('{"font_size": 8, "formats": ["pdf"]}')
    changed = sd.make_plan(tmp_path, tmp_path / "input/manifest.tsv", plot_config=config)
    assert str(config) in changed["input_hashes"]
    assert changed["plot_config"]["formats"] == ["pdf"]
    assert changed["permutation_replicates"] == 100000
    assert changed["inference_version"] == 3 and plan["schema_version"] == 2


@pytest.mark.parametrize("scale", [1e-300, 1e-14, 1, 1e10, 1e300])
def test_sign_flip_is_scale_invariant_at_small_and_large_magnitudes(scale):
    result = sd.sign_flip(np.full(17, .13 * scale), np.random.default_rng(1), 100, 512)
    assert result["p_value"] == 2 / 2 ** 17


def test_integer_null_avoids_float_sum_overflow():
    result = sd.sign_flip(np.full(18, 1e308), np.random.default_rng(1), 100, 512)
    assert result["test_method"] == "exact_integer_dp"
    assert result["p_value"] == 2 / 2 ** 18


def test_enumeration_and_meet_in_middle_match_independent_rational_null():
    # Exact integer arithmetic provides the independent oracle for decimal weights.
    integers = [13, -72, 49, 57, -11, 38, 61, -24, 16, 43, -92, 55, 23, -38, 69, 32, -17]
    for n in (8, 17):
        weights = integers[:n]
        observed = abs(sum(weights))
        expected = sum(abs(sum(s * w for s, w in zip(signs, weights, strict=True))) >= observed
                       for signs in itertools.product((-1, 1), repeat=n)) / 2 ** n
        for scale in (1e-14, .01, 1e12):
            result = sd.sign_flip(np.array(weights) * scale, np.random.default_rng(1), 100, 512)
            assert result["p_value"] == expected


def test_monte_carlo_scaled_null_preserves_draws_and_probability():
    weights = np.linspace(-.48, 1.43, 38)
    results = [sd.sign_flip(weights * scale, np.random.default_rng(4), 2000, 1)
               for scale in (1e-200, 1, 1e200)]
    assert results[0] == results[1] == results[2]


@pytest.mark.parametrize("a,b", [(1e300, 1e-300), (1e308, 1e308), (5e-324, 5e-324)])
def test_expression_ratio_handles_finite_extreme_abundances(tmp_path, a, b):
    import math
    plan = fixture_inputs(tmp_path / "input")
    path = tmp_path / "input/expr.tsv"
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    for row in rows:
        if row["gene_id"] in {"g0A", "g0B"}:
            value = a if row["gene_id"] == "g0A" else b
            row.update({c: str(value) for c in ("leaf1", "leaf1tech", "leaf2")})
    sd.write_table(path, rows, rows[0].keys())
    genes, _, _ = sd.load_mapping(plan["analyses"][0]["mapping_file"])
    expression, coverage = sd.expression_rows(plan["analyses"][0], genes)
    assert expression[0]["log2_ratio"] == pytest.approx(math.log2(a) - math.log2(b))
    assert coverage[0]["positive_samples"] == 2
    assert sd.log_mean_abundance([5e-324, 0]) == -1075


def test_zero_eligible_expression_contrasts_are_retained(tmp_path):
    plan = fixture_inputs(tmp_path / "input")
    path = tmp_path / "input/expr.tsv"
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    for row in rows:
        row.update({c: "0" for c in ("leaf1", "leaf1tech", "leaf2")})
    sd.write_table(path, rows, rows[0].keys())
    result = sd.analyse(plan["analyses"][0], tmp_path / "out", 100, 1,
                        plot_config={"metrics": ["expression_log2_ratio"], "formats": ["svg"]})
    row = next(r for r in result["statistics"] if r["metric"] == "expression_log2_ratio")
    assert row["n_opportunities"] == 12 and row["n_loci"] == row["n_blocks"] == 0
    assert row["effect"] is row["ci_low"] is row["p_value"] is row["q_value"] is None
    assert row["status"].startswith("not_estimable") and row["multiple_testing_n"] == 0


def test_benjamini_hochberg_matches_independent_definition_and_resets_missing_q():
    from fractions import Fraction
    ps = [Fraction(1, 200), Fraction(1, 50), Fraction(1, 50), Fraction(1, 2), Fraction(0), Fraction(1)]
    ordered = sorted(ps)
    expected = {p: min(Fraction(1), *(q * len(ps) / rank for rank, q in enumerate(ordered, 1) if q >= p)) for p in ps}
    rows = [{"p_value": float(p)} for p in ps] + [{"p_value": None, "q_value": .001}]
    sd.adjust_p(rows)
    for row, p in zip(rows, ps, strict=False):
        assert row["q_value"] == pytest.approx(float(expected[p]))
        assert row["multiple_testing_n"] == len(ps)
    assert rows[-1]["q_value"] is None
    assert rows[-1]["multiple_testing_n"] == len(ps)


@pytest.mark.parametrize("p", [float("nan"), float("inf"), -.1, 1.1, True])
def test_invalid_p_values_fail(p):
    with pytest.raises(ValueError, match="p_value"):
        sd.adjust_p([{"p_value": p}])


@pytest.mark.parametrize("parameters", [{"replicates": True}, {"replicates": 100.0}, {"seed": -1},
                                        {"permutation_replicates": 99}, {"exact_max_states": False}])
def test_resampling_parameters_are_validated_before_planning(tmp_path, parameters):
    fixture_inputs(tmp_path / "input")
    with pytest.raises(ValueError, match="must be an integer"):
        sd.make_plan(tmp_path, tmp_path / "input/manifest.tsv", **parameters)


def test_implementation_guard_and_legacy_plan_compatibility(tmp_path):
    plan = fixture_inputs(tmp_path / "input")
    assert len(plan["implementation_hashes"]) == 2
    sd.verify_inputs(plan)
    implementation = tmp_path / "helper.py"
    implementation.write_text("version one")
    plan["implementation_hashes"] = {str(implementation): sd.digest(implementation)}
    implementation.write_text("version two")
    with pytest.raises(ValueError, match="Input changed"):
        sd.verify_inputs(plan)
    plan.pop("implementation_hashes")
    sd.verify_inputs(plan)


def test_short_optional_manifest_row_fails_cleanly(tmp_path):
    fixture_inputs(tmp_path / "input")
    path = tmp_path / "input/manifest.tsv"
    path.write_text(path.read_text().replace("samples_file\n", "samples_file\treference\n"))
    with pytest.raises(ValueError, match="malformed"):
        sd.make_plan(tmp_path, path)


def test_five_subgenomes_use_all_ten_independent_pairwise_contrasts(tmp_path):
    plan = fixture_inputs(tmp_path / "input")
    mapping = ["gene_id\tgroup_id\tsubgenome\tassignment_scope\tassignment_basis\tevidence"]
    pairs = ["pair_id\tblock_id\tgene_id"]
    retention = ["group_id\tblock_id\tlocus_id\tsubgenome\tcallable\tretained"]
    expression = ["gene_id\tleaf1\tleaf1tech\tleaf2"]
    for i in range(12):
        for label in ("D", "R1", "R2", "R3", "R4"):
            gene = f"g{i}{label}"
            mapping.append(f"{gene}\tG\t{label}\tlocal\tcurated\tindependent_phasing")
            pairs.append(f"p{i}\tb{i}\t{gene}")
            retention.append(f"G\tb{i}\tl{i}\t{label}\t1\t{int(label == 'D')}")
            value = 4 if label == "D" else 1
            expression.append(f"{gene}\t{value}\t{value}\t{value}")
    for name, rows in (("mapping", mapping), ("pairs", pairs), ("retention", retention), ("expr", expression)):
        (tmp_path / f"input/{name}.tsv").write_text("\n".join(rows) + "\n")
    result = sd.analyse(plan["analyses"][0], tmp_path / "output", 100, 1,
                        plot_config={"formats": ["svg"], "metrics": list(sd.plotting.METRICS[:2])})
    for metric, magnitude in (("retention_difference", 1), ("expression_log2_ratio", 2)):
        rows = [r for r in result["statistics"] if r["metric"] == metric]
        assert len(rows) == 10 and {r["multiple_testing_n"] for r in rows} == {10}
        assert sum(r["q_value"] < .05 for r in rows) == 4
        assert all(r["effect"] == (magnitude if r["subgenome_a"] == "D" else 0) for r in rows)
        assert all(r["p_value"] == (2 / 4096 if r["subgenome_a"] == "D" else 1) for r in rows)
    assert result["genomewide_identity"] == "unresolved" and result["genomewide_dominance"] == "not_tested"


def test_missing_expression_id_is_not_counted_as_zero(tmp_path):
    plan = fixture_inputs(tmp_path / "input")
    path = tmp_path / "input/expr.tsv"
    lines = path.read_text().splitlines()
    path.write_text("\n".join(line for line in lines if not line.startswith("g0B\t")) + "\n")
    result = sd.analyse(plan["analyses"][0], tmp_path / "output", 100, 1,
                        plot_config={"formats": ["svg"]})
    log_ratio, detection = result["statistics"][1:]
    assert log_ratio["n_opportunities"] == 12 and log_ratio["n_loci"] == 11
    assert detection["n_loci"] == 11 and detection["effect"] == 0
    coverage = list(csv.DictReader((tmp_path / "output/expression_coverage.tsv").open(), delimiter="\t"))
    assert coverage[0]["gene_ids_present"] == "0" and coverage[0]["positive_samples"] == "0"

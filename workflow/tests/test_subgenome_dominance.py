import csv
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

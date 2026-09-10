import sqlite3
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from gene_family_output_store import GeneFamilyOutputStore
from pgls_multiplicity import adjust_associations, write_association_table


def rows():
    return pd.DataFrame([
        dict(tree_id=family, analysis_method=method, aggregation=aggregation,
             analysis_id="a", model_id="m", response=response, source_term="temperature",
             term="temperature", term_test="coefficient", p_value=p,
             inference_status=status, coefficient=2.0)
        for family, method, aggregation, response, p, status in [
            ("OG1", "rsc", "gene_lineage_contrasts", "leaf", .01, "ok"),
            ("OG2", "rsc", "gene_lineage_contrasts", "leaf", .04, "ok"),
            ("OG1", "species_nwkit", "sum", "leaf", .2, "ok"),
            ("OG3", "rsc", "gene_lineage_contrasts", "leaf", .0001, "failed"),
            ("OG1", "rsc", "gene_lineage_contrasts", "root", .02, "ok"),
        ]
    ])


def test_pairs_pool_methods_and_keep_family_adjustments():
    out = adjust_associations(rows())
    assert out.p_value_global_bh.iloc[:3].tolist() == pytest.approx([.04, .08, .2 * 4 / 3])
    assert pd.isna(out.p_value_global_bh.iloc[3])
    assert out.p_value_global_bh.iloc[4] == .02
    assert out.global_n_associations.tolist() == [4, 4, 4, 4, 1]
    assert out.global_n_usable.tolist() == [3, 3, 3, 3, 1]
    assert out.p_value_family_holm.iloc[0] == .02
    assert out.p_value_family_bh.iloc[0] == .02
    assert not any("best" in c or "min_p" in c for c in out)
    shuffled = adjust_associations(rows().sample(frac=1, random_state=2)).sort_values(["tree_id", "response", "analysis_method"])
    expected = out.sort_values(["tree_id", "response", "analysis_method"])
    pd.testing.assert_frame_equal(shuffled.reset_index(drop=True), expected.reset_index(drop=True))


def test_duplicate_identity_and_missing_schema_are_errors():
    with pytest.raises(ValueError, match="Duplicate"):
        adjust_associations(pd.concat([rows(), rows().iloc[:1]]))
    with pytest.raises(ValueError, match="missing columns"):
        adjust_associations(rows().drop(columns="source_term"))


def test_intercepts_invalid_values_and_empty():
    frame = rows()
    frame.loc[0, "term"] = "(Intercept)"
    frame.loc[1, "p_value"] = float("inf")
    out = adjust_associations(frame)
    assert len(out) == 4
    assert pd.isna(out.p_value_global_bh.iloc[0])
    assert out.global_n_associations.iloc[0] == 3
    assert adjust_associations(rows().iloc[:0]).empty


def test_database_rebuild_from_comparison_only(tmp_path):
    sqlalchemy = pytest.importorskip("sqlalchemy")
    directory = tmp_path / "pgls_comparison"
    directory.mkdir()
    rows().to_csv(directory / "OG_comparison.tsv", sep="\t", index=False)
    native = tmp_path / "rsc_regression"
    native.mkdir()
    rows().to_csv(native / "duplicate.tsv", sep="\t", index=False)
    db = tmp_path / "result.db"
    engine = sqlalchemy.create_engine(f"sqlite:///{db}")
    for _ in range(2):
        write_association_table(engine, GeneFamilyOutputStore(tmp_path))
    with sqlite3.connect(db) as conn:
        actual = pd.read_sql_query("SELECT * FROM pgls_association", conn)
    assert len(actual) == 5
    assert actual.p_value_global_bh.iloc[0] == .04
    engine.dispose()


@pytest.mark.parametrize("archived", [False, True])
def test_summary_database_cli_live_and_zip(tmp_path, archived):
    import subprocess

    from gene_family_output_store import archive_completed_outputs, family_context

    root = tmp_path / "orthogroup"
    genecount = tmp_path / "counts.tsv"
    genecount.write_text("Orthogroup\tTotal\nOG1\t1\nOG2\t1\nOG3\t1\n")
    for family, group in rows().groupby("tree_id"):
        for directory, suffix, frame in [
            ("pgls_comparison", "_comparison.tsv", group),
            ("stat_tree", "_stat.tree.tsv", pd.DataFrame([dict(num_branch=3, num_spe=1, num_dup=0, num_sp=2, rsc_best_p_value=.01)])),
            ("stat_branch", "_stat.branch.tsv", pd.DataFrame([dict(branch_id=0, node_name="n0", num_sp=2, so_event="S")])),
        ]:
            path = root / directory / (family + suffix)
            path.parent.mkdir(parents=True, exist_ok=True)
            frame.to_csv(path, sep="\t", index=False)
        plot = root / "tree_plot" / (family + "_tree_plot.pdf")
        plot.parent.mkdir(exist_ok=True)
        plot.write_bytes(b"%PDF-test\n")
    if archived:
        ids, family_from_name = family_context("orthogroup", genecount=genecount)
        archive_completed_outputs(root, "orthogroup", ids, family_from_name)
        assert not list((root / "pgls_comparison").glob("*.tsv"))
    db = root / "result.db"
    script = Path(__file__).resolve().parents[1] / "support" / "generate_orthogroup_database.py"
    subprocess.run([sys.executable, str(script), "--overwrite", "1", "--dbpath", str(db),
                    "--dir_gene_family", str(root), "--dir_stat_tree", str(root / "stat_tree"),
                    "--dir_stat_branch", str(root / "stat_branch")], check=True, cwd=tmp_path, capture_output=True)
    with sqlite3.connect(db) as conn:
        actual = pd.read_sql_query("SELECT * FROM pgls_association", conn)
        assert len(actual) == 5
        assert actual.loc[(actual.tree_id == "OG1") & (actual.response == "leaf") & (actual.analysis_method == "rsc"), "p_value_global_bh"].iloc[0] == .04
        assert "rsc_best_p_value" not in pd.read_sql_query("SELECT * FROM tree", conn).columns


def test_bh_matches_r_reference():
    import subprocess

    raw = [0, .001, .01, .01, .2, 1]
    frame = pd.concat([rows().iloc[[0]]] * len(raw), ignore_index=True)
    frame["tree_id"] = [f"OG{i}" for i in range(len(raw))]
    frame["p_value"] = raw
    expected = subprocess.check_output(["Rscript", "-e", "cat(p.adjust(c(0,.001,.01,.01,.2,1), 'BH'), sep='\\n')"], text=True)
    assert adjust_associations(frame).p_value_global_bh.tolist() == pytest.approx([float(x) for x in expected.split()], abs=1e-12)


def test_predictors_are_separate_but_levels_and_omnibus_share_a_pair():
    frame = pd.concat([rows().iloc[[0]]] * 4, ignore_index=True)
    frame["term"] = ["climate[dry]", "climate[wet]", "climate", "size"]
    frame["source_term"] = ["climate", "climate", "climate", "size"]
    frame["term_test"] = ["coefficient", "coefficient", "omnibus", "coefficient"]
    frame["p_value"] = [.01, .04, .2, .01]
    out = adjust_associations(frame)
    assert out.p_value_global_bh.tolist() == pytest.approx([.03, .06, .2, .01])
    assert out.global_n_associations.tolist() == [3, 3, 3, 1]

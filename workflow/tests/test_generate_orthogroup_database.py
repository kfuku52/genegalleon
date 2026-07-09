from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import math
import sqlite3
import subprocess
import sys

import pandas


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "generate_orthogroup_database.py"


def load_module():
    spec = spec_from_file_location("generate_orthogroup_database", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_read_header_columns_returns_tsv_header_columns(tmp_path):
    mod = load_module()
    infile = tmp_path / "sample.tsv"
    infile.write_text("col1\tcol2\tcol3\n1\t2\t3\n", encoding="utf-8")

    assert mod.read_header_columns(str(infile)) == ["col1", "col2", "col3"]


def test_process_files_uses_single_read_csv_call(tmp_path, monkeypatch):
    mod = load_module()
    infile = tmp_path / "OG0001.test.tsv"
    pandas.DataFrame({"a": [1, 2], "b": [3, 4]}).to_csv(infile, sep="\t", index=False)

    call_count = {"n": 0}
    original_read_csv = mod.pd.read_csv

    def spy_read_csv(*args, **kwargs):
        call_count["n"] += 1
        return original_read_csv(*args, **kwargs)

    monkeypatch.setattr(mod.pd, "read_csv", spy_read_csv)
    out = mod.process_files(str(infile), ["orthogroup", "a", "b"])

    assert call_count["n"] == 1
    assert out.columns.tolist() == ["orthogroup", "a", "b"]
    assert out["orthogroup"].iloc[0] == "OG0001"


def test_process_files_returns_empty_dataframe_when_columns_missing(tmp_path):
    mod = load_module()
    infile = tmp_path / "OG0002.test.tsv"
    pandas.DataFrame({"a": [1], "b": [2]}).to_csv(infile, sep="\t", index=False)

    out = mod.process_files(str(infile), ["orthogroup", "a", "b", "c"])
    assert out.empty


def test_parse_cutoff_stat_parses_valid_tokens_and_ignores_invalid():
    mod = load_module()
    parsed = mod.parse_cutoff_stat("OCNany2spe,0.8|badtoken|,1|ABC,notfloat|XYZ,1.5")
    assert parsed == [("OCNany2spe", 0.8), ("XYZ", 1.5)]


def test_apply_cutoff_accepts_preparsed_cutoff_list():
    mod = load_module()
    df = pandas.DataFrame(
        [
            {"A": "0.9", "B": "0.1"},
            {"A": "0.7", "B": "0.9"},
            {"A": "0.95", "B": "0.95"},
        ]
    )
    out = mod.apply_cutoff(df, [("A", 0.8), ("B", 0.5)])
    assert out.shape[0] == 1
    assert out.iloc[0]["A"] == "0.95"


def test_gene_family_id_from_path_recognizes_csubst_scan_suffixes():
    mod = load_module()

    assert mod.gene_family_id_from_path("/tmp/OG0001_csubst_scan.tsv") == "OG0001"
    assert mod.gene_family_id_from_path("/tmp/OG0001_csubst_scan_units.tsv") == "OG0001"
    assert mod.gene_family_id_from_path("/tmp/HOG0002.csubst_scan.tsv") == "HOG0002"
    assert mod.gene_family_id_from_path("/tmp/HOG0002.csubst_scan_units.tsv") == "HOG0002"


def test_calculate_bh_fdr_preserves_nan_and_original_order():
    mod = load_module()

    out = mod.calculate_bh_fdr([0.01, 0.04, 0.03, float("nan")])

    assert out[:3].tolist() == [0.03, 0.04, 0.04]
    assert math.isnan(out[3])


def test_database_builder_adds_aa_change_tables_and_global_fdr(tmp_path):
    stat_tree = tmp_path / "stat_tree"
    stat_branch = tmp_path / "stat_branch"
    aa_change = tmp_path / "csubst_scan"
    aa_change_unit = tmp_path / "csubst_scan_units"
    for path in [stat_tree, stat_branch, aa_change, aa_change_unit]:
        path.mkdir()

    pandas.DataFrame({"tree_metric": [1.0]}).to_csv(stat_tree / "OG0001_stat.tree.tsv", sep="\t", index=False)
    pandas.DataFrame({"branch_id": [1], "branch_metric": [2.0]}).to_csv(stat_branch / "OG0001_stat.branch.tsv", sep="\t", index=False)
    pandas.DataFrame(
        [
            {"trait": "traitA", "state_change": "10K", "p_rate_enrichment": 0.01, "p_rate_enrichment_empirical": 0.02},
            {"trait": "traitA", "state_change": "12S", "p_rate_enrichment": 0.04, "p_rate_enrichment_empirical": 0.50},
        ]
    ).to_csv(aa_change / "OG0001_csubst_scan.tsv", sep="\t", index=False)
    pandas.DataFrame(
        [
            {"trait": "traitA", "unit_id": 1, "matched_leaf_names": "Species_A"},
            {"trait": "traitA", "unit_id": 2, "matched_leaf_names": "Species_B"},
        ]
    ).to_csv(aa_change_unit / "OG0001_csubst_scan_units.tsv", sep="\t", index=False)

    db_path = tmp_path / "gg_orthogroup.db"
    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--overwrite",
            "1",
            "--dbpath",
            str(db_path),
            "--dir_stat_tree",
            str(stat_tree),
            "--dir_stat_branch",
            str(stat_branch),
            "--dir_csubst_aa_change",
            str(aa_change),
            "--dir_csubst_aa_change_unit",
            str(aa_change_unit),
            "--ncpu",
            "1",
        ],
        cwd=str(tmp_path),
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr

    with sqlite3.connect(db_path) as conn:
        tables = {row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        assert "aa_change" in tables
        assert "aa_change_unit" in tables
        aa_df = pandas.read_sql_query(
            "SELECT orthogroup, state_change, q_rate_enrichment_global, q_rate_enrichment_empirical_global FROM aa_change ORDER BY state_change",
            conn,
        )
        unit_df = pandas.read_sql_query("SELECT orthogroup, unit_id FROM aa_change_unit ORDER BY unit_id", conn)

    assert aa_df["orthogroup"].tolist() == ["OG0001", "OG0001"]
    assert aa_df["q_rate_enrichment_global"].round(4).tolist() == [0.02, 0.04]
    assert aa_df["q_rate_enrichment_empirical_global"].round(4).tolist() == [0.04, 0.5]
    assert unit_df["orthogroup"].tolist() == ["OG0001", "OG0001"]


def test_import_has_no_logfile_side_effect(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    log_path = tmp_path / "generate_orthogroup_database.log"
    assert not log_path.exists()
    load_module()
    assert not log_path.exists()


def test_help_has_no_logfile_side_effect(tmp_path):
    log_path = tmp_path / "generate_orthogroup_database.log"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "--help"],
        cwd=str(tmp_path),
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0
    assert not log_path.exists()


def test_visible_entries_and_has_visible_entries_ignore_hidden_files(tmp_path):
    mod = load_module()
    directory = tmp_path / "d"
    directory.mkdir()
    (directory / ".DS_Store").write_text("", encoding="utf-8")
    (directory / "data.tsv").write_text("a\tb\n1\t2\n", encoding="utf-8")

    entries = mod.visible_entries(str(directory))
    assert entries == ["data.tsv"]
    assert mod.has_visible_entries(str(directory)) is True

    only_hidden = tmp_path / "only_hidden"
    only_hidden.mkdir()
    (only_hidden / ".keep").write_text("", encoding="utf-8")
    assert mod.visible_entries(str(only_hidden)) == []
    assert mod.has_visible_entries(str(only_hidden)) is False


def test_visible_files_returns_only_non_hidden_regular_files(tmp_path):
    mod = load_module()
    directory = tmp_path / "d"
    directory.mkdir()
    (directory / "table.tsv").write_text("a\tb\n1\t2\n", encoding="utf-8")
    (directory / ".DS_Store").write_text("", encoding="utf-8")
    (directory / "nested").mkdir()

    assert mod.visible_files(str(directory)) == ["table.tsv"]

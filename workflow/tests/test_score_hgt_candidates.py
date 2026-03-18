from pathlib import Path
import os
import sqlite3
import subprocess
import sys

import pandas


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "score_hgt_candidates.py"


def write_branch_db(db_path: Path) -> None:
    branch = pandas.DataFrame(
        [
            {
                "orthogroup": "OG0001",
                "branch_id": 3,
                "node_name": "n3",
                "gene_labels": "geneA; geneB",
                "num_leaf": 2,
                "so_event": "S",
                "taxon": "",
                "spnode_coverage": "Viridiplantae",
                "generax_event": "H",
                "generax_transfer": "Y",
                "generax_event_parent": "S",
                "clade_min_expression_pearsoncor": 0.6,
                "num_intron": pandas.NA,
                "intron_present": pandas.NA,
                "synteny_support_score": pandas.NA,
                "sprot_best": pandas.NA,
                "organism": pandas.NA,
                "taxid_y": pandas.NA,
                "expression_leaf": pandas.NA,
            },
            {
                "orthogroup": "OG0001",
                "branch_id": 1,
                "node_name": "geneA",
                "gene_labels": "geneA",
                "num_leaf": 1,
                "so_event": "L",
                "taxon": "Arabidopsis thaliana",
                "spnode_coverage": "Arabidopsis thaliana",
                "generax_event": "L",
                "generax_transfer": "",
                "generax_event_parent": "H",
                "clade_min_expression_pearsoncor": pandas.NA,
                "num_intron": 2,
                "intron_present": 1,
                "synteny_support_score": 0.8,
                "sprot_best": "P00001",
                "organism": "Escherichia coli",
                "taxid_y": pandas.NA,
                "expression_leaf": 10,
            },
            {
                "orthogroup": "OG0001",
                "branch_id": 2,
                "node_name": "geneB",
                "gene_labels": "geneB",
                "num_leaf": 1,
                "so_event": "L",
                "taxon": "Arabidopsis thaliana",
                "spnode_coverage": "Arabidopsis thaliana",
                "generax_event": "L",
                "generax_transfer": "",
                "generax_event_parent": "H",
                "clade_min_expression_pearsoncor": pandas.NA,
                "num_intron": 1,
                "intron_present": 1,
                "synteny_support_score": 0.4,
                "sprot_best": "P00002",
                "organism": "Bacillus subtilis",
                "taxid_y": pandas.NA,
                "expression_leaf": 5,
            },
        ]
    )
    with sqlite3.connect(db_path) as conn:
        branch.to_sql("branch", conn, index=False)


def run_script(tmp_path: Path, contamination_dir: Path | None = None):
    db_path = tmp_path / "gg_orthogroup.db"
    write_branch_db(db_path)
    branch_out = tmp_path / "hgt_branch_candidates.tsv"
    gene_out = tmp_path / "hgt_gene_candidates.tsv"
    orthogroup_out = tmp_path / "hgt_orthogroup_summary.tsv"
    cmd = [
        sys.executable,
        str(SCRIPT_PATH),
        "--dbpath",
        str(db_path),
        "--branch_out",
        str(branch_out),
        "--gene_out",
        str(gene_out),
        "--orthogroup_out",
        str(orthogroup_out),
    ]
    if contamination_dir is not None:
        cmd.extend(["--dir_contamination_tsv", str(contamination_dir)])
    completed = subprocess.run(cmd, capture_output=True, text=True, check=False)
    assert completed.returncode == 0, completed.stderr
    return (
        pandas.read_csv(branch_out, sep="\t"),
        pandas.read_csv(gene_out, sep="\t"),
        pandas.read_csv(orthogroup_out, sep="\t"),
    )


def run_script_with_env(tmp_path: Path, extra_args: list[str], env_overrides: dict[str, str]):
    db_path = tmp_path / "gg_orthogroup.db"
    write_branch_db(db_path)
    branch_out = tmp_path / "hgt_branch_candidates.tsv"
    gene_out = tmp_path / "hgt_gene_candidates.tsv"
    orthogroup_out = tmp_path / "hgt_orthogroup_summary.tsv"
    cmd = [
        sys.executable,
        str(SCRIPT_PATH),
        "--dbpath",
        str(db_path),
        "--branch_out",
        str(branch_out),
        "--gene_out",
        str(gene_out),
        "--orthogroup_out",
        str(orthogroup_out),
        *extra_args,
    ]
    env = os.environ.copy()
    env.update(env_overrides)
    completed = subprocess.run(cmd, capture_output=True, text=True, check=False, env=env)
    assert completed.returncode == 0, completed.stderr
    return pandas.read_csv(branch_out, sep="\t")


def test_score_hgt_candidates_emits_branch_gene_and_orthogroup_outputs(tmp_path):
    branch_out, gene_out, orthogroup_out = run_script(tmp_path)

    assert branch_out.shape[0] == 1
    row = branch_out.iloc[0]
    assert row["orthogroup"] == "OG0001"
    assert row["branch_id"] == 3
    assert row["matched_leaf_count"] == 2
    assert row["besthit_taxonomy_method"] in {"taxonomy_db", "name_heuristic"}
    assert isinstance(row["besthit_lca_rank_mode"], str)
    assert row["besthit_lca_rank_mode"] != ""
    assert "hgt_score" not in branch_out.columns
    assert "hgt_confidence" not in branch_out.columns
    assert "hgt_reason" not in branch_out.columns

    assert set(gene_out["gene_id"].tolist()) == {"geneA", "geneB"}
    assert gene_out["candidate_branch_count"].eq(1).all()
    assert set(gene_out["candidate_branch_ids"].astype(str).tolist()) == {"3"}
    assert "besthit_conflict_score" not in gene_out.columns
    assert orthogroup_out.loc[0, "orthogroup"] == "OG0001"
    assert orthogroup_out.loc[0, "hgt_branch_count"] == 1
    assert orthogroup_out.loc[0, "hgt_gene_count"] == 2
    assert str(orthogroup_out.loc[0, "candidate_branch_ids"]) == "3"
    assert "top_hgt_score" not in orthogroup_out.columns
    assert "top_hgt_confidence" not in orthogroup_out.columns


def test_score_hgt_candidates_applies_contamination_penalty(tmp_path):
    contamination_dir = tmp_path / "contamination"
    contamination_dir.mkdir(parents=True)
    contamination_tsv = contamination_dir / "species.tsv"
    contamination_tsv.write_text(
        "query\tis_compatible_lineage\tlca_taxid\tlca_sciname\n"
        "geneA\tfalse\t562\tEscherichia coli\n"
        "geneB\tfalse\t562\tEscherichia coli\n",
        encoding="utf-8",
    )

    branch_out, gene_out, _ = run_script(tmp_path, contamination_dir=contamination_dir)

    assert branch_out.shape[0] == 1
    row = branch_out.iloc[0]
    assert row["contamination_incompatible_fraction"] == 1.0
    assert str(row["contamination_top_lca_taxid"]) == "562"
    assert row["contamination_top_lca_sciname"] == "Escherichia coli"
    assert gene_out["contamination_is_compatible_lineage"].eq(False).all()
    assert set(gene_out["contamination_lca_sciname"].tolist()) == {"Escherichia coli"}


def test_score_hgt_candidates_explicit_empty_taxonomy_dbfile_overrides_env(tmp_path, monkeypatch):
    monkeypatch.setenv("GG_TAXONOMY_DBFILE", "/nonexistent/taxa.sqlite")
    branch_out = run_script_with_env(
        tmp_path,
        extra_args=["--taxonomy_dbfile", ""],
        env_overrides={"GG_TAXONOMY_DBFILE": "/nonexistent/taxa.sqlite"},
    )

    assert branch_out.shape[0] == 1
    assert branch_out.loc[0, "besthit_taxonomy_method"] == "name_heuristic"

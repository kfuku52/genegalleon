from pathlib import Path
import subprocess
import sys

import pandas


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "annotate_hgt_tree_plot.py"


def test_annotate_hgt_tree_plot_merges_branch_and_gene_evidence(tmp_path: Path):
    stat_branch = tmp_path / "stat.branch.tsv"
    branch_tsv = tmp_path / "hgt_branch_candidates.tsv"
    gene_tsv = tmp_path / "hgt_gene_candidates.tsv"
    outfile = tmp_path / "annotated.tsv"

    pandas.DataFrame(
        [
            {"orthogroup": "OG0001", "branch_id": 3, "node_name": "n3", "so_event": "S"},
            {"orthogroup": "OG0001", "branch_id": 1, "node_name": "geneA", "so_event": "L"},
            {"orthogroup": "OG0001", "branch_id": 2, "node_name": "geneB", "so_event": "L"},
        ]
    ).to_csv(stat_branch, sep="\t", index=False)
    pandas.DataFrame(
        [
            {
                "orthogroup": "OG0001",
                "branch_id": 3,
                "candidate_gene_count": 2,
                "matched_leaf_count": 2,
                "besthit_gene_count": 2,
                "besthit_taxid_count": 2,
                "besthit_taxonomy_method": "taxonomy_db",
                "besthit_same_superkingdom_fraction": 0.0,
                "besthit_lca_rank_mode": "phylum",
                "intron_support_fraction": 1.0,
                "expression_measured_fraction": 1.0,
                "clade_min_expression_pearsoncor": 0.5,
                "synteny_support_fraction": 1.0,
                "synteny_mean_support_score": 0.7,
                "contamination_incompatible_fraction": 0.0,
                "contamination_top_lca_sciname": "",
            }
        ]
    ).to_csv(branch_tsv, sep="\t", index=False)
    pandas.DataFrame(
        [
            {
                "orthogroup": "OG0001",
                "gene_id": "geneA",
                "gene_taxon": "Arabidopsis thaliana",
                "candidate_branch_count": 1,
                "candidate_branch_ids": "3",
                "besthit_accession": "P1",
                "besthit_organism": "Escherichia coli str. K-12",
                "besthit_taxid": "562",
                "besthit_taxonomy_method": "taxonomy_db",
                "besthit_lca_rank": "phylum",
                "besthit_same_superkingdom": 0,
                "intron_supported": True,
                "expression_measured": True,
                "synteny_support_score": 0.8,
                "contamination_lca_taxid": "",
                "contamination_lca_sciname": "",
                "contamination_is_compatible_lineage": True,
            }
        ]
    ).to_csv(gene_tsv, sep="\t", index=False)

    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--stat_branch",
            str(stat_branch),
            "--branch_tsv",
            str(branch_tsv),
            "--gene_tsv",
            str(gene_tsv),
            "--orthogroup",
            "OG0001",
            "--outfile",
            str(outfile),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr

    out = pandas.read_csv(outfile, sep="\t")
    branch_row = out.loc[out["branch_id"] == 3].iloc[0]
    gene_a = out.loc[out["node_name"] == "geneA"].iloc[0]
    gene_b = out.loc[out["node_name"] == "geneB"].iloc[0]

    assert branch_row["hgtbranch_candidate_flag"] == 1
    assert branch_row["hgtbranch_candidate_gene_count"] == 2
    assert gene_a["hgt_Cand"] == 1
    assert gene_a["hgt_SameSK"] == 0
    assert gene_a["hgt_ContamOK"] == 1
    assert gene_a["besthit_lca_rank_short"] == "phylum"
    assert gene_a["besthit_lca_rank_display"] == "phylum"
    assert gene_a["besthit_organism_short"].startswith("E. coli")
    assert gene_a["besthit_organism_display"].startswith("E. coli")
    assert gene_a["contamination_lca_sciname_display"] == "-"
    assert gene_b["hgt_Cand"] == 0

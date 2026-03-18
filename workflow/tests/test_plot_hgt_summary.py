from pathlib import Path
import subprocess
import sys

import pandas


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "plot_hgt_summary.py"


def test_plot_hgt_summary_generates_overview_and_taxonomy_flow_pdfs(tmp_path: Path):
    branch_tsv = tmp_path / "hgt_branch_candidates.tsv"
    gene_tsv = tmp_path / "hgt_gene_candidates.tsv"
    overview_pdf = tmp_path / "hgt_branch_overview.pdf"
    flow_pdf = tmp_path / "hgt_taxonomy_flow.pdf"
    readme_md = tmp_path / "README.md"

    pandas.DataFrame(
        [
            {
                "orthogroup": "OG0001",
                "branch_id": 3,
                "candidate_gene_count": 2,
                "matched_leaf_count": 2,
                "besthit_gene_count": 2,
                "besthit_taxid_count": 0,
                "besthit_taxonomy_method": "name_heuristic",
                "besthit_same_superkingdom_fraction": 0.0,
                "besthit_lca_rank_mode": "genus_mismatch",
                "intron_support_fraction": 1.0,
                "expression_measured_fraction": 1.0,
                "clade_min_expression_pearsoncor": 0.4,
                "synteny_support_fraction": 0.5,
                "synteny_mean_support_score": 0.6,
                "contamination_incompatible_fraction": 0.0,
                "contamination_top_lca_sciname": "",
            },
            {
                "orthogroup": "OG0002",
                "branch_id": 7,
                "candidate_gene_count": 1,
                "matched_leaf_count": 1,
                "besthit_gene_count": 1,
                "besthit_taxid_count": 0,
                "besthit_taxonomy_method": "name_heuristic",
                "besthit_same_superkingdom_fraction": 1.0,
                "besthit_lca_rank_mode": "species",
                "intron_support_fraction": 0.0,
                "expression_measured_fraction": 1.0,
                "clade_min_expression_pearsoncor": 0.0,
                "synteny_support_fraction": 0.0,
                "synteny_mean_support_score": 0.0,
                "contamination_incompatible_fraction": 1.0,
                "contamination_top_lca_sciname": "Escherichia coli",
            },
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
                "besthit_organism": "Escherichia coli",
                "besthit_taxid": "",
                "besthit_taxonomy_method": "name_heuristic",
                "besthit_lca_rank": "genus_mismatch",
                "besthit_same_superkingdom": 0,
                "intron_supported": True,
                "expression_measured": True,
                "synteny_support_score": 0.8,
                "contamination_lca_taxid": "",
                "contamination_lca_sciname": "",
                "contamination_is_compatible_lineage": True,
            },
            {
                "orthogroup": "OG0002",
                "gene_id": "geneB",
                "gene_taxon": "Arabidopsis thaliana",
                "candidate_branch_count": 1,
                "candidate_branch_ids": "7",
                "besthit_accession": "P2",
                "besthit_organism": "Bacillus subtilis",
                "besthit_taxid": "",
                "besthit_taxonomy_method": "name_heuristic",
                "besthit_lca_rank": "genus_mismatch",
                "besthit_same_superkingdom": 0,
                "intron_supported": False,
                "expression_measured": True,
                "synteny_support_score": 0.0,
                "contamination_lca_taxid": "562",
                "contamination_lca_sciname": "Escherichia coli",
                "contamination_is_compatible_lineage": False,
            },
        ]
    ).to_csv(gene_tsv, sep="\t", index=False)

    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--branch_tsv",
            str(branch_tsv),
            "--gene_tsv",
            str(gene_tsv),
            "--overview_pdf",
            str(overview_pdf),
            "--taxonomy_flow_pdf",
            str(flow_pdf),
            "--flow_rank",
            "phylum",
            "--flow_max_categories",
            "8",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert overview_pdf.exists() and overview_pdf.stat().st_size > 0
    assert flow_pdf.exists() and flow_pdf.stat().st_size > 0
    assert readme_md.exists()
    readme_text = readme_md.read_text(encoding="utf-8")
    assert "`Cand`" in readme_text
    assert "`TopContam`" in readme_text

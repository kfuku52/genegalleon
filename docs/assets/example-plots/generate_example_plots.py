#!/usr/bin/env python3
"""Generate documentation example plots from small test inputs."""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve()
REPO_ROOT = SCRIPT.parents[3]
ASSET_DIR = SCRIPT.parent
TEST_DATA_DIR = ASSET_DIR / "test-data"
SUPPORT_DIR = REPO_ROOT / "workflow" / "support"


SPECIES_TREE = (
    "((Amborella_trichopoda:140,"
    "((Arabidopsis_thaliana:40,Oryza_sativa:40)92:35,"
    "Spinacia_oleracea:75)88:65)95:25,"
    "(Cephalotus_follicularis:70,"
    "(Dionaea_muscipula:35,Nepenthes_gracilis:35)90:35)93:95);"
)


def run(cmd: list[str], cwd: Path = REPO_ROOT) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def write_tsv(
    path: Path,
    rows: list[dict[str, object]],
    fieldnames: list[str],
    lineterminator: str = "\r\n",
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator=lineterminator,
        )
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def require_tool(name: str) -> str:
    resolved = shutil.which(name)
    if resolved is None:
        raise RuntimeError(f"Required command was not found on PATH: {name}")
    return resolved


def convert_pdf_first_page_to_png(pdf_path: Path, png_path: Path, dpi: int = 180) -> None:
    pdftoppm = require_tool("pdftoppm")
    png_path.parent.mkdir(parents=True, exist_ok=True)
    prefix = png_path.with_suffix("")
    run([pdftoppm, "-png", "-singlefile", "-r", str(dpi), str(pdf_path), str(prefix)])


def write_query2family_test_data() -> dict[str, Path]:
    base = TEST_DATA_DIR / "query2family"
    query_dir = base / "query_gene"
    stat_dir = base / "query2family" / "stat_branch"
    write_text(base / "species_tree.nwk", SPECIES_TREE + "\n")
    write_text(base / "species_tree_support.nwk", SPECIES_TREE + "\n")

    families = {
        "AHA": {
            "Amborella_trichopoda": 1,
            "Arabidopsis_thaliana": 2,
            "Cephalotus_follicularis": 1,
            "Dionaea_muscipula": 1,
            "Nepenthes_gracilis": 2,
            "Oryza_sativa": 1,
            "Spinacia_oleracea": 1,
        },
        "STRICTCHK": {
            "Amborella_trichopoda": 1,
            "Arabidopsis_thaliana": 1,
            "Cephalotus_follicularis": 0,
            "Dionaea_muscipula": 0,
            "Nepenthes_gracilis": 0,
            "Oryza_sativa": 1,
            "Spinacia_oleracea": 1,
        },
        "YABBY": {
            "Amborella_trichopoda": 1,
            "Arabidopsis_thaliana": 1,
            "Cephalotus_follicularis": 1,
            "Dionaea_muscipula": 2,
            "Nepenthes_gracilis": 2,
            "Oryza_sativa": 0,
            "Spinacia_oleracea": 1,
        },
        "WOX": {
            "Amborella_trichopoda": 1,
            "Arabidopsis_thaliana": 1,
            "Cephalotus_follicularis": 1,
            "Dionaea_muscipula": 1,
            "Nepenthes_gracilis": 1,
            "Oryza_sativa": 1,
            "Spinacia_oleracea": 0,
        },
        "MADS": {
            "Amborella_trichopoda": 1,
            "Arabidopsis_thaliana": 2,
            "Cephalotus_follicularis": 1,
            "Dionaea_muscipula": 1,
            "Nepenthes_gracilis": 1,
            "Oryza_sativa": 2,
            "Spinacia_oleracea": 1,
        },
        "KNOX": {
            "Amborella_trichopoda": 1,
            "Arabidopsis_thaliana": 1,
            "Cephalotus_follicularis": 0,
            "Dionaea_muscipula": 1,
            "Nepenthes_gracilis": 1,
            "Oryza_sativa": 0,
            "Spinacia_oleracea": 1,
        },
        "PPR": {
            "Amborella_trichopoda": 2,
            "Arabidopsis_thaliana": 2,
            "Cephalotus_follicularis": 2,
            "Dionaea_muscipula": 2,
            "Nepenthes_gracilis": 2,
            "Oryza_sativa": 2,
            "Spinacia_oleracea": 2,
        },
        "AP2": {
            "Amborella_trichopoda": 1,
            "Arabidopsis_thaliana": 1,
            "Cephalotus_follicularis": 1,
            "Dionaea_muscipula": 0,
            "Nepenthes_gracilis": 0,
            "Oryza_sativa": 1,
            "Spinacia_oleracea": 1,
        },
        "ARF": {
            "Amborella_trichopoda": 1,
            "Arabidopsis_thaliana": 1,
            "Cephalotus_follicularis": 1,
            "Dionaea_muscipula": 1,
            "Nepenthes_gracilis": 1,
            "Oryza_sativa": 1,
            "Spinacia_oleracea": 1,
        },
        "BHLH": {
            "Amborella_trichopoda": 1,
            "Arabidopsis_thaliana": 2,
            "Cephalotus_follicularis": 0,
            "Dionaea_muscipula": 1,
            "Nepenthes_gracilis": 1,
            "Oryza_sativa": 2,
            "Spinacia_oleracea": 1,
        },
        "GRAS": {
            "Amborella_trichopoda": 1,
            "Arabidopsis_thaliana": 1,
            "Cephalotus_follicularis": 1,
            "Dionaea_muscipula": 1,
            "Nepenthes_gracilis": 0,
            "Oryza_sativa": 1,
            "Spinacia_oleracea": 0,
        },
        "MYB": {
            "Amborella_trichopoda": 1,
            "Arabidopsis_thaliana": 2,
            "Cephalotus_follicularis": 1,
            "Dionaea_muscipula": 2,
            "Nepenthes_gracilis": 2,
            "Oryza_sativa": 1,
            "Spinacia_oleracea": 1,
        },
        "NAC": {
            "Amborella_trichopoda": 0,
            "Arabidopsis_thaliana": 1,
            "Cephalotus_follicularis": 1,
            "Dionaea_muscipula": 1,
            "Nepenthes_gracilis": 1,
            "Oryza_sativa": 0,
            "Spinacia_oleracea": 1,
        },
    }
    for query_id in families:
        write_text(query_dir / query_id, f"{query_id}_seed_gene\n")
    stat_fields = ["branch_id", "node_name", "spnode_coverage", "so_event", "num_leaf"]
    for family_id, counts in families.items():
        rows = []
        branch_id = 1
        for species, count in counts.items():
            for copy_idx in range(1, count + 1):
                rows.append(
                    {
                        "branch_id": branch_id,
                        "node_name": f"{species}_{family_id}_{copy_idx}",
                        "spnode_coverage": species,
                        "so_event": "L",
                        "num_leaf": 1,
                    }
                )
                branch_id += 1
        rows.append(
            {
                "branch_id": branch_id,
                "node_name": f"{family_id}_internal_duplication",
                "spnode_coverage": "Arabidopsis_thaliana",
                "so_event": "D",
                "num_leaf": 2,
            }
        )
        write_tsv(stat_dir / f"{family_id}_stat.branch.tsv", rows, stat_fields)

    busco_rows = [
        {"species": "Amborella_trichopoda", "Single": 915, "Duplicated": 55, "Fragmented": 20, "Missing": 10},
        {"species": "Arabidopsis_thaliana", "Single": 890, "Duplicated": 80, "Fragmented": 18, "Missing": 12},
        {"species": "Cephalotus_follicularis", "Single": 840, "Duplicated": 95, "Fragmented": 42, "Missing": 23},
        {"species": "Dionaea_muscipula", "Single": 802, "Duplicated": 118, "Fragmented": 55, "Missing": 25},
        {"species": "Nepenthes_gracilis", "Single": 815, "Duplicated": 110, "Fragmented": 50, "Missing": 25},
        {"species": "Oryza_sativa", "Single": 902, "Duplicated": 60, "Fragmented": 24, "Missing": 14},
        {"species": "Spinacia_oleracea", "Single": 870, "Duplicated": 75, "Fragmented": 35, "Missing": 20},
    ]
    write_tsv(base / "busco_summary.tsv", busco_rows, ["species", "Single", "Duplicated", "Fragmented", "Missing"])
    busco_full_dir = base / "busco_full"
    busco_statuses = {
        "Amborella_trichopoda": ["Complete", "Complete", "Duplicated", "Fragmented", "Missing"],
        "Arabidopsis_thaliana": ["Complete", "Complete", "Duplicated", "Duplicated", "Missing"],
        "Cephalotus_follicularis": ["Complete", "Duplicated", "Fragmented", "Missing", "Missing"],
        "Dionaea_muscipula": ["Complete", "Duplicated", "Duplicated", "Fragmented", "Missing"],
        "Nepenthes_gracilis": ["Complete", "Complete", "Duplicated", "Fragmented", "Missing"],
        "Oryza_sativa": ["Complete", "Complete", "Complete", "Duplicated", "Missing"],
        "Spinacia_oleracea": ["Complete", "Complete", "Duplicated", "Fragmented", "Missing"],
    }
    for species, statuses in busco_statuses.items():
        lines = ["# Busco id\tStatus\tSequence\tScore\tLength"]
        for idx, status in enumerate(statuses, start=1):
            lines.append(f"BUSCO{idx:04d}\t{status}\t{species}_gene{idx}\t100\t250")
        write_text(busco_full_dir / species, "\n".join(lines) + "\n")
    return {
        "base": base,
        "query_dir": query_dir,
        "gene_family_dir": base / "query2family",
        "tree": base / "species_tree.nwk",
        "support_tree": base / "species_tree_support.nwk",
        "busco": busco_full_dir,
    }


def generate_query2family_presence_absence() -> None:
    paths = write_query2family_test_data()
    out_dir = TEST_DATA_DIR / "query2family" / "gene_summary"
    out_dir.mkdir(parents=True, exist_ok=True)
    run(
        [
            sys.executable,
            str(SUPPORT_DIR / "gene_family_presence_absence.py"),
            "--mode",
            "query2family",
            "--dir_gene_family",
            str(paths["gene_family_dir"]),
            "--dir_query_gene",
            str(paths["query_dir"]),
            "--species_tree",
            str(paths["tree"]),
            "--out_presence",
            str(out_dir / "query2family_presence_absence.tsv"),
            "--out_copy_number",
            str(out_dir / "query2family_copy_number.tsv"),
            "--out_long",
            str(out_dir / "query2family_presence_absence.long.tsv"),
            "--out_plot_presence",
            str(out_dir / "query2family_presence_absence.plot.tsv"),
            "--out_plot_copy_number",
            str(out_dir / "query2family_copy_number.plot.tsv"),
            "--out_plot_long",
            str(out_dir / "query2family_presence_absence.plot.long.tsv"),
            "--out_selection",
            str(out_dir / "query2family_presence_absence.plot_selection.tsv"),
            "--include_incomplete",
            "0",
            "--max_families",
            "all",
        ]
    )
    run(
        [
            require_tool("Rscript"),
            str(SUPPORT_DIR / "plot_query2family_presence_absence.R"),
            f"--species_tree={paths['tree']}",
            f"--support_tree={paths['support_tree']}",
            f"--busco_table={paths['busco']}",
            f"--long_table={out_dir / 'query2family_presence_absence.plot.long.tsv'}",
            f"--out_pdf={out_dir / 'query2family_presence_absence.pdf'}",
            f"--out_svg={ASSET_DIR / 'query2family-presence-absence.svg'}",
            "--value=presence",
            "--width=7.2",
        ]
    )
    convert_pdf_first_page_to_png(
        out_dir / "query2family_presence_absence.pdf",
        ASSET_DIR / "query2family-presence-absence.png",
        dpi=180,
    )


def generate_single_copy_ortholog_decay() -> None:
    base = TEST_DATA_DIR / "single-copy-ortholog-decay"
    genecount = base / "Orthogroups.GeneCount.tsv"
    selected_genecount = base / "Orthogroups.GeneCount.selected.tsv"
    species = [
        "Amborella_trichopoda",
        "Arabidopsis_thaliana",
        "Cephalotus_follicularis",
        "Dionaea_muscipula",
        "Nepenthes_gracilis",
        "Oryza_sativa",
        "Spinacia_oleracea",
    ]
    rows = [
        ["OG0000001", 1, 1, 1, 1, 1, 1, 1],
        ["OG0000002", 1, 2, 1, 1, 1, 1, 1],
        ["OG0000003", 1, 1, 0, 1, 1, 1, 1],
        ["OG0000004", 1, 1, 1, 0, 0, 1, 1],
        ["OG0000005", 0, 1, 1, 2, 2, 0, 1],
        ["OG0000006", 1, 1, 1, 1, 0, 1, 0],
        ["OG0000007", 0, 0, 1, 1, 1, 1, 1],
        ["OG0000008", 1, 3, 1, 0, 1, 1, 1],
        ["OG0000009", 1, 1, 1, 1, 1, 0, 1],
        ["OG0000010", 1, 1, 2, 2, 2, 1, 1],
        ["OG0000011", 0, 1, 0, 1, 1, 1, 0],
        ["OG0000012", 1, 0, 1, 0, 1, 0, 1],
    ]
    dict_rows = []
    for row in rows:
        counts = dict(zip(species, row[1:], strict=True))
        counts["Orthogroup"] = row[0]
        counts["Total"] = sum(row[1:])
        dict_rows.append(counts)
    fields = ["Orthogroup", *species, "Total"]
    write_tsv(genecount, dict_rows, fields)
    selected_ids = {
        "OG0000001",
        "OG0000003",
        "OG0000005",
        "OG0000007",
        "OG0000009",
        "OG0000011",
    }
    write_tsv(
        selected_genecount,
        [row for row in dict_rows if row["Orthogroup"] in selected_ids],
        fields,
        lineterminator="\n",
    )
    run(
        [
            sys.executable,
            str(SUPPORT_DIR / "single_copy_ortholog_decay_plot.py"),
            "--orthogroup-genecount",
            str(genecount),
            "--selected-orthogroup-genecount",
            str(selected_genecount),
            "--outdir",
            str(base / "out"),
            "--replicates",
            "50",
            "--species-counts",
            "1-7",
            "--seed",
            "7",
            "--plot-basename",
            "single-copy-ortholog-decay",
            "--summary-name",
            "single-copy-ortholog-decay.tsv",
            "--formats",
            "svg",
        ]
    )
    asset_path = ASSET_DIR / "single-copy-ortholog-decay.svg"
    shutil.copyfile(base / "out" / "single-copy-ortholog-decay.svg", asset_path)
    svg_text = asset_path.read_text(encoding="utf-8")
    asset_path.write_text(
        "\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n",
        encoding="utf-8",
    )


def generate_hgt_summary_plots() -> None:
    base = TEST_DATA_DIR / "hgt-summary"
    branch_tsv = base / "hgt_branch_candidates.tsv"
    gene_tsv = base / "hgt_gene_candidates.tsv"
    out_dir = base / "out"
    branch_rows = [
        {
            "orthogroup": "OG0000001",
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
            "clade_min_expression_pearsoncor": 0.42,
            "synteny_support_fraction": 0.5,
            "synteny_mean_support_score": 0.6,
            "contamination_incompatible_fraction": 0.0,
            "contamination_top_lca_sciname": "",
        },
        {
            "orthogroup": "OG0000002",
            "branch_id": 7,
            "candidate_gene_count": 3,
            "matched_leaf_count": 3,
            "besthit_gene_count": 3,
            "besthit_taxid_count": 0,
            "besthit_taxonomy_method": "name_heuristic",
            "besthit_same_superkingdom_fraction": 0.0,
            "besthit_lca_rank_mode": "phylum_mismatch",
            "intron_support_fraction": 0.33,
            "expression_measured_fraction": 0.67,
            "clade_min_expression_pearsoncor": 0.15,
            "synteny_support_fraction": 0.0,
            "synteny_mean_support_score": 0.0,
            "contamination_incompatible_fraction": 0.33,
            "contamination_top_lca_sciname": "Pseudomonas fluorescens",
        },
        {
            "orthogroup": "OG0000003",
            "branch_id": 12,
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
        {
            "orthogroup": "OG0000004",
            "branch_id": 19,
            "candidate_gene_count": 4,
            "matched_leaf_count": 4,
            "besthit_gene_count": 3,
            "besthit_taxid_count": 0,
            "besthit_taxonomy_method": "name_heuristic",
            "besthit_same_superkingdom_fraction": 0.25,
            "besthit_lca_rank_mode": "class_mismatch",
            "intron_support_fraction": 0.75,
            "expression_measured_fraction": 0.75,
            "clade_min_expression_pearsoncor": 0.31,
            "synteny_support_fraction": 0.25,
            "synteny_mean_support_score": 0.22,
            "contamination_incompatible_fraction": 0.0,
            "contamination_top_lca_sciname": "",
        },
    ]
    branch_fields = list(branch_rows[0].keys())
    write_tsv(branch_tsv, branch_rows, branch_fields)

    gene_rows = [
        {
            "orthogroup": "OG0000001",
            "gene_id": "Arabidopsis_thaliana_AHA1",
            "gene_taxon": "Arabidopsis thaliana",
            "candidate_branch_count": 1,
            "candidate_branch_ids": "3",
            "besthit_accession": "P001",
            "besthit_organism": "Bacillus subtilis",
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
            "orthogroup": "OG0000001",
            "gene_id": "Spinacia_oleracea_AHA1",
            "gene_taxon": "Spinacia oleracea",
            "candidate_branch_count": 1,
            "candidate_branch_ids": "3",
            "besthit_accession": "P002",
            "besthit_organism": "Bacillus subtilis",
            "besthit_taxid": "",
            "besthit_taxonomy_method": "name_heuristic",
            "besthit_lca_rank": "genus_mismatch",
            "besthit_same_superkingdom": 0,
            "intron_supported": True,
            "expression_measured": True,
            "synteny_support_score": 0.4,
            "contamination_lca_taxid": "",
            "contamination_lca_sciname": "",
            "contamination_is_compatible_lineage": True,
        },
        {
            "orthogroup": "OG0000002",
            "gene_id": "Cephalotus_follicularis_HGT1",
            "gene_taxon": "Cephalotus follicularis",
            "candidate_branch_count": 1,
            "candidate_branch_ids": "7",
            "besthit_accession": "P003",
            "besthit_organism": "Pseudomonas fluorescens",
            "besthit_taxid": "",
            "besthit_taxonomy_method": "name_heuristic",
            "besthit_lca_rank": "phylum_mismatch",
            "besthit_same_superkingdom": 0,
            "intron_supported": False,
            "expression_measured": True,
            "synteny_support_score": 0.0,
            "contamination_lca_taxid": "294",
            "contamination_lca_sciname": "Pseudomonas fluorescens",
            "contamination_is_compatible_lineage": False,
        },
        {
            "orthogroup": "OG0000003",
            "gene_id": "Dionaea_muscipula_HGT2",
            "gene_taxon": "Dionaea muscipula",
            "candidate_branch_count": 1,
            "candidate_branch_ids": "12",
            "besthit_accession": "P004",
            "besthit_organism": "Escherichia coli",
            "besthit_taxid": "",
            "besthit_taxonomy_method": "name_heuristic",
            "besthit_lca_rank": "species",
            "besthit_same_superkingdom": 1,
            "intron_supported": False,
            "expression_measured": True,
            "synteny_support_score": 0.0,
            "contamination_lca_taxid": "562",
            "contamination_lca_sciname": "Escherichia coli",
            "contamination_is_compatible_lineage": False,
        },
        {
            "orthogroup": "OG0000004",
            "gene_id": "Nepenthes_gracilis_HGT3",
            "gene_taxon": "Nepenthes gracilis",
            "candidate_branch_count": 1,
            "candidate_branch_ids": "19",
            "besthit_accession": "P005",
            "besthit_organism": "Streptomyces coelicolor",
            "besthit_taxid": "",
            "besthit_taxonomy_method": "name_heuristic",
            "besthit_lca_rank": "class_mismatch",
            "besthit_same_superkingdom": 0,
            "intron_supported": True,
            "expression_measured": False,
            "synteny_support_score": 0.2,
            "contamination_lca_taxid": "",
            "contamination_lca_sciname": "",
            "contamination_is_compatible_lineage": True,
        },
    ]
    gene_fields = list(gene_rows[0].keys())
    write_tsv(gene_tsv, gene_rows, gene_fields)

    run(
        [
            sys.executable,
            str(SUPPORT_DIR / "plot_hgt_summary.py"),
            "--branch_tsv",
            str(branch_tsv),
            "--gene_tsv",
            str(gene_tsv),
            "--overview_pdf",
            str(out_dir / "hgt_branch_overview.pdf"),
            "--taxonomy_flow_pdf",
            str(out_dir / "hgt_taxonomy_flow.pdf"),
            "--flow_rank",
            "phylum",
            "--flow_max_categories",
            "8",
        ]
    )
    convert_pdf_first_page_to_png(out_dir / "hgt_branch_overview.pdf", ASSET_DIR / "hgt-branch-overview.png")
    convert_pdf_first_page_to_png(out_dir / "hgt_taxonomy_flow.pdf", ASSET_DIR / "hgt-taxonomy-flow.png")


def copy_quickstart_tree_plot() -> None:
    quickstart_pdf = REPO_ROOT / "workspace" / "output" / "query2family" / "tree_plot" / "AHA_tree_plot.pdf"
    out_png = ASSET_DIR / "query2family-tree-plot.png"
    if not quickstart_pdf.is_file():
        print(f"Skipping quick-start tree plot; source PDF was not found: {quickstart_pdf}")
        return
    convert_pdf_first_page_to_png(quickstart_pdf, out_png, dpi=150)


def main() -> int:
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    generate_query2family_presence_absence()
    generate_single_copy_ortholog_decay()
    generate_hgt_summary_plots()
    copy_quickstart_tree_plot()
    print("Example plots generated under:", ASSET_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

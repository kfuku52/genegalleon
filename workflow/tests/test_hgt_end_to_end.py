from pathlib import Path
import gzip
import os
import shutil
import sqlite3
import subprocess

import pandas
import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
HGT_CORE = REPO_ROOT / "workflow" / "core" / "gg_hgt_core.sh"


def _has_r_package(package: str) -> bool:
    completed = subprocess.run(
        ["Rscript", "-e", f'suppressPackageStartupMessages(library("{package}", quietly=TRUE))'],
        capture_output=True,
        text=True,
        check=False,
    )
    return completed.returncode == 0


def _write_gzip_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def _install_stub_ggimage(tmp_path: Path) -> Path:
    lib_dir = tmp_path / "r_stub_lib"
    src_dir = tmp_path / "ggimage_stub"
    r_dir = src_dir / "R"
    r_dir.mkdir(parents=True, exist_ok=True)
    (src_dir / "DESCRIPTION").write_text(
        "\n".join(
            [
                "Package: ggimage",
                "Type: Package",
                "Title: Stub ggimage for workflow integration tests",
                "Version: 0.0.0.9000",
                "Authors@R: person('Codex', 'Stub', email = 'stub@example.com', role = c('aut', 'cre'))",
                "Description: Minimal stub package so requireNamespace('ggimage') succeeds during integration tests.",
                "License: MIT",
                "Encoding: UTF-8",
                "LazyData: true",
            ]
        ),
        encoding="utf-8",
    )
    (src_dir / "NAMESPACE").write_text('exportPattern("^[[:alpha:]]+")\n', encoding="utf-8")
    (r_dir / "stub.R").write_text("stub_ggimage <- function() invisible(TRUE)\n", encoding="utf-8")
    lib_dir.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        ["R", "CMD", "INSTALL", "-l", str(lib_dir), str(src_dir)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    return lib_dir


def _install_stub_ggmsa(tmp_path: Path) -> Path:
    lib_dir = tmp_path / "r_stub_lib"
    src_dir = tmp_path / "ggmsa_stub"
    r_dir = src_dir / "R"
    r_dir.mkdir(parents=True, exist_ok=True)
    (src_dir / "DESCRIPTION").write_text(
        "\n".join(
            [
                "Package: ggmsa",
                "Type: Package",
                "Title: Stub ggmsa for workflow integration tests",
                "Version: 0.0.0.9000",
                "Authors@R: person('Codex', 'Stub', email = 'stub@example.com', role = c('aut', 'cre'))",
                "Description: Minimal stub package providing geom_msa so alignment panels render during integration tests.",
                "License: MIT",
                "Encoding: UTF-8",
                "Imports: ggplot2",
                "LazyData: true",
            ]
        ),
        encoding="utf-8",
    )
    (src_dir / "NAMESPACE").write_text("export(geom_msa)\nimportFrom(ggplot2, geom_blank)\n", encoding="utf-8")
    (r_dir / "stub.R").write_text(
        "\n".join(
            [
                "geom_msa <- function(data = NULL, color = NULL, ...) {",
                "  ggplot2::geom_blank(data = data, ...)",
                "}",
            ]
        ),
        encoding="utf-8",
    )
    lib_dir.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        ["R", "CMD", "INSTALL", "-l", str(lib_dir), str(src_dir)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    return lib_dir


def _install_fake_conda(tmp_path: Path) -> Path:
    bin_dir = tmp_path / "fake_conda_bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    script_path = bin_dir / "conda"
    script_path.write_text(
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "if [[ \"$1\" == \"shell.bash\" && \"$2\" == \"hook\" ]]; then",
                "  cat <<'EOF'",
                "conda() {",
                "  case \"$1\" in",
                "    activate|deactivate)",
                "      return 0",
                "      ;;",
                "    *)",
                "      return 0",
                "      ;;",
                "  esac",
                "}",
                "EOF",
                "  exit 0",
                "fi",
                "exit 0",
            ]
        ),
        encoding="utf-8",
    )
    script_path.chmod(0o755)
    return bin_dir


def _write_branch_db(db_path: Path) -> None:
    branch = pandas.DataFrame(
        [
            {
                "orthogroup": "OG0001",
                "branch_id": 3,
                "node_name": "n3",
                "gene_labels": "Arabidopsis_thaliana_geneA; Arabidopsis_thaliana_geneB",
                "num_leaf": 2,
                "so_event": "S",
                "taxon": "",
                "spnode_coverage": "Arabidopsis_thaliana",
                "generax_event": "H",
                "generax_transfer": "Y",
                "generax_event_parent": "S",
                "clade_min_expression_pearsoncor": 0.75,
                "num_intron": pandas.NA,
                "intron_present": pandas.NA,
                "synteny_support_score": pandas.NA,
                "sprot_best": pandas.NA,
                "organism": pandas.NA,
                "taxid_y": pandas.NA,
                "expression_cond1": pandas.NA,
                "expression_cond2": pandas.NA,
            },
            {
                "orthogroup": "OG0001",
                "branch_id": 1,
                "node_name": "Arabidopsis_thaliana_geneA",
                "gene_labels": "Arabidopsis_thaliana_geneA",
                "num_leaf": 1,
                "so_event": "L",
                "taxon": "Arabidopsis thaliana",
                "spnode_coverage": "Arabidopsis_thaliana",
                "generax_event": "L",
                "generax_transfer": "",
                "generax_event_parent": "H",
                "clade_min_expression_pearsoncor": pandas.NA,
                "num_intron": 2,
                "intron_present": 1,
                "synteny_support_score": 0.8,
                "sprot_best": "P00001",
                "organism": "Escherichia coli str. K-12",
                "taxid_y": pandas.NA,
                "expression_cond1": 12.0,
                "expression_cond2": 9.0,
            },
            {
                "orthogroup": "OG0001",
                "branch_id": 2,
                "node_name": "Arabidopsis_thaliana_geneB",
                "gene_labels": "Arabidopsis_thaliana_geneB",
                "num_leaf": 1,
                "so_event": "L",
                "taxon": "Arabidopsis thaliana",
                "spnode_coverage": "Arabidopsis_thaliana",
                "generax_event": "L",
                "generax_transfer": "",
                "generax_event_parent": "H",
                "clade_min_expression_pearsoncor": pandas.NA,
                "num_intron": 1,
                "intron_present": 1,
                "synteny_support_score": 0.6,
                "sprot_best": "P00002",
                "organism": "Bacillus subtilis subsp. subtilis",
                "taxid_y": pandas.NA,
                "expression_cond1": 8.0,
                "expression_cond2": 7.0,
            },
        ]
    )
    db_path.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(db_path) as conn:
        branch.to_sql("branch", conn, index=False)


def _write_stat_branch(path: Path) -> None:
    df = pandas.DataFrame(
        [
            {
                "orthogroup": "OG0001",
                "branch_id": 3,
                "parent": -999,
                "child1": 1,
                "child2": 2,
                "node_name": "n3",
                "gene_labels": "Arabidopsis_thaliana_geneA; Arabidopsis_thaliana_geneB",
                "num_leaf": 2,
                "so_event": "S",
                "taxon": "",
                "spnode_coverage": "Arabidopsis_thaliana",
                "generax_event": "H",
                "generax_transfer": "Y",
                "generax_event_parent": "S",
                "bl_rooted": 0.25,
                "support_unrooted": 100,
                "clade_min_expression_pearsoncor": 0.75,
                "expression_cond1": pandas.NA,
                "expression_cond2": pandas.NA,
                "start": pandas.NA,
                "end": pandas.NA,
                "chromosome": "",
                "tmhmm_predhel": pandas.NA,
                "num_intron": pandas.NA,
                "intron_positions": "",
                "targetp_noTP": pandas.NA,
                "targetp_SP": pandas.NA,
                "targetp_mTP": pandas.NA,
                "targetp_cTP": pandas.NA,
                "targetp_luTP": pandas.NA,
                "fimo_start": "",
                "fimo_end": "",
                "fimo_strand": "",
                "fimo_alt_id": "",
                "fimo_qvalue": "",
                "promoter_available": "",
                "promoter_N": "",
            },
            {
                "orthogroup": "OG0001",
                "branch_id": 1,
                "parent": 3,
                "child1": pandas.NA,
                "child2": pandas.NA,
                "node_name": "Arabidopsis_thaliana_geneA",
                "gene_labels": "Arabidopsis_thaliana_geneA",
                "num_leaf": 1,
                "so_event": "L",
                "taxon": "Arabidopsis thaliana",
                "spnode_coverage": "Arabidopsis_thaliana",
                "generax_event": "L",
                "generax_transfer": "",
                "generax_event_parent": "H",
                "bl_rooted": 0.10,
                "support_unrooted": 95,
                "clade_min_expression_pearsoncor": pandas.NA,
                "expression_cond1": 12.0,
                "expression_cond2": 9.0,
                "start": 100,
                "end": 220,
                "chromosome": "chr1",
                "tmhmm_predhel": 1,
                "num_intron": 2,
                "intron_positions": "3;9",
                "targetp_noTP": 0.10,
                "targetp_SP": 0.70,
                "targetp_mTP": 0.10,
                "targetp_cTP": 0.05,
                "targetp_luTP": 0.05,
                "fimo_start": "10;40",
                "fimo_end": "14;44",
                "fimo_strand": "+;+",
                "fimo_alt_id": "M1;M2",
                "fimo_qvalue": "0.001;0.003",
                "promoter_available": "Y",
                "promoter_N": "",
            },
            {
                "orthogroup": "OG0001",
                "branch_id": 2,
                "parent": 3,
                "child1": pandas.NA,
                "child2": pandas.NA,
                "node_name": "Arabidopsis_thaliana_geneB",
                "gene_labels": "Arabidopsis_thaliana_geneB",
                "num_leaf": 1,
                "so_event": "L",
                "taxon": "Arabidopsis thaliana",
                "spnode_coverage": "Arabidopsis_thaliana",
                "generax_event": "L",
                "generax_transfer": "",
                "generax_event_parent": "H",
                "bl_rooted": 0.12,
                "support_unrooted": 90,
                "clade_min_expression_pearsoncor": pandas.NA,
                "expression_cond1": 8.0,
                "expression_cond2": 7.0,
                "start": 260,
                "end": 360,
                "chromosome": "chr1",
                "tmhmm_predhel": 0,
                "num_intron": 1,
                "intron_positions": "6",
                "targetp_noTP": 0.75,
                "targetp_SP": 0.05,
                "targetp_mTP": 0.10,
                "targetp_cTP": 0.05,
                "targetp_luTP": 0.05,
                "fimo_start": "16",
                "fimo_end": "20",
                "fimo_strand": "-",
                "fimo_alt_id": "M1",
                "fimo_qvalue": "0.002",
                "promoter_available": "Y",
                "promoter_N": "",
            },
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def _write_workspace_fixture(workspace: Path) -> None:
    output_root = workspace / "output" / "orthogroup"
    _write_branch_db(output_root / "gg_orthogroup.db")
    _write_stat_branch(output_root / "stat_branch" / "OG0001_stat.branch.tsv")

    synteny_df = pandas.DataFrame(
        [
            {"node_name": "Arabidopsis_thaliana_geneA", "offset": -1, "group_id": "grpL"},
            {"node_name": "Arabidopsis_thaliana_geneB", "offset": -1, "group_id": "grpL"},
            {"node_name": "Arabidopsis_thaliana_geneA", "offset": 1, "group_id": "grpR"},
            {"node_name": "Arabidopsis_thaliana_geneB", "offset": 1, "group_id": "grpR"},
        ]
    )
    (output_root / "synteny").mkdir(parents=True, exist_ok=True)
    synteny_df.to_csv(output_root / "synteny" / "OG0001_synteny.tsv", sep="\t", index=False)

    rps_df = pandas.DataFrame(
        [
            {
                "qacc": "Arabidopsis_thaliana_geneA",
                "sacc": "PF00001",
                "stitle": "PF00001,DomainA",
                "qlen": 4,
                "slen": 50,
                "qstart": 1,
                "qend": 3,
            },
            {
                "qacc": "Arabidopsis_thaliana_geneB",
                "sacc": "PF00002",
                "stitle": "PF00002,DomainB",
                "qlen": 4,
                "slen": 40,
                "qstart": 2,
                "qend": 4,
            },
        ]
    )
    (output_root / "rpsblast").mkdir(parents=True, exist_ok=True)
    rps_df.to_csv(output_root / "rpsblast" / "OG0001_rpsblast.tsv", sep="\t", index=False)

    fasta_text = (
        ">Arabidopsis_thaliana_geneA\nATGGCCATGGCC\n"
        ">Arabidopsis_thaliana_geneB\nATGAACATGAAC\n"
    )
    _write_gzip_text(output_root / "clipkit" / "OG0001_cds.clipkit.fa.gz", fasta_text)
    _write_gzip_text(output_root / "orthogroup_extraction_fasta" / "OG0001_orthogroup_extraction.fa.gz", fasta_text)
    _write_gzip_text(output_root / "maxalign" / "OG0001_cds.maxalign.fa.gz", fasta_text)
    _write_gzip_text(output_root / "mafft" / "OG0001_cds.aln.fa.gz", fasta_text)
    _write_gzip_text(output_root / "protein_fasta" / "OG0001_pep.fa.gz", fasta_text)
    _write_gzip_text(output_root / "cds_fasta" / "OG0001_cds.fa.gz", fasta_text)

    (output_root / "dated_tree").mkdir(parents=True, exist_ok=True)
    (output_root / "dated_tree" / "OG0001_dated.nwk").write_text(
        "(Arabidopsis_thaliana_geneA:0.1,Arabidopsis_thaliana_geneB:0.1);\n",
        encoding="utf-8",
    )

    (output_root / "meme").mkdir(parents=True, exist_ok=True)
    (output_root / "meme" / "OG0001_meme.xml").write_text(
        "\n".join(
            [
                "<MEME>",
                "  <motifs>",
                "    <motif id=\"motif_1\" name=\"MotifAlpha\" alt=\"TF_ALPHA\" width=\"10\" sites=\"5\" e_value=\"1e-6\"/>",
                "    <motif id=\"motif_2\" name=\"MotifBeta\" width=\"8\" sites=\"3\" e_value=\"2e-4\"/>",
                "  </motifs>",
                "</MEME>",
            ]
        ),
        encoding="utf-8",
    )


@pytest.mark.skipif(shutil.which("seqkit") is None, reason="seqkit is required for gg_hgt_core end-to-end validation")
def test_hgt_core_end_to_end_generates_tables_and_pdfs(tmp_path: Path):
    if not _has_r_package("rkftools"):
        pytest.skip("rkftools is required for ortholog panel rendering")

    workspace = tmp_path / "workspace"
    _write_workspace_fixture(workspace)
    stub_lib = _install_stub_ggimage(tmp_path)
    _install_stub_ggmsa(tmp_path)
    fake_conda_bin = _install_fake_conda(tmp_path)

    env = {
        key: value
        for key, value in os.environ.items()
        if not (key.startswith("CONDA") or key.startswith("MAMBA"))
    }
    env["gg_workspace_dir"] = str(workspace)
    env["run_hgt_eval"] = "1"
    env["run_hgt_plot"] = "1"
    env["hgt_use_taxonomy_db"] = "0"
    env["MPLCONFIGDIR"] = str(tmp_path / "mplconfig")
    env["PATH"] = f"{fake_conda_bin}{os.pathsep}{env['PATH']}"
    env["R_LIBS_USER"] = (
        f"{stub_lib}{os.pathsep}{env['R_LIBS_USER']}"
        if env.get("R_LIBS_USER")
        else str(stub_lib)
    )

    completed = subprocess.run(
        ["bash", str(HGT_CORE)],
        cwd=str(REPO_ROOT),
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr

    hgt_root = workspace / "output" / "hgt"
    branch_tsv = hgt_root / "hgt_branch_candidates.tsv"
    gene_tsv = hgt_root / "hgt_gene_candidates.tsv"
    orthogroup_tsv = hgt_root / "hgt_orthogroup_summary.tsv"
    overview_pdf = hgt_root / "plots" / "hgt_branch_overview.pdf"
    taxonomy_flow_pdf = hgt_root / "plots" / "hgt_taxonomy_flow.pdf"
    tree_plot_pdf = hgt_root / "tree_plot" / "OG0001_hgt_tree_plot.pdf"

    for path in [branch_tsv, gene_tsv, orthogroup_tsv, overview_pdf, taxonomy_flow_pdf, tree_plot_pdf]:
        assert path.exists(), f"Expected output not found: {path}"
        assert path.stat().st_size > 0, f"Output was empty: {path}"

    branch_df = pandas.read_csv(branch_tsv, sep="\t")
    gene_df = pandas.read_csv(gene_tsv, sep="\t")
    orthogroup_df = pandas.read_csv(orthogroup_tsv, sep="\t")

    assert branch_df.shape[0] == 1
    assert branch_df.loc[0, "orthogroup"] == "OG0001"
    assert branch_df.loc[0, "candidate_gene_count"] == 2
    assert gene_df["gene_id"].tolist() == ["Arabidopsis_thaliana_geneA", "Arabidopsis_thaliana_geneB"]
    assert orthogroup_df.loc[0, "hgt_branch_count"] == 1

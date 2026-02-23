from pathlib import Path
import subprocess
import sys

import pandas


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "parse_grampa.py"


def test_parse_grampa_writes_summary_with_species_gene_columns(tmp_path):
    grampa_det = tmp_path / "grampa.det.tsv"
    grampa_out = tmp_path / "grampa.out.tsv"
    gene_trees = tmp_path / "gene_trees.nwk"
    sorted_names = tmp_path / "sorted_gene_tree_file_names.tsv"

    grampa_det.write_text(
        "# GT/MT combo\tdups\tlosses\tTotal score\tMaps\n"
        "* GT-1 to MT-1\t2\t3\t5\t1\n"
    )
    grampa_out.write_text("MT-1\tH1\tH2\t(sp1_sp1,sp2_sp2);\t7\n")
    gene_trees.write_text("(geneA_sp1_sp1,geneB_sp2_sp2);\n")
    sorted_names.write_text("gene_tree_1.nwk\n")

    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--grampa_det",
            str(grampa_det),
            "--grampa_out",
            str(grampa_out),
            "--gene_trees",
            str(gene_trees),
            "--species_tree",
            "(sp1_sp1,sp2_sp2);",
            "--sorted_gene_tree_file_names",
            str(sorted_names),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr

    out = pandas.read_csv(tmp_path / "grampa_summary.tsv", sep="\t")
    assert out.shape[0] == 1
    assert out.loc[0, "gene_tree"] == "GT-1"
    assert out.loc[0, "mul_tree"] == "MT-1"
    assert out.loc[0, "sp1_sp1"] == "geneA"
    assert out.loc[0, "sp2_sp2"] == "geneB"


def test_parse_grampa_writes_placeholder_when_no_maps(tmp_path):
    grampa_det = tmp_path / "grampa.det.tsv"
    grampa_out = tmp_path / "grampa.out.tsv"
    gene_trees = tmp_path / "gene_trees.nwk"
    sorted_names = tmp_path / "sorted_gene_tree_file_names.tsv"

    grampa_det.write_text(
        "# GT/MT combo\tdups\tlosses\tTotal score\tMaps\n"
        "GT-1 to MT-1\tNo maps found!\t3\t5\t0\n"
    )
    grampa_out.write_text("")
    gene_trees.write_text("(geneA_sp1_sp1,geneB_sp2_sp2);\n")
    sorted_names.write_text("gene_tree_1.nwk\n")

    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--grampa_det",
            str(grampa_det),
            "--grampa_out",
            str(grampa_out),
            "--gene_trees",
            str(gene_trees),
            "--species_tree",
            "(sp1_sp1,sp2_sp2);",
            "--sorted_gene_tree_file_names",
            str(sorted_names),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    placeholder = (tmp_path / "grampa_summary.tsv").read_text()
    assert "placeholder" in placeholder

import subprocess
import sys
from pathlib import Path

import pandas
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "parse_grampa.py"


@pytest.mark.parametrize("ncpu", [1, 2])
@pytest.mark.parametrize("format_name", ["legacy", "modern"])
def test_grampa_repeated_maps_keep_order_and_missing_tree_blanks(tmp_path, ncpu, format_name):
    det = tmp_path / "det.tsv"
    out = tmp_path / "out.tsv"
    trees = tmp_path / "trees.nwk"
    names = tmp_path / "names.tsv"
    species_tree = "(Genus_speciesB,Genus_speciesA);"
    if format_name == "modern":
        det.write_text(
            "mul.tree\tgene.tree\tdups\tlosses\ttotal.score\tmaps\n"
            "2\t2\t1\t2\t3\t1\n1\t1\t2\t3\t5\t1\n1\t2\t1\t2\t3\t1\n1\t3\t0\t0\t0\t1\n")
        out.write_text(
            "mul.tree\th1.node\th2.node\tscore\tlabeled.tree\n"
            f"2\tH3\tH4\t7\t{species_tree}\n1\tH1\tH2\t5\t{species_tree}\n")
    else:
        det.write_text(
            "# GT/MT combo\tdups\tlosses\tTotal score\tMaps\n"
            "* GT-2 to MT-2\t1\t2\t3\t1\n* GT-1 to MT-1\t2\t3\t5\t1\n"
            "* GT-2 to MT-1\t1\t2\t3\t1\n* GT-3 to MT-1\t0\t0\t0\t1\n")
        out.write_text(f"MT-2\tH3\tH4\t{species_tree}\t7\nMT-1\tH1\tH2\t{species_tree}\t5\n")
    trees.write_text(
        "(Genus_speciesB_b1,Genus_speciesA_a2,a1_Genus_speciesA,Unknown_species_x);\n"
        "(b2_Genus_speciesB,Genus_speciesB_b3,a3_Genus_speciesA);\n")
    names.write_text("first.nwk\nsecond.nwk\n")
    original = [path.read_bytes() for path in (det, out, trees, names)]
    completed = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "--grampa_det", str(det), "--grampa_out", str(out),
         "--gene_trees", str(trees), "--species_tree", species_tree,
         "--sorted_gene_tree_file_names", str(names), "--ncpu", str(ncpu)],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert completed.returncode == 0, completed.stderr
    assert completed.stderr == ""
    result = pandas.read_csv(tmp_path / "grampa_summary.tsv", sep="\t", keep_default_na=False)
    assert result["gene_tree"].tolist() == ["GT-2", "GT-1", "GT-2", "GT-3"]
    assert result["mul_tree"].tolist() == ["MT-2", "MT-1", "MT-1", "MT-1"]
    assert result["Genus_speciesA"].tolist() == ["a3", "a2,a1", "a3", ""]
    assert result["Genus_speciesB"].tolist() == ["b2,b3", "b1", "b2,b3", ""]
    assert result["file_name"].tolist() == ["second.nwk", "first.nwk", "second.nwk", ""]
    assert result["multree_score"].tolist() == [7, 5, 5, 5]
    assert result.columns.tolist() == [
        "file_name", "gene_tree", "mul_tree", "dups", "losses", "total_score", "maps",
        "Genus_speciesA", "Genus_speciesB", "H1_node", "H2_node", "mul_tree_string", "multree_score",
    ]
    assert [path.read_bytes() for path in (det, out, trees, names)] == original


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

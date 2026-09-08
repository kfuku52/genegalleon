import subprocess
import sys
from pathlib import Path

import pandas
import pytest

from workflow.support.fasta_sequence_store import build_store
from workflow.support.gff_source_contract import paired_gff_name, source_bound_gff_names


def indexed_source(tmp_path):
    source = tmp_path / "Species_a_TAIR10.cds.all.fa"
    source.write_text(">Species_a_gene1\nATGAAA\n", encoding="utf-8")
    database = tmp_path / "sources.sqlite3"
    build_store(database, tmp_path / "sources.json", [(source, "Species_a")], 0, 0)
    return source, database


def test_source_pair_is_exact_not_longest_prefix():
    names = ["Species_a_TAIR10.62.gff.gz", "Species_a_Araport11.gene.gff3", "Species_a_TAIR100.62.gff.gz"]
    assert paired_gff_name("Species_a_TAIR10.cds.all.fa.gz", names) == "Species_a_TAIR10.62.gff.gz"
    with pytest.raises(ValueError, match="Ambiguous GFF source pair"):
        paired_gff_name("Species_a_TAIR10.cds.all.fa.gz", names + ["Species_a_TAIR10.61.gff.gz"])
    with pytest.raises(ValueError, match="Cannot determine"):
        paired_gff_name("Species_a.unrecognized.fa", names)


def test_source_binding_rejects_changed_sequence_or_input(tmp_path):
    source, database = indexed_source(tmp_path)
    candidates = {"Species_a_gene1": ["Species_a_TAIR10.62.gff", "Species_a_other.gff"]}
    assert source_bound_gff_names(database, {"Species_a_gene1": "ATGAAA"}, candidates) == {
        "Species_a_gene1": "Species_a_TAIR10.62.gff"
    }
    with pytest.raises(ValueError, match="identity/sequence mismatch"):
        source_bound_gff_names(database, {"Species_a_gene1": "ATGTTT"}, candidates)
    source.write_text(">Species_a_gene1\nATGAAAAAA\n", encoding="utf-8")
    with pytest.raises(ValueError, match="FASTA changed"):
        source_bound_gff_names(database, {"Species_a_gene1": "ATGAAA"}, candidates)


@pytest.mark.parametrize("use_store", [True, False])
def test_gff_cli_never_merges_conflicting_source_annotations(tmp_path, use_store):
    source, database = indexed_source(tmp_path)
    gff = tmp_path / "gff"
    gff.mkdir()
    (gff / "Species_a_TAIR10.62.gff").write_text(
        "chr1\tsrc\tCDS\t101\t106\t.\t+\t0\tParent=gene1;\n", encoding="utf-8"
    )
    (gff / "Species_a_Araport11.gene.gff3").write_text(
        "region\tsrc\tCDS\t1\t3\t.\t+\t0\tParent=gene1;\n", encoding="utf-8"
    )
    output = tmp_path / "stats.tsv"
    script = Path(__file__).resolve().parents[1] / "support/gff2genestat.py"
    command = [sys.executable, str(script), "--dir_gff", str(gff), "--seqfile", str(source),
               "--outfile", str(output), "--ncpu", "2"]
    if use_store:
        command += ["--sequence-store", str(database)]
    result = subprocess.run(command, text=True, capture_output=True)
    if use_store:
        assert result.returncode == 0, result.stdout + result.stderr
        rows = pandas.read_csv(output, sep="\t")
        assert len(rows) == 1
        assert rows.iloc[0]["feature_size"] == 6
        assert rows.iloc[0]["chromosome"] == "chr1"
    else:
        assert result.returncode != 0
        assert "Ambiguous GFF source for genes" in result.stderr
        assert not output.exists()

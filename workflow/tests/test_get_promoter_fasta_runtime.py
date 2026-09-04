import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "get_promoter_fasta.py"


@pytest.mark.runtime
@pytest.mark.parametrize("chromosome", ["1", "01", "NA"])
@pytest.mark.parametrize("strand", ["+", "-"])
def test_promoter_cli_preserves_chromosome_identifiers(tmp_path, chromosome, strand):
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "Species_a.fa").write_text(f">{chromosome}\n" + "AACG" * 10 + "\n")
    info = tmp_path / "genes.tsv"
    info.write_text(
        "gene_id\tchromosome\tstart\tend\tstrand\n"
        f"Species_a_gene1\t{chromosome}\t21\t28\t{strand}\n"
    )
    output = tmp_path / "promoter.fa"
    result = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "--dir_genome", str(genomes), "--geneinfo_tsv", str(info),
         "--promoter_bp", "8", "--outfile", str(output)], capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    sequence = "AACGAACG" if strand == "+" else "CGTTCGTT"
    assert output.read_text() == f">Species_a_gene1\n{sequence}\n"

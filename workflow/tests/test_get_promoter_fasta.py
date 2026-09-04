import subprocess
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "get_promoter_fasta.py"


def load_module():
    spec = spec_from_file_location("get_promoter_fasta", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_normalize_ncpu_clamps_non_positive_values():
    mod = load_module()
    assert mod.normalize_ncpu(0) == 1
    assert mod.normalize_ncpu(-5) == 1
    assert mod.normalize_ncpu(3) == 3


def test_normalize_ncpu_rejects_invalid_values():
    mod = load_module()
    with pytest.raises(ValueError):
        mod.normalize_ncpu("abc")


def test_get_genome_file_ignores_hidden_entries_and_directories(tmp_path: Path):
    mod = load_module()
    (tmp_path / ".Species_a.fa").write_text(">a\nAT\n", encoding="utf-8")
    (tmp_path / "Species_a.extra").mkdir()
    target = tmp_path / "Species_a.fa"
    target.write_text(">a\nAT\n", encoding="utf-8")

    detected = mod.get_genome_file(str(tmp_path), "Species_a")
    assert detected == str(target)


def test_get_genome_file_does_not_match_longer_species_label(tmp_path: Path):
    mod = load_module()
    (tmp_path / "Species_a_subsp_x.genome.fa").write_text(">x\nAT\n", encoding="utf-8")

    assert mod.get_genome_file(str(tmp_path), "Species_a") is None


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

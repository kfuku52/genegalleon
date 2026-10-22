import gzip
import importlib.util
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "support" / "fasta_sequence_contract.py"
SPEC = importlib.util.spec_from_file_location("fasta_sequence_contract", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_codon_contract_accepts_plain_and_gzip_alignments(tmp_path):
    plain = tmp_path / "valid.fa"
    plain.write_text(">a\nATG---AAC\n>b\nATGCCCAAC\n", encoding="utf-8")
    summary = MODULE.validate_fasta(plain, "codon")
    assert summary.records == 2
    assert summary.aligned_length == 9

    compressed = tmp_path / "valid.fa.gz"
    with gzip.open(compressed, "wt", encoding="utf-8") as handle:
        handle.write(plain.read_text(encoding="utf-8"))
    assert MODULE.validate_fasta(compressed, "codon") == summary


def test_codon_contract_rejects_protein_artifact(tmp_path):
    alignment = tmp_path / "legacy.fa"
    alignment.write_text(">a\nMPEPTIDE\n>b\nMPEPTIDQ\n", encoding="utf-8")
    with pytest.raises(MODULE.SequenceContractError, match="outside the codon alphabet"):
        MODULE.validate_fasta(alignment, "codon")


def test_codon_contract_counts_mixed_gap_symbols_after_normalization(tmp_path):
    path = tmp_path / "mixed.fa"
    path.write_text(">description retained\n au g - ? . c c g \n>b\nATG\nCCG\n-?.\n")
    assert MODULE.validate_fasta(path, "codon") == MODULE.FastaSummary(2, 18, 9)
    path.write_text(">bad\nA--??.CCG\n")
    with pytest.raises(MODULE.SequenceContractError, match="ungapped sequence lengths"):
        MODULE.validate_fasta(path, "codon")


def test_contract_consumes_gzip_trailer_before_reporting_invalid_alphabet(tmp_path):
    path = tmp_path / "truncated.gz"
    path.write_bytes(gzip.compress(b">bad\nATG!\n")[:-4])
    with pytest.raises(EOFError):
        MODULE.validate_fasta(path, "dna")


@pytest.mark.parametrize(
    "body, message",
    [
        (">a\nATGA\n>b\nATGA\n", "aligned length is not a multiple"),
        (">a\nATG---\n>b\nATGAA-\n", "ungapped sequence lengths"),
        ("ATG\n", "before the first header"),
    ],
)
def test_codon_contract_rejects_structurally_invalid_fasta(tmp_path, body, message):
    alignment = tmp_path / "invalid.fa"
    alignment.write_text(body, encoding="utf-8")
    with pytest.raises(MODULE.SequenceContractError, match=message):
        MODULE.validate_fasta(alignment, "codon")


def test_gene_evolution_enforces_legacy_mafft_and_generax_contracts():
    core = (Path(__file__).resolve().parents[1] / "core" / "gg_gene_evolution_core.sh").read_text(encoding="utf-8")
    assert 'mafft_provenance_args+=(--output-fasta-type "mafft=${mafft_output_fasta_type}")' in core
    assert 'python "${gg_support_dir}/fasta_sequence_contract.py"' in core
    assert '--expected "${generax_input_fasta_type}"' in core


def test_gene_evolution_serializes_only_mafft_iterative_refinement():
    core = (Path(__file__).resolve().parents[1] / "core" / "gg_gene_evolution_core.sh").read_text(encoding="utf-8")
    start = core.index('task="In-frame mafft alignment"')
    end = core.index('task="cdskit stats for original alignment"', start)
    mafft_block = core[start:end]

    assert "mafft_threadit=0" in mafft_block
    assert mafft_block.count('--thread "${GG_TASK_CPUS}"') == 2
    assert mafft_block.count('--threadit "${mafft_threadit}"') == 2
    assert '--parameter "mafft_threads=${GG_TASK_CPUS}"' in mafft_block
    assert '--parameter "mafft_threadit=${mafft_threadit}"' in mafft_block

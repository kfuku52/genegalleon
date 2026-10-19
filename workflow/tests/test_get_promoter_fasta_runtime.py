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


@pytest.mark.parametrize('compressed', [False, True])
@pytest.mark.parametrize('duplicate_contig', [False, True])
def test_batched_promoters_match_legacy_bytes_and_use_custom_seqkit_once(tmp_path, compressed, duplicate_contig):
    import gzip
    import shlex
    from types import SimpleNamespace

    import pandas as pd

    from workflow.support import get_promoter_fasta as promoter

    genome = tmp_path / ('genome.fa.gz' if compressed else 'genome.fa')
    text = '>01 description\n' + 'ACGTacgtnNRYSWKMBDHV'*6 + '\n>chr.2\n' + 'GATTACA'*40 + '\n'
    if duplicate_contig:
        text += '>01 duplicate description\n' + 'TTAACCGG'*15 + '\n'
    if compressed:
        with gzip.open(genome, 'wt') as handle:
            handle.write(text)
    else:
        genome.write_text(text)
    spans = [('01', 1, 12, '+'), ('01', 2, 12, '+'), ('01', 30, 40, '+'),
             ('01', 30, 40, '-'), ('01', 111, 114, '-'), ('01', 1, 12, '?'),
             ('missing', 30, 40, '+'), ('chr.2', 35, 50, '-'), ('01', 1, 200, '-')]
    frame = pd.DataFrame([dict(gene_id=f'Species_a_g{i}', chromosome=c, start=a, end=b, strand=s)
                          for i, (c, a, b, s) in enumerate(spans)])
    args = SimpleNamespace(seqkit_exe='seqkit', ncpu=1, promoter_bp=18)
    expected = promoter.sanitize_fasta([promoter.run_seqkit(args, frame, i, str(genome)) for i in frame.index])
    calls = tmp_path / 'calls'
    executable = tmp_path / 'custom seqkit'
    executable.write_text('#!/bin/sh\nprintf "%s\\n" "$1" >> ' + shlex.quote(str(calls)) + '\nexec seqkit "$@"\n')
    executable.chmod(0o755)
    args.seqkit_exe = str(executable)
    actual = promoter.sanitize_fasta(promoter.run_seqkit_batch(args, frame, frame.index, str(genome)))
    assert actual == expected
    assert calls.read_text() == 'subseq\n'
    assert sorted(p.name for p in tmp_path.iterdir() if p.name.startswith('genome.')) == [genome.name]


def test_batched_promoters_reject_source_replacement_before_accepting_output(tmp_path, monkeypatch):
    from types import SimpleNamespace

    import pandas as pd

    from workflow.support import get_promoter_fasta as promoter

    genome = tmp_path / 'genome.fa'
    genome.write_text('>chr1\nACGTACGT\n')
    frame = pd.DataFrame([dict(gene_id='Species_a_g', chromosome='chr1', start=4, end=6, strand='+')])
    original = promoter.subprocess.run
    def replace(*args, **kwargs):
        result = original(*args, **kwargs)
        replacement = tmp_path / 'replacement'
        replacement.write_bytes(genome.read_bytes())
        replacement.replace(genome)
        return result
    monkeypatch.setattr(promoter.subprocess, 'run', replace)
    with pytest.raises(OSError, match='changed'):
        promoter.run_seqkit_batch(SimpleNamespace(seqkit_exe='seqkit', ncpu=1, promoter_bp=3), frame, frame.index, str(genome))

"""Exercise the actual shell output stage with the runtime's real seqkit."""
import os
import shutil
import subprocess
from pathlib import Path

import pytest

CORE = Path(__file__).resolve().parents[1] / "core/gg_gene_evolution_core.sh"
FIELDS = "qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle"


def completion_script():
    source = CORE.read_text()
    start = source.index("  # Complete query coverage")
    end = source.index('\n  cp_out "${og_id}.rpsblast.tsv"', start)
    return 'set -euo pipefail\nog_id=OG\noutfmt="' + FIELDS + '"\n' + source[start:end]


@pytest.mark.parametrize("hits", [[], ["Species_g1"], ["Species_g10"], ["Species_g1", "Species_g10", "NA", "nan"]])
@pytest.mark.parametrize("wrapped", [False, True])
def test_rpsblast_completion_preserves_hits_query_order_lengths_and_messages(tmp_path, hits, wrapped):
    queries = [("Species_g1", "ACDEF"), ("Species_g10", "MKK"), ("NA", "WW"), ("nan", "C")]
    sequences = [f">{name}\n" + ("\n".join(sequence) if wrapped else sequence) + "\n" for name, sequence in queries]
    (tmp_path / "ungapped_translated_cds.fas").write_text("".join(sequences))
    # Keep multiple domains and original hit order, even when it differs from query order.
    rows = ["\t".join([name, "pfam", "80", "2", "0", "0", "1", "2", "1", "2", "1e-8", "42", "5", "8", "domain"]) + "\n" for name in reversed(hits)]
    if hits:
        rows.append(rows[0].replace("domain\n", "second domain\n"))
    original = "".join(rows)
    (tmp_path / "OG.rpsblast.tmp.tsv").write_text(original)
    result = subprocess.run(["bash", "-c", completion_script()], cwd=tmp_path, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    expected = FIELDS.replace(" ", "\t") + "\n" + original
    messages = []
    for name, sequence in queries:
        if name in hits:
            messages.append(f"{name}: RPS-BLAST hit found.\n")
        else:
            expected += "\t".join([name, *([""] * 11), str(len(sequence)), "", ""]) + "\n"
            messages.append(f"{name}: no hit in RPS-BLAST. Appending qlen to output tsv.\n")
    assert (tmp_path / "OG.rpsblast.tsv").read_text() == expected
    assert result.stdout == "".join(messages)


def test_rpsblast_completion_converts_the_full_fasta_once(tmp_path):
    real_seqkit = shutil.which("seqkit")
    assert real_seqkit, "This workflow check requires the GeneGalleon runtime's seqkit"
    binary = tmp_path / "bin"
    binary.mkdir()
    calls = tmp_path / "calls"
    wrapper = binary / "seqkit"
    wrapper.write_text('#!/bin/sh\nprintf "call\\n" >> "$SEQKIT_CALLS"\nexec "$REAL_SEQKIT" "$@"\n')
    wrapper.chmod(0o755)
    (tmp_path / "ungapped_translated_cds.fas").write_text("".join(f">Species_g{i}\nMACDE\n" for i in range(1000)))
    (tmp_path / "OG.rpsblast.tmp.tsv").touch()
    env = dict(os.environ, PATH=str(binary) + os.pathsep + os.environ["PATH"], REAL_SEQKIT=real_seqkit, SEQKIT_CALLS=str(calls))
    result = subprocess.run(["bash", "-c", completion_script()], cwd=tmp_path, env=env, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    assert calls.read_text() == "call\n"
    assert len((tmp_path / "OG.rpsblast.tsv").read_text().splitlines()) == 1001

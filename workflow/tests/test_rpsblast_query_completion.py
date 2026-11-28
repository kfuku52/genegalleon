"""Check the family RPS-BLAST output contract through the CDSKIT-backed stage."""
import os
import subprocess
from pathlib import Path

import pytest

CORE = Path(__file__).resolve().parents[1] / "core/gg_gene_evolution_core.sh"
FIELDS = "qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle"
QUERIES = [("Species_g1", "ACDEF"), ("Species_g10", "MKK"), ("NA", "WW"), ("nan", "C")]
STUB = """#!/usr/bin/env python3
import os
import sys
from pathlib import Path

fasta = Path(sys.argv[sys.argv.index("-query") + 1])
result = Path(sys.argv[sys.argv.index("-out") + 1])
ids = [line[1:].split()[0] for line in fasta.read_text().splitlines() if line.startswith(">")]
hits = Path(os.environ["FAKE_RPS_HITS"]).read_text().splitlines(keepends=True)
with result.open("w") as handle:
    for row in hits:
        index, fields = row.split("\\t", 1)
        handle.write(ids[int(index)] + "\\t" + fields)
with Path(os.environ["FAKE_RPS_CALLS"]).open("a") as handle:
    handle.write("search\\n")
"""


def stage_script():
    source = CORE.read_text()
    start = source.index("prepare_gene_family_rpsblast_options()")
    end = source.index('task="Gene trait extraction from gff files"', start)
    return r"""
set -euo pipefail
source "${gg_support_dir}/gg_util.sh"
gg_artifact_server_request() { local program=$1; shift; python "${program}" "$@"; }
gg_step_start() { printf 'started=%s\n' "$1"; }
gg_step_skip() { printf 'skipped=%s\n' "$1"; }
""" + source[start:end]


def run_stage(tmp_path, hits, wrapped):
    primary = tmp_path / "input.fa"
    primary.write_text("".join(
        f">{name} description\n" + ("\n".join(seq) if wrapped else seq) + "\n"
        for name, seq in QUERIES
    ))
    binary = tmp_path / "bin"
    binary.mkdir()
    rpsblast = binary / "rpsblast"
    rpsblast.write_text(STUB)
    rpsblast.chmod(0o755)
    seqkit = binary / "seqkit"
    seqkit.write_text("#!/bin/sh\nexit 1\n")
    seqkit.chmod(0o755)
    profile = tmp_path / "downloads/pfam/Pfam_LE"
    profile.mkdir(parents=True)
    for suffix in ("loo", "rps"):
        (profile / f"Pfam.{suffix}").write_text("unused by the stub")
    columns = [
        (name, "pfam", "80", "1", "0", "0", "1", "1", "1", "1",
         "1e-8", "42", str(len(seq)), "8", "domain")
        for name, seq in QUERIES if name in hits
    ][::-1]
    if columns:
        columns.append((*columns[0][:-1], "second domain"))
    ids = [name for name, _ in QUERIES]
    fixture = tmp_path / "blast.tsv"
    fixture.write_text("".join(
        "\t".join((str(ids.index(row[0])), *row[1:])) + "\n" for row in columns
    ))
    calls = tmp_path / "calls"
    output = tmp_path / "output"
    table = output / "rpsblast/OG_rpsblast.tsv"
    env = dict(os.environ, PATH=str(binary) + os.pathsep + os.environ["PATH"],
               FAKE_RPS_HITS=str(fixture), FAKE_RPS_CALLS=str(calls),
               gg_support_dir=str(CORE.parent.parent / "support"),
               gg_workspace_dir=str(tmp_path), dir_output_active=str(output),
               file_og_primary_fasta=str(primary), file_og_rpsblast=str(table),
               og_id="OG", input_sequence_mode="protein", genetic_code="1",
               GG_TASK_CPUS="2", run_rps_blast="1", artifact_stale_policy="stop")
    command = ["bash", "-c", stage_script()]
    result = subprocess.run(command, cwd=tmp_path, env=env, text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr
    return table, columns, calls, command, env


@pytest.mark.parametrize("hits", [[], ["Species_g1"], ["Species_g10"], ["Species_g1", "Species_g10", "NA", "nan"]])
@pytest.mark.parametrize("wrapped", [False, True])
def test_rpsblast_stage_preserves_hit_order_and_reports_no_hit_lengths(tmp_path, hits, wrapped):
    table, columns, calls, _, _ = run_stage(tmp_path, hits, wrapped)
    expected = FIELDS.replace(" ", "\t") + "\n"
    expected += "".join("\t".join(row) + "\n" for row in columns)
    expected += "".join(
        "\t".join([name, *([""] * 11), str(len(seq)), "", ""]) + "\n"
        for name, seq in QUERIES if name not in hits
    )
    assert table.read_text() == expected
    assert calls.read_text() == "search\n"
    assert table.with_name("OG_rpsblast.search.json").is_file()


def test_rpsblast_stage_reuses_verified_results_without_blast(tmp_path):
    table, _, calls, command, env = run_stage(tmp_path, ["Species_g10"], wrapped=False)
    before = table.read_bytes()
    (tmp_path / "bin/rpsblast").unlink()
    result = subprocess.run(command, cwd=tmp_path, env=env, text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "skipped=Protein RPS-BLAST" in result.stdout
    assert table.read_bytes() == before
    assert calls.read_text() == "search\n"

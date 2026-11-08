"""Exercise native primary-FASTA resume with existing CDS outputs and provenance."""

import gzip
import json
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from workflow.support.fasta_sequence_store import build_store

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / "workflow/support"


@pytest.fixture
def family(tmp_path):
    sources = tmp_path / "species_cds"
    sources.mkdir()
    source = sources / "Species_a.fa"
    source.write_text(
        ">Species_a_keep\nATGAAA\n>Species_a_remove\nATGCCC\n>Species_a_other\nATGGGG\n"
    )
    database, index = tmp_path / "sequences.sqlite3", tmp_path / "sequences.json"
    build_store(database, index, [(source, "Species_a")], 0, 0)
    output, work = tmp_path / "output", tmp_path / "work"
    (output / "cds_fasta").mkdir(parents=True)
    work.mkdir()
    catalog = tmp_path / "Orthogroups.selected.tsv"
    fasta = output / "cds_fasta/OGtest_cds.fa.gz"
    manifest = output / "artifact_provenance/OGtest.primary_fasta.json"
    environment = dict(
        os.environ, gg_support_dir=str(SUPPORT), gg_workspace_dir=str(tmp_path),
        gg_workspace_output_dir=str(output), dir_output_active=str(output),
        file_og=str(catalog), file_og_cds_fasta=str(fasta), file_og_primary_fasta=str(fasta),
        file_species_cds_store_manifest=str(index), file_species_cds_store_db=str(database),
        file_species_genetic_code="", dir_sp_cds=str(sources), og_id="OGtest",
        input_sequence_mode="cds", mode_gene_evolution="orthogroup", genetic_code="1",
        query_blast_coverage="0.25", max_num_gene_blast_hit_retrieval="5000",
        GG_TASK_CPUS="1", run_extract_primary_fasta="1", artifact_stale_policy="rebuild",
        artifact_legacy_policy="stop",
    )
    return SimpleNamespace(root=tmp_path, output=output, work=work, catalog=catalog,
                           fasta=fasta, manifest=manifest, index=index, env=environment)


def catalog(family, genes):
    family.catalog.write_text("Orthogroup\tSpecies_a\nOGtest\t" + ", ".join(genes) + "\n")


def existing(family, genes, *, schema):
    with gzip.open(family.fasta, "wt") as stream:
        for gene in genes:
            stream.write(f">{gene}\nATGAAA\n")
    arguments = [
        "--manifest", str(family.manifest), "--step", "primary_fasta", "--family-id", "OGtest",
        "--logical-root", str(family.output), "--workspace-root", str(family.root),
        "--output", f"primary_fasta={family.fasta}", "--input", f"orthogroups={family.catalog}",
        "--input", f"species_cds_index={family.index}",
    ]
    for value in ["mode=orthogroup", "input_sequence_mode=cds", "genetic_code=1",
                  "query_blast_coverage=0.25", "max_query_hits=5000", "fasta_sequence_store_schema=2"]:
        arguments += ["--parameter", value]
    if schema:
        arguments += ["--parameter", "primary_fasta_generation_schema=2"]
    result = subprocess.run([sys.executable, str(SUPPORT / "artifact_provenance.py"),
                             "record", *arguments], text=True, capture_output=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr


def run(family, *, policy="rebuild"):
    core = (ROOT / "workflow/core/gg_gene_evolution_core.sh").read_text()
    start = core.index('task="Fasta generation"')
    stop = core.index('task="Protein RPS-BLAST"', start)
    script = r'''set -euo pipefail
source "${gg_support_dir}/gg_util.sh"
# The fixture prebuilds the real store; use the real provenance CLI per request.
ensure_species_fasta_sequence_store() { :; }
gg_artifact_server_request() { local program=$1; shift; python "${program}" "$@"; }
gg_step_start() { printf 'started=%s\n' "$1"; }
gg_step_skip() { printf 'skipped=%s\n' "$1"; }
''' + core[start:stop]
    return subprocess.run(["bash", "-c", script], env=dict(family.env, artifact_stale_policy=policy),
                          cwd=family.work, text=True, capture_output=True, timeout=45)


def headers(family):
    with gzip.open(family.fasta, "rt") as stream:
        return [line[1:].strip() for line in stream if line.startswith(">")]


def test_rebuild_replaces_existing_cds_when_catalog_membership_changes(family):
    catalog(family, ["Species_a_keep", "Species_a_remove"])
    existing(family, ["Species_a_keep", "Species_a_remove"], schema=True)
    catalog(family, ["Species_a_keep", "Species_a_other"])
    result = run(family)
    assert result.returncode == 0, result.stdout + result.stderr
    assert headers(family) == ["Species_a_keep", "Species_a_other"]
    with gzip.open(family.fasta, "rt") as stream:
        assert "ATGGGG" in stream.read()


def test_previously_restamped_legacy_cds_is_regenerated_once(family):
    catalog(family, ["Species_a_keep"])
    # The old defect recorded this current catalog with a FASTA containing an extra gene.
    existing(family, ["Species_a_keep", "Species_a_remove"], schema=False)
    result = run(family)
    assert result.returncode == 0, result.stdout + result.stderr
    assert headers(family) == ["Species_a_keep"]
    assert json.loads(family.manifest.read_text())["parameters"]["primary_fasta_generation_schema"] == "2"
    before = {p: (p.stat().st_mtime_ns, p.read_bytes()) for p in (family.fasta, family.manifest)}
    again = run(family)
    assert again.returncode == 0, again.stdout + again.stderr
    assert "skipped=Fasta generation" in again.stdout
    assert {p: (p.stat().st_mtime_ns, p.read_bytes()) for p in before} == before


@pytest.mark.parametrize("failure", ["missing-gene", "strict-stale-policy"])
def test_failed_rebuild_preserves_previous_fasta_and_provenance(family, failure):
    catalog(family, ["Species_a_keep", "Species_a_remove"])
    existing(family, ["Species_a_keep", "Species_a_remove"], schema=True)
    catalog(family, ["Species_a_keep", "Species_a_missing" if failure == "missing-gene" else "Species_a_other"])
    before = {p: (p.stat().st_mtime_ns, p.read_bytes()) for p in (family.fasta, family.manifest)}
    result = run(family, policy="stop" if failure == "strict-stale-policy" else "rebuild")
    assert result.returncode != 0
    assert {p: (p.stat().st_mtime_ns, p.read_bytes()) for p in before} == before

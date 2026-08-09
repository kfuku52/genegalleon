from __future__ import annotations

import gzip
import hashlib
import json
import os
import sqlite3
import subprocess
import sys
import time
import zlib
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SUPPORT = REPO_ROOT / "workflow" / "support"


def run_helper(script: str, *args: object) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SUPPORT / script), *(str(arg) for arg in args)],
        check=False,
        capture_output=True,
        text=True,
    )


def run_bash(script: str, *args: object, timeout: float = 10) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["bash", "-c", script, "bash", *(str(arg) for arg in args)],
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def test_background_job_limiter_is_event_driven_and_respects_the_cap():
    script = r"""
source "$1"
max_seen=0
for _index in $(seq 1 18); do
  wait_until_jobn_le 3
  (sleep 0.05) &
  gg_background_register "$!"
  current=$(jobs -pr | wc -l | tr -d '[:space:]')
  if [[ ${current} -gt ${max_seen} ]]; then
    max_seen=${current}
  fi
done
wait_for_background_jobs
printf '%s\n' "${max_seen}"
"""
    started = time.monotonic()
    result = run_bash(script, SUPPORT / "gg_util.sh")
    elapsed = time.monotonic() - started
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "3"
    # The old one-second polling loop takes several seconds for this workload.
    assert elapsed < 2.0


def test_background_job_limiter_reports_an_already_finished_failure():
    script = r"""
source "$1"
(exit 7) &
gg_background_register "$!"
sleep 0.1
if wait_for_background_jobs; then
  exit 1
fi
printf '%s\n' caught
"""
    result = run_bash(script, SUPPORT / "gg_util" / "06_workspace_validation.sh")
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "caught"


def test_memory_parallel_job_cap_has_a_one_job_floor():
    script = r"""
source "$1"
printf '%s %s\n' "$(gg_memory_parallel_job_cap 16 4)" "$(gg_memory_parallel_job_cap 3 4)"
"""
    result = run_bash(script, SUPPORT / "gg_util.sh")
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "4 1"


def test_artifact_provenance_shell_server_handles_repeated_requests(tmp_path: Path):
    workspace = tmp_path / "workspace"
    logical_root = workspace / "output" / "orthogroup"
    input_file = logical_root / "rooted_tree" / "input.txt"
    output_file = logical_root / "stat_branch" / "output.txt"
    manifest = logical_root / "artifact_provenance" / "manifest.json"
    input_file.parent.mkdir(parents=True)
    output_file.parent.mkdir(parents=True)
    input_file.write_text("input\n")
    output_file.write_text("output\n")
    script = r"""
source "$1"
gg_support_dir=$2
args=(
  --manifest "$3"
  --step server_test
  --family-id family
  --logical-root "$4"
  --workspace-root "$5"
  --input "source=$6"
  --output "result=$7"
)
gg_artifact_record "${args[@]}"
if gg_artifact_needs_run "${args[@]}"; then
  exit 10
else
  status=$?
fi
[[ ${status} -eq 1 ]]
"""
    result = run_bash(
        script,
        SUPPORT / "gg_util" / "06_workspace_validation.sh",
        SUPPORT,
        manifest,
        logical_root,
        workspace,
        input_file,
        output_file,
        timeout=5,
    )
    assert result.returncode == 0, result.stderr


def test_fasta_sequence_store_reuses_index_and_extracts_query_variants(tmp_path: Path):
    source_a = tmp_path / "Species_a.fa"
    source_a.write_text(">Species_a_gene1 description\nACGT\n>Species_a_gene2\nTTTT\n")
    source_b = tmp_path / "Species_b.fa.gz"
    with gzip.open(source_b, "wt") as handle:
        handle.write(">gene3\nCCCC\n")
    source_list = tmp_path / "sources.tsv"
    source_list.write_text(f"{source_a}\tSpecies_a\n{source_b}\tSpecies_b\n")
    database = tmp_path / "sequences.sqlite3"
    manifest = tmp_path / "sequences.json"
    digest_cache = tmp_path / "digests.sqlite3"

    first = run_helper(
        "fasta_sequence_store.py",
        "ensure",
        "--database",
        database,
        "--manifest",
        manifest,
        "--source-list",
        source_list,
        "--digest-cache",
        digest_cache,
    )
    assert first.returncode == 0, first.stderr
    first_inode = database.stat().st_ino
    first_digests = {entry["species"]: entry["sha256"] for entry in json.loads(manifest.read_text())["sources"]}
    assert str(tmp_path) not in manifest.read_text()
    second = run_helper(
        "fasta_sequence_store.py",
        "ensure",
        "--database",
        database,
        "--manifest",
        manifest,
        "--source-list",
        source_list,
        "--digest-cache",
        digest_cache,
    )
    assert second.returncode == 0, second.stderr
    assert database.stat().st_ino == first_inode

    os.utime(source_a, None)
    touched = run_helper(
        "fasta_sequence_store.py",
        "ensure",
        "--database",
        database,
        "--manifest",
        manifest,
        "--source-list",
        source_list,
        "--digest-cache",
        digest_cache,
    )
    assert touched.returncode == 0, touched.stderr
    assert database.stat().st_ino == first_inode

    patterns = tmp_path / "patterns.txt"
    patterns.write_text("gene1\ngene3\n")
    output = tmp_path / "query.fa"
    extracted = run_helper(
        "fasta_sequence_store.py",
        "extract",
        "--database",
        database,
        "--pattern-file",
        patterns,
        "--output",
        output,
        "--ignore-case",
        "--query-variants",
        "--prefix-species",
    )
    assert extracted.returncode == 0, extracted.stderr
    assert output.read_text() == ">Species_a_gene1\nACGT\n>Species_b_gene3\nCCCC\n"

    original_stat = source_a.stat()
    source_a.write_text(">Species_a_gene1 description\nTGCA\n>Species_a_gene2\nTTTT\n")
    os.utime(source_a, ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns))
    rebuilt = run_helper(
        "fasta_sequence_store.py",
        "ensure",
        "--database",
        database,
        "--manifest",
        manifest,
        "--source-list",
        source_list,
        "--digest-cache",
        digest_cache,
    )
    assert rebuilt.returncode == 0, rebuilt.stderr
    rebuilt_digests = {entry["species"]: entry["sha256"] for entry in json.loads(manifest.read_text())["sources"]}
    assert rebuilt_digests["Species_a"] != first_digests["Species_a"]
    assert rebuilt_digests["Species_b"] == first_digests["Species_b"]


def test_fasta_sequence_store_rebuilds_a_valid_but_modified_database(tmp_path: Path):
    source = tmp_path / "Species_a.fa"
    source.write_text(">Species_a_gene1\nACGT\n")
    source_list = tmp_path / "sources.tsv"
    source_list.write_text(f"{source}\tSpecies_a\n")
    database = tmp_path / "sequences.sqlite3"
    manifest = tmp_path / "sequences.json"
    ensure_args = (
        "ensure",
        "--database",
        database,
        "--manifest",
        manifest,
        "--source-list",
        source_list,
        "--minimum-free-bytes",
        0,
    )
    first = run_helper("fasta_sequence_store.py", *ensure_args)
    assert first.returncode == 0, first.stderr
    assert database.with_suffix(database.suffix + ".state.json").is_file()
    with sqlite3.connect(database) as connection:
        connection.execute(
            "UPDATE sequences SET sequence=? WHERE identifier=?",
            (sqlite3.Binary(zlib.compress(b"TTTT")), "Species_a_gene1"),
        )

    repaired = run_helper("fasta_sequence_store.py", *ensure_args)
    assert repaired.returncode == 0, repaired.stderr
    patterns = tmp_path / "patterns.txt"
    patterns.write_text("Species_a_gene1\n")
    output = tmp_path / "output.fa"
    extracted = run_helper(
        "fasta_sequence_store.py",
        "extract",
        "--database",
        database,
        "--pattern-file",
        patterns,
        "--output",
        output,
        "--require-all",
    )
    assert extracted.returncode == 0, extracted.stderr
    assert output.read_text() == ">Species_a_gene1\nACGT\n"


def test_fasta_sequence_store_compresses_sequences_and_enforces_a_size_cap(tmp_path: Path):
    source = tmp_path / "Species_a.fa.gz"
    with gzip.open(source, "wt") as handle:
        handle.write(">Species_a_gene1\n" + ("ACGT" * 250_000) + "\n")
    source_list = tmp_path / "sources.tsv"
    source_list.write_text(f"{source}\tSpecies_a\n")
    database = tmp_path / "sequences.sqlite3"
    manifest = tmp_path / "sequences.json"
    built = run_helper(
        "fasta_sequence_store.py",
        "ensure",
        "--database",
        database,
        "--manifest",
        manifest,
        "--source-list",
        source_list,
        "--minimum-free-bytes",
        0,
    )
    assert built.returncode == 0, built.stderr
    assert database.stat().st_size < 100_000

    capped = run_helper(
        "fasta_sequence_store.py",
        "ensure",
        "--database",
        tmp_path / "capped.sqlite3",
        "--manifest",
        tmp_path / "capped.json",
        "--source-list",
        source_list,
        "--max-database-bytes",
        1_000,
        "--minimum-free-bytes",
        0,
    )
    assert capped.returncode == 2
    assert "size limit" in capped.stderr


def test_content_digest_cache_is_shared_safely_across_threads(tmp_path: Path):
    files = []
    for index in range(32):
        path = tmp_path / f"input-{index}.txt"
        path.write_text(f"record {index}\n")
        files.append(path)
    script = r"""
import sqlite3
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, sys.argv[1])
from content_digest_cache import cached_sha256_file, configure

database = Path(sys.argv[2])
paths = [Path(value) for value in sys.argv[3:]]
configure(database)
with ThreadPoolExecutor(max_workers=8) as executor:
    list(executor.map(cached_sha256_file, paths))
with sqlite3.connect(database) as connection:
    print(connection.execute("SELECT COUNT(*) FROM digest_cache").fetchone()[0])
"""
    result = subprocess.run(
        [sys.executable, "-c", script, str(SUPPORT), str(tmp_path / "digests.sqlite3"), *(str(path) for path in files)],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "32"


def test_batch_busco_extraction_streams_inputs_once_for_both_modes(tmp_path: Path):
    fasta_a = tmp_path / "a.fa"
    fasta_b = tmp_path / "b.fa"
    fasta_a.write_text(">A_single desc\nAAAA\n>A_dup1\nCCCC\n")
    fasta_b.write_text(">B_single\nTTTT\n>B_dup1\nGGGG\n>B_dup2\nACAC\n")
    fasta_list = tmp_path / "fastas.txt"
    fasta_list.write_text(f"{fasta_a}\n{fasta_b}\n")
    summary = tmp_path / "summary.tsv"
    summary.write_text(
        "busco_id\tnum_present\tnum_single\tA\tB\n"
        "BUSCO1\t2\t2\tA_single\tB_single\n"
        "BUSCO2\t2\t0\tA_dup1\tB_dup1,B_dup2\n"
    )

    strict_dir = tmp_path / "strict"
    strict = run_helper(
        "batch_extract_busco_fasta.py",
        "--summary",
        summary,
        "--fasta-list",
        fasta_list,
        "--output-dir",
        strict_dir,
        "--suffix",
        ".fa",
        "--report",
        tmp_path / "strict.tsv",
        "--mode",
        "single-copy",
        "--strict",
        "--require-complete",
    )
    assert strict.returncode == 0, strict.stderr
    assert (strict_dir / "BUSCO1.fa").read_text() == ">A_single desc\nAAAA\n>B_single\nTTTT\n"
    assert not (strict_dir / "BUSCO2.fa").exists()

    duplicate_dir = tmp_path / "duplicate"
    duplicate = run_helper(
        "batch_extract_busco_fasta.py",
        "--summary",
        summary,
        "--fasta-list",
        fasta_list,
        "--output-dir",
        duplicate_dir,
        "--suffix",
        ".fa",
        "--report",
        tmp_path / "duplicate.tsv",
        "--mode",
        "duplicate-aware",
        "--require-complete",
    )
    assert duplicate.returncode == 0, duplicate.stderr
    duplicate_text = (duplicate_dir / "BUSCO2.fa").read_text()
    assert duplicate_text.count(">") == 3
    assert ">A_dup1" in duplicate_text
    assert ">B_dup2" in duplicate_text


def test_batch_busco_extraction_rejects_traversal_identifiers(tmp_path: Path):
    fasta = tmp_path / "species.fa"
    fasta.write_text(">gene1\nAAAA\n")
    fasta_list = tmp_path / "fastas.txt"
    fasta_list.write_text(f"{fasta}\n")
    summary = tmp_path / "summary.tsv"
    summary.write_text("busco_id\tnum_present\tnum_single\tA\n../escaped\t1\t1\tgene1\n")
    result = run_helper(
        "batch_extract_busco_fasta.py",
        "--summary",
        summary,
        "--fasta-list",
        fasta_list,
        "--output-dir",
        tmp_path / "outputs",
        "--suffix",
        ".fa",
        "--report",
        tmp_path / "report.tsv",
        "--mode",
        "single-copy",
    )
    assert result.returncode == 2
    assert "Invalid BUSCO identifier" in result.stderr
    assert not (tmp_path / "escaped.fa").exists()


def test_batch_busco_extraction_refuses_existing_output_symlinks(tmp_path: Path):
    fasta = tmp_path / "species.fa"
    fasta.write_text(">gene1\nAAAA\n")
    fasta_list = tmp_path / "fastas.txt"
    fasta_list.write_text(f"{fasta}\n")
    summary = tmp_path / "summary.tsv"
    summary.write_text("busco_id\tnum_present\tnum_single\tA\nBUSCO1\t1\t1\tgene1\n")
    output_dir = tmp_path / "outputs"
    output_dir.mkdir()
    victim = tmp_path / "victim.fa"
    victim.write_text("do not replace\n")
    (output_dir / "BUSCO1.fa").symlink_to(victim)
    result = run_helper(
        "batch_extract_busco_fasta.py",
        "--summary",
        summary,
        "--fasta-list",
        fasta_list,
        "--output-dir",
        output_dir,
        "--suffix",
        ".fa",
        "--report",
        tmp_path / "report.tsv",
        "--mode",
        "single-copy",
    )
    assert result.returncode == 2
    assert "symlink" in result.stderr
    assert victim.read_text() == "do not replace\n"


def test_busco_reuse_requires_matching_content_and_resolved_lineage(tmp_path: Path):
    species_input = tmp_path / "Species_a.fa"
    full_source = tmp_path / "Species_a.busco.full.tsv"
    short_source = tmp_path / "Species_a.busco.short.txt"
    species_input.write_text(">Species_a_gene1\nATG\n")
    full_source.write_text("BUSCO1\tComplete\tSpecies_a_gene1\n")
    short_source.write_text("C:100.0%\n")

    def entry(label: str, path: Path) -> dict[str, object]:
        data = path.read_bytes()
        return {
            "label": label,
            "sha256": hashlib.sha256(data).hexdigest(),
            "size_bytes": len(data),
        }

    manifest = tmp_path / "busco.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "step": "input_generation_species_busco",
                "family_id": "Species_a",
                "inputs": [entry("species_cds", species_input)],
                "outputs": [entry("busco_full", full_source), entry("busco_short", short_source)],
                "parameters": {
                    "busco_mode": "transcriptome",
                    "busco_lineage_resolved": "eukaryota_odb12",
                    "evalue": "1e-03",
                    "limit": "20",
                },
            }
        )
    )
    full_destination = tmp_path / "promoted" / "full.tsv"
    short_destination = tmp_path / "promoted" / "short.txt"
    reused = run_helper(
        "reuse_busco_outputs.py",
        "--manifest",
        manifest,
        "--species",
        "Species_a",
        "--input",
        species_input,
        "--lineage",
        "eukaryota_odb12",
        "--full-source",
        full_source,
        "--short-source",
        short_source,
        "--full-destination",
        full_destination,
        "--short-destination",
        short_destination,
        "--digest-cache",
        tmp_path / "digests.sqlite3",
    )
    assert reused.returncode == 0, reused.stderr
    assert full_destination.read_bytes() == full_source.read_bytes()
    assert short_destination.read_bytes() == short_source.read_bytes()

    wrong_lineage = run_helper(
        "reuse_busco_outputs.py",
        "--manifest",
        manifest,
        "--species",
        "Species_a",
        "--input",
        species_input,
        "--lineage",
        "metazoa_odb12",
        "--full-source",
        full_source,
        "--short-source",
        short_source,
        "--full-destination",
        full_destination,
        "--short-destination",
        short_destination,
    )
    assert wrong_lineage.returncode == 1

    original_stat = full_source.stat()
    full_source.write_text("BUSCO2\tComplete\tSpecies_a_gene1\n")
    os.utime(full_source, ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns))
    changed_content = run_helper(
        "reuse_busco_outputs.py",
        "--manifest",
        manifest,
        "--species",
        "Species_a",
        "--input",
        species_input,
        "--lineage",
        "eukaryota_odb12",
        "--full-source",
        full_source,
        "--short-source",
        short_source,
        "--full-destination",
        full_destination,
        "--short-destination",
        short_destination,
        "--digest-cache",
        tmp_path / "digests.sqlite3",
    )
    assert changed_content.returncode == 1

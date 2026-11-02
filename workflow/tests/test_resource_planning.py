import gzip
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "support"))
from resource_metrics import run  # noqa: E402 - support directory is added above
from resource_profile import profile, scan  # noqa: E402


def test_profile_streams_compressed_input_and_binds_settings(tmp_path):
    first = tmp_path / "species1.fa.gz"
    with gzip.open(first, "wt") as stream:
        stream.write(">a\nACGT\n>b\nAC\n")
    second = tmp_path / "species2.fa"
    second.write_text(">c\nAACCGG\n")
    result = profile("gg_genome_evolution", [first, second], settings={"mode": "protein"})
    assert result["features"]["species"] == 2
    assert result["features"]["sequences"] == 3
    assert result["features"]["residues"] == 12
    assert profile("gg_genome_evolution", [second, first], settings={"mode": "protein"})["input_sha256"] == result["input_sha256"]
    assert profile("gg_genome_evolution", [first, second], settings={"mode": "DNA"})["input_sha256"] != result["input_sha256"]
    with gzip.open(first, "wt") as stream:
        stream.write(">a\nCCCC\n>b\nAC\n")
    assert profile("gg_genome_evolution", [first, second], settings={"mode": "protein"})["input_sha256"] != result["input_sha256"]


def test_gene_species_and_array_input_are_per_member(tmp_path):
    fasta = tmp_path / "gene.fa"
    fasta.write_text(">a\nAA\n>b\nCC\n>c\nGG\n")
    result = profile("gg_gene_evolution", [fasta], task_count=75000)
    assert result["features"]["species"] == 3
    assert result["task_count"] == 75000
    assert result["features"]["sequences"] == 3


def test_bad_fastq_and_duplicate_inputs_are_rejected(tmp_path):
    fastq = tmp_path / "reads.fastq"
    fastq.write_text("@a\nACGT\n+\n!!!\n")
    with pytest.raises(ValueError):
        profile("gg_transcriptome_generation", fastq=[fastq])
    fasta = tmp_path / "a.fa"
    fasta.write_text(">a\nAA\n")
    with pytest.raises(ValueError):
        profile("gg_genome_evolution", [fasta, fasta])


@pytest.mark.parametrize("newline", [b"\n", b"\r\n"])
@pytest.mark.parametrize("final_newline", [False, True])
@pytest.mark.parametrize("compressed", [False, True])
def test_fastq_profile_preserves_exact_decoded_hash_and_record_lengths(
        tmp_path, newline, final_newline, compressed):
    content = newline.join([b"@first", b"ACGT", b"+first", b"!!!!",
                            b"@second", b"AC", b"+", b"II"])
    if final_newline:
        content += newline
    path = tmp_path / ("reads.fastq.gz" if compressed else "reads.fastq")
    path.write_bytes(gzip.compress(content) if compressed else content)
    observed = scan(path, "fastq")
    assert observed["sha256"] == hashlib.sha256(content).hexdigest()
    assert observed["sequences"] == 2
    assert observed["residues"] == 6
    assert observed["max_length"] == 4
    assert observed["bytes"] == path.stat().st_size


@pytest.mark.parametrize("record", [
    b"invalid\nAC\n+\nII\n", b"@a\nAC\ninvalid\nII\n",
    b"@a\n\n+\n\n", b"@a\nAC\n+\n", b"@a\nAC\n",
    b"@a\nAC\n+\nI\n", b"@a\nAC\n+\nIII\n",
])
def test_fastq_profile_rejects_late_invalid_or_truncated_records(tmp_path, record):
    path = tmp_path / "reads.fastq"
    path.write_bytes(b"@valid\nACGT\n+\nIIII\n" + record)
    with pytest.raises(ValueError, match="invalid four-line FASTQ"):
        scan(path, "fastq")


def test_fastq_profile_consumes_compressed_tail(tmp_path):
    path = tmp_path / "reads.fastq.gz"
    path.write_bytes(gzip.compress(b"@a\nAC\n+\nII\n")[:-4])
    with pytest.raises(EOFError):
        scan(path, "fastq")


def test_fastq_profile_rejects_mutation_during_scan(tmp_path, monkeypatch):
    import os

    import resource_profile

    path = tmp_path / "reads.fastq"
    path.write_bytes(b"@a\nAC\n+\nII\n")
    original_mtime = path.stat().st_mtime_ns
    checksum = hashlib.sha256()

    class MutatingChecksum:
        def update(self, value):
            checksum.update(value)
            os.utime(path, ns=(original_mtime, original_mtime + 1000000000))

        def hexdigest(self):
            return checksum.hexdigest()

    monkeypatch.setattr(resource_profile.hashlib, "sha256", MutatingChecksum)
    with pytest.raises(ValueError, match="input changed during profiling"):
        scan(path, "fastq")


def test_metrics_preserve_stdout_exit_and_measure_child_cpu(tmp_path, capfd):
    code = "print('scientific-output'); sum(i*i for i in range(1000000)); raise SystemExit(7)"
    status = run([sys.executable, "-c", code], tmp_path, workflow="gg_gene_evolution",
                 input_sha256="a" * 64, cpus=2, memory_gb=4, runtime_id="test", server_id="test")
    assert status == 7
    assert capfd.readouterr().out == "scientific-output\n"
    record = json.loads(next(tmp_path.glob("attempt-*.json")).read_text())
    assert record["exit_code"] == 7
    assert record["cpu_seconds"] > 0
    assert record["wall_seconds"] > 0
    assert record["peak_process_rss_gb"] > 0


def test_success_publishes_comparable_history_and_telemetry_failure_is_nonfatal(tmp_path, capfd):
    import hashlib
    key = {"workflow": "gg_gene_evolution", "input_sha256": "a" * 64,
           "runtime_id": "runtime-a", "server_id": "test-server", "stage": "workflow"}
    key_hash = hashlib.sha256(json.dumps(key, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    code = "import os,json,time; open(os.environ['GG_RESOURCE_EVENT_FILE'],'w').write(json.dumps(dict(stage='compute',kind='start',at=time.monotonic(),cpu=0))+'\\n')"
    assert run([sys.executable, "-c", code], tmp_path,
               workflow=key["workflow"], input_sha256=key["input_sha256"],
               cpus=4, memory_gb=8, runtime_id="runtime-a", server_id="test-server", expected_stages=['compute']) == 0
    indexed = tmp_path / "history" / key_hash / "cpus-4.json"
    assert json.loads(indexed.read_text())["exit_code"] == 0
    previous = indexed.read_bytes()
    for skipped in ("pass", code.replace("kind='start'", "kind='skip'")):
        assert run([sys.executable, "-c", skipped], tmp_path,
                   workflow=key["workflow"], input_sha256=key["input_sha256"], cpus=4,
                   memory_gb=8, runtime_id="runtime-a", server_id="test-server", expected_stages=['compute']) == 0
        assert indexed.read_bytes() == previous
    blocked = tmp_path / "not-a-directory"
    blocked.write_text("owned-file")
    assert run([sys.executable, "-c", "print('still-runs')"], blocked,
               workflow=key["workflow"], input_sha256=None, cpus=4,
               memory_gb=8, runtime_id="unidentified") == 0
    captured = capfd.readouterr()
    assert "still-runs" in captured.out
    assert "Resource metrics" in captured.err
    assert blocked.read_text() == "owned-file"


def test_core_parallel_limits_execute_from_current_source(tmp_path):
    source = (ROOT / "core/gg_genome_evolution_core.sh").read_text()
    start = source.index("GG_GENOME_PARALLEL_JOBS=${GG_TASK_CPUS}")
    end = source.index("iqtree_full_mem_args=", start)
    block = source[start:end]
    # The budget helper is the real sourced implementation; no scientific tool
    # stand-in is used to infer any end-to-end speedup.
    util = ROOT / "support/gg_util.sh"
    code = f'source "{util}"\nGG_TASK_CPUS=32\nGG_MEM_TOOL_GB=8\ngenome_parallel_jobs=auto\ngenome_parallel_memory_gb_per_job=2\n' + block + '\nprintf "jobs=%s\\n" "$GG_GENOME_PARALLEL_JOBS"\n'
    result = subprocess.run(["bash", "-c", code], text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    assert "jobs=4" in result.stdout

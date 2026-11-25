"""Exercise FASTQ discovery with the runtime's real amalgkit mate and metadata checks."""

import gzip
import hashlib
import json
import subprocess
from pathlib import Path

import pytest
from test_transcriptome_fastq_discovery import discover


@pytest.mark.parametrize("problem", [None, "swapped", "different", "empty"])
def test_submission_preflight_checks_mates_at_published_paths(tmp_path, problem):
    reads = tmp_path / "inputs"
    reads.mkdir()
    for mate in (1, 2):
        read_id = "different" if mate == 2 and problem == "different" else "same"
        header_mate = 3 - mate if problem == "swapped" else mate
        payload = b"" if problem == "empty" else f"@{read_id}/{header_mate}\nAC\n+\nII\n".encode()
        source = tmp_path / f"raw{mate}.fq.gz"
        source.write_bytes(gzip.compress(payload))
        (reads / f"run_{mate}.fq.gz").symlink_to(source)
    result = discover(reads, "paired")
    if problem is None:
        assert result.returncode == 0, result.stderr
        assert len(result.stdout.splitlines()) == 2
    else:
        assert result.returncode != 0
        assert result.stdout == ""
        assert {"swapped": "mate annotation", "different": "IDs/order", "empty": "Empty paired-end"}[problem] in result.stderr


def test_quant_resolves_published_gsa_metadata_against_original_digest(tmp_path):
    source = tmp_path / "selected.tsv"
    source.write_text("run\tdata_source\tgsa_selection_status\nCRR1\tgsa\tpending_input_counts\n")
    reads = tmp_path / "finished_reads" / "Species_one"
    reads.mkdir(parents=True)
    measured = reads / "metadata.tsv"
    measured.write_text("run\tdata_source\tread_count_status\ttotal_spots\nCRR1\tgsa\tmeasured\t12\n")
    (reads / "gsa_metadata_state.json").write_text(json.dumps({
        "source_digest": hashlib.sha256(source.read_bytes()).hexdigest(),
        "snapshot_digest": hashlib.sha256(measured.read_bytes()).hexdigest(),
    }))
    argv = [
        "bash", "-c", 'set -euo pipefail; GG_CORE_SOURCE_ONLY=1 source "$1"; '
        'resolve_transcriptome_quant_metadata "$2" "$3"', "_", str(Path(__file__).resolve().parents[1] / 'core/gg_transcriptome_generation_core.sh'), str(source), str(reads),
    ]
    result = subprocess.run(argv, text=True, capture_output=True, check=False)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == str(measured)
    with source.open("a") as handle:
        handle.write("CRR2\tgsa\tpending_input_counts\n")
    changed = subprocess.run(argv, text=True, capture_output=True, check=False)
    assert changed.returncode != 0
    assert "does not match the selected source" in changed.stderr

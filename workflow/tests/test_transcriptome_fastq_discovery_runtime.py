"""Exercise FASTQ discovery with the runtime's real amalgkit mate checks."""

import gzip

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

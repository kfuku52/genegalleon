import gzip
import hashlib
import json
import re
import sys
import urllib.request
from pathlib import Path

import pytest

CORE_PATH = Path(__file__).resolve().parents[1] / "core" / "gg_transcriptome_generation_core.sh"


class _Response:
    def __init__(self, payload: bytes):
        self.payload = payload

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return False

    def read(self):
        return self.payload


def _fallback_python_source() -> str:
    text = CORE_PATH.read_text(encoding="utf-8")
    function_start = text.index("download_public_original_fastqs_for_metadata() {")
    function_end = text.index("\nrun_amalgkit_getfastq_or_fallback() {", function_start)
    function_body = text[function_start:function_end]
    match = re.search(r"<<'PY'\n(.*?)\nPY\n", function_body, re.DOTALL)
    assert match is not None
    return match.group(1)


def _run_fallback(monkeypatch, metadata_path: Path, output_dir: Path, responses: dict[str, bytes]):
    def fake_urlopen(url, timeout):
        assert timeout == 120
        return _Response(responses[url])

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    old_argv = sys.argv
    sys.argv = ["fallback", str(metadata_path), str(output_dir)]
    try:
        exec(compile(_fallback_python_source(), str(CORE_PATH), "exec"), {"__name__": "__main__"})
    finally:
        sys.argv = old_argv


def _metadata(path: Path, runs: list[str]):
    path.write_text("run\n" + "\n".join(runs) + "\n", encoding="utf-8")


def _xml(filename: str, url: str) -> bytes:
    return (
        '<ROOT><SRAFile semantic_name="fastq" supertype="Original" '
        f'filename="{filename}" url="{url}" /></ROOT>'
    ).encode()


def test_public_fallback_reuses_valid_fastq_and_atomically_completes_missing_run(monkeypatch, tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "getfastq"
    _metadata(metadata_path, ["RUN1", "RUN2"])

    existing = output_dir / "RUN1" / "RUN1.amalgkit.fastq.gz"
    existing.parent.mkdir(parents=True)
    with gzip.open(existing, "wb") as handle:
        handle.write(b"@existing\nACGT\n+\n!!!!\n")
    existing_hash = hashlib.sha256(existing.read_bytes()).hexdigest()

    run1_url = "https://example.invalid/RUN1.fastq.gz"
    run2_url = "https://example.invalid/RUN2.fastq"
    responses = {
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=RUN1": _xml("RUN1.fastq.gz", run1_url),
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=RUN2": _xml("RUN2.fastq", run2_url),
        run2_url: b"@downloaded\nTGCA\n+\n!!!!\n",
    }

    _run_fallback(monkeypatch, metadata_path, output_dir, responses)

    assert hashlib.sha256(existing.read_bytes()).hexdigest() == existing_hash
    downloaded = output_dir / "RUN2" / "RUN2.amalgkit.fastq.gz"
    with gzip.open(downloaded, "rb") as handle:
        assert handle.read() == b"@downloaded\nTGCA\n+\n!!!!\n"
    assert not list(output_dir.rglob("*.part"))
    manifest = json.loads((output_dir / "getfastq_completion.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "complete"
    assert manifest["run_count"] == 2
    assert [entry["run"] for entry in manifest["runs"]] == ["RUN1", "RUN2"]


def test_public_fallback_preserves_existing_outputs_and_prior_manifest_on_failure(monkeypatch, tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "getfastq"
    _metadata(metadata_path, ["RUN1", "RUN2"])

    existing = output_dir / "RUN1" / "RUN1.amalgkit.fastq.gz"
    existing.parent.mkdir(parents=True)
    with gzip.open(existing, "wb") as handle:
        handle.write(b"@existing\nACGT\n+\n!!!!\n")
    existing_hash = hashlib.sha256(existing.read_bytes()).hexdigest()
    previous_manifest = {"status": "incomplete", "run_count": 1, "runs": [{"run": "RUN1"}]}
    (output_dir / "getfastq_completion.json").write_text(
        json.dumps(previous_manifest) + "\n",
        encoding="utf-8",
    )

    run1_url = "https://example.invalid/RUN1.fastq.gz"
    responses = {
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=RUN1": _xml("RUN1.fastq.gz", run1_url),
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=RUN2": b"<ROOT />",
    }

    with pytest.raises(SystemExit, match="No public original FASTQ URLs were found for run: RUN2"):
        _run_fallback(monkeypatch, metadata_path, output_dir, responses)

    assert hashlib.sha256(existing.read_bytes()).hexdigest() == existing_hash
    assert not (output_dir / "getfastq_completion.json").exists()
    preserved = output_dir / "getfastq_completion.pre_public_fallback.json"
    assert json.loads(preserved.read_text(encoding="utf-8")) == previous_manifest
    assert not list(output_dir.rglob("*.part"))

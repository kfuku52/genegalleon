import gzip
import hashlib
import json
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Union

import pytest

CORE_PATH = Path(__file__).resolve().parents[1] / "core" / "gg_transcriptome_generation_core.sh"


class _Response:
    def __init__(self, payload: bytes):
        self.payload = payload
        self.offset = 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return False

    def read(self, size=-1):
        assert size >= 0, "fallback downloads must use bounded reads"
        chunk = self.payload[self.offset : self.offset + size]
        self.offset += len(chunk)
        return chunk


def _fallback_python_source() -> str:
    text = CORE_PATH.read_text(encoding="utf-8")
    function_start = text.index("download_public_original_fastqs_for_metadata() {")
    function_end = text.index("\nrun_amalgkit_getfastq_or_fallback() {", function_start)
    function_body = text[function_start:function_end]
    match = re.search(r"<<'PY'\n(.*?)\nPY\n", function_body, re.DOTALL)
    assert match is not None
    return match.group(1)


def _shell_function_source(name: str) -> str:
    text = CORE_PATH.read_text(encoding="utf-8")
    start_match = re.search(rf"^{re.escape(name)}\(\) \{{\n", text, re.MULTILINE)
    assert start_match is not None
    next_match = re.search(r"^[A-Za-z_][A-Za-z0-9_]*\(\) \{\n", text[start_match.end() :], re.MULTILINE)
    if next_match is None:
        end_match = re.search(r"^\}$", text[start_match.end() :], re.MULTILINE)
        assert end_match is not None
        end = start_match.end() + end_match.end()
    else:
        end = start_match.end() + next_match.start()
    return text[start_match.start() : end].rstrip()


def _run_fallback(
    monkeypatch,
    metadata_path: Path,
    output_dir: Path,
    responses: dict[str, Union[bytes, list[bytes]]],
    attempt_counts=None,
):
    attempts: dict[str, int] = attempt_counts if attempt_counts is not None else {}

    def fake_urlopen(url, timeout):
        assert timeout == 120
        attempts[url] = attempts.get(url, 0) + 1
        response = responses[url]
        if isinstance(response, list):
            assert response, "response sequence was exhausted"
            payload = response.pop(0)
        else:
            payload = response
        return _Response(payload)

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(time, "sleep", lambda _seconds: None)
    old_argv = sys.argv
    sys.argv = ["fallback", str(metadata_path), str(output_dir)]
    try:
        exec(compile(_fallback_python_source(), str(CORE_PATH), "exec"), {"__name__": "__main__"})
    finally:
        sys.argv = old_argv
    return attempts


def _metadata(path: Path, runs: list[str]):
    path.write_text("run\n" + "\n".join(runs) + "\n", encoding="utf-8")


def _partial_files(root: Path) -> list[Path]:
    return [path for path in root.rglob("*") if ".part." in path.name or path.name.endswith(".part")]


def _manifest_paths(manifest: dict, run_index: int = 0) -> list[str]:
    return [item["path"] for item in manifest["runs"][run_index]["files"]]


def _run_manifest_validator(metadata_path: Path, output_dir: Path):
    script = "\n".join(
        [
            "set -u",
            _shell_function_source("validate_amalgkit_getfastq_completion_manifest"),
            'validate_amalgkit_getfastq_completion_manifest "$1/getfastq_completion.json" "$2"',
        ]
    )
    return subprocess.run(
        ["bash", "-s", "--", str(output_dir), str(metadata_path)],
        input=script,
        text=True,
        capture_output=True,
        check=False,
    )


def _run_manifest_index_validator(metadata_path: Path, output_dir: Path):
    script = "\n".join(
        [
            "set -u",
            _shell_function_source("validate_amalgkit_getfastq_completion_manifest_index"),
            'validate_amalgkit_getfastq_completion_manifest_index "$1/getfastq_completion.json" "$2"',
        ]
    )
    return subprocess.run(
        ["bash", "-s", "--", str(output_dir), str(metadata_path)],
        input=script,
        text=True,
        capture_output=True,
        check=False,
    )


def _run_contract_recovery(output_dir: Path):
    script = "\n".join(
        [
            "set -u",
            _shell_function_source("prepare_amalgkit_getfastq_contract_recovery"),
            'prepare_amalgkit_getfastq_contract_recovery "$1/getfastq_completion.json"',
        ]
    )
    return subprocess.run(
        ["bash", "-s", "--", str(output_dir)],
        input=script,
        text=True,
        capture_output=True,
        check=False,
    )


def _run_manifest_binder(metadata_path: Path, output_dir: Path):
    script = "\n".join(
        [
            "set -u",
            _shell_function_source("bind_amalgkit_getfastq_completion_manifest"),
            'bind_amalgkit_getfastq_completion_manifest "$1/getfastq_completion.json" "$2"',
        ]
    )
    return subprocess.run(
        ["bash", "-s", "--", str(output_dir), str(metadata_path)],
        input=script,
        text=True,
        capture_output=True,
        check=False,
    )


def _xml(filename: str, url: str) -> bytes:
    return (
        '<ROOT><SRAFile semantic_name="fastq" supertype="Original" '
        f'filename="{filename}" url="{url}" /></ROOT>'
    ).encode()


def _ena_report_url(run: str) -> str:
    return "https://www.ebi.ac.uk/ena/portal/api/filereport?{}".format(
        urllib.parse.urlencode(
            {
                "accession": run,
                "result": "read_run",
                "fields": "run_accession,fastq_ftp,fastq_md5,submitted_ftp,submitted_md5,submitted_format",
                "format": "json",
            }
        )
    )


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
    assert not _partial_files(output_dir)
    manifest = json.loads((output_dir / "getfastq_completion.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == 3
    assert manifest["status"] == "complete"
    assert manifest["read_source"] == "public-original"
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
        _ena_report_url("RUN2"): b"[]",
    }

    with pytest.raises(SystemExit, match="No supported public FASTQ URLs were found for run: RUN2"):
        _run_fallback(monkeypatch, metadata_path, output_dir, responses)

    assert hashlib.sha256(existing.read_bytes()).hexdigest() == existing_hash
    assert not (output_dir / "getfastq_completion.json").exists()
    preserved = output_dir / "getfastq_completion.pre_public_fallback.json"
    assert json.loads(preserved.read_text(encoding="utf-8")) == previous_manifest
    assert not _partial_files(output_dir)


def test_public_fallback_uses_ena_fastq_when_trace_has_no_original(monkeypatch, tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "getfastq"
    _metadata(metadata_path, ["ERR123456"])

    fastq_payload = gzip.compress(b"@ena\nACGT\n+\n!!!!\n")
    fastq_url = "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR123456.fastq.gz"
    responses = {
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=ERR123456": b"<ROOT />",
        _ena_report_url("ERR123456"): json.dumps(
            [
                {
                    "run_accession": "ERR123456",
                    "fastq_ftp": fastq_url.removeprefix("https://"),
                    "fastq_md5": hashlib.md5(fastq_payload).hexdigest(),
                }
            ]
        ).encode(),
        fastq_url: fastq_payload,
    }

    _run_fallback(monkeypatch, metadata_path, output_dir, responses)

    recovered = output_dir / "ERR123456" / "ERR123456.amalgkit.fastq.gz"
    assert recovered.read_bytes() == fastq_payload
    manifest = json.loads(
        (output_dir / "getfastq_completion.json").read_text(encoding="utf-8")
    )
    assert _manifest_paths(manifest) == ["ERR123456/ERR123456.amalgkit.fastq.gz"]
    assert manifest["runs"][0]["files"][0]["size"] == recovered.stat().st_size
    assert manifest["runs"][0]["files"][0]["sha256"] == hashlib.sha256(fastq_payload).hexdigest()


def test_public_fallback_uses_checksum_validated_ena_submitted_fastq_pair(monkeypatch, tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "getfastq"
    _metadata(metadata_path, ["ERR4643641"])

    first_payload = gzip.compress(b"@ena1\nACGT\n+\n!!!!\n")
    second_payload = gzip.compress(b"@ena2\nTGCA\n+\n!!!!\n")
    first_url = "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR464/001/ERR4643641_1.fastq.gz"
    second_url = "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR464/001/ERR4643641_2.fastq.gz"
    responses = {
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=ERR4643641": b"<ROOT />",
        _ena_report_url("ERR4643641"): json.dumps(
            [
                {
                    "run_accession": "ERR4643641",
                    "fastq_ftp": "",
                    "fastq_md5": "",
                    "submitted_ftp": ";".join(
                        [first_url.removeprefix("https://"), second_url.removeprefix("https://")]
                    ),
                    "submitted_md5": ";".join(
                        [hashlib.md5(first_payload).hexdigest(), hashlib.md5(second_payload).hexdigest()]
                    ),
                    "submitted_format": "FASTQ;FASTQ",
                }
            ]
        ).encode(),
        first_url: first_payload,
        second_url: second_payload,
    }

    _run_fallback(monkeypatch, metadata_path, output_dir, responses)

    assert (output_dir / "ERR4643641" / "ERR4643641_1.amalgkit.fastq.gz").read_bytes() == first_payload
    assert (output_dir / "ERR4643641" / "ERR4643641_2.amalgkit.fastq.gz").read_bytes() == second_payload
    manifest = json.loads((output_dir / "getfastq_completion.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "complete"
    assert _manifest_paths(manifest) == [
        "ERR4643641/ERR4643641_1.amalgkit.fastq.gz",
        "ERR4643641/ERR4643641_2.amalgkit.fastq.gz",
    ]
    assert not _partial_files(output_dir)


@pytest.mark.parametrize(
    ("submitted_format", "submitted_name"),
    [("BAM", "ERR4643641.bam"), ("FASTQ", "ERR4643641.cram")],
)
def test_public_fallback_rejects_non_fastq_ena_submitted_files(
    monkeypatch,
    tmp_path,
    submitted_format,
    submitted_name,
):
    metadata_path = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "getfastq"
    _metadata(metadata_path, ["ERR4643641"])
    submitted_url = "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR464/001/{}".format(submitted_name)
    responses = {
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=ERR4643641": b"<ROOT />",
        _ena_report_url("ERR4643641"): json.dumps(
            [
                {
                    "run_accession": "ERR4643641",
                    "fastq_ftp": "",
                    "submitted_ftp": submitted_url.removeprefix("https://"),
                    "submitted_md5": "0" * 32,
                    "submitted_format": submitted_format,
                }
            ]
        ).encode(),
    }

    with pytest.raises(SystemExit, match="No supported public FASTQ URLs were found"):
        _run_fallback(monkeypatch, metadata_path, output_dir, responses)

    assert not (output_dir / "getfastq_completion.json").exists()
    assert not _partial_files(output_dir)


def test_public_fallback_reuses_one_valid_mate_and_recovers_the_other(monkeypatch, tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "getfastq"
    _metadata(metadata_path, ["ERR123456"])

    run_dir = output_dir / "ERR123456"
    run_dir.mkdir(parents=True)
    existing = run_dir / "ERR123456_1.amalgkit.fastq.gz"
    with gzip.open(existing, "wb") as handle:
        handle.write(b"@existing\nACGT\n+\n!!!!\n")
    existing_hash = hashlib.sha256(existing.read_bytes()).hexdigest()
    first_payload = existing.read_bytes()
    second_payload = gzip.compress(b"@ena2\nTGCA\n+\n!!!!\n")
    first_url = "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR123456_1.fastq.gz"
    second_url = "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR123456_2.fastq.gz"
    responses = {
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=ERR123456": b"<ROOT />",
        _ena_report_url("ERR123456"): json.dumps(
            [
                {
                    "run_accession": "ERR123456",
                    "fastq_ftp": ";".join(
                        [
                            first_url.removeprefix("https://"),
                            second_url.removeprefix("https://"),
                        ]
                    ),
                    "fastq_md5": ";".join(
                        [
                            hashlib.md5(first_payload).hexdigest(),
                            hashlib.md5(second_payload).hexdigest(),
                        ]
                    ),
                }
            ]
        ).encode(),
        second_url: second_payload,
    }

    _run_fallback(monkeypatch, metadata_path, output_dir, responses)

    assert hashlib.sha256(existing.read_bytes()).hexdigest() == existing_hash
    assert (run_dir / "ERR123456_2.amalgkit.fastq.gz").read_bytes() == second_payload
    manifest = json.loads(
        (output_dir / "getfastq_completion.json").read_text(encoding="utf-8")
    )
    assert manifest["run_count"] == 1
    assert _manifest_paths(manifest) == [
        "ERR123456/ERR123456_1.amalgkit.fastq.gz",
        "ERR123456/ERR123456_2.amalgkit.fastq.gz",
    ]


def test_completion_manifest_rejects_fastq_content_changed_after_publication(monkeypatch, tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "getfastq"
    _metadata(metadata_path, ["RUN1"])
    fastq_url = "https://example.invalid/RUN1.fastq"
    responses = {
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=RUN1": _xml(
            "RUN1.fastq", fastq_url
        ),
        fastq_url: b"@read\nACGT\n+\n!!!!\n",
    }
    _run_fallback(monkeypatch, metadata_path, output_dir, responses)

    current = _run_manifest_validator(metadata_path, output_dir)
    assert current.returncode == 0, current.stderr

    fastq_path = output_dir / "RUN1" / "RUN1.amalgkit.fastq.gz"
    fastq_path.write_bytes(fastq_path.read_bytes()[:-8])
    changed = _run_manifest_validator(metadata_path, output_dir)

    assert changed.returncode != 0
    assert "content contract changed" in changed.stderr or "complete gzip/FASTQ" in changed.stderr


@pytest.mark.parametrize("validator", ["full", "index"])
@pytest.mark.parametrize("invalid_contract", ["schema2", "missing_source", "mixed_source"])
def test_completion_manifest_requires_schema3_read_source(
    monkeypatch,
    tmp_path,
    validator,
    invalid_contract,
):
    metadata_path = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "getfastq"
    _metadata(metadata_path, ["RUN1"])
    fastq_url = "https://example.invalid/RUN1.fastq"
    responses = {
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=RUN1": _xml(
            "RUN1.fastq", fastq_url
        ),
        fastq_url: b"@read\nACGT\n+\n!!!!\n",
    }
    _run_fallback(monkeypatch, metadata_path, output_dir, responses)
    manifest_path = output_dir / "getfastq_completion.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if invalid_contract == "schema2":
        manifest["schema_version"] = 2
    elif invalid_contract == "missing_source":
        del manifest["read_source"]
    else:
        manifest["read_source"] = "mixed"
    manifest_path.write_text(json.dumps(manifest) + "\n", encoding="utf-8")

    if validator == "full":
        completed = _run_manifest_validator(metadata_path, output_dir)
    else:
        completed = _run_manifest_index_validator(metadata_path, output_dir)

    assert completed.returncode != 0


@pytest.mark.parametrize("validator", ["full", "index"])
def test_completion_manifest_rejects_unreferenced_fastq(
    monkeypatch,
    tmp_path,
    validator,
):
    metadata_path = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "getfastq"
    _metadata(metadata_path, ["RUN1"])
    fastq_url = "https://example.invalid/RUN1.fastq"
    responses = {
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=RUN1": _xml(
            "RUN1.fastq", fastq_url
        ),
        fastq_url: b"@read\nACGT\n+\n!!!!\n",
    }
    _run_fallback(monkeypatch, metadata_path, output_dir, responses)
    stale = output_dir / "STALE" / "STALE.amalgkit.fastq.gz"
    stale.parent.mkdir()
    with gzip.open(stale, "wb") as handle:
        handle.write(b"@stale\nTGCA\n+\n!!!!\n")

    if validator == "full":
        completed = _run_manifest_validator(metadata_path, output_dir)
    else:
        completed = _run_manifest_index_validator(metadata_path, output_dir)

    assert completed.returncode != 0
    assert "FASTQ set differs from files on disk" in completed.stderr


def test_public_fallback_quarantines_fastq_from_run_absent_in_current_metadata(
    monkeypatch,
    tmp_path,
):
    metadata_path = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "getfastq"
    _metadata(metadata_path, ["RUN1"])
    stale = output_dir / "STALE" / "STALE.amalgkit.fastq.gz"
    stale.parent.mkdir(parents=True)
    with gzip.open(stale, "wb") as handle:
        handle.write(b"@stale\nTGCA\n+\n!!!!\n")
    fastq_url = "https://example.invalid/RUN1.fastq"
    responses = {
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=RUN1": _xml(
            "RUN1.fastq", fastq_url
        ),
        fastq_url: b"@read\nACGT\n+\n!!!!\n",
    }

    _run_fallback(monkeypatch, metadata_path, output_dir, responses)

    assert not stale.exists()
    assert stale.with_name(stale.name + ".not_in_completion_manifest").is_file()
    manifest = json.loads(
        (output_dir / "getfastq_completion.json").read_text(encoding="utf-8")
    )
    assert [entry["run"] for entry in manifest["runs"]] == ["RUN1"]
    completed = _run_manifest_validator(metadata_path, output_dir)
    assert completed.returncode == 0, completed.stderr


def test_skip_index_validator_rejects_truncated_bound_fastq_without_full_content_scan(
    monkeypatch,
    tmp_path,
):
    metadata_path = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "getfastq"
    _metadata(metadata_path, ["RUN1"])
    fastq_url = "https://example.invalid/RUN1.fastq"
    responses = {
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=RUN1": _xml(
            "RUN1.fastq", fastq_url
        ),
        fastq_url: b"@read\nACGT\n+\n!!!!\n",
    }
    _run_fallback(monkeypatch, metadata_path, output_dir, responses)

    current = _run_manifest_index_validator(metadata_path, output_dir)
    assert current.returncode == 0, current.stderr

    fastq_path = output_dir / "RUN1" / "RUN1.amalgkit.fastq.gz"
    fastq_path.write_bytes(fastq_path.read_bytes()[:-8])
    changed = _run_manifest_index_validator(metadata_path, output_dir)

    assert changed.returncode != 0
    assert "FASTQ size changed" in changed.stderr
    source = _shell_function_source("validate_amalgkit_getfastq_completion_manifest_index")
    assert 'path.open("rb")' not in source
    assert "gzip.open" not in source
    assert "hashlib" not in source


def test_drifted_bound_fastq_and_manifest_are_preserved_for_recovery(monkeypatch, tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "getfastq"
    _metadata(metadata_path, ["RUN1"])
    fastq_url = "https://example.invalid/RUN1.fastq"
    responses = {
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=RUN1": _xml(
            "RUN1.fastq", fastq_url
        ),
        fastq_url: b"@read\nACGT\n+\n!!!!\n",
    }
    _run_fallback(monkeypatch, metadata_path, output_dir, responses)
    fastq_path = output_dir / "RUN1" / "RUN1.amalgkit.fastq.gz"
    truncated = fastq_path.read_bytes()[:-8]
    fastq_path.write_bytes(truncated)

    recovered = _run_contract_recovery(output_dir)

    assert recovered.returncode == 0, recovered.stderr
    assert not fastq_path.exists()
    assert (fastq_path.with_name(fastq_path.name + ".pre_integrity_recovery")).read_bytes() == truncated
    assert not (output_dir / "getfastq_completion.json").exists()
    assert (output_dir / "getfastq_completion.json.pre_integrity_recovery").is_file()


def test_amalgkit_manifest_is_bound_to_validated_fastq_bytes(tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "getfastq"
    run_dir = output_dir / "RUN1"
    run_dir.mkdir(parents=True)
    _metadata(metadata_path, ["RUN1"])
    fastq_path = run_dir / "RUN1.amalgkit.fastq.gz"
    with gzip.open(fastq_path, "wb") as handle:
        handle.write(b"@read\nACGT\n+\n!!!!\n")
    stale_path = output_dir / "STALE" / "STALE.amalgkit.fastq.gz"
    stale_path.parent.mkdir()
    with gzip.open(stale_path, "wb") as handle:
        handle.write(b"@stale\nTGCA\n+\n!!!!\n")
    (output_dir / "getfastq_completion.json").write_text(
        json.dumps(
            {
                "status": "complete",
                "run_count": 1,
                "runs": [
                    {
                        "run": "RUN1",
                        "status": "complete",
                        "files": ["RUN1/RUN1.amalgkit.fastq.gz"],
                    }
                ],
            }
        )
        + "\n",
        encoding="utf-8",
    )

    bound = _run_manifest_binder(metadata_path, output_dir)

    assert bound.returncode == 0, bound.stderr
    manifest = json.loads((output_dir / "getfastq_completion.json").read_text(encoding="utf-8"))
    contract = manifest["runs"][0]["files"][0]
    assert manifest["schema_version"] == 3
    assert manifest["read_source"] == "amalgkit"
    assert contract == {
        "path": "RUN1/RUN1.amalgkit.fastq.gz",
        "size": fastq_path.stat().st_size,
        "sha256": hashlib.sha256(fastq_path.read_bytes()).hexdigest(),
    }
    assert not stale_path.exists()
    assert stale_path.with_name(
        stale_path.name + ".not_in_completion_manifest"
    ).is_file()
    validated = _run_manifest_validator(metadata_path, output_dir)
    assert validated.returncode == 0, validated.stderr


def test_downstream_fastq_gates_revalidate_the_completion_contract():
    text = CORE_PATH.read_text(encoding="utf-8")

    assert "Published getfastq completion index changed; routing the exact run through bounded recovery" in text
    assert "Published getfastq completion contract changed before resume staging; preserving drifted bytes" in text
    assert "getfastq completion contract changed before transcriptome assembly" in text
    assert "getfastq completion contract changed before Corset read reuse" in text
    assert "getfastq completion contract changed before quant FASTQ reuse" in text
    once_body = _shell_function_source(
        "validate_amalgkit_getfastq_completion_manifest_once"
    )
    attempt_body = _shell_function_source("run_amalgkit_getfastq_attempt")
    assert "validate_amalgkit_getfastq_completion_manifest_index" in once_body
    assert "validate_amalgkit_getfastq_completion_manifest \\" in once_body
    assert "validate_amalgkit_getfastq_completion_manifest \\" not in attempt_body
    assert "mark_amalgkit_getfastq_content_validated" in attempt_body
    assert text.count("! validate_amalgkit_getfastq_completion_manifest_once \\") == 4
    prepare = text.index("gg_artifact_prepare_stage getfastq_needs_update run_amalgkit_getfastq")
    skip_gate = text.index("validate_amalgkit_getfastq_completion_manifest_index", prepare)
    recovery = text.index("prepare_amalgkit_getfastq_contract_recovery", skip_gate)
    assembly = text.index("task='De novo transcriptome assembly'", recovery)
    assembly_gate = text.index("getfastq completion contract changed before transcriptome assembly", assembly)
    first_assembly_fastq_read = text.index("seqkit sample", assembly)
    assert prepare < skip_gate < recovery < assembly < assembly_gate < first_assembly_fastq_read


def test_content_validator_reuses_stat_identity_and_rescans_after_change(tmp_path):
    trace_path = tmp_path / "trace.txt"
    script = "\n".join(
        [
            "set -u",
            _shell_function_source("mark_amalgkit_getfastq_content_validated"),
            _shell_function_source(
                "validate_amalgkit_getfastq_completion_manifest_once"
            ),
            r'''
trace_path=$1
mock_fingerprint=stable
getfastq_content_validated=0
getfastq_content_validation_fingerprint=""
validate_amalgkit_getfastq_completion_manifest() {
  printf 'full\n' >> "${trace_path}"
}
validate_amalgkit_getfastq_completion_manifest_index() {
  printf 'index\n' >> "${trace_path}"
  printf '%s\n' "${mock_fingerprint}"
}
validate_amalgkit_getfastq_completion_manifest_once manifest metadata
validate_amalgkit_getfastq_completion_manifest_once manifest metadata
mock_fingerprint=changed
validate_amalgkit_getfastq_completion_manifest_once manifest metadata
''',
        ]
    )

    completed = subprocess.run(
        ["bash", "-s", "--", str(trace_path)],
        input=script,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr + completed.stdout
    assert trace_path.read_text(encoding="utf-8").splitlines() == [
        "full",
        "index",
        "index",
        "index",
        "full",
        "index",
    ]
    assert "filesystem identity changed" in completed.stderr


def test_content_validator_rejects_same_size_fastq_mutation(monkeypatch, tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "getfastq"
    _metadata(metadata_path, ["RUN1"])
    fastq_url = "https://example.invalid/RUN1.fastq"
    responses = {
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=RUN1": _xml(
            "RUN1.fastq", fastq_url
        ),
        fastq_url: b"@read\nACGTACGT\n+\n!!!!!!!!\n",
    }
    _run_fallback(monkeypatch, metadata_path, output_dir, responses)
    script = "\n".join(
        [
            "set -u",
            _shell_function_source("validate_amalgkit_getfastq_completion_manifest"),
            _shell_function_source("validate_amalgkit_getfastq_completion_manifest_index"),
            _shell_function_source("mark_amalgkit_getfastq_content_validated"),
            _shell_function_source(
                "validate_amalgkit_getfastq_completion_manifest_once"
            ),
            r'''
getfastq_content_validated=0
getfastq_content_validation_fingerprint=""
validate_amalgkit_getfastq_completion_manifest_once "$1/getfastq_completion.json" "$2"
python - "$1/RUN1/RUN1.amalgkit.fastq.gz" <<'PY'
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
payload = bytearray(path.read_bytes())
payload[len(payload) // 2] ^= 1
path.write_bytes(payload)
PY
validate_amalgkit_getfastq_completion_manifest_once "$1/getfastq_completion.json" "$2"
''',
        ]
    )

    completed = subprocess.run(
        ["bash", "-s", "--", str(output_dir), str(metadata_path)],
        input=script,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode != 0
    assert "filesystem identity changed" in completed.stderr
    assert "content contract changed" in completed.stderr


def test_public_fallback_rejects_ena_checksum_mismatch(monkeypatch, tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "getfastq"
    _metadata(metadata_path, ["ERR123456"])

    fastq_payload = gzip.compress(b"@ena\nACGT\n+\n!!!!\n")
    fastq_url = "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR123456.fastq.gz"
    responses = {
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=ERR123456": b"<ROOT />",
        _ena_report_url("ERR123456"): json.dumps(
            [
                {
                    "run_accession": "ERR123456",
                    "fastq_ftp": fastq_url.removeprefix("https://"),
                    "fastq_md5": "0" * 32,
                }
            ]
        ).encode(),
        fastq_url: fastq_payload,
    }

    attempts = {}
    with pytest.raises(SystemExit, match="failed bounded integrity retries") as exc_info:
        _run_fallback(monkeypatch, metadata_path, output_dir, responses, attempts)

    assert "checksum mismatch expected={}".format("0" * 32) in str(exc_info.value)
    assert "actual={}".format(hashlib.md5(fastq_payload).hexdigest()) in str(exc_info.value)
    assert "bytes={}".format(len(fastq_payload)) in str(exc_info.value)
    assert attempts[fastq_url] == 5
    assert not (output_dir / "getfastq_completion.json").exists()
    assert not _partial_files(output_dir)


@pytest.mark.parametrize("first_failure", ["checksum", "gzip"])
def test_public_fallback_retries_integrity_failure_before_atomic_publication(
    monkeypatch,
    tmp_path,
    first_failure,
):
    metadata_path = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "getfastq"
    _metadata(metadata_path, ["ERR123456"])

    valid_payload = gzip.compress(b"@ena\nACGT\n+\n!!!!\n")
    if first_failure == "checksum":
        first_payload = gzip.compress(b"@wrong\nTGCA\n+\n!!!!\n")
    else:
        first_payload = valid_payload[:-8]
    fastq_url = "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR123456.fastq.gz"
    responses = {
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=ERR123456": b"<ROOT />",
        _ena_report_url("ERR123456"): json.dumps(
            [
                {
                    "run_accession": "ERR123456",
                    "fastq_ftp": fastq_url.removeprefix("https://"),
                    "fastq_md5": hashlib.md5(valid_payload).hexdigest(),
                }
            ]
        ).encode(),
        fastq_url: [first_payload, valid_payload],
    }

    attempts = _run_fallback(monkeypatch, metadata_path, output_dir, responses)

    recovered = output_dir / "ERR123456" / "ERR123456.amalgkit.fastq.gz"
    assert recovered.read_bytes() == valid_payload
    assert attempts[fastq_url] == 2
    assert not _partial_files(output_dir)
    validated = _run_manifest_validator(metadata_path, output_dir)
    assert validated.returncode == 0, validated.stderr


def test_public_fallback_preserves_and_recovers_truncated_existing_fastq(monkeypatch, tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "getfastq"
    _metadata(metadata_path, ["RUN1"])
    existing = output_dir / "RUN1" / "RUN1.amalgkit.fastq.gz"
    existing.parent.mkdir(parents=True)
    truncated = gzip.compress(b"@read\nACGT\n+\n!!!!\n" * 1000)[:-8]
    existing.write_bytes(truncated)
    valid_payload = gzip.compress(b"@replacement\nTGCA\n+\n!!!!\n")
    fastq_url = "https://example.invalid/RUN1.fastq.gz"
    responses = {
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=RUN1": _xml(
            "RUN1.fastq.gz", fastq_url
        ),
        fastq_url: valid_payload,
    }

    attempts = _run_fallback(monkeypatch, metadata_path, output_dir, responses)

    preserved = existing.with_name(existing.name + ".pre_integrity_recovery")
    assert preserved.read_bytes() == truncated
    assert existing.read_bytes() == valid_payload
    assert attempts[fastq_url] == 1
    assert not _partial_files(output_dir)
    validated = _run_manifest_validator(metadata_path, output_dir)
    assert validated.returncode == 0, validated.stderr


def test_public_fallback_rejects_non_fastq_success_response(monkeypatch, tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "getfastq"
    _metadata(metadata_path, ["RUN1"])
    fastq_url = "https://example.invalid/RUN1.fastq"
    responses = {
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=RUN1": _xml(
            "RUN1.fastq", fastq_url
        ),
        fastq_url: b"<html>temporary upstream error</html>",
    }

    attempts = {}
    with pytest.raises(SystemExit, match="failed bounded integrity retries") as exc_info:
        _run_fallback(monkeypatch, metadata_path, output_dir, responses, attempts)

    assert "incomplete FASTQ gzip" in str(exc_info.value)
    assert attempts[fastq_url] == 5
    assert not (output_dir / "RUN1" / "RUN1.amalgkit.fastq.gz").exists()
    assert not (output_dir / "getfastq_completion.json").exists()
    assert not _partial_files(output_dir)


def _run_resume_staging(work_dir: Path, published_dir: Path) -> subprocess.CompletedProcess[str]:
    script = "\n".join(
        [
            "set -u",
            _shell_function_source("stage_getfastq_outputs_for_resume"),
            r'''
dir_tmp=$1
dir_amalgkit_getfastq_sp=$2
ensure_dir() { mkdir -p -- "$1"; }
ensure_parent_dir() { mkdir -p -- "$(dirname -- "$1")"; }
stage_getfastq_outputs_for_resume
''',
        ]
    )
    return subprocess.run(
        ["bash", "-s", "--", str(work_dir), str(published_dir)],
        input=script,
        text=True,
        capture_output=True,
        check=False,
    )


def test_resume_staging_materializes_exact_published_symlink(tmp_path):
    work_dir = tmp_path / "work"
    published_dir = tmp_path / "published"
    run_dir = published_dir / "RUN1"
    run_dir.mkdir(parents=True)
    state = run_dir / "getfastq_run_state.json"
    state.write_text('{"status":"validated"}\n', encoding="utf-8")
    expected_state = state.read_bytes()
    work_dir.mkdir()
    (work_dir / "getfastq").symlink_to(published_dir, target_is_directory=True)

    completed = _run_resume_staging(work_dir, published_dir)

    assert completed.returncode == 0, completed.stderr + completed.stdout
    resume_dir = work_dir / "getfastq"
    assert resume_dir.is_dir()
    assert not resume_dir.is_symlink()
    assert (resume_dir / "RUN1" / state.name).read_bytes() == expected_state
    assert "Materialized the exact published getfastq symlink" in completed.stdout


def test_resume_staging_rejects_symlink_outside_published_root(tmp_path):
    work_dir = tmp_path / "work"
    published_dir = tmp_path / "published"
    outside_dir = tmp_path / "outside"
    published_dir.mkdir()
    outside_dir.mkdir()
    marker = outside_dir / "keep.txt"
    marker.write_text("keep\n", encoding="utf-8")
    work_dir.mkdir()
    resume_link = work_dir / "getfastq"
    resume_link.symlink_to(outside_dir, target_is_directory=True)

    completed = _run_resume_staging(work_dir, published_dir)

    assert completed.returncode != 0
    assert resume_link.is_symlink()
    assert marker.read_text(encoding="utf-8") == "keep\n"
    assert "does not resolve to the exact published output directory" in completed.stderr


def test_fatal_then_incomplete_retry_routes_through_public_fallback(tmp_path):
    trace_path = tmp_path / "trace.txt"
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    script = "\n".join(
        [
            "set -u",
            _shell_function_source("run_amalgkit_getfastq_or_fallback"),
            r'''
trace_path=$1
dir_tmp=$2
file_amalgkit_metadata="${dir_tmp}/metadata.tsv"
dir_amalgkit_getfastq_sp="${dir_tmp}/published"
amalgkit_rrna_filter=yes
attempt_count=0

run_amalgkit_getfastq_attempt() {
  attempt_count=$((attempt_count + 1))
  printf 'attempt:%s:%s\n' "$1" "$2" >> "${trace_path}"
  if [[ ${attempt_count} -eq 1 ]]; then
    return 2
  fi
  return 3
}
has_resumable_getfastq_run_state() {
  printf 'resume-state-present\n' >> "${trace_path}"
  return 0
}
prepare_getfastq_outputs_for_public_fallback() {
  printf 'prepare\n' >> "${trace_path}"
}
download_public_original_fastqs_for_metadata() {
  printf 'download:%s:%s\n' "$1" "$2" >> "${trace_path}"
  return 0
}
validate_amalgkit_getfastq_completion_manifest_index() {
  printf 'prevalidate:%s:%s\n' "$1" "$2" >> "${trace_path}"
  return 0
}
mark_amalgkit_getfastq_content_validated() {
  printf 'validate:%s:%s\n' "$1" "$2" >> "${trace_path}"
  return 0
}
mv_out_replace_dir() {
  printf 'publish:%s:%s\n' "$1" "$2" >> "${trace_path}"
}
rm() {
  printf 'cleanup:%s\n' "$*" >> "${trace_path}"
}

run_amalgkit_getfastq_or_fallback
''',
        ]
    )

    completed = subprocess.run(
        ["bash", "-s", "--", str(trace_path), str(work_dir)],
        input=script,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr + completed.stdout
    assert trace_path.read_text(encoding="utf-8").splitlines() == [
        "attempt:yes:initial",
        "attempt:no:retry_rrna_filter_no",
        "prepare",
        f"download:{work_dir / 'metadata.tsv'}:{work_dir / 'getfastq_public_original'}",
        f"prevalidate:{work_dir / 'getfastq_public_original/getfastq_completion.json'}:{work_dir / 'metadata.tsv'}",
        f"publish:{work_dir / 'getfastq_public_original'}:{work_dir / 'published'}",
        f"validate:{work_dir / 'published/getfastq_completion.json'}:{work_dir / 'metadata.tsv'}",
        f"cleanup:-rf -- {work_dir / 'getfastq_public_original'}",
    ]
    assert "incomplete all-run manifest" in completed.stdout


def _run_fatal_route_fixture(tmp_path: Path, final_log_lines: list[str]) -> subprocess.CompletedProcess[str]:
    trace_path = tmp_path / "trace.txt"
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    script = "\n".join(
        [
            "set -u",
            _shell_function_source("amalgkit_getfastq_log_has_only_download_source_exhaustion"),
            _shell_function_source("run_amalgkit_getfastq_or_fallback"),
            r'''
trace_path=$1
dir_tmp=$2
final_log_lines=$3
file_amalgkit_metadata="${dir_tmp}/metadata.tsv"
dir_amalgkit_getfastq_sp="${dir_tmp}/published"
amalgkit_rrna_filter=yes
attempt_count=0

run_amalgkit_getfastq_attempt() {
  attempt_count=$((attempt_count + 1))
  local log_file="${dir_tmp}/amalgkit_getfastq.$2.log"
  printf 'attempt:%s:%s\n' "$1" "$2" >> "${trace_path}"
  if [[ ${attempt_count} -eq 2 ]]; then
    printf '%s\n' "${final_log_lines}" > "${log_file}"
  fi
  return 2
}
has_resumable_getfastq_run_state() {
  printf 'resume-state-present\n' >> "${trace_path}"
  return 0
}
prepare_getfastq_outputs_for_public_fallback() {
  printf 'prepare\n' >> "${trace_path}"
}
download_public_original_fastqs_for_metadata() {
  printf 'download\n' >> "${trace_path}"
  return 0
}
validate_amalgkit_getfastq_completion_manifest_index() {
  printf 'prevalidate\n' >> "${trace_path}"
  return 0
}
mark_amalgkit_getfastq_content_validated() {
  printf 'validate\n' >> "${trace_path}"
  return 0
}
mv_out_replace_dir() {
  printf 'publish\n' >> "${trace_path}"
}
rm() {
  printf 'cleanup\n' >> "${trace_path}"
}

run_amalgkit_getfastq_or_fallback
''',
        ]
    )
    return subprocess.run(
        ["bash", "-s", "--", str(trace_path), str(work_dir), "\n".join(final_log_lines)],
        input=script,
        text=True,
        capture_output=True,
        check=False,
    )


def test_download_source_exhaustion_only_routes_through_public_fallback(tmp_path):
    completed = _run_fatal_route_fixture(
        tmp_path,
        ["ERROR: Configured download sources were exhausted."],
    )

    assert completed.returncode == 0, completed.stderr + completed.stdout
    trace = (tmp_path / "trace.txt").read_text(encoding="utf-8").splitlines()
    assert trace == [
        "attempt:yes:initial",
        "attempt:no:retry_rrna_filter_no",
        "prepare",
        "download",
        "prevalidate",
        "publish",
        "validate",
        "cleanup",
    ]
    assert "Every fatal condition" in completed.stdout


def test_aggregate_download_source_exhaustion_routes_through_public_fallback(tmp_path):
    completed = _run_fatal_route_fixture(
        tmp_path,
        [
            "ERROR: getfastq failed for 1/56 SRA runs. "
            "ERR4643641: SRA file download failed for ERR4643641. "
            "Expected PATH: /tmp/ERR4643641.sra. "
            "Configured download sources were exhausted."
        ],
    )

    assert completed.returncode == 0, completed.stderr + completed.stdout
    trace = (tmp_path / "trace.txt").read_text(encoding="utf-8").splitlines()
    assert trace == [
        "attempt:yes:initial",
        "attempt:no:retry_rrna_filter_no",
        "prepare",
        "download",
        "prevalidate",
        "publish",
        "validate",
        "cleanup",
    ]
    assert "Every fatal condition" in completed.stdout


def test_download_source_exhaustion_mixed_with_another_fatal_fails_closed(tmp_path):
    completed = _run_fatal_route_fixture(
        tmp_path,
        [
            "ERROR: Configured download sources were exhausted.",
            "ERROR: Metadata validation failed.",
        ],
    )

    assert completed.returncode != 0
    trace = (tmp_path / "trace.txt").read_text(encoding="utf-8").splitlines()
    assert trace == ["attempt:yes:initial", "attempt:no:retry_rrna_filter_no"]
    assert "Exiting without fallback" in completed.stdout


def test_aggregate_exhaustion_mixed_with_another_run_failure_fails_closed(tmp_path):
    completed = _run_fatal_route_fixture(
        tmp_path,
        [
            "ERROR: getfastq failed for 2/56 SRA runs. "
            "ERR4643641: SRA file download failed for ERR4643641. "
            "Expected PATH: /tmp/ERR4643641.sra. "
            "Configured download sources were exhausted. "
            "SRR0000001: Metadata validation failed."
        ],
    )

    assert completed.returncode != 0
    trace = (tmp_path / "trace.txt").read_text(encoding="utf-8").splitlines()
    assert trace == ["attempt:yes:initial", "attempt:no:retry_rrna_filter_no"]
    assert "Exiting without fallback" in completed.stdout


def test_relaxed_metadata_normalization_preserves_old_output_when_sed_fails(tmp_path):
    output = tmp_path / "metadata.tsv"
    output.write_text("old\n", encoding="utf-8")
    metadata = tmp_path / "source.tsv"
    accessions = tmp_path / "accessions.txt"
    metadata.write_text("run\nRUN1\n", encoding="utf-8")
    accessions.write_text("RUN1\n", encoding="utf-8")
    function_source = _shell_function_source(
        "extract_transcriptomic_rows_for_requested_accessions"
    )
    script = "\n".join(
        [
            "set -u",
            function_source,
            r'''
gg_support_dir=.
python() {
  local raw_output=$5
  printf 'new\n' > "${raw_output}"
}
sed() {
  return 7
}
mv_out() {
  command mv -- "$1" "$2"
}
extract_transcriptomic_rows_for_requested_accessions "$1" "$2" "$3"
''',
        ]
    )

    completed = subprocess.run(
        ["bash", "-s", "--", str(metadata), str(accessions), str(output)],
        input=script,
        text=True,
        capture_output=True,
        check=False,
        cwd=tmp_path,
    )

    assert completed.returncode != 0
    assert output.read_text(encoding="utf-8") == "old\n"
    assert list(tmp_path.glob("metadata.tsv.raw.*"))
    assert not list(tmp_path.glob("metadata.tsv.normalized.*"))

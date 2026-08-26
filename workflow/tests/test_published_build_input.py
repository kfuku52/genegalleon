import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "container/scripts/check_published_build_input.py"
AMD64_HASH = "a" * 64
ARM64_HASH = "b" * 64


def run_check(payload, *expected):
    arguments = [sys.executable, str(SCRIPT)]
    for value in expected:
        arguments.extend(("--expected", value))
    return subprocess.run(
        arguments,
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        check=False,
    )


def image_payload(amd64_hash=AMD64_HASH, arm64_hash=ARM64_HASH, revision="c" * 40):
    return {
        "linux/amd64": {
            "config": {
                "Labels": {
                    "io.genegalleon.build-input": amd64_hash,
                    "org.opencontainers.image.revision": revision,
                }
            }
        },
        "linux/arm64": {
            "config": {
                "Labels": {
                    "io.genegalleon.build-input": arm64_hash,
                    "org.opencontainers.image.revision": revision,
                }
            }
        },
    }


def test_published_build_inputs_match_both_platforms():
    completed = run_check(
        image_payload(),
        f"linux/amd64={AMD64_HASH}",
        f"linux/arm64={ARM64_HASH}",
    )

    assert completed.returncode == 0, completed.stderr
    assert "linux/amd64, linux/arm64" in completed.stderr


def test_changed_upstream_build_input_requires_rebuild():
    completed = run_check(
        image_payload(amd64_hash="c" * 64),
        f"linux/amd64={AMD64_HASH}",
        f"linux/arm64={ARM64_HASH}",
    )

    assert completed.returncode == 1
    assert "linux/amd64 build input changed" in completed.stderr
    assert f"published={'c' * 64}" in completed.stderr
    assert f"expected={AMD64_HASH}" in completed.stderr


def test_missing_platform_or_label_requires_rebuild():
    payload = image_payload()
    del payload["linux/arm64"]
    payload["linux/amd64"]["config"]["Labels"].clear()

    completed = run_check(
        payload,
        f"linux/amd64={AMD64_HASH}",
        f"linux/arm64={ARM64_HASH}",
    )

    assert completed.returncode == 1
    assert "lacks a valid io.genegalleon.build-input label" in completed.stderr
    assert "missing platform linux/arm64" in completed.stderr


def test_invalid_published_json_fails_closed():
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--expected",
            f"linux/amd64={AMD64_HASH}",
        ],
        input="not-json",
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 1
    assert "comparison failed" in completed.stderr


def test_prints_common_published_revision_for_both_platforms():
    revision = "d" * 40
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--print-common-revision",
            "linux/amd64",
            "--print-common-revision",
            "linux/arm64",
        ],
        input=json.dumps(image_payload(revision=revision)),
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == revision


def test_rejects_different_published_revisions_across_platforms():
    payload = image_payload(revision="d" * 40)
    payload["linux/arm64"]["config"]["Labels"]["org.opencontainers.image.revision"] = "e" * 40

    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--print-common-revision",
            "linux/amd64",
            "--print-common-revision",
            "linux/arm64",
        ],
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 1
    assert "revisions differ" in completed.stderr


def test_rejects_missing_published_revision():
    payload = image_payload()
    del payload["linux/arm64"]["config"]["Labels"]["org.opencontainers.image.revision"]

    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--print-common-revision",
            "linux/amd64",
            "--print-common-revision",
            "linux/arm64",
        ],
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 1
    assert "lacks a valid org.opencontainers.image.revision label" in completed.stderr

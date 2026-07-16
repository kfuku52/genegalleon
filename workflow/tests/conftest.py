import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


@pytest.fixture(autouse=True)
def isolate_test_runtime(monkeypatch, tmp_path):
    monkeypatch.setenv("GG_TEST_ALLOW_ONLY_LOOPBACK_HTTP", "1")
    monkeypatch.setenv(
        "GG_INPUT_GENERATION_OUTPUT_ROOT",
        str(tmp_path / "input_generation_output"),
    )

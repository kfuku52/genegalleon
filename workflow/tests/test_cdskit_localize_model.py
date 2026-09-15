import hashlib
import io
import json
import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from workflow.support import resolve_cdskit_localize_model as model

PAYLOAD = b"verified checkpoint"


def release(date="2026-09-11", name="new", **kwargs):
    return dict(
        draft=False,
        prerelease=False,
        tag_name=f"localize-{name}",
        published_at=date,
        assets=[
            dict(
                name=f"cdskit-localize-{name}.pt",
                digest="sha256:" + hashlib.sha256(PAYLOAD).hexdigest(),
                browser_download_url=model.DOWNLOAD_PREFIX + name + "/model.pt",
            )
        ],
        **kwargs,
    )


def test_first_download_selects_latest_then_reuses_offline(tmp_path):
    calls = []

    def fetch(url):
        calls.append(url)
        if url.startswith(model.RELEASES_URL):
            return io.BytesIO(json.dumps([release("2026-01-01", "old"), release()]).encode())
        return io.BytesIO(PAYLOAD)

    with patch.object(model, "open_url", fetch):
        path = model.resolve(tmp_path)
    assert path.name == "cdskit-localize-new.pt"
    assert path.read_bytes() == PAYLOAD
    assert len(calls) == 2
    with patch.object(model, "open_url", side_effect=AssertionError("Must not check for updates")):
        assert model.resolve(tmp_path) == path
        assert model.resolve(tmp_path, allow_download=False) == path
    spec = json.loads((path.parent / "selection.json").read_text())
    assert spec["release"] == "localize-new"


def test_failed_download_does_not_pin_selection(tmp_path):
    with patch.object(
        model, "open_url", side_effect=[io.BytesIO(json.dumps([release()]).encode()), io.BytesIO(b"corrupt")]
    ):
        with pytest.raises(ValueError, match="checksum"):
            model.resolve(tmp_path)
    assert not (tmp_path / "genegalleon-latest" / "selection.json").exists()


def test_offline_first_download_fails_without_network(tmp_path, monkeypatch):
    with patch.object(model, "open_url", side_effect=AssertionError("Must not download")):
        with pytest.raises(FileNotFoundError, match="disabled"):
            model.resolve(tmp_path, False)
        monkeypatch.setenv("CDSKIT_OFFLINE", "1")
        with pytest.raises(FileNotFoundError, match="disabled"):
            model.resolve(tmp_path)


def test_corrupt_cache_is_not_automatically_replaced(tmp_path):
    with patch.object(
        model, "open_url", side_effect=[io.BytesIO(json.dumps([release()]).encode()), io.BytesIO(PAYLOAD)]
    ):
        path = model.resolve(tmp_path)
    path.write_bytes(b"corrupt")
    with patch.object(model, "open_url", side_effect=AssertionError("Must not update")):
        with pytest.raises(ValueError, match="checksum"):
            model.resolve(tmp_path)


def test_draft_prerelease_and_non_model_releases_are_excluded():
    excluded = [release("2027", "draft"), release("2028", "preview"), release("2029", "program")]
    excluded[0]["draft"] = True
    excluded[1]["prerelease"] = True
    excluded[2]["tag_name"] = "v9.0.0"
    with patch.object(model, "open_url", return_value=io.BytesIO(json.dumps(excluded + [release()]).encode())):
        assert model.latest_model()["release"] == "localize-new"


def test_shell_helper_uses_saved_path_and_preserves_explicit_model(tmp_path):
    repo = Path(__file__).resolve().parents[2]
    cache = tmp_path / "models"
    with patch.object(
        model, "open_url", side_effect=[io.BytesIO(json.dumps([release()]).encode()), io.BytesIO(PAYLOAD)]
    ):
        path = model.resolve(cache)
    script = """
source "$1/workflow/support/gg_util.sh"
export CDSKIT_MODEL_DIR="$2"
cdskit() { printf '%s\\n' "$@"; }
gg_run_cdskit_localize input protein output "$3" unknown 0 1 2 1
"""
    for request, expected in [("latest", str(path)), ("targeting5", "targeting5")]:
        result = subprocess.run(
            ["bash", "-c", script, "test", str(repo), str(cache), request], check=True, capture_output=True, text=True
        )
        args = result.stdout.splitlines()
        assert args[args.index("--model") + 1] == expected
        assert args[args.index("--model_download") + 1] == "no"


def test_concurrent_first_requests_download_once(tmp_path):
    from concurrent.futures import ThreadPoolExecutor

    calls = []

    def fetch(url):
        calls.append(url)
        if url.startswith(model.RELEASES_URL):
            return io.BytesIO(json.dumps([release()]).encode())
        return io.BytesIO(PAYLOAD)

    with patch.object(model, "open_url", fetch), ThreadPoolExecutor(max_workers=4) as pool:
        paths = list(pool.map(lambda _: model.resolve(tmp_path), range(4)))
    assert len(set(paths)) == 1
    assert len(calls) == 2


def test_release_discovery_reads_all_pages():
    first_page = [release("2026-01-01", "old") for _ in range(100)]
    with patch.object(
        model,
        "open_url",
        side_effect=[
            io.BytesIO(json.dumps(first_page).encode()),
            io.BytesIO(json.dumps([release()]).encode()),
        ],
    ) as fetch:
        assert model.latest_model()["release"] == "localize-new"
    assert fetch.call_count == 2

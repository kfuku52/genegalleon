import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

INTEGRATION_DOWNLOAD_FILES = {
    "test_format_species_inputs_download.py",
    "test_format_species_inputs_download_runtime.py",
    "test_format_species_inputs_provider_pages.py",
    "test_generate_species_trait.py",
}
INTEGRATION_WORKFLOW_FILES = {
    "test_busco_hmmsearch_wrapper.py",
    "test_genome_evolution_protein_mode.py",
    "test_gg_input_generation_end_to_end.py",
    "test_hgt_end_to_end.py",
}
RUNTIME_FILES = {
    "test_csubst_scan_runtime_integration.py",
}
SMOKE_NODE_IDS = {
    "workflow/tests/test_busco_hmmsearch_wrapper.py::test_busco_output_exists_uses_canonical_path_without_directory_rescan",
    "workflow/tests/test_busco_hmmsearch_wrapper.py::test_hmmsearch_wrapper_creates_modified_fas_symlink_when_missing",
    "workflow/tests/test_busco_hmmsearch_wrapper.py::test_gg_run_busco_with_metaeuk_modified_fas_compat_wraps_busco_env",
    "workflow/tests/test_busco_hmmsearch_wrapper.py::test_gg_busco_stderr_matches_known_metaeuk_modified_fas_bug",
    "workflow/tests/test_gg_input_generation_end_to_end.py::test_gg_input_generation_single_mode_end_to_end",
}
SMOKE_FILES = {node_id.split("::", 1)[0].rsplit("/", 1)[-1] for node_id in SMOKE_NODE_IDS}


def _test_lane(filename: str) -> str:
    if filename in RUNTIME_FILES:
        return "runtime"
    if filename in INTEGRATION_DOWNLOAD_FILES:
        return "integration-download"
    if filename in INTEGRATION_WORKFLOW_FILES:
        return "integration-workflow"
    stem = Path(filename).stem
    if "_static" in stem or stem.startswith("test_static") or filename == "test_ci_static.py":
        return "static"
    return "fast"


def pytest_addoption(parser):
    parser.addoption(
        "--gg-suite",
        action="store",
        choices=(
            "fast",
            "integration-download",
            "integration-workflow",
            "runtime",
            "smoke",
            "static",
        ),
        help="Collect one non-overlapping GeneGalleon test suite.",
    )


def pytest_ignore_collect(collection_path: Path, config):
    suite = config.getoption("--gg-suite")
    filename = collection_path.name
    if suite is None or not filename.startswith("test_") or collection_path.suffix != ".py":
        return None
    if suite == "smoke":
        return filename not in SMOKE_FILES
    return _test_lane(filename) != suite


def pytest_collection_modifyitems(config, items):
    suite = config.getoption("--gg-suite")
    selected = []
    deselected = []
    for item in items:
        lane = _test_lane(item.path.name)
        if lane.startswith("integration-"):
            item.add_marker(pytest.mark.integration)
            item.add_marker(pytest.mark.slow)
        else:
            item.add_marker(getattr(pytest.mark, lane))
            if lane == "runtime":
                item.add_marker(pytest.mark.slow)
        if item.nodeid in SMOKE_NODE_IDS:
            item.add_marker(pytest.mark.smoke)

        duplicate_smoke = suite in {"integration-download", "integration-workflow"} and item.nodeid in SMOKE_NODE_IDS
        missing_smoke = suite == "smoke" and item.nodeid not in SMOKE_NODE_IDS
        if duplicate_smoke or missing_smoke:
            deselected.append(item)
        else:
            selected.append(item)

    if deselected:
        items[:] = selected
        config.hook.pytest_deselected(items=deselected)


@pytest.fixture(autouse=True)
def isolate_test_runtime(monkeypatch, tmp_path):
    monkeypatch.setenv("GG_TEST_ALLOW_ONLY_LOOPBACK_HTTP", "1")
    monkeypatch.setenv(
        "GG_INPUT_GENERATION_OUTPUT_ROOT",
        str(tmp_path / "input_generation_output"),
    )

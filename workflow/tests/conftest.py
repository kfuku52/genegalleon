import os
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def pytest_collection_modifyitems(config, items):
    if os.environ.get("GG_RUN_NCBI_DOWNLOAD_TESTS") == "1":
        return

    ncbi_download_skip = pytest.mark.skip(
        reason="NCBI download resolver tests are skipped by default; set GG_RUN_NCBI_DOWNLOAD_TESTS=1 to run them."
    )
    ncbi_name_patterns = (
        "ncbi",
        "refseq_and_genbank_id_only_auto_resolve",
    )

    for item in items:
        item_path = getattr(item, "path", None)
        if item_path is None:
            item_path = Path(str(item.fspath))
        if item_path.name != "test_format_species_inputs_download.py":
            continue
        if any(pattern in item.name.lower() for pattern in ncbi_name_patterns):
            item.add_marker(ncbi_download_skip)

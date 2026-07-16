import csv
import sys
from pathlib import Path

import pytest

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
if str(SUPPORT_DIR) not in sys.path:
    sys.path.insert(0, str(SUPPORT_DIR))

import format_species_inputs as format_species  # noqa: E402
import format_species_manifest as manifest  # noqa: E402


def test_read_download_manifest_detects_tab_delimiter(tmp_path):
    path = tmp_path / "manifest.tsv"
    path.write_text("provider\tid\nncbi\tGCA_000001.1\n", encoding="utf-8")

    assert manifest.read_download_manifest(path) == [{"provider": "ncbi", "id": "GCA_000001.1"}]


def test_resolved_manifest_fieldnames_preserve_fernbase_confidence_order():
    fieldnames = manifest.resolved_manifest_fieldnames(
        [
            {
                "extra": "value",
                "fernbase_confidence_mode": "high-low combined",
            }
        ]
    )

    assert fieldnames.index("fernbase_confidence_mode") < fieldnames.index("extra")


def test_write_resolved_manifest_tsv_fills_missing_fields(tmp_path):
    path = tmp_path / "resolved.tsv"
    fieldnames = ["provider", "id", "species_key", "extra"]

    manifest.write_resolved_manifest_tsv(
        path,
        fieldnames,
        [{"provider": "direct", "id": "row1", "species_key": "Species_one"}],
    )

    with path.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert rows == [{"provider": "direct", "id": "row1", "species_key": "Species_one", "extra": ""}]


def test_direct_catalog_xlsx_defaults_fill_direct_rows(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    path = tmp_path / "manifest.xlsx"
    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.append(["provider", "id", "cds_url", "gff_url", "genome_url"])
    sheet.append(["direct", "Choice A", "", "", ""])
    catalog = workbook.create_sheet(manifest.DIRECT_CATALOG_XLSX_SHEET)
    catalog.append(["choice", "id", "cds_url"])
    catalog.append(["Choice A", "actual-id", "https://example.test/cds.fa.gz"])
    workbook.save(path)

    rows = manifest.read_download_manifest(path)

    assert rows[0]["id"] == "actual-id"
    assert rows[0]["cds_url"] == "https://example.test/cds.fa.gz"


def test_format_species_inputs_delegates_manifest_helpers():
    assert format_species.read_download_manifest is manifest.read_download_manifest
    assert format_species.resolved_manifest_fieldnames is manifest.resolved_manifest_fieldnames
    assert format_species.write_resolved_manifest_tsv is manifest.write_resolved_manifest_tsv

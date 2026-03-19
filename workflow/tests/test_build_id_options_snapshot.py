import importlib.util
import json
from pathlib import Path

from openpyxl import Workbook


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "build_id_options_snapshot.py"


def load_module():
    spec = importlib.util.spec_from_file_location("build_id_options_snapshot", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_build_direct_catalog_from_manifest_rows_reads_xlsx(tmp_path):
    module = load_module()
    manifest_path = tmp_path / "direct_manifest.xlsx"
    workbook = Workbook()
    sheet = workbook.active
    sheet.append(list(module.MANIFEST_FIELDNAMES))
    sheet.append(
        [
            "direct",
            "Fragaria_ananassa_kazusa",
            "Fragaria_ananassa",
            "https://example.org/fragaria.cds.fa.gz",
            "https://example.org/fragaria.gff3.gz",
            "",
            "https://example.org/fragaria.genome.fa.gz",
            "",
            "",
            "",
            "",
            "fragaria.cds.fa.gz",
            "fragaria.gff3.gz",
            "",
            "fragaria.genome.fa.gz",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
        ]
    )
    workbook.save(manifest_path)

    rows = module.load_manifest_rows(manifest_path)
    catalog = module.build_direct_catalog_from_manifest_rows(rows)
    assert len(catalog) == 1
    assert catalog[0]["id"] == "Fragaria_ananassa_kazusa"
    assert catalog[0]["species"] == "Fragaria ananassa"
    assert catalog[0]["species_key"] == "Fragaria_ananassa"
    assert catalog[0]["cds_url"] == "https://example.org/fragaria.cds.fa.gz"
    assert catalog[0]["gff_url"] == "https://example.org/fragaria.gff3.gz"
    assert catalog[0]["genome_url"] == "https://example.org/fragaria.genome.fa.gz"


def test_build_snapshot_preserves_rich_direct_catalog(monkeypatch):
    module = load_module()
    monkeypatch.setattr(module, "FETCH_PROVIDERS", ())
    direct_catalog_entries = [
        {
            "id": "Fragaria_ananassa_kazusa",
            "species": "Fragaria ananassa",
            "species_key": "Fragaria_ananassa",
            "cds_url": "https://example.org/fragaria.cds.fa.gz",
            "gff_url": "https://example.org/fragaria.gff3.gz",
            "gbff_url": "",
            "genome_url": "https://example.org/fragaria.genome.fa.gz",
            "cds_archive_member": "",
            "gff_archive_member": "",
            "gbff_archive_member": "",
            "genome_archive_member": "",
            "cds_filename": "fragaria.cds.fa.gz",
            "gff_filename": "fragaria.gff3.gz",
            "gbff_filename": "",
            "genome_filename": "fragaria.genome.fa.gz",
            "cds_url_template": "",
            "gff_url_template": "",
            "genome_url_template": "",
            "local_cds_path": "",
            "local_gff_path": "",
            "local_gbff_path": "",
            "local_genome_path": "",
        }
    ]

    snapshot, warnings = module.build_snapshot(
        timeout=0.1,
        input_root=None,
        fallback_options={},
        fallback_direct_catalog=[],
        direct_catalog_entries=direct_catalog_entries,
    )

    assert warnings == []
    assert snapshot["providers"]["direct"] == direct_catalog_entries
    assert snapshot["providers"]["direct"][0]["species_key"] == "Fragaria_ananassa"
    assert snapshot["providers"]["direct"][0]["cds_filename"] == "fragaria.cds.fa.gz"


def test_load_snapshot_direct_catalog_preserves_direct_fields(tmp_path):
    module = load_module()
    snapshot_path = tmp_path / "snapshot.json"
    snapshot_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "providers": {
                    "direct": [
                        {
                            "id": "Zoysia_matrella_kazusa",
                            "species": "Zoysia matrella",
                            "species_key": "Zoysia_matrella",
                            "genome_url": "https://example.org/zoysia.genome.fa.gz",
                            "gff_url": "https://example.org/zoysia.gff3.gz",
                            "cds_url": "https://example.org/zoysia.cds.fa.gz",
                        }
                    ]
                },
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    catalog = module.load_snapshot_direct_catalog(snapshot_path)
    assert len(catalog) == 1
    assert catalog[0]["id"] == "Zoysia_matrella_kazusa"
    assert catalog[0]["species_key"] == "Zoysia_matrella"
    assert catalog[0]["cds_url"] == "https://example.org/zoysia.cds.fa.gz"

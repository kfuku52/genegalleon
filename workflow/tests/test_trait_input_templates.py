from pathlib import Path
import csv
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
TRAIT_PLAN_PATH = REPO_ROOT / "workspace" / "input" / "input_generation" / "trait_plan.tsv"
TRAIT_SOURCES_PATH = REPO_ROOT / "workspace" / "input" / "input_generation" / "trait_database_sources.tsv"
SCRIPT_PATH = REPO_ROOT / "workflow" / "support" / "generate_species_trait.py"


def read_tsv_rows(path: Path):
    with open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    return reader.fieldnames or [], rows


def test_trait_plan_template_has_required_columns():
    assert TRAIT_PLAN_PATH.exists(), f"Missing trait plan template: {TRAIT_PLAN_PATH}"
    fieldnames, rows = read_tsv_rows(TRAIT_PLAN_PATH)
    required = {"database", "source_column", "output_trait"}
    assert required.issubset(set(fieldnames))
    assert len(rows) > 0
    for row in rows:
        assert str(row.get("database", "")).strip() != ""
        assert str(row.get("source_column", "")).strip() != ""
        assert str(row.get("output_trait", "")).strip() != ""


def test_trait_database_sources_template_has_required_columns():
    assert TRAIT_SOURCES_PATH.exists(), f"Missing trait source template: {TRAIT_SOURCES_PATH}"
    fieldnames, rows = read_tsv_rows(TRAIT_SOURCES_PATH)
    required = {"database", "acquisition_mode"}
    assert required.issubset(set(fieldnames))
    assert len(rows) > 0
    allowed_modes = {"bulk", "species_api", "gift_api"}
    for row in rows:
        mode = str(row.get("acquisition_mode", "")).strip()
        assert mode in allowed_modes


def test_trait_plan_databases_are_declared_in_sources_template():
    _, plan_rows = read_tsv_rows(TRAIT_PLAN_PATH)
    _, source_rows = read_tsv_rows(TRAIT_SOURCES_PATH)
    plan_dbs = {str(row.get("database", "")).strip().lower() for row in plan_rows}
    source_dbs = {str(row.get("database", "")).strip().lower() for row in source_rows}
    assert plan_dbs.issubset(source_dbs)


def test_trait_templates_use_supported_database_ids():
    completed = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "--print-supported-databases"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    supported = {line.split("\t", 1)[0].strip().lower() for line in lines[1:]}
    _, plan_rows = read_tsv_rows(TRAIT_PLAN_PATH)
    _, source_rows = read_tsv_rows(TRAIT_SOURCES_PATH)
    used = {
        str(row.get("database", "")).strip().lower()
        for row in list(plan_rows) + list(source_rows)
        if str(row.get("database", "")).strip() != ""
    }
    assert used.issubset(supported)


def test_trait_plan_gift_rows_have_trait_key_and_trait_key_column():
    _, plan_rows = read_tsv_rows(TRAIT_PLAN_PATH)
    gift_rows = [row for row in plan_rows if str(row.get("database", "")).strip().lower() == "gift"]
    assert len(gift_rows) > 0
    for row in gift_rows:
        assert str(row.get("trait_key", "")).strip() != ""
        assert str(row.get("trait_key_column", "")).strip() == "trait_ID"

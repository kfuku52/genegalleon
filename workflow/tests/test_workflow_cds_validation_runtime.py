"""Validate identifiers with real seqkit and GNU parallel in the owned runtime."""

import pytest
from test_workflow_audit_regressions import run_shell


@pytest.mark.parametrize("headers", [
    ["Alpha_one_valid", "Beta_two_foreign"],
    ["Beta_two_foreign", "Alpha_one_valid"],
    ["Alpha_one_valid", "Alpha_one"],
    ["Alpha_one_valid", "Alpha_one_"],
    [],
])
def test_cds_validation_rejects_every_foreign_or_incomplete_identifier(tmp_path, headers):
    directory = tmp_path / "cds"
    directory.mkdir()
    (directory / "Alpha_one.cds.fa").write_text("".join(f">{header}\nATGAAA\n" for header in headers))
    result = run_shell('GG_TASK_CPUS=1\ncheck_species_cds_dir "$2"', directory, timeout=30)
    assert result.returncode != 0, result.stdout + result.stderr


def test_cds_validation_accepts_all_valid_identifiers_and_reports_reader_failure(tmp_path):
    directory = tmp_path / "cds"
    directory.mkdir()
    (directory / "Alpha_one.cds.fa").write_text(">Alpha_one_a\nATGAAA\n>Alpha_one_b\nATGCCC\n")
    result = run_shell('GG_TASK_CPUS=1\ncheck_species_cds_dir "$2"', directory, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    result = run_shell('GG_TASK_CPUS=1\nseqkit() { return 23; }\nexport -f seqkit\ncheck_species_cds_dir "$2"', directory, timeout=30)
    assert result.returncode != 0
    assert "Failed to read CDS sequence names" in result.stdout

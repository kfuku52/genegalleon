from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SUPPORT_DIR = REPO_ROOT / "workflow" / "support"


def test_format_species_inputs_delegates_writer_helpers():
    script_text = (SUPPORT_DIR / "format_species_inputs.py").read_text(encoding="utf-8")

    assert "from format_species_writers import (" in script_text
    assert "def write_fasta_records_gzip(" not in script_text
    assert "def write_gff_gzip(" not in script_text
    assert "def write_text_lines(" not in script_text


def test_format_species_inputs_delegates_stage_and_provider_input_helpers():
    script_text = (SUPPORT_DIR / "format_species_inputs.py").read_text(encoding="utf-8")

    assert "from format_species_download_stage import apply_download_input_dir, run_download_stage" in script_text
    assert "from format_species_provider_inputs import (" in script_text
    assert "manifest_declared_providers," in script_text
    assert "manifest_declared_species_keys," in script_text
    assert "resolve_provider_inputs," in script_text
    assert "def resolve_provider_inputs(" not in script_text
    assert "def manifest_declared_providers(" not in script_text
    assert "def manifest_declared_species_keys(" not in script_text

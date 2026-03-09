from pathlib import Path


WORKFLOW_DIR = Path(__file__).resolve().parents[1]


def test_gg_input_generation_core_validates_configured_species_cds_dir():
    core = WORKFLOW_DIR / "core" / "gg_input_generation_core.sh"
    text = core.read_text(encoding="utf-8")
    assert 'check_species_cds_dir "${species_cds_dir}"' in text
    assert 'check_species_cds "${gg_workspace_dir}"' not in text

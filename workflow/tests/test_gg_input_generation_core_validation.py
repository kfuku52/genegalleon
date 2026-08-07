from pathlib import Path

WORKFLOW_DIR = Path(__file__).resolve().parents[1]


def test_gg_input_generation_core_validates_configured_species_cds_dir():
    core = WORKFLOW_DIR / "core" / "gg_input_generation_core.sh"
    text = core.read_text(encoding="utf-8")
    assert 'check_species_cds_dir "${species_cds_dir}"' in text
    assert 'check_species_cds "${gg_workspace_dir}"' not in text
    assert 'gene_grouping_mode="${gene_grouping_mode:-rescue_overlap}"' in text
    assert '--gene-grouping-mode "${gene_grouping_mode}"' in text
    assert 'gff_repair_mode="${gff_repair_mode:-safe}"' in text
    assert '--gff-repair-mode "${gff_repair_mode}"' in text
    assert 'validate_longest_cds_selection.py' in text


def test_array_worker_validates_before_fx2tab_and_busco():
    core = WORKFLOW_DIR / "core" / "gg_input_generation_core.sh"
    text = core.read_text(encoding="utf-8")
    worker_start = text.index("run_array_worker_mode()")
    worker_end = text.index("run_array_finalize_mode()", worker_start)
    worker = text[worker_start:worker_end]

    assert worker.index("run_validate_stage_one_worker") < worker.index("run_cds_fx2tab_stage_one_worker")
    assert worker.index("run_validate_stage_one_worker") < worker.index("run_species_busco_stage_one_worker")

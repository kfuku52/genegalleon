import argparse
import fcntl
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
import input_generation_stage_resume as resume


@pytest.fixture
def checkpoint(tmp_path):
    root = tmp_path / "output/input_generation"
    plan = root / "tmp/task_plan.json"
    plan.parent.mkdir(parents=True)
    plan.write_text(json.dumps({"task_count": 1, "tasks": [{"species_prefix": "Species_one"}]}))
    settings = dict(provider="direct", gene_grouping_mode="rescue_overlap", gff_repair_mode="safe",
                    strict="1", run_validate_inputs="1", busco_lineage="eukaryota_odb12")
    Path(str(plan) + ".settings.json").write_text(json.dumps(settings))
    paths = {name: str(tmp_path / name) for name in ("cds_path", "gff_path", "genome_path",
             "cds_output_path", "gff_output_path", "genome_output_path")}
    for path in paths.values():
        Path(path).write_text("validated contents\n")
    meta = root / "tmp/task_meta_shards/1.json"
    meta.parent.mkdir()
    meta.write_text(json.dumps({"task_index": 1, "species_prefix": "Species_one", **paths}))
    for name in ("task_stats_shards/1.json", "task_stats_shards/1.mapping.json",
                 "task_stats_shards/1.longest.json", "species_summary_shards/1.tsv"):
        path = root / "tmp" / name
        path.parent.mkdir(exist_ok=True)
        path.write_text("validated shard\n")
    for stage in ("format", "validate"):
        resume.record(plan, 1, root, stage, "10")
    return plan, root, paths


@pytest.mark.parametrize("stage", ["format", "validate"])
def test_lineage_change_does_not_invalidate_upstream_stages(checkpoint, stage):
    plan, root, _ = checkpoint
    path = Path(str(plan) + ".settings.json")
    settings = json.loads(path.read_text())
    settings["busco_lineage"] = "embryophyta_odb12"
    path.write_text(json.dumps(settings))
    assert resume.valid(plan, 1, root, stage, "10")


@pytest.mark.parametrize("old_contract", ["1", "2"])
def test_old_validation_checkpoint_cannot_bypass_source_gene_ownership_check(checkpoint, old_contract):
    plan, root, _ = checkpoint
    path = resume.checkpoint_path(root, "Species_one", "validate")
    payload = json.loads(path.read_text())
    payload["parameters"]["validation_contract_version"] = old_contract
    path.write_text(json.dumps(payload))
    assert not resume.valid(plan, 1, root, "validate", "10")


@pytest.mark.parametrize("stage,label", [("format", "cds_path"), ("format", "gff_path"),
                                        ("format", "genome_path"), ("format", "cds_output_path"),
                                        ("validate", "cds_output_path"), ("validate", "gff_output_path")])
def test_content_change_cannot_reuse_checkpoint(checkpoint, stage, label):
    plan, root, paths = checkpoint
    Path(paths[label]).write_text("different contents\n")
    assert not resume.valid(plan, 1, root, stage, "10")


@pytest.mark.parametrize("name", ["task_stats_shards/1.json", "species_summary_shards/1.tsv",
                                 "task_stats_shards/1.mapping.json", "task_stats_shards/1.longest.json"])
def test_missing_validation_evidence_cannot_be_skipped(checkpoint, name):
    plan, root, _ = checkpoint
    (root / "tmp" / name).unlink()
    assert not resume.valid(plan, 1, root, "validate", "10")


def test_stricter_settings_invalidate_upstream_proofs(checkpoint):
    plan, root, _ = checkpoint
    path = Path(str(plan) + ".settings.json")
    settings = json.loads(path.read_text())
    settings["require_genome"] = "1"
    path.write_text(json.dumps(settings))
    assert not resume.valid(plan, 1, root, "format", "10")
    assert not resume.valid(plan, 1, root, "validate", "10")


def test_new_formatter_contract_cannot_reuse_old_checkpoints(checkpoint):
    plan, root, _ = checkpoint
    assert not resume.valid(plan, 1, root, "format", "11")
    assert not resume.valid(plan, 1, root, "validate", "11")


def test_disabled_validation_cannot_be_certified(checkpoint):
    plan, root, _ = checkpoint
    path = Path(str(plan) + ".settings.json")
    settings = json.loads(path.read_text())
    settings["run_validate_inputs"] = "0"
    path.write_text(json.dumps(settings))
    with pytest.raises(ValueError, match="Disabled validation"):
        resume.record(plan, 1, root, "validate", "10")


def test_outputs_alone_do_not_prove_validation(checkpoint):
    plan, root, _ = checkpoint
    resume.checkpoint_path(root, "Species_one", "validate").unlink()
    assert not resume.valid(plan, 1, root, "validate", "10")


def test_import_cannot_copy_from_active_donor(checkpoint, tmp_path):
    plan, root, _ = checkpoint
    args = argparse.Namespace(source_plan=plan, source_root=root, task_plan=tmp_path / "new_plan.json",
                              root=tmp_path / "new_root", source_plan_sha256=resume.digest(plan))
    with (root / ".array-phase.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_SH)
        with pytest.raises(ValueError, match="active workers"):
            resume.import_stages(args)


def test_import_rejects_changed_donor_plan(checkpoint, tmp_path):
    plan, root, _ = checkpoint
    args = argparse.Namespace(source_plan=plan, source_root=root, task_plan=tmp_path / "new_plan.json",
                              root=tmp_path / "new_root", source_plan_sha256="0" * 64)
    with pytest.raises(ValueError, match="sealed resume SHA-256"):
        resume.import_stages(args)

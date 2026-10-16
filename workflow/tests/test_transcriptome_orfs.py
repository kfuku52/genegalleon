import copy
import gzip
import json
import subprocess
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))

import format_species_discovery as formatter  # noqa: E402
from transcriptome_orfs import apply_input_type, extract_orf, input_type, verify_input_type  # noqa: E402


def task(tmp_path, text=">TRINITY_DN1_c0_g1_i1\nCCCATGAAATAACCCC\n"):
    source = tmp_path / "source.fa"
    source.write_text(text)
    return dict(provider="local", species_key="Test_species", species_prefix="Test_species",
                cds_path=source, gff_path=None, gbff_path=None, genome_path=None,
                cds_input_type="transcriptome", genetic_code=1, format_strict=True)


def test_native_cdskit_calls_both_strands_and_records_transcript_coordinates():
    value = dict(cds_input_type="transcriptome", genetic_code=1)
    assert extract_orf(value, "forward", "CCCATGAAATAACCCC") == "ATGAAATAA"
    assert extract_orf(value, "reverse", "GGGGTTATTTCATGGG") == "ATGAAATAA"
    assert [r["strand"] for r in value["_transcriptome_orf_audit"]] == ["+", "-"]
    assert all(r["category"] == "complete" for r in value["_transcriptome_orf_audit"])


def test_short_and_invalid_transcripts_do_not_become_fake_cds():
    value = dict(cds_input_type="transcriptome", genetic_code=1)
    assert extract_orf(value, "short", "AC") is None
    assert value["_transcriptome_orf_audit"] == [dict(source_id="short", transcript_length=2, status="no_orf")]
    with pytest.raises(ValueError, match="Invalid transcript"):
        extract_orf(value, "protein", "MPEPTIDE")


def test_formatter_publishes_orf_audit_and_reuses_only_unchanged_input(tmp_path):
    value = task(tmp_path)
    result = formatter.format_cds(value, tmp_path / "out", False, False)
    with gzip.open(result["output_path"], "rt") as stream:
        lines = stream.read().splitlines()
    assert lines[1] == "ATGAAATAA"
    audit = Path(str(result["output_path"]) + ".orf-extraction.json")
    evidence = json.loads(audit.read_text())
    assert evidence["coordinate_system"] == "source transcript; not genomic"
    assert evidence["records"][0]["start"] == 4
    assert evidence["contract"]["selection"] == "complete-first"
    assert formatter.format_cds(copy.deepcopy(value), tmp_path / "out", False, False)["status"] == "skip"
    value["cds_path"].write_text(">TRINITY_DN1_c0_g1_i1\nATGCCCTAA\n")
    with pytest.raises(ValueError, match="stale transcript ORF audit"):
        formatter.format_cds(value, tmp_path / "out", False, False)


def test_existing_ordinary_cds_cannot_be_reused_as_derived_orf(tmp_path):
    value = task(tmp_path)
    ordinary = {**value, "cds_input_type": "cds"}
    formatter.format_cds(ordinary, tmp_path / "out", False, False)
    with pytest.raises(ValueError, match="ORF audit"):
        formatter.format_cds(value, tmp_path / "out", False, False)


@pytest.mark.parametrize("role", ["gff", "gbff", "genome"])
def test_transcript_mode_rejects_genomic_coordinate_roles(role):
    with pytest.raises(ValueError, match="genomic annotation"):
        apply_input_type({}, dict(cds_input_type="transcriptome", cds_url="file:///raw.fa", **{role + "_url": "file:///x"}))


def test_invalid_type_and_staging_type_drift_are_rejected():
    with pytest.raises(ValueError, match="cds_input_type"):
        input_type(dict(cds_input_type="auto"))
    with pytest.raises(ValueError, match="frozen plan"):
        verify_input_type({}, dict(cds_input_type="transcriptome"))
    assert apply_input_type({}, {}) == {}  # Old tasks retain their original shape.


def test_scalar_manifest_uses_the_same_explicit_orf_route(tmp_path):
    value = task(tmp_path)
    manifest = tmp_path / "input.tsv"
    manifest.write_text("provider\tid\tspecies_key\tcds_url\tcds_input_type\n"
                        + "local\tTest_species\tTest_species\t" + value["cds_path"].as_uri() + "\ttranscriptome\n")
    command = [sys.executable, str(SUPPORT / "format_species_inputs.py"), "--provider", "all",
               "--download-manifest", str(manifest), "--download-dir", str(tmp_path / "raw"),
               "--species-cds-dir", str(tmp_path / "cds"), "--species-gff-dir", str(tmp_path / "gff"),
               "--species-genome-dir", str(tmp_path / "genome"), "--strict"]
    subprocess.run(command, check=True, capture_output=True, text=True)
    with gzip.open(next((tmp_path / "cds").glob("*.gz")), "rt") as stream:
        assert stream.read().splitlines()[1] == "ATGAAATAA"


def test_array_prepare_and_worker_preserve_explicit_transcript_input_type(tmp_path):
    value = task(tmp_path)
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("provider\tid\tspecies_key\tbind_local_sources\tcds_url\tcds_input_type\n"
                        + "local\tTest_species\tTest_species\t1\t" + value["cds_path"].as_uri() + "\ttranscriptome\n")
    plan = tmp_path / "plan.json"
    commands = [
        ["plan_input_generation_tasks.py", "--provider", "all", "--download-manifest", str(manifest),
         "--download-dir", str(tmp_path / "downloads"), "--stage-downloads", "--outfile", str(plan), "--strict"],
        ["stage_input_generation_downloads.py", "--task-plan", str(plan)],
        ["run_input_generation_task.py", "--task-plan", str(plan), "--task-index", "1",
         "--species-cds-dir", str(tmp_path / "cds"), "--species-gff-dir", str(tmp_path / "gff"),
         "--species-genome-dir", str(tmp_path / "genome")],
    ]
    for command in commands:
        subprocess.run([sys.executable, str(SUPPORT / command[0]), *command[1:]], check=True, capture_output=True, text=True)
    saved = json.loads(Path(str(plan) + ".tasks/1.json").read_text())
    assert saved["task"]["cds_input_type"] == "transcriptome"
    assert json.loads(plan.read_text())["tasks"][0]["cds_input_type"] == "transcriptome"
    with gzip.open(next((tmp_path / "cds").glob("*.gz")), "rt") as stream:
        assert stream.read().splitlines()[1] == "ATGAAATAA"
    saved["task"].pop("cds_input_type")
    Path(str(plan) + ".tasks/1.json").write_text(json.dumps(saved))
    result = subprocess.run([sys.executable, str(SUPPORT / commands[2][0]), *commands[2][1:]], capture_output=True, text=True)
    assert result.returncode != 0 and "input type differs" in result.stderr

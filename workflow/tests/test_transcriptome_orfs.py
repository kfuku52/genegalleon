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
from transcriptome_orfs import apply_input_type, extract_orf, input_type, signature, verify_input_type  # noqa: E402


def task(tmp_path, text=">TRINITY_DN1_c0_g1_i1\nCCCATGAAATAACCCC\n"):
    source = tmp_path / "source.fa"
    source.write_text(text)
    (tmp_path / "out").mkdir(exist_ok=True)
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
    assert evidence["contract"]["gene_grouping"] == "declared-trinity-gene-v1"
    assert formatter.format_cds(copy.deepcopy(value), tmp_path / "out", False, False)["status"] == "skip"
    value["cds_path"].write_text(">TRINITY_DN1_c0_g1_i1\nATGCCCTAA\n")
    with pytest.raises(ValueError, match="stale transcript ORF audit"):
        formatter.format_cds(value, tmp_path / "out", False, False)


def test_old_transcript_grouping_audit_is_not_reused(tmp_path):
    value = task(tmp_path)
    result = formatter.format_cds(value, tmp_path / "out", False, False)
    audit = Path(str(result["output_path"]) + ".orf-extraction.json")
    evidence = json.loads(audit.read_text())
    evidence["contract"].pop("gene_grouping")
    audit.write_text(json.dumps(evidence))
    with pytest.raises(ValueError, match="stale transcript ORF audit"):
        formatter.format_cds(value, tmp_path / "out", False, False)


@pytest.mark.parametrize('problem', ['none', 'sequence', 'source', 'policy', 'tool', 'grouping'])
def test_independent_historical_rederivation_keeps_ids_and_requires_exact_sources_tools_and_sequences(tmp_path, problem):
    import validate_longest_cds_selection as validator
    from transcriptome_orfs import contract

    value = task(tmp_path, '>TRINITY_DN1_c0_g1_i1\nATGAAATAA\n>TRINITY_DN1_c0_g1_i2\nATGCCCTAA\n')
    output = tmp_path / 'out/Test_species.fa.gz'
    with gzip.open(output, 'wt') as stream:
        stream.write('>Test_species_TRINITY_DN1_c0_g1_i1\nATGAAATAA\n>Test_species_TRINITY_DN1_c0_g1_i2\nATGCCCTAA\n')
    original = contract(value)
    original.pop('gene_grouping')
    original['implementation'] = 'a' * 64
    if problem == 'sequence':
        with gzip.open(output, 'wt') as stream:
            stream.write('>Test_species_TRINITY_DN1_c0_g1_i1\nATGTTTTAA\n>Test_species_TRINITY_DN1_c0_g1_i2\nATGCCCTAA\n')
    if problem == 'source':
        value['cds_path'].write_text('>TRINITY_DN1_c0_g1_i1\nATGTTTTAA\n')
    if problem == 'policy':
        original['selection'] = 'longest-any'
    if problem == 'tool':
        original['cdskit_version'] = 'different'
    if problem == 'grouping':
        original['gene_grouping'] = 'unknown'
    audit = Path(str(output) + '.orf-extraction.json')
    audit.write_text(json.dumps(dict(contract=original, output=signature(output))))
    before = audit.read_bytes()
    row = dict(provider='local', species_key='Test_species', species_prefix='Test_species',
               cds_input_path=str(value['cds_path']), cds_input_type='transcriptome', cds_output_path=str(output),
               cds_sequences_before='2', cds_sequences_after='2', aggregated_cds_removed='0')
    args = dict(cds_file=output, summary_row=row, index=1)
    assert not validator.validate_single_species(args, 10)['ok']
    result = validator.validate_single_species({**args, 'revalidate_historical_orfs': True}, 10)
    assert result['ok'] == (problem == 'none'), result
    assert audit.read_bytes() == before


def test_existing_ordinary_cds_cannot_be_reused_as_derived_orf(tmp_path):
    value = task(tmp_path)
    ordinary = {**value, "cds_input_type": "cds"}
    formatter.format_cds(ordinary, tmp_path / "out", False, False)
    with pytest.raises(ValueError, match="ORF audit"):
        formatter.format_cds(value, tmp_path / "out", False, False)


def test_declared_trinity_transcripts_select_one_longest_cds_per_gene(tmp_path):
    value = task(tmp_path,
                 ">TRINITY_DN1_c0_g1_i1\nATGAAATAA\n"
                 ">TRINITY_DN1_c0_g1_i2\nATGAAAAAATAA\n"
                 ">TRINITY_DN1_c0_g2_i1\nATGCCCTAA\n")
    result = formatter.format_cds(value, tmp_path / "out", False, False)
    assert result["before_count"] == 3 and result["after_count"] == 2
    with gzip.open(result["output_path"], "rt") as stream:
        assert stream.read() == (
            ">Test_species_TRINITY_DN1_c0_g1\nATGAAAAAATAA\n"
            ">Test_species_TRINITY_DN1_c0_g2\nATGCCCTAA\n")


@pytest.mark.parametrize("source_type,identifier", [
    ("cds", "TRINITY_DN1_c0_g1_i1"),
    ("transcriptome", "other_i1"),
    ("transcriptome", "TRINITY_DN1_c0_g1_i1_extra"),
])
def test_transcript_gene_grouping_leaves_other_identifiers_unchanged(tmp_path, source_type, identifier):
    value = task(tmp_path, f">{identifier}\nATGAAATAA\n")
    value["cds_input_type"] = source_type
    result = formatter.format_cds(value, tmp_path / "out", False, False)
    with gzip.open(result["output_path"], "rt") as stream:
        assert stream.readline().strip() == ">Test_species_" + identifier


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


def test_transcript_summary_validates_rederived_orfs_and_rejects_tampering(tmp_path):
    import csv

    import format_species_summary as summary
    import validate_longest_cds_selection as validator

    value = task(tmp_path, ">forward\nCCCATGAAATAACCCC\n>reverse\nGGGGTTATTTCATGGG\n>short\nAC\n")
    result = formatter.format_cds(value, tmp_path / "out", False, False)
    missing = dict(status="missing", output_path=None)
    row = summary.build_species_summary_row(value, result, missing, missing, "now", False, False)
    output = tmp_path / "summary.tsv"
    summary.write_species_summary_rows(output, {"fixture": row})
    with output.open() as handle:
        saved = next(csv.DictReader(handle, delimiter="\t"))
    assert saved["cds_input_type"] == "transcriptome"
    args = dict(cds_file=result["output_path"], summary_row=saved, index=0)
    assert validator.validate_single_species(args, 10)["ok"]
    changed = dict(saved, cds_input_type="cds")
    assert not validator.validate_single_species({**args, "summary_row": changed}, 10)["ok"]
    audit = Path(str(result["output_path"]) + ".orf-extraction.json")
    original = audit.read_bytes()
    audit.unlink()
    assert not validator.validate_single_species(args, 10)["ok"]
    audit.write_bytes(original)
    output_bytes = result["output_path"].read_bytes()
    with gzip.open(result["output_path"], "wt") as stream:
        stream.write(">Test_species_forward\nATGCCCTAA\n>Test_species_reverse\nATGCCCTAA\n")
    forged = json.loads(original)
    forged["output"] = signature(result["output_path"])
    audit.write_text(json.dumps(forged))
    assert not validator.validate_single_species(args, 10)["ok"]
    result["output_path"].write_bytes(output_bytes)
    audit.write_bytes(original)
    value["cds_path"].write_text(">forward\nATGCCCTAA\n")
    assert not validator.validate_single_species(args, 10)["ok"]


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
         "--species-genome-dir", str(tmp_path / "genome"), "--species-summary-output", str(tmp_path / "summary.tsv")],
        ["validate_longest_cds_selection.py", "--species-cds-dir", str(tmp_path / "cds"),
         "--species-summary", str(tmp_path / "summary.tsv"), "--summary-outputs-only"],
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

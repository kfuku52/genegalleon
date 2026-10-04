import gzip
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
import format_species_inputs as formatter
from curate_paired_species_inputs import curate_pair, fingerprint, validated_formatted_pair
from format_species_summary import build_species_summary_row, write_species_summary_rows
from stage_input_generation_downloads import bound_local_manifest_task, explicit_manifest_task
from validate_longest_cds_selection import build_task_from_summary_row, validate_single_species


def curated_task(tmp_path):
    cds = tmp_path / "cds.fa.gz"
    gff = tmp_path / "annotations.gff.gz"
    genome = tmp_path / "genome.fa.gz"
    with gzip.open(cds, "wt") as handle:
        handle.write(">Test_species_g1 literal description\nATGAAACCCTAA\n>Test_species_g2\nTAATGAAACCCTAAC\n")
    with gzip.open(genome, "wt") as handle:
        handle.write(">chr1\nATGAAACCCTAA\n>chr2\nTAATGAAACCCTAAC\n")
    with gzip.open(gff, "wt") as handle:
        handle.write("##gff-version 3\n"
                     "chr1\ts\tgene\t1\t12\t.\t+\t.\tID=g1\n"
                     "chr1\ts\tmRNA\t1\t12\t.\t+\t.\tID=t1;Parent=g1\n"
                     "chr1\ts\tCDS\t1\t12\t.\t+\t0\tParent=t1\n"
                     "chr2\ts\tgene\t1\t15\t.\t+\t.\tID=g2\n"
                     "chr2\ts\tmRNA\t1\t15\t.\t+\t.\tID=t2;Parent=g2\n"
                     "chr2\ts\tCDS\t1\t4\t.\t+\t1\tParent=t2\n"
                     "chr2\ts\tCDS\t5\t15\t.\t+\t2\tParent=t2\n")
    policy = tmp_path / "decisions.json"
    policy.write_text(json.dumps(dict(schema_version=1, species="Test_species", decision_basis="Keep already adopted records",
                                     input_sha256={key: fingerprint(path)["sha256"] for key, path in dict(cds=cds, gff=gff, genome=genome).items()}, records=[])))
    report = curate_pair("Test_species", cds, gff, genome, policy, tmp_path / "curated")
    return dict(provider="direct", species_key="Test_species", species_prefix="Test_species",
                cds_path=Path(report["cds_output"]["path"]), gff_path=Path(report["gff_output"]["path"]), genome_path=genome,
                gbff_path=None, paired_curation=json.dumps(report), gene_grouping_mode="strict", genetic_code=1)


def test_explicit_curated_import_preserves_partial_codon_CDS_GFF_and_genome(tmp_path):
    task = curated_task(tmp_path)
    ordinary = {key: value for key, value in task.items() if key != "paired_curation"}
    normal = tmp_path / "ordinary"
    normal.mkdir()
    result = formatter.format_cds(ordinary, normal, False, False)
    assert result["cds_normalised_records"] == 1
    assert dict(formatter.iter_fasta_records(result["output_path"]))["Test_species_g2"] == "ATGAAACCCTAA"
    output = tmp_path / "preserved"
    cds = formatter.format_cds(task, output, False, False)
    gff = formatter.format_gff(task, output, False, False, formatted_cds_path=cds["output_path"])
    genome = formatter.format_genome(task, output, False, False)
    for role, result in dict(cds=cds, gff=gff, genome=genome).items():
        assert result["output_path"].read_bytes() == task[role + "_path"].read_bytes()
    assert cds["before_count"] == cds["after_count"] == 2 and cds["duplicates"] == 0
    assert formatter.format_cds(task, output, False, False)["status"] == "skip"
    row = build_species_summary_row(task, cds, gff, genome, "test", False, False)
    summary = tmp_path / "species.tsv"
    write_species_summary_rows(summary, {"test": row})
    import csv
    with summary.open() as handle:
        persisted = next(csv.DictReader(handle, delimiter="\t"))
    assert json.loads(build_task_from_summary_row(persisted)["paired_curation"]) == json.loads(task["paired_curation"])
    checked = validate_single_species(dict(index=1, cds_file=cds["output_path"], summary_row=persisted), 10)
    assert checked["ok"], checked
    assert checked["stats"]["aggregated_cds_removed"] == 0


@pytest.mark.parametrize("role", ["cds", "gff", "genome"])
def test_changed_curated_source_is_rejected_before_output(tmp_path, role):
    task = curated_task(tmp_path)
    with task[role + "_path"].open("ab") as handle:
        handle.write(b"changed")
    with pytest.raises(ValueError, match="input hashes differ"):
        formatter.format_cds(task, tmp_path / "output", False, False)
    assert not (tmp_path / "output").exists()


def test_cached_receipt_does_not_accept_a_replaced_source_path(tmp_path):
    task = curated_task(tmp_path)
    validated_formatted_pair(task)
    other = tmp_path / "other.fa.gz"
    other.write_bytes(b"unapproved")
    task["cds_path"] = other
    with pytest.raises(ValueError, match="input hashes differ"):
        formatter.format_cds(task, tmp_path / "output", False, False)


def test_curated_import_never_overwrites_its_input_archive(tmp_path):
    task = curated_task(tmp_path)
    before = {role: task[role+"_path"].read_bytes() for role in ("cds", "gff", "genome")}
    with pytest.raises(ValueError, match="separate from its immutable source"):
        formatter.format_cds(task, task["cds_path"].parent, True, False)
    assert {role: task[role+"_path"].read_bytes() for role in before} == before


def test_declared_receipt_does_not_waive_coordinate_checks(tmp_path):
    task = curated_task(tmp_path)
    with gzip.open(task["gff_path"], "rt") as handle:
        text = handle.read()
    with gzip.open(task["gff_path"], "wt") as handle:
        handle.write(text.replace("CDS\t5\t15", "CDS\t5\t99"))
    receipt = json.loads(task["paired_curation"])
    receipt["gff_output"]["sha256"] = fingerprint(task["gff_path"])["sha256"]
    task["paired_curation"] = json.dumps(receipt)
    with pytest.raises(ValueError, match="coordinates exceed"):
        formatter.format_cds(task, tmp_path / "output", False, False)


def test_declared_receipt_does_not_waive_duplicate_gene_identity(tmp_path):
    task = curated_task(tmp_path)
    with gzip.open(task["cds_path"], "at") as handle:
        handle.write(">Test_species_g1\nATGAAACCCTAA\n")
    receipt = json.loads(task["paired_curation"])
    receipt["cds_output"]["sha256"] = fingerprint(task["cds_path"])["sha256"]
    task["paired_curation"] = json.dumps(receipt)
    with pytest.raises(ValueError, match="unique"):
        formatter.format_cds(task, tmp_path / "output", False, False)


def test_native_manifest_routes_retain_the_frozen_curation_receipt(tmp_path):
    from format_species_provider_resolvers import provider_raw_dir
    task = curated_task(tmp_path)
    row = dict(paired_curation=task["paired_curation"], bind_local_sources="1")
    for role in ("cds", "gff", "genome"):
        row[role + "_url"] = task[role + "_path"].as_uri()
        row[role + "_filename"] = task[role + "_path"].name
    planned = dict(task, manifest_row=row, input_sha256={str(task[role + "_path"]): fingerprint(task[role + "_path"])["sha256"] for role in ("cds", "gff", "genome")})
    assert bound_local_manifest_task(planned)["paired_curation"] == task["paired_curation"]
    raw = provider_raw_dir("direct", tmp_path / "download", "Test_species")
    raw.mkdir(parents=True)
    for role in ("cds", "gff", "genome"):
        shutil.copyfile(task[role + "_path"], raw / task[role + "_path"].name)
    assert explicit_manifest_task(planned, row, tmp_path / "download")["paired_curation"] == task["paired_curation"]


def test_curator_cli_produces_an_exclusive_native_manifest(tmp_path):
    import csv
    task = curated_task(tmp_path)
    decision = tmp_path / "second-decisions.json"
    decision.write_text(json.dumps(dict(schema_version=1, species="Test_species", decision_basis="Preserve approved formatted input",
                                       input_sha256={role: fingerprint(task[role + "_path"])["sha256"] for role in ("cds", "gff", "genome")}, records=[])))
    manifest = tmp_path / "native.tsv"
    command = [sys.executable, str(Path(__file__).resolve().parents[1]/"support/curate_paired_species_inputs.py"), "curate",
               "--species", "Test_species", "--cds", str(task["cds_path"]), "--gff", str(task["gff_path"]), "--genome", str(task["genome_path"]),
               "--decision-manifest", str(decision), "--output-dir", str(tmp_path/"second-curated"), "--download-manifest", str(manifest)]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    with manifest.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    receipt = json.loads(row["paired_curation"])
    assert receipt["remaining_cds_records"] == 2 and row["bind_local_sources"] == "1"
    assert row["genome_sha256"] == fingerprint(task["genome_path"])["sha256"]
    before = manifest.read_bytes()
    again = subprocess.run(command, capture_output=True, text=True)
    assert again.returncode != 0 and "manifest already exists" in again.stderr
    assert manifest.read_bytes() == before


def test_approved_pair_import_rejects_invalid_attributes_without_editing_source(tmp_path):
    task = curated_task(tmp_path)
    source = task["gff_path"]
    with gzip.open(source, "rt") as handle:
        text = handle.read()
    with gzip.open(source, "wt") as handle:
        handle.write(text.replace("ID=g1\n", "ID=g1;Name=SULTR4;1;\n"))
    receipt = json.loads(task["paired_curation"])
    receipt["gff_output"]["sha256"] = fingerprint(source)["sha256"]
    task["paired_curation"] = json.dumps(receipt)
    before = source.read_bytes()
    with pytest.raises(ValueError, match="regenerate with gg_input_generation"):
        formatter.format_gff(task, tmp_path / "output", False, False)
    assert source.read_bytes() == before and not (tmp_path / "output").exists()

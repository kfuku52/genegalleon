"""Coding bounds require complete exact publisher evidence, including UTRs."""
import gzip
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from format_species_annotation.coding_coordinates import unique_exact_offset
from format_species_annotation.common import parse_gff_attributes, reverse_complement
from format_species_annotation.source_identity import coding_coordinate_audit, task_annotation_path
from format_species_discovery import format_cds, format_gff
from validate_cds_gff_mapping import validate_single_species


def bundle(tmp_path, strand="+"):
    # Two UTR bases, then coding sequence split 4+5 across an intron, then UTR.
    genomic = "CCATGA" + "GGGGGG" + "AATAATT"
    if strand == "-":
        genomic = reverse_complement(genomic)
        intervals = [(1, 7), (14, 19)]
    else:
        intervals = [(1, 6), (13, 19)]
    genome, gff, cds = (tmp_path / name for name in ("genome.fa", "source.gff", "cds.fa"))
    genome.write_text(">chr1\n" + genomic + "\n>unplaced\nACGT\n")
    gff.write_text("##gff-version 3\n" + f"chr1\tauthor\tgene\t1\t19\t.\t{strand}\t.\tID=g\n" + "".join(
        f"chr1\tauthor\texon\t{a}\t{b}\t.\t{strand}\t.\tID=e{i};Parent=g\n"
        for i, (a, b) in enumerate(intervals)))
    cds.write_text(">g\nATGAAATAA\n")
    return dict(provider="direct", species_key="Test_species", species_prefix="Test_species",
                gff_path=gff, genome_path=genome, cds_path=cds, format_strict=True,
                gff_repair_mode="safe", genetic_code=1)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_native_formatter_reconstructs_exact_cds_bounds_and_preserves_exons(tmp_path, strand):
    task = bundle(tmp_path, strand)
    original = {key: task[key].read_bytes() for key in ("gff_path", "genome_path", "cds_path")}
    view = task_annotation_path(task)
    rows = [line.split("\t") for line in view.read_text().splitlines() if not line.startswith("#")]
    coding = [row for row in rows if row[2] == "CDS"]
    assert [(int(row[3]), int(row[4]), int(row[7])) for row in coding] == (
        [(3, 6, 0), (13, 17, 2)] if strand == "+" else [(14, 17, 0), (3, 7, 2)])
    assert all(parse_gff_attributes(row[8])["Parent"] == ("g.gg_published_cds",)
               for row in rows if row[2] in {"exon", "CDS"})
    out = tmp_path / "out"
    out.mkdir()
    cds = format_cds(task, out, False, False)
    gff = format_gff(task, out, False, False, formatted_cds_path=cds["output_path"])
    with gzip.open(cds["output_path"], "rt") as handle:
        assert handle.read() == ">Test_species_g\nATGAAATAA\n"
    audit = json.loads(Path(str(gff["output_path"]) + ".repair.json").read_text())
    assert audit["coding_coordinates"]["status"] == "reconstructed"
    assert audit["coding_coordinates"]["gene_models"] == 1
    assert len(audit["coding_coordinates"]["source_sha256"]) == 3
    result = validate_single_species(dict(index=1, species_prefix="Test_species",
                                         cds_file=Path(cds["output_path"]), gff_file=Path(gff["output_path"]),
                                         genome_file=task["genome_path"], strict=True), 10)
    assert result["ok"], result
    assert {key: task[key].read_bytes() for key in original} == original
    assert task_annotation_path(task) == view


@pytest.mark.parametrize("fault", ["ambiguous", "mismatch", "extra_cds", "duplicate_cds",
                                  "missing_exon", "overlap", "strand", "bounds", "missing_reference",
                                  "collision", "multiple_parents", "ambiguous_bases"])
def test_unproved_coding_coordinates_are_never_exposed(tmp_path, fault):
    task = bundle(tmp_path)
    gff, cds, genome = (task[key] for key in ("gff_path", "cds_path", "genome_path"))
    if fault == "ambiguous":
        genome.write_text(">chr1\nATGAAATAANNNNATGAAATAA\n")
        gff.write_text("chr1\t.\tgene\t1\t22\t.\t+\t.\tID=g\nchr1\t.\texon\t1\t22\t.\t+\t.\tParent=g\n")
    elif fault == "mismatch":
        cds.write_text(">g\nATGCCCTAA\n")
    elif fault in {"extra_cds", "duplicate_cds"}:
        with cds.open("a") as handle:
            handle.write(">" + ("extra" if fault == "extra_cds" else "g") + "\nATGAAATAA\n")
    elif fault == "missing_exon":
        gff.write_text(gff.read_text().replace("Parent=g", "Parent=other"))
    elif fault == "overlap":
        gff.write_text(gff.read_text().replace("exon\t13\t19", "exon\t5\t19"))
    elif fault == "strand":
        gff.write_text(gff.read_text().replace("exon\t1\t6\t.\t+", "exon\t1\t6\t.\t-"))
    elif fault == "bounds":
        genome.write_text(">chr1\nACGT\n")
    elif fault == "missing_reference":
        genome.write_text(genome.read_text().replace(">chr1", ">other"))
    elif fault == "collision":
        gff.write_text(gff.read_text().replace("ID=e0", "ID=g.gg_published_cds"))
    elif fault == "multiple_parents":
        gff.write_text(gff.read_text().replace("Parent=g", "Parent=g,other"))
    elif fault == "ambiguous_bases":
        cds.write_text(">g\nATGNAATAA\n")
    with pytest.raises((ValueError, KeyError)):
        task_annotation_path(task)


def test_repair_off_and_existing_coding_models_keep_original_source(tmp_path):
    task = bundle(tmp_path)
    assert task_annotation_path({**task, "gff_repair_mode": "off"}) == task["gff_path"]
    assert task_annotation_path({**task, "gff_repair_mode": " off "}) == task["gff_path"]
    with task["gff_path"].open("a") as handle:
        handle.write("chr1\tauthor\tCDS\t3\t6\t.\t+\t0\tParent=g\n")
    assert task_annotation_path(task) == task["gff_path"]
    assert coding_coordinate_audit(task)["status"] == "not_applicable"


def test_changed_input_invalidates_cached_view(tmp_path):
    task = bundle(tmp_path)
    task_annotation_path(task)
    task["cds_path"].write_text(">g\nATGCCCTAA\n")
    with pytest.raises(ValueError, match="one exact"):
        task_annotation_path(task)


@pytest.mark.parametrize("spliced,coding,expected", [
    ("CCATGNAATAATT", "ATGNAATAA", 2),
    ("ATGNATGA", "ATGN", None),
    ("ATGAATGN", "ATGA", None),
    ("ATGAAATAA", "ATGNAATAA", None),
    ("ACGT", "NNNN", None),
    ("NATGA", "ATGA", 1),
])
def test_unknown_bases_cannot_establish_false_unique_boundaries(spliced, coding, expected):
    assert unique_exact_offset(spliced, coding) == expected

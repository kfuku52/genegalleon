"""Source ID collisions must never merge paralogs or split ordinary CDS parts."""
import gzip
import json
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))
from format_species_annotation.common import parse_gff_attributes, reverse_complement  # noqa: E402
from format_species_annotation.locus_identity import LocusIdentityError, suffix  # noqa: E402
from format_species_annotation.source_identity import (  # noqa: E402
    locus_identity_audit,
    source_annotation_path,
    task_annotation_path,
)
from format_species_discovery import format_cds, format_gff  # noqa: E402
from validate_cds_gff_mapping import validate_single_species  # noqa: E402
from validate_longest_cds_selection import collect_expected_longest_records  # noqa: E402


def bundle(tmp_path, layout="contigs", provided=False):
    locations = [("chr1", "+", 1), ("chr2", "+", 1)]
    if layout == "same_chromosome":
        locations[1] = ("chr1", "+", 31)
    if layout == "opposite_strands":
        locations[1] = ("chr1", "-", 31)
    sequences = {seqid: list("N" * 60) for seqid, _strand, _start in locations}
    rows, fasta = [], []
    for copy, (seqid, strand, start) in enumerate(locations):
        sequence = "ATGAAATAACCCGGGTTT"
        sequences[seqid][start - 1:start + 17] = sequence if strand == "+" else reverse_complement(sequence)
        def line(kind, left, right, attrs, phase=".", seqid=seqid, strand=strand):
            return f"{seqid}\tsource\t{kind}\t{left}\t{right}\t.\t{strand}\t{phase}\t{attrs}\n"
        rows.append(line("gene", start, start + 17, "ID=g;gene_id=g;Note=original"))
        for isoform, length in ((1, 9), (2, 6)):
            left, right = (start, start + length - 1) if strand == "+" else (start + 18 - length, start + 17)
            rows.append(line("mRNA", left, right, f"ID=r{isoform};Parent=g;gene_id=g"))
            parts = [(left, right)] if isoform == 2 else [(left, left + 2), (left + 3, right)]
            for a, b in parts:
                rows.append(line("CDS", a, b, f"ID=c{isoform};Parent=r{isoform};protein_id=P{copy}_{isoform};gene_id=g", "0"))
            fasta.append(f">P{copy}_{isoform}\n{sequence[:length]}\n")
    genome, gff, cds = (tmp_path / name for name in ("genome.fa", "models.gff", "cds.fa"))
    genome.write_text("".join(">" + seqid + "\n" + "".join(sequence) + "\n" for seqid, sequence in sequences.items()))
    gff.write_text("##gff-version 3\n" + "".join(rows))
    task = dict(provider="direct", species_key="Test_species", species_prefix="Test_species",
                gff_path=gff, genome_path=genome, format_strict=True, gff_repair_mode="safe")
    if provided:
        cds.write_text("".join(fasta))
        task["cds_path"] = cds
    return task


@pytest.mark.parametrize("layout", ["contigs", "same_chromosome", "opposite_strands"])
@pytest.mark.parametrize("provided", [False, True])
def test_collision_keeps_both_loci_and_selects_one_isoform_each(tmp_path, layout, provided):
    task = bundle(tmp_path, layout, provided)
    source = task["gff_path"].read_bytes()
    cds = format_cds(task, tmp_path, False, False)
    assert (cds["before_count"], cds["after_count"]) == (4, 2)
    with gzip.open(cds["output_path"], "rt") as handle:
        records = handle.read().splitlines()
    assert records[1::2] == ["ATGAAATAA"] * 2
    gene_ids = {line[1:].removeprefix("Test_species_") for line in records[::2]}
    assert len(gene_ids) == 2 and all(identifier.startswith("g.locus") for identifier in gene_ids)
    gff = format_gff(task, tmp_path, False, False, formatted_cds_path=cds["output_path"])
    with gzip.open(gff["output_path"], "rt") as handle:
        rows = [line.split("\t") for line in handle if not line.startswith("#")]
    ids = {identifier for row in rows for identifier in parse_gff_attributes(row[8]).get("ID", ())}
    assert {parse_gff_attributes(row[8])["ID"][0] for row in rows if row[2] == "gene"} == gene_ids
    assert all(parent in ids for row in rows for parent in parse_gff_attributes(row[8]).get("Parent", ()))
    assert all(row[7] == "0" for row in rows if row[2] == "CDS")
    audit = json.loads(Path(str(gff["output_path"]) + ".repair.json").read_text())
    assert audit["locus_identity"]["status"] == "repaired"
    assert len(audit["locus_identity"]["mappings"]) == 10
    assert task["gff_path"].read_bytes() == source
    expected, proof = collect_expected_longest_records(task)
    assert set(expected) == {"Test_species_" + gene for gene in gene_ids}
    assert proof["source_gene_complete"]
    mapping = validate_single_species(dict(index=1, species_prefix="Test_species",
        cds_file=cds["output_path"], gff_file=gff["output_path"], genome_file=task["genome_path"], strict=True), 10)
    assert mapping["ok"], mapping.get("error")


def test_ids_are_stable_under_row_reordering_and_repair_is_idempotent(tmp_path):
    task = bundle(tmp_path)
    initial = task_annotation_path(task)
    original = task["gff_path"].read_text().splitlines(True)
    task["gff_path"].write_text(original[0] + "".join(reversed(original[1:])))
    shuffled = task_annotation_path(task)
    def ids(path):
        return {parse_gff_attributes(line.split("\t")[8])["ID"][0]
                for line in path.read_text().splitlines() if not line.startswith("#")}
    assert ids(initial) == ids(shuffled)
    assert source_annotation_path(shuffled, repair_locus_ids=True) == shuffled


def test_off_preserves_source_identifiers(tmp_path):
    task = bundle(tmp_path)
    task["gff_repair_mode"] = "off"
    assert task_annotation_path(task) == task["gff_path"]
    assert locus_identity_audit(task)["status"] == "off"


def test_ambiguous_supplied_cds_id_is_not_assigned_to_a_copy(tmp_path):
    task = bundle(tmp_path, provided=True)
    task["cds_path"].write_text(">r1\nATGAAATAA\n")
    with pytest.raises(ValueError, match="cannot identify one repaired GFF locus"):
        format_cds(task, tmp_path, False, False)
    assert not (tmp_path / "Test_species_cds.fasta.gz").exists()


def test_exact_source_locations_disambiguate_supplied_cds_copies(tmp_path):
    task = bundle(tmp_path, provided=True)
    task["cds_path"].write_text("".join(
        f">lcl|chr{copy}_cds_r1_{copy} [protein_id=r1] [location=join(1..3,4..9)]\nATGAAATAA\n"
        for copy in (1, 2)))
    result = format_cds(task, tmp_path, False, False)
    assert result["after_count"] == 2
    expected, proof = collect_expected_longest_records(task)
    assert len(expected) == 2 and proof["source_gene_complete"]


@pytest.mark.parametrize("declaration", ["", ";number=1", ";exception=trans-splicing;part=1/2"])
def test_multipart_cds_is_preserved(tmp_path, declaration):
    source = tmp_path / "source.gff"
    source.write_text("chr1\ts\tgene\t1\t30\t.\t+\t.\tID=g\n"
                      "chr1\ts\tmRNA\t1\t30\t.\t+\t.\tID=r;Parent=g\n"
                      f"chr1\ts\tCDS\t1\t9\t.\t+\t0\tID=c;Parent=r{declaration}\n"
                      f"chr1\ts\tCDS\t21\t29\t.\t+\t0\tID=c;Parent=r{declaration.replace('1', '2')}\n")
    assert source_annotation_path(source, repair_locus_ids=True) == source


@pytest.mark.parametrize("declaration", ["number=1", "exception=trans-splicing;part=1/2"])
def test_ordered_cross_contig_cds_identity_is_preserved(tmp_path, declaration):
    source = tmp_path / "source.gff"
    source.write_text("chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=r\n"
                      f"chr1\ts\tCDS\t1\t9\t.\t+\t0\tID=c;Parent=r;{declaration}\n"
                      f"chr2\ts\tCDS\t1\t9\t.\t+\t0\tID=c;Parent=r;{declaration.replace('1', '2')}\n")
    assert source_annotation_path(source, repair_locus_ids=True) == source


@pytest.mark.parametrize("defect", ["missing_rna", "outside_rna", "overlapping_definitions", "generated_id"])
def test_unproven_collision_stops_with_source_evidence(tmp_path, defect):
    task = bundle(tmp_path)
    text = task["gff_path"].read_text()
    if defect == "missing_rna":
        text = "\n".join(line for line in text.splitlines() if not (line.startswith("chr2") and "\tmRNA\t" in line)) + "\n"
    elif defect == "outside_rna":
        text = text.replace("chr2\tsource\tmRNA\t1\t9", "chr2\tsource\tmRNA\t1\t2")
    elif defect == "overlapping_definitions":
        text += "chr1\tsource\tgene\t1\t20\t.\t+\t.\tID=g\n"
    else:
        identifier = suffix("g", ("gene", ("chr1", "+", 1, 18), ()))
        text += f"chr3\tsource\tgene\t1\t9\t.\t+\t.\tID={identifier}\n"
    task["gff_path"].write_text(text)
    with pytest.raises(LocusIdentityError) as error:
        format_cds(task, tmp_path, False, False)
    audit = error.value.audit
    assert audit["status"] == "blocked" and audit["problems"]
    assert any(path.name.endswith(".locus-identity.json") for path in tmp_path.iterdir())
    assert task["gff_path"].read_text() == text


def test_pinus_like_cds_parent_reuse_is_not_reinterpreted_as_two_transcripts(tmp_path):
    source = tmp_path / "source.gff"
    source.write_text("chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=r\n"
                      "chr1\ts\tCDS\t1\t9\t.\t+\t0\tID=c;Parent=r\n"
                      "chr2\ts\tCDS\t1\t9\t.\t+\t0\tID=c;Parent=r\n")
    with pytest.raises(LocusIdentityError, match="parent model"):
        source_annotation_path(source, repair_locus_ids=True)


def test_shared_physical_cds_with_multiple_parents_keeps_its_id(tmp_path):
    source = tmp_path / "source.gff"
    source.write_text("chr1\ts\tgene\t1\t9\t.\t+\t.\tID=g1\n"
                      "chr1\ts\tgene\t1\t9\t.\t+\t.\tID=g2\n"
                      "chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=r1;Parent=g1\n"
                      "chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=r2;Parent=g2\n"
                      "chr1\ts\tCDS\t1\t9\t.\t+\t0\tID=c;Parent=r1\n"
                      "chr1\ts\tCDS\t1\t9\t.\t+\t0\tID=c;Parent=r2\n")
    assert source_annotation_path(source, repair_locus_ids=True) == source


def test_collision_does_not_waive_cyclic_parent_ancestry(tmp_path):
    task = bundle(tmp_path)
    task["gff_path"].write_text(task["gff_path"].read_text().replace("ID=g;gene_id=g", "ID=g;Parent=g;gene_id=g"))
    with pytest.raises(LocusIdentityError) as error:
        task_annotation_path(task)
    assert any(item["reason"] == "cyclic parent ancestry" for item in error.value.audit["problems"])


def test_unassignable_identity_attribute_is_rejected_before_any_view_is_written(tmp_path):
    task = bundle(tmp_path)
    task["gff_path"].write_text(task["gff_path"].read_text() +
        "chr3\ts\tregion\t1\t9\t.\t+\t.\tID=outside;gene_id=g\n")
    with pytest.raises(LocusIdentityError) as error:
        task_annotation_path(task)
    assert any(item["reason"] == "identity reference cannot be assigned to one locus"
               for item in error.value.audit["problems"])


def test_source_id_is_part_of_the_suffix_namespace():
    scope = ("gene", ("chr1", "+", 1, 18), ())
    assert suffix("g-a", scope).split(".locus")[1] != suffix("g_a", scope).split(".locus")[1]

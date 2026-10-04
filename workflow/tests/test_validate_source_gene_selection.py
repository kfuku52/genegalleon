import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
import validate_longest_cds_selection as longest
from validate_source_gene_selection import SourceGeneSelection


@pytest.fixture
def annotation(tmp_path):
    path = tmp_path / "source.gff"
    # Children intentionally precede their parents; ancestry must be order-free.
    path.write_text(
        "chr1\ts\tmRNA\t1\t20\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ts\tmRNA\t1\t20\t.\t+\t.\tID=t2;Parent=g1\n"
        "chr1\ts\tgene\t1\t20\t.\t+\t.\tID=g1\n"
        "chr1\ts\tgene\t30\t50\t.\t+\t.\tID=g2\n"
        "chr1\ts\tmRNA\t30\t50\t.\t+\t.\tID=t3;Parent=g2\n")
    return path


def test_independent_source_check_rejects_retained_isoforms(annotation):
    check = SourceGeneSelection(annotation)
    check.observe("t1", "Species_one_gene1", True)
    check.observe("t2", "Species_one_transcript2", True)
    with pytest.raises(ValueError, match="retained_isoform_gene_groups=1"):
        check.validate()


def test_independent_source_check_rejects_distinct_gene_collapse(annotation):
    check = SourceGeneSelection(annotation)
    check.observe("t1", "Species_one_collapsed", True)
    check.observe("t3", "Species_one_collapsed", False)
    with pytest.raises(ValueError, match="distinct_source_gene_merges=1"):
        check.validate()


def test_independent_source_check_only_counts_retained_representatives(annotation):
    check = SourceGeneSelection(annotation)
    check.observe("t1", "one", True)
    check.observe("t2", "unresolved", False)
    check.observe("no_source_identity", "unresolved", True)
    stats = check.validate()
    assert stats["source_gene_resolved_records"] == 2
    assert stats["source_gene_unresolved_records"] == 1
    assert stats["source_gene_check"] == "partial_explicit_parent"
    assert not stats["source_gene_complete"]
    assert stats["source_gene_unresolved_selected_records"] == 1


def test_independent_source_check_accepts_correct_one_gene_selection(annotation):
    check = SourceGeneSelection(annotation)
    check.observe("t1", "one", True)
    check.observe("t2", "one", True)
    check.observe("t3", "two", True)
    assert check.validate()["source_gene_retained_isoform_groups"] == 0


def test_native_longest_validation_does_not_trust_its_own_wrong_grouping(annotation, monkeypatch):
    task = dict(provider="direct", species_key="Species_one", species_prefix="Species_one",
                gff_path=annotation, cds_path=annotation)
    monkeypatch.setattr(longest.formatter, "prepare_cds_identifier_task", lambda task: task)
    monkeypatch.setattr(longest.formatter, "iter_normalised_cds_records", lambda task: iter([("t1", "ATGAAA", {}), ("t2", "ATGAAATTT", {})]))
    monkeypatch.setattr(longest.formatter, "build_formatted_cds_id", lambda task, header: header)
    monkeypatch.setattr(longest.formatter, "build_gene_aggregate_id", lambda task, header, tid: "Species_one_" + tid)
    with pytest.raises(ValueError, match="retained_isoform_gene_groups=1"):
        longest.collect_expected_longest_records(task)


def test_source_check_does_not_guess_missing_parent_from_suffix(tmp_path):
    path = tmp_path / "source.gff"
    path.write_text("chr1\ts\tCDS\t1\t9\t.\t+\t0\tID=g1.2\n")
    check = SourceGeneSelection(path)
    check.observe("g1.2", "gene1", True)
    assert check.validate()["source_gene_unresolved_records"] == 1
    assert check.validate()["source_gene_check"] == "unresolved"
    assert not check.validate()["source_gene_complete"]


def test_missing_gene_features_preserve_explicit_parent_ownership(tmp_path):
    path = tmp_path / "source.gff"
    path.write_text(
        "chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=gene%3Ag1\n"
        "chr1\ts\tmRNA\t30\t38\t.\t+\t.\tID=t2;Parent=gene%3Ag1\n"
        "chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=t3;Parent=g2\n")
    check = SourceGeneSelection(path)
    check.observe("t1", "g1", True)
    check.observe("t2", "g1", False)
    check.observe("t3", "g2", True)
    stats = check.validate()
    assert stats["source_gene_complete"]
    assert stats["source_gene_resolved_records"] == 3
    split = SourceGeneSelection(path)
    split.observe("t1", "first", True)
    split.observe("t2", "second", True)
    with pytest.raises(ValueError, match="retained_isoform_gene_groups=1"):
        split.validate()
    merge = SourceGeneSelection(path)
    merge.observe("t1", "merged", True)
    merge.observe("t3", "merged", False)
    with pytest.raises(ValueError, match="distinct_source_gene_merges=1"):
        merge.validate()


@pytest.mark.parametrize("axis", [("chr2", "+"), ("chr1", "-")])
def test_missing_parent_conflicting_axes_fail_without_guessing(tmp_path, axis):
    path = tmp_path / "source.gff"
    path.write_text("chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n"
                    f"{axis[0]}\ts\tmRNA\t1\t9\t.\t{axis[1]}\t.\tID=t2;Parent=g1\n")
    with pytest.raises(ValueError, match="Conflicting axes"):
        SourceGeneSelection(path)


@pytest.fixture
def reused_accession(tmp_path):
    path = tmp_path / "reused_accession.gff"
    path.write_text(
        "chr1\ts\tgene\t1\t9\t.\t+\t.\tID=g1;Name=LOC1\n"
        "chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ts\tCDS\t1\t9\t.\t+\t0\tID=cds-P1;Parent=t1;protein_id=P1\n"
        "chr2\ts\tgene\t11\t29\t.\t-\t.\tID=g2;Name=LOC1\n"
        "chr2\ts\tmRNA\t11\t29\t.\t-\t.\tID=t2;Parent=g2\n"
        "chr2\ts\tCDS\t11\t13\t.\t-\t0\tID=cds-P1;Parent=t2;protein_id=P1\n"
        "chr2\ts\tCDS\t24\t29\t.\t-\t0\tID=cds-P1;Parent=t2;protein_id=P1\n")
    return path


@pytest.mark.parametrize("location", ["1..9", "<1..>9", "join(1..9)", "order(1..9)"])
def test_exact_source_location_disambiguates_reused_protein_accession(reused_accession, location):
    check = SourceGeneSelection(reused_accession)
    assert check.header_roots(f"lcl|chr1_cds_P1_1 [protein_id=P1] [gene=LOC1] [location={location}]") == {"g1"}
    assert check.header_roots("lcl|chr2_cds_P1_2 [protein_id=P1] [location=complement(join(11..13,24..29))]") == {"g2"}


def test_reused_accession_cannot_hide_a_distinct_gene_merge(reused_accession):
    check = SourceGeneSelection(reused_accession)
    check.observe("lcl|chr1_cds_P1_1 [location=1..9]", "merged", True)
    check.observe("lcl|chr2_cds_P1_2 [location=complement(join(11..13,24..29))]", "merged", False)
    with pytest.raises(ValueError, match="distinct_source_gene_merges=1"):
        check.validate()


def test_exact_locations_still_reject_two_selected_isoforms(reused_accession):
    check = SourceGeneSelection(reused_accession)
    header = "lcl|chr1_cds_P1_1 [location=1..9]"
    check.observe(header, "first", True)
    check.observe(header, "second", True)
    with pytest.raises(ValueError, match="retained_isoform_gene_groups=1"):
        check.validate()


def test_correct_reused_accession_selection_has_complete_ownership(reused_accession):
    check = SourceGeneSelection(reused_accession)
    check.observe("lcl|chr1_cds_P1_1 [location=1..9]", "first", True)
    check.observe("lcl|chr2_cds_P1_2 [location=complement(join(11..13,24..29))]", "second", True)
    assert check.validate()["source_gene_complete"]
    assert check.validate()["source_gene_unresolved_records"] == 0


@pytest.mark.parametrize("header", [
    "lcl|chr1_cds_P1_1", "lcl|chr1_cds_P1_1 [location=1..8]",
    "lcl|chr1_cds_P1_1 [location=complement(1..9)]",
    "lcl|chr3_cds_P1_1 [location=1..9]",
    "lcl|chr1_cds_P1_1 [location=9..1]",
    "lcl|chr1_cds_P1_1 [location=0..9]",
    "lcl|chr1_cds_P1_1 [location=join(1..9,complement(24..29))]",
    "chr1_cds_P1_1 [protein_id=P1] [location=1..9]",
])
def test_inexact_or_unsupported_locations_do_not_guess_an_owner(reused_accession, header):
    check = SourceGeneSelection(reused_accession)
    check.observe(header, "first", True)
    assert not check.validate()["source_gene_complete"]
    assert check.validate()["source_gene_unresolved_records"] == 1


def test_identifier_location_conflict_is_rejected(reused_accession):
    check = SourceGeneSelection(reused_accession)
    with pytest.raises(ValueError, match="identifier and exact GFF location disagree"):
        check.header_roots("lcl|chr2_cds_unused_1 [transcript_id=t1] [location=complement(join(11..13,24..29))]")


def test_exact_shared_coordinates_do_not_define_one_gene(reused_accession):
    with reused_accession.open("a") as handle:
        handle.write("chr1\ts\tgene\t1\t9\t.\t+\t.\tID=g3\n"
                     "chr1\ts\tCDS\t1\t9\t.\t+\t0\tID=cds-P2;Parent=g3;protein_id=P2\n")
    check = SourceGeneSelection(reused_accession)
    check.observe("lcl|chr1_cds_anonymous_1 [location=1..9]", "first", True)
    assert not check.validate()["source_gene_complete"]

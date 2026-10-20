import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
import validate_longest_cds_selection as longest
from validate_source_gene_selection import SourceGeneSelection


def test_source_check_retains_ownership_with_padded_gff_attributes(tmp_path):
    path = tmp_path / "padded.gff"
    path.write_text(
        "chr1\tEVM\tgene\t1\t9\t.\t+\t.\t ID=g1 ; Name=LOC1 \n"
        "chr1\tEVM\tmRNA\t1\t9\t.\t+\t.\tID=t1; Name=t1; Parent= g1 \n"
        "chr1\tEVM\tCDS\t1\t9\t.\t+\t0\tID=c1; Name=c1; Parent= t1 \n"
    )
    check = SourceGeneSelection(path)
    assert check.header_roots("t1") == {"g1"}
    assert check.header_roots("c1") == {"g1"}
    check.observe("t1", "selected", True)
    stats = check.validate()
    assert stats["source_gene_complete"]
    assert stats["source_gene_unresolved_records"] == 0


def test_padded_attributes_do_not_hide_multiple_explicit_parents(tmp_path):
    path = tmp_path / "ambiguous.gff"
    path.write_text(
        "chr1\tEVM\tgene\t1\t9\t.\t+\t.\tID=g1\n"
        "chr1\tEVM\tgene\t1\t9\t.\t+\t.\tID=g2\n"
        "chr1\tEVM\tmRNA\t1\t9\t.\t+\t.\tID=t1; Parent= g1, g2 \n"
    )
    check = SourceGeneSelection(path)
    assert check.header_roots("t1") == {"g1", "g2"}
    check.observe("t1", "selected", True)
    assert not check.validate()["source_gene_complete"]


def test_attribute_padding_preserves_encoded_spaces_and_opaque_id_commas(tmp_path):
    path = tmp_path / "opaque.gff"
    path.write_text(
        "chr1\tEVM\tgene\t1\t9\t.\t+\t.\t ID=gene%20one \n"
        "chr1\tEVM\tmRNA\t1\t9\t.\t+\t.\t ID=t1,description ; Parent=gene%20one \n"
        "chr1\tEVM\tCDS\t1\t9\t.\t+\t0\tID=c1; Parent=t1%2Cdescription\n"
    )
    check = SourceGeneSelection(path)
    assert check.header_roots("t1,description") == {"gene one"}
    assert check.header_roots("c1") == {"gene one"}
    assert "description" not in check.nodes_by_alias


def test_comma_decorated_ids_keep_unique_explicit_name_and_parent_ancestry(tmp_path):
    path = tmp_path / "source.gff"
    path.write_text(
        "chr1\tHelixer\tgene\t1\t9\t.\t+\t.\tID=g1\n"
        "chr1\tHelixer\tgene\t20\t28\t.\t+\t.\tID=g2\n"
        "chr1\tHelixer\tmRNA\t1\t9\t.\t+\t.\tID=t1,shared description;Parent=g1;Name=t1\n"
        "chr1\tHelixer\tmRNA\t20\t28\t.\t+\t.\tID=t2,shared description;Parent=g2;Name=t2\n"
        "chr1\tHelixer\tCDS\t1\t9\t.\t+\t0\tParent=t1\n"
        "chr1\tHelixer\tCDS\t20\t28\t.\t+\t0\tParent=t2\n")
    check = SourceGeneSelection(path)
    assert check.header_roots("t1") == {"g1"}
    assert check.header_roots("t2") == {"g2"}
    assert "shared description" not in check.parents
    check.observe("t1", "one", True)
    check.observe("t2", "two", True)
    assert check.validate()["source_gene_complete"]
    check.observe("t2", "one", False)
    with pytest.raises(ValueError, match="distinct_source_gene_merges=1"):
        check.validate()


@pytest.mark.parametrize("attributes", [
    "ID=t1;Parent=g1,g2",
    "ID=t1%2Cshared;Parent=g1,g2;Name=t1",
])
def test_single_id_parsing_keeps_genuine_multiple_gene_parents_ambiguous(tmp_path, attributes):
    path = tmp_path / "source.gff"
    path.write_text(
        "chr1\ts\tgene\t1\t9\t.\t+\t.\tID=g1\n"
        "chr1\ts\tgene\t1\t9\t.\t+\t.\tID=g2\n"
        "chr1\ts\tmRNA\t1\t9\t.\t+\t.\t" + attributes + "\n")
    check = SourceGeneSelection(path)
    assert check.header_roots("t1") == {"g1", "g2"}
    check.observe("t1", "one", True)
    assert not check.validate()["source_gene_complete"]
    assert check.validate()["source_gene_unresolved_selected_records"] == 1


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


@pytest.fixture
def gtf_annotation(tmp_path):
    path = tmp_path / "source.gtf"
    path.write_text(
        'chr1\ts\tCDS\t1\t9\t.\t+\t0\tgene_id "g1"; transcript_id "t1.1"; protein_id "p1"; note "quoted; metadata";\n'
        'chr1\ts\ttranscript\t1\t9\t.\t+\t.\tgene_id "g1"; transcript_id "t1.1";\n'
        'chr1\ts\ttranscript\t1\t12\t.\t+\t.\tgene_id "g1"; transcript_id "t1.2";\n'
        'chr2\ts\tCDS\t1\t9\t.\t-\t0\tgene_id "g2,opaque%20"; transcript_id "t2.1";\n')
    return path


def test_gtf_explicit_hierarchy_resolves_transcripts_proteins_and_plain_gene_tags(gtf_annotation):
    check = SourceGeneSelection(gtf_annotation)
    for header in ("t1.1", "p1", "unknown gene=g1"):
        assert check.header_roots(header) == {"g1"}
    assert check.header_roots("t2.1") == {"g2,opaque%20"}
    assert check.header_roots("t1") == set()  # Never guess a numeric suffix.
    assert check.header_roots("t1.1 gene=g2,opaque%20") == {"g1"}
    check.observe("t1.1", "one", True)
    check.observe("t1.2", "one", False)
    check.observe("unknown gene=g2,opaque%20", "two", True)
    assert check.validate()["source_gene_complete"]


def test_gtf_independent_hierarchy_rejects_splits_and_merges(gtf_annotation):
    check = SourceGeneSelection(gtf_annotation)
    check.observe("t1.1", "one", True)
    check.observe("t1.2", "two", True)
    with pytest.raises(ValueError, match="retained_isoform_gene_groups=1"):
        check.validate()
    check = SourceGeneSelection(gtf_annotation)
    check.observe("t1.1", "one", True)
    check.observe("t2.1", "one", False)
    with pytest.raises(ValueError, match="distinct_source_gene_merges=1"):
        check.validate()


@pytest.mark.parametrize("attributes", [
    'transcript_id "g1.2";',
    'transcript_id "t1"; gene_name "g1";',
])
def test_gtf_missing_gene_id_stays_unresolved(tmp_path, attributes):
    path = tmp_path / "missing.gtf"
    path.write_text('chr1\ts\tCDS\t1\t9\t.\t+\t0\t' + attributes + '\n')
    check = SourceGeneSelection(path)
    check.observe(attributes.split('"')[1], "one", True)
    assert not check.validate()["source_gene_complete"]


@pytest.mark.parametrize("second", [
    'chr1\ts\texon\t1\t9\t.\t+\t.\tgene_id "g2"; transcript_id "t1";\n',
    'chr2\ts\texon\t1\t9\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n',
    'chr1\ts\texon\t1\t9\t.\t-\t.\tgene_id "g1"; transcript_id "t2";\n',
])
def test_gtf_conflicting_explicit_hierarchy_is_rejected(tmp_path, second):
    path = tmp_path / "conflict.gtf"
    path.write_text('chr1\ts\tCDS\t1\t9\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n' + second)
    with pytest.raises(ValueError, match="Conflicting"):
        SourceGeneSelection(path)


@pytest.mark.parametrize("attributes", [
    'gene_id "g1"; gene_id "g2"; transcript_id "t1";',
    'gene_id "g1"; transcript_id "t1"; transcript_id "t2";',
    'gene_id "g1"; transcript_id "g1";',
    'gene_id "g1" transcript_id "t1";',
])
def test_gtf_malformed_or_ambiguous_identity_is_rejected(tmp_path, attributes):
    path = tmp_path / "invalid.gtf"
    path.write_text('chr1\ts\tCDS\t1\t9\t.\t+\t0\t' + attributes + '\n')
    with pytest.raises(ValueError, match="GTF"):
        SourceGeneSelection(path)


def test_absent_gff_attributes_remain_unresolved(tmp_path):
    path = tmp_path / "empty.gff"
    path.write_text("chr1\ts\tCDS\t1\t9\t.\t+\t0\t.\n")
    check = SourceGeneSelection(path)
    check.observe("unknown", "one", True)
    assert not check.validate()["source_gene_complete"]


def test_gtf_gene_transcript_collision_across_rows_is_rejected(tmp_path):
    path = tmp_path / "collision.gtf"
    path.write_text(
        'chr1\ts\tCDS\t1\t9\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\ts\tgene\t1\t9\t.\t+\t.\tgene_id "t1";\n')
    with pytest.raises(ValueError, match="identity collision"):
        SourceGeneSelection(path)

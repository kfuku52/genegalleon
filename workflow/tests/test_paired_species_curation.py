import gzip
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from curate_paired_species_inputs import curate_pair, fingerprint, inspect_pair
from format_species_annotation.reference import gff_reference_mapping


def pair(tmp_path, *, extras="", missing="absent", mixed=False):
    cds = tmp_path / "cds.fa"
    gff = tmp_path / "annotation.gff"
    genome = tmp_path / "genome.fa"
    cds.write_text(">Test_species_good literal description\natgaaa\n>Test_species_bad\nATGCCC\n" + extras)
    genome.write_text(">chr1\nATGAAAATGCCC\n")
    gff.write_text("##gff-version 3\n##sequence-region chr1 1 12\n##sequence-region absent 1 12\n"
                   "chr1\tx\tgene\t1\t6\t.\t+\t.\tID=good\n"
                   "chr1\tx\tmRNA\t1\t6\t.\t+\t.\tID=good.t;Parent=good\n"
                   "chr1\tx\tCDS\t1\t6\t.\t+\t0\tParent=good.t\n"
                   + f"{missing}\tx\tgene\t1\t6\t.\t+\t.\tID=bad\n"
                   + f"{missing}\tx\tmRNA\t1\t6\t.\t+\t.\tID=bad.t;Parent=bad\n"
                   + f"{missing}\tx\tCDS\t1\t6\t.\t+\t0\tParent=bad.t\n"
                   + ("chr1\tx\tCDS\t7\t12\t.\t+\t0\tParent=bad.t\n" if mixed else ""))
    return cds, gff, genome


def policy(tmp_path, paths, records=None, approve=True):
    path = tmp_path / "decision.json"
    path.write_text(json.dumps(dict(schema_version=1, species="Test_species", decision_basis="Explicit source review",
                                    input_sha256={key: fingerprint(value)["sha256"] for key, value in zip(("cds", "gff", "genome"), paths, strict=True)},
                                    exclude_missing_reference_annotations=approve,
                                    records=records if records is not None else [dict(cds_id="Test_species_bad", action="exclude", reason="missing_genome_reference")])))
    return path


def run(tmp_path, paths, decision):
    return curate_pair("Test_species", *paths, decision, tmp_path / "curated")


def test_approved_curation_preserves_sequences_and_coordinates_and_default_error(tmp_path):
    paths = pair(tmp_path)
    with pytest.raises(ValueError, match="absent from genome"):
        gff_reference_mapping(paths[1], paths[2])
    source_bytes = [path.read_bytes() for path in paths]
    result = run(tmp_path, paths, policy(tmp_path, paths))
    assert result["excluded_cds_ids"] == ["Test_species_bad"]
    assert result["remaining_cds_records"] == 1
    assert result["excluded_gff_features"] == 3
    assert gzip.decompress(Path(result["cds_output"]["path"]).read_bytes()).decode() == ">Test_species_good literal description\natgaaa\n"
    expected = "".join(line for line in paths[1].read_text().splitlines(keepends=True) if "absent" not in line)
    assert gzip.decompress(Path(result["gff_output"]["path"]).read_bytes()).decode() == expected
    assert [path.read_bytes() for path in paths] == source_bytes
    with pytest.raises(FileExistsError):
        run(tmp_path, paths, policy(tmp_path, paths))


@pytest.mark.parametrize("records,approve,match", [
    ([], True, "exclusion IDs differ"),
    ([dict(cds_id="Test_species_good", action="exclude", reason="missing_genome_reference")], True, "exclusion IDs differ"),
    (None, False, "explicit approval"),
    ([dict(cds_id="Test_species_nonexistent", action="exclude", reason="missing_genome_reference")], True, "Missing or duplicate"),
])
def test_no_unapproved_or_wrong_exclusions(tmp_path, records, approve, match):
    paths = pair(tmp_path)
    with pytest.raises(ValueError, match=match):
        run(tmp_path, paths, policy(tmp_path, paths, records, approve))
    assert not (tmp_path / "curated").exists()


@pytest.mark.parametrize("key", [0, 1, 2])
def test_changed_source_invalidates_decision(tmp_path, key):
    paths = pair(tmp_path)
    decision = policy(tmp_path, paths)
    with paths[key].open("a") as handle:
        handle.write("\n")
    with pytest.raises(ValueError, match="exact CDS/GFF/genome hashes"):
        run(tmp_path, paths, decision)
    assert not (tmp_path / "curated").exists()


def test_mixed_parent_is_not_silently_split(tmp_path):
    paths = pair(tmp_path, mixed=True)
    with pytest.raises(ValueError, match="Parent crosses"):
        run(tmp_path, paths, policy(tmp_path, paths))
    assert not (tmp_path / "curated").exists()


def test_unique_reference_alias_is_retained_not_excluded(tmp_path):
    paths = pair(tmp_path, missing="chr1")
    paths[2].write_text(paths[2].read_text().replace(">chr1", ">lcl|chr1"))
    text = paths[1].read_text().replace("##sequence-region absent 1 12\n", "")
    paths[1].write_text(text)
    report = run(tmp_path, paths, policy(tmp_path, paths, []))
    assert not report["excluded_cds_ids"]
    assert report["remaining_cds_records"] == 2
    assert "\nlcl|chr1\t" in gzip.decompress(Path(report["gff_output"]["path"]).read_bytes()).decode()


def test_unmatched_CDS_requires_exact_retention_flag(tmp_path):
    paths = pair(tmp_path, extras=">Test_species_unmapped\nATGTTT\n")
    decision = policy(tmp_path, paths)
    with pytest.raises(ValueError, match="Unapproved CDS without GFF"):
        run(tmp_path, paths, decision)
    data = json.loads(decision.read_text())
    data["records"].append(dict(cds_id="Test_species_unmapped", action="retain_and_flag", reason="missing_gff_counterpart"))
    decision.write_text(json.dumps(data))
    result = run(tmp_path, paths, decision)
    assert result["remaining_cds_records"] == 2
    assert result["retained_source_exceptions"] == data["records"][1:]
    assert "Test_species_unmapped" in gzip.decompress(Path(result["cds_output"]["path"]).read_bytes()).decode()


def test_coding_span_flag_requires_observed_exact_lengths(tmp_path):
    paths = pair(tmp_path)
    paths[0].write_text(paths[0].read_text().replace("atgaaa", "atgaaaatg"))
    decisions = [dict(cds_id="Test_species_bad", action="exclude", reason="missing_genome_reference"),
                 dict(cds_id="Test_species_good", action="retain_and_flag", reason="coding_span_conflict", cds_length=9, gff_coding_span_length=6)]
    result = run(tmp_path, paths, policy(tmp_path, paths, decisions))
    assert result["retained_source_exceptions"] == decisions[1:]


def test_ambiguous_genome_alias_does_not_enable_exclusion(tmp_path):
    paths = pair(tmp_path)
    paths[2].write_text(">chr1 OriSeqID=shared\nATGAAAATGCCC\n>chr2 OriSeqID=shared\nATGAAA\n")
    with pytest.raises(ValueError, match="Ambiguous genome FASTA alias"):
        inspect_pair("Test_species", *paths)


def test_mapped_out_of_bounds_is_not_hidden_by_exclusion(tmp_path):
    paths = pair(tmp_path)
    paths[1].write_text(paths[1].read_text().replace("chr1\tx\tCDS\t1\t6", "chr1\tx\tCDS\t1\t99"))
    with pytest.raises(ValueError, match="coordinates exceed"):
        run(tmp_path, paths, policy(tmp_path, paths))
    assert not (tmp_path / "curated").exists()


def test_duplicate_formatted_CDS_and_cycles_rejected(tmp_path):
    paths = pair(tmp_path, extras=">Test_species_good\nATGAAA\n")
    with pytest.raises(ValueError, match="unique"):
        inspect_pair("Test_species", *paths)
    paths = pair(tmp_path)
    paths[1].write_text(paths[1].read_text().replace("ID=bad\n", "ID=loop;Parent=bad.t\n").replace("Parent=bad\n", "Parent=loop\n"))
    with pytest.raises(ValueError, match="Cyclic"):
        inspect_pair("Test_species", *paths)


def test_curated_pair_remains_a_native_formatting_input(tmp_path):
    import format_species_inputs as fsi

    paths = pair(tmp_path)
    report = run(tmp_path, paths, policy(tmp_path, paths))
    task = dict(species_key="Test_species", species_prefix="Test_species", provider="local",
                cds_path=Path(report["cds_output"]["path"]), gff_path=Path(report["gff_output"]["path"]),
                genome_path=paths[2], gbff_path=None, gene_grouping_mode="rescue_overlap", gff_repair_mode="safe", format_strict=True)
    output = tmp_path / "native"
    output.mkdir()
    cds = fsi.format_cds(task, output, False, False)
    gff = fsi.format_gff(task, output, False, False, formatted_cds_path=cds["output_path"])
    assert cds["after_count"] == 1
    assert list(fsi.iter_fasta_records(cds["output_path"])) == [("Test_species_good", "ATGAAA")]
    assert gff["status"] == "write"

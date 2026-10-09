import contextlib
import gzip
import io
import os
import sys
import tarfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from format_species_annotation import cds_normalisation, genbank
from format_species_annotation.genome_intervals import genome_intervals
from format_species_annotation.reference import genome_fragments, genome_reference_index
from gff_attribute_syntax import normalise_attributes, normalise_line, validate_attributes


@pytest.mark.parametrize("kind", ["plain", "gzip", "archive"])
def test_wrapped_unwrapped_aliases_and_selected_disk_bytes(tmp_path, kind):
    content = ">acc OriSeqID=chr1 Len=12\r\natg AAA\r\ncccTAG\r\n>lcl|chr2\nACGTACGT"
    path = tmp_path / "genome.fa"
    if kind == "gzip":
        path = path.with_suffix(".fa.gz")
        with gzip.open(path, "wt") as handle:
            handle.write(content)
    elif kind == "archive":
        path = path.with_suffix(".fa.tar.gz")
        with tarfile.open(path, "w:gz") as archive:
            member = tarfile.TarInfo("genome.fa")
            member.size = len(content.encode())
            archive.addfile(member, io.BytesIO(content.encode()))
    else:
        path.write_text(content)
    with genome_intervals(path, [("chr1", 0, 3), ("chr1", 2, 6), ("chr2", 3, 5)]) as reader:
        assert reader.fetch("acc", 1, 5) == "TGAA"
        assert reader.fetch("lcl|chr2", 3, 5) == "TA"
        assert reader.lengths == (12, 8)
        assert len(reader["chr1"]) == 12
        assert reader["chr1"][0:3] == "ATG"
        reader.scratch.seek(0, 2)
        assert reader.scratch.tell() == 8
        with pytest.raises(ValueError, match="outside declared"):
            reader.fetch("acc", 5, 9)
        with pytest.raises(ValueError, match="outside genome bounds"):
            reader.fetch("acc", 0, 13)


@pytest.mark.parametrize("chunk", [1, 2, 3, 7, 16, 1024 * 1024])
def test_fragment_boundaries_preserve_headers_and_sequence(tmp_path, chunk):
    path = tmp_path / "genome.fa"
    path.write_text(">chr1 long header\r\nacgt\ntgca\n>chr2\nAa Tt\n>chr3")
    records = []
    for kind, text in genome_fragments(path, chunk):
        if kind == "header":
            records.append([text, ""])
        else:
            records[-1][1] += text
    assert records == [["chr1 long header", "acgttgca"], ["chr2", "AaTt"], ["chr3", ""]]


def test_memory_is_independent_of_long_contig_length(tmp_path):
    path = tmp_path / "long.fa"
    with path.open("w") as handle:
        handle.write(">chr1\n")
        for _ in range(32):
            handle.write("a" * 1024 * 1024)
    with genome_intervals(path, [("chr1", 31 * 1024 * 1024, 31 * 1024 * 1024 + 9)]) as reader:
        assert reader.fetch("chr1", 31 * 1024 * 1024, 31 * 1024 * 1024 + 9) == "A" * 9
        assert reader.index["chr1"] == 32 * 1024 * 1024
        reader.scratch.seek(0, 2)
        assert reader.scratch.tell() == 9
        assert not any(isinstance(value, str) for value in reader.index.values())


@pytest.mark.parametrize("content, error", [
    (">a OriSeqID=x\nATG\n>b OriSeqID=x\nCCC\n", "Ambiguous"),
    (">a\nATG\n>a\nCCC\n", "duplicate"),
    (">a OriSeqID=x Len=4\nATG\n", "length disagrees"),
    ("ATG\n>a\nCCC\n", "precedes"),
])
def test_invalid_reference_evidence_still_fails(tmp_path, content, error):
    path = tmp_path / "bad.fa"
    path.write_text(content)
    with pytest.raises(ValueError, match=error):
        with genome_intervals(path, [("a", 0, 3)]):
            pass
    with pytest.raises(ValueError, match=error):
        genome_reference_index(path)


def test_source_change_during_use_cannot_finish(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_text(">chr1\nATGAAA\n")
    initial = path.stat()
    with pytest.raises(OSError, match="changed during"):
        with genome_intervals(path, [("chr1", 0, 3)]) as reader:
            assert reader.fetch("chr1", 0, 3) == "ATG"
            path.write_text(">chr1\nATGCCC\n")
            # Some HPC scratch mounts have coarse timestamp resolution.
            os.utime(path, ns=(initial.st_atime_ns, initial.st_mtime_ns + 2_000_000_000))
    assert reader.scratch.closed


def test_non_ascii_sequence_cannot_expand_and_shift_interval_offsets(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_text(">chr1\naßtc\n")
    with pytest.raises(ValueError, match="ASCII"):
        with genome_intervals(path, [("chr1", 0, 4)]):
            pass


def test_mixed_strand_exons_match_existing_biological_order(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_text(">chr1\nATGAAACCCGGGTTT\n")
    gff = tmp_path / "models.gff"
    gff.write_text("chr1\tsrc\tgene\t1\t15\t.\t-\t.\tID=g\n"
                   "chr1\tsrc\tmRNA\t1\t15\t.\t-\t.\tID=t;Parent=g\n"
                   "chr1\tsrc\tCDS\t1\t6\t.\t-\t0\tParent=t\n"
                   "chr1\tsrc\tCDS\t10\t15\t.\t-\t0\tParent=t\n")
    task = dict(provider="direct", species_prefix="Example_species", species_key="Example_species",
                gff_path=gff, genome_path=path)
    assert list(genbank.derive_cds_records_from_gff_and_genome(task)) == [("t [gene=g]", "AAACCCTTTCAT")]


def test_normalisation_scratch_excludes_gene_spans_and_introns(tmp_path, monkeypatch):
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nATGAAATAA" + "N" * 1000000 + "\n")
    gff = tmp_path / "models.gff"
    gff.write_text("chr1\tsrc\tgene\t1\t1000009\t.\t+\t.\tID=g\n"
                   "chr1\tsrc\tmRNA\t1\t1000009\t.\t+\t.\tID=t;Parent=g\n"
                   "chr1\tsrc\texon\t1\t9\t.\t+\t.\tParent=t\n"
                   "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tParent=t\n")
    sizes, handles = [], []
    @contextlib.contextmanager
    def observed(*args, **kwargs):
        with genome_intervals(*args, **kwargs) as reader:
            reader.scratch.seek(0, 2)
            sizes.append(reader.scratch.tell())
            handles.append(reader.scratch)
            yield reader
    monkeypatch.setattr(cds_normalisation, "genome_intervals", observed)
    task = dict(provider="direct", species_prefix="Example_species", species_key="Example_species",
                gff_path=gff, genome_path=genome, _normalisation_scratch=tmp_path)
    records = list(cds_normalisation.iter_normalised_cds_records(task))
    assert [(header, sequence) for header, sequence, _decision in records] == [("t [gene=g]", "ATGAAATAA")]
    assert sizes == [9]
    assert all(handle.closed for handle in handles)


@pytest.mark.parametrize("kind", ["gbff", "embl"])
def test_assembly_only_insdc_never_constructs_sequence_records(tmp_path, monkeypatch, kind):
    path = tmp_path / ("assembly." + kind)
    path.write_text("LOCUS       chr1\nFEATURES             Location/Qualifiers\n"
                    "     source          1..6\n                     /note=\"CDS is absent\"\nORIGIN\n        1 atgaaa\n//\n")
    monkeypatch.setattr(genbank, "iter_genbank_records", lambda *a: pytest.fail("Parsed assembly chromosome"))
    assert list(genbank.derive_cds_records_from_gbff(dict(gbff_path=path))) == []
    path.write_text("FT   CDS             1..6\n" if kind == "embl" else "     CDS             1..6\n")
    assert genbank.has_insdc_coding_features(path)


def test_evidenced_gwh_confidence_flag_preserves_identifiers():
    text = "ID=EVM0000018;Accession=GWHGBFXR003635.1;HC;transl_table=1"
    result = normalise_attributes(text, "EVM", "gene")
    assert result == "ID=EVM0000018;Accession=GWHGBFXR003635.1;gwh_gene_confidence=HC;transl_table=1"
    validate_attributes(result)
    for source, feature, invalid in [("other", "gene", text), ("EVM", "mRNA", text),
                                     ("EVM", "gene", text.replace("GWHGBFXR003635.1", "unknown"))]:
        with pytest.raises(ValueError, match="Unrecoverable"):
            normalise_attributes(invalid, source, feature)


def test_fragaria_note_semicolon_is_metadata_not_a_parent_change():
    note = "ID=g1;Note=nucleosome assly protein 1;1-like, partial [Fragaria vesca subsp. vesca]"
    result = normalise_attributes(note, "maker", "gene")
    assert result.startswith("ID=g1;Note=nucleosome assly protein 1%3B1-like%2C partial")
    changes = []
    line = "chr1\tmaker\tgene\t1\t3\t.\t+\t.\t" + note + "\n"
    assert normalise_line(line, "fixture", 1, changes).split("\t")[:8] == line.split("\t")[:8]
    assert changes[0]["reason"] == "escaped_NAP_metadata_semicolon"
    with pytest.raises(ValueError, match="Unrecoverable"):
        normalise_attributes(note.replace("Note=", "Parent="), "maker", "gene")


@pytest.mark.parametrize("note", [
    "nucleosome assly protein 1;4-like [Fragaria vesca subsp. vesca]",
    "nucleosome assly protein 1;2 [Fragaria vesca subsp. vesca]",
    "nucleosome assly protein 1;3-like isoform X2 [Fragaria vesca subsp. vesca]",
    "PROF_FRAAN RecName: Full=Profilin; AltName: Allergen=Fra a 4",
    "URT1_FRAAN RecName: Full=UDP-rhamnosyltransferase; Short=FaRT1; AltName: Full=Glycosyltransferase 4",
    "POLX_TOBAC RecName: Full=Pol polyprotein; Includes: RecName: Full=Protease",
    "TPIC_FRAAN RecName: Full=Triosephosphate isomerase; Short=TIM; Flags: Precursor",
    "endo-1,3;1,4-beta-D-glucanase-like [Fragaria vesca subsp. vesca]",
])
def test_publisher_metadata_continuations_preserve_one_note(note):
    from urllib.parse import unquote
    text = "ID=g;Note=" + note + ";Parent=p"
    result = normalise_attributes(text, "maker", "gene")
    validate_attributes(result)
    assert result.startswith("ID=g;Note=") and result.endswith(";Parent=p")
    assert result.count(";") == 2
    assert unquote(result.split(";", 2)[1][5:]) == note

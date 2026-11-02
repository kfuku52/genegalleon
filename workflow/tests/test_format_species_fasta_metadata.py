"""FASTA metadata readers must retain full-reader boundaries and errors."""

import bz2
import gzip
import io
import sys
import tarfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from format_species_annotation import common  # noqa: E402
from format_species_annotation.tasks import discover_generic_species_dir_tasks  # noqa: E402
from format_species_common import is_probable_genome_filename  # noqa: E402


def write_fasta(tmp_path, raw, transport):
    suffix = {"raw": "", "gz": ".gz", "bz2": ".bz2", "tar-gz": ".tar.gz",
              "tar-bz2": ".tar.bz2", "tar-single": ".tar.gz", "magic": ""}[transport]
    path = tmp_path / ("source.fa" + suffix)
    if transport.startswith("tar"):
        mode = "w:bz2" if transport == "tar-bz2" else "w:gz"
        with tarfile.open(path, mode) as archive:
            if transport != "tar-single":
                ignored = tarfile.TarInfo("ignored.txt")
                ignored.size = 8
                archive.addfile(ignored, io.BytesIO(b"not DNA\n"))
            member = tarfile.TarInfo("single.data" if transport == "tar-single" else "member.fa")
            member.size = len(raw)
            archive.addfile(member, io.BytesIO(raw))
    else:
        encoded = gzip.compress(raw) if transport in {"gz", "magic"} else bz2.compress(raw) if transport == "bz2" else raw
        path.write_bytes(encoded)
    return path


def reference_count(path):
    count, first_name = 0, ""
    for header, _sequence in common.iter_fasta_records(path):
        count += 1
        if first_name == "":
            first_name = common.first_token(header)
    return count, first_name


def reference_overlap(path, identifiers):
    for header, _sequence in common.iter_fasta_records(path):
        identifier = common.first_token(header)
        if identifier in identifiers or (identifier.startswith("lcl|") and identifier[4:] in identifiers):
            return True
    return False


def observe(call):
    try:
        return call(), None
    except Exception as error:
        return None, (type(error), str(error))


@pytest.mark.parametrize("transport", ["raw", "gz", "bz2", "tar-gz", "tar-bz2", "tar-single", "magic"])
@pytest.mark.parametrize("payload", [
    b"", b"sequence before any header\n", b">\nAC\n>\n\n>001 description\nA C\tG\n>NA\nT\n",
    b"before\r\n>lcl|chr1 description\r\na c\tgt\r\n\r\n>chr2\r\nA\r\n>chr2\nT",
    (">あ Unicode\nAC\u2003GT\n>NULL\nTT\n>long " + "header " * 16000 + "\n" + "ac gt\t" * 16000 + "\n").encode(),
])
def test_metadata_matches_complete_record_reader_across_transports(tmp_path, transport, payload):
    path = write_fasta(tmp_path, payload, transport)
    assert common.count_fasta_records(path) == reference_count(path)
    for identifiers in ({"chr1"}, {"001"}, {"あ"}, {"not-present"}):
        gff = tmp_path / "annotation.gff"
        gff.write_text("".join(f"{identifier}\tsrc\tgene\t1\t2\t.\t+\t.\tID=g\n" for identifier in identifiers))
        assert common.fasta_headers_overlap_gff_seqids(path, gff) == reference_overlap(path, identifiers)


@pytest.mark.parametrize("transport", ["raw", "gz", "bz2", "tar-gz"])
@pytest.mark.parametrize("bad_position", ["preamble", "first", "late"])
def test_metadata_preserves_strict_utf8_failures_at_consumed_record_boundaries(tmp_path, transport, bad_position):
    records = [b"before\n", b">match\nAC\n", b">late\n" + b"ACGT\n" * 4000]
    index = {"preamble": 0, "first": 1, "late": 2}[bad_position]
    records[index] += b"\xff\n"
    path = write_fasta(tmp_path, b"".join(records), transport)
    assert observe(lambda: common.count_fasta_records(path)) == observe(lambda: reference_count(path))
    gff = tmp_path / "annotation.gff"
    gff.write_text("match\tsrc\tgene\t1\t2\t.\t+\t.\tID=g\n")
    assert observe(lambda: common.fasta_headers_overlap_gff_seqids(path, gff)) == observe(
        lambda: reference_overlap(path, {"match"}))


@pytest.mark.parametrize("truncated", [True, False])
def test_gzip_members_and_corrupt_tail_keep_existing_errors(tmp_path, truncated):
    first = gzip.compress(b">first\nAC\n>second\nGT\n")
    second = gzip.compress(b">last\n" + b"ACGT\n" * 4000)
    path = tmp_path / "source.fa.gz"
    path.write_bytes(first + (second[:-4] if truncated else second))
    assert observe(lambda: common.count_fasta_records(path)) == observe(lambda: reference_count(path))
    gff = tmp_path / "annotation.gff"
    gff.write_text("first\tsrc\tgene\t1\t2\t.\t+\t.\tID=g\n")
    assert observe(lambda: common.fasta_headers_overlap_gff_seqids(path, gff)) == observe(
        lambda: reference_overlap(path, {"first"}))


@pytest.mark.parametrize("member", ["only", "later"])
@pytest.mark.parametrize("match", ["first", "last"])
def test_crc_corruption_preserves_complete_and_early_read_errors_and_closure(tmp_path, monkeypatch, member, match):
    payload = b">first\n" + b"ACGT\n" * 5000 + b">last\n" + b"TTTT\n" * 5000
    encoded = bytearray(gzip.compress(payload))
    encoded[-8] ^= 1  # Same length and valid deflate stream, invalid CRC32.
    prefix = gzip.compress(b">earlier\n" + b"ACGT\n" * 5000) if member == "later" else b""
    path = tmp_path / "source.fa.gz"
    path.write_bytes(prefix + encoded)
    real_open, handles = common.open_text, []

    def tracked_open(source, *args, **kwargs):
        handle = real_open(source, *args, **kwargs)
        handles.append(handle)
        return handle

    monkeypatch.setattr(common, "open_text", tracked_open)
    expected = observe(lambda: reference_count(path))
    assert expected[1] is not None and expected[1][0] is gzip.BadGzipFile
    assert observe(lambda: common.count_fasta_records(path)) == expected
    target = "earlier" if member == "later" and match == "first" else match
    gff = tmp_path / "annotation.gff"
    gff.write_text(f"{target}\tsrc\tgene\t1\t2\t.\t+\t.\tID=g\n")
    assert observe(lambda: common.fasta_headers_overlap_gff_seqids(path, gff)) == observe(
        lambda: reference_overlap(path, {target}))
    assert all(handle.closed for handle in handles)


def test_declared_compression_and_archive_fallback_remain_strict(tmp_path):
    bad = tmp_path / "source.fa.gz"
    bad.write_bytes(b">first\nAC\n")
    assert observe(lambda: common.count_fasta_records(bad)) == observe(lambda: reference_count(bad))
    path = tmp_path / "source.fa.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        for name in ("first.data", "second.data"):
            payload = b">not-selected\nAC\n"
            member = tarfile.TarInfo(name)
            member.size = len(payload)
            archive.addfile(member, io.BytesIO(payload))
    assert common.count_fasta_records(path) == (0, "")


def test_archive_member_order_and_duplicate_headers_are_retained(tmp_path):
    path = tmp_path / "source.fa.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        for name, payload in (("z.fa", b">z-first\nAC\n>duplicate\nT\n"),
                              ("a.fasta", b">a-second\nGT\n>duplicate\nC\n")):
            member = tarfile.TarInfo(name)
            member.size = len(payload)
            archive.addfile(member, io.BytesIO(payload))
    assert common.count_fasta_records(path) == (4, "z-first")
    assert list(common.iter_fasta_records(path)) == [
        ("z-first", "AC"), ("duplicate", "T"), ("a-second", "GT"), ("duplicate", "C")]


def test_overlap_yields_only_at_next_header_and_closes_reader(tmp_path, monkeypatch):
    fasta = tmp_path / "source.fa"
    gff = tmp_path / "annotation.gff"
    gff.write_text("match\tsrc\tgene\t1\t2\t.\t+\t.\tID=g\n")
    real_open = common.open_text
    observed = []

    class TrackedReader(io.StringIO):
        def __next__(self):
            line = super().__next__()
            observed.append(line)
            return line

    handles = []

    def tracked_open(path, *args, **kwargs):
        if Path(path) != fasta:
            return real_open(path, *args, **kwargs)
        handle = TrackedReader("ignored preamble\n>match\nAA\nCC\n>later\nnot consumed\n")
        handles.append(handle)
        return handle

    monkeypatch.setattr(common, "open_text", tracked_open)
    assert common.fasta_headers_overlap_gff_seqids(fasta, gff)
    assert observed == ["ignored preamble\n", ">match\n", "AA\n", "CC\n", ">later\n"]
    assert all(handle.closed for handle in handles)


@pytest.mark.parametrize("transport", ["raw", "gz", "tar-gz"])
def test_actual_direct_provider_discovery_preserves_generic_genome_classification(tmp_path, transport):
    species = tmp_path / "Species_one"
    species.mkdir()
    path = write_fasta(species, b">lcl|chr1 desc\nACGT\n>chr2\nTT\n", transport)
    gff = species / "annotation.gff"
    gff.write_text("chr1\tsrc\tgene\t1\t4\t.\t+\t.\tID=g1\n")
    tasks, warnings, errors = discover_generic_species_dir_tasks("direct", tmp_path)
    assert not warnings and not errors
    assert len(tasks) == 1
    assert tasks[0]["genome_path"] == path
    assert tasks[0]["cds_path"] is None


@pytest.mark.parametrize("species", ["Cucurbita_pepo", "Euphorbia_peplus", "Peperomia_pellucida"])
def test_species_pep_substring_does_not_turn_genome_into_cds(tmp_path, species):
    directory = tmp_path / species
    directory.mkdir()
    genome = directory / (species + ".genome.fa.gz")
    with gzip.open(genome, "wt") as handle:
        handle.write(">chr1\nATGAAATAA\n")
    (directory / (species + ".gff3")).write_text(
        "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=g1\n")
    tasks, warnings, errors = discover_generic_species_dir_tasks("local", tmp_path)
    assert not warnings and not errors
    assert len(tasks) == 1
    assert tasks[0]["genome_path"] == genome
    assert tasks[0]["cds_path"] is None


@pytest.mark.parametrize("name", ["Plant_one.pep.fa", "Plant_one.protein.fa",
    "Plant_one.cds_from_genomic.fna.gz", "Plant_one.mrna.genome.fa", "Plant_one.rna_from_genomic.fna"])
def test_actual_nongenome_role_tokens_still_exclude_genome(name):
    assert not is_probable_genome_filename("local", name)

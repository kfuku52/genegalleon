"""Compression follows public-file bytes when the filename omits its suffix."""
import bz2
import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from format_species_annotation.common import iter_fasta_records
from format_species_annotation.reference import genome_reference_index
from format_species_common import gff_contains_cds_feature
from format_species_writers import inspect_invalid_utf8, open_binary, open_text


@pytest.mark.parametrize("compress", [gzip.compress, bz2.compress])
@pytest.mark.parametrize("filename", ["genome.fa", "49345861"])
def test_compressed_inputs_without_suffix_reach_all_formatter_readers(tmp_path, compress, filename):
    sequence = b">chr1 description\nATGAAATAA\n"
    source = tmp_path / filename
    source.write_bytes(compress(sequence))
    assert list(iter_fasta_records(source)) == [("chr1 description", "ATGAAATAA")]
    assert genome_reference_index(source)["chr1"] == 9
    with open_binary(source) as stream:
        assert stream.read() == sequence
    assert inspect_invalid_utf8(source)["invalid_utf8_bytes"] == 0
    source.write_bytes(compress(b"chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=c1\n"))
    assert gff_contains_cds_feature(source)


@pytest.mark.parametrize("filename", ["genome.fa", "genome.fa.gz", "genome.fa.GZ"])
@pytest.mark.parametrize("damage", ["truncated", "crc"])
def test_corrupt_gzip_is_rejected_with_or_without_suffix(tmp_path, filename, damage):
    source = tmp_path / filename
    payload = gzip.compress(b">chr1\nATGAAATAA\n")
    if damage == "truncated":
        payload = payload[:-4]
    else:
        payload = payload[:-8] + bytes([payload[-8] ^ 1]) + payload[-7:]
    source.write_bytes(payload)
    with pytest.raises((OSError, EOFError)):
        list(iter_fasta_records(source))
    with pytest.raises((OSError, EOFError)):
        inspect_invalid_utf8(source)


@pytest.mark.parametrize("filename", ["genome.fa.gz", "genome.fa.bz2"])
def test_declared_compression_cannot_silently_accept_plain_text(tmp_path, filename):
    source = tmp_path / filename
    source.write_text(">chr1\nATGAAATAA\n")
    with pytest.raises(OSError):
        list(iter_fasta_records(source))


@pytest.mark.parametrize("filename", ["output.fa", "output.fa.gz", "output.fa.bz2"])
def test_output_compression_still_follows_destination_suffix(tmp_path, filename):
    target = tmp_path / filename
    target.write_bytes(gzip.compress(b"old compressed input"))
    sequence = ">chr1\nATGAAATAA\n"
    with open_text(target, "wt") as stream:
        stream.write(sequence)
    payload = target.read_bytes()
    if filename.endswith(".gz"):
        payload = gzip.decompress(payload)
    elif filename.endswith(".bz2"):
        payload = bz2.decompress(payload)
    assert payload == sequence.encode()

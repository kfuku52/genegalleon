import gzip
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
SCRIPT_PATH = SUPPORT_DIR / "gff_attribute_identifier.py"


def load_module():
    sys.path.insert(0, str(SUPPORT_DIR))
    try:
        spec = spec_from_file_location("gff_attribute_identifier", SCRIPT_PATH)
        module = module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        try:
            sys.path.remove(str(SUPPORT_DIR))
        except ValueError:
            pass


def write_fasta(tmp_path, payload, compressed):
    path = tmp_path / ("input.fa.gz" if compressed else "input.fa")
    path.write_bytes(gzip.compress(payload, mtime=0) if compressed else payload)
    return path


@pytest.mark.parametrize("compressed", [False, True], ids=["plain", "gzip"])
@pytest.mark.parametrize(
    ("payload", "expected"),
    [
        pytest.param(b"", [], id="empty"),
        pytest.param(
            b">gene1 description with spaces\nACGT\nACGT\n>gene1\nATG\n>gene2",
            ["gene1 description with spaces", "gene1", "gene2"],
            id="description-duplicates-and-final-header",
        ),
        pytest.param(
            ">gene_α description 日本語\nACGT\n".encode("utf-8"),
            ["gene_α description 日本語"],
            id="unicode",
        ),
        pytest.param(
            b"ACGT\n>gene1\nAA>BB\n",
            ["gene1", "AA>BB"],
            id="existing-header-substring-rule",
        ),
    ],
)
def test_read_fasta_seqname_preserves_headers(tmp_path, compressed, payload, expected):
    mod = load_module()
    path = write_fasta(tmp_path, payload, compressed)

    result = mod.read_fasta_seqname(str(path))

    assert result.tolist() == expected
    assert result.index.tolist() == list(range(len(expected)))


@pytest.mark.parametrize("compressed", [False, True], ids=["plain", "gzip"])
def test_read_fasta_seqname_preserves_crlf_handling(tmp_path, compressed):
    mod = load_module()
    path = write_fasta(tmp_path, b">gene1\r\nATG\r\n>gene2\r\n", compressed)

    # Plain text reads normalize CRLF; the existing binary gzip reader strips LF only.
    expected = ["gene1\r", "gene2\r"] if compressed else ["gene1", "gene2"]
    assert mod.read_fasta_seqname(str(path)).tolist() == expected


@pytest.mark.parametrize("compressed", [False, True], ids=["plain", "gzip"])
def test_read_fasta_seqname_rejects_invalid_utf8_in_sequence_lines(tmp_path, compressed):
    mod = load_module()
    path = write_fasta(tmp_path, b">gene1\nATG\n\xff\n", compressed)

    with pytest.raises(UnicodeDecodeError):
        mod.read_fasta_seqname(str(path))


def test_read_fasta_seqname_rejects_truncated_gzip(tmp_path):
    mod = load_module()
    path = tmp_path / "input.fa.gz"
    path.write_bytes(gzip.compress(b">gene1\nATG\n", mtime=0)[:-5])

    with pytest.raises(EOFError, match="end-of-stream"):
        mod.read_fasta_seqname(str(path))


def test_read_fasta_seqname_reads_all_gzip_members(tmp_path):
    mod = load_module()
    path = tmp_path / "input.fa.gz"
    path.write_bytes(
        gzip.compress(b">gene1\nATG\n", mtime=0)
        + gzip.compress(b">gene2 description\nAAA\n", mtime=0)
    )

    assert mod.read_fasta_seqname(str(path)).tolist() == ["gene1", "gene2 description"]


@pytest.mark.parametrize(
    "species",
    [
        "Oryza_sativa_subsp_japonica",
        "Cenchrus_americanus_x_Cenchrus_purpureus",
    ],
)
def test_identify_attribute_key_strips_species_prefix(tmp_path, species):
    mod = load_module()
    fasta = tmp_path / f"{species}.cds.fa"
    gff = tmp_path / f"{species}.gff3"
    fasta.write_text(f">{species}_gene1\nATG\n", encoding="utf-8")
    gff.write_text("chr1\tsrc\tmRNA\t1\t3\t.\t+\t.\tID=gene1\n", encoding="utf-8")

    args = type("Args", (), {"fasta_file": str(fasta), "gff_file": str(gff), "gff_feature": "mRNA"})()

    assert mod.identify_attribute_key(args) == "ID"

import gzip
import io

import pandas as pd
import pytest

from workflow.support.line_comment_stream import LineCommentStream


@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize("size", [1, 7, 8192])
def test_comment_stream_preserves_literal_hashes_and_quoted_continuations(tmp_path, compressed, size):
    payload = ('# initial comment\ngene_id\tdescription\n'
               'gene#1\t"first\n# still data\nlast with ""quote"""\n'
               '# comment between rows\n'
               'gene2\tunquoted"quote\n# final comment\n')
    expected = ('\ngene_id\tdescription\ngene#1\t"first\n# still data\nlast with ""quote"""\n'
                '\ngene2\tunquoted"quote\n\n')
    path = tmp_path / ("table.tsv.gz" if compressed else "table.tsv")
    path.write_bytes(gzip.compress(payload.encode(), mtime=0) if compressed else payload.encode())
    with LineCommentStream(path, quoted=True) as stream:
        assert stream.read(0) == ""
        parts = []
        while chunk := stream.read(size):
            parts.append(chunk)
        assert "".join(parts) == expected
        stream.seek(0)
        assert "".join(stream) == expected
        stream.seek(0)
        actual = pd.read_csv(stream, sep="\t")
    assert actual.gene_id.tolist() == ["gene#1", "gene2"]
    assert actual.description.tolist() == ['first\n# still data\nlast with "quote"', 'unquoted"quote']


def test_borrowed_stream_remains_open_and_replayable():
    source = io.StringIO("# comment\nx\tID=gene#1\n")
    with LineCommentStream(source) as stream:
        assert stream.readline() == "\n"
        assert stream.read(2) == "x\t"
        assert stream.readline() == "ID=gene#1\n"
        stream.seek(0)
        assert stream.read() == "\nx\tID=gene#1\n"
    assert not source.closed


def test_gff_fasta_tail_is_consumed_for_integrity(tmp_path):
    path = tmp_path / "annotations.gff.gz"
    payload = b"# comment\nx\ts\tgene\t1\t3\t.\t+\t.\tID=gene#1\n##FASTA\n>x\nATG\n"
    path.write_bytes(gzip.compress(payload, mtime=0)[:-5])
    with LineCommentStream(path, gff=True) as stream, pytest.raises(EOFError):
        stream.read()


def test_comment_filter_preserves_physical_error_row_numbers():
    payload = "# header comment\nid\tvalue\ng1\t1\n# middle comment\ng2\t2\textra\n"
    with pytest.raises(pd.errors.ParserError) as reference:
        pd.read_csv(io.StringIO(payload), sep="\t", comment="#")
    with LineCommentStream(io.StringIO(payload)) as stream, pytest.raises(pd.errors.ParserError) as actual:
        pd.read_csv(stream, sep="\t")
    assert str(actual.value) == str(reference.value)

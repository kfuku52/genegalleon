"""Coding GFF parsing preserves row order without an expanded temporary file."""
import gzip
import io
import sys
from pathlib import Path
from unittest import mock

import pandas
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
from coding_gff_stream import CodingGffStream
from gff2genestat import _parse_gff_attributes, read_gff_table
from input_generation_cohort import qc_directory


@pytest.mark.parametrize('compressed', [False, True])
def test_stream_matches_historical_spool_including_ancestor_order(tmp_path, compressed):
    rows = ['# annotation\n',
            'chr1\ts\tmatch\t1\t9\t.\t+\t.\tID=outer\n',
            'chr1\ts\tmatch_part\t1\t9\t.\t+\t.\tID=inner;Parent=outer\n',
            'chr1\ts\tmatch\t1\t9\t.\t+\t.\tID=unrelated\n',
            'chr1\ts\tCDS\t1\t9\t.\t+\t0\tID=c;Parent=inner\n',
            '##FASTA\n', '>chr1\n', 'ATGAAATAA\n']
    expected = rows[0] + rows[4] + rows[2] + rows[1]
    path = tmp_path / ('input.gff.gz' if compressed else 'input.gff')
    if compressed:
        with gzip.open(path, 'wt') as out:
            out.writelines(rows)
    else:
        path.write_text(''.join(rows))
    with CodingGffStream(path, _parse_gff_attributes) as stream:
        assert stream.read(0) == ''
        fragments = list(iter(lambda: stream.read(7), ''))
        assert ''.join(fragments) == expected
        stream.seek(0)
        actual = read_gff_table(stream)
    baseline = tmp_path / 'historical.gff'
    baseline.write_text(expected)
    pandas.testing.assert_frame_equal(actual, read_gff_table(str(baseline)))


def test_parser_fallback_rewinds_without_spooling(tmp_path):
    path = tmp_path / 'input.gff'
    path.write_text('c\ts\tCDS\t1\t9\t.\t+\t0\tID=c\n')
    with CodingGffStream(path, _parse_gff_attributes) as stream:
        stream.read(5)
        with mock.patch.object(pandas, 'read_csv', side_effect=pandas.errors.ParserError('fixture')):
            result = read_gff_table(stream)
    assert result.shape == (1, 9)
    assert result.iloc[0, 8] == 'ID=c'


def test_qc_attempt_directory_is_bounded_and_preserves_default(tmp_path):
    assert qc_directory(tmp_path) == tmp_path / 'qc'
    assert qc_directory(tmp_path, tmp_path / 'qc-attempt-2') == tmp_path / 'qc-attempt-2'
    with pytest.raises(ValueError, match='inside'):
        qc_directory(tmp_path, tmp_path.parent / 'foreign')
    (tmp_path / 'redirect').symlink_to(tmp_path / 'qc')
    with pytest.raises(ValueError, match='symlink'):
        qc_directory(tmp_path, tmp_path / 'redirect')


def test_only_rewind_is_supported(tmp_path):
    path = tmp_path / 'empty.gff'
    path.write_text('')
    with CodingGffStream(path, _parse_gff_attributes) as stream:
        with pytest.raises(io.UnsupportedOperation):
            stream.seek(1)

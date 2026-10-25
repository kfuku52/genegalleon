"""Execute the real nested promoter annotator without running unrelated stages."""

import ast
import gzip
import re
from pathlib import Path

import pandas
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "support" / "orthogroup_statistics.py"


@pytest.fixture(scope="module")
def get_n_coordinates():
    tree = ast.parse(SCRIPT.read_text(), filename=str(SCRIPT))
    main = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "main")
    functions = [node for node in main.body
                 if isinstance(node, ast.FunctionDef) and node.name in {"open_text", "get_N_coordinate"}]
    namespace = {"gzip": gzip, "pandas": pandas, "re": re}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(SCRIPT), "exec"), namespace)
    return namespace["get_N_coordinate"]


def write_fasta(tmp_path, text, compressed):
    path = tmp_path / ("promoters.fa.gz" if compressed else "promoters.fa")
    if compressed:
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)
    return str(path)


def legacy_n_coordinates(file):
    """Frozen full-input implementation for literal-record compatibility."""
    opener = gzip.open if file.endswith(".gz") else open
    with opener(file, "rt") as handle:
        seqs = handle.read().split(">")
    if seqs[0] == "":
        seqs = seqs[1::]
    if len(seqs) == 0 or all(s.strip() == "" for s in seqs):
        return pandas.DataFrame({"node_name": [], "promoter_N": [], "promoter_available": []})
    seqnames = []
    for i in range(len(seqs)):
        seqnames.append(re.sub("\n.*", "", seqs[i].replace(">", "")))
        seqs[i] = re.sub("^[^\n]*\n", "", seqs[i])
        seqs[i] = re.sub("\n", "", seqs[i])
    maximum = max(len(seq) for seq in seqs)
    out = []
    for seq in seqs:
        slices = []
        if "N" in seq or "n" in seq:
            for match in re.finditer("[Nn]+", seq):
                start, end = match.span()
                slices.append(str(end) if end == start + 1 else f"{start + 1}:{end}")
        if len(seq) < maximum:
            slices.append(str(len(seq)) + ":" + str(maximum))
        out.append(re.sub("^,", "", ",".join(slices)))
    result = pandas.DataFrame({"node_name": seqnames, "promoter_N": out})
    result.loc[:, "promoter_available"] = "Y"
    return result


@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize(("sequence", "coordinates"), [
    ("N", "1"),
    ("n", "1"),
    ("nN", "1:2"),
    ("ACNNnN", "3:6"),
    ("NAnNCn", "1,3:4,6"),
    ("NNACNnN", "1:2,5:7"),
    ("ANAN", "2,4"),
    ("ACGT", ""),
    ("", ""),
    ("Nn" * 5000, "1:10000"),
])
def test_promoter_coordinates_are_one_based_and_keep_terminal_mixed_case_runs(
    tmp_path, get_n_coordinates, compressed, sequence, coordinates
):
    path = write_fasta(tmp_path, f">node\n{sequence}\n", compressed)
    assert get_n_coordinates(path).to_dict("list") == {
        "node_name": ["node"], "promoter_N": [coordinates], "promoter_available": ["Y"],
    }


@pytest.mark.parametrize("compressed", [False, True])
def test_short_promoters_preserve_existing_outside_scaffold_coordinates(tmp_path, get_n_coordinates, compressed):
    path = write_fasta(tmp_path, ">terminal\nACNNnN\n>isolated\nNAnCNN\n>short\nN\n>noNs\nACG\n>empty\n", compressed)
    assert get_n_coordinates(path).to_dict("list") == {
        "node_name": ["terminal", "isolated", "short", "noNs", "empty"],
        "promoter_N": ["3:6", "1,3,5:6", "1,1:6", "3:6", "0:6"],
        "promoter_available": ["Y"] * 5,
    }


@pytest.mark.parametrize("compressed", [False, True])
def test_wrapped_fasta_keeps_literal_headers_and_combines_runs_across_lines(tmp_path, get_n_coordinates, compressed):
    path = write_fasta(tmp_path, ">001 description\r\nAN\r\nnN\r\n>NA\r\nAC\r\nGT\r\n", compressed)
    assert get_n_coordinates(path).to_dict("list") == {
        "node_name": ["001 description", "NA"],
        "promoter_N": ["2:4", ""],
        "promoter_available": ["Y", "Y"],
    }


@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize("text", ["", "\n", " \n\t\n", ">"])
def test_empty_promoter_input_retains_empty_annotation_schema(tmp_path, get_n_coordinates, compressed, text):
    path = write_fasta(tmp_path, text, compressed)
    result = get_n_coordinates(path)
    assert result.empty
    assert result.columns.tolist() == ["node_name", "promoter_N", "promoter_available"]


@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize("text", [
    "", "\n", " \n\t\n", ">", ">>>\n", "> \n\t>\n ",
    "ACNN", "prefix\nNn\n>node\nAN", ">node", ">node\nAN\n>tail",
    ">>N", ">node\nAN>>tail\nN\n>", ">node\nAN>inline\nN\n",
    " >node\nN\n", ">\nNn\n> \nAN\n", ">001 description\r\nN\r\nn\r\n>NA\r\nAC\r\n",
])
def test_literal_promoter_records_preserve_complete_legacy_frame(tmp_path, get_n_coordinates, compressed, text):
    path = write_fasta(tmp_path, text, compressed)
    pandas.testing.assert_frame_equal(get_n_coordinates(path), legacy_n_coordinates(path))


@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize("offset", [-1, 0, 1])
def test_literal_record_delimiters_at_read_boundaries(tmp_path, get_n_coordinates, compressed, offset):
    boundary = 1024 * 1024 + offset
    prefix = ">first\n"
    text = prefix + "A" * (boundary - len(prefix) - 2) + "Nn>second\nnN\n>>last"
    path = write_fasta(tmp_path, text, compressed)
    pandas.testing.assert_frame_equal(get_n_coordinates(path), legacy_n_coordinates(path))


@pytest.mark.parametrize("compressed", [False, True])
def test_promoter_record_and_n_run_cross_multiple_reads(tmp_path, get_n_coordinates, compressed):
    prefix = ">long\n"
    text = prefix + "A" * (1024 * 1024 - len(prefix) - 1) + "N\nnN\n" + "A" * (2 * 1024 * 1024)
    text += ">short\nNN\n"
    path = write_fasta(tmp_path, text, compressed)
    pandas.testing.assert_frame_equal(get_n_coordinates(path), legacy_n_coordinates(path))


@pytest.mark.parametrize("compressed", [False, True])
def test_long_n_record_followed_by_n_free_record(tmp_path, get_n_coordinates, compressed):
    text = ">large_n\n" + "Nn" * (1024 * 1024 // 2 + 17)
    text += "\n>without_n\n" + "A" * (1024 * 1024 + 37) + "\n>last\nN\n"
    path = write_fasta(tmp_path, text, compressed)
    pandas.testing.assert_frame_equal(get_n_coordinates(path), legacy_n_coordinates(path))


def test_promoter_reads_are_bounded(tmp_path, get_n_coordinates, monkeypatch):
    path = write_fasta(tmp_path, ">long\n" + "A" * (2 * 1024 * 1024) + "N\n", False)
    original_open = get_n_coordinates.__globals__["open_text"]
    sizes = []

    class BoundedReader:
        def __enter__(self):
            self.handle = original_open(path)
            return self

        def __exit__(self, *exc):
            self.handle.close()

        def read(self, size=-1):
            assert 0 < size <= 1024 * 1024
            sizes.append(size)
            return self.handle.read(size)

    monkeypatch.setitem(get_n_coordinates.__globals__, "open_text", lambda _: BoundedReader())
    pandas.testing.assert_frame_equal(get_n_coordinates(path), legacy_n_coordinates(path))
    assert len(sizes) >= 3


@pytest.mark.parametrize("damage", ["truncated", "crc"])
def test_gzip_tail_failure_is_not_hidden_after_complete_records(tmp_path, get_n_coordinates, damage):
    first = gzip.compress((">early\nNN\n>filler\n" + "A" * (1024 * 1024) + "\n").encode())
    final = gzip.compress(b">late\nACN\n")
    if damage == "truncated":
        final = final[:-5]
        error = EOFError
    else:
        final = final[:-8] + bytes([final[-8] ^ 1]) + final[-7:]
        error = gzip.BadGzipFile
    path = tmp_path / "promoters.fa.gz"
    path.write_bytes(first + final)
    with pytest.raises(error):
        legacy_n_coordinates(str(path))
    with pytest.raises(error):
        get_n_coordinates(str(path))


def test_all_concatenated_gzip_members_are_processed(tmp_path, get_n_coordinates):
    path = tmp_path / "promoters.fa.gz"
    path.write_bytes(gzip.compress(b">first\nAN\n") + gzip.compress(b"nN\n>last"))
    pandas.testing.assert_frame_equal(get_n_coordinates(str(path)), legacy_n_coordinates(str(path)))

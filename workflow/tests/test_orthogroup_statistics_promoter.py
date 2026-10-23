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

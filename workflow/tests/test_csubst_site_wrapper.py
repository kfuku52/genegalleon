import ast
import os
import re
from pathlib import Path

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "script" / "csubst_site_wrapper.py"


def load_extract_pdb_id():
    source = SCRIPT_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(SCRIPT_PATH))
    extract_node = None
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == "extract_pdb_id":
            extract_node = node
            break
    if extract_node is None:
        raise AssertionError("extract_pdb_id function not found in csubst_site_wrapper.py")
    module_tree = ast.Module(body=[extract_node], type_ignores=[])
    code = compile(module_tree, filename=str(SCRIPT_PATH), mode="exec")
    namespace = {"os": os, "re": re}
    exec(code, namespace)
    return namespace["extract_pdb_id"]


def test_extract_pdb_id_returns_none_for_missing_directory(tmp_path):
    extract_pdb_id = load_extract_pdb_id()
    missing = tmp_path / "missing"
    assert extract_pdb_id(str(missing)) is None


def test_extract_pdb_id_returns_none_when_no_match(tmp_path):
    extract_pdb_id = load_extract_pdb_id()
    (tmp_path / "README.txt").write_text("x", encoding="utf-8")
    (tmp_path / "csubst_site.tsv").write_text("x", encoding="utf-8")
    assert extract_pdb_id(str(tmp_path)) is None


def test_extract_pdb_id_parses_expected_identifier(tmp_path):
    extract_pdb_id = load_extract_pdb_id()
    (tmp_path / "csubst_site.2XYZ.fa").write_text(">x\nAA\n", encoding="utf-8")
    assert extract_pdb_id(str(tmp_path)) == "2XYZ"

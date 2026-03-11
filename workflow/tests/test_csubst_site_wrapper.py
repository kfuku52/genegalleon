import ast
from importlib.util import module_from_spec, spec_from_file_location
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "csubst_site_wrapper.py"


def load_module():
    spec = spec_from_file_location("csubst_site_wrapper", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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


def test_extract_pdb_id_ignores_hidden_files_and_is_deterministic(tmp_path):
    extract_pdb_id = load_extract_pdb_id()
    (tmp_path / ".csubst_site.ZZZZ.fa").write_text(">x\nAA\n", encoding="utf-8")
    (tmp_path / "csubst_site.BBBB.fa").write_text(">x\nAA\n", encoding="utf-8")
    (tmp_path / "csubst_site.AAAA.fa").write_text(">x\nAA\n", encoding="utf-8")
    assert extract_pdb_id(str(tmp_path)) == "AAAA"


def test_get_cb_required_columns_keeps_needed_columns_only():
    mod = load_module()
    cols = [
        "orthogroup",
        "OCNany2spe",
        "ECNany2spe",
        "OCSany2spe",
        "ECSany2spe",
        "omegaCany2spe",
        "OCNCoD",
        "branch_id_1",
        "branch_id_2",
        "is_fg_traitA",
        "branch_num_fg_stem_traitA",
        "unused_col",
    ]

    out = mod.get_cb_required_columns(cols, ["traitA"])

    assert out == [
        "orthogroup",
        "OCNany2spe",
        "ECNany2spe",
        "OCSany2spe",
        "ECSany2spe",
        "omegaCany2spe",
        "OCNCoD",
        "branch_id_1",
        "branch_id_2",
        "is_fg_traitA",
        "branch_num_fg_stem_traitA",
    ]


def test_skip_lower_order_filters_subset_rows_per_orthogroup():
    mod = load_module()
    cb_passed = pandas.DataFrame(
        [
            {"orthogroup": "OG1", "branch_id_1": 1, "branch_id_2": 2},
            {"orthogroup": "OG1", "branch_id_1": 2, "branch_id_2": 3},
            {"orthogroup": "OG2", "branch_id_1": 5, "branch_id_2": 6},
        ]
    )
    already = {
        "traitA": {
            "OG1": [frozenset([1, 2, 9])],
            "OG2": [frozenset([7, 8])],
        }
    }

    out, updated = mod.skip_lower_order(cb_passed, 2, "traitA", already)

    assert out["orthogroup"].tolist() == ["OG1", "OG2"]
    assert out.loc[0, ["branch_id_1", "branch_id_2"]].tolist() == [2, 3]
    assert out.loc[1, ["branch_id_1", "branch_id_2"]].tolist() == [5, 6]
    assert frozenset([2, 3]) in updated["traitA"]["OG1"]
    assert frozenset([5, 6]) in updated["traitA"]["OG2"]


def test_load_annotation_besthits_reads_besthit_columns_only(tmp_path):
    mod = load_module()
    annot_dir = tmp_path / "Orthogroups"
    annot_dir.mkdir(parents=True)
    infile = annot_dir / "Orthogroups.GeneCount.annotated.tsv"
    infile.write_text(
        "Orthogroup\tTotal\tbesthit_0.5\tbesthit_0.95\n"
        "OG1\t4\thitA\thitB\n",
        encoding="utf-8",
    )

    out = mod.load_annotation_besthits(str(tmp_path))

    assert out.columns.tolist() == ["orthogroup", "besthit_0.5", "besthit_0.95"]
    assert out.loc[0, "orthogroup"] == "OG1"


def test_help_has_no_side_effect_files(tmp_path):
    proc = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "--help"],
        cwd=str(tmp_path),
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0
    assert not (tmp_path / "generate_orthogroup_database.log").exists()
    assert not (tmp_path / "mpl").exists()

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "script" / "extract_astral_support_label.py"


def load_module():
    spec = spec_from_file_location("extract_astral_support_label", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_extract_support_labels_handles_reordered_fields_and_extra_keys():
    mod = load_module()
    tree = "(A,(B,C)'[q1=0.75;q2=0.10;q3=0.15;f1=3.0;pp1=0.82;QC=2;EN=3.0]':0.1);"
    converted, missing = mod.extract_support_labels(tree, "q1")
    assert missing == 0
    assert "[q1=" not in converted
    assert "QC=2" not in converted
    assert "0.75" in converted


def test_extract_support_labels_handles_legacy_order():
    mod = load_module()
    tree = "(A,(B,C)'[pp1=0.8;pp2=0.1;pp3=0.1;f1=3;f2=0;f3=0;q1=1.0;q2=0.0;q3=0.0]':0.2);"
    converted, missing = mod.extract_support_labels(tree, "pp1")
    assert missing == 0
    assert "[pp1=" not in converted
    assert "0.8" in converted

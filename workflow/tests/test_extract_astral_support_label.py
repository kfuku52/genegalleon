from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "extract_astral_support_label.py"


def load_module():
    spec = spec_from_file_location("extract_astral_support_label", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize(
    ("tree", "support_key", "expected"),
    [
        (
            "(A,(B,C)'[q1=0.75;q2=0.10;q3=0.15;f1=3.0;pp1=0.82;QC=2;EN=3.0]':0.1);",
            "q1",
            "0.75",
        ),
        (
            "(A,(B,C)'[pp1=0.8;pp2=0.1;pp3=0.1;f1=3;f2=0;f3=0;q1=1.0;q2=0.0;q3=0.0]':0.2);",
            "pp1",
            "0.8",
        ),
    ],
)
def test_extract_support_labels_handles_field_order(tree, support_key, expected):
    mod = load_module()
    converted, missing = mod.extract_support_labels(tree, support_key)
    assert missing == 0
    assert f"[{support_key}=" not in converted
    assert expected in converted

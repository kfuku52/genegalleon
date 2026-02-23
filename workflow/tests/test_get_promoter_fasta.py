from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "get_promoter_fasta.py"


def load_module():
    spec = spec_from_file_location("get_promoter_fasta", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_normalize_ncpu_clamps_non_positive_values():
    mod = load_module()
    assert mod.normalize_ncpu(0) == 1
    assert mod.normalize_ncpu(-5) == 1
    assert mod.normalize_ncpu(3) == 3


def test_normalize_ncpu_rejects_invalid_values():
    mod = load_module()
    with pytest.raises(ValueError):
        mod.normalize_ncpu("abc")


def test_get_genome_file_ignores_hidden_entries_and_directories(tmp_path: Path):
    mod = load_module()
    (tmp_path / ".Species_a.fa").write_text(">a\nAT\n", encoding="utf-8")
    (tmp_path / "Species_a.extra").mkdir()
    target = tmp_path / "Species_a.fa"
    target.write_text(">a\nAT\n", encoding="utf-8")

    detected = mod.get_genome_file(str(tmp_path), "Species_a")
    assert detected == str(target)

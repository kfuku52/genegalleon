from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
import types


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "iqtree2mapnh.py"


def load_module(monkeypatch):
    package = types.ModuleType("kftools")
    kfphylo = types.ModuleType("kftools.kfphylo")
    kfseq = types.ModuleType("kftools.kfseq")
    monkeypatch.setitem(sys.modules, "kftools", package)
    monkeypatch.setitem(sys.modules, "kftools.kfphylo", kfphylo)
    monkeypatch.setitem(sys.modules, "kftools.kfseq", kfseq)
    spec = spec_from_file_location("iqtree2mapnh_module", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_alignment_subset_nuc_freqs_does_not_use_boundaryless_leaf_prefix(tmp_path: Path, monkeypatch):
    mod = load_module(monkeypatch)
    alignment = tmp_path / "aln.fa"
    alignment.write_text(
        ">Species_a\nAAA\n"
        ">Species_a_subsp_x\nCCC\n",
        encoding="utf-8",
    )

    nuc_freqs = mod.alignment_subset_nuc_freqs(str(alignment), "F3X4", leaf_names=["Species_a"])

    assert all(pos_freq["A"] == 1.0 for pos_freq in nuc_freqs)
    assert all(pos_freq["C"] == 0.0 for pos_freq in nuc_freqs)

import sys
import types
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

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


@pytest.mark.parametrize('repeats', [1, 7, 122, 123, 257])
def test_alignment_frequencies_weight_codon_counts_and_ignore_whole_ambiguous_codons(tmp_path, monkeypatch, repeats):
    mod = load_module(monkeypatch)
    path = tmp_path / 'mixed.fa'
    path.write_text('>leaf\n' + 'acg' * repeats + 'uuu' * 2 + 'RTA' + 'A--' + 'éAT' + 'AA\n')
    result = mod.alignment_subset_nuc_freqs(path, 'F3X4+G4')
    expected = []
    for base in 'ACG':
        row = dict.fromkeys('ACGT', 0.0)
        row[base] = repeats / (repeats + 2)
        row['T'] = 2 / (repeats + 2)
        expected.append(row)
    assert result == expected


@pytest.mark.parametrize('sequence', ['NNNNNNNN', 'A--éATAA', '', '??A--T'])
def test_alignment_without_usable_codons_keeps_uniform_frequencies(tmp_path, monkeypatch, sequence):
    mod = load_module(monkeypatch)
    path = tmp_path / 'unknown.fa'
    path.write_text(f'>leaf\n{sequence}\n')
    assert mod.alignment_subset_nuc_freqs(path, 'F3X4') == [dict.fromkeys('ACGT', 0.25)] * 3


def test_alignment_subset_keeps_normalized_collisions_and_last_duplicate_record(tmp_path, monkeypatch):
    mod = load_module(monkeypatch)
    path = tmp_path / 'names.fa'
    path.write_text('>leaf-1\nAAA\n>leaf_1\nGGG\n>leaf-1\nCCC\nCCC\n>unselected\nTTT\n')
    result = mod.alignment_subset_nuc_freqs(path, 'F3X4', ['leaf_1'])
    assert result == [{'A': 0.0, 'C': 2 / 3, 'G': 1 / 3, 'T': 0.0}] * 3
    with pytest.raises(ValueError, match='No matching sequences'):
        mod.alignment_subset_nuc_freqs(path, 'F3X4', ['absent'])
    with pytest.raises(ValueError, match='supports only F3X4'):
        mod.alignment_subset_nuc_freqs(tmp_path / 'absent.fa', 'F1X4')

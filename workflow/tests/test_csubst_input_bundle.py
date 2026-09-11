import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
import csubst_input_bundle as bundle
import csubst_site_wrapper as sites


def make_bundle(root):
    directory = root / bundle.STRUCTURAL_DIR
    directory.mkdir(parents=True)
    for suffix in bundle.FIT_SUFFIXES:
        (directory / f"csubst.{suffix}").write_text("fit\n")
    (directory / "csubst.iqtree").write_text("Model of substitution: GY+FQ\n")
    (directory / "csubst.fasta").write_text(">a\nATG---NNN\n>b\nATGTTA---\n")
    (directory / "csubst_3di_state_cache.npz").write_bytes(b"cache")
    (directory / "inspect").mkdir()
    (directory / "inspect/csubst_alignment_3di.fa").write_text(">a\nA-X\n>b\nCD-\n>Node1\nACD\n")
    bundle.finalize(root, 2)
    return directory


def test_full_cds_filters_tips_but_retains_gap_and_unknown_columns(tmp_path):
    source = tmp_path / "source.fa"
    source.write_text(">extra\nATG---NNN\n>a\nATG---NNN\n>b\nATGTTA---\n")
    tree = tmp_path / "tree.nwk"
    tree.write_text("(b:1,a:1);")
    output = tmp_path / "full.fa"
    bundle.prepare_full_alignment(source, tree, output)
    assert [(r.id, str(r.seq)) for r in bundle.read_alignment(output)] == [("b", "ATGTTA---"), ("a", "ATG---NNN")]
    source.write_text(">a\nATG\n")
    with pytest.raises(ValueError, match="lacks tree tips"):
        bundle.prepare_full_alignment(source, tree, output)


@pytest.mark.parametrize("text", [">a\nAT\n", ">a\nATG\n>a\nATG\n", ">a\nATG\n>b\nATGATG\n"])
def test_full_alignment_rejects_invalid_axes(tmp_path, text):
    path = tmp_path / "bad.fa"
    path.write_text(text)
    with pytest.raises(ValueError):
        bundle.read_alignment(path)


def test_bundle_relocation_and_sites_use_full_fit_and_real_structural_states(tmp_path):
    original = tmp_path / "original"
    make_bundle(original)
    moved = tmp_path / "moved"
    original.rename(moved)
    directory = bundle.structural_directory(moved)
    command = sites.build_csubst_sites_command(str(moved), str(moved), "1,2", 1, "3di20", 2, pdb="none")
    assert "--alignment_file" not in command
    assert command[command.index("--full_cds_alignment_file") + 1] == str(directory / "csubst.fasta")
    for suffix in ("treefile", "state", "rate", "iqtree", "log"):
        assert command[command.index("--iqtree_" + suffix) + 1] == str(directory / ("csubst." + suffix))
    output = tmp_path / "states.fa"
    bundle.write_structural_tip_alignment(directory / "csubst.fasta", output)
    assert output.read_text() == ">a\nA-X\n>b\nCD-\n"
    assert sites.resolve_csubst_genetic_code(moved) == 2
    (directory / "csubst.state").write_text("changed\n")
    with pytest.raises(ValueError, match="Missing or changed"):
        bundle.structural_directory(moved)


def test_old_bundle_requires_full_input_regeneration(tmp_path):
    (tmp_path / "csubst.input.json").write_text(json.dumps({"schema": "genegalleon-csubst-input-v1", "genetic_code": 1}))
    with pytest.raises(ValueError, match="Regenerate the iqtree_anc"):
        bundle.structural_directory(tmp_path)

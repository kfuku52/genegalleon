import hashlib
import json
import os
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


@pytest.mark.parametrize("size", [0, 262143, 524305])
def test_digest_preserves_sha256_for_complete_binary_content(tmp_path, size):
    payload = (bytes(range(256)) * ((size + 255) // 256))[:size]
    path = tmp_path / "state.bin"
    path.write_bytes(payload)
    assert bundle._digest(path) == hashlib.sha256(payload).hexdigest()


def test_finalize_manifest_hashes_every_complete_member(tmp_path):
    directory = make_bundle(tmp_path)
    (directory / "csubst.state").write_bytes(bytes(range(256)) * 4097)
    bundle.finalize(tmp_path, 2)
    expected = {f"csubst.{suffix}" for suffix in bundle.FIT_SUFFIXES} | {
        "csubst_3di_state_cache.npz", "inspect/csubst_alignment_3di.fa",
    }
    metadata = json.loads((tmp_path / "csubst.input.json").read_text())
    assert metadata == {
        "schema": bundle.SCHEMA,
        "genetic_code": 2,
        "structural_directory": bundle.STRUCTURAL_DIR,
        "site_coordinates": "full_cds_codon_1based",
        "structural_sha256": {
            name: hashlib.sha256((directory / name).read_bytes()).hexdigest()
            for name in expected
        },
    }


@pytest.mark.parametrize("name", ["csubst.state", "csubst_3di_state_cache.npz"])
def test_relocated_bundle_rejects_late_same_size_mutation(tmp_path, name):
    original = tmp_path / "original"
    directory = make_bundle(original)
    (directory / name).write_bytes(bytes(range(256)) * 4097)
    bundle.finalize(original, 2)
    moved = tmp_path / "moved"
    original.rename(moved)
    directory = bundle.structural_directory(moved)
    member = directory / name
    before = member.stat()
    with member.open("r+b") as handle:
        handle.seek(-1, os.SEEK_END)
        handle.write(b"\x00")
    os.utime(member, ns=(before.st_atime_ns, before.st_mtime_ns))
    assert member.stat().st_size == before.st_size
    assert member.stat().st_mtime_ns == before.st_mtime_ns
    with pytest.raises(ValueError, match="Missing or changed 3Di bundle input"):
        bundle.structural_directory(moved)


@pytest.mark.parametrize("name", [
    "csubst.state", "csubst_3di_state_cache.npz", "inspect/csubst_alignment_3di.fa",
])
@pytest.mark.parametrize("kind", ["missing", "empty", "directory"])
def test_bundle_rejects_unusable_required_members(tmp_path, name, kind):
    member = make_bundle(tmp_path) / name
    if kind == "empty":
        member.write_bytes(b"")
    else:
        member.unlink()
        if kind == "directory":
            member.mkdir()
    with pytest.raises(ValueError, match="Missing or changed 3Di bundle input"):
        bundle.structural_directory(tmp_path)


@pytest.mark.parametrize("extra", [False, True])
def test_bundle_requires_exact_hash_manifest_members(tmp_path, extra):
    make_bundle(tmp_path)
    path = tmp_path / "csubst.input.json"
    metadata = json.loads(path.read_text())
    if extra:
        metadata["structural_sha256"]["unexpected.state"] = "0" * 64
    else:
        metadata["structural_sha256"].pop("csubst.state")
    path.write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match="Incomplete 3Di bundle metadata"):
        bundle.structural_directory(tmp_path)


def test_finalize_rechecks_content_after_manifest_publication(tmp_path, monkeypatch):
    directory = make_bundle(tmp_path)
    member = directory / "csubst.state"
    member.write_bytes(bytes(range(256)) * 4097)
    write_text = Path.write_text

    def publish_then_mutate(path, *args, **kwargs):
        result = write_text(path, *args, **kwargs)
        if path == tmp_path / "csubst.input.json":
            with member.open("r+b") as handle:
                handle.seek(-1, os.SEEK_END)
                handle.write(b"\x00")
        return result

    monkeypatch.setattr(Path, "write_text", publish_then_mutate)
    with pytest.raises(ValueError, match="Missing or changed 3Di bundle input"):
        bundle.finalize(tmp_path, 2)

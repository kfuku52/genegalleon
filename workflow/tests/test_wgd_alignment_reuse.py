import json

import pytest

from workflow.support import wgd_alignment_reuse as reuse


def test_reuse_validates_prior_hash_and_identity(tmp_path):
    path = tmp_path / "cross_ks_aligned_pairs.code1.tsv"
    path.write_text("data\n")
    summary = tmp_path / "summary.json"
    identity = {"wrapper": "test-wrapper", "helpers": {"tbfast": "test-helper"}}
    summary.write_text(json.dumps({"plan": {"alignment_identity": identity},
                                  "output_hashes": {path.name: reuse.digest(path), "../escape": "ignored"}}))
    sources, inputs = reuse.reuse_sources(summary, identity)
    assert sources == {"1": str(path)}
    assert inputs == {str(summary): reuse.digest(summary), str(path): reuse.digest(path)}
    with pytest.raises(ValueError, match="same verified"):
        reuse.reuse_sources(summary, {**identity, "helpers": {"tbfast": "changed"}})
    path.write_text("tampered\n")
    with pytest.raises(ValueError, match="hash does not match"):
        reuse.reuse_sources(summary, identity)


def test_reuse_keeps_order_and_requires_exact_cds(tmp_path):
    path = tmp_path / "pairs.tsv"
    header = "\t".join(reuse.PAIR_COLUMNS) + "\n"
    content = "a|b\ta\tb\tATG---GCT\tATGAAAGCT\nc|d\tc\td\tATG\tATG\n"
    path.write_text(header + content)
    request = {("a", "b"): ("A", "B"), ("d", "c"): ("D", "C")}
    sequences = {"A": {"a": ("ATGGCT", "MA")}, "B": {"b": ("ATGAAAGCT", "MKA")}}
    matches = reuse.matching_rows(path, request.items(), sequences)
    assert set(matches) == {("a", "b")}
    assert matches["a", "b"]["sequence_1"] == "ATG---GCT"
    sequences["A"]["a"] = ("ATGGCC", "MA")
    assert reuse.matching_rows(path, request.items(), sequences) == {}
    path.write_text(header + content.replace("ATG---GCT", "AT--GGGCT"))
    with pytest.raises(ValueError, match="Malformed"):
        reuse.matching_rows(path, request.items(), sequences)
    path.write_text(header + content + content)
    with pytest.raises(ValueError, match="Duplicate"):
        reuse.matching_rows(path, request.items(), sequences)


def test_alignment_identity_tracks_helpers_and_declines_unknown_layout(tmp_path, monkeypatch):
    wrapper = tmp_path / "mafft"
    helpers = tmp_path / "libexec"
    helpers.mkdir()
    tbfast = helpers / "tbfast"
    tbfast.write_bytes(b"first")
    wrapper.write_text(f"#!/bin/bash\nprefix={helpers}\n")
    monkeypatch.setattr(reuse.shutil, "which", lambda _: str(wrapper))
    monkeypatch.delenv("MAFFT_BINARIES", raising=False)
    identity = {"source_hashes": {}, "mafft_version": "v7.x", "mafft_executable_sha256": reuse.digest(wrapper)}
    before, files = reuse.alignment_identity(identity)
    assert str(tbfast) in dict(files)
    tbfast.write_bytes(b"new helper")
    assert reuse.alignment_identity(identity)[0] != before
    wrapper.write_text("unsupported wrapper\n")
    assert reuse.alignment_identity(identity) == (None, ())

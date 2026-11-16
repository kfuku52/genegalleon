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


def test_reuse_merges_two_verified_bundles_and_rejects_conflicting_overlap(tmp_path):
    identity = {"test": True}
    header = "\t".join(reuse.PAIR_COLUMNS) + "\n"
    a, b = tmp_path / "old", tmp_path / "new"
    for folder in (a, b):
        folder.mkdir()
    files = [p / "cross_ks_aligned_pairs.code1.tsv" for p in (a, b)]
    summaries = [p / "summary.json" for p in (a, b)]
    shared = "u|v\tu\tv\tATG---ATG\tATGAAAATG\n"
    files[0].write_text(header + shared + "a|b\ta\tb\tATG\tATG\n")
    files[1].write_text(header + shared + "c|d\tc\td\tATG\tATG\n")
    for summary, file in zip(summaries, files, strict=True):
        summary.write_text(json.dumps({"plan": {"alignment_identity": identity},
                                       "output_hashes": {file.name: reuse.digest(file)}}))
    bundles, hashes = reuse.reuse_bundle([*summaries, summaries[0]], identity)
    assert bundles == {"1": list(map(str, files))}
    assert set(hashes) == {str(p) for p in (*files, *summaries)}
    requests = [((x, y), ("A", "B")) for x, y in ("ab", "cd", "uv")]
    sequences = {"A": {"a": ("ATG",), "c": ("ATG",), "u": ("ATGATG",)},
                 "B": {"b": ("ATG",), "d": ("ATG",), "v": ("ATGAAAATG",)}}
    assert reuse.matching_cached_rows(bundles["1"], requests, sequences).keys() == {("a", "b"), ("c", "d"), ("u", "v")}
    files[1].write_text(header + shared.replace("ATG---ATG", "ATGATG---"))
    with pytest.raises(ValueError, match="Conflicting verified"):
        reuse.matching_cached_rows(bundles["1"], requests, sequences)
    with pytest.raises(ValueError, match="hash does not match"):
        reuse.reuse_bundle(summaries, identity)


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


def test_likelihood_reuse_requires_exact_aligned_bytes_and_method(tmp_path):
    aligned = tmp_path / "cross_ks_aligned_pairs.code1.tsv"
    aligned.write_text("properly ordered original alignments\n")
    report = tmp_path / "cross_ks_pair_ds.code1.tsv"
    report.write_text("original report\n")
    summary = tmp_path / "summary.json"
    method = {"version": 1, "method": "gy94"}
    summary.write_text(json.dumps({"plan": {"alignment_identity": {"test": True}, "cross_ks_likelihood_identity": method},
                                  "output_hashes": {p.name: reuse.digest(p) for p in (aligned, report)}}))
    _, hashes = reuse.reuse_bundle([summary], {"test": True})
    assert reuse.likelihood_reuse([summary], {**method, "method": "gy94-regularized"}, hashes) == {}
    source = reuse.likelihood_reuse([summary], method, hashes)
    out = tmp_path / "result.tsv"
    assert reuse.reuse_likelihood_report(source, 1, aligned, out)
    assert out.read_bytes() == report.read_bytes()
    assert not reuse.reuse_likelihood_report(source, 2, aligned, out)
    report.write_text("tampered estimate\n")
    with pytest.raises(ValueError, match="changed while copying"):
        reuse.reuse_likelihood_report(source, 1, aligned, out)
    with pytest.raises(ValueError, match="hash does not match"):
        reuse.likelihood_reuse([summary], method, hashes)
    aligned.write_text("different current alignments\n")
    assert not reuse.reuse_likelihood_report(source, 1, aligned, out)
    summary.write_text("{}")
    with pytest.raises(ValueError, match="summary changed"):
        reuse.likelihood_reuse([summary], method, hashes)

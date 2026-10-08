"""Frozen full-record benchmark proofs, ownership, failure and publication."""
import copy
import hashlib
import json
import os
from pathlib import Path

import pytest

from workflow.benchmarks import benchmark_rescue_storage as bench


@pytest.fixture
def publication(tmp_path):
    root = tmp_path / "producer"
    directory = root / "rescued" / "Plant_example"
    directory.mkdir(parents=True)
    plan = root / "plan.json"
    plan.write_text(json.dumps({"request": {"tools": {"implementation": "fixture-frozen"}}}))
    records = []
    for i in range(8):
        status = "accepted" if i in {0, 4} else "accepted_alternative_path" if i == 5 else "unresolved"
        record = {"query": f"relative_gene_{i}", "id": f"prediction_{i}", "seqid": "chr",
                  "strand": "+", "cds": [[20, 32, 0]], "sequence": "ATGAAACCCTAA",
                  "coverage": .99, "identity": .9, "frameshift": False,
                  "paf": f"#PAF\trelative_gene_{i}\t4\t0\t4\t+\tchr\t100\t20\t32\t4\t4\t60",
                  "status": status, "problems": [] if status == "accepted" else ["missing_start"],
                  "partial_evidence": {"partial": i in {1, 3}},
                  "evidence": {"species": "Donor_example", "source": "nearest", "query": f"relative_gene_{i}"},
                  "gene_id": "rescue_same_gene", "locus_support": {"external_donor_species_count": 2},
                  "path_selection": {"representative_ambiguous": True, "canonical_policy": "longest_cds"},
                  "alternative_coding_paths": [{"cds": [[23, 32, 0]], "support": [{"species": "Other_donor"}]}],
                  "unknown_preserved": {"unicode": "β🙂", "large_integer": 2**70, "negative_zero": -0.0}}
        record["raw_prediction"] = {key: copy.deepcopy(record[key]) for key in
                                    ("query", "seqid", "strand", "cds", "paf", "coverage", "identity", "frameshift", "evidence")}
        records.append(record)
    partials = [records[3], records[1], records[3]]
    revisions = [{**copy.deepcopy(records[0]), "status": "healthy_representative_preserved",
                  "candidate_protein": "MKP", "rna_adoption_eligible": False}]
    values = {"models": records, "partial": partials, "revision": revisions}
    files = {}
    for kind, rows in values.items():
        p = directory / bench.MEMBERS[kind]
        p.write_text(json.dumps(rows, sort_keys=True, indent=2) + "\n")
        files[p.name] = bench.digest(p)
    # Execution input is inventoried but must never enter the compact worker.
    (directory / "regions.fa").write_text(">large_execution_input\nACGTACGT\n")
    files["regions.fa"] = bench.digest(directory / "regions.fa")
    receipt = directory / "receipt.json"
    receipt.write_text(json.dumps({"key": {"plan": bench.digest(plan), "species": "Plant_example"}, "files": files}))
    frozen = bench.freeze_inputs(root, "Plant_example", expected_plan=bench.digest(plan), expected_receipt=bench.digest(receipt))
    return root, directory, frozen, values


def config(frozen):
    import rescue_model_store
    result = {"inputs": frozen, "codec": "gzip", "shard_bytes": 512, "restart_after": 2,
              "module_sha256": bench.digest(Path(rescue_model_store.__file__)),
              "driver_sha256": bench.digest(Path(bench.__file__))}
    result["frozen_sha256"] = bench.configuration_sha(result)
    return result


def test_full_lifecycle_preserves_all_fields_orders_and_selected_science(publication, tmp_path):
    _, source, frozen, _ = publication
    before = {p.name: (bench.identity(p), bench.digest(p)) for p in source.iterdir() if p.is_file()}
    configuration = config(frozen)
    baseline = bench.worker(configuration, "legacy_read", tmp_path / "legacy")
    creation = bench.worker(configuration, "restart", tmp_path / "scratch_worker")
    assert creation["interrupted_publication_absent"] is True
    assert creation["interrupted_prefix_records"] == 3
    assert creation["manifest_counts"]["partial"] == 3
    assert not (tmp_path / "scratch_worker" / "regions.fa").exists()
    assert not (tmp_path / "scratch_worker" / "models.json").exists()
    bench.assert_parity(baseline, creation)
    read = bench.worker(configuration, "store_read", tmp_path / "scratch_worker")
    bench.assert_parity(baseline, read)
    assert read["summaries"]["accepted"]["records"] == 2
    assert read["summaries"]["models"]["status_counts"]["accepted_alternative_path"] == 1
    assert read["maximum_owned_disk_sampled"]["bytes"] > 0
    assert (tmp_path / "scratch_worker" / ".benchmark-tmp").is_dir()
    destination = tmp_path / "canonical_worker"
    copied = bench.publish_store(tmp_path / "scratch_worker", destination)
    assert copied["files"] >= 6
    bench.assert_parity(baseline, bench.worker(configuration, "store_read", destination))
    assert before == {p.name: (bench.identity(p), bench.digest(p)) for p in source.iterdir() if p.is_file()}


@pytest.mark.parametrize("field", ["canonical_sha256", "ordered_fields_sha256", "qc_sequence_cds_sha256"])
def test_full_proof_rejects_changed_scientific_values_or_order(publication, field):
    _, _, frozen, _ = publication
    rows = bench.SourceRows(frozen["files"]["models"])
    for _ in rows:
        pass
    baseline = {"summaries": {kind: bench.scan_rows(bench.SourceRows(spec)) for kind, spec in frozen["files"].items()}}
    current = copy.deepcopy(baseline)
    current["summaries"]["models"][field] = "0" * 64
    with pytest.raises(ValueError, match="parity differs"):
        bench.assert_parity(baseline, current)


def test_full_digest_detects_record_order_and_unknown_fields(publication):
    _, _, _, values = publication
    rows = values["models"]
    assert bench.scan_rows(rows)["canonical_sha256"] != bench.scan_rows(list(reversed(rows)))["canonical_sha256"]
    changed = copy.deepcopy(rows)
    changed[0]["unknown_preserved"]["large_integer"] += 1
    assert bench.scan_rows(rows)["canonical_sha256"] != bench.scan_rows(changed)["canonical_sha256"]
    reordered = [{key: row[key] for key in reversed(row)} for row in rows]
    assert bench.scan_rows(rows)["canonical_sha256"] == bench.scan_rows(reordered)["canonical_sha256"]
    assert bench.scan_rows(rows)["ordered_fields_sha256"] != bench.scan_rows(reordered)["ordered_fields_sha256"]


def test_source_raw_checksum_rejects_same_size_mtime_corruption(publication):
    _, _, frozen, _ = publication
    spec = copy.deepcopy(frozen["files"]["models"])
    path = Path(spec["path"])
    before = path.stat()
    original = path.read_bytes()
    corrupted = original.replace(b'"identity": 0.9', b'"identity": 0.8', 1)
    assert corrupted != original and len(corrupted) == len(original)
    path.write_bytes(corrupted)
    os.utime(path, ns=(before.st_atime_ns, before.st_mtime_ns))
    with pytest.raises(OSError, match="Frozen model input changed"):
        list(bench.SourceRows(spec))
    # Exercise the independent byte checksum even with a freshly observed stat.
    spec["identity"] = bench.identity(path)
    with pytest.raises(ValueError, match="checksum differs"):
        list(bench.SourceRows(spec))


def test_freeze_requires_independent_receipt_and_plan(publication):
    root, _, frozen, _ = publication
    with pytest.raises(ValueError, match="digests are required"):
        bench.freeze_inputs(root, "Plant_example", expected_plan=None, expected_receipt=None)
    with pytest.raises(ValueError, match="plan/receipt differs"):
        bench.freeze_inputs(root, "Plant_example", expected_plan="0" * 64, expected_receipt=frozen["receipt"]["sha256"])
    with pytest.raises(ValueError, match="Unsafe species"):
        bench.freeze_inputs(root, "../Plant_example", expected_plan=frozen["plan"]["sha256"], expected_receipt=frozen["receipt"]["sha256"])


def test_missing_bound_member_and_unbound_optional_are_rejected(publication):
    root, source, frozen, _ = publication
    (source / "models.json").unlink()
    with pytest.raises(ValueError, match="Missing producer member"):
        bench.freeze_inputs(root, "Plant_example", expected_plan=frozen["plan"]["sha256"], expected_receipt=frozen["receipt"]["sha256"])


def test_publication_never_overwrites_and_does_not_follow_symlinks(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    (source / "data").write_text("frozen")
    destination = tmp_path / "destination"
    destination.mkdir()
    with pytest.raises(ValueError, match="already exists"):
        bench.publish_store(source, destination)
    (source / "link").symlink_to(source / "data")
    with pytest.raises(ValueError, match="symlink"):
        bench.publish_store(source, tmp_path / "other")


def test_canonical_copy_corruption_keeps_source_and_no_committed_target(tmp_path, monkeypatch):
    source = tmp_path / "source"
    source.mkdir()
    (source / "data").write_text("frozen")
    before = bench.identity(source / "data")
    original_copy = bench.shutil.copytree
    def corrupt_copy(src, dest, **kwargs):
        value = original_copy(src, dest, **kwargs)
        p = Path(dest) / "data"
        token = p.stat()
        p.write_text("broken")
        os.utime(p, ns=(token.st_atime_ns, token.st_mtime_ns))
        return value
    monkeypatch.setattr(bench.shutil, "copytree", corrupt_copy)
    destination = tmp_path / "canonical"
    with pytest.raises(ValueError, match="publication differs"):
        bench.publish_store(source, destination)
    assert not destination.exists()
    assert bench.identity(source / "data") == before
    assert (source / "data").read_text() == "frozen"


def test_inventory_charges_hardlinks_once_and_no_symlink_target(tmp_path):
    (tmp_path / "file").write_bytes(b"abc")
    os.link(tmp_path / "file", tmp_path / "hardlink")
    (tmp_path / "symlink").symlink_to(tmp_path / "file")
    values = bench.inventory(tmp_path)
    assert values["regular_files"] == 2 and values["hardlink_aliases"] == 1
    assert values["bytes"] == 3 and values["symlinks"] == 1


def test_destination_delayed_ctime_requires_matching_every_byte(tmp_path, monkeypatch):
    source = tmp_path / "source"
    source.mkdir()
    (source / "data").write_text("frozen")
    actual = bench.identity
    calls = 0
    def delayed_ctime(path):
        nonlocal calls
        token = actual(path)
        if any(p.name.startswith(".publishing-") for p in Path(path).parents):
            calls += 1
            token[-1] += calls
        return token
    monkeypatch.setattr(bench, "identity", delayed_ctime)
    destination = tmp_path / "canonical"
    bench.publish_store(source, destination)
    assert (destination / "data").read_text() == "frozen" and calls >= 2


def test_source_ctime_changes_remain_strict_even_with_identical_bytes(tmp_path, monkeypatch):
    source = tmp_path / "source"
    source.mkdir()
    data = source / "data"
    data.write_text("frozen")
    actual = bench.identity
    calls = 0
    def changing_source(path):
        nonlocal calls
        token = actual(path)
        if Path(path) == data:
            calls += 1
            token[-1] += calls
        return token
    monkeypatch.setattr(bench, "identity", changing_source)
    with pytest.raises(OSError, match="File changed while hashing"):
        bench.publish_store(source, tmp_path / "canonical")


def test_driver_freeze_and_complete_cli_results_are_canonical(publication, tmp_path):
    root, _, frozen, _ = publication
    output, scratch = tmp_path / "canonical", tmp_path / "owned_scratch"
    args = ["--source-root", str(root), "--species", "Plant_example", "--output-directory", str(output),
            "--scratch-directory", str(scratch), "--codec", "gzip", "--repeats", "2", "--restart-after", "2",
            "--shard-bytes", "512", "--expected-source-plan-sha256", frozen["plan"]["sha256"],
            "--expected-source-receipt-sha256", frozen["receipt"]["sha256"]]
    args.append("--owned-copy-cold-advice")
    bench.main(args + ["--freeze-only"])
    assert (output / "frozen_inputs.json").is_file()
    bench.main(args)
    result = json.loads((output / "result.json").read_text())
    assert result["all_record_values_and_orders_equal"] is True
    assert result["canonical_publication"]["files"] >= 6
    assert result["legacy_model_stream_bytes"] == sum(spec["identity"][2] for spec in frozen["files"].values())
    assert result["maximum_total_owned_disk_sampled"]["bytes"] > 0
    assert result["owned_input_advice"]["advice"] == "POSIX_FADV_DONTNEED"
    assert result["owned_input_copy"]["capacity"]["bytes"] == result["legacy_model_stream_bytes"]
    assert result["legacy_warmed"]["summaries"] == result["legacy"]["summaries"]
    assert "Not performed" in result["consumer_pipeline_and_busco_reexecution"]
    with pytest.raises(ValueError, match="Completed benchmark exists"):
        bench.main(args)


def test_implementation_change_is_not_benchmark_resume(publication, tmp_path):
    _, _, frozen, _ = publication
    configuration = config(frozen)
    configuration["module_sha256"] = hashlib.sha256(b"other version").hexdigest()
    configuration["frozen_sha256"] = bench.configuration_sha(configuration)
    with pytest.raises(ValueError, match="implementation changed"):
        bench.worker(configuration, "legacy_read", tmp_path / "worker")


def test_advice_targets_only_newly_owned_verified_inodes(publication, tmp_path, monkeypatch):
    _, _, frozen, _ = publication
    configuration = config(frozen)
    destination = tmp_path / "owned_copies"
    copied, _ = bench.owned_input_copies(configuration, destination)
    advised = []
    def capture(descriptor, offset, length, advice):
        info = os.fstat(descriptor)
        advised.append((info.st_dev, info.st_ino))
        assert offset == length == 0 and advice == os.POSIX_FADV_DONTNEED
    monkeypatch.setattr(bench.os, "posix_fadvise", capture)
    proof = bench.advise_owned_copies(destination, configuration["frozen_sha256"])
    assert len(advised) == 3 and "not guaranteed" in proof["cache_regime"]
    source_inodes = {tuple(spec["identity"][:2]) for spec in frozen["files"].values()}
    assert not set(advised) & source_inodes
    assert set(advised) == {tuple(spec["identity"][:2]) for spec in copied["files"].values()}
    with pytest.raises(ValueError, match="owner differs"):
        bench.advise_owned_copies(destination, "0" * 64)


def test_advice_rejects_alias_to_original_and_corruption(publication, tmp_path):
    _, _, frozen, _ = publication
    configuration = config(frozen)
    destination = tmp_path / "owned_copies"
    _, receipt = bench.owned_input_copies(configuration, destination)
    original = Path(frozen["files"]["models"]["path"])
    target = Path(receipt["files"]["models"]["destination"])
    target.unlink()
    os.link(original, target)
    with pytest.raises(ValueError, match="ownership/content changed"):
        bench.advise_owned_copies(destination, configuration["frozen_sha256"])

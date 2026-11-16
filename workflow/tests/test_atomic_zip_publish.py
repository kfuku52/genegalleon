from __future__ import annotations

import os
import zipfile
from pathlib import Path

import pytest

from workflow.support import atomic_zip_publish as PUBLISH


def _write_zip(path: Path, payload: bytes = b"payload\n") -> None:
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("result/value.txt", payload)


def test_publish_zip_verifies_and_atomically_replaces_destination(tmp_path: Path):
    source = tmp_path / "source.zip"
    destination = tmp_path / "published" / "result.zip"
    _write_zip(source)
    destination.parent.mkdir()
    destination.write_bytes(b"old")

    result = PUBLISH.publish_zip(source, destination)

    assert result == destination
    assert source.is_file()
    with zipfile.ZipFile(destination) as archive:
        assert archive.read("result/value.txt") == b"payload\n"
    assert not list(destination.parent.glob(".result.zip.partial.*"))


def test_publish_zip_remove_source_happens_after_success(tmp_path: Path):
    source = tmp_path / "source.zip"
    destination = tmp_path / "result.zip"
    _write_zip(source)

    PUBLISH.publish_zip(source, destination, remove_source=True)

    assert destination.is_file()
    assert not source.exists()


def test_remove_source_preserves_destination_reached_through_symlinked_parent(
    tmp_path: Path,
):
    real_parent = tmp_path / "real"
    real_parent.mkdir()
    alias_parent = tmp_path / "alias"
    alias_parent.symlink_to(real_parent, target_is_directory=True)
    source = real_parent / "result.zip"
    destination = alias_parent / "result.zip"
    _write_zip(source)

    PUBLISH.publish_zip(source, destination, remove_source=True)

    assert source.is_file()
    assert destination.is_file()
    with zipfile.ZipFile(destination) as archive:
        assert archive.read("result/value.txt") == b"payload\n"


def test_corrupt_source_does_not_replace_destination(tmp_path: Path):
    source = tmp_path / "source.zip"
    destination = tmp_path / "result.zip"
    _write_zip(source)
    source.write_bytes(source.read_bytes()[:-12])
    destination.write_bytes(b"old")

    with pytest.raises(PUBLISH.ZipPublishError, match="Failed to verify"):
        PUBLISH.publish_zip(source, destination, remove_source=True)

    assert destination.read_bytes() == b"old"
    assert source.is_file()
    assert not list(tmp_path.glob(".result.zip.partial.*"))


def test_copy_failure_does_not_replace_destination(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    source = tmp_path / "source.zip"
    destination = tmp_path / "result.zip"
    _write_zip(source)
    destination.write_bytes(b"old")

    def fail_copy(input_handle, output_handle, length=0):
        output_handle.write(input_handle.read(8))
        raise OSError("injected copy failure")

    monkeypatch.setattr(PUBLISH.shutil, "copyfileobj", fail_copy)

    with pytest.raises(OSError, match="injected copy failure"):
        PUBLISH.publish_zip(source, destination)

    assert destination.read_bytes() == b"old"
    assert not list(tmp_path.glob(".result.zip.partial.*"))


def test_same_size_source_change_with_restored_mtime_is_detected(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    source = tmp_path / "source.zip"
    destination = tmp_path / "result.zip"
    _write_zip(source)
    destination.write_bytes(b"old")
    original_copy = PUBLISH.shutil.copyfileobj

    def copy_then_mutate(input_handle, output_handle, length=0):
        original_copy(input_handle, output_handle, length=length)
        original_mtime_ns = source.stat().st_mtime_ns
        source.write_bytes(b"X" * source.stat().st_size)
        os.utime(source, ns=(original_mtime_ns, original_mtime_ns))

    monkeypatch.setattr(PUBLISH.shutil, "copyfileobj", copy_then_mutate)

    with pytest.raises(PUBLISH.ZipPublishError, match="changed while being published"):
        PUBLISH.publish_zip(source, destination)

    assert destination.read_bytes() == b"old"
    assert not list(tmp_path.glob(".result.zip.partial.*"))


def test_expected_prefix_rejects_unexpected_members_without_publication(tmp_path: Path):
    source = tmp_path / "source.zip"
    destination = tmp_path / "result.zip"
    with zipfile.ZipFile(source, "w") as archive:
        archive.writestr("../escaped.txt", b"unsafe")
    destination.write_bytes(b"old")

    with pytest.raises(PUBLISH.ZipPublishError, match="Failed to verify"):
        PUBLISH.publish_zip(source, destination, expected_prefix="result")

    assert destination.read_bytes() == b"old"


def test_expected_prefix_rejects_empty_archive_without_publication(tmp_path: Path):
    source = tmp_path / "source.zip"
    destination = tmp_path / "result.zip"
    with zipfile.ZipFile(source, "w"):
        pass
    destination.write_bytes(b"old")

    with pytest.raises(PUBLISH.ZipPublishError, match="expected prefix"):
        PUBLISH.publish_zip(source, destination, expected_prefix="result", remove_source=True)

    assert source.is_file()
    assert destination.read_bytes() == b"old"
    assert not list(tmp_path.glob(".result.zip.partial.*"))


@pytest.mark.parametrize("storage", ["raw", "zip"])
def test_diagnostic_archive_omits_published_payload_and_reconstitutes_exact_bytes(tmp_path, storage):
    import json

    from workflow.support import diagnostic_archive
    from workflow.support.gene_family_output_store import convert_storage_to_zip, orthogroup_id_from_name

    stage, logical = tmp_path / "scan", tmp_path / "orthogroup"
    stage.mkdir()
    (logical / "csubst_scan").mkdir(parents=True)
    table = b"site\tvalue\n1\t0.125\n" * 100
    (stage / "sites.tsv").write_bytes(table)
    (stage / "inference.json").write_text('{"calibration":"none"}\n')
    alignment = tmp_path / "source.fa"
    alignment.write_text(">a\nACTG\n")
    member, name = "csubst_scan", "OG0001_csubst_scan.tsv"
    published = logical / member / name
    published.write_bytes(table)
    compact, full = tmp_path / "compact.zip", tmp_path / "full.zip"
    diagnostic_archive.pack(stage, compact, external=[("sites.tsv", f"{member}/{name}")],
                            inputs=[("inputs/alignment.fa", alignment)])
    assert not (stage / "inputs").exists()
    with zipfile.ZipFile(compact) as archive:
        assert "scan/sites.tsv" not in archive.namelist()
        assert len(json.loads(archive.read("scan/external_members.json"))["members"]) == 1
    if storage == "zip":
        convert_storage_to_zip(logical, "orthogroup", ["OG0001"], orthogroup_id_from_name)
        assert not published.exists()

    diagnostic_archive.materialize(compact, logical, full)
    with zipfile.ZipFile(full) as archive:
        assert archive.read("scan/sites.tsv") == table
        assert archive.read("scan/inputs/alignment.fa") == alignment.read_bytes()
        assert archive.read("scan/inference.json") == (stage / "inference.json").read_bytes()
        assert "scan/external_members.json" not in archive.namelist()
    if storage == "raw":
        published.write_bytes(b"stale")
        with pytest.raises(ValueError, match="differs from audit"):
            diagnostic_archive.materialize(compact, logical, tmp_path / "invalid.zip")
        assert not (tmp_path / "invalid.zip").exists()


def test_diagnostic_archive_materializes_self_contained_legacy_zip(tmp_path):
    from workflow.support import diagnostic_archive

    source, result = tmp_path / "source.zip", tmp_path / "result.zip"
    _write_zip(source)
    diagnostic_archive.materialize(source, tmp_path / "not-needed", result)
    with zipfile.ZipFile(result) as archive:
        assert archive.read("result/value.txt") == b"payload\n"
    assert not (tmp_path / "not-needed").exists()


def test_expression_trait_workflow_packs_three_audits_as_one_without_scientific_tables(tmp_path):
    import subprocess

    repo = Path(__file__).resolve().parents[2]
    text = (repo / "workflow/core/gg_gene_evolution_core.sh").read_text()
    start = text.index('  python "${gg_support_dir}/diagnostic_archive.py" pack')
    end = text.index("  mv_out_bundle", start)
    scratch, archive_path = tmp_path / "OGTEST.rsc", tmp_path / "bundle.zip"
    scratch.mkdir()
    names = ("rsc.audit.jsonl", "pgls.method-audit.jsonl", "species_expression.audit.tsv")
    for name in (*names, "regression.tsv"):
        (scratch / name).write_text(name + "\n")
    environment = dict(os.environ, gg_support_dir=str(repo / "workflow/support"), rsc_work_dir=str(scratch),
                       expression_trait_audit_tmp=str(archive_path), rsc_combined_audit=str(scratch / names[0]),
                       pgls_method_audit_tmp=str(scratch / names[1]), species_expression_audit_tmp=str(scratch / names[2]))
    subprocess.run(["bash", "-eu", "-c", text[start:end]], env=environment, check=True)
    with zipfile.ZipFile(archive_path) as archive:
        assert set(archive.namelist()) == {scratch.name + "/" + name for name in names}
        for name in names:
            assert archive.read(scratch.name + "/" + name) == (scratch / name).read_bytes()


def test_expression_trait_workflow_checks_previous_audit_layout_before_reuse(tmp_path):
    import subprocess

    script = (Path(__file__).resolve().parents[2] / "workflow/core/gg_gene_evolution_core.sh").read_text()
    start = script.index('rsc_audit_output="diagnostic_archive=')
    stop = script.index('if [[ ${rsc_needs_update} -eq 1', start)
    previous = tmp_path / "rsc_audit/OGTEST_rsc.audit.jsonl"
    previous.parent.mkdir()
    previous.touch()
    packed = tmp_path / "packed.zip"
    environment = dict(os.environ, dir_output_active=str(tmp_path), og_id="OGTEST",
                       file_og_expression_trait_audit=str(packed))
    harness = ('rsc_provenance_args=(--output trait_origins=a --output method_status=b '
               '--output species_expression_summary=c); '
               'gg_artifact_prepare_stage() { shift 2; printf "%s\\n" "$@"; };\n')
    out = subprocess.check_output(["bash", "-eu", "-c", harness + script[start:stop]], env=environment, text=True)
    assert [row.split("=", 1)[0] for row in out.splitlines()[1::2]] == [
        "trait_origins", "audit", "method_status", "method_audit", "species_expression_summary", "species_expression_audit"]
    assert f"audit={previous}" in out.splitlines()
    packed.touch()
    out = subprocess.check_output(["bash", "-eu", "-c", harness + script[start:stop]], env=environment, text=True)
    assert out.splitlines()[-2:] == ["--output", f"diagnostic_archive={packed}"]

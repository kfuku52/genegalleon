"""Exercise default ZIP lifecycles without changing research workspace data."""

import gzip
import os
import subprocess
import zipfile
from pathlib import Path

import pytest

from workflow.support import artifact_provenance as PROVENANCE
from workflow.support import directory_output_store as DIRECTORY
from workflow.support import orthofinder_output_store as STORE
from workflow.support import workspace_output_storage as WORKSPACE
from workflow.support.directory_storage_policy import DIRECTORY_SCOPES

REPO = Path(__file__).resolve().parents[2]


def _fixture(root, names):
    contents = {}
    for name in names:
        raw = root / name
        (raw / "empty").mkdir(parents=True)
        payloads = {"result.tsv": b"id\tvalue\nBUSCO1\t7\n", "sequence.fa.gz": gzip.compress(b">s_gene\nATG\n")}
        for filename, payload in payloads.items():
            path = raw / filename
            path.write_bytes(payload)
            path.chmod(0o640)
            os.utime(path, ns=(1700000000123456789, 1700000000123456789))
            contents[str(path.relative_to(root))] = (payload, path.stat().st_mode & 0o777, path.stat().st_mtime_ns)
    return contents


@pytest.mark.parametrize("kind", DIRECTORY_SCOPES)
def test_all_retained_directory_scopes_preserve_content_metadata_and_provenance(tmp_path, kind):
    root = tmp_path / kind
    names = DIRECTORY_SCOPES[kind]
    contents = _fixture(root, names)
    before = PROVENANCE.sha256_path(root)
    directory_before = {name: PROVENANCE.sha256_path(root / name) for name in names}
    DIRECTORY.convert(root, kind=kind, workers=2)
    assert PROVENANCE.sha256_path(root) == before
    for name in names:
        assert not (root / name).exists()
        assert PROVENANCE.declared_path_exists(root / name)
        assert PROVENANCE.sha256_path(root / name) == directory_before[name]
        with zipfile.ZipFile(root / (name + ".zip")) as archive:
            assert archive.getinfo(Path(name).name + "/sequence.fa.gz").compress_type == zipfile.ZIP_STORED
    with DIRECTORY.open_binary(root, names[0] + "/result.tsv", kind=kind) as handle:
        assert handle.read() == contents[names[0] + "/result.tsv"][0]
    with DIRECTORY.materialized_directory(root, names[0], tmp_path, kind=kind) as native:
        assert (native / "result.tsv").read_bytes() == contents[names[0] + "/result.tsv"][0]
    assert not native.exists()
    assert (root / (names[0] + ".zip")).is_file()
    DIRECTORY.convert(root, kind=kind, to="raw", workers=2)
    assert PROVENANCE.sha256_path(root) == before
    for relative, expected in contents.items():
        p = root / relative
        assert (p.read_bytes(), p.stat().st_mode & 0o777, p.stat().st_mtime_ns) == expected
    assert all((root / name / "empty").is_dir() for name in names)
    assert not list(root.rglob("*.zip"))


def test_workspace_default_conversion_and_rollback_include_every_new_scope(tmp_path):
    workspace = tmp_path / "workspace"
    output = workspace / "output"
    roots = {"genome_evolution": output / "genome_evolution", "hgt": output / "hgt",
             "query2family_hgt": output / "query2family_hgt", "analyses/hgt_saved": output / "analyses/hgt_saved"}
    for target, root in roots.items():
        kind = "analyses" if target.startswith("analyses/") else "hgt" if target == "query2family_hgt" else target
        _fixture(root, DIRECTORY_SCOPES[kind])
    # Curated inputs, native runtime caches and environments never join a result catalog.
    protected = [workspace / "input/species_cds/s.fa", workspace / "downloads/busco_downloads/lineages/hmms/a.hmm",
                 output / "analyses/hgt_saved/.venv/bin/python"]
    for path in protected:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"preserve")
    for to in ("zip", "raw"):
        args = WORKSPACE.build_parser().parse_args(["convert", "--workspace", str(workspace), "--to", to,
                                                  "--verification", "deep", "--workers", "2"])
        code, report, _, _ = WORKSPACE.run(args)
        assert code == 0, report
        assert {row["name"] for row in report["targets"]} == set(roots)
        assert all(row["status"] == "complete" for row in report["targets"])
        assert all(path.read_bytes() == b"preserve" for path in protected)
    assert not list(output.rglob("*.zip"))


def test_generic_pack_conflict_keeps_both_forms_and_rejects_unmanaged_scope(tmp_path):
    root = tmp_path / "hgt"
    _fixture(root, ["tree_plot"])
    DIRECTORY.convert(root, kind="hgt")
    raw = root / "tree_plot"
    raw.mkdir()
    (raw / "new.tsv").write_bytes(b"new producer")
    with pytest.raises(STORE.OrthoFinderStorageError, match="Unarchived source"):
        DIRECTORY.convert(root, kind="hgt")
    assert (raw / "new.tsv").read_bytes() == b"new producer"
    assert (root / "tree_plot.zip").is_file()
    with pytest.raises(STORE.OrthoFinderStorageError, match="Invalid directory selection"):
        DIRECTORY.convert(root, kind="hgt", directories=["busco_mafft"])


def _function(core, name):
    text = (REPO / "workflow/core" / core).read_text()
    return name + "() {" + text.split(name + "() {", 1)[1].split("\n}\n", 1)[0] + "\n}\n"


@pytest.mark.parametrize("mode", ["zip", "files", "raw"])
def test_genome_core_restores_native_inputs_and_defaults_to_zip_on_finish(tmp_path, mode):
    root = tmp_path / "output/genome_evolution"
    _fixture(root, ["busco_mafft"])
    DIRECTORY.convert(root, kind="genome_evolution")
    script = "set -euo pipefail\n" + _function("gg_genome_evolution_core.sh", "genome_evolution_storage_require_native_paths")
    script += _function("gg_genome_evolution_core.sh", "genome_evolution_storage_prepare")
    script += _function("gg_genome_evolution_core.sh", "genome_evolution_storage_finish")
    script += '''source "${gg_support_dir}/../gg_common_params.sh"
    genome_evolution_storage_mode="${GG_COMMON_GENOME_EVOLUTION_OUTPUT_STORAGE}"
    genome_evolution_storage_lock_token=""
    genome_evolution_storage_prepare
    test -f "${dir_genome_evolution}/busco_mafft/result.tsv"
    printf 'retained\\n' > "${dir_genome_evolution}/busco_mafft/new.tsv"
    genome_evolution_storage_finish
    '''
    env = dict(os.environ, gg_support_dir=str(REPO / "workflow/support"), gg_workspace_output_dir=str(root.parent),
               dir_genome_evolution=str(root), run_busco_dupaware_mafft="1")
    env.pop("GG_COMMON_GENOME_EVOLUTION_OUTPUT_STORAGE", None)
    if mode != "zip":
        env["GG_COMMON_GENOME_EVOLUTION_OUTPUT_STORAGE"] = mode
    completed = subprocess.run(["bash", "-c", script], env=env, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    if mode == "zip":
        assert not (root / "busco_mafft").exists()
        assert {row["path"] for row in STORE.archive_records(root, "busco_mafft", deep=True)["files"]} == {
            "result.tsv", "sequence.fa.gz", "new.tsv"}
    else:
        assert (root / "busco_mafft/new.tsv").read_text() == "retained\n"
        assert not (root / "busco_mafft.zip").exists()


@pytest.mark.parametrize("exit_status", [0, 7])
def test_hgt_exit_waits_for_family_producer_then_packs_default_output(tmp_path, exit_status):
    root = tmp_path / "hgt"
    root.mkdir()
    script = "set -euo pipefail\n" + _function("gg_hgt_core.sh", "hgt_cleanup_on_exit")
    script += '''source "${gg_support_dir}/../gg_common_params.sh"
    hgt_storage_mode="${GG_COMMON_HGT_OUTPUT_STORAGE}"
    hgt_storage_lock_file="$(dirname "${dir_hgt}")/.gg_hgt_directory_storage.lock"
    hgt_storage_lock_token=$(python "${gg_support_dir}/shared_namespace_lock.py" acquire-exclusive "$hgt_storage_lock_file" --owner-pid "$$")
    hgt_cleanup_materialization_run() { :; }
    mkdir -p "${dir_hgt}/tree_plot"
    (sleep 0.1; printf 'complete\\n' > "${dir_hgt}/tree_plot/OG1.pdf") &
    hgt_plot_pids=("$!")
    trap hgt_cleanup_on_exit EXIT
    '''
    env = dict(os.environ, gg_support_dir=str(REPO / "workflow/support"), dir_hgt=str(root))
    env.pop("GG_COMMON_HGT_OUTPUT_STORAGE", None)
    script += f"exit {exit_status}\n"
    completed = subprocess.run(["bash", "-c", script], env=env, capture_output=True, text=True)
    assert completed.returncode == exit_status, completed.stderr
    assert not (root / "tree_plot").exists()
    assert STORE.archive_records(root, "tree_plot", deep=True)["files"][0]["size"] == len(b"complete\n")


def test_restore_never_overwrites_a_file_created_during_publication(tmp_path, monkeypatch):
    root = tmp_path / "hgt"
    _fixture(root, ["tree_plot_input"])
    DIRECTORY.convert(root, kind="hgt")
    target = root / "tree_plot_input/result.tsv"
    def race(operation):
        def publish(source, destination, *args, **kwargs):
            if Path(destination) == target:
                target.write_bytes(b"concurrent producer result")
            return operation(source, destination, *args, **kwargs)
        return publish
    monkeypatch.setattr(os, "replace", race(os.replace))
    monkeypatch.setattr(os, "link", race(os.link))
    with pytest.raises(STORE.OrthoFinderStorageError, match="conflict|appeared"):
        DIRECTORY.convert(root, kind="hgt", to="raw")
    assert target.read_bytes() == b"concurrent producer result"
    assert (root / "tree_plot_input.zip").is_file()


@pytest.mark.parametrize("reader", [PROVENANCE.sha256_path, PROVENANCE.declared_path_exists])
def test_direct_archived_directory_is_fenced_during_observation(tmp_path, reader):
    from workflow_observation import observe_files
    root = tmp_path / "hgt"
    _fixture(root, ["tree_plot_input"])
    DIRECTORY.convert(root, kind="hgt")
    archive = root / "tree_plot_input.zip"
    with pytest.raises(ValueError, match="changed during observation"):
        with observe_files():
            reader(root / "tree_plot_input")
            original = archive.stat()
            os.utime(archive, ns=(original.st_atime_ns, original.st_mtime_ns + 1000000000))


def test_nested_managed_zip_preserves_direct_and_ancestor_provenance(tmp_path):
    root = tmp_path / "hgt"
    nested = root / "trait_focus"
    _fixture(nested, ["tree_plot_input"])
    DIRECTORY.convert(nested, kind="hgt", directories=["tree_plot_input"])
    before_root = PROVENANCE.sha256_path(root)
    before_directory = PROVENANCE.sha256_path(nested)
    DIRECTORY.convert(root, kind="hgt", directories=["trait_focus"])
    assert PROVENANCE.sha256_path(root) == before_root
    assert PROVENANCE.sha256_path(nested) == before_directory
    DIRECTORY.convert(root, kind="hgt", directories=["trait_focus"], to="raw")
    assert PROVENANCE.sha256_path(root) == before_root


def test_invalid_storage_namespace_is_rejected_before_any_path_changes(tmp_path):
    root = tmp_path / "hgt"
    _fixture(root, ["tree_plot_input"])
    with pytest.raises(STORE.OrthoFinderStorageError, match="namespace"):
        DIRECTORY.convert(root, kind="typo")
    assert (root / "tree_plot_input/result.tsv").is_file()


def test_genome_core_keeps_archives_intact_when_busco_stages_are_disabled(tmp_path):
    root = tmp_path / "output/genome_evolution"
    _fixture(root, ["busco_mafft"])
    DIRECTORY.convert(root, kind="genome_evolution")
    archive = root / "busco_mafft.zip"
    before = (archive.read_bytes(), archive.stat().st_ino, archive.stat().st_mtime_ns)
    script = "set -euo pipefail\n" + _function("gg_genome_evolution_core.sh", "genome_evolution_storage_require_native_paths")
    script += _function("gg_genome_evolution_core.sh", "genome_evolution_storage_prepare")
    script += _function("gg_genome_evolution_core.sh", "genome_evolution_storage_finish")
    script += '''genome_evolution_storage_mode=zip
    genome_evolution_storage_lock_token=""
    run_busco_dupaware_mafft=0
    genome_evolution_storage_prepare
    test ! -d "${dir_genome_evolution}/busco_mafft"
    genome_evolution_storage_finish
    '''
    env = {key: value for key, value in os.environ.items() if not key.startswith("run_busco_dupaware_")}
    env.update(gg_support_dir=str(REPO / "workflow/support"), gg_workspace_output_dir=str(root.parent),
               dir_genome_evolution=str(root))
    completed = subprocess.run(["bash", "-c", script], env=env, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    assert (archive.read_bytes(), archive.stat().st_ino, archive.stat().st_mtime_ns) == before


def test_nested_large_unmarked_zip_is_preserved_as_an_ordinary_file(tmp_path):
    root = tmp_path / "hgt"
    raw = root / "trait_focus"
    raw.mkdir(parents=True)
    ordinary = raw / "tree_plot_input.zip"
    with zipfile.ZipFile(ordinary, "w", compression=zipfile.ZIP_STORED) as archive:
        archive.writestr("ordinary.bin", b"sequence-data" * 900000)
    before = PROVENANCE.sha256_path(raw)
    DIRECTORY.convert(root, kind="hgt", directories=["trait_focus"])
    assert PROVENANCE.sha256_path(raw) == before


def test_late_enabled_busco_stage_can_restore_native_paths(tmp_path):
    root = tmp_path / "output/genome_evolution"
    _fixture(root, ["busco_mafft"])
    DIRECTORY.convert(root, kind="genome_evolution")
    script = "set -euo pipefail\n"
    for name in ("genome_evolution_storage_require_native_paths", "genome_evolution_storage_prepare", "genome_evolution_storage_finish"):
        script += _function("gg_genome_evolution_core.sh", name)
    script += '''genome_evolution_storage_mode=zip
    genome_evolution_storage_lock_token=""
    run_busco_dupaware_mafft=0
    genome_evolution_storage_prepare
    test ! -d "${dir_genome_evolution}/busco_mafft"
    run_busco_dupaware_mafft=1
    genome_evolution_storage_require_native_paths
    test -f "${dir_genome_evolution}/busco_mafft/result.tsv"
    genome_evolution_storage_finish
    '''
    env = {key: value for key, value in os.environ.items() if not key.startswith("run_busco_dupaware_")}
    env.update(gg_support_dir=str(REPO / "workflow/support"), gg_workspace_output_dir=str(root.parent),
               dir_genome_evolution=str(root))
    completed = subprocess.run(["bash", "-c", script], env=env, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    assert (root / "busco_mafft.zip").is_file()


def test_restoration_reads_each_zip_payload_once_and_verifies_restored_bytes(tmp_path, monkeypatch):
    from collections import Counter
    root = tmp_path / "hgt"
    _fixture(root, ["tree_plot_input"])
    DIRECTORY.convert(root, kind="hgt")
    reads = Counter()
    original = zipfile.ZipFile.open
    def count(archive, name, *args, **kwargs):
        member = name.filename if isinstance(name, zipfile.ZipInfo) else name
        if member.startswith("tree_plot_input/"):
            reads[member] += 1
        return original(archive, name, *args, **kwargs)
    monkeypatch.setattr(zipfile.ZipFile, "open", count)
    DIRECTORY.convert(root, kind="hgt", to="raw")
    assert reads == {"tree_plot_input/result.tsv": 1, "tree_plot_input/sequence.fa.gz": 1}
    assert (root / "tree_plot_input/result.tsv").read_bytes() == b"id\tvalue\nBUSCO1\t7\n"


def test_crc_corruption_is_not_published_and_recovery_zip_survives(tmp_path):
    root = tmp_path / "hgt"
    _fixture(root, ["tree_plot_input"])
    DIRECTORY.convert(root, kind="hgt")
    archive_path = root / "tree_plot_input.zip"
    with zipfile.ZipFile(archive_path) as archive:
        info = archive.getinfo("tree_plot_input/sequence.fa.gz")
    with archive_path.open("r+b") as handle:
        handle.seek(info.header_offset + 26)
        header = handle.read(4)
        offset = info.header_offset + 30 + int.from_bytes(header[:2], "little") + int.from_bytes(header[2:], "little")
        handle.seek(offset)
        byte = handle.read(1)
        handle.seek(offset)
        handle.write(bytes([byte[0] ^ 1]))
    with pytest.raises(STORE.OrthoFinderStorageError, match="Invalid ZIP payload"):
        DIRECTORY.convert(root, kind="hgt", to="raw")
    assert archive_path.is_file()
    assert not (root / "tree_plot_input/sequence.fa.gz").exists()


@pytest.mark.parametrize("managed", [True, False])
def test_disabled_stage_still_audits_known_archived_outputs(tmp_path, managed):
    root = tmp_path / "output/genome_evolution"
    if managed:
        _fixture(root, ["busco_mafft"])
        DIRECTORY.convert(root, kind="genome_evolution")
        output = root / "busco_mafft"
    else:
        root.mkdir(parents=True)
        output = root / "ordinary_output"
        with zipfile.ZipFile(output.with_suffix(".zip"), "w") as archive:
            archive.writestr("ordinary.txt", "unrelated")
    script = '''set -euo pipefail
    source "${gg_support_dir}/gg_util.sh"
    gg_artifact_set_needs_update() { printf -v "$1" 1; printf audited > "$RESULT"; }
    run_demo=0
    artifact_stale_policy=rebuild
    gg_artifact_prepare_stage needs_update run_demo --output "result=${TARGET}"
    printf '%s' "$run_demo"
    '''
    result = tmp_path / "called"
    completed = subprocess.run(["bash", "-c", script], capture_output=True, text=True,
                               env=dict(os.environ, gg_support_dir=str(REPO / "workflow/support"),
                                        RESULT=str(result), TARGET=str(output)))
    assert completed.returncode == 0, completed.stderr
    assert result.exists() == managed
    assert completed.stdout.endswith("1" if managed else "0")


def test_archived_file_contracts_and_manifests_remain_readable_without_restoration(tmp_path):
    root = tmp_path / "hgt"
    _fixture(root, ["tree_plot_input"])
    manifest = root / "artifact_provenance/stage.json"
    manifest.parent.mkdir()
    manifest.write_text('{"schema_version": 1, "step": "hgt"}')
    file = root / "tree_plot_input/result.tsv"
    before = PROVENANCE.sha256_path(file)
    DIRECTORY.convert(root, kind="hgt")
    assert PROVENANCE.declared_path_exists(file)
    assert PROVENANCE.sha256_path(file) == before
    with PROVENANCE.open_declared_binary(file) as handle:
        assert handle.read() == b"id\tvalue\nBUSCO1\t7\n"
    assert PROVENANCE.load_manifest(manifest)["step"] == "hgt"
    entry = dict(label="result", scope="workspace", path=file.relative_to(tmp_path).as_posix(),
                 artifact_type="file", sha256=before[0], size_bytes=before[1])
    payload = dict(schema_version=1, inputs=[], outputs=[entry], parameters={})
    assert PROVENANCE.audit_manifest(payload, PROVENANCE.GeneFamilyOutputStore(root), root, tmp_path)[0] == "current"
    assert not file.exists() and not manifest.exists()

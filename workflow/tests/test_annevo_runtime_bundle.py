import gzip
import io
import json
import subprocess
import tarfile
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "container/scripts/build_annevo_runtime_bundle.py"


def module():
    spec = spec_from_file_location("build_annevo_runtime_bundle", SCRIPT)
    loaded = module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


def official_fixture(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    for name in ["annotation.py", "prediction.py", "decoding.py", "src/utils.py", "LICENSE",
                 "training/not_runtime.py", "example/not_runtime.fa", "saved_model/ANNEVO_Magnoliopsida.pt"]:
        path = source / name
        path.parent.mkdir(exist_ok=True)
        path.write_bytes(b"test license\n" if name == "LICENSE" else b"x" * 2048)
    for argv in [["git", "init", "-q", str(source)], ["git", "-C", str(source), "remote", "add", "origin", "https://github.com/xjtu-omics/ANNEVO"],
                 ["git", "-C", str(source), "add", "."], ["git", "-C", str(source), "-c", "user.name=fixture", "-c", "user.email=fixture@example.invalid", "commit", "-qm", "offline fixture"]]:
        subprocess.run(argv, check=True, capture_output=True)
    return source


def environment_archive(path, extra=None):
    with tarfile.open(path, "w:gz") as archive:
        for name in ["bin/python", "bin/conda-unpack", "lib/retained-license"]:
            member = tarfile.TarInfo(name)
            data = b"synthetic environment fixture " + name.encode()
            member.size = len(data)
            member.mode = 0o755 if name.startswith("bin/") else 0o644
            archive.addfile(member, io.BytesIO(data))
        if extra:
            archive.addfile(extra)


def shared_fixture(tmp_path):
    loaded = module()
    prefix = tmp_path / "environment"
    (prefix / "lib").mkdir(parents=True)
    (prefix / "conda-meta").mkdir()
    names = sorted(loaded.SHARED_OPENMP_FILES)
    for index, name in enumerate(names):
        path = prefix / name
        path.write_bytes(b"\x7fELF OFFLINE physical native library " + name.encode())
        path.chmod(0o644 if index == 0 else 0o755)
    for owner in ("intel-openmp", "llvm-openmp"):
        paths = [{"_path": name, "path_type": "hardlink" if owner == "intel-openmp" else "softlink",
                  "sha256_in_prefix": loaded.sha256(prefix / name) if owner == "intel-openmp" else "b" * 64}
                 for name in names]
        data = {"name": owner, "version": "OFFLINE", "build": "fixture", "files": names,
                "paths_data": {"paths": paths}}
        (prefix / "conda-meta" / (owner + ".json")).write_text(json.dumps(data))
    return loaded, prefix


def test_only_effective_prefix_library_bytes_are_preserved_once_with_conflict_proof(tmp_path, monkeypatch):
    loaded, prefix = shared_fixture(tmp_path)
    preserved = loaded.preserved_shared_files(prefix)
    assert set(preserved) == loaded.SHARED_OPENMP_FILES
    assert all(value["matches_providers"] == ["intel-openmp"] for value in preserved.values())
    assert all(len(value["providers"]) == 2 for value in preserved.values())
    source = official_fixture(tmp_path)
    source_info, _ = loaded.inspect_source(source, "main")
    calls = []
    monkeypatch.setattr(loaded, "checked_run", lambda argv, **kwargs: calls.append(argv) or "fixture")
    loaded.pack_environment(prefix, tmp_path / "packed.tar.gz", "conda-pack", exclude_files=preserved)
    assert calls[-1].count("--exclude") == 2
    assert not {"--ignore-missing-files", "--ignore-editable-packages"} & set(calls[-1])
    env, payload = tmp_path / "packed.tar.gz", tmp_path / "payload.tar.gz"
    environment_archive(env)
    loaded.write_payload(env, source, source_info, {}, {"effective_prefix_shared_files": preserved}, payload,
                         prefix=prefix, shared_files=preserved)
    with tarfile.open(payload, "r:gz") as archive:
        for name, proof in preserved.items():
            full = "environment/" + name
            assert archive.getnames().count(full) == 1
            assert loaded.hashlib.sha256(archive.extractfile(full).read()).hexdigest() == proof["sha256"]
        pack = json.load(archive.extractfile("provenance/conda-pack.json"))
        assert pack["effective_prefix_shared_files"] == preserved
        assert "environment/lib/retained-license" in archive.getnames() and "source/LICENSE" in archive.getnames()


@pytest.mark.parametrize("failure", ["unmatched_SHA", "metadata_placeholder", "physical_placeholder", "provider", "mode", "type", "extra_overlap"])
def test_shared_library_preservation_rejects_unqualified_prefix_changes(tmp_path, failure):
    loaded, prefix = shared_fixture(tmp_path)
    name = "lib/libiomp5.so"
    physical = prefix / name
    metadata_path = prefix / "conda-meta/intel-openmp.json"
    metadata = json.loads(metadata_path.read_text())
    if failure == "unmatched_SHA":
        physical.write_bytes(b"unknown physical bytes")
    elif failure == "metadata_placeholder":
        metadata["paths_data"]["paths"][0]["prefix_placeholder"] = "/old/prefix"
    elif failure == "physical_placeholder":
        physical.write_bytes(b"\x7fELF " + str(prefix).encode())
        metadata["paths_data"]["paths"][0]["sha256_in_prefix"] = loaded.sha256(physical)
    elif failure == "provider":
        metadata["name"] = "unrecognized-openmp"
    elif failure == "mode":
        physical.chmod(0o4777)
    elif failure == "type":
        physical.unlink()
        physical.symlink_to("libomptarget.so")
    else:
        for file in (prefix / "conda-meta").glob("*.json"):
            data = json.loads(file.read_text())
            data["files"].append("lib/unrecognized_overlap.so")
            file.write_text(json.dumps(data))
        metadata["files"].append("lib/unrecognized_overlap.so")
    metadata_path.write_text(json.dumps(metadata))
    with pytest.raises(loaded.BundleError):
        loaded.preserved_shared_files(prefix)


@pytest.mark.parametrize("failure", ["conflicting_archive_entry", "physical_drift"])
def test_shared_library_preservation_never_accepts_conflicting_archive_entries_or_changed_bytes(tmp_path, failure):
    loaded, prefix = shared_fixture(tmp_path)
    preserved = loaded.preserved_shared_files(prefix)
    source = official_fixture(tmp_path)
    source_info, _ = loaded.inspect_source(source, "main")
    extra = tarfile.TarInfo("lib/libiomp5.so") if failure == "conflicting_archive_entry" else None
    env = tmp_path / "packed.tar.gz"
    environment_archive(env, extra)
    if failure == "physical_drift":
        (prefix / "lib/libiomp5.so").write_bytes(b"changed after qualification")
    with pytest.raises(loaded.BundleError):
        loaded.write_payload(env, source, source_info, {}, {}, tmp_path / "payload.tar.gz", prefix=prefix, shared_files=preserved)


def test_bundle_excludes_models_training_and_examples_and_binds_official_model(tmp_path, monkeypatch):
    loaded = module()
    source = official_fixture(tmp_path)
    abi = {"system": "Linux", "machine": "x86_64", "libc": {"name": "glibc", "version": "2.39"},
           "image_digest": "sha256:" + "a" * 64}
    monkeypatch.setattr(loaded, "inspect_environment", lambda *args: {"machine": "x86_64", "torch_version": "fixture", "compiled_cuda": "fixture"})

    def fixture_pack(prefix, output, executable):
        environment_archive(output)
        return {"version": "offline synthetic fixture", "strict_missing_files": True}

    monkeypatch.setattr(loaded, "pack_environment", fixture_pack)
    manifest_path, manifest = loaded.build_bundle(tmp_path / "env", source, tmp_path / "out", abi, "fixture", part_bytes=193)
    assert json.loads(manifest_path.read_text()) == manifest
    assert manifest["source"]["repository"] == loaded.OFFICIAL_REPOSITORY
    model = manifest["models"]["Magnoliopsida"]
    assert manifest["source"]["commit"] in model["official_url"]
    assert model["sha256"] == loaded.sha256(source / "saved_model/ANNEVO_Magnoliopsida.pt")
    parts = [tmp_path / "out" / p["name"] for p in manifest["parts"]]
    joined = b"".join(p.read_bytes() for p in parts)
    assert loaded.hashlib.sha256(joined).hexdigest() == manifest["bundle"]["sha256"]
    assert len(joined) == manifest["bundle"]["size"]
    assert all(p.stat().st_size < 2 * 1024 ** 3 for p in parts)
    with tarfile.open(fileobj=io.BytesIO(joined), mode="r:gz") as archive:
        names = archive.getnames()
        assert "source/LICENSE" in names and "environment/bin/conda-unpack" in names
        assert not any("saved_model" in n or "/training/" in n or "/example/" in n for n in names)
        payload = json.load(archive.extractfile("provenance/payload-manifest.json"))
        assert {x["path"] for x in payload["entries"]} == set(names) - {"provenance/payload-manifest.json"}
        assert manifest["payload_uncompressed_size"] == sum(member.size for member in archive if member.isfile())
    # Exercise the actual installer-side exact inventory contract, including implicit parents.
    consumer_path = SCRIPT.parents[2] / "workflow/support/annevo_runtime_cache.py"
    monkeypatch.syspath_prepend(str(consumer_path.parent))
    consumer_spec = spec_from_file_location("annevo_runtime_cache", consumer_path)
    consumer = module_from_spec(consumer_spec)
    consumer_spec.loader.exec_module(consumer)
    archive_path = tmp_path / "joined.tar.gz"
    archive_path.write_bytes(joined)
    installed = tmp_path / "unpacked"
    installed.mkdir()
    consumer._safe_extract(archive_path, installed)
    consumer._check_payload(installed, manifest)
    assert manifest["qualification"]["GPU_inference"] == "not_performed"
    assert manifest["relocation"]["rename_after_unpack"] is False


def test_input_pax_long_paths_and_link_targets_roundtrip_through_actual_consumer(tmp_path, monkeypatch):
    loaded = module()
    source = official_fixture(tmp_path)
    long_name = "lib/" + "long_native_library_" * 7 + ".so"
    data = b"OFFLINE fixture: native library bytes must survive PAX relocation"
    abi = {"system": "Linux", "machine": "x86_64", "libc": {"name": "glibc", "version": "2.39"},
           "image_digest": "sha256:" + "a" * 64}
    monkeypatch.setattr(loaded, "inspect_environment", lambda *a: {"machine": "x86_64", "torch_version": "fixture", "compiled_cuda": "fixture"})

    def fixture_pack(prefix, output, executable):
        with tarfile.open(output, "w:gz", format=tarfile.PAX_FORMAT) as archive:
            for name in ["bin/python", "bin/conda-unpack", long_name]:
                member = tarfile.TarInfo(name)
                member.size = len(data)
                member.mode = 0o755
                member.pax_headers = {"path": name, "mtime": "1774882182.0684626", "uid": "123", "uname": "old-owner"}
                archive.addfile(member, io.BytesIO(data))
            for name, kind, link in [("lib/hardlink", tarfile.LNKTYPE, long_name),
                                     ("lib/symlink", tarfile.SYMTYPE, long_name.removeprefix("lib/")),
                                     ("lib/rootlink", tarfile.SYMTYPE, "..")]:
                member = tarfile.TarInfo(name)
                member.type, member.linkname = kind, link
                member.pax_headers = {"path": name, "linkpath": link, "mtime": "1774882182.0684626"}
                archive.addfile(member)
        return {"strict_missing_files": True}

    monkeypatch.setattr(loaded, "pack_environment", fixture_pack)
    _, manifest = loaded.build_bundle(tmp_path / "env", source, tmp_path / "out", abi, "offline-pax-fixture")
    archive_path = tmp_path / "joined.tar.gz"
    archive_path.write_bytes(b"".join((tmp_path / "out" / p["name"]).read_bytes() for p in manifest["parts"]))
    with tarfile.open(archive_path, "r:gz") as archive:
        members = archive.getmembers()
        assert len({member.name for member in members}) == len(members)
        assert all(Path(member.name).parts[0] in {"environment", "source", "provenance"} for member in members)
        for member in members:
            assert (member.uid, member.gid, member.mtime, member.uname, member.gname) == (0, 0, 0, "", "")
        hardlink = archive.getmember("environment/lib/hardlink")
        assert hardlink.linkname == "environment/" + long_name
        assert hardlink.pax_headers["linkpath"] == hardlink.linkname
        assert archive.getmember("environment/" + long_name).pax_headers["path"] == "environment/" + long_name
        assert archive.getmember("environment/lib/symlink").linkname == long_name.removeprefix("lib/")
    consumer_path = SCRIPT.parents[2] / "workflow/support/annevo_runtime_cache.py"
    monkeypatch.syspath_prepend(str(consumer_path.parent))
    consumer_spec = spec_from_file_location("pax_consumer_fixture", consumer_path)
    consumer = module_from_spec(consumer_spec)
    consumer_spec.loader.exec_module(consumer)
    installed = tmp_path / "unpacked"
    installed.mkdir()
    consumer._safe_extract(archive_path, installed, manifest["payload_uncompressed_size"])
    consumer._check_payload(installed, manifest)
    assert (installed / "environment" / long_name).read_bytes() == data
    assert (installed / "environment/lib/hardlink").samefile(installed / "environment" / long_name)
    assert (installed / "environment/lib/symlink").resolve() == installed / "environment" / long_name
    assert (installed / "environment/lib/rootlink").resolve() == installed / "environment"


@pytest.mark.parametrize("pax_name", ["../outside", "/outside"])
def test_input_pax_path_overrides_are_validated_before_rewriting(tmp_path, pax_name):
    loaded = module()
    path = tmp_path / "bad-pax.tar.gz"
    with tarfile.open(path, "w:gz", format=tarfile.PAX_FORMAT) as archive:
        member = tarfile.TarInfo("lib/apparently-safe")
        member.pax_headers = {"path": pax_name}
        archive.addfile(member)
    with tarfile.open(path, "r:gz") as archive, pytest.raises(loaded.BundleError, match="Unsafe"):
        list(loaded.environment_members(archive))


@pytest.mark.parametrize("failure", ["pax_path", "duplicate", "mode", "bytes", "escaping_link", "special_type", "owner", "total"])
def test_serialized_payload_preflight_rejects_header_or_inventory_drift(tmp_path, failure):
    loaded = module()
    source = official_fixture(tmp_path)
    source_info, _ = loaded.inspect_source(source, "main")
    environment = tmp_path / "environment.tar.gz"
    environment_archive(environment)
    payload = tmp_path / "payload.tar.gz"
    _, inventory, total = loaded.write_payload(environment, source, source_info, {}, {}, payload)
    loaded.verify_serialized_payload(payload, inventory, total)
    changed = tmp_path / "changed.tar.gz"
    with tarfile.open(payload, "r:gz") as original, tarfile.open(changed, "w:gz", format=tarfile.PAX_FORMAT) as target:
        for member in original:
            data = original.extractfile(member).read() if member.isfile() else None
            if member.name == "environment/bin/python":
                if failure == "pax_path":
                    member.pax_headers = {"path": "bin/python"}
                elif failure == "mode":
                    member.mode = 0o700
                elif failure == "bytes":
                    data = b"X" + data[1:]
                elif failure == "escaping_link":
                    member.type, member.size, member.linkname, data = tarfile.SYMTYPE, 0, "../../../outside", None
                elif failure == "special_type":
                    member.type, member.size, data = tarfile.CHRTYPE, 0, None
                elif failure == "owner":
                    member.pax_headers = {"uid": "123"}
                target.addfile(member, io.BytesIO(data) if data is not None else None)
                if failure == "duplicate":
                    target.addfile(member, io.BytesIO(data))
            else:
                target.addfile(member, io.BytesIO(data) if data is not None else None)
    with pytest.raises(loaded.BundleError):
        loaded.verify_serialized_payload(changed, inventory, total + 1 if failure == "total" else total)


def test_failed_serialized_payload_preflight_publishes_no_parts(tmp_path, monkeypatch):
    loaded = module()
    source = official_fixture(tmp_path)
    abi = {"system": "Linux", "machine": "x86_64", "libc": {"name": "glibc", "version": "2.39"},
           "image_digest": "sha256:" + "a" * 64}
    monkeypatch.setattr(loaded, "inspect_environment", lambda *a: {"machine": "x86_64", "torch_version": "fixture", "compiled_cuda": "fixture"})

    def pack(prefix, output, executable):
        environment_archive(output)
        return {}

    def fail(*a):
        raise loaded.BundleError("serialized fixture rejected")

    monkeypatch.setattr(loaded, "pack_environment", pack)
    monkeypatch.setattr(loaded, "verify_serialized_payload", fail)
    with pytest.raises(loaded.BundleError, match="serialized fixture rejected"):
        loaded.build_bundle(tmp_path / "env", source, tmp_path / "out", abi, "fixture")
    assert list((tmp_path / "out").iterdir()) == []


@pytest.mark.parametrize("kind", ["absolute", "parent", "symlink", "device", "hardlink"])
def test_packed_environment_rejects_unsafe_members(tmp_path, kind):
    loaded = module()
    bad = tarfile.TarInfo("bad")
    if kind == "absolute":
        bad.name = "/bad"
    elif kind == "parent":
        bad.name = "../bad"
    elif kind == "symlink":
        bad.type = tarfile.SYMTYPE
        bad.linkname = "../../outside"
    elif kind == "device":
        bad.type = tarfile.CHRTYPE
    else:
        bad.type = tarfile.LNKTYPE
        bad.linkname = "../outside"
    path = tmp_path / "environment.tar.gz"
    environment_archive(path, bad)
    with tarfile.open(path, "r:gz") as archive, pytest.raises(loaded.BundleError):
        list(loaded.environment_members(archive))


def test_model_and_native_source_git_binding_rejects_dirty_or_wrong_origin(tmp_path):
    loaded = module()
    source = official_fixture(tmp_path)
    (source / "annotation.py").write_text("modified\n")
    with pytest.raises(loaded.BundleError, match="changes"):
        loaded.inspect_source(source, "main")
    subprocess.run(["git", "-C", str(source), "checkout", "--", "annotation.py"], check=True)
    subprocess.run(["git", "-C", str(source), "remote", "set-url", "origin", "https://example.invalid/ANNEVO"], check=True)
    with pytest.raises(loaded.BundleError, match="official"):
        loaded.inspect_source(source, "main")


def test_conda_pack_never_weakens_validation_or_bakes_install_prefix(tmp_path, monkeypatch):
    loaded = module()
    calls = []
    options = []
    monkeypatch.setattr(loaded, "checked_run", lambda argv, **kwargs: calls.append(argv) or options.append(kwargs) or "fixture")
    loaded.pack_environment(tmp_path / "read-only-env", tmp_path / "env.tar.gz", "conda-pack")
    argv = calls[-1]
    assert "--prefix" in argv and "--format" in argv
    assert not {"--ignore-missing-files", "--ignore-editable-packages", "--dest-prefix"} & set(argv)
    assert options[-1]["env"]["TMPDIR"] == str(tmp_path)


def test_existing_bundle_output_is_not_overwritten(tmp_path):
    loaded = module()
    out = tmp_path / "out"
    out.mkdir()
    existing = out / "owned-by-other-publication"
    existing.write_bytes(b"retain")
    with pytest.raises(loaded.BundleError, match="empty"):
        loaded.build_bundle(tmp_path / "env", tmp_path / "source", out, {}, "profile")
    assert existing.read_bytes() == b"retain"


def test_manifest_and_archive_headers_do_not_contain_model_bytes(tmp_path, monkeypatch):
    loaded = module()
    source = official_fixture(tmp_path)
    info, _ = loaded.inspect_source(source, "main")
    env = tmp_path / "environment.tar.gz"
    environment_archive(env)
    output = tmp_path / "payload.tar.gz"
    original = loaded.tempfile.TemporaryFile

    def require_bound_storage(*args, **kwargs):
        assert kwargs["dir"] == tmp_path  # RAM-backed default /tmp is never used for CUDA libraries.
        return original(*args, **kwargs)

    monkeypatch.setattr(loaded.tempfile, "TemporaryFile", require_bound_storage)
    loaded.write_payload(env, source, info, {}, {}, output)
    assert output.read_bytes()[4:8] == b"\0" * 4  # normalized gzip mtime
    assert b"saved_model" not in gzip.decompress(output.read_bytes())


def test_output_inside_input_is_rejected_without_creating_any_directory(tmp_path):
    loaded = module()
    output = tmp_path / "environment" / "output"
    with pytest.raises(loaded.BundleError, match="outside"):
        loaded.build_bundle(tmp_path / "environment", tmp_path / "source", output, {}, "profile")
    assert not output.exists()


@pytest.mark.parametrize("bad", ["CPU", "OS", "architecture", "libc", "digest"])
def test_environment_metadata_rejects_cpu_fallback_or_wrong_base_abi(tmp_path, monkeypatch, bad):
    loaded = module()
    prefix = tmp_path / "environment"
    (prefix / "bin").mkdir(parents=True)
    (prefix / "conda-meta").mkdir()
    (prefix / "bin/python").write_bytes(b"fixture")
    abi = {"system": "Linux", "machine": "x86_64", "libc": {"name": "glibc", "version": "2.39"},
           "image_digest": "sha256:" + "a" * 64}
    observed = {**abi, "python_prefix": str(prefix), "compiled_cuda": "12.1"}
    if bad == "CPU":
        observed["compiled_cuda"] = None
    elif bad == "OS":
        observed["system"] = "Darwin"
    elif bad == "architecture":
        observed["machine"] = "aarch64"
    elif bad == "libc":
        observed["libc"] = {"name": "glibc", "version": "2.35"}
    else:
        abi["image_digest"] = "latest"
    monkeypatch.setattr(loaded, "checked_run", lambda *args, **kwargs: json.dumps(observed))
    with pytest.raises(loaded.BundleError):
        loaded.inspect_environment(prefix, tmp_path / "source", abi, tmp_path)

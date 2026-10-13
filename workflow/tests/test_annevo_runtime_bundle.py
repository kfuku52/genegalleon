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

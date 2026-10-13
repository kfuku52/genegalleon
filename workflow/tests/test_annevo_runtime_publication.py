import hashlib
import io
import json
import subprocess
from copy import deepcopy
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]


def load(name):
    spec = spec_from_file_location(name, ROOT / "container/scripts" / (name + ".py"))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def candidate(tmp_path):
    content = b"OFFLINE PUBLICATION FIXTURE"
    name = "archive.part0000"
    (tmp_path / name).write_bytes(content)
    sha = hashlib.sha256(content).hexdigest()
    manifest = {"kind": "genegalleon-annevo-runtime", "platform": "linux/amd64",
                "source": {"commit": "b" * 40}, "bundle": {"sha256": sha, "size": len(content)},
                "parts": [{"name": name, "sha256": sha, "size": len(content)}], "models": {}}
    path = tmp_path / "fixture.manifest.json"
    path.write_text(json.dumps(manifest))
    return path, manifest


def anonymous_graph(module, path, manifest):
    raw = path.read_bytes()
    metadata_digest = "sha256:" + hashlib.sha256(raw).hexdigest()
    part = manifest["parts"][0]
    oci = {"artifactType": module.ARTIFACT_TYPE, "layers": [
        {"mediaType": module.MANIFEST_TYPE, "digest": metadata_digest, "size": len(raw)},
        {"mediaType": module.PART_TYPE, "digest": "sha256:" + part["sha256"], "size": part["size"]}]}
    blob = json.dumps(oci).encode()
    digest = "sha256:" + hashlib.sha256(blob).hexdigest()
    values = {("manifests", digest): blob, ("blobs", metadata_digest): raw,
              ("blobs", "sha256:" + part["sha256"]): (path.parent / part["name"]).read_bytes()}
    return digest, values


def qualification(path, manifest):
    """Synthetic offline promotion receipts; never experimental GPU evidence."""
    manifest["models"] = {"Magnoliopsida": {"official_url": "https://example.invalid/OFFLINE.pt",
                                          "sha256": "d" * 64, "size": 123,
                                          "path": "source/saved_model/ANNEVO_Magnoliopsida.pt"}}
    manifest["cuda_profile"] = {"torch_version": "OFFLINE+cu121", "compiled_cuda": "12.1"}
    path.write_text(json.dumps(manifest))
    resolved = deepcopy(manifest)
    resolved["parts"][0]["url"] = "https://example.invalid/OFFLINE.part"
    selected = {"runtime_digest": manifest["bundle"]["sha256"], "manifest": resolved,
                "source_commit": manifest["source"]["commit"], "model_profile": "Magnoliopsida",
                "model_sha256": "d" * 64, "model_size": 123}
    canonical = (json.dumps(resolved, sort_keys=True, separators=(",", ":")) + "\n").encode()
    proof = {"qualification_completed": True, "integration_returncode": 0,
             "cold_first_install_verified": True, "warm_no_download_reuse_verified": True,
             "cold_warm_gff_bytes_equal": True,
             "publisher_manifest_sha256": hashlib.sha256(canonical).hexdigest(),
             "gpu_probe": {"verified": True, **manifest["cuda_profile"]}, "annotation": {"sha256": "e" * 64}}
    proof_path, prepared_path = path.parent / "OFFLINE_gpu_proof.json", path.parent / "OFFLINE_prepared.json"
    proof_path.write_text(json.dumps(proof))
    prepared_path.write_text(json.dumps(selected))
    return {"gpu_qualification_proof": proof_path, "prepared_receipt": prepared_path}


def test_anonymous_complete_manifest_parts_whole_archive_are_really_hashed(tmp_path):
    module = load("publish_annevo_runtime_bundle")
    path, manifest = candidate(tmp_path)
    digest, values = anonymous_graph(module, path, manifest)
    def fetch(kind, sha):
        return io.BytesIO(values[(kind, sha)])
    result = module.verify_anonymous("ghcr.io/kfuku52/genegalleon", digest, path, manifest, fetch=fetch)
    assert result["bundle_sha256"] == manifest["bundle"]["sha256"]
    assert result["anonymous_parts_verified"] == 1
    values[("blobs", "sha256:" + manifest["parts"][0]["sha256"])] = b"altered"
    with pytest.raises(ValueError, match="part bytes"):
        module.verify_anonymous("ghcr.io/kfuku52/genegalleon", digest, path, manifest, fetch=fetch)


@pytest.mark.parametrize("immutable_only", [True, False])
def test_anonymous_failure_never_updates_any_alias(tmp_path, monkeypatch, immutable_only):
    module = load("publish_annevo_runtime_bundle")
    path, manifest = candidate(tmp_path)
    evidence = qualification(path, manifest) if not immutable_only else {}
    calls = []
    monkeypatch.setattr(module, "anonymous_reader", lambda repository: None)

    def command(argv, **kwargs):
        calls.append(argv)
        return subprocess.CompletedProcess(argv, 0, stdout="sha256:" + "c" * 64, stderr="")

    def failure(*args, **kwargs):
        raise ValueError("anonymous bytes not verified")

    monkeypatch.setattr(module, "oras", command)
    monkeypatch.setattr(module, "verify_anonymous", failure)
    with pytest.raises(ValueError, match="anonymous bytes"):
        module.publish(path, "ghcr.io/kfuku52/genegalleon", "annevo-runtime-linux-amd64-latest",
                       promote=not immutable_only, promote_existing=not immutable_only, **evidence)
    assert not any(argv[0] in {"tag", "push"} for argv in calls)


def test_immutable_only_has_tag_locator_and_digest_but_never_promotes(tmp_path, monkeypatch):
    module = load("publish_annevo_runtime_bundle")
    path, manifest = candidate(tmp_path)
    calls = []
    monkeypatch.setattr(module, "anonymous_reader", lambda repository: None)
    monkeypatch.setattr(module, "oras", lambda argv, **kwargs: calls.append(argv) or subprocess.CompletedProcess(argv, 0, "sha256:" + "c" * 64, ""))
    monkeypatch.setattr(module, "verify_anonymous", lambda *args: {"fixture_verified": True})
    result = module.publish(path, "ghcr.io/kfuku52/genegalleon", "annevo-runtime-linux-amd64-latest", promote=False)
    assert result["immutable_ref"] == "oci://ghcr.io/kfuku52/genegalleon:annevo-runtime-linux-amd64-sha256-" + manifest["bundle"]["sha256"]
    assert result["immutable_digest"] == "sha256:" + "c" * 64
    assert result["moving_ref"] is None and result["moving_alias_updated"] is False
    assert not any(argv[0] in {"tag", "push"} for argv in calls)


def test_wrong_local_part_fails_before_any_registry_write(tmp_path, monkeypatch):
    module = load("publish_annevo_runtime_bundle")
    path, manifest = candidate(tmp_path)
    (tmp_path / manifest["parts"][0]["name"]).write_bytes(b"corrupt")
    calls = []
    monkeypatch.setattr(module, "anonymous_reader", lambda repository: None)
    monkeypatch.setattr(module, "oras", lambda argv, **kwargs: calls.append(argv))
    with pytest.raises(ValueError, match="Local part"):
        module.publish(path, "ghcr.io/kfuku52/genegalleon", "annevo-runtime-linux-amd64-latest", promote=False)
    assert calls == []


def test_promotion_requires_existing_candidate_and_auth_failure_is_not_missing(tmp_path, monkeypatch):
    module = load("publish_annevo_runtime_bundle")
    path, manifest = candidate(tmp_path)
    evidence = qualification(path, manifest)
    calls = []
    monkeypatch.setattr(module, "anonymous_reader", lambda repository: None)
    monkeypatch.setattr(module, "oras", lambda argv, **kwargs: calls.append(argv) or subprocess.CompletedProcess(argv, 1, "", "401 unauthorized"))
    with pytest.raises(ValueError, match="already exist"):
        module.publish(path, "ghcr.io/kfuku52/genegalleon", "annevo-runtime-linux-amd64-latest",
                       promote=True, promote_existing=True, **evidence)
    assert not any(argv[0] in {"tag", "push"} for argv in calls)
    calls.clear()
    with pytest.raises(ValueError, match="lookup failed"):
        module.publish(path, "ghcr.io/kfuku52/genegalleon", "annevo-runtime-linux-amd64-latest", promote=False)
    assert not any(argv[0] in {"tag", "push"} for argv in calls)


def test_explicit_promotion_happens_only_after_verification_and_keeps_image_tags(tmp_path, monkeypatch):
    module = load("publish_annevo_runtime_bundle")
    path, manifest = candidate(tmp_path)
    evidence = qualification(path, manifest)
    calls = []
    monkeypatch.setattr(module, "anonymous_reader", lambda repository: None)
    monkeypatch.setattr(module, "oras", lambda argv, **kwargs: calls.append(argv) or subprocess.CompletedProcess(argv, 0, "sha256:" + "c" * 64, ""))
    monkeypatch.setattr(module, "verify_anonymous", lambda *args: calls.append(["verified"]) or {})
    result = module.publish(path, "ghcr.io/kfuku52/genegalleon", "annevo-runtime-linux-amd64-latest",
                            promote=True, promote_existing=True, **evidence)
    assert result["moving_alias_updated"]
    assert [argv[0] for argv in calls] == ["resolve", "verified", "tag", "resolve"]
    assert calls[2][-1] == "annevo-runtime-linux-amd64-latest"
    assert not any(argument.endswith(":latest") or ":sha-" in argument for argv in calls for argument in argv)
    assert result["GPU_qualification"]["qualified_bundle_sha256"] == manifest["bundle"]["sha256"]


def test_promotion_without_exact_gpu_evidence_fails_before_registry_calls(tmp_path, monkeypatch):
    module = load("publish_annevo_runtime_bundle")
    path, _ = candidate(tmp_path)
    monkeypatch.setattr(module, "anonymous_reader", lambda _: pytest.fail("must fail before network"))
    with pytest.raises(ValueError, match="requires an already published"):
        module.publish(path, "ghcr.io/kfuku52/genegalleon", "annevo-runtime-linux-amd64-latest", promote=True)
    with pytest.raises(ValueError, match="requires the GPU"):
        module.publish(path, "ghcr.io/kfuku52/genegalleon", "annevo-runtime-linux-amd64-latest",
                       promote=True, promote_existing=True)


@pytest.mark.parametrize("field,value", [
    ("qualification_completed", False), ("integration_returncode", 1), ("integration_returncode", False),
    ("cold_first_install_verified", False), ("warm_no_download_reuse_verified", False),
    ("cold_warm_gff_bytes_equal", False), ("publisher_manifest_sha256", "f" * 64),
    ("gpu_probe", {"verified": True, "torch_version": "CPU", "compiled_cuda": None}),
    ("annotation", {"sha256": "missing"}),
])
def test_promotion_rejects_incomplete_or_different_measured_gpu_proof(tmp_path, field, value):
    module = load("publish_annevo_runtime_bundle")
    path, manifest = candidate(tmp_path)
    evidence = qualification(path, manifest)
    proof_path = evidence["gpu_qualification_proof"]
    proof = json.loads(proof_path.read_text())
    proof[field] = value
    proof_path.write_text(json.dumps(proof))
    with pytest.raises(ValueError):
        module.verify_gpu_qualification(manifest, proof_path, evidence["prepared_receipt"])


@pytest.mark.parametrize("field,value", [("runtime_digest", "f" * 64), ("source_commit", "f" * 40),
                                         ("model_sha256", "f" * 64), ("model_size", 456)])
def test_promotion_rejects_other_bundle_source_or_model_selection(tmp_path, field, value):
    module = load("publish_annevo_runtime_bundle")
    path, manifest = candidate(tmp_path)
    evidence = qualification(path, manifest)
    receipt_path = evidence["prepared_receipt"]
    selected = json.loads(receipt_path.read_text())
    selected[field] = value
    receipt_path.write_text(json.dumps(selected))
    with pytest.raises(ValueError, match="selection differs"):
        module.verify_gpu_qualification(manifest, evidence["gpu_qualification_proof"], receipt_path)


def test_gpu_receipt_manifest_metadata_must_match_candidate_even_with_self_consistent_proof(tmp_path):
    module = load("publish_annevo_runtime_bundle")
    path, manifest = candidate(tmp_path)
    evidence = qualification(path, manifest)
    receipt_path, proof_path = evidence["prepared_receipt"], evidence["gpu_qualification_proof"]
    selected, proof = json.loads(receipt_path.read_text()), json.loads(proof_path.read_text())
    selected["manifest"]["parts"][0]["sha256"] = "f" * 64
    canonical = (json.dumps(selected["manifest"], sort_keys=True, separators=(",", ":")) + "\n").encode()
    proof["publisher_manifest_sha256"] = hashlib.sha256(canonical).hexdigest()
    receipt_path.write_text(json.dumps(selected))
    proof_path.write_text(json.dumps(proof))
    with pytest.raises(ValueError, match="metadata differs"):
        module.verify_gpu_qualification(manifest, proof_path, receipt_path)


def test_ci_event_and_permissions_parse_to_trusted_candidate_only(tmp_path):
    document = yaml.load((ROOT / ".github/workflows/annevo-runtime-ghcr.yml").read_text(), Loader=yaml.BaseLoader)
    assert set(document["on"]) == {"push", "workflow_dispatch"}
    assert document["on"]["push"]["branches"] == ["main"]
    assert document["on"]["push"]["paths"]
    build, publish = document["jobs"]["build"], document["jobs"]["publish"]
    assert document["permissions"] == {"contents": "read"}
    assert build["permissions"] == {"contents": "read"}
    assert publish["permissions"] == {"contents": "read", "packages": "write"}
    assert build["runs-on"] == publish["runs-on"] == "ubuntu-24.04"
    assert "default_branch" in build["if"] and "default_branch" in publish["if"]
    assert publish["needs"] == "build"
    steps = publish["steps"]
    assert "--immutable-only" in steps[-1]["run"]
    assert "GHCR_TOKEN" not in json.dumps(build)
    assert all(step.get("with", {}).get("persist-credentials") == "false" for job in [build, publish]
               for step in job["steps"] if "actions/checkout@" in step.get("uses", ""))


def test_bootstrap_resolves_a_single_native_amd64_digest_without_changing_source_defaults(monkeypatch):
    module = load("bootstrap_annevo_runtime_distribution")
    payload = {"manifests": [{"platform": {"architecture": "amd64", "os": "linux"}, "digest": "sha256:" + "a" * 64},
                             {"platform": {"architecture": "arm64", "os": "linux"}, "digest": "sha256:" + "b" * 64}]}
    monkeypatch.setattr(module, "run", lambda *args, **kwargs: subprocess.CompletedProcess([], 0, json.dumps(payload)))
    reference, digest = module.immutable_amd64_reference("ghcr.io/kfuku52/genegalleon:latest")
    assert reference == "ghcr.io/kfuku52/genegalleon@" + digest
    assert digest == "sha256:" + "a" * 64
    assert "--ignore-missing-files" not in module.INNER and "--ignore-editable-packages" not in module.INNER
    assert "--file /work/source/ANNEVO.yml" in module.INNER

"""Convert one closed, independently owned Rescue worker to ZIP-backed storage."""

import argparse
import base64
import fcntl
import hashlib
import json
import os
import shutil
import stat
import tempfile
import zipfile
from contextlib import contextmanager
from pathlib import Path

if __package__:
    from . import rescue_model_store as store
else:
    import rescue_model_store as store


RECOVERY = ".model-store-loose-layout.json"
LOCK = ".model-store-archive.lock"


def _sha(raw):
    return hashlib.sha256(raw).hexdigest()


def _replace(path, raw):
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        pending = Path(handle.name)
        try:
            handle.write(raw)
            handle.flush()
            os.fsync(handle.fileno())
            os.replace(pending, path)
        finally:
            pending.unlink(missing_ok=True)
    store._sync_directory(path.parent)


def _file_bytes(path):
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Expected a regular model-store file: {path}")
    return path.read_bytes()


def _build_archive(directory, archive, shards, files):
    identities = {}
    with tempfile.NamedTemporaryFile(dir=directory) as handle:
        with zipfile.ZipFile(handle, "w", compression=zipfile.ZIP_STORED) as output:
            for shard in shards:
                name = shard["path"]
                source = store._member(directory, name)
                identities[name] = store._stat(source)
                store._copy_zip_member(source, output, name, files[name])
                if store._stat(source) != identities[name]:
                    raise ValueError(f"Model-store source changed during packing: {name}")
        handle.flush()
        os.fsync(handle.fileno())
        os.link(handle.name, archive)
    store._sync_directory(directory)
    return identities


def archive_closed_store(worker):
    """Retain exact original metadata for restoration and resume a partial retirement.

    A changed receipt is a new publication generation. The caller must exclude
    both writers and consumers that pin the previous receipt before proceeding.
    """
    with _locked_store(worker) as (root, directory):
        return _archive_locked(root, directory)


@contextmanager
def _locked_store(worker):
    logical = Path(worker)
    if logical.is_symlink():
        raise ValueError("The Rescue worker must be a real directory")
    root = logical.resolve(strict=True)
    directory = root / "model_store"
    if directory.is_symlink() or not directory.is_dir():
        raise ValueError("Expected the worker's real model_store directory")
    lock_path = directory / LOCK
    if lock_path.is_symlink() or (lock_path.exists() and not stat.S_ISREG(lock_path.lstat().st_mode)):
        raise ValueError("Unsafe model-store maintenance lock")
    with lock_path.open("a+b") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        yield root, directory


def _archive_locked(root, directory):
    manifest_path = directory / "manifest.json"
    receipt_path = root / "receipt.json"
    saved = directory / RECOVERY
    archive = directory / store.SHARD_ARCHIVE
    if os.path.lexists(directory / store.LEGACY_SHARD_ARCHIVE):
        raise ValueError("Restore the receipt-preserving archive before replacing its generation")
    current_manifest, current_receipt = map(_file_bytes, (manifest_path, receipt_path))
    if saved.exists() or saved.is_symlink():
        recovery = json.loads(_file_bytes(saved))
        original_manifest = base64.b64decode(recovery["manifest"], validate=True)
        original_receipt = base64.b64decode(recovery["receipt"], validate=True)
        # An interrupted conversion may already have replaced one of the two.
        expected_archive = recovery["archive_sha256"]
    else:
        original_manifest, original_receipt = current_manifest, current_receipt
        expected_archive = None
    manifest = json.loads(original_manifest)
    store._validate_manifest(manifest)
    if "shard_archive" in manifest or store.SHARD_ARCHIVE in manifest["files"]:
        raise ValueError("Expected a loose model store")
    receipt = json.loads(original_receipt)
    prefix = "model_store/"
    if (receipt["files"].get(prefix + "manifest.json") != _sha(original_manifest)
            or any(receipt["files"].get(prefix + name) != sha for name, sha in manifest["files"].items())):
        raise ValueError("The loose store does not match its original Rescue receipt")
    shards = manifest["shards"]
    files = manifest["files"]
    if not shards or any("/" in item["path"] or not item["path"].endswith((".jsonl.gz", ".jsonl.zst"))
                         for item in shards):
        raise ValueError("Expected named compressed JSONL shards in the model_store directory")

    identities = {}
    if not archive.exists() and not archive.is_symlink():
        if expected_archive is not None:
            raise ValueError("The pending migration archive is missing")
        identities = _build_archive(directory, archive, shards, files)
    archive_path = store._member(directory, archive.name)
    archive_identity = store._stat(archive_path)
    archive_sha = store._digest(archive_path)
    if expected_archive is not None and expected_archive != archive_sha:
        raise ValueError("The pending migration archive changed")
    if any(files[name] != digest for name, digest in store._archive_shard_hashes(directory, shards).items()):
        raise ValueError("A ZIP member differs from the original Rescue shard")
    if store._stat(archive_path) != archive_identity:
        raise ValueError("The archive changed during verification")
    manifest["shard_archive"] = archive.name
    manifest["files"][archive.name] = archive_sha
    packed_manifest = store._json(manifest) + b"\n"
    for shard in shards:
        receipt["files"].pop(prefix + shard["path"])
    receipt["files"][prefix + "manifest.json"] = _sha(packed_manifest)
    receipt["files"][prefix + archive.name] = archive_sha
    packed_receipt = store._json(receipt) + b"\n"
    if current_manifest not in (original_manifest, packed_manifest) or current_receipt not in (original_receipt, packed_receipt):
        raise ValueError("Model-store metadata changed since the migration started")
    if expected_archive is None:
        _replace(saved, store._json({"manifest": base64.b64encode(original_manifest).decode(),
                                    "receipt": base64.b64encode(original_receipt).decode(),
                                    "archive_sha256": archive_sha}) + b"\n")
    for path, original, replacement in (
            (manifest_path, current_manifest, packed_manifest), (receipt_path, current_receipt, packed_receipt)):
        if _file_bytes(path) != original:
            raise ValueError(f"Metadata changed during migration: {path}")
        if original != replacement:
            _replace(path, replacement)
    store.frozen_model_store_key(root)

    removed = 0
    for shard in shards:
        name = shard["path"]
        path = directory / name
        if not path.exists() and not path.is_symlink():
            continue
        path = store._member(directory, name)
        before = store._stat(path)
        if identities.get(name) != before and store._digest(path) != files[name]:
            raise ValueError(f"Loose shard changed before retirement: {name}")
        if store._stat(path) != before or store._stat(archive_path) != archive_identity:
            raise ValueError(f"Source or archive changed before retirement: {name}")
        path.unlink()
        removed += 1
    store._sync_directory(directory)
    return {"worker": str(root), "shards": len(shards), "removed": removed, "archive_sha256": archive_sha,
            "original_receipt_sha256": _sha(original_receipt), "receipt_sha256": _sha(packed_receipt)}


def move_receipted_shards(worker, expected_receipt, *, restore=False):
    """Change storage for a pinned worker without changing its producer identity."""
    with _locked_store(worker) as (root, directory):
        receipt, receipt_sha = store._receipt(root)
        manifest, manifest_sha = store._manifest(directory)
        shards, files = manifest["shards"], manifest["files"]
        if (receipt_sha != expected_receipt or "shard_archive" in manifest or not shards
                or receipt["files"].get("model_store/manifest.json") != manifest_sha
                or any(receipt["files"].get(f"model_store/{name}") != sha for name, sha in files.items())
                or any("/" in item["path"] or not item["path"].endswith((".jsonl.gz", ".jsonl.zst"))
                       for item in shards)):
            raise ValueError("Expected the exact receipted loose-shard publication")
        archive = store._legacy_archive(directory, manifest)
        identities = {}
        if archive is None:
            if restore:
                store.verify_model_store_key(root, store.frozen_model_store_key(root))
                return {"worker": str(root), "restored": 0, "receipt_sha256": receipt_sha}
            archive = directory / store.LEGACY_SHARD_ARCHIVE
            identities = _build_archive(directory, archive, shards, files)
        archive_identity = store._stat(archive)

        with store.verified_legacy_shards(root, receipt["files"]), \
                store._open_shards(directory, manifest, legacy_archive=archive) as zipped:
            changed = 0
            for shard in shards:
                name = shard["path"]
                path = directory / name
                if store._stat(store._member(directory, archive.name)) != archive_identity:
                    raise ValueError("Migration ZIP changed")
                if os.path.lexists(path):
                    path = store._member(directory, name)
                    current = store._stat(path)
                    if ((identities.get(name) != current and store._digest(path) != files[name])
                            or store._stat(path) != current):
                        raise ValueError(f"Loose shard changed: {name}")
                    if not restore:
                        path.unlink()
                        changed += 1
                elif restore:
                    with tempfile.NamedTemporaryFile(dir=directory) as handle, (zipped / name).open("rb") as source:
                        shutil.copyfileobj(source, handle, length=1024 * 1024)
                        handle.flush()
                        os.fsync(handle.fileno())
                        if store._digest(handle.name) != files[name]:
                            raise ValueError(f"Restored shard differs: {name}")
                        os.link(handle.name, path)
                    changed += 1
                elif name in identities:
                    raise ValueError(f"Loose shard disappeared during migration: {name}")
            store._sync_directory(directory)
        if restore:
            archive.unlink()
            store._sync_directory(directory)
        return {"worker": str(root), "shards": len(shards), "restored" if restore else "removed": changed,
                "receipt_sha256": receipt_sha, "manifest_sha256": manifest_sha}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("worker", type=Path)
    parser.add_argument("--writers-and-receipt-consumers-stopped", action="store_true", required=True)
    parser.add_argument("--preserve-receipt", metavar="RECEIPT_SHA256")
    parser.add_argument("--restore-loose-shards", action="store_true")
    args = parser.parse_args()
    if args.restore_loose_shards and args.preserve_receipt is None:
        parser.error("--restore-loose-shards requires --preserve-receipt")
    result = (archive_closed_store(args.worker) if args.preserve_receipt is None
              else move_receipted_shards(args.worker, args.preserve_receipt, restore=args.restore_loose_shards))
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

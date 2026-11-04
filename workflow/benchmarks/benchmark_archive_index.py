#!/usr/bin/env python3
"""Synthetic native archive-index publication; run in a GeneGalleon container.

Only index publication is timed. Valid ZIPs and deterministic logical payloads
are prepared beforehand and all payloads are verified afterwards. Use separate
processes/fixtures for repetitions and --module-root for the baseline helper.
"""
import argparse
import hashlib
import importlib
import json
import resource
import sys
import tempfile
import time
import zipfile
import zlib
from dataclasses import asdict
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--module-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
parser.add_argument('--families', type=int, default=2500)
parser.add_argument('--subdirs', type=int, default=29)
parser.add_argument('--changed-families', type=int, default=100)
parser.add_argument('--work-dir', type=Path)
args = parser.parse_args()
if not 1 <= args.changed_families <= args.families or args.subdirs < 1:
    parser.error('Require positive subdirs and 1 <= changed-families <= families')
module_path = args.module_root / 'gene_family_output_store.py'
if not module_path.is_file():
    parser.error(f'Module root has no gene_family_output_store.py: {args.module_root}')
sys.path[:0] = [str(args.module_root), str(Path(__file__).resolve().parents[1] / 'support')]
m = importlib.import_module('gene_family_output_store')
if Path(m.__file__).resolve() != module_path.resolve():
    parser.error('Imported archive helper does not match --module-root')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

with tempfile.TemporaryDirectory(prefix='gg-archive-index-', dir=args.work_dir) as temporary:
    root = Path(temporary) / 'orthogroup'
    root.mkdir()
    store = m.GeneFamilyOutputStore(root)
    artifacts = {}
    zip_paths = []
    expected = hashlib.sha256()
    for j in range(args.subdirs):
        subdir = f'outputs{j:02d}'
        directory = store.payload_root / subdir
        directory.mkdir(parents=True)
        for start in range(0, args.families, 100):
            path = directory / f'{start // 100:08d}.zip'
            zip_paths.append(path)
            members = []
            with zipfile.ZipFile(path, 'w', compression=zipfile.ZIP_STORED) as archive:
                for i in range(start, min(start + 100, args.families)):
                    family = f'OG{i:07d}'
                    logical = f'{subdir}/{family}_data.tsv'
                    data = hashlib.sha256(logical.encode()).digest() * 2
                    digest = hashlib.sha256(data).hexdigest()
                    archive.writestr(logical, data)
                    artifact = m.Artifact(logical_path=logical, subdir=subdir, name=logical.split('/', 1)[1],
                                          generation=1, zip_path=path, member_name=logical, size=len(data),
                                          crc=zlib.crc32(data), sha256=digest, mtime_ns=1700000000000000000,
                                          family_id=family)
                    artifacts[logical] = artifact
                    member = asdict(artifact)
                    member.pop('zip_path')
                    member.pop('name')
                    member.pop('subdir')
                    member.pop('generation')
                    members.append(member)
                    expected.update(logical.encode() + b'\0' + digest.encode() + b'\n')
                archive.writestr(m.MANIFEST_MEMBER, json.dumps(dict(schema_version=1, generation=1,
                                 mode='orthogroup', subdir=subdir, compression='store', compression_level=0,
                                 created_utc='2026-01-01T00:00:00Z', members=members), sort_keys=True))
    before_zip = {str(path.relative_to(root)): sha(path) for path in zip_paths}
    with m.producer_quiescence_lock(store.archive_root), m.archive_lock(store.archive_root):
        store._write_index(artifacts)
    updated = {subdir: {key: value for key, value in artifacts.items() if value.subdir == subdir}
               for subdir in sorted({artifact.subdir for artifact in artifacts.values()})}
    changed = {key for key, value in artifacts.items() if int(value.family_id[2:]) < args.changed_families}
    started = time.monotonic()
    with m.producer_quiescence_lock(store.archive_root), m.archive_lock(store.archive_root):
        with store._index_update():
            store._merge_index_subdirs_uncommitted(updated, changed_paths=changed)
    elapsed = time.monotonic() - started
    paths = (list((store.archive_root / m.INDEX_DIR_NAME).glob('*.json')) +
             list((store.archive_root / m.SUBDIR_INDEX_DIR_NAME).glob('*.json')) +
             [store.archive_root / m.INDEX_CATALOG_FILE])
    manifest = {str(path.relative_to(store.archive_root)): sha(path) for path in paths}
    observed = hashlib.sha256()
    with store.read_snapshot():
        for logical in sorted(artifacts):
            subdir, name = logical.split('/', 1)
            with store.open_binary(subdir, name) as handle:
                data = handle.read()
            digest = hashlib.sha256(data).hexdigest()
            assert digest == artifacts[logical].sha256
            observed.update(logical.encode() + b'\0' + digest.encode() + b'\n')
    assert observed.hexdigest() == expected.hexdigest()
    assert before_zip == {str(path.relative_to(root)): sha(path) for path in zip_paths}
    assert not (store.archive_root / m.INDEX_UPDATE_FILE).exists()
    print(json.dumps(dict(kind='synthetic_native_archive_index_commit', families=args.families,
                          subdirs=args.subdirs, indexed_members=len(artifacts), changed_members=len(changed),
                          elapsed_seconds=elapsed, peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                          module_sha256=sha(Path(m.__file__)), index_files=len(manifest),
                          index_manifest_sha256=hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest(),
                          logical_payload_sha256=observed.hexdigest(), zip_bodies_unchanged=True), sort_keys=True))

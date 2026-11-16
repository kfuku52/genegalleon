#!/usr/bin/env python3
"""Bundle diagnostic files and record checksums of separately published results."""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import shutil
import tempfile
import zipfile
from pathlib import Path, PurePosixPath

try:
    from atomic_zip_publish import publish_zip
    from input_generation_array_state import FreshDigestBatch
    from safe_zip_extract import validated_members
except ImportError:
    from .atomic_zip_publish import publish_zip
    from .input_generation_array_state import FreshDigestBatch
    from .safe_zip_extract import validated_members

REFERENCES = 'external_members.json'


def relative_name(value):
    path = PurePosixPath(value)
    if (not path.parts or path.is_absolute() or '..' in path.parts
            or str(path) != value or '\\' in value):
        raise ValueError('Unsafe archive path: ' + value)
    return path


def selected_paths(root, names):
    for name in names:
        path = root / relative_name(name)
        if not path.is_file() or path.is_symlink() or path.resolve() != path:
            raise ValueError('Expected a regular diagnostic file: ' + str(path))
        yield path


def pack(root, destination, *, members=None, external=(), inputs=()):
    root, destination = Path(root).resolve(strict=True), Path(destination).absolute()
    if not root.is_dir() or destination.resolve().is_relative_to(root):
        raise ValueError('Diagnostic archive must be outside its source directory')
    names = sorted(members if members is not None else (
        path.relative_to(root).as_posix() for path in root.rglob('*') if path.is_symlink() or path.is_file()))
    refs = dict(external)
    additional = {str(relative_name(name)): Path(path).resolve(strict=True) for name, path in inputs}
    if (len(set(names)) != len(names) or len(refs) != len(external) or refs.keys() - set(names)
            or len(additional) != len(inputs) or set(names) & additional.keys()
            or REFERENCES in set(names) | additional.keys()):
        raise ValueError('Invalid diagnostic archive member inventory')
    for location in refs.values():
        if len(relative_name(location).parts) != 2:
            raise ValueError('External diagnostic members require a logical subdirectory and filename')
    paths = list(selected_paths(root, names))
    batch = FreshDigestBatch()
    hashes = batch.read([*(root / name for name in refs), *additional.values()])
    excluded = [{'member': root.name + '/' + name, 'logical_path': location,
                 'sha256': hashes[str(root / name)]} for name, location in sorted(refs.items())]
    with tempfile.TemporaryDirectory(prefix='gg-audit-') as work:
        path = Path(work) / 'audit.zip'
        with zipfile.ZipFile(path, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
            for source in paths:
                name = source.relative_to(root).as_posix()
                if name not in refs:
                    archive.write(source, root.name + '/' + name)
            for name, source in sorted(additional.items()):
                archive.write(source, root.name + '/' + name)
            if excluded:
                archive.writestr(root.name + '/' + REFERENCES, json.dumps({'schema': 1, 'members': excluded}) + '\n')
        batch.check()
        return publish_zip(path, destination, expected_prefix=root.name)


def materialize(source, logical_root, destination):
    """Create a self-contained ZIP on demand, using a raw or archived family store."""
    try:
        from gene_family_output_store import GeneFamilyOutputStore
    except ImportError:
        from .gene_family_output_store import GeneFamilyOutputStore
    source, destination = Path(source).resolve(strict=True), Path(destination).absolute()
    if destination.exists():
        raise FileExistsError(destination)
    with tempfile.TemporaryDirectory(prefix='gg-complete-audit-') as work:
        path = Path(work) / 'audit.zip'
        with zipfile.ZipFile(source) as archive, zipfile.ZipFile(
                path, 'w', zipfile.ZIP_DEFLATED) as complete:
            roots = {relative_name(name.rstrip('/')).parts[0] for name in archive.namelist()}
            if len(roots) != 1:
                raise ValueError('Expected one diagnostic archive root')
            prefix, = roots
            validated_members(archive, source, prefix)
            ref_name = prefix + '/' + REFERENCES
            refs = json.loads(archive.read(ref_name)) if ref_name in archive.namelist() else {'schema': 1, 'members': []}
            if refs['schema'] != 1:
                raise ValueError('Unsupported external diagnostic references')
            for info in archive.infolist():
                if not info.is_dir() and info.filename != ref_name:
                    with archive.open(info) as reader, complete.open(info.filename, 'w', force_zip64=True) as writer:
                        shutil.copyfileobj(reader, writer)
            names = set(archive.namelist())
            store = GeneFamilyOutputStore(Path(logical_root)) if refs['members'] else None
            with store.read_snapshot() if store is not None else contextlib.nullcontext():
                for row in refs['members']:
                    member = relative_name(row['member'])
                    if member.parts[0] != prefix or str(member) in names:
                        raise ValueError('Unexpected external diagnostic member')
                    names.add(str(member))
                    subdir, name = relative_name(row['logical_path']).parts
                    checksum = hashlib.sha256()
                    with store.open_binary(subdir, name) as reader, complete.open(
                            str(member), 'w', force_zip64=True) as writer:
                        for chunk in iter(lambda: reader.read(1024 * 1024), b''):
                            checksum.update(chunk)
                            writer.write(chunk)
                    if checksum.hexdigest() != row['sha256']:
                        raise ValueError('Published output differs from audit: ' + row['logical_path'])
        return publish_zip(path, destination, expected_prefix=prefix)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    create = commands.add_parser('pack')
    create.add_argument('--root', type=Path, required=True)
    create.add_argument('--destination', type=Path, required=True)
    create.add_argument('--member', action='append')
    create.add_argument('--external', nargs=2, action='append', default=[])
    create.add_argument('--input', nargs=2, action='append', default=[], metavar=('MEMBER', 'PATH'))
    extract = commands.add_parser('materialize')
    extract.add_argument('--source', type=Path, required=True)
    extract.add_argument('--logical-root', type=Path, required=True)
    extract.add_argument('--destination', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'pack':
        pack(args.root, args.destination, members=args.member, external=args.external, inputs=args.input)
    else:
        materialize(args.source, args.logical_root, args.destination)


if __name__ == '__main__':
    main()

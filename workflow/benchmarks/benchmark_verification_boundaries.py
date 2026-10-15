#!/usr/bin/env python3
"""Measure runtime, stage-input and directory verification in fresh processes."""
import argparse
import contextlib
import hashlib
import io
import json
import os
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORKLOADS = ('runtime', 'index', 'synteny-plan', 'cds-view', 'logical-directory', 'busco-cache')


def worker(args):
    sys.path.insert(0, str(args.support))
    import annevo_runtime_cache as runtime
    import gene_model_catalog as catalog
    import input_generation_array_state as state

    counters = dict(sha256_reads=0, sha256_bytes=0, directory_scans=0)
    real_sha, real_digest, real_signature, real_scandir = runtime.file_sha, state.digest, catalog.source_signature, os.scandir

    def hashed(function, path):
        counters['sha256_reads'] += 1
        counters['sha256_bytes'] += Path(path).stat().st_size
        return function(path)

    runtime.file_sha = lambda path: hashed(real_sha, path)
    state.digest = lambda path: hashed(real_digest, path)
    catalog.source_signature = lambda path: hashed(real_signature, path)
    if args.workload == 'synteny-plan':
        import pairwise_synteny as synteny
        real_synteny_digest = synteny.digest
        synteny.digest = lambda path: hashed(real_synteny_digest, path)
    elif args.workload == 'cds-view':
        import cds_resolution as resolution
        real_resolution_digest = resolution.digest
        resolution.digest = lambda path: hashed(real_resolution_digest, path)
    elif args.workload == 'logical-directory':
        import artifact_provenance as provenance
        real_provenance_digest = provenance.cached_sha256_file
        provenance.cached_sha256_file = lambda path: hashed(real_provenance_digest, path)
        provenance.configure_digest_cache(None)
    elif args.workload == 'busco-cache':
        import inspect

        import busco_guide_tree as guide
        guide.digest = state.digest
        proof_supported = 'proof' in inspect.signature(guide.cache_matches).parameters

    def scandir(*a, **kw):
        counters['directory_scans'] += 1
        return real_scandir(*a, **kw)

    os.scandir = scandir
    with tempfile.TemporaryDirectory(prefix='gg-verification-') as temporary:
        root = Path(temporary)
        block = b'ACGT' * (1024 * 256)
        if args.workload in ('runtime', 'index'):
            model = root / 'model.pt'
            with model.open('wb') as stream:
                for _ in range(args.mib):
                    stream.write(block)
        if args.workload in ('runtime', 'logical-directory'):
            tree = root / 'runtime'
            tree.mkdir()
            for index in range(args.files):
                directory = tree / f'lib{index % 128:03d}' / 'package'
                directory.mkdir(parents=True, exist_ok=True)
                (directory / f'module{index:05d}.py').write_text('runtime fixture\n')
            if args.workload == 'runtime':
                (tree / 'module-alias').symlink_to('lib000/package/module00000.py')
                recorded = runtime.inspect_tree(tree)
        elif args.workload == 'index':
            genome = root / 'genome.fa'
            with genome.open('wb') as stream:
                stream.write(b'>chr\n')
                for _ in range(args.mib):
                    stream.write(block + b'\n')
            os.environ['GG_GENOME_INDEX_CACHE'] = str(root / 'index-cache')
            with contextlib.redirect_stdout(io.StringIO()), catalog.indexed_genome(genome):
                pass
        elif args.workload == 'synteny-plan':
            import csv
            sources = []
            for index in range(4):
                name = f'Species_{index}'
                fasta, gff = root / (name + '.faa'), root / (name + '.gff3')
                with fasta.open('wb') as stream:
                    stream.write(('>' + name + '_g\n').encode())
                    for _ in range(args.mib):
                        stream.write(block + b'\n')
                gff.write_text('chr\ts\tgene\t1\t9\t.\t+\t.\tID=g\n')
                sources.append((name, fasta.name, gff.name))
            pairs = root / 'pairs.tsv'
            with pairs.open('w') as stream:
                writer = csv.writer(stream, delimiter='\t')
                writer.writerow(['analysis_id', 'target_species', 'target_fasta', 'target_gff',
                                 'query_species', 'query_fasta', 'query_gff'])
                for index in range(args.boundaries * 2):
                    writer.writerow([f'pair{index}', *sources[index % 4], *sources[(index + 1) % 4]])
            planning = argparse.Namespace(workspace=root, pairs=pairs, sequence_mode='protein', genetic_code=1,
                cscore=.7, min_anchors=4, distance=20, minimum_mapping_fraction=1,
                formats='svg', karyotype_sort='both_length', dotplot_min_length=0)
            tools = synteny.tool_identity()
            synteny.tool_identity = lambda: tools
        elif args.workload == 'cds-view':
            sources, resolutions, views = root / 'sources', root / 'resolutions', root / 'views'
            sources.mkdir()
            resolutions.mkdir()
            view_inputs = {}
            for index in range(2):
                source = sources / f'Species_{index}.fa'
                with source.open('wb') as stream:
                    stream.write(b'>gene\n')
                    for _ in range(args.mib):
                        stream.write(block + b'\n')
                target = resolutions / source.name
                import shutil
                shutil.copyfile(source, target)
                sha = real_resolution_digest(source)
                meta = {'contract': {'sources': {'cds': {'path': str(source), 'sha256': sha}}},
                        'output_sha256': sha}
                (resolutions / (source.name + '.resolution.meta.json')).write_text(json.dumps(meta))
                view_inputs[source.name] = {'path': str(target), 'sha256': sha}
            view_identity = hashlib.sha256(json.dumps(view_inputs, sort_keys=True).encode()).hexdigest()
            resolution.resolved_view(sources, resolutions, views)
        elif args.workload == 'busco-cache':
            sketches = []
            for index in range(args.boundaries):
                payload, metadata = root / f'{index}.bin', root / f'{index}.json'
                with payload.open('wb') as stream:
                    for _ in range(args.mib):
                        stream.write(block)
                metadata.write_text(json.dumps({'sha256': real_digest(payload)}))
                sketches.append((payload, metadata))
        counters.update({key: 0 for key in counters})
        started = time.perf_counter()
        logical = []
        with contextlib.redirect_stdout(io.StringIO()):
            if args.workload == 'runtime':
                verification = runtime.RuntimeVerification() if hasattr(runtime, 'RuntimeVerification') else None
                for _ in range(args.boundaries):
                    if verification is None:
                        runtime.verify_tree(tree, recorded)
                        logical.append(runtime.file_sha(model))
                    else:
                        verification.verify_tree(tree, recorded)
                        logical.append(verification.file_sha(model))
                if verification is not None:
                    verification.check()
                logical.append({name: {k: v for k, v in item.items() if k != 'signature'}
                                for name, item in recorded.items()})
            elif args.workload == 'index':
                for _ in range(args.boundaries):
                    with catalog.indexed_genome(genome) as indexed:
                        logical.append(indexed.fetch('chr', 13, 113))
                        receipt = json.loads((Path(os.fsdecode(indexed.filename)).parent / 'receipt.json').read_text())
                        receipt['key']['source']['path'] = '<genome>'
                        logical.append(receipt)
            elif args.workload == 'synteny-plan':
                plan = synteny.build_plan(planning)
                # GG implementation identities change with the code, while all
                # scientific sources, parameters, pairs and tool identities must match.
                plan['input_hashes'] = {key: value for key, value in plan['input_hashes'].items()
                                        if not Path(key).is_relative_to(args.support)}
                logical.append(json.loads(json.dumps(plan).replace(str(root), '<workspace>')))
            elif args.workload == 'cds-view':
                for _ in range(args.boundaries):
                    view = resolution.resolved_view(sources, resolutions, views)
                    assert view.name == view_identity
                    logical.append(sorted(path.name for path in view.iterdir()))
            elif args.workload == 'logical-directory':
                for _ in range(args.boundaries):
                    with provenance.digest_observation():
                        logical.append(provenance.raw_or_zip_directory_digest(tree))
            else:
                paths, expected = [], []
                for payload, metadata in sketches:
                    proof = state.FreshDigestBatch()
                    if proof_supported:
                        assert guide.cache_matches(payload, metadata, proof=proof)
                        observed = proof.read([payload, metadata])[str(payload)]
                    else:
                        assert guide.cache_matches(payload, metadata)
                        observed = guide.digest(payload)
                    assert observed == json.loads(metadata.read_text())['sha256']
                    proof.check()
                    paths.append(payload)
                    expected.append(observed)
                fresh = state.FreshDigestBatch()
                actual = fresh.read(paths)
                assert [actual[str(path)] for path in paths] == expected
                fresh.check()
                logical.append(expected)
        elapsed = time.perf_counter() - started
        if args.workload == 'cds-view':
            logical.append({path.name: real_resolution_digest(path) for path in view.iterdir()})
        print(json.dumps(dict(seconds=elapsed, peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                              logical_sha256=hashlib.sha256(json.dumps(logical, sort_keys=True).encode()).hexdigest(),
                              **counters)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-support', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mib', type=int, default=256)
    parser.add_argument('--files', type=int, default=4096)
    parser.add_argument('--boundaries', type=int, default=7)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--cases', choices=WORKLOADS, nargs='+', default=['runtime', 'index'])
    parser.add_argument('--support', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--workload', choices=WORKLOADS, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if min(args.mib, args.files, args.boundaries, args.repeats) < 1:
        parser.error('workload sizes and repeats must be positive')
    if args.support:
        worker(args)
        return
    supports = {'current': ROOT / 'workflow/support'}
    if args.baseline_support:
        supports = {'baseline': args.baseline_support, **supports}
    result = dict(environment=dict(platform=platform.platform(), python=platform.python_version()),
                  configuration=dict(mib=args.mib, files=args.files, boundaries=args.boundaries,
                                     repeats=args.repeats, warmups=1), workloads={},
                  source_sha256={label: {name: hashlib.sha256((support / name).read_bytes()).hexdigest()
                                         for name in ('annevo_runtime_cache.py', 'annevo_runtime.py', 'annevo_pipeline.py',
                                                      'oriongeno_runtime.py', 'gene_model_catalog.py',
                                                      'pairwise_synteny.py', 'cds_resolution.py', 'artifact_provenance.py',
                                                      'busco_guide_tree.py')}
                                 for label, support in supports.items()})
    for workload in args.cases:
        rows = {label: [] for label in supports}
        for trial in range(args.repeats + 1):
            items = list(supports.items())
            if trial % 2:
                items.reverse()
            for label, support in items:
                command = [sys.executable, str(Path(__file__).resolve()), '--support', str(support),
                           '--workload', workload, '--output', str(args.output), '--mib', str(args.mib),
                           '--files', str(args.files), '--boundaries', str(args.boundaries)]
                completed = subprocess.run(command, check=True, capture_output=True, text=True)
                row = json.loads(completed.stdout)
                if trial:
                    rows[label].append(row)
        fingerprints = {row['logical_sha256'] for values in rows.values() for row in values}
        if len(fingerprints) != 1:
            raise ValueError('Verification outputs differ: ' + workload)
        result['workloads'][workload] = {label: dict(trials=values, medians={key: statistics.median(row[key] for row in values)
                                             for key in values[0] if key != 'logical_sha256'})
                                           for label, values in rows.items()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

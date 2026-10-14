#!/usr/bin/env python3
"""One allocated GPU producer with bounded, independent native CPU decoders.

Invoke ANNEVO's unmodified prediction.py and decoding.py instead of serializing
both through annotation.py. Complete contigs, native thresholds and receipts
retain their meaning. A completed source-bound GFF is reusable by rescue.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import annevo_runtime as native
from annevo_runtime_cache import runtime_environment
from annevo_runtime_download import atomic_json, file_sha
from genome_sequence_identity import genome_sequence_digest
from shared_namespace_lock import namespace_lock


def read_inputs(path, output):
    path, output = Path(path).absolute(), Path(output).absolute()
    rows, names, protected = [], set(), {path}
    before = file_sha(path)
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if reader.fieldnames not in (['species', 'genome', 'receipt'],
                                    ['species', 'genome', 'receipt', 'selection_receipt']):
            raise ValueError('Pipeline TSV requires species, genome, receipt and optional selection_receipt')
        for row in reader:
            name = row['species']
            if (None in row or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', name or '')
                    or name in names or any(not row.get(k) for k in ('genome', 'receipt'))):
                raise ValueError('Invalid, incomplete or duplicate pipeline species')
            names.add(name)
            for key in ('genome', 'receipt', 'selection_receipt'):
                if row.get(key):
                    p = Path(row[key])
                    p = (p if p.is_absolute() else path.parent / p).absolute()
                    if p.is_symlink() or not p.is_file():
                        raise ValueError('Pipeline input is missing or redirected: ' + str(p))
                    row[key] = str(p)
                    protected.add(p.resolve())
            rows.append(row)
    if not rows or file_sha(path) != before:
        raise ValueError('Pipeline inputs are empty or changed during reading')
    # Publication and intermediate directories must not contain an input.
    if any(p == output.resolve() or p.is_relative_to(output.resolve()) for p in protected):
        raise ValueError('Pipeline output must be separate from its inputs')
    return rows, before


def genome_bases(path):
    with Path(path).open('rb') as stream:
        compressed = stream.read(2) == b'\x1f\x8b'
    opener = gzip.open if compressed else open
    with opener(path, 'rt') as stream:
        return sum(len(line.strip()) for line in stream if not line.startswith('>'))


def parameters(args):
    return dict(threads=args.decode_threads, batch_size=args.batch_size,
                sequence_threshold=args.sequence_threshold,
                num_workers=min(2, args.predict_threads), min_intron_length=20, min_prot_length=100)


def existing_pair(output, receipt, expected):
    if output.is_symlink() or receipt.is_symlink():
        raise ValueError('Redirected pipeline prediction pair')
    if not output.exists() and not receipt.exists():
        return None
    if output.is_symlink() or receipt.is_symlink() or not output.is_file() or not receipt.is_file():
        raise ValueError('Incomplete or redirected pipeline prediction pair')
    value = json.loads(receipt.read_text())
    if (not isinstance(value, dict) or any(value.get(k) != v for k, v in expected.items())
            or value.get('kind') != 'genegalleon-annevo-native-prediction'
            or value.get('native_returncode') != 0 or value.get('gpu_probe', {}).get('gpu_verified') is not True
            or value.get('GFF_sha256') != file_sha(output) or value.get('GFF_bytes') != output.stat().st_size):
        raise ValueError('Stale pipeline prediction; use a new output directory')
    native.gff_counts(output)
    return value


def command(receipt, script, arguments):
    return [receipt['python'], '-E', '-s', '-B',
            str(Path(receipt['native_entrypoint']).with_name(script)), *arguments]


def execute(argv, env, work, logs, stage):
    start = time.monotonic()
    with (logs / (stage + '.stdout.txt')).open('w') as stdout, (logs / (stage + '.stderr.txt')).open('w') as stderr:
        result = subprocess.run(argv, env=env, cwd=work, stdout=stdout, stderr=stderr)
    if result.returncode:
        raise RuntimeError('ANNEVO native ' + stage + ' failed; retained logs: ' + str(logs))
    return dict(stage=stage, argv=argv, seconds=time.monotonic() - start,
                started_monotonic=start, finished_monotonic=time.monotonic(), returncode=result.returncode)


def run(args):
    if not math.isfinite(args.max_tmp_gib) or min(args.predict_threads, args.decode_threads, args.decoder_workers,
           args.batch_size, args.sequence_threshold, args.max_tmp_gib) <= 0:
        raise ValueError('Pipeline thread, worker, batch and storage bounds must be positive')
    cpu_limit = args.predict_threads + min(2, args.predict_threads) + args.decoder_workers * args.decode_threads
    allocation = os.environ.get('SLURM_CPUS_PER_TASK')
    if allocation and cpu_limit > int(allocation):
        raise ValueError('Pipeline CPU budget exceeds the scheduler allocation')
    rows, table_sha = read_inputs(args.inputs, args.output)
    root = Path(args.output).absolute()
    root.mkdir(mode=0o700, parents=True, exist_ok=True)
    tmp_root = Path(args.tmp_root).absolute()
    if not tmp_root.is_dir() or tmp_root.is_symlink():
        raise ValueError('Pipeline temporary root is missing or redirected')
    budget = int(args.max_tmp_gib * 1024**3)
    completed, pending, reserved = {}, [], 0
    started = time.monotonic()
    runtime_receipts = {}
    stats = dict(gpu_workers=1, decoder_workers=args.decoder_workers,
                 predict_threads=args.predict_threads, decode_threads=args.decode_threads,
                 cpu_budget=cpu_limit, max_tmp_bytes=budget, peak_reserved_bytes=0,
                 max_inflight_species=0, input_sha256=table_sha)

    def reap(wait=False):
        nonlocal reserved
        if wait and pending:
            pending[0][0].result()
        for item in list(pending):
            future, size = item
            if future.done():
                name, value = future.result()
                completed[name] = value
                reserved -= size
                pending.remove(item)

    def process(row, executor):
        nonlocal reserved
        name = row['species']
        directory = root / name
        if directory.is_symlink():
            raise ValueError('Redirected species output')
        directory.mkdir(mode=0o700, exist_ok=True)
        output, receipt_path = directory / (name + '.gff3'), directory / 'native_prediction.json'
        selected, genome = Path(row['receipt']), Path(row['genome'])
        runtime = native.read_selected(selected)
        if root.resolve().is_relative_to(Path(runtime['generation']).resolve()):
            raise ValueError('Pipeline must not publish inside its runtime')
        runtime_receipts[str(selected)] = file_sha(selected)
        profile = runtime['model_profile']
        selection_args = argparse.Namespace(model_profile=profile, species=name,
                                            selection_receipt=row.get('selection_receipt', ''))
        selection, selection_sha = native.selection_record(selection_args, genome)
        original_sha = file_sha(genome)
        settings = native.inference_settings(runtime, profile)
        expected = dict(runtime_digest=runtime['runtime_digest'], cache_key=runtime['cache_key'],
                        source_commit=runtime['source_commit'], model_sha256=runtime['model_sha256'],
                        model_profile=profile, genome_sha256=original_sha,
                        prepared_receipt_sha256=file_sha(selected), parameters=parameters(args),
                        inference_settings=settings)
        if selection is not None:
            expected.update(model_selection=selection, model_selection_sha256=selection_sha, species=name)
        # Keep the same publication lock throughout native prediction and decode.
        lock = namespace_lock(output.with_name(output.name + '.publication.lock'), exclusive=True, timeout=86400)
        lock.__enter__()
        temporary = None
        try:
            old = existing_pair(output, receipt_path, expected)
            if old:
                completed[name] = dict(reused=True, gff=str(output), receipt=str(receipt_path))
                lock.__exit__(None, None, None)
                return
            sequence_identity = genome_sequence_digest(genome)
            bp = genome_bases(genome)
            if file_sha(genome) != original_sha:
                raise RuntimeError('Pipeline genome changed during identity reading')
            # Native float16 probabilities: 15 classes, two strands = 60 B/base.
            # Include uncompressed FASTA, HDF5 metadata and ten percent headroom.
            required = int(bp * 68 + 16 * 1024**2)
            if required > budget:
                raise ValueError('Genome intermediate exceeds pipeline storage budget: ' + name)
            while pending and (len(pending) >= args.decoder_workers + 1 or reserved + required > budget):
                reap(wait=True)
            if shutil.disk_usage(tmp_root).free < required:
                raise ValueError('Insufficient local temporary space for: ' + name)
            reserved += required
            stats['peak_reserved_bytes'] = max(stats['peak_reserved_bytes'], reserved)
            stats['max_inflight_species'] = max(stats['max_inflight_species'], len(pending) + 1)
            temporary = tempfile.TemporaryDirectory(prefix='gg-ap-', dir=tmp_root)
            work = Path(temporary.name)
            native_genome = genome
            compressed = genome.name.endswith('.gz')
            if compressed:
                native_genome = work / 'genome.fa'
                with gzip.open(genome, 'rb') as source, native_genome.open('xb') as target:
                    shutil.copyfileobj(source, target, length=8 * 1024 * 1024)
            prepared_sha = file_sha(native_genome)
            logs = directory / 'logs'
            logs.mkdir(exist_ok=True)
            env = runtime_environment(runtime)
            env.update(OMP_NUM_THREADS=str(args.predict_threads), MKL_NUM_THREADS=str(args.predict_threads),
                       OPENBLAS_NUM_THREADS=str(args.predict_threads), TMPDIR=str(work),
                       NUMBA_CACHE_DIR=str(work / 'numba-cache'), TORCH_HOME=str(work / 'torch-cache'),
                       XDG_CACHE_HOME=str(work / 'xdg-cache'))
            gpu = native.probe_gpu(runtime, env)
            predicted = work / 'probabilities.h5'
            pred_argv = command(runtime, 'prediction.py', ['-g', str(native_genome), '-m', runtime['model_path'],
                '-l', profile, '-s', str(args.sequence_threshold), '-p', str(predicted),
                '--batch_size', str(args.batch_size), '--num_workers', str(min(2, args.predict_threads))])
            if settings['overlap_pred']:
                pred_argv.append('--overlap_pred')
            phase = execute(pred_argv, env, work, logs, 'prediction')
            if not predicted.is_file() or not predicted.stat().st_size:
                raise RuntimeError('Native prediction produced no probabilities')
            if file_sha(native_genome) != prepared_sha:
                raise RuntimeError('Prepared genome changed during GPU prediction')

            def decode():
                try:
                    env_decode = {**env, 'CUDA_VISIBLE_DEVICES': '', 'OMP_NUM_THREADS': str(args.decode_threads),
                                  'MKL_NUM_THREADS': str(args.decode_threads), 'OPENBLAS_NUM_THREADS': str(args.decode_threads)}
                    gff = work / 'prediction.gff3'
                    argv = command(runtime, 'decoding.py', ['-g', str(native_genome), '-p', str(predicted),
                        '-o', str(gff), '-t', str(args.decode_threads), '--min_intron_length', '20',
                        '--min_prot_length', '100', '--show_log'])
                    decoded = execute(argv, env_decode, work, logs, 'decoding')
                    counts = native.gff_counts(gff)
                    native.read_selected(selected)
                    if file_sha(selected) != expected['prepared_receipt_sha256'] or file_sha(genome) != original_sha:
                        raise RuntimeError('ANNEVO inputs changed during pipeline execution')
                    if file_sha(native_genome) != prepared_sha:
                        raise RuntimeError('Prepared genome changed during CPU decoding')
                    if selection_sha is not None and native.selection_record(selection_args, genome)[1] != selection_sha:
                        raise RuntimeError('ANNEVO model selection changed during pipeline execution')
                    value = {**expected, 'schema_version': 1, 'kind': 'genegalleon-annevo-native-prediction',
                             'native_returncode': 0, 'gpu_probe': gpu, 'argv': pred_argv,
                             'genome_sequence_sha256': sequence_identity, 'native_genome_sha256': prepared_sha,
                             'genome_gzip': compressed, 'GFF_sha256': file_sha(gff), 'GFF_bytes': gff.stat().st_size,
                             'feature_counts': counts, 'native_wall_seconds': decoded['finished_monotonic'] - phase['started_monotonic'],
                             'pipeline': {'gpu_workers': 1, 'phases': [phase, decoded], 'predict_threads': args.predict_threads},
                             'scientific_completion_verified': False}
                    native.publish_pair(gff, output, receipt_path, value)
                    return name, dict(reused=False, gff=str(output), receipt=str(receipt_path))
                finally:
                    temporary.cleanup()
                    lock.__exit__(None, None, None)

            pending.append((executor.submit(decode), required))
        except BaseException:
            if temporary is not None:
                temporary.cleanup()
            lock.__exit__(*sys.exc_info())
            raise

    with namespace_lock(root / '.pipeline.lock', exclusive=True, timeout=86400):
        try:
            with ThreadPoolExecutor(max_workers=args.decoder_workers) as executor:
                for row in rows:
                    reap()
                    process(row, executor)
                while pending:
                    reap(wait=True)
            if file_sha(args.inputs) != table_sha:
                raise RuntimeError('Pipeline input table changed during execution')
            for selected, expected_sha in runtime_receipts.items():
                native.read_selected(selected)
                if file_sha(selected) != expected_sha:
                    raise RuntimeError('Pipeline runtime receipt changed during execution')
            result = {**stats, 'species': completed, 'wall_seconds': time.monotonic() - started,
                      'status': 'complete', 'scientific_completion_verified': False}
            atomic_json(root / 'pipeline.json', result)
            # Published only after all predictions are verified. Existing rescue
            # consumers accept exactly these three columns and native receipts.
            fd, stage = tempfile.mkstemp(prefix='.candidates-', dir=root)
            with os.fdopen(fd, 'w', newline='') as stream:
                writer = csv.writer(stream, delimiter='\t')
                writer.writerow(['species', 'gff', 'receipt'])
                for row in rows:
                    value = completed[row['species']]
                    writer.writerow([row['species'], value['gff'], value['receipt']])
            os.replace(stage, root / 'candidate_inputs.tsv')
            return result
        except BaseException as error:
            atomic_json(root / 'pipeline.failed.json', {**stats, 'status': 'failed',
                        'error': str(error), 'wall_seconds': time.monotonic() - started})
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('inputs', 'output'):
        parser.add_argument('--' + name, required=True, type=Path)
    parser.add_argument('--decoder-workers', type=int, default=2)
    parser.add_argument('--decode-threads', type=int, default=4)
    parser.add_argument('--predict-threads', type=int, default=4)
    parser.add_argument('--batch-size', type=int, default=32)
    parser.add_argument('--sequence-threshold', type=int, default=25)
    parser.add_argument('--tmp-root', type=Path, default=Path('/tmp'))
    parser.add_argument('--max-tmp-gib', type=float, default=128)
    try:
        result = run(parser.parse_args())
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print('ANNEVO pipeline: ' + str(error), file=sys.stderr)
        return 1
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

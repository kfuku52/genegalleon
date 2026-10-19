#!/usr/bin/env python3
"""Matched alternating trials for ordering, LZF and native Numba-cache reuse.

Use an existing one-GPU allocation. Native models/windows/thresholds stay fixed.
Memory peaks are sampled process-tree RSS/PSS, not whole-job cgroup memory.
"""
import argparse
import json
import os
import statistics
import sys
import tempfile
import threading
import time
from pathlib import Path

from benchmark_annevo_pipeline import GPUSamples

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
import annevo_pipeline as pipeline
import annevo_runtime as native
from annevo_prediction_equivalence import feature_digest
from annevo_runtime_cache import runtime_environment
from annevo_runtime_download import atomic_json, file_sha

CONFIGURATIONS = {
    'baseline': {'order': 'input', 'compression': 'none', 'shared_jit': False},
    'optimized': {'order': 'largest-first', 'compression': 'lzf', 'shared_jit': True},
    'order_only': {'order': 'largest-first', 'compression': 'none', 'shared_jit': False},
    'order_lzf': {'order': 'largest-first', 'compression': 'lzf', 'shared_jit': False},
}


def process_memory(pid):
    rows = {}
    for path in Path('/proc').glob('[0-9]*'):
        try:
            status = dict(line.split(':', 1) for line in (path/'status').read_text().splitlines() if ':' in line)
            rows[int(path.name)] = (int(status['PPid']), int(status.get('VmRSS', '0 kB').split()[0]))
        except (OSError, ValueError, KeyError):
            continue  # Processes can exit during a sample.
    selected = {pid}
    while True:
        children = {p for p, (parent, _) in rows.items() if parent in selected}
        if children <= selected:
            break
        selected |= children
    rss, pss, sampled = 0, 0, 0
    for child in selected:
        if child not in rows:
            continue
        try:
            rollup = dict(line.split(':', 1) for line in (Path('/proc')/str(child)/'smaps_rollup').read_text().splitlines() if ':' in line)
            child_pss = int(rollup['Pss'].split()[0])
        except (FileNotFoundError, ProcessLookupError):
            # Linux can return ESRCH, as well as ENOENT, after a child exits
            # between the ancestry snapshot and its smaps_rollup read.
            continue
        rss += rows[child][1]
        pss += child_pss
        sampled += 1
    if not sampled:
        raise RuntimeError('No process-tree memory sample available')
    return dict(time=time.monotonic(), RSS_KiB=rss, PSS_KiB=pss, processes=sampled)


class MemorySamples:
    def __init__(self):
        self.stop = threading.Event()
        self.rows, self.errors = [], []
        self.thread = threading.Thread(target=self.sample, daemon=True)

    def sample(self):
        while not self.stop.is_set():
            try:
                self.rows.append(process_memory(os.getpid()))
            except (OSError, RuntimeError, KeyError, ValueError) as error:
                self.errors.append(type(error).__name__+': '+str(error))
            self.stop.wait(1)

    def __enter__(self):
        self.thread.start()
        return self

    def __exit__(self, *args):
        self.stop.set()
        self.thread.join(timeout=15)


def probability_proof(row, genome, root, tmp_root):
    """Untimed, allocated native GPU predictions on one qualified small fixture."""
    root.mkdir()
    runtime = native.read_selected(row['receipt'])
    source_sha, receipt_sha = file_sha(genome), file_sha(row['receipt'])
    settings = native.inference_settings(runtime, runtime['model_profile'])
    signatures, exits = {}, []
    with tempfile.TemporaryDirectory(prefix='gg-prob-proof-', dir=tmp_root) as temporary:
        work = Path(temporary)
        env = runtime_environment(runtime)
        env.update(OMP_NUM_THREADS='4', MKL_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4', TMPDIR=str(work),
                   NUMBA_CACHE_DIR=str(work/'numba-cache'), TORCH_HOME=str(work/'torch-cache'), XDG_CACHE_HOME=str(work/'xdg-cache'))
        gpu = native.probe_gpu(runtime, env)
        for mode in ('none', 'lzf'):
            predicted = work/(mode+'.h5')
            argv = pipeline.command(runtime, 'prediction.py', ['-g', str(genome), '-m', runtime['model_path'],
                '-l', runtime['model_profile'], '-s', '25', '-p', str(predicted), '--batch_size', '32', '--num_workers', '2'])
            if settings['overlap_pred']:
                argv.append('--overlap_pred')
            if mode == 'lzf':
                argv.append('--comp')
            logs = root/mode
            logs.mkdir()
            exits.append(pipeline.execute(argv, env, work, logs, 'prediction'))
            digest_tool = Path(pipeline.__file__).with_name('annevo_probability_digest.py')
            output = root/(mode+'.signature.json')
            command = [runtime['python'], '-E', '-s', '-B', str(digest_tool), '--input', str(predicted), '--output', str(output)]
            exits.append(pipeline.execute(command, env, work, logs, 'probability_signature'))
            signatures[mode] = json.loads(output.read_text())
        native.read_selected(row['receipt'])
        if file_sha(genome) != source_sha or file_sha(row['receipt']) != receipt_sha:
            raise RuntimeError('Probability qualification inputs changed')
    if signatures['none']['semantic_sha256'] != signatures['lzf']['semantic_sha256']:
        raise RuntimeError('Native compression changed probability values, dtype, shape or metadata')
    result = dict(kind='native-annevo-probability-compression-proof', genome_sha256=source_sha,
                  receipt_sha256=receipt_sha, model_sha256=runtime['model_sha256'], runtime_digest=runtime['runtime_digest'],
                  inference_settings=settings, GPU=gpu, signatures=signatures, exits=exits,
                  all_probability_bits_and_metadata_identical=True, scope='qualified small fixture; outside performance timing')
    atomic_json(root/'comparison.json', result)
    return result


def sampled_pipeline(inputs, directory, options, config, shared_cache):
    directory.mkdir()
    started = time.monotonic()
    gpu, memory = GPUSamples(), MemorySamples()
    try:
        with gpu, memory:
            result = pipeline.run(argparse.Namespace(inputs=inputs, output=directory, **options,
                order=config['order'], compression=config['compression'],
                jit_cache_root=shared_cache if config['shared_jit'] else None))
            wall = time.monotonic()-started
    finally:
        atomic_json(directory/'resource_samples.json', dict(
            GPU_samples=gpu.rows, GPU_sample_errors=gpu.errors,
            memory_samples=memory.rows, memory_sample_errors=memory.errors))
    if gpu.errors or memory.errors or not gpu.rows or not memory.rows:
        raise RuntimeError('Performance/resource sampling was incomplete; inspect '+str(directory/'resource_samples.json'))
    return result, wall, gpu, memory


def shared_cache_preflight(row, genome, root, cache_root, options):
    """Actual cold/warm GPU+decode probes, outside timing and its JIT cache."""
    root.mkdir()
    inputs = root/'inputs.tsv'
    inputs.write_text('species\tgenome\treceipt\nfixture_probe\t'+str(genome)+'\t'+row['receipt']+'\n')
    config = CONFIGURATIONS['optimized']
    expected, results = None, []
    for repeat in range(2):
        directory = root/str(repeat)
        print('PREFLIGHT_STARTED', repeat, flush=True)
        result, wall, _, _ = sampled_pipeline(inputs, directory, options, config, cache_root)
        gff = result['species']['fixture_probe']['gff']
        signature = feature_digest(gff)
        if expected is not None and signature != expected:
            raise RuntimeError('Cold/warm native cache preflight changed GFF features')
        expected = signature
        logs = ''.join(p.read_text() for p in directory.glob('*/logs/decoding.stdout.txt'))
        compiled = len(list(cache_root.rglob('*.nbc')))
        loaded = logs.count('[cache] data loaded')
        if not compiled or (repeat and not loaded):
            raise RuntimeError('Native JIT cache creation/reuse was not demonstrated in preflight')
        results.append(dict(repeat=repeat, wall_seconds=wall, features_sha256=signature,
                            compiled_files=compiled, cache_load_events=loaded))
        print('PREFLIGHT_COMPLETE', repeat, flush=True)
    proof = dict(kind='native-annevo-shared-cache-preflight', results=results,
                 cold_warm_GFF_features_identical=True, outside_performance_timing=True,
                 separate_from_timed_cache=True, cache_directory=str(cache_root))
    atomic_json(root/'comparison.json', proof)
    return proof


def benchmark(args):
    if args.repeats < 2:
        raise ValueError('Use at least two alternating paired trials')
    root = args.output.absolute()
    root.mkdir(parents=True, exist_ok=False)
    rows, input_sha = pipeline.read_inputs(args.inputs, root)
    hashes = {str(args.inputs): input_sha}
    for row in rows:
        for key in ('genome', 'receipt', 'selection_receipt'):
            if row.get(key):
                hashes[row[key]] = file_sha(row[key])
    hashes[str(args.probability_proof_genome)] = file_sha(args.probability_proof_genome)
    proof = probability_proof(rows[0], args.probability_proof_genome, root/'probability_proof', args.tmp_root)
    options = dict(predict_threads=4, decode_threads=4, decoder_workers=2, batch_size=32,
                   sequence_threshold=25, tmp_root=args.tmp_root, max_tmp_gib=args.max_tmp_gib)
    # Shared compiled kernels are executable cache data. Use the allocation's
    # explicitly selected temporary filesystem, which must preserve private
    # POSIX permissions, rather than assuming the GFF output share does so.
    old_debug = os.environ.get('NUMBA_DEBUG_CACHE')
    os.environ['NUMBA_DEBUG_CACHE'] = '1'
    try:
        with tempfile.TemporaryDirectory(prefix='gg-benchmark-jit-', dir=args.tmp_root) as work:
            work = Path(work)
            preflight = shared_cache_preflight(rows[0], args.probability_proof_genome,
                                               root/'cache_preflight', work/'preflight', options)
            return benchmark_trials(args, root, hashes, proof, options, work/'timed', preflight)
    finally:
        if old_debug is None:
            os.environ.pop('NUMBA_DEBUG_CACHE', None)
        else:
            os.environ['NUMBA_DEBUG_CACHE'] = old_debug


def benchmark_trials(args, root, hashes, proof, options, shared_cache, preflight):
    trials, expected = [], {}
    for repeat in range(args.repeats):
        names = list(CONFIGURATIONS) if repeat % 2 == 0 else list(reversed(CONFIGURATIONS))
        for name in names:
            config = CONFIGURATIONS[name]
            directory = root/(str(repeat)+'_'+name)
            cache_before = len(list(shared_cache.rglob('*.nbc'))) if shared_cache.exists() else 0
            print('TRIAL_STARTED', repeat, name, flush=True)
            result, wall, gpu, memory = sampled_pipeline(args.inputs, directory, options, config, shared_cache)
            signatures = {}
            for species, row in result['species'].items():
                signatures[species] = dict(features_sha256=feature_digest(row['gff']),
                    raw_sha256=file_sha(row['gff']), counts=native.gff_counts(row['gff']))
                if species in expected and signatures[species]['features_sha256'] != expected[species]['features_sha256']:
                    raise RuntimeError('Native GFF annotation fields changed: '+species)
                expected.setdefault(species, signatures[species])
            cache_logs = ''.join(p.read_text() for p in directory.glob('*/logs/decoding.stdout.txt'))
            trial = dict(repeat=repeat, configuration=name, options=config, wall_seconds=wall,
                species=result['species'], signatures=signatures, pipeline=result,
                GPU_samples=gpu.rows, memory_samples=memory.rows,
                peak_GPU_MiB=max(x['memory_MiB'] for x in gpu.rows),
                peak_tree_RSS_KiB=max(x['RSS_KiB'] for x in memory.rows),
                peak_tree_PSS_KiB=max(x['PSS_KiB'] for x in memory.rows),
                jit_compiled_files_before=cache_before,
                jit_compiled_files_after=len(list(shared_cache.rglob('*.nbc'))) if shared_cache.exists() else 0,
                jit_cache_load_events=cache_logs.count('[cache] data loaded'),
                jit_cache_save_events=cache_logs.count('[cache] data saved'))
            trials.append(trial)
            atomic_json(root/'partial_results.json', trials)
            print('TRIAL_COMPLETE', repeat, name, wall, flush=True)
    if any(file_sha(path) != value for path, value in hashes.items()):
        raise RuntimeError('Benchmark inputs changed')
    medians = {name: statistics.median(t['wall_seconds'] for t in trials if t['configuration']==name) for name in CONFIGURATIONS}
    report = dict(schema=1, trials=trials, inputs_sha256=hashes, median_seconds=medians,
        baseline_median_seconds=medians['baseline'], optimized_median_seconds=medians['optimized'],
        speedup=medians['baseline']/medians['optimized'], wall_time_reduction_fraction=1-medians['optimized']/medians['baseline'],
        all_GFF_features_identical=True, probability_proof=proof, cache_preflight=preflight,
        cache_storage='private temporary filesystem selected by tmp-root; shared across timed trials, separate preflight cache',
        all_GFF_bytes_identical=all(t['signatures']==expected for t in trials),
        scientific_completion_verified=False,
        scope='Native full-genome GPU/decode equivalence and sampled resources; admission requires separate actual assays')
    atomic_json(root/'comparison.json', report)
    print(json.dumps({k:v for k,v in report.items() if k not in {'trials','inputs_sha256','probability_proof'}}), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=2)
    parser.add_argument('--tmp-root', type=Path, default=Path('/tmp'))
    parser.add_argument('--max-tmp-gib', type=float, default=128)
    parser.add_argument('--probability-proof-genome', type=Path, required=True)
    benchmark(parser.parse_args())


if __name__ == '__main__':
    main()

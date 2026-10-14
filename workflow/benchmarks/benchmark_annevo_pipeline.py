#!/usr/bin/env python3
"""Compare serial annotation and one-GPU/two-decoder execution on identical inputs.

Run in one site's allocated GPU-enabled GG container. Assets must already be
qualified; this harness does not submit jobs or silently acquire a model.
"""
import argparse
import json
import statistics
import subprocess
import sys
import threading
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
import annevo_pipeline as pipeline
import annevo_runtime as native
from annevo_prediction_equivalence import feature_digest
from annevo_runtime_download import atomic_json, file_sha


class GPUSamples:
    def __init__(self):
        self.stop = threading.Event()
        self.rows, self.errors = [], []
        self.thread = threading.Thread(target=self.sample, daemon=True)

    def sample(self):
        while not self.stop.is_set():
            try:
                result = subprocess.run(['nvidia-smi', '--query-gpu=uuid,utilization.gpu,memory.used',
                    '--format=csv,noheader,nounits'], text=True, capture_output=True, timeout=10, check=True)
                for line in result.stdout.strip().splitlines():
                    uuid, utilization, memory = [x.strip() for x in line.split(',')]
                    self.rows.append(dict(time=time.monotonic(), uuid=uuid,
                                          utilization_percent=float(utilization), memory_MiB=float(memory)))
            except (OSError, ValueError, subprocess.SubprocessError) as error:
                self.errors.append(str(error))
            self.stop.wait(1)

    def __enter__(self):
        self.thread.start()
        return self

    def __exit__(self, *args):
        self.stop.set()
        self.thread.join(timeout=15)


def benchmark(args):
    if args.repeats < 2:
        raise ValueError('Use at least two paired trials to assess noise')
    root = args.output.absolute()
    root.mkdir(parents=True, exist_ok=False)
    rows, input_sha = pipeline.read_inputs(args.inputs, root)
    hashes = {str(args.inputs): input_sha}
    for row in rows:
        for key in ('genome', 'receipt', 'selection_receipt'):
            if row.get(key):
                hashes[row[key]] = file_sha(row[key])
    options = dict(predict_threads=4, decode_threads=4, decoder_workers=2,
                   batch_size=32, sequence_threshold=25, tmp_root=args.tmp_root, max_tmp_gib=args.max_tmp_gib)
    if args.warmup_genome:
        warm = rows[0]
        runtime = native.read_selected(warm['receipt'])
        native.run(argparse.Namespace(receipt=warm['receipt'], genome=str(args.warmup_genome),
            output=str(root/'warmup/prediction.gff3'), output_receipt=str(root/'warmup/native_prediction.json'),
            model_profile=runtime['model_profile'], threads=4, batch_size=32, sequence_threshold=25))
        hashes[str(args.warmup_genome)] = file_sha(args.warmup_genome)
    trials, expected = [], {}
    for repeat in range(args.repeats):
        # Alternating AB/BA order gives both routes the earlier trial position.
        for mode in (('serial', 'pipeline') if repeat % 2 == 0 else ('pipeline', 'serial')):
            print('TRIAL_STARTED', repeat, mode, flush=True)
            directory = root / (str(repeat) + '_' + mode)
            started = time.monotonic()
            with GPUSamples() as samples:
                if mode == 'pipeline':
                    result = pipeline.run(argparse.Namespace(inputs=args.inputs, output=directory, **options))
                    species = result['species']
                else:
                    species = {}
                    for row in rows:
                        name = row['species']
                        runtime = native.read_selected(row['receipt'])
                        output = directory / name / (name+'.gff3')
                        receipt = output.parent/'native_prediction.json'
                        native.run(argparse.Namespace(receipt=row['receipt'], genome=row['genome'],
                            output=str(output), output_receipt=str(receipt), model_profile=runtime['model_profile'],
                            selection_receipt=row.get('selection_receipt', ''), species=name,
                            threads=4, batch_size=32, sequence_threshold=25))
                        species[name] = dict(gff=str(output), receipt=str(receipt))
                wall = time.monotonic()-started
            signatures = {}
            for name, row in species.items():
                signatures[name] = dict(features_sha256=feature_digest(row['gff']),
                                       raw_sha256=file_sha(row['gff']), counts=native.gff_counts(row['gff']))
                if name in expected and signatures[name]['features_sha256'] != expected[name]['features_sha256']:
                    raise RuntimeError('GFF features, IDs, phase or scores changed: ' + name)
                expected.setdefault(name, signatures[name])
            trial = dict(repeat=repeat, mode=mode, wall_seconds=wall, species=species,
                         signatures=signatures, GPU_samples=samples.rows, GPU_sample_errors=samples.errors)
            trials.append(trial)
            atomic_json(root/'partial_results.json', trials)
            print('TRIAL_COMPLETE', repeat, mode, wall, flush=True)
    if any(file_sha(path) != sha for path, sha in hashes.items()):
        raise RuntimeError('Benchmark input changed')
    serial = statistics.median(t['wall_seconds'] for t in trials if t['mode']=='serial')
    concurrent = statistics.median(t['wall_seconds'] for t in trials if t['mode']=='pipeline')
    report = dict(schema=1, trials=trials, inputs_sha256=hashes,
                  serial_median_seconds=serial, pipeline_median_seconds=concurrent,
                  wall_time_reduction_fraction=1-concurrent/serial, speedup=serial/concurrent,
                  all_GFF_features_identical=True,
                  all_GFF_bytes_identical=all(t['signatures']==expected for t in trials),
                  scientific_completion_verified=False,
                  scope='Native GPU inference and CPU decode equivalence/performance; rescue acceptance needs its separate actual comparison')
    atomic_json(root/'comparison.json', report)
    print(json.dumps({k:v for k,v in report.items() if k not in ('trials', 'inputs_sha256')}), flush=True)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--repeats', type=int, default=2)
    p.add_argument('--tmp-root', type=Path, default=Path('/tmp'))
    p.add_argument('--max-tmp-gib', type=float, default=128)
    p.add_argument('--warmup-genome', type=Path, help='Qualified small matched-profile genome, excluded from trial timing')
    benchmark(p.parse_args())


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Measure candidate input-state annotation with complete cache-key fingerprints."""
import argparse
import hashlib
import json
import platform
import resource
import statistics
import subprocess
import sys
import time
from pathlib import Path


def fingerprint(frame):
    payload = {'table': frame.to_json(orient='split'),
               'dtypes': [str(dtype) for dtype in frame.dtypes],
               'index_name': frame.index.name}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def worker(args):
    sys.path.insert(0, str(args.support_root))
    import csubst_scan_candidate_sites as scan
    import pandas as pd

    family_ids = [f'OG{i:07}' for i in range(args.families)]
    frame = pd.DataFrame({
        'orthogroup': [family_ids[i % args.families] for i in range(args.candidates)],
        '_analysis_key': [hashlib.sha256(f'analysis{i}'.encode()).hexdigest() for i in range(args.candidates)],
        '_candidate_id': [f'candidate{i}_site{i % 1024 + 1}' for i in range(args.candidates)],
        'trait': ['aquatic' if i % 2 else 'terrestrial' for i in range(args.candidates)],
        'codon_site_alignment': [i % 1024 + 1 for i in range(args.candidates)],
    })
    for i in range(32):
        frame[f'probability_{i}'] = [float('nan') if j % 17 == 0 else ((i + j) % 100) / 100
                                    for j in range(args.candidates)]
    frame.index = pd.Index([i * 3 + 7 for i in range(args.candidates)], name='source_row')
    states = {family: {'missing_required_inputs': [f'stat/{family}.tsv', f'tree/{family}.nwk'] if i % 7 == 0 else [],
                       'required_input_signature': hashlib.sha256(f'input{i}'.encode()).hexdigest()}
              for i, family in enumerate(family_ids)}
    original = fingerprint(frame)
    started = time.perf_counter()
    annotated = scan.annotate_candidate_input_state(frame, states)
    elapsed = time.perf_counter() - started
    assert fingerprint(frame) == original
    sample = {'seconds': elapsed, 'output_sha256': fingerprint(annotated),
              'missing_candidates': int(annotated['_missing_required_inputs'].ne('').sum()),
              'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    print(json.dumps(sample))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--candidates', type=int, default=8192)
    parser.add_argument('--families', type=int, default=256)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 32 <= args.candidates <= 32768 or not 16 <= args.families <= min(2048, args.candidates):
        parser.error('Require 32..32768 candidates and 16..min(2048, candidates) families')
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error('--output is required')
    results = {'support_root': str(args.support_root), 'candidates': args.candidates, 'families': args.families,
               'python': sys.version, 'platform': platform.platform(), 'samples': []}
    for trial in range(4):
        raw = subprocess.check_output([sys.executable, str(Path(__file__).resolve()), '--worker',
            '--support-root', str(args.support_root), '--candidates', str(args.candidates),
            '--families', str(args.families)], text=True)
        sample = json.loads(raw)
        if trial:
            results['samples'].append(sample)
        print('warmup' if trial == 0 else trial, sample['seconds'], flush=True)
    assert len({sample['output_sha256'] for sample in results['samples']}) == 1
    for metric in ('seconds', 'peak_rss_kib_linux'):
        results['median_' + metric] = statistics.median(sample[metric] for sample in results['samples'])
    args.output.write_text(json.dumps(results, indent=2) + '\n')
    print('medians', {key: value for key, value in results.items() if key.startswith('median_')})


if __name__ == '__main__':
    main()

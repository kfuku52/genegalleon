"""Compare full rewriting with verified sparse admission publication updates.

Synthetic 35,000-candidate/5,000-body input, 350 changed decisions, two warmups
and five fresh-process repetitions per mode. This measures publication after
decisions, not complete rescue speed or cold filesystem throughput. Both paths
verify original bytes; output parity covers all ordered model-store projections.
"""
import argparse
import json
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
import rescue_model_store as store  # noqa: E402
from benchmark_rescue_model_selection import digest_records, records  # noqa: E402
from rescue_model_store_patch import patch_model_store  # noqa: E402


def receipt(directory):
    files = {str(p.relative_to(directory)): store._digest(p) for p in directory.rglob('*') if p.is_file()}
    (directory / 'receipt.json').write_text(json.dumps({'key': {'species': 'synthetic'}, 'files': files}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--mode', choices=('baseline', 'patch'), help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        source = args.worker / 'source'
        changes = {int(i): m for i, m in json.loads((args.worker / 'changes.json').read_text()).items()}
        with tempfile.TemporaryDirectory(dir=args.worker, prefix='publication-') as temporary:
            destination = Path(temporary)
            start = time.perf_counter()
            if args.mode == 'baseline':
                rows = (changes.get(i, m) for i, m in enumerate(store.iter_models(source)))
                store.write_model_store(destination, rows, codec='gzip', shard_bytes=1024 * 1024)
            else:
                patch_model_store(source, destination, changes)
            seconds = time.perf_counter() - start
            receipt(destination)
            proof = {kind: digest_records(store._iter(destination, kind)) for kind in sorted(store.KINDS)}
            print(json.dumps({'seconds': seconds, 'proof': proof,
                              'maxrss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))
        return
    with tempfile.TemporaryDirectory(prefix='gg-admission-publication-') as temporary:
        root = Path(temporary)
        source = root / 'source'
        manifest = store.write_model_store(source, records(35000), codec='gzip', shard_bytes=1024 * 1024)
        receipt(source)
        changes = {i: {**m, 'status': 'accepted', 'decision': {'basis': 'synthetic-current-admission'}}
                   for i, m in enumerate(records(35000)) if i % 100 == 0}
        (root / 'changes.json').write_text(json.dumps(changes))
        modes, expected = {}, None
        for mode in ('baseline', 'patch'):
            runs = []
            for i in range(7):
                result = json.loads(subprocess.check_output([sys.executable, str(Path(__file__).resolve()),
                                    '--worker', str(root), '--mode', mode], text=True))
                if expected is not None and result['proof'] != expected:
                    raise ValueError('Admission publication output parity failed')
                expected = result['proof']
                if i >= 2:
                    runs.append(result)
            modes[mode] = {'median_seconds': statistics.median(r['seconds'] for r in runs),
                           'median_maxrss_kib': statistics.median(r['maxrss_kib'] for r in runs), 'runs': runs}
        print(json.dumps({'records': 35000, 'bodies': manifest['counts']['bodies'], 'changes': len(changes),
                          'warmups': 2, 'repeats': 5, 'modes': modes, 'python': sys.version,
                          'platform': platform.platform(), 'scope': 'publication update only; ordered projections verified'}))


if __name__ == '__main__':
    main()

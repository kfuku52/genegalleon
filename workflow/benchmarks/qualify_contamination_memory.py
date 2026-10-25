#!/usr/bin/env python3
"""Qualify FASTA filtering under an externally applied Docker/cgroup RAM limit.

Use a freshness-checked GeneGalleon image with --memory=256m --memory-swap=256m.
The controller streams the fixture before starting one measured worker. Its
memory, the worker and filesystem cache all share the same cgroup limit.
"""

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

from benchmark_contamination_memory import create_fixture


def memory_events():
    return dict((key, int(value)) for key, value in
                (line.split() for line in Path('/sys/fs/cgroup/memory.events').read_text().splitlines()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--expect', choices=('oom', 'fits'), required=True)
    parser.add_argument('--memory-mib', type=int, default=256)
    parser.add_argument('--records', type=int, default=10000)
    parser.add_argument('--length', type=int, default=40960)
    parser.add_argument('--reference-report', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or min(args.memory_mib, args.records, args.length) < 1:
        parser.error('Use a new output file and positive memory/fixture sizes')
    limit = Path('/sys/fs/cgroup/memory.max').read_text().strip()
    swap_limit = Path('/sys/fs/cgroup/memory.swap.max').read_text().strip()
    if limit != str(args.memory_mib * 2**20) or swap_limit != '0':
        parser.error('The requested cgroup memory limit and disabled swap must be applied externally')
    source = args.source.resolve()
    benchmark = Path(__file__).with_name('benchmark_contamination_memory.py').resolve()
    with tempfile.TemporaryDirectory(prefix='gg-ram-limit-') as temporary:
        root = Path(temporary)
        fixture = create_fixture(root, args.records, args.length, 'retain-most', 'plain')
        before = memory_events()
        result = subprocess.run([sys.executable, str(benchmark), '--source', str(source),
                                 '--worker-root', str(root), '--codec', 'plain'],
                                capture_output=True, text=True)
        after = memory_events()
        payload = {'source': str(source), 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                   'returncode': result.returncode, 'memory_max': limit, 'memory_swap_max': swap_limit,
                   'memory_peak': Path('/sys/fs/cgroup/memory.peak').read_text().strip(),
                   'memory_events_delta': {key: after[key] - before[key] for key in before},
                   'tempdir': tempfile.gettempdir(), 'fixture': fixture}
        if result.returncode == 0:
            payload['worker'] = json.loads(result.stdout)
        else:
            payload['stderr'] = result.stderr[-2000:]
        if args.reference_report:
            reference = json.loads(args.reference_report.read_text())
            payload['reference_outputs_match'] = (fixture == reference['fixture'] and result.returncode == 0
                and payload['worker']['outputs'] == reference['trials']['baseline'][0]['outputs'])
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2) + '\n')
        print(json.dumps(payload))
        if args.expect == 'oom':
            assert result.returncode == -9 and payload['memory_events_delta']['oom_kill'] >= 1
        else:
            assert result.returncode == 0 and payload['memory_events_delta']['oom_kill'] == 0
        if args.reference_report:
            assert payload['reference_outputs_match'], 'Limited-memory outputs differ from the reference'


if __name__ == '__main__':
    main()

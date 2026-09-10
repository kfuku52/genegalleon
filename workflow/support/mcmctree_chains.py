#!/usr/bin/env python3
"""Run independent MCMCtree chains, retain evidence, and diagnose without gating.

Completed chains can be reused after interruption. Incomplete chains are rerun
from their original seed in a new attempt; binary checkpoints are never spliced.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import signal
import threading
import time
import tempfile


STOP = threading.Event()
PROCESSES = set()


def interrupt(signum, frame):
    STOP.set()
    for process in list(PROCESSES):
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass


def save_json(path, value):
    temporary = path.with_suffix('.pending')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def control(text, **values):
    text = text.replace('FossilErrprint', 'FossilErr\nprint')
    for key, value in values.items():
        text = re.sub(r'^\s*' + re.escape(key) + r'\s*=.*\n?', '', text, flags=re.M)
        text = text.rstrip() + f'\n{key} = {value}\n'
    return text


def sample_stream(path):
    """Validate while streaming: wide MCMC files need not fit in Python memory."""
    with path.open() as stream:
        header = stream.readline().split()
        if len(header) < 2 or header[0] != 'Gen' or len(set(header)) != len(header):
            raise ValueError('Invalid sample header')
        yield header
        previous = -math.inf
        count = 0
        for line in stream:
            if not line.strip():
                continue
            row = line.split()
            if len(row) != len(header) or not all(math.isfinite(float(x)) for x in row):
                raise ValueError('Incomplete or non-finite sample row')
            if float(row[0]) <= previous:
                raise ValueError('Non-increasing sample index')
            previous = float(row[0])
            count += 1
            yield row
        if not count:
            raise ValueError('Empty sample file')


def samples(path):
    rows = sample_stream(path)
    return next(rows), list(rows)


def positive_control_integer(template, name):
    values = re.findall(r'^\s*' + re.escape(name) + r'\s*=\s*([^\s*#]+)', template, re.M)
    if len(values) != 1 or not re.fullmatch(r'[0-9]+', values[0]) or int(values[0]) < 1:
        raise ValueError(f'Expected one positive integer {name} in control file')
    return int(values[0])


def wait_process(argv, directory, out, err):
    if STOP.is_set():
        return -signal.SIGTERM
    process = subprocess.Popen(argv, cwd=directory, stdout=out, stderr=err,
                               start_new_session=True,
                               env={**os.environ, 'OMP_NUM_THREADS': '1'})
    PROCESSES.add(process)
    stopped_at = None
    try:
        while True:
            if STOP.is_set():
                if stopped_at is None:
                    stopped_at = time.monotonic()
                    try:
                        os.killpg(process.pid, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                elif time.monotonic() - stopped_at > 3:
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
            try:
                return process.wait(timeout=0.2)
            except subprocess.TimeoutExpired:
                continue
    finally:
        PROCESSES.discard(process)


def execute(binary, directory, ctl):
    with (directory / 'stdout.log').open('w') as out, (directory / 'stderr.log').open('w') as err:
        return wait_process([binary, ctl], directory, out, err)


def reusable_attempt(chain, receipt, index, seed):
    required = {'run.ctl', 'mcmc.txt', 'chain.out', 'stdout.log', 'stderr.log'}
    try:
        prior = json.loads(receipt.read_text())
        name = prior['attempt']
        if not isinstance(name, str) or Path(name).name != name or not name.startswith('attempt-'):
            return None
        attempt = chain / name
        if attempt.is_symlink() or prior['state'] != 'completed' or prior['seed'] != seed or prior['chain'] != index:
            return None
        if set(prior['hashes']) != required:
            return None
        if all((attempt / name).is_file() and not (attempt / name).is_symlink()
               and digest(attempt / name) == value for name, value in prior['hashes'].items()):
            return attempt
    except (OSError, ValueError, KeyError, TypeError):
        pass
    return None


def run_chain(root, template, files, binary, index, seed, expected, frequency):
    chain = root / f'chain-{index:02d}'
    chain.mkdir(exist_ok=True)
    receipt = chain / 'completed.json'
    previous = reusable_attempt(chain, receipt, index, seed)
    if previous is not None:
        return previous
    attempt = Path(tempfile.mkdtemp(prefix='attempt-', dir=chain))
    state = {'chain': index, 'seed': seed, 'state': 'running', 'attempt': attempt.name, 'started_at': time.time()}
    save_json(attempt / 'state.json', state)
    try:
        for path in files:
            shutil.copy2(path, attempt / path.name)
        (attempt / 'run.ctl').write_text(control(template, seed=seed, print=1,
                                                outfile='chain.out', mcmcfile='mcmc.txt', checkpoint=1))
        print(f'MCMCtree chain {index}: seed={seed}, evidence={attempt}', flush=True)
        code = execute(binary, attempt, 'run.ctl')
        state['exit_code'] = code
        if code != 0:
            raise ValueError(f'MCMCtree exited {code}')
        rows = sample_stream(attempt / 'mcmc.txt')
        next(rows)  # header
        count, last = 0, None
        for row in rows:
            count += 1
            last = float(row[0])
        if count < expected or last != expected * frequency or not (attempt / 'chain.out').stat().st_size:
            raise ValueError('Missing output or unexpected sample count')
        state.update(state='completed', finished_at=time.time(), hashes={name: digest(attempt / name)
                     for name in ('run.ctl', 'mcmc.txt', 'chain.out', 'stdout.log', 'stderr.log')})
        save_json(attempt / 'state.json', state)
        save_json(receipt, state)
        return attempt
    except (OSError, ValueError) as error:
        state.update(state='interrupted' if STOP.is_set() else 'failed', finished_at=time.time(), error=str(error))
        save_json(attempt / 'state.json', state)
        return None


def diagnose(paths, root):
    output = root / 'diagnostics.tsv'
    try:
        output.unlink(missing_ok=True)
        with (root / 'diagnostics.log').open('w') as log:
            code = wait_process(['Rscript', str(Path(__file__).with_name('mcmctree_diagnostics.R')),
                                 str(output), *map(str, paths)], root, log, log)
        if code:
            raise ValueError('Diagnostics unavailable; see diagnostics.log')
        with output.open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        if not rows:
            raise ValueError('No diagnostic parameters')
        failed = []
        for row in rows:
            try:
                rhat, bulk, tail = [float(row[k]) for k in ('rhat', 'ess_bulk', 'ess_tail')]
            except ValueError:
                raise ValueError('Undefined diagnostics (including constant parameters)')
            if not all(map(math.isfinite, (rhat, bulk, tail))):
                raise ValueError('Non-finite diagnostics')
            if rhat >= 1.01 or bulk < 400 or tail < 400:
                failed.append(row['parameter'])
        return {'state': 'unconverged' if failed else 'passed', 'failed_parameters': failed}
    except (OSError, ValueError, KeyError, TypeError) as error:
        return {'state': 'diagnostic_failed', 'reason': str(error)}


def run(args):
    if not 0 < args.seed < 2147483643 or not 1 <= args.jobs <= 4:
        raise ValueError('seed must be 1..2147483642 and jobs must be 1..4')
    template_dir = args.template.resolve()
    files = sorted(p for p in template_dir.iterdir() if p.is_file())
    template = (template_dir / args.control).read_text()
    expected = positive_control_integer(template, 'nsample')
    frequency = positive_control_integer(template, 'sampfreq')
    binary = Path(shutil.which('mcmctree') or '').resolve()
    if not binary.is_file():
        raise ValueError('mcmctree not found')
    identity = {'inputs': {p.name: digest(p) for p in files}, 'binary_sha256': digest(binary),
                'seed': args.seed, 'chains': 4, 'time_factor': args.time_factor,
                'runner_sha256': digest(Path(__file__)),
                'diagnostics_sha256': digest(Path(__file__).with_name('mcmctree_diagnostics.R'))}
    run_id = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    root = args.store.resolve() / run_id
    root.mkdir(parents=True, exist_ok=True)
    with (root / '.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            save_json(root / 'manifest.json', identity)
            running = {'state': 'running', 'run': str(root), 'policy': 'warn',
                       'started_at': time.time(), 'parallel_jobs': args.jobs}
            save_json(root / 'status.json', running)
            save_json(args.status, running)
            with concurrent.futures.ThreadPoolExecutor(max_workers=min(4, args.jobs)) as pool:
                tasks = [pool.submit(run_chain, root, template, files, str(binary), i + 1,
                                     args.seed + i, expected, frequency) for i in range(4)]
                attempts = [task.result() for task in tasks]
            if STOP.is_set():
                save_json(args.status, {'state': 'interrupted', 'run': str(root), 'policy': 'warn'})
                save_json(root / 'status.json', {'state': 'interrupted', 'policy': 'warn'})
                raise ValueError('MCMCtree run interrupted; evidence retained')
            good = [p for p in attempts if p is not None]
            diagnostic_dir = Path(tempfile.mkdtemp(prefix='diagnostics-', dir=root))
            status = diagnose([p / 'mcmc.txt' for p in good], diagnostic_dir) if len(good) >= 2 else {
                'state': 'diagnostic_failed', 'reason': 'Fewer than two complete chains'}
            if len(good) != 4:
                status.update(diagnostic_state=status['state'], state='incomplete',
                              reason=f'{len(good)} of 4 chains completed')
            status['diagnostics'] = str(diagnostic_dir)
            status.update(run=str(root), completed_chains=len(good), expected_chains=4,
                          policy='warn', rhat_max=1.01, ess_min=400, time_factor=args.time_factor)
            interim = {**status, 'state': 'summarizing', 'diagnostic_state': status['state']}
            save_json(args.status, interim)
            save_json(root / 'status.json', interim)
            if status['state'] != 'passed':
                print(f"WARNING: MCMCtree {status['state']}; downstream continues with provisional ages. See {args.status}", file=sys.stderr)
            if not good:
                raise ValueError('No complete MCMCtree chains; no usable tree to publish')
            # PAML's documented print=-1 mode summarizes retained draws without MCMC.
            summary = Path(tempfile.mkdtemp(prefix='summary-', dir=root))
            for path in files:
                shutil.copy2(path, summary / path.name)
            header = None
            count = 0
            with (summary / 'mcmc.txt').open('w') as output:
                for path in good:
                    rows = sample_stream(path / 'mcmc.txt')
                    names = next(rows)
                    if header is None:
                        header = names
                        output.write('\t'.join(names) + '\n')
                    elif names != header:
                        raise ValueError('Cannot pool chains with different parameter headers')
                    for row in rows:
                        count += 1
                        output.write('\t'.join([str(count), *row[1:]]) + '\n')
            (summary / 'run.ctl').write_text(control(template, print=-1, seed=args.seed,
                                                    outfile='summary.out', mcmcfile='mcmc.txt',
                                                    nsample=count, checkpoint=0))
            if execute(str(binary), summary, 'run.ctl'):
                raise ValueError(f'MCMCtree summary failed; see {summary}')
            if not (summary / 'summary.out').is_file() or not (summary / 'summary.out').stat().st_size:
                raise ValueError(f'MCMCtree produced no summary; see {summary}')
            pending = args.output.with_suffix(args.output.suffix + '.pending')
            shutil.copy2(summary / 'summary.out', pending)
            pending.replace(args.output)
            status.update(summary=str(summary), pooled_samples=count)
            save_json(root / 'status.json', status)
            save_json(args.status, status)
        except (OSError, ValueError, KeyError, TypeError) as error:
            failed = {'state': 'interrupted' if STOP.is_set() else 'failed',
                      'run': str(root), 'reason': str(error), 'policy': 'warn'}
            save_json(root / 'status.json', failed)
            save_json(args.status, failed)
            raise ValueError(str(error)) from error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--template', type=Path, required=True)
    parser.add_argument('--control', required=True)
    parser.add_argument('--store', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--status', type=Path, required=True)
    parser.add_argument('--seed', type=int, default=1729)
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--time-factor', default='1')
    args = parser.parse_args()
    signal.signal(signal.SIGTERM, interrupt)
    signal.signal(signal.SIGINT, interrupt)
    try:
        run(args)
    except (OSError, ValueError) as error:
        print(f'Error: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())

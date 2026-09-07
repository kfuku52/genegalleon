#!/usr/bin/env python3
"""Supervise a workflow's disposable scratch space; durable state stays in workspace."""
import argparse
import fcntl
import hashlib
import json
import math
import re
import os
from pathlib import Path
import shutil
import signal
import socket
import subprocess
import sys
import tempfile
import time


def owned_directory(path):
    if path.is_symlink():
        raise ValueError(f"Refusing symlinked scratch directory: {path}")
    path.mkdir(mode=0o700, exist_ok=True)
    if not path.is_dir() or path.stat().st_uid != os.getuid():
        raise ValueError(f"Scratch directory is not owned by this user: {path}")
    if path.stat().st_mode & 0o077:
        raise ValueError(f"Scratch directory must be private (mode 700): {path}")
    return path


def write_record(path, record):
    staged = path.with_suffix('.new')
    staged.write_text(json.dumps(record, sort_keys=True) + '\n')
    os.replace(staged, path)


def lock_run(path):
    fd = os.open(path / 'lock', os.O_RDWR | os.O_CREAT | os.O_NOFOLLOW, 0o600)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        os.close(fd)
        return None
    return fd


def retained_runs(scope):
    """Return only recognized, idle runs. Caller holds the scope lock."""
    result = []
    for path in scope.glob('run-*'):
        if path.is_symlink() or not path.is_dir() or path.stat().st_uid != os.getuid():
            continue
        record_path = path / 'record.json'
        if record_path.is_symlink() or not record_path.is_file():
            continue
        try:
            record = json.loads(record_path.read_text())
            if (not isinstance(record, dict)
                    or record.get('schema') != 1 or record.get('run') != path.name
                    or not isinstance(record.get('task'), str)
                    or not isinstance(record.get('updated'), (int, float))
                    or not math.isfinite(record['updated']) or record['updated'] < 0):
                continue
            fd = lock_run(path)
            if fd is not None:
                result.append((path, record, fd))
        except (OSError, ValueError):
            continue
    return result


def usage(path):
    size = count = 0
    for parent, dirs, files in os.walk(path, followlinks=False):
        for name in files:
            size += (Path(parent) / name).lstat().st_size
            count += 1
    return size, count


def prune(scope, task, reuse, limits):
    if not task and not any(limits):
        return None
    runs = retained_runs(scope)
    selected = None
    try:
        # A preserved task may be resumed only while no process holds its lock.
        matching = sorted((r for r in runs if r[1].get('identity', r[1]['task']) == task),
                          key=lambda r: r[1]['updated'], reverse=True)
        if reuse and matching:
            selected = matching[0]
        elif not reuse:
            for path, record, fd in matching:
                shutil.rmtree(path)
        survivors = sorted((r for r in runs if r != selected and r[0].exists()),
                           key=lambda r: r[1]['updated'])
        days, max_dirs, max_bytes, max_files = limits
        sizes = {p: usage(p) for p, _, _ in survivors} if max_bytes or max_files else {}
        while survivors:
            p, record, _ = survivors[0]
            too_old = days > 0 and time.time() - record['updated'] > days * 86400
            too_many = max_dirs > 0 and len(survivors) > max_dirs
            too_big = max_bytes > 0 and sum(sizes[r[0]][0] for r in survivors) > max_bytes
            too_full = max_files > 0 and sum(sizes[r[0]][1] for r in survivors) > max_files
            if not any((too_old, too_many, too_big, too_full)):
                break
            shutil.rmtree(p)
            survivors.pop(0)
        return selected
    finally:
        for run in runs:
            if run != selected:
                os.close(run[2])


def gene_task_identity(task, mode, workspace):
    """Do not treat reassigned array slots as the previous gene family.

    For query mode hash the filename inventory and sorting locale, avoiding a
    second implementation of the core's locale-sensitive `find | sort` choice.
    For orthogroups the table row gives the family directly; stop at that row.
    Missing input is left to core validation, and cannot match resolved input.
    """
    workspace = Path(workspace)
    try:
        if mode == 'orthogroup':
            table = workspace / 'output/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv'
            with table.open(encoding='utf-8', errors='surrogateescape', newline='') as handle:
                next(handle, None)
                family = next((line.split('\t', 1)[0].rstrip('\n')
                               for index, line in enumerate(handle, 1) if index == int(task)), None)
            if not family:
                return task + ':unresolved'
            identity = family
        else:
            names = sorted(p.name for p in (workspace / 'input/query_gene').iterdir()
                           if not p.name.startswith('.') and not p.is_symlink() and p.is_file())
            identity = json.dumps([names, [os.environ.get(k, '')
                                  for k in ('LC_ALL', 'LC_COLLATE', 'LANG')]])
    except FileNotFoundError:
        return task + ':unresolved'
    return task + ':' + hashlib.sha256(identity.encode('utf-8', errors='surrogateescape')).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workflow', required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command:
        parser.error('a command is required')
    task = os.environ.get('GG_ARRAY_TASK_ID', '1')
    if not re.fullmatch(r'[0-9]+', task) or int(task) < 1:
        raise ValueError('GG_ARRAY_TASK_ID must be a positive integer')
    is_gene = args.workflow == 'gg_gene_evolution'
    mode = (os.environ.get('mode_gene_evolution') or 'query2family').lower() if is_gene else ''
    if is_gene and mode not in {'orthogroup', 'query2family'}:
        raise ValueError(f'Invalid mode_gene_evolution: {mode}')
    delete = os.environ.get('delete_tmp_dir', '1')
    preexisting = os.environ.get('delete_preexisting_tmp_dir', '1')
    if delete not in {'0', '1'} or (is_gene and preexisting not in {'0', '1'}):
        raise ValueError('Scratch deletion flags must be 0 or 1')
    reuse = is_gene and preexisting == '0'
    identity = gene_task_identity(task, mode, os.environ.get('gg_workspace_dir', '/workspace')) if is_gene else task
    defaults = (7, 100, 107374182400, 100000) if is_gene else (0, 0, 0, 0)
    names = ('retention_days', 'max_dirs', 'max_bytes', 'max_files')
    limits = tuple(int(os.environ.get(f'gene_family_tmp_{name}', str(default)))
                   if is_gene else 0 for name, default in zip(names, defaults))
    if any(value < 0 for value in limits):
        raise ValueError('Temporary retention limits must be nonnegative')
    root = Path(os.environ['GG_TMP_MOUNT'])
    workspace = os.environ.get('GG_TMP_WORKSPACE_ID', os.environ.get('gg_workspace_dir', '/workspace'))
    key = hashlib.sha256(workspace.encode()).hexdigest()[:24]
    user_root = owned_directory(root / f'genegalleon-{os.getuid()}')
    workspace_root = owned_directory(user_root / key)
    workflow_key = hashlib.sha256((args.workflow + '\0' + mode).encode()).hexdigest()[:24]
    scope = owned_directory(workspace_root / workflow_key)
    scope_fd = os.open(scope / 'scope.lock', os.O_RDWR | os.O_CREAT | os.O_NOFOLLOW, 0o600)
    fcntl.flock(scope_fd, fcntl.LOCK_EX)
    # Only gene evolution has an explicit preexisting-task deletion policy.
    selected = prune(scope, identity if is_gene else '', reuse, limits)
    if selected:
        run, record, run_fd = selected
    else:
        run = Path(tempfile.mkdtemp(prefix='run-', dir=scope))
        run_fd = lock_run(run)
        record = {'schema': 1, 'run': run.name, 'task': task, 'workflow': args.workflow,
                  'workspace': workspace, 'mode': mode, 'identity': identity, 'host': socket.gethostname()}
    record.update(updated=time.time(), state='running')
    write_record(run / 'record.json', record)
    fcntl.flock(scope_fd, fcntl.LOCK_UN)
    work = owned_directory(run / 'work')
    runtime = owned_directory(run / 'runtime')
    env = os.environ.copy()
    env.update(GG_TMP_TASK_ROOT=str(work), GG_TMP_LOCK_FD=str(run_fd),
               TMPDIR=str(runtime), TMP=str(runtime), TEMP=str(runtime))
    host_path = str(Path(os.environ.get('GG_TMP_HOST_ROOT', str(root))) / run.relative_to(root))
    print(f'GeneGalleon scratch: {host_path} (container: {run}); free bytes: {shutil.disk_usage(run).free}', flush=True)
    child = None
    interrupted = 0

    def forward(signum, _frame):
        nonlocal interrupted
        interrupted = signum
        if child is None:
            return
        try:
            os.killpg(child.pid, signum)
        except ProcessLookupError:
            pass
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(sig, forward)
    child = subprocess.Popen(command, env=env, stdin=sys.stdin, start_new_session=True, pass_fds=(run_fd,))
    if interrupted:
        forward(interrupted, None)
    status = child.wait()
    if interrupted:
        status = 128 + interrupted
    record.update(updated=time.time(), state='finished', exit_code=status)
    write_record(run / 'record.json', record)
    # Core exit 8 means an already-completed gene family.
    success = status == 0 or (is_gene and status == 8)
    if success and delete == '1':
        shutil.rmtree(run)
    else:
        print(f'GeneGalleon scratch retained: {host_path}', flush=True)
    os.close(run_fd)
    fcntl.flock(scope_fd, fcntl.LOCK_EX)
    # Apply retention without treating completion as a new same-task invocation.
    prune(scope, '', False, limits)
    os.close(scope_fd)
    return status if status >= 0 else 128 - status


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (OSError, ValueError) as exc:
        print(f'GeneGalleon scratch error: {exc}', file=sys.stderr)
        sys.exit(1)

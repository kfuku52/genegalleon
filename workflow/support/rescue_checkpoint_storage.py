"""Share byte-identical, immutable Miniprot checkpoint payloads."""
import argparse
import errno
import json
import os
import stat
import tempfile
from pathlib import Path

try:
    from fasta_sequence_store import exclusive_lock
    from input_generation_array_state import digest
except ImportError:
    from .fasta_sequence_store import exclusive_lock
    from .input_generation_array_state import digest


# Postprocessing only reads these predictor outputs. Files it updates are copied.
SHARED = frozenset({
    'genome.unique.gff', 'candidates.json', 'genome_prediction_query_mapping.tsv',
    'intervals.zip', 'local_search_inputs.jsonl.gz', 'genome_search_nomination.json',
    'genome_search_inputs.jsonl.gz', 'model_repair_search/intervals.zip',
    'model_repair_search/local_search_inputs.jsonl.gz',
})
UNSUPPORTED = {errno.EXDEV, errno.EPERM, errno.EACCES, errno.EMLINK, errno.ENOSYS, errno.EOPNOTSUPP}


def link_member(source, target):
    if not stat.S_ISREG(source.lstat().st_mode):
        raise ValueError('Miniprot checkpoint member must be a regular file')
    try:
        os.link(source, target, follow_symlinks=False)
    except OSError as exc:
        if exc.errno not in UNSUPPORTED:
            raise
        return False
    return True


def file_identity(info):
    return info.st_dev, info.st_ino, info.st_mode, info.st_nlink, info.st_size, info.st_mtime_ns, info.st_ctime_ns


def retire_copy(source, target, checksum):
    """Atomically replace one duplicate directory entry with a verified hardlink."""
    old = target.lstat()
    if not stat.S_ISREG(old.st_mode) or not stat.S_ISREG(source.lstat().st_mode):
        raise ValueError('Expected regular immutable Miniprot files')
    original = source.lstat()
    if (stat.S_IMODE(old.st_mode), old.st_uid, old.st_gid) != (
            stat.S_IMODE(original.st_mode), original.st_uid, original.st_gid):
        return 0
    if os.path.samefile(source, target):
        return 0
    if digest(target) != checksum or digest(source) != checksum:
        raise ValueError('Miniprot output differs from its receipt: ' + str(target))
    with tempfile.TemporaryDirectory(prefix='.gg-checkpoint-', dir=target.parent) as tmp:
        linked = Path(tmp) / 'verified'
        if not link_member(source, linked):
            return 0
        if file_identity(target.lstat()) != file_identity(old) or not os.path.samefile(source, linked) or digest(linked) != checksum:
            raise ValueError('Miniprot output changed during retirement: ' + str(target))
        os.replace(linked, target)
        return old.st_size if old.st_nlink == 1 else 0


def deduplicate(root, *, apply=False, task_index=None):
    """Retire duplicate payloads for a closed, same-plan checkpoint/worker pair."""
    checkpoints = sorted(path for path in (root / 'miniprot_predictions').glob('*') if not path.name.startswith('.')
                         and (path / 'receipt.json').is_file())
    if task_index is not None:
        if not 1 <= task_index <= len(checkpoints):
            raise ValueError('Checkpoint task index out of range')
        checkpoints = [checkpoints[task_index - 1]]
    for checkpoint in checkpoints:
        worker = root / 'rescued' / checkpoint.name
        checkpoint_path, worker_path = checkpoint / 'receipt.json', worker / 'receipt.json'
        if not checkpoint_path.is_file() or not worker_path.is_file():
            continue
        checkpoint_sha, worker_sha = digest(checkpoint_path), digest(worker_path)
        before, after = json.loads(checkpoint_path.read_text()), json.loads(worker_path.read_text())
        marker = worker / 'miniprot_checkpoint.json'
        if (before['key'] != after['key'] or 'miniprot_checkpoint.json' not in after['files']
                or digest(marker) != after['files']['miniprot_checkpoint.json']
                or json.loads(marker.read_text()) != {'receipt_sha256': checkpoint_sha}):
            continue
        candidate_bytes = freed = linked = 0
        lock = root / '.locks' / ('rescued__' + checkpoint.name + '.lock')
        with exclusive_lock(lock):
            for name in sorted(SHARED & before['files'].keys() & after['files'].keys()):
                checksum = before['files'][name]
                source, target = checkpoint / name, worker / name
                if checksum != after['files'][name] or os.path.samefile(source, target):
                    continue
                candidate_bytes += target.stat().st_size
                if apply:
                    if digest(checkpoint_path) != checkpoint_sha or digest(worker_path) != worker_sha:
                        raise ValueError('Miniprot stage receipt changed during retirement')
                    size = retire_copy(source, target, checksum)
                    freed += size
                    linked += os.path.samefile(source, target)
        yield {'species': checkpoint.name, 'candidate_bytes': candidate_bytes,
               'freed_bytes': freed, 'linked_files': linked}


def main():
    parser = argparse.ArgumentParser(description=deduplicate.__doc__)
    parser.add_argument('--output', type=Path, required=True, help='Closed split-Miniprot rescue directory')
    parser.add_argument('--apply', action='store_true', help='Retire verified separate copies in a closed run')
    parser.add_argument('--task-index', type=int, help='One-based checkpoint directory index in sorted order')
    args = parser.parse_args()
    for row in deduplicate(args.output.resolve(strict=True), apply=args.apply, task_index=args.task_index):
        print(json.dumps(row, sort_keys=True), flush=True)


if __name__ == '__main__':
    main()

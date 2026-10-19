"""Scratch retirement closes its own lock before NFS removes that inode."""
import fcntl
import os
import subprocess
import time
from pathlib import Path

CORE = Path(__file__).resolve().parents[1] / 'core/gg_transcriptome_generation_core.sh'


def run(script, *args, env=None):
    return subprocess.run(['bash', '-c', 'source "$(dirname "$1")/../support/gg_util.sh"; '
        'GG_CORE_SOURCE_ONLY=1 source "$1"; ' + script,
        'scratch-test', str(CORE), *map(str, args)], text=True, capture_output=True, env=env, check=True)


def test_current_tree_is_retired_while_locked_then_closed_before_removal(tmp_path):
    root = tmp_path / 'tmp'
    scratch = root / '2_Plant_one'
    scratch.mkdir(parents=True)
    (scratch / 'result').write_text('temporary data')
    script = '''set -euo pipefail
dir_tmp=$2; transcriptome_tmp_root=$3; delete_tmp_dir=1
acquire_transcriptome_tmp_lock
rm() {
 python - "${@: -1}" <<'PY'
import os, pathlib, sys
target = str(pathlib.Path(sys.argv[1])) + '/'
for path in pathlib.Path('/proc/self/fd').iterdir():
 try:
  link = os.readlink(path)
 except OSError:
  continue
 if link.startswith(target) and link.endswith('/.gg_active.lock'):
  raise SystemExit('NFS would retain the open lock inode')
PY
 command rm "$@"
}
'''
    source = CORE.read_text()
    ending = source[source.index('if [[ ${delete_tmp_dir} -eq 1 ]]; then', source.index('### Modify below')):]
    ending = ending[:ending.index('\ncleanup_transcriptome_tmp "${transcriptome_tmp_root}" ""')]
    result = run(script + ending, scratch, root)
    assert 'Normal completion' in result.stdout
    assert not scratch.exists()
    assert list(root.iterdir()) == []


def test_retention_closes_its_lock_before_rmtree_and_preserves_active_or_other_data(tmp_path):
    root = tmp_path / 'tmp'
    expired, active, current = (root / name for name in ('1_old', '2_active', '3_current'))
    for path in (expired, active, current):
        path.mkdir(parents=True)
        (path / 'data').write_text('generated scratch')
        (path / '.gg_active.lock').touch()
        old = time.time() - 3 * 86400
        for item in (*path.iterdir(), path):
            os.utime(item, (old, old))
    unrelated = root / 'curated_input'
    unrelated.write_text('preserve')
    injection = tmp_path / 'injection'
    injection.mkdir()
    (injection / 'sitecustomize.py').write_text('''import os, pathlib, shutil
original = shutil.rmtree
def nfs_rmtree(path, *args, **kwargs):
 target = str(path) + '/'
 for fd in pathlib.Path('/proc/self/fd').iterdir():
  try:
   link = os.readlink(fd)
  except OSError:
   continue
  if link.startswith(target) and link.endswith('/.gg_active.lock'):
   raise OSError('NFS open-lock unlink would leave directory nonempty')
 return original(path, *args, **kwargs)
shutil.rmtree = nfs_rmtree
''')
    with (active / '.gg_active.lock').open('r+') as held:
        fcntl.flock(held, fcntl.LOCK_EX | fcntl.LOCK_NB)
        run('transcriptome_tmp_retention_days=1; transcriptome_tmp_max_dirs=0; '
            'transcriptome_tmp_max_bytes=0; transcriptome_tmp_max_files=0; '
            'cleanup_transcriptome_tmp "$2" 3_current', root,
            env={**os.environ, 'PYTHONPATH': str(injection)})
    assert not expired.exists()
    assert (active / 'data').read_text() == 'generated scratch'
    assert (current / 'data').read_text() == 'generated scratch'
    assert unrelated.read_text() == 'preserve'
    assert not list(root.glob('.gg-tmp-retire-*'))


def test_unsafe_current_scratch_is_preserved(tmp_path):
    outside = tmp_path / 'input'
    outside.mkdir()
    (outside / 'keep').write_text('curated')
    scratch = tmp_path / '2_Plant_one'
    scratch.symlink_to(outside, target_is_directory=True)
    result = subprocess.run(['bash', '-c', 'GG_CORE_SOURCE_ONLY=1 source "$1"; '
        'dir_tmp=$2; transcriptome_tmp_root=$3; retire_transcriptome_tmp',
        'scratch-test', str(CORE), str(scratch), str(tmp_path)], text=True, capture_output=True)
    assert result.returncode != 0
    assert (outside / 'keep').read_text() == 'curated'

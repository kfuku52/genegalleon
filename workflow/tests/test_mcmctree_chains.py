"""Real PAML and advisory-diagnostic integration tests; run in GeneGalleon."""
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / 'support' / 'mcmctree_chains.py'
spec = importlib.util.spec_from_file_location('chains', SCRIPT)
chains = importlib.util.module_from_spec(spec)
spec.loader.exec_module(chains)


def fixture(directory, count=30):
    directory.mkdir()
    (directory / 'tree.nwk').write_text("3 1\n((a,b)'B(0.1,0.2,0.025,0.025)',c)'U(0.3,0.025)';\n")
    (directory / 'dummy.phy').write_text('3 4\na    ACGT\nb    ACGT\nc    ACGT\n')
    (directory / 'input.ctl').write_text(f'''seed = 1
seqfile = dummy.phy
treefile = tree.nwk
outfile = out
mcmcfile = mcmc.txt
ndata = 1
seqtype = 0
usedata = 0
clock = 2
RootAge = <0.3
model = 4
alpha = 0.5
ncatG = 5
cleandata = 1
BDparas = 1 1 0.5 M
kappa_gamma = 6 2
alpha_gamma = 1 10
rgene_gamma = 2 20 1
sigma2_gamma = 1 10 1
burnin = 20
sampfreq = 2
nsample = {count}
''')


def test_samples_reject_corruption(tmp_path):
    path = tmp_path / 'samples'
    for text in ('Gen t\n1 0.1\n2\n', 'Gen t\n1 nan\n', 'Gen t\n1 1\n1 2\n'):
        path.write_text(text)
        with pytest.raises(ValueError):
            chains.samples(path)


def test_ctl_replaces_duplicate_values():
    text = chains.control(' seed = -1\nseed = 8\nprint = 0\n', seed=12, print=1)
    assert text.count('seed =') == 1
    assert 'seed = 12' in text
    assert 'print = 1' in text


@pytest.mark.skipif(shutil.which('mcmctree') is None, reason='GeneGalleon runtime required')
def test_parallel_paml_pools_samples_and_resumes_completed_chains(tmp_path):
    template = tmp_path / 'template'
    fixture(template)
    status = tmp_path / 'status.json'
    output = tmp_path / 'summary.out'
    command = [sys.executable, str(SCRIPT), '--template', str(template), '--control', 'input.ctl',
               '--store', str(tmp_path / 'runs'), '--output', str(output), '--status', str(status)]
    result = subprocess.run(command, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads(status.read_text())
    assert report['completed_chains'] == 4
    assert report['pooled_samples'] >= 120
    assert report['state'] != 'passed'  # deliberately far too short
    assert 'WARNING' in result.stderr
    assert 'FigTree' in output.read_text()
    from nwkit.convert import convert_tree_text
    public = convert_tree_text(output.read_text(), target='figtree', time_factor=1000)
    assert public.startswith('#NEXUS') and '95%HPD' in public
    run = Path(report['run'])
    attempts = sorted(run.glob('chain-*/attempt-*'))
    assert len(attempts) == 4
    assert len({(p / 'mcmc.txt').read_bytes() for p in attempts}) == 4
    for i, path in enumerate(attempts):
        assert f'seed = {1729 + i}' in (path / 'run.ctl').read_text()
        assert (path / 'stdout.log').is_file()
        assert (path / 'stderr.log').is_file()
    result = subprocess.run(command, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert sorted(run.glob('chain-*/attempt-*')) == attempts
    # A corrupted completed chain is rerun, without overwriting old evidence.
    (attempts[0] / 'mcmc.txt').write_text('corrupt')
    result = subprocess.run(command, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert len(list(run.glob('chain-*/attempt-*'))) == 5
    assert (attempts[0] / 'mcmc.txt').read_text() == 'corrupt'


def test_diagnostics_detect_shift_and_constant(tmp_path):
    import numpy as np
    rng = np.random.default_rng(871)
    paths = []
    for i in range(4):
        path = tmp_path / f'chain{i}.txt'
        with path.open('w') as out:
            out.write('Gen stable shifted\n')
            for j in range(1500):
                out.write(f'{j+1} {rng.normal()} {rng.normal()+i*4}\n')
        paths.append(path)
    report = chains.diagnose(paths, tmp_path)
    assert report['state'] == 'unconverged', (tmp_path / 'diagnostics.log').read_text()
    assert 'shifted' in report['failed_parameters']
    assert 'stable' not in report['failed_parameters']
    for path in paths:
        with path.open('w') as out:
            out.write('Gen stable\n')
            for j in range(1500):
                out.write(f'{j+1} {rng.normal()}\n')
    assert chains.diagnose(paths, tmp_path)['state'] == 'passed'
    for path in paths:
        path.write_text('Gen constant\n' + ''.join(f'{j+1} 1\n' for j in range(20)))
    assert chains.diagnose(paths, tmp_path)['state'] == 'diagnostic_failed'


@pytest.mark.parametrize('fail_one', [False, True])
def test_four_processes_overlap_and_partial_failure_only_warns(tmp_path, monkeypatch, fail_one):
    template = tmp_path / 'template'
    fixture(template, count=10)
    binary_dir = tmp_path / 'bin'
    binary_dir.mkdir()
    marker = tmp_path / 'started'
    marker.mkdir()
    binary = binary_dir / 'mcmctree'
    binary.write_text('''#!/usr/bin/env python
import os, pathlib, re, sys, time
ctl = pathlib.Path(sys.argv[1]).read_text()
def value(key):
    return re.search(r'^' + key + r' = (.+)$', ctl, re.M)[1]
if value('print') != '-1':
    marker = pathlib.Path(os.environ['MCMC_TEST_MARKERS'])
    (marker / value('seed')).touch()
    deadline = time.monotonic() + 5
    while len(list(marker.iterdir())) < 4:
        if time.monotonic() > deadline:
            sys.exit(9)  # sequential execution would deadlock at this barrier
        time.sleep(0.01)
    if os.environ['MCMC_TEST_FAIL'] == '1' and value('seed') == '1729':
        sys.exit(8)
    pathlib.Path('mcmc.txt').write_text('Gen t\\n' + ''.join(f'{i*2} {i/100}\\n' for i in range(1,11)))
pathlib.Path(value('outfile')).write_text('Species tree for FigTree\\n(a:1,b:1);\\n')
''')
    binary.chmod(0o755)
    monkeypatch.setenv('PATH', str(binary_dir) + ':' + __import__('os').environ['PATH'])
    monkeypatch.setenv('MCMC_TEST_MARKERS', str(marker))
    monkeypatch.setenv('MCMC_TEST_FAIL', '1' if fail_one else '0')
    status = tmp_path / 'status.json'
    result = subprocess.run([sys.executable, str(SCRIPT), '--template', str(template),
                             '--control', 'input.ctl', '--store', str(tmp_path / 'runs'),
                             '--output', str(tmp_path / 'out'), '--status', str(status)],
                            capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, result.stderr
    assert len(list(marker.iterdir())) == 4
    report = json.loads(status.read_text())
    assert report['completed_chains'] == (3 if fail_one else 4)
    assert 'WARNING' in result.stderr
    if fail_one:
        assert report['state'] == 'incomplete'
        assert report['pooled_samples'] == 30


def test_tree_status_rejects_replaced_tree_and_warns_without_failure(tmp_path):
    tool = SCRIPT.with_name('mcmctree_status.py')
    tree = tmp_path / 'tree.nwk'
    tree.write_text('(a:1,b:1);')
    status = tmp_path / 'status.json'
    status.write_text(json.dumps({'state': 'passed'}))
    subprocess.run([sys.executable, str(tool), 'bind', '--tree', str(tree), '--status', str(status)], check=True)
    command = [sys.executable, str(tool), 'warn', '--tree', str(tree), '--status', str(status)]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0 and not result.stderr
    tree.write_text('(a:2,b:2);')
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0
    assert 'legacy_unverified' in result.stderr
    status.unlink()
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0 and 'WARNING' in result.stderr


@pytest.mark.parametrize('content', ['not json', '{}', '{"attempt":"../../escape","hashes":{}}',
                                   '{"attempt":"attempt-old","state":"completed","seed":1729,"chain":1,"hashes":{}}'])
def test_corrupt_receipt_cannot_reuse_an_attempt(tmp_path, content):
    receipt = tmp_path / 'completed.json'
    receipt.write_text(content)
    assert chains.reusable_attempt(tmp_path, receipt, 1, 1729) is None


@pytest.mark.parametrize('text', ['nsample=1.5', 'nsample=-1', 'nsample=0', 'nsample=1\nnsample=2', ''])
def test_control_integer_rejects_missing_ambiguous_or_invalid_values(text):
    with pytest.raises(ValueError):
        chains.positive_control_integer(text, 'nsample')


def fake_engine(tmp_path, monkeypatch, *, hang=False):
    import os
    binary_dir = tmp_path / 'bin'
    binary_dir.mkdir()
    marker = tmp_path / 'started'
    marker.mkdir()
    binary = binary_dir / 'mcmctree'
    binary.write_text('''#!/usr/bin/env python
import os, pathlib, re, signal, sys, time
ctl = pathlib.Path(sys.argv[1]).read_text()
def value(key):
    return re.search(r'^' + key + r' = (.+)$', ctl, re.M)[1]
if value('print') == '-1':
    sys.exit(9)
if os.environ.get('MCMC_HANG') == '1':
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    (pathlib.Path(os.environ['MCMC_MARKER']) / value('seed')).write_text(str(os.getpid()))
    time.sleep(60)
pathlib.Path('mcmc.txt').write_text('Gen t\\n' + ''.join(f'{i*2} {i/100}\\n' for i in range(1,11)))
pathlib.Path(value('outfile')).write_text('Species tree for FigTree\\n(a:1,b:1);\\n')
''')
    binary.chmod(0o755)
    monkeypatch.setenv('PATH', str(binary_dir) + ':' + os.environ['PATH'])
    monkeypatch.setenv('MCMC_HANG', '1' if hang else '0')
    monkeypatch.setenv('MCMC_MARKER', str(marker))
    return marker


def test_summary_failure_cannot_leave_passed_state_or_replace_output(tmp_path, monkeypatch):
    import argparse
    template = tmp_path / 'template'
    fixture(template, count=10)
    fake_engine(tmp_path, monkeypatch)
    monkeypatch.setattr(chains, 'diagnose', lambda *a: {'state': 'passed'})
    status = tmp_path / 'status.json'
    output = tmp_path / 'out'
    output.write_text('previous result')
    args = argparse.Namespace(template=template, control='input.ctl', seed=1729, jobs=4,
                              time_factor='1', store=tmp_path / 'runs', status=status, output=output)
    with pytest.raises(ValueError, match='summary failed'):
        chains.run(args)
    report = json.loads(status.read_text())
    assert report['state'] == 'failed'
    assert json.loads((Path(report['run']) / 'status.json').read_text())['state'] == 'failed'
    assert output.read_text() == 'previous result'
    assert len(list(args.store.glob('*/chain-*/attempt-*/mcmc.txt'))) == 4


def test_termination_stops_children_and_retains_interrupted_state(tmp_path, monkeypatch):
    import signal
    import time
    template = tmp_path / 'template'
    fixture(template, count=10)
    marker = fake_engine(tmp_path, monkeypatch, hang=True)
    status = tmp_path / 'status.json'
    with (tmp_path / 'log').open('w') as log:
        process = subprocess.Popen([sys.executable, str(SCRIPT), '--template', str(template),
                                    '--control', 'input.ctl', '--store', str(tmp_path / 'runs'),
                                    '--output', str(tmp_path / 'out'), '--status', str(status)], stdout=log, stderr=log)
        try:
            deadline = time.monotonic() + 10
            while len(list(marker.iterdir())) < 4 and time.monotonic() < deadline:
                time.sleep(0.05)
            assert len(list(marker.iterdir())) == 4
            process.send_signal(signal.SIGTERM)
            assert process.wait(timeout=10) != 0
            report = json.loads(status.read_text())
            assert report['state'] == 'interrupted'
            states = list(Path(report['run']).glob('chain-*/attempt-*/state.json'))
            assert len(states) == 4
            assert all(json.loads(p.read_text())['state'] == 'interrupted' for p in states)
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()
            for pid_file in marker.iterdir():
                try:
                    __import__('os').kill(int(pid_file.read_text()), signal.SIGKILL)
                except ProcessLookupError:
                    pass


def test_cached_conversion_refreshes_advisory_state_without_rewriting_tree(tmp_path):
    tool = SCRIPT.with_name('mcmctree_status.py')
    source, target, status = [tmp_path / name for name in ('source', 'target', 'status')]
    source.write_text('(a:1,b:1);')
    target.write_text('(a:1,b:1)root;')
    def command(action, *extra):
        subprocess.run([sys.executable, str(tool), action, '--tree', str(source),
                        '--status', str(status), *extra], check=True)
    status.write_text('{"state":"passed"}')
    command('bind')
    command('copy', '--target', str(target))
    status.write_text('{"state":"unconverged"}')
    command('bind')
    command('copy', '--target', str(target), '--cached')
    sidecar = Path(str(target) + '.convergence.json')
    assert json.loads(sidecar.read_text())['state'] == 'unconverged'
    assert target.read_text() == '(a:1,b:1)root;'
    target.write_text('(x:1,y:1);')
    previous = sidecar.read_bytes()
    command('copy', '--target', str(target), '--cached')
    assert sidecar.read_bytes() == previous  # do not certify an unrelated replacement

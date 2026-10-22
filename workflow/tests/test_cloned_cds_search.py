import importlib
import json
import subprocess
import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
resources = importlib.import_module('cloned_cds_resources')
search = importlib.import_module('cloned_cds_search')


def test_allocation_planner_uses_smallest_limit_and_thread_overhead():
    plan = resources.memory_plan(4, environ={'GG_MEM_TOTAL_GB': '32', 'GG_MEM_TOOL_GB': '28',
        'SLURM_MEM_PER_CPU': '4096'}, limits={'cgroup': 12 * resources.GIB})
    assert plan['allocation_bytes'] == 12 * resources.GIB
    assert plan['budget_bytes'] == 9 * resources.GIB
    assert plan['target_bases'] == (9 * resources.GIB - (512 + 4 * 64) * 1024 ** 2) // 16
    capped = resources.memory_plan(4, 1000, environ={'SLURM_MEM_PER_NODE': '32768'}, limits={})
    assert capped['target_bases'] == 1000
    assert resources.memory_plan(1, environ={}, limits={})['target_bases'] == 250_000_000


@pytest.mark.parametrize('value', ['-1', 'nan', 'Infinity', '0', 'nonsense'])
def test_invalid_explicit_memory_never_falls_back_to_host_ram(value):
    with pytest.raises(ValueError, match='Invalid allocated memory'):
        resources.memory_plan(2, environ={'GG_MEM_TOTAL_GB': value}, limits={})


@pytest.mark.parametrize('controller,base,filename', [('', '', 'memory.max'),
    ('cpu,memory', 'memory', 'memory.limit_in_bytes')])
def test_cgroup_hierarchy_keeps_parent_and_namespace_root_limits(tmp_path, controller, base, filename):
    membership = tmp_path / 'membership'
    membership.write_text(f'0:{controller}:/slurm/task\n')
    root = tmp_path / 'cgroup'
    leaf = root / base / 'slurm/task'
    leaf.mkdir(parents=True)
    (root / base / filename).write_text(str(32 * resources.GIB))
    (leaf.parent / filename).write_text(str(8 * resources.GIB))
    (leaf / filename).write_text('max' if not controller else str(2 ** 63 - 4096))
    limits = resources.cgroup_limits(root, membership)
    assert min(limits.values()) == 8 * resources.GIB
    membership.write_text(f'0:{controller}:/../../outside\n')
    assert min(resources.cgroup_limits(root, membership).values()) == 32 * resources.GIB


def test_whole_contig_cannot_silently_exceed_known_budget():
    plan = resources.memory_plan(1, environ={'GG_MEM_TOTAL_GB': '1'}, limits={})
    with pytest.raises(ValueError, match='Contig huge.*larger allocation'):
        resources.groups_for([['huge', 100_000_000]], plan)
    # A target size is soft at contig boundaries, but the estimated budget is hard.
    plan['target_bases'] = 10
    assert resources.groups_for([['chr', 100]], plan) == [[['chr', 100]]]
    with pytest.raises(ValueError, match='too small'):
        resources.memory_plan(4, environ={'GG_MEM_TOTAL_GB': '.5'}, limits={})


def test_interrupted_stage_resumes_verified_parts_and_publishes_same_cds_gff(tmp_path, monkeypatch):
    import rescue_gene_models as rescue
    from cloned_cds_models import apply_to_task, attach
    from test_cloned_cds_models import fixture
    task, _sequence, paths, row = fixture(tmp_path / 'inputs')
    with paths['genome'].open('a') as handle:
        handle.write('>chr2\n' + 'C' * 1500 + '\n>chr3\n' + 'C' * 1500 + '\n')
    attach(task, {**row, 'cloned_cds_index_batch_bases': '1000'}, tmp_path)
    baseline = dict(task)
    apply_to_task(baseline, tmp_path / 'baseline')
    run = rescue.run
    calls = []

    def interrupted(*args, **kwargs):
        calls.append(args[0])
        if len(calls) == 2:
            raise RuntimeError('simulated interruption')
        return run(*args, **kwargs)
    monkeypatch.setattr(rescue, 'run', interrupted)
    root = tmp_path / 'retry'
    with pytest.raises(RuntimeError, match='interruption'):
        apply_to_task(dict(task), root)
    assert len(list((root / '.cloned_search').glob('*/parts/*/receipt.json'))) == 1
    assert not list(root.glob('Test_species/*/receipt.json'))
    assert task['cds_path'] == paths['cds']
    calls.clear()

    def counted(*args, **kwargs):
        calls.append(args[0])
        return run(*args, **kwargs)
    monkeypatch.setattr(rescue, 'run', counted)
    resumed = dict(task)
    apply_to_task(resumed, root)
    assert len(calls) == 2
    for role in ('cds', 'gff'):
        assert resumed[role + '_path'].read_bytes() == baseline[role + '_path'].read_bytes()
    batches = json.loads((Path(resumed['cloned_model_audit']).parent / 'index_batches.json').read_text())
    assert batches['complete'] is True and len(batches['batches']) == 3
    assert all(row['receipt_sha256'] for row in batches['batches'])
    assert not list((root / '.cloned_search').rglob('genome.fa'))


def test_smaller_allocation_refines_only_unfinished_groups_and_corrupt_parts_rerun(tmp_path, monkeypatch):
    import rescue_gene_models as rescue
    from gene_model_catalog import indexed_genome
    from test_cloned_cds_models import fixture
    _task, _sequence, paths, _row = fixture(tmp_path / 'inputs')
    with paths['genome'].open('a') as handle:
        for i in (2, 3, 4):
            handle.write(f'>chr{i}\n' + 'C' * 1500 + '\n')
    query = tmp_path / 'protein.fa'
    query.write_text('>q\nM' + 'A' * 200 + '\n')
    plan = resources.memory_plan(1, environ={}, limits={})
    plan['target_bases'] = 3000
    monkeypatch.setattr(resources, 'memory_plan', lambda *_: dict(plan))
    root = tmp_path / 'checkpoints'
    out = tmp_path / 'out'
    out.mkdir()
    run = rescue.run
    calls = []

    def counted(*args, **kwargs):
        calls.append(args[0])
        return run(*args, **kwargs)
    monkeypatch.setattr(rescue, 'run', counted)
    with indexed_genome(paths['genome']) as genome:
        generator = search.alignments(genome, query, out, 1, 1, checkpoint_root=root)
        next(generator)
        generator.close()  # abrupt stop after the first verified group
    plan.update(target_bases=1500, budget_bytes=plan['overhead_bytes'] + 2000 * 16)
    calls.clear()
    with indexed_genome(paths['genome']) as genome:
        completed = list(search.alignments(genome, query, out, 1, 1, checkpoint_root=root))
    assert len(completed) == 3 and len(calls) == 2
    state = json.loads((root / 'plan.json').read_text())
    assert len(state['groups'][0]) == 2  # completed group survives the lower budget
    part = next(root.glob('parts/*/models.gff'))
    part.write_text('corrupt\n')
    calls.clear()
    # Raise the budget so the corrupt old multi-contig part can be rebuilt intact.
    plan.update(target_bases=3000, budget_bytes=None)
    with indexed_genome(paths['genome']) as genome:
        list(search.alignments(genome, query, out, 1, 1, checkpoint_root=root))
    assert len(calls) == 1
    query.write_text('>q\nM' + 'A' * 201 + '\n')
    with indexed_genome(paths['genome']) as genome, pytest.raises(ValueError, match='complete reference'):
        list(search.alignments(genome, query, out, 1, 1, checkpoint_root=root))


def test_sigkill_releases_locks_and_retains_completed_batch_outside_temporary_stage(tmp_path):
    from cloned_cds_models import apply_to_task, attach
    from test_cloned_cds_models import fixture
    task, _sequence, paths, row = fixture(tmp_path / 'inputs')
    with paths['genome'].open('a') as handle:
        handle.write('>chr2\n' + 'C' * 1500 + '\n')
    attach(task, {**row, 'cloned_cds_index_batch_bases': '1000'}, tmp_path)
    request = tmp_path / 'task.json'
    request.write_text(json.dumps(task, default=str))
    root, ready = tmp_path / 'retry', tmp_path / 'ready'
    code = '''
import json,sys,time
from pathlib import Path
sys.path.insert(0,sys.argv[1])
import rescue_gene_models as rescue
from cloned_cds_models import apply_to_task
task=json.loads(Path(sys.argv[2]).read_text())
run=rescue.run
calls=0
def interrupted(*args,**kwargs):
    global calls
    calls+=1
    if calls==2:
        Path(sys.argv[4]).write_text('ready')
        time.sleep(60)
    return run(*args,**kwargs)
rescue.run=interrupted
apply_to_task(task,Path(sys.argv[3]))
'''
    with (tmp_path / 'child.log').open('w') as log:
        child = subprocess.Popen([sys.executable, '-c', code, str(Path(__file__).resolve().parents[1] / 'support'),
                                  str(request), str(root), str(ready)], stdout=log, stderr=log)
        try:
            deadline = time.monotonic() + 20
            while not ready.exists() and child.poll() is None and time.monotonic() < deadline:
                time.sleep(.05)
            assert ready.exists(), (tmp_path / 'child.log').read_text()
        finally:
            child.kill()
            child.wait(timeout=5)
    assert len(list((root / '.cloned_search').glob('*/parts/*/receipt.json'))) == 1
    assert not list(root.glob('Test_species/*/receipt.json'))
    apply_to_task(task, root)
    assert len(list((root / '.cloned_search').glob('*/parts/*/receipt.json'))) == 2
    state = json.loads((Path(task['cloned_model_audit']).parent / 'index_batches.json').read_text())
    assert state['complete'] is True

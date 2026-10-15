import errno
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmarks'))
spec=importlib.util.spec_from_file_location('pipeline_optimization_benchmark',Path(__file__).resolve().parents[1]/'benchmarks/benchmark_annevo_pipeline_optimizations.py')
bm=importlib.util.module_from_spec(spec)
spec.loader.exec_module(bm)


def test_all_configurations_keep_plain_baseline_and_isolate_incremental_changes():
    assert bm.CONFIGURATIONS['baseline']==dict(order='input',compression='none',shared_jit=False)
    assert bm.CONFIGURATIONS['order_only']==dict(order='largest-first',compression='none',shared_jit=False)
    assert bm.CONFIGURATIONS['order_lzf']==dict(order='largest-first',compression='lzf',shared_jit=False)
    assert bm.CONFIGURATIONS['optimized']==dict(order='largest-first',compression='lzf',shared_jit=True)


def test_process_memory_counts_only_descendants_and_uses_proportional_memory(tmp_path,monkeypatch):
    for pid,parent,rss,pss in [(10,1,100,40),(20,10,200,60),(30,20,300,80),(40,1,10000,10000)]:
        directory=tmp_path/str(pid)
        directory.mkdir()
        (directory/'status').write_text('PPid: '+str(parent)+'\nVmRSS: '+str(rss)+' kB\n')
        (directory/'smaps_rollup').write_text('Pss: '+str(pss)+' kB\n')
    actual_path=bm.Path
    def redirected(value):
        value=str(value)
        return actual_path(str(tmp_path)+value[5:]) if value=='/proc' or value.startswith('/proc/') else actual_path(value)
    monkeypatch.setattr(bm,'Path',redirected)
    sample=bm.process_memory(10)
    assert sample['RSS_KiB']==600 and sample['PSS_KiB']==180 and sample['processes']==3


@pytest.mark.parametrize('error_number', [errno.ENOENT, errno.ESRCH, errno.EACCES])
def test_process_exit_race_is_skipped_but_permission_errors_fail(tmp_path, monkeypatch, error_number):
    for pid, parent in [(10, 1), (20, 10)]:
        directory = tmp_path/str(pid)
        directory.mkdir()
        (directory/'status').write_text('PPid: '+str(parent)+'\nVmRSS: 100 kB\n')
        (directory/'smaps_rollup').write_text('Pss: 40 kB\n')
    actual_path, original_read = bm.Path, Path.read_text
    def redirected(value):
        value = str(value)
        return actual_path(str(tmp_path)+value[5:]) if value == '/proc' or value.startswith('/proc/') else actual_path(value)
    def read(path, *a, **kw):
        if path == tmp_path/'20/smaps_rollup':
            raise OSError(error_number, 'controlled proc read failure')
        return original_read(path, *a, **kw)
    monkeypatch.setattr(bm, 'Path', redirected)
    monkeypatch.setattr(Path, 'read_text', read)
    if error_number == errno.EACCES:
        with pytest.raises(PermissionError):
            bm.process_memory(10)
    else:
        sample = bm.process_memory(10)
        assert sample['RSS_KiB'] == 100 and sample['PSS_KiB'] == 40 and sample['processes'] == 1


@pytest.mark.parametrize('failure', ['GPU', 'memory', 'empty', 'prediction'])
def test_failed_trial_retains_sampling_diagnostics_and_never_reports_success(tmp_path, monkeypatch, failure):
    fixture = tmp_path/'genome.fa'
    fixture.write_text('>chr1\nACGT\n')
    inputs = tmp_path/'input.tsv'
    inputs.write_text('controlled fixture\n')
    monkeypatch.setattr(bm.pipeline, 'read_inputs', lambda *a: ([{'genome': str(fixture), 'receipt': str(fixture)}], bm.file_sha(inputs)))
    monkeypatch.setattr(bm, 'probability_proof', lambda *a: {})
    monkeypatch.setattr(bm, 'shared_cache_preflight', lambda *a: {})
    class GPU:
        rows = [] if failure == 'empty' else [{'memory_MiB': 1}]
        errors = ['GPU query failed'] if failure == 'GPU' else []
        def __enter__(self):
            return self
        def __exit__(self, *a):
            pass
    class Memory(GPU):
        rows = [{'RSS_KiB': 100, 'PSS_KiB': 40}]
        errors = ['PermissionError: controlled failure'] if failure == 'memory' else []
    def run(args):
        if failure == 'prediction':
            raise RuntimeError('controlled native prediction failure')
        return {'species': {}}
    monkeypatch.setattr(bm, 'GPUSamples', GPU)
    monkeypatch.setattr(bm, 'MemorySamples', Memory)
    monkeypatch.setattr(bm.pipeline, 'run', run)
    args = SimpleNamespace(repeats=2, output=tmp_path/'output', inputs=inputs,
                           probability_proof_genome=fixture, tmp_root=tmp_path, max_tmp_gib=1)
    with pytest.raises(RuntimeError, match='controlled native prediction failure' if failure == 'prediction' else 'resource_samples.json'):
        bm.benchmark(args)
    saved = json.loads((args.output/'0_baseline/resource_samples.json').read_text())
    assert saved['GPU_sample_errors'] == GPU.errors
    assert saved['memory_sample_errors'] == Memory.errors
    assert saved['GPU_samples'] == GPU.rows and saved['memory_samples'] == Memory.rows
    assert not (args.output/'comparison.json').exists()


@pytest.mark.parametrize('problem', [None, 'no_reuse', 'changed_GFF'])
def test_preflight_requires_actual_cache_reuse_and_same_GFF(tmp_path, monkeypatch, problem):
    cache = tmp_path/'local-cache'
    calls = []
    def sampled(inputs, directory, options, config, shared_cache):
        assert shared_cache == cache and config == bm.CONFIGURATIONS['optimized']
        directory.mkdir()
        cache.mkdir(exist_ok=True, mode=0o700)
        (cache/'kernel.nbc').write_bytes(b'controlled compiled-cache fixture')
        logs = directory/'fixture_probe/logs'
        logs.mkdir(parents=True)
        loaded = bool(calls) and problem != 'no_reuse'
        (logs/'decoding.stdout.txt').write_text('[cache] data loaded\n' if loaded else '[cache] data saved\n')
        gff = directory/'probe.gff3'
        end = 20 if calls and problem == 'changed_GFF' else 10
        gff.write_text('chr1\tfixture\tgene\t1\t'+str(end)+'\t.\t+\t.\tID=gene1\n')
        calls.append(dict(inputs=inputs, cache=str(cache)))
        return {'species': {'fixture_probe': {'gff': str(gff)}}}, 1.0, None, None
    monkeypatch.setattr(bm, 'sampled_pipeline', sampled)
    argv = ({'receipt': '/fixture/qualified.json'}, tmp_path/'genome.fa', tmp_path/'proof', cache, {})
    if problem:
        with pytest.raises(RuntimeError, match='reuse|GFF features'):
            bm.shared_cache_preflight(*argv)
        assert not (tmp_path/'proof/comparison.json').exists()
    else:
        proof = bm.shared_cache_preflight(*argv)
        assert proof['cold_warm_GFF_features_identical'] is True
        assert proof['outside_performance_timing'] and proof['separate_from_timed_cache']
        assert [r['cache_load_events'] for r in proof['results']] == [0, 1]
    assert len(calls) == 2


def test_failed_preflight_never_starts_timed_trials_and_cleans_local_cache(tmp_path, monkeypatch):
    inputs, genome = tmp_path/'input.tsv', tmp_path/'fixture.fa'
    inputs.write_text('controlled inputs\n')
    genome.write_text('>chr1\nACGT\n')
    monkeypatch.setattr(bm.pipeline, 'read_inputs', lambda *a: ([{'genome': str(genome), 'receipt': str(genome)}], bm.file_sha(inputs)))
    monkeypatch.setattr(bm, 'probability_proof', lambda *a: {})
    seen = []
    def preflight(row, genome, root, cache, options):
        assert cache.parent.parent == tmp_path and cache.parent.is_dir()
        assert not cache.is_relative_to(tmp_path/'output')
        seen.append(cache.parent)
        raise ValueError('controlled private-cache incompatibility')
    monkeypatch.setattr(bm, 'shared_cache_preflight', preflight)
    monkeypatch.setattr(bm, 'benchmark_trials', lambda *a: pytest.fail('Timed trials must follow successful preflight'))
    monkeypatch.setenv('NUMBA_DEBUG_CACHE', 'original-value')
    args = SimpleNamespace(repeats=2, output=tmp_path/'output', inputs=inputs,
                           probability_proof_genome=genome, tmp_root=tmp_path, max_tmp_gib=1)
    with pytest.raises(ValueError, match='private-cache incompatibility'):
        bm.benchmark(args)
    assert seen and not seen[0].exists()
    assert bm.os.environ['NUMBA_DEBUG_CACHE'] == 'original-value'
    assert not (args.output/'0_baseline').exists()

import importlib.util
import sys
from pathlib import Path

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

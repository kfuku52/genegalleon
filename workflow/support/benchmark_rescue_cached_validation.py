"""Compare frozen cached-record validation with an archived serial caller.

This is a read-only bounded benchmark, not a rescue or admission publication.
Every source model shard is verified; only a declared prefix is materialised.
Run inside the normal GeneGalleon Linux runtime.
"""
import argparse
import ast
import copy
import hashlib
import json
import os
import statistics
import threading
import time
from pathlib import Path

import rescue_cached_validation as parallel
import rescue_gene_models as rescue
import rescue_model_store as store
from gene_model_catalog import indexed_genome
from input_generation_array_state import FreshDigestBatch
from rescue_prediction_cache import _raw_prediction


def baseline_caller(source,verification):
    """Keep the previous caller, and require its decision kernels to match."""
    support = Path(source)/'workflow/support'
    current = Path(rescue.__file__).parent
    modules = ('rescue_raw_validation.py','rescue_terminal_completion.py','rescue_boundary_revalidation.py',
               'rescue_homology_policy.py','rescue_prediction_cache.py','gene_model_selection.py')
    files=[support/'rescue_gene_models.py',current/'rescue_gene_models.py',
           *[directory/name for directory in (support,current) for name in modules]]
    observed=verification.read(files)
    for name in modules:
        if observed[str(support/name)] != observed[str(current/name)]:
            raise ValueError('Benchmark decision kernel changed: '+name)
    original = (support/'rescue_gene_models.py').read_bytes()
    if hashlib.sha256(original).hexdigest()!=observed[str(support/'rescue_gene_models.py')]:
        raise ValueError('Benchmark archived source changed during reading')
    tree = ast.parse(original)
    present = ast.parse(Path(rescue.__file__).read_bytes())
    for name in ('validate_model','check_interval'):
        old = next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name)
        new = next(n for n in present.body if isinstance(n,ast.FunctionDef) and n.name==name)
        if ast.dump(old,include_attributes=False) != ast.dump(new,include_attributes=False):
            raise ValueError('Benchmark decision function changed: '+name)
    callers = [n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='checked_prediction']
    if len(callers)!=1:
        raise ValueError('Benchmark requires one archived serial caller')
    code = compile(ast.fix_missing_locations(ast.Module(body=callers,type_ignores=[])),str(support/'rescue_gene_models.py'),'exec')
    verification.check()
    return code,hashlib.sha256(original).hexdigest()


def rss_tree():
    todo,seen,total=[os.getpid()],set(),0
    while todo:
        pid=todo.pop()
        if pid in seen:
            continue
        seen.add(pid)
        p=Path('/proc')/str(pid)
        try:
            for line in (p/'status').read_text().splitlines():
                if line.startswith('VmRSS:'):
                    total+=int(line.split()[1])*1024
            for task in (p/'task').iterdir():
                todo.extend(map(int,(task/'children').read_text().split()))
        except (FileNotFoundError,ProcessLookupError):
            pass
    return total


def trial(rows,genome,code,params,proteins,cpus,baseline):
    stop=threading.Event()
    peak=[rss_tree()]
    def sample():
        while not stop.wait(.02):
            peak[0]=max(peak[0],rss_tree())
    observer=threading.Thread(target=sample,daemon=True)
    memo=rescue.RawValidationMemo(rescue.validate_model,genome,code,params,
                                 candidate_homology=rescue.homology_problems)
    context={**vars(rescue),'raw_validation':memo,'genome':genome,'source':{'genetic_code':code},
             'params':params,'proteins':proteins}
    exec(baseline,context)
    serial=context['checked_prediction']
    audit={}
    iterator=(map(serial,copy.deepcopy(rows)) if cpus==1 else
              parallel.validate_cached(iter(copy.deepcopy(rows)),genome,code,params,proteins,cpus,serial,audit))
    digest=hashlib.sha256()
    count=0
    observer.start()
    started=time.perf_counter()
    try:
        for row in iterator:
            digest.update(json.dumps(row,ensure_ascii=False,allow_nan=False,separators=(',',':')).encode())
            digest.update(b'\n')
            count+=1
    finally:
        elapsed=time.perf_counter()-started
        stop.set()
        observer.join()
    return dict(cpus=cpus,records=count,wall_seconds=elapsed,peak_rss_sum_bytes=peak[0],
                ordered_full_record_sha256=digest.hexdigest(),parallel_audit=audit)


def benchmark(args):
    if not Path('/proc/self/status').is_file() or args.cpus<2 or args.records<512 or args.repeats<1:
        raise ValueError('Benchmark requires a Linux runtime, >=2 CPUs, >=512 records and positive repeats')
    if args.output.exists():
        raise ValueError('Benchmark output already exists')
    root=args.rescue_output.resolve(strict=True)
    if args.output.resolve().is_relative_to(root) or args.output.resolve().is_relative_to(args.baseline_source.resolve()):
        raise ValueError('Benchmark output must be separate from frozen scientific inputs/source')
    source_proof=FreshDigestBatch()
    baseline,baseline_sha=baseline_caller(args.baseline_source,source_proof)
    plan_sha=source_proof.read([root/'plan.json'])[str(root/'plan.json')]
    plan=json.loads((root/'plan.json').read_text())
    name=args.species
    if name not in plan['species']:
        raise ValueError('Benchmark species absent from frozen plan')
    directory=root/'rescued'/name
    key=store.frozen_model_store_key(directory)
    receipt=json.loads((directory/'receipt.json').read_text())
    if receipt['key']['plan']!=plan_sha or receipt['key']['species']!=name:
        raise ValueError('Benchmark worker is not bound to the producer plan')
    candidate=directory/'candidates.json'
    if source_proof.read([candidate]) != {str(candidate):receipt['files']['candidates.json']}:
        raise ValueError('Benchmark candidate metadata changed')
    regions={r['id']:r for r in json.loads(candidate.read_text())}
    seen=[0]
    def prefix(context):
        seen[0]+=1
        return seen[0]<=args.records
    rows=[]
    for _,model in store.iter_model_items(directory,context_filter=prefix,frozen_key=key,verification=source_proof):
        raw=_raw_prediction(model,regions)
        if raw is not None:
            rows.append(raw)
    if len(rows)<512:
        raise ValueError('Benchmark prefix has fewer than 512 reusable raw records')
    source=plan['request']['sources'][name]
    genome_path=source['genome']
    if source_proof.read([genome_path]) != {genome_path:plan['request']['files'][genome_path]}:
        raise ValueError('Benchmark genome changed')
    proteins={}
    for donor in [name,*plan['donors'][name]]:
        prepared=root/'prepared'/donor
        prepared_receipt=prepared/'receipt.json'
        if source_proof.read([prepared_receipt])!={str(prepared_receipt):receipt['key']['prepared'][donor]}:
            raise ValueError('Benchmark prepared receipt differs from producer worker')
        proof=json.loads((prepared/'receipt.json').read_text())
        if proof['key'] != {'plan':plan_sha,'species':donor}:
            raise ValueError('Benchmark prepared generation changed')
        path=prepared/'genes.pep'
        if source_proof.read([path]) != {str(path):proof['files']['genes.pep']}:
            raise ValueError('Benchmark donor proteins changed')
        proteins[donor]={i:s for i,_,s in rescue.fasta_records(path)}
    params=rescue.parameters_for(plan['request'],name)
    runs=[]
    with indexed_genome(genome_path) as genome:
        for repeat in range(args.repeats):
            for cpus in ([1,args.cpus] if repeat%2==0 else [args.cpus,1]):
                runs.append(trial(rows,genome,source['genetic_code'],params,proteins,cpus,baseline))
                source_proof.check()
    if len({(r['records'],r['ordered_full_record_sha256']) for r in runs})!=1:
        raise ValueError('Parallel cached validation differs from the archived serial caller')
    medians={str(cpus):statistics.median(r['wall_seconds'] for r in runs if r['cpus']==cpus) for cpus in [1,args.cpus]}
    result=dict(schema=1,species=name,requested_prefix_records=args.records,raw_records=len(rows),
        source_records=seen[0],source_worker_receipt_sha256=key['producer_receipt_sha256'],
        source_plan_sha256=plan_sha,baseline_source_sha256=baseline_sha,
        parallel_source_sha256=rescue.digest(parallel.__file__),runs=runs,median_wall_seconds=medians,
        speedup=medians['1']/medians[str(args.cpus)],full_records_and_order_identical=True,
        memory_definition='Sampled sum of main and descendant RSS, including shared pages once per process',
        scope='Declared prefix of frozen raw candidates; all source model-store bytes verified; no admission publication')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as out:
        json.dump(result,out,sort_keys=True,indent=2)
        out.write('\n')
    args.output.chmod(0o400)
    print(json.dumps({k:result[k] for k in ['raw_records','median_wall_seconds','speedup','full_records_and_order_identical']}))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rescue-output',required=True,type=Path)
    p.add_argument('--species',required=True)
    p.add_argument('--baseline-source',required=True,type=Path)
    p.add_argument('--records',type=int,default=2048)
    p.add_argument('--cpus',type=int,default=8)
    p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--output',required=True,type=Path)
    benchmark(p.parse_args())


if __name__=='__main__':
    main()

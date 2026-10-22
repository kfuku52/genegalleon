"""Synthetic admission-update memory benchmark; no research data are used."""
import argparse
import hashlib
import importlib.util
import json
import resource
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO=Path(__file__).resolve().parents[2]
def modules():
    sys.path.insert(0, str(REPO))
    sys.path.insert(0, str(REPO / 'workflow/support'))
    from workflow.support import rescue_model_store as store
    from workflow.tests.test_rescue_model_store import model, receipt
    return store, model, receipt


def record(index, body_group_size, model):
    row = model(index)
    row['sequence'] = 'ATG' + 'AAA' * 4000 + 'TAA'
    if body_group_size:
        start = 3 * (1 + index // body_group_size)
        row['cds'] = [[start, start + 12006, 0]]
    return row


def worker(args):
    store, model, _ = modules()
    source=Path(args.source)
    spec=importlib.util.spec_from_file_location('benchmark_patch',args.implementation)
    patch=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(patch)
    changes={i:{**record,'status':'unresolved','problems':['ambiguous_locus']}
             for i,record in enumerate(store.iter_models(source)) if i % args.change_every == 0}
    appended=[{**record(i,args.body_group_size,model),'status':'unresolved','problems':['ambiguous_locus']}
              for i in range(args.records,args.records+args.appended_records)]
    started=time.perf_counter()
    result=patch.patch_model_store(source,Path(args.destination),changes,appended=appended,
                                  append_partial_evidence=bool(appended))
    maximum=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(json.dumps({'wall_seconds':time.perf_counter()-started,'peak_rss_bytes':maximum if sys.platform=='darwin' else maximum*1024,
                      'output_sha256':hashlib.sha256(store._json(result)).hexdigest(),'output_files':len(result['files'])}))


def main():
    store, model, receipt = modules()
    p=argparse.ArgumentParser()
    p.add_argument('--implementation',required=True)
    p.add_argument('--baseline',help='Optional prior implementation; alternate fresh before/after workers.')
    p.add_argument('--change-every',type=int,default=1,help='Change one in every N candidates (default: all).')
    p.add_argument('--appended-records',type=int,default=0,help='Additional independently ordered candidates.')
    p.add_argument('--body-group-size',type=int,default=0,help='Consecutive records per distinct coding body; 0 shares one body.')
    p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--records',type=int,default=20000)
    p.add_argument('--source')
    p.add_argument('--destination')
    p.add_argument('--codec',choices=['auto','gzip','zstd'],default='auto')
    args=p.parse_args()
    if args.change_every < 1 or args.repeats < 1 or args.records < 1 or min(args.appended_records,args.body_group_size)<0:
        p.error('records, repeats and change-every must be positive; append/body counts must be nonnegative')
    if args.source:
        worker(args)
        return
    with tempfile.TemporaryDirectory(prefix='gg-patch-memory-') as name:
        base=Path(name)
        source=base/'source'
        def rows():
            for i in range(args.records):
                yield record(i,args.body_group_size,model)
        before=store.write_model_store(source,rows(),codec=args.codec,shard_bytes=65536)
        receipt(source)
        implementations={'after':args.implementation}
        if args.baseline:
            implementations={'before':args.baseline,**implementations}
        measurements={label:[] for label in implementations}
        for i in range(args.repeats):
            order=list(implementations)
            if i % 2:
                order.reverse()
            for label in order:
                command=[sys.executable,'-B',str(Path(__file__).resolve()),'--implementation',implementations[label],
                         '--change-every',str(args.change_every),
                         '--records',str(args.records),'--appended-records',str(args.appended_records),
                         '--body-group-size',str(args.body_group_size),
                         '--source',str(source),'--destination',str(base/(label+'-'+str(i)))]
                result=subprocess.run(command,capture_output=True,text=True,check=True)
                measurements[label].append(json.loads(result.stdout))
        assert len({x['output_sha256'] for runs in measurements.values() for x in runs})==1
        report={'records':args.records,'changed_records':len(range(0,args.records,args.change_every)),
                          'appended_records':args.appended_records,'body_group_size':args.body_group_size,
                          'sequence_nt':12006,'shards':len(before['shards']),
                          'python':sys.version,'codec':before['codec'],'implementation_sha256':hashlib.sha256(Path(args.implementation).read_bytes()).hexdigest(),
                          'measurements':measurements['after']}
        if args.baseline:
            report['baseline_sha256']=hashlib.sha256(Path(args.baseline).read_bytes()).hexdigest()
            report['comparisons']=measurements
        print(json.dumps(report))


if __name__ == '__main__':
    main()

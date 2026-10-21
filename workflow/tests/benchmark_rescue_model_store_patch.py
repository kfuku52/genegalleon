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


def worker(args):
    store, _, _ = modules()
    source=Path(args.source)
    spec=importlib.util.spec_from_file_location('benchmark_patch',args.implementation)
    patch=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(patch)
    changes={i:{**record,'status':'unresolved','problems':['ambiguous_locus']} for i,record in enumerate(store.iter_models(source))}
    started=time.perf_counter()
    result=patch.patch_model_store(source,Path(args.destination),changes)
    maximum=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(json.dumps({'wall_seconds':time.perf_counter()-started,'peak_rss_bytes':maximum if sys.platform=='darwin' else maximum*1024,
                      'output_sha256':hashlib.sha256(store._json(result)).hexdigest(),'output_files':len(result['files'])}))


def main():
    store, model, receipt = modules()
    p=argparse.ArgumentParser()
    p.add_argument('--implementation',required=True)
    p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--records',type=int,default=20000)
    p.add_argument('--source')
    p.add_argument('--destination')
    p.add_argument('--codec',choices=['auto','gzip','zstd'],default='auto')
    args=p.parse_args()
    if args.source:
        worker(args)
        return
    with tempfile.TemporaryDirectory(prefix='gg-patch-memory-') as name:
        base=Path(name)
        source=base/'source'
        def rows():
            for i in range(args.records):
                m=model(i)
                m['sequence']='ATG'+'AAA'*4000+'TAA'
                yield m
        before=store.write_model_store(source,rows(),codec=args.codec,shard_bytes=65536)
        receipt(source)
        measurements=[]
        for i in range(args.repeats):
            command=[sys.executable,'-B',str(Path(__file__).resolve()),'--implementation',args.implementation,
                     '--source',str(source),'--destination',str(base/('output-'+str(i)))]
            result=subprocess.run(command,capture_output=True,text=True,check=True)
            measurements.append(json.loads(result.stdout))
        assert len({x['output_sha256'] for x in measurements})==1
        print(json.dumps({'records':args.records,'sequence_nt':12006,'shards':len(before['shards']),
                          'python':sys.version,'codec':before['codec'],'implementation_sha256':hashlib.sha256(Path(args.implementation).read_bytes()).hexdigest(),
                          'measurements':measurements}))


if __name__ == '__main__':
    main()

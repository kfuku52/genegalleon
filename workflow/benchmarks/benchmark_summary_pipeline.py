"""Compare audited database publication on identical private family outputs."""
import argparse
import contextlib
import io
import json
import os
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from benchmark_remaining_io import db_fingerprint


def worker(args):
    sys.path.insert(0, str(args.support_root))
    import artifact_provenance as p
    root = args.worker / 'output/orthogroup'
    os.chdir(args.worker)
    for index in range(128):
        family = f'OG{index:07d}'
        alignment = root / 'alignment' / (family + '.fa')
        branch = root / 'stat_branch' / (family + '_stat.branch.tsv')
        tree = root / 'stat_tree' / (family + '_stat.tree.tsv')
        for path in (alignment, branch, tree):
            path.parent.mkdir(parents=True, exist_ok=True)
        alignment.write_bytes(b'>gene\n' + b'ACGT' * 65536)
        branch.write_text('branch_id\tnode_name\tnum_sp\tso_event\n0\tg1\t1\tL\n')
        tree.write_text('num_branch\tnum_spe\tnum_dup\tnum_sp\n1\t2\t0\t1\n')
        argv = ['--logical-root',str(root),'--workspace-root',str(args.worker),
                '--manifest',str(root / 'artifact_provenance' / (family + '.summary_statistics.json')),
                '--family-id',family,'--step','summary_statistics','--input','alignment='+str(alignment),
                '--output','branch='+str(branch),'--output','tree='+str(tree)]
        with contextlib.redirect_stdout(io.StringIO()):
            assert p.dispatch(['record',*argv]) == 0
    database = args.worker / 'result.db'
    audit = ['--logical-root',str(root),'--workspace-root',str(args.worker),'--output-tsv',str(args.worker/'audit.tsv'),
             '--mode','orthogroup','--check-csubst-branches','--progress-interval','0']
    db = ['--dbpath',str(database),'--dir_gene_family',str(root),'--dir_stat_tree',str(root/'stat_tree'),
          '--dir_stat_branch',str(root/'stat_branch'),'--ncpu','2','--row_threshold','4096']
    global_root = args.worker / 'output/.gg_global_artifacts'
    global_root.mkdir()
    record = ['--logical-root',str(global_root),'--workspace-root',str(args.worker),'--family-id','orthogroup',
              '--step','gene_family_database','--manifest',str(args.worker/'database.json'),
              '--input-gene-family-store','families='+str(root),'--output','database='+str(database)]
    begin = time.perf_counter()
    if (args.support_root/'gene_family_database_pipeline.py').exists():
        command=[sys.executable,str(args.support_root/'gene_family_database_pipeline.py')]
        command.extend('--audit='+item for item in audit)
        command.extend('--database='+item for item in db)
        command.extend('--record='+item for item in record)
        subprocess.run(command,check=True)
    else:
        assert p.dispatch(['audit',*audit]) == 0
        subprocess.run([sys.executable,str(args.support_root/'generate_orthogroup_database.py'),*db],check=True)
        assert p.dispatch(['record',*record]) == 0
    seconds=time.perf_counter()-begin
    contract=json.loads((args.worker/'database.json').read_text())
    assert contract['outputs'][0]['sha256'] == __import__('hashlib').sha256(database.read_bytes()).hexdigest()
    (args.worker/'measurement.json').write_text(json.dumps({'seconds':seconds,'output_sha256':db_fingerprint(database),
        'audit':(args.worker/'audit.tsv').read_text(),
        'contract':p.contract_comparison_payload(json.loads((args.worker/'database.json').read_text())),
        'peak_rss_kib_linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'peak_child_rss_kib_linux':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root',type=Path,default=Path(__file__).resolve().parents[1]/'support')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--worker',type=Path)
    args=parser.parse_args()
    args.support_root=args.support_root.resolve()
    if args.worker:
        return worker(args)
    samples=[]
    for i in range(4):
        with tempfile.TemporaryDirectory() as tmp, args.output.with_suffix('.log').open('a') as log:
            subprocess.run([sys.executable,str(Path(__file__).resolve()),'--support-root',str(args.support_root),
                            '--worker',tmp],check=True,stdout=log,stderr=subprocess.STDOUT)
            sample=json.loads((Path(tmp)/'measurement.json').read_text())
        if i:
            samples.append(sample)
    assert len({item['output_sha256'] for item in samples})==1
    args.output.write_text(json.dumps({'samples':samples,'median_seconds':statistics.median(s['seconds'] for s in samples)},indent=2)+'\n')
    print(statistics.median(s['seconds'] for s in samples),flush=True)


if __name__=='__main__':
    main()

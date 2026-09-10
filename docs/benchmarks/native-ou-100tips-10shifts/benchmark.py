"""Paired development benchmark; execute inside the frozen GeneGalleon image."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import sys
import time

ROOT = Path(__file__).resolve().parent

def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')

def generate():
    import numpy as np
    # Independent branch-recursion simulator: no inference-package simulation code.
    for seed in [27101, 27102, 27103]:
        nodes = []
        def build(start, count):
            node = {'name': f't{start}' if count == 1 else '', 'children': []}
            if count > 1:
                node['children'] = [build(start, count//2), build(start+count//2, count-count//2)]
            node['height'] = max((c['height'] for c in node['children']), default=-1)+1
            node['tips'] = sorted(sum((c['tips'] for c in node['children']), [])) if count > 1 else [node['name']]
            nodes.append(node)
            return node
        tree = build(0, 100)
        H = tree['height']
        def newick(node):
            if not node['children']: return node['name']
            return '(' + ','.join(newick(c)+':'+format((node['height']-c['height'])/H,'.17g') for c in node['children']) + ')'
        candidates = [n for n in nodes if 5 <= len(n['tips']) <= 7]
        rng = np.random.default_rng(seed)
        chosen = [candidates[i] for i in rng.choice(len(candidates),10,replace=False)]
        assert len(set(sum([n['tips'] for n in chosen],[]))) == sum(len(n['tips']) for n in chosen)
        multipliers = rng.uniform(.8,1.2,10)*rng.choice([-1,1],10)
        noise = rng.normal(size=len(nodes))
        node_index = {id(n):i for i,n in enumerate(nodes)}
        for strength in [2,6]:
            directory = ROOT/'data'/f'effect{strength}-seed{seed}'
            directory.mkdir(parents=True, exist_ok=False)
            optima = {id(n):float(strength*m) for n,m in zip(chosen,multipliers)}
            rows = []
            alpha = 3.0
            stationary_variance = 1 / (-math.expm1(-2*alpha))
            def simulate(n, state=0.0, mean=0.0, theta=0.0, parent_height=None):
                theta = optima.get(id(n),theta)
                if parent_height is not None:
                    length = (parent_height-n['height'])/H
                    decay = math.exp(-alpha*length)
                    state = decay*state + math.sqrt(stationary_variance*(-math.expm1(-2*alpha*length)))*noise[node_index[id(n)]]
                    mean = decay*mean+(1-decay)*theta
                for c in n['children']: simulate(c,state,mean,theta,n['height'])
                if not n['children']: rows.append([n['name'],state+mean,mean])
            simulate(tree)
            (directory/'tree.nwk').write_text(newick(tree)+';\n')
            with (directory/'traits.tsv').open('w') as f:
                writer=csv.writer(f,delimiter='\t'); writer.writerow(['taxon','x']); writer.writerows([r[:2] for r in rows])
            dump(directory/'truth.json',{'seed':seed,'effect':strength,'alpha_height':alpha,'process_tip_variance':1,'root_model':'OUfixedRoot','shift_clades':[n['tips'] for n in chosen],'optima':list(optima.values()),'tip_mean':{r[0]:r[2] for r in rows},'sha256':{name:hashlib.sha256((directory/name).read_bytes()).hexdigest() for name in ['tree.nwk','traits.tsv']}})

def measure(command, log, limit):
    start=time.perf_counter()
    pid=os.fork()
    if pid==0:
        os.setsid()
        fd=os.open(log,os.O_WRONLY|os.O_CREAT|os.O_TRUNC,0o644)
        os.dup2(fd,1);os.dup2(fd,2);os.close(fd)
        os.execvp(command[0],command)
    timed_out=False
    while True:
        found,status,usage=os.wait4(pid,os.WNOHANG)
        if found:break
        if time.perf_counter()-start > limit:
            timed_out=True
            os.killpg(pid,signal.SIGKILL)
            _,status,usage=os.wait4(pid,0)
            break
        time.sleep(.1)
    return {'wall_seconds':time.perf_counter()-start,'peak_rss_mib':usage.ru_maxrss/1024,'user_seconds':usage.ru_utime,'system_seconds':usage.ru_stime,'timed_out':timed_out,'exit_code':os.waitstatus_to_exitcode(status),'command':command}

def run():
    results=[]
    directories=sorted((ROOT/'data').iterdir())
    for index,directory in enumerate(directories):
        methods=['kfl1ou','nwkit'] if index%2==0 else ['nwkit','kfl1ou']
        for method in methods:
            output=directory/(method+'.json')
            cmd=['Rscript',str(ROOT/'kfl1ou.R'),str(directory),str(output)] if method=='kfl1ou' else ['python',str(ROOT/'native.py'),str(directory),str(output)]
            print('START',directory.name,method,flush=True)
            measured=measure(cmd,str(directory/(method+'.log')),300)
            row={'dataset':directory.name,'method':method,**measured}
            if output.exists():row['result']=json.loads(output.read_text())
            results.append(row);dump(ROOT/'results.json',results)
            print('FINISH',directory.name,method,json.dumps({k:v for k,v in measured.items() if k!='command'}),flush=True)
    dump(ROOT/'complete.json',{'runs':len(results)})

if __name__=='__main__':
    {'generate':generate,'run':run}[sys.argv[1]]()

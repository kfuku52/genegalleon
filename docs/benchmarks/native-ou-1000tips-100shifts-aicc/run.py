"""Sequential paired large-tree AICc search timing; no resampling or fixed alpha."""
import hashlib,json
from benchmark import ROOT,measure,dump
assert not (ROOT/'results.json').exists(), 'Refuse to overwrite completed timings'
results=[]
for i,directory in enumerate(sorted((ROOT/'data').iterdir())):
    for method in (['path','beam'] if i%2==0 else ['beam','path']):
        output=directory/(method+'.json')
        command=['python',str(ROOT/'native.py'),str(directory),str(output),method,'AICc']
        print('START',directory.name,method,flush=True)
        measured=measure(command,str(directory/(method+'.log')),1800)
        row={'dataset':directory.name,'method':method,'input_sha256':{p:hashlib.sha256((directory/p).read_bytes()).hexdigest() for p in ['tree.nwk','traits.tsv','truth.json']},**measured}
        if output.exists():row['result']=json.loads(output.read_text())
        results.append(row);dump(ROOT/'results.json',results)
        print('FINISH',len(results),directory.name,method,measured['exit_code'],round(measured['wall_seconds'],2),flush=True)
dump(ROOT/'complete.json',{'runs':len(results),'complete':sum(r['exit_code']==0 for r in results),'timeouts':sum(r['timed_out'] for r in results)})

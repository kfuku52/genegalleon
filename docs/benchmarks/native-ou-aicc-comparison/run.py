"""Prespecified AICc addendum on unchanged benchmark inputs, sequential runs."""
import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
PREVIOUS=ROOT.parent/'native-ou-100tips-10shifts-ic'
sys.path.insert(0,str(PREVIOUS))
from benchmark import measure
jobs=[]
nonnull=sorted((PREVIOUS/'data').iterdir())
null=sorted((ROOT.parent/'native-ou-aic-improvement/data').glob('effect0-*'))
for i,d in enumerate(nonnull):
    methods=['nwkit-beam-AICc','nwkit-path-AIC','nwkit-path-AICc','kfl1ou-AICc','kfl1ou-AIC-replay']
    methods=methods[i%5:]+methods[:i%5]
    jobs.extend((d,m) for m in methods)
for i,d in enumerate(null):
    methods=['nwkit-beam-AICc','nwkit-path-AICc','kfl1ou-AICc']
    methods=methods[i%3:]+methods[:i%3]
    jobs.extend((d,m) for m in methods)
assert not (ROOT/'results.json').exists(), 'Refuse to overwrite results'
results=[]
for d,method in jobs:
    destination=ROOT/'runs'/d.name;destination.mkdir(parents=True,exist_ok=True)
    output=destination/(method+'.json')
    if method.startswith('kfl1ou'):
        criterion='AIC' if method.endswith('replay') else 'AICc'
        command=['Rscript',str(ROOT/'kfl1ou.R'),str(d),str(output),criterion]
    else:
        _,strategy,criterion=method.split('-')
        command=['python',str(ROOT/'native.py'),str(d),str(output),strategy,criterion]
    print('START',d.name,method,flush=True)
    measured=measure(command,str(destination/(method+'.log')),1800)
    row={'dataset':d.name,'method':method,'inputs':{p:hashlib.sha256((d/p).read_bytes()).hexdigest() for p in ['tree.nwk','traits.tsv','truth.json']},**measured}
    if output.exists():row['result']=json.loads(output.read_text())
    results.append(row)
    (ROOT/'results.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
    print('FINISH',len(results),len(jobs),d.name,method,measured['exit_code'],flush=True)
(ROOT/'complete.json').write_text(json.dumps({'runs':len(results),'complete':sum(r['exit_code']==0 and r['result']['status']=='complete' for r in results)},indent=2)+'\n')

"""Repair the documented harness preflight; extend every censored run equally.

No input, inference configuration, seed or stopping rule is changed. The sole
extension is a 20-minute resource ceiling, allowing completed accuracy scoring.
"""
import json
from pathlib import Path
from benchmark import ROOT,dump,measure

initial=json.loads((ROOT/'results.json').read_text())
archive=ROOT/'first-pass-results.json'
assert not archive.exists(), 'Refuse to overwrite the initial record.'
dump(archive,initial)
results=list(initial)
existing={(r['dataset'],r['method']) for r in initial}
for index,directory in enumerate(sorted((ROOT/'data').iterdir())):
    order=['kfl1ou','nwkit'] if index%2==0 else ['nwkit','kfl1ou']
    for method in order:
        if (directory.name,method) not in existing:
            initial.append({'dataset':directory.name,'method':method,'timed_out':False,'pending':True})
            results.append(initial[-1])
for index,run in enumerate(initial):
    preflight=(run['dataset']=='effect2-seed27101' and run['method']=='kfl1ou' and 'postorder' in run.get('result',{}).get('error',''))
    if not preflight and not run['timed_out'] and not run.get('pending'):continue
    directory=ROOT/'data'/run['dataset']
    method=run['method']
    output=directory/(method+'-followup.json')
    command=['Rscript',str(ROOT/'kfl1ou.R'),str(directory),str(output)] if method=='kfl1ou' else ['python',str(ROOT/'native.py'),str(directory),str(output)]
    print('FOLLOWUP START',run['dataset'],method,flush=True)
    measured=measure(command,str(directory/(method+'-followup.log')),1200)
    reason='tree-order harness preflight correction' if preflight else ('pending run under extended ceiling' if run.get('pending') else 'unchanged inference with longer resource ceiling')
    row={'dataset':run['dataset'],'method':method,**measured,'resource_limit_seconds':1200,'followup_reason':reason}
    if output.exists():row['result']=json.loads(output.read_text())
    results[index]=row;dump(ROOT/'results.json',results)
    print('FOLLOWUP FINISH',run['dataset'],method,json.dumps({k:v for k,v in measured.items() if k!='command'}),flush=True)
dump(ROOT/'followup-complete.json',{'rows':len(results),'initial_timeouts':sum(r['timed_out'] for r in initial),'final_timeouts':sum(r['timed_out'] for r in results)})

"""Sequential paired validation; do not tune from this output."""
import os,json,sys
from benchmark import ROOT,dump,measure
methods=['baseline','candidate','kfl1ou'];results=[]
run_root=ROOT/sys.argv[1] if len(sys.argv)>1 else ROOT
assert not (run_root/'results.json').exists(),'Refuse to overwrite results'
for index,directory in enumerate(sorted((run_root/'data').iterdir())):
 for method in methods[index%3:]+methods[:index%3]:
  output=directory/(method+'.json')
  if method=='kfl1ou':command=['Rscript',str(ROOT/'kfl1ou.R'),str(directory),str(output),'AIC']
  else:command=['env','PYTHONPATH='+str(ROOT/method),'python',str(ROOT/'native.py'),str(directory),str(output),method]
  measured=measure(command,str(directory/(method+'.log')),1800)
  row={'dataset':directory.name,'method':method,**measured}
  if output.exists():row['result']=json.loads(output.read_text())
  results.append(row);dump(run_root/'results.json',results)
  print('FINISH',len(results),'runs',directory.name,method,'exit',measured['exit_code'],round(measured['wall_seconds'],3),flush=True)
dump(run_root/'complete.json',{'runs':len(results),'success':sum(r['exit_code']==0 for r in results)})

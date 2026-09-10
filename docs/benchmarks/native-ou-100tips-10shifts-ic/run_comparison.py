"""Seven configured methods on unchanged paired inputs; run inside verified image."""
import json
from benchmark import ROOT,dump,measure

methods=['nwkit-bootstrap','nwkit-pBIC','nwkit-AIC','nwkit-BIC','kfl1ou-pBIC','kfl1ou-AIC','kfl1ou-BIC']
results=[]
assert not (ROOT/'results.json').exists(), 'Refuse to overwrite completed measurements.'
for index,directory in enumerate(sorted((ROOT/'data').iterdir())):
    order=methods[index:]+methods[:index]
    for method in order:
        package,criterion=method.split('-')
        output=directory/(method+'.json')
        command=['Rscript' if package=='kfl1ou' else 'python',str(ROOT/('kfl1ou.R' if package=='kfl1ou' else 'native.py')),str(directory),str(output),criterion]
        print('START',directory.name,method,flush=True)
        measured=measure(command,str(directory/(method+'.log')),1800)
        row={'dataset':directory.name,'method':method,**measured,'resource_limit_seconds':1800}
        if output.exists():row['result']=json.loads(output.read_text())
        results.append(row);dump(ROOT/'results.json',results)
        print('FINISH',directory.name,method,json.dumps({k:v for k,v in measured.items() if k!='command'}),flush=True)
dump(ROOT/'complete.json',{'runs':len(results)})

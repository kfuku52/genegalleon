import os,json,subprocess
from pathlib import Path
root=Path('/bench');checks=[]
for directory in sorted(Path('/old/data').iterdir()):
 output=root/'validation'/('baseline-'+directory.name+'.json')
 subprocess.run(['python',str(root/'native.py'),str(directory),str(output),'baseline'],env={**os.environ,'PYTHONPATH':str(root/'baseline')},check=True)
 actual=json.loads(output.read_text());expected=json.loads((directory/'nwkit-AIC.json').read_text())
 assert actual['selected']['shift_clades']==expected['selected']['shift_clades']
 assert abs(actual['selected']['log_likelihood']-expected['selected']['log_likelihood'])<1e-10
 assert max(abs(v-expected['selected']['predicted'][n]) for n,v in actual['selected']['predicted'].items())<1e-10
 checks.append({'dataset':directory.name,'status':'passed'})
(root/'validation'/'baseline-replay.json').write_text(json.dumps(checks,indent=2)+'\n')
print('Baseline: six original layouts, likelihoods and tip means reproduced.')

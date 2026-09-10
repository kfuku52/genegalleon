import csv,json,math
from pathlib import Path
import numpy as np
from nwkit.util import read_tree
from nwkit.shift_native_model import ShiftData,ShiftLayout
from nwkit.shift_native_fit import fit_native_layout,NativeFitOptions
from nwkit.shift_native_ic import native_information_criterion
root=Path('/bench'); old=Path('/old'); records=[]
for row in json.loads((root/'crossfit-r.json').read_text()):
 directory=old/'data'/row['dataset']
 tree=read_tree(str(directory/'tree.nwk'),1,True,quiet=True)
 tab=list(csv.DictReader((directory/'traits.tsv').open(),delimiter='\t')); lookup={r['taxon']:float(r['x']) for r in tab}
 data=ShiftData.build(tree,np.array([[lookup[n]] for n in tree.leaf_names()]),['x'])
 nodes=list(tree.traverse('levelorder')); mapping={tuple(sorted(n.leaf_names())):i for i,n in enumerate(nodes)}
 saved=json.loads((directory/(row['origin']+'.json')).read_text()); clades=saved['shift_clades'] if row['origin']=='truth' else saved['selected']['shift_clades']
 layout=ShiftLayout.build(data.tree,[mapping[tuple(sorted([c] if isinstance(c,str) else c))] for c in clades])
 pool=json.loads((directory/'nwkit-AIC.json').read_text())['search_metadata']['screening']['pool']
 result={**row,'pool_coverage':sum(b in pool for b in layout.shifts),'shifts':len(layout.shifts)}
 for mode,kwargs in [('native',{}),('r_alpha',{'alpha_height':row['alpha']})]:
  fit=fit_native_layout(data,layout,**kwargs); f=fit['fits'][0]
  result[mode]={'alpha':str(f.alpha_height),'log_likelihood':fit['log_likelihood'],'score':native_information_criterion(data,fit,'AIC')['score'],'predicted':dict(zip(data.tree.leaf_names,(data.centers[0]+data.scales[0]*f.predicted).tolist()))}
 assert abs(result['r_alpha']['log_likelihood']-row['log_likelihood'])<2e-5,result
 records.append(result)
 print(row['dataset'],row['origin'],'coverage',result['pool_coverage'],result['shifts'],'alpha',row['alpha'],result['native']['alpha'],'LL',round(row['log_likelihood'],3),round(result['native']['log_likelihood'],3),flush=True)
(root/'crossfit.json').write_text(json.dumps(records,indent=2)+'\n')

import csv,json,time,sys
from pathlib import Path
import numpy as np
from nwkit.util import read_tree
from nwkit.shift_native_model import ShiftData
from nwkit.shift_native_heuristic import NativeSearchOptions
from variants import search
root=Path('/bench'); old=Path('/old'); records=[]
for variant in sys.argv[1:] or ['scale','refine','adaptive','scale-refine','scale-adaptive']:
 for directory in sorted((old/'data').iterdir()):
  started=time.perf_counter()
  tree=read_tree(str(directory/'tree.nwk'),1,True,quiet=True)
  rows=list(csv.DictReader((directory/'traits.tsv').open(),delimiter='\t')); lookup={r['taxon']:float(r['x']) for r in rows}
  data=ShiftData.build(tree,np.array([[lookup[n]] for n in tree.leaf_names()]),['x'])
  options=NativeSearchOptions(max_shifts=10)
  implementation=__import__('path_variant').search if variant=='path' else search
  if variant=='production':
   from nwkit.shift_native_path import sparse_native_search
   result=sparse_native_search(data,options=options)
  else: result=implementation(data,options,variant)
  selected=result.best_information; fit=selected['fits'][0]
  nodes=list(tree.traverse('levelorder'));clades=[sorted(nodes[i].leaf_names()) for i in selected['layout'].shifts]
  truth=json.loads((directory/'truth.json').read_text());tp=len(set(map(tuple,clades))&set(map(tuple,truth['shift_clades'])))
  pred=data.centers[0]+data.scales[0]*fit.predicted
  rmse=float(np.sqrt(np.mean((pred-np.array([truth['tip_mean'][n] for n in data.tree.leaf_names]))**2)))
  record={'variant':variant,'dataset':directory.name,'seconds':time.perf_counter()-started,'f1':2*tp/(len(clades)+10),'rmse':rmse,'aic':selected['information_criterion']['score'],'alpha':str(fit.alpha_height),'clades':clades,'metadata':result.metadata}
  records.append(record);(root/('pilot-'+variant+'.json')).write_text(json.dumps([r for r in records if r['variant']==variant],indent=2)+'\n')
  print(variant,directory.name,round(record['seconds'],2),round(record['f1'],3),round(rmse,3),round(record['aic'],2),flush=True)

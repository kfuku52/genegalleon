"""Fresh-process baseline/candidate adapter, selected via PYTHONPATH."""
import csv,json,sys,time
from pathlib import Path
import numpy as np
from nwkit.util import read_tree
from nwkit.shift_native_model import ShiftData
from nwkit.shift_native_heuristic import NativeSearchOptions,heuristic_native_search
started=time.perf_counter()
directory,output=map(Path,sys.argv[1:3]);method=sys.argv[3]
rows=list(csv.DictReader((directory/'traits.tsv').open(),delimiter='\t'));lookup={r['taxon']:float(r['x']) for r in rows}
tree=read_tree(str(directory/'tree.nwk'),1,True,quiet=True)
data=ShiftData.build(tree,np.array([[lookup[n]] for n in tree.leaf_names()]),['x'])
options=NativeSearchOptions(max_shifts=10,candidate_pool=24,refit_budget=48,screening_budget=2000,beam_width=2,lasso_iterations=150)
if method=='candidate':
 from nwkit.shift_native_path import sparse_native_search
 search=sparse_native_search
else:search=heuristic_native_search
result=search(data,options=options,criterion='AIC');selected=result.best_information;fit=selected['fits'][0]
nodes=list(tree.traverse('levelorder'));predicted=data.centers[0]+data.scales[0]*fit.predicted
record={'status':'complete','selection':'AIC','method':method,'score':selected['information_criterion']['score'],'selected':{'shift_clades':[sorted(nodes[i].leaf_names()) for i in selected['layout'].shifts],'log_likelihood':selected['log_likelihood'],'predicted':dict(zip(data.tree.leaf_names,predicted.tolist()))},'alpha':str(fit.alpha_height),'metadata':result.metadata,'candidate_scores':result.records,'elapsed_seconds':time.perf_counter()-started}
output.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')

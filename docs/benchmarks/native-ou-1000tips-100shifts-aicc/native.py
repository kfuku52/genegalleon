"""NWKIT beam/path adapter with explicit information criterion."""
import csv,json,sys,time
from pathlib import Path
import numpy as np
from nwkit.util import read_tree
from nwkit.shift_native_model import ShiftData
from nwkit.shift_native_heuristic import NativeSearchOptions,heuristic_native_search
from nwkit.shift_native_path import sparse_native_search
started=time.perf_counter()
directory,output=map(Path,sys.argv[1:3]);strategy,criterion=sys.argv[3:5]
rows=list(csv.DictReader((directory/'traits.tsv').open(),delimiter='\t'));lookup={r['taxon']:float(r['x']) for r in rows}
tree=read_tree(str(directory/'tree.nwk'),1,True,quiet=True)
data=ShiftData.build(tree,np.array([[lookup[n]] for n in tree.leaf_names()]),['x'])
options=NativeSearchOptions(max_shifts=100,candidate_pool=128,refit_budget=220,screening_budget=20000,beam_width=1,lasso_iterations=150)
search=sparse_native_search if strategy=='path' else heuristic_native_search
result=search(data,options=options,criterion=criterion);selected=result.best_information;fit=selected['fits'][0]
nodes=list(tree.traverse('levelorder'));predicted=data.centers[0]+data.scales[0]*fit.predicted
record={'status':'complete','selection':criterion,'strategy':strategy,'score':selected['information_criterion']['score'],'information_criterion':selected['information_criterion'],'selected':{'shift_clades':[sorted(nodes[i].leaf_names()) for i in selected['layout'].shifts],'log_likelihood':selected['log_likelihood'],'predicted':dict(zip(data.tree.leaf_names,predicted.tolist()))},'alpha':str(fit.alpha_height),'metadata':result.metadata,'candidate_scores':result.records,'elapsed_seconds':time.perf_counter()-started}
output.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')

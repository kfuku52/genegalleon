"""Delivered default/explicit routing against frozen independent outputs."""
import csv,json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from nwkit.util import read_tree
from nwkit.shift_native_model import ShiftData
from nwkit.shift_native_selection import NativeSearchRunner
root=Path('/bench');checks=[]
for dataset in ['effect2-seed28101','effect6-seed28101']:
 d=root/'data'/dataset;tree=read_tree(str(d/'tree.nwk'),1,True,quiet=True);rows=list(csv.DictReader((d/'traits.tsv').open(),delimiter='\t'));lookup={r['taxon']:float(r['x']) for r in rows};data=ShiftData.build(tree,np.array([[lookup[n]] for n in tree.leaf_names()]),['x']);nodes=list(tree.traverse('levelorder'))
 for strategy,method in [('auto','baseline'),('native-path','candidate')]:
  args=SimpleNamespace(search_strategy=strategy,criterion='AIC',max_shifts=10,convergence=False,candidate_pool=24,refit_budget=48,screening_budget=2000,beam_width=2,lasso_iterations=150,search_memory_mb=512,exhaustive_max_configurations=5000)
  fit=NativeSearchRunner(data,args,{})(data).best_information;expected=json.loads((d/(method+'.json')).read_text());pred=data.centers[0]+data.scales[0]*fit['fits'][0].predicted
  assert [sorted(nodes[i].leaf_names()) for i in fit['layout'].shifts]==expected['selected']['shift_clades']
  error=abs(fit['log_likelihood']-expected['selected']['log_likelihood']);prediction_error=max(abs(v-expected['selected']['predicted'][n]) for n,v in zip(data.tree.leaf_names,pred))
  assert max(error,prediction_error)<1e-10
  checks.append({'dataset':dataset,'strategy':strategy,'reference':method,'likelihood_error':error,'prediction_error':prediction_error})
(root/'validation'/'delivered-routing-replay.json').write_text(json.dumps({'status':'passed','checks':checks},indent=2)+'\n');print('Delivered routing: four layouts, likelihoods and means reproduced.')

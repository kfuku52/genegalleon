import csv
import json
import sys
import time
from pathlib import Path
import numpy as np
from nwkit.util import read_tree
from nwkit.shift_native_model import ShiftData
from nwkit.shift_native_fit import NativeFitOptions
from nwkit.shift_native_heuristic import NativeSearchOptions, heuristic_native_search
from nwkit.shift_native_bootstrap import calibrate_native_search
from nwkit.shift_native_provenance import native_implementation

directory,output=map(Path,sys.argv[1:])
started=time.perf_counter()
with (directory/'traits.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
tree=read_tree(str(directory/'tree.nwk'),1,True,quiet=True)
lookup={r['taxon']:float(r['x']) for r in rows}
names=list(tree.leaf_names())
data=ShiftData.build(tree,np.array([[lookup[n]] for n in names]),['x'])
options=NativeSearchOptions(max_shifts=10,convergence=False,candidate_pool=24,refit_budget=48,screening_budget=2000,beam_width=2)
record={'status':'running','implementation':native_implementation(),'selection':'full-search plug-in bootstrap B=19 level=.05','search_calls':0}
def save():output.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
def summarize(result):
    nodes=list(tree.traverse('levelorder'))
    fit=result['fits'][0]
    predicted=data.centers[0]+data.scales[0]*fit.predicted
    return {'shift_clades':[sorted(nodes[i].leaf_names()) for i in result['layout'].shifts],'log_likelihood':result['log_likelihood'],'predicted':dict(zip(data.tree.leaf_names,predicted.tolist()))}
def search(sample):
    result=heuristic_native_search(sample,options=options,fit_arguments={'options':NativeFitOptions(estimate_measurement_error=False)})
    record['search_calls']+=1
    record['elapsed_seconds']=time.perf_counter()-started
    save()
    return result
try:
    searched=search(data)
    record['search_seconds']=time.perf_counter()-started
    record['search_metadata']=searched.metadata
    record['uncalibrated_best']=summarize(searched.families()[-1][1])
    record['largest_fitted_shift_count']=max(len(r['shift_branch_ids']) for r in searched.records)
    save()
    selected,calibration=calibrate_native_search(data,searched,search,replicates=19,seed=json.loads((directory/'truth.json').read_text())['seed']+100000,level=.05)
    record.update(status='complete',selected=summarize(selected),calibration=calibration)
except Exception as exc:
    record.update(status='failed',error=f'{type(exc).__name__}: {exc}')
record['elapsed_seconds']=time.perf_counter()-started
save()
if record['status']=='failed':sys.exit(1)

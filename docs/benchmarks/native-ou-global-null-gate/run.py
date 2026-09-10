"""Prespecified paired pilot using independently simulated frozen AIC fixtures."""
import csv
import hashlib
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import time
import numpy as np
from nwkit.util import read_tree
from nwkit.shift_native_model import ShiftData
from nwkit.shift_native_path import sparse_native_search
from nwkit.shift_native_heuristic import NativeSearchOptions
from nwkit.shift_native_bootstrap import gate_native_aic
ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'native-ou-aic-improvement'

def run(job):
    directory, mode, replicates = job
    name = ('extension-' if directory.parent.parent.name == 'extension' else '') + directory.name + '-' + mode
    if replicates != 19:
        name += '-B' + str(replicates)
    output = ROOT / (name + '.json')
    if output.exists():
        raise ValueError('Refuse overwrite: ' + str(output))
    rows = list(csv.DictReader((directory/'traits.tsv').open(), delimiter='\t'))
    lookup = {r['taxon']:float(r['x']) for r in rows}
    tree = read_tree(str(directory/'tree.nwk'), 1, True, quiet=True)
    data = ShiftData.build(tree, np.array([[lookup[n]] for n in tree.leaf_names()]), ['x'])
    truth = json.loads((directory/'truth.json').read_text())
    args = {'alpha_height':3.0} if mode == 'fixed-alpha' else {}
    options = NativeSearchOptions(max_shifts=10, candidate_pool=24, refit_budget=48, screening_budget=2000, beam_width=2, lasso_iterations=150)
    def search(sample):
        return sparse_native_search(sample, options=options, criterion='AIC', fit_arguments=args)
    start = time.perf_counter()
    searched = search(data)
    search_seconds = time.perf_counter()-start
    def metrics(selected):
        nodes = list(tree.traverse('levelorder'))
        chosen = {tuple(sorted(nodes[i].leaf_names())) for i in selected['layout'].shifts}
        expected = {tuple(sorted(c)) for c in truth['shift_clades']}
        tp = len(chosen & expected)
        predicted = data.centers[0]+data.scales[0]*selected['fits'][0].predicted
        return {'k':len(chosen), 'tp':tp, 'fp':len(chosen-expected), 'precision':tp/len(chosen) if chosen else None, 'recall':tp/len(expected) if expected else None, 'rmse':float(np.sqrt(np.mean((predicted-np.array([truth['tip_mean'][n] for n in data.tree.leaf_names]))**2))), 'alpha':str(selected['fits'][0].alpha_height)}
    record = {'dataset':directory.name,'mode':mode,'true_k':len(truth['shift_clades']),'effect':truth['effect'],'search_seconds':search_seconds,'ungated':metrics(searched.best_information),'null_alpha':str(searched.best_by_complexity[0]['fits'][0].alpha_height),'inputs':{p:hashlib.sha256((directory/p).read_bytes()).hexdigest() for p in ['tree.nwk','traits.tsv','truth.json']}}
    if mode == 'estimated-alpha':
        start = time.perf_counter()
        selected, calibration = gate_native_aic(data, searched, search, replicates=replicates, seed=truth['seed']+71000, level=.05)
        record.update(gate_seconds=time.perf_counter()-start, gated=metrics(selected), calibration=calibration)
    output.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    return name

if __name__ == '__main__':
    null = sorted((SOURCE/'data').glob('effect0-*'))
    nonnull = [SOURCE/'data'/f'effect{e}-seed{s}' for e in [2,6] for s in range(28101,28104)]
    nonnull += sorted((SOURCE/'extension/data').iterdir())
    jobs = [(d,'fixed-alpha',19) for d in null]+[(d,'estimated-alpha',19) for d in null+nonnull]
    if '--default-check' in sys.argv:
        jobs = [(SOURCE/'data'/f'effect{e}-seed28101','estimated-alpha',199) for e in [0,2,6]]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for name in pool.map(run,jobs):
            print('FINISH',name,flush=True)

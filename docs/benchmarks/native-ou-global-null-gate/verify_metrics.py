"""Audit outcomes against previous raw candidate predictions and newly fitted nulls."""
import csv
import json
from pathlib import Path
import numpy as np
from nwkit.util import read_tree
from nwkit.shift_native_model import ShiftData, ShiftLayout
from nwkit.shift_native_fit import fit_native_layout
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parent/'native-ou-aic-improvement'
checked=[]
for path in sorted(ROOT.glob('*-estimated-alpha*.json')):
    r=json.loads(path.read_text())
    directory=SOURCE/('extension/data' if path.name.startswith('extension-') else 'data')/r['dataset']
    truth=json.loads((directory/'truth.json').read_text())
    previous=json.loads((directory/'candidate.json').read_text())['selected']
    expected={tuple(sorted(c)) for c in truth['shift_clades']}
    def metrics(clades,predicted):
        chosen={tuple(sorted(c)) for c in clades};tp=len(chosen & expected)
        return {'k':len(chosen),'tp':tp,'fp':len(chosen-expected),'precision':tp/len(chosen) if chosen else None,'recall':tp/len(expected) if expected else None,'rmse':float(np.sqrt(np.mean([(predicted[n]-v)**2 for n,v in truth['tip_mean'].items()])))}
    raw={'ungated':previous}
    if r['calibration']['rejected']:
        raw['gated']=previous
    else:
        lookup={x['taxon']:float(x['x']) for x in csv.DictReader((directory/'traits.tsv').open(),delimiter='\t')}
        tree=read_tree(str(directory/'tree.nwk'),1,True,quiet=True)
        data=ShiftData.build(tree,np.array([[lookup[n]] for n in tree.leaf_names()]),['x'])
        fitted=fit_native_layout(data,ShiftLayout.build(data.tree))
        predicted=data.centers[0]+data.scales[0]*fitted['fits'][0].predicted
        raw['gated']={'shift_clades':[],'predicted':dict(zip(data.tree.leaf_names,predicted.tolist()))}
    for field in ['ungated','gated']:
        actual=metrics(raw[field]['shift_clades'],raw[field]['predicted'])
        for key,value in actual.items():
            if value is None:assert r[field][key] is None
            else:assert np.isclose(value,r[field][key],rtol=1e-9,atol=1e-10),(path.name,field,key,value,r[field][key])
    checked.append(path.name)
(ROOT/'independent-metric-qa.json').write_text(json.dumps({'checked':checked,'result':'all clade metrics and prediction errors agree with previous raw AIC outputs and independent null refits','tolerance':'rtol=1e-9, atol=1e-10'},indent=2)+'\n')
print('Verified',len(checked),'paired results')

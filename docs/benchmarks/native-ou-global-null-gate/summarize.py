"""Recompute paired pilot outcomes; no inference from individual selected-branch p-values."""
import json
import math
from pathlib import Path
from statistics import mean
ROOT = Path(__file__).resolve().parent

def interval(k,n):
    # Two-sided 95% Wilson interval for independent outer replicates.
    z=1.959963984540054
    p=k/n; den=1+z*z/n
    center=(p+z*z/(2*n))/den
    half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [center-half,center+half]

def summary(rows,field):
    results=[r[field] for r in rows]
    selected=sum(r['k']>0 for r in results)
    tp=sum(r['tp'] for r in results); fp=sum(r['fp'] for r in results)
    true=sum(r['true_k'] for r in rows)
    return {'n':len(rows),'any_shift_count':selected,'any_shift_rate':selected/len(rows),'any_shift_wilson95':interval(selected,len(rows)), 'mean_k':mean(r['k'] for r in results),'mean_fp':mean(r['fp'] for r in results),'micro_precision':tp/(tp+fp) if tp+fp else None,'micro_recall':tp/true if true else None,'mean_rmse':mean(r['rmse'] for r in results)}

if __name__=='__main__':
    rows=[json.loads(p.read_text()) for p in sorted(ROOT.glob('*-estimated-alpha.json'))]
    assert len(rows)==48, len(rows)
    groups={}
    for k,e in sorted({(r['true_k'],r['effect']) for r in rows}):
        selected=[r for r in rows if r['true_k']==k and r['effect']==e]
        groups[f'K{k}-effect{e}']={f:summary(selected,f) for f in ['ungated','gated']}
    fixed=[json.loads(p.read_text()) for p in sorted(ROOT.glob('*-fixed-alpha.json'))]
    assert len(fixed)==30
    result={'pilot_B':19,'groups':groups,'fixed_alpha_null':summary(fixed,'ungated'),'mean_search_seconds':mean(r['search_seconds'] for r in rows),'mean_gate_seconds':mean(r['gate_seconds'] for r in rows),'paired_outcomes':len(rows)}
    checks=[]
    for r in rows:
        c=r['calibration']; s=c['bootstrap_statistics']; t=c['statistic']
        exceed=sum(v>=t-1e-10*max(1,t) for v in s)
        assert len(s)==19 and exceed==c['exceedances'] and c['p_value']==(1+exceed)/20
        assert c['rejected']==(c['p_value']<=.05)
        assert r['gated']==r['ungated'] if c['rejected'] else r['gated']['k']==0
        checks.append(r['dataset'])
    (ROOT/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    (ROOT/'metric-qa.json').write_text(json.dumps({'validated_pairs':len(checks),'p_values_and_gate_decisions':'passed'},indent=2)+'\n')
    print(json.dumps(result,indent=2))

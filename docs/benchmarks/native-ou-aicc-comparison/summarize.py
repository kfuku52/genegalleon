"""Combine frozen seven-method accuracy with AICc and current path measurements."""
import csv,json,math
from pathlib import Path
from statistics import mean,median
ROOT=Path(__file__).resolve().parent
PREVIOUS=ROOT.parent/'native-ou-100tips-10shifts-ic'
PATH_STUDY=ROOT.parent/'native-ou-aic-improvement'
METHODS=['nwkit-bootstrap','nwkit-pBIC','nwkit-BIC','nwkit-AIC','nwkit-beam-AICc','nwkit-path-AIC','nwkit-path-AICc','kfl1ou-pBIC','kfl1ou-BIC','kfl1ou-AIC','kfl1ou-AICc']
LABELS={'nwkit-bootstrap':'NWKIT beam bootstrap (B=19)','nwkit-pBIC':'NWKIT beam pBIC','nwkit-BIC':'NWKIT beam BIC','nwkit-AIC':'NWKIT beam AIC','nwkit-beam-AICc':'NWKIT beam AICc','nwkit-path-AIC':'NWKIT path AIC','nwkit-path-AICc':'NWKIT path AICc','kfl1ou-pBIC':'kfl1ou pBIC','kfl1ou-BIC':'kfl1ou BIC','kfl1ou-AIC':'kfl1ou AIC','kfl1ou-AICc':'kfl1ou AICc','nwkit-path-AIC-gate':'NWKIT path AIC + gate (B=19)'}
def metrics(model,truth):
    chosen={tuple(sorted([c] if isinstance(c,str) else c)) for c in model['shift_clades']}
    actual={tuple(sorted(c)) for c in truth['shift_clades']}
    assert set(model['predicted'])==set(truth['tip_mean'])
    tp=len(chosen & actual);fp=len(chosen-actual);fn=len(actual-chosen)
    return {'k':len(chosen),'tp':tp,'fp':fp,'fn':fn,'precision':tp/len(chosen) if chosen else 0.,'recall':tp/len(actual) if actual else None,'f1':2*tp/(2*tp+fp+fn) if actual else None,'rmse':math.sqrt(mean((model['predicted'][n]-v)**2 for n,v in truth['tip_mean'].items()))}
def wilson(k,n):
    z=1.959963984540054;p=k/n;den=1+z*z/n;c=(p+z*z/(2*n))/den;h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [c-h,c+h]
raw=json.loads((ROOT/'results.json').read_text())
assert len(raw)==120
assert len({(r['dataset'],r['method']) for r in raw})==120
old=json.loads((PREVIOUS/'results.json').read_text())
old_lookup={(r['dataset'],r['method']):r for r in old}
rows=[];replay=[]
for run in old+raw:
    is_null=run['dataset'].startswith('effect0-')
    directory=(PATH_STUDY if is_null else PREVIOUS)/'data'/run['dataset']
    truth=json.loads((directory/'truth.json').read_text())
    result=run.get('result',{})
    assert run['exit_code']==0 and not run['timed_out'] and result.get('status')=='complete',run['dataset']
    values=metrics(result['selected'],truth)
    row={'dataset':run['dataset'],'method':run['method'],'effect':truth['effect'],**values,'wall_seconds':run['wall_seconds'],'peak_rss_mib':run['peak_rss_mib']}
    if run['method']=='kfl1ou-AIC-replay':
        baseline=old_lookup[run['dataset'],'kfl1ou-AIC']['result']['selected']
        expected=metrics(baseline,truth)
        assert all((expected[k] is None and v is None) or (v is not None and math.isclose(expected[k],v,rel_tol=1e-9,abs_tol=1e-10)) for k,v in values.items())
        assert abs(baseline['log_likelihood']-result['selected']['log_likelihood'])<1e-9
        assert max(abs(baseline['predicted'][n]-v) for n,v in result['selected']['predicted'].items())<1e-9
        replay.append(run['dataset'])
    else:rows.append(row)
# Three previously measured AIC procedures on the same null fixtures.
for directory in sorted((PATH_STUDY/'data').glob('effect0-*')):
    truth=json.loads((directory/'truth.json').read_text())
    for method,filename in [('nwkit-AIC','baseline.json'),('nwkit-path-AIC','candidate.json'),('kfl1ou-AIC','kfl1ou.json')]:
        result=json.loads((directory/filename).read_text())
        rows.append({'dataset':directory.name,'method':method,'effect':0,**metrics(result['selected'],truth),'wall_seconds':None,'peak_rss_mib':None})
    gate=json.loads((ROOT.parent/'native-ou-global-null-gate'/(directory.name+'-estimated-alpha.json')).read_text())['gated']
    rows.append({'dataset':directory.name,'method':'nwkit-path-AIC-gate','effect':0,'k':gate['k'],'tp':0,'fp':gate['fp'],'fn':0,'precision':0.,'recall':None,'f1':None,'rmse':gate['rmse'],'wall_seconds':None,'peak_rss_mib':None})
summary=[]
for effect in [2,6]:
    for method in METHODS:
        group=[r for r in rows if r['effect']==effect and r['method']==method]
        assert len(group)==3,(effect,method,len(group))
        summary.append({'effect':effect,'method':method,'n':len(group),**{field:mean(r[field] for r in group) for field in ['k','fp','precision','recall','f1','rmse']},'wall_seconds_median':median(r['wall_seconds'] for r in group)})
null=[]
for method in ['nwkit-AIC','nwkit-beam-AICc','nwkit-path-AIC','nwkit-path-AICc','kfl1ou-AIC','kfl1ou-AICc','nwkit-path-AIC-gate']:
    group=[r for r in rows if r['effect']==0 and r['method']==method]
    assert len(group)==30
    k=sum(r['k']>0 for r in group)
    null.append({'method':method,'n':30,'any_shift':k,'fpr':k/30,'wilson95':wilson(k,30),'mean_k':mean(r['k'] for r in group),'rmse':mean(r['rmse'] for r in group)})
with (ROOT/'metrics.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
(ROOT/'summary.json').write_text(json.dumps({'nonnull':summary,'null':null},indent=2)+'\n')
(ROOT/'metric-qa.json').write_text(json.dumps({'fresh_runs':len(raw),'all_completed':True,'nonnull_comparison_outcomes':66,'null_comparison_outcomes':210,'kfl1ou_AIC_replay_verified':replay},indent=2)+'\n')
lookup={(r['effect'],r['method']):r for r in summary}
lines=['# AICc added to the OU shift comparison','','Same 100-tip/10-shift nonnull fixtures as the original seven-method comparison; three seeds per effect. Precision, recall, F1 and RMSE are averages across datasets. Higher F1 and lower RMSE are preferable. Native beam and path searches are distinct methods.','', '| Method | Weak F1 | Weak RMSE | Strong F1 | Strong RMSE |','|---|---:|---:|---:|---:|']
for method in METHODS:
    a=lookup[2,method];b=lookup[6,method]
    lines.append(f'| {LABELS[method]} | {a["f1"]:.3f} | {a["rmse"]:.3f} | {b["f1"]:.3f} | {b["rmse"]:.3f} |')
lines+=['','## No-shift cases','','Same thirty no-shift fixtures for every row. The global-null gate row uses 19 draws from the preceding study. No null data were used to tune the AICc formula or search.','', '| Method | Any false shift | Mean false branches | Mean RMSE |','|---|---:|---:|---:|']
for r in null:lines.append(f'| {LABELS[r["method"]]} | {r["any_shift"]}/30 | {r["mean_k"]:.2f} | {r["rmse"]:.3f} |')
lines += ['', 'AICc reduces the average number of false branches and null prediction error, but all three AICc procedures still select at least one shift in 30/30 null datasets (95% Wilson interval for the rate: 88.65–100%). Thus it does not solve the false-positive problem in this setting. The gate reduces null detections, but its substantial nonnull power cost was measured on different nonnull fixtures; see the [gate study](../native-ou-global-null-gate/RESULTS.md).']
lines+=['','## Full nonnull recovery table','', '| Effect | Method | Precision | Recall | F1 | Mean false branches | Mean K | RMSE |','|---|---|---:|---:|---:|---:|---:|---:|']
for r in summary:lines.append(f'| {r["effect"]} | {LABELS[r["method"]]} | '+ ' | '.join(f'{r[f]:.3f}' for f in ['precision','recall','f1','fp','k','rmse'])+' |')
lines+=['','All 120 newly measured runs completed. The six repeated kfl1ou AIC fits reproduce the original predictions, likelihoods and recovery metrics. The seven original procedures retain their frozen accuracy results; the AICc/path additions use the same inputs and runtime image. Runtime measurements come from separate batches, some under concurrent validation load; no speed ranking is inferred here. Individual wall/CPU/RSS records remain in the raw runs.','','The nonnull comparison has only three seeds per condition, the cap equals the true shift count, and these fixtures were used during native-path development. Null rates use thirty seeds; see the Wilson intervals in `summary.json`. These results do not establish general superiority, a calibrated significance test, or production adoption. See [protocol](README.md), [raw metrics](metrics.csv), and [score validation](validation/native-kfl-score-agreement.json).']
(ROOT/'RESULTS.md').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines[:20]))

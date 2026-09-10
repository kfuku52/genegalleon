"""Summarize completed configured-procedure measurements, retaining all failures."""
import csv,hashlib,json,math
from pathlib import Path
from statistics import mean,median
ROOT=Path(__file__).resolve().parent
raw=json.loads((ROOT/'results.json').read_text())
assert len(raw)==12,'Wait for the entire prespecified measurement batch'
assert len({(r['dataset'],r['method']) for r in raw})==12
rows=[]
for run in raw:
    directory=ROOT/'data'/run['dataset']
    truth=json.loads((directory/'truth.json').read_text())
    assert all(hashlib.sha256((directory/name).read_bytes()).hexdigest()==sha for name,sha in run['input_sha256'].items())
    row={k:v for k,v in run.items() if k not in ['input_sha256','result','command']}
    row['effect']=truth['effect'];row['seed']=truth['seed']
    row['cpu_seconds']=run['user_seconds']+run['system_seconds']
    result=run.get('result',{})
    row['status']='timeout' if run['timed_out'] else result.get('status','failed')
    if row['status']=='complete':
        assert run['exit_code']==0 and result['selection']=='AICc'
        actual={tuple(sorted(c)) for c in truth['shift_clades']}
        selected=result['selected'];chosen={tuple(sorted(c)) for c in selected['shift_clades']}
        assert len(actual)==100 and set(selected['predicted'])==set(truth['tip_mean'])
        tp=len(actual&chosen);fp=len(chosen-actual);fn=len(actual-chosen)
        row.update(k=len(chosen),tp=tp,fp=fp,precision=tp/len(chosen) if chosen else 0.,recall=tp/100,f1=2*tp/(2*tp+fp+fn),rmse=math.sqrt(mean((selected['predicted'][n]-v)**2 for n,v in truth['tip_mean'].items())))
        candidates=result['candidate_scores'];meta=result['metadata']
        row.update(refits=len(candidates),largest_fitted_shifts=max(len(c['shift_branch_ids']) for c in candidates),screening_evaluations=meta.get('screening_evaluations',meta.get('quick_evaluations')),budget_exhausted=meta['budget_exhausted'],all_paths_converged=meta.get('all_paths_converged'),alpha=result['alpha'])
        assert all(c['complete_covariance_modes'] for c in candidates)
        assert selected['log_likelihood']==next(c['log_likelihood'] for c in candidates if c['information_criterion']['score']==result['score'])
        score=result['information_criterion'];p=score['parameter_count'];n=score['sample_size']
        assert n==1000 and p==2*len(chosen)+3
        expected=-2*selected['log_likelihood']+2*p+2*p*(p+1)/(n-p-1)
        assert abs(expected-result['score'])<1e-7
    rows.append(row)
keys=list(dict.fromkeys(k for row in rows for k in row))
with (ROOT/'metrics.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
summary=[]
for effect in [2,6]:
    for method in ['path','beam']:
        group=[r for r in rows if r['effect']==effect and r['method']==method]
        assert len(group)==3
        complete=[r for r in group if r['status']=='complete']
        item={'effect':effect,'method':method,'complete':len(complete),'runs':3,'timeouts':sum(r['status']=='timeout' for r in group)}
        if complete:
            for field in ['wall_seconds','cpu_seconds','peak_rss_mib']:
                item[field+'_median_completed']=median(r[field] for r in complete)
            item['wall_range_completed']=[min(r['wall_seconds'] for r in complete),max(r['wall_seconds'] for r in complete)]
            for field in ['k','tp','fp','precision','recall','f1','rmse','refits','largest_fitted_shifts']:
                item[field+'_mean_completed']=mean(r[field] for r in complete)
        summary.append(item)
ratios=[]
for effect in [2,6]:
    a,b=[next(r for r in summary if r['effect']==effect and r['method']==m) for m in ['path','beam']]
    if a['complete']==b['complete']==3:
        paired=[next(r['wall_seconds'] for r in rows if r['effect']==effect and r['seed']==seed and r['method']=='beam')/next(r['wall_seconds'] for r in rows if r['effect']==effect and r['seed']==seed and r['method']=='path') for seed in [30101,30102,30103]]
        ratios.append({'effect':effect,'ratio_of_median_beam_to_path_wall_time':b['wall_seconds_median_completed']/a['wall_seconds_median_completed'],'paired_ratios':paired,'median_paired_ratio':median(paired)})
(ROOT/'summary.json').write_text(json.dumps({'groups':summary,'ratios':ratios},indent=2)+'\n')
lines=['# 1,000-tip / 100-shift AICc timing comparison','','One trait; balanced tree; unknown alpha and process variance refitted; max 100 shifts; three datasets per effect. Cache warmups excluded. The two procedures can select different models.','', '| Effect | Method | Complete | Wall median (s) | Wall range (s) | CPU median (s) | Peak RSS median (MiB) |','|---|---|---:|---:|---:|---:|---:|']
for r in summary:
    if r['complete']:
        low,high=r['wall_range_completed'];values=f'{r["wall_seconds_median_completed"]:.2f} | {low:.2f}–{high:.2f} | {r["cpu_seconds_median_completed"]:.2f} | {r["peak_rss_mib_median_completed"]:.1f}'
    else:values='— | — | — | —'
    lines.append(f'| {r["effect"]} | {r["method"]} AICc | {r["complete"]}/3 | {values} |')
if len(ratios)==2:
    lines += ['', 'In these workloads, path used about one fifth of the beam search wall time: ratios of beam/path medians were '+ ' and '.join(f"{r['ratio_of_median_beam_to_path_wall_time']:.2f}" for r in ratios) + ' at effect scales 2 and 6, respectively. The procedures produced different models; this is a configured-procedure timing comparison.']
lines+=['','Medians are among completed runs; interpret them with completion/timeout counts. No censored run is discarded. Ratios are reported only when both methods completed all three cases.','', '| Effect | Method | Mean selected K | Mean true shifts recovered | Precision | Recall | F1 | Mean RMSE | Mean refits | Largest fitted K (mean) |','|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
for r in summary:
    if r['complete']:lines.append(f'| {r["effect"]} | {r["method"]} AICc | '+' | '.join(f'{r[f+"_mean_completed"]:.3f}' for f in ['k','tp','precision','recall','f1','rmse','refits','largest_fitted_shifts'])+' |')
limited=sum(r['status']=='complete' and r['method']=='path' and r['all_paths_converged'] is False for r in rows)
lines += ['', f'{limited} completed path runs contain candidate-path points that did not meet convergence within the configured iteration budget. All retained final covariance-fit modes passed the search completion checks; the candidate generator remains approximate.']
lines+=['','These are configured search costs, not times to guarantee recovery of all 100 true shifts. The cap equals truth. Candidate paths may hit their 150-iteration limit; final selected coefficients are unpenalized fits. Per-run candidate/convergence/budget metadata is preserved. Shared regimes, sampling error, multiple traits, null calibration and support bootstrap are excluded.','','See [protocol](README.md), [raw metrics](metrics.csv), [summary](summary.json), and [input validation](input-qa.json). This is GeneGalleon Docker validation on a shared host, with no concurrent assistant benchmark or test suite. No SIF claim is made.']
(ROOT/'RESULTS.md').write_text('\n'.join(lines)+'\n')
(ROOT/'metric-qa.json').write_text(json.dumps({'rows':len(rows),'complete':sum(r['status']=='complete' for r in rows),'input_hashes':'passed','aicc_formula':'passed for completed fits','clade_and_prediction_counts':'passed'},indent=2)+'\n')
print(json.dumps({'summary':summary,'ratios':ratios},indent=2))

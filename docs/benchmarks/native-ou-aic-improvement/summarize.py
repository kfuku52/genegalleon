"""Independent paired validation metrics and seed-bootstrap uncertainty."""
import csv,json,statistics
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
METHODS=['baseline','candidate','kfl1ou']
LABELS={'baseline':'NWKIT original','candidate':'NWKIT updated path','kfl1ou':'kfl1ou AIC'}
raw=json.loads((ROOT/'results.json').read_text())
assert len(raw)==270 and len({(r['dataset'],r['method']) for r in raw})==270
rows=[]
for run in raw:
 truth=json.loads((ROOT/'data'/run['dataset']/'truth.json').read_text())
 result=run.get('result',{});status='timeout' if run['timed_out'] else result.get('status','failed')
 if run['exit_code']!=0 and status=='complete':status='failed'
 row={k:v for k,v in run.items() if k not in ['command','result']}
 row.update(effect=truth['effect'],seed=truth['seed'],status=status,cpu_seconds=run['user_seconds']+run['system_seconds'])
 if status=='complete':
  selected=result['selected'];pred={tuple(sorted([c] if isinstance(c,str) else c)) for c in selected['shift_clades']};actual=set(map(tuple,truth['shift_clades']))
  assert set(selected['predicted'])==set(truth['tip_mean'])
  tp=len(pred&actual);fp=len(pred-actual);fn=len(actual-pred)
  rmse=float(np.sqrt(np.mean([(selected['predicted'][n]-x)**2 for n,x in truth['tip_mean'].items()])))
  row.update(tp=tp,fp=fp,fn=fn,selected_shifts=len(pred),any_shift=int(bool(pred)),precision=tp/len(pred) if pred else 0,recall=tp/len(actual) if actual else None,f1=2*tp/(len(pred)+len(actual)) if actual else None,mean_rmse=rmse,aic=result['score'])
  assert len(pred)<=10 and tp+fn==(10 if truth['effect'] else 0)
  if actual:assert 0<=row['f1']<=1
 rows.append(row)
keys=list(dict.fromkeys(k for r in rows for k in r))
with (ROOT/'metrics.csv').open('w') as f:
 writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)
summary=[]
for effect in [0,2,6]:
 for method in METHODS:
  runs=[r for r in rows if r['effect']==effect and r['method']==method];complete=[r for r in runs if r['status']=='complete']
  r={'effect':effect,'method':method,'complete':len(complete),'runs':len(runs),'failed':sum(x['status']=='failed' for x in runs),'timeout':sum(x['status']=='timeout' for x in runs)}
  for key in ['wall_seconds','cpu_seconds','peak_rss_mib']:r[key+'_median']=statistics.median(x[key] for x in runs)
  for key in ['f1','mean_rmse','any_shift','selected_shifts','precision','recall']:
   values=[x[key] for x in complete if x[key] is not None];r[key+'_mean']=statistics.mean(values) if values else None
  summary.append(r)
(ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
# Resample matched seeds, never methods independently. The same resample indices
# are used for each metric/comparator within an effect condition.
rng=np.random.default_rng(904731);comparisons=[]
for effect in [0,2,6]:
 for comparator in ['baseline','kfl1ou']:
  groups={m:{r['seed']:r for r in rows if r['effect']==effect and r['method']==m and r['status']=='complete'} for m in ['candidate',comparator]}
  seeds=sorted(set(groups['candidate'])&set(groups[comparator]));n=len(seeds)
  if not n:continue
  indices=rng.integers(0,n,size=(10000,n))
  for key in ['mean_rmse','any_shift'] if effect==0 else ['f1','mean_rmse']:
   a=np.array([groups['candidate'][s][key] for s in seeds]);b=np.array([groups[comparator][s][key] for s in seeds]);difference=a-b
   sample=difference[indices].mean(axis=1)
   entry={'effect':effect,'comparator':comparator,'metric':key,'paired_seeds':n,'candidate_mean':float(a.mean()),'comparator_mean':float(b.mean()),'mean_difference':float(difference.mean()),'difference_ci95':np.quantile(sample,[.025,.975]).tolist()}
   if key=='mean_rmse':
    entry['ratio_of_means']=float(a.mean()/b.mean());entry['ratio_ci95']=np.quantile(a[indices].mean(axis=1)/b[indices].mean(axis=1),[.025,.975]).tolist()
   comparisons.append(entry)
(ROOT/'paired-comparisons.json').write_text(json.dumps({'bootstrap_replicates':10000,'seed':904731,'interval':'paired seed percentile bootstrap, pointwise 95%, not simultaneous','comparisons':comparisons},indent=2)+'\n')
(ROOT/'metric-qa.json').write_text(json.dumps({'status':'passed','unique_runs':len(rows),'complete':sum(r['status']=='complete' for r in rows),'null_f1':'undefined; not scored','predicted_tip_sets':'all equal to truth'},indent=2)+'\n')
lines=['# Independent native AIC comparison','','[Adoption decision and limitations](DECISION.md): the path is available explicitly; the original auto default is retained.','','Frozen development candidate evaluated on 30 new seeds at each of three effect strengths (90 datasets; 270 configured runs). All methods use a 10-shift cap. No tuning used these results.','','| Effect | Method | Complete | F1 | Mean RMSE | Any shift | Mean shift count | Median seconds | Median RSS MiB |','| --- | --- | --- | --- | --- | --- | --- | --- | --- |']
for r in summary:
 f1='—' if r['f1_mean'] is None else f"{r['f1_mean']:.3f}"
 lines.append(f"| {r['effect']} | {LABELS[r['method']]} | {r['complete']}/{r['runs']} | {f1} | {r['mean_rmse_mean']:.3f} | {r['any_shift_mean']:.3f} | {r['selected_shifts_mean']:.2f} | {r['wall_seconds_median']:.2f} | {r['peak_rss_mib_median']:.1f} |")
lines+=['','F1 measures exact descendant-clade recovery. RMSE compares fitted tip means with noise-free simulated means. Any shift is the global-null false-positive frequency only at effect 0. Null F1/recall are undefined and omitted.','','## Paired differences','','Positive F1 differences favor the candidate; negative RMSE differences favor it. Intervals resample paired seeds (10,000 replicates, pointwise 95% percentile intervals). These are not simultaneous intervals or proof of equivalence.','','| Effect | Compared with | Metric | Candidate minus comparator | 95% interval | RMSE ratio |','| --- | --- | --- | --- | --- | --- |']
for r in comparisons:
 lo,hi=r['difference_ci95'];ratio=f"{r['ratio_of_means']:.3f}" if 'ratio_of_means' in r else '—'
 lines.append(f"| {r['effect']} | {LABELS[r['comparator']]} | {r['metric']} | {r['mean_difference']:.3f} | [{lo:.3f}, {hi:.3f}] | {ratio} |")
lines+=['','![Comparison](comparison.png)','','One balanced tree, one trait, no observation errors or convergence. Nonnull shift count equals the search cap; null count is below it. Do not generalize these estimates to arbitrary trees or traits. Time includes startup and serialization; comparisons are between configured procedures with different candidate sets, not equivalent-output kernel benchmarks. All failed and timed-out runs remain in the denominators and raw records.','','[Protocol](PROTOCOL.md) · [Individual metrics](metrics.csv) · [Raw results](results.json) · [Paired intervals](paired-comparisons.json) · [SVG figure](comparison.svg)']
(ROOT/'RESULTS.md').write_text('\n'.join(lines)+'\n')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
fig,axes=plt.subplots(2,2,figsize=(12,9));colors={2:'#0072B2',6:'#D55E00'}
for ax,key,title,xlabel in [(axes[0,0],'f1','Exact shift recovery','Mean F1 (higher is better)'),(axes[0,1],'mean_rmse','Expected tip-mean recovery','Mean RMSE (lower is better)')]:
 for i,m in enumerate(METHODS):
  for effect,off in [(2,-.16),(6,.16)]:
   values=np.array([r[key] for r in rows if r['method']==m and r['effect']==effect and r['status']=='complete'])
   ci=np.quantile(values[rng.integers(0,len(values),size=(10000,len(values)))].mean(axis=1),[.025,.975]);mean=float(values.mean())
   ax.errorbar(mean,i+off,xerr=[[mean-ci[0]],[ci[1]-mean]],fmt='o' if effect==2 else '^',color=colors[effect],capsize=3,label=f'Effect {effect}' if i==0 else None)
 ax.set_yticks(range(3),[LABELS[m] for m in METHODS]);ax.set_ylim(2.6,-.6);ax.set_xlabel(xlabel);ax.set_title(title,loc='left',weight='bold');ax.set_xlim(left=0);ax.grid(axis='x',alpha=.15)
 if key=='f1':ax.set_xlim(0,1);ax.legend(loc='lower right',frameon=False)
ax=axes[1,0]
for i,m in enumerate(METHODS):
 for effect,off in [(0,-.22),(2,0),(6,.22)]:
  values=[r['wall_seconds'] for r in rows if r['method']==m and r['effect']==effect]
  ax.scatter(values,np.full(len(values),i+off),s=10,color=colors.get(effect,'#888888'),alpha=.3)
  ax.plot(statistics.median(values),i+off,'|',color='black',ms=12,mew=2)
ax.set_yticks(range(3),[LABELS[m] for m in METHODS]);ax.set_ylim(2.6,-.6);ax.set_xlim(left=0);ax.set_xlabel('Seconds; dots = runs, black marks = medians');ax.set_title('Time including selection',loc='left',weight='bold');ax.grid(axis='x',alpha=.15)
ax=axes[1,1]
for i,m in enumerate(METHODS):
 values=[r['any_shift'] for r in rows if r['method']==m and r['effect']==0 and r['status']=='complete'];count=sum(values);n=len(values);p=count/n
 # Wilson binomial interval for the unpaired descriptive null proportions.
 z=1.95996398454;den=1+z*z/n;center=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
 ax.errorbar(p,i,xerr=[[max(0.,p-(center-half))],[max(0.,(center+half)-p)]],fmt='o',color='#555555',capsize=3);ax.text(min(p+.05,.86),i-.1,f'{count}/{n}',fontsize=10)
ax.set_yticks(range(3),[LABELS[m] for m in METHODS]);ax.set_ylim(2.6,-.6);ax.set_xlim(0,1.02);ax.set_xlabel('Fraction selecting any shift (lower is better)');ax.set_title('Global-null false positives',loc='left',weight='bold');ax.grid(axis='x',alpha=.15)
fig.suptitle('Native AIC: independent 100-tip comparison',fontsize=17,weight='bold',y=.98)
fig.text(.5,.93,'30 new paired seeds per effect • 90 datasets • original NWKIT / updated path / kfl1ou',ha='center',fontsize=11)
fig.text(.06,.035,'Accuracy bars: seed-bootstrap 95% intervals. Null bars: Wilson 95% intervals.\nTiming: gray = null, blue = weak, orange = strong. One balanced tree; no broad equivalence claim.',fontsize=10)
fig.tight_layout(rect=[0,.09,1,.9]);fig.savefig(ROOT/'comparison.png',dpi=180);fig.savefig(ROOT/'comparison.svg');plt.close(fig)

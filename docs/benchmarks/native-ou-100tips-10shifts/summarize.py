"""Summarize censored paired runs; no failed selection is scored as zero."""
import csv
import json
from pathlib import Path
import statistics
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parent
audit_path=ROOT/'pbic-postrun-audit.json'
invalid_pbic=audit_path.exists() and json.loads(audit_path.read_text())['status']=='failed'

def metrics(model, truth):
    predicted={tuple(sorted([c] if isinstance(c,str) else c)) for c in model['shift_clades']}
    actual={tuple(sorted(c)) for c in truth['shift_clades']}
    tp=len(predicted & actual); fp=len(predicted-actual); fn=len(actual-predicted)
    # Empty positive prediction: precision defined as zero for this benchmark.
    precision=tp/len(predicted) if predicted else 0.
    recall=tp/len(actual)
    f1=2*tp/(2*tp+fp+fn)
    assert set(model['predicted'])==set(truth['tip_mean'])
    rmse=float(np.sqrt(np.mean([(model['predicted'][n]-value)**2 for n,value in truth['tip_mean'].items()])))
    # Tip partitions help distinguish a nearby branch error from a wholly wrong
    # trait clustering. Nested shifts assign tips to the smallest shifted clade.
    def labels(clades):
        ordered=sorted(clades,key=lambda c:(len(c),c))
        return [next((i+1 for i,c in enumerate(ordered) if tip in c),0) for tip in sorted(truth['tip_mean'])]
    from collections import Counter
    a,b=labels(actual),labels(predicted)
    choose2=lambda n:n*(n-1)/2
    pairs=sum(choose2(n) for n in Counter(zip(a,b)).values())
    aa=sum(choose2(n) for n in Counter(a).values());bb=sum(choose2(n) for n in Counter(b).values())
    expected=aa*bb/choose2(len(a));maximum=(aa+bb)/2
    ari=(pairs-expected)/(maximum-expected) if maximum!=expected else 1.
    return dict(tp=tp,fp=fp,fn=fn,selected_shifts=len(predicted),precision=precision,recall=recall,f1=f1,mean_rmse=rmse,tip_partition_ari=ari)

raw=json.loads((ROOT/'results.json').read_text())
rows=[];diagnostics=[]
for run in raw:
    truth=json.loads((ROOT/'data'/run['dataset']/'truth.json').read_text())
    result=run.get('result',{})
    status='timeout' if run['timed_out'] else result.get('status','failed')
    row={k:v for k,v in run.items() if k not in ['result','command']}
    row.setdefault('resource_limit_seconds',300)
    row.update(effect=truth['effect'],seed=truth['seed'],status=status,search_seconds=result.get('search_seconds'),search_calls=result.get('search_calls'))
    if status=='complete':row.update(metrics(result['selected'],truth))
    rows.append(row)
    if 'uncalibrated_best' in result:
        diagnostics.append(dict(dataset=run['dataset'],label='NWKIT uncalibrated maximum-likelihood candidate; not completed selection',**metrics(result['uncalibrated_best'],truth)))
keys=list(dict.fromkeys(k for r in rows for k in r))
with (ROOT/'metrics.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
(ROOT/'diagnostic-search-accuracy.json').write_text(json.dumps(diagnostics,indent=2)+'\n')
summary=[]
for effect in [2,6]:
    for method in ['kfl1ou','nwkit']:
        selected=[r for r in rows if r['method']==method and r['effect']==effect]
        complete=[r for r in selected if r['status']=='complete']
        item={'effect':effect,'method':method,'runs':len(selected),'complete':len(complete),'timeout':sum(r['status']=='timeout' for r in selected),'failed':sum(r['status']=='failed' for r in selected)}
        for key in ['wall_seconds','peak_rss_mib']:
            item[key+'_median_observed']=statistics.median(r[key] for r in selected)
        search_times=[r['search_seconds'] for r in selected if r.get('search_seconds') is not None]
        item['native_initial_search_seconds_median']=statistics.median(search_times) if search_times else None
        for key in ['precision','recall','f1','mean_rmse','selected_shifts','tip_partition_ari']:
            item[key+'_mean_completed']=statistics.mean(r[key] for r in complete) if complete else None
        summary.append(item)
(ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
by_group={(r['effect'],r['method']):r for r in summary}
ratios={effect:by_group[effect,'nwkit']['wall_seconds_median_observed']/by_group[effect,'kfl1ou']['wall_seconds_median_observed'] for effect in [2,6]}
full={method:sum(r['effect']==6 and r['method']==method and r.get('tp')==10 and r.get('fp')==0 for r in rows) for method in ['kfl1ou','nwkit']}
search_seconds=[r['search_seconds'] for r in rows if r.get('search_seconds') is not None]
search_calls=[r['search_calls'] for r in rows if r.get('search_calls') is not None]
lines=['# Measured results', '',
       'Six paired datasets: three seeds per effect scale, one trait, 100 tips and ten planted shifts. '
       'kfl1ou uses pBIC; NWKIT native includes full-search bootstrap selection (B=19).', '',
       f'Configured kfl1ou was faster: native median total time was {ratios[2]:.1f} times higher at effect scale 2 '
       f'and {ratios[6]:.1f} times higher at scale 6. Peak RSS was comparable (roughly 106–122 MiB). '
       f'At scale 6, exact recovery of all ten branches occurred in {full["kfl1ou"]}/3 kfl1ou runs and {full["nwkit"]}/3 native runs. '
       'Native recovery improved in one strong-effect dataset, but another had more false branch selections; '
       'neither method consistently recovered the planted configuration.', '',
       f'Native initial search took {min(search_seconds):.2f}–{max(search_seconds):.2f} seconds. '
       f'Including calibration, each dataset required {min(search_calls)}–{max(search_calls)} complete searches. '
       'Thus these total-time ratios compare the configured selection procedures; they are not a kernel or language speed comparison. '
       'B=19 is a coarse development calibration setting, below the workflow default of 199.', '',
       '![Paired benchmark](comparison.png)', '',
       '| Effect | Method | Complete | Median time (s) | Median RSS (MiB) | Precision | Recall | F1 | Mean RMSE |',
       '| --- | --- | --- | --- | --- | --- | --- | --- | --- |']
def number(value):return 'NA' if value is None else f'{value:.3f}'
for r in summary:
    lines.append(f'| {r["effect"]} | {r["method"]} | {r["complete"]}/{r["runs"]} | '
                 f'{r["wall_seconds_median_observed"]:.1f} | {r["peak_rss_mib_median_observed"]:.1f} | '+
                 ' | '.join(number(r[k+'_mean_completed']) for k in ['precision','recall','f1','mean_rmse'])+' |')
lines.extend(['', 'Accuracy columns are means over completed runs; precision is defined as zero for empty predictions. '
              'Resource columns summarize the final observed attempt, with any censoring identified in metrics.csv. '
              'No population confidence interval is estimated from three replicates.', '',
              'The first 300-second native timeout was rerun without changing inference settings under a 1,200-second ceiling; '
              'pending runs also used the longer ceiling. Initial and follow-up records are retained separately. '
              'The documented kfl1ou tree-order preflight error is excluded after correction and rerun.', '',
              '[Protocol and reproduction](README.md) · [Point-level metrics](metrics.csv) · [Raw final runs](results.json) · [Vector figure](comparison.svg)', '',
              'This experiment does not establish false-positive control, convergence-group accuracy, missing/error robustness, '
              'multivariate performance or production adoption. The methods have different selection rules and search/optimizer defaults.'])
if invalid_pbic:
    lines[2:2]=['**CORRECTION: the kfl1ou baseline failed the pBIC capability probe.** '
                'These are historical measurements of the uncorrected backend, not a valid comparison '
                'against corrected pBIC. Baseline rerunning is required. '
                '[Post-run audit](pbic-postrun-audit.json).', '']
(ROOT/'RESULTS.md').write_text('\n'.join(lines)+'\n')

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','svg.fonttype':'none'})
fig,axes=plt.subplots(3,2,figsize=(10.5,12.4))
colors={'kfl1ou':'#0072B2','nwkit':'#D55E00'}
markers={'kfl1ou':'o','nwkit':'s'}
fields=[('wall_seconds','Wall time (seconds; log scale)','Time incl. model selection'),('peak_rss_mib','Peak process RSS (MiB)','RAM high-water mark'),('precision','Exact branch precision','Shift precision'),('recall','Exact branch recall','Shift recall'),('f1','Exact branch F1','Shift F1'),('mean_rmse','RMSE (trait units)','True tip-mean recovery')]
for ax,(field,ylabel,title) in zip(axes.flat,fields):
    for method,offset in [('kfl1ou',-.15),('nwkit',.15)]:
        for effect,x in [(2,0),(6,1)]:
            group=[r for r in rows if r['effect']==effect and r['method']==method]
            for j,r in enumerate(group):
                if field not in r:continue
                xx=x+offset+(j-1)*.06
                censored=r['status']=='timeout'
                marker='^' if censored else markers[method]
                ax.scatter(xx,r[field],s=65,marker=marker,facecolors='none' if censored else colors[method],edgecolors=colors[method],linewidths=1.5,zorder=3)
            available=[r[field] for r in group if field in r and r['status']=='complete']
            if available:
                center=np.median(available) if field in ['wall_seconds','peak_rss_mib'] else np.mean(available)
                ax.plot([x+offset-.065,x+offset+.065],[center]*2,color=colors[method],lw=2.2,zorder=2)
    ax.set_xticks([0,1],['Effect scale 2','Effect scale 6'])
    ax.set_xlim(-.5,1.5);ax.set_ylabel(ylabel);ax.set_title(title,loc='left',fontsize=12)
    ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
    if field in ['precision','recall','f1']:
        ax.set_ylim(-.04,1.06);ax.set_yticks([0,.25,.5,.75,1])
    elif field=='wall_seconds':
        ax.set_yscale('log')
        values=[r[field] for r in rows if field in r]
        ax.set_ylim(min(values)*.8,max(values)*1.3)
    else:ax.set_ylim(0,max(r[field] for r in rows if field in r)*1.15)
    if field not in ['wall_seconds','peak_rss_mib'] and any(r['status']!='complete' for r in rows):
        for effect,x in [(2,0),(6,1)]:
            label='\n'.join(f'{m}: {sum(r["status"]=="complete" for r in rows if r["method"]==m and r["effect"]==effect)}/3' for m in ['kfl1ou','nwkit'])
            ax.text(x,-.25,label,transform=ax.get_xaxis_transform(),ha='center',fontsize=10)
handles=[Line2D([],[],marker=markers[m],linestyle='none',color=colors[m],label='kfl1ou / pBIC' if m=='kfl1ou' else 'NWKIT native / bootstrap B=19',markersize=7) for m in ['kfl1ou','nwkit']]
if any(r['status']=='timeout' for r in rows):
    handles.append(Line2D([],[],marker='^',linestyle='none',color='#555555',markerfacecolor='none',label='Timeout: censored time / RAM',markersize=7))
fig.suptitle('100 tips, 10 planted shifts',x=.09,ha='left',fontsize=18,fontweight='bold',y=.975)
fig.text(.09,.942,'Paired development benchmark | 1 trait | 3 seeds per effect | 1 thread',fontsize=12)
fig.text(.09,.919,f'{sum(r["status"]=="complete" for r in rows)}/12 runs completed; fixed-root OU; alpha and variance estimated.',fontsize=11)
if invalid_pbic:
    fig.text(.09,.84,'CORRECTION: legacy pBIC failed validation; baseline must be rerun.',fontsize=11,fontweight='bold')
fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.08,.9),ncol=1,frameon=False,fontsize=11)
fig.subplots_adjust(top=.80,bottom=.10,left=.10,right=.98,hspace=.50,wspace=.40)
fig.text(.09,.025,'Points = datasets; bars = median time/RAM, mean accuracy (completed runs only).\nInitial limit 5 min; extended to 20 min for the timed-out case and pending runs.\nNo convergence, missingness or measurement error. No global-null false-positive validation.',fontsize=10)
fig.savefig(ROOT/'comparison.png',dpi=180)
fig.savefig(ROOT/'comparison.svg')
print(json.dumps(summary,indent=2))

"""Point-level scientific comparison with explicit completion denominators."""
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
METHODS=['nwkit-bootstrap','nwkit-pBIC','nwkit-AIC','nwkit-BIC','kfl1ou-pBIC','kfl1ou-AIC','kfl1ou-BIC']
LABELS={m:m.replace('nwkit-','NWKIT ').replace('kfl1ou-','kfl1ou ') for m in METHODS}

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
assert len(raw)==42, 'Only summarize the complete planned experiment.'
assert len({(r['dataset'],r['method']) for r in raw})==42
assert json.loads((ROOT/'validation/pbic-attestation.json').read_text())['status']=='passed'
rows=[]
for run in raw:
    truth=json.loads((ROOT/'data'/run['dataset']/'truth.json').read_text())
    result=run.get('result',{})
    status='timeout' if run['timed_out'] else result.get('status','failed')
    row={k:v for k,v in run.items() if k not in ['result','command']}
    row.update(effect=truth['effect'],seed=truth['seed'],status=status,cpu_seconds=run['user_seconds']+run['system_seconds'],search_seconds=result.get('search_seconds'),search_calls=result.get('search_calls'))
    if status=='complete':
        assert run['exit_code']==0
        row.update(metrics(result['selected'],truth))
        for key in ['precision','recall','f1']:
            assert 0<=row[key]<=1
        assert row['tp']+row['fn']==10
        assert row['tp']+row['fp']==row['selected_shifts']<=10
    rows.append(row)
keys=list(dict.fromkeys(k for r in rows for k in r))
with (ROOT/'metrics.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)
summary=[]
for effect in [2,6]:
    for method in METHODS:
        group=[r for r in rows if r['effect']==effect and r['method']==method]
        complete=[r for r in group if r['status']=='complete']
        item={'effect':effect,'method':method,'runs':len(group),'complete':len(complete),'timeout':sum(r['status']=='timeout' for r in group),'failed':sum(r['status']=='failed' for r in group)}
        for key in ['wall_seconds','cpu_seconds','peak_rss_mib']:
            item[key+'_median_observed']=statistics.median(r[key] for r in group)
        for key in ['precision','recall','f1','mean_rmse','selected_shifts','tip_partition_ari']:
            item[key+'_mean_completed']=statistics.mean(r[key] for r in complete) if complete else None
        summary.append(item)
(ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
(ROOT/'metric-qa.json').write_text(json.dumps({'status':'passed','unique_runs':len(rows),'complete':sum(r['status']=='complete' for r in rows),'truth_shift_denominator':10,'input_prediction_tip_sets':'all equal','precision_recall_f1_bounds':'passed'},indent=2)+'\n')
by_group={(r['effect'],r['method']):r for r in summary}
def val(effect,method,key):return by_group[effect,method][key]
native_times=[r['wall_seconds_median_observed'] for r in summary if r['method'].startswith('nwkit-') and r['method']!='nwkit-bootstrap']
r_times=[r['wall_seconds_median_observed'] for r in summary if r['method'].startswith('kfl1ou-')]
ratios=[val(e,'kfl1ou-'+c,'wall_seconds_median_observed')/val(e,'nwkit-'+c,'wall_seconds_median_observed') for e in [2,6] for c in ['AIC','BIC','pBIC']]
f1=lambda e,m:val(e,m,'f1_mean_completed')
rmse=lambda e,m:val(e,m,'mean_rmse_mean_completed')
lines=['# Measured comparison', '',
       'Seven configured selection procedures on six paired datasets: 100 tips, ten planted shifts, one trait and three seeds per effect scale. Corrected kfl1ou pBIC passed its capability probe before timing.', '',
       f"Native IC median runtime was {min(native_times):.1f}–{max(native_times):.1f} seconds, versus {min(r_times):.1f}–{max(r_times):.1f} seconds for kfl1ou IC. Matching-criterion ratios of median time were {min(ratios):.2f}–{max(ratios):.2f} in favor of native. Bootstrap medians were {val(2,'nwkit-bootstrap','wall_seconds_median_observed'):.1f} seconds at effect scale 2 and {val(6,'nwkit-bootstrap','wall_seconds_median_observed'):.1f} seconds at scale 6. Peak process memory was similar in scale ({min(r['peak_rss_mib'] for r in rows):.0f}–{max(r['peak_rss_mib'] for r in rows):.0f} MiB across individual runs).", '',
       f"Accuracy varied by method and dataset. At weak effect scale 2, mean exact-branch F1 was {f1(2,'kfl1ou-AIC'):.3f} for kfl1ou AIC and {f1(2,'nwkit-AIC'):.3f} for native AIC. At scale 6, native BIC/AIC had mean F1 {f1(6,'nwkit-BIC'):.3f}/{f1(6,'nwkit-AIC'):.3f} and kfl1ou AIC {f1(6,'kfl1ou-AIC'):.3f}, with large differences between seeds. Expected-tip-mean RMSE for kfl1ou AIC was {rmse(2,'kfl1ou-AIC'):.3f}/{rmse(6,'kfl1ou-AIC'):.3f} across the two strengths, versus native AIC {rmse(2,'nwkit-AIC'):.3f}/{rmse(6,'nwkit-AIC'):.3f}. These small-sample results do not establish a general ranking.", '',
       'For every dataset, all four native methods evaluated the same 48 candidate layouts with identical likelihoods. The default bootstrap selected exactly the same branches and tip means as the previous implementation; only its timing was measured afresh. [Search and default-behavior audit](validation/native-search-equivalence.json).', '',
       '![Seven-method comparison](comparison.png)', '',
       '| Effect | Method | Complete | Wall median (s) | CPU median (s) | RSS median (MiB) | Precision | Recall | F1 | Mean RMSE |',
       '| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |']
def fmt(value):return 'NA' if value is None else f'{value:.3f}'
for r in summary:
    lines.append(f'| {r["effect"]} | {LABELS[r["method"]]} | {r["complete"]}/{r["runs"]} | '+
                 ' | '.join(f'{r[k+"_median_observed"]:.2f}' for k in ['wall_seconds','cpu_seconds','peak_rss_mib'])+' | '+
                 ' | '.join(fmt(r[k+'_mean_completed']) for k in ['precision','recall','f1','mean_rmse'])+' |')
lines.extend(['', 'Resource summaries are medians of observed final runs; accuracy summaries are means over completed runs. All individual measurements, failures and timeouts remain in [metrics.csv](metrics.csv). Empty predictions have precision zero. CPU time is user plus system CPU. The [summary JSON](summary.json) also contains selected shift counts and tip-partition ARI.', '',
              'Bootstrap B=19 is a coarse development calibration; the workflow default is 199. IC methods perform no calibration. Search and optimizer defaults differ between packages. These are configured-procedure comparisons, not equivalent-kernel speedups. Three replicates on one tree do not establish population confidence intervals or general reliability.', '',
              '[Protocol](README.md) · [Raw runs](results.json) · [Vector figure](comparison.svg) · [Score validation](validation/native-kfl-score-agreement.json)', '',
              'No convergence, missingness, observation error or multivariate traits are tested. No global-null false-positive or production-adoption claim is made. The 10-shift search cap equals the truth.'])
(ROOT/'RESULTS.md').write_text('\n'.join(lines)+'\n')

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','svg.fonttype':'none'})
fig,axes=plt.subplots(3,2,figsize=(15,12))
colors={2:'#0072B2',6:'#D55E00'};markers={2:'o',6:'^'}
fields=[('wall_seconds','Seconds (log scale)','Time including model selection'),('peak_rss_mib','MiB','Peak process memory'),('precision','Exact branch precision','Shift precision'),('recall','Exact branch recall','Shift recall'),('f1','Exact branch F1','Shift F1'),('mean_rmse','RMSE (trait units)','Expected tip-mean recovery')]
for ax,(field,xlabel,title) in zip(axes.flat,fields):
    for i,method in enumerate(METHODS):
        for effect,offset in [(2,-.25),(6,.25)]:
            group=sorted([r for r in rows if r['method']==method and r['effect']==effect],key=lambda r:r['seed'])
            available=[]
            for j,r in enumerate(group):
                if field not in r:continue
                censored=r['status']=='timeout'
                ax.scatter(r[field],i+offset+(j-1)*.14,s=24,marker='>' if censored else markers[effect],facecolors='none' if censored else colors[effect],edgecolors=colors[effect],linewidths=1,zorder=3)
                if r['status']=='complete':available.append(r[field])
            if available:
                center=statistics.median(available) if field in ['wall_seconds','peak_rss_mib'] else statistics.mean(available)
                ax.plot([center,center],[i+offset-.08,i+offset+.08],color='#333333',lw=1.8,zorder=4)
    ax.set_yticks(range(len(METHODS)),[LABELS[m] for m in METHODS]);ax.set_ylim(6.55,-.55)
    ax.set_xlabel(xlabel);ax.set_title(title,loc='left',fontsize=12,pad=10)
    ax.grid(axis='x',color='#dddddd',lw=.7);ax.set_axisbelow(True)
    if field=='wall_seconds':
        ax.set_xscale('log');ax.set_xlim(min(r[field] for r in rows)*.75,max(r[field] for r in rows)*1.35)
    elif field in ['precision','recall','f1']:
        ax.set_xlim(-.035,1.035);ax.set_xticks([0,.25,.5,.75,1])
    else:ax.set_xlim(0,max(r[field] for r in rows if field in r)*1.12)
handles=[Line2D([],[],marker=markers[e],linestyle='none',color=colors[e],label=f'Effect scale {e}',markersize=7) for e in [2,6]]
handles.append(Line2D([],[],marker='|',linestyle='none',color='#333333',label='Median time/RSS; mean accuracy',markersize=12))
fig.suptitle('OU shift selection: seven-method comparison',x=.13,ha='left',fontsize=19,fontweight='bold',y=.975)
fig.text(.13,.944,f'100 tips · 10 planted shifts · 1 trait · 3 paired seeds per effect · {sum(r["status"]=="complete" for r in rows)}/42 completed runs',fontsize=12)
fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.124,.931),ncol=3,frameon=False,fontsize=11)
fig.subplots_adjust(top=.855,bottom=.12,left=.13,right=.98,hspace=.4,wspace=.52)
fig.text(.13,.029,'Each point is one dataset. Bootstrap B=19; IC methods have no calibration. One thread; fixed-root OU; estimated alpha/variance.\nOne tree, 10-shift search cap; no convergence, missingness or measurement error. Exploratory comparison, not production validation.',fontsize=10)
fig.savefig(ROOT/'comparison.png',dpi=180)
fig.savefig(ROOT/'comparison.svg')
print(json.dumps(summary,indent=2))

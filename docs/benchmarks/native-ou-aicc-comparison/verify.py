"""Independent scalar score, clade metric, input and source audit."""
import csv,hashlib,json,math,tarfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
raw=json.loads((ROOT/'results.json').read_text())
assert len(raw)==120
max_score_error=0.
for row in raw:
    assert row['exit_code']==0 and row['result']['status']=='complete' and not row['timed_out']
    directory=ROOT.parent/('native-ou-aic-improvement' if row['dataset'].startswith('effect0-') else 'native-ou-100tips-10shifts-ic')/'data'/row['dataset']
    assert all(hashlib.sha256((directory/name).read_bytes()).hexdigest()==value for name,value in row['inputs'].items())
    result=row['result'];selected=result['selected'];k=len(selected['shift_clades']);p=2*k+3
    expected=-2*selected['log_likelihood']+2*p
    if result['selection']=='AICc':expected+=2*p*(p+1)/(100-p-1)
    else:assert result['selection']=='AIC'
    error=abs(result['score']-expected);max_score_error=max(max_score_error,error)
    assert error<1e-7,(row['dataset'],row['method'],error)
# Independently count sets, then use a direct sum for the RMSE.
checked=0
for row in csv.DictReader((ROOT/'metrics.csv').open()):
    if row['method']=='nwkit-path-AIC-gate':continue
    directory=ROOT.parent/('native-ou-aic-improvement' if row['effect']=='0' else 'native-ou-100tips-10shifts-ic')/'data'/row['dataset']
    truth=json.loads((directory/'truth.json').read_text())
    candidates=[r for r in raw if r['dataset']==row['dataset'] and r['method']==row['method']]
    if not candidates:
        continue  # Frozen rows were audited in their original study.
    fit=candidates[0]['result']['selected']
    actual={frozenset(x) for x in truth['shift_clades']}
    found={frozenset([x] if isinstance(x,str) else x) for x in fit['shift_clades']}
    tp=len(actual&found);fp=len(found-actual);fn=len(actual-found)
    assert (tp,fp,fn,len(found))==tuple(int(row[n]) for n in ['tp','fp','fn','k'])
    expected=math.sqrt(sum((fit['predicted'][t]-m)**2 for t,m in truth['tip_mean'].items())/len(truth['tip_mean']))
    assert math.isclose(expected,float(row['rmse']),rel_tol=1e-12,abs_tol=1e-12)
    checked+=1
manifest=json.loads((ROOT/'source-manifest.json').read_text())
with tarfile.open(ROOT/'source.tar.gz') as archive:
    assert all(hashlib.sha256(archive.extractfile(name).read()).hexdigest()==sha for name,sha in manifest.items())
(ROOT/'independent-qa.json').write_text(json.dumps({'fresh_runs':len(raw),'all_input_hashes_match':True,'score_formula_max_error':max_score_error,'fresh_metric_rows_checked':checked,'source_archive_modules':len(manifest),'status':'passed'},indent=2)+'\n')
print('Independent audit passed:',len(raw),'scores;',checked,'metrics; max score error',max_score_error)

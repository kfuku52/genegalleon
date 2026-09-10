"""Independent direct-formula and clade/prediction audit of large-search results."""
import csv,hashlib,json,math,tarfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
runs=json.loads((ROOT/'results.json').read_text())
assert len(runs)==12
metric_rows={(r['dataset'],r['method']):r for r in csv.DictReader((ROOT/'metrics.csv').open())}
max_score_error=0.;complete=0
for run in runs:
    if run['exit_code']!=0:continue
    complete+=1
    d=ROOT/'data'/run['dataset'];truth=json.loads((d/'truth.json').read_text());fit=run['result']['selected'];metrics=metric_rows[run['dataset'],run['method']]
    a={frozenset(c) for c in truth['shift_clades']};b={frozenset(c) for c in fit['shift_clades']};tp=len(a&b)
    assert len(a)==100
    assert (len(b),tp,len(b-a))==tuple(int(metrics[f]) for f in ['k','tp','fp'])
    mse=sum((fit['predicted'][n]-v)**2 for n,v in truth['tip_mean'].items())/1000
    assert math.isclose(math.sqrt(mse),float(metrics['rmse']),rel_tol=1e-12,abs_tol=1e-12)
    p=2*len(b)+3
    score=-2*fit['log_likelihood']+2*p+2*p*(p+1)/(1000-p-1)
    max_score_error=max(max_score_error,abs(score-run['result']['score']))
    assert max_score_error<1e-7
    assert all(hashlib.sha256((d/name).read_bytes()).hexdigest()==sha for name,sha in run['input_sha256'].items())
manifest=json.loads((ROOT/'source-manifest.json').read_text())
with tarfile.open(ROOT/'source.tar.gz') as archive:
    assert all(hashlib.sha256(archive.extractfile(name).read()).hexdigest()==sha for name,sha in manifest.items())
(ROOT/'independent-qa.json').write_text(json.dumps({'status':'passed','completed_runs_verified':complete,'total_runs_including_failures':len(runs),'max_aicc_score_error':max_score_error,'source_modules_hash_verified':len(manifest)},indent=2)+'\n')
print('Verified',complete,'completed large-search runs')

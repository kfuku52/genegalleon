"""Summarize the paired experiment; run in the same container as benchmark.py."""
import gzip
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).parent
rows = json.loads((root / 'results.json').read_text())
if (root / 'raw.json').exists():
    raw = json.loads((root / 'raw.json').read_text())
else:
    with gzip.open(root / 'raw.json.gz', 'rt') as source:
        raw = json.load(source)
environment = json.loads((root / 'environment.json').read_text())
index = {(r['p'], r['rho'], r['scenario'], r['mode']): r for r in rows}
fig, axes = plt.subplots(2, 2, figsize=(9, 6), sharex=True, sharey=True)
for i, scenario in enumerate(('aligned', 'opposed')):
    for j, rho in enumerate((0., .8)):
        ax = axes[i, j]
        for mode, color in [('diagonal', '#888888'), ('full', '#006b8f')]:
            values = [index[p, rho, scenario, mode]['correct_branch'] for p in (2, 5, 10)]
            errors = [1.96 * np.sqrt(v*(1-v)/environment['evaluation_reps']) for v in values]
            ax.errorbar([2, 5, 10], np.array(values)*100, yerr=np.array(errors)*100,
                        marker='o', color=color, label=mode, capsize=3)
        ax.set_title(f'{scenario} shift; correlation = {rho}')
        ax.set_xticks([2, 5, 10])
        ax.set_ylim(-3, 103)
        ax.spines[['top', 'right']].set_visible(False)
        if j == 0:
            ax.set_ylabel('Correct branch detected (%)')
        if i == 1:
            ax.set_xlabel('Number of traits')
axes[0, 0].legend(frameon=False)
fig.suptitle('100-tip OU prototype: diagonal vs full evolutionary covariance')
fig.text(.5, .005, 'Matched generating-null calibration; shared grid-estimated alpha; no observation error. Not NWKIT runtime.', ha='center', fontsize=8)
fig.tight_layout(rect=[0, .025, 1, .96])
fig.savefig(root / 'comparison.png', dpi=180)
fig.savefig(root / 'comparison.svg')

lines = ['# 100-tip evolutionary trait-covariance experiment', '',
         'This is a controlled matrix-normal OU prototype, **not a benchmark of an implemented NWKIT full-covariance shift mode**. The inspected local NWKIT native implementation supports diagonal trait covariance only. No production workflow settings were changed.', '',
         '## Protocol', '',
         '- Fixed nearly balanced 100-tip ultrametric tree, root height 1, OU fixed root. All 198 branches are searched, with zero or one shift. Exact branches, not nearby branches, count as correct.',
         '- Shared alpha is profiled over 0.1, 0.3, 1, 3, 10; generating alpha is 1. Both methods use the identical grid and candidate set. Each layout refits its diagonal or full trait covariance by ML.',
         '- 2, 5, or 10 traits; equal pairwise evolutionary correlation 0 or 0.8; marginal process variance 1. Complete observations without sampling or additional measurement errors.',
         '- Shift clades contain 10–30 tips. Tip mean displacement has total Euclidean length 2, either equally spread in the same direction across all traits, or opposite directions in the first two traits with other traits unchanged.',
         f"- {environment['calibration_reps']} separate null simulations per trait-count/correlation cell calibrate each method’s maximum likelihood-ratio statistic to nominal 5%; {environment['evaluation_reps']} new paired simulations per cell/scenario evaluate it. **Calibration uses the true generating covariance, including correlation for the diagonal method: this is an oracle comparison, not a deployable bootstrap procedure.**",
         '- Timing is one exhaustive search including covariance refits, after warmup; it excludes imports, simulation, reusable tree/design preparation, and null calibration. Method order alternates within each paired replicate. Both methods use the same batched algebra; neither calls NWKIT.',
         '- One container process on an Apple M2 Max; BLAS/OpenMP threads fixed to 1. Image local/genegalleon:nwkit-ou-auto-dev, arm64, image ID sha256:40c8c0c545ca39445a57faac6ba8879e6e1075553426c9ee6d7a8205381b56f8.', '',
         '## Detection', '',
         'Correct-branch detection percentages (diagonal → full); intervals are approximate paired 95% Monte Carlo intervals for the change in percentage points, conditional on the calibrated thresholds. They exclude uncertainty in the calibration thresholds.', '',
         '| Traits | Correlation | Shift | Correct branch, % | Difference, percentage points (95% MC interval) |',
         '|---:|---:|---|---:|---:|']
for p in (2,5,10):
    for rho in (0.,.8):
        for scenario in ('aligned','opposed'):
            a, b = [index[p,rho,scenario,m] for m in ('diagonal','full')]
            pairs = {}
            for e in raw:
                if (e['p'],e['rho'],e['scenario']) == (p,rho,scenario):
                    pairs.setdefault(e['rep'],{})[e['mode']] = int(e['selected']==e['true'])
            diff = np.array([v['full']-v['diagonal'] for v in pairs.values()]) * 100
            margin = 1.96 * diff.std(ddof=1)/np.sqrt(len(diff))
            lines.append(f"| {p} | {rho} | {scenario} | {a['correct_branch']*100:.1f} → {b['correct_branch']*100:.1f} | {diff.mean():+.1f} ({diff.mean()-margin:+.1f}, {diff.mean()+margin:+.1f}) |")
lines += ['', '![Detection comparison](comparison.png)', '', '## False positives and search time', '',
          '| Traits | Correlation | Null false positive %, diagonal → full | Median search ms, diagonal → full | Ratio |',
          '|---:|---:|---:|---:|---:|']
for p in (2,5,10):
    for rho in (0.,.8):
        a,b = [index[p,rho,'null',m] for m in ('diagonal','full')]
        med = {m: float(np.median([e['seconds'] for e in raw if (e['p'],e['rho'],e['mode'])==(p,rho,m)]))*1000 for m in ('diagonal','full')}
        lines.append(f"| {p} | {rho} | {a['detection']*100:.1f} → {b['detection']*100:.1f} | {med['diagonal']:.3f} → {med['full']:.3f} | {med['full']/med['diagonal']:.2f}× |")
lines += ['', '## Interpretation and limits', '',
          'Full covariance greatly helps this experiment’s shifts along low-variance contrasts between strongly correlated traits. It does not generally help changes aligned with their common high-variance direction. With zero correlation, estimating extra covariance entries has no consistent advantage. The effects depend on signal size and direction, tree shape, trait count and noise; these are not general accuracy estimates for GeneGalleon.', '',
          'The full model estimates p(p+1)/2 covariance entries instead of p (3 vs 2; 15 vs 5; 55 vs 10). Shared alpha and no observation error allow covariance ML in closed form. Trait-specific alpha, missing data, observation variances, multiple shifts and GeneGalleon search budgets could materially change both accuracy and runtime. The measured ratios must not be applied directly to production NWKIT wall time.', '',
          'Raw results also include an exploratory AICc score with n=100 and mean/covariance/alpha parameter counts. It does not correct for searching 198 branch locations and is not equivalent to NWKIT’s criterion. It is not the basis of the detection comparison. A change in the production covariance model would also require validation of model selection and calibration.', '',
          f"Independent dense Kronecker Gaussian likelihood check: absolute log-likelihood discrepancy {environment['oracle_absolute_error']:.3g}. Peak RSS for the whole experiment including accumulated results: {environment['max_rss_kib']/1024:.1f} MiB; this does not isolate memory by method. Docker validation only; no SIF validation.", '',
          '## Reproduce', '', '```sh',
          'docker run --rm --entrypoint python -e REPS=1000 -e CAL=2000 \\',
          '  -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 -e MKL_NUM_THREADS=1 \\',
          '  -v "$PWD/docs/benchmarks/ou-trait-covariance-100tips:/bench" \\',
          '  local/genegalleon:nwkit-ou-auto-dev /bench/benchmark.py',
          'docker run --rm --entrypoint python \\',
          '  -v "$PWD/docs/benchmarks/ou-trait-covariance-100tips:/bench" \\',
          '  local/genegalleon:nwkit-ou-auto-dev /bench/summarize.py', '```', '',
          'Files: `benchmark.py`, `summarize.py`, `results.json`, compressed paired `raw.json.gz`, `environment.json`, and plots. Seed and library versions are in environment.json.', '']
(root / 'README.md').write_text('\n'.join(lines))
with gzip.open(root / 'raw.json.gz', 'wt') as out:
    json.dump(raw, out)

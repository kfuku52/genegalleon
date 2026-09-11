"""Paired diagonal/full matrix-normal OU prototype; not NWKIT shift runtime.

100 ultrametric tips, shared alpha estimated on a declared grid, ML covariance
refit for every layout, exhaustive zero/one-shift search. No measurement error.
Calibration draws and evaluation draws are disjoint. Timing includes search,
not simulation, imports or reusable tree/design preparation.
"""
import json
import os
import platform
import resource
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import solve_triangular

N = 100
GRID = [0.1, 0.3, 1.0, 3.0, 10.0]
REPS = int(os.environ.get("REPS", "200"))
CAL = int(os.environ.get("CAL", "400"))
OUT = Path(__file__).parent
rng = np.random.default_rng(20260911)
mrca = np.zeros((N, N))
clades = []


def split(ids, depth=0.0):
    if len(ids) == 1:
        mrca[ids[0], ids[0]] = 1.0
        return
    mid = len(ids) // 2
    left, right = ids[:mid], ids[mid:]
    mrca[np.ix_(left, right)] = depth
    mrca[np.ix_(right, left)] = depth
    for child in (left, right):
        clades.append(child)
        split(child, 1 - np.log2(len(child)) / np.log2(N))


split(list(range(N)))
C = np.zeros((N, len(clades)))
for j, ids in enumerate(clades):
    C[ids, j] = 1
truth_candidates = [j for j, ids in enumerate(clades) if 10 <= len(ids) <= 30]


def kernel(alpha):
    # Fixed root, normalized to unit marginal process variance.
    return (np.exp(-2 * alpha * (1 - mrca)) - np.exp(-2 * alpha)) / (-np.expm1(-2 * alpha))


prepared = []
for alpha in GRID:
    K = kernel(alpha)
    L = np.linalg.cholesky(K)
    W = solve_triangular(L, np.eye(N), lower=True)
    x = W @ np.ones(N)
    x /= np.linalg.norm(x)
    z = W @ C
    z -= np.outer(x, x @ z)
    z /= np.linalg.norm(z, axis=0)
    prepared.append((W, x, z, 2 * np.log(np.diag(L)).sum()))
Ltrue = np.linalg.cholesky(kernel(1.0))


def fit(y, mode):
    p = y.shape[1]
    ll = np.full(len(clades) + 1, -np.inf)
    for W, x, z, logdetK in prepared:
        wy = W @ y
        residual = wy - np.outer(x, x @ wy)
        base = residual.T @ residual
        b = z.T @ residual
        S = np.concatenate([base[None], base[None] - b[:, :, None] * b[:, None, :]]) / N
        if mode == "diagonal":
            logdet = np.log(np.diagonal(S, axis1=1, axis2=2)).sum(axis=1)
        else:
            sign, logdet = np.linalg.slogdet(S)
            assert (sign > 0).all()
        current = -0.5 * (N * logdet + p * logdetK + N * p * (1 + np.log(2 * np.pi)))
        ll = np.maximum(ll, current)
    winner = 1 + int(np.argmax(ll[1:]))
    covariance_k = p if mode == "diagonal" else p * (p + 1) // 2
    k = np.full(len(ll), 2 * p + covariance_k + 1)
    k[0] = p + covariance_k + 1
    aicc = -2 * ll + 2 * k + 2 * k * (k + 1) / (N - k - 1)
    selected = int(np.argmin(aicc))
    return float(2 * (ll[winner] - ll[0])), winner - 1, selected - 1


def simulate(p, rho, scenario):
    sigma = (1 - rho) * np.eye(p) + rho * np.ones((p, p))
    y = Ltrue @ rng.normal(size=(N, p)) @ np.linalg.cholesky(sigma).T
    true = -1
    if scenario != "null":
        true = int(rng.choice(truth_candidates))
        effect = np.ones(p) / np.sqrt(p)
        if scenario == "opposed":
            effect = np.zeros(p)
            effect[:2] = [1 / np.sqrt(2), -1 / np.sqrt(2)]
        # Realized tip mean displacement, total Euclidean length 2 SD.
        y += C[:, true, None] * (2.0 * effect)
    return y, true


def timed_fit(y, mode):
    start = time.perf_counter()
    result = fit(y, mode)
    return result, time.perf_counter() - start


def oracle_check():
    # Independent dense Kronecker likelihood for one fixed-alpha fixed layout.
    p = 2
    y, _ = simulate(p, 0.8, "aligned")
    W, x, z, logdetK = prepared[2]
    wy = W @ y
    residual = wy - np.outer(x, x @ wy)
    residual -= np.outer(z[:, 10], z[:, 10] @ residual)
    sigma = residual.T @ residual / N
    V = np.kron(kernel(1.0), sigma)
    raw = np.linalg.solve(W, residual).ravel()
    dense = -0.5 * (N*p*np.log(2*np.pi) + np.linalg.slogdet(V)[1] + raw @ np.linalg.solve(V, raw))
    compact = -0.5 * (N*np.linalg.slogdet(sigma)[1] + p*logdetK + N*p*(1+np.log(2*np.pi)))
    assert abs(dense-compact) < 1e-7, (dense, compact)
    return abs(float(dense-compact))


oracle_error = oracle_check()
rows = []
raw_rows = []
for p in (2, 5, 10):
    for rho in (0.0, 0.8):
        calibration = {m: [] for m in ("diagonal", "full")}
        warm, _ = simulate(p, rho, "null")
        for mode in calibration:
            fit(warm, mode)
        for rep in range(CAL):
            y, _ = simulate(p, rho, "null")
            for mode in calibration:
                calibration[mode].append(fit(y, mode)[0])
        thresholds = {m: float(np.quantile(v, .95, method="higher")) for m, v in calibration.items()}
        for scenario in ("null", "aligned", "opposed"):
            data = {m: [] for m in calibration}
            for rep in range(REPS):
                y, true = simulate(p, rho, scenario)
                modes = list(calibration) if rep % 2 == 0 else list(reversed(calibration))
                for mode in modes:
                    (score, best, aicc), seconds = timed_fit(y, mode)
                    selected = best if score > thresholds[mode] else -1
                    entry = dict(p=p, rho=rho, scenario=scenario, rep=rep, mode=mode, true=true,
                                 selected=selected, aicc_selected=aicc, seconds=seconds, score=score)
                    raw_rows.append(entry)
                    data[mode].append(entry)
            for mode, entries in data.items():
                row = dict(p=p, rho=rho, scenario=scenario, mode=mode, replicates=REPS,
                           threshold=thresholds[mode],
                           detection=float(np.mean([e['selected'] >= 0 for e in entries])),
                           correct_branch=float(np.mean([e['selected'] == e['true'] for e in entries])),
                           wrong_branch=float(np.mean([e['selected'] >= 0 and e['selected'] != e['true'] for e in entries])),
                           aicc_detection=float(np.mean([e['aicc_selected'] >= 0 for e in entries])),
                           median_seconds=float(np.median([e['seconds'] for e in entries])))
                rows.append(row)
            print(json.dumps(rows[-2:]), flush=True)
        (OUT / "results.json").write_text(json.dumps(rows, indent=2) + "\n")
(OUT / "raw.json").write_text(json.dumps(raw_rows) + "\n")
environment = dict(numpy=np.__version__, scipy=scipy.__version__, platform=platform.platform(),
                   processor=platform.processor(), seed=20260911, tips=N, alpha_grid=GRID,
                   calibration_reps=CAL, evaluation_reps=REPS, oracle_absolute_error=oracle_error,
                   max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                   thread_env={k: os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')})
(OUT / "environment.json").write_text(json.dumps(environment, indent=2) + "\n")

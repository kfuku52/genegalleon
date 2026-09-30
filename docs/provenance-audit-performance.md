# Full provenance audit performance and progress

`artifact_provenance.py audit` preserves the TSV columns, sorting, scientific
hashes, stale-policy classifications and CSUBST branch checks. It hashes shared
declared sources once during classification, then rereads every unique source
and manifest before publishing the report. Persistent digests from earlier runs
cannot substitute for this audit's content reads. Final content checks are needed
because shared-filesystem attributes alone did not detect a rapid same-size
rewrite in the actual SIF regression fixture.

One producer read lock and an initial/final metadata fence replace repeated
archive-index scans. Raw/ZIP source selection, pending-index, symlink, tombstone
and generation checks remain enforced. Thread pools and their submission windows
are bounded; the default is up to four workers, capped by `GG_TASK_CPUS`.
`--workers 1` selects serial operation. Run-scoped caches hold at most 250,000
sources/guards; a guard-capacity overflow fails instead of dropping a check.
The audit holds up to eight ZIP readers per thread and closes them at the fence.

## Evidence, 2026-09-30

Synthetic fixtures repeat a 24 MiB shared table across families and a 16 KiB
family table across 15 manifests. Each run uses Linux Python 3.12.14 in the same
`local/genegalleon:dev` Docker runtime, limited to 4 CPUs/4 GiB. The baseline is
GeneGalleon `63efea916bee99eeae98d2ae30ead62b374b91f4`. Separate warmups precede
three measurements per version; baseline trials use an ABBAAB sequence.
Run-scoped RAM caches start empty; the baseline retains its persistent digest
cache warmed by prior trials. The new full audit intentionally ignores those
persistent digests. OS caches are not forcibly cleared.

| Layout | Families / manifests | Baseline median | New median | Speedup |
|---|---:|---:|---:|---:|
| Raw | 1,000 / 15,000 | 45.54 s | 8.33 s | 5.47× |
| ZIP | 300 / 4,500 | 90.91 s | 2.78 s | 32.66× |
| Mixed raw/ZIP | 300 / 4,500 | 63.57 s | 3.17 s | 20.04× |

Every baseline/candidate TSV SHA256 matches within its fixture. These are
complete fixture audits; the fixture has no CSUBST branch data. Existing branch
identity regressions exercise those unchanged scientific checks separately.
These numbers are not an ETA for a running scientific job.

Six 6,815-family / 102,225-manifest runs took 60.30–65.18 seconds with a maximum
262.65 MiB peak RSS. Every TSV SHA256 matched. Instrumented atomic progress
publication consumed at most 0.27% of wall time. Progress-on/off medians were
64.33/61.77 seconds, with overlapping ranges: this jitter does not independently
establish an end-to-end overhead bound. The 0.27% measurement covers snapshot
publication rather than total filesystem traffic or all audit bookkeeping.

Reproduce private fixtures outside any scientific workspace:

```bash
python -B workflow/benchmarks/benchmark_provenance_audit.py \
  --workspace /tmp/gg-audit-benchmark-raw --create
python -B workflow/benchmarks/benchmark_provenance_audit.py \
  --workspace /tmp/gg-audit-benchmark-zip --create --families 300 \
  --layout zip --support-root "$PWD/workflow/support"
python -B workflow/benchmarks/benchmark_provenance_audit.py \
  --workspace /tmp/gg-audit-benchmark-raw --support-root workflow/support \
  --workers 4 --label trial-1
```

Fixture creation refuses an existing destination. For baseline measurements,
extract its immutable `workflow/support` tree and pass that directory; omit
`--workers` because the old CLI has no such option. Keep the runtime, CPU/memory
limits, fixture, warmups and cache conditions comparable.

## Machine-readable observations

The default `--progress-interval 10` writes atomic owner-private
`REPORT.tsv.progress.json`; zero disables progress snapshots. Phases are
inventory, manifest hashing, legacy inventory, branch identity, source
revalidation and report publication. Counters and ETA describe only the current
phase; heartbeat time and last advancement time are separate.

`REPORT.tsv.result.json` records source-closure, manifest-inventory and report
SHA256, status counts, exit code, elapsed time and bounded cache counters.
`digest_bytes_read` counts declared-source hash stream bytes, including the final
reread; it excludes manifest JSON and physical/compressed filesystem traffic.
`progress_write_ns` measures progress publication time.

Observed jobs additionally write `progress-artifact_audit.json` and
`result-artifact_audit.json` beside their exact `run.json`.
The optional [workflow API](workflow-inspection-api.md) progress capability
returns these as attempt evidence, with strict bounds and identity validation.
Neither progress completion nor a successful audit result proves that the full
workflow or database generation finished. Existing jobs retain their immutable
source and observation coverage.

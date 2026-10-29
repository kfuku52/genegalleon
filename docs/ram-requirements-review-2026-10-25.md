# RAM requirements review (2026-10-25)

Twelve separate subagents reviewed the following proposals against baseline
`610adefe2d432ca74085496d8db82a9c51f62b50` (0.8.408). Reviews covered callers,
validation boundaries, output ordering, and measurement design before selecting
changes for implementation. Curated inputs and scientific settings are outside
this change.

## Twelve proposals

| # | Proposal | Review and decision |
| --- | --- | --- |
| 1 | Read the GFF attribute identifier in bounded chunks | Rejected: pandas 3.0.6 silently accepted a tenth field at the start of the next chunk. The baseline rejects it; consuming all chunks is insufficient. |
| 2 | Aggregate the best BLAST hit per gene while reading chunks | Deferred: whole-file identifier inference, missing/infinite e-values, and stable tie ordering need a dedicated equivalence study. |
| 3 | Stream Orthogroup conversion and gene counts | Rejected: the same parser defect accepted a malformed late row and published output. A dependency fix is required before reconsideration. |
| 4 | Spool contamination-filtered FASTA records to disk | Accepted: 71% lower peak RSS, with about 8% more time and additional temporary disk I/O. Preserve last-wins duplicate IDs and repeated-record renaming; the existing sequence store is not interchangeable. |
| 5 | Retain CDS identifiers/lengths rather than sequence bodies for structural-tip export | Accepted: 46–54% lower peak RSS in the large fixtures, with complete output equivalence. Keep the public alignment reader and complete bundle validation. |
| 6 | Bound SQL input buffers and/or externally sort global BH-FDR | Deferred: chunking SQL alone leaves BH arrays proportional to all tests; external sorting trades RAM for I/O and must preserve the global test family, rowids, numeric validation and rollback. |
| 7 | Stream scaffold-taxonomy context construction | Deferred: global gene/rank/locus consistency checks need compact cross-chunk state or disk spill. Dense output contexts remain large. |
| 8 | Project only five synteny gene-info columns | Deferred: plain `usecols` can suppress malformed-row errors. Full validation plus projection needs a RAM/time comparison of its own. |
| 9 | Retain only family-relevant WGD classification evidence | Accepted: 87% lower peak RSS for a 100-gene family. Validate every row and all gene/event identities; retain invalid/unmatched anchors used in conflict detection. The dense control has a time cost. |
| 10 | Retain selected copy-number families and species columns | Deferred: useful for small explicit selections; the default all-family analysis still needs a large matrix. Preserve the current validation scope for unselected numeric cells. |
| 11 | Decode trait-index SQL rows without `fetchall()` | Deferred: sparse family queries already select few rows and wide payloads are compressed. Dense-query measurements are needed before prioritizing this smaller duplicate-buffer saving. |
| 12 | Scan the OrthoFinder completion log without `read_text().splitlines()` | Deferred: no representative large production log was established. A streaming reader must still consume malformed UTF-8 after a marker and preserve Unicode line separators. |

## Measurement boundaries

Fixtures are synthetic and created outside fresh measured workers. Comparisons
use the same GeneGalleon Docker image, alternate baseline/current order, discard
one warmup, and compare complete outputs as well as wall time, CPU time and Linux
process peak RSS. This measures the changed helper or evidence-loading stage,
not total workflow RAM, scheduler requests, production storage throughput, or
SIF compatibility. The WGD dense control retains the entire family scope.

The Orthogroup benchmark now creates CLI fixtures in the controller and hashes
outputs in bounded blocks. Its former worker retained the fixture DataFrame,
which could obscure the formatter's actual RAM requirements. Kernel profiles
retain their original measurement scope.

## Rejected chunk readers

Before production edits, all 131 selected fast tests passed, including the new
input-validation guards. After the GFF chunk implementation, both plain and gzip
cases failed to raise `ParserError` for the ninth row with a chunk size of eight.
After the Orthogroup implementation, a malformed seventh data row with chunk size two
was accepted and published. The production changes were withdrawn. These are
instances of the dependency's [chunk-boundary parsing issue](https://github.com/pandas-dev/pandas/issues/61973),
not justification for a weaker input contract or a local parser workaround.
The guards and fixture-generation corrections remain useful for future reviews.

## Environment and final measurements

Freshness verification passed through the normal runtime wrapper, with no
freshness bypass. Docker Linux/aarch64 used Python 3.12.15, pandas 3.0.6 and NumPy
1.26.4, image
`sha256:340685eba10a8bf74ad03b85b1315107bdfb8dba84a37dc921b96b3e18298339`.
Docker Desktop provided 12 CPUs and 8,214,585,344 bytes of RAM. Benchmarks had no
additional container limit except the explicit memory-limit qualification below.
Our measured benchmarks ran serially, separately from our builds and tests.
Other work also ran on this shared desktop during some samples; peak RSS was
stable across fresh workers. Timings expose local tradeoffs and are not
whole-workflow speed claims. Raw sources, logs, samples and hashes remain in
`tmp/ram-review-20261025-7Fud2j/`.

Each table reports three measured trials after one discarded warmup, with
alternating baseline/current order. RSS is MiB (2^20 bytes). GFF timings include
imports; the other operation timers exclude imports, fixture creation and output
hashing. RSS includes worker imports. Comparisons use complete output hashes,
not summaries of selected numerical values.

### Structural-tip export

The large fixtures contain 4,096 tips with 5,120 codons each (60 MiB of CDS
sequence), with zero or 4,095 additional structural ancestral sequences. The
small control has 16 tips and 100 codons. Complete bundle/member verification
remains in the measured call. Every exported FASTA byte matched.

| Profile | Peak RSS before → after | Reduction | Wall before → after | CPU before → after |
| --- | ---: | ---: | ---: | ---: |
| Small | 52.54 → 52.54 MiB | 0% | 0.0006 → 0.0006 s | 0.0006 → 0.0006 s |
| Tips only | 133.73 → 71.61 MiB | 46.4% | 0.3081 → 0.3105 s | 0.3081 → 0.3105 s |
| Tips + ancestors | 157.01 → 72.26 MiB | 54.0% | 0.3801 → 0.3772 s | 0.3801 → 0.3772 s |
| Tips + ancestors, gzip CDS | 157.02 → 72.30 MiB | 54.0% | 0.5201 → 0.4390 s | 0.5201 → 0.4390 s |

No speed improvement is claimed for the sub-millisecond control or the plain
large profiles: their paired after/before timing ranges cross unity. Necessary
tip-state sequences and the largest current FASTA record still occupy RAM.

### WGD family evidence

The fixture contains 200,000 gene positions, 500,000 anchors and 32 events. The
sparse family selects 100 gene identities and 515 anchor rows (105 pairs).
Every retained field, row order and duplicate anchor matched. Real classifier
tests also compare unchanged classification TSVs, pair evidence and NHX output
after adding unrelated evidence, and retain invalid-anchor conflict behavior.

| Scope | Peak RSS before → after | Reduction | Wall before → after | CPU before → after |
| --- | ---: | ---: | ---: | ---: |
| 100-gene family | 603.21 → 75.79 MiB | 87.4% | 1.7838 → 1.1358 s | 1.7836 → 1.1356 s |
| All 200,000 genes | 633.24 → 621.85 MiB | 1.8% | 1.9519 → 2.3659 s | 1.9518 → 2.3643 s |

The dense control takes about 21% longer: scoped retention is useful for ordinary
family-sized queries, not a time optimization for a tree covering every input
gene. Global identity validation remains proportional to gene count during the
position scan; that set is released before anchor loading. Complete pre/post
input-hash fences and all scientific evidence rules remain in place.

### Contamination-filtered FASTA

The large input has 10,000 records and 409,120,396 sequence bases; 9,500 records
(388,664,315 bases) are retained. A fixed two-class taxonomy replaces database
lookup for the benchmark; the real `main()` still performs parsing, joins,
sorting, lineage classification, renaming and both output writes. Every FASTA
and lineage-compatibility TSV byte matched. The reject-most control retains ten
records; the small control has 32 records of at most 2,048 bases.

| Profile | Peak RSS before → after | Reduction | Wall before → after | CPU before → after |
| --- | ---: | ---: | ---: | ---: |
| Retain-most, plain | 529.42 → 154.76 MiB | 70.8% | 1.5717 → 1.6925 s | 1.5611 → 1.6925 s |
| Retain-most, gzip | 529.56 → 154.81 MiB | 70.8% | 1.8424 → 1.9797 s | 1.8409 → 1.9784 s |
| Reject-most | 154.25 → 154.12 MiB | 0.1% | 0.4485 → 0.4470 s | 0.4485 → 0.4470 s |
| Small | 148.42 → 148.71 MiB | −0.2% | 0.0116 → 0.0095 s | 0.0116 → 0.0095 s |

This is a RAM/I/O tradeoff. The retained plain fixture takes about 8% longer.
Logical process reads increase from 391.14 to 836.44 MiB and writes from 377.62
to 750.44 MiB, measured with `/proc/self/io` outside the operation timer. These
counters include buffered/cache I/O and are not physical-device throughput.
Temporary space holds every compatible FASTA occurrence plus pickle metadata,
including records later superseded by a duplicate ID. Only self-generated
SeqRecords are deserialized, within one private temporary-file context.

The remaining RAM is the taxonomy/metadata tables, ID indexes and work buffers
proportional to the largest FASTA record. **Use disk-backed `TMPDIR` storage**:
tmpfs would charge the stored sequences to node/cgroup RAM instead. The tested
Docker `/tmp` was disk-backed; no workflow temporary-path setting was changed.

For an actual allocation check, the fixture creator, worker and filesystem cache
shared a 256 MiB cgroup limit with swap disabled. The baseline worker was killed
by the memory controller (`SIGKILL`, `oom_kill=1`). The current worker completed
with `oom_kill=0` and identical FASTA/TSV hashes to the unlimited baseline. Both
cgroups reached the limit because filesystem cache also counts; success under
the limit, rather than a lower cache-inclusive peak, is the qualification.
This does not guarantee 256 MiB for arbitrary genomes with larger records.

## Reproduce and validate

Create saved baseline helpers in a new ignored directory, then use a runtime
that passes the normal freshness check:

```bash
mkdir -p tmp
review_dir="$(mktemp -d tmp/ram-review.XXXXXX)"
review_dir="$(cd "${review_dir}" && pwd)"
for helper in csubst_input_bundle.py wgd_ssd.py remove_contaminated_sequences.py; do
  git show 610adefe2d432ca74085496d8db82a9c51f62b50:workflow/support/"${helper}" > "${review_dir}/${helper}"
done
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_csubst_bundle_memory.py \
  --baseline "${review_dir}/csubst_input_bundle.py" --output "${review_dir}/structural.json"
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_wgd_evidence_memory.py \
  --baseline "${review_dir}/wgd_ssd.py" --output "${review_dir}/wgd.json"
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_contamination_memory.py \
  --baseline-source "${review_dir}/remove_contaminated_sequences.py" \
  --scenario retain-most --codec plain --output "${review_dir}/contamination.json"
```

Use `--profiles ancestors --gzip-cds` for compressed structural inputs. For FASTA,
repeat with `--codec gzip`, `--scenario reject-most` and `--scenario small`, each
with a new output path. All comparison scripts support `--baseline-only` for a
pre-implementation baseline. `benchmark_gff_attribute_memory.py --include-errors`
and `benchmark_gene_counts.py --profiles hog2og sp2og --cli-rows 100000` retain the
rejected readers' qualification workloads.

The bounded-memory experiment uses `qualify_contamination_memory.py` in the same
freshness-checked image, externally constrained with Docker `--memory=256m
--memory-swap=256m` and the checkout mounted at its absolute path. Pass the saved
baseline with `--expect oom`; pass the current source with `--expect fits
--reference-report <contamination.json>`. Each run needs a new `--output` path.
The script verifies the actual cgroup limit, disabled swap, OOM counters and,
when given the reference report, the complete outputs.

## Verification record

- `bash ./dev check fast workflow/tests/test_gff_attribute_identifier.py workflow/tests/test_orthogroup_table_formatter.py workflow/tests/test_csubst_input_bundle.py workflow/tests/test_csubst_site_wrapper.py workflow/tests/test_wgd_evidence.py workflow/tests/test_remove_contaminated_sequences.py -x`: **209 passed**, no skips, Docker.
- `GG_PYTEST_WORKERS=1 bash ./dev check runtime-python workflow/tests/test_wgd_ssd_contracts.py workflow/tests/test_wgd_ssd_runtime.py workflow/tests/test_wgd_ssd_gene_stage.py -x`: **174 passed**, no skips, Docker. Includes real native count/synteny/Ds/classification and wrapper/cache checks.
- `bash ./dev lint`: passed with the documented Docker Bash parser wrapper selected by `GG_LINT_BASH`; Ruff and configuration validation ran on the host.
- `git diff --check`: passed.
- Bounded-memory baseline OOM and current success/output equivalence: passed with the published qualification script, swap disabled.

Runtime-wrapper checks selected the immutable image above with
`GG_TEST_RUNTIME=docker` and `GG_CONTAINER_DOCKER_IMAGE`. The allocation check
used that same image directly after freshness verification, with Docker memory
and swap limits. The two deliberately
rejected candidate test failures are recorded separately and their production
changes are absent. SIF and full scientific workflow validation were not run;
these results do not establish SIF compatibility or production genome accuracy.

## Follow-up: literal identifiers and GFF scanning (2026-10-29)

GFF attribute identification now scans tab-delimited records directly. Unlike
the rejected pandas chunk experiment, it checks every annotation row's nine
columns even after finding a key and consumes gzip/UTF-8 input through EOF.
Only whole-line comments are skipped, preserving literal hashes in attributes.

Using the same qualified Docker image, Python 3.12.15 / pandas 3.0.6 on Linux
ARM64, `benchmark_gff_attribute_memory.py --rows 500000 --small-rows 32
--fasta-records 10000 --positions first last missing --codecs plain gzip
--repeats 3` compared baseline `b055605a` with the direct reader. Each variant
had one discarded warmup, alternating fresh workers and identical selected-key
results. The plain GFF was about 58.7 MB. Across the six large cases, median
peak RSS fell from 328–331 MiB to 76–78 MiB (76–77% lower), and wall time
fell by 32–49%. For example, a final-row match in gzip took 1.021 → 0.652 s
with 329.8 → 78.3 MiB. Small inputs remained around 0.17–0.20 s.

The error qualification also passed for extra fields and truncated gzip. It
compares ParserError class, row and expected/observed column counts while
retaining complete diagnostic messages; other outcomes compare exactly.

The shared literal-hash reader adds a small annotation-loading cost. In the
existing BUSCO aggregation benchmark (three repeats), 50,000 single-row genes
took 0.075 → 0.089 s; the 10,000-gene merge CLI took 0.099 → 0.109 s. Complete
typed frames, TSV/CLI bytes, errors and logs matched for coordinate-suffix
fixtures. Literal colon namespaces, hashes and multiline comments have separate
regression tests because their corrected outputs intentionally differ.
These measurements concern the named helpers, not complete GPU workflows.

# Copy-number and CAFE validation

Gene-family copy numbers are finite, nonnegative integers. A zero is an observed
absence; missing, fractional, Boolean or infinite values are not converted to
zero. Species columns and family identifiers must be unambiguous, and TSV rows
must have the declared width. Literal identifiers such as `001` are preserved.

Validate selected gene counts without creating plots or a prepared matrix:

```bash
python workflow/support/prepare_orthogroup_copy_number.py \
  --genecount Orthogroups.GeneCount.selected.tsv \
  --dated_species_tree dated_species_tree.nwk --validate-only
```

The workflow performs this check before preparing or reusing a copy-number
matrix. This matters because an old cached matrix must not bypass corrected
input validation. In particular, pandas can rename duplicate species headers
or infer an index from extra cells, and native CAFE can truncate fractional
counts. These inputs are rejected before producing output.

## Reusing a saved CAFE fit

```bash
python workflow/support/validate_cafe_outputs.py \
  --copy-number orthogroup_copy_number.tsv \
  --tree dated_species_tree.nwk --output-dir cafe_output
```

Before adopting a manifest-free fit, before recording a new fit, and when reusing
a recorded fit, the workflow checks species and topology, branch/family
identities, and agreement of reconstructed extant counts with the supplied
matrix. Foreign or duplicate families, malformed changes and ambiguous headers
are rejected. Input families legitimately omitted by CAFE, such as root-absent
families, are reported as omissions and are not imputed as zeros. Native
`FamilyID` and `#FamilyID` headers are supported.

This semantic check cannot establish the gamma categories, optimisation settings
or other historical parameters of an unrecorded fit. Preserve its original
provenance. Explicit `artifact_legacy_policy=rebuild` still requests a new fit;
validation does not prevent that intentional rebuild.

## GO analysis

Both event methods reject duplicate or foreign probability IDs, inconsistent
branch columns, probabilities outside [0,1], malformed reported probabilities,
and nonfinite or fractional reconstructed changes. CAFE's absent probability
rows remain unreported, not zero-probability events. The legacy all-missing
zero-change root-column convention remains supported. Failed invocations remove
their own prior GO success summaries, retaining other evidence.

GO readers preserve family and gene identifiers literally and reject shifted
TSV rows. The frozen GO input remains an input, not an output to overwrite.
See [GO enrichment](go-enrichment.md) for the distinction between event counts,
family flags, BH FDR and Tarone FWER. These validation fixes do not change valid
inputs' event model, cutoffs, raw P values or correction scope.

## Summary aggregation performance

CAFE summary counts now use one reference-tree geometry, indexed change rows,
column sums and significance labels. They no longer fortify each family's tree
or scan all change rows and merge an aggregate for every family. All native
branch counts and reference geometry are retained; inconsistent branch
identities fail explicitly.

The reproducible synthetic workload is
`workflow/tests/benchmark_cafe_summary.R`: 24 species, 1,500 families, fixed seed,
one warmup and three timed aggregations. Run it from the repository root in a
qualified GeneGalleon runtime, passing either source revision:

```bash
Rscript workflow/tests/benchmark_cafe_summary.R saved_baseline_summary.r before.tsv
Rscript workflow/tests/benchmark_cafe_summary.R workflow/support/cafe_plot_summary.r after.tsv
cmp before.tsv after.tsv
```

In the same Docker image (`sha256:00e94be099f4ff2e0816ef02644d4187b57d08df79722476d81cfdbf382a3a27`),
the original aggregation took 10.831/11.868/11.627 seconds and the indexed
aggregation took 0.103/0.079/0.076 seconds. The medians are 11.627 and 0.079
seconds, about 147 times faster for this aggregation workload. A separate
full-workload comparison verified exact equality of all output columns,
including coordinates, with the original function. This is not a speedup
measurement for CAFE inference, individual-family rendering or an entire study.

The targeted tests include a real, small Gamma CAFE invocation and read-only
reuse checks; Docker validation is not evidence of SIF or project-data validation.

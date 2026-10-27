# Read k-mer profiles with GenomeScope

Enable `run_genomescope=1` in `workflow/gg_genome_annotation_entrypoint.sh`.
Provide canonical species FASTQ directories as for [Smudgeplot](smudgeplot.md).
The standard runtime installs KMC, GenomeScope and Smudgeplot. Enabling
GenomeScope automatically prepares Smudgeplot from the same canonical reads
and k before fitting, even when `run_smudgeplot=0`. The default initial lambda
is Smudgeplot's native `haploid_coverage` (H), rounded to the nearest positive
integer (half values round up, matching GenomeScope's integer CLI).

| Parameter | Default | Meaning |
| --- | --- | --- |
| `genomescope_kmer_length` | 21 | Canonical k-mer length |
| `genomescope_kmer_lower` | 1 | Minimum count retained |
| `genomescope_kmer_upper` | 1000 | Counter and histogram ceiling |
| `genomescope_ploidy` | 2 | Model ploidy, 1 to 6 |
| `genomescope_initial_coverage` | -1 | Positive manual lambda takes precedence; -1 uses the selected source |
| `genomescope_coverage_source` | `smudgeplot` | Native H, or `auto` for GenomeScope's own initialization |
| `genomescope_max_kmercov` | -1 | Model coverage limit, or complete spectrum with -1 |
| `genomescope_histogram` | empty | Existing two-column histogram, bypassing KMC |
| `genomescope_program` | `/usr/local/bin/genomescope2.0/genomescope.R` | Runtime executable |

Scoped overrides use the `GG_GENOME_ANNOTATION_` prefix and upper-case variable
names. Guided runs retain H/2, H, 2H and automatic initialization as sensitivity fits
(duplicate rounded starts are omitted). The primary fit always uses H; the
pipeline does not choose a result by reference assembly size or fit score.
Model ploidy remains the explicit `genomescope_ploidy` setting, independent of
Smudgeplot's inferred haplotype structure.
For example, set `GG_GENOME_ANNOTATION_GENOMESCOPE_KMER_UPPER=10000000` to retain
the high-frequency tail. KMC uses canonical counting without homopolymer
compression. Its maximum-count exclusion is disabled; counts above the ceiling
are retained as saturated boundary counts. `histogram_metadata.json` reports
that boundary. A zero boundary count is required before treating the new KMC
histogram as uncensored. Histogram-file inputs retain their original counter's
limitations, which cannot be inferred from the two columns alone.

Set both k-mer lengths alike if changing k. Histogram-only input cannot verify
matching Smudgeplot reads: supply a positive manual lambda or set
`genomescope_coverage_source="auto"` for this use case. Unusable Smudgeplot H or
failed inference stops guided fitting; inspect the histogram and choose an
explicit Smudgeplot cutoff when necessary (see [limitations](smudgeplot.md)).
There is no automatic cutoff adjustment or fallback to another initial lambda.

The archive at `output/species_dnaseq_genomescope/<species>_genomescope.zip`
retains the original `kmc.histo.tsv`, a dense `model.histo.tsv` with a final
empty bin, histogram mass/boundary metadata, the commands, native JSON reports and model outputs. `coverage_selection.json`
records raw H, rounded lambda, source archive SHA-256, input identities, cutoff,
fixed model ploidy and each fit's status. The primary outputs retain their
existing archive-root paths; alternatives live under `sensitivity/`. Failed
sensitivity fits retain their logs and status without replacing the primary
fit. Missing outputs and malformed native reports include a failure reason in
`coverage_selection.json`. The fit helper refuses pre-existing model outputs
before running any fit, preventing stale results from being mistaken for a new
successful run. A failed primary fit leaves diagnostics in the task directory
and does not replace the published archive or provenance.
The final empty bin prevents GenomeScope's terminal-row removal from losing the
highest observed multiplicity. Existing archives are replaced atomically only
after a primary model, JSON report and plot exist. Changes to the Smudgeplot
archive invalidate the GenomeScope artifact through the normal provenance and
stale-artifact policy. Existing artifacts may require an explicit rebuild when
upgrading to this initialization protocol. Successful convergence alone does not establish
identifiable genome size, heterozygosity or ploidy; inspect sensitivity and fit
diagnostics before reporting biological estimates.

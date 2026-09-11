# Alignment statistics

Gene-evolution workflows generate original and cleaned alignment summaries with
`cdskit stats --mode alignment`. CDSKIT 0.31.0 or newer is required; rebuild older
GeneGalleon containers before using this workflow. Container source defaults
continue to follow the moving branches in `container/source_branches.env`.

Both DNA (`--seq_type dna`) and protein (`--seq_type aa`) inputs are supported.
seqkit converts the compressed alignment to a temporary FASTA, and cdskit writes
an AMAS-compatible TSV before the workflow publishes it to the family output.
AMAS is no longer a runtime dependency. Other seqkit operations remain in use.

For existing workspaces, the `run_amas_original` and `run_amas_cleaned` switches,
`amas_original` / `amas_cleaned` directories, TSV filenames, and summary columns
retain their names. The orthogroup and query2family summary readers also retain
support for existing ZIP-backed and legacy dot-named outputs. Provenance now
records `statistics_engine=cdskit-stats-alignment`, invalidating older manifests
when the stage runs again.

`GC_content` retains the AMAS-compatible 0–1 scale and per-taxon averaging;
it is not a pooled GC percentage. Protein summaries report this field as `NA`.
See the [CDSKIT stats documentation](https://github.com/kfuku52/cdskit/wiki/cdskit-stats)
for missing-character and informative-site definitions.

## Migration validation

Validated with a Docker-backed GeneGalleon image derived from
`local/genegalleon:nwkit-ou-auto-dev`, installing CDSKIT 0.31.0 from its committed
source and removing the `AMAS.py` executable. The modified workflow source was
mounted from the checkout.

- Four runtime cases execute the actual original/cleaned stage blocks on gzip
  DNA/protein FASTA and read each TSV with both summary readers.
- 419 related tests passed: orthogroup/query2family summaries, family output
  storage and ZIP compatibility, FASTA contracts, shell static safety, and
  container build entrypoint checks.
- Bash syntax checks passed for the core and entrypoint scripts.
- Wrapper dry-run covered eight configurations; the orthogroup configuration
  was skipped because this workspace lacked the selected gene-count table.

This was focused Docker validation, not a complete workflow run, a clean
multi-platform container build, or SIF/Apptainer validation.

# BUSCO-filtered input export plans

`plan_input_cohort_export.py` creates English metadata for a cohort from a
finalized native species summary. It selects species using exact
`Complete + Duplicated` BUSCO group counts, not rounded displayed percentages
or the number of duplicate gene hits. Short/full summaries must agree and
the native BUSCO provenance must match the current CDS and QC file hashes.
All species must use one named lineage and CDS transcriptome BUSCO mode.
This cutoff describes CDS completeness, not genome assembly completeness.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/support/plan_input_cohort_export.py \
  --workspace-root /data/project/workspace \
  --species-summary /data/project/workspace/output/input_generation/gg_input_generation_species.tsv \
  --busco-short-dir /data/project/workspace/output/input_generation/species_cds_busco_short \
  --busco-full-dir /data/project/workspace/output/input_generation/species_cds_busco_full \
  --busco-provenance-dir /data/project/workspace/output/input_generation/artifact_provenance \
  --minimum-complete-busco 50 --target-part-gb 20 \
  --metadata-table traits=/data/project/workspace/input/species_trait/species_trait.tsv \
  --output /data/project/cohort-export-plan
```

The output directory must be new. Selected species need a nonempty gzip CDS,
GFF and genome. Each triplet stays in one planned part. Payload totals use
decimal GB and exclude tar headers/padding; a species larger than the target
gets its own part. Source files are freshly hashed with mutation fences and
are not copied, rewritten, recompressed or uploaded. Planning is an additional
handoff step and does not replace the source workflow's scientific QC.

Outputs:

- `cohort.tsv`: species, TaxID, exact BUSCO counts/percentage, lineage/version,
  planned part, triplet size and relative CDS/GFF/genome paths.
- `files.tsv`: recipient-relative paths, roles, part, compressed size and SHA256.
- Optional metadata TSVs: only selected species, preserving English headers
  and source values. Each source needs exactly one `species` or `species_prefix`
  column, unique labels and every selected species. Add taxonomy with another
  `--metadata-table taxonomy=PATH` if the native taxonomy table fits this contract.
- `export_plan.json`: internal absolute paths and hashes for future packaging,
  part membership/size and excluded species. Keep this internal when sharing
  the recipient tables.
- `README.txt` and `SHA256SUMS`: explanation and metadata integrity checks.

The planner requires exact native BUSCO provenance. Legacy results with renamed
IDs or different compressed CDS representations need a separately reviewed
adoption step; sequence-only agreement is not silently accepted here. Existing
plans are never overwritten. Archives and transfer-service uploads remain
separate actions.

# Aggregate verified formatted input cohorts

`workflow/support/input_generation_cohort.py` combines selected successful native
workers from multiple input-generation workspaces. It preserves compressed CDS,
GFF and genome bytes, sequence IDs, transcript ORF audits, experimental projection
proofs and cloned-model history. It never reruns the formatter or collapses
experimental sequence groups. GFF and genome remain optional for species whose
reviewed source has no such output.

The selection TSV has `species`, `source_workspace`, `plan_path`, `plan_sha256`,
`cds_path`, `gff_path`, `genome_path` and optional `cloned_model_audit_path`.
Each selected species must have exactly one matching native plan task and a
completion receipt binding its summary shard and all selected outputs/evidence.
A published multispecies summary is not required for unrelated donor species.

```bash
python workflow/support/input_generation_cohort.py plan \
  --selection selected_sources.tsv --expected-species 383 --plan cohort_plan.json
python workflow/support/input_generation_cohort.py prepare \
  --plan cohort_plan.json --workspace /data/project/cohort/workspace
```

Prepare takes native donor locks, freshly verifies every receipt-bound file,
and creates regular-file views using hardlinks where supported and bounded
streaming copies across filesystems. It verifies content again after publication.
Donor metadata and evidence remain required: retain those workspaces and treat
all linked inputs as read-only. No stat-only or persisted checksum establishes
success. A partially prepared identical view can resume; conflicting files or
changed selection/proofs require investigation.

The view provides standard `input/species_cds`, `input/species_gff` and
`input/species_genome` directories, native BUSCO/fx2tab artifacts under
`output/input_generation`, the selected original summary rows, byte-identical
cloned audit snapshots and a `cohort_view.json` provenance index. Original summary
paths remain in their donor namespaces so source-gene and ORF checks retain their
meaning. BUSCO artifact contracts remain usable against byte-identical files.

Run `qc --task-index N --plan cohort_plan.json --workspace ...` once per species,
using the standard longest-CDS/source-ownership and strict CDS/GFF/reference
validators. When saved rows contain `/workspace` paths, mount that exact donor
workspace at `/workspace` inside the same qualified container. Other project
paths must remain visible at their original absolute locations. An absent or
wrong donor mount is an error; ORF audits are never rewritten to waive checks.

Run `finalize --plan ... --workspace ...` only after every QC worker succeeds.
It verifies the exact cohort/view and each species result, writes `cohort_qc.json`
and a BUSCO `marker_occupancy.tsv` (S/D/F/M, `?` for ambiguous Complete rows,
`.` for absent marker evidence). Marker inventories/lineages must agree; taxa
with no single-copy markers are reported without automatic exclusion. Scheduler
acceptance and a prepared view do not certify final QC completion.

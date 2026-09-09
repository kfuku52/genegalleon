# Intron ancestral-state reconstruction

`gg_gene_evolution_entrypoint.sh` uses `nwkit asr` to infer the probability
of intron presence at every node of a dated gene tree. Enable it with
`run_asr_intron=1`. The default is `0`; `gene_evolution_profile="hgt"`
enables it together with GFF collection. The stage requires the family's
GFF statistics and dated analysis tree, after any tree pruning.

## Model and observations

The GFF table supplies `gene_id` and `num_intron`. Zero means absent, a positive
integer means present, and a missing count or absent GFF row means unknown.
Unknown states are inferred by ASR; they are never converted to observed zeros.
Duplicate gene IDs and invalid counts are errors.

The fixed generator, in absent/present order, is:

| From / to | Absent | Present |
|---|---:|---:|
| Absent | −g | g |
| Present | l | −l |

- `g = intron_gain_rate`, default `0.0001`.
- `l = retrotransposition_rate`, default `0.001`.
- Both rates are per unit of the input dated-tree branch lengths.
- `nwkit asr --model CUSTOM --root-prior equal` fixes the root prior at
  `(0.5, 0.5)` and calculates marginal posterior probabilities. Rates are
  not re-estimated, and no stochastic simulation count or random seed is needed.
- The dated tree is used as supplied. ASR does not pad short edges, rescale
  lengths, or force an ultrametric tree. Upstream tree dating owns the time scale.

Even when every observed gene has the same state, ancestral and unknown tip
states can remain uncertain. Their probabilities are not forced to 0 or 1.

## Outputs

Paths below are relative to the active `query2family` or `orthogroup` output root.

| Directory | Family file | Contents |
|---|---|---|
| `asr_intron_summary` | `<family>_asr.intron.tsv` | Native NWKIT node probabilities plus observed `num_intron` |
| `asr_intron_model` | `<family>_asr.intron.model.tsv` | Fixed Q, root prior, likelihood, and model metadata |
| `asr_intron_tree` | `<family>_asr.intron.nhx` | All native ASR annotations, plus `Absent`/`Present` probability aliases for the plot |
| `asr_intron_plot` | `<family>_asr.intron.pdf` | NWKIT probability pies at all nodes; light gray is absent, dark is present |

The summary retains `branch_id`, `parent`, `node_class`, `name`,
`observed_state`, `is_imputed`, `map_state`, `map_probability`,
`p_intron_absent`, and `p_intron_present`. `num_intron` is attached to leaves
by gene ID and is missing for ancestors and unknown tips. The PDF is written
for constant-state and all-missing inputs too.

NWKIT branch IDs are level-order indices. GeneGalleon's `stat_branch` IDs are
clade ranks, so `orthogroup_statistics.py --asr_intron` verifies the NWKIT
node IDs, names, parents, and classes against the dated tree before translating
them. It also rejects invalid counts, counts attached to ancestors, inconsistent
imputation flags, and probabilities that contradict observed tip states. It retains
`intron_present`, `num_intron`, and `intron_is_imputed`, and
computes `delta_intron_present` for existing retrotransposition visualizations.
The default intron-number panel uses observed `num_intron` from the GFF
summary even when ancestral reconstruction is disabled. Zero remains zero;
tips without GFF counts remain missing. When ASR is available, its observed
counts must agree with GFF counts.
HGT intron-support counts use observed introns, not imputed tip probabilities.

The stage's `artifact_provenance/<family>.asr_intron.json` covers all four
outputs, the input tree and GFF table, the adapter, rates, root prior, branch-length
policy, and NWKIT identity. Existing stale-artifact policies still apply.
Computation failures preserve previous results. The adapter uses NWKIT's output
transaction, and the core stage uses GeneGalleon's bundle publisher to preserve
output inventory tracking. Both restore the previous set on handled publication
failures; neither guarantees multi-file atomicity across a process crash or power
loss. On restart, the provenance check rejects incomplete or mixed results.
`artifact_stale_policy="rebuild"` regenerates the complete set; the default `stop`
policy reports the mismatch before downstream analysis. Output aliases of adapter
inputs are rejected before computation.

## Migration from SCM

Replace `run_scm_intron` with `run_asr_intron` in saved overrides, including
`GG_GENE_EVOLUTION_RUN_ASR_INTRON` for environment overrides. The entrypoint
rejects retired SCM overrides with a migration instruction. The old
R stochastic-mapping helper, `--scm_intron` summary argument, and `scm_intron_*`
output names are retired. Old SCM files are not inputs to the new stage and
are not automatically deleted. Rerun ASR and downstream summaries with the
normal stale-artifact policy; use `artifact_stale_policy="rebuild"` when
intentionally regenerating dependent results.

This migration deliberately removes the old all-nodes 0/1 override and corrects
the misplaced `num_intron` values. Results can therefore differ from historical
SCM summaries, especially for missing observations and constant-state families.

## Validation

`workflow/tests/test_asr_intron_runtime.py` executes the real NWKIT commands
inside a GeneGalleon runtime. It compares probabilities with independent
two-state CTMC enumeration, checks counts and imputation, exercises output
publication failures, and runs the actual core stage through cache reuse,
missing-output recovery, rate-change invalidation, and forced process-group death
mid-publication followed by stale-result rejection and rebuild. It also passes
ASR results through the real branch-statistics, SQLite, and HGT scoring commands.
The runtime validation
manifest includes this test; run it in the container when changing this stage.

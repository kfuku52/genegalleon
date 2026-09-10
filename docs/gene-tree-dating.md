# Sequence dating and OU shifts

Gene-family dating (`run_tree_dating=1`) uses `nwkit radte --backend native`.
The native sequence engine remains the default; IQ-TREE is also selectable.
Native RADTE is an **experimental, exploratory estimator**. A nominal 95%
interval is not a demonstrated 95% guarantee across gene families. This also
applies when IQ-TREE supplies the sequence likelihood. Report ages as
"exploratory RADTE estimates with conditional, nominal intervals" and retain
the actual estimator and unavailable/calibration-limited interval diagnostics.
Dating is disabled by default (`run_tree_dating=0`). See the
[independent native validation record](native-dating-validation.md) for measured
coverage, numerical fixes and validation limits.
The species tree must be a dated, rooted, ultrametric chronogram with positive
branch lengths. An undated substitution tree is not a time calibration.

The default uses the trimmed CDS alignment with **GY94 + F3x4**, estimates one
shared κ and ω on an unclocked tree, and holds substitution parameters fixed for
subsequent dating. Protein mode defaults to LG. Both modes use four gamma rate
categories, native `inference=auto`, `likelihood=auto`, and **95% profile intervals**.
Auto may use validated quadratic marginal inference or exact conditional joint
MAP; the actual estimator and diagnostics are recorded. Profile intervals are
conditional on the species ages, reconciliation, and fitted substitution model;
they are not MCMC posterior intervals. Failed or calibration-limited intervals
are reported by NWKIT, without silently changing uncertainty methods. In auto
mode, a failed quadratic check during profile exploration triggers a recorded
refit of both point estimates and profile intervals with exact likelihood.

Species ages remain fixed. `radte_species_intervals_tsv` optionally supplies
external uncertainty for display only, with columns `node` (or stable
`species_event_id`), `lower`, `upper`, `level`, `kind`, `source`. These intervals
are neither hard calibration bounds nor propagated samples. Species not included
in the family's pruned tree must be removed from this table before use.

## Configuration

Edit `workflow/gg_gene_evolution_entrypoint.sh`, or use the corresponding
`GG_GENE_EVOLUTION_` environment overrides:

```sh
export GG_GENE_EVOLUTION_RUN_TREE_DATING=1
export GG_GENE_EVOLUTION_RADTE_SUBSTITUTION_MODEL=ecmk07
# Empty uses ECM's published equilibrium frequencies; f uses smoothed codon counts.
export GG_GENE_EVOLUTION_RADTE_CODON_FREQUENCIES=f
```

| Parameter | Default | Meaning |
|---|---|---|
| `radte_sequence_engine` | `native` | `native` or `iqtree`; both use NWKIT dating constraints and conditional intervals |
| `radte_iqtree_interface` | `auto` | Use an installed external library worker when available; `cli` selects standard exports, `library` requires the worker |
| `radte_iqtree_worker` | empty | Optional path to an externally built `nwkit-iqtree-worker` |
| `radte_iqtree_model` | empty | Complete IQ-TREE model such as `GY+F3X4+R4`; overrides the substitution model and gamma category setting |
| `radte_substitution_model` | `auto` | GY94 for CDS, LG for protein; `gy94`, `ecmk07`, `ecmrest`, and NWKIT NT/AA models can be selected |
| `radte_codon_frequencies` | empty | GY94: `f3x4`; ECM: published frequencies. Alternatives: `f`, `f1x4`, `f3x4`, `fq` |
| `radte_kappa`, `radte_omega` | empty | Estimate supported parameters; numbers fix them. ω is only for GY94, κ for HKY/GY94 |
| `radte_gamma_shape` | empty | Estimate gamma shape |
| `radte_gamma_categories` | `4` | Site-rate categories; `1` disables gamma variation |
| `radte_inference` | `auto` | `auto`, `marginal`, `joint-map` |
| `radte_likelihood` | `auto` | `auto`, `exact`, `quadratic` |
| `radte_uncertainty` | `profile` | `profile`, `laplace`, `bootstrap`, `none` |
| `radte_interval_level` | `0.95` | Nominal conditional interval level; general coverage is unvalidated |
| `radte_rate_sd` | empty | Estimate log-rate SD; a number fixes it |
| `radte_max_age` | `1000` | Upper age limit for above-root duplications, in the species tree's time units |
| `radte_maxiter`, `radte_seed` | `1000`, `1` | Optimizer iterations and random seed |
| `radte_species_intervals_tsv` | empty | Workspace-relative or absolute external interval table |

Codon models currently require standard genetic code 1 for every represented
species, including per-species overrides. Complete in-frame codon gaps and IUPAC
ambiguities are accepted; partial gaps, frames not divisible by three, and stop
codons are errors. GY94 permits only single-base changes. ECMrest uses published
single-base exchangeabilities; ECMK07 also permits multi-base changes. ECM has no
free κ/ω multiplier. The native engine records codon rates in expected nucleotide changes per
codon site per time unit. Its `RADTE.md` documents normalization and frequency
pseudocounts.

## Selecting IQ-TREE

```sh
export GG_GENE_EVOLUTION_RADTE_SEQUENCE_ENGINE=iqtree
# Optional: FreeRate instead of the default four-category gamma model.
export GG_GENE_EVOLUTION_RADTE_IQTREE_MODEL='GY+F3X4+R4'
```

Leave `radte_iqtree_model` empty to translate the existing model settings into
IQ-TREE syntax. When supplying a complete model, leave `radte_kappa`,
`radte_omega`, `radte_codon_frequencies`, and `radte_gamma_shape` empty; put those
controls in the model string instead. `radte_gamma_categories` is then unused.
Switch `radte_sequence_engine` back to `native` and clear `radte_iqtree_model` to
use the internal likelihood implementation.

IQ-TREE receives the reconciled topology and the same alignment used for dating.
NWKIT retains the species-age constraints, clock optimization and CI computation.
IQ-TREE model parameters are fitted once and frozen during dating. NWKIT uses
unmodified IQ-TREE 3 through its standard IQ2MC CLI (`iqtree3`) for the initial
model fit. For repeated branch evaluations, `auto` uses an installed external
library worker when available, otherwise the CLI. `library` explicitly requires
the worker. No analysis downloads or builds software. NWKIT computes full observed
curvature from branch scores, caches repeated evaluations and retains shared
marginal-profile matrices. Threads follow `GG_TASK_CPUS`. No MCMCTree process is
launched and no upstream source extension is needed.

The adapter supports reversible NT models JC/HKY/GTR/F81, AA models
Poisson/LG/WAG/JTT, and codon models GY/MG/ECMK07/ECMrest, with compatible frequency,
`+I`, `+Gk`, and `+Rk` options. It requires at least three sequences and standard
code 1 for codons. ModelFinder, partitions and mixture models are not yet exposed.
Failed or unsupported IQ-TREE evaluations stop the stage and preserve the previous
published results; there is no implicit change of engine.

IQ-TREE uses its own frequency estimator and branch-length units. For ECMK07,
its codon-event rate differs from the native nucleotide-change rate. Compare
absolute rates only after converting units. The manifest records the engine,
requested/frozen model, fitted parameters and IQ-TREE executable identity. The
normal downstream filenames and fixed-species-age/display-only-CI policy apply
to both engines. The executable must be IQ-TREE 3 or later with IQ2MC support.

New standard GeneGalleon container builds include the official IQ-TREE 3 CLI
and external library worker. With `radte_sequence_engine=iqtree`, the default
`radte_iqtree_interface=auto` selects it automatically. Check availability with
`python -m nwkit.iqtree_library check --interface library` inside the container.

For development against a local NWKIT checkout, an overlay is also available:

```sh
BASE_IMAGE=local/genegalleon:dev \
IMAGE=local/genegalleon:iqtree3-dev \
bash container/build_iqtree3_overlay.sh /path/to/nwkit
```

The wrapper clones https://github.com/iqtree/iqtree3 recursively without source
changes. It builds the CLI and library from the same snapshot and installs the
library worker as a separate runtime executable, outside NWKIT’s Python package. The default branch is declared in `container/source_branches.env`;
`IQTREE_REPO_SHA` may select a one-off reproduction. Resolved revisions, source
checksums and binary identity are build metadata under `/opt/pg/logs`.
The worker identity is recorded in `iqtree3_library_worker.json`. For manual
setup outside this overlay, see the
[NWKIT library guide](https://github.com/kfuku52/nwkit/blob/master/IQTREE_LIBRARY.md).
Use `GG_CONTAINER_RUNTIME=docker` and
`GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:iqtree3-dev` for this image.
Docker validation does not establish SIF compatibility. Historical custom-session
benchmark reports do not describe this official CLI path.

## Reconciliation, results and reruns

GeneRax NHX annotations or Notung parsable reconciliation identify species events.
GeneRax's species-tree node names are mapped by descendant clades onto the dated
species tree, requiring identical rooted topology. GeneRax substitution branches
do not supply species ages. Transfer events fail explicitly.

Existing downstream paths remain `dated_tree/<family>_dated.nwk` and
`dated_tree_log/<family>_dated.log.txt`. The latter now contains the NWKIT run
manifest. Full native tables, likelihood summary, manifest, and PDF are stored
under `dated_tree_native/<family>_radte.*`. The orthogroup summary reads the actual
estimator, sequence model, uncertainty status, and interval level from this
manifest. Native results also retain `dating_interpretation` (exploratory,
conditional, with general coverage unestablished) and `dating_diagnostics` in
the orthogroup statistics. An unavailable interval remains unavailable even
when the point estimate completed. Historical text logs remain readable for
historical results.

The cache includes the alignment, reconciliation, species tree, model settings,
external interval table, genetic-code overrides, NWKIT identity, and the selected
sequence engine, IQ-TREE model/executable identity, selected interface and
external library/worker identity. Old RADTE
results do not satisfy this provenance. To recompute stale results, use the
workflow's explicit `artifact_stale_policy=rerun` setting. Output publication is
transactional: an inference or rendering failure preserves the previous bundle
and does not delete downstream analyses.

OU shift detection already uses **kfl1ou**, through
`workflow/support/detect_OU_shift_kfl1ou.r`. The stage and provenance identify that
engine and its version, while `run_l1ou`, the existing options, and `l1ou_*` output
names continue to connect to downstream readers. Its dated-tree input participates
in cache validation, so changed dates invalidate an old fit. The standalone legacy
RADTE R program is no longer installed in new containers. PAML remains available
for the separate species-tree workflow.

### IQ-TREE numerical limits

IQ-TREE can produce nonfinite scores near zero-length codon branches. NWKIT
rejects nonfinite exports and preserves the previous published bundle. It does
not patch IQ-TREE, floor branches or change the fitted model to conceal numerical
failures. Upstream numerical defects require an upstream fix. Validation of
ordinary Gamma and FreeRate inputs does not cover every numerical boundary.

The unmodified official IQ-TREE 3.1.4 validation also compares sibling codon
branches near `8e-11` and `6e-11` with independent SciPy matrix-exponential
pruning. This case passes at the existing IQ2MC text-export tolerances
(NLL absolute `2e-5`; score relative `1e-5`, absolute `3e-4`). It requires no
upstream patch and does not establish accuracy at every numerical boundary.

See the [official IQ-TREE 3 validation record](iqtree3-validation.md) for the
container build and completed checks.

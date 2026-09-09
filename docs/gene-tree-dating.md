# Sequence dating and OU shifts

Gene-family dating (`run_tree_dating=1`) uses `nwkit radte --backend native`.
The native sequence engine remains the default; IQ-TREE is also selectable.
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
| `radte_iqtree_mode` | `persistent` | Keep the frozen model/alignment loaded; `subprocess` explicitly selects the original CLI route |
| `radte_iqtree_model` | empty | Complete IQ-TREE model such as `GY+F3X4+R4`; overrides the substitution model and gamma category setting |
| `radte_substitution_model` | `auto` | GY94 for CDS, LG for protein; `gy94`, `ecmk07`, `ecmrest`, and NWKIT NT/AA models can be selected |
| `radte_codon_frequencies` | empty | GY94: `f3x4`; ECM: published frequencies. Alternatives: `f`, `f1x4`, `f3x4`, `fq` |
| `radte_kappa`, `radte_omega` | empty | Estimate supported parameters; numbers fix them. ω is only for GY94, κ for HKY/GY94 |
| `radte_gamma_shape` | empty | Estimate gamma shape |
| `radte_gamma_categories` | `4` | Site-rate categories; `1` disables gamma variation |
| `radte_inference` | `auto` | `auto`, `marginal`, `joint-map` |
| `radte_likelihood` | `auto` | `auto`, `exact`, `quadratic` |
| `radte_uncertainty` | `profile` | `profile`, `laplace`, `bootstrap`, `none` |
| `radte_interval_level` | `0.95` | Conditional interval coverage |
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
IQ-TREE model parameters are fitted once and frozen during dating. After this
separate prefit, one persistent worker keeps the alignment, topology and model
loaded and accepts branch-length vectors over pipes. It returns likelihoods and
scores directly without repeatedly writing a full Hessian. Full curvature is
computed from those scores, and the same exact-validation/profile-refit rules
apply. Worker threads follow `GG_TASK_CPUS`. No MCMCTree process is launched.
Use `GG_GENE_EVOLUTION_RADTE_IQTREE_MODE=subprocess` to explicitly select the
original per-evaluation CLI route. There is no automatic fallback when the
persistent protocol is absent or fails.

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
to both engines. Persistent mode requires the updated NWKIT package and the local
IQ-TREE source extension providing `--likelihood-session`, as well as IQ2MC
support for the prefit. Stock released IQ-TREE binaries are not assumed to
provide the session extension.

Build a local GeneGalleon overlay from the two updated source checkouts:

```sh
BASE_IMAGE=local/genegalleon:dev \
IMAGE=local/genegalleon:iqtree-session-dev \
bash container/build_iqtree_session_overlay.sh /path/to/iqtree /path/to/nwkit
```

Initialize IQ-TREE's submodules before building. The overlay keeps the IQ-TREE
source extension in its owning repository, installs the matching NWKIT package,
and records the source-content checksums and compiled binary identity under
`/opt/pg/logs`. It does not embed upstream version or commit defaults. Use this
image with `GG_CONTAINER_RUNTIME=docker` and
`GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:iqtree-session-dev`. This is Docker
validation/build support; it does not establish SIF compatibility.

See the [persistent-session validation measurements](iqtree-session-benchmark.md)
for controlled timing comparisons, numerical checks and their limits.

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
manifest. Historical text logs remain readable for historical results.

The cache includes the alignment, reconciliation, species tree, model settings,
external interval table, genetic-code overrides, NWKIT identity, and the selected
sequence engine and IQ-TREE model/executable identity. Old RADTE
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

The locally extended persistent worker checks for spectral cancellation and
recomputes ill-conditioned requests inside IQ-TREE using original-state pruning
and a stable matrix exponential of the same fitted generator. A regression with
two sibling codon branches near `8e-11` and `6e-11` checks likelihood and scores
against an independent SciPy matrix exponential. Branch lengths are not floored,
and the sequence engine and fitted model remain unchanged. Such requests take
more computation than ordinary spectral evaluations.

The subprocess/IQ2MC route can still return nonfinite branch scores for this
boundary case (observed in IQ-TREE 3.1.3); a finite reported likelihood alone does
not establish accuracy. NWKIT rejects nonfinite results and preserves the
previous published bundle. The default Gamma and FreeRate end-to-end tests,
and the numerical regression, do not establish support for every alignment and
boundary case.

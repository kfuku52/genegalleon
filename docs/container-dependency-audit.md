> Update: GeneGalleon now uses NWKIT for OU inference; kfl1ou has been removed
> from the workflow and container dependency set. The earlier audit observations
> below describe the environment before that migration. Current OU integration
> coverage is in `workflow/tests/test_native_ou_shift.py`.

# Container dependencies after the NWKIT migrations

The audit covers the standard amd64/arm64 Conda manifests, source artifacts,
CRAN installers, required-command checks, workflow calls and relevant R package
dependency declarations. A filename mentioning an old program is not by itself
evidence that its executable is still used.

## Changes

| Dependency | Decision | Evidence |
| --- | --- | --- |
| Explicit Conda `iqtree=3.*` request | Remove the direct request; build official IQ-TREE 3 CLI and its library worker together | NWKIT requires the matching library adapter. Conda's CLI package does not provide it. Both `iqtree` and `iqtree3` resolve to `/opt/pg/iqtree3/iqtree3`. |
| R `magic` | Remove from the required R manifest | No current GeneGalleon call, no current kfl1ou runtime import and no installed Conda package depends on it in the audited environment. |
| R `grplasso` | Remove from the required R manifest | Current kfl1ou lists it under `Suggests`; only its optional `package` backend requires it. GeneGalleon's adapter uses the default internal backend and exposes no package-backend selector. |
| IQ-TREE library worker | Add as a standard, separately installed executable | Build from the official unmodified `BUILD_LIB=ON` library and the adapter extracted from the exact NWKIT wheel installed in the image. |

OrthoFinder still declares a Conda dependency on IQ-TREE, so removing the direct
request does **not** remove that transitive package from the solved environment.
Do not force-remove it or bypass Conda's dependency checks. Runtime entrypoints
use the official source artifact, whose revision is recorded separately.

## Dependencies retained

| Dependency | Active consumer or reason |
| --- | --- |
| IQ2MC-compatible PAML/MCMCtree source build | `gg_genome_evolution_core.sh` still prepares IQ2MC input and runs the species-tree MCMCtree analysis. NWKIT now handles gene-tree dating, not this separate species-tree fit. |
| Conda PAML | `gg_gene_evolution_core.sh` still runs `codeml` for two-ratio models. The source artifact installs MCMCtree only. |
| kfl1ou | Removed after migration to NWKIT native OU inference. |
| `ape`, `Rcpp`, `nlme` | kfl1ou imports ape/Rcpp; ape itself imports nlme/Rcpp. Tree and plotting helpers also use ape. |
| `phytools`, `phangorn`, `igraph` | rkftools imports phytools; phytools imports phangorn; phangorn imports igraph. Annotation plots also load phytools. |
| `missMDA` | `multispecies_transcriptome_summary.r` still calls `estim_ncpPCA` and `imputePCA`. |
| ggplot2/ggtree, cowplot, svglite and related plotting packages | The packaged `genegalleon.treevis` API and other R plotting helpers remain active despite migration of some individual figures to NWKIT. |
| NOTUNG, GeneRax and GRAMPA | Reconciliation and polyploidy-related workflow stages remain active. NWKIT rooting and tree conversion do not replace those analyses. |

## Already removed by earlier migrations

- Standalone R RADTE: gene-tree dating now invokes `nwkit radte`.
- Rphylopars: neither the current CRAN installer nor its required runtime checks
  request it; regression moved to NWKIT.
- Standalone MAD and legacy tree conversion/plotting dependencies are not
  explicit standard-container installs. Retained historical output filenames
  such as `mcmctree_95CI.nhx` do not imply a separate plotting package is needed.

## Verification

The standard image build checks the worker's protocol/version and compares
three resident likelihood evaluations, gradients and diagonal second
derivatives against the installed official CLI. It removes the worker's input
alignment after initialization to verify reuse of loaded data. A missing worker
or numerical mismatch fails the build.

Dependency-removal validation also exercises a multivariate kfl1ou fit and
loads `genegalleon.treevis` without `magic` or `grplasso`. The intended deployment
paths remain Docker and native Apptainer; Docker execution alone does not
establish SIF compatibility. Build results and platform coverage should be
reported with the actual validation run.

Validated on 2026-09-10 using the standard Linux ARM64 Docker build, with
IQ-TREE 3.1.4 and NWKIT 0.43.15:

- Complete standard image build and all 56 required runtime checks passed.
- The resident-worker/CLI numerical comparison passed; source archive,
  submodule contents and NWKIT's MIT metadata were verified in the image.
- Both Gamma and FreeRate gene-dating integration cases passed (102.14 seconds).
- A multivariate kfl1ou fit and treevis imports passed with both removed R
  packages absent from the actual standard image.
- 366 container, CI, source-policy and shell tests passed. Ruff, ShellCheck
  and workflow syntax checks also passed.

Native amd64 and SIF execution were not performed in this local validation.

## Alignment statistics migration

AMAS is removed from both container environments and command inventories. Gene-evolution
alignment summaries use `cdskit stats --mode alignment` (CDSKIT 0.31.0 or newer),
with `--seq_type dna` or `aa` according to the input mode. The source build
continues to follow the moving branch in `container/source_branches.env`.

Existing `run_amas_original` / `run_amas_cleaned` switches, `amas_*` output
paths, provenance step names, and summary columns are retained for existing
workspaces and ZIP archives. Provenance records the new statistics engine so
previous manifests invalidate on the next enabled stage run. Input FASTA is
still decompressed with seqkit. Protein summaries include `GC_content=NA`,
allowing the shared summary readers to handle both sequence types.

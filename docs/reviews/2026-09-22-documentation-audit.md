# Documentation/implementation audit — 2026-09-22

## Scope and baseline

Baseline: `main`, `43e65cef7de7ff0485f88c3d850be34bef62e9fa`, version
0.8.1. The worktree was clean; `git fetch origin` confirmed the same remote
baseline. Applied root `AGENTS.md`, development/runtime validation guidance,
and the `validate-change` and `prepare-github-push` skills. No nested tracked
`AGENTS.md` was present. Publication advances only `VERSION` to 0.8.2 alongside
these documentation changes; analytical behavior and dependencies are unchanged.

Prioritized installation, the bundled query2family entrypoint, configuration
forwarding, input identifiers, family outputs, progress/status commands and
reruns. Read their argument definitions, processing code and relevant tests,
not only help. Historical reviews, benchmarks and release descriptions were
not rewritten. This is a bounded audit, not certification of every stage.

## A — corrected documentation

| Location and previous statement/omission | Implementation/test evidence | Correction |
| --- | --- | --- |
| [README](../../README.md#quick-start), [recipe 1](../common-workflow-recipes.md#1-run-the-bundled-query2family-example): missing checkout prerequisites; bundled analysis called a smoke test without task/resource/download caveats | `gg_gene_evolution_entrypoint.sh` enables RPS-BLAST/localization; core selects one sorted visible query file; `gg_util/02_container_scheduler.sh` defaults local task/CPU to 1; tracked `query_gene/AHA` exists | Separate runtime preparation from execution; identify task 1, first-use downloads, and the distinction from `dev check smoke`; link rerun policy |
| [Container guide](../container-build-and-runtime.md#host-prerequisites): host/architecture and writable-cache prerequisites were scattered | `container/gg_container_build_impl.sh`, `container/env/base.arm64.required.txt`, runtime dispatcher and reference helpers | State Linux SIF versus Docker Linux/macOS paths, arm64 assembly limitations, network/cache prerequisites |
| [Container guide](../container-build-and-runtime.md): IQ-TREE described as conda-pinned to `3.*` | `container/Dockerfile` IQ-TREE build/worker stages, `source_branches.env`, `scripts/install_source_artifacts.sh`; conda manifests have no IQ-TREE entry | Describe the source-built official CLI and worker; container README introduction no longer implies a required overlay |
| [Container README](../../container/README.md#important-caveats): references must be populated manually under `/usr/local/db` | `gg_util/09_sequence_databases.sh` and `10_reference_databases.sh` reuse system references or download/build workspace caches | Explain automatic cache preparation and link troubleshooting |
| [Container guide](../container-build-and-runtime.md): large release SIF kept only for 90 days and historical SIF should be rebuilt | `.github/workflows/release-sif.yml`: durable OCI publication, one-day handoff, `.oci.txt`/qualification/checksum assets, 2 GiB upload threshold | Replace stale storage description and refer to container README's existing canonical retention policy |
| Container pull examples resemble real immutable/release tags | The illustrated owner/hash/version are not established published artifacts | Label them as templates requiring an existing published tag; no registry availability claim |
| [Alignment migration](../alignment-statistics.md#migrating-an-existing-workspace): example writes `workspace/output/orthogroup_summary.tsv` | Progress core passes `--out orthogroup_summary.tsv` after workspace bootstrap; orthogroup summary tests preserve selected input | Align migration example with `workspace/orthogroup_summary.tsv`, avoiding an aggregate the normal progress wrapper does not update |
| [Progress guide](../gene-family-outputs-and-progress-monitoring.md#gene-family-progress-summaries): implied unconditional table generation and completion by visible subdirectory, omitted mutation | `gg_progress_summary_core.sh` checks prerequisites, drains queues and cleans retained temporary data; both summary readers use logical live/ZIP stores and filename presence | Explain conditional outputs, overwrite of aggregates, augmented gene-count sidecar, ZIP/shared-folder handling, presence versus verified completion, and read-only inspection alternative |

Configuration precedence was cross-checked against
`gg_apply_registered_env_overrides`, `gg_apply_env_override_to_config_var`,
entrypoint activation and the forwarding registry. The documented registered
scoped override > editable block/shared value > core fallback order agrees;
no change was needed. Generated metadata was checked rather than hand-edited.
`Parsimony_informative_sites_clean` is still the actual aggregate column name;
its spelling was deliberately retained.

## B: protein input validation

Follow-up in 0.8.3: fixed the first-header-only restriction to honor the
existing external-protein-ID contract. All records may use unprefixed IDs;
empty-file, duplicate-ID and prohibited-character rejection remain. Regression
coverage includes both mixed-prefix orders and preserved external IDs through
the genome-evolution protein input path. The evidence below describes the
0.8.2 baseline, not the corrected behavior.

Follow-up validation used the freshness-checked `local/genegalleon:dev` Docker
runtime: 9 focused helper tests, 2 protein-input preparation cases, and
`bash ./dev check static` (291 tests) passed without skips. Six temporary-input
cases also passed using real seqkit, including rejection of empty, duplicate-ID
and prohibited-character inputs. Host `dev lint` was blocked by Bash 3.2;
container Bash syntax checks passed for all 74 tracked shell entrypoints,
while host Ruff and `dev config-check` passed separately. SIF was not tested.

[Input conventions](../input-conventions.md#workspaceinputspecies_protein)
recommend, but do not require, species-prefixed protein IDs. In contrast,
`workflow/support/gg_util/06_workspace_validation.sh`, function
`check_species_protein_dir`, checks the **first** header against the filename's
species prefix. `prepare_species_protein_tmp` in the genome-evolution core calls
this validator for provided proteins. The remaining headers are checked for
uniqueness/prohibited characters but not for that species prefix.

Executed in the verified Docker runtime with temporary `Species_name.fa`:

| FASTA headers (sequence `MPEPTIDE` for each) | Validator exit/result |
| --- | --- |
| `gene1` | 1: inconsistent with species parsed from filename |
| `Species_name_gene1`, then `gene2` | 0: reported all per-species protein files valid |

Minimal reproduction, from repository root (no workflow or persistent workspace):

```bash
bash workflow/tests/run_in_runtime.sh bash -s <<'SH'
set -u
audit_dir=$(mktemp -d)
trap 'rm -rf -- "$audit_dir"' EXIT
source workflow/support/gg_util.sh
export GG_TASK_CPUS=1
printf '>gene1\nMPEPTIDE\n' > "$audit_dir/Species_name.fa"
if (check_species_protein_dir "$audit_dir"); then
  printf 'first-header exit: 0\n'
else
  printf 'first-header exit: %s\n' "$?"
fi
printf '>Species_name_gene1\nMPEPTIDE\n>gene2\nMPEPTIDE\n' > "$audit_dir/Species_name.fa"
if (check_species_protein_dir "$audit_dir"); then
  printf 'later-header exit: 0\n'
else
  printf 'later-header exit: %s\n' "$?"
fi
SH
```

The executed reproduction used Python `TemporaryDirectory` plus the same shell
function, recording and asserting both exit codes. The shell form above is also
checked separately. Existing protein-mode fixtures use prefixed headers and do
not resolve whether unprefixed IDs should be supported. This requires a contract
decision: either support external IDs consistently, or explicitly require and
validate prefixes for every record. No implementation or recommendation was
silently changed; the input guide now flags the discrepancy.

## C: scientific interpretation held for review

The progress guide says rows with `Parsimony_informative_sites_clean == 0`
“cannot produce normal IQ-TREE-based downstream outputs.” The core's IQ-TREE
stage does not gate execution on that column; the value is an alignment
statistic. No tested definition of “normal” or general scientific impossibility
was established here. Whether to replace this with an identifiability/quality
caveat needs model-specific evidence and maintainer review. The statement is
marked as pending review in the guide and left unresolved, rather than inferring scientific validity from a file flag or
successful tool exit. No scientific threshold/model was changed.

## Executed validation

Host: macOS arm64. Docker: `local/genegalleon:dev`, image ID prefix
`509ccb364688`, Python 3.12.14. The normal runtime wrapper confirmed current
container inputs and the daily owned-upstream snapshot; freshness was not
bypassed. Docker evidence does not establish SIF compatibility.

- `bash workflow/tests/run_in_runtime.sh python -c 'import sys; print(sys.version)'`:
  exit 0, runtime/freshness available.
- `bash ./dev check fast workflow/tests/test_query2family_output_summary.py workflow/tests/test_orthogroup_output_summary.py workflow/tests/test_docker_runtime_shim.py -x`:
  **21 passed**, no skips. Includes live/ZIP outputs, alignment statistics,
  stable task rows and Docker dispatch. Docker dispatch tests use command stubs;
  the real progress wrapper was additionally executed below.
- `bash ./dev config-check`: exit 0, **8 entrypoints / 23 common parameters**.
- `bash ./dev config-schema json > /tmp/genegalleon-doc-audit-schema.json`:
  exit 0. These two are the documented commands themselves, executed on host.
- Through `run_in_runtime.sh`, `bash ./dev --help`,
  `python workflow/support/query2family_output_summary.py --help`,
  `bash workflow/gg_gene_family_archive.sh status --help`, and
  `python workflow/migrations/migrate_alignment_statistics.py --help`: all exit 0.
- Scaled alternative: invoke `query2family_output_summary.py` with required
  `--dir_query2family`, `--dir_query_gene`, `--out` pointing into one container
  `TemporaryDirectory`. Two query files `AHA`, `YABBY`, one empty
  `tree_plot/AHA_tree_plot.pdf`: exit 0, TSV columns
  `query`, `GG_ARRAY_TASK_ID`, `tree_plot`; rows `AHA,1,1` and `YABBY,2,0`.
  `gg_gene_family_archive.sh status --root ... --mode query2family --query-dir ...`
  on the same input: exit 0, `present=1 missing=1 expected=2`.
- Scaled wrapper execution: `gg_workspace_dir=<host TemporaryDirectory>
  GG_CONTAINER_RUNTIME=docker GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:dev
  bash workflow/gg_progress_summary_entrypoint.sh`: exit 0. One synthetic family
  produced `query2family_summary.tsv` at the workspace root with the expected
  row/columns. Queue reported idle; runtime/version/resource records were also
  confined to that disposable workspace. This does **not** validate a real PDF
  or execute the README's gene-evolution analysis.
- `bash ./dev bump patch`: 0.8.1 → 0.8.2, only `VERSION` changed.
- `bash ./dev check fast workflow/tests/test_development_tooling.py -k 'generated_config_schema or version_helper' -x`:
  **2 passed**, no skips; generated-config and version-helper verification.
- Local Markdown target scan of README, container README and 58 top-level docs:
  **60 files, no missing local file targets** before edits; **61 files, zero
  missing targets** after edits including this report.
  This lightweight scan does not verify external URLs or every heading slug.
  No dedicated documentation build/link-check command was found in `dev`/CI.
- `git diff --check`: passed. Final diff reviewed; only Markdown and the required
  `VERSION` bump are included. Host lint/full scientific suites were not run
  because the documented documentation-only check lane does not require them.

## Unverified boundaries

No clean image build, registry pull, model/database download, complete bundled
scientific run, migration of existing outputs, or scheduler submission was
performed: these require substantial downloads/computation or mutate persistent
state. SIF and HPC scheduler execution were unavailable on this macOS host.
External website links and availability of illustrative registry tags were not
verified. Deep statistical semantics (dating, OU/RSC, CSUBST/FDR, HGT), all
provider-specific input schemas, every helper CLI/docstring, and Windows support
were not comprehensively audited. Existing scientific and historical validation
claims in other documents remain their original evidence, not results of this
audit. No new test framework or production behavior change was introduced.

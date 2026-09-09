# Workflow inspection API and compatible migration

GeneGalleon provides an optional JSON interface for monitoring, inspecting a
declared restart, and verifying artifacts. It works without kfauto, a scheduler,
an approval database, or a background service. The API never submits a job,
repairs an artifact, drains an archive queue, or authorizes execution.

Existing entrypoints, configuration blocks, output formats, and progress-summary
behavior remain supported. In particular, progress-summary still drains archive
requests. Inspection is a separate operation.

## Query capabilities

Run from the repository root, using the GeneGalleon runtime for workflow data:

```bash
python -B workflow/support/workflow_api.py capabilities
```

All query responses use `schema=genegalleon-api-v1`, include `observed_at_ns` and
`read_only=true`, and write one JSON object to stdout. Successful inspection
returns process status 0 even when it reports an artifact problem. Query errors
return status 2 with `error_code=query_unavailable` (or `cursor_invalid` for an
invalid paged cursor); they are not job failures.
CLI usage errors return argparse's usual status 2 on stderr. Consumers must
check the schema and command-specific coverage, and reject unsupported schemas.
Additional object fields may be added without changing the schema identifier.

## Record new attempts, optionally

```bash
GG_OBSERVABILITY=1 bash workflow/gg_gene_evolution_entrypoint.sh
```

With an `exec` container adapter, this adds per-attempt JSON under
`workspace/output/observations/<attempt-id>/`. The default is off. Enabling it
does not convert existing projects or create records for historical executions.
Recording failure emits a warning and preserves the analytical command's result.

The observer preserves stdin/stdout/stderr and the exit status, forwards
termination to the child process group, and keeps inherited lock descriptors.
The gene-evolution bridge recognizes its existing status 8 (outputs already
present) as accepted, while still returning 8 to the existing entrypoint.
Other failures without structured evidence remain `unknown`; status 137 alone
does not prove an out-of-memory failure.

`run.json` separates execution acceptance from scientific completion. A started
record without a final record means that termination is unobserved; it is not
proof that a job is still running. SIGKILL, node loss, and container startup
failure can leave incomplete or absent observations. Consult the scheduler for
liveness and allocation errors.

During an observed attempt, provenance helpers publish `contract-*.json` receipts
with the actual declared inputs, outputs, parameters, policy, operation, and
manifest digest. Newly recorded artifact manifests additionally retain the
attempt ID as a diagnostic; this does not affect artifact freshness. Reusing an
older manifest records a receipt without rewriting its historical diagnostics.
Receipt paths are made absolute at recording time so inspection does not depend
on the original process's working directory. A later failed contract check
replaces an earlier successful receipt for the same family and step.

Known provenance boundaries also publish `error-*.json`, containing `step`,
`error_code`, `affected_input` (nullable), `retryability`, and diagnostic detail.
Stable codes include `inputs_changed`, `outputs_changed`,
`optional_outputs_changed`, `parameters_changed`, `identity_changed`,
`schema_changed`, `artifact_stale`, and `provenance_error`. Arbitrary dependency
logs are not parsed into confident diagnoses. Error evidence can describe an
earlier handled problem, so consumers must not automatically treat the latest
error record as the cause of the process exit. A retryability hint is never
permission to retry.

## Read status and changes

```bash
python -B workflow/support/workflow_api.py status \
  --directory workspace/output/observations

# NEXT_CURSOR is the opaque token from the last complete response.
python -B workflow/support/workflow_api.py status \
  --directory workspace/output/observations --since "$NEXT_CURSOR"
```

For large cursors, save the token as a JSON string and use `--since-file` instead
of a command-line argument. The reader returns changed records, removed attempt
IDs, errors, and `next_cursor`. An incomplete response has no next cursor; keep
the previous cursor and retry. A missing directory is unavailable, not an empty
successful observation. An existing empty directory covers no recorded attempts.
An attempt directory with a missing or malformed `run.json` makes the response
incomplete; it is not reported as a removed attempt. JSON observation files must
be regular files with unambiguous keys and finite values. Each read checks that
the file was not replaced or modified while being read.

Cursors bind to the observation directory and use content fingerprints, so
late publication and preserved timestamps do not hide updates. Changes are
also sensitive to receipt declarations: each compact contract summary includes
`receipt_sha256`, covering the complete receipt, including its argument list.
Responses are deltas; the implementation still inventories and hashes observation files. It
does not scan analytical output trees, and does not claim O(changes) I/O. The
inventory defaults to a 100,000-record bound. Snapshots are per-file observations,
not a globally atomic view. Old runs without observations are outside coverage.

### Bounded pages for large inventories

Consumers of large inventories should negotiate `status_pages=attempt-pages-v1`
and use `status --directory PATH --page-size 512`. Continue with the opaque
`--page-cursor` until `snapshot_complete=true`. Page size is 1–512; records are
bounded to 8 MiB, responses to 12 MiB, and cursors to 4 KiB. A page can contain
fewer records than requested to respect the byte limit. The legacy `--since`
interface above remains available, but its full inventory cursor is unsuitable
for bounded large-inventory consumers.

Each page reports `root_identity`, `inventory_sha256`, `inventory_count`,
`page_after`, `page_last`, and ordered `[attempt_id, sha256]` pairs in
`fingerprints`. Optional `--known-records-file` accepts an object with
`schema=genegalleon-status-known-v1`, the same `root_identity`, and a `records`
map of at most 512 attempt IDs to fingerprints. Unchanged hinted records are
omitted from `records`; fingerprints still cover every visited attempt. Consumers
must resolve omitted records from verified local copies and detect removals by
comparing the completed inventory with the previous generation.

Compact continuation cursors bind page size, position, physical directory
identity (resolved path/device/inode), and the directory-membership digest.
Membership changes or a stale/wrong-root cursor produce `cursor_invalid`. Discard
the incomplete generation and retry once from page one without hints. Ordinary
I/O errors remain `query_unavailable` and must not trigger a cache reset fallback.
Publish consumer state only after validating the complete inventory and all
record digests. Records remain per-file observations; pages do not create a
globally atomic snapshot. Each page inventories directory names, so pagination
bounds payload and record memory, not total directory-enumeration I/O.

## Inspect a proposed restart without writing

```bash
python -B workflow/support/workflow_api.py preflight \
  --attempt workspace/output/observations/ATTEMPT_ID
```

This rechecks the contracts actually encountered by that attempt. Steps not yet
reached, scheduler admission, external services, and undeclared workflow inputs
are outside its coverage. It does not infer the complete workflow dependency
graph or promise that a subsequent execution will finish.

For new settings, an old project without observations, or additional steps,
provide an explicit plan. Each contract uses the existing
`artifact_provenance.py needs-run` arguments, omitting the command and `--dry-run`:

```json
{
  "schema": "genegalleon-preflight-plan-v1",
  "contracts": [[
    "--workspace-root", "/workspace",
    "--logical-root", "/workspace/output/orthogroup",
    "--manifest", "/workspace/output/orthogroup/artifact_provenance/OG0001.summary_statistics.json",
    "--family-id", "OG0001",
    "--step", "summary_statistics",
    "--input", "tree=/workspace/output/orthogroup/rooted_tree/OG0001_root.nwk",
    "--output", "table=/workspace/output/orthogroup/stat_branch/OG0001_stat.branch.tsv",
    "--parameter", "mode=a",
    "--stale-policy", "stop"
  ]]
}
```

```bash
python -B workflow/support/workflow_api.py preflight --plan restart-plan.json
```

Paths are resolved in the inspecting process's environment. Inspect observed
container contracts using the same `/workspace` mapping. If the workspace has
moved, pass `--workspace-root /new/workspace`: only declared paths beneath the
old workspace are relocated. External absolute paths and parameter values are
preserved. The resulting binding digest covers the relocated declarations.
Pass `--stale-policy stop|rebuild|reuse` only when explicitly previewing a new
policy; omission preserves each recorded policy. This override does not infer
other changed entrypoint settings or authorize a restart.

Preflight calls the same provenance decision engine used by execution, with
dry-run enforced. It disables the persistent digest cache and legacy adoption
writes. Dry-run can read logical raw/ZIP artifacts without materializing them,
including validating archive-backed FASTA through the same sequence validator.
Recovery recipes are checked but never executed.

| Result | Meaning |
| --- | --- |
| `verified_current` | Declared contract matches recorded inputs, outputs, and parameters |
| `legacy_reusable` | Existing legacy policy permits reuse; historical generation is unproven |
| `policy_reuse_or_adoption` | Explicit stale reuse or compatible adoption would occur |
| `needs_run` | Runtime would generate or restore this declared artifact |
| `blocked` | Runtime's stale policy stops this contract |
| `unavailable` | The contract cannot be inspected reliably |

Missing historical evidence does not itself require a rebuild. A manifest
backfilled by an earlier version remains an adoption-time baseline, not proof of
the original generating inputs. The preview includes runtime exit codes,
diagnostics, and a binding digest; runtime validation must still run immediately
before execution. Preflight never changes the configured stale policy.

## Verify declared completion, including ZIP storage

```bash
python -B workflow/support/workflow_api.py verify \
  --root workspace/output/orthogroup --workspace-root workspace \
  --family-id OG0001 --require-step summary_statistics --require-step tree_plot

# Additionally require receipts from this exact accepted execution:
python -B workflow/support/workflow_api.py verify \
  --root workspace/output/orthogroup --workspace-root workspace \
  --family-id OG0001 --require-step summary_statistics \
  --attempt workspace/output/observations/ATTEMPT_ID
```

The result reports the family generation/run token, per-step provenance checks,
recorded parameters, and validator digest. Add `--include-queue` to also inspect
the root-wide archive queue; ordinary family verification does not inventory
unrelated families. Manifests use the workflow's `FAMILY.STEP.json` naming rule
within `artifact_provenance`. For a custom historical name, provide
`--manifest STEP=FILENAME`; its embedded family/step must still match. The
completion value `verified_declared_steps` covers only the required steps, not
the whole workflow or a scheduler job. Without `--attempt`, it does not attribute
artifacts to any particular historical execution. With `--attempt`, every step
must also have a matching manifest-bound receipt from an accepted attempt.
Its declarations are revalidated against current artifacts and parameters;
explicit stale reuse cannot establish completion. The receipt must belong to
the requested workspace, logical root, family, manifest, and step, and fall
within the attempt's recorded time interval. By default `verify --attempt`
requires the original path mapping. To inspect a container's workspace from its
host, explicitly pass `--recorded-workspace-root /workspace` together with the
actual `--workspace-root`. Every receipt must name that original workspace;
only its workspace-owned declarations are relocated and revalidated. The
caller must still establish which scheduler job/project owns the attempt.

For the gene-evolution terminal contract, additionally pass
`--profile gene-evolution-terminal-v1` and require both `summary_statistics` and
`tree_plot`. The API checks nonempty, consistently shaped branch/tree TSVs,
unique branch IDs, one tree-statistics row and PDF header/end framing, reading
up to 16 MiB per terminal artifact. `terminal_validation` reports this limited
structural coverage separately; it is not a full PDF renderer or a scientific
correctness proof. This profile is opt-in because disabled analyses may require
a different set of terminal steps.

Missing manifests return `evidence_missing`/`unverified`. Adopted manifests remain
`legacy_reusable`. A family published as running or failed, ambiguous manifests,
or an observed family/index generation change cannot establish completion.
Scientific `not_estimable`/disabled outcomes are not inferred from missing files:
their producer diagnostics remain available, and callers must choose the steps
appropriate to the configured analysis. Pending archival is not an analytical
failure.
Legacy adoption of optional outputs also remains an adoption-time baseline.
An optional output recorded as absent becoming a directory invalidates the
contract just as an unexpected file does.

Inspection uses optimistic archive reads instead of creating namespace-lock
files. It rejects index changes or pending index updates; ordinary workflow and
maintenance operations retain their existing locks. This is observational
evidence, not a transaction protecting a later write. Manual modification of
files or a concurrent run always requires fresh validation before reuse.

The query remembers filesystem signatures for accessed raw artifacts,
provenance and observation records, and ZIP payloads, then checks them again
before returning. Detected changes make the query unavailable, including when
an earlier input changes while later outputs are being checked. This check uses
filesystem metadata and does not claim a globally atomic snapshot or protection
against changes after the final check. Inspected store inputs are hashed from
their actual bytes, including ZIP members, rather than accepting indexed hashes
alone. This adds I/O proportional to the declared inputs being verified.

## Inspect the effective runtime record

```bash
python -B workflow/support/workflow_api.py runtime \
  --attempt workspace/output/observations/ATTEMPT_ID
```

The record includes registered entrypoint settings forwarded after overrides,
unset keys, configuration digest, actual streamed core-script digest, available
workflow-source fingerprints, version, scheduler metadata, and runtime identity.
Core-local derived defaults are explicitly outside the configuration snapshot.
The full environment is never dumped. Container digest is nullable and explicitly
caller-supplied when present; it is not a scientific capability qualification or
an independently verified image digest. Missing historical values remain unknown.

## Migration boundaries and verification

1. Keep code/container bindings of running and queued jobs intact. Stage a new
   runtime for later jobs; do not replace a runtime still used by an older job.
2. Existing projects may continue through their original entrypoints without
   enabling observations. No bulk conversion, deletion, or reanalysis is needed
   merely to add this API.
3. Inspect old artifacts and the intended restart contracts. Resolve only actual
   incompatibilities according to the project's existing policy.
4. Enable observations for new attempts when desired. Old and new records coexist;
   the API does not synthesize records for old executions.

Compatibility tests cover untracked and adopted outputs, provenance schema 1,
partial output sets, changed content and parameters, raw/ZIP verification,
workspace relocation, exact-attempt isolation, unknown schemas, missing records,
signal/exit-status preservation, and absence of workspace writes during queries.
They establish compatibility with those representations, not every historical
GeneGalleon release or every scientific toolchain. Production SIF and individual
project migrations still require their appropriate runtime validation.

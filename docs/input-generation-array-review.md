# Input-generation array review — 2026-09-09

This records the September 9 review. For the subsequent single download/prepare
job, staged local workers, and atomic namespace locks, see the current
[input-generation array guide](input-generation-arrays.md).

The review found correctness and isolation problems in the initial array changes.
The items below have been fixed and covered by regression tests. These are
software checks, not evidence that a particular HPC shared filesystem honors
cross-node locks.

| Finding | Correction and regression evidence |
| --- | --- |
| Different plans could use the same workspace shard indices, or different workspaces could write the same custom species directories. | Bind the workspace and canonical output directories to one plan. Hold locks for the actual output locations, including single-mode writers. Tests reject conflicting ownership without partially claiming a second workspace. |
| `01` and `1` could use different lock/shard names for the same task. | Validate and normalize task indices before acquiring worker locks. An integration test executes `01` and verifies only receipt `1.json` exists. |
| Resolved download JSON was reusable without proving which plan/task/species it belonged to. | Bind cached resolutions to plan hash, task index, provider and species; reject mismatches before formatting. |
| Local files referenced by manifests could change after planning; inputs could also change while a worker was running. | Hash local manifest sources at prepare, verify before use and before completion, and retain the hashes in completion receipts. Tests change sources before download and before receipt creation. Nonregular files are rejected without blocking. |
| Finalization could reread a changed source manifest or accept changed trait settings/shared lineage. | Reconstruct trait species from the frozen plan, skip source rediscovery in workers/finalize, freeze trait settings/files, and verify the prepared lineage hash. |
| Canonical tables were published before the last shared stages; staging in workspace scratch could require a non-atomic cross-filesystem move. | Stage beside each destination and rename only after all enabled final stages succeed. A forced report-stage failure preserves both previous canonical tables; retry then publishes successfully. |
| CDN redirects could switch a database request into an unrelated host bucket, losing the shared database limit/cooldown. | Carry a logical provider/database through redirects and isolate provider context between threads. Known database destinations override hints. Local HTTP tests cover redirects and 429/503 cooldowns. |
| Waiting for the limiter metadata lock could exceed the configured admission timeout indefinitely. | Use bounded nonblocking lock acquisition. Test a held metadata lock and verify timeout/recovery. Nonfinite Retry-After values use the normal backoff. |
| Malformed completion receipts could abort retry selection instead of selecting the affected task. | Validate receipt structure and treat invalid receipts as unfinished. |
| Large contiguous arrays produced unnecessarily long argument lists, and captured Slurm rejection details were lost. | Compress consecutive indices into ranges and retain scheduler errors/job IDs. Fake-scheduler tests cover successful dependencies and failed preparation/submission. |
| Species keys/download filenames could contain path traversal components. | Reject unsafe or hidden path components during planning; validate loaded plan species prefixes before execution. |
| Standard authentication headers could follow a redirect to a different origin. | Strip Authorization, Cookie and Proxy-Authorization on cross-origin redirects in the limited opener. Preserve them on same-origin redirects; test both paths with dummy credentials. |

An independently failing existing test was also investigated:
`test_format_species_inputs_uses_locus_tag_for_genbank_style_ncbi_cds` used
inconsistent negative-strand CDS phases. Its block lengths of 214, 136 and 79
require phases 0, 2 and 1. Only the fixture values were corrected; production GFF
validation was not relaxed.

## Verification

Validation used `local/genegalleon:dev` with the working checkout mounted into
Docker. The broad regression run passed 474 tests. After the final path and
redirect-credential guards, the affected 55-test selection also passed (these
runs overlap). Ruff, Bash syntax, and `git diff --check` passed. Coverage includes actual formatting/validation helpers, local HTTP
servers, multiple processes, SIGKILL, incomplete and corrupted state, frozen
manifest behavior, single/array biological-output equivalence, shell static
checks, and the existing provider download suites. Integration orchestration
uses fake Slurm/BUSCO/seqkit/report executables; it does not establish live Slurm
or real BUSCO performance. Ruff, Bash syntax and whitespace checks are separate
from the runtime tests.

## Remaining deployment limits

- Cross-node `flock` behavior and SIF execution on the target HPC have not been
  tested. Node-local scratch cannot coordinate database limits between nodes.
- No real scheduler jobs or representative large-data performance benchmarks
  were run. Slurm's configured `MaxArraySize` still applies; the helper reports
  rejection rather than silently splitting or renumbering a plan.
- Direct URLs on different unrecognized hosts cannot automatically be inferred
  to represent one database. Use the appropriate provider in the manifest when
  that logical relationship is known.
- Canonical table renames are atomic individually. Optional trait/report stages
  can leave their own partial files on failure, and publishing multiple files
  is not one filesystem-wide transaction. The workflow returns failure and
  retains worker receipts for retry.
- Input hashing and retaining raw download evidence consume I/O/storage. No
  speedup claim is made without representative before/after measurements.

See [operation and restart instructions](input-generation-arrays.md).

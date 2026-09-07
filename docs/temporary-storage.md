# Temporary computation storage

`GG_COMMON_TMP_ROOT` selects disposable workflow scratch storage. The default is
`workspace` on every site, including NIG, Shirokane and audrey1. No hostname-based
scratch override is applied. With the default repository workspace, paths are
below `<repository>/workspace/output`; an explicit `gg_workspace_dir` changes
that workspace root.

| Value | Location |
| --- | --- |
| `workspace` (default) | Existing workflow-specific output `tmp` directories |
| `/tmp` | Private GeneGalleon directories inside host `/tmp` |
| `/scratch/user` | Private GeneGalleon directories inside that existing directory |
| `env` | Execution host's `TMPDIR`, resolved immediately before container launch |

For example, in a job script running on the compute node:

```bash
GG_COMMON_TMP_ROOT=env bash workflow/gg_gene_evolution_entrypoint.sh
# Or choose an existing writable directory explicitly:
GG_COMMON_TMP_ROOT=/scratch/user bash workflow/gg_genome_evolution_entrypoint.sh
```

`env` fails when `TMPDIR` is unset. An explicit root must be absolute, existing,
writable and searchable. Commas, colons and newlines are rejected because the
container bind syntax cannot represent them safely. Invalid roots fail without
falling back to workspace or `/tmp`. The selected root is mounted at `/gg_tmp`
with Apptainer, Singularity, and the Docker-backed runtime. Custom adapters must
use `exec`; `shell` adapters are rejected for external scratch. Use the entrypoints;
direct core execution does not set up or supervise external scratch.

Each external run uses a private directory scoped by user, workspace, workflow,
gene-evolution mode (where applicable), and execution. The host path, container path and available bytes are logged.
Child tools receive `TMPDIR`, `TMP` and `TEMP` pointing inside that run. Main
computation directories in input generation, genome annotation, transcriptome
generation, gene evolution, genome evolution and fractionation bias are routed
there too, including the derived protein staging directories used by genome
evolution. Tools that insist on writing next to their output may still write
there; this option is not a promise to relocate every temporary byte.

Input-generation task plans, array summary shards, shared download caches,
validation stamps, locks and ZIP materialization receipts remain in the
workspace. Output publication staging also remains beside the destination so
atomic rename works across filesystems. A completed scratch computation is not
a durable checkpoint: only published outputs are authoritative.

Successful runs remove their external run directory unless `delete_tmp_dir=0`.
Failures and interrupted runs retain scratch. A forwarded termination signal
remains a failure even if a child's signal handler returns exit code zero.
Gene evolution also honors `delete_preexisting_tmp_dir`: 1 removes idle scratch
for the same task number, mode and input assignment before a new run; 0 reuses
the most recent matching idle directory. Query mode also requires an unchanged
query filename inventory and sorting locale; orthogroup mode checks the family
assigned to that table row. Active runs are protected with
file locks inherited by the computation, including through the metrics monitor.
Changing workspace, workflow, scratch root or compute node does not
recover node-local files from a different location.

Gene evolution's existing `gene_family_tmp_retention_days`,
`gene_family_tmp_max_dirs`, `gene_family_tmp_max_bytes` and
`gene_family_tmp_max_files` limits also apply to idle external runs, scoped by
workspace, workflow and gene-evolution mode at that scratch location. Cleanup
occurs on invocation and completion, not as a background service. Other workflows retain failed runs
until manually removed or removed by the site's policy. Never remove an active
run. `delete_tmp_dir=0` prevents successful-run cleanup by GeneGalleon but does
not disable retention limits or site cleanup.

The default workspace paths are:

| Workflow | Default computation directory (relative to workspace) |
| --- | --- |
| Input generation | `output/input_generation/tmp` |
| Genome annotation | `output/tmp/<task>_<species>` |
| Transcriptome generation | `output/transcriptome_assembly/tmp/<task>_<species>` |
| Gene evolution | `output/{orthogroup,query2family}/tmp/<task>_<family>` |
| Genome evolution | `output/species_tree/tmp`, `downloads/tmp/species_protein` and its staging variants |
| Fractionation bias | `output/tmp/kffractbias/<task>_<analysis>.<unique>` |

Node-local storage may disappear at job termination even when GeneGalleon keeps
it. Choose workspace storage when retained intermediate files must remain
accessible after the job; choose site-approved scratch explicitly when disposable
local computation is desired. No performance improvement is claimed without a
representative benchmark.

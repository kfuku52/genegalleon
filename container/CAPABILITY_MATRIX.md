# Capability Matrix

This matrix summarizes intended behavior for the multi-arch `GeneGalleon` runtime.

| Tier | Scope | amd64 | arm64 |
|---|---|---|---|
| Core | Single `base` env install for required commands (Git-sourced programs follow the moving branches in `source_branches.env`) | Expected pass | Expected pass with arm64 profile (excludes `Trinity`, `jellyfish`) |
| Advanced optional | none (strict profile: tools referenced by pipeline scripts are promoted to required) | N/A | N/A |
| External DB assets | `/usr/local/db/Pfam_LE`, `/usr/local/db/uniprot_sprot.pep`, `/usr/local/db/jaspar` | Manual population | Manual population |

Validation artifacts produced at build time:
- `/opt/pg/logs/runtime_validation_amd64.tsv` or `/opt/pg/logs/runtime_validation_arm64.tsv`
- `/opt/pg/logs/failed_optional_*.txt`

# Capability Matrix

This matrix summarizes intended behavior for the multi-arch `GeneGalleon` scaffold.

| Tier | Scope | amd64 | arm64 |
|---|---|---|---|
| Core | Single `base` env install for required commands (GitHub installs for `amalgkit`, `cdskit`, `csubst`, `nwkit`, `kftools`, `rkftools`; `amalgkit` auto-selects newer commit between `master` and `kfdevel`) | Expected pass | Expected pass with arm64 profile (excludes `Trinity`, `jellyfish`) |
| Advanced optional | none (strict profile: tools referenced by pipeline scripts are promoted to required) | N/A | N/A |
| External DB assets | `/usr/local/db/Pfam_LE`, `/usr/local/db/uniprot_sprot.pep`, `/usr/local/db/jaspar` | Manual population | Manual population |

Validation artifacts produced at build time:
- `/opt/pg/logs/runtime_validation_amd64.tsv` or `/opt/pg/logs/runtime_validation_arm64.tsv`
- `/opt/pg/logs/failed_optional_*.txt`

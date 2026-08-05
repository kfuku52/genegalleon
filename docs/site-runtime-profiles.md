# Site Runtime Profiles

GeneGalleon applies small host-specific runtime adaptations from
`workflow/support/gg_site_runtime.sh` before starting the container. Site
profiles do not change scientific parameters or workflow stages.

## Profile selection

The profile is selected in this order:

1. explicit `GG_SITE_PROFILE`,
2. known scheduler or hostname patterns,
3. `default`.

Use an explicit override when automatic hostname detection is not appropriate:

```bash
GG_SITE_PROFILE=nig qsub workflow/gg_gene_evolution_entrypoint.sh
```

The startup summary prints the selected profile and container runtime.

## Available profiles

### `shirokane`

- detected from the SHIROKANE AGE environment,
- initializes Environment Modules and loads `apptainer` on compute nodes,
- uses the standard Apptainer execution command.

See the [SHIROKANE AGE Guide](shirokane-age.md) for transfer, submission, and
resource examples.

### `nig`

- detected from known National Institute of Genetics hostnames, or selected
  with `GG_SITE_PROFILE=nig`,
- changes to `PBS_O_WORKDIR` for scheduler-spooled jobs,
- discovers versioned `apptainer` or `singularity` installations under
  `/opt/pkg`, with the legacy Singularity path as a fallback,
- binds available UGE/AGE spool and package paths needed by the compute-node
  runtime.

Submit the normal `workflow/gg_*_entrypoint.sh` directly with the site's
`qsub`; separate NIG-specific wrappers are not required.

### `nhr-fau`

- detected from `*.nhr.fau.de` hostnames,
- starts the container with `--contain`.

### `default`

- changes to `PBS_O_WORKDIR` when present,
- otherwise uses the normal Apptainer/Singularity or Docker-backed runtime
  without additional site flags.

## Troubleshooting

If the wrong profile is selected, set `GG_SITE_PROFILE` explicitly and inspect
the startup runtime summary. If the right profile is selected but no container
runtime is found, check the compute-node environment rather than the login-node
environment and confirm that the SIF path is visible from the job.

Site-specific behavior is implementation-level runtime setup. Keep project
parameters in the entrypoint config blocks or their scoped environment
overrides rather than adding site-specific scientific defaults here.

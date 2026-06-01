# Agent / Developer Validation Policy

Use the repository `genegalleon.sif` runtime for workflow validation, integration tests, R helper checks, and toolchain-dependent behavior.

Do not treat host-local Python, R, or command-line tool behavior as authoritative for GeneGalleon runtime compatibility. Host-local checks are acceptable for quick syntax or narrow static checks, but container-backed checks are the source of truth when dependencies matter.

Preferred validation entrypoint:

```bash
bash workflow/tests/run_in_sif.sh python -m pytest -q workflow/tests/test_hgt_end_to_end.py
bash workflow/tests/run_in_sif.sh Rscript workflow/tests/test_treevis_main.R
```

If `genegalleon.sif` or an Apptainer/Singularity runtime is unavailable, report that clearly and do not conclude SIF compatibility from host-local results.

Do not add backward-compatibility workarounds for older dependency behavior.

If the root cause is in a dependency program or package, do not patch GeneGalleon to absorb it. Report it as dependency-side so the dependency can be fixed or updated instead.

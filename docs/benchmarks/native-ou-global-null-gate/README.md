# Native AIC global-null bootstrap gate

This study evaluates the optional `--global-null-gate` in NWKIT. It reruns the
whole AIC search under a fitted no-shift null and returns zero shifts after
non-rejection. It leaves AIC's choice unchanged after rejection.

- [Protocol and commands](PROTOCOL.md)
- [Results and limitations](RESULTS.md)
- [Validation evidence](VALIDATION.md)
- [Machine-readable pilot summary](summary.json)
- [Implementation delta](implementation.patch) from the previous delivered AIC path
- `source.tar.gz`: complete measured NWKIT Python source, verified against
  `source-manifest.json`; `image-id.txt`: GeneGalleon Docker image identity
- `*-estimated-alpha.json` / `*-fixed-alpha.json`: paired observed metrics,
  input hashes, full bootstrap statistics and decision records
- `*-B199.json`: separate default-count smoke checks
- `run.py` and `summarize.py`: execution and independent p-value/decision checks

For replay in a clean directory, extract `source.tar.gz` to a new source folder
and substitute that folder for the NWKIT mount in the protocol command. Use a
copy of this study directory without generated result JSON files: the runner
refuses to overwrite results. Keep its sibling `native-ou-aic-improvement`
fixtures at the same relative location. Run `summarize.py` after the 48 paired
estimated-alpha runs finish. This is a development pilot using previously
studied fixtures, not an independent held-out calibration study.

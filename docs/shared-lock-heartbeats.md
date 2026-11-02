# Shared-lock heartbeat shutdown

Shell publishers start an ownership-aware Python heartbeat and wait for its
startup acknowledgement. A private `gg-lock-heartbeat.*` file in the job's
`TMPDIR` defines that helper's lifetime. Stopping removes this inode before
sending TERM and reaping the child. The helper checks its parent and the open
control inode every 0.1 seconds, so a lost shutdown signal cannot leave the
publisher waiting indefinitely while its output stripe remains occupied.

The helper removes its own control file on exit. It does not release or steal
the publication lease when the workflow parent disappears; normal token-aware
recovery remains responsible for that lease. No persistent workspace control
files, storage formats, stale thresholds, scientific parameters or artifact
recipes change. Failure to acknowledge startup aborts bundle publication and
rolls back its staged sources and previous outputs.

The regression suite covers immediate shutdown, a deliberately dropped TERM,
parent death, replacement control files, unsafe control paths, ownership
handoff, live lease renewal, and publication rollback:

```bash
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_shared_lock.py \
  workflow/tests/test_gg_shared_lock_shell.py \
  workflow/tests/test_shared_namespace_lock.py \
  workflow/tests/test_gg_util_paths.py
```

A synthetic Docker comparison ran ten sequential lease-protected publications
per trial, with three trials per source and signal condition. Successful
outputs had identical SHA-256 in all nine completed trials. With normal TERM
delivery the median changed from 0.963 to 1.570 seconds per ten publications
(about 0.061 additional seconds per publication). With TERM deliberately
dropped, all three previous-implementation trials reached the three-second
timeout; all three revised trials finished, with a median of 2.487 seconds.
These measurements isolate lease lifecycle overhead, not whole-workflow
throughput. A Python helper remains resident for each active heartbeat, in
exchange for removing repeated shell sleep/helper process creation.

The comparison used the same fresh Docker runtime and temporary fixtures. It
does not establish performance or compatibility of an HPC SIF; deployments
into a frozen project also qualify the affected tests against that exact SIF.

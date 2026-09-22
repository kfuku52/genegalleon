# Runtime validation

Use a GeneGalleon container for workflow integration, R helpers, and
toolchain-dependent behavior. Host syntax and narrow static checks are useful
feedback but do not establish runtime compatibility.

From the repository root, use `bash ./dev check <suite>` or
`bash workflow/tests/run_in_runtime.sh <command>`. These select a usable SIF on
Linux/HPC or Docker on macOS and check runtime freshness before execution.
For prerequisites, focused commands, build setup, freshness policy, and expensive
checks, use [Development and Tests](development-and-tests.md#choose-checks-for-a-change).
Keep executable procedures there rather than maintaining a second copy here.

Report the runtime actually used. Docker results do not establish SIF
compatibility. If no suitable runtime is available, report the missing
prerequisite and the unverified checks. Do not silently use host tools, bypass
freshness, or claim a skipped check as successful.

Do not add backward-compatibility workarounds for older dependency behavior.
When the root cause belongs to a dependency, fix/update the owning dependency
within the authorized scope or report the blocker; do not absorb it into
GeneGalleon.

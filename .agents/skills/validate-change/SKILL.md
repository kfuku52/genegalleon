---
name: validate-change
description: Select and execute focused GeneGalleon checks after a repository change, distinguishing static, Docker, and SIF evidence. Use for change validation, not full scientific analyses or performance benchmarking.
---

# Validate a GeneGalleon Change

Inputs: the changed files or intended diff, the behavior to verify, and available
runtime/resource constraints. Work from the repository root. If the scope is not
supplied, inspect `git status --short` and the relevant diff before choosing tests.

1. Read [the change-to-check table](../../../docs/development-and-tests.md#choose-checks-for-a-change)
   and only the linked runtime/domain guidance needed for this change. Follow
   affected callers and output consumers when a shared contract changes.
2. Find candidate tests with targeted `rg` searches for the changed helper or
   behavior. Read their fixtures and resource requirements. Check file membership
   in `workflow/tests/conftest.py` and `validation_manifest.json`; function markers
   alone do not select a lane.
3. Preview unfamiliar selection with `run_checks.py SUITE --list` and, when
   useful, `dev check SUITE FILE --collect-only`. An explicit file in the wrong
   lane may select no tests. Use `run_in_runtime.sh python -m pytest` for a focused
   Python check without lane filtering. Runtime/full lanes also run R checks.
4. Execute the relevant bounded commands from the guide. Use existing temporary
   fixtures; do not run a real entrypoint against curated data as a quick check.
   Preserve failures and inspect skipped tests. Builds and model downloads need
   resource scope appropriate to the user's task, not automatic escalation from
   a small check.
5. Review `git diff --check` and status for unintended outputs. Report the exact
   commands, selected behavior, runtime, results, and remaining verification.

For a validation-tooling change, a worked starting point is:

```bash
python3 workflow/tests/run_checks.py fast --list workflow/tests/test_validation_runner.py workflow/tests/test_development_tooling.py -x
bash ./dev check fast workflow/tests/test_validation_runner.py workflow/tests/test_development_tooling.py -x
```

The expected outcome is that both files are selected and their checks pass;
test counts can change. This example exercises routing and failure propagation,
not scientific tool compatibility. Choose different tests for different changes.

Deliver a concise validation report, not a new permanent log or test framework.
Classify evidence as executed behavior (with SIF/Docker identified), static-only,
or blocked/unverified. For missing tools, stale images, dependency failures, or
unaffordable runtime needs, record the actionable prerequisite and continue
independent checks. Never weaken strict checks, disable freshness implicitly,
or turn host results into container evidence. Publication remains governed by
the existing prepare-github-push skill.

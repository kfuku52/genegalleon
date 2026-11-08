"""Measure version inventory reuse across synthetic scheduler scratch contexts.

The inventory and runtime are fixtures, not scientific inputs or an actual SIF.
Native shell coordination, image hashing, cache verification and SQLite packing
are exercised. A fixed probe delay models inventory startup separately from
the software record's byte-for-byte equivalence check.
"""

import argparse
import hashlib
import json
import os
import resource
import shlex
import stat
import statistics
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def file_sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def shell_fixture(source, work, delay, runtime):
    quote = shlex.quote
    return f"""
set -euo pipefail
source {quote(str(source / 'workflow/support/gg_util.sh'))}
gg_support_dir={quote(str(source / 'workflow/support'))}
gg_workflow_dir={quote(str(source / 'workflow'))}
gg_workspace_dir={quote(str(work / 'workspace'))}
gg_container_image_path={quote(str(work / 'fixture.sif'))}
singularity_command=({quote(str(runtime))} exec --cleanenv
  --bind "${{GG_TMP_HOST_ROOT}}:${{GG_TMP_MOUNT}}"
  --bind "${{gg_workspace_dir}}:/workspace" --bind "${{gg_workflow_dir}}:/script")
gg_shared_lock_python() {{ printf '%s\\n' {quote(sys.executable)}; }}
gg_container_shell_command_is_set() {{ return 0; }}
gg_container_shell_command_runtime_bin() {{ printf '%s\\n' {quote(str(runtime))}; }}
gg_container_bind_destination_exists() {{ return 0; }}
gg_print_version_summary() {{ printf '%s\\n' 'synthetic runtime inventory'; }}
date() {{ if [[ $# -eq 0 ]]; then printf '%s\\n' 'fixed synthetic benchmark time'; else command date "$@"; fi; }}
gg_run_container_shell_script() {{
  printf '%s\\n' called > {quote(str(work / 'calls'))}/"${{BASHPID:-$$}}"
  sleep {delay}
  cat {quote(str(work / 'payload'))}
}}
gg_require_versions_dump synthetic-benchmark
"""


def trial(args):
    source = args.source.resolve()
    work = args.workdir.resolve()
    work.mkdir(mode=0o700)
    (work / "workspace").mkdir(mode=0o700)
    (work / "fixture.sif").write_bytes(b"synthetic immutable runtime\n" * 64)
    runtime = args.runtime_fixture.resolve()
    unit = b"synthetic package 1.0\n"
    payload = (unit * ((args.payload_bytes + len(unit) - 1) // len(unit)))[:args.payload_bytes]
    (work / "payload").write_bytes(payload)
    (work / "calls").mkdir(mode=0o700)
    support = source / "workflow/support"
    sys.path.insert(0, str(support))
    import metadata_store
    import version_report_identity

    cache = work / "workspace/downloads/version_inventory"
    version_report_identity.private_cache(cache)
    metadata_store.pack(cache, writers_stopped=True)
    if args.cache_limit is not None:
        metadata_store.pack(cache, writers_stopped=True, cache_limit=args.cache_limit)
    script = shell_fixture(source, work, args.probe_delay, runtime)
    started = time.perf_counter()
    def run_context(task):
        scratch = work / "scratch" / ("job.1.task." + str(task))
        scratch.mkdir(mode=0o700, parents=True)
        env = dict(os.environ, TMPDIR=str(scratch), TMP=str(scratch), TEMP=str(scratch),
                   JOB_ID="900001", SGE_TASK_ID=str(task), GG_JOB_ID="900001",
                   GG_ARRAY_TASK_ID=str(task), GG_COMMON_TMP_ROOT=str(scratch),
                   GG_TMP_HOST_ROOT=str(scratch), GG_TMP_MOUNT="/gg_tmp",
                   GG_TMP_WORKSPACE_ID=str(work / "workspace"),
                   SINGULARITYENV_GG_VERSION="0.0.0-synthetic-benchmark",
                   APPTAINERENV_GG_VERSION="0.0.0-synthetic-benchmark")
        if args.distinct_contexts:
            env["GG_BENCHMARK_VERSION_CONTEXT"] = str(task)
        for prefix in ("SINGULARITYENV_", "APPTAINERENV_"):
            for key in ("GG_COMMON_TMP_ROOT", "GG_TMP_HOST_ROOT", "GG_TMP_MOUNT"):
                env[prefix + key] = env[key]
            env[prefix + "TMPDIR"] = "/gg_tmp/synthetic-runtime/tmp"
        result = subprocess.run([args.bash, "-c", script], env=env, cwd=source,
                                capture_output=True, text=True, timeout=120)
        if result.returncode:
            raise RuntimeError(result.stdout[-2000:] + result.stderr[-4000:])
        reports = [line.split(" ", 3)[-1] for line in result.stdout.splitlines()
                   if line.startswith(("gg_trigger_versions_dump: wrote ",
                                       "gg_trigger_versions_dump: reused verified "))]
        if len(reports) != 1:
            raise RuntimeError("Missing native version inventory receipt: " + result.stdout[-2000:])
        report = Path(reports[0])
        if not report.is_file():
            raise RuntimeError("Native inventory receipt does not name a regular report")
        return file_sha(report)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        output_sha = list(pool.map(run_context, range(args.contexts)))
    elapsed = time.perf_counter() - started
    invocations = list((work / "calls").iterdir())
    assert all(path.read_bytes() == b"called\n" for path in invocations)
    run_provenance = list((work / "workspace/output/versions").glob("run.*.log"))
    if len(run_provenance) != args.contexts or len(set(output_sha)) != 1:
        raise AssertionError("Inventory changed or per-invocation provenance was lost")
    scale = 1 if sys.platform == "darwin" else 1024
    maximum = max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss) * scale
    return {"seconds": elapsed, "probe_invocations": len(invocations),
            "output_sha256": output_sha, "run_provenance_count": len(run_provenance),
            "max_single_process_rss_bytes": maximum}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--runtime-fixture", type=Path, required=True,
                        help="Use the same private synthetic runtime path for paired comparisons")
    parser.add_argument("--bash", default="bash")
    parser.add_argument("--contexts", type=int, default=8)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--distinct-contexts", action="store_true",
                        help="Retain a different software context for each task")
    parser.add_argument("--cache-limit", type=int)
    parser.add_argument("--payload-bytes", type=int, default=213000)
    parser.add_argument("--probe-delay", type=float, default=0.1)
    parser.add_argument("--trial", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if (args.contexts < 2 or args.repeats < 1 or args.payload_bytes < 1 or args.probe_delay < 0
            or not 1 <= args.workers <= 32 or args.cache_limit is not None and args.cache_limit < 1):
        parser.error("Require multiple contexts, positive repeats/payload and nonnegative delay")
    fixture = args.runtime_fixture.absolute()
    expected = b"#!/bin/sh\nprintf '%s\\n' 'synthetic runtime inspect identity'\n"
    if not fixture.exists():
        with fixture.open("xb") as handle:
            handle.write(expected)
            os.fchmod(handle.fileno(), 0o500)
    if (fixture.resolve() != fixture or not stat.S_ISREG(fixture.stat().st_mode)
            or fixture.stat().st_uid != os.getuid()
            or fixture.stat().st_nlink != 1 or fixture.stat().st_mode & 0o022
            or fixture.read_bytes() != expected):
        parser.error("Synthetic runtime must be an owned regular, unchanged fixture")
    if args.trial:
        print(json.dumps(trial(args)))
        return
    root = args.workdir.resolve()
    root.mkdir(mode=0o700)
    records = []
    for index in range(args.repeats + 1):
        command = [sys.executable, "-B", str(Path(__file__).resolve()), "--trial",
                   "--source", str(args.source.resolve()), "--workdir", str(root / str(index)),
                   "--runtime-fixture", str(fixture),
                   "--bash", args.bash, "--contexts", str(args.contexts), "--workers", str(args.workers),
                   "--payload-bytes", str(args.payload_bytes), "--probe-delay", str(args.probe_delay)]
        if args.distinct_contexts:
            command.append("--distinct-contexts")
        if args.cache_limit is not None:
            command += ["--cache-limit", str(args.cache_limit)]
        result = subprocess.run(command, capture_output=True, text=True, timeout=600)
        if result.returncode:
            raise RuntimeError(result.stderr[-6000:])
        records.append(json.loads(result.stdout))
    measured = records[1:]
    if any(row["output_sha256"] != records[0]["output_sha256"] for row in measured):
        raise AssertionError("Trials produced different logical inventories")
    report = {"synthetic_workload": True, "command": sys.argv,
              "python": sys.version, "contexts": args.contexts, "workers": args.workers,
              "distinct_contexts": args.distinct_contexts, "cache_limit": args.cache_limit,
              "runtime_fixture": str(fixture), "runtime_fixture_sha256": file_sha(fixture),
              "probe_delay_seconds": args.probe_delay, "payload_bytes": args.payload_bytes,
              "warmup": records[0], "measured": measured,
              "median_seconds": statistics.median(row["seconds"] for row in measured),
              "median_max_single_process_rss_bytes": statistics.median(
                  row["max_single_process_rss_bytes"] for row in measured),
              "source_files": {str(path.relative_to(args.source.resolve())): file_sha(path)
                               for path in (args.source.resolve() / "workflow/support/version_report_identity.py",
                                            args.source.resolve() / "workflow/support/gg_util/08_execution_reporting.sh",
                                            Path(__file__).resolve())}}
    path = root / "report.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Preview or submit input-generation prepare -> species array -> finalize on Slurm."""
import argparse
import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "support"))
from input_generation_array_state import digest, load_plan, prepared, verify_receipt  # noqa: E402


def array_expression(indices):
    ranges = []
    for index in indices:
        if ranges and index == ranges[-1][1] + 1:
            ranges[-1][1] = index
        else:
            ranges.append([index, index])
    return ",".join(str(first) if first == last else f"{first}-{last}" for first, last in ranges)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-plan", required=True)
    parser.add_argument("--entrypoint", default=str(Path(__file__).resolve().with_name("gg_input_generation_entrypoint.sh")))
    parser.add_argument("--max-running", type=int, default=8)
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--memory", default="32G", help="Total memory per worker")
    parser.add_argument("--time", default="3-00:00:00")
    parser.add_argument("--partition", default="", help="Slurm partition (default: scheduler default)")
    parser.add_argument("--retry", action="store_true", help="Skip prepare and submit only workers without verified receipts")
    parser.add_argument("--submit", action="store_true", help="Submit jobs; default prints a dry-run preview")
    args = parser.parse_args()
    if min(args.cpus, args.max_running) < 1:
        parser.error("CPU and concurrency counts must be positive")
    plan_path = Path(args.task_plan).expanduser().resolve()
    base = ["sbatch", "--parsable", "--cpus-per-task=" + str(args.cpus), "--mem=" + args.memory,
            "--time=" + args.time]
    if args.partition:
        base += ["--partition=" + args.partition]
    env = os.environ.copy()
    env["GG_INPUT_TASK_PLAN_OUTPUT"] = str(plan_path)

    def command(mode, extra):
        # Values are inherited via the process environment, avoiding comma/space
        # escaping problems in Slurm's --export parser.
        return base + ["--export=ALL", "--job-name=gg_input_" + mode] + extra + ["--wrap=" + shlex.join(["exec", "bash", str(Path(args.entrypoint).resolve())])]

    def dispatch(mode, extra):
        cmd = command(mode, extra)
        print("GG_INPUT_TASK_PLAN_OUTPUT=" + shlex.quote(str(plan_path)) + " GG_INPUT_INPUT_GENERATION_MODE=" + mode + " " + shlex.join(cmd), flush=True)
        if not args.submit:
            return "WORKER_JOB_ID"
        result = subprocess.run(cmd, env={**env, "GG_INPUT_INPUT_GENERATION_MODE": mode}, capture_output=True, text=True, check=False)
        job_id = result.stdout.strip().split(";")[0]
        valid_job_id = re.fullmatch(r"\d+", job_id)
        if valid_job_id:
            print("Submitted " + mode + ": " + job_id, flush=True)
        if result.returncode:
            raise RuntimeError(f"{mode} failed (exit {result.returncode}, job {job_id or 'not submitted'}): {result.stderr.strip()}")
        if not valid_job_id:
            raise RuntimeError("Unrecognized sbatch job ID: " + result.stdout)
        return job_id

    if not args.retry:
        # The task count is available only after prepare. --wait requires a
        # successful prepare before anything dependent is submitted.
        dispatch("array_prepare", ["--wait"])
    if not plan_path.exists():
        if args.submit or args.retry:
            parser.error("Task plan not found: " + str(plan_path))
        print("After successful prepare, read " + str(plan_path) + " for N species tasks.")
        worker_id = dispatch("array_worker", ["--array=1-N%" + str(args.max_running)])
    else:
        plan = load_plan(plan_path)
        if (args.submit or args.retry) and not prepared(plan_path):
            parser.error("Prepare has not completed for this plan/settings; run without --retry first")
        plan_sha256 = digest(plan_path)
        pending = [i for i in range(1, plan["task_count"] + 1) if not args.retry or not verify_receipt(plan_path, i, plan_sha256)]
        print(json.dumps({"species_count": plan["task_count"], "selected_tasks": pending, "cpus_per_task": args.cpus,
                          "memory_per_task": args.memory, "max_running": args.max_running}))
        worker_id = dispatch("array_worker", ["--array=" + array_expression(pending) + "%" + str(args.max_running)]) if pending else ""
    dispatch("array_finalize", (["--dependency=afterok:" + worker_id] if worker_id else []))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Preview or submit input-generation download/prepare -> compute array -> finalize on Slurm."""
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


def ensure_no_active_legacy_worker_array():
    """Fail closed before submitting a retry while any legacy-named array is active."""
    result = subprocess.run(
        ["squeue", "--noheader", "--me", "--name", "gg_input_array_worker", "--format", "%i"],
        capture_output=True, text=True, check=False,
    )
    if result.returncode:
        raise RuntimeError("Cannot verify active input-generation arrays: " + result.stderr.strip())
    active = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    if active:
        raise RuntimeError(
            "An input-generation worker array is still active or pending ({}); retry after it exits.".format(
                ", ".join(active[:8])
            )
        )


def ensure_no_active_rescue_array():
    result = subprocess.run(["squeue", "--noheader", "--me", "--name",
                             "gg_input_rescue_synteny,gg_input_rescue_models,gg_input_rescue_finalize", "--format", "%i"],
                            capture_output=True, text=True, check=False)
    if result.returncode or result.stdout.strip():
        raise RuntimeError("Cannot submit rescue retry while rescue jobs are active or scheduler status is unavailable")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-plan", required=True)
    parser.add_argument("--entrypoint", default=str(Path(__file__).resolve().with_name("gg_input_generation_entrypoint.sh")))
    parser.add_argument("--max-running", type=int, default=8)
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--memory", default="32G", help="Total memory per worker")
    parser.add_argument("--prepare-cpus", type=int, help="Download/prepare CPUs and download workers (default: --cpus)")
    parser.add_argument("--prepare-memory", help="Total download/prepare memory (default: --memory)")
    parser.add_argument("--prepare-partition", help="Network-enabled download/prepare partition (default: --partition)")
    parser.add_argument("--time", default="3-00:00:00")
    parser.add_argument("--partition", default="", help="Slurm partition (default: scheduler default)")
    parser.add_argument("--retry", action="store_true", help="Skip prepare and submit only workers without verified receipts")
    parser.add_argument("--refinement", action="store_true", help="After initial finalize, preserve all isoforms, build correspondence, predict and finalize effective representatives")
    parser.add_argument("--refinement-output", type=Path, help="Refinement directory; default: TASK_PLAN parent/../gene_model_refinement")
    parser.add_argument("--rescue", action="store_true", help="After initial finalize, submit sparse synteny -> species rescue -> rescue finalize arrays")
    parser.add_argument("--rescue-output", type=Path, help="Rescue directory; default: TASK_PLAN parent/../gene_model_rescue")
    parser.add_argument("--rescue-max-running", type=int, help="Concurrent species rescue workers (default: --max-running); independent of synteny/formatting")
    parser.add_argument("--rescue-cpus", type=int, help="CPUs per species rescue/finalizer (default: --cpus)")
    parser.add_argument("--rescue-memory", help="Total memory per species rescue/finalizer (default: --memory)")
    parser.add_argument("--submit", action="store_true", help="Submit jobs; default prints a dry-run preview")
    args = parser.parse_args()
    if min(args.cpus, args.max_running, args.prepare_cpus if args.prepare_cpus is not None else args.cpus,
           args.rescue_cpus if args.rescue_cpus is not None else args.cpus,
           args.rescue_max_running if args.rescue_max_running is not None else args.max_running) < 1:
        parser.error("CPU and concurrency counts must be positive")
    plan_path = Path(args.task_plan).expanduser().resolve()
    env = os.environ.copy()
    env["GG_INPUT_TASK_PLAN_OUTPUT"] = str(plan_path)
    rescue_output = (args.rescue_output or plan_path.parent.parent / "gene_model_rescue").resolve()
    refinement_output = (args.refinement_output or plan_path.parent.parent / "gene_model_refinement").resolve()
    if args.refinement:
        env["GG_INPUT_RUN_GENE_MODEL_REFINEMENT"] = "1"
        env["GG_INPUT_GENE_MODEL_REFINEMENT_DIR"] = str(refinement_output)
        if args.retry and args.submit:
            result = subprocess.run(["squeue", "--noheader", "--me", "--name", "gg_input_rescue_synteny,gg_input_refinement_prepare,gg_input_refinement_catalog,gg_input_refinement_correspondence,gg_input_refinement_predict,gg_input_refinement_finalize", "--format", "%i"], capture_output=True, text=True, check=False)
            if result.returncode or result.stdout.strip():
                raise RuntimeError("Cannot retry while refinement workers are active or scheduler status is unavailable")
    if args.rescue:
        env["GG_INPUT_RUN_GENE_MODEL_RESCUE"] = "1"
        env["GG_INPUT_GENE_MODEL_RESCUE_DIR"] = str(rescue_output)
        if args.retry and args.submit:
            ensure_no_active_rescue_array()

    def command(mode, extra):
        preparing = mode == "array_prepare"
        model_worker = mode in {"rescue_models", "rescue_finalize", "refinement_predict", "refinement_finalize"}
        cpus = (args.prepare_cpus or args.cpus) if preparing else args.cpus
        memory = (args.prepare_memory or args.memory) if preparing else args.memory
        if model_worker:
            cpus = args.rescue_cpus or cpus
            memory = args.rescue_memory or memory
        partition = (args.prepare_partition if args.prepare_partition is not None else args.partition) if preparing else args.partition
        base = ["sbatch", "--parsable", "--cpus-per-task=" + str(cpus), "--mem=" + memory, "--time=" + args.time]
        if partition:
            base += ["--partition=" + partition]
        # Values are inherited via the process environment, avoiding comma/space
        # escaping problems in Slurm's --export parser.
        return base + ["--export=ALL", "--job-name=gg_input_" + mode] + extra + ["--wrap=" + shlex.join(["exec", "bash", str(Path(args.entrypoint).resolve())])]

    def dispatch(mode, extra):
        cmd = command(mode, extra)
        rescue_env = ("GG_INPUT_RUN_GENE_MODEL_RESCUE=1 GG_INPUT_GENE_MODEL_RESCUE_DIR=" + shlex.quote(env.get("GG_INPUT_GENE_MODEL_RESCUE_DIR", str(rescue_output))) + " ") if args.rescue else ""
        refinement_env = ("GG_INPUT_RUN_GENE_MODEL_REFINEMENT=1 GG_INPUT_GENE_MODEL_REFINEMENT_DIR=" + shlex.quote(str(refinement_output)) + " ") if args.refinement else ""
        anchor_env = ("GG_INPUT_GENE_MODEL_RESCUE_DIR=" + shlex.quote(env["GG_INPUT_GENE_MODEL_RESCUE_DIR"]) + " ") if args.refinement and not args.rescue and env.get("GG_INPUT_GENE_MODEL_RESCUE_DIR") else ""
        print(refinement_env + rescue_env + anchor_env + "GG_INPUT_TASK_PLAN_OUTPUT=" + shlex.quote(str(plan_path)) + " GG_INPUT_INPUT_GENERATION_MODE=" + mode + " " + shlex.join(cmd), flush=True)
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
        if args.retry and args.submit and pending:
            ensure_no_active_legacy_worker_array()
        worker_id = dispatch("array_worker", ["--array=" + array_expression(pending) + "%" + str(args.max_running)]) if pending else ""
    dispatch("array_finalize", (["--dependency=afterok:" + worker_id] if worker_id else []) + (["--wait"] if args.rescue or args.refinement else []))
    if args.rescue:
        rescue_plan = rescue_output / "plan.json"
        if not rescue_plan.exists():
            if args.submit:
                parser.error("Rescue plan missing after initial finalize: " + str(rescue_plan))
            print("After initial finalize, read " + str(rescue_plan) + " for P comparisons and S species.")
            pair_id = dispatch("rescue_synteny", ["--array=1-P%" + str(args.max_running)])
            species_id = dispatch("rescue_models", ["--dependency=afterok:" + pair_id, "--array=1-S%" + str(args.rescue_max_running or args.max_running)])
        else:
            rescue = json.loads(rescue_plan.read_text())
            plan_hash = digest(rescue_plan)
            def done(directory, key):
                try:
                    receipt = json.loads((directory / "receipt.json").read_text())
                    return isinstance(receipt, dict) and receipt.get("key") == key and isinstance(receipt.get("files"), dict) and bool(receipt["files"]) and all(
                        isinstance(p, str) and isinstance(value, str) and
                        (directory / p).is_file() and digest(directory / p) == value for p, value in receipt["files"].items())
                except (OSError, ValueError):
                    return False
            pending_prepared = {n for n in rescue["species"] if not done(rescue_output / "prepared" / n, {"plan": plan_hash, "species": n})}
            def comparison_done(job):
                try:
                    key = {"plan": plan_hash, "job": job,
                           "prepared": {n: digest(rescue_output / "prepared" / n / "receipt.json")
                                        for n in sorted({job["a"], job["b"]})}}
                    return done(rescue_output / "synteny" / job["id"], key)
                except (OSError, ValueError, KeyError, TypeError):
                    return False
            pairs = [j["index"] for j in rescue["synteny_jobs"] if not args.retry
                     or not comparison_done(j)
                     or j.get("a") in pending_prepared or j.get("b") in pending_prepared]
            def worker_done(name):
                try:
                    donors = rescue["donors"][name]
                    jobs = [j for j in rescue["synteny_jobs"] if name in {j["a"], j["b"]}
                            and (j["a"] == j["b"] or (j["b"] if name == j["a"] else j["a"]) in donors)]
                    key = {"plan": plan_hash, "species": name,
                           "comparisons": {j["id"]: digest(rescue_output / "synteny" / j["id"] / "receipt.json") for j in jobs},
                           "prepared": {n: digest(rescue_output / "prepared" / n / "receipt.json") for n in [name, *donors]}}
                    if not done(rescue_output / "rescued" / name, key):
                        return False
                    effective_key = {"plan": plan_hash, "rescue_receipts": {
                        name: digest(rescue_output / "rescued" / name / "receipt.json")}}
                    return (done(rescue_output / "effective" / name, effective_key)
                            and done(rescue_output / "workers" / name, {"plan": plan_hash, "species": name}))
                except (OSError, ValueError, KeyError, TypeError):
                    return False
            # A repaired comparison can change candidate evidence. Revisit all
            # species after comparison retries; workers verify their own caches.
            species = [i for i, n in enumerate(rescue["species"], 1)
                       if pairs or not args.retry or not worker_done(n)]
            pair_id = dispatch("rescue_synteny", ["--array=" + array_expression(pairs) + "%" + str(args.max_running)]) if pairs else ""
            species_id = dispatch("rescue_models", (["--dependency=afterok:" + pair_id] if pair_id else []) +
                                  ["--array=" + array_expression(species) + "%" + str(args.rescue_max_running or args.max_running)]) if species else ""
        rescue_final_id = dispatch("rescue_finalize", (["--dependency=afterok:" + species_id] if species_id else []) + (["--wait"] if args.refinement else []))
        if args.refinement:
            dispatch("refinement_prepare", ["--dependency=afterok:" + rescue_final_id, "--wait"])

    if args.refinement:
        frozen_path = refinement_output / "plan.json"
        if frozen_path.exists():
            refinement = json.loads(frozen_path.read_text())
            species_count = len(refinement["species"])
            anchor_output = refinement["request"].get("rescue_output")
            pair_count = 0
            if anchor_output and not refinement["request"].get("edges"):
                anchor_plan = json.loads((Path(anchor_output) / "plan.json").read_text())
                pairs = [j["index"] for j in anchor_plan["synteny_jobs"] if j["a"] != j["b"]]
                pair_count = len(pairs)
                env["GG_INPUT_GENE_MODEL_RESCUE_DIR"] = anchor_output
                pair_id = dispatch("rescue_synteny", ["--array=" + array_expression(pairs) + "%" + str(args.max_running)]) if pairs else ""
            else:
                pair_id = ""
            indices = list(range(1, species_count + 1))
            # Receipts are checked again by each idempotent worker. Retrying all
            # indices also revisits dependencies repaired by comparison retries.
            catalog_id = dispatch("refinement_catalog", ["--array=" + array_expression(indices) + "%" + str(args.max_running)])
            dependencies = ":".join(job for job in (pair_id, catalog_id) if job)
            print(json.dumps({"refinement_species": species_count, "refinement_pairs": pair_count}))
        else:
            if args.submit:
                parser.error("Refinement plan missing after preparation: " + str(frozen_path))
            print("After refinement preparation, read " + str(frozen_path) + " for sparse comparisons and species tasks.")
            pair_id = dispatch("rescue_synteny", ["--array=1-P%" + str(args.max_running)])
            catalog_id = dispatch("refinement_catalog", ["--array=1-S%" + str(args.max_running)])
            dependencies = pair_id + ":" + catalog_id
            indices = None
        graph_id = dispatch("refinement_correspondence", ["--dependency=afterok:" + dependencies])
        worker_id = dispatch("refinement_predict", ["--dependency=afterok:" + graph_id, "--array=" + (array_expression(indices) if indices else "1-S") + "%" + str(args.rescue_max_running or args.max_running)])
        dispatch("refinement_finalize", ["--dependency=afterok:" + worker_id])


if __name__ == "__main__":
    main()

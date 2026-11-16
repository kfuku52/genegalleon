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
from input_generation_array_state import FreshDigestBatch, digest, load_plan, prepared, verify_receipt  # noqa: E402


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
                             "gg_input_rescue_synteny,gg_input_rescue_models,gg_input_rescue_predict,"
                             "gg_input_rescue_postprocess,gg_input_rescue_finalize", "--format", "%i"],
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
    parser.add_argument("--rescue-split-miniprot", action="store_true",
                        help="Use separate miniprot and postprocessing jobs with verified prediction checkpoints")
    parser.add_argument("--rescue-predict-cpus", type=int, help="CPUs for the split prediction jobs (default: rescue CPUs)")
    parser.add_argument("--rescue-predict-memory", help="Memory for the split prediction jobs (default: rescue memory)")
    parser.add_argument("--rescue-predict-max-running", type=int,
                        help="Maximum concurrent split prediction jobs (default: rescue concurrency)")
    parser.add_argument("--submit", action="store_true", help="Submit jobs; default prints a dry-run preview")
    args = parser.parse_args()
    if any(value is not None and value < 1 for value in (
            args.cpus, args.max_running, args.prepare_cpus, args.rescue_cpus, args.rescue_max_running,
            args.rescue_predict_cpus, args.rescue_predict_max_running)):
        parser.error("CPU and concurrency counts must be positive")
    if args.rescue_split_miniprot and not args.rescue:
        parser.error("--rescue-split-miniprot requires --rescue")
    if not args.rescue_split_miniprot and any(value is not None for value in (
            args.rescue_predict_cpus, args.rescue_predict_memory, args.rescue_predict_max_running)):
        parser.error("Prediction resources require --rescue-split-miniprot")
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
        model_worker = mode in {"rescue_models", "rescue_predict", "rescue_postprocess", "rescue_finalize",
                                "refinement_predict", "refinement_finalize"}
        cpus = (args.prepare_cpus or args.cpus) if preparing else args.cpus
        memory = (args.prepare_memory or args.memory) if preparing else args.memory
        if model_worker:
            cpus = args.rescue_cpus or cpus
            memory = args.rescue_memory or memory
        if mode == "rescue_predict":
            cpus = args.rescue_predict_cpus or cpus
            memory = args.rescue_predict_memory or memory
        partition = (args.prepare_partition if args.prepare_partition is not None else args.partition) if preparing else args.partition
        base = ["sbatch", "--parsable", "--cpus-per-task=" + str(cpus), "--mem=" + memory, "--time=" + args.time]
        if partition:
            base += ["--partition=" + partition]
        # Values are inherited via the process environment, avoiding comma/space
        # escaping problems in Slurm's --export parser.
        return base + ["--export=ALL", "--job-name=gg_input_" + mode] + extra + ["--wrap=" + shlex.join(["exec", "bash", str(Path(args.entrypoint).resolve())])]

    def dispatch(mode, extra):
        if (args.rescue_split_miniprot and mode in {"rescue_predict", "rescue_postprocess", "rescue_finalize"}
                and any(argument.startswith("--dependency=") for argument in extra)):
            extra = [*extra, "--kill-on-invalid-dep=yes"]
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

    def dispatch_rescue_species(species, prediction, pair_id, ready=""):
        pair_dependency = ["--dependency=afterok:" + pair_id] if pair_id else []
        limit = args.rescue_max_running or args.max_running
        if not args.rescue_split_miniprot:
            return dispatch("rescue_models", pair_dependency + [f"--array={species}%{limit}"]) if species else ""
        post_jobs = []
        ready_limit = limit if not prediction else max(1, limit // 2)
        if ready:
            post_jobs.append(dispatch("rescue_postprocess", [f"--array={ready}%{ready_limit}"]))
        if prediction:
            prediction_limit = args.rescue_predict_max_running or limit
            prediction_id = dispatch("rescue_predict", pair_dependency + [f"--array={prediction}%{prediction_limit}"])
            dependency = "aftercorr:" + prediction_id
            post_limit = limit if not ready else max(1, limit - ready_limit)
            if ready and limit == 1:
                dependency += ",afterany:" + post_jobs[0]
            post_jobs.append(dispatch("rescue_postprocess", ["--dependency=" + dependency,
                                                             f"--array={prediction}%{post_limit}"]))
        return ":".join(post_jobs)

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
            species_id = dispatch_rescue_species("1-S", "1-S", pair_id)
        else:
            rescue = json.loads(rescue_plan.read_text())
            plan_hash = digest(rescue_plan)
            def done(directory, key, proof=None):
                try:
                    receipt = json.loads((directory / "receipt.json").read_text())
                    if not (isinstance(receipt, dict) and receipt.get("key") == key
                            and isinstance(receipt.get("files"), dict) and receipt["files"]
                            and all(isinstance(p, str) and isinstance(value, str)
                                    for p, value in receipt["files"].items())):
                        return False
                    if proof is None:
                        return all((directory / p).is_file() and digest(directory / p) == value
                                   for p, value in receipt["files"].items())
                    expected = {str(directory / p): value for p, value in receipt["files"].items()}
                    proof.read(expected, expected=expected)
                    return True
                except (OSError, ValueError):
                    return False
            pending_prepared = {n for n in rescue["species"] if args.retry
                                and not done(rescue_output / "prepared" / n, {"plan": plan_hash, "species": n})}
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
            def species_jobs(name):
                donors = rescue["donors"][name]
                return [j for j in rescue["synteny_jobs"] if name in {j["a"], j["b"]}
                        and (j["a"] == j["b"] or (j["b"] if name == j["a"] else j["a"]) in donors)]

            def species_key(name):
                donors = rescue["donors"][name]
                return {"plan": plan_hash, "species": name,
                        "comparisons": {j["id"]: digest(rescue_output / "synteny" / j["id"] / "receipt.json")
                                        for j in species_jobs(name)},
                        "prepared": {n: digest(rescue_output / "prepared" / n / "receipt.json") for n in [name, *donors]}}

            def prediction_done(name):
                try:
                    return done(rescue_output / "miniprot_predictions" / name, species_key(name))
                except (OSError, ValueError, KeyError, TypeError):
                    return False

            def worker_done(name):
                try:
                    key = species_key(name)
                    proof = FreshDigestBatch()
                    if not done(rescue_output / "rescued" / name, key, proof):
                        return False
                    effective_key = {"plan": plan_hash, "rescue_receipts": {
                        name: digest(rescue_output / "rescued" / name / "receipt.json")}}
                    return (done(rescue_output / "effective" / name, effective_key, proof)
                            and done(rescue_output / "workers" / name, {"plan": plan_hash, "species": name}, proof))
                except (OSError, ValueError, KeyError, TypeError):
                    return False
            species, predictions, ready = [], [], []
            pending_pairs = set(pairs)
            for i, name in enumerate(rescue["species"], 1):
                needs_dependencies = bool(pairs)
                if args.retry and args.rescue_split_miniprot:
                    try:
                        needs_dependencies = (bool(pending_prepared.intersection([name, *rescue["donors"][name]]))
                                              or any(job["index"] in pending_pairs for job in species_jobs(name)))
                    except (KeyError, TypeError):
                        needs_dependencies = True
                if args.retry and not needs_dependencies and worker_done(name):
                    continue
                species.append(i)
                if args.rescue_split_miniprot:
                    target = (predictions if needs_dependencies or not args.retry or not prediction_done(name) else ready)
                    target.append(i)
            pair_id = dispatch("rescue_synteny", ["--array=" + array_expression(pairs) + "%" + str(args.max_running)]) if pairs else ""
            species_id = dispatch_rescue_species(array_expression(species), array_expression(predictions), pair_id,
                                                 array_expression(ready))
        rescue_final_id = dispatch("rescue_finalize", (["--dependency=afterok:" + species_id] if species_id else []) + (["--wait"] if args.refinement else []))
        if args.refinement:
            dispatch("refinement_prepare", ["--dependency=afterok:" + rescue_final_id, "--wait"])

    if args.refinement:
        frozen_path = refinement_output / "plan.json"
        if frozen_path.exists():
            refinement = json.loads(frozen_path.read_text())
            targets = refinement["request"].get("target_species")
            species_count = len(targets if targets is not None else refinement["species"])
            anchor_output = refinement["request"].get("rescue_output")
            pair_count = 0
            if anchor_output and targets is None and (not refinement["request"].get("edges") or refinement["request"].get("transcriptome_sources")):
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
            # Correspondence verifies/builds the shared donor catalogs and the
            # target catalog. Keep donor indices out of fixed-panel arrays.
            catalog_id = (dispatch("refinement_catalog", ["--array=" + array_expression(indices) + "%" + str(args.max_running)])
                          if targets is None else "")
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
        graph_id = dispatch("refinement_correspondence", ["--dependency=afterok:" + dependencies] if dependencies else [])
        worker_id = dispatch("refinement_predict", ["--dependency=afterok:" + graph_id, "--array=" + (array_expression(indices) if indices else "1-S") + "%" + str(args.rescue_max_running or args.max_running)])
        dispatch("refinement_finalize", ["--dependency=afterok:" + worker_id])


if __name__ == "__main__":
    main()

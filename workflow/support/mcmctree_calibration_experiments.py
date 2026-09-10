#!/usr/bin/env python3
"""Run isolated, provenance-tracked prior/posterior calibration sensitivity experiments."""

from __future__ import annotations

import argparse
import csv
import importlib.metadata
import json
import math
import re
import shutil
import subprocess
import time
from decimal import Decimal
from pathlib import Path

import numpy as np
import pandas as pd
from mcmctree_calibration_audit import (
    atomic_text,
    digest,
    json_text,
    load_tree,
    node_id,
    reviewed_rows,
    tree_text,
)
from mcmctree_time_scale import scale_calibration_labels, scale_ctl_rootage_text


def file_hash(path):
    import hashlib

    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def control_values(text):
    values = {}
    for line in text.splitlines():
        line = re.split(r"[\*#]", line, maxsplit=1)[0].strip()
        if not line:
            continue
        match = re.fullmatch(r"([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.+)", line)
        if not match or match[1] in values:
            raise ValueError(f"Malformed or duplicate control setting: {line}")
        values[match[1]] = match[2].strip()
    return values


def scenarios(tree, manifest=None):
    from nwkit.time_tree import parse_mcmctree_calibration

    labels = {node_id(node): node.name for node in tree.traverse()
              if not node.is_leaf and parse_mcmctree_calibration(node.name)}
    if not labels:
        raise ValueError("No calibrations in the supplied tree.")
    result = [{"id": "baseline", "omit_nodes": [], "reason": "full supplied calibration set"}]
    for key in sorted(labels):
        result.append({"id": f"omit_node_{key}", "omit_nodes": [key], "reason": "leave one calibration out"})
    if manifest:
        _, accepted = reviewed_rows(manifest, tree)

        def normalize(text):
            return {k: v for k, v in parse_mcmctree_calibration(text).items() if k != "raw"}
        expected = {row["node_id"]: normalize(row["calibration"]) for row in accepted}
        if expected != {key: normalize(label) for key, label in labels.items()}:
            raise ValueError("Reviewed manifest does not describe the supplied calibration tree.")
        groups = {}
        for row in accepted:
            for group in json.loads(row["dependency_groups"]):
                groups.setdefault(group, []).append(row["node_id"])
        for group, keys in sorted(groups.items()):
            result.append({"id": "omit_group_" + digest(group.encode()),
                           "omit_nodes": sorted(keys), "reason": "leave shared-evidence group out",
                           "group": group})
    for scenario in result:
        # Do not silently replace the last absolute-age anchor with IQ2MC's
        # generic RootAge bound. Such experiments require a separately justified prior.
        scenario["status"] = ("skipped_no_remaining_calibration"
                              if len(scenario["omit_nodes"]) == len(labels) else "pending")
    return result


def write_tsv(path, rows):
    if not rows:
        return
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_chain(path, tree, expected_draws, sampfreq):
    from nwkit.time_tree import read_mcmctree_posterior

    # Delegate PAML node mapping, duplicate headers and age ordering to NWKIT.
    read_mcmctree_posterior(str(path), tree)
    table = pd.read_csv(path, sep=r"\s+")
    values = table.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("MCMC samples contain nonfinite values.")
    scheduled = np.arange(1, expected_draws + 1, dtype=np.int64) * sampfreq
    # PAML also writes iteration 1 before the regular sampling grid (see
    # mcmctree.c: ir == 0 || (ir + 1) % sampfreq == 0). Retain that raw row,
    # but use equally spaced draws for autocorrelation-based diagnostics.
    expected = np.unique(np.concatenate(([1], scheduled)))
    if "Gen" not in table or not np.array_equal(table["Gen"].to_numpy(), expected):
        raise ValueError("Incomplete or irregular MCMC generation grid.")
    return table.loc[table["Gen"].isin(scheduled)].copy()


def summarize(target_dir, tree, chain_paths, scale):
    from nwkit.time_tree import paml_node_mapping, parse_mcmctree_calibration

    raw = target_dir / "diagnostics_internal.tsv"
    script = Path(__file__).with_name("mcmctree_calibration_diagnostics.R")
    with (target_dir / "diagnostics.log").open("w") as log:
        subprocess.run(["Rscript", str(script), str(raw), *map(str, chain_paths)],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    table = pd.read_csv(raw, sep="\t")
    mapping = {f"t_n{key}": node for key, node in paml_node_mapping(tree).items() if not node.is_leaf}
    samples = pd.concat([pd.read_csv(path, sep=r"\s+") for path in chain_paths], ignore_index=True)
    rows = []
    for row in table.to_dict("records"):
        parameter = row["parameter"]
        node = mapping.get(parameter)
        row.update(node_id=node_id(node) if node is not None else "",
                   unit="Ma" if node is not None else "internal_control_unit",
                   mass_below_bound="", mass_above_bound="")
        if node is not None:
            for field in ("mean", "median", "eti_low", "eti_high", "sd", "mcse_mean"):
                row[field] *= float(scale)
            calibration = parse_mcmctree_calibration(node.name)
            if calibration:
                values = samples[parameter].to_numpy() * float(scale)
                if "lower" in calibration:
                    row["mass_below_bound"] = float(np.mean(values < calibration["lower"]))
                if "upper" in calibration:
                    row["mass_above_bound"] = float(np.mean(values > calibration["upper"]))
        # lnL is constant by construction in prior-only runs. Do not exempt
        # constant ages or rate parameters: their diagnostics remain undefined.
        fixed_likelihood = parameter == "lnL" and target_dir.name == "prior"
        metrics = [row[key] for key in ("rhat", "ess_bulk", "ess_tail", "mcse_mean", "sd")]
        passed = (all(math.isfinite(value) for value in metrics)
                  and row["rhat"] <= 1.01 and min(row["ess_bulk"], row["ess_tail"]) >= 400
                  and row["sd"] > 0 and row["mcse_mean"] / row["sd"] <= 0.05)
        row["diagnostic_status"] = "not_applicable" if fixed_likelihood else ("pass" if passed else "not_converged")
        rows.append(row)
    write_tsv(target_dir / "summary.tsv", rows)
    correlation = samples[list(mapping)].corr()
    correlation.index = [node_id(mapping[key]) for key in correlation.index]
    correlation.columns = [node_id(mapping[key]) for key in correlation.columns]
    correlation.to_csv(target_dir / "node_correlations.tsv", sep="\t", index_label="node_id")
    return rows


def compare_results(run_dir, all_rows):
    by_key = {(row["scenario"], row["target"], row["node_id"]): row
              for row in all_rows if row["node_id"]}
    comparisons = []
    for row in all_rows:
        if not row["node_id"] or row["target"] != "posterior":
            continue
        for comparison, key in (
            ("posterior_minus_prior", (row["scenario"], "prior", row["node_id"])),
            ("posterior_minus_baseline", ("baseline", "posterior", row["node_id"])),
        ):
            reference = by_key.get(key)
            if reference is None:
                continue
            width = row["eti_high"] - row["eti_low"]
            reference_width = reference["eti_high"] - reference["eti_low"]
            comparisons.append({"scenario": row["scenario"], "node_id": row["node_id"],
                                "comparison": comparison, "median_difference_Ma": row["median"] - reference["median"],
                                "eti_width_ratio": width / reference_width if reference_width > 0 else "NA",
                                "diagnostic_status": "pass" if row["diagnostic_status"] == reference["diagnostic_status"] == "pass" else "not_converged",
                                "scientific_robustness": "requires_prespecified_biological_tolerance"})
    write_tsv(run_dir / "comparisons.tsv", comparisons)


def run_experiments(args):
    if args.chains < 2 or args.seed < 1 or args.seed + args.chains * 100000 >= 2147483647:
        raise ValueError("Use at least two chains and a positive seed safely below 2^31.")
    scale = Decimal(args.time_scale)
    if not scale.is_finite() or scale <= 0:
        raise ValueError("Time scale must be finite and positive.")
    control_text = args.control.read_text(encoding="utf-8")
    values = control_values(scale_ctl_rootage_text(control_text, scale, "down"))
    if values.get("usedata", "").split()[0:1] != ["2"]:
        raise ValueError("Experiments require IQ2MC approximate-likelihood inputs (usedata=2).")
    for key in ("burnin", "sampfreq", "nsample"):
        if not re.fullmatch(r"\d+", values.get(key, "")):
            raise ValueError(f"Missing or invalid {key}.")
    if int(values["nsample"]) < 8 or int(values["sampfreq"]) < 1:
        raise ValueError("Use at least 8 samples and a positive sampling frequency.")
    tree = load_tree(args.tree)
    planned = scenarios(tree, args.manifest)
    if args.seed + len(planned) * 2 * args.chains >= 2147483647:
        raise ValueError("Too many scenarios for the requested seed.")
    sources = {"control": args.control, "tree": args.tree, "alignment": args.alignment, "hessian": args.hessian}
    if args.manifest:
        sources["manifest"] = args.manifest
    executable = shutil.which("mcmctree")
    if not executable:
        raise ValueError("MCMCtree is unavailable; use the GeneGalleon container.")
    r_version = subprocess.run(["Rscript", "-e", "stopifnot(requireNamespace('posterior', quietly=TRUE)); cat(as.character(packageVersion('posterior')))"],
                               capture_output=True, text=True, check=True).stdout.strip()
    contract = {"schema_version": 1, "inputs": {key: file_hash(path) for key, path in sources.items()},
                "chains": args.chains, "seed": args.seed, "time_scale": str(scale),
                "mcmctree_sha256": file_hash(executable), "nwkit_version": importlib.metadata.version("nwkit"),
                "posterior_version": r_version, "python_adapter_sha256": file_hash(__file__),
                "audit_adapter_sha256": file_hash(Path(__file__).with_name("mcmctree_calibration_audit.py")),
                "time_scale_adapter_sha256": file_hash(Path(__file__).with_name("mcmctree_time_scale.py")),
                "r_adapter_sha256": file_hash(Path(__file__).with_name("mcmctree_calibration_diagnostics.R"))}
    run_dir = args.outdir.resolve() / digest(json_text(contract).encode())
    status_path = run_dir / "status.json"
    if run_dir.exists():
        if not status_path.is_file():
            raise ValueError(f"Incomplete or active experiment directory, retained for inspection: {run_dir}")
        status = json.loads(status_path.read_text())
        if status.get("status") not in {"diagnostics_pass", "not_converged"}:
            raise ValueError(f"Previous experiment is incomplete or failed; choose a new --outdir: {run_dir}")
        for relative, expected in status["output_sha256"].items():
            path = run_dir / relative
            if not path.is_file() or file_hash(path) != expected:
                raise ValueError(f"Experiment output changed or missing: {path}")
        print(f"Reusing intact calibration diagnostics: {run_dir} ({status['status']})")
        return run_dir
    run_dir.mkdir(parents=True, exist_ok=False)
    status = {"status": "running", "contract": contract, "scenarios": planned,
              "main_dated_tree": "not_modified", "scientific_validation": "not_established",
              "root_age_control_internal": values.get("RootAge", "not_specified"),
              "group_sensitivity": "reviewed_groups" if args.manifest else "unavailable_without_manifest",
              "interval": "equal_tail_95_percent", "time_unit": "Ma for node ages; internal units otherwise"}
    atomic_text(status_path, json_text(status) + "\n")
    all_rows = []
    try:
        inputs = run_dir / "inputs"
        inputs.mkdir()
        for key, source in sources.items():
            shutil.copyfile(source, inputs / key)
            if file_hash(inputs / key) != contract["inputs"][key]:
                raise ValueError(f"Input changed while creating experiment snapshot: {source}")
        for scenario_index, scenario in enumerate(planned):
            if scenario["status"].startswith("skipped"):
                continue
            scenario_tree = load_tree(inputs / "tree")
            for node in scenario_tree.traverse():
                if not node.is_leaf and node_id(node) in scenario["omit_nodes"]:
                    node.name = ""
            public_tree = tree_text(scenario_tree)
            scaled_tree = scale_calibration_labels(public_tree, scale, "down")
            for target_index, target in enumerate(("prior", "posterior")):
                target_dir = run_dir / scenario["id"] / target
                chain_paths = []
                for chain_index in range(args.chains):
                    chain_dir = target_dir / f"chain_{chain_index + 1}"
                    chain_dir.mkdir(parents=True)
                    chain_values = dict(values)
                    chain_values.update(seed=str(args.seed + (scenario_index * 2 + target_index) * args.chains + chain_index),
                                        seqfile="../../../inputs/alignment", treefile="tree.nwk",
                                        outfile="mcmctree.out", mcmcfile="mcmc.txt", print="1",
                                        hessianfile="../../../inputs/hessian", ckpfile="mcmctree.ckp", checkpoint="0",
                                        usedata="0" if target == "prior" else "2 ../../../inputs/hessian")
                    (chain_dir / "tree.nwk").write_text(f"{len(scenario_tree)} 1\n" + scaled_tree, encoding="utf-8")
                    (chain_dir / "mcmctree.ctl").write_text("\n".join(f"{key} = {value}" for key, value in chain_values.items()) + "\n", encoding="utf-8")
                    started = time.monotonic()
                    with (chain_dir / "stdout.log").open("w") as out, (chain_dir / "stderr.log").open("w") as err:
                        subprocess.run([executable, "mcmctree.ctl"], cwd=chain_dir, stdout=out, stderr=err, check=True)
                    # Use the scaled topology for validating samples; labels do not affect node mapping.
                    path = chain_dir / "mcmc.txt"
                    samples = read_chain(path, load_tree(chain_dir / "tree.nwk"),
                                         int(values["nsample"]), int(values["sampfreq"]))
                    diagnostic_path = chain_dir / "diagnostic_samples.tsv"
                    samples.to_csv(diagnostic_path, sep="\t", index=False)
                    atomic_text(chain_dir / "timing.json", json_text({"elapsed_seconds": time.monotonic() - started}) + "\n")
                    chain_paths.append(diagnostic_path)
                    print(f"Calibration diagnostic {scenario['id']} {target} chain {chain_index + 1}/{args.chains} complete", flush=True)
                rows = summarize(target_dir, scenario_tree, chain_paths, scale)
                for row in rows:
                    row.update(scenario=scenario["id"], target=target)
                all_rows.extend(rows)
            scenario["status"] = "complete"
            atomic_text(status_path, json_text(status) + "\n")
        write_tsv(run_dir / "summary.tsv", all_rows)
        compare_results(run_dir, all_rows)
        status["status"] = "diagnostics_pass" if all(row["diagnostic_status"] in {"pass", "not_applicable"} for row in all_rows) else "not_converged"
        status["output_sha256"] = {str(path.relative_to(run_dir)): file_hash(path)
                                   for path in sorted(run_dir.rglob("*")) if path.is_file() and path != status_path}
    except BaseException as exc:
        status.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        atomic_text(status_path, json_text(status) + "\n")
    print(f"Calibration diagnostics: {run_dir} ({status['status']}); main dated tree unchanged")
    return run_dir


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("control", "tree", "alignment", "hessian", "outdir"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--chains", type=int, default=4)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--time-scale", default="1")
    run_experiments(parser.parse_args())


if __name__ == "__main__":
    main()

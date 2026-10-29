#!/usr/bin/env python3
"""Delegate profiles to kfFractBias and compare paired, disjoint reference blocks."""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import os
import tempfile
from pathlib import Path

import numpy as np

if __package__:
    from .subgenome_dominance import adjust_p, inference, validate_resampling
else:
    from subgenome_dominance import adjust_p, inference, validate_resampling


def paired_contrasts(loci, *, block_size=100, replicates=2000, permutations=100000, seed=1):
    if type(block_size) is not int or block_size < 1 or type(seed) is not int or seed < 0:
        raise ValueError("Require a positive block size and a nonnegative integer seed")
    profiles = {}
    for r in loci:
        key = (r["target_seqid"], r["target_gene"], int(r["target_rank"]))
        profile = profiles.setdefault(r["analysis_id"], {})
        if key in profile or r["common_callable"] not in {"0", "1"}:
            raise ValueError("Duplicate locus or invalid callable value")
        count, depth = int(r["retained_copies"]), int(r["expected_copies"])
        if count < 0 or depth < 1 or count > depth:
            raise ValueError("Invalid copy count for masked comparison")
        profile[key] = (r["common_callable"] == "1", count / depth)
    keys = next(iter(profiles.values())).keys() if profiles else []
    if len(profiles) < 2 or any(p.keys() != keys for p in profiles.values()):
        raise ValueError("Require at least two profiles on identical reference loci")
    if any(len({p[k][0] for p in profiles.values()}) != 1 for k in keys):
        raise ValueError("Shared callable mask differs among profiles")
    selected = [k for k in keys if next(iter(profiles.values()))[k][0]]
    blocks = [(k[0], (k[2] - 1) // block_size) for k in selected]
    rows = []
    for a, b in itertools.combinations(profiles, 2):
        stream = int.from_bytes(hashlib.sha256(json.dumps([seed, a, b]).encode()).digest()[:16], "big")
        result = inference([profiles[a][k][1] - profiles[b][k][1] for k in selected], blocks,
                           replicates, np.random.default_rng(stream),
                           permutation_replicates=permutations)
        rows.append({"analysis_a": a, "analysis_b": b, "block_size": block_size, **result})
    adjust_p(rows)
    return rows


def run(config_path, output_dir):
    from kffractbias.comparison import comparison_report

    config_path = Path(config_path).resolve()
    raw_config = config_path.read_bytes()
    config = json.loads(raw_config)
    allowed = {"manifest", "prefix", "window_size", "step_size", "denominator",
               "minimum_callable_fraction", "minimum_anchors", "trim", "raw", "target_seqids",
               "plot_config", "statistics", "comparison_scope"}
    if not isinstance(config, dict) or "manifest" not in config or set(config) - allowed:
        raise ValueError("Require a report configuration object with supported keys and manifest")
    manifest = (config_path.parent / config.pop("manifest")).resolve()
    statistics = config.pop("statistics", None)
    if statistics is not None and config.get("comparison_scope", "shared") == "pairwise":
        raise ValueError("Paired block inference requires comparison_scope=shared")
    if statistics is not None:
        if not isinstance(statistics, dict) or set(statistics) - {"block_size", "replicates", "permutations", "seed"}:
            raise ValueError("Unsupported block statistics configuration")
        validate_resampling(statistics.get("replicates", 2000), statistics.get("permutations", 100000),
                            262144, statistics.get("seed", 1))
        block_size = statistics.get("block_size", 100)
        if type(block_size) is not int or block_size < 1:
            raise ValueError("Require a positive integer block size")
    if statistics is not None and config.get("raw", False):
        raise ValueError("Block inference requires masked local-copy profiles")
    prefix = config.get("prefix", "comparison")
    output_dir = Path(output_dir)
    if config.get("comparison_scope", "shared") == "pairwise" and (output_dir / (prefix + ".contrasts.json")).exists():
        raise ValueError("Use a new prefix for pairwise reports: previous shared contrasts exist")
    summary = comparison_report(manifest, output_dir, **config)
    if config_path.read_bytes() != raw_config:
        raise RuntimeError("Report configuration changed during calculation")
    if statistics is not None:
        loci_path = output_dir / (prefix + ".loci.tsv")
        raw_loci = loci_path.read_bytes()
        if hashlib.sha256(raw_loci).hexdigest() != summary["output_sha256"]["loci"]:
            raise RuntimeError("Closed profile table changed before inference")
        with loci_path.open() as stream:
            rows = paired_contrasts(list(csv.DictReader(stream, delimiter="\t")), **statistics)
        if config_path.read_bytes() != raw_config or loci_path.read_bytes() != raw_loci:
            raise RuntimeError("Report inputs changed during inference")
        result = {"method": "Locus-weighted paired differences; non-overlapping reference blocks; cluster bootstrap and block sign flips; BH across all species pairs",
                  "effect_units": "fraction of expected copies (analysis_a minus analysis_b)",
                  "config_sha256": hashlib.sha256(raw_config).hexdigest(),
                  "loci_sha256": hashlib.sha256(raw_loci).hexdigest(), "parameters": statistics, "contrasts": rows}
        path = output_dir / (prefix + ".contrasts.json")
        descriptor, name = tempfile.mkstemp(prefix="." + prefix + ".contrasts-", dir=output_dir)
        temporary = Path(name)
        try:
            with os.fdopen(descriptor, "w") as stream:
                json.dump(result, stream, indent=2, sort_keys=True)
                stream.write("\n")
            temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    run(args.config, args.output_dir)


if __name__ == "__main__":
    main()

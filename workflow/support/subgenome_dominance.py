#!/usr/bin/env python3
"""Analyse independently assigned homoeologs, with explicit callable retention loci.

Neither expression nor retained-gene counts are used to assign subgenomes.
Local group labels are never pooled into a genome-wide subgenome identity.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

if __package__:
    from . import subgenome_dominance_plot as plotting
else:
    import subgenome_dominance_plot as plotting

INFERENCE_VERSION = 3
DEFAULT_PERMUTATIONS = 100000
DEFAULT_EXACT_STATES = 262144


def contrast_seed(seed, analysis_id, metric, group, a, b, tissue, purpose):
    """Stable streams independent of other contrasts and bootstrap draw counts."""
    key = ["subgenome-contrast-v1", seed, analysis_id, metric, group, a, b, tissue, purpose]
    return int.from_bytes(hashlib.sha256(json.dumps(key, ensure_ascii=True).encode()).digest()[:16], "big")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def table(path, required):
    with Path(path).open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f"Missing or duplicate column headers: {path}")
        missing = set(required) - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing columns {sorted(missing)}: {path}")
        rows = list(reader)
    if any(None in row or any(v is None for v in row.values()) or
           any(not row[k].strip() for k in required) for row in rows):
        raise ValueError(f"Incomplete required cells or malformed rows: {path}")
    return rows


def write_table(path, rows, fields):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def validate_resampling(replicates, permutation_replicates, exact_max_states, seed=0):
    for name, value, minimum in (("replicates", replicates, 100),
                                 ("permutation_replicates", permutation_replicates, 100),
                                 ("exact_max_states", exact_max_states, 1), ("seed", seed, 0)):
        if type(value) is not int or value < minimum:
            raise ValueError(f"{name} must be an integer >= {minimum}")


def make_plan(workspace, manifest, replicates=2000, seed=1,
              permutation_replicates=DEFAULT_PERMUTATIONS, exact_max_states=DEFAULT_EXACT_STATES,
              plot_config=None):
    workspace = Path(workspace).resolve()
    manifest = Path(manifest).resolve()
    validate_resampling(replicates, permutation_replicates, exact_max_states, seed)
    analyses = table(manifest, ("analysis_id", "species", "mapping_file"))
    if not analyses or len({r["analysis_id"] for r in analyses}) != len(analyses):
        raise ValueError("Require unique, nonempty analysis IDs")
    inputs = {str(manifest): digest(manifest)}
    for row in analyses:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", row["analysis_id"]):
            raise ValueError("Unsafe analysis_id")
        for field in ("mapping_file", "retention_file", "homoeolog_file", "expression_file", "samples_file"):
            value = row.get(field, "").strip()
            if value:
                path = Path(value)
                if not path.is_absolute():
                    path = manifest.parent / path
                row[field] = str(path.resolve())
                inputs[row[field]] = digest(row[field])
            else:
                row[field] = ""
        expression_inputs = [bool(row[field]) for field in ("homoeolog_file", "expression_file", "samples_file")]
        row["expression_unit"] = row.get("expression_unit", "").strip() or "TPM"
        if row["expression_unit"] not in {"TPM", "FPKM"}:
            raise ValueError("expression_unit must be TPM or FPKM; counts and transformed values are unsupported")
        if any(expression_inputs) and not all(expression_inputs):
            raise ValueError("Expression requires homoeolog_file, expression_file and samples_file together")
        for flag in ("mapping_validated", "retention_validated"):
            if row.get(flag, "0") not in {"0", "1", ""}:
                raise ValueError(f"{flag} must be 0 or 1")
    config, plot_inputs = plotting.load_config(plot_config)
    inputs.update(plot_inputs)
    return {"schema_version": 2, "workspace": str(workspace), "analyses": analyses,
            "replicates": replicates, "seed": seed, "input_hashes": inputs,
            "permutation_replicates": permutation_replicates, "exact_max_states": exact_max_states,
            "plot_config": config, "inference_version": INFERENCE_VERSION,
            "implementation_hashes": {str(path): digest(path) for path in
                                      (Path(__file__).resolve(), Path(plotting.__file__).resolve())},
            "numpy_version": np.__version__}


def verify_inputs(plan):
    for path, expected in {**plan["input_hashes"], **plan.get("implementation_hashes", {})}.items():
        if digest(path) != expected:
            raise ValueError(f"Input changed during subgenome analysis: {path}")


def contract_args(plan):
    workspace = Path(plan["workspace"])
    args = ["--manifest", str(workspace / "output/artifact_provenance/genome_evolution/subgenome_dominance.json"),
            "--step", "subgenome_dominance", "--family-id", "all_analyses",
            "--logical-root", str(workspace / "output/.gg_global_artifacts"), "--workspace-root", str(workspace),
            "--output", f"results={workspace / 'output/genome_evolution/subgenome_dominance'}",
            "--input", f"implementation={Path(__file__).resolve()}",
            "--input", f"plot_implementation={Path(__file__).with_name('subgenome_dominance_plot.py').resolve()}"]
    for index, path in enumerate(sorted(plan["input_hashes"])):
        args.extend(("--input", f"source{index}={path}"))
    args.extend(("--parameter", "configuration=" + json.dumps(plan, sort_keys=True)))
    return args


def load_mapping(path):
    rows = table(path, ("gene_id", "group_id", "subgenome", "assignment_scope", "assignment_basis", "evidence"))
    genes, groups = {}, defaultdict(set)
    scopes = set()
    for row in rows:
        if row["group_id"] == "ALL_GROUPS":
            raise ValueError("ALL_GROUPS is reserved for independently phased global contrasts")
        if row["gene_id"] in genes:
            raise ValueError(f"Duplicate mapped gene: {row['gene_id']}")
        if row["assignment_scope"] not in {"local", "global"}:
            raise ValueError("assignment_scope must be local or global")
        if row["assignment_basis"] not in {"synteny", "phylogeny", "kmer", "curated"}:
            raise ValueError("Assignments must be independent of expression and retention")
        genes[row["gene_id"]] = row
        groups[row["group_id"]].add(row["subgenome"])
        scopes.add(row["assignment_scope"])
    if not genes or any(len(sg) < 2 for sg in groups.values()) or len(scopes) != 1:
        raise ValueError("Require at least two subgenomes per group and one consistent assignment scope")
    if scopes == {"global"} and len({tuple(sorted(v)) for v in groups.values()}) != 1:
        raise ValueError("Global assignments require the same subgenome identities in every group")
    return genes, groups, scopes.pop()


def sign_sums(values):
    sums = np.zeros(1)
    for value in values:
        sums = np.concatenate((sums - value, sums + value))
    return sums


def sign_flip(sums, rng, replicates, exact_max_states):
    """The same two-sided block null, evaluated exactly when bounded."""
    validate_resampling(100, replicates, exact_max_states)
    sums = np.asarray(sums, dtype=float)
    if sums.ndim != 1 or not np.isfinite(sums).all():
        raise ValueError("Sign-flip test requires finite one-dimensional block sums")
    n = len(sums)
    # Integral block sums allow exact dynamic programming without 2**n draws.
    if all(float(s).is_integer() for s in sums):
        weights = [abs(int(s)) for s in sums if s]
        divisor = math.gcd(*weights) if weights else 1
        weights = [w // divisor for w in weights]
        if sum(weights) + 1 <= exact_max_states:
            observed = abs(sum(int(s) for s in sums))
            counts = {0: 1}
            for weight in weights:
                updated = counts.copy()
                for subtotal, count in counts.items():
                    updated[subtotal + weight] = updated.get(subtotal + weight, 0) + count
                counts = updated
            total = sum(weights)
            extreme = sum(count for subtotal, count in counts.items()
                          if abs(2 * subtotal - total) * divisor >= observed)
            return {"p_value": extreme / (1 << len(weights)), "test_method": "exact_integer_dp",
                    "null_draws": 1 << n}
    # Scale before summation to avoid overflow; relative tolerance preserves
    # near-ties without turning all small effects into null statistics of zero.
    scale = max(abs(sums), default=0)
    if scale:
        sums = sums / scale
    observed = abs(math.fsum(sums))
    threshold = max(0, observed * (1 - 100 * np.finfo(float).eps))
    if n <= 16:
        null = abs(sign_sums(sums))
        return {"p_value": float(np.mean(null >= threshold)), "test_method": "exact_enumeration",
                "null_draws": len(null)}
    half = n // 2
    if (1 << (n - half)) <= exact_max_states:
        left, right = sign_sums(sums[:half]), np.sort(sign_sums(sums[half:]))
        if threshold == 0:
            extreme = 1 << n
        else:
            extreme = int(np.sum(np.searchsorted(right, -threshold - left, side="right"), dtype=np.int64))
            extreme += int(np.sum(len(right) - np.searchsorted(right, threshold - left, side="left"), dtype=np.int64))
        return {"p_value": extreme / (1 << n), "test_method": "exact_meet_in_middle", "null_draws": 1 << n}
    extreme = 0
    batch_size = max(1, min(4096, 1000000 // n))
    for start in range(0, replicates, batch_size):
        null = abs((rng.choice([-1, 1], (min(batch_size, replicates - start), n)) * sums).sum(axis=1))
        extreme += int(np.count_nonzero(null >= threshold))
    # Wilson interval for the sampled null tail probability, not biological uncertainty.
    tail, z = extreme / replicates, 1.959963984540054
    denominator = 1 + z * z / replicates
    centre = (tail + z * z / (2 * replicates)) / denominator
    radius = z * math.sqrt(tail * (1 - tail) / replicates + z * z / (4 * replicates ** 2)) / denominator
    return {"p_value": (1 + extreme) / (replicates + 1), "test_method": "monte_carlo",
            "null_draws": replicates, "mc_p_ci_low": max(0, centre - radius),
            "mc_p_ci_high": min(1, centre + radius)}


def inference(values, blocks, replicates, rng, *, null_rng=None,
              permutation_replicates=DEFAULT_PERMUTATIONS, exact_max_states=DEFAULT_EXACT_STATES):
    """Locus-weighted mean, cluster bootstrap CI and block sign-flip null test.

    Input blocks must partition loci, and must be defined independently of bias.
    The CI describes loci conditional on the sampled tissues and annotation.
    """
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not np.isfinite(values).all():
        raise ValueError("Inference requires finite one-dimensional values")
    validate_resampling(replicates, permutation_replicates, exact_max_states)
    if null_rng is None:
        state = json.dumps(rng.bit_generator.state, sort_keys=True, default=lambda x: x.tolist())
        null_rng = np.random.default_rng(int.from_bytes(hashlib.sha256(state.encode()).digest()[:16], "big"))
    grouped = defaultdict(list)
    for value, block in zip(values, blocks, strict=True):
        grouped[block].append(value)
    result = {"n_loci": len(values), "n_blocks": len(grouped),
              "effect": math.fsum(values) / len(values) if len(values) else None,
              "ci_low": None, "ci_high": None, "p_value": None,
              "status": "insufficient_blocks" if len(values) else "not_estimable_no_loci",
              "test_method": "not_tested", "null_draws": 0, "bootstrap_replicates": 0,
              "mc_p_ci_low": None, "mc_p_ci_high": None,
              "n_nonzero_blocks": sum(math.fsum(v) != 0 for v in grouped.values())}
    if len(grouped) < 3:
        return result
    keys = sorted(grouped, key=lambda k: json.dumps(k, sort_keys=True))
    sums = np.array([math.fsum(grouped[k]) for k in keys])
    counts = np.array([len(grouped[k]) for k in keys])
    draws = rng.integers(0, len(sums), (replicates, len(sums)))
    estimates = sums[draws].sum(axis=1) / counts[draws].sum(axis=1)
    result["ci_low"], result["ci_high"] = map(float, np.quantile(estimates, [0.025, 0.975]))
    result.update(sign_flip(sums, null_rng, permutation_replicates, exact_max_states),
                  status="estimated", bootstrap_replicates=replicates)
    return result


def retention_rows(path, groups):
    loci = defaultdict(dict)
    for row in table(path, ("group_id", "block_id", "locus_id", "subgenome", "callable", "retained")):
        group = row["group_id"]
        if group not in groups or row["subgenome"] not in groups[group]:
            raise ValueError("Retention refers to an unmapped group/subgenome")
        if row["callable"] not in {"0", "1"} or row["retained"] not in {"0", "1"}:
            raise ValueError("callable and retained must be 0 or 1")
        key = (group, row["locus_id"])
        if row["subgenome"] in loci[key]:
            raise ValueError(f"Duplicate retention opportunity: {key}")
        if loci[key] and row["block_id"] != next(iter(loci[key].values()))["block_id"]:
            raise ValueError("Each locus must belong to exactly one resampling block")
        loci[key][row["subgenome"]] = row
    for (group, _), rows in loci.items():
        if set(rows) != groups[group]:
            raise ValueError("Each retention locus needs an explicit row for every expected subgenome")
    return loci


def expression_rows(analysis, genes):
    samples = table(analysis["samples_file"], ("column", "tissue", "biological_id"))
    if not samples or len({s["column"] for s in samples}) != len(samples):
        raise ValueError("Sample columns must be unique and nonempty")
    if len({s["tissue"] for s in samples}) > 1:
        for biological_id in {s["biological_id"] for s in samples}:
            if len({s["tissue"] for s in samples if s["biological_id"] == biological_id}) > 1:
                raise ValueError("Use tissue-specific biological IDs; cross-tissue pairing is not inferred")
    data = table(analysis["expression_file"], ("gene_id", *(s["column"] for s in samples)))
    expression = {}
    for row in data:
        if row["gene_id"] in expression:
            raise ValueError("Duplicate expression gene ID")
        values = {s["column"]: float(row[s["column"]]) for s in samples}
        if any(not math.isfinite(v) or v < 0 for v in values.values()):
            raise ValueError("Expression must contain finite, nonnegative, untransformed TPM/FPKM values")
        expression[row["gene_id"]] = values
    pairs, used = defaultdict(dict), set()
    for row in table(analysis["homoeolog_file"], ("pair_id", "block_id", "gene_id")):
        gene = row["gene_id"]
        if gene not in genes or gene in used:
            raise ValueError("Homoeolog genes must be mapped and occur in exactly one pair")
        used.add(gene)
        mapping = genes[gene]
        key = (mapping["group_id"], row["pair_id"])
        sg = mapping["subgenome"]
        if sg in pairs[key] or (pairs[key] and row["block_id"] != next(iter(pairs[key].values()))[1]):
            raise ValueError("Homoeolog sets require one gene per subgenome and a single block")
        pairs[key][sg] = (gene, row["block_id"])
    if any(len(p) < 2 for p in pairs.values()):
        raise ValueError("Homoeolog sets need at least two subgenomes")
    if pairs and not any(sum(gene in expression for gene, _ in pair.values()) >= 2 for pair in pairs.values()):
        raise ValueError("No homoeolog gene pairs match the expression matrix identifiers")
    by_tissue = defaultdict(lambda: defaultdict(list))
    for sample in samples:
        by_tissue[sample["tissue"]][sample["biological_id"]].append(sample["column"])
    rows, coverage = [], []
    for (group, pair_id), pair in sorted(pairs.items()):
        for a, b in itertools.combinations(sorted(pair), 2):
            ga, block = pair[a]
            gb = pair[b][0]
            for tissue, biological in by_tissue.items():
                ratios, positive_a, positive_b = [], 0, 0
                for columns in biological.values():
                    va = log_mean_abundance([expression[ga][c] for c in columns]) if ga in expression else None
                    vb = log_mean_abundance([expression[gb][c] for c in columns]) if gb in expression else None
                    positive_a += int(va is not None)
                    positive_b += int(vb is not None)
                    if ga not in expression or gb not in expression:
                        continue
                    if va is not None and vb is not None:
                        ratios.append(va - vb)
                coverage.append({"group_id": group, "pair_id": pair_id, "subgenome_a": a, "subgenome_b": b,
                                 "tissue": tissue, "biological_samples": len(biological), "positive_samples": len(ratios),
                                 "gene_a": ga, "gene_b": gb, "block_id": block,
                                 "positive_a_samples": positive_a, "positive_b_samples": positive_b,
                                 "gene_ids_present": int(ga in expression and gb in expression)})
                # Use the same pair set across biological samples. Zeros are audited,
                # never converted to arbitrary pseudocounts or treated as loss.
                if len(ratios) == len(biological):
                    rows.append({"group_id": group, "pair_id": pair_id, "block_id": block,
                                 "subgenome_a": a, "subgenome_b": b, "tissue": tissue,
                                 "log2_ratio": math.fsum(ratios) / len(ratios), "biological_samples": len(ratios)})
    return rows, coverage


def log_mean_abundance(values):
    """Log2 of an arithmetic mean, without overflowing or underflowing the mean."""
    scale = max(values)
    if scale == 0:
        return None
    return math.log2(scale) + math.log2(math.fsum(v / scale for v in values) / len(values))


def adjust_p(rows):
    for row in rows:
        p = row["p_value"]
        if p is not None and (isinstance(p, bool) or not math.isfinite(p) or not 0 <= p <= 1):
            raise ValueError("p_value must be finite and between zero and one")
        row["q_value"] = None
    selected = [(i, row["p_value"]) for i, row in enumerate(rows) if row["p_value"] is not None]
    running = 1.0
    for rank, (i, p) in reversed(list(enumerate(sorted(selected, key=lambda item: item[1]), 1))):
        running = min(running, p * len(selected) / rank)
        rows[i]["q_value"] = running
    for row in rows:
        row["multiple_testing_n"] = len(selected)


def analyse(analysis, output, replicates, seed, permutation_replicates=DEFAULT_PERMUTATIONS,
            exact_max_states=DEFAULT_EXACT_STATES, plot_config=None):
    validate_resampling(replicates, permutation_replicates, exact_max_states, seed)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    genes, groups, scope = load_mapping(analysis["mapping_file"])
    def infer(values, blocks, metric, group, a, b, tissue=""):
        seeds = {purpose: contrast_seed(seed, analysis["analysis_id"], metric, group, a, b, tissue, purpose)
                 for purpose in ("bootstrap", "permutation")}
        result = inference(values, blocks, replicates, np.random.default_rng(seeds["bootstrap"]),
                           null_rng=np.random.default_rng(seeds["permutation"]),
                           permutation_replicates=permutation_replicates, exact_max_states=exact_max_states)
        return {**result, "bootstrap_seed": str(seeds["bootstrap"]), "permutation_seed": str(seeds["permutation"]),
                "inference_version": INFERENCE_VERSION}
    results = []
    if analysis["retention_file"]:
        loci = retention_rows(analysis["retention_file"], groups)
        contrasts = list(sorted(groups.items()))
        if scope == "global":
            ancestral_groups = defaultdict(set)
            for group, locus in loci:
                ancestral_groups[locus].add(group)
            if any(len(g) > 1 for g in ancestral_groups.values()):
                raise ValueError("Global retention requires unique ancestral loci across groups")
            contrasts.append(("ALL_GROUPS", next(iter(groups.values()))))
        for group, sgs in contrasts:
            for a, b in itertools.combinations(sorted(sgs), 2):
                selected = [rows for (g, _), rows in loci.items() if g == group or group == "ALL_GROUPS"]
                called = [r for r in selected if r[a]["callable"] == r[b]["callable"] == "1"]
                va = [int(r[a]["retained"]) for r in called]
                vb = [int(r[b]["retained"]) for r in called]
                result = infer(np.subtract(va, vb), [(r[a]["group_id"], r[a]["block_id"]) for r in called],
                               "retention_difference", group, a, b)
                results.append({"metric": "retention_difference", "group_id": group, "subgenome_a": a,
                                "subgenome_b": b, "tissue": "", "n_opportunities": len(selected),
                                "retained_a": sum(va), "retained_b": sum(vb), **result})
    expression, coverage = [], []
    if analysis["expression_file"]:
        expression, coverage = expression_rows(analysis, genes)
        buckets = defaultdict(list)
        for row in coverage:
            buckets[(row["group_id"], row["subgenome_a"], row["subgenome_b"], row["tissue"])]
            if scope == "global":
                buckets[("ALL_GROUPS", row["subgenome_a"], row["subgenome_b"], row["tissue"])]
        for row in expression:
            buckets[(row["group_id"], row["subgenome_a"], row["subgenome_b"], row["tissue"])].append(row)
            if scope == "global":
                buckets[("ALL_GROUPS", row["subgenome_a"], row["subgenome_b"], row["tissue"])].append(row)
        for (group, a, b, tissue), rows in sorted(buckets.items()):
            results.append({"metric": "expression_log2_ratio", "group_id": group, "subgenome_a": a,
                            "subgenome_b": b, "tissue": tissue, "n_opportunities": len([r for r in coverage if
                            (r["group_id"] == group or group == "ALL_GROUPS") and
                            (r["subgenome_a"], r["subgenome_b"], r["tissue"]) == (a, b, tissue)]),
                            "retained_a": None, "retained_b": None,
                            **infer([r["log2_ratio"] for r in rows], [(r["group_id"], r["block_id"]) for r in rows],
                                    "expression_log2_ratio", group, a, b, tissue)})
        detection = defaultdict(list)
        for row in coverage:
            if not row["gene_ids_present"]:
                continue
            detection[row["group_id"], row["subgenome_a"], row["subgenome_b"], row["tissue"]].append(row)
            if scope == "global":
                detection["ALL_GROUPS", row["subgenome_a"], row["subgenome_b"], row["tissue"]].append(row)
        for (group, a, b, tissue), rows in sorted(detection.items()):
            values = [(r["positive_a_samples"] - r["positive_b_samples"]) / r["biological_samples"] for r in rows]
            results.append({"metric": "expression_detection_difference", "group_id": group, "subgenome_a": a,
                            "subgenome_b": b, "tissue": tissue, "n_opportunities": len(rows),
                            "retained_a": None, "retained_b": None,
                            **infer(values, [(r["group_id"], r["block_id"]) for r in rows],
                                    "expression_detection_difference", group, a, b, tissue)})
    for metric in {r["metric"] for r in results}:
        adjust_p([r for r in results if r["metric"] == metric])
    fields = ("metric", "group_id", "subgenome_a", "subgenome_b", "tissue", "n_opportunities",
              "retained_a", "retained_b", "n_loci", "n_blocks", "effect", "ci_low", "ci_high", "p_value", "q_value", "status",
              "test_method", "null_draws", "bootstrap_replicates", "n_nonzero_blocks", "mc_p_ci_low", "mc_p_ci_high",
              "bootstrap_seed", "permutation_seed", "inference_version", "multiple_testing_n")
    write_table(output / "statistics.tsv", results, fields)
    write_table(output / "expression_pairs.tsv", expression, ("group_id", "pair_id", "block_id", "subgenome_a", "subgenome_b", "tissue", "log2_ratio", "biological_samples"))
    write_table(output / "expression_coverage.tsv", coverage, ("group_id", "pair_id", "block_id", "gene_a", "gene_b", "subgenome_a", "subgenome_b", "tissue", "biological_samples", "positive_samples", "positive_a_samples", "positive_b_samples", "gene_ids_present"))
    summary = {"schema_version": 2, "analysis_id": analysis["analysis_id"], "species": analysis["species"],
               "reference": analysis.get("reference", ""), "pair_set": analysis.get("pair_set", ""),
               "assignment_scope": scope, "mapped_genes": len(genes), "groups": len(groups),
               "genomewide_identity": "declared_independent" if scope == "global" else "unresolved",
               "retention_status": ("validated_callable_loci" if analysis.get("retention_validated") == "1" else "exploratory_syntelog_detection") if analysis["retention_file"] else "not_estimable_missing_callable_outgroup_loci",
               "expression_status": ("mapping_validated" if analysis.get("mapping_validated") == "1" else "exploratory_mapping_not_validated") if analysis["expression_file"] else "not_estimable_missing_expression",
               "expression_unit": analysis.get("expression_unit", "TPM"),
               "expression_ratio_status": "estimated" if expression else "not_estimable_no_positive_complete_pairs",
               "genomewide_dominance": "contrasts_available" if scope == "global" and results else "not_tested", "statistics": results,
               "inference": {"resampling": "nonoverlapping_block_bootstrap", "null_test": "block_sign_flip",
                             "version": INFERENCE_VERSION, "replicates": replicates, "seed": seed,
                             "floating_null_comparison": "scaled block sums; 100 machine eps relative tolerance",
                             "permutation_replicates": permutation_replicates, "exact_max_states": exact_max_states,
                             "rng_scheme": "sha256-contrast-v1; separate bootstrap/permutation streams",
                             "point_unit": "group-specific subgenome contrast", "resampling_unit": "nonoverlapping block",
                             "effect": "A minus B retention/detection fraction or mean log2(normalised_abundance_A/normalised_abundance_B)",
                             "ci_scope": "loci conditional on sampled tissues, mapping and annotation",
                             "multiple_testing": "Benjamini–Hochberg within each metric and analysis across groups, contrasts and tissues"}}
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    if any(r["effect"] is not None and r["metric"] in (plot_config or {}).get("metrics", plotting.METRICS) for r in results):
        plot(results, output, analysis["species"], analysis.get("expression_unit", "TPM"), config=plot_config, metadata=summary)
        plot(results, output, analysis["species"], analysis.get("expression_unit", "TPM"), absolute=True,
             config=plot_config, metadata=summary)
    return summary


def plot(rows, output, species, unit="TPM", absolute=False, *, config=None, metadata=None):
    metadata = metadata or {"analysis_id": "contrast", "species": species, "expression_unit": unit}
    context = {key: metadata[key] for key in ("analysis_id", "species", "expression_unit", "assignment_scope",
                                             "reference", "pair_set", "retention_status") if key in metadata}
    individual_config = dict(config or {})
    individual_config.update(width_pt=individual_config.get("individual_width_pt"),
                             height_pt=individual_config.get("individual_height_pt"))
    return plotting.comparison_plot([{**context, "statistics": rows}], output, config=individual_config, absolute=absolute,
                                    stem="contrasts_absolute" if absolute else "contrasts", apply_filters=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("plan")
    p.add_argument("--workspace", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--replicates", type=int, default=2000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--permutation-replicates", type=int, default=DEFAULT_PERMUTATIONS)
    p.add_argument("--exact-max-states", type=int, default=DEFAULT_EXACT_STATES)
    p.add_argument("--plot-config")
    p.add_argument("--outfile", required=True)
    p = sub.add_parser("report", help="Render existing statistics without repeating inference")
    p.add_argument("--manifest", required=True)
    p.add_argument("--plot-config")
    p.add_argument("--output", required=True)
    for command in ("contract", "verify", "run"):
        p = sub.add_parser(command)
        p.add_argument("--plan", required=True)
        if command == "run":
            p.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.command == "report":
        plotting.report(args.manifest, args.output, args.plot_config)
        return
    if args.command == "plan":
        plan = make_plan(args.workspace, args.manifest, args.replicates, args.seed,
                         args.permutation_replicates, args.exact_max_states, args.plot_config)
        Path(args.outfile).write_text(json.dumps(plan, indent=2) + "\n")
        return
    plan = json.loads(Path(args.plan).read_text())
    verify_inputs(plan)
    if args.command == "contract":
        print("\0".join(contract_args(plan)) + "\0", end="")
    elif args.command == "run":
        output = Path(args.output)
        output.mkdir(parents=True, exist_ok=False)
        summaries = [analyse(row, output / row["analysis_id"], plan["replicates"], plan["seed"],
                             plan.get("permutation_replicates", DEFAULT_PERMUTATIONS),
                             plan.get("exact_max_states", DEFAULT_EXACT_STATES), plan.get("plot_config"))
                     for row in plan["analyses"]]
        available = [{**{k: s[k] for k in ("analysis_id", "species", "expression_unit", "assignment_scope",
                                            "reference", "pair_set", "retention_status")},
                      "statistics": s["statistics"]} for s in summaries]
        if available:
            for absolute in (False, True):
                plotting.comparison_plot(available, output, config=plan.get("plot_config"), absolute=absolute,
                                         input_hashes=plan["input_hashes"])
        verify_inputs(plan)
        (output / "run.json").write_text(json.dumps({"plan": plan, "summaries": summaries}, indent=2) + "\n")


if __name__ == "__main__":
    main()

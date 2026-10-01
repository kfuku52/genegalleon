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
    if any(None in row or any(row[k] is None or not row[k].strip() for k in required) for row in rows):
        raise ValueError(f"Incomplete required cells or malformed rows: {path}")
    return rows


def write_table(path, rows, fields):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def make_plan(workspace, manifest, replicates=2000, seed=1):
    workspace = Path(workspace).resolve()
    manifest = Path(manifest).resolve()
    if replicates < 100 or seed < 0:
        raise ValueError("Require at least 100 bootstrap replicates and a nonnegative seed")
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
    return {"schema_version": 1, "workspace": str(workspace), "analyses": analyses,
            "replicates": replicates, "seed": seed, "input_hashes": inputs,
            "numpy_version": np.__version__}


def verify_inputs(plan):
    for path, expected in plan["input_hashes"].items():
        if digest(path) != expected:
            raise ValueError(f"Input changed during subgenome analysis: {path}")


def contract_args(plan):
    workspace = Path(plan["workspace"])
    args = ["--manifest", str(workspace / "output/artifact_provenance/genome_evolution/subgenome_dominance.json"),
            "--step", "subgenome_dominance", "--family-id", "all_analyses",
            "--logical-root", str(workspace / "output/.gg_global_artifacts"), "--workspace-root", str(workspace),
            "--output", f"results={workspace / 'output/genome_evolution/subgenome_dominance'}",
            "--input", f"implementation={Path(__file__).resolve()}"]
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


def inference(values, blocks, replicates, rng):
    """Locus-weighted mean, cluster bootstrap CI and block sign-flip null test.

    Input blocks must partition loci, and must be defined independently of bias.
    The CI describes loci conditional on the sampled tissues and annotation.
    """
    values = np.asarray(values, dtype=float)
    grouped = defaultdict(list)
    for value, block in zip(values, blocks, strict=True):
        grouped[block].append(value)
    result = {"n_loci": len(values), "n_blocks": len(grouped), "effect": float(values.mean()) if len(values) else None,
              "ci_low": None, "ci_high": None, "p_value": None, "status": "insufficient_blocks"}
    if len(grouped) < 3:
        return result
    sums = np.array([sum(v) for v in grouped.values()])
    counts = np.array([len(v) for v in grouped.values()])
    draws = rng.integers(0, len(sums), (replicates, len(sums)))
    estimates = sums[draws].sum(axis=1) / counts[draws].sum(axis=1)
    result["ci_low"], result["ci_high"] = map(float, np.quantile(estimates, [0.025, 0.975]))
    observed = abs(sums.sum())
    if len(sums) <= 16:
        null = np.array([abs(np.dot(signs, sums)) for signs in itertools.product((-1, 1), repeat=len(sums))])
        p = float(np.mean(null >= observed - 1e-12))
    else:
        null = abs((rng.choice([-1, 1], (replicates, len(sums))) * sums).sum(axis=1))
        p = float((1 + np.sum(null >= observed - 1e-12)) / (replicates + 1))
    result.update(p_value=p, status="estimated")
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
                    va = np.mean([expression[ga][c] for c in columns]) if ga in expression else None
                    vb = np.mean([expression[gb][c] for c in columns]) if gb in expression else None
                    positive_a += int(va is not None and va > 0)
                    positive_b += int(vb is not None and vb > 0)
                    if ga not in expression or gb not in expression:
                        continue
                    if va > 0 and vb > 0:
                        ratios.append(float(np.log2(va / vb)))
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
                                 "log2_ratio": float(np.mean(ratios)), "biological_samples": len(ratios)})
    return rows, coverage


def adjust_p(rows):
    selected = [(i, row["p_value"]) for i, row in enumerate(rows) if row["p_value"] is not None]
    running = 1.0
    for rank, (i, p) in reversed(list(enumerate(sorted(selected, key=lambda item: item[1]), 1))):
        running = min(running, p * len(selected) / rank)
        rows[i]["q_value"] = running
    for row in rows:
        row.setdefault("q_value", None)


def analyse(analysis, output, replicates, seed):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    genes, groups, scope = load_mapping(analysis["mapping_file"])
    rng = np.random.default_rng(seed)
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
                result = inference(np.subtract(va, vb), [(r[a]["group_id"], r[a]["block_id"]) for r in called], replicates, rng)
                results.append({"metric": "retention_difference", "group_id": group, "subgenome_a": a,
                                "subgenome_b": b, "tissue": "", "n_opportunities": len(selected),
                                "retained_a": sum(va), "retained_b": sum(vb), **result})
    expression, coverage = [], []
    if analysis["expression_file"]:
        expression, coverage = expression_rows(analysis, genes)
        buckets = defaultdict(list)
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
                            **inference([r["log2_ratio"] for r in rows], [(r["group_id"], r["block_id"]) for r in rows], replicates, rng)})
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
                            **inference(values, [(r["group_id"], r["block_id"]) for r in rows], replicates, rng)})
    for metric in {r["metric"] for r in results}:
        adjust_p([r for r in results if r["metric"] == metric])
    fields = ("metric", "group_id", "subgenome_a", "subgenome_b", "tissue", "n_opportunities",
              "retained_a", "retained_b", "n_loci", "n_blocks", "effect", "ci_low", "ci_high", "p_value", "q_value", "status")
    write_table(output / "statistics.tsv", results, fields)
    write_table(output / "expression_pairs.tsv", expression, ("group_id", "pair_id", "block_id", "subgenome_a", "subgenome_b", "tissue", "log2_ratio", "biological_samples"))
    write_table(output / "expression_coverage.tsv", coverage, ("group_id", "pair_id", "block_id", "gene_a", "gene_b", "subgenome_a", "subgenome_b", "tissue", "biological_samples", "positive_samples", "positive_a_samples", "positive_b_samples", "gene_ids_present"))
    summary = {"schema_version": 1, "analysis_id": analysis["analysis_id"], "species": analysis["species"],
               "assignment_scope": scope, "mapped_genes": len(genes), "groups": len(groups),
               "genomewide_identity": "declared_independent" if scope == "global" else "unresolved",
               "retention_status": ("validated_callable_loci" if analysis.get("retention_validated") == "1" else "exploratory_syntelog_detection") if analysis["retention_file"] else "not_estimable_missing_callable_outgroup_loci",
               "expression_status": ("mapping_validated" if analysis.get("mapping_validated") == "1" else "exploratory_mapping_not_validated") if analysis["expression_file"] else "not_estimable_missing_expression",
               "expression_unit": analysis.get("expression_unit", "TPM"),
               "expression_ratio_status": "estimated" if expression else "not_estimable_no_positive_complete_pairs",
               "genomewide_dominance": "contrasts_available" if scope == "global" and results else "not_tested", "statistics": results,
               "inference": {"resampling": "nonoverlapping_block_bootstrap", "null_test": "block_sign_flip",
                             "replicates": replicates, "seed": seed, "effect": "A minus B retention/detection fraction or mean log2(normalised_abundance_A/normalised_abundance_B)",
                             "ci_scope": "loci conditional on sampled tissues, mapping and annotation",
                             "multiple_testing": "BH within each metric and analysis across groups, contrasts and tissues"}}
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    if results:
        plot(results, output, analysis["species"], analysis.get("expression_unit", "TPM"))
        plot(results, output, analysis["species"], analysis.get("expression_unit", "TPM"), absolute=True)
    return summary


def plot(rows, output, species, unit="TPM", absolute=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    metrics = sorted({r["metric"] for r in rows})
    fig, axes = plt.subplots(len(metrics), 1, figsize=(7.2, max(3, len(rows) * 0.22 + 1.4)), squeeze=False)
    for axis, metric in zip(axes[:, 0], metrics, strict=True):
        selected = [r for r in rows if r["metric"] == metric]
        for y, row in enumerate(selected):
            if row["effect"] is None:
                continue
            effect = abs(row["effect"]) if absolute else row["effect"]
            axis.plot(effect, y, "o", color="#365d8d", markersize=4)
            if row["ci_low"] is not None:
                low, high = row["ci_low"], row["ci_high"]
                if absolute:
                    low, high = (0 if low <= 0 <= high else min(abs(low), abs(high))), max(abs(low), abs(high))
                axis.plot([low, high], [y, y], color="#365d8d")
        axis.set_yticks(range(len(selected)), [f"{r['group_id']} {r['subgenome_a']}/{r['subgenome_b']} {r['tissue']}".strip() for r in selected], fontsize=8)
        axis.axvline(0, color="0.6", lw=0.8)
        label = {"retention_difference": "Retention fraction difference",
                         "expression_detection_difference": f"Expression detection fraction difference ({unit} > 0)",
                         "expression_log2_ratio": f"Mean log2 {unit} ratio (A/B)"}[metric]
        axis.set_xlabel("Absolute value: " + label if absolute else label)
        axis.invert_yaxis()
        axis.spines[["top", "right"]].set_visible(False)
    subtitle = "Bias magnitudes; signed 95% intervals mapped under abs" if absolute else "Contrasts; 95% block bootstrap intervals"
    fig.suptitle(species.replace("_", " ") + "\n" + subtitle, fontsize=11)
    fig.tight_layout()
    for extension in ("png", "svg"):
        stem = "contrasts_absolute" if absolute else "contrasts"
        fig.savefig(output / f"{stem}.{extension}", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("plan")
    p.add_argument("--workspace", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--replicates", type=int, default=2000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--outfile", required=True)
    for command in ("contract", "verify", "run"):
        p = sub.add_parser(command)
        p.add_argument("--plan", required=True)
        if command == "run":
            p.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.command == "plan":
        plan = make_plan(args.workspace, args.manifest, args.replicates, args.seed)
        Path(args.outfile).write_text(json.dumps(plan, indent=2) + "\n")
        return
    plan = json.loads(Path(args.plan).read_text())
    verify_inputs(plan)
    if args.command == "contract":
        print("\0".join(contract_args(plan)) + "\0", end="")
    elif args.command == "run":
        output = Path(args.output)
        output.mkdir(parents=True, exist_ok=False)
        summaries = [analyse(row, output / row["analysis_id"], plan["replicates"], plan["seed"]) for row in plan["analyses"]]
        verify_inputs(plan)
        (output / "run.json").write_text(json.dumps({"plan": plan, "summaries": summaries}, indent=2) + "\n")


if __name__ == "__main__":
    main()

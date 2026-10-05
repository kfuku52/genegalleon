#!/usr/bin/env python3
"""Evaluate paired representative DNA CDS sets with one frozen BUSCO contract."""

from __future__ import annotations

import argparse
import csv
import fcntl
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from busco_quality_metadata import parse_short_summary
from gene_model_refinement import verify_inputs
from input_generation_array_state import FreshDigestBatch, atomic_json, digest
from species_labeling import extract_species_label


def input_pairs(root, cds_dir=None):
    """Native sources are frozen; extra CDS-only species are explicit passthroughs."""
    plan = json.loads((root / "plan.json").read_text())
    effective = {r["species"]: r for r in verify_inputs(root / "effective" / "inputs.tsv")}
    if set(effective) != set(plan["request"]["sources"]):
        raise ValueError("Refinement species membership differs from the frozen plan")
    pairs = {}
    for name, source in plan["request"]["sources"].items():
        before = Path(source["fasta"])
        if digest(before) != plan["request"]["files"][str(before)]:
            raise ValueError("Refinement source CDS changed: " + str(before))
        pairs[name] = {"species": name, "before": str(before), "after": effective[name]["cds"],
                       "refinement_status": "analysed", "reason": ""}
    seen = set()
    if cds_dir:
        for path in sorted(cds_dir.iterdir()):
            if path.name.startswith(".") or not path.is_file() or not re.search(r"\.(fa|fasta|fna)(\.gz)?$", path.name):
                continue
            name = extract_species_label(path.name, strip_extension=True)
            if not name:
                raise ValueError("Cannot identify CDS species: " + path.name)
            if name in seen:
                raise ValueError("Ambiguous CDS files for species: " + name)
            seen.add(name)
            if name not in pairs:
                pairs[name] = {"species": name, "before": str(path.resolve()), "after": str(path.resolve()),
                               "refinement_status": "not_analysed",
                               "reason": "No matching genome/GFF in the refinement plan; CDS retained unchanged"}
    return [pairs[n] for n in sorted(pairs)]


def read_result(path):
    text = path.read_text()
    row = parse_short_summary(text, path)
    fields = {
        "complete": r"Complete BUSCOs \(C\)", "single": r"Complete and single-copy BUSCOs \(S\)",
        "duplicated": r"Complete and duplicated BUSCOs \(D\)", "fragmented": r"Fragmented BUSCOs \(F\)",
        "missing": r"Missing BUSCOs \(M\)", "total": r"Total BUSCO groups searched",
    }
    for name, label in fields.items():
        hits = re.findall(r"^\s*(\d+)\s+" + label, text, re.M)
        if len(hits) != 1:
            raise ValueError("Missing/ambiguous BUSCO count: " + name)
        row[name] = int(hits[0])
    if (row["total"] <= 0 or row["complete"] != row["single"] + row["duplicated"]
            or row["complete"] + row["fragmented"] + row["missing"] != row["total"]):
        raise ValueError("Inconsistent BUSCO counts")
    match = re.search(r"Creation date:\s*([^,)]+)", text)
    row["lineage_creation_date"] = match.group(1).strip() if match else ""
    row["dependencies"] = dict(re.findall(r"^\s*(hmmsearch|metaeuk):\s*(\S+)", text, re.M))
    for key in ("lineage", "mode", "busco_version", "lineage_creation_date"):
        if not row[key]:
            raise ValueError("Missing BUSCO comparability metadata: " + key)
    return row


def paired_result(pair, before, after):
    keys = ("lineage", "lineage_creation_date", "total", "mode", "busco_version", "dependencies")
    if any(before[k] != after[k] for k in keys):
        raise ValueError("Noncomparable BUSCO pair: " + pair["species"])
    return dict(pair, before_result=before, after_result=after,
                delta_complete=after["complete"] - before["complete"],
                delta_complete_pp=100 * (after["complete"] - before["complete"]) / before["total"],
                delta_duplicated=after["duplicated"] - before["duplicated"])


def run_one(pair, phase, report, contract, cpus):
    source = Path(pair[phase])
    source_hash = digest(source)
    key = {"contract": contract, "source_sha256": source_hash}
    directory = report / "runs" / pair["species"] / phase
    directory.mkdir(parents=True, exist_ok=True)
    summary, receipt = directory / "summary.txt", directory / "receipt.json"
    with (directory / ".lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if receipt.exists():
            saved = json.loads(receipt.read_text())
            if saved["key"] != key or digest(summary) != saved["summary_sha256"]:
                raise ValueError("BUSCO cache changed; use a new report directory")
            return read_result(summary)
        # Identical bytes and an identical complete contract permit reuse, including excluded species.
        previous = report / "runs" / pair["species"] / "before"
        if phase == "after" and (previous / "receipt.json").exists():
            saved = json.loads((previous / "receipt.json").read_text())
            if saved["key"] == key and digest(previous / "summary.txt") == saved["summary_sha256"]:
                shutil.copyfile(previous / "summary.txt", summary)
                atomic_json(receipt, {"key": key, "summary_sha256": digest(summary), "reused_identical_before": True})
                return read_result(summary)
        # Scratch contains large temporary predictor files; only the bound summary and log persist.
        with tempfile.TemporaryDirectory(prefix="gg-refinement-busco-") as scratch:
            scratch = Path(scratch)
            opener = gzip.open if source.suffix == ".gz" else open
            with opener(source, "rb") as src, (scratch / "input.fasta").open("wb") as dst:
                shutil.copyfileobj(src, dst)
            support = Path(__file__).resolve().parent
            command = ["bash", "-c", 'source "$1"; gg_run_busco_with_metaeuk_modified_fas_compat "${@:2}"',
                       "gg-busco", str(support / "gg_busco.sh"), "--in", str(scratch / "input.fasta"),
                       "--mode", "transcriptome", "--out", "busco", "--out_path", str(scratch),
                       "--cpu", str(cpus), "--evalue", "1e-03", "--limit", "20", "--offline",
                       "--lineage_dataset", contract["lineage_path"], "--download_path", contract["download_path"]]
            with (directory / "busco.log").open("w") as log:
                subprocess.run(command, cwd=scratch, stdout=log, stderr=subprocess.STDOUT, check=True)
            summaries = sorted((scratch / "busco").glob("short_summary.specific.*.txt"))
            if len(summaries) != 1:
                raise ValueError("Expected one BUSCO specific summary")
            shutil.copyfile(summaries[0], summary)
        if digest(source) != source_hash:
            raise OSError("BUSCO input mutated during evaluation")
        result = read_result(summary)
        atomic_json(receipt, {"key": key, "summary_sha256": digest(summary), "reused_identical_before": False})
        return result


def plot_comparison(rows, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 3, figsize=(19, max(6, .43 * len(rows) + 2)),
                             gridspec_kw={"width_ratios": [1, 1, .65], "wspace": .12}, sharey=True)
    colors = ["#2a8d75", "#8bc9ab", "#e6b757", "#d6dbe2"]
    for i, row in enumerate(rows):
        if row["refinement_status"] == "not_analysed":
            for ax in axes:
                ax.axhspan(i - .48, i + .48, color="#eef0f3", zorder=0)
        for j, phase in enumerate(("before", "after")):
            result, left = row[phase + "_result"], 0
            for category, color in zip(("single", "duplicated", "fragmented", "missing"), colors, strict=True):
                value = 100 * result[category] / result["total"]
                axes[j].barh(i, value, left=left, color=color, height=.7)
                left += value
            axes[j].text(101, i, f'{100 * result["complete"] / result["total"]:.2f}', va="center", fontsize=9)
        delta = row["delta_complete_pp"]
        axes[2].barh(i, delta, color="#2a8d75" if delta >= 0 else "#b85250", height=.65)
        axes[2].text(.98, i, f'{row["delta_complete"]:+d} ({delta:+.2f} pp)',
                     transform=axes[2].get_yaxis_transform(), ha="right", va="center", fontsize=9)
    labels = [r["species"].replace("_", " ") + ("  [not analysed]" if r["refinement_status"] == "not_analysed" else "") for r in rows]
    axes[0].set_yticks(range(len(rows)), labels)
    axes[0].invert_yaxis()
    for ax, title in zip(axes, ("Before refinement", "After refinement", "Change in complete BUSCOs"), strict=True):
        ax.set_title(title, fontweight="bold", pad=16)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.grid(axis="x", alpha=.15)
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=0)
    for ax in axes[:2]:
        ax.set_xlim(0, 116)
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.set_xlabel("BUSCO groups (%)     C (%) at right")
    limit = max(1, max(abs(r["delta_complete_pp"]) for r in rows) * 1.2)
    axes[2].set_xlim(-limit, limit * 2.5)
    axes[2].axvline(0, color="#647383", linewidth=.7)
    axes[2].set_xlabel("Percentage points (pp)")
    fig.subplots_adjust(left=.24, right=.98, top=.90, bottom=.14)
    fig.suptitle("Representative CDS completeness before and after refinement", x=.24, ha="left", y=.98, fontsize=17, fontweight="bold")
    identity = rows[0]["before_result"]
    fig.text(.24, .94, f'BUSCO {identity["busco_version"]}; {identity["lineage"]} ({identity["lineage_creation_date"]}); '
             f'n = {identity["total"]}; transcriptome mode; one representative per locus', fontsize=10)
    fig.legend([Patch(facecolor=c) for c in colors], ["Single-copy", "Duplicated", "Fragmented", "Missing"],
               loc="lower left", bbox_to_anchor=(.24, .065), ncol=4, frameon=False)
    fig.text(.24, .025, "Grey rows: excluded from structural refinement; unchanged CDS are still evaluated by BUSCO.\n"
             "Before = refinement source CDS (including earlier rescued genes); after = selected DNA CDS, not all isoforms.", fontsize=10)
    for suffix in ("png", "svg"):
        fig.savefig(output / ("busco_comparison." + suffix), dpi=180, facecolor="white")
    plt.close(fig)


def evaluate(pairs, report, lineage, download_path, cpus=4, jobs=1):
    if not pairs:
        raise ValueError("Empty BUSCO comparison")
    report.mkdir(parents=True, exist_ok=True)
    files = sorted(p for p in lineage.rglob("*") if p.is_file() and not p.name.startswith("."))
    if not (lineage / "dataset.cfg").is_file() or not files:
        raise ValueError("A complete local BUSCO lineage is required")
    batch = FreshDigestBatch()
    hashes = batch.read(files)
    batch.read([pair[phase] for pair in pairs for phase in ("before", "after")])
    lineage_hash = hashlib.sha256(json.dumps({str(p.relative_to(lineage)): hashes[str(p)] for p in files}, sort_keys=True).encode()).hexdigest()
    version = subprocess.check_output(["busco", "--version"], text=True).strip()
    support = Path(__file__).resolve().parent
    contract = {"schema": 1, "busco_version": version, "lineage_sha256": lineage_hash,
                "lineage_path": str(lineage), "download_path": str(download_path),
                "mode": "transcriptome", "evalue": "1e-03", "limit": 20,
                "cpus_per_job": cpus,
                "tool_sha256": {tool: digest(shutil.which(tool)) for tool in ("busco", "metaeuk", "hmmsearch")},
                "implementation": digest(Path(__file__)), "wrapper_sha256": digest(support / "gg_busco.sh"),
                "hmmsearch_wrapper_sha256": digest(support / "gg_wrapper_bin/hmmsearch")}
    atomic_json(report / "contract.json", {"contract": contract, "pairs": pairs}, immutable=True)
    def one(pair):
        before = run_one(pair, "before", report, contract, cpus)
        after = run_one(pair, "after", report, contract, cpus)
        row = paired_result(pair, before, after)
        print(json.dumps({"species": pair["species"], "delta_complete": row["delta_complete"]}), flush=True)
        return row
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        rows = list(pool.map(one, pairs))
    batch.check()
    atomic_json(report / "busco_comparison.json", {"contract": contract, "species": rows})
    with (report / "busco_comparison.tsv").open("w", newline="") as handle:
        columns = ("species", "refinement_status", "reason", "before_complete", "after_complete", "delta_complete", "delta_complete_pp", "delta_duplicated")
        writer = csv.DictWriter(handle, columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in columns if k not in ("before_complete", "after_complete")}
                            | {"before_complete": row["before_result"]["complete"], "after_complete": row["after_result"]["complete"]})
    plot_comparison(rows, report)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--cds-dir", type=Path, help="All dataset species, including unrefined CDS-only species")
    parser.add_argument("--lineage", required=True, type=Path, help="Frozen local lineage directory")
    parser.add_argument("--download-path", required=True, type=Path)
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--jobs", type=int, default=1, help="Total CPU budget = jobs times cpus")
    args = parser.parse_args()
    root, report = args.output.resolve(), args.report.resolve()
    if args.cpus < 1 or args.jobs < 1:
        parser.error("CPU and job counts must be positive")
    if root == report or root in report.parents or report in root.parents:
        parser.error("Report must be separate from the immutable refinement tree")
    evaluate(input_pairs(root, args.cds_dir), report, args.lineage.resolve(), args.download_path.resolve(), args.cpus, args.jobs)


if __name__ == "__main__":
    main()

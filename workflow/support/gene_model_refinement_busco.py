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
from dated_tree_presentation import STATUS, STATUS_COLOURS, STATUS_LABELS
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
    full_table = directory / "full_table.tsv"
    with (directory / ".lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if receipt.exists():
            saved = json.loads(receipt.read_text())
            if (saved["key"] != key or digest(summary) != saved["summary_sha256"]
                    or not full_table.is_file() or digest(full_table) != saved.get("full_table_sha256")):
                raise ValueError("BUSCO cache changed; use a new report directory")
            return read_result(summary)
        # Identical bytes and an identical complete contract permit reuse, including excluded species.
        previous = report / "runs" / pair["species"] / "before"
        if phase == "after" and (previous / "receipt.json").exists():
            saved = json.loads((previous / "receipt.json").read_text())
            if (saved["key"] == key and digest(previous / "summary.txt") == saved["summary_sha256"]
                    and (previous / "full_table.tsv").is_file()
                    and digest(previous / "full_table.tsv") == saved.get("full_table_sha256")):
                shutil.copyfile(previous / "summary.txt", summary)
                shutil.copyfile(previous / "full_table.tsv", full_table)
                atomic_json(receipt, {"key": key, "summary_sha256": digest(summary),
                                      "full_table_sha256": digest(full_table), "reused_identical_before": True})
                return read_result(summary)
        # Preserve small QC tables and logs; discard large temporary predictor files.
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
            tables = list((scratch / "busco").glob("run_*/full_table.tsv")) + list((scratch / "busco").glob("run_*/full_table.tsv.gz"))
            if len(tables) != 1:
                raise ValueError("Expected one BUSCO full table")
            raw_table = tables[0]
            if raw_table.suffix != ".gz":
                shutil.copyfile(raw_table, full_table)
            else:
                with gzip.open(raw_table, "rb") as src, full_table.open("wb") as dst:
                    shutil.copyfileobj(src, dst)
        if digest(source) != source_hash:
            raise OSError("BUSCO input mutated during evaluation")
        result = read_result(summary)
        atomic_json(receipt, {"key": key, "summary_sha256": digest(summary),
                              "full_table_sha256": digest(full_table), "reused_identical_before": False})
        return result


def collect_model_changes(root, pairs):
    """Bind gene rescue and accepted coding-path counts to the same publication."""
    from format_species_annotation.common import parse_gff_attributes
    from plot_gene_model_refinement import verified_json

    plan_hash = digest(root / "plan.json")
    request = json.loads((root / "plan.json").read_text())["request"]
    analysed = {r["species"] for r in pairs if r["refinement_status"] == "analysed"}
    if analysed != set(request["sources"]):
        raise ValueError("Model-count species differ from the refinement plan")
    stats, evidence = {}, {}
    batch = FreshDigestBatch()
    for pair in pairs:
        name = pair["species"]
        if pair["refinement_status"] == "not_analysed":
            stats[name] = {"refinement_status": "not_analysed", "prior_rescued_loci": None,
                           "accepted_repair_paths": None, "accepted_isoform_paths": None}
            continue
        source = request["sources"][name]
        if (Path(pair["before"]).resolve() != Path(source["fasta"]).resolve()
                or Path(pair["after"]).resolve() != (root / "effective/species_cds" / (name + ".fa")).resolve()):
            raise ValueError("Model counts belong to different BUSCO inputs: " + name)
        metadata = verified_json(root / "catalog" / name, "catalog_metadata.json", plan_hash)
        models = verified_json(root / "predictions" / name, "predictions.json", plan_hash)
        gff = Path(source["gff"])
        source_hash = batch.read([gff])[str(gff)]
        if source_hash != metadata["sources"]["gff"]["sha256"] or source_hash != request["files"][str(gff)]:
            raise ValueError("Rescue source annotation changed: " + name)
        rescued = set()
        opener = gzip.open if gff.suffix == ".gz" else open
        with opener(gff, "rt") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) == 9 and fields[1:3] == ["genegalleon_rescue", "gene"]:
                    identifiers = parse_gff_attributes(fields[8]).get("ID", [])
                    if len(identifiers) != 1:
                        raise ValueError("Rescued gene lacks a unique ID: " + name)
                    rescued.add(identifiers[0])
        accepted = [m for m in models if m["status"] == "accepted"]
        stats[name] = {
            "refinement_status": "analysed", "prior_rescued_loci": len(rescued),
            "accepted_repair_paths": sum(m["change_type"] == "model_revision" for m in accepted),
            "accepted_isoform_paths": sum(m["change_type"] == "isoform_addition" for m in accepted),
        }
        evidence[name] = {"source_gff_sha256": source_hash,
                          "catalog_receipt_sha256": digest(root / "catalog" / name / "receipt.json"),
                          "prediction_receipt_sha256": digest(root / "predictions" / name / "receipt.json")}
    batch.check()
    if digest(root / "plan.json") != plan_hash:
        raise ValueError("Refinement plan changed while collecting model counts")
    return {"schema": 1, "plan_sha256": plan_hash,
            "effective_receipt_sha256": digest(root / "effective/receipt.json"),
            "species": stats, "evidence": evidence,
            "units": "Prior missing-gene rescue counts unique source gene loci already present before refinement; repairs and additional isoforms count accepted coding paths, not unique loci."}


def plot_comparison(rows, output, model_changes=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.ticker import MaxNLocator

    if model_changes is not None:
        stats = model_changes["species"]
        if len(rows) != len(stats) or set(stats) != {r["species"] for r in rows}:
            raise ValueError("Model-count and BUSCO species membership differ")
        for row in rows:
            value = stats[row["species"]]
            if value["refinement_status"] != row["refinement_status"]:
                raise ValueError("Model-count analysis status differs from BUSCO")
            for key in ("prior_rescued_loci", "accepted_repair_paths", "accepted_isoform_paths"):
                count = value[key]
                if row["refinement_status"] == "not_analysed":
                    if count is not None:
                        raise ValueError("Unanalysed model counts must be unavailable, not zero")
                elif type(count) is not int or count < 0:
                    raise ValueError("Model counts must be nonnegative integers")
    extra = model_changes is not None
    margin_left = .20 if extra else .24
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11 if extra else 10, "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 5 if extra else 3, figsize=(25 if extra else 19, max(6, .43 * len(rows) + 2)),
                             gridspec_kw={"width_ratios": [1, 1, .65, .75, .9] if extra else [1, 1, .65],
                                          "wspace": .16 if extra else .12}, sharey=True)
    colors = STATUS_COLOURS
    for i, row in enumerate(rows):
        if row["refinement_status"] == "not_analysed":
            for ax in axes:
                ax.axhspan(i - .48, i + .48, color="#eef0f3", zorder=0)
        for j, phase in enumerate(("before", "after")):
            result, left = row[phase + "_result"], 0
            for category, color in zip(STATUS, colors, strict=True):
                value = 100 * result[category] / result["total"]
                axes[j].barh(i, value, left=left, color=color, height=.7)
                left += value
            axes[j].text(101, i, f'{100 * result["complete"] / result["total"]:.2f}', va="center", fontsize=9)
        delta = row["delta_complete_pp"]
        axes[2].barh(i, delta, color="#2a8d75" if delta >= 0 else "#b85250", height=.65)
        axes[2].text(.98, i, f'{row["delta_complete"]:+d} ({delta:+.2f} pp)',
                     transform=axes[2].get_yaxis_transform(), ha="right", va="center", fontsize=9)
        if extra:
            value = stats[row["species"]]
            if row["refinement_status"] == "not_analysed":
                for ax in axes[3:]:
                    ax.text(.03, i, "Not analysed", transform=ax.get_yaxis_transform(), va="center", color="#657585", fontsize=9)
            else:
                rescued = value["prior_rescued_loci"]
                repair, isoform = value["accepted_repair_paths"], value["accepted_isoform_paths"]
                axes[3].barh(i, rescued, color="#5275b5", height=.7)
                axes[3].annotate(str(rescued), (rescued, i), xytext=(4, 0), textcoords="offset points", va="center", fontsize=9)
                axes[4].barh(i, repair, color="#187d97", height=.7)
                axes[4].barh(i, isoform, left=repair, color="#d38b21", height=.7)
                axes[4].annotate(f"{repair} / {isoform}", (repair + isoform, i), xytext=(4, 0),
                                 textcoords="offset points", va="center", fontsize=9)
    labels = [r["species"].replace("_", " ") + ("  [not analysed]" if r["refinement_status"] == "not_analysed" else "") for r in rows]
    axes[0].set_yticks(range(len(rows)), labels)
    axes[0].invert_yaxis()
    titles = ["Before refinement", "After refinement", "Change in complete BUSCOs"]
    if extra:
        titles += ["Missing-gene rescue", "Accepted coding paths"]
    for ax, title in zip(axes, titles, strict=True):
        ax.set_title(title, fontweight="bold", pad=16, fontsize=11 if extra else 12)
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
    if extra:
        rescued_max = max((v["prior_rescued_loci"] or 0) for v in stats.values())
        paths_max = max((v["accepted_repair_paths"] or 0) + (v["accepted_isoform_paths"] or 0) for v in stats.values())
        axes[3].set_xlim(0, max(1, rescued_max) * 1.25)
        axes[4].set_xlim(0, max(1, paths_max) * 1.65)
        for ax in axes[3:]:
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
        axes[3].set_xlabel("Previously added gene loci")
        axes[4].set_xlabel("Accepted paths; labels: repair / isoform")
    fig.subplots_adjust(left=margin_left, right=.98, top=.90, bottom=.18 if extra else .14)
    fig.suptitle("Representative CDS completeness and gene-model improvement" if extra else
                 "Representative CDS completeness before and after refinement", x=margin_left, ha="left", y=.98, fontsize=17, fontweight="bold")
    identity = rows[0]["before_result"]
    fig.text(margin_left, .94, f'BUSCO {identity["busco_version"]}; {identity["lineage"]} ({identity["lineage_creation_date"]}); '
             f'n = {identity["total"]}; transcriptome mode; one representative per locus', fontsize=10)
    fig.legend([Patch(facecolor=c) for c in colors], STATUS_LABELS,
               loc="lower left", bbox_to_anchor=(margin_left, .09 if extra else .065), ncol=4, frameon=False)
    if extra:
        fig.legend([Patch(facecolor=c) for c in ("#5275b5", "#187d97", "#d38b21")],
                   ["Previously rescued gene loci", "Repair coding paths", "Additional isoform paths"],
                   loc="lower left", bbox_to_anchor=(.59, .09), ncol=3, frameon=False)
    note = "Grey rows: excluded from structural refinement; unchanged CDS are still evaluated by BUSCO.\n"
    note += "Before = refinement source CDS (including earlier rescued genes); after = selected DNA CDS, not all isoforms."
    if extra:
        note += "\nRescue counts are gene loci already in Before; repair / isoform counts are accepted paths and may share a locus."
    fig.text(margin_left, .025, note, fontsize=10)
    for suffix in ("png", "svg"):
        fig.savefig(output / ("busco_comparison." + suffix), dpi=180, facecolor="white")
    plt.close(fig)


def render_existing(report, root=None):
    """Redraw a historical evaluation without executing its predictor again."""
    value = json.loads((report / "busco_comparison.json").read_text())
    frozen = json.loads((report / "contract.json").read_text())
    pairs = {p["species"]: p for p in frozen["pairs"]}
    rows = value["species"]
    if (value["contract"] != frozen["contract"] or len(pairs) != len(frozen["pairs"])
            or len(rows) != len(pairs) or {r["species"] for r in rows} != set(pairs)):
        raise ValueError("Comparison contract or species membership changed")
    verified_tables = 0
    for row in rows:
        pair = pairs[row["species"]]
        for phase in ("before", "after"):
            directory = report / "runs" / pair["species"] / phase
            saved = json.loads((directory / "receipt.json").read_text())
            if (saved["key"]["contract"] != frozen["contract"]
                    or saved["key"]["source_sha256"] != digest(pair[phase])
                    or saved["summary_sha256"] != digest(directory / "summary.txt")
                    or read_result(directory / "summary.txt") != row[phase + "_result"]):
                raise ValueError("Comparison input or score changed")
            # Summary-only historical publications remain renderable as such.
            if "full_table_sha256" in saved:
                if digest(directory / "full_table.tsv") != saved["full_table_sha256"]:
                    raise ValueError("Comparison full table changed")
                verified_tables += 1
        if paired_result(pair, row["before_result"], row["after_result"]) != row:
            raise ValueError("Comparison delta changed")
    changes = collect_model_changes(root, rows) if root is not None else None
    if changes is not None:
        atomic_json(report / "model_change_summary.json", changes)
    plot_comparison(rows, report, changes)
    atomic_json(report / "rendering_provenance.json", {
        "comparison_sha256": digest(report / "busco_comparison.json"),
        "evaluation_contract": frozen["contract"], "renderer_sha256": digest(Path(__file__)),
        "verified_full_tables": verified_tables, "predictor_executed": False,
        "model_change_summary_sha256": digest(report / "model_change_summary.json") if changes is not None else None,
    })
    return rows


def evaluate(pairs, report, lineage, download_path, cpus=4, jobs=1, model_changes=None):
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
    if model_changes is not None:
        atomic_json(report / "model_change_summary.json", model_changes)
    plot_comparison(rows, report, model_changes)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Refinement publication; optional with --plot-only to include model-change counts")
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--cds-dir", type=Path, help="All dataset species, including unrefined CDS-only species")
    parser.add_argument("--lineage", type=Path, help="Frozen local lineage directory")
    parser.add_argument("--download-path", type=Path)
    parser.add_argument("--plot-only", action="store_true", help="Validate and redraw an existing comparison using its original evaluation contract")
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--jobs", type=int, default=1, help="Total CPU budget = jobs times cpus")
    args = parser.parse_args()
    if args.plot_only:
        if args.lineage or args.download_path or args.cds_dir:
            parser.error("--plot-only uses the saved evaluation; do not supply new inputs or lineage settings")
        render_existing(args.report.resolve(), args.output.resolve() if args.output else None)
        return
    if not args.output or not args.lineage or not args.download_path:
        parser.error("--output, --lineage and --download-path are required for evaluation")
    root, report = args.output.resolve(), args.report.resolve()
    if args.cpus < 1 or args.jobs < 1:
        parser.error("CPU and job counts must be positive")
    if root == report or root in report.parents or report in root.parents:
        parser.error("Report must be separate from the immutable refinement tree")
    pairs = input_pairs(root, args.cds_dir)
    evaluate(pairs, report, args.lineage.resolve(), args.download_path.resolve(), args.cpus, args.jobs,
             model_changes=collect_model_changes(root, pairs))


if __name__ == "__main__":
    main()

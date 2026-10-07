#!/usr/bin/env python3
"""Create a self-contained review of a verified refinement publication."""

from __future__ import annotations

import argparse
import base64
import csv
import html
import json
from collections import Counter, defaultdict
from pathlib import Path

from gene_model_refinement import verify_inputs
from gene_model_store import _connection as store_connection
from gene_model_store import load_locus
from input_generation_array_state import atomic_json, digest


def verified_json(directory, name, plan_hash):
    receipt = json.loads((directory / "receipt.json").read_text())
    if receipt["key"].get("plan") != plan_hash:
        raise ValueError("Review artifact belongs to another frozen plan")
    raw = (directory / name).read_bytes()
    import hashlib
    if hashlib.sha256(raw).hexdigest() != receipt["files"].get(name):
        raise ValueError("Review artifact changed: " + str(directory / name))
    return json.loads(raw)


def table(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def collect(root, max_loci=200, preferred_species="", cds_dir=None):
    """Count the entire publication; bound only the detailed locus gallery."""
    rows = verify_inputs(root / "effective" / "inputs.tsv")
    plan_hash = digest(root / "plan.json")
    changes = verified_json(root / "effective", "changes.json", plan_hash)
    selections = {(r["species"], r["gene_id"]): r for r in changes}
    stats, predictions = {}, {}
    detail_keys = set()
    phase_detail_keys = set()
    for row in rows:
        name = row["species"]
        metadata = verified_json(root / "catalog" / name, "catalog_metadata.json", plan_hash)
        models = verified_json(root / "predictions" / name, "predictions.json", plan_hash)
        predictions[name] = models
        selected = [r for r in changes if r["species"] == name]
        accepted = [r for r in models if r["status"] == "accepted"]
        stats[name] = {
            "refinement_status": "analysed", "refinement_reason": "",
            "source_loci": metadata["summary"]["loci"],
            "source_candidates": metadata["summary"]["candidates"],
            "accepted_repair_paths": sum(r["change_type"] == "model_revision" for r in accepted),
            "accepted_isoform_paths": sum(r["change_type"] == "isoform_addition" for r in accepted),
            "accepted_loci": len({r["gene_id"] for r in accepted}),
            "changed_representatives": sum(r["status"] == "conserved" for r in selected),
            "predicted_representatives": sum(r["selected_origin"] == "predicted" for r in selected),
            "proposed_paths": len(models) - len(accepted),
            "selection_status": dict(Counter(r["status"] for r in selected)),
            "proposal_reasons": dict(Counter(p for r in models if r["status"] != "accepted" for p in r["problems"])),
            "source_fasta_mapping": metadata["summary"]["fasta_mapping"],
        }
        detail_keys.update((name, r["gene_id"]) for r in accepted)
        detail_keys.update((name, r["gene_id"]) for r in selected if r["status"] == "conserved")
    for r in table(root / "effective" / "effective_exclusions.tsv"):
        stats[r["species"]].setdefault("effective_exclusions", 0)
        stats[r["species"]]["effective_exclusions"] += 1
    for r in table(root / "effective" / "translation_admission.tsv"):
        quality = json.loads(r["quality"])
        if r["status"] == "included" and quality.get("phase_inferred"):
            stats[r["species"]].setdefault("phase_resolved_representatives", 0)
            stats[r["species"]]["phase_resolved_representatives"] += 1
            if quality.get("phase_inference_evidence") == "complete_genomic_cds_and_unique_source_coding_path":
                stats[r["species"]].setdefault("coding_path_phase_resolved_representatives", 0)
                stats[r["species"]]["coding_path_phase_resolved_representatives"] += 1
                phase_detail_keys.add((r["species"], r["gene_id"]))
        if r["status"] == "excluded":
            stats[r["species"]].setdefault("translation_withheld", 0)
            stats[r["species"]]["translation_withheld"] += 1
    # Include informative rejected proposals when the changed-locus gallery fits.
    ordered = sorted(detail_keys, key=lambda key: (key[0] != preferred_species, key))
    for key in sorted(phase_detail_keys, key=lambda key: (key[0] != preferred_species, key)):
        if len(ordered) >= max_loci:
            break
        if key not in ordered:
            ordered.append(key)
    proposals = sorted(
        {(n, r["gene_id"]) for n, models in predictions.items() for r in models},
        key=lambda key: (key[0] != preferred_species, key),
    )
    for key in proposals:
        if len(ordered) >= max_loci:
            break
        if key not in ordered:
            ordered.append(key)
    included = set(ordered[:max_loci])
    grouped = defaultdict(list)
    for name, models in predictions.items():
        for r in models:
            if (name, r["gene_id"]) in included:
                grouped[name, r["gene_id"]].append(r)
    details = []
    db = root / "catalog_index_final" / "loci.sqlite3"
    db_receipt = json.loads((db.parent / "receipt.json").read_text())
    if db_receipt["key"].get("plan") != plan_hash or digest(db) != db_receipt["files"].get(db.name):
        raise ValueError("Review locus database changed")
    before_db = db.stat()
    with store_connection(db) as connection:
        for key in sorted(included):
            name, gene_id = key
            locus = load_locus(connection, name, gene_id)
            chosen = selections[key]
            candidates = [
                {k: c.get(k) for k in ("candidate_id", "source_transcript_id", "origin", "blocks", "source_blocks", "quality", "support")}
                | {"cds_length": len(c["cds"]), "protein_length": len(c["protein"])}
                for c in locus["candidates"]
            ]
            proposals_at_locus = [
                {k: r.get(k) for k in ("status", "change_type", "evidence_class", "donors", "rna_paths", "problems", "alignments")}
                | {"candidate_id": r["candidate"]["candidate_id"], "blocks": r["candidate"]["blocks"],
                   "cds_length": len(r["candidate"]["cds"]), "protein_length": len(r["candidate"]["protein"])}
                for r in grouped[key]
            ]
            details.append({
                "species": name, "gene_id": locus["gene_id"], "seqid": locus["seqid"], "strand": locus["strand"],
                "baseline_id": locus.get("source_baseline_candidate_id") or locus.get("source_baseline_coding_candidate_id"),
                "baseline_basis": locus.get("source_baseline_basis", "source_transcript_identity"), "selection": chosen,
                "candidates": candidates, "predictions": proposals_at_locus,
            })
    after_db = db.stat()
    if (before_db.st_size, before_db.st_mtime_ns, before_db.st_ctime_ns) != (
            after_db.st_size, after_db.st_mtime_ns, after_db.st_ctime_ns):
        raise OSError("Review locus database mutated while reading")
    performance = {}
    for path in sorted(root.glob("*/performance.json")) + sorted(root.glob("*/*/performance.json")):
        performance[str(path.parent.relative_to(root))] = verified_json(path.parent, path.name, plan_hash)
    if cds_dir:
        from gene_model_refinement_busco import input_pairs
        for pair in input_pairs(root, cds_dir):
            if pair["species"] not in stats:
                stats[pair["species"]] = {"refinement_status": "not_analysed", "refinement_reason": pair["reason"]}
    stats = dict(sorted(stats.items()))
    return {
        "schema": 1, "plan_sha256": plan_hash,
        "effective_receipt_sha256": digest(root / "effective" / "receipt.json"),
        "species": stats, "details": details, "changed_loci_available": len(detail_keys),
        "coding_path_phase_loci_available": len(phase_detail_keys),
        "gallery_limit": max_loci, "performance": performance,
        "limits": "Counts cover all species and paths. The locus gallery is bounded. Ranking margins are not probabilities. RNA chain support does not establish translation initiation. Accepted paths and selected representatives are distinct decisions.",
    }


def plot_summary(data, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    names = list(data["species"])
    stats = [data["species"][n] for n in names]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 4, figsize=(18, max(6, len(names) * .32)),
                             sharey=True, gridspec_kw={"wspace": .15})
    ys = list(range(len(names)))
    repairs = [s.get("accepted_repair_paths", 0) for s in stats]
    isoforms = [s.get("accepted_isoform_paths", 0) for s in stats]
    axes[0].barh(ys, repairs, color="#187d97", label="Repair coding paths")
    axes[0].barh(ys, isoforms, left=repairs, color="#d38b21", label="Additional isoform paths")
    axes[1].barh(ys, [s.get("changed_representatives", 0) for s in stats], color="#5275b5")
    axes[2].barh(ys, [s.get("phase_resolved_representatives", 0) for s in stats], color="#32856b")
    axes[3].barh(ys, [s.get("effective_exclusions", 0) for s in stats], color="#a14c57", label="Source/structure mismatch")
    axes[3].barh(ys, [s.get("translation_withheld", 0) - s.get("effective_exclusions", 0) for s in stats],
                 left=[s.get("effective_exclusions", 0) for s in stats], color="#bbb1bc", label="Other translation exclusions")
    for ax, title in zip(axes, ("Accepted coding paths", "Changed representatives", "Phase-resolved representatives", "Withheld from coding analysis"), strict=True):
        ax.set_title(title, fontweight="bold")
        ax.set_xlabel("Count")
        ax.xaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", alpha=.15)
        ax.set_axisbelow(True)
        ax.set_xlim(left=0)
        for i, s in enumerate(stats):
            if s.get("refinement_status") == "not_analysed":
                ax.axhspan(i - .5, i + .5, color="#eef0f3", zorder=0)
                ax.text(.02, i, "Not analysed", transform=ax.get_yaxis_transform(), va="center", color="#647383", fontsize=9)
    axes[0].set_yticks(ys, [n.replace("_", " ") + (" [not analysed]" if stats[i].get("refinement_status") == "not_analysed" else "") for i, n in enumerate(names)], fontstyle="italic")
    axes[0].invert_yaxis()
    handles, labels = axes[0].get_legend_handles_labels()
    other_handles, other_labels = axes[3].get_legend_handles_labels()
    fig.legend(handles + other_handles, labels + other_labels, loc="lower left", bbox_to_anchor=(.25, .06), ncol=2, frameon=False, fontsize=9)
    fig.suptitle("Synteny-guided gene-model refinement", x=.32, y=.995, ha="left", fontweight="bold", fontsize=16)
    fig.subplots_adjust(left=.25, right=.99, top=.92, bottom=.20)
    fig.text(.25, .02, "Grey rows: no matching genome/GFF in the plan; structural refinement not analysed, CDS retained unchanged.\n"
             "Coding paths and genes are counted separately. Ranking scores are not confidence probabilities.", fontsize=9)
    for suffix in ("png", "svg"):
        fig.savefig(output / ("summary." + suffix), dpi=180, facecolor="white")
    plt.close(fig)


def locus_svg(locus):
    """Draw genomic coding blocks with true genomic spacing and phase labels."""
    paths = [(c["candidate_id"], c["blocks"], c["origin"], c["cds_length"]) for c in locus["candidates"]]
    seen = {r[0] for r in paths}
    paths.extend((p["candidate_id"], p["blocks"], p["status"], p["cds_length"])
                 for p in locus["predictions"] if p["candidate_id"] not in seen)
    blocks = [b for _, bs, _, _ in paths for b in bs]
    if not blocks:
        return "<p>No reconstructible coding blocks.</p>"
    lo, hi = min(b[0] for b in blocks), max(b[1] for b in blocks)
    left, width = 335, 650
    def scale(x):
        return left + width * (x - lo) / max(1, hi - lo)
    chosen = locus["selection"]["candidate_id"]
    phase_resolved = {c["candidate_id"] for c in locus["candidates"] if c.get("quality", {}).get("phase_inferred")}
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1020 {75 + 38 * len(paths)}" role="img" aria-label="Coding exon structure">']
    for i, (identifier, bs, origin, length) in enumerate(paths):
        y = 35 + i * 38
        color = "#5275b5" if identifier == chosen else "#187d97" if origin == "original" else "#d38b21" if origin != "proposal" else "#a14c57"
        label = identifier.removeprefix(locus["species"] + "_")
        short_label = label if len(label) <= 38 else label[:35] + "…"
        tags = [origin, str(length) + " nt"]
        if identifier in phase_resolved:
            tags.append("phase resolved")
        if identifier == chosen:
            tags.append("selected")
        if identifier == locus.get("baseline_id"):
            tags.append("source coding path" if locus.get("baseline_basis") == "matching_coding_path_transcript_identity_ambiguous" else "source")
        parts.append(f'<text x="5" y="{y}" font-size="11"><title>{html.escape(identifier)}</title>{html.escape(short_label)}</text>')
        parts.append(f'<text x="5" y="{y + 14}" font-size="9">{html.escape("; ".join(tags))}</text>')
        if bs:
            a, b = min(v[0] for v in bs), max(v[1] for v in bs)
            parts.append(f'<line x1="{scale(a):.2f}" x2="{scale(b):.2f}" y1="{y}" y2="{y}" stroke="{color}"/>')
            tip = scale(b) + 8 if locus["strand"] == "+" else scale(a) - 8
            back = tip - 6 if locus["strand"] == "+" else tip + 6
            parts.append(f'<path d="M {back:.2f},{y-4} L {tip:.2f},{y} L {back:.2f},{y+4}" fill="none" stroke="{color}"/>')
        for a, b, phase in bs:
            parts.append(f'<rect x="{scale(a):.2f}" y="{y - 8}" width="{max(.8, scale(b)-scale(a)):.2f}" height="16" fill="{color}"><title>{a+1:,}–{b:,}; phase {phase}</title></rect>')
    y = 40 + len(paths) * 38
    parts.append(f'<text x="{left}" y="{y}" font-size="11">{html.escape(locus["seqid"])}:{lo+1:,}–{hi:,} ({locus["strand"]})</text></svg>')
    return "".join(parts)


def write_review(data, output):
    image = base64.b64encode((output / "summary.png").read_bytes()).decode()
    cards = []
    for locus in data["details"]:
        choice = locus["selection"]
        rows = []
        for p in locus["predictions"]:
            rows.append("<tr>" + "".join("<td>" + html.escape(str(v)) + "</td>" for v in (
                p["candidate_id"], p["status"], p["change_type"], p["evidence_class"],
                len(p["donors"]), len(p["rna_paths"]), ", ".join(p["problems"]))) + "</tr>")
        details = json.dumps({"selection": choice, "candidate_quality": [
            {"candidate_id": c["candidate_id"], "quality": c["quality"], "support": c["support"], "source_blocks": c.get("source_blocks")}
            for c in locus["candidates"]], "prediction_alignments": locus["predictions"]}, indent=2)
        cards.append(
            f'<article data-search="{html.escape((locus["species"]+" "+locus["gene_id"]+" "+choice["status"]).lower(), quote=True)}">'
            f'<h2>{html.escape(locus["gene_id"])}</h2><p>{html.escape(choice["status"])}: {html.escape(choice["reason"])}; '
            f'margin {html.escape(str(choice.get("margin")))}</p>{locus_svg(locus)}'
            '<div class="scroll"><table><thead><tr><th>Candidate</th><th>Decision</th><th>Change</th><th>Evidence</th>'
            '<th>Donor species</th><th>RNA path records</th><th>Reasons</th></tr></thead><tbody>'
            + "".join(rows) + '</tbody></table></div><details><summary>Quality, selection and alignment evidence</summary>'
            f'<pre>{html.escape(details)}</pre></details></article>')
    totals = Counter()
    for stats in data["species"].values():
        for key, value in stats.items():
            if isinstance(value, int):
                totals[key] += value
    html_text = """<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>Gene-model refinement review</title><style>
body{font:16px system-ui,sans-serif;margin:30px auto;padding:0 24px;max-width:1450px;color:#182b38;background:#fafafa}
h1{font-size:30px}h2{font-size:19px;overflow-wrap:anywhere}img,svg{width:100%;height:auto}article{background:white;padding:20px;margin:24px 0;border:1px solid #ccd6df;border-radius:8px}
input{padding:12px;width:min(95%,700px);font:inherit}table{border-collapse:collapse;font-size:12px}td,th{padding:9px;border:1px solid #ddd;text-align:left;overflow-wrap:anywhere}pre{white-space:pre-wrap;font-size:12px}.scroll{overflow:auto}.muted{color:#506674}
.metrics{display:flex;flex-wrap:wrap;gap:20px;margin:24px 0}.metrics div{padding:16px 22px;background:white;border:1px solid #ccd6df;border-radius:8px}.metrics strong{font-size:28px;color:#187d97}
</style><h1>Synteny-guided gene-model refinement</h1>"""
    html_text += "<p>" + html.escape(data["limits"]) + "</p>"
    metrics = [("Repair coding paths", "accepted_repair_paths"), ("Additional isoform paths", "accepted_isoform_paths"),
               ("Changed representatives", "changed_representatives"), ("Phase-resolved representatives", "phase_resolved_representatives"),
               ("Accepted loci", "accepted_loci")]
    html_text += '<div class="metrics">' + "".join(
        f'<div><strong>{totals[key]:,}</strong><br>{label}</div>' for label, key in metrics) + '</div>'
    html_text += f'<img alt="Species-wide coding path additions, representative changes and exclusions" src="data:image/png;base64,{image}">'
    html_text += '<details><summary>Complete species counts</summary><div class="scroll"><table><thead><tr>'
    columns = [("Species", ""), ("Source loci", "source_loci"), ("Source coding paths", "source_candidates"),
               ("Repair paths", "accepted_repair_paths"), ("Added isoform paths", "accepted_isoform_paths"),
               ("Accepted loci", "accepted_loci"), ("Changed representatives", "changed_representatives"),
               ("Predicted representatives", "predicted_representatives"), ("Proposal paths", "proposed_paths"),
               ("Phase-resolved representatives", "phase_resolved_representatives"),
               ("Sequence/GFF exclusions", "effective_exclusions"), ("Translations withheld", "translation_withheld")]
    html_text += '<th>Refinement status / reason</th>' + "".join(f'<th>{label}</th>' for label, _ in columns) + '</tr></thead><tbody>'
    for name, stats in data["species"].items():
        status = stats.get("refinement_status", "analysed")
        html_text += '<tr><td>' + html.escape(status + (": " + stats["refinement_reason"] if stats.get("refinement_reason") else "")) + '</td>'
        html_text += "".join('<td>' + (html.escape(name.replace("_", " ")) if not key else "Not analysed" if status == "not_analysed" else f'{stats.get(key, 0):,}') + '</td>'
                                    for _, key in columns) + '</tr>'
    html_text += '</tbody></table></div></details>'
    html_text += f'<p>Gallery: {len(cards)} loci; all changed/accepted loci available: {data["changed_loci_available"]}. '
    html_text += f'Source coding-path phase resolutions available: {data["coding_path_phase_loci_available"]}. '
    html_text += "Phase resolutions and rejected proposals fill unused gallery slots. Exon plots use genomic spacing and 1-based display coordinates.</p>"
    html_text += '<input id="search" aria-label="Filter gallery" placeholder="Filter by species, gene ID or selection decision"><span id="count"></span>'
    html_text += "".join(cards) + """<script>
const q=document.getElementById('search'),items=[...document.querySelectorAll('article')];
function filter(){let n=0;for(const a of items){a.hidden=!a.dataset.search.includes(q.value.toLowerCase());if(!a.hidden)n++}
document.getElementById('count').textContent=' '+n+' visible loci'}q.addEventListener('input',filter);filter();
</script>"""
    html_text += '<p class="muted">Plan SHA256: ' + html.escape(data["plan_sha256"]) + "; effective receipt SHA256: " + html.escape(data["effective_receipt_sha256"]) + "</p></html>"
    (output / "review.html").write_text(html_text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="Completed native refinement directory")
    parser.add_argument("--report", type=Path, required=True, help="Separate review output directory")
    parser.add_argument("--max-loci", type=int, default=200)
    parser.add_argument("--preferred-species", default="")
    parser.add_argument("--cds-dir", type=Path, help="Include all dataset CDS species; mark species absent from the native plan as not analysed")
    args = parser.parse_args()
    if args.max_loci < 1:
        parser.error("--max-loci must be positive")
    root, report = args.output.resolve(), args.report.resolve()
    if report == root or report in root.parents or root in report.parents:
        parser.error("--report must be separate from the immutable refinement tree")
    data = collect(root, args.max_loci, args.preferred_species, args.cds_dir)
    report.mkdir(parents=True, exist_ok=True)
    plot_summary(data, report)
    write_review(data, report)
    atomic_json(report / "review_data.json", data)
    print(json.dumps({"species": len(data["species"]), "gallery_loci": len(data["details"]), "report": str(report)}))


if __name__ == "__main__":
    main()

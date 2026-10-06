"""Render native gene trees for category-1 events with bilateral scaffold background.

This presentation subset does not modify the input focused event cohort. Background
support is evaluated per retained gene, never from family/scaffold averages.
"""

import csv
import hashlib
import io
import json
import logging
import math
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

from gene_family_output_store import GeneFamilyOutputStore, read_only_observation

PROFILE = {"rank": "class", "minimum_classified_units": 10,
           "minimum_classified_fraction": 0.5, "minimum_compatible_fraction": 0.9,
           "minimum_ufboot": 90, "minimum_passing_genes_per_side": 1}
INDEX_FIELDS = ["orthogroup", "status", "reason", "event_count", "node_count",
                "pdf", "annotated_stat_branch", "source_stat_branch_sha256"]
EVENT_FIELDS = ["event_id", "orthogroup", "gene_tree_branch_id", "gene_tree_node",
                "generax_transfer", "status", "reason", "support_generax_ufboot",
                "donor_passing_gene_count", "recipient_passing_gene_count", "plot_label"]


def number(value):
    try:
        result = float(value)
    except (ValueError, TypeError):
        return None
    return result if math.isfinite(result) else None


def background_supported(row):
    if str(row.get("eligible_for_context", "")).lower() not in {"true", "1"}:
        return False
    if row.get("host_scaffold_status") != "measured" or not row.get("host_scaffold_id"):
        return False
    prefix = "host_scaffold_background_class_"
    total, compatible, incompatible, unresolved = [number(row.get(prefix + suffix)) for suffix in
                                                   ("total_count", "compatible_count",
                                                    "incompatible_count", "unresolved_count")]
    if any(value is None for value in (total, compatible, incompatible, unresolved)):
        return False
    if any(value < 0 or not value.is_integer() for value in (total, compatible, incompatible, unresolved)):
        raise ValueError("Invalid per-gene scaffold background counts")
    if compatible + incompatible + unresolved != total:
        raise ValueError("Inconsistent per-gene scaffold background counts")
    classified = compatible + incompatible
    coverage = classified / total if total else None
    fraction = compatible / classified if classified else None
    for suffix, expected in (("classified_fraction", coverage), ("compatible_fraction", fraction)):
        measured = number(row.get(prefix + suffix))
        if (measured is None) != (expected is None) or (
                measured is not None and abs(measured - expected) > 1e-9):
            raise ValueError("Inconsistent per-gene scaffold background fraction")
    return (classified >= PROFILE["minimum_classified_units"]
            and coverage >= PROFILE["minimum_classified_fraction"]
            and fraction >= PROFILE["minimum_compatible_fraction"])


def identity(event, project_name, native_name):
    values = {str(event[name]) for name in (project_name, native_name) if event.get(name) not in (None, "")}
    if len(values) != 1:
        raise ValueError(f"Missing or conflicting event identity: {project_name}/{native_name}")
    return values.pop()


def annotate(stat_rows, events, links):
    """Match exact family/branch/node/token and audit every requested event."""
    branches = {row["branch_id"]: row for row in stat_rows}
    if len(branches) != len(stat_rows) or len({row["node_name"] for row in stat_rows}) != len(stat_rows):
        raise ValueError("Duplicate gene-tree branch or node identity")
    linked = defaultdict(list)
    tips = {row['node_name'] for row in stat_rows
            if all(row.get(name) == '-999' for name in ('child1', 'child2'))}
    for row in links:
        linked[(row["event_id"], row["side"])].append(row)
    audits, selected, genes = [], defaultdict(list), set()
    for event in sorted(events, key=lambda row: row["event_id"]):
        branch_id = identity(event, "gene_tree_branch_id", "branch_id")
        node = identity(event, "gene_tree_node", "node_name")
        row = branches.get(branch_id)
        reason = ""
        support = None
        sides = {}
        if row is None or row["node_name"] != node:
            reason = "gene_tree_branch_node_unmapped"
        elif any(row.get(name) in ("", "-999", None) for name in ("child1", "child2")):
            reason = "terminal_gene_tree_branch"
        else:
            tokens = re.split(r"\s*[;,|]\s*", row.get("generax_transfer", ""))
            position = number(event.get("event_index"))
            if (position is None or not position.is_integer() or not 1 <= position <= len(tokens)
                    or tokens[int(position) - 1] != event["generax_transfer"]):
                reason = "transfer_token_or_event_index_unmapped"
            support = number(row.get("support_generax_ufboot"))
            if not reason and support is None:
                reason = "generax_ufboot_unavailable"
            elif not reason and not 0 <= support <= 100:
                raise ValueError("Gene-tree UFBoot must be in [0, 100]")
            elif not reason and support < PROFILE["minimum_ufboot"]:
                reason = "generax_ufboot_below_threshold"
            if not reason:
                recorded = number(event.get("support_used", event.get("support_generax_ufboot")))
                if recorded is not None and abs(recorded - support) > 1e-9:
                    raise ValueError("Focused event support disagrees with its exact gene-tree branch")
        if not reason:
            for side in ("donor", "recipient"):
                sides[side] = [link for link in linked[(event["event_id"], side)] if background_supported(link)]
            if any(not sides[side] for side in ("donor", "recipient")):
                reason = "bilateral_class_background_not_supported_or_unavailable"
            elif any(link['gene_id'] not in tips or link['orthogroup'] != event['orthogroup']
                     for side in ('donor', 'recipient') for link in sides[side]):
                reason = "supported_context_gene_unmapped_to_family_tip"
        label = ""
        if not reason:
            label = f"HGT{sum(len(group) for group in selected.values()) + 1}"
            selected[branch_id].append((event, label, support))
            genes.update(link["gene_id"] for link in sides["recipient"])
        audits.append(dict(event_id=event["event_id"], orthogroup=event["orthogroup"],
                           gene_tree_branch_id=branch_id, gene_tree_node=node,
                           generax_transfer=event["generax_transfer"], status="selected" if not reason else "withheld",
                           reason=reason, support_generax_ufboot="" if support is None else support,
                           donor_passing_gene_count=len(sides['donor']) if 'donor' in sides else '',
                           recipient_passing_gene_count=len(sides['recipient']) if 'recipient' in sides else '', plot_label=label))
    output = []
    for row in stat_rows:
        matched = selected.get(row["branch_id"], [])
        output.append(dict(row, hgtfocus_event_count=len(matched),
                           hgtfocus_event_ids="; ".join(e["event_id"] for e, _, _ in matched),
                           hgtfocus_node_label="; ".join(f"{label} UF={support:g}" for _, label, support in matched),
                           hgtfocus_recipient_flag=int(row["node_name"] in genes),
                           hgtfocus_tip_status="Supported recipient" if row["node_name"] in genes else ""))
    return output, audits


def write(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def export_gene_trees(directory, events, links, family_root, renderer=None, gff_root=''):
    """Export one native PDF per family; unavailable mappings remain in the audit."""
    csv.field_size_limit(100_000_000)
    directory.mkdir(parents=True)
    helper = Path(__file__).resolve().parent
    families = defaultdict(list)
    for event in events:
        family = event["orthogroup"]
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", family) or family in {".", ".."}:
            raise ValueError("Unsafe orthogroup identifier")
        families[family].append(event)
    from focus_hgt_context import GenomeCoordinates, render_context
    from gene_tree_plot_config import replay

    coordinates = GenomeCoordinates(gff_root)
    index, audit, sources, context_audit, configurations = [], [], {}, [], {}
    with read_only_observation():
        store = GeneFamilyOutputStore(family_root)
        with store.read_snapshot():
            for family, group in sorted(families.items()):
                try:
                    with store.open_binary("stat_branch", family + "_stat.branch.tsv") as handle:
                        raw = handle.read()
                except FileNotFoundError:
                    for event in group:
                        audit.append(dict(event_id=event["event_id"], orthogroup=family,
                                          gene_tree_branch_id=identity(event, "gene_tree_branch_id", "branch_id"),
                                          gene_tree_node=identity(event, "gene_tree_node", "node_name"),
                                          generax_transfer=event["generax_transfer"], status="withheld",
                                          reason="stat_branch_unavailable", support_generax_ufboot="",
                                          donor_passing_gene_count="", recipient_passing_gene_count="", plot_label=""))
                    index.append(dict(orthogroup=family, status="withheld", reason="stat_branch_unavailable",
                                      event_count=0, node_count=0, pdf="", annotated_stat_branch="",
                                      source_stat_branch_sha256=""))
                    continue
                sha = hashlib.sha256(raw).hexdigest()
                sources["stat_branch/" + family + "_stat.branch.tsv"] = sha
                rows = list(csv.DictReader(io.StringIO(raw.decode()), delimiter="\t"))
                annotated, checks = annotate(rows, group, [link for link in links if link["orthogroup"] == family])
                audit.extend(checks)
                passed = [row for row in checks if row["status"] == "selected"]
                table = directory / "tree_plot_input" / (family + "_focused_stat.branch.tsv")
                write(table, list(annotated[0]), annotated)
                pdf = directory / (family + "_focused_hgt_tree_plot.pdf")
                if passed:
                    logging.info('Focused gene tree %s: %d supported events', family, len(passed))
                    if renderer is not None:
                        renderer(table, pdf)
                    else:
                        with tempfile.TemporaryDirectory(prefix="hgt-focus-tree-") as tmp:
                            materialized = Path(tmp) / 'family_inputs'
                            spec = replay(store, family, rows, materialized, sources)
                            configurations[family] = spec
                            command = [
                                "Rscript", str(helper / "stat_branch2tree_plot.r"), f"--stat_branch={table.resolve()}",
                                *spec['arguments'],
                            ]
                            try:
                                import os
                                environment = dict(os.environ, TREEVIS_SPECIES_PARSER=spec['species_label_parser'])
                                subprocess.run(command, cwd=tmp, check=True, stdout=subprocess.PIPE,
                                               stderr=subprocess.STDOUT, text=True, env=environment)
                            except subprocess.CalledProcessError as exc:
                                raise RuntimeError(f"Native focused gene-tree rendering failed for {family}:\n{exc.stdout}") from exc
                            source = Path(tmp) / "stat_branch2tree_plot.pdf"
                            if not source.is_file() or not source.read_bytes().startswith(b"%PDF"):
                                raise ValueError("Native focused gene-tree renderer did not produce a PDF")
                            context = Path(tmp) / 'context.pdf'
                            selected_ids = {r['event_id'] for r in passed}
                            context_audit += render_context(context, rows, [e for e in group if e['event_id'] in selected_ids],
                                                            links, coordinates)
                            from pypdf import PdfReader, PdfWriter
                            if len(PdfReader(source).pages) != 1 or len(PdfReader(context).pages) != 1:
                                raise ValueError('Focused gene-tree PDF must have exactly one tree and one context page')
                            writer = PdfWriter()
                            writer.append(str(source))
                            writer.append(str(context))
                            with pdf.open('wb') as handle:
                                writer.write(handle)
                index.append(dict(orthogroup=family, status="rendered" if passed else "withheld",
                                  reason="" if passed else "no_qualifying_mapped_event", event_count=len(passed),
                                  node_count=len({row["gene_tree_branch_id"] for row in passed}),
                                  pdf=pdf.name if passed else "", annotated_stat_branch=str(table.relative_to(directory)),
                                  source_stat_branch_sha256=sha))
            for logical, expected in sources.items():
                subdir, name = logical.split("/", 1)
                with store.open_binary(subdir, name) as handle:
                    if hashlib.sha256(handle.read()).hexdigest() != expected:
                        raise ValueError("Gene-tree input changed during focused rendering")
    write(directory / "index.tsv", INDEX_FIELDS, index)
    write(directory / "event_node_audit.tsv", EVENT_FIELDS, audit)
    if context_audit:
        write(directory / 'context_gene_audit.tsv', list(context_audit[0]), context_audit)
    (directory / 'renderer_settings.json').write_text(json.dumps(configurations, indent=2) + '\n')
    coordinates.verify()
    (directory / "README.txt").write_text(
        "Native GeneGalleon gene trees for observed category-1 recipients\n\n"
        "Orange diamonds and HGT labels mark exact gene-tree transfer nodes, including internal nodes.\n"
        "UF labels are the matched branch's support_generax_ufboot (>=90 inclusive).\n"
        "At least one retained event-linked gene on each side must have candidate-free class background\n"
        "with >=10 classified units, >=50% classification coverage and >=90% host compatibility.\n"
        "Orange recipient tips are the genes that individually pass that background check.\n"
        "Page 1 replays gg_gene_evolution panels and saved settings, including domain, gene structure and alignment.\n"
        "Missing optional measurements are not invented. Renderer settings and input availability are recorded.\n"
        "Page 2 gives exact gene-tree paths and representative donor/recipient GFF neighborhoods.\n"
        "All genomic tracks share a linear kb axis centered on their focal-gene midpoint, without intron compression.\n"
        "CDS blocks are coding exons; UTR blocks are shown when recorded; unavailable structures stay unconfirmed.\n"
        "This is whole-scaffold context, not conserved gene order or proof of physical integration.\n"
        "All event IDs, branch/node IDs and selection/withholding reasons are in event_node_audit.tsv.\n"
        "No sequence or phylogenetic analysis is run. The parent focused event tables are unchanged.\n")
    return dict(profile=PROFILE, family_source_sha256=sources, gff_source_sha256=coordinates.sources,
                rendered_family_count=sum(row["status"] == "rendered" for row in index),
                selected_event_count=sum(row["status"] == "selected" for row in audit))

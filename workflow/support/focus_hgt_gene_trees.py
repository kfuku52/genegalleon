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
           "minimum_ufboot": None, "minimum_passing_genes_per_side": 1}
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


def validate_link_identity(event, link):
    """Match native and project branch aliases without borrowing another event."""
    for field in ("event_id", "orthogroup", "generax_transfer", "event_index",
                  "generax_donor_node", "generax_recipient_node"):
        if field in link and field in event and link[field] != event[field]:
            raise ValueError("Event-gene identity mismatch: " + field)
    for aliases in (("branch_id", "gene_tree_branch_id"), ("node_name", "gene_tree_node")):
        values = [{str(row[name]) for name in aliases if row.get(name) not in (None, "")}
                  for row in (event, link)]
        if any(len(group) > 1 for group in values) or all(values) and values[0] != values[1]:
            raise ValueError("Event-gene identity mismatch: " + "/".join(aliases))
    if not str(link.get("gene_id", "")).strip():
        raise ValueError("Empty event-gene identity")


def validated_event_links(events, links):
    """Validate the requested cohort and its exact links, ignoring other events."""
    ids = {event['event_id']: event for event in events}
    if len(ids) != len(events) or any(not str(value).strip() for value in ids):
        raise ValueError('Duplicate or empty transfer event_id')
    selected, identities = [], set()
    for link in links:
        if link['event_id'] not in ids:
            continue
        key = link['event_id'], link['side'], link['gene_id']
        if link['side'] not in {'donor', 'recipient'} or key in identities:
            raise ValueError('Duplicate or invalid event-gene link')
        identities.add(key)
        validate_link_identity(ids[link['event_id']], link)
        if str(link.get('eligible_for_context', '')).lower() in {'true', '1'} \
                and link.get('lineage_status', 'retained') != 'retained':
            raise ValueError('Eligible context gene does not have a retained transfer lineage')
        selected.append(link)
    return selected


def annotate(stat_rows, events, links, minimum_ufboot=None):
    """Match exact family/branch/node/token and audit every requested event."""
    links = validated_event_links(events, links)
    if len({event['orthogroup'] for event in events}) > 1:
        raise ValueError('Gene-tree annotation requires events from a single orthogroup')
    if minimum_ufboot is not None:
        minimum_ufboot = number(minimum_ufboot)
        if minimum_ufboot is None or not 0 <= minimum_ufboot <= 100:
            raise ValueError("Gene-tree UFB threshold must be a finite number in [0, 100]")
    branches = {row["branch_id"]: row for row in stat_rows}
    if len(branches) != len(stat_rows) or len({row["node_name"] for row in stat_rows}) != len(stat_rows):
        raise ValueError("Duplicate gene-tree branch or node identity")
    linked = defaultdict(list)
    tips = {row['node_name'] for row in stat_rows
            if all(row.get(name) == '-999' for name in ('child1', 'child2'))}
    for row in links:
        linked[(row["event_id"], row["side"])].append(row)
    audits, selected, genes = [], defaultdict(list), set()
    donor_genes, roles = set(), defaultdict(set)
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
            if not reason and support is None and minimum_ufboot is not None:
                reason = "generax_ufboot_unavailable"
            elif not reason and support is not None and not 0 <= support <= 100:
                raise ValueError("Gene-tree UFBoot must be in [0, 100]")
            elif not reason and minimum_ufboot is not None and support < minimum_ufboot:
                reason = "generax_ufboot_below_threshold"
            if not reason:
                recorded = number(event.get("support_used", event.get("support_generax_ufboot")))
                if recorded is not None and (support is None or abs(recorded - support) > 1e-9):
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
            donor_genes.update(link['gene_id'] for link in sides['donor'])
            for side in ('donor', 'recipient'):
                roles[side].update(link['gene_id'] for link in linked[event['event_id'], side]
                                   if str(link.get('eligible_for_context', '')).lower() in {'true', '1'}
                                   and link['gene_id'] in tips and link['orthogroup'] == event['orthogroup'])
        audits.append(dict(event_id=event["event_id"], orthogroup=event["orthogroup"],
                           gene_tree_branch_id=branch_id, gene_tree_node=node,
                           generax_transfer=event["generax_transfer"], status="selected" if not reason else "withheld",
                           reason=reason, support_generax_ufboot="" if support is None else support,
                           donor_passing_gene_count=len(sides['donor']) if 'donor' in sides else '',
                           recipient_passing_gene_count=len(sides['recipient']) if 'recipient' in sides else '', plot_label=label))
    output = []
    for row in stat_rows:
        matched = selected.get(row["branch_id"], [])
        name = row['node_name']
        status = []
        for side, passing in [('donor', donor_genes), ('recipient', genes)]:
            if name in roles[side]:
                status.append(('Scaffold-supported ' if name in passing else 'Scaffold-unconfirmed ') + side + ' descendant')
        output.append(dict(row, hgtfocus_event_count=len(matched),
                           hgtfocus_event_ids="; ".join(e["event_id"] for e, _, _ in matched),
                           hgtfocus_node_label="; ".join(f"{label} UFB=" + (f"{support:g}" if support is not None else "NA")
                                                         for _, label, support in matched),
                           hgtfocus_recipient_flag=int(row["node_name"] in genes),
                           hgtfocus_donor_flag=int(name in donor_genes),
                           hgtfocus_tip_status='; '.join(status)))
    return output, audits


def write(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def export_gene_trees(directory, events, links, family_root, renderer=None, gff_root='', context_annotations='',
                      mmseqs2_taxonomy_dir='', scaffold_taxonomy_dir='', taxonomy_dbfile='', minimum_ufboot=None):
    """Export one native PDF per family; unavailable mappings remain in the audit."""
    links = validated_event_links(events, links)
    csv.field_size_limit(100_000_000)
    directory.mkdir(parents=True)
    helper = Path(__file__).resolve().parent
    families = defaultdict(list)
    for event in events:
        family = event["orthogroup"]
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", family) or family in {".", ".."}:
            raise ValueError("Unsafe orthogroup identifier")
        families[family].append(event)
    from focus_hgt_context import CONTEXT_MAX_GENES_PER_SIDE, GenomeCoordinates, render_context
    from focus_hgt_context_annotations import ContextAnnotations
    from gene_tree_plot_config import replay

    coordinates = GenomeCoordinates(gff_root)
    annotations = ContextAnnotations(context_annotations, mmseqs2_taxonomy_dir=mmseqs2_taxonomy_dir,
                                     scaffold_taxonomy_dir=scaffold_taxonomy_dir, taxonomy_dbfile=taxonomy_dbfile)
    index, audit, sources, context_audit, configurations = [], [], {}, [], {}
    with read_only_observation():
        store = GeneFamilyOutputStore(family_root)
        annotations.store = store
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
                annotated, checks = annotate(rows, group, [link for link in links if link["orthogroup"] == family],
                                             minimum_ufboot=minimum_ufboot)
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
                                                            links, coordinates, gene_tree_panel=False,
                                                            max_genes_per_side=CONTEXT_MAX_GENES_PER_SIDE,
                                                            annotations=annotations)
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
            annotations.verify()
    write(directory / "index.tsv", INDEX_FIELDS, index)
    write(directory / "event_node_audit.tsv", EVENT_FIELDS, audit)
    if context_audit:
        write(directory / 'context_gene_audit.tsv', list(context_audit[0]), context_audit)
    if annotations.display_audit:
        write(directory / 'context_annotation_audit.tsv', list(annotations.display_audit[0]), annotations.display_audit)
    (directory / 'renderer_settings.json').write_text(json.dumps(configurations, indent=2) + '\n')
    coordinates.verify()
    (directory / "README.txt").write_text(
        "Native GeneGalleon gene trees for observed category-1 recipients\n\n"
        "Orange diamonds and HGT labels mark exact gene-tree transfer nodes, including internal nodes.\n"
        "UFB = Ultrafast bootstrap; labels retain the matched branch's measured support_generax_ufboot.\n"
        + ("No UFB threshold is applied; missing support is NA, never imputed.\n" if minimum_ufboot is None
           else f"Gene-tree UFB threshold: >= {minimum_ufboot} inclusive.\n") +
        "At least one retained event-linked gene on each side must have candidate-free class background\n"
        "with >=10 classified units, >=50% classification coverage and >=90% host compatibility.\n"
        "Orange recipient and blue donor tips individually pass that background check.\n"
        "The shared descendants column also records eligible genes with unconfirmed scaffold support.\n"
        "Page 1 replays gg_gene_evolution panels and saved settings, including domain, gene structure and alignment.\n"
        "Syntenic similarity and Sequence identity are disabled in focused replay; the synteny neighborhood remains.\n"
        "Missing optional measurements are not invented. Renderer settings and input availability are recorded.\n"
        "Page 2 separates donor descendants (blue, left) and recipient descendants (orange, right).\n"
        "It shows at most three distinct genes per side on one page, including eligible genes with unconfirmed or failing scaffold evidence.\n"
        "Shown/total/omitted gene counts and individual scaffold status remain explicit; omitted genes stay in context_gene_audit.tsv.\n"
        "Display priority is passing scaffold support, available GFF, coverage, compatibility, then gene ID.\n"
        "Repeated links for one side/gene are drawn once and preserve every event in the audit. No extra gene-tree inset is drawn.\n"
        "Tracks share an exon/UTR kb scale centered on the focal midpoint; noncoding gaps >5 kb are capped at 2 kb.\n"
        "Numbered // marks and titles identify intergenic/intronic omissions; recorded exon/UTR blocks remain uncompressed.\n"
        "CDS and UTR blocks are distinct; exon-only blocks have unknown CDS/UTR identity; missing structures stay unconfirmed.\n"
        "Each displayed focal/neighbor gene has its own product, best-hit organism/accession and kingdom-to-genus ranks.\n"
        "Protein products always use Swiss-Prot best-hit predictions; missing names/ranks stay unavailable. GFF products are not displayed.\n"
        "MMseqs2 query classification (LCA name/rank/taxid), kingdom-to-genus names and per-gene host labels are separate from Swiss-Prot.\n"
        "Query ranks reuse saved lineage taxids and an existing read-only database; lower unresolved ranks are not filled from the host or best hit.\n"
        "Column order: Track label, Protein product, Swiss-Prot best hit, MMseqs2 classification.\n"
        "Unresolved classification is not a host match. These display fields do not change candidate selection.\n"
        "Every coordinate-bearing model intersecting the display range is drawn; overlapping models have separate lanes.\n"
        "Only the focal, two nearest loci on each coordinate side and up to two overlaps are numbered/listed in the annotation table.\n"
        "Insufficient scaffold annotations are explicit; overlapping loci never substitute for missing flanks.\n"
        "Best-hit taxonomy does not identify the modeled donor or establish host background for a neighbor.\n"
        "context_annotation_audit.tsv retains the per-gene input fields, sources and exact event/context mapping.\n"
        "This is whole-scaffold context, not conserved gene order or proof of physical integration.\n"
        "All event IDs, branch/node IDs and selection/withholding reasons are in event_node_audit.tsv.\n"
        "No sequence or phylogenetic analysis is run. The parent focused event tables are unchanged.\n")
    return dict(profile=dict(PROFILE, minimum_ufboot=minimum_ufboot), family_source_sha256=sources, gff_source_sha256=coordinates.sources,
                context_annotation_source_sha256=annotations.sources,
                context_neighbor_family_source_sha256=annotations.family_sources,
                rendered_family_count=sum(row["status"] == "rendered" for row in index),
                selected_event_count=sum(row["status"] == "selected" for row in audit))

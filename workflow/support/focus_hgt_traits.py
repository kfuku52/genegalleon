#!/usr/bin/env python3
"""Return an existing HGT event cohort to observed category-1 recipient targets."""

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import tempfile
from collections import defaultdict
from pathlib import Path

import pandas
from hgt_species_tree import read_species_tree
from plot_hgt_summary import DEFAULT_TRANSFER_ARROW_ALPHA, validate_transfer_arrow_alpha
from species_trait_contract import select_foreground_traits
from species_trait_schema import schema_path, schema_payload, trait_value_types

VERSION = 1
MISSING = {"", ".", "na", "nan", "none", "null"}
EVENT_REQUIRED = {"event_id", "orthogroup", "generax_transfer", "generax_donor_node", "generax_recipient_node"}
LINK_REQUIRED = {"event_id", "orthogroup", "gene_id", "gene_species", "side", "eligible_for_context"}
FOCUS_FIELDS = ["focus_trait", "focus_category", "focus_recipient_basis", "focus_ancestral_state_inferred"]
PAIR_FIELDS = ["generax_donor_node", "generax_recipient_node", "event_count", "orthogroup_count",
               "event_ids", "orthogroups", "donor_clade_tip_labels", "recipient_clade_tip_labels"]


def digest(path):
    sha = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def key(value):
    return str(value).strip().replace(" ", "_")


def category_one(value):
    text = str(value).strip()
    if text.lower() in MISSING:
        return None
    try:
        numeric = float(text)
    except ValueError:
        return 0
    return int(math.isfinite(numeric) and numeric == 1)


def binary_value(value):
    text = str(value).strip()
    if text.lower() in MISSING:
        return None
    try:
        numeric = float(text)
    except ValueError:
        raise ValueError(f"Invalid binary trait value: {value!r}") from None
    if numeric not in (0, 1):
        raise ValueError(f"Invalid binary trait value: {value!r}")
    return int(numeric)


def focused_traits(path):
    kinds = trait_value_types(path)
    frame, audit = select_foreground_traits(path, allow_empty=True)
    identifiers = frame.iloc[:, 0].map(key)
    if identifiers.eq("").any() or identifiers.duplicated().any():
        raise ValueError("Species traits require unique nonempty species IDs, including aliases")
    selected, reports = {}, []
    for name in frame.columns[1:]:
        kind = kinds[name]
        values = frame[name].tolist()
        reason = ""
        if kind in {"numeric", "text"}:
            reason = "declared_" + kind
        elif kind == "unspecified":
            try:
                parsed = [binary_value(value) for value in values]
                if all(value is None for value in parsed):
                    reason = "no_observed_binary_values"
            except ValueError:
                reason = "not_an_implicit_binary_trait"
        elif kind == "binary":
            parsed = [binary_value(value) for value in values]
        else:
            parsed = [category_one(value) for value in values]
        if not reason:
            selected[name] = dict(zip(identifiers, parsed, strict=True))
        reports.append(dict(trait=name, value_type=kind, status="excluded" if reason else "selected", reason=reason))
    return selected, reports, audit


def read_tsv(path, required=()):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if not fields or len(fields) != len(set(fields)) or not set(required) <= set(fields):
            raise ValueError(f"Missing or duplicate required TSV columns: {path}")
        rows = []
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed TSV row: {path}")
            rows.append(row)
    return fields, rows


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({name: row.get(name, "") for name in fields} for row in rows)


def safe_names(names):
    result = {}
    used = set()
    for name in sorted(names):
        token = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._")[:100] or "target"
        if token in used:
            token += "_" + hashlib.sha256(name.encode()).hexdigest()[:12]
        if token in used:
            raise ValueError("Ambiguous output target names")
        used.add(token)
        result[name] = token
    return result


def add_clade_labels(events, fields, nodes):
    for side in ("donor", "recipient"):
        for suffix in ("clade_tip_labels", "clade_tip_count"):
            name = f"{side}_{suffix}"
            if name not in fields:
                fields.append(name)
        for event in events:
            tips = nodes.get(key(event[f"generax_{side}_node"]), ())
            if not tips:
                continue
            if f"{side}_clade_tip_labels" in event:
                existing = {key(tip) for tip in event[f"{side}_clade_tip_labels"].split(";") if tip.strip()}
                if existing != set(tips):
                    raise ValueError("Event clade tips disagree with the supplied analysis species tree")
                if f"{side}_clade_tip_count" in event and event[f"{side}_clade_tip_count"] != str(len(tips)):
                    raise ValueError("Event clade tip count disagrees with the analysis species tree")
            else:
                event[f"{side}_clade_tip_labels"] = "; ".join(tips)
                event[f"{side}_clade_tip_count"] = str(len(tips))


def export_bundle(directory, events, fields, links, link_fields, trait, trait_table, tree_path, plots, tip=None,
                  arrow_alpha=DEFAULT_TRANSFER_ARROW_ALPHA, render_tree=True):
    directory.mkdir(parents=True)
    write_tsv(directory / "events.tsv", fields, events)
    direct = [row for row in events if row["focus_recipient_basis"] == "observed_tip_category1"]
    ancestral = [row for row in events if row not in direct]
    write_tsv(directory / "direct_events.tsv", fields, direct)
    write_tsv(directory / "ancestral_recipient_events.tsv", fields, ancestral)
    event_ids = {row["event_id"] for row in events}
    selected_links = [row for row in links if row["event_id"] in event_ids]
    write_tsv(directory / "event_gene_links.tsv", link_fields, selected_links)
    counts = {}
    for side in ("donor", "recipient"):
        genes = defaultdict(list)
        for row in selected_links:
            if row["side"] == side and str(row["eligible_for_context"]).lower() in {"true", "1"}:
                if tip and side == "recipient" and key(row["gene_species"]) != tip:
                    continue
                genes[(row["orthogroup"], row["gene_id"])].append(row)
        rows = []
        event_identity = {"event_id", "branch_id", "node_name", "event_index", "generax_transfer",
                          "generax_donor_node", "generax_recipient_node"}
        gene_fields = [name for name in link_fields if name not in event_identity]
        for name in ("focus_event_ids", "focus_event_link_count"):
            if name in gene_fields:
                raise ValueError("Reserved focus gene column already exists")
        for _, group in sorted(genes.items()):
            row = {name: group[0][name] for name in gene_fields}
            row["focus_event_ids"] = "; ".join(sorted({g["event_id"] for g in group}))
            row["focus_event_link_count"] = str(len({g["event_id"] for g in group}))
            rows.append(row)
        write_tsv(directory / f"{side}_genes.tsv", gene_fields + ["focus_event_ids", "focus_event_link_count"], rows)
        counts[f"{side}_gene_count"] = len({gene_id for _, gene_id in genes})
    pairs = defaultdict(list)
    for row in events:
        pairs[(row["generax_donor_node"], row["generax_recipient_node"])].append(row)
    pair_rows = []
    for (donor, recipient), group in sorted(pairs.items()):
        pair_rows.append(dict(generax_donor_node=donor, generax_recipient_node=recipient,
                              event_count=len(group), orthogroup_count=len({r["orthogroup"] for r in group}),
                              event_ids="; ".join(sorted(r["event_id"] for r in group)),
                              orthogroups="; ".join(sorted({r["orthogroup"] for r in group})),
                              donor_clade_tip_labels=group[0]["donor_clade_tip_labels"],
                              recipient_clade_tip_labels=group[0]["recipient_clade_tip_labels"]))
    write_tsv(directory / "branch_pairs.tsv", PAIR_FIELDS, pair_rows)
    summary = dict(event_count=len(events), direct_event_count=len(direct), ancestral_recipient_event_count=len(ancestral),
                   orthogroup_count=len({row["orthogroup"] for row in events}), branch_pair_count=len(pairs), **counts)
    (directory / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if plots:
        from plot_hgt_summary import plot_transfer_tree

        plot_transfer_tree(pandas.DataFrame(events, columns=fields),
                           str(directory / "transfer_tree.pdf") if render_tree else "",
                           species_tree_path=str(tree_path), edges_tsv=str(directory / "transfer_edges.tsv"),
                           max_edges=0, species_trait_path=str(trait_table), highlight_trait=trait,
                           arrow_alpha=arrow_alpha)
    return summary


def build_focus(stage, event_path, link_path, tree_path, trait_path, plots=True,
                arrow_alpha=DEFAULT_TRANSFER_ARROW_ALPHA, gene_family_root="", gff_root='', filter_audit='',
                context_annotations=''):
    csv.field_size_limit(100_000_000)
    fields, events = read_tsv(event_path, EVENT_REQUIRED)
    link_fields, links = read_tsv(link_path, LINK_REQUIRED)
    ids = [event["event_id"] for event in events]
    if any(not value for value in ids) or len(ids) != len(set(ids)):
        raise ValueError("Duplicate or empty transfer event_id")
    events_by_id = {event["event_id"]: event for event in events}
    links = [row for row in links if row["event_id"] in events_by_id]
    identities = set()
    for row in links:
        identity = (row["event_id"], row["side"], row["gene_id"])
        if row["side"] not in {"donor", "recipient"} or identity in identities:
            raise ValueError("Duplicate or invalid event-gene link")
        identities.add(identity)
        event = events_by_id[row["event_id"]]
        for name in EVENT_REQUIRED & set(link_fields):
            if row[name] != event[name]:
                raise ValueError(f"Event-gene identity mismatch: {name}")
    tree = read_species_tree(tree_path)
    nodes = {key(node.name): tuple(sorted(key(tip.name) for tip in node.get_terminals()))
             for node in tree.find_clades() if node.name is not None}
    terminals = {key(tip.name) for tip in tree.get_terminals()}
    for row in links:
        if str(row.get('eligible_for_context', '')).lower() in {'true', '1'}:
            event = events_by_id[row['event_id']]
            branch = key(event[f"generax_{row['side']}_node"])
            species = key(row.get('gene_species', ''))
            if branch in nodes and species not in nodes[branch]:
                raise ValueError('Eligible context gene is outside its species branch: ' + row['gene_id'])
            if row.get('lineage_status', 'retained') != 'retained':
                raise ValueError('Eligible context gene does not have a retained transfer lineage')
    add_clade_labels(events, fields, nodes)
    if set(FOCUS_FIELDS) & set(fields):
        raise ValueError("Reserved focus columns already exist in event input")
    fields += FOCUS_FIELDS
    traits, reports, audit = focused_traits(trait_path)
    trait_dirs = safe_names(traits)
    write_tsv(stage / "trait_selection.tsv", ["trait", "value_type", "status", "reason"], reports)
    index, gene_trees, figures = [], {}, {}
    for trait, values in traits.items():
        positive = {species for species, value in values.items() if value == 1}
        if positive - terminals:
            raise ValueError(f"Category-1 species absent from analysis species tree: {sorted(positive - terminals)}")
        selected_nodes = {name for name, tips in nodes.items() if tips and set(tips) <= positive}
        root = stage / "traits" / trait_dirs[trait]
        root.mkdir(parents=True)
        indicator = root / "species_trait_category1.tsv"
        write_tsv(indicator, ["species", trait], [dict(species=species, **{trait: "" if value is None else value})
                                                  for species, value in values.items()])
        schema_path(indicator).write_bytes(schema_payload(indicator.read_bytes(), {trait: "binary"}))
        branch_rows = [dict(species_branch=name, branch_type="terminal" if name in terminals else "internal",
                            clade_tip_labels="; ".join(tips), descendant_tip_count=len(tips),
                            category1_tip_count=sum(values.get(tip) == 1 for tip in tips),
                            unknown_tip_count=sum(values.get(tip) is None for tip in tips),
                            selected=int(name in selected_nodes), ancestral_state_inferred=0)
                       for name, tips in nodes.items()]
        write_tsv(root / "species_branches.tsv", list(branch_rows[0]), branch_rows)
        selected, withheld = [], []
        for event in events:
            donor, recipient = key(event["generax_donor_node"]), key(event["generax_recipient_node"])
            if donor not in nodes or recipient not in nodes:
                reason = "species_branch_unmapped"
            elif event["generax_transfer"] != f"Y@{event['generax_donor_node']}@{event['generax_recipient_node']}":
                raise ValueError("Event transfer token disagrees with donor/recipient branch IDs")
            elif event.get("mapping_status", "matched") != "matched":
                reason = "event_mapping_unresolved"
            elif recipient not in selected_nodes:
                reason = "recipient_clade_not_all_observed_category1"
            else:
                row = dict(event, focus_trait=trait, focus_category="1", focus_ancestral_state_inferred="0",
                           focus_recipient_basis="observed_tip_category1" if recipient in terminals
                           else "all_descendant_tips_category1_no_ancestral_reconstruction")
                selected.append(row)
                continue
            withheld.append(dict(event_id=event["event_id"], reason=reason))
        write_tsv(root / "events_not_focused.tsv", ["event_id", "reason"], withheld)
        aggregate = export_bundle(root / "all_category1", selected, fields, links, link_fields,
                                  trait, indicator, tree_path, plots, arrow_alpha=arrow_alpha, render_tree=False)
        if plots and gene_family_root:
            from focus_hgt_gene_trees import export_gene_trees

            gene_trees[trait] = export_gene_trees(root / "all_category1/tree_plot", selected, links, gene_family_root,
                                                 gff_root=gff_root, context_annotations=context_annotations)
            from focus_hgt_figures import export_figures
            checks = read_tsv(root / 'all_category1/tree_plot/event_node_audit.tsv')[1]
            ids = {r['event_id'] for r in checks if r['status'] == 'selected'}
            figure_events = [r for r in selected if r['event_id'] in ids]
            figures[trait] = export_figures(root / 'all_category1/plots', events, figure_events,
                                            [r for r in links if r['event_id'] in ids], tree, values,
                                            gene_family_root, trait, filter_audit=filter_audit,
                                            context_annotations=context_annotations)
        index.append(dict(trait=trait, target="ALL_CATEGORY1", target_type="aggregate",
                          relative_path=str((root / "all_category1").relative_to(stage)), **aggregate))
        target_dirs = safe_names(selected_nodes)
        for target in sorted(selected_nodes):
            terminal = target in terminals
            selected_events = [row for row in selected if key(row["generax_recipient_node"]) == target
                               or (terminal and target in nodes[key(row["generax_recipient_node"])])]
            destination = root / ("tips" if terminal else "internal_branches") / target_dirs[target]
            summary = export_bundle(destination, selected_events, fields, links, link_fields,
                                    trait, indicator, tree_path, plots, tip=target if terminal else None,
                                    arrow_alpha=arrow_alpha, render_tree=False)
            index.append(dict(trait=trait, target=target, target_type="tip" if terminal else "internal_branch",
                              relative_path=str(destination.relative_to(stage)), **summary))
    index_fields = ["trait", "target", "target_type", "relative_path", "event_count", "direct_event_count",
                    "ancestral_recipient_event_count", "orthogroup_count", "branch_pair_count", "donor_gene_count",
                    "recipient_gene_count"]
    write_tsv(stage / "index.tsv", index_fields, index)
    (stage / "README.txt").write_text(
        "HGT candidates focused on observed trait category 1\n\n"
        "The input event cohort and its existing support, direction, scaffold, quality and annotation fields are preserved.\n"
        "No additional HGT threshold is applied. UFBoot is gene-tree branch support, not an HGT probability.\n"
        "Binary traits, categorical category 1, and schema-free observed 0/1 traits are selected automatically.\n"
        "Declared numeric/text traits and observation/quality columns are excluded. Missing traits stay unknown.\n"
        "Internal recipient branches qualify only when every descendant tip has observed category 1. This is not ancestral reconstruction.\n"
        "Per-tip direct_events.tsv contains transfers into that tip. ancestral_recipient_events.tsv contains qualifying ancestral-branch context.\n"
        "The same ancestral event can appear in multiple tip reports; aggregate event IDs count it once. Do not sum per-tip totals.\n"
        "Recipient gene tables in tip reports contain that tip's eligible event-linked genes. Donor gene tables preserve event-linked donor homologs.\n"
        "An empty table means no selected result in this input cohort, not biological absence of HGT.\n"
        "Aggregate figures are three single-page PDFs: filtering flow, annotated family/species distribution and donor/recipient counts.\n"
        "Individual qualifying orthogroup PDFs have gg_gene_evolution panels on page 1 and genomic context on page 2.\n"
        "Per-recipient tips and internal branches retain tables and directed edge TSVs but have no separate tree PDF.\n"
        "When gene-family inputs are supplied, each trait aggregate's tree_plot/ contains native per-orthogroup gene-tree PDFs.\n"
        "Only exact transfer branches with UFBoot >=90 and individually passing scaffold background genes on both sides are marked.\n"
        "See tree_plot/event_node_audit.tsv for every selected/withheld event and tree_plot/README.txt for the evidence profile.\n"
        "Scaffold background and shared-neighbor synteny are distinct evidence. Neither proves physical integration.\n")
    return dict(schema_version=VERSION, source_event_count=len(events), trait_contract=audit,
                trait_selection=reports, result_index=index, plots=plots, transfer_arrow_alpha=arrow_alpha,
                plot_scope="trait_aggregate_only", gene_tree_plots=gene_trees, summary_figures=figures)


def generate(event_path, link_path, tree_path, trait_path, output, plots=True,
             arrow_alpha=DEFAULT_TRANSFER_ARROW_ALPHA, gene_family_root="", gff_root='', filter_audit='',
             context_annotations=''):
    arrow_alpha = validate_transfer_arrow_alpha(arrow_alpha)
    inputs = [Path(path).resolve() for path in (event_path, link_path, tree_path, trait_path)]
    if filter_audit and plots:
        inputs.append(Path(filter_audit).resolve())
    if context_annotations and plots and gene_family_root:
        inputs.append(Path(context_annotations).resolve())
    for suffix in (".schema.json", ".metadata.json"):
        sidecar = Path(str(trait_path) + suffix)
        if sidecar.exists():
            inputs.append(sidecar.resolve())
    helper_root = Path(__file__).resolve().parent
    code = [helper_root / name for name in ("focus_hgt_traits.py", "plot_hgt_summary.py", "hgt_species_tree.py",
                                            "species_trait_contract.py", "species_trait_schema.py")]
    if gene_family_root and plots:
        code += [helper_root / name for name in ('focus_hgt_gene_trees.py', 'stat_branch2tree_plot.r',
                                                'focus_hgt_context.py', 'focus_hgt_context_annotations.py',
                                                'focus_hgt_figures.py', 'gene_tree_plot_config.py')]
        code += sorted((helper_root / "treevis/R").glob("*.R"))
    output = Path(output).absolute()
    # These directories are read inputs too. Replacing a managed report must
    # never remove curated family/GFF sources nested underneath it, even via aliases.
    protected = inputs + [Path(path).resolve() for path in (gene_family_root, gff_root) if path]
    if output.is_symlink() or any(path == output.resolve() or output.resolve() in path.parents for path in protected):
        raise ValueError("Output must not replace or contain an input")
    if output.exists():
        manifest = output / "manifest.json"
        if not manifest.is_file() or json.loads(manifest.read_text()).get("schema_version") != VERSION:
            raise ValueError("Refusing to replace an unmanaged focused-results directory")
    before = {str(path): digest(path) for path in inputs}
    code_before = {path.name: digest(path) for path in code}
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".hgt-trait-focus-", dir=output.parent))
    backup = None
    try:
        manifest = build_focus(stage, *inputs[:4], plots=plots, arrow_alpha=arrow_alpha,
                               gene_family_root=gene_family_root, gff_root=gff_root, filter_audit=filter_audit,
                               context_annotations=context_annotations)
        if any(digest(path) != before[str(path)] for path in inputs):
            raise ValueError("Focused-analysis inputs changed during generation")
        if any(digest(path) != code_before[path.name] for path in code):
            raise ValueError("Focused-analysis code changed during generation")
        manifest["inputs_sha256"] = before
        manifest["code_sha256"] = code_before
        manifest["outputs_sha256"] = {str(path.relative_to(stage)): digest(path)
                                      for path in sorted(stage.rglob("*")) if path.is_file()}
        (stage / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        if output.exists():
            backup = Path(tempfile.mkdtemp(prefix=".hgt-trait-focus-backup-", dir=output.parent))
            backup.rmdir()
            os.replace(output, backup)
        try:
            os.replace(stage, output)
        except BaseException:
            if backup is not None:
                os.replace(backup, output)
                backup = None
            raise
        if backup is not None:
            shutil.rmtree(backup)
        return manifest
    finally:
        if stage.exists():
            shutil.rmtree(stage)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event_tsv", required=True)
    parser.add_argument("--event_gene_tsv", required=True)
    parser.add_argument("--species_tree", required=True)
    parser.add_argument("--species_trait", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--plots", choices=("0", "1"), default="1")
    parser.add_argument("--gene_family_root", default="",
                        help="Existing raw/ZIP-backed family outputs for scaffold-supported category-1 gene-tree PDFs")
    parser.add_argument('--gff_info_root', default='', help='Existing per-species gff_info TSVs for context-page structures')
    parser.add_argument('--context_annotations_tsv', default='',
                        help='Existing per-gene product and best-hit taxonomy annotations for context pages')
    parser.add_argument('--filter_audit_tsv', default='', help='Optional project event-level direction/support filtering audit (TSV or TSV.gz)')
    parser.add_argument("--transfer_arrow_alpha", type=validate_transfer_arrow_alpha,
                        default=DEFAULT_TRANSFER_ARROW_ALPHA)
    args = parser.parse_args()
    manifest = generate(args.event_tsv, args.event_gene_tsv, args.species_tree, args.species_trait,
                        args.output_dir, plots=args.plots == "1", arrow_alpha=args.transfer_arrow_alpha,
                        gene_family_root=args.gene_family_root, gff_root=args.gff_info_root, filter_audit=args.filter_audit_tsv,
                        context_annotations=args.context_annotations_tsv)
    print(json.dumps(dict(output_dir=str(Path(args.output_dir).resolve()), source_event_count=manifest["source_event_count"],
                          result_sets=len(manifest["result_index"])), indent=2))


if __name__ == "__main__":
    main()

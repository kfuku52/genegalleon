#!/usr/bin/env python3
"""Resolve input-species taxonomy and map ranks onto a species or NCBI tree."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import math
import os
import re
import shutil
import sqlite3
import sys
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import quote

from format_species_taxonomy import build_species_key_from_tokens, tokenize_taxonomic_name

DEFAULT_RANKS = "all"
FILES = (
    "species_taxonomy.tsv", "species_lineage.tsv", "taxonomy_mapping.tsv",
    "taxonomy_tree.nwk", "taxonomy_tree.nhx", "taxonomy_tree.pdf",
    "taxonomy_tree.svg", "taxonomy_tree.png", "taxonomy_columns.tsv", "provenance.json",
)
FASTA_RE = re.compile(r"\.(?:fa|fas|fasta|fna|faa)(?:\.gz)?$", re.I)


def rank_column(rank):
    return "species_rank" if rank == "species" else rank


def rank_property(rank):
    return "gg_" + rank_column(rank).replace(" ", "_")


def pack_values(values):
    """Keep singleton cells simple; repeated ranks retain every ancestor in order."""
    values = [str(value) for value in values]
    return json.dumps(values, ensure_ascii=False) if len(values) > 1 else (values[0] if values else "")


def display_value(value):
    if value.startswith("["):
        try:
            values = json.loads(value)
        except ValueError:
            return value  # Scientific names can themselves start with brackets.
        if isinstance(values, list) and all(isinstance(item, str) for item in values):
            return "\n".join(values)
    return value


def available_ranks(lineages):
    depths, predecessors = {}, defaultdict(set)
    for lineage in lineages.values():
        seen, previous = set(), None
        for depth, node in enumerate(lineage):
            rank = node["rank"]
            if rank:
                depths[rank] = min(depth, depths.get(rank, depth))
                if rank not in seen:
                    if previous is not None:
                        predecessors[rank].add(previous)
                    seen.add(rank)
                    previous = rank
    # Shared ancestry orders sparse columns (e.g. subclass before order even
    # when other inputs lack subclass). Repeated ranks use their first position.
    pending, ordered = set(depths), []
    while pending:
        ready = [rank for rank in pending if not predecessors[rank] & pending]
        # Inconsistent rank ordering across lineages has no unique solution;
        # retain every column with the same deterministic depth/name tie-break.
        rank = min(ready or pending, key=lambda value: (depths[value], value))
        ordered.append(rank)
        pending.remove(rank)
    return ordered


def aligned_columns(rows, lineages, ranks, tree):
    """Align clades by TaxID, interleaving their columns with named ranks."""
    columns = {
        rank_column(rank): dict(column=rank_column(rank), taxid_column=rank + "_taxid",
                                rank=rank, taxid="", label=rank.capitalize())
        for rank in ranks if rank != "clade"
    }
    depths, predecessors = {}, defaultdict(set)
    for lineage in lineages.values():
        previous, seen = None, set()
        for depth, node in enumerate(lineage):
            rank = node["rank"]
            if rank not in ranks:
                continue
            key = rank_column(rank)
            if rank == "clade":
                key = f"clade_{node['taxid']}"
                descriptor = dict(column=key, taxid_column=key + "_taxid", rank=rank,
                                  taxid=str(node["taxid"]), label=node["spname"])
                if key in columns and columns[key] != descriptor:
                    raise ValueError(f"Conflicting taxonomy column identity: {key}")
                columns[key] = descriptor
            depths[key] = min(depth, depths.get(key, depth))
            if key not in seen:
                if previous is not None:
                    predecessors[key].add(previous)
                seen.add(key)
                previous = key
    # Equal-depth siblings are independent columns, never a shared rank slot.
    pending, ordered = set(columns), []
    while pending:
        ready = [key for key in pending if not predecessors[key] & pending]
        key = min(ready or pending, key=lambda value: (depths.get(value, math.inf), columns[value]["label"], value))
        descriptor = dict(columns[key], position=len(ordered) + 1)
        descriptor["unmet_predecessors"] = json.dumps(sorted(predecessors[key] & pending))
        ordered.append(descriptor)
        pending.remove(key)
    clade_columns = [column for column in ordered if column["rank"] == "clade"]
    tips = {} if tree is None else {leaf.name: leaf for leaf in tree.leaves()}
    for row in rows:
        members = {str(node["taxid"]) for node in lineages[row["species"]] if node["rank"] == "clade"}
        for column in clade_columns:
            member = column["taxid"] in members
            row[column["column"]] = column["label"] if member else ""
            row[column["taxid_column"]] = column["taxid"] if member else ""
            if member and row["tree_tip"] in tips:
                tips[row["tree_tip"]].add_prop("gg_" + column["column"], quote(column["label"], safe=" _-."))
    return ordered


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_table(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f"Missing or duplicate TSV columns: {path}")
        rows = []
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed TSV row: {path}")
            rows.append({key: value.strip() for key, value in row.items()})
        return reader.fieldnames, rows


def species_rows(path, *, merge_summary=False):
    fields, rows = read_table(path)
    key = next((key for key in ("species", "species_key", "leaf_name") if key in fields), None)
    if key is None:
        raise ValueError(f"TSV needs species, species_key, or leaf_name: {path}")
    result = {}
    for row in rows:
        label = row[key]
        if not label or (label in result and not merge_summary):
            raise ValueError(f"Empty or duplicate species identifier in {path}: {label!r}")
        if label in result:
            previous = result[label]
            candidates = set(previous.get("taxid_candidates", [previous["taxid"]])) | {row.get("taxid", "")}
            candidates.discard("")
            previous["taxid"] = next(iter(candidates)) if len(candidates) == 1 else ""
            if len(candidates) > 1:
                previous["taxid_candidates"] = sorted(candidates)
            continue
        result[label] = {"species": label, "taxid": row.get("taxid", "")}
    return result


def collect_species(args):
    """An explicit table is authoritative; otherwise use current FASTA inputs."""
    sources = []
    fasta_inputs = []
    if args.species_table:
        sources.append(Path(args.species_table))
        rows = species_rows(args.species_table)
    else:
        summary = Path(args.species_summary) if args.species_summary else (
            args.workspace / "output/input_generation/gg_input_generation_species.tsv")
        metadata = {}
        if args.species_summary or summary.is_file():
            metadata = species_rows(summary, merge_summary=True)
            sources.append(summary)
        directories = [Path(p) for p in args.species_dir] or [
            args.workspace / "input" / "species_cds",
            args.workspace / "input" / "species_protein",
        ]
        rows = {}
        for directory in directories:
            if not directory.exists():
                continue
            for path in sorted(directory.iterdir()):
                if path.name.startswith(".") or not path.is_file() or not FASTA_RE.search(path.name):
                    continue
                fasta_inputs.append(path)
                stem = FASTA_RE.sub("", path.name)
                matches = [key for key in metadata if stem == key or stem.startswith(key + "_")]
                label = max(matches, key=len) if matches else build_species_key_from_tokens(tokenize_taxonomic_name(stem))
                if label:
                    rows[label] = dict(metadata.get(label, {"species": label, "taxid": ""}))
    if not rows:
        raise ValueError("No input species found; provide --species-table or --species-dir.")
    if args.taxid_map:
        sources.append(Path(args.taxid_map))
        overrides = species_rows(args.taxid_map)
        unknown = set(overrides) - set(rows)
        if unknown:
            raise ValueError(f"TaxID overrides contain unknown input species: {sorted(unknown)}")
        for label, row in overrides.items():
            if not re.fullmatch(r"[1-9][0-9]*", row["taxid"]):
                raise ValueError(f"Invalid override TaxID for {label}: {row['taxid']!r}")
            rows[label].update(taxid=row["taxid"], taxid_source="override")
            rows[label].pop("taxid_candidates", None)
    inline = getattr(args, "taxid_override", "")
    if inline:
        match = re.fullmatch(r"([A-Za-z][A-Za-z0-9_.-]*):([1-9][0-9]*)", inline)
        if match is None:
            raise ValueError("TaxID override must be species:positive_taxid")
        label, taxid = match.groups()
        if label not in rows:
            raise ValueError(f"TaxID override contains unknown input species: {label}")
        if args.taxid_map and label in overrides and overrides[label]["taxid"] != taxid:
            raise ValueError(f"Conflicting TaxID overrides for {label}")
        rows[label].update(taxid=taxid, taxid_source="inline_override")
        rows[label].pop("taxid_candidates", None)
    return [rows[label] for label in sorted(rows)], sources, fasta_inputs


class Taxonomy:
    """Read the shared ETE database without refreshing it or choosing ambiguous hits."""

    def __init__(self, path):
        self.conn = sqlite3.connect(Path(path).resolve().as_uri() + "?mode=ro", uri=True)
        self.conn.row_factory = sqlite3.Row
        self.tables = {row[0] for row in self.conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        self.cache = {}

    def node(self, taxid):
        taxid = int(taxid)
        if taxid not in self.cache:
            self.cache[taxid] = self.conn.execute(
                "SELECT taxid, spname, rank, track FROM species WHERE taxid=?", (taxid,)).fetchone()
        return self.cache[taxid]

    def resolve(self, source):
        row = dict(source)
        row.update(input_taxid=source["taxid"], scientific_name="", resolution_status="unresolved", resolution_source="",
                   query=source["species"].replace("_", " "))
        if source.get("taxid_candidates"):
            row.update(input_taxid=";".join(source["taxid_candidates"]), resolution_status="conflicting_taxids",
                       resolution_source="input_summary")
            return row, []
        taxid = source["taxid"]
        if taxid:
            if not re.fullmatch(r"[1-9][0-9]*", taxid):
                raise ValueError(f"Invalid TaxID for {source['species']}: {taxid!r}")
            row["resolution_source"] = source.get("taxid_source", "input_taxid")
            node = self.node(taxid)
            if node is None and "merged" in self.tables:
                merged = self.conn.execute("SELECT taxid_new FROM merged WHERE taxid_old=?", (int(taxid),)).fetchone()
                if merged:
                    node = self.node(merged[0])
                    row["resolution_source"] += ":merged"
        else:
            # Query the entire identifier. Never degrade qualified/infraspecific
            # labels to a binomial or genus to make an unresolved name succeed.
            hits = self.conn.execute(
                "SELECT taxid FROM species WHERE spname=? COLLATE NOCASE", (row["query"],)).fetchall()
            source_name = "scientific_name"
            if not hits and "synonym" in self.tables:
                hits = self.conn.execute(
                    "SELECT taxid FROM synonym WHERE spname=? COLLATE NOCASE", (row["query"],)).fetchall()
                source_name = "synonym"
            ids = {hit[0] for hit in hits}
            if len(ids) > 1:
                row.update(taxid="", resolution_status="ambiguous", resolution_source=source_name)
                return row, []
            node = self.node(next(iter(ids))) if ids else None
            row["resolution_source"] = source_name
        if node is None:
            row.update(taxid="", resolution_status="invalid_taxid" if taxid else "unresolved")
            return row, []
        row.update(taxid=str(node["taxid"]), scientific_name=node["spname"], resolution_status="resolved")
        # ETE tracks run from the terminal taxon towards the root.
        ids = [int(value) for value in node["track"].split(",") if value]
        if not ids or ids[0] != node["taxid"] or len(ids) != len(set(ids)):
            raise ValueError(f"Invalid ETE lineage for TaxID {node['taxid']}")
        lineage = [self.node(value) for value in reversed(ids)]
        if any(value is None for value in lineage):
            raise ValueError(f"Incomplete taxonomy database lineage for TaxID {node['taxid']}")
        return row, [dict(value) for value in lineage]


def choose_tree(args):
    if args.species_tree != "auto":
        path = Path(args.species_tree)
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"Configured species tree is missing or empty: {path}")
        return path, "species_tree"
    output = args.workspace / "output"
    directories = [output / "species_tree/species_tree_summary", output / "species_tree",
                   output / "query2family/parameters", output / "orthogroup/parameters"]
    for directory in directories:
        for name in ("dated_species_tree.nwk", "undated_species_tree.nwk",
                     "dated_species_tree.pruned.nwk", "undated_species_tree.pruned.nwk"):
            path = directory / name
            if path.exists() or path.is_symlink():
                if not path.is_file() or not path.stat().st_size:
                    raise ValueError(f"Discovered species tree is empty or invalid: {path}")
                return path, "species_tree"
    return None, "ncbi_taxonomy"


def load_or_build_tree(path, lineages):
    from nwkit.util import read_tree, validate_unique_named_leaves

    if path:
        tree = read_tree(str(path), "auto", True, quiet=True)
    else:
        from nwkit.constrain import get_taxid_counts, taxid2tree
        from nwkit.util import remove_singleton

        resolved = {label: [node["taxid"] for node in lineage]
                    for label, lineage in lineages.items() if lineage}
        if not resolved:
            return None
        # The native constrain backend receives exactly the lineage snapshot
        # used for the tables, avoiding a second independently refreshed cache.
        tree = taxid2tree(resolved, get_taxid_counts(resolved))
        tree = remove_singleton(tree, verbose=False, preserve_branch_length=False)
        for node in tree.traverse():
            node.dist = 0 if node.is_root else 1
            node.props.pop("ancestors", None)
        from nwkit.rooting_state import set_rooting_info

        set_rooting_info(tree, True, source="ncbi_taxonomy")
    validate_unique_named_leaves(tree, option_name="species tree")
    return tree


def map_taxonomy(tree, rows, lineages, ranks):
    """Clade membership is evaluated on all tree tips, including unclassified tips."""
    by_species = {row["species"]: row for row in rows}
    aliases = defaultdict(list)
    for label in by_species:
        aliases[label.replace(" ", "_")].append(label)
    tips, tip_keys = {}, {}
    for leaf in ([] if tree is None else tree.leaves()):
        matches = [leaf.name] if leaf.name in by_species else aliases.get(leaf.name.replace(" ", "_"), [])
        if len(matches) > 1:
            raise ValueError(f"Ambiguous species-tree label: {leaf.name}")
        key = matches[0] if matches else leaf.name
        if key in tips:
            raise ValueError(f"Multiple tree tips map to the same input species: {key}")
        tips[key], tip_keys[leaf] = leaf, key
    if tree is not None and not (tips.keys() & by_species.keys()):
        raise ValueError("Species tree has no tips matching the current input species.")
    for row in rows:
        row["tree_status"] = "mapped" if row["species"] in tips else "missing_from_tree"
        row["tree_tip"] = tips[row["species"]].name if row["species"] in tips else ""
        by_rank = defaultdict(list)
        for node in lineages[row["species"]]:
            by_rank[node["rank"]].append(node)
        for rank in ranks:
            row[rank_column(rank)] = pack_values(node["spname"] for node in by_rank[rank])
            row[rank + "_taxid"] = pack_values(node["taxid"] for node in by_rank[rank])
    if tree is None:
        return [], {}
    descendants, node_ids, id_counts = {}, {}, defaultdict(int)
    for node in tree.traverse("postorder"):
        members = {tip_keys[node]} if node.is_leaf else set().union(*(descendants[ch] for ch in node.children))
        descendants[node] = members
        clade_id = "clade_" + hashlib.sha256(json.dumps(sorted(members)).encode()).hexdigest()[:16]
        id_counts[clade_id] += 1
        node_ids[node] = clade_id + (f"_u{id_counts[clade_id]}" if id_counts[clade_id] > 1 else "")
        node.add_prop("gg_clade_id", node_ids[node])
    for label, leaf in tips.items():
        row = by_species.get(label)
        if row is None:
            leaf.add_prop("gg_taxonomy_status", "not_input_species")
            continue
        leaf.add_prop("gg_taxonomy_status", row["resolution_status"])
        leaf.add_prop("gg_species_key", quote(row["species"], safe=" _-."))
        for key in ("taxid", *(rank_column(rank) for rank in ranks), *(rank + "_taxid" for rank in ranks)):
            if row[key]:
                leaf.add_prop("gg_" + key.replace(" ", "_"), quote(row[key], safe=" _-."))
    groups = defaultdict(set)
    names, depths = {}, {}
    for row in rows:
        for depth, taxon in enumerate(lineages[row["species"]]):
            rank = taxon["rank"]
            if rank in ranks:
                key = (rank, str(taxon["taxid"]))
                groups[key].add(row["species"])
                names[key] = taxon["spname"]
                depths[key] = min(depth, depths.get(key, depth))
    mappings, labels = [], defaultdict(list)
    annotations = defaultdict(lambda: defaultdict(list))
    for (rank, taxid), members in sorted(groups.items(), key=lambda item: (ranks.index(item[0][0]), depths[item[0]], item[0][1])):
        present = members & tips.keys()
        missing = members - tips.keys()
        node = None
        extra = set()
        if present:
            node = tips[next(iter(present))]
            while not present <= descendants[node]:
                node = node.up
            extra = descendants[node] - members
        unclassified = {label for label in extra if not by_species.get(label, {}).get(rank + "_taxid")}
        status = "missing_from_tree" if node is None else (
            ("non_monophyletic" if extra - unclassified else "unresolved_membership") if extra else (
                "singleton" if len(present) == 1 else "monophyletic"))
        mappings.append(dict(rank=rank, taxid=taxid, name=names[(rank, taxid)], status=status,
                             node_id=node_ids.get(node, ""), mapped_count=len(present),
                             missing_count=len(missing), members=json.dumps(sorted(present)),
                             missing_members=json.dumps(sorted(missing)), other_descendants=json.dumps(sorted(extra)),
                             unclassified_descendants=json.dumps(sorted(unclassified))))
        if node is not None and not extra:
            # Tips already contain the complete lineage, including repeated ranks.
            if not node.is_leaf:
                annotations[node][rank].append((names[(rank, taxid)], taxid))
                labels[node].append((rank, names[(rank, taxid)]))
    for node, by_rank in annotations.items():
        for rank, values in by_rank.items():
            node.add_prop(rank_property(rank), quote(pack_values(name for name, _ in values), safe=" _-."))
            node.add_prop("gg_" + rank.replace(" ", "_") + "_taxid", quote(pack_values(taxid for _, taxid in values), safe=" _-."))
    return mappings, labels


def write_tsv(path, rows, fields):
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def plot(tree, rows, columns, labels, source, directory, plot_clades=False):
    os.environ.setdefault("MPLCONFIGDIR", str(directory / ".matplotlib"))
    logging.getLogger("fontTools.subset").setLevel(logging.WARNING)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgb
    from matplotlib.patches import Rectangle

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "svg.fonttype": "none",
                         "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
    by_species = {row["species"]: row for row in rows}
    by_species.update({row["tree_tip"]: row for row in rows if row["tree_tip"]})
    leaves = [] if tree is None else list(tree.leaves())
    absent = [row["species"] for row in rows if row["tree_status"] == "missing_from_tree"]
    display = [leaf.name for leaf in leaves] + absent
    display_columns = [column for column in columns
                       if (plot_clades or column["rank"] != "clade")
                       and any(row[column["column"]] for row in rows)]
    # Label space grows with actual names, instead of clipping long species/ranks.
    name_width = max(2.6, max(map(len, display), default=1) * 0.075)
    band_widths = [max(1.1, len(column["label"]) * 0.075,
                      max((len(line) for row in rows for line in display_value(row[column["column"]]).split("\n")), default=1) * 0.065)
                   for column in display_columns]
    widths = [4.2, name_width, *band_widths] if display_columns else [4.2, name_width]
    max_lines = max((display_value(row[column["column"]]).count("\n") + 1 for row in rows for column in display_columns), default=1)
    height = max(4.8, len(display) * max(0.38, max_lines * 0.19) + 2.65)
    plot_height = height - 2.45
    fig = plt.figure(figsize=(sum(widths) + 0.8, height), facecolor="white")
    grid = fig.add_gridspec(1, len(widths), width_ratios=widths, wspace=0.04,
                            left=0.035, right=0.985, top=1 - 1.4 / height, bottom=1.05 / height)
    axes = [fig.add_subplot(grid[0, index]) for index in range(len(widths))]
    y = {label: len(display) - index - 1 for index, label in enumerate(display)}
    for ax in axes:
        ax.set_ylim(-0.7, max(len(display) - 0.3, 0.7))
        ax.set_axis_off()
    tree_ax, text_ax = axes[:2]
    if tree is not None:
        xpos, ypos = {tree: 0.0}, {}
        distances = [node.dist for node in tree.traverse() if not node.is_root]
        use_lengths = source == "species_tree" and all(
            value is not None and math.isfinite(value) and value >= 0 for value in distances) and any(distances)
        for node in tree.traverse("preorder"):
            if not node.is_root:
                xpos[node] = xpos[node.up] + (node.dist if use_lengths else 1.0)
        for node in tree.traverse("postorder"):
            ypos[node] = y[node.name] if node.is_leaf else sum(ypos[ch] for ch in node.children) / len(node.children)
            if not node.is_leaf:
                tree_ax.plot([xpos[node]] * 2, [min(ypos[ch] for ch in node.children), max(ypos[ch] for ch in node.children)],
                             color="#334155", lw=1.1)
            if not node.is_root:
                tree_ax.plot([xpos[node.up], xpos[node]], [ypos[node]] * 2, color="#334155", lw=1.1)
        extent = max(xpos.values()) or 1
        for leaf in leaves:
            tree_ax.plot([xpos[leaf], extent * 1.025], [ypos[leaf]] * 2, color="#cbd5e1", lw=0.6, ls=":")
        # Restrict inline labels to one useful rank per node; all ranks remain
        # available in the bands, mapping table, and NHX properties.
        for node, annotations in labels.items():
            preferred = next(((rank, name) for rank, name in reversed(annotations)
                              if rank in {"phylum", "class", "order", "family"}), None)
            if preferred and not node.is_root:
                tree_ax.annotate(preferred[1], (xpos[node], ypos[node]), xytext=(3, 5), textcoords="offset points",
                                 fontsize=7, color="#475569", bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 0.5})
        tree_ax.set_xlim(-extent * 0.04, extent * 1.05)
        if use_lengths:
            scale = extent / 5
            tree_ax.plot([0, scale], [-0.18 / plot_height] * 2, transform=tree_ax.get_xaxis_transform(), clip_on=False, color="#334155")
            tree_ax.text(scale / 2, -0.38 / plot_height, f"{scale:.3g} branch-length units", transform=tree_ax.get_xaxis_transform(), ha="center", fontsize=8)
    else:
        tree_ax.text(0.5, 0.5, "No resolved TaxIDs\nTree unavailable", transform=tree_ax.transAxes, ha="center", color="#64748b")
    text_ax.set_xlim(0, 1)
    text_ax.set_title("Input species", loc="left", fontsize=10, fontweight="bold", pad=12)
    for label in display:
        row = by_species.get(label)
        suffix = " †" if label in absent else (" ‡" if row is None else "")
        text_ax.text(0.01, y[label], label.replace("_", " ") + suffix, va="center", fontsize=9,
                     color="#0f172a" if row and row["resolution_status"] == "resolved" else "#64748b")
    palette = ["#b7d9d0", "#f3c6a5", "#bfcbed", "#d9c5e7", "#e8dba4", "#b6d9ea", "#e8bdcb", "#c5dbad"]
    for column, ax in zip(display_columns, axes[2:], strict=True):
        categories = sorted({row[column["taxid_column"]] for row in rows if row[column["taxid_column"]]})
        colors = {key: palette[index % len(palette)] for index, key in enumerate(categories)}
        ax.set_xlim(0, 1)
        heading = column["label"]
        if column["rank"] == "clade":
            heading += f"\nclade ({column['taxid']})"
        ax.set_title(heading, fontsize=9, fontweight="bold", pad=12)
        for label in display:
            row = by_species.get(label, {})
            value = display_value(row.get(column["column"], ""))
            color = colors.get(row.get(column["taxid_column"], ""), "#f1f5f9")
            ax.add_patch(Rectangle((0.015, y[label] - 0.43), 0.97, 0.86, facecolor=color, edgecolor="white", linewidth=0.6))
            # Pale colors keep every taxon label readable without a separate legend.
            ink = "#253448" if sum(to_rgb(color)) > 1.5 else "white"
            ax.text(0.5, y[label], value or "—", ha="center", va="center", fontsize=8, color=ink)
    title = "Input species · taxonomic ranks"
    subtitle = "NCBI taxonomy backbone · topology only" if source == "ncbi_taxonomy" else "Species tree · NCBI taxonomic annotation"
    resolved_count = sum(row["resolution_status"] == "resolved" for row in rows)
    fig.text(0.035, 1 - 0.33 / height, title, fontsize=19, fontweight="bold", color="#0f172a")
    fig.text(0.035, 1 - 0.70 / height, f"{subtitle}   |   {len(rows)} inputs · {resolved_count} resolved", fontsize=10, color="#475569")
    fig.text(0.035, 0.24 / height, "† Missing from tree   ‡ Tree tip outside input set   — Rank unavailable\n"
             "Clade labels require identical sampled membership; classification bands retain non-monophyletic groups.",
             fontsize=8, color="#64748b", linespacing=1.6)
    for extension in ("pdf", "svg", "png"):
        fig.savefig(directory / f"taxonomy_tree.{extension}", dpi=160, facecolor="white")
    plt.close(fig)


def run(args):
    from nwkit.output_transaction import output_transaction, validate_output_targets
    from nwkit.util import get_tree_property_names, write_tree

    args.workspace = Path(args.workspace).expanduser().resolve()
    output = Path(args.output_dir).expanduser().absolute() if args.output_dir else args.workspace / "output/species_taxonomy"
    db = Path(args.taxonomy_db or os.environ.get("GG_TAXONOMY_DBFILE", "") or
              args.workspace / "downloads/ete_taxonomy/taxa.sqlite").expanduser().resolve()
    if not db.is_file():
        raise ValueError(f"Taxonomy database not found: {db}. Prepare it through a GeneGalleon entrypoint.")
    ranks = list(dict.fromkeys(rank.strip() for rank in args.ranks.split(",") if rank.strip()))
    if not ranks or any(not re.fullmatch(r"[a-z][a-z_ ]*", rank) for rank in ranks) or ("all" in ranks and ranks != ["all"]):
        raise ValueError("--ranks must be all or comma-separated rank names (letters, spaces, underscores).")
    if any(rank in {"taxid", "input_taxid", "species_key", "species_rank", "scientific_name", "query",
                    "resolution_status", "resolution_source", "tree_status", "tree_tip"} or rank.endswith("_taxid") for rank in ranks):
        raise ValueError("Rank name conflicts with a metadata column.")
    inputs, sources, fasta_inputs = collect_species(args)
    tree_path, source = choose_tree(args)
    sources += [db, Path(__file__), Path(__file__).with_name("format_species_taxonomy.py"),
                Path(__file__).with_name("format_species_constants.py")] + ([tree_path] if tree_path else [])
    targets = [output / name for name in FILES]
    validate_output_targets(targets, follow_symlinks=False)
    for target in targets:
        for path in [*sources, *fasta_inputs]:
            if target.resolve() == path.resolve() or (target.exists() and os.path.samefile(target, path)):
                raise ValueError(f"Output would overwrite an input: {target}")
    fingerprint = dict(inputs=inputs, ranks=ranks, plot_clades=args.plot_clades, tree_source=source,
                       files={str(path.resolve()): sha256(path) for path in sources},
                       versions={name: version(name) for name in ("nwkit", "ete4", "matplotlib")})
    provenance_path = output / "provenance.json"
    if provenance_path.is_file():
        try:
            previous = json.loads(provenance_path.read_text())
            # gg-cache-guard: audited - fingerprint covers inputs/settings/tool versions; every output is SHA-256 checked.
            if previous["fingerprint"] == fingerprint and all(
                (output / name).is_file() and sha256(output / name) == previous["outputs"][name]
                for name in FILES if name != "provenance.json"
            ):
                print(f"Species taxonomy is current: {output}")
                return
        except (ValueError, KeyError, TypeError, OSError):
            pass
    taxonomy = Taxonomy(db)
    try:
        rows, lineages = [], {}
        for item in inputs:
            row, lineage = taxonomy.resolve(item)
            rows.append(row)
            lineages[row["species"]] = lineage
    finally:
        taxonomy.conn.close()
    if ranks == ["all"]:
        ranks = available_ranks(lineages)
    tree = load_or_build_tree(tree_path, lineages)
    mappings, labels = map_taxonomy(tree, rows, lineages, ranks)
    columns = aligned_columns(rows, lineages, ranks, tree)
    # Build and validate every artifact before replacing any public output.
    with tempfile.TemporaryDirectory(prefix="gg-species-taxonomy-") as temporary:
        staged = Path(temporary)
        fields = ["species", "input_taxid", "taxid", "scientific_name", "query", "resolution_status", "resolution_source", "tree_status", "tree_tip"]
        # 'species' is the persistent input key; the species rank uses a separate name.
        write_tsv(staged / FILES[0], rows, fields + [column["column"] for column in columns] + [column["taxid_column"] for column in columns])
        write_tsv(staged / "taxonomy_columns.tsv", columns,
                  ["position", "column", "taxid_column", "rank", "taxid", "label", "unmet_predecessors"])
        lineage_rows = [dict(species=label, lineage_index=index, taxid=node["taxid"], name=node["spname"], rank=node["rank"])
                        for label, lineage in lineages.items() for index, node in enumerate(lineage)]
        write_tsv(staged / FILES[1], lineage_rows, ["species", "lineage_index", "taxid", "name", "rank"])
        write_tsv(staged / FILES[2], mappings, ["rank", "taxid", "name", "status", "node_id", "mapped_count", "missing_count", "members", "missing_members", "other_descendants", "unclassified_descendants"])
        if tree is None:
            (staged / FILES[3]).write_text("")
            (staged / FILES[4]).write_text("")
        else:
            if tree_path:
                shutil.copyfile(tree_path, staged / FILES[3])
            else:
                write_tree(tree, SimpleNamespace(outfile=str(staged / FILES[3])), format=9, quiet=True, props=[])
            write_tree(tree, SimpleNamespace(outfile=str(staged / FILES[4])), format="auto" if tree_path else 9,
                       quiet=True, props=get_tree_property_names(tree))
        plot(tree, rows, columns, labels, source, staged, args.plot_clades)
        provenance = dict(created_utc=datetime.now(timezone.utc).isoformat(), fingerprint=fingerprint,
                          resolved_ranks=ranks,
                          tree_path=str(tree_path.resolve()) if tree_path else None, tree_source=source,
                          tree_available=tree is not None, input_count=len(inputs),
                          taxonomy_db_mtime_utc=datetime.fromtimestamp(db.stat().st_mtime, timezone.utc).isoformat(),
                          unresolved=[row["query"] for row in rows if row["resolution_status"] != "resolved"],
                          outputs={name: sha256(staged / name) for name in FILES if name != "provenance.json"})
        (staged / FILES[-1]).write_text(json.dumps(provenance, indent=2) + "\n")
        if any(sha256(path) != fingerprint["files"][str(path.resolve())] for path in sources):
            raise ValueError("An input or taxonomy database changed during generation; rerun the stage.")
        with output_transaction(targets, create_parents=True) as pending:
            for target in targets:
                shutil.copyfile(staged / target.name, pending[target])
    print(f"Species taxonomy ({source}): {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", default="workspace")
    parser.add_argument("--species-table", default="", help="Authoritative species/species_key/leaf_name TSV, optional taxid.")
    parser.add_argument("--species-summary", default="", help="Metadata for current input FASTAs; never adds historical species.")
    parser.add_argument("--species-dir", action="append", default=[], help="Current input FASTA directory; repeatable.")
    parser.add_argument("--taxid-map", default="", help="Explicit species-to-TaxID override TSV.")
    parser.add_argument("--taxid-override", default="", help="One explicit species:TaxID correction for a scheduled task.")
    parser.add_argument("--taxonomy-db", default="")
    parser.add_argument("--species-tree", default="auto")
    parser.add_argument("--ranks", default=DEFAULT_RANKS, help="all (default) includes every observed lineage rank, or specify a comma-separated list.")
    parser.add_argument("--plot-clades", type=int, choices=(0, 1), default=0,
                        help="Draw clade columns (1); hidden by default (0). Tables and NHX retain clades.")
    parser.add_argument("--output-dir", default="")
    args = parser.parse_args()
    try:
        run(args)
    except (ValueError, OSError, sqlite3.Error) as exc:
        print(f"Species taxonomy: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

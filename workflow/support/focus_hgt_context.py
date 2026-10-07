"""Existing-coordinate donor/recipient context pages; no sequence reanalysis."""

import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path

ORANGE = "#b34d00"
BLUE = "#2b6ca3"
CONTEXT_MAX_GENES_PER_SIDE = 3
CONTEXT_MIN_FLANK_GENES = 2
CONTEXT_MAX_OVERLAPPING_GENES = 2
NONCODING_GAP_THRESHOLD_BP = 5000
NONCODING_GAP_DISPLAY_BP = 2000


def support_label(branch):
    from focus_hgt_gene_trees import number
    support = number(branch.get('support_generax_ufboot'))
    return f'{support:g}' if support is not None else 'NA'


def blocks(text):
    if not text or text.lower() in {"na", "nan"}:
        return []
    pairs = []
    for item in text.split(";"):
        match = re.fullmatch(r"(\d+)-(\d+)", item)
        if not match:
            raise ValueError("Invalid existing GFF block: " + item)
        start, end = map(int, match.groups())
        if start < 1 or start > end:
            raise ValueError("Invalid existing GFF block bounds")
        pairs.append((start, end))
    return sorted(pairs)


def structure(row):
    """Keep coding blocks and UTRs distinct; introns use only cis coordinates."""
    coding = blocks(row.get("feature_blocks", ""))
    utr = blocks(row.get("utr_blocks", ""))
    if row.get("splice_mode", "cis") != "cis":
        return dict(status="trans_splicing_not_drawn_on_single_scaffold", coding=[], utr=[], exon=[], introns=[])
    if row.get("feature_type") not in {"CDS", "exon"} or not coding:
        return dict(status="exon_coordinates_unavailable", coding=[], utr=[], exon=[], introns=[])
    start, end = int(row["start"]), int(row["end"])
    # gff_info start/end span the recorded CDS/exon features. UTR coordinates
    # legitimately extend beyond a CDS span and must retain their own bounds.
    if any(a < start or b > end for a, b in coding):
        raise ValueError("GFF feature block outside focal feature span")
    merged = []
    for a, b in sorted(coding + utr):
        if merged and a <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(b, merged[-1][1]))
        else:
            merged.append((a, b))
    introns = [
        (left[1] + 1, right[0] - 1) for left, right in zip(merged, merged[1:], strict=False) if right[0] > left[1] + 1
    ]
    return dict(
        status=("annotated_exons_CDS_UTR_unavailable" if row['feature_type'] == 'exon'
                else "coding_exons_with_utr" if utr else "coding_exons_utr_unavailable"),
        coding=coding if row["feature_type"] == "CDS" else [],
        utr=utr if row["feature_type"] == "CDS" else [],
        exon=coding if row['feature_type'] == 'exon' else [],
        introns=introns,
    )


def neighbor_relation(row, focal):
    if focal is None:
        return 'unavailable'
    if row['gene_id'] == focal['gene_id']:
        return 'focal'
    if int(row['end']) < int(focal['start']):
        return 'left'
    if int(row['start']) > int(focal['end']):
        return 'right'
    return 'overlapping'


def neighbor_counts(focal, neighbors):
    counts = {side: sum(neighbor_relation(row, focal) == side for row in neighbors)
              for side in ('left', 'right', 'overlapping')} if focal else {}
    return {f'neighbor_{side}_{field}': value
            for side in ('left', 'right', 'overlapping')
            for field, value in (
                ('gene_count', counts.get(side, '')),
                ('status', 'gff_unavailable' if not focal else 'shown' if side == 'overlapping'
                 else 'minimum_met' if counts[side] >= CONTEXT_MIN_FLANK_GENES else 'insufficient_annotated_loci'))}


def model_span(row):
    utr = blocks(row.get('utr_blocks', ''))
    return (min([int(row['start'])] + [a for a, _ in utr]),
            max([int(row['end']) + 1] + [b + 1 for _, b in utr]))


class GapCompressedCoordinates:
    """Monotone map; recorded exon/UTR blocks and unknown spans never shrink."""

    def __init__(self, focal, neighbors):
        self.gaps, self.center, self.anchor = [], 0, 0
        if focal is None:
            return
        spans, protected = [], []
        for row in neighbors:
            spans.append(model_span(row))
            info = structure(row)
            recorded = info['coding'] + info['exon'] + info['utr']
            protected.extend([(a, b + 1) for a, b in recorded] if recorded else [spans[-1]])
        merged = []
        for start, end in sorted(protected):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
            else:
                merged.append((start, end))
        for left, right in zip(merged, merged[1:], strict=False):
            length = right[0] - left[1]
            if length > NONCODING_GAP_THRESHOLD_BP:
                within_gene = any(a <= left[1] and right[0] <= b for a, b in spans)
                mixed = any(a < right[0] and b > left[1] for a, b in spans)
                self.gaps.append(dict(gap_type='intronic' if within_gene else 'mixed_noncoding' if mixed else 'intergenic',
                                      genomic_start_bp=left[1], genomic_end_exclusive_bp=right[0],
                                      original_gap_bp=length, display_gap_bp=NONCODING_GAP_DISPLAY_BP,
                                      omitted_bp=length - NONCODING_GAP_DISPLAY_BP))
        self.center = (int(focal['start']) + int(focal['end'])) / 2
        self.anchor = self.transform(self.center)
        self.bounds = (min(start for start, _ in spans), max(end for _, end in spans))

    def transform(self, position):
        reduction = 0
        for gap in self.gaps:
            covered = min(max(position - gap['genomic_start_bp'], 0), gap['original_gap_bp'])
            reduction += covered * gap['omitted_bp'] / gap['original_gap_bp']
        return position - reduction

    def point(self, position):
        return (self.transform(position) - self.anchor) / 1000

    def audit(self, extent=None):
        gaps = [gap for gap in self.gaps if extent is None or
                self.point(gap['genomic_start_bp']) <= extent and self.point(gap['genomic_end_exclusive_bp']) >= -extent]
        return [dict(gap, break_label=str(i), display_start_kb=self.point(gap['genomic_start_bp']),
                     display_end_kb=self.point(gap['genomic_end_exclusive_bp']))
                for i, gap in enumerate(gaps, 1)]


class GenomeCoordinates:
    def __init__(self, root):
        self.root = Path(root) if root else None
        self.cache = {}
        self.sources = {}

    def load(self, species):
        if species in self.cache:
            return self.cache[species]
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", species):
            raise ValueError("Unsafe GFF species identifier")
        path = self.root / (species + ".gff_info.tsv") if self.root else None
        rows = []
        if path and path.is_file():
            raw = path.read_bytes()
            self.sources[str(path.resolve())] = hashlib.sha256(raw).hexdigest()
            import io

            reader = csv.DictReader(io.StringIO(raw.decode('utf-8-sig')), delimiter="\t")
            fields = reader.fieldnames or []
            if len(set(fields)) != len(fields) or not {'gene_id', 'chromosome', 'start', 'end'} <= set(fields):
                raise ValueError('Missing or duplicate GFF coordinate columns: ' + species)
            rows = list(reader)
            if any(None in row or any(v is None for v in row.values()) for row in rows):
                raise ValueError('Malformed GFF coordinate row: ' + species)
        by_gene = {r["gene_id"]: r for r in rows}
        if len(by_gene) != len(rows) or '' in by_gene:
            raise ValueError("Duplicate GFF gene identity: " + species)
        by_scaffold = defaultdict(list)
        from focus_hgt_context_annotations import available
        for row in rows:
            if all(available(row.get(k)) for k in ('chromosome', 'start', 'end')):
                if int(row['start']) < 1 or int(row['end']) < int(row['start']):
                    raise ValueError('Invalid GFF coordinate span: ' + row['gene_id'])
                by_scaffold[row["chromosome"]].append(row)
        self.cache[species] = by_gene, by_scaffold
        return self.cache[species]

    def neighborhood(self, link):
        by_gene, scaffolds = self.load(link.get("gene_species", ""))
        row = by_gene.get(link["gene_id"])
        if row is None:
            return None, [], "gff_gene_unavailable"
        from focus_hgt_context_annotations import available
        if not all(available(row.get(k)) for k in ('chromosome', 'start', 'end')):
            return None, [], 'gff_coordinates_unavailable'
        if link.get("host_scaffold_id") and row["chromosome"] != link["host_scaffold_id"]:
            raise ValueError("Event-gene scaffold and GFF scaffold disagree: " + link["gene_id"])
        candidates = [r for r in scaffolds[row["chromosome"]] if r["gene_id"] != row["gene_id"]]
        left = sorted(
            [r for r in candidates if int(r["end"]) < int(row["start"])],
            key=lambda r: (-int(r["end"]), r["gene_id"]),
        )[:CONTEXT_MIN_FLANK_GENES]
        right = sorted(
            [r for r in candidates if int(r["start"]) > int(row["end"])],
            key=lambda r: (int(r["start"]), r["gene_id"]),
        )[:CONTEXT_MIN_FLANK_GENES]
        center = (int(row['start']) + int(row['end'])) / 2
        # Intronic/nested loci are additional context, never substitutes for
        # the required left/right flanks. Missing scaffold annotations stay missing.
        overlapping = sorted([r for r in candidates if int(r['start']) <= int(row['end'])
                              and int(r['end']) >= int(row['start'])],
                             key=lambda r: (abs((int(r['start']) + int(r['end'])) / 2 - center), r['gene_id']))[:CONTEXT_MAX_OVERLAPPING_GENES]
        selected = left + right + overlapping
        return row, sorted(selected + [row], key=lambda r: (int(r['start']), r['gene_id'])), structure(row)["status"]

    def scaffold_models(self, link, focal):
        """All coordinate-bearing models, including outer and nested loci."""
        if focal is None:
            return []
        _, scaffolds = self.load(link.get('gene_species', ''))
        return sorted(scaffolds[focal['chromosome']], key=lambda r: (int(r['start']), r['gene_id']))

    def verify(self):
        for path, expected in self.sources.items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != expected:
                raise ValueError("GFF coordinate input changed during focused rendering")


def choose_representatives(events, links, coordinates):
    from focus_hgt_gene_trees import background_supported

    selected, audit = [], []
    for event in events:
        for side in ("donor", "recipient"):
            passing = [
                r for r in links if r["event_id"] == event["event_id"] and r["side"] == side and background_supported(r)
            ]
            ranked = sorted(
                passing,
                key=lambda r: (
                    coordinates.neighborhood(r)[0] is None,
                    -float(r["host_scaffold_background_class_classified_fraction"]),
                    -float(r["host_scaffold_background_class_compatible_fraction"]),
                    r["gene_id"],
                ),
            )
            if not ranked:
                raise ValueError("Context page requires passing event-linked gene on each side")
            link = ranked[0]
            focal, neighbors, status = coordinates.neighborhood(link)
            selected.append((event, side, link, focal, neighbors))
            audit.append(
                dict(
                    event_id=event["event_id"],
                    side=side,
                    gene_id=link["gene_id"],
                    gene_species=link.get("gene_species", ""),
                    scaffold=link["host_scaffold_id"],
                    passing_gene_count=len(passing),
                    representative_rule="available_gff_then_coverage_compatibility_gene_id",
                    structure_status=status,
                    feature_blocks=focal.get("feature_blocks", "") if focal else "",
                    utr_blocks=focal.get("utr_blocks", "") if focal else "",
                    intron_count=focal.get("num_intron", "") if focal else "",
                    neighbor_gene_ids="; ".join(
                        r["gene_id"]
                        for r in sorted(neighbors, key=lambda r: int(r["start"]))
                        if r["gene_id"] != link["gene_id"]
                    ),
                    coordinate_unit="genomic_bp_1_based_inclusive",
                )
            )
    return selected, audit


def selected_gene_clade(ax, rows, event, donor, recipient):
    """A compact view of existing paths, preserving the exact HGT node."""
    by_id = {r["branch_id"]: r for r in rows}
    by_name = {r["node_name"]: r["branch_id"] for r in rows}
    node = event.get("gene_tree_branch_id", event.get("branch_id"))

    def ancestors(bid):
        path = []
        while bid in by_id:
            if bid in path:
                raise ValueError("Cycle in gene-tree topology")
            path.append(bid)
            bid = by_id[bid].get("parent")
        return path

    paths = [ancestors(by_name[name]) for name in (donor, recipient)] + [ancestors(node)]
    root = next(bid for bid in paths[0] if all(bid in path for path in paths))
    kept = {bid for path in paths for bid in path[: path.index(root) + 1]}
    children = {bid: [x for x in (by_id[bid].get("child1"), by_id[bid].get("child2")) if x in kept] for bid in kept}
    ys, xs, edges = {}, {}, []
    tip_index = [0]

    def layout(bid, x):
        xs[bid] = x
        child = children[bid]
        if not child:
            ys[bid] = tip_index[0]
            tip_index[0] += 1
        else:
            for c in child:
                value = float(by_id[c].get("bl_rooted") or 0)
                if not math.isfinite(value) or value < 0:
                    raise ValueError("Invalid existing gene-tree branch length")
                layout(c, x + value)
                edges.append((bid, c))
            ys[bid] = sum(ys[c] for c in child) / len(child)

    layout(root, 0)
    for parent, child in edges:
        ax.plot([xs[parent], xs[child]], [ys[child], ys[child]], color=BLUE, lw=1)
        ax.plot([xs[parent], xs[parent]], [ys[parent], ys[child]], color=BLUE, lw=1)
    maximum = max(xs.values()) or 1
    for name, color, side in [(donor, BLUE, "donor"), (recipient, ORANGE, "recipient")]:
        bid = by_name[name]
        ax.text(xs[bid] + maximum * 0.02, ys[bid], name + " [" + side + "]", color=color, fontsize=8, va="center")
    ax.scatter([xs[node]], [ys[node]], marker="D", s=25, color=ORANGE, zorder=5)
    ax.text(
        0.01,
        1.02,
        f"{event['event_id']} | HGT node {by_id[node]['node_name']} | UFB {support_label(by_id[node])}",
        transform=ax.transAxes,
        fontsize=9,
        color=ORANGE,
    )
    ax.set_xlim(-maximum * 0.02, maximum * 1.9)
    ax.set_ylim(-0.5, tip_index[0] - 0.5)
    ax.set_yticks([])
    ax.set_xlabel("Existing gene-tree branch length (substitutions/site)", fontsize=8)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(labelsize=7)


def choose_context_genes(events, links, coordinates, max_genes_per_side=CONTEXT_MAX_GENES_PER_SIDE):
    """Cap distinct genes per modeled side, retaining every event/gene in the audit."""
    from focus_hgt_gene_trees import background_supported, number

    if (type(max_genes_per_side) is not int or not 1 <= max_genes_per_side <= CONTEXT_MAX_GENES_PER_SIDE):
        raise ValueError(f"Context display limit must be an integer in [1, {CONTEXT_MAX_GENES_PER_SIDE}]")
    by_event = {e['event_id']: e for e in events}
    if not by_event or len(by_event) != len(events):
        raise ValueError('Context requires nonempty, distinct event IDs')
    grouped = {}
    passing_by_event = defaultdict(set)
    prefix = 'host_scaffold_background_class_'
    for link in links:
        event_id = link['event_id']
        if event_id not in by_event or str(link.get('eligible_for_context', '')).lower() not in {'true', '1'}:
            continue
        side = link['side']
        if side not in {'donor', 'recipient'} or link['orthogroup'] != by_event[event_id]['orthogroup']:
            raise ValueError('Context event/gene side or family is inconsistent')
        supported = background_supported(link)
        counts = [number(link.get(prefix + key)) for key in
                  ('total_count', 'compatible_count', 'incompatible_count', 'unresolved_count')]
        measured = (link.get('host_scaffold_status') == 'measured' and bool(link.get('host_scaffold_id'))
                    and all(value is not None for value in counts))
        status = ('scaffold_supported' if supported else 'scaffold_thresholds_not_met' if measured
                  else 'scaffold_evidence_unavailable')
        signature = {key: value for key, value in link.items()
                     if key.startswith('host_scaffold_') or key == 'gene_species'}
        key = side, link['gene_id']
        if key in grouped and grouped[key]['signature'] != signature:
            raise ValueError('Conflicting scaffold evidence for the same context gene')
        if key not in grouped:
            focal, neighbors, structure_status = coordinates.neighborhood(link)
            grouped[key] = dict(side=side, link=link, event_ids=set(), signature=signature,
                                supported=supported, status=status, focal=focal,
                                neighbors=neighbors, structure_status=structure_status)
        grouped[key]['event_ids'].add(event_id)
        if supported:
            passing_by_event[event_id].add(side)
    if any(passing_by_event[event_id] != {'donor', 'recipient'} for event_id in by_event):
        raise ValueError('Context page requires a passing event-linked gene on each side of every event')
    selected, audit, totals = [], [], {}
    for side in ('donor', 'recipient'):
        ordered = sorted([entry for entry in grouped.values() if entry['side'] == side], key=lambda entry: (
            not entry['supported'], entry['focal'] is None,
            -(number(entry['link'].get(prefix + 'classified_fraction')) or 0),
            -(number(entry['link'].get(prefix + 'compatible_fraction')) or 0), entry['link']['gene_id']))
        totals[side] = dict(total=len(ordered), shown=min(len(ordered), max_genes_per_side),
                            omitted=max(0, len(ordered) - max_genes_per_side),
                            supported=sum(entry['supported'] for entry in ordered))
        for rank, entry in enumerate(ordered, 1):
            link, focal, neighbors = entry['link'], entry['focal'], entry['neighbors']
            displayed = rank <= max_genes_per_side
            if displayed:
                selected.append(entry)
            for event_id in sorted(entry['event_ids']):
                audit.append(dict(
                    event_id=event_id, side=side, gene_id=link['gene_id'], gene_species=link.get('gene_species', ''),
                    scaffold=focal['chromosome'] if focal else link.get('host_scaffold_id', ''),
                    scaffold_basis='existing_gff' if focal else 'event_gene_summary' if link.get('host_scaffold_id') else '',
                    passing_gene_count=sum(g['supported'] and g['side'] == side and event_id in g['event_ids']
                                           for g in grouped.values()),
                    representative_rule='supported_then_available_gff_then_coverage_compatibility_gene_id',
                    structure_status=entry['structure_status'],
                    feature_blocks=focal.get('feature_blocks', '') if focal else '',
                    utr_blocks=focal.get('utr_blocks', '') if focal else '',
                    intron_count=focal.get('num_intron', '') if focal else '',
                    neighbor_gene_ids='; '.join(r['gene_id'] for r in sorted(neighbors, key=lambda r: int(r['start']))
                                               if r['gene_id'] != link['gene_id']),
                    **neighbor_counts(focal, neighbors),
                    coordinate_unit='genomic_bp_1_based_inclusive', scaffold_support_status=entry['status'],
                    scaffold_count_unit=link.get('host_scaffold_count_unit', ''),
                    **{prefix + suffix: link.get(prefix + suffix, '') for suffix in
                       ('total_count', 'compatible_count', 'incompatible_count', 'unresolved_count',
                        'classified_fraction', 'compatible_fraction')},
                    displayed=int(displayed), display_rank=rank, max_genes_per_side=max_genes_per_side,
                    display_reason='within_side_limit' if displayed else 'side_display_limit',
                    side_total_gene_count=totals[side]['total'], side_shown_gene_count=totals[side]['shown'],
                    side_omitted_gene_count=totals[side]['omitted'], side_supported_gene_count=totals[side]['supported']))
    return selected, audit, totals


def model_lanes(models, focal, display):
    """Pack full model spans; overlapping models never share a vertical lane."""
    if focal is None:
        return {}
    occupied, result = {0: [tuple(display.point(x) for x in model_span(focal))]}, {focal['gene_id']: 0}
    for row in models:
        if row['gene_id'] == focal['gene_id']:
            continue
        start, end = (display.point(x) for x in model_span(row))
        for index in range(2 * len(models) + 1):
            lane = 0 if index == 0 else (index + 1) // 2 * (1 if index % 2 else -1)
            if all(end <= a or start >= b for a, b in occupied.get(lane, [])):
                occupied.setdefault(lane, []).append((start, end))
                result[row['gene_id']] = lane * .48
                break
    return result


def prepare_context_models(selected, coordinates):
    """Use all scaffold exons for compression, then draw every model in view."""
    for entry in selected:
        models = coordinates.scaffold_models(entry['link'], entry['focal'])
        entry['scaffold_models'] = models
        entry['display_coordinates'] = GapCompressedCoordinates(entry['focal'], models)
    extent = math.ceil(max([20.0] + [max(abs(entry['display_coordinates'].point(x))
                     for row in entry['neighbors'] for x in model_span(row)) + 2
                     for entry in selected if entry['focal']]) / 5) * 5
    for entry in selected:
        display = entry['display_coordinates']
        entry['models'] = [row for row in entry['scaffold_models']
                           if display.point(model_span(row)[0]) <= extent and display.point(model_span(row)[1]) >= -extent]
        entry['model_lanes'] = model_lanes(entry['models'], entry['focal'], display)
        entry['model_height_pt'] = max(62, 24 + 18 * len(set(entry['model_lanes'].values())))
    return extent


def draw_context_neighborhood(ax, entry, extent, color):
    """Shared genic scale with marked intergenic omissions and separate overlap lanes."""
    from matplotlib.patches import Rectangle

    link, focal, neighbors = entry['link'], entry['focal'], entry['neighbors']
    ax.set_xlim(-extent, extent)
    levels = list(entry.get('model_lanes', {}).values()) or [0]
    ax.set_ylim(min(levels) - 1, max(levels) + .6)
    ax.set_yticks([])
    ax.axvline(0, color='#cccccc', lw=0.5, zorder=0)
    if focal is None:
        ax.text(0, 0, 'GFF coordinates unavailable', ha='center', color='#777777', fontsize=8)
    else:
        display = entry.get('display_coordinates') or GapCompressedCoordinates(focal, neighbors)
        numbered = [r['gene_id'] for r in sorted(neighbors, key=lambda r: int(r['start'])) if r['gene_id'] != link['gene_id']]
        numbers = {gene: str(i) for i, gene in enumerate(numbered, 1)}
        models = entry.get('models', neighbors)
        lanes = entry.get('model_lanes') or model_lanes(models, focal, display)
        ordered = enumerate(sorted(models, key=lambda r: (int(r['start']), r['gene_id'])))
        # Keep genomic label positions, but paint the focal structure last so
        # an overlapping neighbor cannot obscure its donor/recipient color.
        for j, row in sorted(ordered, key=lambda item: item[1]['gene_id'] == link['gene_id']):
            focal_flag = row['gene_id'] == link['gene_id']
            level = lanes.get(row['gene_id'], 0)
            edge = color if focal_flag else '#858585'
            face = edge if not focal_flag or entry['supported'] else '#f6e8df' if entry['side'] == 'recipient' else '#e1edf5'
            hatch = '///' if focal_flag and not entry['supported'] else None
            a, b = display.point(int(row['start'])), display.point(int(row['end']) + 1)
            info = structure(row)
            for x, y in info['introns']:
                ax.plot([display.point(x), display.point(y + 1)], [level, level], color='#444444', lw=0.7)
            for kind, height in [('coding', 0.20), ('exon', 0.20), ('utr', 0.10)]:
                for x, y in info[kind]:
                    ax.add_patch(Rectangle((display.point(x), level - height / 2), display.point(y + 1) - display.point(x),
                                           height, facecolor=face if kind == 'coding' else '#f0f0f0' if kind == 'exon' else '#a6adb2',
                                           edgecolor=edge, lw=0.6, hatch=('..' + (hatch or '')) if kind == 'exon' else hatch))
            if not any(info[kind] for kind in ('coding', 'utr', 'exon')):
                ax.add_patch(Rectangle((a, level - 0.1), b - a, 0.2, fill=False, ec=edge, ls=':', lw=0.8))
            direction = 1 if row.get('strand') == '+' else -1 if row.get('strand') == '-' else 0
            if direction:
                endpoint = b if direction == 1 else a
                ax.annotate('', xy=(endpoint, level), xytext=(endpoint - direction * min(0.4, max(0.05, b - a)), level),
                            arrowprops=dict(arrowstyle='->', color=edge, lw=0.7))
            if focal_flag:
                label = 'Focal'
            else:
                label = numbers.get(row['gene_id'], '')
            if label and (a < -extent or b > extent):
                label += '*'
            label_level = level + (.25 if level >= 0 else -.25) if level else (.27 if j % 2 == 0 else -.27)
            if label:
                ax.text((max(a, -extent) + min(b, extent)) / 2, label_level,
                        label, ha='center', va='center', fontsize=7, color=edge,
                        weight='bold' if focal_flag else 'normal')
        for gap in display.audit(extent):
            ax.text((max(-extent, gap['display_start_kb']) + min(extent, gap['display_end_kb'])) / 2, min(levels) - .86,
                    '//' + gap['break_label'], ha='center', va='center', fontsize=7, color='#666666')
        ax.text(0.01, 0.98, structure(focal)['status'].replace('_', ' '), transform=ax.transAxes,
                fontsize=7, color='#777777', va='top')
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.tick_params(labelsize=7)
    ax.set_xlabel('Compressed display position (kb-equivalent; focal-feature midpoint = 0)', fontsize=8)


def render_bounded_context(path, rows, events, links, coordinates, max_genes_per_side=CONTEXT_MAX_GENES_PER_SIDE,
                           annotations=None):
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from focus_hgt_context_annotations import ContextAnnotations, context_annotation_rows, draw_annotation_table
    from focus_hgt_gene_trees import number

    selected, audit, totals = choose_context_genes(events, links, coordinates, max_genes_per_side)
    tips = {r['node_name'] for r in rows if r['child1'] == r['child2'] == '-999'}
    if any(row['gene_id'] not in tips for row in audit):
        raise ValueError('Context gene is not a tip of its orthogroup gene tree')
    branches = {r['branch_id']: r for r in rows}
    ordered_events = sorted(events, key=lambda e: e['event_id'])
    labels, references = {}, []
    for i, event in enumerate(ordered_events, 1):
        branch = branches[event.get('gene_tree_branch_id', event.get('branch_id'))]
        if branch['node_name'] != event.get('gene_tree_node', event.get('node_name')):
            raise ValueError('Context transfer branch/node mapping is inconsistent')
        labels[event['event_id']] = f'HGT{i}'
        references.append(f"HGT{i}: node {branch['node_name']} | UFB {support_label(branch)}")
    extent = prepare_context_models(selected, coordinates)
    annotations = annotations if annotations is not None else ContextAnnotations()
    family = ordered_events[0]['orthogroup']
    leaves = {r['node_name']: r for r in rows if r['child1'] == r['child2'] == '-999'}
    entries = {side: [e for e in selected if e['side'] == side] for side in totals}
    tables = {id(e): context_annotation_rows(e, annotations, family, leaves) for e in selected}
    for entry in selected:
        annotations.display_audit.extend({k: v for k, v in r.items() if k not in {'cells', 'height_pt'}}
                                          for r in tables[id(entry)])
    nrows = max(totals[side]['shown'] for side in totals)
    plot_heights = [max(entries[side][i]['model_height_pt'] for side in entries if i < len(entries[side])) for i in range(nrows)]
    row_heights = [83 + plot_heights[i] + max(25 + sum(r['height_pt'] for r in tables[id(entries[side][i])])
                            for side in entries if i < len(entries[side])) + 30 for i in range(nrows)]
    page_height = 165 + sum(row_heights) + 115
    fig = plt.figure(figsize=(22, page_height / 72))
    def y(point):
        return 1 - point / page_height
    fig.suptitle('A traceable candidate with donor and recipient context', fontsize=17, x=0.05, ha='left', y=y(23))
    fig.text(0.05, y(55), f"{family} | {len(events)} modeled transfer event(s) | "
             f"at most {max_genes_per_side} distinct genes per side; existing GFF coordinates and annotations", fontsize=10)
    reference_text = '; '.join(references[:4])
    if len(references) > 4:
        reference_text += f'; +{len(references) - 4} events (see context audit)'
    fig.text(0.05, y(78), reference_text, fontsize=9)
    for column, (side, color, title) in enumerate([('donor', BLUE, 'DONOR DESCENDANTS'), ('recipient', ORANGE, 'RECIPIENT DESCENDANTS')]):
        left = 0.05 if column == 0 else 0.535
        fig.text(left, y(115), title, color=color, fontsize=15, weight='bold')
        count = totals[side]
        fig.text(left, y(139), f"Shown {count['shown']} of {count['total']} genes | {count['omitted']} omitted | "
                 f"{count['supported']} scaffold-supported in total", fontsize=9, color=color)
        top = 165
        for index, entry in enumerate(entries[side]):
            plot_height = plot_heights[index]
            ax = fig.add_axes([left, y(top+53+plot_height), .415, plot_height / page_height])
            draw_context_neighborhood(ax, entry, extent, color)
            link = entry['link']
            prefix = 'host_scaffold_background_class_'
            coverage, compatible = [number(link.get(prefix + suffix)) for suffix in ('classified_fraction', 'compatible_fraction')]
            total, host, other = [number(link.get(prefix + suffix)) for suffix in
                                  ('total_count', 'compatible_count', 'incompatible_count')]
            classified = host + other if host is not None and other is not None else None
            def ratio(value, numerator, denominator):
                fraction = f'{value:.1%}' if value is not None else 'undefined' if denominator == 0 else 'unavailable'
                return (f'{numerator:g}/{denominator:g} ({fraction})'
                        if numerator is not None and denominator is not None else fraction)
            unit = {'gff_locus': 'GFF loci', 'cds_id': 'CDS IDs'}.get(link.get('host_scaffold_count_unit'),
                                                                  link.get('host_scaffold_count_unit') or 'count unit unavailable')
            measured = (f'coverage {ratio(coverage, classified, total)} | '
                        f'compatible {ratio(compatible, host, classified)} | {unit}')
            event_tags = [labels[eid] for eid in sorted(entry['event_ids'])]
            tags = ', '.join(event_tags[:3]) + (f', +{len(event_tags)-3} events' if len(event_tags) > 3 else '')
            status = entry['status'].replace('_', ' ')
            scaffold = entry['focal']['chromosome'] if entry['focal'] else link.get('host_scaffold_id') or 'unavailable'
            counts = neighbor_counts(entry['focal'], entry['neighbors'])
            flank_label = ('Flanks unavailable' if entry['focal'] is None else
                           f"Left {counts['neighbor_left_gene_count']}/2 | Right {counts['neighbor_right_gene_count']}/2 | "
                           f"Overlapping {counts['neighbor_overlapping_gene_count']}")
            gaps = entry['display_coordinates'].audit(extent)
            omitted = ', '.join(f"{kind} {sum(g['omitted_bp'] for g in gaps if g['gap_type'] == kind)/1000:.3g}kb"
                                for kind in ('intergenic', 'intronic', 'mixed_noncoding') if any(g['gap_type'] == kind for g in gaps))
            if omitted:
                flank_label += ' | Omitted gaps: ' + omitted
            ax.set_title(f"{tags} | {link['gene_id']}\nScaffold {scaffold} | {status}\n{measured}\n"
                         f"All models in view: {len(entry['models'])} | {flank_label}",
                         loc='left', fontsize=9, color=color, pad=10)
            table_rows = tables[id(entry)]
            table_height = 25 + sum(r['height_pt'] for r in table_rows)
            table_ax = fig.add_axes([left, y(top+83+plot_height+table_height), .435, table_height / page_height])
            draw_annotation_table(table_ax, table_rows, color)
            top += row_heights[index]
    fig.text(0.05, y(page_height-23),
             'Blue: donor descendant focal gene; orange: recipient descendant focal gene; gray: nearby annotated loci. Pale hatched focal blocks: scaffold support not established.\n'
             'Thick blocks: CDS; thin gray blocks: recorded UTR; gray dotted blocks: exons with unknown CDS/UTR identity; lines: introns. Overlapping loci use separate vertical lanes.\n'
             'Every annotated model intersecting the display range is drawn. Only numbered neighbors (nearest two left/right plus up to two overlaps) appear in the annotation table.\n'
             'Shared exon/UTR kb scale: recorded blocks remain uncompressed. Noncoding gaps >5 kb (intergenic or intronic) are capped at 2 kb; numbered // marks identify omissions.\n'
             'Display priority: scaffold-supported, available GFF, background coverage, host compatibility, gene ID. Counts are distinct genes per side, not acquisitions.\n'
             'UFB = Ultrafast bootstrap. Candidate-free class background: at least 10 classified units, 50% coverage, 90% host compatibility. Best-hit taxonomy is annotation, not the modeled transfer donor.\n'
             'MMseqs2 shows each focal/neighbor query classification and its saved host-class/species match; unavailable and unresolved are distinct. These labels do not alter candidate selection.\n'
             'Swiss-Prot best hits provide predicted products and hit taxonomy, separately from MMseqs2 query classification. Full sources and gene/event mappings are in the annotation audit.\n'
             'Distances across // marks and long introns are compressed, not physical genomic distances. Titles give total omitted lengths; audits retain each original interval and length. Neighbor order does not establish conserved synteny.',
             fontsize=8, color='#666666')
    fig.savefig(path, format='pdf')
    plt.close(fig)
    annotations.verify()
    by_id = {e['event_id']: e for e in events}
    displayed = {(e['side'], e['link']['gene_id']): e for e in selected}
    for row in audit:
        event = by_id[row['event_id']]
        branch_id = event.get('gene_tree_branch_id', event.get('branch_id'))
        row.update(shared_axis_min_kb=-extent, shared_axis_max_kb=extent,
                   shared_axis_unit='kb_equivalent_after_noncoding_gap_compression',
                   noncoding_gap_threshold_bp=NONCODING_GAP_THRESHOLD_BP,
                   noncoding_gap_display_bp=NONCODING_GAP_DISPLAY_BP,
                   drawn_gene_model_ids='; '.join(r['gene_id'] for r in displayed[row['side'], row['gene_id']]['models'])
                   if (row['side'], row['gene_id']) in displayed else '',
                   drawn_gene_model_count=len(displayed[row['side'], row['gene_id']]['models'])
                   if (row['side'], row['gene_id']) in displayed else '',
                   compressed_gaps_json=json.dumps(displayed[row['side'], row['gene_id']]['display_coordinates'].audit(extent))
                   if (row['side'], row['gene_id']) in displayed else '',
                   plot_label=labels[row['event_id']], gene_tree_branch_id=branch_id,
                   gene_tree_node=event.get('gene_tree_node', event.get('node_name')),
                   support_generax_ufboot=branches[branch_id]['support_generax_ufboot'])
    return audit


def render_context(path, rows, events, links, coordinates, *, gene_tree_panel=True, max_genes_per_side=None,
                   annotations=None):
    if max_genes_per_side is not None:
        if gene_tree_panel:
            raise ValueError("Bounded context pages require gene_tree_panel=False")
        return render_bounded_context(path, rows, events, links, coordinates, max_genes_per_side, annotations)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    selected, audit = choose_representatives(events, links, coordinates)
    # Every genomic track on this page uses exactly the same limits and scale.
    extent = max(
        [20.0]
        + [
            max(
                abs(a - (int(f["start"]) + int(f["end"])) / 2)
                for a in [int(f["start"]), int(f["end"])]
                + [x for pair in blocks(f.get("utr_blocks", "")) for x in pair]
            )
            / 1000
            + 20
            for _, _, _, f, near in selected
            if f
        ]
    )
    extent = math.ceil(extent / 5) * 5
    tracks_per_event = 3 if gene_tree_panel else 2
    fig = plt.figure(figsize=(15, max(7, (5 if gene_tree_panel else 4) * len(events) + 3)))
    fig.suptitle("A traceable candidate with donor and recipient context", fontsize=16, x=0.07, ha="left", y=0.98)
    fig.text(
        0.07,
        0.947 if gene_tree_panel else 0.93,
        "One passing gene per side per event; "
        + ("existing gene-tree paths and GFF coordinates." if gene_tree_panel else "existing GFF coordinates."),
        fontsize=10,
        va="baseline" if gene_tree_panel else "top",
    )
    grid = fig.add_gridspec(
        len(events) * tracks_per_event,
        1,
        left=0.09,
        right=0.96,
        top=0.89 if gene_tree_panel else 0.78,
        bottom=0.15 if gene_tree_panel else 0.22,
        hspace=1.12 if gene_tree_panel else 1.0,
        height_ratios=([1.5, 1, 1] if gene_tree_panel else [1, 1]) * len(events),
    )
    for i, event in enumerate(events):
        entries = selected[i * 2 : i * 2 + 2]
        if gene_tree_panel:
            ax = fig.add_subplot(grid[i * tracks_per_event])
            selected_gene_clade(ax, rows, event, entries[0][2]["gene_id"], entries[1][2]["gene_id"])
        for offset, (_, side, link, focal, neighbors) in enumerate(entries, int(gene_tree_panel)):
            ax = fig.add_subplot(grid[i * tracks_per_event + offset])
            ax.set_xlim(-extent, extent)
            ax.set_ylim(-0.65, 0.85)
            ax.set_yticks([])
            ax.axvline(0, color="#cccccc", lw=0.5, zorder=0)
            coverage = 100 * float(link["host_scaffold_background_class_classified_fraction"])
            compatible = 100 * float(link["host_scaffold_background_class_compatible_fraction"])
            title = (
                f"{event['event_id']} | {side}: {link.get('gene_species', '')} | {link['host_scaffold_id']}"
                f" | background coverage {coverage:.1f}%, compatible {compatible:.1f}%"
            )
            if not gene_tree_panel and side == "donor":
                branch = next(
                    r for r in rows if r["branch_id"] == event.get("gene_tree_branch_id", event.get("branch_id"))
                )
                title = f"HGT node {branch['node_name']} | UFB {support_label(branch)}\n" + title
            neighbor_key = []
            if focal is None:
                ax.text(0, 0, "GFF coordinates unavailable", ha="center", color="#777777")
            else:
                center = (int(focal["start"]) + int(focal["end"])) / 2
                neighbor_number = 0
                for j, row in enumerate(sorted(neighbors, key=lambda r: int(r["start"]))):
                    a, b = (int(row["start"]) - center) / 1000, (int(row["end"]) - center) / 1000
                    focal_flag = row["gene_id"] == link["gene_id"]
                    color = ORANGE if focal_flag else BLUE
                    info = structure(row)
                    if any(info[kind] for kind in ('coding', 'utr', 'exon')):
                        for x, y in info["introns"]:
                            ax.plot([(x - center) / 1000, (y - center) / 1000], [0, 0], color="#444444", lw=0.7)
                        for kind, height in [("coding", 0.20), ("exon", 0.20), ("utr", 0.10)]:
                            for x, y in info[kind]:
                                ax.add_patch(
                                    Rectangle(
                                        ((x - center) / 1000, -height / 2),
                                        (y - x + 1) / 1000,
                                        height,
                                        facecolor=color if kind == "coding" else "#f0f0f0" if kind == 'exon' else "#a6adb2",
                                        edgecolor=color,
                                        lw=0.6,
                                        hatch='..' if kind == 'exon' else None,
                                    )
                                )
                    else:
                        ax.add_patch(Rectangle((a, -0.1), b - a, 0.2, fill=False, ec=color, ls=":", lw=0.8))
                    direction = 1 if row.get("strand") == "+" else -1 if row.get("strand") == "-" else 0
                    if direction:
                        endpoint = b if direction == 1 else a
                        ax.annotate(
                            "",
                            xy=(endpoint, 0),
                            xytext=(endpoint - direction * min(0.4, max(0.05, b - a)), 0),
                            arrowprops=dict(arrowstyle="->", color=color, lw=0.7),
                        )
                    label = row["gene_id"].removeprefix(link.get("gene_species", "") + "_").replace("GeneID", "GID")
                    if not focal_flag:
                        neighbor_number += 1
                        neighbor_key.append(f"{neighbor_number}={label}")
                        label = str(neighbor_number)
                    if a < -extent or b > extent:
                        label += "*"
                    ax.text(
                        (max(a, -extent) + min(b, extent)) / 2,
                        0.34 if j % 2 == 0 else -0.34,
                        label,
                        ha="center",
                        va="center",
                        fontsize=7,
                        color=color,
                        weight="bold" if focal_flag else "normal",
                    )
                ax.text(
                    0.01,
                    0.98,
                    structure(focal)["status"].replace("_", " "),
                    transform=ax.transAxes,
                    fontsize=7,
                    color="#777777",
                )
            if neighbor_key:
                title += "\nNearby annotation IDs: " + ", ".join(neighbor_key)
            ax.set_title(title, loc="left", fontsize=8, pad=17)
            ax.spines[["top", "right", "left"]].set_visible(False)
            ax.tick_params(labelsize=8)
            ax.set_xlabel("Genomic position relative to recorded focal-feature midpoint (kb)", fontsize=8)
    fig.text(
        0.07,
        0.055,
        "Orange: focal gene; blue: nearby annotations; thick blocks: coding exons; thin gray blocks: UTR; lines: introns.\n"
        "UFB = Ultrafast bootstrap. All genomic tracks share one uncompressed kb axis."
        + (" Gene-tree paths use their own substitution/site axis.\n" if gene_tree_panel else "\n")
        + "The shared window includes each focal feature plus 20 kb flanks; * = neighboring feature extends beyond the display window.\n"
        "Neighbors are not asserted to be host-classified or conserved in order. CDS-only records do not establish complete exon/UTR structure.\n"
        "Representative selection and all event/gene identities are exported in the context audit; no sequence analysis was run.",
        fontsize=8,
        color="#666666",
    )
    fig.savefig(path, format="pdf")
    plt.close(fig)
    for row in audit:
        row["shared_axis_min_kb"] = -extent
        row["shared_axis_max_kb"] = extent
    return audit

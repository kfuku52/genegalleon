"""Combine native pairwise dS dotplots with shared species layouts and one key."""

import json
import math
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

from kffractbias.io import read_bed

try:
    from busco_reference_quality import safe_token
    from input_generation_array_state import FreshDigestBatch
    from pairwise_synteny_ds import anchor_pairs, read_ds
    from pairwise_synteny_style import MISSING_COLOR, STYLE, ds_colormap
    from synteny_multitrack import read_json
except ImportError:
    from .busco_reference_quality import safe_token
    from .input_generation_array_state import FreshDigestBatch
    from .pairwise_synteny_ds import anchor_pairs, read_ds
    from .pairwise_synteny_style import MISSING_COLOR, STYLE, ds_colormap
    from .synteny_multitrack import read_json

TRACK_FIELDS = {"species", "seqids", "orientation", "labels"}


def load_panels(path, pairs):
    """Select existing analyses; all input paths remain owned by the workflow."""
    payload = read_json(path)
    if (not isinstance(payload, dict) or set(payload) != {"panels"}
            or not isinstance(payload["panels"], list) or not payload["panels"]):
        raise ValueError("Dotplot panels file requires a nonempty panels list")
    analyses = {p["analysis_id"]: p for p in pairs}
    result, identifiers = [], set()
    for panel in payload["panels"]:
        if not isinstance(panel, dict) or set(panel) - {"panel_id", "comparisons", "tracks", "columns"}:
            raise ValueError("Unknown dotplot panel options")
        identifier = safe_token(panel.get("panel_id", ""), "panel_id")
        if identifier in identifiers:
            raise ValueError("Dotplot panel_id values must be unique")
        identifiers.add(identifier)
        comparisons, tracks = panel.get("comparisons"), panel.get("tracks")
        if (not isinstance(comparisons, list) or not comparisons
                or any(not isinstance(v, str) for v in comparisons)
                or len(set(comparisons)) != len(comparisons) or set(comparisons) - set(analyses)):
            raise ValueError("Dotplot comparisons must name unique analysis_id values in the pair table")
        if (not isinstance(tracks, list) or len(tracks) < 2
                or any(not isinstance(t, dict) or set(t) - TRACK_FIELDS for t in tracks)):
            raise ValueError("Dotplot panels require tracks with documented species/layout options")
        species = [t.get("species") for t in tracks]
        if (any(not isinstance(s, str) or not s for s in species)
                or len(set(species)) != len(species)):
            raise ValueError("Dotplot tracks require unique species")
        edges = [frozenset(analyses[c][s]["species"] for s in ("target", "query")) for c in comparisons]
        if set().union(*edges) != set(species) or len(set(edges)) != len(edges):
            raise ValueError("Dotplot comparisons must connect its tracks once per species pair")
        for track in tracks:
            seqids = track.get("seqids")
            if seqids is not None and (not isinstance(seqids, list) or not seqids
                    or any(not isinstance(s, str) or not s or any(c in s for c in "\n\r\0") for s in seqids)
                    or len(set(seqids)) != len(seqids)):
                raise ValueError("Dotplot seqids must be unique nonempty sequence identifiers")
            for field in ("orientation", "labels"):
                values = track.get(field, {})
                if (not isinstance(values, dict) or any(not isinstance(k, str) or not k for k in values)
                        or (seqids is not None and set(values) - set(seqids))):
                    raise ValueError(f"Dotplot {field} must map selected seqids")
            if any(v not in ("+", "-") for v in track.get("orientation", {}).values()):
                raise ValueError("Dotplot orientation values must be + or -")
            labels = track.get("labels", {}).values()
            if (any(not isinstance(v, str) or not v.strip() or any(c in v for c in "\n\r\0") for v in labels)
                    or len(set(labels)) != len(labels)):
                raise ValueError("Dotplot labels must be unique nonempty single-line text")
        columns = panel.get("columns", min(2, len(comparisons)))
        if isinstance(columns, bool) or not isinstance(columns, int) or not 1 <= columns <= 8:
            raise ValueError("Dotplot columns must be an integer from 1 to 8")
        result.append({**panel, "columns": columns})
    return result


def gene_ranks(genes, track):
    """Retain each chromosome's native BED ranks, with explicit display reversal."""
    by_sequence = defaultdict(list)
    for gene in genes:
        by_sequence[gene.seqid].append(gene)
    order = track["seqids"]
    if set(order) - set(by_sequence):
        raise ValueError("Dotplot track requests an unknown or unannotated sequence")
    if any(set(track.get(field, {})) - set(order) for field in ("orientation", "labels")):
        raise ValueError("Dotplot labels/orientation must map selected seqids")
    labels = {s: track.get("labels", {}).get(s, s) for s in order}
    if len(set(labels.values())) != len(labels):
        raise ValueError("Dotplot display labels must be unique within each species")
    bounds, ranks = [0], {}
    for seqid in order:
        rows = by_sequence[seqid]
        if track.get("orientation", {}).get(seqid, "+") == "-":
            rows = rows[::-1]
        for i, gene in enumerate(rows):
            ranks[gene.gene_id] = bounds[-1] + i + .5
        bounds.append(bounds[-1] + len(rows))
    return ranks, bounds, labels


def ordered_bed(path):
    """Validate with kfFractBias while preserving the native JCVI display ranks."""
    from jcvi.formats.bed import Bed

    genes = {gene.gene_id: gene for gene in read_bed(path)}
    return tuple(genes[row.accn] for row in Bed(str(path), sorted=False))


def _data(plan, panel, plots, ds, proof):
    analyses = {p["analysis_id"]: p for p in plan["pairs"]}
    beds, hashes, eligible, defaults, tables, source_hashes = {}, {}, {}, {}, {}, {}
    for comparison in panel["comparisons"]:
        pair = analyses[comparison]
        directory = plots / comparison
        paths = [directory / (side + ".bed") for side in ("target", "query")]
        paths += [directory / ("dotplot." + side + ".bed") for side in ("target", "query")]
        paths += [directory / "target.query.lifted.anchors", ds / comparison / "ds.tsv"]
        paths += [Path(pair["ds"][side]["fasta"]) for side in ("target", "query")]
        observed = proof.read(paths)
        source_hashes.update(observed)
        for side in ("target", "query"):
            species = pair[side]["species"]
            full, filtered = directory / (side + ".bed"), directory / ("dotplot." + side + ".bed")
            identity = (observed[str(full.resolve())], observed[str(Path(pair["ds"][side]["fasta"]).resolve())])
            genes = ordered_bed(full)
            selected = list(dict.fromkeys(g.seqid for g in ordered_bed(filtered)))
            if species in hashes and (identity != hashes[species] or set(selected) != eligible[species]):
                raise ValueError(f"Shared dotplot species has different BED, CDS or eligible sequences: {species}")
            beds.setdefault(species, genes)
            hashes[species], eligible[species] = identity, set(selected)
            defaults.setdefault(species, selected)
        pairs = anchor_pairs(directory / "target.query.lifted.anchors")
        _fields, rows, lookup = read_ds(ds / comparison / "ds.tsv")
        if set(lookup) != {a + "|" + b for a, b in pairs}:
            raise ValueError("Dotplot dS results must cover exactly every unique lifted anchor")
        all_genes = [{g.gene_id for g in beds[pair[s]["species"]]} for s in ("target", "query")]
        if any(genes[i] not in all_genes[i] for genes in pairs for i in (0, 1)):
            raise ValueError("Dotplot anchor has no source BED coordinate")
        tables[comparison] = pair, pairs, lookup, {r["pair_id"]: r["status"] for r in rows}, paths[:6]
    tracks, layouts = [], {}
    for source in panel["tracks"]:
        species = source["species"]
        track = {**source, "seqids": source.get("seqids", defaults[species])}
        if not track["seqids"] or set(track["seqids"]) - eligible[species]:
            raise ValueError("Dotplot track seqids must pass the pairwise physical-length filter")
        layouts[species] = gene_ranks(beds[species], track)
        tracks.append(track)
    species_order = [t["species"] for t in tracks]
    points, reports, paths = [], [], {}
    for comparison, (pair, pairs, lookup, statuses, inputs) in tables.items():
        source_species = [pair[s]["species"] for s in ("target", "query")]
        display_species = sorted(source_species, key=species_order.index)
        ranks = [layouts[s][0] for s in source_species]
        plotted = [(a, b) for a, b in pairs if a in ranks[0] and b in ranks[1]]
        if not plotted:
            raise ValueError(f"No anchors connect the selected dotplot tracks: {comparison}")
        swap = source_species != display_species
        values = [(ranks[0][a], ranks[1][b], lookup[a + "|" + b]) for a, b in plotted]
        if swap:
            values = [(y, x, value) for x, y, value in values]
        counts = Counter(statuses[a + "|" + b] for a, b in plotted)
        reports.append({"analysis_id": comparison, "source_species": source_species,
                        "display_species": display_species, "axes_swapped": swap,
                        "unique_anchor_pairs": len(pairs), "displayed_pairs": len(plotted),
                        "excluded_pairs": len(pairs) - len(plotted), "status_counts": dict(counts),
                        "missing_ds_pairs": sum(not math.isfinite(v) for _, _, v in values),
                        "above_color_limit": sum(v > plan["ds_color_max"] for _, _, v in values)})
        points.append(values)
        paths[comparison] = inputs
    proof.check()
    return tracks, layouts, points, reports, paths, source_hashes


def _spread(values, extent, width):
    """Place crowded labels with leader lines, retaining their exact boundaries."""
    gap = 8 / width
    if (len(values) - 1) * gap > .98:
        raise ValueError("Too many sequence labels at 8 pt; select fewer seqids or shorten labels")
    result = [max(v / extent, .01 + i * gap) for i, v in enumerate(values)]
    for i in range(1, len(result)):
        result[i] = max(result[i], result[i - 1] + gap)
    if result[-1] > .99:
        result[-1] = .99
        for i in range(len(result) - 2, -1, -1):
            result[i] = min(result[i], result[i + 1] - gap)
    return result


def _draw(plan, panel, tracks, layouts, points, reports, directory):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.backends.backend_pdf import RendererPdf
    from matplotlib.colors import Normalize
    from matplotlib.font_manager import FontProperties
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    columns = min(panel["columns"], len(points))
    rows = math.ceil(len(points) / columns)
    renderer = RendererPdf(None, 300, 8, 7.2)
    font = FontProperties(family="Helvetica", size=8)
    label_rows = {t["species"]: {s: layouts[t["species"]][2][s] +
                  (" (rev)" if t.get("orientation", {}).get(s) == "-" else "") for s in t["seqids"]} for t in tracks}
    label_width = max(renderer.get_text_width_height_descent(v, font, False)[0]
                      for labels in label_rows.values() for v in labels.values())
    cell_width, left, right, top, bottom = 3.6 * 72, 44, label_width + 9, label_width + 22, 33
    square = cell_width - left - right
    if square < 72:
        raise ValueError("Dotplot sequence labels leave insufficient room at 8 pt; use shorter labels")
    for species, (_ranks, bounds, _labels) in layouts.items():
        _spread([(a + b) / 2 for a, b in zip(bounds[:-1], bounds[1:], strict=True)], bounds[-1], square)
    cell_height = top + square + bottom
    empty_cell = len(points) < rows * columns
    legend_height = 0 if empty_cell else 62
    width, height = columns * cell_width, rows * cell_height + legend_height
    fig = plt.figure(figsize=(width / 72, height / 72))
    cmap, norm = ds_colormap(), Normalize(0, plan["ds_color_max"], clip=True)
    try:
        for i, (values, report) in enumerate(zip(points, reports, strict=True)):
            col, row = i % columns, i // columns
            x = col * cell_width + left
            y = height - (row + 1) * cell_height + bottom
            ax = fig.add_axes((x / width, y / height, square / width, square / height))
            data = np.asarray(values, dtype=float)
            scatter = ax.scatter(data[:, 0], data[:, 1], c=data[:, 2],
                                 cmap=cmap, norm=norm, s=2, edgecolors="none", linewidths=0, plotnonfinite=True)
            if np.ma.count(scatter.get_offsets()) != 2 * len(values):
                raise ValueError("The combined renderer discarded an anchor")
            for dim, species in enumerate(report["display_species"]):
                _ranks, bounds, _labels = layouts[species]
                extent = bounds[-1]
                mids = [(a + b) / 2 for a, b in zip(bounds[:-1], bounds[1:], strict=True)]
                labels = label_rows[species]
                positions = _spread(mids, extent, square)
                for seqid, mid, pos in zip(labels, mids, positions, strict=True):
                    xy, text_xy = ((mid / extent, 1), (pos, 1 + 5 / square)) if dim == 0 else ((1, 1 - mid / extent), (1 + 5 / square, 1 - pos))
                    ax.annotate(labels[seqid], xy=xy, xycoords='axes fraction', xytext=text_xy,
                                textcoords='axes fraction', ha='center' if dim == 0 else 'left',
                                va='bottom' if dim == 0 else 'center', rotation=90 if dim == 0 else 0,
                                annotation_clip=False, arrowprops={'arrowstyle': '-', 'lw': .35,
                                                                 'color': '#777777', 'shrinkA': 1, 'shrinkB': 1})
                for boundary in bounds[1:-1]:
                    (ax.axvline if dim == 0 else ax.axhline)(boundary, color='#dddddd', lw=.4, zorder=0)
                axis = ax.xaxis if dim == 0 else ax.yaxis
                axis.set_major_locator(MaxNLocator(nbins=3, integer=True))
                axis.set_major_formatter(FuncFormatter(lambda v, _pos: f'{v / 1000:g}k' if v >= 1000 else f'{v:g}'))
            x_species, y_species = report["display_species"]
            ax.set_xlim(0, layouts[x_species][1][-1])
            ax.set_ylim(layouts[y_species][1][-1], 0)
            ax.set_xlabel(x_species.replace('_', ' '), fontstyle='italic', labelpad=4)
            ax.set_ylabel(y_species.replace('_', ' '), fontstyle='italic', labelpad=4)
            ax.tick_params(direction='out', length=2.5, width=.5, pad=2)
            letter = chr(97 + i) if i < 26 else str(i + 1)
            report['panel_label'] = letter
            fig.text((col * cell_width + 6) / width, (y + square + top - 12) / height,
                     letter, fontweight='bold')
            fig.text((x + square / 2) / width, (y + square + top - 12) / height,
                     f"{report['displayed_pairs']:,} pairs; NA: {report['missing_ds_pairs']:,}", ha='center')
        if empty_cell:
            lx = (columns - 1) * cell_width + 36
            ly = cell_height / 2
        else:
            lx, ly = max(12, (width - 180) / 2), 32
        cax = fig.add_axes((lx / width, ly / height, 180 / width, 6 / height))
        bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                          orientation='horizontal', extend='max' if any(r['above_color_limit'] for r in reports) else 'neither')
        bar.outline.set_linewidth(.5)
        bar.ax.tick_params(length=2, width=.5, pad=2)
        fig.text((lx + 90) / width, (ly + 17) / height, 'dS', fontstyle='italic', ha='center')
        key_y = ly - 25 if empty_cell else ly - 26
        fig.patches.append(Rectangle((lx / width, key_y / height), 6 / width, 6 / height,
                                     transform=fig.transFigure, facecolor=MISSING_COLOR, edgecolor='none'))
        fig.text((lx + 10) / width, (key_y + 3) / height, 'Unestimable or saturated', va='center')
        if empty_cell:
            fig.text(lx / width, (key_y - 20) / height, 'Axes: cumulative gene rank')
            fig.text(lx / width, (key_y - 33) / height, 'Grid lines: sequence boundaries')
            if any('-' in t.get('orientation', {}).values() for t in tracks):
                fig.text(lx / width, (key_y - 46) / height, '(rev): display reversal')
        original_dpi = fig.dpi
        fig.set_dpi(72)
        try:
            box = fig.get_tightbbox(renderer)
            if box.x0 < 0 or box.y0 < 0 or box.x1 > width / 72 or box.y1 > height / 72:
                raise ValueError('Dotplot labels do not fit at 8 pt; use shorter species/sequence labels')
        finally:
            fig.set_dpi(original_dpi)
        for fmt in plan['formats']:
            fig.savefig(directory / ('dotplot.' + fmt), format=fmt, dpi=300, facecolor='white')
    finally:
        plt.close(fig)
    return width, height


def render_workflow_panels(plan, plots, ds):
    """Publish new, self-contained bundles after every source fence passes."""
    if not plan.get('dotplot_panels'):
        return
    if plan.get('dotplot_color') != 'ds' or not math.isfinite(plan['ds_color_max']) or plan['ds_color_max'] <= 0:
        raise ValueError('Combined dotplot panels require dS coloring and a positive finite color limit')
    import matplotlib
    matplotlib.use('Agg')

    for panel in plan['dotplot_panels']:
        destination = plots / 'dotplot_panels' / panel['panel_id']
        if destination.exists():
            raise ValueError('Combined dotplot output directory must be new')
        proof = FreshDigestBatch()
        tracks, layouts, points, reports, paths, source_hashes = _data(plan, panel, plots, ds, proof)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='.gg-dotplots-', dir=destination.parent) as temporary:
            staged = Path(temporary) / 'bundle'
            staged.mkdir()
            copied = FreshDigestBatch()
            for comparison, inputs in paths.items():
                folder = staged / 'inputs' / comparison
                folder.mkdir(parents=True)
                for source in inputs:
                    target = folder / source.name
                    shutil.copyfile(source, target)
                    if copied.read([target])[str(target.resolve())] != source_hashes[str(source.resolve())]:
                        raise OSError('Combined dotplot input copy changed: ' + str(source))
            with matplotlib.rc_context(STYLE):
                width, height = _draw(plan, panel, tracks, layouts, points, reports, staged)
            proof.check()
            report = {'schema_version': 1, 'coordinate_system': 'gene_rank', 'gene_rank_centres': 'zero-based rank + 0.5',
                      'species_order': [t['species'] for t in tracks], 'tracks': tracks, 'comparisons': reports,
                      'shared_ds_scale': [0, plan['ds_color_max']], 'missing_color': MISSING_COLOR,
                      'minimum_length_bp': plan.get('dotplot_min_length', 1_000_000),
                      'downsampled': False, 'filtered_by_ds': False, 'figure_size_pt': [width, height],
                      'interpretation': 'Interspecies gene-pair divergence; not a WGD date or test',
                      'source_hashes': source_hashes}
            (staged / 'dotplot_panels.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n', encoding='utf-8')
            (staged / 'manifest.json').write_text(json.dumps(panel, indent=2, sort_keys=True) + '\n', encoding='utf-8')
            proof.check()
            copied.check()
            if destination.exists():
                raise ValueError('Combined dotplot destination appeared during rendering')
            os.rename(staged, destination)

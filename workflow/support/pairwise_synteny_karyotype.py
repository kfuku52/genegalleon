"""Compact JCVI ribbons with explicit colors and shared or independent scales."""

import argparse
import hashlib
import json
import math
import os
import re
import tempfile
from contextlib import ExitStack
from pathlib import Path

try:
    from pairwise_synteny_dotplot import homoeolog_order
    from pairwise_synteny_style import SOFT_COLORS, STYLE, style_text
except ImportError:
    from .pairwise_synteny_dotplot import homoeolog_order
    from .pairwise_synteny_style import SOFT_COLORS, STYLE, style_text


def chromosome_colors(selected, genomes, anchors, mode="chromosome"):
    from kffractbias.io import natural_key

    if mode not in {"chromosome", "homoeolog"}:
        raise ValueError("karyotype-color must be chromosome or homoeolog")
    groups = []
    if mode == "homoeolog":
        _, metadata = homoeolog_order(selected, genomes, anchors)
        groups = metadata["groups"]
    maps = [{}, {}]
    for index, group in enumerate(groups):
        for side, track in enumerate(("target", "query")):
            maps[side].update((seqid, SOFT_COLORS[index % len(SOFT_COLORS)]) for seqid in group[track])
    for side, genome in enumerate(genomes):
        next_color = len(groups)
        for seqid in sorted({gene.seqid for gene in genome}, key=natural_key):
            if seqid not in maps[side]:
                maps[side][seqid] = SOFT_COLORS[next_color % len(SOFT_COLORS)]
                next_color += 1
    return {"mode": mode, "palette": list(SOFT_COLORS), "groups": groups,
            "chromosomes": dict(zip(("target", "query"), maps, strict=True)),
            "interpretation": "display colors only; homoeolog groups are not a biological test"}


def gene_scale(totals):
    """An integral 1/2/5 scale small enough for both tracks."""
    limit = max(1, min(totals) / 6)
    power = 10 ** math.floor(math.log10(limit))
    return int(max(x * power for x in (1, 2, 5) if x * power <= limit))


def scale_tracks(tracks, mode="shared"):
    """Keep native gaps and rank mapping; left-align with a common gene width."""
    if mode not in {"shared", "independent"}:
        raise ValueError("karyotype-scale must be shared or independent")
    if mode == "shared":
        ratio = min(track.ratio for track in tracks)
        for track in tracks:
            track.ratio = ratio
            track.xend = track.xstart + ratio * track.total + track.gap * (len(track.seqids) - 1)
            track.update_offsets()


def connection_criteria(analysis):
    """Use recorded analysis settings, never current/default plotting settings."""
    analysis = Path(analysis)
    summary_file, log_file = analysis / "summary.json", analysis / "logs/01.mcscan.log"
    parameters = json.loads(summary_file.read_text(encoding="utf-8"))["parameters"]
    log = log_file.read_text(encoding="utf-8")

    def recorded(pattern, name, convert):
        values = {convert(value) for value in re.findall(pattern, log)}
        if len(values) != 1:
            raise ValueError(f"Cannot annotate connections: missing or ambiguous recorded {name}")
        return values.pop()

    evalue = recorded(r"--evalue\s+([0-9.eE+-]+)", "DIAMOND E-value", float)
    tandem = recorded(r"tandem_Nmax=(\d+)", "tandem distance", int)
    lifted = recorded(r"new pairs found \(dist=(\d+)\)", "liftover distance", int)
    cscore, minimum, distance = parameters["cscore"], parameters["min_anchors"], parameters["distance"]
    if not math.isfinite(evalue) or evalue <= 0 or not 0 < cscore <= 1 or minimum < 2 or distance < 1:
        raise ValueError("Cannot annotate connections: invalid recorded thresholds")
    return {"protein_search": "DIAMOND blastp", "protein_evalue": evalue,
            "seed_cscore": cscore, "seed_min_unique_genes_per_genome": minimum,
            "chaining_max_gene_rank_gap_per_genome": distance, "tandem_gene_rank_distance": tandem,
            "liftover_gene_rank_distance": lifted, "liftover_metric": "Manhattan", "liftover_bound": "strict",
            "quota": parameters["quota"], "ds_filter": False,
            "source_hashes": {str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in (summary_file, log_file)}}


def connection_legend_text(criteria):
    evalue = f"{criteria['protein_evalue']:g}".replace("e-0", "e-").replace("e+0", "e+")
    quota = "no quota" if criteria["quota"] is None else f"quota = {criteria['quota']}"
    return ("Connections: syntenic blocks (MCscan + liftover)\n"
            f"Protein hits: E <= {evalue}; seed C-score >= {criteria['seed_cscore']:g}; "
            f"tandem distance = {criteria['tandem_gene_rank_distance']} gene ranks\n"
            f"Seed: >= {criteria['seed_min_unique_genes_per_genome']} unique genes / genome; "
            f"chaining gap <= {criteria['chaining_max_gene_rank_gap_per_genome']} gene ranks / genome\n"
            f"Liftover: |dx| + |dy| < {criteria['liftover_gene_rank_distance']} gene ranks; {quota}; no dS filter")


def resolved_layout(directory):
    """Resolve JCVI's BED/block paths beside the saved layout, without chdir."""
    source = directory / "layout"
    inputs, lines = [("layout", source), ("seqids", directory / "seqids")], []
    for line in source.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            lines.append(line)
            continue
        fields = [field.strip() for field in line.split(",")]
        edge = fields[0] == "e"
        if (edge and len(fields) not in {4, 5}) or (not edge and len(fields) not in {8, 9}):
            raise ValueError("Invalid pairwise karyotype layout row.")
        index = 3 if edge else 7
        if not fields[index]:
            raise ValueError("Karyotype layout requires a BED or block file path.")
        path = Path(fields[index])
        path = path if path.is_absolute() else directory / path
        fields[index] = str(path.resolve())
        inputs.append(("blocks" if edge else "BED", path))
        lines.append(",".join(fields))
    return "\n".join(lines) + "\n", inputs


def render_karyotype(directory, pair, fmt, colors, scale_mode="shared", analysis=None, track_order="target-query",
                     *, layout_report=None):
    directory = Path(directory).resolve()
    if fmt not in {"pdf", "svg", "png"}:
        raise ValueError("Karyotype format must be pdf, svg or png")
    if scale_mode not in {"shared", "independent"}:
        raise ValueError("karyotype-scale must be shared or independent")
    if track_order not in {"target-query", "query-target"}:
        raise ValueError("karyotype-track-order must be target-query or query-target")
    sides = ("target", "query")
    display_sides = sides if track_order == "target-query" else sides[::-1]
    os.environ.setdefault("MPLCONFIGDIR", str(directory / ".mplconfig"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from jcvi.graphics.chromosome import HorizontalChromosome
    from jcvi.graphics.karyotype import Layout, ShadeManager, Track
    from matplotlib.backends.backend_pdf import RendererPdf
    from matplotlib.offsetbox import AnnotationBbox, DrawingArea, HPacker, TextArea
    from matplotlib.patches import PathPatch, Rectangle
    from matplotlib.path import Path as MplPath
    from matplotlib.transforms import Bbox
    from nwkit.file_paths import validate_distinct_output_paths, validate_outputs_do_not_replace_inputs
    from nwkit.output_transaction import output_transaction

    layout_text, inputs = resolved_layout(directory)
    colors_path = directory / "karyotype_colors.json"
    if colors_path.exists():
        inputs.append(("colors", colors_path))
    if analysis is not None:
        inputs.extend(("analysis metadata", Path(analysis) / path) for path in ("summary.json", "logs/01.mcscan.log"))
    outfile = directory / f"karyotype.{fmt}"
    outputs = [("plot", outfile)] + ([("layout_report", layout_report)] if layout_report is not None else [])
    validate_distinct_output_paths(outputs)
    validate_outputs_do_not_replace_inputs(inputs, outputs)
    criteria = connection_criteria(analysis) if analysis is not None else None
    legend_text = connection_legend_text(criteria) if criteria is not None else None
    with matplotlib.rc_context(STYLE), ExitStack() as cleanup:
        fig = plt.figure(figsize=(9, 3.6), dpi=72)
        cleanup.callback(plt.close, fig)
        root = fig.add_axes((0, 0, 1, 1))
        root.set_xlim(0, 1)
        root.set_ylim(0, 1)
        root.set_axis_off()
        with tempfile.TemporaryDirectory(prefix="gg-karyotype-") as temporary:
            absolute_layout = Path(temporary) / "layout"
            absolute_layout.write_text(layout_text, encoding="utf-8")
            layout = Layout(str(absolute_layout), generank=True, seed=1)
        selected = (directory / "seqids").read_text(encoding="utf-8").splitlines()
        if len(layout) != 2 or len(selected) != 2 or any(entry.empty for entry in layout):
            raise ValueError("Pairwise karyotype requires exactly two tracks")
        for entry in layout:
            if (entry.rotation != 0 or not all(math.isfinite(value) for value in (entry.y, entry.xstart, entry.xend))
                    or not 0 < entry.y < 1 or not 0 <= entry.xstart < entry.xend <= 1):
                raise ValueError("Compact karyotype requires horizontal tracks with valid normalized coordinates.")
        upper, lower = sorted((layout[0].y, layout[1].y), reverse=True)
        if upper == lower:
            raise ValueError("Pairwise karyotype needs distinct vertical track positions")
        # Keep BEDs, seqids and edge indices in analysis order. Moving only
        # the vertical coordinates preserves every ribbon endpoint.
        layout[0].y, layout[1].y = (upper, lower) if track_order == "target-query" else (lower, upper)
        # Use native JCVI ranks/geometry and ribbon primitives; never infer bp
        # length from the last gene or label a gene-rank scale as Mbp.
        for entry, line in zip(layout, selected, strict=True):
            entry.seqids = [sid.strip() for sid in line.split(",")]
            if not all(entry.seqids) or len(entry.seqids) != len(set(entry.seqids)):
                raise ValueError("Each karyotype track requires unique nonempty chromosome IDs.")
            identifiers = [gene.accn for gene in entry.bed]
            if len(identifiers) != len(set(identifiers)):
                raise ValueError("Karyotype BED gene IDs must be unique within each track.")
            entry.rev = set()
            entry.sizes = {sid: len(list(entry.bed.sub_bed(sid))) for sid in entry.seqids}
            if not all(entry.sizes.values()):
                raise ValueError("Every selected karyotype chromosome must have genes in its BED.")
        for first, second, blocks, _ in layout.edges:
            if {first, second} != {0, 1}:
                raise ValueError("Pairwise ribbon edges must connect the two different tracks.")
            for a, b, c, d, score, orientation, _ in blocks:
                if score <= 0 or orientation not in {"+", "-"}:
                    raise ValueError("Invalid karyotype block score or orientation.")
                for side, endpoints in ((first, (a, b)), (second, (c, d))):
                    coordinates = layout[side].order_in_chr
                    if any(gene not in coordinates for gene in endpoints):
                        raise ValueError("Karyotype block endpoints must occur in their track BED.")
                    if coordinates[endpoints[0]][0] != coordinates[endpoints[1]][0]:
                        raise ValueError("A karyotype ribbon block cannot span multiple chromosomes.")
        tracks = [Track(root, entry, draw=False) for entry in layout]
        scale_tracks(tracks, scale_mode)
        if any(not math.isfinite(track.ratio) or track.ratio <= 0 for track in tracks):
            raise ValueError("Karyotype chromosome gaps leave no positive space for genes; select fewer chromosomes.")
        ShadeManager(root, tracks, layout)
        scale = gene_scale([track.total for track in tracks])
        scale_lines, species_labels, unit_labels, chromosome_labels = [], [], [], []
        for side, track in zip(sides, tracks, strict=True):
            direction = 1 if side == display_sides[0] else -1
            track_labels = []
            for sid in track.seqids:
                start, length = track.offsets[sid], track.ratio * track.sizes[sid]
                chromosome = HorizontalChromosome(root, start, start + length, track.y, height=0.012,
                                                  lw=0.5, ec="#707070", fc=colors["chromosomes"][side][sid],
                                                  style="rect")
                chromosome.set_transform(track.tr)
                label = root.annotate(sid, (start + length / 2, track.y), xytext=(0, direction * 5),
                                      textcoords="offset points", rotation=90, rotation_mode="anchor",
                                      ha="left" if direction == 1 else "right", va="center")
                track_labels.append(label)
            chromosome_labels.append(track_labels)
            species = root.annotate(pair[f"{side}_species"].replace("_", " "),
                                    ((track.xstart + track.xend) / 2, track.y), xytext=(0, 0),
                                    textcoords="offset points", ha="center",
                                    va="bottom" if direction == 1 else "top")
            units = root.annotate(" (gene rank)", (1, 0), xycoords=species, xytext=(0, 0),
                                  textcoords="offset points", ha="left", va="bottom")
            species_labels.append(species)
            unit_labels.append(units)
            # One bottom bar is valid for both shared-scale tracks. Independent
            # normalization needs a separate, correctly sized bar per track.
            if scale_mode == "independent" or side == display_sides[-1]:
                line, = root.plot([track.xstart, track.xstart + track.ratio * scale], [track.y, track.y],
                                  color="black", linewidth=0.6, solid_capstyle="butt", clip_on=False)
                text = root.annotate(f"{scale:,} {'gene' if scale == 1 else 'genes'}", (track.xstart + track.ratio * scale / 2, track.y),
                                     xytext=(0, 0), textcoords="offset points", ha="center",
                                     va="bottom", annotation_clip=False)
                scale_lines.append((track, direction, line, text))
        legend = None
        if legend_text is not None:
            symbol = DrawingArea(26, 26, 0, 0)
            vertices = [(2, 23), (10, 23), (10, 15), (24, 11), (24, 3),
                        (16, 3), (16, 11), (2, 15), (2, 23), (2, 23)]
            path = MplPath(vertices, [MplPath.MOVETO, MplPath.LINETO,
                                     *([MplPath.CURVE4] * 3), MplPath.LINETO,
                                     *([MplPath.CURVE4] * 3), MplPath.CLOSEPOLY])
            symbol.add_artist(PathPatch(path, facecolor=SOFT_COLORS[0], edgecolor="none", alpha=0.6))
            for x, y in ((0, 23), (14, 1)):
                symbol.add_artist(Rectangle((x, y), 12, 2, facecolor=SOFT_COLORS[0], edgecolor="none"))
            label = TextArea(legend_text, textprops={"fontfamily": "Helvetica", "fontsize": 8, "color": "black"})
            legend = AnnotationBbox(HPacker(children=[symbol, label], align="center", pad=0, sep=8),
                                    (0, 0), xycoords=root.transAxes, box_alignment=(0, 1),
                                    frameon=False, pad=0, annotation_clip=False)
            root.add_artist(legend)
        style_text(fig)
        for species in species_labels:
            species.set_fontstyle("italic")
        renderer = RendererPdf(None, 300, 3.6, 9)
        props = unit_labels[0].get_fontproperties()
        space = (renderer.get_text_width_height_descent("x x", props, False)[0]
                 - renderer.get_text_width_height_descent("xx", props, False)[0])
        suffix = space + renderer.get_text_width_height_descent("(gene rank)", props, False)[0]
        label_gap = 4
        label_layout = []

        def content(width):
            fig.set_size_inches(width, width * 0.4)
            label_layout.clear()
            for index, track in enumerate(tracks):
                direction = 1 if sides[index] == display_sides[0] else -1
                labels = chromosome_labels[index]
                boxes = [label.get_window_extent(renderer) for label in labels]
                anchor_y = root.transData.transform((0, track.y))[1]
                edge = max(box.y1 for box in boxes) if direction == 1 else min(box.y0 for box in boxes)
                extent = direction * (edge - anchor_y)
                species_labels[index].set_position((-suffix / 2, direction * (extent + label_gap)))
                unit_labels[index].set_position((space, 0))
                longest = max(labels, key=lambda label: renderer.get_text_width_height_descent(
                    label.get_text(), label.get_fontproperties(), False)[0])
                label_layout.append({"side": "target" if index == 0 else "query",
                                     "longest_chromosome_label": longest.get_text(),
                                     "chromosome_label_extent_pt": extent,
                                     "species_offset_pt": direction * (extent + label_gap),
                                     "species_label_gap_pt": label_gap})
            for track, direction, line, text in scale_lines:
                index = tracks.index(track)
                species_box = Bbox.union([species_labels[index].get_window_extent(renderer),
                                          unit_labels[index].get_window_extent(renderer)])
                edge = species_box.y1 if direction == 1 else species_box.y0
                # Keep the bottom scale label above its bar and below the
                # species label, with room for the label's actual font height.
                clearance = 12 if direction == 1 else max(12, text.get_window_extent(renderer).height + 5 + label_gap)
                y = (edge + direction * clearance) / fig.bbox.height
                line.set_ydata((y, y))
                text.xy = (text.xy[0], y)
                text.set_position((0, 5))
            artists = [*root.patches, *root.lines, *root.texts]
            box = Bbox.union([artist.get_window_extent(renderer) for artist in artists if artist.get_visible()])
            if legend is not None:
                legend.xy = (box.x0 / fig.bbox.width, (box.y0 - 12) / fig.bbox.height)
                legend.xybox = legend.xy
                box = Bbox.union([box, legend.get_window_extent(renderer)])
            return box

        width, padding = 7.2, 0.06
        low, high = 1.0, 12.0
        if content(low).width > (width - 2 * padding) * 72:
            raise ValueError("Karyotype labels cannot fit a 7.2-inch page at 8 pt")
        for _ in range(32):
            middle = (low + high) / 2
            if content(middle).width <= (width - 2 * padding) * 72:
                low = middle
            else:
                high = middle
        box = content(low).transformed(fig.dpi_scale_trans.inverted())
        page = Bbox.from_bounds((box.x0 + box.x1 - width) / 2,
                                box.y0 - padding, width, box.height + 2 * padding)
        report = {"coordinate_system": "gene_rank", "scale_mode": scale_mode, "scale_unit": "genes", "scale_value": scale,
            "track_order": list(display_sides),
            "display_species_order": [pair[f"{side}_species"] for side in display_sides],
            "track_metadata_order": list(sides),
            "scale_bar_count": len(scale_lines), "track_gene_counts": [track.total for track in tracks],
            "scale_label_position": "above",
            "track_extents": [[track.xstart, track.xend] for track in tracks], "track_alignment": "left",
            "track_ratios": [track.ratio for track in tracks], "pdf_width_inches": 7.2,
            "font": "Helvetica", "font_size_pt": 8, "text_color": "black",
            "species_font_style": "italic", "chromosome_label_rotation": 90,
            "species_label_layout": label_layout,
            "chromosome_style": "rectangular",
            "connection_criteria": criteria, "connection_legend_text": legend_text}
        with output_transaction([Path(path) for _, path in outputs]) as staged:
            fig.savefig(staged[outfile], format=fmt, dpi=300, bbox_inches=page)
            if layout_report is not None:
                Path(staged[Path(layout_report)]).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                                                           encoding="utf-8")
        return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", required=True, type=Path)
    parser.add_argument("--target-species", required=True)
    parser.add_argument("--query-species", required=True)
    parser.add_argument("--format", choices=("pdf", "png", "svg"), required=True)
    parser.add_argument("--scale", choices=("shared", "independent"), default="shared")
    parser.add_argument("--track-order", choices=("target-query", "query-target"), default="target-query",
                        help="Top-to-bottom display order; analysis inputs and ribbon endpoints are unchanged.")
    parser.add_argument("--analysis", type=Path, help="Recorded analysis directory for the connection legend")
    args = parser.parse_args()
    colors = json.loads((args.directory / "karyotype_colors.json").read_text(encoding="utf-8"))
    render_karyotype(args.directory, vars(args), args.format, colors, args.scale, args.analysis, args.track_order,
                     layout_report=args.directory / "karyotype_style.json")


if __name__ == "__main__":
    main()

"""Physical-length selection and reproducible chromosome ordering for dotplots."""

import csv
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path

try:
    from fasta_sequence_store import open_text
    from pairwise_synteny_style import SOFT_COLORS, STYLE, style_text
except ImportError:
    from .fasta_sequence_store import open_text
    from .pairwise_synteny_style import SOFT_COLORS, STYLE, style_text


def save_pdf_dotplot(fig, root, ax, pair, output, colorbar_ax=None, fmt="pdf", missing_legend=None,
                     rich_title=None, rich_colorbar_label=None):
    """3.6-inch PDF page; square plot box and unscaled 8-point typography."""
    import matplotlib
    from matplotlib.backends.backend_pdf import RendererPdf
    from matplotlib.text import Text
    from matplotlib.transforms import Bbox

    with matplotlib.rc_context(STYLE):
        # JCVI's chromosome labels use the fixed 0.1..0.9 figure coordinates.
        # A square canvas keeps them aligned with this square 0.8 x 0.8 box.
        fig.set_size_inches(6, 6)
        ax.set_box_aspect(1)
        # JCVI hides numeric tick marks. Restore them on the bottom/left,
        # pointing away from the plot even when the y axis is inverted.
        ax.tick_params(axis="both", which="major", direction="out", length=3, width=0.8,
                       colors="black", bottom=True, left=True, top=False, right=False,
                       labelbottom=True, labelleft=True, labeltop=False, labelright=False)
        for text in root.texts:
            if text.get_rotation() == 45:
                text.set_rotation(90)
        ax.set_xlabel(pair["target_species"].replace("_", " "))
        ax.set_ylabel(pair["query_species"].replace("_", " "))
        style_text(fig)
        ax.xaxis.label.set_fontstyle("italic")
        ax.yaxis.label.set_fontstyle("italic")
        for label in (rich_title, rich_colorbar_label):
            if label is not None:
                for text in label.findobj(Text):
                    if text.get_text() == "dS":
                        text.set_fontstyle("italic")
        # Separate plain-text artists retain real Helvetica/Helvetica-Oblique;
        # mathtext would substitute other fonts. Artist-relative annotations
        # follow the labels even when savefig translates the cropped page.
        units_x = ax.annotate(" (gene rank)", (1, 0), xycoords=ax.xaxis.label,
                              xytext=(0, 0), textcoords="offset points",
                              ha="left", va="bottom", fontstyle="normal", annotation_clip=False)
        units_y = ax.annotate(" (gene rank)", (0, 1), xycoords=ax.yaxis.label,
                              xytext=(0, 0), textcoords="offset points", rotation=90,
                              ha="left", va="bottom", fontstyle="normal", annotation_clip=False)
        # Agg uses substituted fonts; measure with the PDF backend's Helvetica
        # metrics at 72 dpi instead. Crop the page, never scale the PDF/text.
        original_dpi = fig.dpi
        fig.set_dpi(72)
        renderer = RendererPdf(None, 300, 6, 6)
        # AFM reports zero ink width for an isolated space. Measure its advance
        # between glyphs; the actual space is included in each upright suffix.
        props = units_x.get_fontproperties()
        spacing = (renderer.get_text_width_height_descent("x x", props, False)[0]
                   - renderer.get_text_width_height_descent("xx", props, False)[0])
        for label in (rich_title, rich_colorbar_label):
            if label is not None:
                label.offsetbox.sep = spacing
        suffix_width = spacing + renderer.get_text_width_height_descent(
            "(gene rank)", units_x.get_fontproperties(), False)[0]
        units_x.set_position((spacing, 0))
        units_y.set_position((0, spacing))

        def layout(size):
            fig.set_size_inches(size, size)
            ax.get_tightbbox(renderer)  # Update axis-label positions/aspect.
            # Center each compound label, not just its italic species portion.
            ax.xaxis.label.set_x(0.5 - suffix_width / (2 * ax.bbox.width))
            ax.yaxis.label.set_y(0.5 - suffix_width / (2 * ax.bbox.height))
            top = max([ax.bbox.y1] + [text.get_window_extent(renderer).y1
                      for text in root.texts if text.get_rotation() == 90])
            legend = root.get_legend()
            if legend is not None:
                legend.set_bbox_to_anchor((0.09, 1))
                shift = (top + 6 - legend.get_window_extent(renderer).y0) / (size * 72)
                legend.set_bbox_to_anchor((0.09, 1 + shift))
                top = legend.get_window_extent(renderer).y1
            title_y = (top + 6) / (size * 72)
            if rich_title is not None:
                rich_title.xy = (0.5, title_y)
            if root.title.get_text():
                root.set_title(root.title.get_text(), y=title_y)
            for text in root.texts:
                if text.get_text().startswith("Pairwise synteny"):
                    text.set_y(title_y)
                    text.set_verticalalignment("bottom")
            if colorbar_ax is not None:
                bottom = min(ax.xaxis.label.get_window_extent(renderer).y0,
                             units_x.get_window_extent(renderer).y0)
                legend_width = missing_legend.get_window_extent(renderer).width if missing_legend is not None else 0
                available = ax.bbox.width - legend_width - 8
                bar_width = max(36, available) if missing_legend is not None else ax.bbox.width * 0.5
                # The gray square and the quantitative bar share one centerline.
                colorbar_ax.set_position((ax.bbox.x0 / (size * 72), (bottom - 23) / (size * 72),
                                          bar_width / (size * 72), 6 / (size * 72)))
                # The compound label follows the native empty axis label's
                # updated anchor in every output backend, including PNG at 300 dpi.
                colorbar_ax.get_tightbbox(renderer)
            return fig.get_tightbbox(renderer)

        try:
            width, padding = 3.6, 0.1
            low, high = 0.5, width
            if layout(low).width > width - 2 * padding:
                raise ValueError("Dotplot labels cannot fit a 3.6-inch PDF at 8 pt")
            # Find the largest square plot canvas fitting labels plus margins.
            for _ in range(32):
                size = (low + high) / 2
                if layout(size).width <= width - 2 * padding:
                    low = size
                else:
                    high = size
            box = layout(low)
            page = Bbox.from_bounds((box.x0 + box.x1 - width) / 2,
                                    box.y0 - padding, width, box.height + 2 * padding)
            fig.savefig(output, format=fmt, dpi=300, bbox_inches=page)
        finally:
            fig.set_dpi(original_dpi)


def render_orientation_pdf(directory, pair, filtered=True, fmt="pdf"):
    """Keep JCVI's orientation palette/anchors while applying PDF typography."""
    import os

    os.environ.setdefault("MPLCONFIGDIR", str(directory / ".mplconfig"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from jcvi.formats.bed import Bed
    from jcvi.graphics.dotplot import Palette, dotplot

    prefix = "dotplot." if filtered else ""
    beds = [Bed(str(directory / f"{prefix}{side}.bed"), sorted=False) for side in ("target", "query")]
    orders = [bed.order for bed in beds]
    anchors = directory / ("dotplot.anchors" if filtered else "target.query.lifted.anchors")
    count = 0
    with anchors.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2 or any(fields[i] not in orders[i] for i in (0, 1)):
                raise ValueError("Invalid dotplot anchor")
            count += 1
    if not count:
        raise ValueError("No anchors to render")
    palette = Palette.from_block_orientation(str(anchors), *beds,
                                             forward_color=SOFT_COLORS[0], reverse_color=SOFT_COLORS[3])
    fig = plt.figure(figsize=(6, 6))
    root = fig.add_axes((0, 0, 1, 1), frameon=False)
    ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))
    try:
        dotplot(str(anchors), *beds, fig, root, ax, palette=palette, sample_number=count,
                minfont=4, chpf=False, usetex=False, sepcolor="#b0b0b0",
                title="Pairwise synteny (gene rank)")
        if len(ax.collections[0].get_offsets()) != count:
            raise ValueError("The renderer discarded an anchor")
        save_pdf_dotplot(fig, root, ax, pair, directory / f"dotplot.{fmt}", fmt=fmt)
    finally:
        plt.close(fig)


def chromosome_lengths(source, genes, minimum):
    """Use explicit assembly lengths, otherwise declared whole GFF regions."""
    lengths = {}
    kind = source.get("kind", "gff")
    path = Path(source["path"])

    def add(seqid, length):
        if seqid in lengths or length < 1:
            raise ValueError(f"Duplicate or invalid chromosome length: {seqid}")
        lengths[seqid] = length

    with open_text(path) as handle:
        if kind == "genome":
            seqid, length = None, 0
            for line in handle:
                if line.startswith(">"):
                    if seqid is not None:
                        add(seqid, length)
                    fields = line[1:].split()
                    if not fields:
                        raise ValueError("Empty genome FASTA identifier")
                    seqid, length = fields[0], 0
                elif line.strip():
                    if seqid is None:
                        raise ValueError("Genome sequence precedes FASTA header")
                    sequence = line.strip()
                    if not sequence.isascii() or set(sequence.upper()) - set("ACGTURYSWKMBDHVN"):
                        raise ValueError(f"Invalid genome FASTA sequence: {seqid}")
                    length += len(sequence)
            if seqid is not None:
                add(seqid, length)
        elif kind == "sizes":
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.split()
                if fields[:2] in (["seqid", "length"], ["seqid", "length_bp"]):
                    continue
                if len(fields) < 2:
                    raise ValueError("Chromosome sizes require seqid and length columns")
                add(fields[0], int(fields[1]))
        else:
            for line in handle:
                if line.startswith("##FASTA"):
                    break
                if line.startswith("##sequence-region"):
                    fields = line.split()
                    if len(fields) != 4 or int(fields[2]) < 1 or int(fields[3]) < int(fields[2]):
                        raise ValueError("Invalid GFF sequence-region directive")
                    if int(fields[2]) == 1:
                        add(fields[1], int(fields[3]))
    available = {gene.seqid for gene in genes}
    missing = sorted(available - lengths.keys())
    if minimum and missing:
        raise ValueError(f"Missing physical chromosome lengths for {missing[:10]}; provide target_genome/query_genome "
                         "or target_sizes/query_sizes in the pair table (gene extents are not assembly lengths)")
    if any(gene.end > lengths[gene.seqid] for gene in genes if gene.seqid in lengths):
        raise ValueError("An annotation coordinate exceeds the declared chromosome length")
    return lengths, {"kind": kind, "path": str(path), "unknown_seqids": missing}


def homoeolog_order(selected, genomes, anchors):
    """Greedily group supported 2x2 chromosome tiles; no ancestry inference."""
    lookup = [{gene.gene_id: gene.seqid for gene in genome} for genome in genomes]
    allowed = [set(seqids) for seqids in selected]
    counts, seen = Counter(), set()
    with Path(anchors).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2 or any(fields[i] not in lookup[i] for i in (0, 1)):
                raise ValueError(f"Invalid dotplot anchor: {line.rstrip()}")
            pair = tuple(fields[:2])
            if pair in seen:
                continue
            seen.add(pair)
            seqids = tuple(lookup[i][pair[i]] for i in (0, 1))
            if all(seqids[i] in allowed[i] for i in (0, 1)):
                counts[seqids] += 1
    ranks = [{seqid: i for i, seqid in enumerate(track)} for track in selected]
    # Enumerate only pairs sharing partners, not all combinations of unplaced scaffolds.
    partners = {a: {b for (x, b), count in counts.items() if x == a and count} for a in selected[0]}
    candidates = []
    for a, b in combinations(selected[0], 2):
        common = sorted(partners[a] & partners[b], key=ranks[1].__getitem__)
        for c, d in combinations(common, 2):
            support = [counts[(x, y)] for x in (a, b) for y in (c, d)]
            score = 4 / sum((Fraction(1, n) for n in support), Fraction())
            candidates.append((score, sum(support), (a, b), (c, d), support))
    candidates.sort(key=lambda item: (-item[0], -item[1],
                                     tuple(ranks[0][x] for x in item[2]), tuple(ranks[1][x] for x in item[3])))
    used = [set(), set()]
    groups = []
    for score, total, target, query, support in candidates:
        if any(used[i].intersection(track) for i, track in enumerate((target, query))):
            continue
        for i, track in enumerate((target, query)):
            used[i].update(track)
        groups.append({"target": list(target), "query": list(query), "anchor_counts_2x2": support,
                       "harmonic_support": float(score), "unique_anchor_count": total})
    groups.sort(key=lambda group: min(ranks[0][x] for x in group["target"]))
    ordered = [[x for group in groups for x in group[side]] + [x for x in selected[i] if x not in used[i]]
               for i, side in enumerate(("target", "query"))]
    return ordered, {"method": "greedy_2x2_harmonic_unique_anchor_support", "groups": groups,
                     "interpretation": "display grouping only; not a test of homoeology or whole-genome duplication",
                     "globally_optimal": False}


def prepare_dotplot(directory, pair, genomes, karyotype_order, minimum, mode):
    lengths, sources, selected = [], [], []
    for i, side in enumerate(("target", "query")):
        values, metadata = chromosome_lengths(pair.get("dotplot_lengths", {}).get(side, {
            "kind": "gff", "path": pair[side]["gff"]}), genomes[i], minimum)
        available = list(dict.fromkeys(gene.seqid for gene in genomes[i]))
        reference = karyotype_order[i] + [x for x in available if x not in karyotype_order[i]]
        selected.append([x for x in reference if not minimum or values[x] >= minimum])
        lengths.append(values)
        sources.append(metadata)
    metadata = {"mode": mode, "minimum_length_bp": minimum, "length_sources": sources,
                "orientation_changed": False, "coordinate_system": "gene_rank"}
    if mode == "homoeolog":
        selected, grouping = homoeolog_order(selected, genomes, directory / "target.query.lifted.anchors")
        metadata.update(grouping)
    elif mode == "none":
        selected = [[x for x in dict.fromkeys(gene.seqid for gene in genome) if x in set(track)]
                    for genome, track in zip(genomes, selected, strict=True)]
    elif mode != "karyotype":
        raise ValueError("dotplot-sort must be karyotype, homoeolog or none")
    metadata["display_order"] = dict(zip(("target", "query"), selected, strict=True))
    ids = []
    for side, genome, track in zip(("target", "query"), genomes, selected, strict=True):
        rank = {seqid: i for i, seqid in enumerate(track)}
        # Stable chromosome-only sorting retains every within-chromosome BED rank.
        genes = sorted((gene for gene in genome if gene.seqid in rank), key=lambda gene: rank[gene.seqid])
        ids.append({gene.gene_id for gene in genes})
        with (directory / f"dotplot.{side}.bed").open("w", encoding="utf-8") as handle:
            for gene in genes:
                handle.write(f"{gene.seqid}\t{gene.start}\t{gene.end}\t{gene.gene_id}\n")
    count, original, blocks = 0, 0, 0
    pending_separator = True
    all_ids = [{gene.gene_id for gene in genome} for genome in genomes]
    with (directory / "target.query.lifted.anchors").open(encoding="utf-8") as source, \
            (directory / "dotplot.anchors").open("w", encoding="utf-8") as dest:
        for line in source:
            if not line.strip():
                continue
            if line.startswith("#"):
                pending_separator = True
                continue
            fields = line.split()
            if len(fields) < 2 or any(fields[i] not in all_ids[i] for i in (0, 1)):
                raise ValueError("Invalid dotplot anchor")
            original += 1
            if all(fields[i] in ids[i] for i in (0, 1)):
                if pending_separator:
                    dest.write("###\n")
                    pending_separator = False
                    blocks += 1
                dest.write(line)
                count += 1
    if not count:
        raise ValueError(f"No anchors connect chromosomes meeting the dotplot minimum length ({minimum} bp)")
    metadata.update(anchor_count=count, original_anchor_count=original, excluded_anchor_count=original - count,
                    displayed_block_count=blocks, downsampled=False, filtered_by_dS=False)
    with (directory / "dotplot_display.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("species", "seqid", "length_bp", "selected_for_dotplot", "display_rank"))
        for i, side in enumerate(("target", "query")):
            ranks = {x: j + 1 for j, x in enumerate(selected[i])}
            for seqid in dict.fromkeys(gene.seqid for gene in genomes[i]):
                writer.writerow((pair[f"{side}_species"], seqid, lengths[i].get(seqid, ""), int(seqid in ranks), ranks.get(seqid, "")))
    return metadata

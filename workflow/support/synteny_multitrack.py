#!/usr/bin/env python3
"""Optimize and draw a saved multi-track JCVI layout without rerunning synteny."""

import argparse
import hashlib
import json
import os
import shutil
import tempfile
from collections import Counter
from pathlib import Path

from kffractbias.io import read_bed

try:
    from busco_reference_quality import safe_token
    from input_generation_array_state import FreshDigestBatch
    from pairwise_synteny_style import SOFT_COLORS, STYLE, style_text
    from synteny_layout import OPTIONS, TRACK_FIELDS, optimize_layout, validate_options
except ImportError:
    from .busco_reference_quality import safe_token
    from .input_generation_array_state import FreshDigestBatch
    from .pairwise_synteny_style import SOFT_COLORS, STYLE, style_text
    from .synteny_layout import OPTIONS, TRACK_FIELDS, optimize_layout, validate_options


def read_json(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate layout JSON key: {key}")
            result[key] = value
        return result
    with Path(path).open(encoding="utf-8") as handle:
        return json.load(handle, object_pairs_hook=unique)


def load_panels(path, pairs):
    """Workflow schema: select completed comparisons, never arbitrary new data."""
    payload = read_json(path)
    if not isinstance(payload, dict) or set(payload) != {"panels"} or not isinstance(payload["panels"], list) or not payload["panels"]:
        raise ValueError("Karyotype panels file requires a nonempty panels list")
    analyses = {p["analysis_id"]: p for p in pairs}
    identifiers, panels = set(), []
    for panel in payload["panels"]:
        if not isinstance(panel, dict) or set(panel) - ({"panel_id", "tracks", "comparisons"} | OPTIONS):
            raise ValueError("Unknown karyotype panel options")
        panel_id = panel.get("panel_id", "")
        safe_token(panel_id, "panel_id")
        if panel_id in identifiers:
            raise ValueError("Karyotype panel_id values must be unique")
        identifiers.add(panel_id)
        tracks, comparisons = panel.get("tracks"), panel.get("comparisons")
        if (not isinstance(tracks, list) or len(tracks) < 2
                or any(not isinstance(t, dict) or set(t) - (TRACK_FIELDS - {"bed"}) for t in tracks)):
            raise ValueError("A panel requires at least two tracks with documented options")
        species = [t.get("species") for t in tracks]
        if any(not isinstance(s, str) or not s for s in species) or len(set(species)) != len(species):
            raise ValueError("Panel tracks require unique species")
        if (not isinstance(comparisons, list) or not comparisons or any(not isinstance(s, str) for s in comparisons)
                or len(set(comparisons)) != len(comparisons) or set(comparisons) - set(analyses)):
            raise ValueError("Panel comparisons must name unique analysis_id values in the pair table")
        connected, edges = set(), set()
        for comparison in comparisons:
            pair = analyses[comparison]
            ends = (pair["target"]["species"], pair["query"]["species"])
            if set(ends) - set(species) or frozenset(ends) in edges:
                raise ValueError("Panel comparisons must connect its tracks once per species pair")
            connected.update(ends)
            edges.add(frozenset(ends))
        if connected != set(species):
            raise ValueError("Every panel track must participate in a comparison")
        validate_options(panel)
        panels.append(panel)
    return panels


def render_panel(config, directory, formats=("pdf", "svg", "png"), scale_mode="shared"):
    """Create a new self-contained plot bundle; never replace source inputs."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from jcvi.graphics.chromosome import HorizontalChromosome
    from jcvi.graphics.karyotype import Layout, ShadeManager, Track
    from matplotlib.backends.backend_pdf import RendererPdf
    from matplotlib.colors import is_color_like
    from matplotlib.font_manager import FontProperties
    from matplotlib.patches import PathPatch

    if (not isinstance(config, dict) or set(config) - ({"tracks", "edges"} | OPTIONS)
            or not isinstance(config.get("tracks"), list) or not isinstance(config.get("edges"), list)):
        raise ValueError("Saved layout requires tracks, edges and documented optimizer options")
    validate_options(config)
    formats = tuple(formats)
    if not formats or len(set(formats)) != len(formats) or set(formats) - {"pdf", "svg", "png"}:
        raise ValueError("Layout formats must be a unique nonempty subset of pdf,svg,png")
    directory = Path(directory).resolve()
    if directory.exists():
        raise ValueError("Multi-track output directory must be new")
    tracks, edges = config["tracks"], config["edges"]
    if (len(tracks) < 2 or any(not isinstance(t, dict) or "bed" not in t or "species" not in t for t in tracks)
            or any(not isinstance(e, dict) or set(e) != {"tracks", "blocks"} for e in edges)):
        raise ValueError("Layout tracks require species/bed and edges require tracks/blocks")
    proof = FreshDigestBatch()
    inputs = [Path(t["bed"]).resolve() for t in tracks] + [Path(e["blocks"]).resolve() for e in edges]
    if any(directory == p or directory in p.parents for p in inputs):
        raise ValueError("Multi-track output cannot contain a source input")
    proof.read(inputs)
    if any(any(c in str(p) for c in ",\n\r\0") for p in inputs):
        raise ValueError("JCVI layout paths cannot contain commas or control characters")
    genomes = [read_bed(t["bed"]) for t in tracks]
    labels = []
    for track, genes in zip(tracks, genomes, strict=True):
        seqids = track.get("seqids", list(dict.fromkeys(g.seqid for g in genes)))
        values = track.get("labels", {})
        if (not isinstance(values, dict) or set(values) - set(seqids)
                or any(not isinstance(v, str) or not v.strip() or any(c in v for c in "\n\r\0") for v in values.values())):
            raise ValueError("Track labels must map selected seqids to nonempty single-line text")
        row = {s: values.get(s, s) for s in seqids}
        if len(set(row.values())) != len(row):
            raise ValueError("Track display labels must be unique within each species")
        labels.append(row)
    width = 7.2 * 72
    renderer = RendererPdf(None, 300, 4, 7.2)
    with matplotlib.rc_context(STYLE):
        font = FontProperties(family="Helvetica", size=8)
        italic = FontProperties(family="Helvetica", size=8, style="italic")
        label_width = max(renderer.get_text_width_height_descent(t.get("species", "").replace("_", " "), italic, False)[0] for t in tracks)
        identifiers = [label + " (rev)" for row in labels for label in row.values()]
        label_height = max((renderer.get_text_width_height_descent(s, font, False)[0] for s in identifiers), default=0) + 6
        # Leave additional clearance for the raster font substitute as well as
        # native Helvetica PDF metrics; species labels must not touch bars.
        left, right = (12 + label_width + 16) / width, (width - 6) / width
        if right - left < .3:
            raise ValueError("Species labels leave insufficient space at 8 pt; use shorter species labels")
        step, top, bottom = max(70., label_height + 18), label_height + 10, 42.
        height = top + step * (len(tracks) - 1) + bottom
        ys = [(height - top - i * step) / height for i in range(len(tracks))]
        ordered, reversed_by_track, report = optimize_layout(tracks, genomes, edges, scale_mode=scale_mode,
            xstart=left, xend=right, figsize=(width / 72, height / 72), ys=ys,
            gaps=[8 / width] * len(tracks),
            **{k: config[k] for k in OPTIONS if k in config})
        colors = []
        for t, selected in zip(tracks, ordered, strict=True):
            values = t.get("colors", {})
            if not isinstance(values, dict) or set(values) - set(selected) or any(not is_color_like(c) for c in values.values()):
                raise ValueError("Track colors must map selected seqids to valid matplotlib colors")
            colors.append({s: values.get(s, SOFT_COLORS[i % len(SOFT_COLORS)])
                           for i, s in enumerate(report["input_order"][t["species"]])})
        proof.check()
        directory.mkdir(parents=True)
        try:
            copied = directory / "inputs"
            copied.mkdir()
            lines = ["# y,xstart,xend,rotation,color,label,va,bed"]
            for i, t in enumerate(tracks):
                shutil.copyfile(t["bed"], copied / f"track{i}.bed")
                lines.append(f"{ys[i]},{left},{right},0,gainsboro,{t['species'].replace('_', ' ')},top,inputs/track{i}.bed")
            for i, edge in enumerate(edges):
                shutil.copyfile(edge["blocks"], copied / f"edge{i}.source.simple")
                a, b = edge["tracks"]
                lookup = {g.gene_id: g.seqid for g in genomes[a]}
                rows = []
                for line in Path(edge["blocks"]).read_text(encoding="utf-8").splitlines():
                    if not line.strip() or (line.startswith("#") and "*" not in line.split()[0]) or line.startswith("StartGeneA"):
                        continue
                    gene = line.split()[0].split("*", 1)[-1]
                    rows.append(line if "*" in line.split()[0] else colors[a].get(lookup[gene], "#cccccc") + "*" + line)
                (copied / f"edge{i}.simple").write_text("\n".join(rows) + "\n", encoding="utf-8")
                lines.append(f"e,{a},{b},inputs/edge{i}.simple")
            (directory / "layout").write_text("\n".join(lines) + "\n", encoding="utf-8")
            (directory / "seqids").write_text("\n".join(",".join(s + ("-" if reverse[s] else "") for s in row)
                for row, reverse in zip(ordered, reversed_by_track, strict=True)) + "\n", encoding="utf-8")
            saved = {k: config[k] for k in OPTIONS if k in config}
            saved["tracks"] = [{**t, "bed": f"inputs/track{i}.bed"} for i, t in enumerate(tracks)]
            saved["edges"] = [{"tracks": e["tracks"], "blocks": f"inputs/edge{i}.source.simple"} for i, e in enumerate(edges)]
            (directory / "manifest.json").write_text(json.dumps(saved, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            absolute = []
            for line in lines:
                fields = line.split(",")
                if line.startswith("#"):
                    absolute.append(line)
                    continue
                index = 3 if fields[0] == "e" else 7
                fields[index] = str(directory / fields[index])
                absolute.append(",".join(fields))
            with tempfile.TemporaryDirectory(prefix="gg-multitrack-") as temporary:
                path = Path(temporary) / "layout"
                path.write_text("\n".join(absolute) + "\n", encoding="utf-8")
                layout = Layout(str(path), generank=True, seed=config.get("seed", 1))
            fig = plt.figure(figsize=(width / 72, height / 72), dpi=72)
            try:
                ax = fig.add_axes([0, 0, 1, 1])
                ax.set(xlim=(0, 1), ylim=(0, 1))
                ax.axis("off")
                for i, entry in enumerate(layout):
                    entry.seqids = ordered[i]
                    entry.rev = {s for s, v in reversed_by_track[i].items() if v}
                    counts = Counter(g.seqid for g in genomes[i])
                    entry.sizes = {s: counts[s] for s in entry.seqids}
                native = [Track(ax, t, draw=False) for t in layout]
                for i, track in enumerate(native):
                    track.ratio, track.gap = report["track_ratios"][i], report["gaps"][i]
                    track.xend = track.xstart + track.total * track.ratio + (len(track.seqids) - 1) * track.gap
                    track.update_offsets()
                ShadeManager(ax, native, layout, heightpad=1.5 / height)
                ribbons = [p for p in ax.patches if isinstance(p, PathPatch)]
                if len(ribbons) != sum(report["displayed_block_counts"]):
                    raise ValueError("Native JCVI ribbon count differs from optimizer inputs")
                for patch in ribbons:
                    patch.set_alpha(.52)
                    patch.set_linewidth(0)
                for i, track in enumerate(native):
                    ax.text(12 / width, track.y, tracks[i]["species"].replace("_", " "), va="center", fontstyle="italic")
                    for sid in track.seqids:
                        start, span = track.offsets[sid], track.sizes[sid] * track.ratio
                        glyph = HorizontalChromosome(ax, start, start + span, track.y, height=3 / height,
                            lw=.3, ec="#555555", fc=colors[i][sid], style="rect")
                        glyph.set_transform(track.tr)
                        label = labels[i][sid] + (" (rev)" if sid in track.rev else "")
                        ax.annotate(label, (start + span / 2, track.y), xytext=(0, 6), textcoords="offset points",
                                    rotation=90, rotation_mode="anchor", ha="left", va="center")
                try:
                    from pairwise_synteny_karyotype import gene_scale
                except ImportError:
                    from .pairwise_synteny_karyotype import gene_scale
                scale = gene_scale([t.total for t in native])
                scale_tracks = [len(native) - 1] if scale_mode == "shared" else list(range(len(native)))
                # Independent scales are drawn beside the species labels, not
                # incorrectly represented by a single shared bar.
                for i in scale_tracks:
                    y = 18 / height if scale_mode == "shared" else native[i].y - 18 / height
                    ax.plot([left, left + scale * native[i].ratio], [y, y], color="black", lw=.6)
                    ax.text(left, y + 4 / height, f"{scale:,} genes", va="bottom")
                ax.text(12 / width, 5 / height, "Ribbons: syntenic blocks; (rev): reversed for display", va="bottom")
                style_text(fig)
                for text in ax.texts:
                    if text.get_text() in {t["species"].replace("_", " ") for t in tracks}:
                        text.set_fontstyle("italic")
                for fmt in formats:
                    proof.check()
                    fig.savefig(directory / ("karyotype." + fmt), format=fmt, dpi=300)
            finally:
                plt.close(fig)
            report.update(font="Helvetica", font_size_pt=8, pdf_width_inches=7.2,
                          scale_value=scale, scale_bar_count=len(scale_tracks),
                          display_labels=dict(zip(report["species_order"], labels, strict=True)),
                          source_hashes={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in set(inputs)})
            (directory / "karyotype_layout.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            (directory / "karyotype_colors.json").write_text(json.dumps(dict(zip(report["species_order"], colors, strict=True)), indent=2) + "\n")
            proof.check()
        except BaseException:
            shutil.rmtree(directory)
            raise
    return report


def render_workflow_panels(plan, source, output):
    pairs = {p["analysis_id"]: p for p in plan["pairs"]}
    proof = FreshDigestBatch()
    proof.read([source / comparison / name for panel in plan.get("karyotype_panels", [])
                for comparison in panel["comparisons"] for name in ("target.bed", "query.bed", "target.query.screened.simple")])
    for panel in plan.get("karyotype_panels", []):
        tracks = [dict(t) for t in panel["tracks"]]
        indices = {t["species"]: i for i, t in enumerate(tracks)}
        beds, edges = {}, []
        for comparison in panel["comparisons"]:
            pair = pairs[comparison]
            for side in ("target", "query"):
                species = pair[side]["species"]
                path = source / comparison / (side + ".bed")
                if species in beds and path.read_bytes() != beds[species].read_bytes():
                    raise ValueError(f"Shared panel species has different gene coordinates/identifiers: {species}")
                beds[species] = path
            edges.append({"tracks": [indices[pair[s]["species"]] for s in ("target", "query")],
                          "blocks": str(source / comparison / "target.query.screened.simple")})
        for track in tracks:
            track["bed"] = str(beds[track["species"]])
        config = {k: panel[k] for k in OPTIONS if k in panel}
        config.update(tracks=tracks, edges=edges)
        render_panel(config, output / "panels" / panel["panel_id"], plan["formats"], plan.get("karyotype_scale", "shared"))
        proof.check()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True, help="Saved layout JSON; relative BED/block paths resolve beside it")
    parser.add_argument("--output", type=Path, required=True, help="New output directory")
    parser.add_argument("--formats", default="pdf,svg,png")
    parser.add_argument("--scale", choices=("shared", "independent"), default="shared")
    args = parser.parse_args()
    proof = FreshDigestBatch()
    proof.read([args.manifest])
    config = read_json(args.manifest)
    for row, key in [*((t, "bed") for t in config["tracks"]), *((e, "blocks") for e in config["edges"])]:
        path = Path(row[key])
        row[key] = str(path if path.is_absolute() else args.manifest.resolve().parent / path)
    output = args.output.resolve()
    if output.exists():
        raise ValueError("Multi-track output directory must be new")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".gg-layout-", dir=output.parent) as temporary:
        staged = Path(temporary) / "panel"
        render_panel(config, staged, args.formats.split(","), args.scale)
        proof.check()
        if output.exists():
            raise ValueError("Multi-track destination appeared during rendering")
        os.rename(staged, output)


if __name__ == "__main__":
    main()

import csv
import json
import os
import random
import re
import subprocess
import sys
from pathlib import Path

import pytest
from Bio.Data import CodonTable
from PIL import Image

from workflow.tests.test_genome_evolution_protein_mode import _run_core

REPO_ROOT = Path(__file__).resolve().parents[2]
CORE = REPO_ROOT / "workflow/core/gg_genome_evolution_core.sh"


def test_pairwise_synteny_help_uses_real_runtime_dependency(tmp_path):
    script = REPO_ROOT / "workflow/support/pairwise_synteny.py"
    result = subprocess.run([sys.executable, str(script), "--help"], cwd=tmp_path,
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "usage:" in result.stdout.lower()
    assert list(tmp_path.iterdir()) == []


@pytest.fixture
def saved_karyotype(tmp_path):
    directory = tmp_path / "saved"
    directory.mkdir()
    for side, prefix in (("target", "t"), ("query", "q")):
        (directory / f"{side}.bed").write_text(
            "".join(f"chr1\t{i * 10}\t{i * 10 + 5}\t{prefix}{i}\n" for i in range(6)))
    (directory / "seqids").write_text("chr1\nchr1\n")
    (directory / "layout").write_text(
        "0.7,0.12,0.92,0,,Target species,top,target.bed,top\n"
        "0.3,0.12,0.92,0,,Query species,bottom,query.bed,bottom\n"
        "e,0,1,colored.simple\n")
    (directory / "colored.simple").write_text("#88afc4*t0 t2 q1 q4 5 +\n")
    colors = {"chromosomes": {side: {"chr1": "#88afc4"} for side in ("target", "query")}}
    pair = {"target_species": "Target_species", "query_species": "Query_species"}
    return directory, pair, colors


def test_karyotype_saved_relative_inputs_resolve_beside_layout(saved_karyotype, tmp_path, monkeypatch):
    from workflow.support.pairwise_synteny_karyotype import render_karyotype

    directory, pair, colors = saved_karyotype
    # A caller's working directory can contain different files with the same names.
    for side, prefix in (("target", "t"), ("query", "q")):
        (tmp_path / f"{side}.bed").write_text(
            "".join(f"chr1\t{i * 10}\t{i * 10 + 5}\t{prefix}{i}\n" for i in range(12)))
    (tmp_path / "colored.simple").write_text("#88afc4*t0 t2 q1 q4 5 +\n")
    monkeypatch.chdir(tmp_path)
    original = (directory / "layout").read_bytes()
    style = render_karyotype(directory, pair, "svg", colors)
    assert style["track_gene_counts"] == [6, 6]
    assert (directory / "layout").read_bytes() == original


@pytest.mark.parametrize("invalid", ["duplicate-seqid", "unknown-seqid", "duplicate-gene"])
def test_karyotype_invalid_saved_coordinates_preserve_plot_and_close_figure(saved_karyotype, monkeypatch, invalid):
    import matplotlib.pyplot as plt

    from workflow.support.pairwise_synteny_karyotype import render_karyotype

    directory, pair, colors = saved_karyotype
    monkeypatch.chdir(directory)
    if invalid == "duplicate-seqid":
        (directory / "seqids").write_text("chr1,chr1\nchr1\n")
    elif invalid == "unknown-seqid":
        (directory / "seqids").write_text("absent\nchr1\n")
    else:
        with (directory / "target.bed").open("a") as handle:
            handle.write("chr1\t100\t105\tt0\n")
    plot = directory / "karyotype.svg"
    plot.write_text("previous plot")
    figures = plt.get_fignums()
    with pytest.raises(ValueError):
        render_karyotype(directory, pair, "svg", colors)
    assert plot.read_text() == "previous plot"
    assert plt.get_fignums() == figures


def test_karyotype_failed_save_preserves_previous_plot(saved_karyotype, monkeypatch):
    from matplotlib.figure import Figure

    from workflow.support.pairwise_synteny_karyotype import render_karyotype

    directory, pair, colors = saved_karyotype
    monkeypatch.chdir(directory)
    plot = directory / "karyotype.svg"
    plot.write_text("previous plot")

    def fail_save(figure, path, **kwargs):
        Path(path).write_text("partial plot")
        raise OSError("simulated disk failure")

    monkeypatch.setattr(Figure, "savefig", fail_save)
    with pytest.raises(OSError, match="simulated disk failure"):
        render_karyotype(directory, pair, "svg", colors)
    assert plot.read_text() == "previous plot"


def test_karyotype_failed_style_publication_rolls_back_plot_and_style(saved_karyotype, monkeypatch):
    import nwkit.output_transaction as transaction

    from workflow.support.pairwise_synteny_karyotype import render_karyotype

    directory, pair, colors = saved_karyotype
    plot, report = directory / "karyotype.svg", directory / "karyotype_style.json"
    plot.write_text("previous plot")
    report.write_text("previous style")
    rename = transaction.os.replace

    def fail_install(source, destination):
        if Path(destination) == report and ".stage." in Path(source).name:
            raise OSError("simulated style installation failure")
        return rename(source, destination)

    monkeypatch.setattr(transaction.os, "replace", fail_install)
    with pytest.raises(OSError, match="simulated style installation failure"):
        render_karyotype(directory, pair, "svg", colors, layout_report=report)
    assert plot.read_text() == "previous plot" and report.read_text() == "previous style"


@pytest.mark.parametrize("defect", ["cross-chromosome", "orientation"])
def test_karyotype_rejects_malformed_saved_blocks(saved_karyotype, monkeypatch, defect):
    from workflow.support.pairwise_synteny_karyotype import render_karyotype

    directory, pair, colors = saved_karyotype
    monkeypatch.chdir(directory)
    if defect == "orientation":
        (directory / "colored.simple").write_text("#88afc4*t0 t2 q1 q4 5 sideways\n")
    else:
        (directory / "query.bed").write_text(
            "".join(f"chr{1 if i < 3 else 2}\t{i * 10}\t{i * 10 + 5}\tq{i}\n" for i in range(6)))
        (directory / "seqids").write_text("chr1\nchr1,chr2\n")
        colors["chromosomes"]["query"]["chr2"] = "#88afc4"
    plot = directory / "karyotype.svg"
    plot.write_text("previous plot")
    with pytest.raises(ValueError, match="block|ribbon"):
        render_karyotype(directory, pair, "svg", colors)
    assert plot.read_text() == "previous plot"


def write_genome(workspace, species, mode, chromosomes, seed_offsets=None):
    sequence_dir = workspace / "input" / f"species_{mode}"
    annotation_dir = workspace / "input/species_gff"
    sequence_dir.mkdir(parents=True, exist_ok=True)
    annotation_dir.mkdir(parents=True, exist_ok=True)
    fasta = []
    gff = ["##gff-version 3"]
    codons = {}
    for codon, aa in CodonTable.standard_dna_table.forward_table.items():
        codons.setdefault(aa, codon)
    for chromosome_index, (chromosome, reverse) in enumerate(chromosomes):
        gff.append(f"##sequence-region {chromosome} 1 2000000")
        indices = list(range(8))
        if reverse:
            indices.reverse()
        for position, index in enumerate(indices):
            offset = 0 if seed_offsets is None else seed_offsets[chromosome_index]
            rng = random.Random(1729 + offset + index)
            protein = "M" + "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(160))
            identifier = f"{chromosome}_g{index}"
            fasta_identifier = species + "_" + identifier if mode == "protein" else identifier
            sequence = protein if mode == "protein" else "".join(codons[aa] for aa in protein)
            fasta.extend((f">{fasta_identifier}", sequence))
            start = position * 1000 + 1
            gff.extend((f"{chromosome}\ttest\tgene\t{start}\t{start + 600}\t.\t+\t.\tID=locus_{identifier}",
                        f"{chromosome}\ttest\tmRNA\t{start}\t{start + 600}\t.\t+\t.\tID={identifier};Parent=locus_{identifier}"))
    (sequence_dir / f"{species}.{mode}.fa").write_text("\n".join(fasta) + "\n")
    (annotation_dir / f"{species}.gff3").write_text("\n".join(gff) + "\n")


def run_core(workspace, **overrides):
    env = {**os.environ, "gg_workspace_dir": str(workspace), "GG_COMMON_TMP_ROOT": "workspace",
           "GG_ARRAY_TASK_ID": "1", "GG_JOB_ID": "synteny_test", "GG_TASK_CPUS": "1", "GG_MEM_PER_CPU_GB": "8",
           "genome_evolution_mode": "synteny", "artifact_stale_policy": "stop", **overrides}
    return subprocess.run(["bash", str(CORE)], cwd=REPO_ROOT, env=env, capture_output=True, text=True, timeout=180)


def test_ds_stage_uses_cdskit_reuses_analysis_and_has_independent_cache(tmp_path):
    workspace = tmp_path / "workspace"
    for species in ("Triphyophyllum_peltatum", "Ancistrocladus_abbreviatus"):
        for mode in ("protein", "cds"):
            write_genome(workspace, species, mode, (("Chr1", False),))
    (workspace / "input/synteny_pairs.tsv").write_text(
        "analysis_id\ttarget_species\tquery_species\n"
        "pair\tTriphyophyllum_peltatum\tAncistrocladus_abbreviatus\n")
    result = run_core(workspace, synteny_karyotype_sort="none", synteny_plot_formats="png")
    assert result.returncode == 0, result.stdout + result.stderr
    root = workspace / "output/genome_evolution/synteny"
    analysis = root / "analysis/pair"
    analysis_mtime = (analysis / "commands.json").stat().st_mtime_ns
    anchors = (analysis / "target.query.lifted.anchors").read_bytes()
    result = run_core(workspace, synteny_dotplot_color="ds", synteny_plot_only="1", artifact_stale_policy="rebuild")
    assert result.returncode != 0
    assert "current ds" in result.stderr
    assert not (root / "ds").exists()
    result = run_core(workspace, synteny_dotplot_color="ds", synteny_plot_formats="png", artifact_stale_policy="rebuild")
    assert result.returncode == 0, result.stdout + result.stderr
    ds = root / "ds/pair"
    metadata = json.loads((root / "plots/pair/dotplot_ds.json").read_text())
    summary = json.loads((ds / "summary.json").read_text())
    assert summary["tools"]["method"] == "YN00_weighting0_F3x4"
    assert metadata["anchor_count"] == json.loads((analysis / "summary.json").read_text())["anchor_count"]
    assert metadata["downsampled"] is metadata["filtered_by_dS"] is False
    assert (analysis / "commands.json").stat().st_mtime_ns == analysis_mtime
    ds_mtime = (ds / "ds.tsv").stat().st_mtime_ns
    result = run_core(workspace, synteny_dotplot_color="ds", synteny_ds_color_max="1",
                      synteny_plot_only="1", synteny_plot_formats="png", artifact_stale_policy="rebuild")
    assert result.returncode == 0, result.stdout + result.stderr
    assert (ds / "ds.tsv").stat().st_mtime_ns == ds_mtime
    assert json.loads((root / "plots/pair/dotplot_ds.json").read_text())["color_range"] == [0, 1]
    style_file = root / "plots/pair/karyotype_style.json"
    style = json.loads(style_file.read_text())
    assert style["scale_mode"] == "shared" and style["scale_bar_count"] == 1
    assert style["track_ratios"][0] == style["track_ratios"][1]
    result = run_core(workspace, synteny_dotplot_color="ds", synteny_ds_color_max="1",
                      synteny_plot_only="1", synteny_plot_formats="png",
                      synteny_karyotype_scale="independent", artifact_stale_policy="rebuild")
    assert result.returncode == 0, result.stdout + result.stderr
    assert (analysis / "commands.json").stat().st_mtime_ns == analysis_mtime
    assert (ds / "ds.tsv").stat().st_mtime_ns == ds_mtime
    style = json.loads(style_file.read_text())
    assert style["scale_mode"] == "independent" and style["scale_bar_count"] == 2
    cds = workspace / "input/species_cds/Triphyophyllum_peltatum.cds.fa"
    original = cds.read_text()
    assert "GCT" in original
    cds.write_text(original.replace("GCT", "GCC", 1))
    result = run_core(workspace, synteny_dotplot_color="ds", synteny_plot_only="1", artifact_stale_policy="rebuild")
    assert result.returncode != 0
    assert (ds / "ds.tsv").stat().st_mtime_ns == ds_mtime
    result = run_core(workspace, synteny_dotplot_color="ds", synteny_plot_formats="png", artifact_stale_policy="rebuild")
    assert result.returncode == 0, result.stdout + result.stderr
    assert (ds / "ds.tsv").stat().st_mtime_ns != ds_mtime
    assert (analysis / "commands.json").stat().st_mtime_ns == analysis_mtime
    assert (analysis / "target.query.lifted.anchors").read_bytes() == anchors


def test_ds_renderer_retains_missing_and_over_limit_anchors(tmp_path):
    from workflow.support.pairwise_synteny_ds import render_ds_dotplot
    plots, ds = tmp_path / "plots", tmp_path / "ds"
    plots.mkdir()
    ds.mkdir()
    (plots / "target.bed").write_text("chr1\t0\t10\ta\nchr1\t10\t20\tb\nchr1\t20\t30\tc\n")
    (plots / "query.bed").write_text("chr2\t0\t10\td\nchr2\t10\t20\te\nchr2\t20\t30\tf\n")
    (plots / "target.query.lifted.anchors").write_text("###\na d 10\nb e 10\nc f 10\n")
    (ds / "ds.tsv").write_text("pair_id\tdS\tstatus\na|d\t0.5\tok\nb|e\t4\tok\nc|f\t\tsaturated\n")
    render_ds_dotplot(plots, ds, {"target_species": "Target_species", "query_species": "Query_species"}, "png", 2)
    metadata = json.loads((plots / "dotplot_ds.json").read_text())
    assert metadata["anchor_count"] == 3
    assert metadata["above_color_max_count"] == metadata["missing_count"] == 1
    assert metadata["filtered_by_dS"] is False
    with Image.open(plots / "dotplot.png") as image:
        image.verify()


@pytest.mark.parametrize("color", ["ds", "orientation"])
@pytest.mark.parametrize("fmt", ["pdf", "svg", "png"])
def test_dotplot_pdf_3p6inch_square_helvetica8_vertical_chromosomes_and_italic_species(tmp_path, monkeypatch, color, fmt):
    from matplotlib.colors import to_rgba
    from matplotlib.figure import Figure
    from matplotlib.lines import TICKDOWN, TICKLEFT
    from matplotlib.offsetbox import AnnotationBbox, DrawingArea
    from matplotlib.patches import Rectangle
    from matplotlib.text import Text

    from workflow.support.pairwise_synteny_ds import render_ds_dotplot

    plots, ds = tmp_path / "plots", tmp_path / "ds"
    plots.mkdir()
    ds.mkdir()
    (plots / "target.bed").write_text("chr1\t0\t10\ta\nchr1\t10\t20\tb\nchr1\t20\t30\tc\n")
    (plots / "query.bed").write_text("chr2\t0\t10\td\nchr2\t10\t20\te\nchr2\t20\t30\tf\nchr2\t30\t40\tg\n")
    (plots / "target.query.lifted.anchors").write_text("###\na d 10\nb e 10\nc f 10\n")
    (ds / "ds.tsv").write_text("pair_id\tdS\tstatus\na|d\t0.5\tok\nb|e\t4\tok\nc|f\t\tsaturated\n")
    pair = {"target_species": "Target_species", "query_species": "Query_species"}
    saved = Figure.savefig
    observations = []

    def inspect_save(figure, *args, **kwargs):
        figure.canvas.draw()
        axes = figure.axes[1]
        box = axes.get_window_extent()
        assert box.width == pytest.approx(box.height, rel=1e-10)
        assert axes.get_xlim() == (0, 3)
        assert axes.get_ylim() == (4, 0)
        assert len(axes.collections[0].get_offsets()) == 3
        for axis, marker in ((axes.xaxis, TICKDOWN), (axes.yaxis, TICKLEFT)):
            ticks = axis.get_major_ticks()
            assert ticks
            for tick in ticks:
                assert tick.tick1line.get_visible() and not tick.tick2line.get_visible()
                assert tick.tick1line.get_marker() == marker
                assert tick.tick1line.get_markersize() == 3
                assert to_rgba(tick.tick1line.get_markeredgecolor()) == to_rgba("black")
                label_box = tick.label1.get_window_extent()
                if axis is axes.xaxis:
                    assert label_box.y1 < box.y0
                else:
                    assert label_box.x1 < box.x0
        texts = [text for text in figure.findobj(Text) if text.get_visible() and text.get_text()]
        assert texts
        assert all(text.get_fontsize() == 8 and text.get_fontfamily() == ["Helvetica"] for text in texts)
        assert all(to_rgba(text.get_color()) == to_rgba("black") for text in texts)
        assert axes.get_xlabel() == "Target species"
        assert axes.get_ylabel() == "Query species"
        assert axes.xaxis.label.get_fontstyle() == axes.yaxis.label.get_fontstyle() == "italic"
        italic = [text for text in texts if text.get_fontstyle() == "italic"]
        ds_tokens = [text for text in texts if text.get_text() == "dS"]
        assert len(ds_tokens) == (2 if color == "ds" else 0)
        assert len(italic) == 2 + len(ds_tokens)
        assert all(text.get_fontstyle() == "italic" for text in ds_tokens)
        units = [text for text in texts if text.get_text().strip() == "(gene rank)"]
        assert len(units) == 2 and all(text.get_fontstyle() == "normal" for text in units)
        assert all(text.get_text().startswith(" ") for text in units)
        assert all(sum(text.get_position()) > 2 for text in units)
        labels = [text for text in figure.axes[0].texts if text.get_text() == "chr1"]
        assert len(labels) == 1 and labels[0].get_rotation() == 90
        if color == "ds":
            from workflow.support.pairwise_synteny_style import MISSING_COLOR, ds_colormap

            cax = figure.axes[2]
            legend = next(artist for artist in cax.artists
                          if isinstance(artist, AnnotationBbox) and artist.xycoords is cax.transAxes)
            assert legend.xy == (1, 0.5) and legend.xycoords is cax.transAxes
            square = next(child for child in legend.offsetbox.get_children() if isinstance(child, DrawingArea))
            patch = next(child for child in square.get_children() if isinstance(child, Rectangle))
            assert patch.get_width() == patch.get_height() == 6
            assert patch.get_facecolor() == to_rgba(MISSING_COLOR)
            assert axes.collections[0].cmap.name == "genegalleon_ds"
            assert axes.collections[0].get_facecolors() == pytest.approx(
                ds_colormap()([0.25, 1.0, float("nan")]))
        observations.append(box.width / box.height)
        return saved(figure, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", inspect_save)
    if color == "ds":
        render_ds_dotplot(plots, ds, pair, fmt, 2)
    else:
        from workflow.support.pairwise_synteny_dotplot import render_orientation_pdf

        render_orientation_pdf(plots, pair, filtered=False, fmt=fmt)
    assert observations == pytest.approx([1])
    if fmt == "png":
        with Image.open(plots / "dotplot.png") as image:
            image.verify()
        return
    if fmt == "svg":
        svg = (plots / "dotplot.svg").read_text()
        assert svg.count(">dS</text>") == (2 if color == "ds" else 0)
        return
    pdf = (plots / "dotplot.pdf").read_bytes()
    media_box = re.search(rb"/MediaBox\s*\[([^]]+)\]", pdf)
    assert media_box is not None
    x0, _, x1, _ = map(float, media_box.group(1).split())
    assert x1 - x0 == pytest.approx(3.6 * 72, abs=1e-7)
    assert b"/BaseFont /Helvetica" in pdf
    assert b"/BaseFont /Helvetica-Oblique" in pdf
    assert b"/Subtype /Type3" not in pdf
    if color == "ds":
        from pypdf import PdfReader
        from pypdf.generic import ContentStream

        reader = PdfReader(plots / "dotplot.pdf")
        page = reader.pages[0]
        fonts = page["/Resources"]["/Font"]
        tokens, current_font = [], None
        for operands, operation in ContentStream(page.get_contents(), reader).operations:
            if operation == b"Tf":
                current_font = fonts[operands[0]]["/BaseFont"]
                assert float(operands[1]) == 8
            elif operation == b"Tj" and str(operands[0]) == "dS":
                tokens.append(current_font)
        assert tokens == ["/Helvetica-Oblique", "/Helvetica-Oblique"]


@pytest.mark.parametrize("fmt", ["pdf", "svg", "png"])
@pytest.mark.parametrize("scale_mode", ["shared", "independent"])
@pytest.mark.parametrize("with_legend", [False, True])
def test_compact_jcvi_karyotype_has_black_helvetica8_species_only_italic_and_track_scales(tmp_path, monkeypatch, fmt, scale_mode, with_legend):
    from jcvi.formats.bed import Bed
    from kffractbias.io import read_bed
    from matplotlib.backends.backend_pdf import RendererPdf
    from matplotlib.colors import to_rgba
    from matplotlib.figure import Figure
    from matplotlib.offsetbox import AnnotationBbox, DrawingArea
    from matplotlib.patches import Polygon
    from matplotlib.text import Text

    from workflow.support.pairwise_synteny_karyotype import chromosome_colors, render_karyotype

    seqids = ["chromosome_with_a_long_label1", "chromosome_with_a_long_label2"]
    (tmp_path / "target.bed").write_text("".join(f"{seqids[i // 3]}\t{i * 100}\t{i * 100 + 10}\tt{i}\n" for i in range(6)))
    (tmp_path / "query.bed").write_text("".join(f"scaffold{3 if i < 5 else 16}\t{i * 100}\t{i * 100 + 10}\tq{i}\n" for i in range(10)))
    (tmp_path / "seqids").write_text(",".join(seqids) + "\nscaffold3,scaffold16\n")
    (tmp_path / "layout").write_text(
        f"0.7,0.12,0.92,0,,Target species (gene rank),top,{tmp_path / 'target.bed'},top\n"
        f"0.3,0.12,0.92,0,,Query species (gene rank),bottom,{tmp_path / 'query.bed'},bottom\n"
        f"e,0,1,{tmp_path / 'colored.simple'}\n")
    (tmp_path / "colored.simple").write_text("#88afc4*t0 t2 q0 q4 5 +\n")
    genomes = [read_bed(tmp_path / f"{side}.bed") for side in ("target", "query")]
    colors = chromosome_colors([seqids, ["scaffold3", "scaffold16"]], genomes, tmp_path / "unused")
    analysis = None
    if with_legend:
        analysis = tmp_path / "analysis"
        (analysis / "logs").mkdir(parents=True)
        (analysis / "summary.json").write_text(json.dumps({"parameters": {
            "cscore": 0.7, "min_anchors": 4, "distance": 20, "quota": None}}))
        (analysis / "logs/01.mcscan.log").write_text(
            "diamond blastp --evalue 1e-5 --outfmt 6\n"
            "local dups filter (tandem_Nmax=10)\n0 new pairs found (dist=10).\n")
    pair = {"target_species": "Target_species", "query_species": "Query_species"}
    bar_count = 1 if scale_mode == "shared" else 2
    saved = Figure.savefig
    observed = []

    def inspect(figure, *args, **kwargs):
        root = figure.axes[0]
        figure.canvas.draw()
        chromosomes = [patch for patch in root.patches if isinstance(patch, Polygon)]
        assert len(chromosomes) == 8  # Four outlines and four coloured interiors.
        for chromosome in chromosomes:
            vertices = chromosome.get_xy()
            assert len(vertices) == 5 and tuple(vertices[0]) == tuple(vertices[-1])
            assert len({x for x, _ in vertices}) == len({y for _, y in vertices}) == 2
        texts = [text for text in figure.findobj(Text) if text.get_visible() and text.get_text()]
        assert texts and all(text.get_fontsize() == 8 and text.get_fontfamily() == ["Helvetica"] for text in texts)
        assert all(to_rgba(text.get_color()) == to_rgba("black") for text in texts)
        assert {text.get_text() for text in texts if text.get_fontstyle() == "italic"} == {"Target species", "Query species"}
        assert sum(text.get_text() == "1 gene" for text in texts) == bar_count
        assert len(root.lines) == bar_count
        # Each species follows only its own track's longest rotated label.
        # Measure with the same vector font metrics used for page fitting.
        renderer = RendererPdf(None, 300, figure.get_figheight(), figure.get_figwidth())
        species_boxes = []
        for direction, species_name, chromosome_ids in (
                (1, "Target species", seqids), (-1, "Query species", ["scaffold3", "scaffold16"])):
            species = next(text for text in root.texts if text.get_text() == species_name)
            species_box = species.get_window_extent(renderer)
            label_boxes = [text.get_window_extent(renderer) for text in root.texts
                           if text.get_text() in chromosome_ids]
            gap = (species_box.y0 - max(box.y1 for box in label_boxes) if direction == 1
                   else min(box.y0 for box in label_boxes) - species_box.y1)
            assert gap == pytest.approx(4, abs=1e-7)
            species_boxes.append(species_box)
        assert root.lines[-1].get_window_extent(renderer).y0 < species_boxes[1].y0
        for line in root.lines:
            scale_label = next(text for text in root.texts if text.get_text() == "1 gene"
                               and text.xy[1] == line.get_ydata()[0])
            assert scale_label.get_window_extent(renderer).y0 > line.get_window_extent(renderer).y1
        bottom_label = next(text for text in root.texts if text.get_text() == "1 gene"
                            and text.xy[1] == root.lines[-1].get_ydata()[0])
        assert species_boxes[1].y0 - bottom_label.get_window_extent(renderer).y1 >= 4 - 1e-7
        if scale_mode == "independent":
            assert root.lines[0].get_window_extent(renderer).y0 > species_boxes[0].y1
        for line, genes in zip(root.lines, (10,) if bar_count == 1 else (6, 10), strict=True):
            assert line.get_xdata()[1] - line.get_xdata()[0] == pytest.approx(0.79 / genes)
            assert line.get_clip_on() is False
        assert all(text.get_annotation_clip() is False for text in root.texts if text.get_text() == "1 gene")
        legends = [artist for artist in root.artists if isinstance(artist, AnnotationBbox)]
        assert len(legends) == int(with_legend)
        if with_legend:
            symbol = next(child for child in legends[0].offsetbox.get_children() if isinstance(child, DrawingArea))
            assert len(symbol.get_children()) == 3
            assert "Connections: syntenic blocks (MCscan + liftover)" in "\n".join(text.get_text() for text in texts)
            assert legends[0].get_window_extent().y1 < min(line.get_window_extent().y0 for line in root.lines)
        assert kwargs["bbox_inches"].width == pytest.approx(7.2)
        observed.append(True)
        return saved(figure, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", inspect)
    style = render_karyotype(tmp_path, pair, fmt, colors, scale_mode, analysis)
    assert observed == [True]
    assert style["scale_unit"] == "genes" and style["scale_value"] == 1
    assert style["scale_mode"] == scale_mode and style["scale_bar_count"] == bar_count
    assert style["scale_label_position"] == "above"
    assert style["chromosome_style"] == "rectangular"
    placement = style["species_label_layout"]
    assert [item["longest_chromosome_label"] for item in placement] == [seqids[0], "scaffold16"]
    assert [item["species_label_gap_pt"] for item in placement] == [4, 4]
    assert abs(placement[1]["species_offset_pt"]) < abs(placement[0]["species_offset_pt"])
    assert style["track_ratios"] == pytest.approx([0.79 / 10, 0.79 / 10] if bar_count == 1 else [0.79 / 6, 0.79 / 10])
    assert [len(Bed(str(tmp_path / f"{side}.bed"))) for side in ("target", "query")] == [6, 10]
    assert (style["connection_criteria"] is not None) == with_legend
    if with_legend:
        assert "|dx| + |dy| < 10 gene ranks" in style["connection_legend_text"]
    output = tmp_path / f"karyotype.{fmt}"
    if fmt == "pdf":
        from pypdf import PdfReader

        data = output.read_bytes()
        media_box = re.search(rb"/MediaBox\s*\[([^]]+)\]", data)
        x0, _, x1, _ = map(float, media_box.group(1).split())
        assert x1 - x0 == pytest.approx(7.2 * 72, abs=1e-7)
        assert b"/BaseFont /Helvetica-Oblique" in data and b"/Subtype /Type3" not in data
        text = PdfReader(output).pages[0].extract_text()
        assert text.count("1 gene") == bar_count
        if with_legend:
            assert "E <= 1e-5" in text and "C-score >= 0.7" in text
            assert "gap <= 20 gene ranks / genome" in text and "no quota; no dS filter" in text
    elif fmt == "svg":
        assert "Helvetica" in output.read_text()
    else:
        with Image.open(output) as image:
            image.verify()


@pytest.mark.parametrize("scale_mode", ["shared", "independent"])
def test_karyotype_track_swap_preserves_inputs_ribbon_endpoints_and_scale(tmp_path, monkeypatch, scale_mode):
    import hashlib

    from kffractbias.io import read_bed
    from matplotlib.figure import Figure
    from matplotlib.patches import PathPatch

    from workflow.support.pairwise_synteny_karyotype import chromosome_colors, render_karyotype

    (tmp_path / "target.bed").write_text("".join(f"chr1\t{i * 10}\t{i * 10 + 5}\tt{i}\n" for i in range(6)))
    (tmp_path / "query.bed").write_text("".join(f"scaffold_long\t{i * 10}\t{i * 10 + 5}\tq{i}\n" for i in range(10)))
    (tmp_path / "seqids").write_text("chr1\nscaffold_long\n")
    (tmp_path / "layout").write_text(
        f"0.7,0.12,0.92,0,,Target species,top,{tmp_path / 'target.bed'},top\n"
        f"0.3,0.12,0.92,0,,Query species,bottom,{tmp_path / 'query.bed'},bottom\n"
        f"e,0,1,{tmp_path / 'colored.simple'}\n")
    (tmp_path / "colored.simple").write_text("#88afc4*t0 t2 q1 q4 5 +\n")
    originals = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in tmp_path.iterdir()}
    genomes = [read_bed(tmp_path / f"{side}.bed") for side in ("target", "query")]
    colors = chromosome_colors([["chr1"], ["scaffold_long"]], genomes, tmp_path / "unused")
    pair = {"target_species": "Target_species", "query_species": "Query_species"}
    observed, saved = [], Figure.savefig

    def inspect(figure, *args, **kwargs):
        root = figure.axes[0]
        observed.append({"ribbons": [patch.get_path().vertices.copy() for patch in root.patches
                                     if isinstance(patch, PathPatch)],
                         "species": {text.get_text(): text.xy[1] for text in root.texts
                                     if text.get_text() in {"Target species", "Query species"}},
                         "bars": [(line.get_xdata()[1] - line.get_xdata()[0], line.get_ydata()[0])
                                  for line in root.lines]})
        return saved(figure, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", inspect)
    styles = [render_karyotype(tmp_path, pair, "svg", colors, scale_mode, track_order=order)
              for order in ("target-query", "query-target")]
    assert observed[0]["species"] == {"Target species": 0.7, "Query species": 0.3}
    assert observed[1]["species"] == {"Target species": 0.3, "Query species": 0.7}
    assert len(observed[0]["ribbons"]) == len(observed[1]["ribbons"]) == 1
    for before, after in zip(observed[0]["ribbons"], observed[1]["ribbons"], strict=True):
        assert before[:, 0] == pytest.approx(after[:, 0])
        assert before[:, 1] == pytest.approx(1 - after[:, 1])
    assert [style["track_gene_counts"] for style in styles] == [[6, 10], [6, 10]]
    assert styles[0]["track_ratios"] == styles[1]["track_ratios"]
    assert styles[1]["track_order"] == ["query", "target"]
    assert styles[1]["display_species_order"] == ["Query_species", "Target_species"]
    assert styles[1]["track_metadata_order"] == ["target", "query"]
    assert styles[0]["scale_value"] == styles[1]["scale_value"]
    if scale_mode == "shared":
        assert len(observed[1]["bars"]) == 1
        assert observed[0]["bars"][0][0] == pytest.approx(observed[1]["bars"][0][0])
        assert observed[1]["bars"][0][1] < 0.3
    else:
        assert [bar[0] for bar in observed[0]["bars"]] == pytest.approx([bar[0] for bar in observed[1]["bars"]])
    assert {name: hashlib.sha256((tmp_path / name).read_bytes()).hexdigest() for name in originals} == originals


def test_connection_legend_bounds_match_native_mcscan_and_liftover():
    from jcvi.compara.synteny import synteny_liftover, synteny_scan

    # Single linkage can span much more than 20 ranks, but each link is bounded
    # inclusively on both axes. Four hits to one subject are not four NR anchors.
    assert len(synteny_scan([(i * 20, i * 20, 1) for i in range(4)], 20, 20, 4)) == 1
    assert not synteny_scan([(i * 21, i * 20, 1) for i in range(4)], 20, 20, 4)
    assert not synteny_scan([(i * 20, 0, 1) for i in range(4)], 20, 20, 4)
    lifted = synteny_liftover([(9, 0, 1), (10, 0, 1), (0, 9, 1), (6, 4, 1)], [(0, 0)], 10)
    assert {tuple(point[:2]) for point, _ in lifted} == {(9, 0), (0, 9)}


def test_ds_cds_mapping_rejects_translation_mismatch_and_ambiguous_alias(tmp_path):
    from workflow.support.pairwise_synteny_ds import load_cds

    (tmp_path / "target.pep").write_text(">gene\nMA\n")
    (tmp_path / "target.id_map.tsv").write_text(
        "original_id\tjcvi_id\tstatus\nSpecies_gene\tgene\tselected\n")
    cds = tmp_path / "cds.fa"
    cds.write_text(">gene\nATGGCTTAA\n")
    assert load_cds(cds, tmp_path, "target", "Species", 1)["gene"] == ("ATGGCT", "MA")
    cds.write_text(">gene\nATGGTT\n")
    with pytest.raises(ValueError, match="translation differs"):
        load_cds(cds, tmp_path, "target", "Species", 1)
    cds.write_text(">gene\nATGGCT\n>Species_gene\nATGGCT\n")
    with pytest.raises(ValueError, match="one matching CDS"):
        load_cds(cds, tmp_path, "target", "Species", 1)
    cds.write_text(">gene\nATGTAAGCT\n")
    with pytest.raises(ValueError, match="internal stop"):
        load_cds(cds, tmp_path, "target", "Species", 1)


@pytest.mark.parametrize("defect", ["duplicate_header", "repeated_original"])
def test_dS_id_map_cannot_silently_assign_another_genes_synonymous_variant(tmp_path, defect):
    from workflow.support.pairwise_synteny_ds import load_cds

    analysis, sources = write_small_ds_analysis(tmp_path)
    # Equal translations cannot reveal an accidental swap of synonymous CDS.
    sources[0].write_text(">gene1\nATGGCT\n>gene2\nATGGCC\n")
    if defect == "duplicate_header":
        text = ("original_id\tjcvi_id\tstatus\toriginal_id\n"
                "gene1\tTarget_species_gene1\tselected\tgene2\n"
                "gene2\tTarget_species_gene2\tselected\tgene2\n")
    else:
        text = ("original_id\tjcvi_id\tstatus\n"
                "gene1\tTarget_species_gene1\tselected\n"
                "gene1\tTarget_species_gene2\tselected\n")
    (analysis / "target.id_map.tsv").write_text(text)
    with pytest.raises(ValueError, match="ID.map|[Dd]uplicate"):
        load_cds(sources[0], analysis, "target", "Target_species", 1)


def test_synteny_only_generates_real_plots_preserves_other_stages_and_reuses_analysis(tmp_path):
    workspace = tmp_path / "workspace"
    write_genome(workspace, "Triphyophyllum_peltatum", "protein", (("Chr1", False),))
    write_genome(workspace, "Ancistrocladus_abbreviatus", "cds", (("Chr2", False), ("Chr10", True)))
    pair_table = workspace / "input/synteny_pairs.tsv"
    pair_table.write_text("analysis_id\ttarget_species\tquery_species\ntriphyophyllum_ancistrocladus\tTriphyophyllum_peltatum\tAncistrocladus_abbreviatus\n")
    for stage in ("species_tree", "orthofinder"):
        path = workspace / "output" / stage
        path.mkdir(parents=True)
        (path / "user-output").write_text("preserve this\n")
    result = run_core(workspace)
    assert result.returncode == 0, result.stdout + result.stderr
    root = workspace / "output/genome_evolution/synteny"
    analysis = root / "analysis/triphyophyllum_ancistrocladus"
    plots = root / "plots/triphyophyllum_ancistrocladus"
    summary = json.loads((analysis / "summary.json").read_text())
    assert summary["syntenic_genes"] == [8, 16]
    assert summary["parameters"]["quota"] is None
    style = json.loads((plots / "karyotype_style.json").read_text())
    assert style["track_order"] == ["target", "query"]
    assert style["connection_criteria"]["seed_cscore"] == summary["parameters"]["cscore"]
    assert style["connection_criteria"]["protein_evalue"] == 1e-5
    assert "Connections: syntenic blocks (MCscan + liftover)" in style["connection_legend_text"]
    assert summary["target"]["source"]["mode"] == "protein"
    assert summary["query"]["source"]["mode"] == "cds"
    with (analysis / "blocks.tsv").open() as handle:
        blocks = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["orientation"] for row in blocks} == {"+", "-"}
    assert {row["query_seqid"] for row in blocks} == {"Chr2", "Chr10"}
    for name in ("dotplot", "karyotype"):
        assert (plots / f"{name}.pdf").read_bytes().startswith(b"%PDF")
        svg = (plots / f"{name}.svg").read_text()
        assert "<svg" in svg
        import xml.etree.ElementTree as ET

        labels = {element.text for element in ET.fromstring(svg).iter() if element.tag.endswith("}text")}
        assert {"Chr1", "Chr2", "Chr10"} <= labels
        with Image.open(plots / f"{name}.png") as image:
            image.verify()
    assert (plots / "seqids").read_text() == "Chr1\nChr2,Chr10\n"
    for stage in ("species_tree", "orthofinder"):
        assert sorted(p.name for p in (workspace / "output" / stage).iterdir()) == ["user-output"]
    analysis_mtime = (analysis / "commands.json").stat().st_mtime_ns
    image_mtime = (plots / "karyotype.png").stat().st_mtime_ns
    result = run_core(workspace)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (analysis / "commands.json").stat().st_mtime_ns == analysis_mtime
    assert (plots / "karyotype.png").stat().st_mtime_ns == image_mtime
    pair_table.write_text("analysis_id\ttarget_species\tquery_species\tquery_seqids\ntriphyophyllum_ancistrocladus\tTriphyophyllum_peltatum\tAncistrocladus_abbreviatus\tChr10,Chr2\n")
    result = run_core(workspace, synteny_plot_only="1", synteny_plot_formats="png", artifact_stale_policy="rebuild")
    assert result.returncode == 0, result.stdout + result.stderr
    assert (analysis / "commands.json").stat().st_mtime_ns == analysis_mtime
    assert (plots / "seqids").read_text() == "Chr1\nChr10,Chr2\n"
    assert (plots / "karyotype.png").stat().st_mtime_ns != image_mtime
    original_inputs = {name: (plots / name).read_bytes() for name in
                       ("target.bed", "query.bed", "seqids", "layout", "colored.simple", "karyotype_colors.json")}
    result = run_core(workspace, synteny_plot_only="1", synteny_plot_formats="png",
                      synteny_karyotype_track_order="query-target", artifact_stale_policy="rebuild")
    assert result.returncode == 0, result.stdout + result.stderr
    assert (analysis / "commands.json").stat().st_mtime_ns == analysis_mtime
    assert {name: (plots / name).read_bytes() for name in original_inputs} == original_inputs
    style = json.loads((plots / "karyotype_style.json").read_text())
    assert style["display_species_order"] == ["Ancistrocladus_abbreviatus", "Triphyophyllum_peltatum"]
    assert style["track_gene_counts"] == [8, 16]
    source = workspace / "input/species_protein/Triphyophyllum_peltatum.protein.fa"
    source.write_text(source.read_text().replace("\nM", "\nA", 1))
    result = run_core(workspace)
    assert result.returncode != 0
    assert (analysis / "commands.json").stat().st_mtime_ns == analysis_mtime
    result = run_core(workspace, synteny_plot_only="1", artifact_stale_policy="rebuild")
    assert result.returncode != 0
    assert "plot-only requires" in result.stderr
    assert (analysis / "commands.json").stat().st_mtime_ns == analysis_mtime
    result = run_core(workspace, synteny_plot_only="1", artifact_stale_policy="reuse")
    assert result.returncode != 0
    assert (analysis / "commands.json").stat().st_mtime_ns == analysis_mtime


def test_plot_only_without_completed_analysis_does_not_publish(tmp_path):
    workspace = tmp_path / "workspace"
    write_genome(workspace, "Target_species", "protein", (("Chr1", False),))
    write_genome(workspace, "Query_species", "protein", (("Chr2", False),))
    (workspace / "input/synteny_pairs.tsv").write_text("analysis_id\ttarget_species\tquery_species\npair\tTarget_species\tQuery_species\n")
    result = run_core(workspace, synteny_plot_only="1")
    assert result.returncode != 0
    assert "plot-only requires" in result.stderr
    assert not (workspace / "output/genome_evolution/synteny/analysis").exists()


def test_sort_redraws_one_track_without_changing_analysis_or_dotplot_inputs(tmp_path):
    workspace = tmp_path / "workspace"
    write_genome(workspace, "Target_species", "protein", (("T1", False), ("T2", False)), (0, 100))
    write_genome(workspace, "Query_species", "protein", (("Q1", False), ("Q2", False)), (100, 0))
    (workspace / "input/synteny_pairs.tsv").write_text("analysis_id\ttarget_species\tquery_species\npair\tTarget_species\tQuery_species\n")
    result = run_core(workspace, synteny_plot_formats="png", synteny_karyotype_sort="none")
    assert result.returncode == 0, result.stdout + result.stderr
    root = workspace / "output/genome_evolution/synteny"
    analysis = root / "analysis/pair"
    plots = root / "plots/pair"
    analysis_mtime = (analysis / "commands.json").stat().st_mtime_ns
    inputs = {name: (plots / name).read_bytes() for name in ("target.bed", "query.bed", "target.query.lifted.anchors")}
    assert (plots / "seqids").read_text() == "T1,T2\nQ1,Q2\n"
    for mode, expected in (("target", "T2,T1\nQ1,Q2\n"), ("query", "T1,T2\nQ2,Q1\n"),
                           ("target_length", "T2,T1\nQ1,Q2\n"), ("query_length", "T1,T2\nQ2,Q1\n")):
        result = run_core(workspace, synteny_plot_only="1", synteny_plot_formats="png",
                          synteny_karyotype_sort=mode, artifact_stale_policy="rebuild")
        assert result.returncode == 0, result.stdout + result.stderr
        assert (analysis / "commands.json").stat().st_mtime_ns == analysis_mtime
        assert (plots / "seqids").read_text() == expected
        assert all((plots / name).read_bytes() == value for name, value in inputs.items())
        ordering = json.loads((plots / "karyotype_order.json").read_text())
        assert ordering["mode"] == mode
        assert ordering["orientation_changed"] is False


def test_normal_genome_evolution_can_opt_in_to_pairwise_stage(tmp_path):
    workspace = tmp_path / "workspace"
    write_genome(workspace, "Target_species", "protein", (("Chr1", False),))
    write_genome(workspace, "Query_species", "protein", (("Chr2", False),))
    (workspace / "input/synteny_pairs.tsv").write_text("analysis_id\ttarget_species\tquery_species\npair\tTarget_species\tQuery_species\n")
    result = _run_core(tmp_path, {"genome_evolution_mode": "all", "run_pairwise_synteny": "1",
                                 "synteny_plot_formats": "png", "run_orthofinder": "0",
                                 # The plot uses real NWKIT transactions, not this
                                 # fixture's species-parser-only Python stub.
                                 "PYTHONPATH": os.environ.get("PYTHONPATH", "")})
    assert result.returncode == 0, result.stdout + result.stderr
    root = workspace / "output/genome_evolution/synteny"
    assert (root / "analysis/pair/summary.json").is_file()
    assert (root / "plots/pair/karyotype.png").is_file()


def test_weighted_layout_matches_real_jcvi_track_coordinates_including_shared_starts(tmp_path):
    import matplotlib.pyplot as plt
    from jcvi.graphics.karyotype import Karyotype
    from kffractbias.io import read_bed

    from workflow.support.pairwise_synteny_layout import offsets, track_geometry

    bed = tmp_path / "genes.bed"
    selected = [f"chr{i}" for i in range(18)]
    bed.write_text("".join(f"{sid}\t0\t10\t{sid}_g2\n{sid}\t0\t20\t{sid}_g10\n" for sid in selected))
    (tmp_path / "seqids").write_text(",".join(selected) + "\n")
    (tmp_path / "layout").write_text(f"0.7,0.12,0.92,0,,Example,top,{bed},top\n")
    widths, ranks, ratio, gap = track_geometry(selected, read_bed(bed))
    starts = offsets(selected, widths, gap)
    figure, axes = plt.subplots(figsize=(20, 8))
    native = Karyotype(axes, str(tmp_path / "seqids"), str(tmp_path / "layout"), plot_label=False).tracks[0]
    assert native.gap == pytest.approx(gap)
    assert native.ratio == pytest.approx(ratio)
    for sid in selected:
        for suffix in ("g2", "g10"):
            gene = f"{sid}_{suffix}"
            assert native.get_coords(gene)[0] == pytest.approx(starts[sid] + ratio * ranks[gene])
    plt.close(figure)


@pytest.mark.parametrize("scale_mode", ["shared", "independent"])
def test_pair_layout_matches_native_jcvi_coordinates_with_unequal_tracks(tmp_path, scale_mode):
    import matplotlib.pyplot as plt
    from jcvi.graphics.karyotype import Layout, Track
    from kffractbias.io import read_bed

    from workflow.support.pairwise_synteny_karyotype import scale_tracks
    from workflow.support.pairwise_synteny_layout import offsets, pair_track_geometry

    selected = [[f"t{i}" for i in range(18)], [f"q{i}" for i in range(3)]]
    for side, seqids in zip(("target", "query"), selected, strict=True):
        (tmp_path / f"{side}.bed").write_text("".join(
            f"{sid}\t0\t10\t{sid}_g2\n{sid}\t0\t20\t{sid}_g10\n" for sid in seqids))
    (tmp_path / "layout").write_text(
        f"0.7,0.12,0.92,0,,Target,top,{tmp_path / 'target.bed'},top\n"
        f"0.3,0.12,0.92,0,,Query,bottom,{tmp_path / 'query.bed'},bottom\n")
    genomes = [read_bed(tmp_path / f"{side}.bed") for side in ("target", "query")]
    geometry = pair_track_geometry(selected, genomes, scale_mode)
    figure, axes = plt.subplots(figsize=(20, 8))
    try:
        layout = Layout(str(tmp_path / "layout"), generank=True, seed=1)
        for entry, seqids in zip(layout, selected, strict=True):
            entry.seqids, entry.rev = seqids, set()
            entry.sizes = {sid: 2 for sid in seqids}
        tracks = [Track(axes, entry, draw=False) for entry in layout]
        scale_tracks(tracks, scale_mode)
        for native, seqids, (widths, ranks, ratio, gap) in zip(tracks, selected, geometry, strict=True):
            starts = offsets(seqids, widths, gap)
            assert native.gap == pytest.approx(gap)
            assert native.ratio == pytest.approx(ratio)
            assert native.xend == pytest.approx(0.12 + sum(widths.values()) + (len(seqids) - 1) * gap)
            for gene, rank in ranks.items():
                sid = gene.rsplit("_", 1)[0]
                assert native.get_coords(gene)[0] == pytest.approx(starts[sid] + ratio * rank)
        if scale_mode == "shared":
            assert tracks[0].ratio == tracks[1].ratio
            assert tracks[1].xend < tracks[0].xend
        else:
            assert tracks[0].xend == tracks[1].xend == pytest.approx(0.92)
        with pytest.raises(ValueError, match="karyotype-scale"):
            scale_tracks(tracks, "bad")
    finally:
        plt.close(figure)


def write_small_ds_analysis(tmp_path):
    analysis = tmp_path / "analysis"
    analysis.mkdir()
    sources = []
    for side, species in (("target", "Target_species"), ("query", "Query_species")):
        sources.append(tmp_path / f"{side}.cds.fa")
        sources[-1].write_text(">gene1\nATGGCT\n>gene2\nATGGCT\n")
        (analysis / f"{side}.pep").write_text(f">{species}_gene1\nMA\n>{species}_gene2\nMA\n")
        (analysis / f"{side}.id_map.tsv").write_text(
            "original_id\tjcvi_id\tstatus\n" + "".join(f"gene{i}\t{species}_gene{i}\tselected\n" for i in (1, 2)))
    (analysis / "target.query.lifted.anchors").write_text(
        "###\nTarget_species_gene1 Query_species_gene1 10\nTarget_species_gene2 Query_species_gene2 10\n")
    return analysis, sources


@pytest.mark.parametrize("collision", ["cds", "symlink", "hardlink", "analysis"])
def test_codon_alignment_output_cannot_overwrite_inputs(tmp_path, collision):
    from workflow.support.pairwise_synteny_ds import prepare_pairs

    analysis, sources = write_small_ds_analysis(tmp_path)
    original = sources[0].read_bytes()
    output = tmp_path / "output.tsv"
    if collision == "cds":
        output = sources[0]
    elif collision == "symlink":
        output.symlink_to(sources[0])
    elif collision == "hardlink":
        os.link(sources[0], output)
    else:
        output = analysis / "new-output.tsv"
    with pytest.raises(ValueError, match="Input and output"):
        prepare_pairs(analysis, sources, ("Target_species", "Query_species"), 1, output)
    assert sources[0].read_bytes() == original


def test_failed_codon_alignment_preserves_existing_output(tmp_path, monkeypatch):
    from workflow.support import pairwise_synteny_ds as ds

    analysis, sources = write_small_ds_analysis(tmp_path)
    output = tmp_path / "pairs.tsv"
    output.write_text("previous complete result\n")
    calls = []
    def fail_second(pair, sequences, code):
        calls.append(pair)
        if len(calls) == 2:
            raise RuntimeError("alignment failed")
        return dict(zip(ds.PAIR_COLUMNS, ("a|b", "a", "b", "ATGGCT", "ATGGCT"), strict=True))
    monkeypatch.setattr(ds, "align_pair", fail_second)
    with pytest.raises(RuntimeError, match="alignment failed"):
        ds.prepare_pairs(analysis, sources, ("Target_species", "Query_species"), 1, output)
    assert output.read_text() == "previous complete result\n"
    assert list(tmp_path.glob(".pairs.tsv.*.tmp")) == []


@pytest.mark.parametrize("stdout,message", [
    (">target\nX\n>target\nX\n>query\nX\n", "Invalid MAFFT"),
    (">target\nA\n>query\nX\n", "changed a protein"),
])
def test_alignment_rejects_duplicate_ids_and_changed_unknown_residues(tmp_path, monkeypatch, stdout, message):
    from types import SimpleNamespace

    from workflow.support import pairwise_synteny_ds as ds

    monkeypatch.setattr(ds.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(stdout=stdout))
    with pytest.raises(ValueError, match=message):
        ds.align_pair(("a", "b"), [{"a": ("NNN", "X")}, {"b": ("NNN", "X")}], 1)


@pytest.mark.parametrize("sequence", ["ATGſCT", "ATG---", "ATG?CT"])
def test_cds_loader_rejects_nonascii_and_aligned_or_invalid_raw_cds(tmp_path, sequence):
    from workflow.support.pairwise_synteny_ds import load_cds

    path = tmp_path / "cds.fa"
    path.write_text(">gene\n" + sequence + "\n")
    with pytest.raises(ValueError, match="alphabet"):
        load_cds(path, tmp_path, "target", "Species", 1)


@pytest.mark.parametrize("code,sequence,protein", [(27, "ATGTGA", "MW"), (28, "ATGTAA", "MQ"), (31, "ATGTAG", "ME")])
def test_dS_loader_preserves_dual_meaning_code_sense_codons(tmp_path, code, sequence, protein):
    from workflow.support.pairwise_synteny_ds import load_cds

    (tmp_path / "target.pep").write_text(">gene\n" + protein + "\n")
    (tmp_path / "target.id_map.tsv").write_text("original_id\tjcvi_id\tstatus\ngene\tgene\tselected\n")
    cds = tmp_path / "cds.fa"
    cds.write_text(">gene\n" + sequence + "\n")
    assert load_cds(cds, tmp_path, "target", "Species", code)["gene"] == (sequence, protein)


@pytest.mark.parametrize("row", ["a|b\t\tok\n", "a|b\tNaN\tok\n", "a|b\t1\tsaturated\n", "a|b\t-1\tok\n", "a|b\t1\n"])
def test_dS_report_rejects_inconsistent_values_and_malformed_rows(tmp_path, row):
    from workflow.support.pairwise_synteny_ds import read_ds

    path = tmp_path / "ds.tsv"
    path.write_text("pair_id\tdS\tstatus\n" + row)
    with pytest.raises(ValueError):
        read_ds(path)


@pytest.mark.parametrize("text", [
    "pair_id\tdS\tstatus\na|b\t\tunknown\n",
    "pair_id\tdS\tstatus\n \t\tsaturated\n",
    "pair_id\tdS\tstatus\t\na|b\t0\tok\t\n",
    'pair_id\tdS\tstatus\n"a|b\t0\tok\n',
])
def test_dS_report_does_not_hide_corrupt_evidence_as_missing(tmp_path, text):
    from workflow.support.pairwise_synteny_ds import read_ds

    path = tmp_path / "ds.tsv"
    path.write_text(text)
    with pytest.raises((ValueError, csv.Error)):
        read_ds(path)


def test_length_source_change_redraws_dS_only_and_failed_filter_preserves_plots(tmp_path):
    workspace = tmp_path / "workspace"
    for species in ("Target_species", "Query_species"):
        for mode in ("protein", "cds"):
            write_genome(workspace, species, mode, (("Chr1", False), ("Chr2", False)), seed_offsets=(0, 1000))
    sizes = workspace / "input/target.sizes"
    sizes.write_text("Chr1\t999999\nChr2\t1000000\n")
    (workspace / "input/synteny_pairs.tsv").write_text(
        "analysis_id\ttarget_species\tquery_species\ttarget_sizes\n"
        "pair\tTarget_species\tQuery_species\tinput/target.sizes\n")
    options = {"synteny_plot_formats": "png", "synteny_dotplot_color": "ds", "synteny_karyotype_sort": "none"}
    result = run_core(workspace, **options)
    assert result.returncode == 0, result.stdout + result.stderr
    root = workspace / "output/genome_evolution/synteny"
    analysis_mtime = (root / "analysis/pair/commands.json").stat().st_mtime_ns
    ds_mtime = (root / "ds/pair/ds.tsv").stat().st_mtime_ns
    plot = root / "plots/pair"
    before = json.loads((plot / "dotplot_ds.json").read_text())
    assert before["anchor_count"] == 8
    assert json.loads((root / "ds/pair/summary.json").read_text())["unique_anchor_pairs"] == 16
    sizes.write_text("Chr1\t1000000\nChr2\t1000000\n")
    result = run_core(workspace, **options, synteny_plot_only="1", artifact_stale_policy="rebuild")
    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads((plot / "dotplot_ds.json").read_text())["anchor_count"] == 16
    assert (root / "analysis/pair/commands.json").stat().st_mtime_ns == analysis_mtime
    assert (root / "ds/pair/ds.tsv").stat().st_mtime_ns == ds_mtime
    previous = (plot / "dotplot.png").read_bytes()
    result = run_core(workspace, **options, synteny_plot_only="1", synteny_dotplot_min_length="3000000", artifact_stale_policy="rebuild")
    assert result.returncode != 0
    assert "No anchors connect" in result.stderr
    assert (plot / "dotplot.png").read_bytes() == previous
    assert (root / "ds/pair/ds.tsv").stat().st_mtime_ns == ds_mtime


@pytest.mark.parametrize("defect", ["wrong_pair", "duplicate_pair", "wrong_code", "wrong_method", "missing_column",
                                    "wrong_schema", "wrong_semantics", "unknown_status"])
def test_dS_stage_checks_report_identity_not_just_row_count(tmp_path, monkeypatch, defect):
    import cdskit.dnds
    from cdskit.codonutil import CODON_SEMANTICS_VERSION
    from cdskit.tsvio import TSV_REPORT_SCHEMA_VERSION

    from workflow.support import pairwise_synteny_ds as ds

    analysis = tmp_path / "output/genome_evolution/synteny/analysis/pair"
    analysis.mkdir(parents=True)
    (analysis / "target.query.lifted.anchors").write_text("###\na b 10\nc d 10\n")
    sources = {side: {"fasta": str(tmp_path / f"{side}.fa"), "genetic_code": 1} for side in ("target", "query")}
    plan = {"workspace": str(tmp_path), "ds_tools": {}, "pairs": [
        {"analysis_id": "pair", "target_species": "Target_species", "query_species": "Query_species", "ds": sources}]}
    def fake_prepare(analysis, paths, species, code, output, cpus):
        output.write_text("audit\n")
        return 2
    def bad_report(args):
        columns = list(cdskit.dnds.REPORT_COLUMNS)
        if defect == "missing_column":
            columns.remove("method")
        with Path(args.outfile).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            for index, pair in enumerate(("a|b", "c|d")):
                row = dict.fromkeys(columns, "")
                row.update(pair_id="unexpected" if defect == "wrong_pair" and index else "a|b" if defect == "duplicate_pair" else pair,
                           dS="" if defect == "unknown_status" else "0",
                           status="unknown" if defect == "unknown_status" else "ok",
                           codon_table="2" if defect == "wrong_code" else "1",
                           schema_version="invalid" if defect == "wrong_schema" else TSV_REPORT_SCHEMA_VERSION,
                           codon_semantics_version="invalid" if defect == "wrong_semantics" else CODON_SEMANTICS_VERSION)
                if "method" in columns:
                    row["method"] = "unknown" if defect == "wrong_method" else cdskit.dnds.METHOD
                writer.writerow(row)
    monkeypatch.setattr(ds, "prepare_pairs", fake_prepare)
    monkeypatch.setattr(cdskit.dnds, "dnds_main", bad_report)
    with pytest.raises(ValueError):
        ds.estimate_ds(plan, tmp_path / "ds", 1)
    assert not (tmp_path / "ds/pair/summary.json").exists()


def test_real_2x2_dotplots_filter_short_scaffolds_and_share_color_mode_order(tmp_path):
    workspace = tmp_path / "workspace"
    chromosomes = tuple((f"Chr{i}", i == 4) for i in range(1, 6))
    for species, offsets in (("Target_species", (0, 1000, 0, 1000, 5000)),
                             ("Query_species", (1000, 0, 1000, 0, 5000))):
        for mode in ("protein", "cds"):
            write_genome(workspace, species, mode, chromosomes, offsets)
        gff = workspace / f"input/species_gff/{species}.gff3"
        gff.write_text(gff.read_text().replace("##sequence-region Chr5 1 2000000", "##sequence-region Chr5 1 999999"))
    (workspace / "input/synteny_pairs.tsv").write_text(
        "analysis_id\ttarget_species\tquery_species\n"
        "pair\tTarget_species\tQuery_species\n")
    options = {"synteny_plot_formats": "png,svg", "synteny_karyotype_sort": "none"}
    result = run_core(workspace, **options)
    assert result.returncode == 0, result.stdout + result.stderr
    root = workspace / "output/genome_evolution/synteny"
    plots = root / "plots/pair"
    metadata = json.loads((plots / "dotplot_order.json").read_text())
    assert metadata["display_order"] == {"target": ["Chr1", "Chr3", "Chr2", "Chr4"],
                                         "query": ["Chr2", "Chr4", "Chr1", "Chr3"]}
    assert metadata["anchor_count"] == 64
    assert metadata["excluded_anchor_count"] == 8
    assert len(metadata["groups"]) == 2
    assert "Chr5" not in (plots / "dotplot.target.bed").read_text()
    assert "Chr5" not in (plots / "dotplot.query.bed").read_text()
    assert "Chr5" in (plots / "target.bed").read_text()
    assert "Chr5" in (plots / "seqids").read_text()
    orientation = (plots / "dotplot.png").read_bytes()
    bed_order = [(plots / f"dotplot.{side}.bed").read_bytes() for side in ("target", "query")]
    result = run_core(workspace, **options, synteny_dotplot_color="ds", artifact_stale_policy="rebuild")
    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads((plots / "dotplot_ds.json").read_text())["anchor_count"] == 64
    assert json.loads((root / "ds/pair/summary.json").read_text())["unique_anchor_pairs"] == 72
    assert [(plots / f"dotplot.{side}.bed").read_bytes() for side in ("target", "query")] == bed_order
    assert (plots / "dotplot.png").read_bytes() != orientation
    with Image.open(plots / "dotplot.png") as image:
        image.verify()


def test_missing_chromosome_lengths_stop_before_dS_alignment(tmp_path):
    workspace = tmp_path / "workspace"
    for species in ("Target_species", "Query_species"):
        for mode in ("protein", "cds"):
            write_genome(workspace, species, mode, (("Chr1", False),))
        gff = workspace / f"input/species_gff/{species}.gff3"
        gff.write_text(gff.read_text().replace("##sequence-region Chr1 1 2000000\n", ""))
    (workspace / "input/synteny_pairs.tsv").write_text(
        "analysis_id\ttarget_species\tquery_species\npair\tTarget_species\tQuery_species\n")
    result = run_core(workspace, synteny_dotplot_color="ds", synteny_plot_formats="png")
    assert result.returncode != 0
    assert "Missing physical chromosome lengths" in result.stderr
    root = workspace / "output/genome_evolution/synteny"
    assert (root / "analysis/pair/summary.json").exists()
    assert not (root / "ds").exists()
    assert not (root / "plots").exists()


def test_no_eligible_anchors_stop_before_dS_alignment(tmp_path):
    workspace = tmp_path / "workspace"
    for species in ("Target_species", "Query_species"):
        for mode in ("protein", "cds"):
            write_genome(workspace, species, mode, (("Chr1", False),))
    (workspace / "input/synteny_pairs.tsv").write_text(
        "analysis_id\ttarget_species\tquery_species\npair\tTarget_species\tQuery_species\n")
    result = run_core(workspace, synteny_dotplot_color="ds", synteny_plot_formats="png",
                      synteny_dotplot_min_length="3000000")
    assert result.returncode != 0
    assert "No anchors connect" in result.stderr
    root = workspace / "output/genome_evolution/synteny"
    assert (root / "analysis/pair/summary.json").exists()
    assert not (root / "ds").exists()
    assert not (root / "plots").exists()


def test_dS_workers_do_not_multiply_the_allocated_cpu_budget_with_blas_threads(tmp_path):
    workspace = tmp_path / "workspace"
    for species in ("Target_species", "Query_species"):
        for mode in ("protein", "cds"):
            write_genome(workspace, species, mode, (("Chr1", False),))
    (workspace / "input/synteny_pairs.tsv").write_text(
        "analysis_id\ttarget_species\tquery_species\npair\tTarget_species\tQuery_species\n")
    report = tmp_path / "thread_environment.tsv"
    startup = tmp_path / "sitecustomize.py"
    startup.write_text('''import os
import sys
from pathlib import Path

if sys.argv[0].endswith("pairwise_synteny.py") and sys.argv[1:2] == ["ds"]:
    names = ("GG_TASK_CPUS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")
    Path(os.environ["GG_TEST_THREAD_REPORT"]).write_text("\\t".join(os.environ.get(name, "") for name in names) + "\\n")
''')
    result = run_core(workspace, synteny_dotplot_color="ds", synteny_plot_formats="png",
                      GG_TASK_CPUS="2", OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2",
                      MKL_NUM_THREADS="2", NUMEXPR_NUM_THREADS="2",
                      PYTHONPATH=str(tmp_path), GG_TEST_THREAD_REPORT=str(report))
    assert result.returncode == 0, result.stdout + result.stderr
    assert report.read_text().strip().split("\t") == ["2", "1", "1", "1", "1"]

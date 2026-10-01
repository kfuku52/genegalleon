"""Shared, restrained typography and color encodings for synteny figures."""

SOFT_COLORS = ("#88afc4", "#b1a6c6", "#8ab7ad", "#ceaaaa", "#d0c097",
               "#a9b2bc", "#d4b39e", "#bba5b7", "#b0bb8b", "#aec5b7")
# Higher dS is darker; even the lightest point has >= 3:1 contrast on white.
DS_COLORS = ("#6f98b2", "#5883a4", "#416e95", "#2e5881", "#20446c")
MISSING_COLOR = "#b0b0b0"
STYLE = {"pdf.use14corefonts": True, "text.usetex": False, "font.family": "Helvetica",
         "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
         "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
         "text.color": "black", "axes.labelcolor": "black", "axes.titlecolor": "black",
         "xtick.color": "black", "ytick.color": "black", "svg.fonttype": "none"}


def style_text(fig):
    from matplotlib.text import Text

    for text in fig.findobj(Text):
        text.set_fontfamily("Helvetica")
        text.set_fontsize(8)
        text.set_fontweight("normal")
        text.set_fontstyle("normal")
        text.set_color("black")
        text.set_usetex(False)
        text.set_parse_math(False)


def ds_colormap():
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("genegalleon_ds", DS_COLORS)
    cmap.set_bad(MISSING_COLOR)
    return cmap

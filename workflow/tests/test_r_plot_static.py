from pathlib import Path
import re


REPO_ROOT = Path(__file__).resolve().parents[2]
SUPPORT_DIR = REPO_ROOT / "workflow" / "support"
NUMERIC_LEGEND_POSITION_RE = re.compile(r"legend\.position\s*=\s*c\(")


def test_r_scripts_do_not_use_deprecated_numeric_legend_position():
    offending = []
    for path in SUPPORT_DIR.rglob("*"):
        if path.suffix not in {".r", ".R"}:
            continue
        text = path.read_text(encoding="utf-8")
        if NUMERIC_LEGEND_POSITION_RE.search(text):
            offending.append(path.relative_to(REPO_ROOT).as_posix())

    assert offending == [], (
        "Found deprecated numeric legend.position usage in: "
        + ", ".join(offending)
    )


def test_amino_acid_site_panel_noops_on_empty_site_list():
    stat_branch2tree = (SUPPORT_DIR / "stat_branch2tree_plot.r").read_text(encoding="utf-8")
    fimo_motif = (SUPPORT_DIR / "treevis" / "R" / "05_fimo_motif.R").read_text(encoding="utf-8")

    assert "parse_amino_acid_site_list = function(site_text)" in stat_branch2tree
    assert "return(integer(0))" in stat_branch2tree
    assert "if (length(selected_amino_acid_sites) == 0)" in stat_branch2tree
    assert "amino_acid_site column will not be added" in stat_branch2tree
    assert "return(g)" in fimo_motif


def test_query2family_busco_plot_uses_shared_species_label_parser():
    text = (SUPPORT_DIR / "plot_query2family_presence_absence.R").read_text(encoding="utf-8")

    assert 'source(file.path(script_dir, "species_label_utils.r"), local = TRUE)' in text
    assert "gg_species_label_from_filename(x)" in text
    assert 'sub("_sp_[^._-]+$", "_sp", x)' not in text

"""The default family figure exposes the independent localization probability."""
from pathlib import Path


def test_default_family_plot_includes_peroxisome_after_targeting():
    workflow = Path(__file__).resolve().parents[1]
    core = (workflow / 'core/gg_gene_evolution_core.sh').read_text()
    targeting = core.index('"--panel${panel_index}=signal_peptide"')
    peroxisome = core.index('"--panel${panel_index}=peroxisome"')
    membrane = core.index('"--panel${panel_index}=transmembrane_domain"')
    assert targeting < peroxisome < membrane
    assert '--cdskit_localize "${file_og_cdskit_localize}"' in core
    driver = (workflow / 'support/stat_branch2tree_plot.r').read_text()
    assert "g = add_peroxisome_column(g, args)" in driver

"""The default family figure exposes the independent localization probability."""
from pathlib import Path


def test_default_family_plot_combines_localization():
    workflow = Path(__file__).resolve().parents[1]
    core = (workflow / 'core/gg_gene_evolution_core.sh').read_text()
    localization = core.index('"--panel${panel_index}=localization"')
    membrane = core.index('"--panel${panel_index}=transmembrane_domain"')
    assert localization < membrane
    assert '"--panel${panel_index}=signal_peptide"' not in core
    assert '"--panel${panel_index}=peroxisome"' not in core
    assert '--cdskit_localize "${file_og_cdskit_localize}"' in core
    driver = (workflow / 'support/stat_branch2tree_plot.r').read_text()
    assert "g = add_localization_column(g, args)" in driver

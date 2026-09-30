import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "support/orthogroup_method_comparison.py"
spec = importlib.util.spec_from_file_location("orthogroup_method_comparison", SCRIPT)
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


def test_native_flat_memberships_match_legacy_totals(tmp_path):
    native = tmp_path / "Orthogroups.txt"
    native.write_text("OG0000000: species_a_g1 species_a_g2 species_b_g1\n\nOG0000001: species_b_g2\n")
    legacy = tmp_path / "Orthogroups.GeneCount.tsv"
    legacy.write_text("Orthogroup\tspecies_a\tspecies_b\tTotal\nOG0000000\t2\t1\t3\nOG0000001\t0\t1\t1\n")
    actual = comparison.read_gene_counts(native)
    expected = comparison.read_gene_counts(legacy)[["Orthogroup", "Total"]]
    pandas.testing.assert_frame_equal(actual, expected)


@pytest.mark.parametrize("text", [
    "", "OG0000000:\n", "HOG0000000: a\n", "OG0000000 a\n",
    "OG0000000: a a\n", "OG0000000: a\nOG0000000: b\n",
    "OG0000000: a\nOG0000001: a\n",
])
def test_invalid_memberships_fail_closed(tmp_path, text):
    path = tmp_path / "Orthogroups.txt"
    path.write_text(text)
    with pytest.raises(ValueError):
        comparison.read_gene_counts(path)


def test_core_binds_the_actual_comparison_input():
    source = (SCRIPT.parents[1] / "core/gg_genome_evolution_core.sh").read_text()
    block = source.split('task="Orthogroup method comparison"', 1)[1].split('task="', 1)[0]
    assert '--input "orthogroup_counts=${file_orthogroup_comparison_input}"' in block
    assert '--orthofinder_og_genecount "${file_orthogroup_comparison_input}"' in block


def test_native_comparison_generates_both_plots(tmp_path):
    native = tmp_path / "Orthogroups.txt"
    native.write_text("OG0000000: a b c\nOG0000001: d e\n")
    hog = tmp_path / "HOG.GeneCount.tsv"
    hog.write_text("Orthogroup\tTotal\nHOG0000000\t2\nHOG0000001\t3\n")
    result = subprocess.run([
        sys.executable, str(SCRIPT), "--orthofinder_og_genecount", str(native),
        "--orthofinder_hog_genecount", str(hog),
    ], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "orthogroup_histogram.pdf").read_bytes().startswith(b"%PDF")
    assert "<svg" in (tmp_path / "orthogroup_histogram.svg").read_text()

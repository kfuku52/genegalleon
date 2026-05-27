from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace

import pandas
import pytest


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "single_copy_ortholog_decay_plot.py"


def load_module():
    spec = spec_from_file_location("single_copy_ortholog_decay_plot", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_gene_count(path: Path):
    pandas.DataFrame(
        [
            {"Orthogroup": "OG1", "spA": 1, "spB": 1, "spC": 1, "Total": 3},
            {"Orthogroup": "OG2", "spA": 1, "spB": 2, "spC": 1, "Total": 4},
            {"Orthogroup": "OG3", "spA": 1, "spB": 0, "spC": 1, "Total": 2},
            {"Orthogroup": "OG4", "spA": 0, "spB": 0, "spC": 2, "Total": 2},
        ]
    ).to_csv(path, sep="\t", index=False)


def test_species_count_parser_accepts_auto_lists_and_ranges():
    mod = load_module()

    assert mod.parse_species_counts("auto", 4) == [1, 2, 3, 4]
    assert mod.parse_species_counts("1,3-4,6-8:2", 8) == [1, 3, 4, 6, 8]

    with pytest.raises(ValueError, match="between 1 and 4"):
        mod.parse_species_counts("1,5", 4)


def test_load_gene_count_table_ignores_metadata_columns(tmp_path):
    mod = load_module()
    path = tmp_path / "Orthogroups.GeneCount.tsv"
    pandas.DataFrame(
        [
            {
                "Orthogroup": "OG1",
                "spA": 1,
                "spB": 2,
                "Total": 3,
                "geneid_0.5": "g1",
                "besthit_0.5": "hit",
            }
        ]
    ).to_csv(path, sep="\t", index=False)

    species_cols, counts = mod.load_gene_count_table(path)

    assert species_cols == ["spA", "spB"]
    assert counts.tolist() == [[1, 2]]


def test_decay_metrics_are_nested_for_all_species_count(tmp_path):
    mod = load_module()
    path = tmp_path / "Orthogroups.GeneCount.tsv"
    _write_gene_count(path)
    species_cols, counts = mod.load_gene_count_table(path)
    species_counts = [len(species_cols)]

    values = mod.calculate_decay(counts, species_counts, replicates=5, seed=7)
    summary = mod.summarize_decay(species_counts, values)
    observed = {
        row["metric"]: (row["mean"], row["sd"])
        for row in summary.to_dict("records")
    }

    assert observed["strict_single_copy"] == (1.0, 0.0)
    assert observed["non_missing"] == (2.0, 0.0)
    assert observed["all_observed"] == (4.0, 0.0)


def test_run_writes_summary_and_plot(tmp_path):
    mod = load_module()
    path = tmp_path / "Orthogroups.GeneCount.tsv"
    outdir = tmp_path / "decay"
    _write_gene_count(path)

    mod.run(
        SimpleNamespace(
            orthogroup_genecount=str(path),
            outdir=str(outdir),
            replicates=3,
            species_counts="3",
            seed=1,
            plot_basename="decay",
            summary_name="summary.tsv",
            formats="svg",
        )
    )

    summary = pandas.read_csv(outdir / "summary.tsv", sep="\t")
    assert set(summary["metric"]) == {
        "strict_single_copy",
        "non_missing",
        "all_observed",
    }
    assert (outdir / "decay.svg").exists()


def test_plot_uses_zero_baseline_integer_x_ticks_and_reversed_legend():
    text = SCRIPT_PATH.read_text(encoding="utf-8")

    assert 'ax.set_yscale("log", base=2)' not in text
    assert "ax.set_ylim(0, y_upper)" in text
    assert "MaxNLocator(integer=True)" in text
    assert "ax.legend(handles[::-1], labels[::-1]" in text

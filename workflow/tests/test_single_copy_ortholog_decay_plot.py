import re
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


def _write_selected_gene_count(path: Path):
    pandas.DataFrame(
        [
            {"Orthogroup": "OG1", "spA": 1, "spB": 1, "spC": 1, "Total": 3},
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
    selected_path = tmp_path / "Orthogroups.GeneCount.selected.tsv"
    _write_gene_count(path)
    _write_selected_gene_count(selected_path)
    species_cols, counts = mod.load_gene_count_table(path)
    selected_species_cols, selected_counts = mod.load_gene_count_table(selected_path)
    assert species_cols == selected_species_cols
    species_counts = [len(species_cols)]

    values = mod.calculate_decay(
        counts,
        species_counts,
        replicates=5,
        seed=7,
        selected_counts=selected_counts,
    )
    summary = mod.summarize_decay(species_counts, values)
    observed = {
        row["metric"]: (row["mean"], row["sd"])
        for row in summary.to_dict("records")
    }

    assert observed["strict_single_copy"] == (1.0, 0.0)
    assert observed["non_missing"] == (2.0, 0.0)
    assert observed["selected_observed"] == (2.0, 0.0)
    assert observed["all_observed"] == (4.0, 0.0)


def test_run_writes_summary_and_plot(tmp_path):
    mod = load_module()
    path = tmp_path / "Orthogroups.GeneCount.tsv"
    selected_path = tmp_path / "Orthogroups.GeneCount.selected.tsv"
    outdir = tmp_path / "decay"
    _write_gene_count(path)
    _write_selected_gene_count(selected_path)

    mod.run(
        SimpleNamespace(
            orthogroup_genecount=str(path),
            selected_orthogroup_genecount=str(selected_path),
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
        "selected_observed",
        "all_observed",
    }
    svg_text = (outdir / "decay.svg").read_text(encoding="utf-8")
    assert 'width="259.2pt" height="316.8pt"' in svg_text
    assert "Helvetica" in svg_text
    assert "Selected orthogroups" in svg_text
    assert svg_text.count('id="FillBetweenPolyCollection_') == 4
    assert set(re.findall(r"font-size: ([0-9.]+)px", svg_text)) == {"8"}


def test_run_rejects_mismatched_selected_species_columns(tmp_path):
    mod = load_module()
    path = tmp_path / "Orthogroups.GeneCount.tsv"
    selected_path = tmp_path / "Orthogroups.GeneCount.selected.tsv"
    _write_gene_count(path)
    pandas.DataFrame(
        [{"Orthogroup": "OG1", "spA": 1, "spC": 1, "spB": 1, "Total": 3}]
    ).to_csv(selected_path, sep="\t", index=False)

    with pytest.raises(ValueError, match="identical species columns in the same order"):
        mod.run(
            SimpleNamespace(
                orthogroup_genecount=str(path),
                selected_orthogroup_genecount=str(selected_path),
                outdir=str(tmp_path / "out"),
                replicates=1,
                species_counts="auto",
                seed=1,
                plot_basename="decay",
                summary_name="summary.tsv",
                formats="svg",
            )
        )


def test_truncate_at_mean_floor_interpolates_crossing_on_log_scale():
    mod = load_module()
    x = mod.numpy.array([34.0, 35.0, 36.0])
    mean = mod.numpy.array([1.454, 1.101, 0.932])
    sd = mod.numpy.array([0.20, 0.15, 0.10])

    x_out, mean_out, sd_out = mod.truncate_at_mean_floor(x, mean, sd, 1.0)

    assert 35.0 < x_out[-1] < 36.0
    assert mean_out[-1] == 1.0
    assert len(sd_out) == len(x_out)


def test_plot_uses_requested_log_scale_dimensions_and_type():
    text = SCRIPT_PATH.read_text(encoding="utf-8")

    assert 'ax.set_yscale("log", base=10)' in text
    assert "figsize=(3.6, 4.4)" in text
    assert '"font.family": "Helvetica"' in text
    assert '"font.size": 8' in text
    assert "MaxNLocator(integer=True)" in text
    assert "truncate_at_mean_floor(x, mean, sd, log_floor)" in text
    assert '"{:,.0f}".format(value)' in text

from types import SimpleNamespace

import matplotlib
import matplotlib.pyplot as plt
import numpy
import pytest
from matplotlib.figure import Figure

from workflow.support import scaffold_size_histogram


@pytest.mark.parametrize("lengths,large,small,bars", [
    ([50, 100, 200, 99], [1, 2], [0, 3], [100, 200, 149]),
    ([50, 70], [], [0, 1], [120]),
    ([100, 200], [0, 1], [], [100, 200]),
])
def test_scaffold_plot_preserves_size_partition_jitter_and_coordinates(
        tmp_path, monkeypatch, lengths, large, small, bars):
    path = tmp_path / "scaffolds.tsv"
    path.write_text("#id\tlength\tGC\tGC-Skew\n" + "".join(
        f"s{i}\t{length}\t{40 + i}\t{(i - 1) / 10}\n"
        for i, length in enumerate(lengths)))
    observed = []

    def inspect(figure, output, **kwargs):
        assert [bar.get_height() for bar in figure.axes[0].patches] == bars
        expected_rng = numpy.random.default_rng(11)
        values = (lengths, [40 + i for i in range(len(lengths))],
                  [(i - 1) / 10 for i in range(len(lengths))])
        for axis, column in zip(figure.axes[1:4], values, strict=True):
            groups = [(positions, offset) for positions, offset in ((large, -0.035), (small, 0.035))
                      if positions]
            assert len(axis.collections) == len(groups)
            for collection, (positions, offset) in zip(axis.collections, groups, strict=True):
                expected_x = offset + expected_rng.normal(0, 0.01, len(positions))
                numpy.testing.assert_allclose(collection.get_offsets()[:, 0], expected_x)
                numpy.testing.assert_allclose(collection.get_offsets()[:, 1], [column[i] for i in positions])
        groups = [positions for positions in (large, small) if positions]
        assert len(figure.axes[4].collections) == len(groups)
        for collection, positions in zip(figure.axes[4].collections, groups, strict=True):
            numpy.testing.assert_allclose(collection.get_offsets()[:, 0],
                                          numpy.log10([lengths[i] for i in positions]))
            numpy.testing.assert_allclose(collection.get_offsets()[:, 1], [40 + i for i in positions])
        observed.append(kwargs["format"])

    monkeypatch.setattr(Figure, "savefig", inspect)
    monkeypatch.setattr(scaffold_size_histogram, "get_matplotlib", lambda: matplotlib)
    with matplotlib.rc_context():
        scaffold_size_histogram.run(SimpleNamespace(fx2tab=path, min_scaffold_size=100))
    plt.close("all")
    assert observed == ["svg", "pdf"]

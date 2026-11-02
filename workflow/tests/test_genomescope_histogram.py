import argparse
import importlib.util
import json
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "support/genomescope_histogram.py"
spec = importlib.util.spec_from_file_location("genomescope_histogram", SCRIPT)
histogram = importlib.util.module_from_spec(spec)
spec.loader.exec_module(histogram)


def test_sparse_histogram_preserves_tail_mass_and_last_observation(tmp_path):
    source, output, metadata = [tmp_path / name for name in ("raw", "model", "metadata")]
    source.write_text("1\t5\n3\t4\n8\t2\n9\t0\n10\t0\n")
    result = histogram.prepare(source, output, metadata, 8)
    rows = [tuple(map(int, line.split())) for line in output.read_text().splitlines()]
    assert rows == [(1, 5), (2, 0), (3, 4), (4, 0), (5, 0), (6, 0), (7, 0), (8, 2), (9, 0)]
    assert sum(x * y for x, y in rows[:-1]) == result["total_kmer_observations"] == 33
    assert result["ceiling_bin_distinct_kmers"] == 2
    assert json.loads(metadata.read_text()) == result


@pytest.mark.parametrize("text", ["1 2\n1 3\n", "0 2\n", "1 -2\n", "1 2 3\n", "1 0\n2 0\n"])
def test_rejects_invalid_or_empty_spectrum(tmp_path, text):
    source = tmp_path / "raw"
    source.write_text(text)
    with pytest.raises(ValueError):
        histogram.prepare(source, tmp_path / "model", tmp_path / "metadata", 1000)


@pytest.mark.parametrize(("name", "value"), [("ploidy", 7), ("lower", 0), ("upper", 0),
                                            ("initial_coverage", 0), ("max_kmercov", -2)])
def test_invalid_settings_fail_before_counting(name, value):
    settings = dict(kmer_length=21, lower=1, upper=10000000, ploidy=2,
                    initial_coverage=-1, max_kmercov=-1)
    settings[name] = value
    with pytest.raises(ValueError):
        histogram.validate(argparse.Namespace(**settings))

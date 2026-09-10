import importlib.util
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))
spec = importlib.util.spec_from_file_location("busco_quality_metadata", SUPPORT / "busco_quality_metadata.py")
quality = importlib.util.module_from_spec(spec)
spec.loader.exec_module(quality)


def summary(path, value="C:95.0%[S:80.0%,D:15.0%],F:3.0%,M:2.0%,n:1000"):
    path.write_text("# BUSCO version is: 6.0.0\n# The lineage dataset is: eukaryota_odb12 (Creation date: today)\n"
                    "# BUSCO was run in mode: proteins\n" + value + "\n")
    return path


def test_completeness_includes_duplicates_and_preserves_species(tmp_path):
    path = summary(tmp_path / "Species_with_strain_A.busco.short.txt")
    row = quality.read_short_summary(path)
    assert row["busco_complete_pct"] == 95
    assert row["busco_duplicated_pct"] == 15
    assert row["species"] == "Species_with_strain_A"
    assert row["lineage"] == "eukaryota_odb12"
    assert row["mode"] == "proteins"
    assert row["source"] == str(path)


def test_missing_is_not_zero(tmp_path):
    assert quality.collect(tmp_path / "absent") == []


@pytest.mark.parametrize("value", ["C:101%[S:99%,D:2%]", "C:50%[S:40%,D:60%]", "not a BUSCO report", "C:90.0%[S:500.0%,D:10.0%]", "C:90.0%[S:40.0%,D:10.0%]"])
def test_invalid_summary_rejected(tmp_path, value):
    with pytest.raises(ValueError):
        quality.read_short_summary(summary(tmp_path / "sp1.busco.short.txt", value))


def test_ambiguous_species_rejected(tmp_path):
    summary(tmp_path / "sp1.busco.short.txt")
    summary(tmp_path / "sp1_busco.short.txt")
    with pytest.raises(ValueError, match="Ambiguous"):
        quality.collect(tmp_path)


def test_rounded_percentages_and_output_input_protection(tmp_path, monkeypatch):
    path = summary(tmp_path / "sp1.busco.short.txt", "C:95.0%[S:80.0%,D:15.1%]")
    assert quality.read_short_summary(path)["busco_complete_pct"] == 95
    before = path.read_bytes()
    monkeypatch.setattr(sys, "argv", ["busco_quality_metadata", "--directory", str(tmp_path), "--output", str(path)])
    with pytest.raises(ValueError, match="must not replace"):
        quality.main()
    assert path.read_bytes() == before
    with pytest.raises(ValueError, match="not a directory"):
        quality.collect(path)

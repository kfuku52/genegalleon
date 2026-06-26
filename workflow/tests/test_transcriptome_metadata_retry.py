import csv
import sys
from pathlib import Path

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
if str(SUPPORT_DIR) not in sys.path:
    sys.path.insert(0, str(SUPPORT_DIR))

import amalgkit_metadata_accessions as ama  # noqa: E402


def test_merge_metadata_tables_by_run_appends_new_columns_from_relaxed_metadata(tmp_path):
    primary_path = tmp_path / "metadata.tsv"
    extra_path = tmp_path / "metadata.relaxed.tsv"
    output_path = tmp_path / "metadata.merged.tsv"

    primary_path.write_text(
        "\n".join(
            [
                "run\tscientific_name\tlib_layout",
                "SRR1\tSpecies one\tsingle",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    extra_path.write_text(
        "\n".join(
            [
                "run\tscientific_name\tlib_layout\tsample_attribute_tissue\tdev_stage",
                "SRR1\tSpecies one\tsingle\tleaf\tadult",
                "SRR2\tSpecies two\tpaired\troot\tjuvenile",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    ama.merge_metadata_tables_by_run(primary_path, extra_path, output_path)

    with output_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    assert reader.fieldnames == [
        "run",
        "scientific_name",
        "lib_layout",
        "sample_attribute_tissue",
        "dev_stage",
    ]
    assert [row["run"] for row in rows] == ["SRR1", "SRR2"]
    assert rows[0]["scientific_name"] == "Species one"
    assert rows[0]["sample_attribute_tissue"] == ""
    assert rows[0]["dev_stage"] == ""
    assert rows[1]["scientific_name"] == "Species two"
    assert rows[1]["sample_attribute_tissue"] == "root"
    assert rows[1]["dev_stage"] == "juvenile"

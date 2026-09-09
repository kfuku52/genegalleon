import importlib.util
from pathlib import Path

import pandas as pd
import pytest

from workflow.support.score_hgt_candidates import intron_observation_masks, intron_support_from_leaf_rows

spec = importlib.util.spec_from_file_location(
    "asr_intron_evolution", Path(__file__).resolve().parents[1] / "support" / "asr_intron_evolution.py"
)
asr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(asr)


def test_gff_counts_preserve_gene_identifiers_and_missing_observations(tmp_path):
    path = tmp_path / "gff.tsv"
    path.write_text("gene_id\tnum_intron\nNA\t2\n001\t0\nmissing\tNA\n")
    counts = asr.read_intron_counts(path)
    assert counts["NA"] == 2
    assert counts["001"] == 0
    assert pd.isna(counts["missing"])


@pytest.mark.parametrize("rows", ["A\t-1\n", "A\t1.5\n", "A\tinf\n", "A\tnone\n", "A\t1\nA\t2\n", "\t1\n"])
def test_gff_counts_reject_invalid_measurements_and_ambiguous_ids(tmp_path, rows):
    path = tmp_path / "gff.tsv"
    path.write_text("gene_id\tnum_intron\n" + rows)
    with pytest.raises(ValueError):
        asr.read_intron_counts(path)


def test_hgt_intron_support_does_not_count_imputed_probabilities_as_measurements():
    leaves = pd.DataFrame({
        "num_intron": [2, 0, None],
        "intron_present": [1, 0, 0.904],
        "intron_is_imputed": [False, False, True],
    })
    observed, supported = intron_observation_masks(leaves)
    assert observed.tolist() == [True, True, False]
    assert supported.tolist() == [True, False, False]
    assert intron_support_from_leaf_rows(leaves) == {
        "measured_count": 2, "supported_count": 1, "support_fraction": 0.5,
    }


def test_hgt_intron_support_respects_imputed_flag_without_count_column():
    leaves = pd.DataFrame({"intron_present": [1, 0.904], "intron_is_imputed": [False, True]})
    assert intron_observation_masks(leaves)[1].tolist() == [True, False]


@pytest.mark.parametrize("name", ["run_scm_intron", "GG_GENE_EVOLUTION_RUN_SCM_INTRON"])
def test_retired_scm_override_fails_with_migration_instruction(name):
    import os
    import subprocess

    root = Path(__file__).resolve().parents[2]
    result = subprocess.run(["bash", str(root / "workflow" / "gg_gene_evolution_entrypoint.sh")],
                            cwd=root, env={**os.environ, name: "1"}, capture_output=True, text=True)
    assert result.returncode != 0
    assert "run_scm_intron is retired; use run_asr_intron" in result.stderr

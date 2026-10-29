"""A saved fit needs scientific identity checks before generic legacy adoption."""
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))
spec = importlib.util.spec_from_file_location("validate_cafe_outputs", SUPPORT / "validate_cafe_outputs.py")
validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)


def fixture(tmp_path):
    copies = tmp_path / "copies.tsv"
    copies.write_text("Desc\tFamilyID\tA\tB\tC\n(null)\tOG1\t2\t4\t1\n(null)\tOG2\t0\t0\t1\n")
    tree = tmp_path / "tree.nwk"
    tree.write_text("((A:1,B:1):1,C:2);\n")
    out = tmp_path / "native"
    out.mkdir()
    (out / "Gamma_asr.tre").write_text("#NEXUS\nBEGIN TREES;\nTREE OG1 = ((A<1>_2:1,B<2>_4*:1)<4>_6:1,C<3>_1:2)<5>_7;\nEND;\n")
    header = "#FamilyID\tA<1>\tB<2>\tC<3>\t<4>\t<5>\n"
    (out / "Gamma_count.tab").write_text(header + "OG1\t2\t4\t1\t6\t7\n")
    (out / "Gamma_change.tab").write_text(header + "OG1\t-4\t-2\t-6\t-1\t0\n")
    (out / "Gamma_branch_probabilities.tab").write_text(header + "OG1\t.1\t.01\t.1\t.1\tN/A\n")
    (out / "Gamma_family_results.txt").write_text("#FamilyID\tpvalue\tSignificant at 0.05\nOG1\t.01\ty\n")
    return copies, tree, out


def test_matching_fit_allows_documented_input_family_omissions(tmp_path):
    result = validator.validate(*fixture(tmp_path))
    assert result["validated"] and result["species"] == 3
    assert result["reconstructed_families"] == 1
    assert result["input_families_without_reconstruction"] == 1


@pytest.mark.parametrize("changed", ["input_species", "topology", "extant_count", "foreign_family", "duplicate_family", "duplicate_header", "branch_header", "foreign_probability", "foreign_tree_family", "invalid_change", "fractional_change"])
def test_incompatible_fit_is_rejected_without_modifying_native_files(tmp_path, changed):
    copies, tree, out = fixture(tmp_path)
    if changed == "input_species":
        copies.write_text(copies.read_text().replace("\tC\n", "\tD\n"))
    elif changed == "topology":
        tree.write_text("((A:1,C:1):1,B:2);\n")
    elif changed == "extant_count":
        copies.write_text(copies.read_text().replace("OG1\t2\t4", "OG1\t3\t4"))
    elif changed == "foreign_family":
        for name in ["Gamma_count.tab", "Gamma_change.tab"]:
            path = out / name
            path.write_text(path.read_text().replace("OG1", "foreign"))
    elif changed == "duplicate_family":
        path = out / "Gamma_count.tab"
        path.write_text(path.read_text() + path.read_text().splitlines(True)[1])
    elif changed == "duplicate_header":
        copies.write_text(copies.read_text().replace("\tB\tC", "\tA\tC"))
    elif changed == "branch_header":
        path = out / "Gamma_change.tab"
        path.write_text(path.read_text().replace("<4>", "<99>"))
    elif changed in {"invalid_change", "fractional_change"}:
        path = out / "Gamma_change.tab"
        path.write_text(path.read_text().replace("-4", "N/A" if changed == "invalid_change" else "0.5"))
    elif changed == "foreign_tree_family":
        path = out / "Gamma_asr.tre"
        path.write_text(path.read_text().replace("TREE OG1 =", "TREE foreign ="))
    else:
        path = out / "Gamma_branch_probabilities.tab"
        path.write_text(path.read_text().replace("OG1", "foreign"))
    before = {p: p.read_bytes() for p in out.iterdir()}
    with pytest.raises(ValueError):
        validator.validate(copies, tree, out)
    assert all(p.read_bytes() == data for p, data in before.items())


@pytest.mark.parametrize("value", ["nan", "inf", "0.5", "-1"])
def test_count_domain_is_finite_nonnegative_integral(value):
    with pytest.raises(ValueError):
        validator.count(value)


def test_count_identity_is_not_rounded_through_binary_float():
    assert validator.count("9007199254740993") == 9007199254740993


@pytest.mark.parametrize("active", ["0", "1"])
def test_raw_validation_runs_before_cached_matrix_adoption(tmp_path, active):
    raw = tmp_path / "counts.tsv"
    raw.write_text("besthit_0.95\tOrthogroup\tA\tA\tB\nY\tOG1\t1\t9\t2\n")
    tree = tmp_path / "tree.nwk"
    tree.write_text("(A:1,B:1);\n")
    cache = tmp_path / "already-prepared.tsv"
    cache.write_text("existing scientific matrix must remain unchanged\n")
    core = (SUPPORT.parent / "core/gg_genome_evolution_core.sh").read_text()
    block = core.split('task="Orthogroup copy-number matrix preparation"', 1)[1]
    block = 'if [[ ${run_orthogroup_copy_number_stage} -eq 1 ]]; then' + block.split('if [[ ${run_orthogroup_copy_number_stage} -eq 1 ]]; then', 1)[1].split('copy_number_needs_update=0', 1)[0]
    env = {**os.environ, "run_orthogroup_copy_number_stage": active, "gg_support_dir": str(SUPPORT),
           "file_orthogroup_genecount_selected": str(raw), "file_dated_species_tree": str(tree)}
    result = subprocess.run(["bash", "-c", "set -euo pipefail\n" + block], env=env, capture_output=True, text=True, timeout=30)
    assert (result.returncode == 0) == (active == "0"), result.stdout + result.stderr
    assert cache.read_text() == "existing scientific matrix must remain unchanged\n"


def test_ambiguous_family_id_aliases_are_rejected(tmp_path):
    copies, tree, out = fixture(tmp_path)
    copies.write_text("Desc\tFamilyID\t#FamilyID\tA\tB\tC\n(null)\tOG1\tOG1\t2\t4\t1\n")
    with pytest.raises(ValueError, match="ambiguous"):
        validator.validate(copies, tree, out)


def test_real_native_gamma_fit_satisfies_semantic_guard(tmp_path):
    copies, tree, _ = fixture(tmp_path)
    output = tmp_path / "native-real"
    command = ["cafe5", "--infile", str(copies), "--tree", str(tree),
               "--fixed_lambda", "0.1", "--fixed_alpha", "1", "--n_gamma_cats", "2",
               "--cores", "1", "--output_prefix", str(output)]
    p = subprocess.run(command, capture_output=True, text=True, timeout=30)
    (tmp_path / "native-fit.log").write_text(p.stdout + p.stderr)
    assert p.returncode == 0, p.stdout + p.stderr
    before = {p.name: p.read_bytes() for p in output.iterdir() if p.is_file()}
    result = validator.validate(copies, tree, output)
    assert result["validated"] and result["species"] == 3
    assert result["reconstructed_families"] == 1
    assert result["input_families_without_reconstruction"] == 1
    assert {p.name: p.read_bytes() for p in output.iterdir() if p.is_file()} == before
    print("native Gamma fit:", result,
          "count header:", (output / "Gamma_count.tab").read_text().splitlines()[0],
          "probability header:", (output / "Gamma_branch_probabilities.tab").read_text().splitlines()[0])

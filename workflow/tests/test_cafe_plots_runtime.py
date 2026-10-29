"""Render saved CAFE reconstructions without repeating family-size inference."""
import hashlib
import json
import os
import shlex
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / "workflow/support"


def saved_cafe(root):
    root.mkdir(parents=True, exist_ok=True)
    tree = "((Dioncophyllum_thollonii<1>_2:1,Dionaea_muscipula<2>_4*:1)<4>_6:1,Amborella_trichopoda<3>_1:2)<5>_7;"
    (root / "Gamma_asr.tre").write_text(
        "#nexus\nBEGIN TREES;\n" + "".join(f" TREE OG{i:06d} = {tree}\n" for i in (1, 2)) + "END;\n")
    fields = "FamilyID\tDioncophyllum_thollonii<1>\tDionaea_muscipula<2>\tAmborella_trichopoda<3>\t<4>\t<5>\n"
    (root / "Gamma_count.tab").write_text(fields + "OG000001\t2\t4\t1\t6\t7\nOG000002\t2\t4\t1\t6\t7\n")
    (root / "Gamma_change.tab").write_text(fields + "OG000001\t-4\t-2\t-6\t-1\t0\nOG000002\t1\t2\t0\t-1\t0\n")
    (root / "Gamma_branch_probabilities.tab").write_text(fields + "OG000001\t0.1\t0.01\t0.1\t0.1\tNA\n")
    (root / "Gamma_family_results.txt").write_text("#FamilyID\tpvalue\tSignificant at 0.05\nOG000001\t0.01\ty\n")
    return root


def hashes(root):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in root.iterdir() if p.is_file()}


@pytest.mark.parametrize("cores", [1, 2])
def test_saved_multispecies_cafe_plots_in_psock_workers(tmp_path, cores):
    native = saved_cafe(tmp_path / "native")
    before = hashes(native)
    out = tmp_path / "plots"
    result = subprocess.run(["Rscript", str(SUPPORT / "cafe_plot_each_family.r"),
        str(native / "Gamma_asr.tre"), str(native / "Gamma_count.tab"),
        str(native / "Gamma_change.tab"), str(out), str(cores)],
        text=True, capture_output=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    for name in ("OG000001", "OG000002"):
        assert (out / f"{name}.pdf").read_bytes().startswith(b"%PDF")
        svg = (out / f"{name}.svg").read_text()
        assert "Dioncophyllum thollonii" in svg and "Amborella trichopoda" in svg
    for script, argv in [("cafe_plot_summary.r", ["Gamma_asr.tre", "Gamma_change.tab"]),
                         ("cafe_plot_branch_id.r", ["Gamma_asr.tre"])]:
        result = subprocess.run(["Rscript", str(SUPPORT / script),
            *[str(native / name) for name in argv], str(out / "summary")],
            text=True, capture_output=True, timeout=120)
        assert result.returncode == 0, result.stdout + result.stderr
    for name in ("summary_all", "summary_significant", "branch_id"):
        assert (out / "summary" / f"{name}.pdf").read_bytes().startswith(b"%PDF")
        assert "<svg" in (out / "summary" / f"{name}.svg").read_text()
    assert hashes(native) == before


def test_cafe_labels_handle_vectors_internal_nodes_and_qualified_species():
    code = f'''
source({json.dumps(str(SUPPORT / "tree_label_formatting.R"))})
source({json.dumps(str(SUPPORT / "cafe_plot_labels.R"))})
x <- c("Dioncophyllum_thollonii<1>_2*", "Oryza_sativa_Indica_Group<2>_3", "<5>_7", NA_character_)
stopifnot(identical(norm_key(x), c("Dioncophyllum_thollonii<1>", "Oryza_sativa_Indica_Group<2>", "<5>", NA_character_)))
stopifnot(identical(get_species(x[1:2]), c("Dioncophyllum thollonii", "Oryza sativa Indica Group")))
stopifnot(identical(norm_key(character()), character()), identical(get_species(character()), character()))
'''
    p = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, timeout=30)
    assert p.returncode == 0, p.stdout + p.stderr


def run_cafe_stage(tmp_path):
    core = (ROOT / "workflow/core/gg_genome_evolution_core.sh").read_text()
    stage = "task='CAFE analysis'" + core.split("task='CAFE analysis'", 1)[1].split('task="Orthogroup copy-number trait PGLS"', 1)[0]
    cafe = tmp_path / "output/cafe"
    env = {**os.environ, "gg_workspace_dir": str(tmp_path), "gg_workspace_output_dir": str(tmp_path / "output"),
           "gg_support_dir": str(SUPPORT), "genome_evolution_provenance_dir": str(tmp_path / "provenance"),
           "file_orthogroup_copy_number": str(tmp_path / "copies.tsv"), "file_dated_species_tree": str(tmp_path / "tree.nwk"),
           "dir_cafe": str(cafe), "dir_cafe_output": str(cafe / "cafe_output"),
           "file_cafe_summary_all_pdf": str(cafe / "summary_plot/summary_all.pdf"),
           "file_cafe_summary_significant_pdf": str(cafe / "summary_plot/summary_significant.pdf"),
           "n_gamma_cats_cafe": "4", "GG_TASK_CPUS": "2", "run_cafe": "1"}
    prelude = f'''set -euo pipefail
source {shlex.quote(str(SUPPORT / "gg_util/06_workspace_validation.sh"))}
# Only transport/logging are isolated; use the real contract builder and validator.
gg_artifact_server_request() {{ local script=$1; shift; python "$script" "$@"; }}
gg_step_start() {{ :; }}
gg_step_skip() {{ :; }}
ensure_dir() {{ mkdir -p "$1"; }}
cafe5() {{ echo "unexpected optimizer invocation" >&2; exit 93; }}
'''
    return subprocess.run(["bash", "-c", prelude + stage], env=env, capture_output=True, text=True, timeout=120)


def test_native_stage_reuses_completed_fit_and_checks_contract_on_replay(tmp_path):
    native = saved_cafe(tmp_path / "output/cafe/cafe_output")
    (tmp_path / "copies.tsv").write_text("Desc\tFamilyID\tDioncophyllum_thollonii\tDionaea_muscipula\tAmborella_trichopoda\n(null)\tOG000001\t2\t4\t1\n(null)\tOG000002\t2\t4\t1\n")
    (tmp_path / "tree.nwk").write_text("((Dioncophyllum_thollonii:1,Dionaea_muscipula:1):1,Amborella_trichopoda:2);\n")
    before = hashes(native)
    first = run_cafe_stage(tmp_path)
    assert first.returncode == 0, first.stdout + first.stderr
    assert "Skipping CAFE run" in first.stdout
    assert hashes(native) == before
    inference = json.loads((tmp_path / "provenance/cafe.inference.json").read_text())
    assert inference["parameters"] == {"gamma_categories": "4", "pvalue": "0.05"}
    assert inference["diagnostics"]["provenance_state"] == "adopted_legacy_output_without_rebuild"
    saved = {p: p.read_bytes() for p in tmp_path.rglob("*") if p.is_file() and not p.name.endswith(".lock")}
    replay = run_cafe_stage(tmp_path)
    assert replay.returncode == 0, replay.stdout + replay.stderr
    assert "Skipping CAFE run" in replay.stdout and "Skipping rendering" in replay.stdout
    assert all(p.read_bytes() == data for p, data in saved.items())

    # A changed input remains a hard stop, with the saved inference untouched.
    (tmp_path / "copies.tsv").write_text("changed scientific input\n")
    rejected = run_cafe_stage(tmp_path)
    assert rejected.returncode != 0
    assert "Skipping CAFE run" not in rejected.stdout
    assert hashes(native) == before


def test_native_stage_rejects_foreign_legacy_fit_before_adoption(tmp_path):
    native = saved_cafe(tmp_path / "output/cafe/cafe_output")
    before = hashes(native)
    (tmp_path / "copies.tsv").write_text("Desc\tFamilyID\tDioncophyllum_thollonii\n(null)\tOG000001\t2\n")
    (tmp_path / "tree.nwk").write_text("(Dioncophyllum_thollonii:1,Dionaea_muscipula:1);\n")
    result = run_cafe_stage(tmp_path)
    assert result.returncode != 0
    assert "CAFE semantic validation failed" in result.stderr
    assert not (tmp_path / "provenance/cafe.inference.json").exists()
    assert hashes(native) == before


def test_native_stage_allows_explicit_legacy_rebuild(tmp_path, monkeypatch):
    saved_cafe(tmp_path / "output/cafe/cafe_output")
    (tmp_path / "copies.tsv").write_text("Desc\tFamilyID\tDioncophyllum_thollonii\tDionaea_muscipula\n(null)\tOG000001\t2\t4\n")
    (tmp_path / "tree.nwk").write_text("(Dioncophyllum_thollonii:1,Dionaea_muscipula:1);\n")
    monkeypatch.setenv("artifact_legacy_policy", "rebuild")
    result = run_cafe_stage(tmp_path)
    assert result.returncode == 93
    assert "unexpected optimizer invocation" in result.stderr
    assert "CAFE semantic validation failed" not in result.stderr

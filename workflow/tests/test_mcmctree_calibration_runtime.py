"""Compatibility and real-engine checks for opt-in calibration auditing."""

import csv
import importlib
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))
audit = importlib.import_module("mcmctree_calibration_audit")
experiments = importlib.import_module("mcmctree_calibration_experiments")


def fixture_manifest(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("((a:0.1,b:0.1)'B(10,20,0.025,0.025)':0.1,c:0.2)'B(30,40,0.025,0.025)';\n")
    snapshot = audit.inventory(tree, tmp_path / "audit")
    with (snapshot / "candidates.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        row.update(decision="accept", source_type="secondary", source_references="synthetic-study-A",
                   interval_kind="study_posterior", dependency_groups='["shared-fossil"]',
                   node_basis="Synthetic known clade", distribution_basis="Synthetic soft bounds",
                   reviewer="test", review_note="Synthetic fixture, not biological evidence")
    manifest = tmp_path / "review.tsv"
    manifest.write_text(audit.table_text(rows))
    return tree, manifest, rows


def test_legacy_inventory_is_additive_idempotent_and_never_infers_sources(tmp_path):
    tree = tmp_path / "tree.nwk"
    original = "((a,b)'B(10,20,0.025,0.025)',c);\n"
    tree.write_text(original)
    first = audit.inventory(tree, tmp_path / "audit")
    before = {path.name: path.read_bytes() for path in first.iterdir()}
    assert audit.inventory(tree, tmp_path / "audit") == first
    assert {path.name: path.read_bytes() for path in first.iterdir()} == before
    assert tree.read_text() == original
    status = json.loads((first / "status.json").read_text())
    assert status["status"] == "unreviewed"
    assert status["historical_mcmc_diagnostics"] == "unknown"


def test_reviewed_selection_uses_exact_clades_preserves_inputs_and_soft_tails(tmp_path):
    tree, manifest, rows = fixture_manifest(tmp_path)
    inputs = tree.read_bytes(), manifest.read_bytes()
    rows[0]["decision"] = "exclude"
    manifest.write_text(audit.table_text(rows))
    output, report = tmp_path / "selected.nwk", tmp_path / "review.json"
    audit.apply_manifest(tree, manifest, output, report)
    assert output.read_text().count("B(") == 1
    assert "0.025,0.025" in output.read_text()
    assert tree.read_bytes() == inputs[0]
    assert json.loads(report.read_text())["status"] == "reviewed_inputs"
    # Sibling order is not node identity; a different topology is.
    tree.write_text("(c:0.2,(b:0.1,a:0.1):0.1);\n")
    audit.apply_manifest(tree, manifest, output, report)
    old = output.read_bytes()
    tree.write_text("((a:0.1,c:0.1):0.1,b:0.2);\n")
    with pytest.raises(ValueError, match="changed topology"):
        audit.apply_manifest(tree, manifest, output, report)
    assert output.read_bytes() == old


def test_current_status_does_not_apply_review_to_changed_calibrations(tmp_path):
    tree, manifest, _ = fixture_manifest(tmp_path)
    output = tmp_path / "selected.nwk"
    audit.apply_manifest(tree, manifest, output, tmp_path / "reviewed_calibrations.json")
    audit.inventory(output, tmp_path / "audit", manifest)
    status = tmp_path / "audit/current_status.json"
    assert json.loads(status.read_text())["status"] == "reviewed_inputs"
    audit.inventory(output, tmp_path / "audit")
    assert json.loads(status.read_text())["status"] == "unreviewed"
    output.write_text(output.read_text().replace("10,20", "11,20"))
    audit.inventory(output, tmp_path / "audit", manifest)
    assert json.loads(status.read_text())["status"] == "unreviewed"


@pytest.mark.parametrize("field,value", [
    ("reviewer", ""), ("interval_kind", "unknown"), ("time_unit", "years"),
    ("dependency_groups", "[]"), ("calibration", "B(10,10,0.025,0.025)"),
    ("calibration", "B(10,20,0.8,0.8)"), ("calibration", "@10"),
    ("calibration", "B(nan,20,0.025,0.025)"),
])
def test_bad_review_does_not_replace_output(tmp_path, field, value):
    tree, manifest, rows = fixture_manifest(tmp_path)
    rows[0][field] = value
    manifest.write_text(audit.table_text(rows))
    output = tmp_path / "out.nwk"
    output.write_text("previous result")
    with pytest.raises(ValueError):
        audit.apply_manifest(tree, manifest, output, tmp_path / "status.json")
    assert output.read_text() == "previous result"


def test_shared_group_omission_does_not_run_without_any_anchor(tmp_path):
    tree, manifest, _ = fixture_manifest(tmp_path)
    scenarios = experiments.scenarios(audit.load_tree(tree), manifest)
    assert len(scenarios) == 4
    assert scenarios[-1]["status"] == "skipped_no_remaining_calibration"
    assert len(scenarios[-1]["omit_nodes"]) == 2
    assert all(row["status"] == "pending" for row in scenarios[:-1])


def test_archive_preserves_old_public_and_raw_outputs(tmp_path):
    source = tmp_path / "mcmctree_main"
    source.mkdir()
    (source / "FigTree.tre").write_text("old dated tree")
    audit.archive_existing(tmp_path)
    archived = list((tmp_path / "calibration_history").glob("*/mcmctree_main/FigTree.tre"))
    assert len(archived) == 1
    assert archived[0].read_bytes() == (source / "FigTree.tre").read_bytes()


def test_diagnostics_distinguishes_stationary_shifted_and_constant_chains(tmp_path):
    rng = np.random.default_rng(419)
    paths = []
    for i in range(4):
        path = tmp_path / f"chain{i}.tsv"
        pd.DataFrame({"Gen": np.arange(2000), "stationary": rng.normal(size=2000),
                      "shifted": rng.normal(i * 3, size=2000), "stuck": np.zeros(2000)}).to_csv(path, sep="\t", index=False)
        paths.append(path)
    output = tmp_path / "diagnostics.tsv"
    subprocess.run(["Rscript", str(SUPPORT / "mcmctree_calibration_diagnostics.R"),
                    str(output), *map(str, paths)], check=True, capture_output=True)
    result = pd.read_csv(output, sep="\t").set_index("parameter")
    assert result.loc["stationary", "rhat"] < 1.01
    assert result.loc["stationary", "ess_bulk"] > 400
    assert result.loc["stationary", "ess_tail"] > 400
    assert result.loc["shifted", "rhat"] > 1.1
    assert np.isnan(result.loc["stuck", "rhat"])


@pytest.mark.parametrize("time_scale", ["1", "100"])
def test_real_iq2mc_prior_and_posterior_keep_inputs_and_detect_tampered_cache(tmp_path, time_scale):
    # Real IQ-TREE Hessian, not a fabricated derivative matrix. This is a small
    # execution test; its short chains intentionally cannot establish convergence.
    rng = np.random.default_rng(15)
    ancestral = rng.choice(list("ACGT"), size=500)
    alignment = tmp_path / "alignment.fa"
    sequences = []
    for name in "abc":
        sequence = ancestral.copy()
        mask = rng.random(500) < 0.15
        sequence[mask] = rng.choice(list("ACGT"), size=mask.sum())
        sequences.append(f">{name}\n{''.join(sequence)}\n")
    alignment.write_text("".join(sequences))
    tree = tmp_path / "input.nwk"
    tree.write_text("((a:0.1,b:0.1)'B(0.1,0.2,0.025,0.025)':0.1,c:0.2)'B(0.3,0.4,0.025,0.025)';\n")
    subprocess.run(["iqtree3", "-s", str(alignment), "-te", str(tree), "-m", "JC", "-T", "1",
                    "--dating", "mcmctree", "--mcmc-bds", "1,1,0.5", "--mcmc-iter", "20,2,40",
                    "--prefix", str(tmp_path / "iq2mc"), "--seed", "37"],
                   check=True, capture_output=True, timeout=60)
    control = tmp_path / "iq2mc.mcmctree.ctl"
    # Use the existing GeneGalleon normalizer, keeping any PAML adaptation in
    # the same code path as a normal workflow execution.
    core = (SUPPORT.parent / "core/gg_genome_evolution_core.sh").read_text()
    functions = core[core.index("mcmctree_requires_bdparas_flag() {"):core.index("iq2mc_option_supported() {")]
    control.write_text(control.read_text().replace("FossilErrprint", "FossilErr\nprint"))
    subprocess.run(["bash", "-c", functions + '\n dir_tmp="$1"; normalize_mcmctree_ctl_for_installed_paml "$2"',
                    "test", str(tmp_path), str(control)], check=True, capture_output=True)
    args = SimpleNamespace(control=control, tree=tmp_path / "iq2mc.rooted.nwk",
                           alignment=tmp_path / "iq2mc.dummy.phy", hessian=tmp_path / "iq2mc.mcmctree.hessian",
                           outdir=tmp_path / "experiments", manifest=None, chains=2, seed=1729, time_scale=time_scale)
    from decimal import Decimal

    from mcmctree_time_scale import scale_calibration_labels, scale_ctl_rootage_text

    args.tree.write_text(scale_calibration_labels(args.tree.read_text(), Decimal(time_scale), "up"))
    control.write_text(scale_ctl_rootage_text(control.read_text(), Decimal(time_scale), "up"))
    originals = {path: path.read_bytes() for path in (args.control, args.tree, args.alignment, args.hessian)}
    directory = experiments.run_experiments(args)
    status = json.loads((directory / "status.json").read_text())
    assert status["status"] == "not_converged"
    assert len(list(directory.glob("*/prior/chain_*/mcmc.txt"))) == 6
    assert len(list(directory.glob("*/posterior/chain_*/mcmc.txt"))) == 6
    assert status["group_sensitivity"] == "unavailable_without_manifest"
    assert (directory / "comparisons.tsv").is_file()
    summary = pd.read_csv(directory / "summary.tsv", sep="\t")
    public_root = summary.loc[(summary["scenario"] == "baseline") & (summary["target"] == "prior") & (summary["parameter"] == "t_n4")].iloc[0]
    assert 0.25 * float(time_scale) < public_root["median"] < 0.5 * float(time_scale)
    assert public_root["unit"] == "Ma"
    assert "checkpoint = 0" in next(directory.glob("*/prior/chain_*/mcmctree.ctl")).read_text()
    assert all(path.read_bytes() == before for path, before in originals.items())
    assert experiments.run_experiments(args) == directory
    sample = next(directory.glob("*/prior/chain_*/mcmc.txt"))
    sample.write_text("corrupt")
    with pytest.raises(ValueError, match="changed or missing"):
        experiments.run_experiments(args)


def test_cached_legacy_core_stage_adds_inventory_without_running_nwkit_or_replacing_tree(tmp_path):
    import shlex

    core = (SUPPORT.parent / "core/gg_genome_evolution_core.sh").read_text()
    section = core[core.index('task="Time-constrained tree preparation"'):core.index('task="Constrained range plotting"')]
    tree, _, _ = fixture_manifest(tmp_path)
    before = tree.read_bytes()
    variables = {"run_mcmctree_calibration_diagnostics": "0", "mcmctree_calibration_manifest": "",
                 "run_constrained_tree": "1", "artifact_stale_policy": "stop",
                 "file_undated_species_tree": str(tree), "file_constrained_tree": str(tree),
                 "genome_evolution_provenance_dir": str(tmp_path / "provenance"),
                 "timetree_constraint": "1", "mcmctree_divergence_time_constraints_str": "",
                 "species_label_parser": "taxonomic", "species_label_regex": "", "species_label_map_tsv": "",
                 "gg_support_dir": str(SUPPORT), "dir_tmp": str(tmp_path), "captured": str(tmp_path / "contract")}
    setup = "set -euo pipefail\n" + "\n".join(f"{key}={shlex.quote(value)}" for key, value in variables.items())
    setup += '''
disable_if_no_input_file() { :; }
gg_artifact_contract_init() { constrained_tree_provenance_args=(); }
gg_artifact_add_input_if_present() { :; }
gg_artifact_prepare_stage() { printf -v "$1" 0; printf '%s\\n' "${@:3}" > "$captured"; }
gg_step_skip() { :; }
nwkit() { echo "Unexpected NWKIT execution" >&2; exit 99; }
'''
    subprocess.run(["bash", "-c", setup + "\n" + section], check=True, capture_output=True)
    assert tree.read_bytes() == before
    assert (tree.parent / "calibration_audit/current_status.json").is_file()
    assert "reviewed_manifest=" not in (tmp_path / "contract").read_text()


def test_legacy_configuration_contract_does_not_gain_required_audit_input():
    core = (SUPPORT.parent / "core/gg_genome_evolution_core.sh").read_text()
    section = core[core.index('task="Time-constrained tree preparation"'):core.index('task="Constrained range plotting"')]
    guard = section.index('if [[ -n "${mcmctree_calibration_manifest}" ]]')
    assert section.index('--input "reviewed_manifest=') > guard
    assert '--parameter "timetree_statistic=ci"' in section[:guard]
    assert '--parameter "min_clade_proportion=0.2"' in section[:guard]
    assert section.index('mcmctree_calibration_audit.py" archive') < section.index('elif [[ ${timetree_constraint} -eq 1 ]]')

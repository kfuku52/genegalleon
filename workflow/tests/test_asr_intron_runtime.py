"""Fixed-Q intron probabilities, publication, and real core-stage invalidation."""

import itertools
import json
import math
import os
import shlex
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
from ete4 import Tree

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / "workflow" / "support"
CORE = ROOT / "workflow" / "core" / "gg_gene_evolution_core.sh"


def exact_probabilities(observed, length):
    # Independent enumeration, analytic two-state CTMC; no NWKIT likelihood code.
    edges = [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)]
    stationary = (10 / 11, 1 / 11)
    decay = math.exp(-0.0011 * length)
    total, sums = 0.0, [0.0] * 7
    for assignment in itertools.product((0, 1), repeat=7):
        if any(value is not None and assignment[i + 3] != int(value > 0) for i, value in enumerate(observed)):
            continue
        weight = 0.5
        for parent, child in edges:
            i, j = assignment[parent], assignment[child]
            weight *= stationary[j] + ((i == j) - stationary[j]) * decay
        total += weight
        for i, value in enumerate(assignment):
            sums[i] += weight * value
    return [value / total for value in sums]


def inputs(tmp_path, values, length=100):
    tree = tmp_path / "dated.nwk"
    gff = tmp_path / "gff.tsv"
    tree.write_text(f"((A:{length},B:{length})N1:{length},(C:{length},D:{length})N2:{length})Root;\n")
    # Reversed rows and a gene absent from the pruned tree exercise the join.
    gff.write_text("gene_id\tnum_intron\n" + "".join(
        f"{gene}\t{value if value is not None else 'NA'}\n"
        for gene, value in reversed(list(zip("ABCD", values, strict=True)))
    ) + "pruned_gene\t7\n")
    return tree, gff


def command(tree, gff, prefix):
    return [sys.executable, str(SUPPORT / "asr_intron_evolution.py"),
            "--tree-file", str(tree), "--trait-file", str(gff),
            "--intron-gain-rate", "0.0001", "--retrotransposition-rate", "0.001",
            "--output-prefix", str(prefix)]


@pytest.mark.parametrize("values,length", [
    ([2, 1, 0, 0], 100), ([1, 2, 3, 4], 100), ([0, 0, 0, 0], 1000),
    ([1, None, 2, None], 1000), ([1, None, 0, None], 100),
    ([None, None, None, None], 100),
])
def test_intron_asr_matches_exact_enumeration_and_preserves_counts(tmp_path, values, length):
    tree, gff = inputs(tmp_path, values, length)
    prefix = tmp_path / "intron"
    completed = subprocess.run(command(tree, gff, prefix), capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    summary = pd.read_csv(prefix.with_suffix(".tsv"), sep="\t").sort_values("branch_id")
    assert list(summary.name) == ["Root", "N1", "N2", "A", "B", "C", "D"]
    assert summary.p_intron_present.tolist() == pytest.approx(exact_probabilities(values, length), abs=1e-12)
    assert summary.loc[summary.node_class.ne("leaf"), "num_intron"].isna().all()
    tips = summary[summary.node_class.eq("leaf")].set_index("name")
    for gene, value in zip("ABCD", values, strict=True):
        assert bool(tips.loc[gene, "is_imputed"]) == (value is None)
        if value is None:
            assert pd.isna(tips.loc[gene, "num_intron"])
        else:
            assert tips.loc[gene, "num_intron"] == value
    model = pd.read_csv(str(prefix) + ".model.tsv", sep="\t").iloc[0]
    assert model["model"] == "CUSTOM"
    assert not model.rate_estimated
    assert model.root_prior == "equal"
    annotated = prefix.with_suffix(".nhx").read_text()
    assert "asr_p_intron_present=" in annotated
    assert "Present=" in annotated and "Absent=" in annotated
    assert prefix.with_suffix(".pdf").read_bytes().startswith(b"%PDF-")


def test_intron_asr_handles_absent_gff_rows_and_internal_gene_name_collision(tmp_path):
    tree, gff = inputs(tmp_path, [1, 0, 2, None])
    tree.write_text(tree.read_text().replace("N1", "A"))
    gff.write_text(gff.read_text().replace("D\tNA\n", ""))
    prefix = tmp_path / "intron"
    completed = subprocess.run(command(tree, gff, prefix), capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    summary = pd.read_csv(prefix.with_suffix(".tsv"), sep="\t")
    assert summary.loc[summary.node_class.ne("leaf"), "num_intron"].isna().all()
    assert summary.loc[summary.name.eq("D"), "is_imputed"].item()


def test_failed_intron_asr_does_not_publish_outputs(tmp_path):
    tree, gff = inputs(tmp_path, [1, 0, 1, 0])
    tree.write_text(tree.read_text().replace("A:100", "A:-100"))
    prefix = tmp_path / "intron"
    completed = subprocess.run(command(tree, gff, prefix), capture_output=True, text=True)
    assert completed.returncode != 0
    assert not list(tmp_path.glob("intron.*"))


def test_intron_asr_preserves_zero_short_and_non_ultrametric_lengths(tmp_path):
    tree, gff = inputs(tmp_path, [1, 1, 0, 0])
    tree.write_text("((A:0,B:0.00001)N1:0.1,(C:100,D:99)N2:0.01)Root;")
    prefix = tmp_path / "intron"
    result = subprocess.run(command(tree, gff, prefix), capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    before = Tree(tree.read_text(), parser=1)
    after = Tree(prefix.with_suffix(".nhx").read_text(), parser=1)
    original = {node.name: node.dist for node in before.traverse() if not node.is_root}
    result = {node.name: node.dist for node in after.traverse() if not node.is_root}
    assert result == original


def test_intron_asr_flows_into_real_branch_statistics_and_parent_deltas(tmp_path):
    tree, gff = inputs(tmp_path, [2, None, 0, None])
    names = dict(zip("ABCD", ["Alpha_one_g1", "Beta_two_g1", "Gamma_three_g1", "Delta_four_g1"], strict=True))
    text = tree.read_text()
    traits = gff.read_text()
    for letter, name in names.items():
        text = text.replace(f"{letter}:", f"{name}:")
        traits = traits.replace(f"{letter}\t", f"{name}\t")
    tree.write_text(text)
    gff.write_text(traits)
    species = tmp_path / "species.nwk"
    species.write_text(text.replace("_g1", ""))
    # A synthetic GeneRax reconciliation supplies one root transfer candidate.
    reconciled = Tree(text, parser=1)
    for node in reconciled.traverse():
        node.add_props(S="0", H="Y" if node.is_root else "N", D="N")
    generax = tmp_path / "generax.nhx"
    generax.write_text(reconciled.write(props=["S", "H", "D"], parser=1, format_root_node=True))
    prefix = tmp_path / "intron"
    result = subprocess.run(command(tree, gff, prefix), capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    result = subprocess.run([
        sys.executable, str(SUPPORT / "orthogroup_statistics.py"),
        "--rooted_tree", str(tree), "--dated_tree", str(tree), "--species_tree", str(species),
        "--character_gff", str(gff), "--asr_intron", str(prefix.with_suffix(".tsv")),
        "--generax_nhx", str(generax), "--ncpu", "1",
    ], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    branch = pd.read_csv(tmp_path / "orthogroup.branch.tsv", sep="\t").set_index("branch_id")
    native = pd.read_csv(prefix.with_suffix(".tsv"), sep="\t").set_index("name")
    assert len(branch) == len(native) == 7
    for _, row in branch.iterrows():
        assert row.intron_present == pytest.approx(native.loc[row.node_name, "p_intron_present"])
        if row.parent in branch.index:
            expected = row.intron_present - branch.loc[row.parent, "intron_present"]
            assert row.delta_intron_present == pytest.approx(expected)
    by_name = branch.set_index("node_name")
    assert by_name.loc[names["A"], "num_intron"] == 2
    assert by_name.loc[names["C"], "num_intron"] == 0
    assert by_name.loc[["Root", "N1", "N2", names["B"], names["D"]], "num_intron"].isna().all()
    assert by_name.loc[[names["B"], names["D"]], "intron_is_imputed"].all()
    branch_dir = tmp_path / "stat_branch"
    tree_dir = tmp_path / "stat_tree"
    branch_dir.mkdir()
    tree_dir.mkdir()
    (tmp_path / "orthogroup.branch.tsv").rename(branch_dir / "OG1_branch.tsv")
    (tmp_path / "orthogroup.tree.tsv").rename(tree_dir / "OG1_tree.tsv")
    database = tmp_path / "families.sqlite"
    result = subprocess.run([
        sys.executable, str(SUPPORT / "generate_orthogroup_database.py"),
        "--dir_stat_branch", str(branch_dir), "--dir_stat_tree", str(tree_dir),
        "--dbpath", str(database), "--ncpu", "1",
    ], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    result = subprocess.run([
        sys.executable, str(SUPPORT / "score_hgt_candidates.py"), "--dbpath", str(database),
        "--branch_out", str(tmp_path / "hgt.branch.tsv"),
        "--gene_out", str(tmp_path / "hgt.gene.tsv"),
        "--orthogroup_out", str(tmp_path / "hgt.family.tsv"),
    ], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    candidates = pd.read_csv(tmp_path / "hgt.branch.tsv", sep="\t")
    assert len(candidates) == 1
    assert candidates.intron_measured_gene_count.item() == 2
    assert candidates.intron_supported_gene_count.item() == 1
    assert candidates.intron_support_fraction.item() == 0.5


def test_intron_core_stage_caches_and_rebuilds_model_tree_and_plot(tmp_path):
    tree, gff = inputs(tmp_path, [1, None, 0, None])
    core = CORE.read_text()
    begin = core.index('task="Ancestral reconstruction of intron presence with NWKIT"')
    stage = core[begin:core.index('\ntask="NWKIT OU shift detection"', begin)]
    output = tmp_path / "output" / "orthogroup"
    declarations = "\n".join(line for line in core.splitlines() if line.startswith("file_og_asr_intron_"))
    script = tmp_path / "stage.sh"
    script.write_text(f'''set -euo pipefail
gg_support_dir={shlex.quote(str(SUPPORT))}
source "${{gg_support_dir}}/gg_util.sh"
gg_workspace_dir={shlex.quote(str(tmp_path))}
dir_output_active={shlex.quote(str(output))}
GG_FAMILY_OUTPUT_ROOT="${{dir_output_active}}"
GG_FAMILY_OUTPUT_INVENTORY="${{gg_workspace_dir}}/inventory"
mkdir -p "${{GG_FAMILY_OUTPUT_INVENTORY}}"
og_id=OG1
file_og_gff_info={shlex.quote(str(gff))}
file_og_dated_tree_analysis={shlex.quote(str(tree))}
run_asr_intron=1
intron_gain_rate=${{1:-0.0001}}
retrotransposition_rate=0.001
artifact_stale_policy=${{GG_TEST_STALE_POLICY:-rebuild}}
mv() {{
  command mv "$@" || return $?
  if [[ ${{GG_TEST_INTERRUPT:-0}} == 1 && ${{2:-}} == *".gg-stage."* && ${{3:-}} == "${{file_og_asr_intron_tree}}" ]]; then
    kill -KILL -- "-$$"
  fi
}}
gg_step_start() {{ echo run >> stage-runs.txt; }}
gg_step_skip() {{ echo skip >> stage-runs.txt; }}
{declarations}
{stage}
''')

    def run(rate="0.0001"):
        result = subprocess.run(["bash", str(script), rate], cwd=tmp_path, capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr

    run()
    run()
    assert (tmp_path / "stage-runs.txt").read_text().splitlines() == ["run", "skip"]
    plot = output / "asr_intron_plot" / "OG1_asr.intron.pdf"
    plot.unlink()
    run()
    assert plot.is_file()
    before = pd.read_csv(output / "asr_intron_summary" / "OG1_asr.intron.tsv", sep="\t")
    run("0.01")
    after = pd.read_csv(output / "asr_intron_summary" / "OG1_asr.intron.tsv", sep="\t")
    assert not before.p_intron_present.equals(after.p_intron_present)
    assert (tmp_path / "stage-runs.txt").read_text().splitlines() == ["run", "skip", "run", "run"]
    manifest = json.loads((output / "artifact_provenance" / "OG1.asr_intron.json").read_text())
    assert manifest["step"] == "asr_intron"
    saved = {path: path.read_bytes() for path in output.glob("asr_intron*/OG1*")}
    failed = subprocess.run(["bash", str(script), "-1"], cwd=tmp_path, capture_output=True, text=True)
    assert failed.returncode != 0
    assert {path: path.read_bytes() for path in saved} == saved
    journal = b"".join(path.read_bytes() for path in (tmp_path / "inventory").glob("*.paths"))
    recorded = set(journal.decode().split("\0"))
    assert {str(path.relative_to(output)) for path in saved} <= recorded
    # Kill the real bundle publisher after installing two new artifacts.
    interrupted = subprocess.run(["bash", str(script), "0.02"], cwd=tmp_path,
                                 env={**os.environ, "GG_TEST_INTERRUPT": "1"}, start_new_session=True,
                                 capture_output=True, text=True, timeout=90)
    assert interrupted.returncode != 0
    assert json.loads((output / "artifact_provenance" / "OG1.asr_intron.json").read_text()) == manifest
    stopped = subprocess.run(["bash", str(script), "0.01"], cwd=tmp_path,
                             env={**os.environ, "GG_TEST_STALE_POLICY": "stop"},
                             capture_output=True, text=True, timeout=90)
    assert stopped.returncode != 0, "Interrupted output must not be reused as a complete result"
    run("0.01")
    recovered = pd.read_csv(output / "asr_intron_summary" / "OG1_asr.intron.tsv", sep="\t")
    pd.testing.assert_frame_equal(recovered, after)


@pytest.mark.parametrize("alias", ["direct", "symlink", "hardlink"])
def test_intron_asr_rejects_output_alias_of_input(tmp_path, alias):
    tree, gff = inputs(tmp_path, [1, 0, 1, 0])
    prefix = gff.with_suffix("") if alias == "direct" else tmp_path / "intron"
    target = prefix.with_suffix(".tsv")
    if alias == "symlink":
        target.symlink_to(gff)
    elif alias == "hardlink":
        target.hardlink_to(gff)
    before = gff.read_bytes()
    result = subprocess.run(command(tree, gff, prefix), capture_output=True, text=True)
    assert result.returncode != 0
    assert "must not overwrite input" in result.stderr
    assert gff.read_bytes() == before


def test_intron_asr_rejects_invalid_output_before_replacing_previous_results(tmp_path):
    tree, gff = inputs(tmp_path, [1, 0, 1, 0])
    prefix = tmp_path / "intron"
    prefix.with_suffix(".nhx").mkdir()
    model = tmp_path / "intron.model.tsv"
    model.write_text("previous model\n")
    result = subprocess.run(command(tree, gff, prefix), capture_output=True, text=True)
    assert result.returncode != 0
    assert model.read_text() == "previous model\n"
    assert not prefix.with_suffix(".pdf").exists()


def test_intron_asr_rolls_back_failed_output_installation(tmp_path, monkeypatch):
    import importlib.util
    from types import SimpleNamespace

    import nwkit.output_transaction as transaction

    spec = importlib.util.spec_from_file_location("intron_adapter", SUPPORT / "asr_intron_evolution.py")
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    tree, gff = inputs(tmp_path, [1, 0, 1, 0])
    prefix = tmp_path / "intron"
    paths = [Path(str(prefix) + suffix) for suffix in (".model.tsv", ".nhx", ".pdf", ".tsv")]
    for path in paths:
        path.write_text("previous " + path.name)
    before = {path: path.read_bytes() for path in paths}
    original = transaction.replace_output

    def fail_second(source, target):
        if str(target) == str(paths[1]) and ".stage." in str(source):
            raise OSError("simulated installation failure")
        original(source, target)

    monkeypatch.setattr(transaction, "replace_output", fail_second)
    with pytest.raises(OSError, match="simulated installation failure"):
        adapter.run_asr(SimpleNamespace(tree_file=str(tree), trait_file=str(gff),
                                       output_prefix=str(prefix), intron_gain_rate=0.0001,
                                       retrotransposition_rate=0.001))
    assert {path: path.read_bytes() for path in paths} == before


@pytest.mark.parametrize("gain,loss,values,expected", [
    (0, 0, [None] * 4, 0.5),
    (0, 0, [0] * 4, 0.0),
    (0, 0, [1] * 4, 1.0),
    (0, 0, [1, 0, 1, 0], None),
])
def test_intron_asr_zero_rate_boundary(tmp_path, gain, loss, values, expected):
    tree, gff = inputs(tmp_path, values)
    prefix = tmp_path / "intron"
    args = command(tree, gff, prefix)
    args[args.index("--intron-gain-rate") + 1] = str(gain)
    args[args.index("--retrotransposition-rate") + 1] = str(loss)
    result = subprocess.run(args, capture_output=True, text=True)
    if expected is None:
        assert result.returncode != 0
        assert not prefix.with_suffix(".tsv").exists()
    else:
        assert result.returncode == 0, result.stdout + result.stderr
        summary = pd.read_csv(prefix.with_suffix(".tsv"), sep="\t")
        assert summary.p_intron_present.tolist() == pytest.approx([expected] * len(summary))

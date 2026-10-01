"""Species-root propagation and plain GeneRax inputs in the real runtime."""

import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / "workflow/support"
CORE = ROOT / "workflow/core/gg_gene_evolution_core.sh"
ANNOTATED = "(A_one:1.2345678901234567,(B_two:0.5,C_three:0.5)s2:0.5)s1[&&NHX:nwkit_rooted=yes];"


def run(*args, **kwargs):
    return subprocess.run([str(arg) for arg in args], capture_output=True, text=True, **kwargs)


@pytest.mark.parametrize("source_text", [ANNOTATED, "[&R]" + ANNOTATED.split("[&&NHX:")[0] + ";"])
def test_plain_tool_input_preserves_root_and_lengths_and_is_read_by_generax(tmp_path, source_text):
    from nwkit.util import read_tree

    source, plain = tmp_path / "source.nwk", tmp_path / "species.nwk"
    source.write_text(source_text)
    converted = run(sys.executable, SUPPORT / "prepare_generax_species_tree.py",
                    "--input", source, "--output", plain)
    assert converted.returncode == 0, converted.stdout + converted.stderr
    assert source.read_text() == source_text
    assert "[" not in plain.read_text()
    before, after = (read_tree(str(path), "auto", True) for path in (source, plain))
    assert [set(child.leaf_names()) for child in before.children] == [set(child.leaf_names()) for child in after.children]
    for first in before.leaf_names():
        for second in before.leaf_names():
            assert before.get_distance(first, second) == after.get_distance(first, second)
    assert "1.2345678901234567" in plain.read_text()

    (tmp_path / "genes.nwk").write_text("(A_one_g:0.1,(B_two_g:0.1,C_three_g:0.1):0.1);")
    (tmp_path / "alignment.fa").write_text("".join(
        f">{species}_g\nMKTLLILAVAAAQGGKKSTVAAAGGLLSLAAAMKTLAAAVAAAVGGKKSTV\n"
        for species in ("A_one", "B_two", "C_three")))
    (tmp_path / "map.txt").write_text("".join(f"{species}_g\t{species}\n" for species in ("A_one", "B_two", "C_three")))
    (tmp_path / "families.txt").write_text(
        "[FAMILIES]\n- family_1\nstarting_gene_tree = genes.nwk\n"
        "alignment = alignment.fa\nmapping = map.txt\nsubst_model = LG\n")
    env = dict(os.environ, OMPI_MCA_ras="^gridengine", OMPI_MCA_plm="isolated",
               OMPI_MCA_plm_rsh_agent="/bin/false", OMPI_MCA_btl="^openib")
    mpi = ["mpiexec", "--oversubscribe", "-np", "1"]
    if os.getuid() == 0:
        mpi.append("--allow-run-as-root")
    result = run(*mpi, "generax", "--species-tree", plain, "--families", "families.txt",
                 "--strategy", "SKIP", "--rec-model", "UndatedDTL", "--prefix", "result",
                 "--seed", "12345", cwd=tmp_path, env=env, timeout=45)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "result/species_trees/starting_species_tree.newick").is_file()


@pytest.mark.parametrize("text", ["", "(A:1,A:1);", "(A:1,B:1,C:1);", "[&U](A:1,(B:1,C:1):1);"])
def test_invalid_species_input_does_not_replace_previous_plain_tree(tmp_path, text):
    source, output = tmp_path / "source.nwk", tmp_path / "plain.nwk"
    source.write_text(text)
    output.write_text("previous output")
    result = run(sys.executable, SUPPORT / "prepare_generax_species_tree.py", "--input", source, "--output", output)
    assert result.returncode != 0
    assert source.read_text() == text
    assert output.read_text() == "previous output"


def test_species_serializer_refuses_to_overwrite_annotated_source(tmp_path):
    source = tmp_path / "source.nwk"
    source.write_text(ANNOTATED)
    result = run(sys.executable, SUPPORT / "prepare_generax_species_tree.py", "--input", source, "--output", source)
    assert result.returncode != 0
    assert source.read_text() == ANNOTATED


def pruning_harness(tmp_path):
    text = CORE.read_text()
    begin = text.index("prepare_species_tree_pruned() (")
    function = text[begin:text.index("\ncleanup_tmp_dir_on_normal_exit()", begin)]
    (tmp_path / "sequences").mkdir()
    (tmp_path / "parameters").mkdir()
    for species in ("A_one", "B_two", "C_three"):
        (tmp_path / "sequences" / f"{species}_cds.fasta").write_text(f">{species}_g\nATG\n")
    (tmp_path / "source").mkdir()
    source = tmp_path / "source/tree.nwk"
    source.write_text(ANNOTATED)
    script = tmp_path / "prune.sh"
    script.write_text(f'''set -euo pipefail
gg_support_dir={shlex.quote(str(SUPPORT))}
source "${{gg_support_dir}}/gg_util.sh"
gg_workspace_dir={shlex.quote(str(tmp_path))}
dir_output_active="${{gg_workspace_dir}}"
file_og_parameters_dir="${{gg_workspace_dir}}/parameters"
species_tree="${{gg_workspace_dir}}/source/tree.nwk"
species_tree_pruned="${{file_og_parameters_dir}}/species.pruned.nwk"
dir_sp_cds="${{gg_workspace_dir}}/sequences"
input_sequence_mode=cds
artifact_stale_policy=${{1:-rebuild}}
{function}
prepare_species_tree_pruned
''')
    return script, source, tmp_path / "parameters/species.pruned.nwk"


def test_pruned_cache_reuses_current_tree_and_refreshes_root_and_tip_selection(tmp_path):
    from nwkit.util import read_tree

    script, source, pruned = pruning_harness(tmp_path)
    pruned.write_text("(Wrong:1,Old:1);")  # Legacy first-writer-only cache.
    for _ in range(2):
        result = run("bash", script)
        assert result.returncode == 0, result.stdout + result.stderr
        if _ == 0:
            stamp = pruned.stat().st_mtime_ns
    assert pruned.stat().st_mtime_ns == stamp
    manifest = Path(str(pruned) + ".provenance.json")
    assert json.loads(manifest.read_text())["step"] == "species_tree_pruned"
    source.write_text("(C_three:0.5,(A_one:1.2345678901234567,B_two:0.5)s3:0.5)s1[&&NHX:nwkit_rooted=yes];")
    stopped = run("bash", script, "stop")
    assert stopped.returncode != 0  # A policy stop never accepts an old root.
    result = run("bash", script)
    assert result.returncode == 0, result.stdout + result.stderr
    assert set(read_tree(str(pruned), "auto", True).children[0].leaf_names()) == {"C_three"}
    (tmp_path / "sequences/B_two_cds.fasta").unlink()
    result = run("bash", script)
    assert result.returncode == 0, result.stdout + result.stderr
    assert set(read_tree(str(pruned), "auto", True).leaf_names()) == {"A_one", "C_three"}
    assert not Path(str(pruned) + ".lock").exists()


def test_pruned_cache_failure_preserves_previous_output_and_manifest(tmp_path):
    script, source, pruned = pruning_harness(tmp_path)
    result = run("bash", script)
    assert result.returncode == 0, result.stdout + result.stderr
    manifest = Path(str(pruned) + ".provenance.json")
    saved = pruned.read_bytes(), manifest.read_bytes()
    source.write_text("not a tree")
    result = run("bash", script)
    assert result.returncode != 0
    assert (pruned.read_bytes(), manifest.read_bytes()) == saved
    source.unlink()
    result = run("bash", script)
    assert result.returncode != 0
    assert (pruned.read_bytes(), manifest.read_bytes()) == saved
    assert not Path(str(pruned) + ".lock").exists()


def test_concurrent_pruning_uses_shared_source_contract(tmp_path):
    script, _, pruned = pruning_harness(tmp_path)
    processes = [subprocess.Popen(["bash", str(script)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                 for _ in range(2)]
    for process in processes:
        stdout, stderr = process.communicate(timeout=30)
        assert process.returncode == 0, stdout + stderr
    assert pruned.is_file()
    assert json.loads(Path(str(pruned) + ".provenance.json").read_text())["step"] == "species_tree_pruned"
    assert not Path(str(pruned) + ".lock").exists()


def test_generax_shared_species_cache_refreshes_when_source_root_changes(tmp_path):
    text = CORE.read_text()
    begin = text.index('  generax_out_sptree="')
    stage = text[begin:text.index('  echo "copying GeneRax output gene tree."', begin)]
    parameters = tmp_path / "parameters"
    parameters.mkdir()
    source = parameters / "source.nwk"
    source.write_text(ANNOTATED)
    generated = tmp_path / "generax_OG1/species_trees/starting_species_tree.newick"
    generated.parent.mkdir(parents=True)
    generated.write_text("(A_one:1,(B_two:1,C_three:1)s2:1)s1;")
    cached = parameters / "species.generax.nwk"
    cached.write_text("(Wrong:1,Old:1);")
    script = f'''set -euo pipefail
gg_support_dir={shlex.quote(str(SUPPORT))}
source "${{gg_support_dir}}/gg_util.sh"
gg_workspace_dir={shlex.quote(str(tmp_path))}
dir_output_active="${{gg_workspace_dir}}"
file_og_parameters_dir="${{gg_workspace_dir}}/parameters"
species_tree_pruned="${{file_og_parameters_dir}}/source.nwk"
species_tree_generax="${{file_og_parameters_dir}}/species.generax.nwk"
og_id=OG1
artifact_stale_policy=rebuild
{stage}
'''
    for _ in range(2):
        result = run("bash", "-s", input=script, cwd=tmp_path)
        assert result.returncode == 0, result.stdout + result.stderr
    assert cached.read_bytes() == generated.read_bytes()
    source.write_text("(C_three:1,(A_one:1,B_two:1)s2:1)s1;")
    generated.write_text(source.read_text())
    result = run("bash", "-s", input=script, cwd=tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    assert cached.read_bytes() == generated.read_bytes()
    generated.unlink()
    result = run("bash", "-s", input=script, cwd=tmp_path)
    assert result.returncode != 0  # Never use an old shared cache to hide a missing tool output.

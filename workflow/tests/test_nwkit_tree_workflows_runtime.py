"""Real NWKIT rendering and NOTUNG-preserving root selection in the GG runtime."""

import json
import shlex
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / "workflow" / "support"


def run(*args):
    return subprocess.run([str(value) for value in args], capture_output=True, text=True)


@pytest.mark.parametrize("text", [
    "(A:1,B:1)", "(A:1,B:1);garbage", "(A:1,B:1);(C:1,D:1);", "();",
    "(A:1,B:1));", "#NEXUS\nBEGIN TREES;\nEND;", "(A:-1,B:1);",
])
def test_dated_tree_plot_rejects_invalid_inputs_without_replacing_plot(tmp_path, text):
    source, plot = tmp_path / "tree.nwk", tmp_path / "tree.pdf"
    source.write_text(text)
    plot.write_bytes(b"previous plot")
    result = run(sys.executable, SUPPORT / "plot_dated_tree.py", "--infile", source, "--outfile", plot)
    assert result.returncode != 0
    assert plot.read_bytes() == b"previous plot"


@pytest.mark.parametrize("text", [
    "(A:10,B:10);",
    "((A:10,B:10)named:20[&&NHX:age=10:age_ci_low=8:age_ci_high=12:age_ci_kind=HPD:age_ci_level=0.95],C:30);",
    "#NEXUS\nBEGIN TREES;\nTREE dated = [&R] (A:10,B:10)[&95%HPD={8,12}];\nEND;",
])
def test_dated_tree_plot_accepts_public_unit_nhx_and_figtree(tmp_path, text):
    source, plot = tmp_path / "tree.nwk", tmp_path / "tree.pdf"
    source.write_text(text)
    result = run(sys.executable, SUPPORT / "plot_dated_tree.py", "--infile", source, "--outfile", plot)
    assert result.returncode == 0, result.stdout + result.stderr
    assert plot.read_bytes().startswith(b"%PDF-")
    assert b"NWKIT" in plot.read_bytes()
    assert source.read_text() == text


@pytest.mark.parametrize("candidate_method,reason", [
    (None, "mad_without_notung_candidates"),
    ("mad", "mad_compatible_with_notung"),
    ("midpoint", "midpoint_compatible_with_notung"),
    ("outgroup", "first_notung_candidate"),
])
def test_nwkit_root_selection_keeps_notung_priority_and_branch_lengths(tmp_path, candidate_method, reason):
    from nwkit.clade_mapping import projected_root_split
    from nwkit.util import read_tree

    source = tmp_path / "tree.nwk"
    source.write_text("((A:7.8,B:0.949):0.244,(C:7.8,D:1.439):2.87);")
    candidates = tmp_path / "notung"
    candidates.mkdir()
    candidate = candidates / "tree.nwk.rooting.0"
    if candidate_method:
        args = ["nwkit", "root", "--method", candidate_method, "--infile", source, "--outfile", candidate]
        if candidate_method == "outgroup":
            args += ["--outgroup", "D"]
        result = run(*args)
        assert result.returncode == 0, result.stderr
    output, table, plot = tmp_path / "selected.nwk", tmp_path / "roots.tsv", tmp_path / "roots.pdf"
    result = run(sys.executable, SUPPORT / "species_tree_guided_gene_tree_rooting.py",
                 "--in-tree", source, "--notung-root-dir", candidates, "--out-tree", output,
                 "--comparison-table", table, "--comparison-plot", plot)
    assert result.returncode == 0, result.stdout + result.stderr
    assert f"Selected root: {reason}" in result.stdout
    before, after = read_tree(str(source), "auto", True), read_tree(str(output), "auto", True)
    assert set(before.leaf_names()) == set(after.leaf_names()) == set("ABCD")
    for a in "ABCD":
        for b in "ABCD":
            assert before.get_distance(a, b) == pytest.approx(after.get_distance(a, b))
    if candidate_method:
        candidate_tree = read_tree(str(candidate), "auto", True)
        assert projected_root_split(after, frozenset("ABCD")) == projected_root_split(candidate_tree, frozenset("ABCD"))
    assert table.read_text().startswith("method\tstatus\t")
    assert plot.read_bytes().startswith(b"%PDF-")


@pytest.mark.parametrize("mode", ["dna", "pep"])
def test_root_core_stage_invalidates_engine_and_preserves_previous_bundle_on_failure(tmp_path, mode):
    core = (ROOT / "workflow/core/gg_genome_evolution_core.sh").read_text()
    begin = core.index("busco_species_tree_assisted_gene_tree_rooting() {")
    function = core[begin:core.index("\nbusco_grampa()", begin)]
    title = "DNA" if mode == "dna" else "protein"
    begin = core.index(f'task="Species-tree-guided gene tree rooting of duplicate-containing BUSCO {title} trees"')
    stage = core[begin:core.index('\ntask=', begin + 1)]
    original = "((A:7.8,B:0.949):0.244,(C:7.8,D:1.439):2.87);"
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "OG1.busco.nwk").write_text(original)
    archives = tmp_path / "archives"
    archives.mkdir()
    archive = archives / "OG1.busco.notung.root.zip"

    def write_candidate(text):
        with zipfile.ZipFile(archive, "w") as bundle:
            bundle.writestr("OG1.notung.root/OG1.busco.nwk.rooting.0", text)

    write_candidate(original)
    output = tmp_path / "output"
    script = tmp_path / "stage.sh"
    script.write_text(f'''set -euo pipefail
gg_support_dir={shlex.quote(str(SUPPORT))}
source "${{gg_support_dir}}/gg_util.sh"
gg_workspace_dir={shlex.quote(str(tmp_path))}
gg_workspace_output_dir={shlex.quote(str(output))}
genome_evolution_provenance_dir="${{gg_workspace_output_dir}}/provenance"
genome_nwkit_identity="${{1:-test-engine-one}}"
file_dated_species_tree={shlex.quote(str(inputs / "OG1.busco.nwk"))}
dir_busco_iqtree_{mode}={shlex.quote(str(inputs))}
dir_busco_notung_{mode}={shlex.quote(str(archives))}
dir_busco_rooted_txt_{mode}="${{gg_workspace_output_dir}}/reports"
dir_busco_rooted_nwk_{mode}="${{gg_workspace_output_dir}}/trees"
run_busco_dupaware_root_{mode}=1
species_label_parser=taxonomic
GG_GENOME_PARALLEL_JOBS=1
artifact_stale_policy=rebuild
gg_step_start() {{ echo run >> stage-runs.txt; }}
gg_step_skip() {{ echo skip >> stage-runs.txt; }}
{function}
{stage}
''')

    def execute(identity="test-engine-one"):
        return subprocess.run(["bash", str(script), identity], cwd=tmp_path, capture_output=True, text=True)

    for _ in range(2):
        result = execute()
        assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "stage-runs.txt").read_text().splitlines() == ["run", "skip"]
    result = execute("test-engine-two")
    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "stage-runs.txt").read_text().splitlines() == ["run", "skip", "run"]
    saved = {p: p.read_bytes() for folder in (output / "reports", output / "trees") for p in folder.glob("OG1*")}
    assert len(saved) == 4
    manifest = output / "provenance" / f"busco.root_{mode}.json"
    assert json.loads(manifest.read_text())["step"] == f"genome_evolution_busco_root_{mode}"
    saved_manifest = manifest.read_bytes()
    write_candidate("((A:1,B:1):1,(C:1,WrongTip:1):1);")
    failed = execute("test-engine-two")
    assert failed.returncode != 0
    assert {p: p.read_bytes() for p in saved} == saved
    assert manifest.read_bytes() == saved_manifest


@pytest.mark.parametrize("invalid", ["missing_directory", "file_as_directory", "duplicate_tips", "different_tips", "unresolved_root"])
def test_root_selection_rejects_bad_candidates_and_preserves_outputs(tmp_path, invalid):
    source = tmp_path / "tree.nwk"
    source.write_text("((A:7.8,B:0.949):0.244,(C:7.8,D:1.439):2.87);")
    candidates = tmp_path / "notung"
    if invalid == "file_as_directory":
        candidates.write_text("not a directory")
    elif invalid != "missing_directory":
        candidates.mkdir()
        texts = {
            "duplicate_tips": "((A:1,A:1,B:1):1,(C:1,D:1):1);",
            "different_tips": "((A:1,B:1):1,(C:1,WrongTip:1):1);",
            "unresolved_root": "(A:1,B:1,C:1,D:1);",
        }
        (candidates / "tree.nwk.rooting.0").write_text(texts[invalid])
    outputs = [tmp_path / name for name in ("selected.nwk", "roots.tsv", "roots.pdf")]
    for path in outputs:
        path.write_bytes(b"previous output")
    result = run(sys.executable, SUPPORT / "species_tree_guided_gene_tree_rooting.py",
                 "--in-tree", source, "--notung-root-dir", candidates, "--out-tree", outputs[0],
                 "--comparison-table", outputs[1], "--comparison-plot", outputs[2])
    assert result.returncode != 0, result.stdout + result.stderr
    assert all(path.read_bytes() == b"previous output" for path in outputs)


def test_dated_tree_stage_keeps_pdf_pair_on_failed_install_and_tracks_engine(tmp_path):
    core = (ROOT / "workflow/core/gg_genome_evolution_core.sh").read_text()
    begin = core.index('task="Dated species tree plotting"')
    stage = core[begin:core.index('\nremove_empty_subdirs', begin)]
    source = tmp_path / "mcmctree_95CI.nhx"
    source.write_text("(A:10,B:10)Root[&&NHX:age=10:age_ci_low=8:age_ci_high=12:age_ci_kind=HPD:age_ci_level=0.95];")
    output = tmp_path / "output"
    script = tmp_path / "stage.sh"
    script.write_text(f'''set -euo pipefail
gg_support_dir={shlex.quote(str(SUPPORT))}
source "${{gg_support_dir}}/gg_util.sh"
gg_workspace_dir={shlex.quote(str(tmp_path))}
gg_workspace_output_dir={shlex.quote(str(output))}
genome_evolution_provenance_dir="${{gg_workspace_output_dir}}/provenance"
genome_nwkit_identity="${{1:-engine-one}}"
dir_mcmctree2={shlex.quote(str(tmp_path))}
file_mcmctree_dated_nwk={shlex.quote(str(tmp_path / "absent.nwk"))}
file_dated_species_tree={shlex.quote(str(source))}
file_plot_mcmctree_pdf="${{gg_workspace_output_dir}}/dated.pdf"
file_dated_species_tree_pdf="${{gg_workspace_output_dir}}/summary.pdf"
run_plot_mcmctreer=1
artifact_stale_policy=rebuild
gg_step_start() {{ echo run >> runs.txt; }}
gg_step_skip() {{ echo skip >> runs.txt; }}
mv() {{
  if [[ ${{2:-}} == *.gg-stage.* && ${{3:-}} == "${{file_dated_species_tree_pdf}}" && ${{FAIL_INSTALL:-0}} == 1 ]]; then
    return 1
  fi
  command mv "$@"
}}
{stage}
''')
    import os

    def execute(identity="engine-one", fail=False):
        return subprocess.run(["bash", str(script), identity], cwd=tmp_path,
                              env=dict(os.environ, FAIL_INSTALL="1" if fail else "0"), capture_output=True, text=True)

    for _ in range(2):
        result = execute()
        assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "runs.txt").read_text().splitlines() == ["run", "skip"]
    paths = [output / "dated.pdf", output / "summary.pdf"]
    old = [p.read_bytes() for p in paths]
    assert old[0] == old[1]
    failed = execute("engine-two", fail=True)
    assert failed.returncode != 0
    assert [p.read_bytes() for p in paths] == old
    result = execute("engine-two")
    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "runs.txt").read_text().splitlines() == ["run", "skip", "run", "run"]
    assert paths[0].read_bytes() == paths[1].read_bytes()


def test_root_adapter_consumes_real_notung_candidates(tmp_path):
    from nwkit.util import read_tree

    jar = Path("/usr/local/bin/Notung.jar")
    assert jar.is_file(), "Run this test in the GeneGalleon container."
    species = tmp_path / "species.nwk"
    species.write_text("((Genus_alpha:1,Genus_beta:1):1,(Genus_gamma:1,Genus_delta:1):1);")
    source = tmp_path / "gene.nwk"
    source.write_text("((Genus_alpha_g1:0.2,Genus_gamma_g1:0.8):0.3,((Genus_alpha_g2:0.4,Genus_beta_g1:0.1):0.4,Genus_delta_g1:0.2):0.5);")
    candidate_dir = tmp_path / "candidates"
    result = run("java", "-jar", "-Xmx1g", jar, "-s", species, "-g", source,
                 "--root", "--infertransfers", "false", "--treeoutput", "newick", "--log", "--treestats",
                 "--events", "--parsable", "--speciestag", "prefix", "--allopt", "--maxtrees", "1000",
                 "--nolosses", "--outputdir", candidate_dir)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (candidate_dir / "gene.nwk.rooting.0").is_file()
    output, table, plot = tmp_path / "selected.nwk", tmp_path / "roots.tsv", tmp_path / "roots.pdf"
    result = run(sys.executable, SUPPORT / "species_tree_guided_gene_tree_rooting.py",
                 "--in-tree", source, "--notung-root-dir", candidate_dir, "--out-tree", output,
                 "--comparison-table", table, "--comparison-plot", plot)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Selected root: mad_compatible_with_notung" in result.stdout
    before, after = read_tree(str(source), "auto", True), read_tree(str(output), "auto", True)
    for a in before.leaf_names():
        for b in before.leaf_names():
            assert before.get_distance(a, b) == pytest.approx(after.get_distance(a, b))
    assert plot.read_bytes().startswith(b"%PDF-")

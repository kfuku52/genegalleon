"""Execute GeneGalleon's native reconciliation stages with the real NWKIT."""

import importlib.util
import json
import shlex
import subprocess
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / "workflow/support"


def stage(text, title):
    start = text.index(f'task="{title}"')
    return text[start:text.index('\ntask=', start + 1)]


def test_gene_root_reconciliation_and_statistics_use_native_outputs(tmp_path):
    core = (ROOT / "workflow/core/gg_gene_evolution_core.sh").read_text()
    species, gene = tmp_path / "species.nwk", tmp_path / "gene.nwk"
    species.write_text("(A_a:1,(B_b:1,C_c:1):1);")
    gene.write_text("((A_a_1:1,B_b_1:1):1,A_a_2:2);")
    script = tmp_path / "stages.sh"
    script.write_text(f'''set -euo pipefail
source {shlex.quote(str(SUPPORT / "gg_util.sh"))}
gg_workspace_dir={shlex.quote(str(tmp_path))}
dir_output_active="${{gg_workspace_dir}}/output"
gg_workspace_output_dir="${{dir_output_active}}"
file_og_unrooted_tree_analysis={shlex.quote(str(gene))}
species_tree_pruned={shlex.quote(str(species))}
file_og_rooted_tree="${{dir_output_active}}/rooted_tree/OG1_root.nwk"
file_og_rooted_tree_analysis="${{file_og_rooted_tree}}"
file_og_rooted_log="${{dir_output_active}}/root_log/OG1.txt"
file_og_root_candidates="${{dir_output_active}}/root_candidates/OG1_roots.nwk"
file_og_reconciliation="${{dir_output_active}}/reconciliation/OG1_reconciliation.tsv"
og_id=OG1
species_label_parser=taxonomic
species_label_regex=""
species_label_map_tsv=""
reconciliation_duplication_cost=1.5
reconciliation_loss_cost=1
gene_nwkit_identity=native-test
run_tree_root=1
run_reconciliation=1
tree_rooting_method=reconciliation
artifact_stale_policy=rebuild
gg_step_start() {{ echo "$1" >> ran.txt; }}
gg_step_skip() {{ :; }}
java() {{ echo 'Java must not be called' >&2; return 99; }}
{stage(core, "Gene tree rooting")}
{stage(core, "NWKIT reconciliation")}
''')
    for _ in range(2):
        result = subprocess.run(["bash", str(script)], cwd=tmp_path, text=True, capture_output=True)
        assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "ran.txt").read_text().splitlines() == ["Gene tree rooting", "NWKIT reconciliation"]
    table_path = tmp_path / "output/reconciliation/OG1_reconciliation.tsv"
    table = pd.read_csv(table_path, sep="\t")
    assert table.event_source.eq("lca").all()
    assert table.implied_losses.sum() == 1
    spec = importlib.util.spec_from_file_location("orthogroup_stats_native", SUPPORT / "orthogroup_statistics.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.read_reconciliation_stats(table_path) == {
        "reconciliation_num_dup": 1,
        "reconciliation_num_speciation": 1,
        "reconciliation_num_loss": 1,
    }
    assert not list((tmp_path / "output").rglob("*notung*"))


@pytest.mark.parametrize("mode", ["dna", "pep"])
def test_busco_candidate_stage_preserves_all_ties_and_previous_results_on_failure(tmp_path, mode):
    from nwkit.util import read_tree_strings

    core = (ROOT / "workflow/core/gg_genome_evolution_core.sh").read_text()
    start = core.index("busco_reconciliation() {")
    function = core[start:core.index("\nbusco_species_tree_assisted_gene_tree_rooting()", start)]
    inputs = tmp_path / "input"
    inputs.mkdir()
    gene = inputs / "BUSCO1.busco.nwk"
    gene.write_text("((A_a_1:1,A_a_2:2):3,(A_a_3:4,A_a_4:5):6);")
    species = tmp_path / "species.nwk"
    species.write_text("(A_a:1,B_b:1);")
    title = "DNA" if mode == "dna" else "protein"
    script = tmp_path / "candidates.sh"
    script.write_text(f'''set -euo pipefail
source {shlex.quote(str(SUPPORT / "gg_util.sh"))}
gg_workspace_dir={shlex.quote(str(tmp_path))}
gg_workspace_output_dir="${{gg_workspace_dir}}/output"
genome_evolution_provenance_dir="${{gg_workspace_output_dir}}/provenance"
dir_busco_iqtree_{mode}={shlex.quote(str(inputs))}
dir_busco_reconciliation_{mode}="${{gg_workspace_output_dir}}/candidates"
file_dated_species_tree={shlex.quote(str(species))}
run_busco_dupaware_reconciliation_root_{mode}=1
species_label_parser=taxonomic
species_label_regex=""
species_label_map_tsv=""
reconciliation_duplication_cost=1.5
reconciliation_loss_cost=1
genome_nwkit_identity=native-test
GG_GENOME_PARALLEL_JOBS=1
artifact_stale_policy=rebuild
gg_step_start() {{ echo run >> ran.txt; }}
gg_step_skip() {{ echo skip >> ran.txt; }}
java() {{ echo 'Java must not be called' >&2; return 99; }}
{function}
{stage(core, f"NWKIT reconciliation rooting of duplicate-containing BUSCO {title} trees")}
''')
    for _ in range(2):
        result = subprocess.run(["bash", str(script)], cwd=tmp_path, text=True, capture_output=True)
        assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "ran.txt").read_text().splitlines() == ["run", "skip"]
    candidates = tmp_path / "output/candidates/BUSCO1.busco.roots.nwk"
    assert len(read_tree_strings(str(candidates))) == 5
    manifest = tmp_path / f"output/provenance/busco.reconciliation_{mode}.json"
    assert json.loads(manifest.read_text())["step"] == f"genome_evolution_busco_reconciliation_{mode}"
    previous = candidates.read_bytes(), manifest.read_bytes()
    gene.write_text("(Unknown_x:1,A_a_1:1);")
    failed = subprocess.run(["bash", str(script)], cwd=tmp_path, text=True, capture_output=True)
    assert failed.returncode != 0
    assert (candidates.read_bytes(), manifest.read_bytes()) == previous


REFERENCE = json.loads((ROOT / "workflow/tests/data/nwkit_reconciliation_reference.json").read_text())


@pytest.mark.parametrize("case", REFERENCE["results"], ids=lambda case: str(case["fixture"]))
def test_native_results_match_recorded_notung_reference(tmp_path, case):
    from nwkit.reconcile import build_reconciliation_table
    from nwkit.root import reconciliation_rooting
    from nwkit.util import read_tree

    gene_path, species_path = tmp_path / "gene.nwk", tmp_path / "species.nwk"
    gene_path.write_text(case["gene_tree"])
    species_path.write_text(case["species_tree"])
    gene, species = read_tree(str(gene_path), "auto", True), read_tree(str(species_path), "auto", True)
    mapping = {name: "_".join(name.split("_")[:2]) for name in gene.leaf_names()}
    table = build_reconciliation_table(gene, species, mapping)
    assert int(table.event_type.eq("duplication").sum()) == case["duplications"]
    assert int(table.implied_losses.sum()) == case["losses"]
    _, evaluation = reconciliation_rooting(gene, species, mapping, duplication_cost=1.5, _return_evaluation=True)
    assert sorted(sorted(candidate.split[0]) for candidate in evaluation.candidates) == case["optimal_root_splits"]


@pytest.mark.parametrize("invalid_candidates", [None, "", " \n"])
def test_query_extraction_reuses_all_roots_and_reconciles_extracted_tree(tmp_path, invalid_candidates):
    from nwkit.util import read_tree

    core = (ROOT / "workflow/core/gg_gene_evolution_core.sh").read_text()
    gene, species = tmp_path / "gene.nwk", tmp_path / "species.nwk"
    gene.write_text("((A_a_1:1,A_a_2:2):3,(A_a_3:4,A_a_4:5):6);")
    species.write_text("(A_a:1,B_b:1);")
    query, alignment = tmp_path / "query.txt", tmp_path / "alignment.fa"
    query.write_text("A_a_1\nA_a_2\n")
    alignment.write_text("".join(f">A_a_{i}\nMAAAA\n" for i in range(1, 5)))
    selected, candidates = tmp_path / "selected.nwk", tmp_path / "candidates.nwk"
    subprocess.run(["nwkit", "root", "--method", "reconciliation", "--infile", str(gene),
                    "--species-tree", str(species), "--outfile", str(selected),
                    "--candidates-out", str(candidates)], check=True, capture_output=True)
    start = core.index("set_analysis_file() {")
    setter = core[start:core.index("\nset_default_analysis_files()", start)]
    script = tmp_path / "extract.sh"
    script.write_text(f'''set -euo pipefail
source {shlex.quote(str(SUPPORT / "gg_util.sh"))}
gg_workspace_dir={shlex.quote(str(tmp_path))}
dir_output_active="${{gg_workspace_dir}}/output"
gg_workspace_output_dir="${{dir_output_active}}"
file_query_gene={shlex.quote(str(query))}
file_og_trimmed_aln_analysis={shlex.quote(str(alignment))}
file_og_iqtree_tree={shlex.quote(str(gene))}
file_og_rooted_tree={shlex.quote(str(selected))}
file_og_rooted_tree_analysis="${{file_og_rooted_tree}}"
file_og_root_candidates={shlex.quote(str(candidates))}
species_tree_pruned={shlex.quote(str(species))}
file_og_query_blast="${{gg_workspace_dir}}/absent.tsv"
file_og_orthogroup_extraction_nwk="${{dir_output_active}}/extracted/OG1.nwk"
file_og_orthogroup_extraction_rooted_nwk="${{dir_output_active}}/extracted_root/OG1.nwk"
file_og_orthogroup_extraction_fasta="${{dir_output_active}}/extracted_fasta/OG1.fa.gz"
file_og_reconciliation="${{dir_output_active}}/reconciliation/OG1_reconciliation.tsv"
mode_gene_evolution=query2family
run_orthogroup_extraction=1
run_reconciliation=1
tree_rooting_method=reconciliation
og_id=OG1
input_sequence_mode=protein
genetic_code=1
GG_TASK_CPUS=1
species_label_parser=taxonomic
species_label_regex=""
species_label_map_tsv=""
gene_nwkit_identity=test-query
artifact_stale_policy=rebuild
gg_step_start() {{ echo "$1" >> ran.txt; }}
gg_step_skip() {{ :; }}
{setter}
{stage(core, "Orthogroup extraction with NWKIT")}
{stage(core, "NWKIT reconciliation")}
''')
    for _ in range(2):
        result = subprocess.run(["bash", str(script)], cwd=tmp_path, text=True, capture_output=True)
        assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "ran.txt").read_text().splitlines() == ["Orthogroup extraction with NWKIT", "NWKIT reconciliation"]
    assert len(pd.read_csv(tmp_path / "tmp_num_leaf.tsv", sep="\t")) == 5
    extracted = read_tree(str(tmp_path / "output/extracted_root/OG1.nwk"), "auto", True)
    assert set(extracted.leaf_names()) == {"A_a_1", "A_a_2"}
    table = pd.read_csv(tmp_path / "output/reconciliation/OG1_reconciliation.tsv", sep="\t")
    assert set(table.loc[table.event_type.eq("leaf"), "gene_name"]) == set(extracted.leaf_names())

    previous = {path: path.read_bytes() for path in (tmp_path / "output").rglob("*") if path.is_file()}
    if invalid_candidates is None:
        candidates.unlink()
    else:
        candidates.write_text(invalid_candidates)
    failed = subprocess.run(["bash", str(script)], cwd=tmp_path, text=True, capture_output=True)
    assert failed.returncode != 0
    assert {path: path.read_bytes() for path in previous} == previous

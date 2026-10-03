import json
import os
import random
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio.Data import CodonTable

from workflow.support.wgd_evidence import branch_for_ks, read_table
from workflow.support.wgd_ssd import classify, normalized_counts, species_tree

ROOT = Path(__file__).resolve().parents[2]
CORE = ROOT / "workflow/core/gg_genome_evolution_core.sh"
HELPER = ROOT / "workflow/support/wgd_ssd.py"
NAMES = ("Alpha_one", "Beta_two", "Gamma_three")


def fixture(workspace, genomes=True):
    tree_dir = workspace / "output/species_tree/species_tree_summary"
    tree_dir.mkdir(parents=True)
    (tree_dir / "undated_species_tree.nwk").write_text("((Alpha_one:0.5,Beta_two:0.5):0.5,Gamma_three:1);\n")
    orthogroups = workspace / "output/orthofinder/Orthogroups"
    orthogroups.mkdir(parents=True)
    count_rows, member_rows = [], []
    for i in range(8):
        count_rows.append(f"OG{i:03d}\t{2 if i < 4 else 1}\t1\t1\t{4 if i < 4 else 3}\n")
        a = f"Alpha_one_ChrA_g{i}" + (f", Alpha_one_ChrB_g{i}" if i < 4 else "")
        member_rows.append(f"OG{i:03d}\t{a}\tBeta_two_ChrA_g{i}\tGamma_three_ChrA_g{i}\n")
    (orthogroups / "Orthogroups.GeneCount.tsv").write_text("Orthogroup\t" + "\t".join(NAMES) + "\tTotal\n" + "".join(count_rows))
    (orthogroups / "Orthogroups.tsv").write_text("Orthogroup\t" + "\t".join(NAMES) + "\n" + "".join(member_rows))
    if not genomes:
        return
    codons = {}
    for codon, aa in CodonTable.standard_dna_table.forward_table.items():
        codons.setdefault(aa, []).append(codon)
    manifest = ["species\tmode\n"]
    for species_index, name in enumerate(NAMES):
        fasta, gff = [], ["##gff-version 3"]
        for chromosome in (("ChrA", "ChrB") if species_index == 0 else ("ChrA",)):
            for i in range(4 if chromosome == "ChrB" else 8):
                rng = random.Random(2130 + i)
                protein = list("M" + "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(180)))
                if chromosome == "ChrB":
                    for j in range(12, len(protein), 12):
                        protein[j] = "A" if protein[j] != "A" else "V"
                dna = "".join(codons[aa][(j + species_index + (chromosome == "ChrB")) % len(codons[aa])]
                              for j, aa in enumerate(protein))
                identifier = f"{chromosome}_g{i}"
                fasta.append(f">{identifier}\n{dna}\n")
                start = 1000 * i + 1
                gff.append(f"{chromosome}\ttest\tmRNA\t{start}\t{start + 600}\t.\t+\t.\tID={identifier};Parent=locus_{identifier}")
        for directory in ("species_cds", "species_gff"):
            (workspace / "input" / directory).mkdir(parents=True, exist_ok=True)
        (workspace / "input/species_cds" / f"{name}.cds.fa").write_text("".join(fasta))
        (workspace / "input/species_gff" / f"{name}.gff3").write_text("\n".join(gff) + "\n")
        manifest.append(f"{name}\tcds\n")
    (workspace / "input/wgd_genomes.tsv").write_text("".join(manifest))


def run_core(workspace, **overrides):
    env = {**os.environ, "gg_workspace_dir": str(workspace), "GG_COMMON_TMP_ROOT": "workspace",
           "GG_ARRAY_TASK_ID": "1", "GG_JOB_ID": "wgd_test", "GG_TASK_CPUS": "1", "GG_MEM_PER_CPU_GB": "8",
           "genome_evolution_mode": "wgd", "wgd_count_bootstrap": "0", "wgd_ks_bootstrap": "20",
           "wgd_max_states": "128", "artifact_stale_policy": "stop", **overrides}
    return subprocess.run(["bash", str(CORE)], cwd=ROOT, env=env, capture_output=True, text=True, timeout=900)


def test_real_native_count_synteny_ds_and_node_classification_with_cache(tmp_path):
    workspace = tmp_path / "workspace"
    fixture(workspace)
    sentinel = workspace / "output/orthofinder/keep.txt"
    sentinel.write_text("keep existing orthofinder outputs")
    result = run_core(workspace)
    assert result.returncode == 0, result.stdout + result.stderr
    root = workspace / "output/genome_evolution/wgd_ssd"
    summary = json.loads((root / "summary.json").read_text())
    assert summary["num_raw_anchor_rows"] >= 4
    assert summary["num_ks_contrasts"] > 0
    assert summary["ks_boundary_ci_methods"] == ["pair-median-bonferroni"]
    boundaries = read_table(root / "ks_boundaries.tsv")
    assert all(row["ci_method"] == "pair-median-bonferroni" for row in boundaries)
    # Four complete single-copy families cannot give finite simultaneous
    # 95% marginal-median bounds. Successful estimation is not supported age.
    assert all(row["interval_status"] != "ok" for row in boundaries)
    assert all(row["placement_status"] != "interval_supported" for row in read_table(root / "anchor_evidence.tsv"))
    assert {"ks_boundaries.tsv", "family_ks_contrasts.tsv", "pair_selection.tsv"}.issubset(summary["output_hashes"])
    assert summary["species_synteny"]["Beta_two"]["num_blocks"] == 0
    assert all(row["event_support"] == "unresolved" for row in read_table(root / "wgd_events.tsv"))
    assert (root / "wgd_candidates.pdf").stat().st_size > 1000
    assert sentinel.read_text() == "keep existing orthofinder outputs"
    mtime = (root / "summary.json").stat().st_mtime_ns
    result = run_core(workspace)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (root / "summary.json").stat().st_mtime_ns == mtime
    tree = tmp_path / "gene.nwk"
    tree.write_text("(Alpha_one_ChrA_g0:0.1,Alpha_one_ChrA_g1:0.1)[&&NHX:D=Y:annotation=retained];\n")
    command = [sys.executable, str(HELPER), "classify", "--gene-tree", str(tree), "--species-tree",
               str(workspace / "output/species_tree/species_tree_summary/undated_species_tree.nwk"),
               "--evidence", str(root), "--family-id", "OGtest", "--output", str(tmp_path / "classification"),
               "--native-tree-likelihood", "1"]
    result = subprocess.run(command, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    rows = read_table(tmp_path / "classification/duplication_origins.tsv")
    assert len(rows) == 1 and rows[0]["classification"] == "SSD-supported"
    assert rows[0]["event_source"] == "nhx"
    assert 0 <= float(rows[0]["native_conditional_wgd_probability"]) <= 1
    assert (tmp_path / "classification/native_tree_likelihood.json").is_file()
    assert (tmp_path / "classification/duplication_origins.pdf").stat().st_size > 1000
    assert "duplication_origin=SSD-supported" in (tmp_path / "classification/classified_gene_tree.nhx").read_text()
    assert "D=Y" in (tmp_path / "classification/classified_gene_tree.nhx").read_text()
    assert "annotation=retained" in (tmp_path / "classification/classified_gene_tree.nhx").read_text()
    (workspace / "input/species_cds/Alpha_one.cds.fa").write_text("changed")
    result = run_core(workspace)
    assert result.returncode != 0
    assert (root / "summary.json").stat().st_mtime_ns == mtime


def test_count_selection_keeps_missing_separate_from_zero(tmp_path):
    workspace = tmp_path / "workspace"
    fixture(workspace, genomes=False)
    counts = workspace / "output/orthofinder/Orthogroups/Orthogroups.GeneCount.tsv"
    counts.write_text("Orthogroup\tAlpha_one\tBeta_two\tGamma_three\n"
                      "keep\tNA\t1\t1\nexclude\t0\t0\t1\nkeep2\t1\t0\t1\n")
    output = tmp_path / "results"
    output.mkdir()
    normalized_counts({"counts": str(counts), "species": list(NAMES), "species_tree": str(workspace / "output/species_tree/species_tree_summary/undated_species_tree.nwk"),
                       "parameters": {"max_count_families": 100, "seed": 1}}, output)
    rows = read_table(output / "counts.tsv")
    assert rows[0]["Alpha_one"] == "NA"
    assert len(rows) == 2
    assert read_table(output / "count_family_selection.tsv")[1]["selection"] == "excluded_root_clade_absence"


def test_boundary_uncertainty_root_and_nonmonotone_values_remain_unresolved(tmp_path):
    from nwkit.clade_index import CladeIndex

    path = tmp_path / "tree.nwk"
    path.write_text("((Alpha_one:0.5,Beta_two:0.5):0.5,Gamma_three:1);")
    tree = species_tree(path)
    index = CladeIndex(tree)
    boundary = index.clade_id_for_node(tree.children[0])
    row = {"status": "ok", "monotone_from_younger_node": "not_comparable", "interval_status": "ok", "ci_lower": "0.9", "ci_upper": "1.1"}
    bounds = {("Alpha_one", boundary): row}
    assert branch_for_ks("Alpha_one", 0.5, tree, bounds)[1] == "interval_supported"
    assert branch_for_ks("Alpha_one", 1.0, tree, bounds)[1] == "boundary_overlap"
    assert branch_for_ks("Alpha_one", 1.5, tree, bounds)[0] is None
    assert branch_for_ks("Alpha_one", 0, tree, bounds)[0] is None
    assert branch_for_ks("Alpha_one", 0.5, tree, {("Alpha_one", boundary): {**row, "monotone_from_younger_node": "no"}})[0] is None


def test_evidence_from_another_full_tree_is_rejected(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("(A:1,B:1);")
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    (evidence / "summary.json").write_text(json.dumps({"plan": {"species_tree": "old.nwk", "input_hashes": {"old.nwk": "wrong"}}}))
    with pytest.raises(ValueError, match="species trees differ"):
        classify(SimpleNamespace(evidence=evidence, output=tmp_path / "results", species_tree=tree))


def test_wgd_ssd_cli_help_in_owned_runtime(tmp_path):
    result = subprocess.run([sys.executable, str(HELPER), "--help"], cwd=tmp_path,
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "usage:" in result.stdout.lower()
    assert list(tmp_path.iterdir()) == []

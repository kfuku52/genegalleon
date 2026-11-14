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
    assert (root / "wgd_species_tree.pdf").stat().st_size > 1000
    assert (root / "wgd_species_tree.png").read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert {"wgd_species_tree.pdf", "wgd_species_tree.png"}.issubset(summary["plot_hashes"])
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


def test_wgd_stage_forwards_observed_and_bootstrap_worker_counts(tmp_path):
    workspace = tmp_path / "workspace"
    fixture(workspace, genomes=False)
    result = run_core(workspace, GG_TASK_CPUS="2", wgd_count_scan_workers="2",
                      wgd_count_bootstrap="2", wgd_count_bootstrap_workers="2",
                      wgd_count_event_fractions="0.5", wgd_count_family_gamma_shape="none",
                      wgd_count_bootstrap_checkpoint="bootstrap.json")
    assert result.returncode == 0, result.stdout + result.stderr
    root = workspace / "output/genome_evolution/wgd_ssd"
    summary = json.loads((root / "summary.json").read_text())
    assert summary["plan"]["parameters"]["count_scan_workers"] == 2
    assert summary["plan"]["parameters"]["count_bootstrap_workers"] == 2
    model = json.loads((root / "count_model.json").read_text())
    checkpoint = workspace / "bootstrap.json"
    assert json.loads(checkpoint.read_text())["statistics"] == model["bootstrap_max_statistics"]
    assert summary["plan"]["count_bootstrap_checkpoint"] == str(checkpoint)
    assert "bootstrap.json" not in summary["output_hashes"]
    assert len(model["bootstrap_max_statistics"]) == 2
    assert len(model["candidates"]) == 4
    assert model["event_fractions"] == [0.5] and model["family_rate_scales"] == [1.0]
    assert all(row["event_fraction"] == "0.5" for row in read_table(root / "count_candidates.tsv"))


def test_wgd_stage_forwards_known_observation_probabilities_and_rate_groups(tmp_path):
    workspace = tmp_path / "workspace"
    fixture(workspace, genomes=False)
    detection = workspace / "detection.tsv"
    detection.write_text("leaf_name\tdetection_probability\nGamma_three\t0.91\nAlpha_one\t0.97\nBeta_two\t1\n")
    rates = workspace / "regimes.tsv"
    rates.write_text("branch_id\tregime\n1\tancestral\n2\tterminal\n3\tterminal\n4\tterminal\n")
    options = {"wgd_count_detection_tsv": detection.name, "wgd_count_rate_groups_tsv": rates.name}
    result = run_core(workspace, **options)
    assert result.returncode == 0, result.stdout + result.stderr
    root = workspace / "output/genome_evolution/wgd_ssd"
    initial = (root / "summary.json").read_bytes()
    plan = json.loads(initial)["plan"]
    assert {str(detection), str(rates)} <= plan["input_hashes"].keys()
    model = json.loads((root / "count_model.json").read_text())
    assert dict(zip(model["tip_names"], model["detection_probabilities"], strict=True)) == dict(zip(NAMES, (.97, 1, .91), strict=True))
    assert model["branch_rate_group_names"] == ["ancestral", "terminal"]
    detection.write_text(detection.read_text().replace("0.91", "0.92"))
    result = run_core(workspace, **options)
    assert result.returncode != 0
    assert (root / "summary.json").read_bytes() == initial


@pytest.mark.parametrize("holdout", [False, True])
def test_wgd_stage_can_plot_duplication_score_without_genomes(tmp_path, holdout):
    from zipfile import ZipFile

    workspace = tmp_path / "workspace"
    fixture(workspace, genomes=False)
    orthogroups = workspace / "output/orthofinder/Orthogroups"
    (orthogroups / "Orthogroups.GeneCount.tsv").write_text(
        "Orthogroup\t" + "\t".join(NAMES) + "\tTotal\nOG000\t2\t2\t1\t5\nOG001\t1\t1\t1\t3\n")
    (orthogroups / "Orthogroups.tsv").write_text(
        "Orthogroup\t" + "\t".join(NAMES) + "\n"
        "OG000\tAlpha_one_a, Alpha_one_b\tBeta_two_a, Beta_two_b\tGamma_three_a\n"
        "OG001\tAlpha_one_c\tBeta_two_c\tGamma_three_b\n")
    archive = workspace / "rooted_tree.zip"
    with ZipFile(archive, "w") as stream:
        stream.writestr("rooted_tree/OG000_root.nwk",
                        "(((Alpha_one_a,Beta_two_a),(Alpha_one_b,Beta_two_b)),Gamma_three_a);")
        stream.writestr("rooted_tree/OG001_root.nwk", "((Alpha_one_c,Beta_two_c),Gamma_three_b);")
    result = run_core(workspace, wgd_gene_trees_file="rooted_tree.zip", wgd_min_duplication_score="0.5",
                      wgd_duplication_taxon_holdout=str(int(holdout)), wgd_duplication_coverage_audit=str(int(holdout)))
    assert result.returncode == 0, result.stdout + result.stderr
    root = workspace / "output/genome_evolution/wgd_ssd"
    summary = json.loads((root / "summary.json").read_text())
    assert summary["candidate_model"] == "rooted_gene_tree_duplication_score"
    assert summary["num_input_gene_trees"] == 2
    row, = read_table(root / "wgd_events.tsv")
    assert row["event_support"] == "duplication-only-candidate"
    assert float(row["duplication_score"]) == 0.5
    assert not (root / "count_model.json").exists()
    assert (root / "wgd_species_tree.png").stat().st_size > 1000
    holdout_path = root / "duplication_taxon_holdout.tsv"
    assert holdout_path.is_file() == holdout
    assert (holdout_path.name in summary["output_hashes"]) == holdout
    coverage = root / "duplication_coverage.tsv"
    assert coverage.is_file() == holdout
    if holdout:
        audited, = read_table(coverage)
        assert audited["eligible_complete"] == "2" and audited["eligible_partial"] == "0"
        assert audited["projected_score"] == audited["complete_case_score"] == "0.5"
        assert audited["partial_case_score"] == ""
        assert coverage.name in summary["output_hashes"]
        estimates = read_table(holdout_path)
        assert {r["excluded_taxon"] for r in estimates} == set(NAMES)
        assert all(r["species_event_id"] == row["species_event_id"] for r in estimates)
        assert all(r["branch_identifiable"] == "0" and r["duplication_score"] == "NA" for r in estimates)


def test_wgd_stage_uses_cross_anchors_in_candidate_plot(tmp_path):
    import shutil
    from zipfile import ZipFile

    workspace = tmp_path / "workspace"
    fixture(workspace, genomes=False)
    manifest = ["species_a\tspecies_b\tbed_a\tbed_b\tanchors\n"]
    a, b, c = NAMES
    for name in (a, b):
        (workspace / f"{name}.bed").write_text("".join(
            f"chr{copy}\t{i * 10}\t{i * 10 + 5}\t{name}_{copy}{i}\n"
            for copy in ("a", "b") for i in range(5)))
        (workspace / f"{name}.anchors").write_text("###\n" + "".join(
            f"{name}_a{i}\t{name}_b{i}\t1\n" for i in range(5)))
        manifest.append(f"{name}\t{name}\t{name}.bed\t{name}.bed\t{name}.anchors\n")
    (workspace / "cross.anchors").write_text("".join("###\n" + "".join(
        f"{a}_{copy}{i}\t{b}_{copy}{i}\t1\n" for i in range(5)) for copy in ("a", "b")))
    manifest.append(f"{a}\t{b}\t{a}.bed\t{b}.bed\tcross.anchors\n")
    (workspace / "cross.tsv").write_text("".join(manifest))
    orthogroups = workspace / "output/orthofinder/Orthogroups"
    (orthogroups / "Orthogroups.GeneCount.tsv").write_text(
        "Orthogroup\t" + "\t".join(NAMES) + "\tTotal\n" + "".join(f"OG{i:03d}\t2\t2\t1\t5\n" for i in range(5)))
    (orthogroups / "Orthogroups.tsv").write_text("Orthogroup\t" + "\t".join(NAMES) + "\n" + "".join(
        f"OG{i:03d}\t{a}_a{i}, {a}_b{i}\t{b}_a{i}, {b}_b{i}\t{c}_a{i}\n" for i in range(5)))
    (workspace / "cross_cds.tsv").write_text(f"species\tcds\n{a}\t{a}.fa\n{b}\t{b}.fa\n")
    for name in (a, b):
        (workspace / f"{name}.fa").write_text("".join(
            f">{name}_{copy}{i}\nATG" + ("GCT" if copy == "a" else "GCA") * 90 + "\n"
            for copy in ("a", "b") for i in range(5)))
    with ZipFile(workspace / "rooted_tree.zip", "w") as stream:
        for i in range(5):
            stream.writestr(f"OG{i:03d}_root.nwk", f"((({a}_a{i},{b}_a{i}),({a}_b{i},{b}_b{i})),{c}_a{i});")
    with ZipFile(workspace / "boot.zip", "w") as stream:
        stream.writestr("OG000_generax.ufboot.nwk", f"(({a}_a0,{b}_a0)100,({a}_b0,{b}_b0)88,{c}_a0);")
    result = run_core(workspace, wgd_gene_trees_file="rooted_tree.zip", wgd_cross_synteny_file="cross.tsv",
                      wgd_min_cross_windows="1", wgd_cross_ks_cds_file="cross_cds.tsv",
                      wgd_cross_tree_support_file="boot.zip", wgd_max_ks_quartets="3", wgd_cross_orientation_audit="1",
                      wgd_max_ks_outgroup_families="11")
    assert result.returncode == 0, result.stdout + result.stderr
    root = workspace / "output/genome_evolution/wgd_ssd"
    summary = json.loads((root / "summary.json").read_text())
    row, = read_table(root / "wgd_events.tsv")
    pair, = read_table(root / "cross_synteny_pairs.tsv")
    assert row["event_support"] == "duplication-and-synteny-candidate"
    assert row["species_event_id"] == pair["species_event_id"]
    assert pair["matched_paired_windows"] == "1" and pair["matched_families"] == "5"
    assert {"cross_synteny_pairs.tsv", "cross_synteny_windows.tsv", "cross_synteny_copy_phase.tsv"} <= summary["output_hashes"].keys()
    phase, = read_table(root / "cross_synteny_copy_phase.tsv")
    assert phase["concordant_only_quartets"] == "5" and phase["concordant_only_windows"] == "1"
    assert phase["ambiguous_quartets"] == phase["discordant_only_quartets"] == "0"
    coverage, = read_table(root / "cross_synteny_archive_coverage.tsv")
    assert coverage["num_cross_anchor_pairs"] == coverage["num_archive_family_anchor_pairs"] == "10"
    assert coverage["num_observable_windows"] == coverage["num_archive_family_windows"] == "2"
    orientation, = read_table(root / "cross_synteny_orientation.tsv")
    assert orientation["status"] == "supports_in_both"
    assert orientation["a_matched_paired_windows"] == orientation["b_matched_paired_windows"] == "1"
    nonreused, = read_table(root / "cross_synteny_anchor_nonreuse.tsv")
    assert nonreused["num_nonreused_anchor_windows"] == nonreused["num_nonreused_archive_family_windows"] == "2"
    region = read_table(root / "cross_synteny_region_summary.tsv")[0]
    assert region["mode"] == "paired" and region["regions"] == "1"
    assert region["num_holdouts_losing_support"] == "1"
    assert "cross_synteny_regions.tsv" in summary["output_hashes"]
    assert read_table(root / "cross_synteny_nested_copies.tsv") == []
    support = {row["family_id"]: row for row in read_table(root / "cross_synteny_tree_support.tsv")}
    assert support["OG000"]["first_copy_support"] == "100.0"
    assert support["OG000"]["second_copy_support"] == "88.0"
    assert support["OG000"]["duplication_support"] == ""
    assert support["OG001"]["status"] == "archive_missing"
    assert str(workspace / "boot.zip") in summary["plan"]["input_hashes"]
    assert row["cross_ks_status"] == "inconclusive"
    assert len(read_table(root / "cross_ks_paralogs.tsv")) == 10
    quartets = read_table(root / "cross_ks_quartets.tsv")
    assert len(quartets) == 3 and {r["species_event_id"] for r in quartets} == {row["species_event_id"]}
    assert "cross_ks_quartets.tsv" in summary["output_hashes"]
    assert summary["plan"]["parameters"]["max_ks_outgroup_families"] == 11
    assert {r["status"] for r in read_table(root / "cross_ks_common_outgroups.tsv")} == {"missing_trio"}
    assert read_table(root / "cross_ks_common_outgroup_pairs.tsv") == []
    assert {"cross_ks_common_outgroups.tsv", "cross_ks_common_outgroup_pairs.tsv"} <= summary["output_hashes"].keys()
    assert all(row["status"] == "missing_trio" for row in read_table(root / "cross_ks_matched_boundaries.tsv"))
    assert all(row["status"] == "insufficient_window_components" for row in read_table(root / "cross_ks_focals.tsv"))
    assert "wgd_cross_ks.png" in summary["plot_hashes"]
    assert str(workspace / "cross.anchors") in summary["plan"]["input_hashes"]
    assert (root / "wgd_species_tree.png").stat().st_size > 1000
    saved = workspace / "prior_wgd"
    saved.mkdir()
    shutil.copy(root / "summary.json", saved)
    shutil.copy(root / "cross_ks_aligned_pairs.code1.tsv", saved)
    result = run_core(workspace, wgd_gene_trees_file="rooted_tree.zip", wgd_cross_synteny_file="cross.tsv",
                      wgd_min_cross_windows="1", wgd_cross_ks_cds_file="cross_cds.tsv",
                      wgd_cross_tree_support_file="boot.zip", wgd_max_ks_quartets="3", wgd_cross_orientation_audit="1",
                      wgd_max_ks_outgroup_families="11", artifact_stale_policy="rebuild",
                      wgd_cross_ks_reuse_summary="prior_wgd/summary.json")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "reusing 22/22 codon alignments" in result.stdout
    reused = json.loads((root / "summary.json").read_text())
    assert reused["output_hashes"] == summary["output_hashes"]
    assert str(saved / "cross_ks_aligned_pairs.code1.tsv") in reused["plan"]["input_hashes"]


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

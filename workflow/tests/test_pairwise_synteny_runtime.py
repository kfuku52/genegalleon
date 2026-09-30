import csv
import json
import os
import random
import subprocess
import sys
from pathlib import Path

from Bio.Data import CodonTable
from PIL import Image

from workflow.tests.test_genome_evolution_protein_mode import _run_core

REPO_ROOT = Path(__file__).resolve().parents[2]
CORE = REPO_ROOT / "workflow/core/gg_genome_evolution_core.sh"


def test_pairwise_synteny_help_uses_real_runtime_dependency(tmp_path):
    script = REPO_ROOT / "workflow/support/pairwise_synteny.py"
    result = subprocess.run([sys.executable, str(script), "--help"], cwd=tmp_path,
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "usage:" in result.stdout.lower()
    assert list(tmp_path.iterdir()) == []


def write_genome(workspace, species, mode, chromosomes):
    sequence_dir = workspace / "input" / f"species_{mode}"
    annotation_dir = workspace / "input/species_gff"
    sequence_dir.mkdir(parents=True, exist_ok=True)
    annotation_dir.mkdir(parents=True, exist_ok=True)
    fasta = []
    gff = ["##gff-version 3"]
    codons = {}
    for codon, aa in CodonTable.standard_dna_table.forward_table.items():
        codons.setdefault(aa, codon)
    for chromosome, reverse in chromosomes:
        indices = list(range(8))
        if reverse:
            indices.reverse()
        for position, index in enumerate(indices):
            rng = random.Random(1729 + index)
            protein = "M" + "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(160))
            identifier = f"{chromosome}_g{index}"
            fasta_identifier = species + "_" + identifier if mode == "protein" else identifier
            sequence = protein if mode == "protein" else "".join(codons[aa] for aa in protein)
            fasta.extend((f">{fasta_identifier}", sequence))
            start = position * 1000 + 1
            gff.extend((f"{chromosome}\ttest\tgene\t{start}\t{start + 600}\t.\t+\t.\tID=locus_{identifier}",
                        f"{chromosome}\ttest\tmRNA\t{start}\t{start + 600}\t.\t+\t.\tID={identifier};Parent=locus_{identifier}"))
    (sequence_dir / f"{species}.{mode}.fa").write_text("\n".join(fasta) + "\n")
    (annotation_dir / f"{species}.gff3").write_text("\n".join(gff) + "\n")


def run_core(workspace, **overrides):
    env = {**os.environ, "gg_workspace_dir": str(workspace), "GG_COMMON_TMP_ROOT": "workspace",
           "GG_ARRAY_TASK_ID": "1", "GG_JOB_ID": "synteny_test", "GG_TASK_CPUS": "1", "GG_MEM_PER_CPU_GB": "8",
           "genome_evolution_mode": "synteny", "artifact_stale_policy": "stop", **overrides}
    return subprocess.run(["bash", str(CORE)], cwd=REPO_ROOT, env=env, capture_output=True, text=True, timeout=180)


def test_synteny_only_generates_real_plots_preserves_other_stages_and_reuses_analysis(tmp_path):
    workspace = tmp_path / "workspace"
    write_genome(workspace, "Triphyophyllum_peltatum", "protein", (("Chr1", False),))
    write_genome(workspace, "Ancistrocladus_abbreviatus", "cds", (("Chr2", False), ("Chr10", True)))
    pair_table = workspace / "input/synteny_pairs.tsv"
    pair_table.write_text("analysis_id\ttarget_species\tquery_species\ntriphyophyllum_ancistrocladus\tTriphyophyllum_peltatum\tAncistrocladus_abbreviatus\n")
    for stage in ("species_tree", "orthofinder"):
        path = workspace / "output" / stage
        path.mkdir(parents=True)
        (path / "user-output").write_text("preserve this\n")
    result = run_core(workspace)
    assert result.returncode == 0, result.stdout + result.stderr
    root = workspace / "output/genome_evolution/synteny"
    analysis = root / "analysis/triphyophyllum_ancistrocladus"
    plots = root / "plots/triphyophyllum_ancistrocladus"
    summary = json.loads((analysis / "summary.json").read_text())
    assert summary["syntenic_genes"] == [8, 16]
    assert summary["parameters"]["quota"] is None
    assert summary["target"]["source"]["mode"] == "protein"
    assert summary["query"]["source"]["mode"] == "cds"
    with (analysis / "blocks.tsv").open() as handle:
        blocks = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["orientation"] for row in blocks} == {"+", "-"}
    assert {row["query_seqid"] for row in blocks} == {"Chr2", "Chr10"}
    for name in ("dotplot", "karyotype"):
        assert (plots / f"{name}.pdf").read_bytes().startswith(b"%PDF")
        svg = (plots / f"{name}.svg").read_text()
        assert "<svg" in svg
        assert all(f"<!-- {seqid} -->" in svg for seqid in ("Chr1", "Chr2", "Chr10"))
        with Image.open(plots / f"{name}.png") as image:
            image.verify()
    assert (plots / "seqids").read_text() == "Chr1\nChr2,Chr10\n"
    for stage in ("species_tree", "orthofinder"):
        assert sorted(p.name for p in (workspace / "output" / stage).iterdir()) == ["user-output"]
    analysis_mtime = (analysis / "commands.json").stat().st_mtime_ns
    image_mtime = (plots / "karyotype.png").stat().st_mtime_ns
    result = run_core(workspace)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (analysis / "commands.json").stat().st_mtime_ns == analysis_mtime
    assert (plots / "karyotype.png").stat().st_mtime_ns == image_mtime
    pair_table.write_text("analysis_id\ttarget_species\tquery_species\tquery_seqids\ntriphyophyllum_ancistrocladus\tTriphyophyllum_peltatum\tAncistrocladus_abbreviatus\tChr10,Chr2\n")
    result = run_core(workspace, synteny_plot_only="1", synteny_plot_formats="png", artifact_stale_policy="rebuild")
    assert result.returncode == 0, result.stdout + result.stderr
    assert (analysis / "commands.json").stat().st_mtime_ns == analysis_mtime
    assert (plots / "seqids").read_text() == "Chr1\nChr10,Chr2\n"
    assert (plots / "karyotype.png").stat().st_mtime_ns != image_mtime
    source = workspace / "input/species_protein/Triphyophyllum_peltatum.protein.fa"
    source.write_text(source.read_text().replace("\nM", "\nA", 1))
    result = run_core(workspace)
    assert result.returncode != 0
    assert (analysis / "commands.json").stat().st_mtime_ns == analysis_mtime
    result = run_core(workspace, synteny_plot_only="1", artifact_stale_policy="rebuild")
    assert result.returncode != 0
    assert "plot-only requires" in result.stderr
    assert (analysis / "commands.json").stat().st_mtime_ns == analysis_mtime
    result = run_core(workspace, synteny_plot_only="1", artifact_stale_policy="reuse")
    assert result.returncode != 0
    assert (analysis / "commands.json").stat().st_mtime_ns == analysis_mtime


def test_plot_only_without_completed_analysis_does_not_publish(tmp_path):
    workspace = tmp_path / "workspace"
    write_genome(workspace, "Target_species", "protein", (("Chr1", False),))
    write_genome(workspace, "Query_species", "protein", (("Chr2", False),))
    (workspace / "input/synteny_pairs.tsv").write_text("analysis_id\ttarget_species\tquery_species\npair\tTarget_species\tQuery_species\n")
    result = run_core(workspace, synteny_plot_only="1")
    assert result.returncode != 0
    assert "plot-only requires" in result.stderr
    assert not (workspace / "output/genome_evolution/synteny/analysis").exists()


def test_normal_genome_evolution_can_opt_in_to_pairwise_stage(tmp_path):
    workspace = tmp_path / "workspace"
    write_genome(workspace, "Target_species", "protein", (("Chr1", False),))
    write_genome(workspace, "Query_species", "protein", (("Chr2", False),))
    (workspace / "input/synteny_pairs.tsv").write_text("analysis_id\ttarget_species\tquery_species\npair\tTarget_species\tQuery_species\n")
    result = _run_core(tmp_path, {"genome_evolution_mode": "all", "run_pairwise_synteny": "1",
                                 "synteny_plot_formats": "png", "run_orthofinder": "0"})
    assert result.returncode == 0, result.stdout + result.stderr
    root = workspace / "output/genome_evolution/synteny"
    assert (root / "analysis/pair/summary.json").is_file()
    assert (root / "plots/pair/karyotype.png").is_file()

import json
import os
import random
import shutil
import subprocess
import zipfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
CORE = REPO_ROOT / "workflow" / "core" / "gg_fractionation_bias_core.sh"


def _gene_sequence(index: int) -> str:
    generator = random.Random(index + 1729)
    return "ATG" + "".join(generator.choice("ACGT") for _ in range(300)) + "TAA"


def _write_genome(directory: Path, species: str, prefix: str, seqids: tuple[str, ...]) -> None:
    cds_path = directory / "species_cds" / f"{species}.cds.fa"
    gff_path = directory / "species_gff" / f"{species}.gff3"
    cds_path.parent.mkdir(parents=True, exist_ok=True)
    gff_path.parent.mkdir(parents=True, exist_ok=True)
    fasta_lines: list[str] = []
    gff_lines = ["##gff-version 3"]
    for copy_index, seqid in enumerate(seqids):
        for gene_index in range(8):
            gene_id = f"{prefix}{copy_index + 1}_{gene_index + 1}"
            sequence = _gene_sequence(gene_index)
            fasta_lines.extend((f">{gene_id}", sequence))
            start = gene_index * 1000 + 1
            end = start + len(sequence) - 1
            gff_lines.append(f"{seqid}\ttest\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gene_id}")
    cds_path.write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")
    gff_path.write_text("\n".join(gff_lines) + "\n", encoding="utf-8")


def _write_self_genome(directory: Path, species: str) -> None:
    cds_path = directory / "species_cds" / f"{species}.cds.fa"
    gff_path = directory / "species_gff" / f"{species}.gff3"
    cds_path.parent.mkdir(parents=True, exist_ok=True)
    gff_path.parent.mkdir(parents=True, exist_ok=True)
    fasta_lines: list[str] = []
    gff_lines = ["##gff-version 3"]
    for copy_index, seqid in enumerate(("homeolog_a", "homeolog_b")):
        for gene_index in range(8):
            gene_id = f"self{copy_index + 1}_{gene_index + 1}"
            sequence = _gene_sequence(gene_index)
            if copy_index == 1:
                bases = list(sequence)
                for position in range(12 + gene_index % 3, len(bases) - 3, 20):
                    bases[position] = {"A": "C", "C": "G", "G": "T", "T": "A"}[bases[position]]
                sequence = "".join(bases)
            fasta_lines.extend((f">{gene_id}", sequence))
            start = gene_index * 1000 + 1
            end = start + len(sequence) - 1
            gff_lines.append(f"{seqid}\ttest\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gene_id}")
    cds_path.write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")
    gff_path.write_text("\n".join(gff_lines) + "\n", encoding="utf-8")


@pytest.mark.skipif(
    os.environ.get("KFFRACTBIAS_RUN_INTEGRATION") != "1",
    reason="set KFFRACTBIAS_RUN_INTEGRATION=1 to run the JCVI/LAST workflow integration",
)
def test_fractionation_bias_core_runs_pair_and_publishes_bundle(tmp_path: Path) -> None:
    assert shutil.which("kffractbias"), "kffractbias must be installed in the GeneGalleon runtime"
    workspace = tmp_path / "workspace"
    input_dir = workspace / "input"
    _write_genome(input_dir, "Target_species", "t", ("target_chr",))
    _write_genome(input_dir, "Query_species", "q", ("query_a", "query_b"))
    pair_table = input_dir / "fractionation_bias_pairs.tsv"
    pair_table.write_text(
        "analysis_id\ttarget_species\tquery_species\tquota\twindow_size\n"
        "synthetic\tTarget_species\tQuery_species\t1:2\t4\n",
        encoding="utf-8",
    )

    env = os.environ.copy()
    env.update(
        {
            "gg_workspace_dir": str(workspace),
            "GG_ARRAY_TASK_ID": "1",
            "GG_JOB_ID": "test",
            "GG_TASK_CPUS": "1",
            "GG_MEM_PER_CPU_GB": "8",
            "run_kffractbias": "1",
            "delete_tmp_dir": "1",
            "artifact_stale_policy": "rebuild",
        }
    )
    completed = subprocess.run(
        ["bash", str(CORE)],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
        timeout=180,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr

    result_dir = workspace / "output" / "kffractbias" / "synthetic"
    expected = {
        "synthetic.genes.tsv",
        "synthetic.windows.tsv",
        "synthetic.summary.json",
        "synthetic.plot.pdf",
        "synthetic.plot.png",
        "synthetic.synteny.zip",
    }
    # Persistent publication guards are coordination state, not tool outputs.
    assert {
        path.name for path in result_dir.iterdir()
        if not path.name.endswith((".gg-bundle.lock.guard", ".gg-bundle.lock.guard.namespace-v1"))
    } == expected
    assert all((result_dir / filename).stat().st_size > 0 for filename in expected)

    summary = json.loads((result_dir / "synthetic.summary.json").read_text(encoding="utf-8"))
    assert summary["program"] == "kfFractBias"
    assert summary["metadata"]["synteny_generation"]["quota"] == "1:2"
    assert summary["outputs"]["genes"] == str(result_dir / "synthetic.genes.tsv")
    assert summary["outputs"]["synteny_archive"] == str(result_dir / "synthetic.synteny.zip")
    for label in ("synteny", "target_bed", "query_bed"):
        assert summary["inputs"][label]["path"].startswith(f"{result_dir / 'synthetic.synteny.zip'}::")
    assert summary["inputs"]["source_target_cds"]["path"] == str(
        input_dir / "species_cds" / "Target_species.cds.fa"
    )
    assert summary["inputs"]["source_target_gff"]["path"] == str(
        input_dir / "species_gff" / "Target_species.gff3"
    )
    assert summary["inputs"]["source_query_cds"]["path"] == str(
        input_dir / "species_cds" / "Query_species.cds.fa"
    )
    assert summary["inputs"]["source_query_gff"]["path"] == str(
        input_dir / "species_gff" / "Query_species.gff3"
    )
    with zipfile.ZipFile(result_dir / "synthetic.synteny.zip") as archive:
        assert "synthetic.synteny/target.query.lifted.1x2.anchors" in archive.namelist()

    manifest = workspace / "output" / "artifact_provenance" / "kffractbias" / "synthetic.json"
    assert manifest.is_file()
    assert not list((workspace / "output" / "tmp" / "kffractbias").glob("1_synthetic.*"))


def test_fractionation_bias_pair_example_has_required_columns() -> None:
    header = (REPO_ROOT / "workspace" / "input" / "fractionation_bias_pairs.example.tsv").read_text(
        encoding="utf-8"
    ).splitlines()[0].split("\t")
    assert header[:4] == ["analysis_id", "target_species", "query_species", "quota"]
    assert {"mode", "window_size", "denominator", "aligner", "minimum_mapping_fraction"} <= set(header)


def test_self_mode_is_array_executed_and_genome_evolution_aggregated() -> None:
    fractionation_core = CORE.read_text(encoding="utf-8")
    genome_evolution_core = (
        REPO_ROOT / "workflow" / "core" / "gg_genome_evolution_core.sh"
    ).read_text(encoding="utf-8")
    genome_evolution_entrypoint = (
        REPO_ROOT / "workflow" / "gg_genome_evolution_entrypoint.sh"
    ).read_text(encoding="utf-8")
    assert 'if values["mode"] not in {"compare", "self"}' in fractionation_core
    assert "selfcompare" in fractionation_core
    assert 'dir_result="${gg_workspace_output_dir}/genome_evolution/self_fractionation_bias/${analysis_id}"' in fractionation_core
    assert "run_self_fractionation_bias=0" in genome_evolution_entrypoint
    assert 'task="Summarizing kfFractBias self-synteny retention across species"' in genome_evolution_core
    assert 'python "${gg_support_dir}/summarize_self_fractionation_bias.py"' in genome_evolution_core


@pytest.mark.skipif(
    os.environ.get("KFFRACTBIAS_RUN_INTEGRATION") != "1",
    reason="set KFFRACTBIAS_RUN_INTEGRATION=1 to run the JCVI/LAST workflow integration",
)
def test_fractionation_bias_core_runs_self_mode_into_genome_evolution(tmp_path: Path) -> None:
    assert shutil.which("kffractbias"), "kffractbias must be installed in the GeneGalleon runtime"
    workspace = tmp_path / "workspace"
    input_dir = workspace / "input"
    _write_self_genome(input_dir, "Self_species")
    pair_table = input_dir / "fractionation_bias_pairs.tsv"
    pair_table.write_text(
        "analysis_id\ttarget_species\tquery_species\tquota\tmode\twindow_size\tdiagonal_bound\n"
        "synthetic_self\tSelf_species\tSelf_species\t1:1\tself\t4\t1\n",
        encoding="utf-8",
    )
    env = os.environ.copy()
    env.update(
        {
            "gg_workspace_dir": str(workspace),
            "GG_ARRAY_TASK_ID": "1",
            "GG_JOB_ID": "test",
            "GG_TASK_CPUS": "1",
            "GG_MEM_PER_CPU_GB": "8",
            "run_kffractbias": "1",
            "delete_tmp_dir": "1",
            "artifact_stale_policy": "rebuild",
        }
    )
    completed = subprocess.run(
        ["bash", str(CORE)],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
        timeout=180,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    result_dir = workspace / "output" / "genome_evolution" / "self_fractionation_bias" / "synthetic_self"
    summary = json.loads((result_dir / "synthetic_self.summary.json").read_text(encoding="utf-8"))
    assert summary["analysis_mode"] == "self_synteny_retention"
    assert summary["genegalleon"]["analysis_mode"] == "self"
    assert summary["genegalleon"]["quota"] == "1:1"
    assert summary["counts"]["removed_identity_pair_count"] == 0
    with zipfile.ZipFile(result_dir / "synthetic_self.synteny.zip") as archive:
        assert "synthetic_self.synteny/self.self.lifted.1x1.anchors" in archive.namelist()

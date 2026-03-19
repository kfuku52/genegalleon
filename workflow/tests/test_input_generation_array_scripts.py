from pathlib import Path
import csv
import gzip
import json
import subprocess
import sys


SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
PLAN_SCRIPT = SUPPORT_DIR / "plan_input_generation_tasks.py"
RUN_TASK_SCRIPT = SUPPORT_DIR / "run_input_generation_task.py"
MERGE_SCRIPT = SUPPORT_DIR / "merge_input_generation_shards.py"


def run_python(script: Path, *args):
    return subprocess.run(
        [sys.executable, str(script), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def write_direct_species_fixture(root: Path, species_name: str) -> None:
    species_dir = root / species_name
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / f"{species_name}.cds.fa").write_text(
        f">{species_name}.gene1.t1\nATGAAATTT\n",
        encoding="utf-8",
    )
    (species_dir / f"{species_name}.gff").write_text(
        "\n".join(
            [
                "##gff-version 3",
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=gene1.t1;Parent=gene1",
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=gene1.t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    (species_dir / f"{species_name}.genome.fa").write_text(
        ">chr1\nATGAAATTT\n",
        encoding="utf-8",
    )


def test_plan_input_generation_tasks_discovers_direct_species(tmp_path: Path):
    input_root = tmp_path / "Direct" / "species_wise_original"
    write_direct_species_fixture(input_root, "Arabidopsis_thaliana")
    write_direct_species_fixture(input_root, "Oryza_sativa")
    task_plan = tmp_path / "task_plan.json"

    completed = run_python(
        PLAN_SCRIPT,
        "--provider",
        "direct",
        "--input-dir",
        str(input_root),
        "--outfile",
        str(task_plan),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    payload = json.loads(task_plan.read_text(encoding="utf-8"))
    assert payload["task_count"] == 2
    assert payload["species"] == ["Arabidopsis_thaliana", "Oryza_sativa"]
    assert payload["tasks"][0]["provider"] == "direct"
    assert payload["tasks"][0]["gene_grouping_mode"] == "rescue_overlap"


def test_run_input_generation_task_and_merge_shards(tmp_path: Path):
    input_root = tmp_path / "Direct" / "species_wise_original"
    write_direct_species_fixture(input_root, "Arabidopsis_thaliana")
    write_direct_species_fixture(input_root, "Oryza_sativa")
    task_plan = tmp_path / "task_plan.json"
    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    shard_dir = tmp_path / "species_summary_shards"
    stats_dir = tmp_path / "task_stats_shards"
    meta_dir = tmp_path / "task_meta_shards"
    shard_dir.mkdir(parents=True, exist_ok=True)
    stats_dir.mkdir(parents=True, exist_ok=True)
    meta_dir.mkdir(parents=True, exist_ok=True)

    completed = run_python(
        PLAN_SCRIPT,
        "--provider",
        "direct",
        "--input-dir",
        str(input_root),
        "--outfile",
        str(task_plan),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    for task_index in (1, 2):
        completed = run_python(
            RUN_TASK_SCRIPT,
            "--task-plan",
            str(task_plan),
            "--task-index",
            str(task_index),
            "--species-cds-dir",
            str(out_cds),
            "--species-gff-dir",
            str(out_gff),
            "--species-genome-dir",
            str(out_genome),
            "--species-summary-output",
            str(shard_dir / f"{task_index}.tsv"),
            "--stats-output",
            str(stats_dir / f"{task_index}.json"),
            "--task-meta-output",
            str(meta_dir / f"{task_index}.json"),
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    aggregate_stats = tmp_path / "aggregate_stats.json"
    merged_species_summary = tmp_path / "gg_input_generation_species.tsv"
    completed = run_python(
        MERGE_SCRIPT,
        "--species-summary-shard-dir",
        str(shard_dir),
        "--species-summary-output",
        str(merged_species_summary),
        "--task-stats-dir",
        str(stats_dir),
        "--aggregate-stats-output",
        str(aggregate_stats),
        "--expected-task-count",
        "2",
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    cds_outputs = sorted(out_cds.glob("*.fa.gz"))
    assert len(cds_outputs) == 2
    with gzip.open(cds_outputs[0], "rt", encoding="utf-8") as handle:
        assert handle.read().startswith(">")

    with open(merged_species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 2
    assert {row["species_prefix"] for row in rows} == {"Arabidopsis_thaliana", "Oryza_sativa"}

    payload = json.loads(aggregate_stats.read_text(encoding="utf-8"))
    assert payload["task_stats_files"] == 2
    assert payload["num_species_cds_files"] == 2
    assert payload["num_species_gff_files"] == 2
    assert payload["num_species_genome_files"] == 2

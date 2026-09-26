import csv
import gzip
import json
import subprocess
import sys
from pathlib import Path

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
PLAN_SCRIPT = SUPPORT_DIR / "plan_input_generation_tasks.py"
RUN_TASK_SCRIPT = SUPPORT_DIR / "run_input_generation_task.py"
MERGE_SCRIPT = SUPPORT_DIR / "merge_input_generation_shards.py"
STAGE_SCRIPT = SUPPORT_DIR / "stage_input_generation_downloads.py"
REQUIRE_GENOMES_SCRIPT = SUPPORT_DIR / "validate_required_genomes.py"


def test_staged_http_inputs_run_without_server_and_reject_missing_or_changed_cache(tmp_path):
    import functools
    import threading
    from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

    raw = tmp_path / "raw"
    species = "Arabidopsis_thaliana"
    write_direct_species_fixture(raw, species)
    server = ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(SimpleHTTPRequestHandler, directory=str(raw)))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    manifest = tmp_path / "source.tsv"
    base = f"http://127.0.0.1:{server.server_port}/{species}/{species}"
    manifest.write_text("provider\tid\tspecies_key\tcds_url\tgff_url\tgenome_url\n"
                        f"direct\tfixture\t{species}\t{base}.cds.fa\t{base}.gff\t{base}.genome.fa\n")
    plan = tmp_path / "plan.json"
    planned = run_python(PLAN_SCRIPT, "--provider", "all", "--download-manifest", str(manifest),
                         "--download-dir", str(tmp_path / "downloads"), "--stage-downloads", "--outfile", str(plan))
    assert planned.returncode == 0, planned.stderr
    args = ("--task-plan", str(plan), "--task-index", "1", "--species-cds-dir", str(tmp_path / "cds"),
            "--species-gff-dir", str(tmp_path / "gff"), "--species-genome-dir", str(tmp_path / "genome"))
    try:
        missing = run_python(RUN_TASK_SCRIPT, *args)
        assert missing.returncode != 0 and "Staged download receipt is missing" in missing.stderr
        staged = run_python(STAGE_SCRIPT, "--task-plan", str(plan), "--jobs", "3")
        assert staged.returncode == 0, staged.stdout + staged.stderr
    finally:
        server.shutdown()
        server.server_close()
        thread.join(3)
    manifest.write_text("no longer valid\n")
    completed = run_python(RUN_TASK_SCRIPT, *args)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert list((tmp_path / "cds").glob("*.fa.gz"))
    restaged = run_python(STAGE_SCRIPT, "--task-plan", str(plan))
    assert restaged.returncode == 0 and "no downloads needed" in restaged.stdout, restaged.stderr
    cached = json.loads(Path(str(plan) + ".tasks/1.json").read_text())["task"]
    Path(cached["cds_path"]).write_text(">changed\nATG\n")
    rejected = run_python(RUN_TASK_SCRIPT, *args)
    assert rejected.returncode != 0 and "Raw input changed" in rejected.stderr
    rejected_stage = run_python(STAGE_SCRIPT, "--task-plan", str(plan))
    assert rejected_stage.returncode != 0 and "Staged raw input changed" in rejected_stage.stderr


def test_prepare_resources_are_separate_from_compute_array(tmp_path):
    helper = SUPPORT_DIR.parent / "gg_input_generation_array.py"
    result = run_python(helper, "--task-plan", str(tmp_path / "plan.json"), "--cpus", "4", "--memory", "32G",
                        "--prepare-cpus", "8", "--prepare-memory", "8G", "--partition", "compute",
                        "--prepare-partition", "network")
    assert result.returncode == 0, result.stderr
    prepare = next(line for line in result.stdout.splitlines() if "MODE=array_prepare " in line)
    worker = next(line for line in result.stdout.splitlines() if "MODE=array_worker " in line)
    assert "--cpus-per-task=8" in prepare and "--mem=8G" in prepare and "--partition=network" in prepare
    assert "--cpus-per-task=4" in worker and "--mem=32G" in worker and "--partition=compute" in worker


def test_required_genome_is_opt_in_for_local_array_planning(tmp_path):
    source = tmp_path / "Direct" / "species_wise_original"
    write_direct_species_fixture(source, "Arabidopsis_thaliana")
    (source / "Arabidopsis_thaliana" / "Arabidopsis_thaliana.genome.fa").unlink()
    args = ("--provider", "direct", "--input-dir", str(source), "--outfile", str(tmp_path / "plan.json"))
    default = run_python(PLAN_SCRIPT, *args)
    assert default.returncode == 0, default.stderr
    required = run_python(PLAN_SCRIPT, *args, "--require-genome")
    assert required.returncode != 0
    assert "Required genome input is missing for Arabidopsis_thaliana" in required.stderr


def test_required_genome_rejects_partial_staged_download_without_receipt(tmp_path):
    import functools
    import threading
    from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

    raw = tmp_path / "raw"
    species = "Arabidopsis_thaliana"
    write_direct_species_fixture(raw, species)
    (raw / species / f"{species}.genome.fa").unlink()
    server = ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(SimpleHTTPRequestHandler, directory=str(raw)))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        manifest = tmp_path / "manifest.tsv"
        base = f"http://127.0.0.1:{server.server_port}/{species}/{species}"
        manifest.write_text("provider\tid\tspecies_key\tcds_url\tgff_url\tgenome_url\n"
                            f"direct\tfixture\t{species}\t{base}.cds.fa\t{base}.gff\t{base}.genome.fa\n")
        plan = tmp_path / "plan.json"
        planned = run_python(PLAN_SCRIPT, "--provider", "all", "--download-manifest", str(manifest),
                             "--download-dir", str(tmp_path / "downloads"), "--stage-downloads", "--outfile", str(plan))
        assert planned.returncode == 0, planned.stderr
        required = run_python(STAGE_SCRIPT, "--task-plan", str(plan), "--require-genome")
        assert required.returncode != 0
        assert "Required genome input is missing for Arabidopsis_thaliana" in required.stderr
        assert not Path(str(plan) + ".tasks/1.json").exists()
        optional = run_python(STAGE_SCRIPT, "--task-plan", str(plan))
        assert optional.returncode == 0, optional.stderr
        assert json.loads(Path(str(plan) + ".tasks/1.json").read_text())["task"]["genome_path"] is None
        required_cached = run_python(STAGE_SCRIPT, "--task-plan", str(plan), "--require-genome")
        assert required_cached.returncode != 0
        assert "Required genome input is missing for Arabidopsis_thaliana" in required_cached.stderr
    finally:
        server.shutdown()
        server.server_close()
        thread.join(3)


def test_required_genome_summary_rejects_missing_or_empty_output(tmp_path):
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nATG\n")
    summary = tmp_path / "species.tsv"
    summary.write_text("species_prefix\tgenome_output_path\nArabidopsis_thaliana\t" + str(genome) + "\n")
    ok = run_python(REQUIRE_GENOMES_SCRIPT, "--species-summary", str(summary), "--expected-task-count", "1")
    assert ok.returncode == 0, ok.stderr
    genome.write_text("")
    missing = run_python(REQUIRE_GENOMES_SCRIPT, "--species-summary", str(summary), "--expected-task-count", "1")
    assert missing.returncode != 0
    assert "Required formatted genome is missing" in missing.stderr
    genome.write_text("not a FASTA\nATG\n")
    invalid = run_python(REQUIRE_GENOMES_SCRIPT, "--species-summary", str(summary))
    assert invalid.returncode != 0
    genome_gz = tmp_path / "genome.fa.gz"
    with gzip.open(genome_gz, "wt") as handle:
        handle.write(">chr1\nATG\n")
    summary.write_text("species_prefix\tgenome_output_path\nArabidopsis_thaliana\t" + str(genome_gz) + "\n")
    compressed = run_python(REQUIRE_GENOMES_SCRIPT, "--species-summary", str(summary))
    assert compressed.returncode == 0, compressed.stderr


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
    assert payload["tasks"][0]["gff_repair_mode"] == "safe"
    assert payload["tasks"][0]["format_strict"] is False


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

    (stats_dir / "1.mapping.json").write_text(json.dumps({
        "phase_conflicts_total": 3,
        "utr_conflicts_total": 1,
    }), encoding="utf-8")

    aggregate_stats = tmp_path / "aggregate_stats.json"
    mapping_qc = tmp_path / "species_mapping_qc.tsv"
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
        "--mapping-qc-output", str(mapping_qc),
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
    assert payload["mapping_qc_files"] == 1
    assert payload["phase_conflicts_total"] == 3
    assert payload["utr_conflicts_total"] == 1
    assert payload["num_species_cds_files"] == 2
    assert payload["num_species_gff_files"] == 2
    assert payload["num_species_genome_files"] == 2
    with mapping_qc.open(newline="") as handle:
        qc_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [(row["species_prefix"], row["qc_status"]) for row in qc_rows] == [
        ("Arabidopsis_thaliana", "available"), ("Oryza_sativa", "not_recorded")
    ]
    assert qc_rows[0]["phase_conflicts_total"] == "3"
    assert payload["cds_gff_records_mapped"] == 2
    assert payload["cds_gff_records_unmapped"] == 0


def test_merge_requires_mapping_qc_to_be_bound_to_worker_receipt(tmp_path: Path):
    raw = tmp_path / "Direct" / "species_wise_original"
    write_direct_species_fixture(raw, "Arabidopsis_thaliana")
    plan = tmp_path / "plan.json"
    state = SUPPORT_DIR / "input_generation_array_state.py"
    assert run_python(PLAN_SCRIPT, "--provider", "direct", "--input-dir", str(raw),
                      "--outfile", str(plan)).returncode == 0
    shard_dir = tmp_path / "species_summary_shards"
    stats_dir = tmp_path / "task_stats_shards"
    shard_dir.mkdir()
    stats_dir.mkdir()
    summary = shard_dir / "1.tsv"
    stats = stats_dir / "1.json"
    generated = run_python(
        RUN_TASK_SCRIPT, "--task-plan", str(plan), "--task-index", "1",
        "--species-cds-dir", str(tmp_path / "cds"),
        "--species-gff-dir", str(tmp_path / "gff"),
        "--species-genome-dir", str(tmp_path / "genome"),
        "--species-summary-output", str(summary), "--stats-output", str(stats),
    )
    assert generated.returncode == 0, generated.stderr
    receipt_args = ("complete", "--task-plan", str(plan), "--task-index", "1",
                    "--file", str(summary), "--file", str(stats))
    assert run_python(state, *receipt_args).returncode == 0
    qc = stats_dir / "1.mapping.json"
    qc.write_text(json.dumps({"phase_conflicts_total": 2, "utr_conflicts_total": 1}), encoding="utf-8")
    merge_args = (
        "--species-summary-shard-dir", str(shard_dir),
        "--species-summary-output", str(tmp_path / "merged.tsv"),
        "--task-stats-dir", str(stats_dir),
        "--aggregate-stats-output", str(tmp_path / "aggregate.json"),
        "--expected-task-count", "1", "--task-plan", str(plan),
    )
    unbound = run_python(MERGE_SCRIPT, *merge_args)
    assert unbound.returncode != 0
    assert "Mapping QC is not bound" in unbound.stderr
    assert run_python(state, *receipt_args, "--file", str(qc)).returncode == 0
    merged = run_python(MERGE_SCRIPT, *merge_args)
    assert merged.returncode == 0, merged.stderr
    qc.write_text(json.dumps({"phase_conflicts_total": 0}), encoding="utf-8")
    assert run_python(MERGE_SCRIPT, *merge_args).returncode != 0


def test_manifest_planning_defers_downloads_and_freezes_inputs(tmp_path):
    source = tmp_path / "sources"
    write_direct_species_fixture(source, "Arabidopsis_thaliana")
    manifest = tmp_path / "manifest.tsv"
    species = "Arabidopsis_thaliana"
    raw = source / species
    fields = ["provider", "id", "species_key", "cds_url", "gff_url", "genome_url"]
    row = ["direct", "fixture", species, (raw / (species + ".cds.fa")).as_uri(),
           (raw / (species + ".gff")).as_uri(), (raw / (species + ".genome.fa")).as_uri()]
    manifest.write_text("\t".join(fields) + "\n" + "\t".join(row) + "\n")
    plan = tmp_path / "plan.json"
    downloads = tmp_path / "downloads"
    args = ["--provider", "all", "--download-manifest", str(manifest), "--download-dir", str(downloads), "--outfile", str(plan)]
    prepared = run_python(PLAN_SCRIPT, *args)
    assert prepared.returncode == 0, prepared.stderr
    assert not downloads.exists()
    frozen = plan.read_bytes()
    manifest.unlink()  # The worker must consume the frozen row, not mutable input.
    worker = run_python(RUN_TASK_SCRIPT, "--task-plan", str(plan), "--task-index", "1",
                        "--species-cds-dir", str(tmp_path / "cds"), "--species-gff-dir", str(tmp_path / "gff"),
                        "--species-genome-dir", str(tmp_path / "genome"), "--describe-only",
                        "--task-meta-output", str(tmp_path / "meta.json"))
    assert worker.returncode == 0, worker.stderr + worker.stdout
    assert downloads.exists()
    meta = json.loads((tmp_path / "meta.json").read_text())
    assert meta["species_prefix"] == species
    assert Path(meta["cds_path"]).is_file()
    assert plan.read_bytes() == frozen


def test_duplicate_species_and_changed_plan_are_rejected(tmp_path):
    raw = tmp_path / "raw"
    write_direct_species_fixture(raw, "Arabidopsis_thaliana")
    plan = tmp_path / "plan.json"
    args = ["--provider", "direct", "--input-dir", str(raw), "--outfile", str(plan)]
    assert run_python(PLAN_SCRIPT, *args).returncode == 0
    frozen = plan.read_bytes()
    assert run_python(PLAN_SCRIPT, *args).returncode == 0
    write_direct_species_fixture(raw, "Oryza_sativa")
    assert run_python(PLAN_SCRIPT, *args).returncode != 0
    assert plan.read_bytes() == frozen
    manifest = tmp_path / "duplicate.tsv"
    manifest.write_text("provider\tid\tspecies_key\n" + "direct\tx\tArabidopsis_thaliana\n" * 2)
    result = run_python(PLAN_SCRIPT, "--provider", "all", "--download-manifest", str(manifest),
                        "--download-dir", str(tmp_path / "downloads"), "--outfile", str(tmp_path / "bad.json"))
    assert result.returncode != 0
    assert "Duplicate species" in result.stderr


def test_receipts_detect_changed_outputs_and_retry_selects_only_pending(tmp_path):
    state = SUPPORT_DIR / "input_generation_array_state.py"
    raw = tmp_path / "raw"
    write_direct_species_fixture(raw, "Arabidopsis_thaliana")
    write_direct_species_fixture(raw, "Oryza_sativa")
    plan = tmp_path / "plan.json"
    assert run_python(PLAN_SCRIPT, "--provider", "direct", "--input-dir", str(raw), "--outfile", str(plan)).returncode == 0
    assert run_python(state, "configure", "--task-plan", str(plan), "--prepare").returncode == 0
    assert run_python(state, "prepared", "--task-plan", str(plan)).returncode == 0
    output = tmp_path / "result"
    output.write_text("valid")
    assert run_python(state, "complete", "--task-plan", str(plan), "--task-index", "1", "--file", str(output)).returncode == 0
    pending = run_python(state, "pending", "--task-plan", str(plan))
    assert pending.stdout.strip() == "2"
    helper = SUPPORT_DIR.parent / "gg_input_generation_array.py"
    preview = run_python(helper, "--task-plan", str(plan), "--retry", "--cpus", "3", "--max-running", "2")
    assert preview.returncode == 0, preview.stderr
    assert "--array=2%2" in preview.stdout
    assert "--dependency=afterok:WORKER_JOB_ID" in preview.stdout
    assert "--cpus-per-task=3" in preview.stdout
    output.write_text("corrupt")
    assert run_python(state, "pending", "--task-plan", str(plan)).stdout.strip() == "1,2"


def test_retry_submission_refuses_any_active_legacy_worker_array(tmp_path, monkeypatch):
    import os
    state = SUPPORT_DIR / "input_generation_array_state.py"
    raw = tmp_path / "raw"
    write_direct_species_fixture(raw, "Arabidopsis_thaliana")
    plan = tmp_path / "plan.json"
    assert run_python(PLAN_SCRIPT, "--provider", "direct", "--input-dir", str(raw), "--outfile", str(plan)).returncode == 0
    assert run_python(state, "configure", "--task-plan", str(plan), "--prepare").returncode == 0
    assert run_python(state, "prepared", "--task-plan", str(plan)).returncode == 0
    fake = tmp_path / "squeue"
    fake.write_text("#!/bin/sh\necho 40457_29\n")
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path) + os.pathsep + os.environ["PATH"])
    helper = SUPPORT_DIR.parent / "gg_input_generation_array.py"
    result = run_python(helper, "--task-plan", str(plan), "--retry", "--submit")
    assert result.returncode != 0
    assert "40457_29" in result.stderr


def test_slurm_helper_waits_for_prepare_and_submits_afterok(tmp_path, monkeypatch):
    import os
    plan = tmp_path / "plan.json"
    calls = tmp_path / "calls.jsonl"
    fake = tmp_path / "sbatch"
    fake.write_text("#!" + sys.executable + "\n" + '''import hashlib, json, os, sys
from pathlib import Path
log = Path(os.environ["TEST_CALLS"])
mode = os.environ["GG_INPUT_INPUT_GENERATION_MODE"]
with log.open("a") as handle:
    handle.write(json.dumps({"mode": mode, "argv": sys.argv[1:]}) + "\\n")
if mode == "array_prepare":
    if os.environ.get("TEST_PREPARE_FAIL"):
        print("simulated Slurm rejection", file=sys.stderr)
        sys.exit(1)
    Path(os.environ["GG_INPUT_TASK_PLAN_OUTPUT"]).write_text(json.dumps({"task_count": 2, "tasks": [{"species_prefix": "A_b"}, {"species_prefix": "C_d"}]}))
    plan = Path(os.environ["GG_INPUT_TASK_PLAN_OUTPUT"])
    Path(str(plan) + ".settings.json").write_text("{}")
    Path(str(plan) + ".prepared.json").write_text(json.dumps({"plan_sha256": hashlib.sha256(plan.read_bytes()).hexdigest(), "settings_sha256": hashlib.sha256(b"{}").hexdigest()}))
print({"array_prepare": "101", "array_worker": "102", "array_finalize": "103"}[mode])
''')
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path) + os.pathsep + os.environ["PATH"])
    monkeypatch.setenv("TEST_CALLS", str(calls))
    helper = SUPPORT_DIR.parent / "gg_input_generation_array.py"
    result = run_python(helper, "--task-plan", str(plan), "--submit", "--max-running", "5")
    assert result.returncode == 0, result.stderr
    submissions = [json.loads(line) for line in calls.read_text().splitlines()]
    assert [call["mode"] for call in submissions] == ["array_prepare", "array_worker", "array_finalize"]
    assert "--wait" in submissions[0]["argv"]
    assert "--array=1-2%5" in submissions[1]["argv"]
    assert "--dependency=afterok:102" in submissions[2]["argv"]
    assert not any("--mem-per-cpu" in arg for call in submissions for arg in call["argv"])
    assert all(any(arg.startswith("--wrap=") for arg in call["argv"]) for call in submissions)
    calls.unlink()
    monkeypatch.setenv("TEST_PREPARE_FAIL", "1")
    result = run_python(helper, "--task-plan", str(plan), "--submit")
    assert result.returncode != 0
    assert len(calls.read_text().splitlines()) == 1
    assert "simulated Slurm rejection" in result.stderr


def test_workspace_rejects_second_plan_and_malformed_receipt_is_pending(tmp_path):
    state = SUPPORT_DIR / "input_generation_array_state.py"
    raw = tmp_path / "raw"
    write_direct_species_fixture(raw, "Arabidopsis_thaliana")
    first = tmp_path / "first.json"
    second = tmp_path / "second.json"
    for plan in (first, second):
        assert run_python(PLAN_SCRIPT, "--provider", "direct", "--input-dir", str(raw), "--outfile", str(plan)).returncode == 0
    workspace = tmp_path / "workspace"
    claim = ["claim-workspace", "--workspace", str(workspace), "--prepare", "--task-plan"]
    assert run_python(state, *claim, str(first)).returncode == 0
    assert run_python(state, *claim, str(second)).returncode != 0
    assert run_python(state, "index", "--task-plan", str(first), "--task-index", "01").stdout.strip() == "1"
    output = tmp_path / "output"
    output.write_text("ok")
    assert run_python(state, "complete", "--task-plan", str(first), "--task-index", "1", "--file", str(output)).returncode == 0
    receipt = Path(str(first) + ".completed") / "1.json"
    content = json.loads(receipt.read_text())
    content["files"] = ["bad structure"]
    receipt.write_text(json.dumps(content))
    pending = run_python(state, "pending", "--task-plan", str(first))
    assert pending.returncode == 0, pending.stderr
    assert pending.stdout.strip() == "1"


def test_manifest_source_change_and_foreign_resolved_cache_are_rejected(tmp_path):
    species = "Arabidopsis_thaliana"
    raw_root = tmp_path / "raw"
    write_direct_species_fixture(raw_root, species)
    raw = raw_root / species
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("provider\tid\tspecies_key\nlocal\t" + str(raw) + "\t" + species + "\n")
    plan = tmp_path / "plan.json"
    downloads = tmp_path / "downloads"
    assert run_python(PLAN_SCRIPT, "--provider", "local", "--download-manifest", str(manifest),
                      "--download-dir", str(downloads), "--outfile", str(plan)).returncode == 0
    worker_args = ["--task-plan", str(plan), "--task-index", "1", "--species-cds-dir", str(tmp_path / "cds"),
                   "--species-gff-dir", str(tmp_path / "gff"), "--species-genome-dir", str(tmp_path / "genome"), "--describe-only"]
    source = raw / (species + ".cds.fa")
    original = source.read_bytes()
    source.write_bytes(original + b"AAA\n")
    changed = run_python(RUN_TASK_SCRIPT, *worker_args)
    assert changed.returncode != 0 and "changed after planning" in changed.stderr
    assert not downloads.exists()
    source.write_bytes(original)
    assert run_python(RUN_TASK_SCRIPT, *worker_args).returncode == 0
    cache = Path(str(plan) + ".tasks") / "1.json"
    content = json.loads(cache.read_text())
    content["plan_sha256"] = "belongs to another plan"
    cache.write_text(json.dumps(content))
    foreign = run_python(RUN_TASK_SCRIPT, *worker_args)
    assert foreign.returncode != 0 and "another plan/task" in foreign.stderr
    assert not (tmp_path / "cds").exists()


def test_completion_rejects_raw_changes_during_worker(tmp_path):
    state = SUPPORT_DIR / "input_generation_array_state.py"
    raw = tmp_path / "raw"
    write_direct_species_fixture(raw, "Arabidopsis_thaliana")
    plan = tmp_path / "plan.json"
    assert run_python(PLAN_SCRIPT, "--provider", "direct", "--input-dir", str(raw), "--outfile", str(plan)).returncode == 0
    source = next(raw.glob("*/*.cds.fa"))
    source.write_text(source.read_text() + "AAA\n")
    output = tmp_path / "output"
    output.write_text("would be stale")
    completed = run_python(state, "complete", "--task-plan", str(plan), "--task-index", "1", "--file", str(output))
    assert completed.returncode != 0
    assert not (Path(str(plan) + ".completed") / "1.json").exists()


def test_custom_output_directory_cannot_be_claimed_by_two_workspaces(tmp_path):
    state = SUPPORT_DIR / "input_generation_array_state.py"
    raw = tmp_path / "raw"
    write_direct_species_fixture(raw, "Arabidopsis_thaliana")
    first = tmp_path / "first.json"
    second = tmp_path / "second.json"
    shared = tmp_path / "custom_species_cds"
    for plan in (first, second):
        assert run_python(PLAN_SCRIPT, "--provider", "direct", "--input-dir", str(raw), "--outfile", str(plan)).returncode == 0
    assert run_python(state, "claim-workspace", "--prepare", "--task-plan", str(first), "--workspace", str(tmp_path / "ws1"), "--file", str(shared)).returncode == 0
    second_claim = run_python(state, "claim-workspace", "--prepare", "--task-plan", str(second), "--workspace", str(tmp_path / "ws2"), "--file", str(shared))
    assert second_claim.returncode != 0 and "another array plan" in second_claim.stderr
    assert not (tmp_path / "ws2" / ".array-plan.json").exists()


def test_prepared_marker_rejects_changed_shared_lineage(tmp_path):
    state = SUPPORT_DIR / "input_generation_array_state.py"
    raw = tmp_path / "raw"
    write_direct_species_fixture(raw, "Arabidopsis_thaliana")
    plan = tmp_path / "plan.json"
    assert run_python(PLAN_SCRIPT, "--provider", "direct", "--input-dir", str(raw), "--outfile", str(plan)).returncode == 0
    assert run_python(state, "configure", "--prepare", "--task-plan", str(plan)).returncode == 0
    lineage = tmp_path / "lineage.txt"
    lineage.write_text("eukaryota_odb12")
    assert run_python(state, "prepared", "--task-plan", str(plan), "--file", str(lineage)).returncode == 0
    assert run_python(state, "check-prepared", "--task-plan", str(plan)).returncode == 0
    lineage.write_text("different_lineage")
    assert run_python(state, "check-prepared", "--task-plan", str(plan)).returncode != 0


def test_nonregular_receipt_input_fails_without_blocking(tmp_path):
    import os
    raw = tmp_path / "raw"
    write_direct_species_fixture(raw, "Arabidopsis_thaliana")
    plan = tmp_path / "plan.json"
    assert run_python(PLAN_SCRIPT, "--provider", "direct", "--input-dir", str(raw), "--outfile", str(plan)).returncode == 0
    fifo = tmp_path / "pipe"
    os.mkfifo(fifo)
    result = subprocess.run([sys.executable, str(SUPPORT_DIR / "input_generation_array_state.py"), "complete",
                             "--task-plan", str(plan), "--task-index", "1", "--file", str(fifo)],
                            capture_output=True, text=True, timeout=3)
    assert result.returncode != 0


def test_array_plan_rejects_escaping_species_and_download_filenames(tmp_path):
    for index, (species, filename) in enumerate((("../Arabidopsis_thaliana", "cds.fa"),
                                               ("Arabidopsis_thaliana", "../../outside.fa"),
                                               (".hidden_species", "cds.fa"))):
        manifest = tmp_path / f"bad{index}.tsv"
        manifest.write_text(f"provider\tid\tspecies_key\tcds_filename\ndirect\tx\t{species}\t{filename}\n")
        output = tmp_path / f"bad{index}.json"
        result = run_python(PLAN_SCRIPT, "--provider", "all", "--download-manifest", str(manifest),
                            "--download-dir", str(tmp_path / "downloads"), "--outfile", str(output))
        assert result.returncode != 0 and "filename components" in result.stderr
        assert not output.exists()

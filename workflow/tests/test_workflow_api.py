"""Compatibility and read-only contracts for the scheduler-neutral public API."""
import hashlib
import json
import os
import shlex
import shutil
import signal
import subprocess
import sys
import time
import zipfile
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
API = SUPPORT / "workflow_api.py"
OBSERVER = SUPPORT / "workflow_observation.py"
PROVENANCE = SUPPORT / "artifact_provenance.py"


def cli(script, *args, env=None, input=None):
    return subprocess.run([sys.executable, str(script), *map(str, args)], text=True,
                          capture_output=True, env=env, input=input, check=False)


def query(*args, env=None):
    result = cli(API, *args, env=env)
    assert result.returncode == 0, result.stdout + result.stderr
    return json.loads(result.stdout)


def snapshot(root):
    return {str(path.relative_to(root)): (path.stat().st_mtime_ns, path.read_bytes() if path.is_file() else None)
            for path in root.rglob("*")}


@pytest.fixture
def project(tmp_path):
    workspace = tmp_path / "workspace"
    root = workspace / "output" / "orthogroup"
    source = root / "rooted_tree" / "OG0001_root.nwk"
    output = root / "stat_branch" / "OG0001_stat.branch.tsv"
    manifest = root / "artifact_provenance" / "OG0001.summary_statistics.json"
    source.parent.mkdir(parents=True)
    output.parent.mkdir(parents=True)
    source.write_text("(A,B);\n")
    output.write_text("branch_id\tnode_name\n0\tA\n")
    argv = ["--workspace-root", str(workspace), "--logical-root", str(root), "--manifest", str(manifest),
            "--step", "summary_statistics", "--family-id", "OG0001", "--input", f"tree={source}",
            "--output", f"table={output}", "--parameter", "mode=a"]
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({"schema": "genegalleon-preflight-plan-v1", "contracts": [argv]}))
    return workspace, root, source, output, manifest, argv, plan


def verify(project, *extra):
    workspace, root, *_ = project
    return query("verify", "--root", root, "--workspace-root", workspace,
                 "--family-id", "OG0001", "--require-step", "summary_statistics", *extra)


def test_capabilities_requires_no_project_or_controller():
    response = query("capabilities")
    assert response["schema"] == "genegalleon-api-v1"
    assert response["requires_kfauto"] is False


def test_old_outputs_remain_reusable_without_adoption_or_cache_writes(project, tmp_path):
    workspace, _, _, _, manifest, _, plan = project
    before = snapshot(workspace)
    env = {**os.environ, "GG_CONTENT_DIGEST_CACHE": str(tmp_path / "forbidden.sqlite3")}
    response = query("preflight", "--plan", plan, env=env)
    assert response["contracts"][0]["state"] == "legacy_reusable"
    assert not response["contracts"][0]["historical_generation_verified"]
    assert not response["blocked"]
    assert snapshot(workspace) == before
    assert not manifest.exists()
    assert not (tmp_path / "forbidden.sqlite3").exists()
    result = verify(project)
    assert result["completion_state"] == "unverified"
    assert result["contracts"][0]["error_code"] == "evidence_missing"
    assert snapshot(workspace) == before


def test_runtime_and_preflight_agree_and_do_not_hide_changed_inputs(project):
    workspace, _, source, _, _, argv, plan = project
    assert cli(PROVENANCE, "record", *argv).returncode == 0
    assert query("preflight", "--plan", plan)["contracts"][0]["state"] == "verified_current"
    source.write_text("((A,B),C);\n")
    before = snapshot(workspace)
    response = query("preflight", "--plan", plan)
    assert response["blocked"]
    assert response["contracts"][0]["runtime_exit_code"] == 3
    assert response["contracts"][0]["error_code"] == "inputs_changed"
    assert snapshot(workspace) == before
    assert cli(PROVENANCE, "needs-run", *argv).returncode == 3
    assert verify(project)["contracts"][0]["error_code"] == "changed_input"


def test_proposed_parameter_changes_and_explicit_rebuild_policy(project):
    _, _, _, _, _, argv, plan = project
    assert cli(PROVENANCE, "record", *argv).returncode == 0
    changed = ["mode=b" if value == "mode=a" else value for value in argv]
    for policy, expected in [("stop", "blocked"), ("rebuild", "needs_run"), ("reuse", "policy_reuse_or_adoption")]:
        plan.write_text(json.dumps({"schema": "genegalleon-preflight-plan-v1",
                                    "contracts": [changed + ["--stale-policy", policy]]}))
        row = query("preflight", "--plan", plan)["contracts"][0]
        assert row["state"] == expected
        assert not row["historical_generation_verified"]
        assert cli(PROVENANCE, "needs-run", *changed, "--stale-policy", policy).returncode == row["runtime_exit_code"]


def test_adoption_is_not_relabelled_as_historical_generation_proof(project):
    _, _, _, _, _, argv, plan = project
    assert cli(PROVENANCE, "needs-run", *argv).returncode == 1
    assert query("preflight", "--plan", plan)["contracts"][0]["state"] == "legacy_reusable"
    assert verify(project)["contracts"][0]["state"] == "legacy_reusable"


def test_verify_survives_archive_and_workspace_relocation(project, tmp_path):
    workspace, root, _, output, _, argv, plan = project
    assert cli(PROVENANCE, "record", *argv).returncode == 0
    assert verify(project)["completion_state"] == "verified_declared_steps"
    count = workspace / "genecount.tsv"
    count.write_text("Orthogroup\tsp1\tTotal\nOG0001\t1\t1\n")
    archived = cli(SUPPORT / "gene_family_output_store.py", "archive-family", "--root", root,
                   "--mode", "orthogroup", "--genecount", count, "--family-id", "OG0001")
    assert archived.returncode == 0, archived.stderr
    assert not output.exists()
    before = snapshot(workspace)
    assert verify(project)["completion_state"] == "verified_declared_steps"
    assert query("preflight", "--plan", plan)["contracts"][0]["state"] == "verified_current"
    assert snapshot(workspace) == before
    moved = tmp_path / "moved"
    shutil.move(str(workspace), moved)
    result = query("verify", "--root", moved / "output/orthogroup", "--workspace-root", moved,
                   "--family-id", "OG0001", "--require-step", "summary_statistics")
    assert result["completion_state"] == "verified_declared_steps"
    before = snapshot(moved)
    assert query("preflight", "--plan", plan, "--workspace-root", moved)["contracts"][0]["state"] == "verified_current"
    assert snapshot(moved) == before


def test_unknown_schema_and_partial_legacy_outputs_are_not_success(project):
    _, _, _, output, manifest, argv, plan = project
    assert cli(PROVENANCE, "record", *argv).returncode == 0
    payload = json.loads(manifest.read_text())
    payload["schema_version"] = 9000
    manifest.write_text(json.dumps(payload))
    assert verify(project)["contracts"][0]["error_code"] == "invalid_manifest"
    manifest.unlink()
    second = output.with_name("second.tsv")
    plan.write_text(json.dumps({"schema": "genegalleon-preflight-plan-v1",
                               "contracts": [argv + ["--output", f"second={second}"]]}))
    assert query("preflight", "--plan", plan)["blocked"]


def test_invalid_plan_returns_one_json_error(tmp_path):
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({"schema": "genegalleon-preflight-plan-v1", "contracts": [["--unknown"]]}))
    result = cli(API, "preflight", "--plan", plan)
    assert result.returncode == 2
    assert json.loads(result.stdout)["error_code"] == "query_unavailable"


def test_status_detects_receipt_declaration_changes(project, tmp_path):
    *_, argv, _plan = project
    directory = tmp_path / "observations"
    assert cli(OBSERVER, "--directory", directory, "--workflow", "gg_gene_evolution", "--",
               sys.executable, PROVENANCE, "record", *argv).returncode == 0
    first = query("status", "--directory", directory)
    receipt = next(directory.glob("*/contract-*.json"))
    payload = json.loads(receipt.read_text())
    payload["argv"] = ["mode=b" if value == "mode=a" else value for value in payload["argv"]]
    receipt.write_text(json.dumps(payload))
    second = query("status", "--directory", directory, "--since", first["next_cursor"])
    assert len(second["records"]) == 1


def test_boolean_receipt_exit_code_cannot_bind_attempt(project, tmp_path):
    *_, argv, _plan = project
    directory = tmp_path / "observations"
    assert cli(OBSERVER, "--directory", directory, "--workflow", "gg_gene_evolution", "--",
               sys.executable, PROVENANCE, "record", *argv).returncode == 0
    receipt = next(directory.glob("*/contract-*.json"))
    payload = json.loads(receipt.read_text())
    payload["exit_code"] = False
    receipt.write_text(json.dumps(payload))
    result = cli(API, "verify", "--root", project[1], "--workspace-root", project[0],
                 "--family-id", "OG0001", "--require-step", "summary_statistics",
                 "--attempt", receipt.parent)
    assert json.loads(result.stdout).get("completion_state") != "verified_declared_steps"


@pytest.mark.parametrize("changed_file", ["input", "output", "manifest"])
def test_verify_rejects_input_changed_after_it_was_hashed(project, changed_file):
    workspace, root, source, output, manifest, argv, _plan = project
    target = {"input": source, "output": output, "manifest": manifest}[changed_file]
    assert cli(PROVENANCE, "record", *argv).returncode == 0
    code = f'''import sys
from pathlib import Path
sys.path.insert(0, {str(SUPPORT)!r})
import workflow_api as api
original = api.provenance.audit_entry_digest
def raced(entry, *args):
    result = original(entry, *args)
    if entry.get("label") == "table":
        Path({str(target)!r}).write_text("changed after input inspection")
    return result
api.provenance.audit_entry_digest = raced
sys.exit(api.main(["verify", "--root", {str(root)!r}, "--workspace-root", {str(workspace)!r},
                   "--family-id", "OG0001", "--require-step", "summary_statistics"]))
'''
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 2
    assert json.loads(result.stdout)["error_code"] == "query_unavailable"


def test_archived_input_verification_reads_bytes_not_only_index(project):
    workspace, root, source, _output, _manifest, argv, plan = project
    archive_root = workspace / "source_store"
    (archive_root / "rooted_tree").mkdir(parents=True)
    shutil.copy2(source, archive_root / "rooted_tree" / source.name)
    index = argv.index("--input")
    argv[index:index + 2] = ["--input-gene-family-artifact", f"tree={archive_root}::rooted_tree::{source.name}"]
    assert cli(PROVENANCE, "record", *argv).returncode == 0
    count = workspace / "genecount.tsv"
    count.write_text("Orthogroup\tsp1\tTotal\nOG0001\t1\t1\n")
    assert cli(SUPPORT / "gene_family_output_store.py", "archive-family", "--root", archive_root,
               "--mode", "orthogroup", "--genecount", count, "--family-id", "OG0001").returncode == 0
    changed = False
    for path in archive_root.rglob("*.zip"):
        with zipfile.ZipFile(path) as archive:
            contents = {name: archive.read(name) for name in archive.namelist()}
        members = [name for name in contents if name.endswith(source.name)]
        if members:
            for name in members:
                contents[name] = b"unexpected replacement bytes\n"
            with zipfile.ZipFile(path, "w") as archive:
                for name, data in contents.items():
                    archive.writestr(name, data)
            changed = True
    assert changed
    assert verify(project)["completion_state"] == "unverified"
    plan.write_text(json.dumps({"schema": "genegalleon-preflight-plan-v1", "contracts": [argv]}))
    assert query("preflight", "--plan", plan)["contracts"][0]["state"] != "verified_current"


def test_optional_fifo_cannot_be_classified_as_absent(project):
    *_, output, _manifest, argv, _plan = project
    optional = output.with_name("optional.tsv")
    assert cli(PROVENANCE, "record", *argv, "--optional-output", f"optional={optional}").returncode == 0
    os.mkfifo(optional)
    assert verify(project)["completion_state"] == "unverified"


def test_reversed_attempt_interval_is_unavailable(tmp_path):
    directory = tmp_path / "observations"
    assert cli(OBSERVER, "--directory", directory, "--workflow", "gg_test", "--", "true").returncode == 0
    path = next(directory.glob("*/run.json"))
    record = json.loads(path.read_text())
    record["finished_at_ns"] = record["started_at_ns"] - 1
    path.write_text(json.dumps(record))
    response = query("status", "--directory", directory)
    assert response["complete"] is False
    assert response["next_cursor"] is None


def test_explicit_workspace_mapping_preserves_exact_attempt_validation(project, tmp_path):
    workspace, _root, _source, _output, _manifest, argv, _plan = project
    directory = tmp_path / "observations"
    assert cli(OBSERVER, "--directory", directory, "--workflow", "gg_gene_evolution", "--",
               sys.executable, PROVENANCE, "record", *argv).returncode == 0
    target = tmp_path / "relocated"
    shutil.copytree(workspace, target)
    base = ["verify", "--root", target / "output/orthogroup", "--workspace-root", target,
            "--family-id", "OG0001", "--require-step", "summary_statistics",
            "--attempt", next(directory.iterdir())]
    assert query(*base)["completion_state"] == "unverified"
    assert query(*base, "--recorded-workspace-root", workspace)["completion_state"] == "verified_declared_steps"
    assert query(*base, "--recorded-workspace-root", tmp_path / "wrong")["completion_state"] == "unverified"
    (target / "output/orthogroup/rooted_tree/OG0001_root.nwk").write_text("changed")
    assert query(*base, "--recorded-workspace-root", workspace)["completion_state"] == "unverified"


def test_explicit_preflight_policy_uses_proposed_policy_without_writes(project):
    workspace, _, source, _, _, argv, plan = project
    assert cli(PROVENANCE, "record", *argv).returncode == 0
    source.write_text("changed")
    before = snapshot(workspace)
    result = query("preflight", "--plan", plan, "--stale-policy", "rebuild")
    assert result["contracts"][0]["state"] == "needs_run"
    assert result["contracts"][0]["stale_policy"] == "rebuild"
    assert snapshot(workspace) == before


@pytest.mark.parametrize("pdf", [b"%PDF-1.7\nfixture\n%%EOF\n", b"broken PDF"])
def test_terminal_profile_validates_even_manifest_matching_content(project, pdf):
    workspace, root, _source, output, _manifest, argv, _plan = project
    tree = root / "stat_tree/OG0001_stat.tree.tsv"
    tree.parent.mkdir()
    tree.write_text("num_leaf\ttotal_branch_length\n2\t0.1\n")
    plot = root / "tree_plot/OG0001_tree_plot.pdf"
    plot.parent.mkdir()
    plot.write_bytes(pdf)
    assert cli(PROVENANCE, "record", *argv, "--output", f"stat_tree={tree}").returncode == 0
    assert cli(PROVENANCE, "record", "--logical-root", root, "--workspace-root", workspace,
               "--manifest", root / "artifact_provenance/OG0001.tree_plot.json", "--step", "tree_plot",
               "--family-id", "OG0001", "--input", f"branch={output}", "--output", f"tree_plot={plot}").returncode == 0
    result = verify(project, "--require-step", "tree_plot", "--profile", "gene-evolution-terminal-v1")
    assert (result["completion_state"] == "verified_declared_steps") is pdf.startswith(b"%PDF-")
    assert (result["terminal_validation"]["state"] == "verified") is pdf.startswith(b"%PDF-")


def test_observer_preserves_streams_status_and_records_only_selected_config(tmp_path):
    directory = tmp_path / "observations"
    script = "echo scientific-output; echo diagnostics >&2; exit 7\n"
    env = {**os.environ, "mode_gene_evolution": "orthogroup", "UNRELATED_SECRET": "never-record-me"}
    result = cli(OBSERVER, "--directory", directory, "--workflow", "gg_gene_evolution",
                 "--config-key", "mode_gene_evolution", "--stdin-script", "--", "bash", "-s", input=script, env=env)
    assert result.returncode == 7
    assert result.stdout == "scientific-output\n"
    assert result.stderr == "diagnostics\n"
    records = query("status", "--directory", directory)
    record = records["records"][0]
    assert record["error_code"] == "unknown"
    assert record["completion_state"] == "unverified"
    runtime = query("runtime", "--attempt", directory / record["attempt_id"])["runtime"]
    assert runtime["effective_config"] == {"mode_gene_evolution": "orthogroup"}
    assert runtime["executed_core_sha256"] == hashlib.sha256(script.encode()).hexdigest()
    assert "never-record-me" not in json.dumps(record)
    before = snapshot(directory)
    response = query("status", "--directory", directory, "--since", records["next_cursor"])
    assert response["records"] == []
    assert snapshot(directory) == before


def test_attempts_are_distinct_and_completion_cannot_cross_attempts(project, tmp_path):
    _, _, _, _, _, argv, _ = project
    directory = tmp_path / "observations"
    first = cli(OBSERVER, "--directory", directory, "--workflow", "gg_gene_evolution", "--",
                sys.executable, PROVENANCE, "record", *argv)
    assert first.returncode == 0, first.stderr
    one = next(directory.iterdir())
    assert verify(project, "--attempt", one)["completion_state"] == "verified_declared_steps"
    assert not query("preflight", "--attempt", one)["blocked"]
    second = cli(OBSERVER, "--directory", directory, "--workflow", "gg_gene_evolution", "--",
                 sys.executable, "-c", "pass")
    assert second.returncode == 0
    two = next(path for path in directory.iterdir() if path != one)
    assert verify(project, "--attempt", two)["completion_state"] == "unverified"


def test_owned_error_boundary_publishes_stable_code(project, tmp_path):
    _, _, source, _, _, argv, _ = project
    assert cli(PROVENANCE, "record", *argv).returncode == 0
    source.write_text("(C,D);\n")
    directory = tmp_path / "observations"
    result = cli(OBSERVER, "--directory", directory, "--workflow", "gg_gene_evolution", "--",
                 sys.executable, PROVENANCE, "needs-run", *argv)
    assert result.returncode == 3
    errors = query("status", "--directory", directory)["records"][0]["error_evidence"]
    assert errors[0]["error_code"] == "inputs_changed"
    assert errors[0]["step"] == "summary_statistics"


def test_status_cursor_detects_same_timestamp_content_changes_and_bad_records(tmp_path):
    directory = tmp_path / "observations"
    assert cli(OBSERVER, "--directory", directory, "--workflow", "gg_test", "--", "true").returncode == 0
    first = query("status", "--directory", directory)
    cursor_file = tmp_path / "cursor.json"
    cursor_file.write_text(json.dumps(first["next_cursor"]))
    assert query("status", "--directory", directory, "--since-file", cursor_file)["records"] == []
    path = next(directory.glob("*/run.json"))
    original_times = path.stat()
    record = json.loads(path.read_text())
    record["exit_code"] = 8
    record["execution_accepted"] = False
    path.write_text(json.dumps(record))
    os.utime(path, ns=(original_times.st_atime_ns, original_times.st_mtime_ns))
    assert len(query("status", "--directory", directory, "--since", first["next_cursor"])["records"]) == 1
    path.write_text("invalid")
    response = query("status", "--directory", directory, "--since", first["next_cursor"])
    assert not response["complete"]
    assert response["next_cursor"] is None
    assert response["removed"] == []


def test_unwritable_observation_does_not_change_child_result(tmp_path):
    directory = tmp_path / "not-directory"
    directory.write_text("preserve")
    result = cli(OBSERVER, "--directory", directory, "--workflow", "gg_test", "--", "bash", "-c", "exit 9")
    assert result.returncode == 9
    assert "observation unavailable" in result.stderr
    assert directory.read_text() == "preserve"


def test_termination_is_forwarded_to_child_subtree(tmp_path):
    directory = tmp_path / "observations"
    ready = tmp_path / "ready"
    child_code = f"from pathlib import Path; import time; Path({str(ready)!r}).write_text('ready'); time.sleep(60)"
    process = subprocess.Popen([sys.executable, str(OBSERVER), "--directory", str(directory),
                                "--workflow", "gg_test", "--", sys.executable, "-c", child_code])
    try:
        deadline = time.monotonic() + 10
        while not ready.exists() and time.monotonic() < deadline:
            time.sleep(0.02)
        assert ready.exists()
        process.send_signal(signal.SIGTERM)
        assert process.wait(timeout=10) == 143
        assert query("status", "--directory", directory)["records"][0]["exit_code"] == 143
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()


def test_shell_bridge_opt_in_and_resolved_setting_capture(tmp_path):
    core = tmp_path / "gg_gene_evolution_core.sh"
    core.write_text("echo scientific-output; exit 7\n")
    directory = tmp_path / "observations"
    adapter = tmp_path / "adapter.py"
    adapter.write_text('''import os,sys
args=sys.argv[3:]
args=[sys.argv[1] if x=="/script/support/workflow_observation.py" else sys.argv[2] if x=="/workspace/output/observations" else x for x in args]
os.execvp(args[0],args)
''')
    code = f'''source {shlex.quote(str(SUPPORT / 'gg_util.sh'))}
container_adapter() {{ shift 2; python {shlex.quote(str(adapter))} {shlex.quote(str(OBSERVER))} {shlex.quote(str(directory))} "$@"; }}
singularity_command=(container_adapter exec)
GG_RESOURCE_METRICS=0
export mode_gene_evolution=orthogroup
GG_OBSERVABILITY=0
gg_run_container_shell_script image.sif {shlex.quote(str(core))}
'''
    assert subprocess.run(["bash", "-c", code], capture_output=True).returncode == 7
    assert not directory.exists()
    result = subprocess.run(["bash", "-c", code.replace("GG_OBSERVABILITY=0", "GG_OBSERVABILITY=1")],
                            text=True, capture_output=True)
    assert result.returncode == 7, result.stderr
    assert result.stdout == "scientific-output\n"
    record = query("status", "--directory", directory)["records"][0]
    assert record["runtime"]["effective_config"]["mode_gene_evolution"] == "orthogroup"


def test_archive_fasta_checks_use_the_runtime_sequence_validator(project):
    workspace, root, _, output, _, argv, plan = project
    output.write_text(">A\nMELK\n")
    argv += ["--output-fasta-type", "table=protein"]
    assert cli(PROVENANCE, "record", *argv).returncode == 0
    count = workspace / "genecount.tsv"
    count.write_text("Orthogroup\tsp1\tTotal\nOG0001\t1\t1\n")
    assert cli(SUPPORT / "gene_family_output_store.py", "archive-family", "--root", root,
               "--mode", "orthogroup", "--genecount", count, "--family-id", "OG0001").returncode == 0
    plan.write_text(json.dumps({"schema": "genegalleon-preflight-plan-v1", "contracts": [argv]}))
    assert query("preflight", "--plan", plan)["contracts"][0]["state"] == "verified_current"
    argv[-1] = "table=codon"
    plan.write_text(json.dumps({"schema": "genegalleon-preflight-plan-v1", "contracts": [argv]}))
    response = query("preflight", "--plan", plan)
    assert response["blocked"]
    assert not response["contracts"][0]["historical_generation_verified"]


def test_observation_fences_reject_concurrent_archive_changes_without_lock_writes(project):
    _, root, *_ = project
    code = f'''import sys
from pathlib import Path
sys.path.insert(0, {str(SUPPORT)!r})
from gene_family_output_store import GeneFamilyOutputStore, read_only_observation, ArchiveStoreError
root = Path({str(root)!r})
try:
    with read_only_observation():
        store = GeneFamilyOutputStore(root)
        store.artifacts("stat_branch")
        store.archive_root.mkdir()
        (store.archive_root / "index.epoch").write_text("changed")
except ArchiveStoreError:
    sys.exit(0)
sys.exit(1)
'''
    result = subprocess.run([sys.executable, "-c", code], text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    assert not (root / ".gg_store_locks").exists()


def test_directory_contract_verification_matches_execution(project):
    _, _, _, output, _, argv, plan = project
    output.unlink()
    output.mkdir()
    (output / "nested.tsv").write_text("nested output\n")
    assert cli(PROVENANCE, "record", *argv).returncode == 0
    assert verify(project)["completion_state"] == "verified_declared_steps"
    assert query("preflight", "--plan", plan)["contracts"][0]["state"] == "verified_current"


def test_success_skip_exit_code_remains_unchanged(tmp_path):
    directory = tmp_path / "observations"
    result = cli(OBSERVER, "--directory", directory, "--workflow", "gg_gene_evolution",
                 "--accept-exit-code", "8", "--", "bash", "-c", "exit 8")
    assert result.returncode == 8
    record = query("status", "--directory", directory)["records"][0]
    assert record["execution_accepted"]
    assert record["error_code"] is None
    assert record["completion_state"] == "unverified"


def test_unrelated_malformed_manifest_does_not_block_exact_family_verification(project):
    _, _, _, _, manifest, argv, _ = project
    assert cli(PROVENANCE, "record", *argv).returncode == 0
    manifest.with_name("OG9999.summary_statistics.json").write_text("malformed unrelated record")
    assert verify(project)["completion_state"] == "verified_declared_steps"


def test_running_family_cannot_inherit_old_completed_outputs(project):
    _, root, _, _, _, argv, _ = project
    assert cli(PROVENANCE, "record", *argv).returncode == 0
    code = f'''import sys
sys.path.insert(0, {str(SUPPORT)!r})
from gene_family_output_store import GeneFamilyOutputStore
GeneFamilyOutputStore({str(root)!r}).mark_family_state("OG0001", "running", "new-attempt")
'''
    assert subprocess.run([sys.executable, "-c", code]).returncode == 0
    response = verify(project)
    assert response["analytical_state"] == "running"
    assert response["completion_state"] == "unverified"


def test_status_reports_new_contract_receipt_without_a_run_record_change(project, tmp_path):
    _, _, _, _, _, argv, _ = project
    directory = tmp_path / "observations"
    assert cli(OBSERVER, "--directory", directory, "--workflow", "gg_gene_evolution", "--", "true").returncode == 0
    first = query("status", "--directory", directory)
    attempt = directory / first["records"][0]["attempt_id"]
    run_bytes = (attempt / "run.json").read_bytes()
    result = cli(PROVENANCE, "record", *argv, env={**os.environ, "GG_OBSERVATION_ATTEMPT_DIR": str(attempt)})
    assert result.returncode == 0
    assert (attempt / "run.json").read_bytes() == run_bytes
    response = query("status", "--directory", directory, "--since", first["next_cursor"])
    assert response["records"][0]["contract_evidence"][0]["step"] == "summary_statistics"


def test_archived_optional_only_legacy_artifact_is_reusable(project):
    workspace, root, _, output, manifest, argv, plan = project
    index = argv.index("--output")
    argv[index] = "--optional-output"
    count = workspace / "genecount.tsv"
    count.write_text("Orthogroup\tsp1\tTotal\nOG0001\t1\t1\n")
    assert cli(SUPPORT / "gene_family_output_store.py", "archive-family", "--root", root,
               "--mode", "orthogroup", "--genecount", count, "--family-id", "OG0001").returncode == 0
    assert not output.exists() and not manifest.exists()
    plan.write_text(json.dumps({"schema": "genegalleon-preflight-plan-v1", "contracts": [argv]}))
    before = snapshot(workspace)
    assert query("preflight", "--plan", plan)["contracts"][0]["state"] == "legacy_reusable"
    assert snapshot(workspace) == before


def test_audit_stale_policy_reuse_is_not_exact_attempt_verification(project, tmp_path):
    _, _, _, _, _, argv, _ = project
    assert cli(PROVENANCE, "record", *argv).returncode == 0
    changed = ["mode=b" if value == "mode=a" else value for value in argv]
    directory = tmp_path / "observations"
    result = cli(OBSERVER, "--directory", directory, "--workflow", "gg_gene_evolution",
                 "--accept-exit-code", "1", "--", sys.executable, PROVENANCE,
                 "needs-run", *changed, "--stale-policy", "reuse")
    assert result.returncode == 1
    attempt = next(directory.iterdir())
    assert verify(project, "--attempt", attempt)["completion_state"] == "unverified"


def test_audit_identical_other_project_is_not_bound_to_attempt(project, tmp_path):
    workspace, _, _, _, _, argv, _ = project
    directory = tmp_path / "observations"
    assert cli(OBSERVER, "--directory", directory, "--workflow", "gg_gene_evolution", "--",
               sys.executable, PROVENANCE, "record", *argv).returncode == 0
    other = tmp_path / "other-workspace"
    shutil.copytree(workspace, other)
    result = query("verify", "--root", other / "output/orthogroup", "--workspace-root", other,
                   "--family-id", "OG0001", "--require-step", "summary_statistics",
                   "--attempt", next(directory.iterdir()))
    assert result["completion_state"] == "unverified"


def test_audit_optional_adoption_is_not_historical_proof(project):
    _, _, _, output, manifest, argv, plan = project
    assert cli(PROVENANCE, "record", *argv).returncode == 0
    payload = json.loads(manifest.read_text())
    payload.pop("optional_outputs")
    manifest.write_text(json.dumps(payload))
    argv += ["--optional-output", "optional=" + str(output.with_name("optional.tsv"))]
    assert cli(PROVENANCE, "needs-run", *argv).returncode == 1
    plan.write_text(json.dumps({"schema": "genegalleon-preflight-plan-v1", "contracts": [argv]}))
    assert query("preflight", "--plan", plan)["contracts"][0]["historical_generation_verified"] is False
    assert verify(project)["completion_state"] == "unverified"


def test_audit_absent_optional_directory_appearing_is_a_change(project):
    _, _, _, output, _, argv, _ = project
    optional = output.with_name("optional-directory")
    assert cli(PROVENANCE, "record", *argv, "--optional-output", f"optional={optional}").returncode == 0
    optional.mkdir()
    (optional / "data").write_text("new output")
    response = verify(project)
    assert response["completion_state"] == "unverified"
    assert response["contracts"][0]["error_code"] == "unexpected_optional_output"


def test_audit_missing_run_file_does_not_advance_cursor(tmp_path):
    directory = tmp_path / "observations"
    assert cli(OBSERVER, "--directory", directory, "--workflow", "gg_test", "--", "true").returncode == 0
    first = query("status", "--directory", directory)
    next(directory.glob("*/run.json")).unlink()
    response = query("status", "--directory", directory, "--since", first["next_cursor"])
    assert response["complete"] is False
    assert response["next_cursor"] is None
    assert response["removed"] == []


def test_audit_help_abbreviation_cannot_pollute_json(tmp_path):
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({"schema": "genegalleon-preflight-plan-v1", "contracts": [["--he"]]}))
    result = cli(API, "preflight", "--plan", plan)
    assert result.returncode == 2
    assert json.loads(result.stdout)["error_code"] == "query_unavailable"


def test_audit_named_pipe_json_does_not_hang(tmp_path):
    fifo = tmp_path / "plan.json"
    os.mkfifo(fifo)
    result = subprocess.run([sys.executable, str(API), "preflight", "--plan", str(fifo)],
                            capture_output=True, text=True, timeout=3)
    assert result.returncode == 2
    assert json.loads(result.stdout)["error_code"] == "query_unavailable"


def test_audit_relative_contracts_can_be_inspected_from_another_cwd(project, tmp_path):
    workspace, _, _, _, _, argv, _ = project
    relative = [value.replace(str(workspace), "workspace") for value in argv]
    directory = tmp_path / "observations"
    result = subprocess.run([sys.executable, str(OBSERVER), "--directory", str(directory),
                             "--workflow", "gg_gene_evolution", "--", sys.executable,
                             str(PROVENANCE), "record", *relative], cwd=tmp_path, capture_output=True)
    assert result.returncode == 0, result.stderr
    assert query("preflight", "--attempt", next(directory.iterdir()))["contracts"][0]["state"] == "verified_current"


def test_audit_failed_later_validation_replaces_success_receipt(project, tmp_path):
    _, _, _, _, manifest, argv, _ = project
    attempt = tmp_path / "observations" / ("a" * 32)
    env = {**os.environ, "GG_OBSERVATION_ATTEMPT_DIR": str(attempt)}
    assert cli(PROVENANCE, "record", *argv, env=env).returncode == 0
    manifest.write_text("malformed")
    assert cli(PROVENANCE, "needs-run", *argv, env=env).returncode == 2
    receipt = json.loads(next(attempt.glob("contract-*.json")).read_text())
    assert receipt["exit_code"] == 2


def test_audit_named_pipe_manifest_does_not_hang(project):
    _, _, _, _, manifest, _, plan = project
    manifest.parent.mkdir(parents=True)
    os.mkfifo(manifest)
    result = subprocess.run([sys.executable, str(API), "preflight", "--plan", str(plan)],
                            capture_output=True, text=True, timeout=3)
    response = json.loads(result.stdout)
    assert response.get("error_code") == "query_unavailable" or response["contracts"][0]["state"] == "unavailable"


@pytest.mark.parametrize("extra", ['"schema":"wrong",', '"extra":NaN,'])
def test_audit_ambiguous_json_is_rejected(project, extra):
    *_, plan = project
    plan.write_text("{" + extra + plan.read_text()[1:])
    result = cli(API, "preflight", "--plan", plan)
    assert result.returncode == 2
    assert json.loads(result.stdout)["error_code"] == "query_unavailable"


def test_audit_boolean_manifest_schema_cannot_verify(project):
    *_, manifest, argv, plan = project
    assert cli(PROVENANCE, "record", *argv).returncode == 0
    payload = json.loads(manifest.read_text())
    payload["schema_version"] = True
    manifest.write_text(json.dumps(payload))
    assert verify(project)["completion_state"] == "unverified"
    result = cli(API, "preflight", "--plan", plan)
    assert result.returncode == 2
    assert json.loads(result.stdout)["error_code"] == "query_unavailable"


def paged_observations(tmp_path, count=3):
    directory = tmp_path / "observations"
    directory.mkdir()
    for index in range(count):
        identity = f"{index + 1:032x}"
        attempt = directory / identity
        attempt.mkdir()
        (attempt / "run.json").write_text(json.dumps({
            "schema": "genegalleon-observation-v1", "attempt_id": identity, "workflow": "gg_gene_evolution",
            "started_at_ns": 10, "finished_at_ns": 20, "accepted_exit_codes": [0],
            "execution_state": "exited", "execution_accepted": True, "exit_code": 0,
            "runtime": {"schema": "genegalleon-runtime-v1", "workflow": "gg_gene_evolution", "scheduler": {
                "GG_JOB_ID": "123", "GG_ARRAY_JOB_ID": "123", "GG_ARRAY_TASK_ID": str(index + 1)}}}))
    return directory


def test_paged_status_preserves_inventory_and_bounds_cursors(tmp_path):
    directory = paged_observations(tmp_path)
    before = snapshot(directory)
    first = query("status", "--directory", directory, "--page-size", 2)
    assert len(first["records"]) == 2
    assert first["snapshot_complete"] is False
    assert len(first["next_cursor"]) < 1024
    second = query("status", "--directory", directory, "--page-size", 2, "--page-cursor", first["next_cursor"])
    assert second["inventory_sha256"] == first["inventory_sha256"]
    assert second["page_after"] == first["page_last"]
    assert second["inventory_count"] == 3
    assert second["snapshot_complete"] is True and second["next_cursor"] is None
    assert len(second["records"]) == 1
    assert snapshot(directory) == before


def test_paged_status_deltas_include_fingerprints_for_unchanged_records(tmp_path):
    directory = paged_observations(tmp_path)
    first = query("status", "--directory", directory, "--page-size", 512)
    known = tmp_path / "known.json"
    known.write_text(json.dumps({"schema": "genegalleon-status-known-v1", "root_identity": first["root_identity"],
                                 "records": dict(first["fingerprints"])}))
    unchanged = query("status", "--directory", directory, "--page-size", 512, "--known-records-file", known)
    assert unchanged["records"] == [] and unchanged["fingerprints"] == first["fingerprints"]
    target = directory / first["fingerprints"][0][0] / "run.json"
    record = json.loads(target.read_text())
    record["runtime"]["genegalleon_version"] = "changed"
    target.write_text(json.dumps(record))
    changed = query("status", "--directory", directory, "--page-size", 512, "--known-records-file", known)
    assert len(changed["records"]) == 1
    assert changed["records"][0]["attempt_id"] == record["attempt_id"]


def test_paged_cursor_reports_invalidation_for_changed_membership_and_root(tmp_path):
    directory = paged_observations(tmp_path)
    first = query("status", "--directory", directory, "--page-size", 1)
    other = tmp_path / "other"
    shutil.copytree(directory, other)
    for root, cursor in ((other, first["next_cursor"]), (directory, "broken")):
        result = cli(API, "status", "--directory", root, "--page-size", 1, "--page-cursor", cursor)
        assert result.returncode == 2 and json.loads(result.stdout)["error_code"] == "cursor_invalid"
    shutil.rmtree(directory / first["fingerprints"][0][0])
    result = cli(API, "status", "--directory", directory, "--page-size", 1, "--page-cursor", first["next_cursor"])
    assert json.loads(result.stdout)["error_code"] == "cursor_invalid"


def test_paged_status_does_not_treat_evidence_errors_as_cursor_errors(tmp_path):
    directory = paged_observations(tmp_path)
    target = next(directory.glob("*/run.json"))
    target.write_text("invalid evidence")
    result = cli(API, "status", "--directory", directory, "--page-size", 512)
    value = json.loads(result.stdout)
    assert result.returncode == 2 and value["complete"] is False
    assert value["error_code"] == "query_unavailable"


def test_paged_status_has_a_payload_bound_independent_of_page_count(tmp_path):
    directory = paged_observations(tmp_path, 4)
    for path in directory.glob("*/run.json"):
        row = json.loads(path.read_text())
        row["runtime"]["fixture_setting"] = "x" * (4 * 1024 * 1024)
        path.write_text(json.dumps(row))
    first = query("status", "--directory", directory, "--page-size", 512)
    assert 0 < len(first["records"]) < 4
    assert len(json.dumps(first).encode()) < 16 * 1024 * 1024
    assert first["snapshot_complete"] is False
    second = query("status", "--directory", directory, "--page-size", 512, "--page-cursor", first["next_cursor"])
    assert len(first["records"]) + len(second["records"]) == 4
    assert second["snapshot_complete"] is True

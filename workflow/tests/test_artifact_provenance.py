import json
import os
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "support" / "artifact_provenance.py"
STORE_SCRIPT = Path(__file__).resolve().parents[1] / "support" / "gene_family_output_store.py"
GG_UTIL = Path(__file__).resolve().parents[1] / "support" / "gg_util.sh"


def bash_major_version():
    bash = shutil.which("bash")
    if bash is None:
        return 0
    completed = subprocess.run(
        [bash, "-c", 'printf "%s" "${BASH_VERSINFO[0]}"'],
        text=True,
        capture_output=True,
        check=False,
    )
    return int(completed.stdout or 0)


def run_cli(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *map(str, args)],
        text=True,
        capture_output=True,
        check=False,
    )


def contract_args(
    workspace,
    manifest,
    input_path,
    output_path,
    parameter="mode=a",
    optional_output=None,
):
    logical_root = workspace / "output" / "orthogroup"
    args = [
        "--manifest",
        manifest,
        "--step",
        "summary_statistics",
        "--family-id",
        "OG0001",
        "--logical-root",
        logical_root,
        "--workspace-root",
        workspace,
        "--input",
        f"rooted_tree={input_path}",
        "--output",
        f"stat_branch={output_path}",
        "--parameter",
        parameter,
    ]
    if optional_output is not None:
        args.extend(["--optional-output", f"optional={optional_output}"])
    return args


def test_manifest_detects_content_and_parameter_changes_but_ignores_versions(tmp_path):
    workspace = tmp_path / "workspace"
    logical_root = workspace / "output" / "orthogroup"
    input_path = logical_root / "rooted_tree" / "OG0001_root.nwk"
    output_path = logical_root / "stat_branch" / "OG0001_stat.branch.tsv"
    manifest = logical_root / "artifact_provenance" / "OG0001.summary_statistics.json"
    input_path.parent.mkdir(parents=True)
    output_path.parent.mkdir(parents=True)
    input_path.write_text("(A,B);\n", encoding="utf-8")
    output_path.write_text("branch_id\tnode_name\n0\tA\n", encoding="utf-8")
    args = contract_args(workspace, manifest, input_path, output_path)

    adopted = run_cli("needs-run", *args)
    assert adopted.returncode == 1, adopted.stderr
    adopted_payload = json.loads(manifest.read_text(encoding="utf-8"))
    assert adopted_payload["diagnostics"]["provenance_state"] == "adopted_legacy_output_without_rebuild"
    changed_after_adoption = contract_args(workspace, manifest, input_path, output_path, parameter="mode=b")
    assert run_cli("needs-run", *changed_after_adoption).returncode == 3
    recorded = run_cli("record", *args, "--diagnostic", "tool_version=1.0")
    assert recorded.returncode == 0, recorded.stderr
    assert run_cli("needs-run", *args).returncode == 1

    payload = json.loads(manifest.read_text(encoding="utf-8"))
    payload["diagnostics"]["tool_version"] = "9.9"
    payload["diagnostics"]["container_version"] = "new"
    manifest.write_text(json.dumps(payload), encoding="utf-8")
    assert run_cli("needs-run", *args).returncode == 1

    changed_parameter_args = contract_args(workspace, manifest, input_path, output_path, parameter="mode=b")
    assert run_cli("needs-run", *changed_parameter_args).returncode == 3

    input_path.write_text("((A,B),C);\n", encoding="utf-8")
    assert run_cli("needs-run", *args).returncode == 3

    input_path.write_text("(A,B);\n", encoding="utf-8")
    output_path.write_text("changed\n", encoding="utf-8")
    assert run_cli("needs-run", *args).returncode == 3


def test_stale_policy_stop_reuse_and_rebuild(tmp_path):
    workspace = tmp_path / "workspace"
    logical_root = workspace / "output" / "orthogroup"
    input_path = logical_root / "rooted_tree" / "OG0001_root.nwk"
    output_path = logical_root / "stat_branch" / "OG0001_stat.branch.tsv"
    manifest = logical_root / "artifact_provenance" / "OG0001.summary_statistics.json"
    input_path.parent.mkdir(parents=True)
    output_path.parent.mkdir(parents=True)
    input_path.write_text("(A,B);\n", encoding="utf-8")
    output_path.write_text("branch_id\tnode_name\n0\tA\n", encoding="utf-8")
    args = contract_args(workspace, manifest, input_path, output_path)
    assert run_cli("record", *args).returncode == 0
    original_manifest = manifest.read_bytes()
    changed_args = contract_args(workspace, manifest, input_path, output_path, parameter="mode=b")

    stopped = run_cli("needs-run", *changed_args)
    assert stopped.returncode == 3
    assert "No artifact files were modified" in stopped.stderr
    assert manifest.read_bytes() == original_manifest

    reused = run_cli("needs-run", *changed_args, "--stale-policy", "reuse")
    assert reused.returncode == 1
    assert "reuse stale output" in reused.stderr
    assert manifest.read_bytes() == original_manifest

    rebuilt = run_cli("needs-run", *changed_args, "--stale-policy", "rebuild")
    assert rebuilt.returncode == 0
    assert "regenerate without confirmation" in rebuilt.stderr
    assert manifest.read_bytes() == original_manifest


@pytest.mark.skipif(bash_major_version() < 4, reason="gg_util.sh requires Bash 4+")
@pytest.mark.parametrize(
    ("policy", "expected_status", "expected_needs_update", "expected_run"),
    [
        ("stop", 3, 0, 0),
        ("reuse", 0, 0, 0),
        ("rebuild", 0, 1, 1),
    ],
)
def test_shell_stage_policy_can_stop_reuse_or_enable_disabled_stage(
    tmp_path, policy, expected_status, expected_needs_update, expected_run
):
    workspace = tmp_path / "workspace"
    logical_root = workspace / "output" / "orthogroup"
    input_path = logical_root / "rooted_tree" / "OG0001_root.nwk"
    output_path = logical_root / "stat_branch" / "OG0001_stat.branch.tsv"
    manifest = logical_root / "artifact_provenance" / "OG0001.summary_statistics.json"
    input_path.parent.mkdir(parents=True)
    output_path.parent.mkdir(parents=True)
    input_path.write_text("(A,B);\n", encoding="utf-8")
    output_path.write_text("branch_id\n0\n", encoding="utf-8")
    original_args = contract_args(workspace, manifest, input_path, output_path, parameter="mode=a")
    assert run_cli("record", *original_args).returncode == 0

    script = f'''
gg_support_dir="{GG_UTIL.parent}"
gg_workspace_dir="{workspace}"
gg_workspace_output_dir="{workspace / 'output'}"
artifact_stale_policy="{policy}"
source "{GG_UTIL}"
run_demo=0
needs_update=0
args=(
  --manifest "{manifest}"
  --step summary_statistics
  --family-id OG0001
  --logical-root "{logical_root}"
  --workspace-root "{workspace}"
  --input "rooted_tree={input_path}"
  --output "stat_branch={output_path}"
  --parameter mode=b
)
if gg_artifact_prepare_stage needs_update run_demo "${{args[@]}}"; then
  stage_status=0
else
  stage_status=$?
fi
printf 'status=%s needs=%s run=%s\n' "${{stage_status}}" "${{needs_update}}" "${{run_demo}}"
'''
    completed = subprocess.run(
        [shutil.which("bash"), "-c", script],
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert (
        f"status={expected_status} needs={expected_needs_update} run={expected_run}"
        in completed.stdout
    )

def test_gene_family_store_input_is_content_based_and_storage_layout_independent(tmp_path):
    workspace = tmp_path / "workspace"
    logical_root = workspace / "output" / "orthogroup"
    stat_path = logical_root / "stat_branch" / "OG0001_stat.branch.tsv"
    database = workspace / "output" / "gene_summary" / "orthogroup.db"
    manifest = workspace / "output" / "gene_summary" / "database.json"
    stat_path.parent.mkdir(parents=True)
    database.parent.mkdir(parents=True)
    stat_path.write_text("branch_id\n1\n", encoding="utf-8")
    database.write_text("database\n", encoding="utf-8")
    args = [
        "--manifest",
        manifest,
        "--step",
        "gene_family_database",
        "--family-id",
        "orthogroup",
        "--logical-root",
        workspace / "output" / ".gg_global_artifacts",
        "--workspace-root",
        workspace,
        "--input-gene-family-store",
        f"gene_family_outputs={logical_root}",
        "--output",
        f"database={database}",
    ]
    assert run_cli("record", *args).returncode == 0
    assert run_cli("needs-run", *args).returncode == 1

    # Store bookkeeping and family provenance are not logical analysis inputs.
    (logical_root / ".gg_store").mkdir()
    (logical_root / ".gg_store" / "archive.lock").write_text("lock\n", encoding="utf-8")
    family_manifest = logical_root / "artifact_provenance" / "OG0001.test.json"
    family_manifest.parent.mkdir()
    family_manifest.write_text("{}\n", encoding="utf-8")
    assert run_cli("needs-run", *args).returncode == 1

    genecount = workspace / "output" / "orthofinder" / "Orthogroups.GeneCount.selected.tsv"
    genecount.parent.mkdir(parents=True)
    genecount.write_text("Orthogroup\tsp1\tTotal\nOG0001\t1\t1\n", encoding="utf-8")
    archived = subprocess.run(
        [
            sys.executable,
            str(STORE_SCRIPT),
            "archive-family",
            "--root",
            str(logical_root),
            "--mode",
            "orthogroup",
            "--genecount",
            str(genecount),
            "--family-id",
            "OG0001",
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    assert archived.returncode == 0, archived.stderr
    assert not stat_path.exists()
    assert run_cli("needs-run", *args).returncode == 1


def test_gene_family_store_input_detects_same_size_same_mtime_content_change(tmp_path):
    workspace = tmp_path / "workspace"
    logical_root = workspace / "output" / "orthogroup"
    stat_path = logical_root / "stat_branch" / "OG0001_stat.branch.tsv"
    database = workspace / "output" / "gene_summary" / "orthogroup.db"
    manifest = workspace / "output" / "gene_summary" / "database.json"
    stat_path.parent.mkdir(parents=True)
    database.parent.mkdir(parents=True)
    stat_path.write_text("AAAA\n", encoding="utf-8")
    original_stat = stat_path.stat()
    database.write_text("database\n", encoding="utf-8")
    args = [
        "--manifest",
        manifest,
        "--step",
        "gene_family_database",
        "--family-id",
        "orthogroup",
        "--logical-root",
        workspace / "output" / ".gg_global_artifacts",
        "--workspace-root",
        workspace,
        "--input-gene-family-store",
        f"gene_family_outputs={logical_root}",
        "--output",
        f"database={database}",
    ]
    assert run_cli("record", *args).returncode == 0
    stat_path.write_text("BBBB\n", encoding="utf-8")
    os.utime(stat_path, ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns))
    assert stat_path.stat().st_size == original_stat.st_size
    assert stat_path.stat().st_mtime_ns == original_stat.st_mtime_ns
    assert run_cli("needs-run", *args).returncode == 3


def test_gene_family_subdir_input_ignores_other_subdirs_and_raw_zip_layout(tmp_path):
    workspace = tmp_path / "workspace"
    logical_root = workspace / "output" / "orthogroup"
    rooted_tree = logical_root / "rooted_tree" / "OG0001_root.nwk"
    stat_path = logical_root / "stat_branch" / "OG0001_stat.branch.tsv"
    output_path = workspace / "output" / "genome_evolution" / "grampa.tsv"
    manifest = workspace / "output" / "artifact_provenance" / "grampa.json"
    rooted_tree.parent.mkdir(parents=True)
    stat_path.parent.mkdir(parents=True)
    output_path.parent.mkdir(parents=True)
    rooted_tree.write_text("(A,B);\n", encoding="utf-8")
    stat_path.write_text("branch_id\n1\n", encoding="utf-8")
    output_path.write_text("model\nH0\n", encoding="utf-8")
    args = [
        "--manifest",
        manifest,
        "--step",
        "orthogroup_grampa",
        "--family-id",
        "selected_orthogroups",
        "--logical-root",
        workspace / "output" / ".gg_global_artifacts",
        "--workspace-root",
        workspace,
        "--input-gene-family-subdir",
        f"rooted_trees={logical_root}::rooted_tree",
        "--output",
        f"summary={output_path}",
    ]
    assert run_cli("record", *args).returncode == 0

    stat_path.write_text("branch_id\n2\n", encoding="utf-8")
    assert run_cli("needs-run", *args).returncode == 1

    genecount = workspace / "output" / "orthofinder" / "Orthogroups.GeneCount.selected.tsv"
    genecount.parent.mkdir(parents=True)
    genecount.write_text("Orthogroup\tsp1\tTotal\nOG0001\t1\t1\n", encoding="utf-8")
    archived = subprocess.run(
        [
            sys.executable,
            str(STORE_SCRIPT),
            "archive-family",
            "--root",
            str(logical_root),
            "--mode",
            "orthogroup",
            "--genecount",
            str(genecount),
            "--family-id",
            "OG0001",
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    assert archived.returncode == 0, archived.stderr
    assert not rooted_tree.exists()
    assert run_cli("needs-run", *args).returncode == 1


def test_gene_family_artifact_input_is_precise_and_raw_zip_layout_independent(tmp_path):
    workspace = tmp_path / "workspace"
    logical_root = workspace / "output" / "orthogroup"
    selected = logical_root / "stat_branch" / "OG0001_stat.branch.tsv"
    unrelated = logical_root / "stat_branch" / "OG0002_stat.branch.tsv"
    output_path = workspace / "output" / "hgt" / "tree_plot" / "OG0001.pdf"
    manifest = workspace / "output" / "hgt" / "artifact_provenance" / "OG0001.json"
    selected.parent.mkdir(parents=True)
    output_path.parent.mkdir(parents=True)
    selected.write_text("branch_id\n1\n", encoding="utf-8")
    unrelated.write_text("branch_id\n2\n", encoding="utf-8")
    output_path.write_text("pdf\n", encoding="utf-8")
    args = [
        "--manifest",
        manifest,
        "--step",
        "hgt_tree_plot",
        "--family-id",
        "OG0001",
        "--logical-root",
        workspace / "output" / ".gg_global_artifacts",
        "--workspace-root",
        workspace,
        "--input-gene-family-artifact",
        f"stat_branch={logical_root}::stat_branch::{selected.name}",
        "--output",
        f"plot={output_path}",
    ]
    assert run_cli("record", *args).returncode == 0

    unrelated.write_text("branch_id\n3\n", encoding="utf-8")
    assert run_cli("needs-run", *args).returncode == 1

    genecount = workspace / "output" / "orthofinder" / "Orthogroups.GeneCount.selected.tsv"
    genecount.parent.mkdir(parents=True)
    genecount.write_text(
        "Orthogroup\tsp1\tTotal\nOG0001\t1\t1\nOG0002\t1\t1\n",
        encoding="utf-8",
    )
    archived = subprocess.run(
        [
            sys.executable,
            str(STORE_SCRIPT),
            "archive-family",
            "--root",
            str(logical_root),
            "--mode",
            "orthogroup",
            "--genecount",
            str(genecount),
            "--family-id",
            "OG0001",
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    assert archived.returncode == 0, archived.stderr
    assert not selected.exists()
    assert run_cli("needs-run", *args).returncode == 1

    materialized = subprocess.run(
        [
            sys.executable,
            str(STORE_SCRIPT),
            "materialize-family",
            "--root",
            str(logical_root),
            "--mode",
            "orthogroup",
            "--family-id",
            "OG0001",
            "--destination-root",
            str(workspace / "materialized"),
            "--subdirs",
            "stat_branch",
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    assert materialized.returncode == 0, materialized.stderr
    selected_copy = workspace / "materialized" / "stat_branch" / selected.name
    selected_copy.write_text("branch_id\n9\n", encoding="utf-8")
    # Changing a temporary materialization cannot affect the logical input.
    assert run_cli("needs-run", *args).returncode == 1


def test_logical_directory_input_is_raw_zip_layout_independent(tmp_path):
    workspace = tmp_path / "workspace"
    managed = workspace / "output" / "species_tree" / "single_copy_alignment"
    member = managed / "OG0001.fa"
    output_path = workspace / "output" / "species_tree" / "concat.fa"
    manifest = workspace / "output" / "artifact_provenance" / "concat.json"
    managed.mkdir(parents=True)
    member.write_text(">A\nAA\n", encoding="utf-8")
    output_path.write_text(">A\nAA\n", encoding="utf-8")
    args = [
        "--manifest",
        manifest,
        "--step",
        "concat",
        "--family-id",
        "all_buscos",
        "--logical-root",
        workspace / "output" / ".gg_global_artifacts",
        "--workspace-root",
        workspace,
        "--input-logical-directory",
        f"alignments={managed}",
        "--output",
        f"concat={output_path}",
    ]
    assert run_cli("record", *args).returncode == 0
    archive_path = managed.with_name(f"{managed.name}.zip")
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.write(member, f"{managed.name}/{member.name}")
    member.unlink()
    managed.rmdir()
    assert run_cli("needs-run", *args).returncode == 1

    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr(f"{managed.name}/{member.name}", ">A\nBB\n")
    assert run_cli("needs-run", *args).returncode == 3


def test_optional_output_records_present_and_absent_states(tmp_path):
    workspace = tmp_path / "workspace"
    logical_root = workspace / "output" / "orthogroup"
    input_path = logical_root / "iqtree_anc" / "OG0001.zip"
    output_path = logical_root / "csubst_b" / "OG0001.tsv"
    optional_path = logical_root / "csubst_cb_2" / "OG0001.tsv"
    manifest = logical_root / "artifact_provenance" / "OG0001.csubst.json"
    input_path.parent.mkdir(parents=True)
    output_path.parent.mkdir(parents=True)
    input_path.write_text("input\n", encoding="utf-8")
    output_path.write_text("branch\n", encoding="utf-8")
    args = contract_args(
        workspace,
        manifest,
        input_path,
        output_path,
        optional_output=optional_path,
    )

    assert run_cli("record", *args).returncode == 0
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    assert payload["optional_outputs"][0]["state"] == "absent"
    assert run_cli("needs-run", *args).returncode == 1

    optional_path.parent.mkdir(parents=True)
    optional_path.write_text("combination\n", encoding="utf-8")
    assert run_cli("needs-run", *args).returncode == 3
    assert run_cli("needs-run", *args, "--stale-policy", "reuse").returncode == 1
    assert run_cli("needs-run", *args, "--stale-policy", "rebuild").returncode == 0


def test_optional_only_absence_is_a_completed_recorded_result(tmp_path):
    workspace = tmp_path / "workspace"
    input_path = workspace / "input.tsv"
    optional_path = workspace / "significant.tsv"
    manifest = workspace / "artifact_provenance" / "optional.json"
    workspace.mkdir()
    input_path.write_text("value\n0\n", encoding="utf-8")
    args = [
        "--manifest",
        manifest,
        "--step",
        "significance_filter",
        "--family-id",
        "all",
        "--logical-root",
        workspace / "output" / ".gg_global_artifacts",
        "--workspace-root",
        workspace,
        "--input",
        f"table={input_path}",
        "--optional-output",
        f"significant={optional_path}",
    ]
    assert run_cli("needs-run", *args).returncode == 0
    assert run_cli("record", *args).returncode == 0
    assert run_cli("needs-run", *args).returncode == 1
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    assert payload["optional_outputs"] == [
        {
            "label": "significant",
            "path": "significant.tsv",
            "scope": "workspace",
            "state": "absent",
        }
    ]


def test_existing_manifest_adopts_new_optional_output_contract_once(tmp_path):
    workspace = tmp_path / "workspace"
    logical_root = workspace / "output" / "orthogroup"
    input_path = logical_root / "iqtree_anc" / "OG0001.zip"
    output_path = logical_root / "csubst_b" / "OG0001.tsv"
    optional_path = logical_root / "csubst_cb_2" / "OG0001.tsv"
    manifest = logical_root / "artifact_provenance" / "OG0001.csubst.json"
    input_path.parent.mkdir(parents=True)
    output_path.parent.mkdir(parents=True)
    optional_path.parent.mkdir(parents=True)
    input_path.write_text("input\n", encoding="utf-8")
    output_path.write_text("branch\n", encoding="utf-8")
    optional_path.write_text("combination\n", encoding="utf-8")
    legacy_args = contract_args(workspace, manifest, input_path, output_path)
    assert run_cli("record", *legacy_args).returncode == 0
    legacy_payload = json.loads(manifest.read_text(encoding="utf-8"))
    legacy_payload.pop("optional_outputs")
    manifest.write_text(json.dumps(legacy_payload), encoding="utf-8")

    args = contract_args(
        workspace,
        manifest,
        input_path,
        output_path,
        optional_output=optional_path,
    )
    adopted = run_cli("needs-run", *args)
    assert adopted.returncode == 1, adopted.stderr
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    assert payload["optional_outputs"][0]["state"] == "present"
    assert payload["diagnostics"]["optional_output_provenance_state"] == "adopted_legacy_optional_outputs"


def test_missing_manifest_still_requests_run_when_declared_output_is_missing(tmp_path):
    workspace = tmp_path / "workspace"
    logical_root = workspace / "output" / "orthogroup"
    input_path = logical_root / "rooted_tree" / "OG0001_root.nwk"
    output_path = logical_root / "stat_branch" / "OG0001_stat.branch.tsv"
    manifest = logical_root / "artifact_provenance" / "OG0001.summary_statistics.json"
    input_path.parent.mkdir(parents=True)
    input_path.write_text("(A,B);\n", encoding="utf-8")

    requested = run_cli("needs-run", *contract_args(workspace, manifest, input_path, output_path))
    assert requested.returncode == 0
    assert not manifest.exists()


def test_audit_reports_legacy_untracked_output_without_failing(tmp_path):
    workspace = tmp_path / "workspace"
    logical_root = workspace / "output" / "orthogroup"
    input_path = logical_root / "rooted_tree" / "OG0001_root.nwk"
    output_path = logical_root / "stat_branch" / "OG0001_stat.branch.tsv"
    manifest = logical_root / "artifact_provenance" / "OG0001.summary_statistics.json"
    input_path.parent.mkdir(parents=True)
    output_path.parent.mkdir(parents=True)
    input_path.write_text("(A,B);\n", encoding="utf-8")
    output_path.write_text("branch_id\tnode_name\n0\tA\n", encoding="utf-8")
    args = contract_args(workspace, manifest, input_path, output_path)
    assert run_cli("record", *args, "--diagnostic", "tool_version=1").returncode == 0

    unverified = logical_root / "csubst_scan" / "OG0002_csubst_scan.tsv"
    unverified.parent.mkdir(parents=True)
    unverified.write_text("site\n1\n", encoding="utf-8")
    report = workspace / "output" / "audit.tsv"
    audited = run_cli(
        "audit",
        "--logical-root",
        logical_root,
        "--workspace-root",
        workspace,
        "--output-tsv",
        report,
        "--mode",
        "orthogroup",
    )
    assert audited.returncode == 0, audited.stderr
    frame = pd.read_csv(report, sep="\t", keep_default_na=False)
    legacy = frame.loc[frame["status"] == "legacy_untracked", :]
    assert legacy[["family_id", "step"]].to_dict(orient="records") == [{"family_id": "OG0002", "step": "csubst_scan"}]


def test_audit_infers_raw_query2family_ids_from_query_directory(tmp_path):
    workspace = tmp_path / "workspace"
    logical_root = workspace / "output" / "query2family"
    query_dir = workspace / "input" / "query_gene"
    query_dir.mkdir(parents=True)
    (query_dir / "gene.with.dot").write_text("query\n", encoding="utf-8")
    output = logical_root / "csubst_scan" / "gene.with.dot_csubst_scan.tsv"
    output.parent.mkdir(parents=True)
    output.write_text("site\n1\n", encoding="utf-8")
    report = workspace / "output" / "audit.tsv"

    audited = run_cli(
        "audit",
        "--logical-root",
        logical_root,
        "--workspace-root",
        workspace,
        "--output-tsv",
        report,
        "--mode",
        "query2family",
        "--query-dir",
        query_dir,
    )
    assert audited.returncode == 0, audited.stderr
    frame = pd.read_csv(report, sep="\t", keep_default_na=False)
    row = frame.loc[frame["status"] == "legacy_untracked", :].iloc[0]
    assert row["family_id"] == "gene.with.dot"
    assert row["step"] == "csubst_scan"


def test_audit_detects_changed_input_from_recorded_contract(tmp_path):
    workspace = tmp_path / "workspace"
    logical_root = workspace / "output" / "orthogroup"
    input_path = logical_root / "rooted_tree" / "OG0001_root.nwk"
    output_path = logical_root / "stat_branch" / "OG0001_stat.branch.tsv"
    manifest = logical_root / "artifact_provenance" / "OG0001.summary_statistics.json"
    input_path.parent.mkdir(parents=True)
    output_path.parent.mkdir(parents=True)
    input_path.write_text("(A,B);\n", encoding="utf-8")
    output_path.write_text("branch_id\tnode_name\n0\tA\n", encoding="utf-8")
    args = contract_args(workspace, manifest, input_path, output_path)
    assert run_cli("record", *args).returncode == 0
    input_path.write_text("((A,B),C);\n", encoding="utf-8")

    report = workspace / "output" / "audit.tsv"
    audited = run_cli(
        "audit",
        "--logical-root",
        logical_root,
        "--workspace-root",
        workspace,
        "--output-tsv",
        report,
        "--mode",
        "orthogroup",
    )
    assert audited.returncode == 1
    frame = pd.read_csv(report, sep="\t", keep_default_na=False)
    row = frame.loc[frame["step"] == "summary_statistics", :].iloc[0]
    assert row["status"] == "changed_input"
    assert row["reason"] == "rooted_tree"

    reused = run_cli(
        "audit",
        "--logical-root",
        logical_root,
        "--workspace-root",
        workspace,
        "--output-tsv",
        report,
        "--mode",
        "orthogroup",
        "--stale-policy",
        "reuse",
    )
    assert reused.returncode == 0, reused.stderr


def test_audit_reads_manifest_inputs_and_outputs_from_zip_store(tmp_path):
    workspace = tmp_path / "workspace"
    logical_root = workspace / "output" / "orthogroup"
    input_path = logical_root / "rooted_tree" / "OG0001_root.nwk"
    output_path = logical_root / "stat_branch" / "OG0001_stat.branch.tsv"
    manifest = logical_root / "artifact_provenance" / "OG0001.summary_statistics.json"
    input_path.parent.mkdir(parents=True)
    output_path.parent.mkdir(parents=True)
    input_path.write_text("(A,B);\n", encoding="utf-8")
    output_path.write_text("branch_id\tnode_name\n0\tA\n", encoding="utf-8")
    args = contract_args(workspace, manifest, input_path, output_path)
    assert run_cli("record", *args).returncode == 0

    genecount = workspace / "output" / "orthofinder" / "Orthogroups.GeneCount.selected.tsv"
    genecount.parent.mkdir(parents=True)
    genecount.write_text("Orthogroup\tsp1\tTotal\nOG0001\t1\t1\n", encoding="utf-8")
    archived = subprocess.run(
        [
            sys.executable,
            str(STORE_SCRIPT),
            "archive-family",
            "--root",
            str(logical_root),
            "--mode",
            "orthogroup",
            "--genecount",
            str(genecount),
            "--family-id",
            "OG0001",
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    assert archived.returncode == 0, archived.stderr
    assert not manifest.exists()
    assert not input_path.exists()
    assert not output_path.exists()

    report = workspace / "output" / "audit.tsv"
    audited = run_cli(
        "audit",
        "--logical-root",
        logical_root,
        "--workspace-root",
        workspace,
        "--output-tsv",
        report,
        "--mode",
        "orthogroup",
    )
    assert audited.returncode == 0, audited.stderr
    frame = pd.read_csv(report, sep="\t", keep_default_na=False)
    assert frame[["family_id", "step", "status"]].to_dict(orient="records") == [
        {
            "family_id": "OG0001",
            "step": "summary_statistics",
            "status": "current",
        }
    ]


def test_audit_checks_all_csubst_and_stat_branch_clades(tmp_path):
    pytest.importorskip("csubst")
    workspace = tmp_path / "workspace"
    logical_root = workspace / "output" / "orthogroup"
    iqtree_path = logical_root / "iqtree_anc" / "OG0001_iqtree.anc.zip"
    stat_path = logical_root / "stat_branch" / "OG0001_stat.branch.tsv"
    iqtree_path.parent.mkdir(parents=True)
    stat_path.parent.mkdir(parents=True)
    with zipfile.ZipFile(iqtree_path, "w") as archive:
        archive.writestr(
            "OG0001.iqtree.anc/csubst.nwk",
            "(sp1:0.1,(sp2:0.1,sp3:0.1)RootedBC:0.1)RootedRoot;\n",
        )
        archive.writestr(
            "OG0001.iqtree.anc/csubst.treefile",
            "(sp1:0.1,(sp2:0.1,sp3:0.1)IqtreeBC:0.1)IqtreeRoot;\n",
        )
    stat_path.write_text(
        "branch_id\tnode_name\tgene_labels\n"
        "0\tsp1\tsp1\n"
        "1\tsp2\tsp2\n"
        "2\tsp3\tsp3\n"
        "3\tN1\tsp2; sp3\n"
        "4\tRoot\tsp1; sp2; sp3\n",
        encoding="utf-8",
    )
    report = workspace / "output" / "audit.tsv"
    command = [
        "audit",
        "--logical-root",
        logical_root,
        "--workspace-root",
        workspace,
        "--output-tsv",
        report,
        "--mode",
        "orthogroup",
        "--require-step-for-subdir",
        "unused=unused",
        "--check-csubst-branches",
    ]

    consistent = run_cli(*command)
    assert consistent.returncode == 0, consistent.stderr
    frame = pd.read_csv(report, sep="\t", keep_default_na=False)
    assert frame.loc[0, "status"] == "current"

    stat_path.write_text(
        "branch_id\tnode_name\tgene_labels\n"
        "0\tsp1\tsp2\n"
        "1\tsp2\tsp1\n"
        "2\tsp3\tsp3\n"
        "3\tN1\tsp2; sp3\n"
        "4\tRoot\tsp1; sp2; sp3\n",
        encoding="utf-8",
    )
    mismatched = run_cli(*command)
    assert mismatched.returncode == 1
    frame = pd.read_csv(report, sep="\t", keep_default_na=False)
    assert frame.loc[0, "status"] == "semantic_mismatch"
    assert "inconsistent rooted trees" in frame.loc[0, "reason"]

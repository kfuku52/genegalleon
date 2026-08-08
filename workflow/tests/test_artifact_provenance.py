import json
import subprocess
import sys
import zipfile
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "support" / "artifact_provenance.py"
STORE_SCRIPT = Path(__file__).resolve().parents[1] / "support" / "gene_family_output_store.py"


def run_cli(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *map(str, args)],
        text=True,
        capture_output=True,
        check=False,
    )


def contract_args(workspace, manifest, input_path, output_path, parameter="mode=a"):
    logical_root = workspace / "output" / "orthogroup"
    return [
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
    assert run_cli("needs-run", *changed_after_adoption).returncode == 0
    recorded = run_cli("record", *args, "--diagnostic", "tool_version=1.0")
    assert recorded.returncode == 0, recorded.stderr
    assert run_cli("needs-run", *args).returncode == 1

    payload = json.loads(manifest.read_text(encoding="utf-8"))
    payload["diagnostics"]["tool_version"] = "9.9"
    payload["diagnostics"]["container_version"] = "new"
    manifest.write_text(json.dumps(payload), encoding="utf-8")
    assert run_cli("needs-run", *args).returncode == 1

    changed_parameter_args = contract_args(workspace, manifest, input_path, output_path, parameter="mode=b")
    assert run_cli("needs-run", *changed_parameter_args).returncode == 0

    input_path.write_text("((A,B),C);\n", encoding="utf-8")
    assert run_cli("needs-run", *args).returncode == 0

    input_path.write_text("(A,B);\n", encoding="utf-8")
    output_path.write_text("changed\n", encoding="utf-8")
    assert run_cli("needs-run", *args).returncode == 0


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

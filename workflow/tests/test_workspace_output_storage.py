import json
from pathlib import Path

import workflow.support.gene_family_output_store as FAMILY
import workflow.support.workspace_output_storage as WORKSPACE


def _workspace(tmp_path: Path, *, include_query: bool = True) -> Path:
    workspace = tmp_path / "workspace"
    output = workspace / "output"
    genecount = (
        output
        / "orthofinder"
        / "Orthogroups_filtered"
        / "Orthogroups.GeneCount.selected.tsv"
    )
    genecount.parent.mkdir(parents=True)
    genecount.write_text(
        "Orthogroup\tTotal\nOG0000001\t1\n",
        encoding="utf-8",
    )
    orthogroup = output / "orthogroup" / "mafft"
    orthogroup.mkdir(parents=True)
    (orthogroup / "OG0000001_cds.aln.fasta").write_text(
        ">a\nATG\n",
        encoding="utf-8",
    )
    if include_query:
        query_dir = workspace / "input" / "query_gene"
        query_dir.mkdir(parents=True)
        (query_dir / "queryA").write_text("queryA\n", encoding="utf-8")
        query_output = output / "query2family" / "mafft"
        query_output.mkdir(parents=True)
        (query_output / "queryA_cds.aln.fasta").write_text(
            ">a\nATG\n",
            encoding="utf-8",
        )
    species = output / "species_tree" / "single_copy_iqtree_dna"
    species.mkdir(parents=True)
    (species / "BUSCO1.dna.nwk").write_text("(A,B);\n", encoding="utf-8")
    return workspace


def _parse(*values: str):
    return WORKSPACE.build_parser().parse_args(list(values))


def test_workspace_dry_run_reports_all_targets_without_conversion(tmp_path: Path):
    workspace = _workspace(tmp_path)
    args = _parse(
        "convert",
        "--workspace",
        str(workspace),
        "--to",
        "zip",
        "--dry-run",
        "--report-prefix",
        "dry-run",
    )

    exit_code, report, json_path, tsv_path = WORKSPACE.run(args)

    assert exit_code == 0
    assert report["status"] == "complete"
    assert {target["name"] for target in report["targets"]} == {
        "query2family",
        "orthogroup",
        "species_tree",
    }
    assert all(target["status"] == "dry-run" for target in report["targets"])
    assert not list((workspace / "output").glob("**/*.zip"))
    assert json.loads(json_path.read_text(encoding="utf-8"))["dry_run"] is True
    assert tsv_path.read_text(encoding="utf-8").splitlines()[0].startswith("target\t")
    assert tsv_path.read_text(encoding="utf-8").splitlines()[-1].startswith("TOTAL\t")


def test_project_root_discovers_nested_workspace(tmp_path: Path):
    workspace = _workspace(tmp_path)
    resolved_workspace, output_root, input_root = WORKSPACE.resolve_workspace_layout(
        tmp_path
    )

    assert resolved_workspace == workspace.resolve()
    assert output_root == (workspace / "output").resolve()
    assert input_root == (workspace / "input").resolve()


def test_workspace_convert_zip_and_raw_roundtrip(tmp_path: Path):
    workspace = _workspace(tmp_path)
    genecount = (
        workspace
        / "output"
        / "orthofinder"
        / "Orthogroups_filtered"
        / "Orthogroups.GeneCount.selected.tsv"
    )
    original_genecount = genecount.read_bytes()
    to_zip = _parse(
        "convert",
        "--workspace",
        str(workspace),
        "--to",
        "zip",
        "--verification",
        "quick",
        "--report-prefix",
        "to-zip",
    )

    exit_code, report, _, _ = WORKSPACE.run(to_zip)

    assert exit_code == 0
    assert report["status"] == "complete"
    assert report["totals"]["before_managed_logical_files"] == 3
    assert report["totals"]["after_managed_physical_files"] > 0
    assert (workspace / "output" / "orthogroup" / "mafft.zip").is_file()
    assert (workspace / "output" / "query2family" / "mafft.zip").is_file()
    assert (
        workspace / "output" / "species_tree" / "single_copy_iqtree_dna.zip"
    ).is_file()
    assert not (workspace / "output" / "orthogroup" / "mafft").exists()
    assert genecount.read_bytes() == original_genecount
    assert str((workspace / "output" / "orthofinder").resolve()) in report["excluded_paths"]
    orthogroup_report = next(
        target for target in report["targets"] if target["name"] == "orthogroup"
    )
    assert orthogroup_report["preflight"]["managed_logical_files"] == 1
    assert orthogroup_report["postflight"]["physical_managed_files"] > 0

    to_raw = _parse(
        "convert",
        "--workspace",
        str(workspace),
        "--to",
        "raw",
        "--pure-raw",
        "--report-prefix",
        "to-raw",
    )
    exit_code, report, _, _ = WORKSPACE.run(to_raw)

    assert exit_code == 0
    assert report["status"] == "complete"
    assert (workspace / "output" / "orthogroup" / "mafft").is_dir()
    assert not (workspace / "output" / "orthogroup" / "mafft.zip").exists()
    assert (
        workspace / "output" / "species_tree" / "single_copy_iqtree_dna"
    ).is_dir()


def test_workspace_preflight_failure_prevents_all_conversions(tmp_path: Path):
    workspace = _workspace(tmp_path, include_query=False)
    query_output = workspace / "output" / "query2family" / "mafft"
    query_output.mkdir(parents=True)
    (query_output / "queryA_cds.aln.fasta").write_text(
        ">a\nATG\n",
        encoding="utf-8",
    )
    args = _parse(
        "convert",
        "--workspace",
        str(workspace),
        "--to",
        "zip",
        "--report-prefix",
        "failed-preflight",
    )

    exit_code, report, _, _ = WORKSPACE.run(args)

    assert exit_code == 1
    assert report["status"] == "failed"
    query = next(target for target in report["targets"] if target["name"] == "query2family")
    orthogroup = next(target for target in report["targets"] if target["name"] == "orthogroup")
    assert query["status"] == "failed"
    assert orthogroup["status"] == "not-run"
    assert not (workspace / "output" / "orthogroup" / "mafft.zip").exists()


def test_workspace_large_subdirectory_can_retain_meaningful_parts(tmp_path: Path):
    workspace = _workspace(tmp_path, include_query=False)
    args = _parse(
        "convert",
        "--workspace",
        str(workspace),
        "--target",
        "orthogroup",
        "--to",
        "zip",
        "--large-zip-warning-bytes",
        "1",
        "--max-final-zip-bytes",
        "1",
        "--report-prefix",
        "split",
    )

    exit_code, report, _, _ = WORKSPACE.run(args)

    assert exit_code == 0
    root = workspace / "output" / "orthogroup"
    assert not (root / "mafft.zip").exists()
    parts = list((root / "archives" / "mafft").glob("mafft.part-*.zip"))
    assert len(parts) == 1
    target = report["targets"][0]
    assert target["conversion"]["large_zip_subdirs"] == "mafft"
    assert target["conversion"]["split_zip_subdirs"] == "mafft"
    FAMILY.GeneFamilyOutputStore(root).verify(deep=True)


def test_workspace_raw_dry_run_is_quota_aware(tmp_path: Path):
    workspace = _workspace(tmp_path)
    to_zip = _parse(
        "convert",
        "--workspace",
        str(workspace),
        "--to",
        "zip",
        "--report-prefix",
        "quota-to-zip",
    )
    assert WORKSPACE.run(to_zip)[0] == 0

    dry_raw = _parse(
        "convert",
        "--workspace",
        str(workspace),
        "--to",
        "raw",
        "--dry-run",
        "--available-bytes",
        "0",
        "--report-prefix",
        "quota-raw-dry-run",
    )
    exit_code, report, _, _ = WORKSPACE.run(dry_raw)

    assert exit_code == 1
    assert report["status"] == "failed"
    assert all(target["status"] == "failed" for target in report["targets"])
    assert (
        workspace / "output" / "species_tree" / "single_copy_iqtree_dna.zip"
    ).is_file()
    assert (workspace / "output" / "orthogroup" / "mafft.zip").is_file()


def test_workspace_quota_preflight_accounts_for_sequential_growth(tmp_path: Path):
    workspace = _workspace(tmp_path, include_query=False)
    second_species = (
        workspace / "output" / "species_tree" / "single_copy_iqtree_pep"
    )
    second_species.mkdir(parents=True)
    (second_species / "BUSCO1.pep.nwk").write_text("(A,B);\n", encoding="utf-8")
    to_zip = _parse(
        "convert",
        "--workspace",
        str(workspace),
        "--target",
        "species_tree",
        "--to",
        "zip",
        "--report-prefix",
        "sequential-to-zip",
    )
    assert WORKSPACE.run(to_zip)[0] == 0

    unrestricted = _parse(
        "convert",
        "--workspace",
        str(workspace),
        "--target",
        "species_tree",
        "--to",
        "raw",
        "--dry-run",
        "--report-prefix",
        "sequential-unrestricted",
    )
    _, unrestricted_report, _, _ = WORKSPACE.run(unrestricted)
    operations = [
        row
        for row in unrestricted_report["quota_plan"]["operations"]
        if row["temporary_bytes"] > 0
    ]
    one_operation_headroom = max(row["temporary_bytes"] for row in operations)
    assert unrestricted_report["quota_plan"]["required_from_initial_bytes"] > (
        one_operation_headroom
    )

    constrained = _parse(
        "convert",
        "--workspace",
        str(workspace),
        "--target",
        "species_tree",
        "--to",
        "raw",
        "--dry-run",
        "--available-bytes",
        str(one_operation_headroom),
        "--report-prefix",
        "sequential-constrained",
    )
    exit_code, report, _, _ = WORKSPACE.run(constrained)

    assert exit_code == 1
    assert report["quota_plan"]["sufficient"] is False
    assert "insufficient-workspace-temporary-space" in (
        report["targets"][0]["preflight"]["issues"]
    )

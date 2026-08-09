import argparse
import io
import json
import os
import subprocess
import sys
import zipfile
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas as pd
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "csubst_scan_candidate_sites.py"


def load_module():
    spec = spec_from_file_location("csubst_scan_candidate_sites", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def candidate_rows():
    return pd.DataFrame(
        [
            {
                "orthogroup": "OG0002",
                "trait": "aquatic",
                "state_change": "A>V",
                "codon_site_alignment": 9,
                "support_unit_count": 5,
                "support_unit_ids": "1,2,3,4,5",
                "support_branch_ids": "7, 3,7",
                "p_rate_enrichment": 0.001,
                "q_rate_enrichment_global": 0.02,
                "besthit_0.05": "protein B",
            },
            {
                "orthogroup": "OG0001",
                "trait": "aquatic",
                "state_change": "G>D",
                "codon_site_alignment": 4,
                "support_unit_count": 6,
                "support_unit_ids": "1,2,3,4,5,6",
                "support_branch_ids": "8,9",
                "p_rate_enrichment": 0.002,
                "q_rate_enrichment_global": 0.01,
                "besthit_0.05": "protein A",
            },
            {
                "orthogroup": "OG0003",
                "trait": "aquatic",
                "state_change": "L>F",
                "codon_site_alignment": 12,
                "support_unit_count": 7,
                "support_unit_ids": "1,2,3,4,5,6,7",
                "support_branch_ids": "10,11",
                "p_rate_enrichment": 0.02,
                "q_rate_enrichment_global": 0.2,
                "besthit_0.05": "protein C",
            },
        ]
    )


def write_summary(path, frame=None):
    if frame is None:
        frame = candidate_rows()
    frame.to_csv(path, sep="\t", index=False)


def test_discover_summary_tables_returns_contiguous_thresholds_in_descending_order(tmp_path):
    mod = load_module()
    prefix = tmp_path / "orthogroup_csubst_aa_change"
    for threshold in (5, 6, 7):
        write_summary(tmp_path / f"{prefix.name}_min_support_{threshold}_summary.tsv")

    discovered = mod.discover_summary_tables(prefix, 5)

    assert list(discovered) == [7, 6, 5]


def test_discover_summary_tables_rejects_a_gap(tmp_path):
    mod = load_module()
    prefix = tmp_path / "orthogroup_csubst_aa_change"
    for threshold in (5, 7):
        write_summary(tmp_path / f"{prefix.name}_min_support_{threshold}_summary.tsv")

    with pytest.raises(FileNotFoundError, match="Missing threshold.*6"):
        mod.discover_summary_tables(prefix, 5)


def test_discover_summary_tables_returns_empty_when_observed_max_is_below_minimum(tmp_path):
    mod = load_module()
    prefix = tmp_path / "orthogroup_csubst_aa_change"
    for threshold in (2, 3, 4):
        write_summary(tmp_path / f"{prefix.name}_min_support_{threshold}_summary.tsv")

    assert mod.discover_summary_tables(prefix, 5) == {}


def test_load_threshold_candidates_filters_q_and_canonicalizes_branches(tmp_path):
    mod = load_module()
    summary = tmp_path / "summary.tsv"
    write_summary(summary)

    selected = mod.load_threshold_candidates(
        summary_path=summary,
        minimum_support=5,
        q_column="q_rate_enrichment_global",
        q_threshold=0.05,
        max_candidates=0,
        csubst_nonsyn_recode="no",
        pdb="none",
    )

    assert selected["orthogroup"].tolist() == ["OG0001", "OG0002"]
    assert selected["_canonical_support_branch_ids"].tolist() == ["8,9", "3,7"]
    assert selected["_candidate_rank"].tolist() == [1, 2]
    assert selected["_candidate_id"].str.contains(r"_[0-9a-f]{16}$").all()


def test_candidate_analysis_identity_is_stable_across_recalculated_q_values(tmp_path):
    mod = load_module()
    first_path = tmp_path / "first.tsv"
    second_path = tmp_path / "second.tsv"
    first = candidate_rows().iloc[[0]].copy()
    second = first.copy()
    second["q_rate_enrichment_global"] = 0.03
    write_summary(first_path, first)
    write_summary(second_path, second)

    loaded = [
        mod.load_threshold_candidates(
            summary_path=path,
            minimum_support=5,
            q_column="q_rate_enrichment_global",
            q_threshold=0.05,
            max_candidates=0,
            csubst_nonsyn_recode="no",
            pdb="none",
        )
        for path in (first_path, second_path)
    ]

    assert loaded[0].loc[0, "_analysis_key"] == loaded[1].loc[0, "_analysis_key"]
    assert loaded[0].loc[0, "_candidate_id"] == loaded[1].loc[0, "_candidate_id"]


def test_candidate_analysis_identity_ignores_tool_versions_but_tracks_parameters(
    monkeypatch, tmp_path
):
    mod = load_module()
    summary = tmp_path / "summary.tsv"
    write_summary(summary, candidate_rows().iloc[[0]].copy())

    first = mod.load_threshold_candidates(
        summary_path=summary,
        minimum_support=5,
        q_column="q_rate_enrichment_global",
        q_threshold=0.05,
        max_candidates=0,
        csubst_nonsyn_recode="no",
        pdb="none",
    )
    monkeypatch.setattr(mod, "analysis_engine_signature", lambda: "new-engine")
    monkeypatch.setattr(mod, "csubst_version", lambda: "new-version")
    second = mod.load_threshold_candidates(
        summary_path=summary,
        minimum_support=5,
        q_column="q_rate_enrichment_global",
        q_threshold=0.05,
        max_candidates=0,
        csubst_nonsyn_recode="no",
        pdb="none",
    )
    changed_parameter = mod.load_threshold_candidates(
        summary_path=summary,
        minimum_support=5,
        q_column="q_rate_enrichment_global",
        q_threshold=0.05,
        max_candidates=0,
        csubst_nonsyn_recode="dayhoff6",
        pdb="none",
    )

    assert first.loc[0, "_analysis_key"] == second.loc[0, "_analysis_key"]
    assert first.loc[0, "_analysis_key"] != changed_parameter.loc[0, "_analysis_key"]


def test_trait_color_paths_do_not_collide_after_filename_sanitizing(tmp_path):
    mod = load_module()
    trait_file = tmp_path / "traits.tsv"
    pd.DataFrame(
        {
            "species": ["sp1", "sp2"],
            "trait/a": [1, 0],
            "trait?a": [0, 1],
        }
    ).to_csv(trait_file, sep="\t", index=False)

    paths = mod.write_trait_color_tables(
        trait_file,
        ["trait/a", "trait?a"],
        tmp_path / "colors",
    )

    assert paths["trait/a"] != paths["trait?a"]
    first = pd.read_csv(paths["trait/a"], sep="\t")
    second = pd.read_csv(paths["trait?a"], sep="\t")
    assert first["color"].tolist() == ["firebrick", "black"]
    assert second["color"].tolist() == ["black", "firebrick"]


def test_exclusive_run_lock_recovers_a_dead_same_host_owner(tmp_path):
    mod = load_module()
    lock_path = tmp_path / ".candidate.lock"
    lock_directory = Path(f"{lock_path}.d")
    lock_directory.mkdir()
    (lock_directory / "owner.json").write_text(
        json.dumps(
            {
                "hostname": mod.socket.gethostname(),
                "pid": 2_000_000_000,
                "token": "dead-owner",
                "created_unix": 0,
            }
        ),
        encoding="utf-8",
    )

    with mod.exclusive_run_lock(lock_path):
        current_owner = json.loads(
            (lock_directory / "owner.json").read_text(encoding="utf-8")
        )
        assert current_owner["pid"] == os.getpid()
        assert current_owner["token"] != "dead-owner"

    assert not lock_directory.exists()


def test_run_lock_uses_heartbeat_to_detect_a_stale_foreign_owner(tmp_path):
    mod = load_module()
    lock_directory = tmp_path / ".candidate.lock.d"
    lock_directory.mkdir()
    owner_path = lock_directory / "owner.json"
    owner_path.write_text(
        json.dumps(
            {
                "hostname": "another-container",
                "pid": 1,
                "token": "foreign-owner",
                "created_unix": mod.time.time(),
            }
        ),
        encoding="utf-8",
    )

    assert not mod.run_lock_is_stale(lock_directory)
    stale_time = mod.time.time() - mod.RUN_LOCK_STALE_SECONDS - 1
    os.utime(owner_path, (stale_time, stale_time))
    assert mod.run_lock_is_stale(lock_directory)


def test_load_threshold_candidates_rejects_invalid_branch_ids(tmp_path):
    mod = load_module()
    summary = tmp_path / "summary.tsv"
    frame = candidate_rows().iloc[[0]].copy()
    frame["support_branch_ids"] = "3,bad"
    write_summary(summary, frame)

    with pytest.raises(ValueError, match="Invalid branch ID"):
        mod.load_threshold_candidates(
            summary_path=summary,
            minimum_support=5,
            q_column="q_rate_enrichment_global",
            q_threshold=0.05,
            max_candidates=0,
            csubst_nonsyn_recode="no",
            pdb="none",
        )


def test_required_report_input_preflight_reports_and_recovers_missing_stat_branch(
    tmp_path,
):
    mod = load_module()
    family_dir = tmp_path / "orthogroup"
    for subdir in ("iqtree_anc", "stat_branch", "clipkit"):
        (family_dir / subdir).mkdir(parents=True)
    for orthogroup in ("OG0001", "OG0002"):
        (family_dir / "iqtree_anc" / f"{orthogroup}_iqtree.anc.zip").write_bytes(
            b"placeholder"
        )
        (family_dir / "clipkit" / f"{orthogroup}_cds.clipkit.fa").write_text(
            ">sp1\nAAA\n", encoding="utf-8"
        )
    (family_dir / "stat_branch" / "OG0001_stat.branch.tsv").write_text(
        "branch_id\n0\n", encoding="utf-8"
    )
    missing_stat = family_dir / "stat_branch" / "OG0002_stat.branch.tsv"
    missing_stat.write_bytes(b"")

    first = mod.inspect_required_report_inputs(
        family_dir, ["OG0001", "OG0002"]
    )

    assert first["OG0001"]["missing_required_inputs"] == []
    assert first["OG0002"]["missing_required_inputs"] == [
        "stat_branch/OG0002_stat.branch.tsv"
    ]
    first_signature = first["OG0002"]["required_input_signature"]

    missing_stat.write_text("branch_id\n0\n", encoding="utf-8")
    second = mod.inspect_required_report_inputs(family_dir, ["OG0002"])

    assert second["OG0002"]["missing_required_inputs"] == []
    assert second["OG0002"]["required_input_signature"] != first_signature


def make_candidate_cache(mod, cache_root, row):
    cache_dir = cache_root / row["_cache_name"]
    cache_dir.mkdir(parents=True)
    focused = cache_dir / f"{row['_candidate_id']}.focused_tree_site.pdf"
    mod.site_wrapper.create_pdf("Focused site tree", str(focused))
    site_dir = cache_dir / "csubst_sites" / "csubst.branch_id8,9"
    site_dir.mkdir(parents=True)
    (site_dir / "csubst.tsv").write_text(
        "codon_site_alignment\tOCNany2spe\n4\t2.0\n",
        encoding="utf-8",
    )
    mod.site_wrapper.create_pdf("Raw CSUBST sites summary", str(site_dir / "csubst.pdf"))
    pd.DataFrame(
        [
            {
                "output_kind": "site_table_tsv",
                "output_file": "csubst.tsv",
                "output_path": str((site_dir / "csubst.tsv").resolve()),
                "file_exists": "Y",
                "file_size_bytes": (site_dir / "csubst.tsv").stat().st_size,
            },
            {
                "output_kind": "site_summary_pdf",
                "output_file": "csubst.pdf",
                "output_path": str((site_dir / "csubst.pdf").resolve()),
                "file_exists": "Y",
                "file_size_bytes": (site_dir / "csubst.pdf").stat().st_size,
            },
            {
                "output_kind": "output_manifest",
                "output_file": "csubst.outputs.tsv",
                "output_path": str((site_dir / "csubst.outputs.tsv").resolve()),
                "file_exists": "Y",
                "file_size_bytes": 0,
            },
        ]
    ).to_csv(site_dir / "csubst.outputs.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "analysis_key": row["_analysis_key"],
                "required_input_signature": row["_required_input_signature"],
                "candidate_id": row["_candidate_id"],
            }
        ]
    ).to_csv(cache_dir / "analysis.complete.tsv", sep="\t", index=False)


def test_package_threshold_writes_self_contained_zip(monkeypatch, tmp_path):
    mod = load_module()
    summary = tmp_path / "summary.tsv"
    write_summary(summary, candidate_rows().iloc[[1]].copy())
    candidates = mod.load_threshold_candidates(
        summary_path=summary,
        minimum_support=5,
        q_column="q_rate_enrichment_global",
        q_threshold=0.05,
        max_candidates=0,
        csubst_nonsyn_recode="no",
        pdb="none",
    )
    candidates["_required_input_signature"] = "test-input-signature"
    row = candidates.iloc[0]
    cache_root = tmp_path / "cache"
    make_candidate_cache(mod, cache_root, row)
    archive = tmp_path / "candidate_sites_min_support_5.zip"

    mod.package_threshold(
        candidates=candidates,
        threshold=5,
        source_summary=summary,
        archive_path=archive,
        packages_root=tmp_path / "packages",
        cache_root=cache_root,
        q_column="q_rate_enrichment_global",
        q_threshold=0.05,
    )

    assert archive.is_file()
    with zipfile.ZipFile(archive) as zipped:
        names = zipped.namelist()
        roots = {name.split("/", 1)[0] for name in names if name}
        assert roots == {archive.stem}
        assert f"{archive.stem}/candidate_manifest.tsv" in names
        assert f"{archive.stem}/skipped_candidates.tsv" in names
        assert f"{archive.stem}/package_metadata.tsv" in names
        candidate_prefix = f"{archive.stem}/candidate_0001_{row['_candidate_id']}"
        assert f"{candidate_prefix}/candidate.tsv" in names
        assert f"{candidate_prefix}/{row['_candidate_id']}.focused_tree_site.pdf" in names
        assert f"{candidate_prefix}/{row['_candidate_id']}.report.pdf" in names
        candidate_table = pd.read_csv(zipped.open(f"{candidate_prefix}/candidate.tsv"), sep="\t")
        assert candidate_table.loc[0, "selection_min_support"] == 5
        assert candidate_table.loc[0, "selection_q_column"] == "q_rate_enrichment_global"
        assert candidate_table.loc[0, "besthit_0.05"] == "protein A"
        output_manifest = pd.read_csv(
            zipped.open(f"{candidate_prefix}/csubst_sites/csubst.branch_id8,9/csubst.outputs.tsv"),
            sep="\t",
        )
        assert output_manifest["output_path"].tolist() == [
            "csubst.tsv",
            "csubst.pdf",
            "csubst.outputs.tsv",
        ]
        assert output_manifest["file_exists"].tolist() == ["Y", "Y", "Y"]
        self_row = output_manifest.loc[output_manifest["output_kind"] == "output_manifest"].iloc[0]
        zipped_manifest_info = zipped.getinfo(
            f"{candidate_prefix}/csubst_sites/csubst.branch_id8,9/csubst.outputs.tsv"
        )
        assert int(self_row["file_size_bytes"]) == zipped_manifest_info.file_size
    assert mod.archive_matches_source(archive, summary)

    with monkeypatch.context() as context:
        context.setattr(mod, "analysis_engine_signature", lambda: "changed-engine")
        context.setattr(mod, "csubst_version", lambda: "changed-csubst")
        context.setattr(mod, "runtime_dependency_versions", lambda: "changed-runtime")
        assert mod.archive_matches_source(archive, summary)

    original_summary = summary.read_text(encoding="utf-8")
    summary.write_text(original_summary + "\n", encoding="utf-8")
    assert not mod.archive_matches_source(archive, summary)
    summary.write_text(original_summary, encoding="utf-8")
    assert mod.archive_matches_source(archive, summary)

    original_archive = archive.read_bytes()
    raw_damaged = tmp_path / "raw-damaged.zip"
    with zipfile.ZipFile(archive, "r") as source_zip, zipfile.ZipFile(
        raw_damaged, "w"
    ) as target_zip:
        for member in source_zip.infolist():
            if member.filename.endswith("/csubst.tsv"):
                continue
            target_zip.writestr(member, source_zip.read(member.filename))
    raw_damaged.replace(archive)
    assert not mod.archive_matches_source(archive, summary)
    archive.write_bytes(original_archive)
    assert mod.archive_matches_source(archive, summary)

    damaged = tmp_path / "damaged.zip"
    with zipfile.ZipFile(archive, "r") as source_zip, zipfile.ZipFile(damaged, "w") as target_zip:
        for member in source_zip.infolist():
            if member.filename.endswith(f"{row['_candidate_id']}.report.pdf"):
                continue
            target_zip.writestr(member, source_zip.read(member.filename))
    damaged.replace(archive)
    assert not mod.archive_matches_source(archive, summary)


def test_package_threshold_records_candidates_skipped_for_missing_inputs(tmp_path):
    mod = load_module()
    summary = tmp_path / "summary.tsv"
    write_summary(summary, candidate_rows().iloc[[0]].copy())
    selected = mod.load_threshold_candidates(
        summary_path=summary,
        minimum_support=5,
        q_column="q_rate_enrichment_global",
        q_threshold=0.05,
        max_candidates=0,
        csubst_nonsyn_recode="no",
        pdb="none",
    )
    selected["_missing_required_inputs"] = (
        "stat_branch/OG0002_stat.branch.tsv"
    )
    selected["_required_input_signature"] = "missing-stat-branch"
    skipped = mod.skipped_candidate_frame(selected)
    eligible = selected.iloc[0:0].copy()
    archive = tmp_path / "skipped.zip"

    mod.package_threshold(
        candidates=eligible,
        threshold=5,
        source_summary=summary,
        archive_path=archive,
        packages_root=tmp_path / "packages",
        cache_root=tmp_path / "cache",
        q_column="q_rate_enrichment_global",
        q_threshold=0.05,
        skipped_candidates=skipped,
    )

    assert mod.archive_matches_source(
        archive,
        summary,
        expected_candidates=eligible,
        expected_skipped_candidates=skipped,
    )
    with zipfile.ZipFile(archive) as zipped:
        root = archive.stem
        candidate_manifest = pd.read_csv(
            zipped.open(f"{root}/candidate_manifest.tsv"), sep="\t"
        )
        skipped_manifest = pd.read_csv(
            zipped.open(f"{root}/skipped_candidates.tsv"), sep="\t"
        )
        metadata = pd.read_csv(
            zipped.open(f"{root}/package_metadata.tsv"), sep="\t"
        )
        assert candidate_manifest.empty
        assert skipped_manifest.loc[0, "orthogroup"] == "OG0002"
        assert skipped_manifest.loc[0, "reason_code"] == "missing_required_input"
        assert (
            skipped_manifest.loc[0, "missing_required_inputs"]
            == "stat_branch/OG0002_stat.branch.tsv"
        )
        assert metadata.loc[0, "selected_candidate_count"] == 1
        assert metadata.loc[0, "packaged_candidate_count"] == 0
        assert metadata.loc[0, "skipped_candidate_count"] == 1
        assert metadata.loc[0, "skipped_gene_family_count"] == 1


def test_archive_names_record_selection_and_optional_analysis_modes(tmp_path):
    mod = load_module()

    path = mod.archive_path_for_threshold(
        summary_prefix=tmp_path / "orthogroup_csubst_aa_change",
        out_dir=tmp_path,
        threshold=7,
        q_column="q_rate_enrichment_global",
        q_threshold=0.05,
        max_candidates=20,
        nonsyn_recode="dayhoff6",
        pdb="besthit",
    )

    assert path.name == (
        "orthogroup_csubst_aa_change_candidate_sites_min_support_7_"
        "q_rate_enrichment_global_le_0.05_top20_nonsynRecode-dayhoff6_pdb-besthit.zip"
    )


def test_validate_args_rejects_invalid_q_threshold(tmp_path):
    mod = load_module()
    trait = tmp_path / "trait.tsv"
    trait.write_text("species\taquatic\nsp1\t1\n", encoding="utf-8")
    families = tmp_path / "orthogroup"
    families.mkdir()
    args = argparse.Namespace(
        min_support=5,
        q_threshold=1.1,
        max_candidates=0,
        q_column="q_rate_enrichment_global",
        ncpu=1,
        dir_orthogroup=str(families),
        file_trait=str(trait),
        summary_prefix=str(tmp_path / "summary"),
        out_dir=str(tmp_path / "out"),
    )

    with pytest.raises(ValueError, match="between 0 and 1"):
        mod.validate_args(args)


def test_safe_extract_zip_rejects_parent_path_traversal(tmp_path):
    mod = load_module()
    archive = tmp_path / "unsafe.zip"
    with zipfile.ZipFile(archive, "w") as zipped:
        zipped.writestr("../escaped.txt", "unsafe")

    with pytest.raises(ValueError, match="Unsafe ZIP archive"):
        mod.safe_extract_zip(archive, tmp_path / "extract", "OG0001.iqtree.anc")

    assert not (tmp_path / "escaped.txt").exists()


def test_make_csubst_manifests_portable_rejects_external_files(tmp_path):
    mod = load_module()
    candidate_dir = tmp_path / "candidate"
    site_dir = candidate_dir / "csubst_sites" / "csubst.branch_id1"
    site_dir.mkdir(parents=True)
    outside = tmp_path / "outside.tsv"
    outside.write_text("x\n", encoding="utf-8")
    pd.DataFrame(
        [
            {
                "output_kind": "site_table_tsv",
                "output_file": "../../../outside.tsv",
                "output_path": str(outside),
                "file_exists": "Y",
                "file_size_bytes": outside.stat().st_size,
            }
        ]
    ).to_csv(site_dir / "csubst.outputs.tsv", sep="\t", index=False)

    mod.make_csubst_manifests_portable(candidate_dir)

    manifest = pd.read_csv(site_dir / "csubst.outputs.tsv", sep="\t", dtype=str)
    assert pd.isna(manifest.loc[0, "output_path"])
    assert manifest.loc[0, "file_exists"] == "N"
    assert manifest.loc[0, "file_size_bytes"] == "-1"


def test_package_threshold_writes_valid_empty_zip(tmp_path):
    mod = load_module()
    summary = tmp_path / "summary.tsv"
    write_summary(summary)
    archive = tmp_path / "empty.zip"
    empty = mod.load_threshold_candidates(
        summary_path=summary,
        minimum_support=5,
        q_column="q_rate_enrichment_global",
        q_threshold=0.0,
        max_candidates=0,
        csubst_nonsyn_recode="no",
        pdb="none",
    )

    mod.package_threshold(
        candidates=empty,
        threshold=5,
        source_summary=summary,
        archive_path=archive,
        packages_root=tmp_path / "packages",
        cache_root=tmp_path / "cache",
        q_column="q_rate_enrichment_global",
        q_threshold=0.0,
    )

    assert mod.archive_matches_source(archive, summary)
    with zipfile.ZipFile(archive) as zipped:
        manifest = pd.read_csv(zipped.open(f"{archive.stem}/candidate_manifest.tsv"), sep="\t")
        assert manifest.empty
        assert manifest.columns.tolist() == mod.CANDIDATE_MANIFEST_COLUMNS


def make_run_args(tmp_path):
    family_dir = tmp_path / "orthogroup"
    family_dir.mkdir(exist_ok=True)
    trait_file = tmp_path / "species_trait.tsv"
    trait_file.write_text("species\taquatic\nsp1\t1\n", encoding="utf-8")
    return argparse.Namespace(
        summary_prefix=str(tmp_path / "orthogroup_csubst_aa_change"),
        dir_orthogroup=str(family_dir),
        file_trait=str(trait_file),
        out_dir=str(tmp_path / "out"),
        min_support=5,
        q_column="q_rate_enrichment_global",
        q_threshold=0.05,
        max_candidates=0,
        ncpu=2,
        csubst_nonsyn_recode="no",
        pdb="none",
    )


def write_run_summaries(tmp_path):
    frame = candidate_rows().copy()
    frame["q_rate_enrichment_global"] = 0.01
    for threshold in (5, 6, 7):
        write_summary(
            tmp_path
            / f"orthogroup_csubst_aa_change_min_support_{threshold}_summary.tsv",
            frame.loc[frame["support_unit_count"] >= threshold].copy(),
        )


def available_input_states(_dir_orthogroup, orthogroups):
    return {
        orthogroup: {
            "missing_required_inputs": [],
            "required_input_signature": f"available-{orthogroup}",
        }
        for orthogroup in orthogroups
    }


def test_run_processes_thresholds_in_descending_order(monkeypatch, tmp_path):
    mod = load_module()
    write_run_summaries(tmp_path)
    args = make_run_args(tmp_path)
    analyzed = []
    packaged = []

    def fake_ensure(candidates, **kwargs):
        analyzed.append(int(candidates["_selection_min_support"].iloc[0]))
        return []

    def fake_package(candidates, threshold, archive_path, **kwargs):
        packaged.append((threshold, candidates.shape[0]))
        Path(archive_path).parent.mkdir(parents=True, exist_ok=True)
        Path(archive_path).write_bytes(b"placeholder")

    monkeypatch.setattr(mod, "ensure_candidate_analyses", fake_ensure)
    monkeypatch.setattr(mod, "package_threshold", fake_package)
    monkeypatch.setattr(mod, "archive_matches_source", lambda *args, **kwargs: False)
    monkeypatch.setattr(mod, "inspect_required_report_inputs", available_input_states)

    manifest_path = mod.run(args)

    manifest = pd.read_csv(manifest_path, sep="\t")
    assert analyzed == [7, 6, 5]
    assert packaged == [(7, 1), (6, 2), (5, 3)]
    assert manifest["min_support"].tolist() == [7, 6, 5]
    assert manifest["status"].tolist() == ["completed", "completed", "completed"]


def test_run_skips_missing_inputs_and_packages_remaining_candidates(
    monkeypatch, tmp_path
):
    mod = load_module()
    write_run_summaries(tmp_path)
    args = make_run_args(tmp_path)
    analyzed = []
    packaged = []

    def input_states(_dir_orthogroup, orthogroups):
        return {
            orthogroup: {
                "missing_required_inputs": (
                    ["stat_branch/OG0002_stat.branch.tsv"]
                    if orthogroup == "OG0002"
                    else []
                ),
                "required_input_signature": f"state-{orthogroup}",
            }
            for orthogroup in orthogroups
        }

    def fake_ensure(candidates, **kwargs):
        analyzed.extend(candidates["orthogroup"].tolist())
        return []

    def fake_package(candidates, threshold, archive_path, skipped_candidates, **kwargs):
        packaged.append(
            (
                threshold,
                candidates["orthogroup"].tolist(),
                skipped_candidates["orthogroup"].tolist(),
            )
        )
        Path(archive_path).parent.mkdir(parents=True, exist_ok=True)
        Path(archive_path).write_bytes(b"placeholder")

    monkeypatch.setattr(mod, "inspect_required_report_inputs", input_states)
    monkeypatch.setattr(mod, "ensure_candidate_analyses", fake_ensure)
    monkeypatch.setattr(mod, "package_threshold", fake_package)
    monkeypatch.setattr(mod, "archive_matches_source", lambda *args, **kwargs: False)

    manifest_path = mod.run(args)

    manifest = pd.read_csv(manifest_path, sep="\t")
    assert packaged == [
        (7, ["OG0003"], []),
        (6, ["OG0001", "OG0003"], []),
        (5, ["OG0001", "OG0003"], ["OG0002"]),
    ]
    assert "OG0002" not in analyzed
    assert manifest["candidate_count"].tolist() == [1, 2, 3]
    assert manifest["packaged_candidate_count"].tolist() == [1, 2, 2]
    assert manifest["skipped_candidate_count"].tolist() == [0, 0, 1]
    assert manifest["skipped_gene_family_count"].tolist() == [0, 0, 1]
    assert manifest["status"].tolist() == [
        "completed",
        "completed",
        "completed_with_skips",
    ]
    skipped_path = tmp_path / "out" / manifest.loc[0, "skipped_candidates_tsv"]
    skipped = pd.read_csv(skipped_path, sep="\t")
    assert skipped["orthogroup"].tolist() == ["OG0002"]
    assert skipped["reason_code"].tolist() == ["missing_required_input"]
    assert skipped["missing_required_inputs"].tolist() == [
        "stat_branch/OG0002_stat.branch.tsv"
    ]


def test_run_records_batch_setup_failure_in_manifest(monkeypatch, tmp_path):
    mod = load_module()
    write_run_summaries(tmp_path)
    args = make_run_args(tmp_path)

    def fail_analysis(**kwargs):
        raise RuntimeError("materialization failed")

    monkeypatch.setattr(mod, "ensure_candidate_analyses", fail_analysis)
    monkeypatch.setattr(mod, "archive_matches_source", lambda *args, **kwargs: False)
    monkeypatch.setattr(mod, "inspect_required_report_inputs", available_input_states)

    with pytest.raises(RuntimeError, match="packaging failed"):
        mod.run(args)

    manifests = list((tmp_path / "out").glob("*_manifest.tsv"))
    assert len(manifests) == 1
    manifest = pd.read_csv(manifests[0], sep="\t")
    assert manifest.loc[0, "min_support"] == 7
    assert manifest.loc[0, "status"] == "failed"
    assert manifest.loc[0, "error"] == "materialization failed"
    assert manifest.loc[1:, "status"].tolist() == ["pending", "pending"]


def test_run_keeps_candidate_analysis_errors_as_hard_failures(monkeypatch, tmp_path):
    mod = load_module()
    write_run_summaries(tmp_path)
    args = make_run_args(tmp_path)
    packaged = []

    def failed_analysis(candidates, **kwargs):
        return [
            {
                "candidate_id": candidates.iloc[0]["_candidate_id"],
                "status": "failed",
                "error": "CSUBST/stat_branch branch identity mismatch",
            }
        ]

    monkeypatch.setattr(mod, "inspect_required_report_inputs", available_input_states)
    monkeypatch.setattr(mod, "ensure_candidate_analyses", failed_analysis)
    monkeypatch.setattr(
        mod, "package_threshold", lambda **kwargs: packaged.append(kwargs)
    )
    monkeypatch.setattr(mod, "archive_matches_source", lambda *args, **kwargs: False)

    with pytest.raises(RuntimeError, match="packaging failed"):
        mod.run(args)

    assert packaged == []
    manifest_path = next((tmp_path / "out").glob("*_manifest.tsv"))
    manifest = pd.read_csv(manifest_path, sep="\t")
    assert manifest.loc[0, "status"] == "failed"
    assert "branch identity mismatch" in manifest.loc[0, "error"]


def test_run_writes_empty_manifest_when_no_threshold_reaches_minimum(tmp_path):
    mod = load_module()
    for threshold in (2, 3, 4):
        write_summary(
            tmp_path
            / f"orthogroup_csubst_aa_change_min_support_{threshold}_summary.tsv"
        )
    args = make_run_args(tmp_path)

    manifest_path = mod.run(args)

    manifest = pd.read_csv(manifest_path, sep="\t")
    assert manifest.empty
    assert manifest.columns.tolist() == mod.ARCHIVE_MANIFEST_COLUMNS
    assert not list((tmp_path / "out").glob("*.zip"))
    assert not any((tmp_path / "out").glob(".*.work"))


def test_archived_family_inputs_are_materialized_once_per_orthogroup(
    monkeypatch, tmp_path
):
    mod = load_module()
    family_dir = tmp_path / "orthogroup"
    (family_dir / ".gg_store").mkdir(parents=True)
    materialized_dir = tmp_path / "materialized" / "OG0001"
    materialized_dir.mkdir(parents=True)
    events = []

    class FakeMaterializationDirectory:
        def __init__(self, parent, orthogroup):
            events.append(("lock", str(parent), orthogroup))
            self.name = str(materialized_dir)

        def cleanup(self):
            events.append(("cleanup",))

    def fake_materialize(dir_og, og, destination_root):
        events.append(("materialize", dir_og, og, destination_root))
        return []

    def fake_analyze(record, effective_dir_orthogroup, **kwargs):
        events.append(("analyze", record["_candidate_id"], effective_dir_orthogroup))
        return {"candidate_id": record["_candidate_id"], "status": "completed", "error": ""}

    monkeypatch.setattr(
        mod.site_wrapper,
        "LockedMaterializationDirectory",
        FakeMaterializationDirectory,
    )
    monkeypatch.setattr(mod.site_wrapper, "materialize_csubst_site_inputs", fake_materialize)
    monkeypatch.setattr(mod, "analyze_candidate", fake_analyze)
    records = [
        {"_candidate_id": "candidate1"},
        {"_candidate_id": "candidate2"},
    ]

    results = mod.analyze_orthogroup_batch(
        orthogroup="OG0001",
        records=records,
        cache_root=tmp_path / "cache",
        dir_orthogroup=str(family_dir),
        materialization_parent=tmp_path / "materialization_parent",
        trait_color_paths={},
        nonsyn_recode="no",
        pdb="none",
    )

    assert [result["candidate_id"] for result in results] == ["candidate1", "candidate2"]
    assert sum(event[0] == "materialize" for event in events) == 1
    assert sum(event[0] == "analyze" for event in events) == 2
    assert events[-1] == ("cleanup",)


def test_cli_end_to_end_reuses_analysis_across_threshold_zips(tmp_path):
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    csubst_log = tmp_path / "csubst.log"
    rscript_log = tmp_path / "rscript.log"
    fake_csubst = fake_bin / "csubst"
    fake_csubst.write_text(
        """#!/usr/bin/env python3
import os
import sys
from pathlib import Path
from reportlab.pdfgen import canvas

args = sys.argv[1:]
def value(flag):
    return args[args.index(flag) + 1]

branch_ids = value("--branch_id")
outdir = Path(value("--outdir")) / (value("--output_prefix") + ".branch_id" + branch_ids)
outdir.mkdir(parents=True, exist_ok=True)
site_tsv = outdir / "csubst.tsv"
site_pdf = outdir / "csubst.pdf"
site_tsv.write_text("codon_site_alignment\\tOCNany2spe\\n1\\t0.0\\n2\\t2.0\\n3\\t0.0\\n", encoding="utf-8")
pdf = canvas.Canvas(str(site_pdf))
pdf.drawString(72, 720, "Fake csubst sites summary")
pdf.save()
manifest = outdir / "csubst.outputs.tsv"
manifest.write_text(
    "output_kind\\toutput_file\\toutput_path\\tfile_exists\\tfile_size_bytes\\n"
    + f"site_table_tsv\\tcsubst.tsv\\t{site_tsv.resolve()}\\tY\\t{site_tsv.stat().st_size}\\n"
    + f"site_summary_pdf\\tcsubst.pdf\\t{site_pdf.resolve()}\\tY\\t{site_pdf.stat().st_size}\\n"
    + f"output_manifest\\tcsubst.outputs.tsv\\t{manifest.resolve()}\\tY\\t0\\n",
    encoding="utf-8",
)
with open(os.environ["FAKE_CSUBST_LOG"], "a", encoding="utf-8") as handle:
    handle.write(" ".join(args) + "\\n")
""",
        encoding="utf-8",
    )
    fake_rscript = fake_bin / "Rscript"
    fake_rscript.write_text(
        """#!/usr/bin/env python3
import os
import sys
from reportlab.pdfgen import canvas

pdf = canvas.Canvas("stat_branch2tree_plot.pdf")
pdf.drawString(72, 720, "Fake focused tree")
pdf.save()
with open(os.environ["FAKE_RSCRIPT_LOG"], "a", encoding="utf-8") as handle:
    handle.write(" ".join(sys.argv[1:]) + "\\n")
""",
        encoding="utf-8",
    )
    fake_csubst.chmod(0o755)
    fake_rscript.chmod(0o755)

    family_dir = tmp_path / "orthogroup"
    for subdir in ("iqtree_anc", "clipkit", "stat_branch", "rpsblast"):
        (family_dir / subdir).mkdir(parents=True)
    for og in ("OG0001", "OG0002"):
        iqtree_zip = family_dir / "iqtree_anc" / f"{og}_iqtree.anc.zip"
        with zipfile.ZipFile(iqtree_zip, "w") as archive:
            for filename, content in {
                "csubst.fasta": ">sp1\nGCTGTTGAT\n>sp2\nGCTGTTGAT\n>sp3\nGCTGTTGAT\n",
                "csubst.nwk": "(sp1:0.1,(sp2:0.1,sp3:0.1)RootedBC:0.1)RootedRoot;\n",
                "csubst.treefile": "(sp1:0.1,(sp2:0.1,sp3:0.1)IqtreeBC:0.1)IqtreeRoot;\n",
                "csubst.state": "placeholder\n",
                "csubst.rate": "placeholder\n",
                "csubst.iqtree": "placeholder\n",
                "csubst.log": "placeholder\n",
            }.items():
                archive.writestr(f"{og}.iqtree.anc/{filename}", content)
        (family_dir / "clipkit" / f"{og}_cds.clipkit.fa").write_text(
            ">sp1\nGCTGTTGAT\n>sp2\nGCTGTTGAT\n>sp3\nGCTGTTGAT\n",
            encoding="utf-8",
        )
        (family_dir / "stat_branch" / f"{og}_stat.branch.tsv").write_text(
            "branch_id\tnode_name\tgene_labels\n"
            "0\tsp1\tsp1\n"
            "1\tsp2\tsp2\n"
            "2\tsp3\tsp3\n"
            "3\tN1\tsp2; sp3\n"
            "4\tRoot\tsp1; sp2; sp3\n",
            encoding="utf-8",
        )
        (family_dir / "rpsblast" / f"{og}_rpsblast.tsv").write_text(
            "query\nsp1\n",
            encoding="utf-8",
        )

    summary_row = pd.DataFrame(
        [
            {
                "orthogroup": "OG0001",
                "trait": "aquatic",
                "state_change": "2V",
                "codon_site_alignment": 2,
                "support_unit_count": 6,
                "support_unit_ids": "1,2,3,4,5,6",
                "support_branch_ids": "0,1",
                "p_rate_enrichment": 0.001,
                "q_rate_enrichment_global": 0.01,
                "besthit_0.05": "annotated protein",
            },
            {
                "orthogroup": "OG0002",
                "trait": "aquatic",
                "state_change": "3D",
                "codon_site_alignment": 3,
                "support_unit_count": 6,
                "support_unit_ids": "1,2,3,4,5,6",
                "support_branch_ids": "2,3",
                "p_rate_enrichment": 0.002,
                "q_rate_enrichment_global": 0.02,
                "besthit_0.05": "second annotated protein",
            },
        ]
    )
    prefix = tmp_path / "orthogroup_csubst_aa_change"
    for threshold in (5, 6):
        summary_row.to_csv(
            tmp_path / f"{prefix.name}_min_support_{threshold}_summary.tsv",
            sep="\t",
            index=False,
        )
    trait_file = tmp_path / "species_trait.tsv"
    trait_file.write_text("species\taquatic\nsp1\t1\nsp2\t0\nsp3\t0\n", encoding="utf-8")
    output_dir = tmp_path / "out"
    env = os.environ.copy()
    env["PATH"] = str(fake_bin) + os.pathsep + env["PATH"]
    env["FAKE_CSUBST_LOG"] = str(csubst_log)
    env["FAKE_RSCRIPT_LOG"] = str(rscript_log)

    command = [
        sys.executable,
        str(SCRIPT_PATH),
        "--summary_prefix",
        str(prefix),
        "--dir_orthogroup",
        str(family_dir),
        "--file_trait",
        str(trait_file),
        "--out_dir",
        str(output_dir),
        "--min_support",
        "5",
        "--ncpu",
        "2",
        "--pdb",
        "none",
    ]
    processes = [
        subprocess.Popen(
            command,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        for _ in range(2)
    ]
    results = []
    for process in processes:
        stdout, stderr = process.communicate(timeout=60)
        results.append((process.returncode, stdout, stderr))

    assert all(returncode == 0 for returncode, _, _ in results), results
    combined_stdout = "\n".join(stdout for _, stdout, _ in results)
    assert "candidates packaged" in combined_stdout
    assert "existing ZIP retained" in combined_stdout
    assert len(csubst_log.read_text(encoding="utf-8").splitlines()) == 2
    assert "--pdb" not in csubst_log.read_text(encoding="utf-8")
    rscript_lines = rscript_log.read_text(encoding="utf-8").splitlines()
    assert len(rscript_lines) == 2
    assert any("amino_acid_site,1,2," in line for line in rscript_lines)
    assert any("amino_acid_site,1,3," in line for line in rscript_lines)
    run_manifest = pd.read_csv(next(output_dir.glob("*_manifest.tsv")), sep="\t")
    assert run_manifest["min_support"].tolist() == [6, 5]
    assert run_manifest["status"].tolist() == ["existing", "existing"]
    archives = sorted(output_dir.glob("*.zip"))
    assert len(archives) == 2
    from pypdf import PdfReader

    for archive_path in archives:
        assert archive_path.is_file()
        with zipfile.ZipFile(archive_path) as archive:
            assert archive.testzip() is None
            candidate_manifest_path = f"{archive_path.stem}/candidate_manifest.tsv"
            candidate_manifest = pd.read_csv(archive.open(candidate_manifest_path), sep="\t")
            assert candidate_manifest.shape[0] == 2
            for _, candidate in candidate_manifest.iterrows():
                report_path = f"{archive_path.stem}/{candidate['report_pdf']}"
                assert len(PdfReader(io.BytesIO(archive.read(report_path))).pages) == 3
                raw_manifest_path = (
                    f"{archive_path.stem}/{candidate['csubst_sites_dir']}"
                    f"/csubst.branch_id{candidate['support_branch_ids']}/csubst.outputs.tsv"
                )
                raw_manifest = pd.read_csv(archive.open(raw_manifest_path), sep="\t")
                assert raw_manifest["output_path"].tolist() == [
                    "csubst.tsv",
                    "csubst.pdf",
                    "csubst.outputs.tsv",
                ]
                self_row = raw_manifest.loc[
                    raw_manifest["output_kind"] == "output_manifest"
                ].iloc[0]
                assert int(self_row["file_size_bytes"]) == archive.getinfo(
                    raw_manifest_path
                ).file_size
    assert not any(output_dir.glob(".*.work"))
    assert not any(output_dir.glob(".*.lock*"))

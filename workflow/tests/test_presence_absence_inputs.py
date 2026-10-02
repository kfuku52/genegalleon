"""Saved-family and conservative query-selection contracts (no external tools)."""
import csv
from pathlib import Path

import pytest

from workflow.support import query_gene_orthologs as mod
from workflow.support.gene_family_output_store import archive_completed_outputs, family_context


def row(node, parent, left, right, event, name, species="", marker=""):
    return dict(branch_id=node, parent=parent, child1=left, child2=right,
                so_event=event, node_name=name, spnode_coverage=species,
                query_marker_source=marker, support_unrooted="")


def query(query_id, source="", **extra):
    return dict(query_id=query_id, query_label=query_id, source_species=source, **extra)


def write_tsv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def gene_rows():
    return [
        row(0, -1, 1, 2, "D", "root"),
        row(1, 0, 3, 4, "S", "lineage1"),
        row(2, 0, 5, 6, "S", "lineage2"),
        row(3, 1, -1, -1, "L", "Anchor_one_a", "Anchor_one", "best:near|best:unknown|direct:near_alias"),
        row(4, 1, -1, -1, "L", "Anchor_two_a", "Anchor_two", "best:far|best:far_alias"),
        row(5, 2, -1, -1, "L", "Anchor_one_b", "Anchor_one", "best:paralog"),
        row(6, 2, -1, -1, "L", "Anchor_two_b", "Anchor_two"),
    ]


def species_tree(tmp_path):
    path = tmp_path / "species.nwk"
    path.write_text("((Target_species,Near_species),Far_species);\n")
    return str(path)


def test_source_metadata_is_explicit_and_preserved(tmp_path):
    path = tmp_path / "query.fa"
    path.write_text(">q1 | Display | species=Near species evidence=curated\nAA\n"
                    ">q2 | Unknown | Near species\nAA\n")
    definitions = mod.read_query_definitions(path)
    assert definitions[0]["source_species"] == "Near_species"
    assert definitions[0]["query_label"] == "Display"
    assert definitions[1]["source_species"] == ""
    metadata = tmp_path / "meta.tsv"
    metadata.write_text("family_id\tquery_id\tsource_species\ttree_species\n"
                        "F\tq1\tNear relative\tNear species\n")
    assert mod.read_query_metadata(metadata)[("F", "q1")] == {
        "source_species": "Near_relative", "tree_species": "Near_species"}


def test_closest_uses_original_source_not_best_hit_and_keeps_paralogs_and_unknown(tmp_path):
    definitions = [query("near", "Near_species"), query("far", "Far_species"),
                   query("paralog", "Far_species"), query("unknown")]
    retained, audit = mod.select_closest_queries(
        gene_rows(), definitions, "F", species_tree(tmp_path), "Target_species", policy="closest")
    assert [d["query_id"] for d in retained] == ["near", "paralog", "unknown"]
    assert audit[1]["decision"] == "excluded"
    assert audit[1]["replaced_by"] == "near"
    assert audit[0]["anchor_cds_fasta_id"] == "Anchor_one_a"
    assert audit[0]["distance"] == 1
    assert audit[-1]["reason"] == "unknown_source_or_tree_species"


def test_equal_distances_and_all_policy_do_not_drop_records(tmp_path):
    definitions = [query("near", "Near_species"), query("far", "Near_species")]
    retained, _audit = mod.select_closest_queries(
        gene_rows(), definitions, "F", species_tree(tmp_path), "Target_species", policy="closest")
    assert retained == definitions
    retained, audit = mod.select_closest_queries(gene_rows(), definitions, "F", "", "")
    assert retained == definitions
    assert {a["reason"] for a in audit} == {"all_queries"}


def test_nearness_depends_on_shared_ancestry_not_source_clade_sampling(tmp_path):
    path = tmp_path / "unbalanced.nwk"
    path.write_text("((Target_species,((Near_species,Other_one),Other_two)),Far_species);\n")
    definitions = [query("near", "Near_species"), query("far", "Far_species")]
    retained, audit = mod.select_closest_queries(
        gene_rows(), definitions, "F", str(path), "Target_species", policy="closest")
    assert [d["query_id"] for d in retained] == ["near"]
    assert [a["distance"] for a in audit] == [1, 2]


def test_farther_preduplication_query_is_kept_when_nearer_query_loses_coverage(tmp_path):
    rows = [row(0, -1, 1, 2, "S", "root"),
            row(1, 0, -1, -1, "L", "Anchor_two_a", "Anchor_two", "best:far"),
            row(2, 0, 3, 4, "D", "duplication"),
            row(3, 2, -1, -1, "L", "Anchor_one_a", "Anchor_one", "best:near"),
            row(4, 2, -1, -1, "L", "Anchor_one_b", "Anchor_one")]
    definitions = [query("near", "Near_species"), query("far", "Far_species")]
    retained, _audit = mod.select_closest_queries(
        rows, definitions, "F", species_tree(tmp_path), "Target_species", policy="closest")
    assert retained == definitions  # Orthology is not transitive across the duplication.


def test_unknown_tree_alias_kept_and_explicit_alias_ranked(tmp_path):
    definitions = [query("near", "Unsequenced_relative", tree_species="Near_species"),
                   query("far", "Far_species"), query("unknown", "No_tree_tip")]
    retained, audit = mod.select_closest_queries(
        gene_rows(), definitions, "F", species_tree(tmp_path), "Target_species", policy="closest")
    assert [d["query_id"] for d in retained] == ["near", "unknown"]
    assert audit[0]["source_species"] == "Unsequenced_relative"
    assert audit[0]["tree_species"] == "Near_species"


def test_selection_species_may_have_no_family_members(tmp_path):
    definitions = [query("near", "Near_species"), query("far", "Far_species")]
    retained, _audit = mod.select_closest_queries(
        gene_rows(), definitions, "F", species_tree(tmp_path), "Target_species",
        selection_species="Target_species", policy="closest")
    assert [definition["query_id"] for definition in retained] == ["near"]
    with pytest.raises(ValueError, match="absent from species tree"):
        mod.select_closest_queries(gene_rows(), definitions, "F", species_tree(tmp_path),
                                  "Target_species", selection_species="Typo_species", policy="closest")


def saved_manifest(tmp_path):
    store = tmp_path / "saved"
    (store / "stat_branch").mkdir(parents=True)
    (store / "cds_fasta").mkdir()
    rows = gene_rows()
    write_tsv(store / "stat_branch/OG1_stat.branch.tsv", rows)
    (store / "cds_fasta/OG1_cds.fasta").write_text("".join(
        f">{r['node_name']}\nATG\n" for r in rows if r["so_event"] == "L"))
    hog = tmp_path / "N0.tsv"
    hog.write_text("HOG\tOG\tGene Tree Parent Clade\tAnchor_one\tAnchor_two\n"
                   "HOG1\tOG1\tn1\ta\tb\nHOG2\tOG1\tn2\tb\ta\n")
    manifest = tmp_path / "families.tsv"
    manifest.write_text("family_id\tsource_dir\tsource_family_id\tanchor_species\thog_table\thog_ids\n"
                        "Combined\tsaved\tOG1\tAnchor_one\tN0.tsv\tHOG1;HOG2\n")
    return store, manifest


def test_manifest_uses_full_saved_tree_and_all_anchor_species_hog_members(tmp_path):
    store, path = saved_manifest(tmp_path)
    records = mod.read_family_manifest(path)
    columns, glyphs, tree, mapping = mod.collect_query_anchor_orthologs(
        store, "", manifest_records=records, query_label="label")
    assert len(columns) == 2
    assert {m["hog_ids"] for m in mapping} == {"HOG1", "HOG2"}
    assert {m["source_species"] for m in mapping} == {"Anchor_one"}
    assert {g["species"] for g in glyphs} == {"Anchor_one", "Anchor_two"}
    assert any(t["event"] == "D" for t in tree)
    assert {c["family_id"] for c in columns} == {"Combined"}
    # Crossed HOG membership must not replace reconciled S/D orthology.
    other = [g for g in glyphs if g["species"] == "Anchor_two"]
    assert {g["reference_cds_fasta_ids"]: g["gene_ids"] for g in other} == {
        "Anchor_one_a": "Anchor_two_a", "Anchor_one_b": "Anchor_two_b"}


@pytest.mark.parametrize("zipped", [False, True])
def test_manifest_combines_query_file_and_hog_sources_including_zip(tmp_path, zipped):
    store, path = saved_manifest(tmp_path)
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    query_file = query_dir / "OG1"
    query_file.write_text(">near | Native query | species=Near_species\nAA\n")
    if zipped:
        for subdir, suffix in (("stat_tree", "_stat.tree.tsv"), ("tree_plot", "_tree_plot.pdf")):
            folder = store / subdir
            folder.mkdir()
            (folder / f"OG1{suffix}").write_text("fixture\n")
        families, from_name = family_context("query2family", query_dir=query_dir)
        archive_completed_outputs(store, "query2family", families, from_name)
        assert not (store / "stat_branch").exists()
    # Both artifact bases are read-only, independent of the plotting block names.
    path.write_text("family_id\tsource_dir\tsource_family_id\tquery_file\tanchor_species\thog_table\thog_ids\n"
                    "Native\tsaved\tOG1\tquery_gene/OG1\t\t\t\n"
                    "HOGs\tsaved\tOG1\t\tAnchor_one\tN0.tsv\tHOG1;HOG2\n")
    records = mod.read_family_manifest(path)
    columns, glyphs, _tree, mapping = mod.collect_query_anchor_orthologs(
        store, "", manifest_records=records, query_label="label")
    assert [c["family_id"] for c in columns] == ["Native", "HOGs", "HOGs"]
    assert columns[0]["plot_label"] == "Native query"
    assert mapping[0]["source_species"] == "Near_species"
    assert {g["family_id"] for g in glyphs} == {"Native", "HOGs"}


@pytest.mark.parametrize("metadata", [
    {("F", "near"): {"source_species": "Far_species"}},
    {("F", "typo"): {"source_species": "Near_species"}},
])
def test_conflicting_or_stale_selected_query_metadata_fails(tmp_path, metadata):
    store, _manifest = saved_manifest(tmp_path)
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "F").write_text(">near | Native query | species=Near_species\nAA\n")
    records = [dict(family_id="F", source_dir=str(store), source_family_id="OG1",
                    query_file=str(query_dir / "F"))]
    with pytest.raises(ValueError, match="metadata"):
        mod.collect_query_anchor_orthologs(store, "", manifest_records=records,
                                         query_metadata=metadata)


@pytest.mark.parametrize("error", ["wrong_og", "missing_hog", "missing_gene", "duplicate_block",
                                   "wrong_species", "duplicate_tip", "missing_cds"])
def test_manifest_membership_and_source_errors_fail_closed(tmp_path, error):
    store, path = saved_manifest(tmp_path)
    if error == "wrong_og":
        hog = tmp_path / "N0.tsv"
        hog.write_text(hog.read_text().replace("HOG1\tOG1", "HOG1\tOG2"))
    elif error == "missing_hog":
        path.write_text(path.read_text().replace("HOG1;HOG2", "HOG1;HOG99"))
    elif error == "missing_gene":
        hog = tmp_path / "N0.tsv"
        hog.write_text(hog.read_text().replace("n1\ta", "n1\tnot_in_tree"))
    elif error == "duplicate_block":
        path.write_text(path.read_text() + path.read_text().splitlines()[1] + "\n")
    elif error == "wrong_species":
        hog = tmp_path / "N0.tsv"
        hog.write_text(hog.read_text().replace("n1\ta", "n1\tAnchor_two_b"))
    elif error == "duplicate_tip":
        stat = store / "stat_branch/OG1_stat.branch.tsv"
        stat.write_text(stat.read_text().replace("Anchor_one_b", "Anchor_one_a"))
    else:
        cds = store / "cds_fasta/OG1_cds.fasta"
        cds.write_text(cds.read_text().replace(">Anchor_two_b\nATG\n", ""))
    with pytest.raises(ValueError):
        records = mod.read_family_manifest(path)
        mod.collect_query_anchor_orthologs(store, "", manifest_records=records)


def test_closest_cli_requires_audit_output_and_records_every_original_query(tmp_path):
    store, _path = saved_manifest(tmp_path)
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / "OG1").write_text(">near | Near | species=Near_species\nAA\n"
                                   ">far | Far | species=Far_species\nAA\n")
    args = ["--basis=query_gene", f"--dir_gene_family={store}", f"--dir_query_gene={query_dir}",
            "--query_selection=closest", f"--species_tree={species_tree(tmp_path)}",
            "--target_species=Target_species"]
    for kind in ("columns", "glyphs", "tree", "synteny", "ufboot", "query_map"):
        args.append(f"--out_{kind}={tmp_path / (kind + '.tsv')}")
    with pytest.raises(ValueError, match="--out_selection"):
        mod.run(mod.build_arg_parser().parse_args(args))
    args.append(f"--out_selection={tmp_path / 'selection.tsv'}")
    mod.run(mod.build_arg_parser().parse_args(args))
    with (tmp_path / "selection.tsv").open() as handle:
        audit = list(csv.DictReader(handle, delimiter="\t"))
    assert [r["decision"] for r in audit] == ["retained", "excluded"]
    assert [r["query_id"] for r in audit] == ["near", "far"]
    with (tmp_path / "query_map.tsv").open() as handle:
        mapping = list(csv.DictReader(handle, delimiter="\t"))
    assert [r["query_id"] for r in mapping] == ["near"]
    assert mapping[0]["source_species"] == "Near_species"


def test_manifest_run_writes_selection_overlap_and_complete_long_table(tmp_path):
    store, path = saved_manifest(tmp_path)
    # An explicitly repeated tree block is legal but overlapping genes must be audited.
    path.write_text(path.read_text() + path.read_text().splitlines()[1].replace("Combined", "Second") + "\n")
    tree = tmp_path / "species.nwk"
    tree.write_text("((Anchor_one,Anchor_two),Missing_species);\n")
    args = ["--basis=query_gene", f"--dir_gene_family={store}", f"--family_manifest={path}",
            f"--species_tree={tree}"]
    for kind in ("columns", "glyphs", "tree", "synteny", "ufboot", "query_map", "selection", "overlap", "long"):
        args.append(f"--out_{kind}={tmp_path / (kind + '.tsv')}")
    with pytest.warns(UserWarning, match="multiple plot blocks"):
        mod.run(mod.build_arg_parser().parse_args(args))
    with (tmp_path / "overlap.tsv").open() as handle:
        overlaps = list(csv.DictReader(handle, delimiter="\t"))
    assert len(overlaps) == 4  # Count unique gene identities, not anchor relationships.
    with (tmp_path / "long.tsv").open() as handle:
        long_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(long_rows) == 6
    assert {r["copy_number"] for r in long_rows if r["species"] == "Missing_species"} == {"0"}
    with (tmp_path / "synteny.tsv").open() as handle:
        evidence = list(csv.DictReader(handle, delimiter="\t"))
    assert {r["synteny_status"] for r in evidence} == {"not_evaluable", "anchor_self"}


def test_new_gene_summary_parameters_are_forwarded():
    root = Path(__file__).resolve().parents[1]
    entrypoint = (root / "gg_gene_summary_entrypoint.sh").read_text()
    core = (root / "core/gg_gene_summary_core.sh").read_text()
    registry = (root / "support/gg_entrypoint_config_vars.sh").read_text()
    for name in ("family_manifest", "query_metadata", "query_selection", "target_species",
                 "selection_species", "query_label", "label_map", "focus_species", "legend_columns"):
        parameter = "presence_absence_" + name
        assert parameter in entrypoint and parameter in core and parameter in registry

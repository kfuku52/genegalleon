"""Native LCA/GeneRax bundles retain details, identities, missingness and ZIP reuse."""

import csv
import hashlib
import json
import sqlite3
import subprocess
from pathlib import Path

import pytest
import sqlalchemy
from ete4 import Tree
from nwkit.reconcile import build_reconciliation_table

from workflow.support import species_branch_statistics as stats
from workflow.support.gene_family_output_store import GeneFamilyOutputStore, convert_storage_to_zip
from workflow.support.reconciliation_database import branch_id, write_tables
from workflow.support.reconciliation_summary import build_results
from workflow.tests.test_reconciliation_database import xml


def fixture(tmp_path, *, hidden_loss=False):
    family = "OG0001"
    species = tmp_path / "species.nwk"
    species.write_text("(A:1,(B:1,C:1):1);")
    gene = tmp_path / "gene.nwk"
    gene.write_text("(OG0001_A:1,OG0001_C:1);")
    tree = Tree(gene.read_text(), parser=1)
    table = build_reconciliation_table(
        tree,
        Tree(species.read_text(), parser=1),
        {str(n.name): str(n.name).rsplit("_", 1)[1] for n in tree.leaves()},
        tree_id=family,
    )
    lca = tmp_path / "lca.tsv"
    table.to_csv(lca, sep="\t", index=False)
    root = tmp_path / "families"
    for name in ("generax_xml", "artifact_provenance", "stat_species_branch", "reconciliation_events"):
        (root / name).mkdir(parents=True)
    raw = xml(family)
    if hidden_loss:
        leaf = f'<clade><name>{family}_A</name><eventsRec><leaf speciesLocation="A"/></eventsRec></clade>'
        replacement = (
            '<clade><name>dup_A</name><eventsRec><duplication speciesLocation="A"/></eventsRec>'
            + leaf
            + '<clade><name>lost_A</name><eventsRec><loss speciesLocation="A"/></eventsRec></clade></clade>'
        )
        raw = raw.replace(leaf.encode(), replacement.encode())
    xp = root / "generax_xml" / (family + "_generax.xml")
    xp.write_bytes(raw)
    producer = root / "artifact_provenance" / (family + ".generax.json")
    producer.write_text(
        json.dumps(
            dict(
                family_id=family,
                step="generax",
                parameters=dict(reconciliation_model="UndatedDL"),
                outputs=[dict(label="generax_xml", sha256=hashlib.sha256(raw).hexdigest())],
            )
        )
    )
    return family, root, species, gene, lca, xp, producer


def publish_fixture(tmp_path, hidden_loss=False):
    family, root, species, gene, lca, xp, producer = fixture(tmp_path, hidden_loss=hidden_loss)
    bundle_path = root / "reconciliation_events" / (family + "_reconciliation.json")
    stats_path = root / "stat_species_branch" / (family + "_stat.species_branch.tsv")
    stats.main(
        [
            "--family-id",
            family,
            "--species-tree",
            str(species),
            "--gene-tree",
            str(gene),
            "--lca",
            str(lca),
            "--generax-xml",
            str(xp),
            "--generax-manifest",
            str(producer),
            "--generax-model",
            "UndatedDL",
            "--bundle-out",
            str(bundle_path),
            "--stats-out",
            str(stats_path),
        ]
    )
    metadata = dict(
        family_id=family,
        step="summary_statistics",
        outputs=[
            dict(label=n, sha256=hashlib.sha256(p.read_bytes()).hexdigest())
            for n, p in [("reconciliation_events", bundle_path), ("stat_species_branch", stats_path)]
        ],
    )
    (root / "artifact_provenance" / (family + ".summary_statistics.json")).write_text(json.dumps(metadata))
    return family, root, species, bundle_path, stats_path


def test_details_and_aggregates_preserve_method_difference(tmp_path):
    _, _, _, bp, sp = publish_fixture(tmp_path, hidden_loss=True)
    bundle = json.loads(bp.read_text())
    rows = list(csv.DictReader(sp.open(), delimiter="\t"))
    a = {r["method"]: r for r in rows if r["species_branch_id"] == branch_id(["A"])}
    assert a["generax"]["loss_count"] == "1" and a["lca"]["loss_count"] == "0"
    assert all(r["transfer_in_count"] == "" for r in rows)
    assert all(r["extant_gene_count"] == "1" for r in a.values())
    sources = {s["method"]: s for s in bundle["sources"]}
    assert len([e for e in sources["generax"]["events"] if e[3] == "loss"]) == 2
    assert len([e for e in sources["lca"]["events"] if e[3] == "loss"]) == 1


def test_missing_method_is_null_not_zero(tmp_path):
    family, _, species, gene, lca, _, _ = fixture(tmp_path)
    bundle = stats.build_bundle(family, species, gene, lca)
    missing = next(s for s in bundle["sources"] if s["method"] == "generax")
    assert all(
        r["loss_count"] is None and r["status"] == "not_available"
        for r in stats.branch_rows(missing, bundle["reference"])
    )


def test_lca_producer_binds_original_parameters_versions_and_exact_trees(tmp_path):
    family, _, species, gene, lca, _, _ = fixture(tmp_path)
    producer = tmp_path / "lca-producer.json"
    record = dict(
        family_id=family,
        step="reconciliation",
        parameters=dict(nwkit_identity="qualified-revision"),
        inputs=[
            dict(label=k, sha256=hashlib.sha256(p.read_bytes()).hexdigest())
            for k, p in [("rooted_tree", gene), ("species_tree", species)]
        ],
        outputs=[dict(label="reconciliation", sha256=hashlib.sha256(lca.read_bytes()).hexdigest())],
    )
    producer.write_text(json.dumps(record))
    bundle = stats.build_bundle(family, species, gene, lca, lca_manifest=producer)
    method = next(s for s in bundle["sources"] if s["method"] == "lca")
    assert method["metadata"]["producer"] == record
    gene.write_text(gene.read_text().replace(":1", ":2"))
    with pytest.raises(ValueError, match="LCA input differs"):
        stats.build_bundle(family, species, gene, lca, lca_manifest=producer)


def test_lca_scope_does_not_invent_losses_above_gene_root(tmp_path):
    family, _, species, gene, lca, _, _ = fixture(tmp_path)
    gene.write_text("(OG0001_A1:1,OG0001_A2:1);")
    tree = Tree(gene.read_text(), parser=1)
    build_reconciliation_table(
        tree, Tree(species.read_text(), parser=1), {str(n.name): "A" for n in tree.leaves()}, tree_id=family
    ).to_csv(lca, sep="\t", index=False)
    bundle = stats.build_bundle(family, species, gene, lca)
    source = next(s for s in bundle["sources"] if s["method"] == "lca")
    rows = {r["species_branch_id"]: r for r in stats.branch_rows(source, bundle["reference"])}
    assert rows[branch_id(["A"])]["loss_count"] == 0
    assert rows[branch_id(["B"])]["loss_count"] is None
    assert rows[branch_id(["B"])]["status"] == "outside_gene_root_scope"


@pytest.mark.parametrize("fault", ["locations", "gene_parent", "tree_id", "species_mapping"])
def test_invalid_lca_preserves_existing_outputs(tmp_path, fault):
    family, _, species, gene, lca, _, _ = fixture(tmp_path)
    rows = list(csv.DictReader(lca.open(), delimiter="\t"))
    if fault == "locations":
        rows[0]["implied_loss_species_clade_ids"] = "[]"
    elif fault == "gene_parent":
        rows[-1]["parent_gene_clade_id"] = "unknown"
    elif fault == "tree_id":
        rows[0]["tree_id"] = "foreign"
    else:
        rows[0]["species_event_id"] = "foreign"
    with lca.open("w") as f:
        writer = csv.DictWriter(f, rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    out, table = tmp_path / "old.json", tmp_path / "old.tsv"
    out.write_text("old")
    table.write_text("old")
    with pytest.raises(ValueError):
        stats.main(
            [
                "--family-id",
                family,
                "--species-tree",
                str(species),
                "--gene-tree",
                str(gene),
                "--lca",
                str(lca),
                "--bundle-out",
                str(out),
                "--stats-out",
                str(table),
            ]
        )
    assert out.read_text() == table.read_text() == "old"


@pytest.mark.parametrize("layout", ["raw", "zip"])
def test_db_import_and_downstream_use_common_source_and_keep_legacy_xml(tmp_path, layout):
    family, root, species, _, _ = publish_fixture(tmp_path, hidden_loss=True)
    if layout == "zip":
        convert_storage_to_zip(root, "orthogroup", [family], lambda name: name.split("_")[0])
    dbpath = tmp_path / "families.db"
    engine = sqlalchemy.create_engine("sqlite:///" + str(dbpath))
    with engine.begin() as conn:
        conn.exec_driver_sql("CREATE TABLE tree (orthogroup TEXT)")
        conn.exec_driver_sql("INSERT INTO tree VALUES (?)", (family,))
    metadata = tmp_path / "methods.json"
    metadata.write_text(json.dumps(dict(default=dict(program="GeneRax"), gene_model_contract_sha256="a" * 64)))
    store = GeneFamilyOutputStore(root)
    with store.read_snapshot():
        write_tables(engine, store, species, metadata)
    engine.dispose()
    with sqlite3.connect(dbpath) as db:
        assert db.execute("SELECT COUNT(*) FROM species_reconciliation_source").fetchone()[0] == 2
        assert db.execute("SELECT COUNT(*) FROM reconciliation_event WHERE event_type='loss'").fetchone()[0] == 2
        meta = dict(db.execute("SELECT key,value FROM reconciliation_meta"))
    branches = tmp_path / "branches.json"
    branches.write_text(json.dumps([dict(branch_id="A", descendants=["A"])]))
    fits = tmp_path / "fits.tsv"
    fits.write_text("orthogroup\tstatus\tlogFC\tp_value\nOG0001\tsuccess\t2\t0.02\n")
    config = dict(
        schema="reconciliation_summary_v1",
        family_catalog_sha256=meta["family_catalog_sha256"],
        gene_model_contract_sha256="a" * 64,
        require_complete_reconciliation=True,
        expression_p_value_semantics="one_sided_trap_enrichment_positive_logfc_v1",
        inputs={
            k: dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())
            for k, p in [("branches", branches), ("expression_transcript", fits)]
        },
    )
    for method in ("generax", "lca"):
        config["event_source"] = method
        cfg = tmp_path / (method + ".json")
        cfg.write_text(json.dumps(config))
        output = tmp_path / (method + "_summary")
        build_results(dbpath, cfg, output)
        with sqlite3.connect(output / "reconciliation_summary.db") as db:
            loss, q, extinction = db.execute(
                "SELECT loss_events,transcript_q_branch_loss,lineage_extinction_events FROM branch_family"
            ).fetchone()
        assert (loss, q, extinction) == ((1, 0.02, 1) if method == "generax" else (0, None, None))


def test_statistics_only_profile_disables_inference_and_keeps_selected_tree(tmp_path):
    core = (Path(__file__).parents[1] / "core/gg_gene_evolution_core.sh").read_text()
    block = core[core.index("gene_evolution_plot_only=${") : core.index('case "${gene_evolution_plot_only}" in')]
    script = (
        """set -euo pipefail
mode_gene_evolution=orthogroup;gg_debug_mode=0;gene_evolution_plot_only=0;gene_evolution_statistics_only=1
run_generax=1;run_iqtree=1;run_tree_root=1;run_reconciliation=0;run_summary=0;run_tree_plot=1
gg_print_entrypoint_config_vars() { printf '%s\\n' run_generax run_iqtree run_tree_root run_reconciliation run_summary run_tree_plot; }
"""
        + block
        + """
test "$statistics_reuse_generax" = 1
test "$run_generax$run_iqtree$run_tree_root$run_tree_plot" = 0000
test "$run_reconciliation$run_summary" = 11
"""
    )
    result = subprocess.run(["bash", "-c", script], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr

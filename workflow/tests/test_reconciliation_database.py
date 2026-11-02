import hashlib
import json
import sqlite3

import pytest
import sqlalchemy

from workflow.support import generate_orthogroup_database as generator
from workflow.support.gene_family_output_store import GeneFamilyOutputStore, convert_storage_to_zip
from workflow.support.reconciliation_database import branch_id, digest, parse_family, reference_branches, write_tables
from workflow.support.reconciliation_summary import bh, build_results, sha256


def xml(family="OG0001", loss=True, no_a=False, namespace=True, internal="renumbered"):
    leaves = (f'<clade><name>{family}_B</name><eventsRec><leaf speciesLocation="B"/></eventsRec></clade>'
              if not loss else '<clade><name>loss</name><eventsRec><loss speciesLocation="B"/></eventsRec></clade>')
    a = (f'<clade><name>{family}_A</name><eventsRec><leaf speciesLocation="A"/></eventsRec></clade>'
         if not no_a else '')
    ns = ' xmlns="http://www.recg.org"' if namespace else ''
    return (f'<recPhylo{ns}><spTree><phylogeny><clade><name>root_other</name>'
            f'<clade><name>A</name></clade><clade><name>{internal}</name>'
            '<clade><name>B</name></clade><clade><name>C</name></clade></clade></clade></phylogeny></spTree>'
            '<recGeneTree><phylogeny><clade><name>NULL</name><eventsRec><speciation speciesLocation="root_other"/></eventsRec>'
            + a + f'<clade><name>NULL</name><eventsRec><speciation speciesLocation="{internal}"/></eventsRec>'
            + leaves + f'<clade><name>{family}_C</name><eventsRec><leaf speciesLocation="C"/></eventsRec></clade>'
            '</clade></clade></phylogeny></recGeneTree></recPhylo>').encode()


def build_fixture(tmp_path, layout="raw"):
    root = tmp_path / "families"
    (root / "generax_xml").mkdir(parents=True)
    (root / "artifact_provenance").mkdir()
    for family, lost, absent in (("OG0001", True, False), ("OG0002", True, False), ("OG0003", False, True)):
        data = xml(family, lost, absent)
        if family == "OG0002":
            leaf = '<clade><name>OG0002_A</name><eventsRec><leaf speciesLocation="A"/></eventsRec></clade>'
            duplicated = '<clade><name>dup</name><eventsRec><duplication speciesLocation="A"/></eventsRec>' + leaf + leaf.replace('OG0002_A', 'OG0002_A2') + '</clade>'
            data = data.replace(leaf.encode(), duplicated.encode())
        (root / "generax_xml" / (family + "_generax.xml")).write_bytes(data)
        record = {"family_id": family, "step": "generax", "parameters": {"reconciliation_model": "UndatedDL", "rooting_mode": "mad"},
                  "inputs": [], "outputs": [{"label": "generax_xml", "sha256": hashlib.sha256(data).hexdigest()}],
                  "diagnostics": {"genegalleon_version": "fixture"}}
        (root / "artifact_provenance" / (family + ".generax.json")).write_text(json.dumps(record))
    if layout == "zip":
        convert_storage_to_zip(root, "orthogroup", ["OG0001", "OG0002", "OG0003"], lambda name: name.split("_")[0])
    tree = tmp_path / "tree.nwk"
    tree.write_text("(A,(C,B)original)root;")
    manifest = tmp_path / "methods.json"
    manifest.write_text(json.dumps({"default": {"program": "GeneRax", "MAD": True},
                                    "gene_model_contract_sha256": "a" * 64,
                                    "families": {"OG0002": {"candidate_core": "candidate", "adopted": True, "upstream_convention_validated": False}}}))
    db = tmp_path / "families.db"
    engine = sqlalchemy.create_engine("sqlite:///" + str(db))
    with engine.begin() as conn:
        conn.exec_driver_sql("CREATE TABLE tree (orthogroup TEXT)")
        conn.exec_driver_sql("INSERT INTO tree VALUES (?)", [("OG0001",), ("OG0002",), ("OG0003",), ("OG0004",)])
    store = GeneFamilyOutputStore(root)
    with store.read_snapshot():
        write_tables(engine, store, tree, manifest)
    engine.dispose()
    return db, root, tree, manifest


def tsv(tmp_path, name, columns, rows):
    path = tmp_path / name
    path.write_text("\t".join(columns) + "\n" + "".join("\t".join(map(str, row)) + "\n" for row in rows))
    return path


def config_fixture(tmp_path, db):
    branches = tmp_path / "branches.json"
    branches.write_text(json.dumps([{"branch_id": "loss_B", "descendants": ["B"]},
                                    {"branch_id": "ancestor", "descendants": ["C", "B"]}]))
    expression = tsv(tmp_path, "fits.tsv", ["orthogroup", "status", "logFC", "p_value"],
                     [("OG0001", "success", 2, .02), ("OG0002", "success", 3, .04),
                      ("OG0003", "success", 4, .0001), ("OG0004", "filtered", "NA", "NA")])
    members = tsv(tmp_path, "members.tsv", ["gene_id", "orthogroup"], [("OG0001_A", "OG0001"), ("OG0002_A", "OG0002"), ("OG0002_A2", "OG0002")])
    proteome = tsv(tmp_path, "proteome.tsv", ["gene_id", "orthogroup", "unique_unmodified_peptides_IL_equivalent", "best_peptide_q", "best_protein_q"],
                   [("OG0001_A", "OG0001", 1, .001, .001)])
    shared = tsv(tmp_path, "shared.tsv", ["orthogroup", "peptide", "MS_file", "spectrum", "proteins", "peptide_q", "protein_q"],
                 [("OG0002", "PEPTIDE", "file1", "scan1", "OG0002_A;OG0002_A2", .001, .001)])
    cafe = tsv(tmp_path, "cafe.tsv", ["branch_id", "orthogroup", "copy_change", "significant_decrease"],
               [(branch, family, -1 if family == "OG0001" else 0, int(family == "OG0001"))
                for branch in ("loss_B", "ancestor") for family in ("OG0001", "OG0002", "OG0003", "OG0004")])
    with sqlite3.connect(db) as conn:
        meta = dict(conn.execute("SELECT key,value FROM reconciliation_meta"))
    config = {"schema": "reconciliation_summary_v1", "family_catalog_sha256": meta["family_catalog_sha256"],
              "gene_model_contract_sha256": meta["gene_model_contract_sha256"], "proteome_species": "A",
              "expression_p_value_semantics": "one_sided_trap_enrichment_positive_logfc_v1",
              "inputs": {key: {"path": path.name, "sha256": sha256(path)} for key, path in
                         {"branches": branches, "expression_transcript": expression, "proteome_members": members,
                          "proteome": proteome, "shared_peptides": shared, "cafe": cafe}.items()}}
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    return path, config


@pytest.mark.parametrize("layout", ["raw", "zip"])
def test_db_import_separates_missing_from_zero_and_records_adopted_method(tmp_path, layout):
    db, *_ = build_fixture(tmp_path, layout)
    with sqlite3.connect(db) as conn:
        assert conn.execute("SELECT status,xml_sha256 FROM reconciliation_family WHERE orthogroup='OG0004'").fetchone() == ("missing_xml", None)
        assert conn.execute("SELECT count(*) FROM reconciliation_event WHERE event_type='loss'").fetchone()[0] == 2
        assert conn.execute("SELECT branch_id FROM reconciliation_event WHERE event_type='loss' LIMIT 1").fetchone()[0] == branch_id(["B"])
        assert conn.execute("SELECT count(*) FROM reconciliation_gene").fetchone()[0] == 7
        method = json.loads(conn.execute("SELECT method_json FROM reconciliation_family WHERE orthogroup='OG0002'").fetchone()[0])
        assert method["adopted"] and method["MAD"] and not method["upstream_convention_validated"]
        source_sha, source_json = conn.execute("SELECT generax_manifest_sha256,generax_manifest_json FROM reconciliation_family WHERE orthogroup='OG0002'").fetchone()
        assert len(source_sha) == 64 and json.loads(source_json)["parameters"]["reconciliation_model"] == "UndatedDL"
        assert conn.execute("SELECT count(*) FROM reconciliation_gene WHERE gene_id='loss'").fetchone()[0] == 0


@pytest.mark.parametrize("namespace", [True, False])
def test_internal_renumbering_and_tip_order_do_not_change_mapping(namespace):
    reference, topology = reference_branches(b"(A,(B,C)123)root;")
    assert topology == reference_branches(b"((C,B)xyz,A)abc;")[1]
    clades, events, genes = parse_family(xml(namespace=namespace), "OG0001", reference)
    assert len(clades) == 5 and len(events) == 5 and len(genes) == 2


@pytest.mark.parametrize("change,match", [
    ((b'speciesLocation="B"', b'speciesLocation="unknown"'), "Unmapped"),
    ((b'<loss speciesLocation', b'<unsupported speciesLocation'), "Unsupported"),
    ((b'<name>C</name>', b'<name>D</name>'), "topology/taxa"),
    ((b'<name>OG0001_C</name>', b'<name>OG0001_A</name>'), "duplicate extant"),
    ((b'<name>renumbered</name>', b'<name>A</name>'), "Duplicate XML species label"),
])
def test_invalid_xml_fails_instead_of_inventing_counts(change, match):
    reference, _ = reference_branches(b"(A,(B,C)n)r;")
    with pytest.raises(ValueError, match=match):
        parse_family(xml().replace(*change), "OG0001", reference)


def test_multievent_order_and_unary_gene_chains_are_preserved():
    reference, _ = reference_branches(b"(A,(B,C)n)r;")
    data = xml().replace(b'<leaf speciesLocation="A"/>', b'<transferBack destinationSpecies="A"/><leaf speciesLocation="A"/>')
    _, events, _ = parse_family(data, "OG0001", reference)
    assert [(row[2], row[3]) for row in events if row[1] == 1] == [(0, "transferBack"), (1, "leaf")]


def test_deep_gene_tree_does_not_use_python_recursion():
    reference, _ = reference_branches(b"(A,(B,C)n)r;")
    opening = b'<clade><name>NULL</name><eventsRec><duplication speciesLocation="root_other"/></eventsRec>'
    data = xml().replace(b'<recGeneTree><phylogeny>', b'<recGeneTree><phylogeny>' + opening * 1500)
    data = data.replace(b'</phylogeny></recGeneTree>', b'</clade>' * 1500 + b'</phylogeny></recGeneTree>')
    clades, events, genes = parse_family(data, "OG0001", reference)
    assert len(clades) == len(events) == 1505 and len(genes) == 2


def test_xml_import_is_inside_standard_private_database_builder(tmp_path):
    import pandas as pd
    db, root, tree, method = build_fixture(tmp_path)
    for folder in ("stat_tree", "stat_branch"):
        (root / folder).mkdir()
    for family in ("OG0001", "OG0002", "OG0003", "OG0004"):
        pd.DataFrame([{"num_branch": 2, "num_spe": 1, "num_dup": 0, "num_sp": 2}]).to_csv(root / "stat_tree" / (family + "_stat.tree.tsv"), sep="\t", index=False)
        pd.DataFrame([{"branch_id": 0, "node_name": "n0", "num_sp": 2, "so_event": "S"}]).to_csv(root / "stat_branch" / (family + "_stat.branch.tsv"), sep="\t", index=False)
    output = tmp_path / "full.db"
    generator.main(["--dbpath", str(output), "--dir_gene_family", str(root), "--dir_stat_tree", str(root / "stat_tree"),
                    "--dir_stat_branch", str(root / "stat_branch"), "--reconciliation_species_tree", str(tree),
                    "--reconciliation_method_manifest", str(method)])
    with sqlite3.connect(output) as conn:
        assert conn.execute("SELECT count(*) FROM reconciliation_family").fetchone()[0] == 4
        assert conn.execute("SELECT count(*) FROM tree").fetchone()[0] == 4


def test_branch_local_bh_and_distinct_proteome_evidence_leave_input_db_unchanged(tmp_path):
    db, *_ = build_fixture(tmp_path)
    config, _ = config_fixture(tmp_path, db)
    original = db.read_bytes()
    out = tmp_path / "results"
    build_results(db, config, out)
    assert db.read_bytes() == original
    assert sorted(path.name for path in out.iterdir()) == ["branch_family.tsv", "branch_summary.tsv", "manifest.json", "reconciliation_summary.db"]
    with sqlite3.connect(out / "reconciliation_summary.db") as conn:
        conn.row_factory = sqlite3.Row
        rows = {row["orthogroup"]: dict(row) for row in conn.execute("SELECT * FROM branch_family WHERE branch_id='loss_B'")}
        assert rows["OG0001"]["transcript_q_branch_loss"] == .04
        assert rows["OG0002"]["transcript_q_branch_loss"] == .04
        assert rows["OG0003"]["transcript_q_branch_loss"] is None
        assert rows["OG0004"]["lineage_extinction_events"] is None
        assert rows["OG0001"]["proteome_status"] == "unique_gene_support"
        assert rows["OG0002"]["proteome_status"] == "shared_family_support"
        assert rows["OG0003"]["proteome_status"] == "not_assessable_no_gene"
        assert json.loads(rows["OG0002"]["unique_proteome_genes"]) == []
        assert rows["OG0001"]["cafe_significant_decrease"] == 1
        assert rows["OG0001"]["transcript_trap_enriched"] == 1
        assert rows["OG0001"]["cds_trap_enriched"] is None
    # Re-running a changed downstream view never appends rows to the input DB.
    build_results(db, config, out)
    assert db.read_bytes() == original


@pytest.mark.parametrize("mutation,match", [
    ("family", "catalogue contract"), ("model", "gene-model contract"), ("sha", "SHA mismatch"),
    ("branch", "Unknown or duplicate"), ("member", "membership differs"),
    ("expression", "all families"), ("cafe", "every target branch/family"),
])
def test_bad_evidence_preserves_previous_results(tmp_path, mutation, match):
    db, *_ = build_fixture(tmp_path)
    path, config = config_fixture(tmp_path, db)
    if mutation in ("family", "model"):
        config["family_catalog_sha256" if mutation == "family" else "gene_model_contract_sha256"] = "wrong"
    elif mutation == "sha":
        config["inputs"]["branches"]["sha256"] = "wrong"
    else:
        key = {"branch": "branches", "member": "proteome_members", "expression": "expression_transcript", "cafe": "cafe"}[mutation]
        source = tmp_path / config["inputs"][key]["path"]
        if mutation == "branch":
            source.write_text('[{"branch_id":"unknown","descendants":["A","B"]}]')
        elif mutation == "member":
            source.write_text(source.read_text().replace("OG0001_A", "wrong"))
        else:
            source.write_text("\n".join(source.read_text().splitlines()[:-1]) + "\n")
        config["inputs"][key]["sha256"] = sha256(source)
    path.write_text(json.dumps(config))
    out = tmp_path / "results"
    out.mkdir()
    (out / "previous").write_text("keep")
    with pytest.raises(ValueError, match=match):
        build_results(db, path, out)
    assert (out / "previous").read_text() == "keep"


def test_missing_both_external_evidence_and_xml_remain_unanalysed(tmp_path):
    db, *_ = build_fixture(tmp_path)
    path, config = config_fixture(tmp_path, db)
    config["inputs"] = {"branches": config["inputs"]["branches"]}
    path.write_text(json.dumps(config))
    out = tmp_path / "results"
    build_results(db, path, out)
    with sqlite3.connect(out / "reconciliation_summary.db") as conn:
        assert conn.execute("SELECT proteome_status,shared_family_psms,transcript_status,cafe_copy_change FROM branch_family LIMIT 1").fetchone() == ("not_analysed", None, "not_analysed", None)


def test_assay_proteins_without_a_selected_family_are_counted_separately(tmp_path):
    db, *_ = build_fixture(tmp_path)
    path, config = config_fixture(tmp_path, db)
    protein_file = tmp_path / config["inputs"]["proteome"]["path"]
    with protein_file.open("a") as handle:
        handle.write("outside_gene\t\t1\t0.001\t0.001\n")
    config["inputs"]["proteome"]["sha256"] = sha256(protein_file)
    path.write_text(json.dumps(config))
    out = tmp_path / "results"
    build_results(db, path, out)
    manifest = json.loads((out / "manifest.json").read_bytes())
    assert manifest["evidence_audit"]["identified_proteome_genes_without_selected_family"] == 1


def test_bh_validates_and_handles_empty_scope():
    assert bh({}) == {}
    assert bh({"x": 1, "y": .01}) == {"x": 1, "y": .02}
    with pytest.raises(ValueError):
        bh({"x": -1})


def test_input_identity_fence_detects_reference_changes(tmp_path):
    tree = tmp_path / "tree"
    tree.write_text("(A,B);")
    args = generator.build_parser().parse_args(["--reconciliation_species_tree", str(tree)])
    before = generator.reconciliation_input_identities(args)
    tree.write_text("(A,C);")
    assert generator.reconciliation_input_identities(args) != before


def test_stable_catalogue_digest_uses_exact_strings():
    assert digest(["001", "002"]) != digest([1, 2])
    assert len(hashlib.sha256(b"fixture").hexdigest()) == 64


def test_incomplete_reconciliation_can_be_required_before_publication(tmp_path):
    db, *_ = build_fixture(tmp_path)
    path, config = config_fixture(tmp_path, db)
    config["require_complete_reconciliation"] = True
    path.write_text(json.dumps(config))
    with pytest.raises(ValueError, match="Complete reconciliation"):
        build_results(db, path, tmp_path / "results")
    assert not (tmp_path / "results").exists()


def test_changed_evidence_aborts_before_replacing_results(tmp_path, monkeypatch):
    import workflow.support.reconciliation_summary as summary
    db, *_ = build_fixture(tmp_path)
    path, config = config_fixture(tmp_path, db)
    out = tmp_path / "results"
    out.mkdir()
    (out / "previous").write_text("keep")
    original = summary.collect_evidence
    def mutate(*args):
        result = original(*args)
        (tmp_path / config["inputs"]["expression_transcript"]["path"]).write_text("changed")
        return result
    monkeypatch.setattr(summary, "collect_evidence", mutate)
    with pytest.raises(ValueError, match="input changed"):
        summary.build_results(db, path, out)
    assert (out / "previous").read_text() == "keep"


def test_append_and_out_of_scope_method_overrides_are_rejected(tmp_path):
    db, root, tree, method = build_fixture(tmp_path)
    engine = sqlalchemy.create_engine("sqlite:///" + str(db))
    store = GeneFamilyOutputStore(root)
    with pytest.raises(ValueError, match="append is unsupported"):
        write_tables(engine, store, tree, method)
    engine.dispose()
    payload = json.loads(method.read_bytes())
    payload["families"]["OGoutside"] = {"adopted": True}
    method.write_text(json.dumps(payload))
    other = sqlalchemy.create_engine("sqlite:///" + str(tmp_path / "other.db"))
    with other.begin() as conn:
        conn.exec_driver_sql("CREATE TABLE tree (orthogroup TEXT)")
        conn.exec_driver_sql("INSERT INTO tree VALUES ('OG0001')")
    with pytest.raises(ValueError, match="outside the database"):
        write_tables(other, store, tree, method)
    other.dispose()
    before = db.read_bytes()
    with pytest.raises(SystemExit):
        generator.main(["--dbpath", str(db), "--mode", "append", "--reconciliation_species_tree", str(tree)])
    assert db.read_bytes() == before


def test_xml_is_bound_to_original_generax_manifest(tmp_path):
    _, root, tree, method = build_fixture(tmp_path)
    producer = root / "artifact_provenance/OG0001.generax.json"
    payload = json.loads(producer.read_bytes())
    payload["outputs"][0]["sha256"] = "wrong"
    producer.write_text(json.dumps(payload))
    engine = sqlalchemy.create_engine("sqlite:///" + str(tmp_path / "bad.db"))
    with engine.begin() as conn:
        conn.exec_driver_sql("CREATE TABLE tree (orthogroup TEXT)")
        conn.exec_driver_sql("INSERT INTO tree VALUES (?)", [(family,) for family in ("OG0001", "OG0002")])
    with pytest.raises(ValueError, match="XML differs from its producer"):
        write_tables(engine, GeneFamilyOutputStore(root), tree, method)
    engine.dispose()


def test_audited_pipeline_requires_bound_reference_and_fences_content(tmp_path, monkeypatch):
    import artifact_provenance as provenance
    import gene_family_database_pipeline as pipeline
    root = tmp_path / "output/orthogroup"
    root.mkdir(parents=True)
    global_root = tmp_path / "output/.gg_global_artifacts"
    global_root.mkdir()
    paths = {}
    for subdir, name, body in (
        ("alignment", "OG0001.fa", ">OG0001_A\nACGT\n"),
        ("stat_branch", "OG0001_stat.branch.tsv", "branch_id\tnode_name\tnum_sp\tso_event\n0\tn0\t1\tL\n"),
        ("stat_tree", "OG0001_stat.tree.tsv", "num_branch\tnum_spe\tnum_dup\tnum_sp\n1\t2\t0\t1\n"),
        ("generax_xml", "OG0001_generax.xml", xml().decode()),
    ):
        path = root / subdir / name
        path.parent.mkdir()
        path.write_text(body)
        paths[subdir] = path
    common = ["--workspace-root", str(tmp_path), "--logical-root", str(root)]
    assert provenance.dispatch(["record", *common, "--family-id", "OG0001", "--step", "summary_statistics",
        "--manifest", str(root / "artifact_provenance/OG0001.summary_statistics.json"),
        "--input", "alignment=" + str(paths["alignment"]),
        "--output", "branch=" + str(paths["stat_branch"]), "--output", "tree=" + str(paths["stat_tree"]),
        "--output", "xml=" + str(paths["generax_xml"])]) == 0
    tree = tmp_path / "reference.nwk"
    tree.write_text("(A,(B,C)n)r;")
    output, record_path = tmp_path / "db.sqlite3", tmp_path / "db.json"
    audit = [*common, "--output-tsv", str(tmp_path / "audit.tsv"), "--mode", "orthogroup", "--progress-interval", "0"]
    db_args = ["--overwrite", "1", "--dbpath", str(output), "--dir_gene_family", str(root),
               "--dir_stat_tree", str(root / "stat_tree"), "--dir_stat_branch", str(root / "stat_branch"),
               "--reconciliation_species_tree", str(tree)]
    record = ["--workspace-root", str(tmp_path), "--logical-root", str(global_root), "--family-id", "orthogroup",
              "--step", "gene_family_database", "--manifest", str(record_path),
              "--input-gene-family-store", "families=" + str(root), "--output", "database=" + str(output)]
    with pytest.raises(ValueError, match="missing from provenance"):
        pipeline.run(audit, db_args, record)
    record.extend(["--input", "reference=" + str(tree)])
    assert pipeline.run(audit, db_args, record) == 0
    saved = output.read_bytes(), record_path.read_bytes()
    original = pipeline.database._populate_database
    def mutate(*args):
        result = original(*args)
        tree.write_text("(A,(C,B)different)r;")
        return result
    monkeypatch.setattr(pipeline.database, "_populate_database", mutate)
    with pytest.raises(ValueError, match="changed before publication"):
        pipeline.run(audit, db_args, record)
    assert (output.read_bytes(), record_path.read_bytes()) == saved


def test_downstream_captures_exact_provenance_before_publishing(tmp_path):
    db, *_ = build_fixture(tmp_path)
    config_path, config = config_fixture(tmp_path, db)
    out = tmp_path / "results"
    record_path = tmp_path / "provenance.json"
    argv = ["--manifest", str(record_path), "--workspace-root", str(tmp_path), "--logical-root", str(tmp_path / "global"),
            "--step", "reconciliation_summary", "--family-id", "orthogroup", "--input", "database=" + str(db),
            "--input", "config=" + str(config_path)]
    for name in ("reconciliation_summary.db", "branch_family.tsv", "branch_summary.tsv", "manifest.json"):
        argv += ["--output", name + "=" + str(out / name)]
    with pytest.raises(ValueError, match="exact input/output"):
        build_results(db, config_path, out, argv)
    assert not out.exists() and not record_path.exists()
    for key, entry in config["inputs"].items():
        argv += ["--input", key + "=" + str(tmp_path / entry["path"])]
    build_results(db, config_path, out, argv)
    record = json.loads(record_path.read_bytes())
    assert {row["label"] for row in record["outputs"]} == {path.name for path in out.iterdir()}
    for entry in record["outputs"]:
        assert entry["sha256"] == sha256(out / entry["label"])
    assert next(entry for entry in record["inputs"] if entry["label"] == "database")["sha256"] == sha256(db)

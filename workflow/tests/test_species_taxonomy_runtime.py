"""Taxonomy fixtures stay offline; tree construction/serialization use real NWKIT."""

import csv
import json
import os
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO = Path(__file__).resolve().parents[2]
SUPPORT = REPO / "workflow/support"
sys.path.insert(0, str(SUPPORT))
import species_taxonomy as subject  # noqa: E402

FIXTURE = Path(__file__).with_name("fixtures") / "species_taxonomy"


@pytest.fixture
def args(tmp_path):
    db = tmp_path / "downloads/ete_taxonomy/taxa.sqlite"
    db.parent.mkdir(parents=True)
    with sqlite3.connect(db) as conn:
        conn.execute("CREATE TABLE species(taxid INTEGER PRIMARY KEY, spname TEXT, rank TEXT, track TEXT)")
        conn.execute("CREATE TABLE synonym(taxid INTEGER, spname TEXT)")
        conn.execute("CREATE TABLE merged(taxid_old INTEGER, taxid_new INTEGER)")
        for node in json.loads((FIXTURE / "lineages.json").read_text())["nodes"]:
            conn.execute("INSERT INTO species VALUES(:taxid,:spname,:rank,:track)", node)
        conn.executemany("INSERT INTO synonym VALUES(?,?)", [(3702, "Thale cress"), (3702, "Shared name"), (4530, "Shared name")])
        conn.execute("INSERT INTO merged VALUES(9999999,3702)")
    table = tmp_path / "species.tsv"
    shutil.copyfile(FIXTURE / "species.tsv", table)
    return SimpleNamespace(workspace=tmp_path, species_table=str(table), species_summary="", species_dir=[],
                           taxid_map="", taxonomy_db=str(db), species_tree="auto", ranks=subject.DEFAULT_RANKS, plot_clades=0,
                           output_dir=str(tmp_path / "output/species_taxonomy"))


def table(args, filename="species_taxonomy.tsv"):
    return subject.read_table(Path(args.output_dir) / filename)[1]


def hashes(args):
    return {name: subject.sha256(Path(args.output_dir) / name) for name in subject.FILES}


def test_ncbi_default_ranks_real_tree_and_resume(args, monkeypatch):
    from nwkit.util import read_tree

    subject.run(args)
    output = Path(args.output_dir)
    rows = {row["species"]: row for row in table(args)}
    assert len(rows) == 8
    assert rows["Dionaea_muscipula"]["species_rank"] == "Dionaea muscipula"
    assert rows["Dionaea_muscipula"]["family"] == "Droseraceae"
    assert rows["Dionaea_muscipula"]["domain"] == "Eukaryota"
    assert "superkingdom" not in rows["Dionaea_muscipula"]
    provenance = json.loads((output / "provenance.json").read_text())
    expected_ranks = {node["rank"] for node in json.loads((FIXTURE / "lineages.json").read_text())["nodes"]}
    assert set(provenance["resolved_ranks"]) == expected_ranks
    assert rows["Oryza_sativa"]["subtribe"] == "Oryzinae"
    assert rows["Dionaea_muscipula"]["subtribe"] == ""
    assert rows["Dionaea_muscipula"]["no rank"] == "root"
    assert rows["Dionaea_muscipula"]["cellular root"] == "cellular organisms"
    assert provenance["resolved_ranks"].index("subclass") < provenance["resolved_ranks"].index("order")
    assert provenance["resolved_ranks"].index("subtribe") < provenance["resolved_ranks"].index("genus")
    lineage_rows = table(args, "species_lineage.tsv")
    expected_clades = [row for row in lineage_rows if row["species"] == "Dionaea_muscipula" and row["rank"] == "clade"]
    assert "clade" not in rows["Dionaea_muscipula"]
    assert "clade_taxid" not in rows["Dionaea_muscipula"]
    for clade in expected_clades:
        assert rows["Dionaea_muscipula"][f"clade_{clade['taxid']}"] == clade["name"]
        assert rows["Dionaea_muscipula"][f"clade_{clade['taxid']}_taxid"] == clade["taxid"]
    assert all(row["tree_status"] == "mapped" for row in rows.values())
    assert ":" not in (output / "taxonomy_tree.nwk").read_text()  # No invented evolutionary lengths.
    tree = read_tree(str(output / "taxonomy_tree.nhx"), "auto", True, quiet=True)
    assert set(tree.leaf_names()) == set(rows)
    assert next(leaf for leaf in tree.leaves() if leaf.name == "Dionaea_muscipula").props["gg_family"] == "Droseraceae"
    from urllib.parse import unquote

    tip = next(leaf for leaf in tree.leaves() if leaf.name == "Dionaea_muscipula")
    assert json.loads(unquote(tip.props["gg_clade"])) == [row["name"] for row in expected_clades]
    assert unquote(tip.props[f"gg_clade_{expected_clades[0]['taxid']}"]) == expected_clades[0]["name"]
    mapped = table(args, "taxonomy_mapping.tsv")
    assert {row["taxid"] for row in mapped if row["rank"] == "clade"} >= {row["taxid"] for row in expected_clades}
    assert json.loads(unquote(tree.props["gg_clade"])) == ["Embryophyta", "Tracheophyta", "Euphyllophyta", "Spermatophyta"]
    assert any(row["name"] == "Droseraceae" and row["status"] == "monophyletic"
               for row in table(args, "taxonomy_mapping.tsv"))
    before = hashes(args)
    monkeypatch.setattr(subject, "plot", lambda *a: pytest.fail("Unchanged bundle should be reused"))
    subject.run(args)
    assert hashes(args) == before


def test_species_tree_preserves_lengths_support_and_original_bytes(args):
    from nwkit.util import read_tree

    args.species_tree = str(FIXTURE / "species_tree.nwk")
    subject.run(args)
    output = Path(args.output_dir)
    assert (output / "taxonomy_tree.nwk").read_bytes() == Path(args.species_tree).read_bytes()
    original = read_tree(args.species_tree, "auto", True, quiet=True)
    annotated = read_tree(str(output / "taxonomy_tree.nhx"), "auto", True, quiet=True)
    def edges(tree):
        return {tuple(sorted(node.leaf_names())): (node.dist, node.support, node.name) for node in tree.traverse()}
    assert edges(original) == edges(annotated)
    assert json.loads((output / "provenance.json").read_text())["tree_source"] == "species_tree"


def test_nonmonophyly_missing_and_extra_tree_tips(args):
    tree = args.workspace / "conflict.nwk"
    tree.write_text("((Dionaea_muscipula:1,Arabidopsis_thaliana:1):1,(Drosera_capensis:1,Unknown_species:1):1);")
    args.species_tree = str(tree)
    subject.run(args)
    mappings = table(args, "taxonomy_mapping.tsv")
    group = next(row for row in mappings if row["name"] == "Droseraceae")
    assert group["status"] == "non_monophyletic"
    assert json.loads(group["other_descendants"]) == ["Arabidopsis_thaliana", "Unknown_species"]
    from nwkit.util import read_tree

    annotated = read_tree(str(Path(args.output_dir) / "taxonomy_tree.nhx"), "auto", True, quiet=True)
    assert "gg_family" not in annotated.props
    assert next(leaf for leaf in annotated.leaves() if leaf.name == "Unknown_species").props["gg_taxonomy_status"] == "not_input_species"
    assert sum(row["tree_status"] == "missing_from_tree" for row in table(args)) == 5
    assert any(row["status"] == "missing_from_tree" for row in mappings)


def test_name_ambiguity_qualifiers_merged_ids_and_duplicate_taxids(args):
    Path(args.species_table).write_text(
        "species\ttaxid\nThale_cress\t\nShared_name\t\nArabidopsis_thaliana_cf\t\n"
        "Arabidopsis_thaliana_subsp_unknown\t\nArabidopsis_sp_unknown\t\n"
        "Arabidopsis_thaliana_accession1\t9999999\nArabidopsis_thaliana_accession2\t3702\n"
        "Unknown_species\t777777777\n")
    subject.run(args)
    rows = {row["species"]: row for row in table(args)}
    assert rows["Thale_cress"]["resolution_source"] == "synonym"
    assert rows["Shared_name"]["resolution_status"] == "ambiguous"
    for name in ("Arabidopsis_thaliana_cf", "Arabidopsis_thaliana_subsp_unknown", "Arabidopsis_sp_unknown"):
        assert rows[name]["resolution_status"] == "unresolved"
    assert rows["Unknown_species"]["resolution_status"] == "invalid_taxid"
    assert rows["Unknown_species"]["input_taxid"] == "777777777"
    assert rows["Arabidopsis_thaliana_accession1"]["taxid"] == "3702"
    assert rows["Arabidopsis_thaliana_accession1"]["resolution_source"] == "input_taxid:merged"
    assert sum(row["tree_status"] == "mapped" for row in rows.values()) == 3
    assert all(rows[name]["tree_status"] == "mapped" for name in
               ("Arabidopsis_thaliana_accession1", "Arabidopsis_thaliana_accession2"))


@pytest.mark.parametrize("species", ["Unknown_species", "Arabidopsis_thaliana"])
def test_zero_or_one_resolved_species_has_reviewable_outputs(args, species):
    Path(args.species_table).write_text("species\n" + species + "\n")
    subject.run(args)
    provenance = json.loads((Path(args.output_dir) / "provenance.json").read_text())
    assert provenance["tree_available"] == (species != "Unknown_species")
    assert (Path(args.output_dir) / "taxonomy_tree.png").stat().st_size > 1000
    assert len(table(args)) == 1
    if species == "Unknown_species":
        assert provenance["resolved_ranks"] == []


def test_explicit_ranks_include_repeated_and_whitespace_names(args):
    assert subject.display_value("[Eubacterium] example") == "[Eubacterium] example"
    args.ranks = "subtribe,clade,no rank,superkingdom"
    subject.run(args)
    rows = {row["species"]: row for row in table(args)}
    assert rows["Oryza_sativa"]["subtribe"] == "Oryzinae"
    assert rows["Dionaea_muscipula"]["subtribe"] == ""
    assert rows["Dionaea_muscipula"]["superkingdom"] == ""
    assert sum(bool(value) for key, value in rows["Dionaea_muscipula"].items()
               if key.startswith("clade_") and not key.endswith("_taxid")) > 1
    assert "family" not in rows["Dionaea_muscipula"]


def test_clade_columns_align_taxids_and_interleave_named_ranks(args):
    subject.run(args)
    columns = table(args, "taxonomy_columns.tsv")
    clades = {column["label"]: column for column in columns if column["rank"] == "clade"}
    fixture = json.loads((FIXTURE / "lineages.json").read_text())["nodes"]
    assert {column["taxid"] for column in clades.values()} == {str(node["taxid"]) for node in fixture if node["rank"] == "clade"}
    assert len(clades) == 14
    assert all(json.loads(column["unmet_predecessors"]) == [] for column in columns)
    positions = {column["label"]: int(column["position"]) for column in columns}
    assert positions["Subphylum"] < positions["Embryophyta"] < positions["Tracheophyta"] < positions["Class"]
    assert positions["Class"] < positions["Mesangiospermae"] < positions["rosids"] < positions["malvids"] < positions["Order"]
    assert positions["Family"] < positions["BOP clade"] < positions["Subfamily"]
    rows = {row["species"]: row for row in table(args)}
    for column in clades.values():
        for row in rows.values():
            assert row[column["column"]] in ("", column["label"])
    assert rows["Arabidopsis_thaliana"][clades["rosids"]["column"]] == "rosids"
    assert rows["Oryza_sativa"][clades["rosids"]["column"]] == ""
    assert rows["Amborella_trichopoda"][clades["Mesangiospermae"]["column"]] == ""
    assert rows["Amborella_trichopoda"][clades["Tracheophyta"]["column"]] == "Tracheophyta"


def test_same_named_clades_keep_distinct_columns_and_order_is_input_independent():
    def node(taxid, name, rank="clade"):
        return dict(taxid=taxid, spname=name, rank=rank)
    lineages = {"A": [node(1, "Class A", "class"), node(10, "Ancestor"), node(11, "Same name"), node(21, "Order A", "order")],
                "B": [node(1, "Class A", "class"), node(12, "Same name"), node(22, "Order B", "order")]}
    rows = [dict(species=key, tree_tip="") for key in lineages]
    columns = subject.aligned_columns(rows, lineages, ["class", "clade", "order"], None)
    assert {column["column"] for column in columns if column["label"] == "Same name"} == {"clade_11", "clade_12"}
    assert rows[0]["clade_11"] == "Same name" and rows[0]["clade_12"] == ""
    assert rows[1]["clade_12"] == "Same name" and rows[1]["clade_11"] == ""
    assert subject.aligned_columns(list(reversed(rows)), dict(reversed(list(lineages.items()))),
                                   ["class", "clade", "order"], None) == columns


def test_override_is_authoritative_and_invalidates_cache(args):
    Path(args.species_table).write_text("species\nUnknown_species\n")
    subject.run(args)
    assert table(args)[0]["resolution_status"] == "unresolved"
    mapping = args.workspace / "overrides.tsv"
    mapping.write_text("species\ttaxid\nUnknown_species\t3702\n")
    args.taxid_map = str(mapping)
    subject.run(args)
    assert table(args)[0]["taxid"] == "3702"
    assert table(args)[0]["resolution_source"] == "override"
    before = hashes(args)
    mapping.write_text("species\ttaxid\nWrong_species\t3702\n")
    with pytest.raises(ValueError, match="unknown input species"):
        subject.run(args)
    assert hashes(args) == before


def test_current_inputs_ignore_historical_summary_rows_and_preserve_strain_keys(args):
    args.species_table = ""
    directory = args.workspace / "input/species_cds"
    directory.mkdir(parents=True)
    (directory / "Arabidopsis_thaliana_strainA_assembly.cds.fa").write_text(">gene\nATG\n")
    (directory / "Dionaea_muscipula_assembly.fa.gz").touch()
    (directory / ".Oryza_sativa.fa").touch()
    (directory / "Oryza_sativa.txt").touch()
    summary = args.workspace / "summary.tsv"
    summary.write_text("species_key\ttaxid\nArabidopsis_thaliana_strainA\t3702\nOryza_sativa\t4530\n")
    args.species_summary = str(summary)
    subject.run(args)
    rows = {row["species"]: row for row in table(args)}
    assert set(rows) == {"Arabidopsis_thaliana_strainA", "Dionaea_muscipula"}
    assert rows["Arabidopsis_thaliana_strainA"]["taxid"] == "3702"


def test_auto_tree_discovery_replaces_ncbi_result_and_rejects_broken_tree(args):
    subject.run(args)
    species_dir = args.workspace / "output/species_tree/species_tree_summary"
    species_dir.mkdir(parents=True)
    selected = species_dir / "undated_species_tree.nwk"
    shutil.copyfile(FIXTURE / "species_tree.nwk", selected)
    subject.run(args)
    assert (Path(args.output_dir) / "taxonomy_tree.nwk").read_bytes() == selected.read_bytes()
    before = hashes(args)
    selected.write_text("broken;")
    # A single valid tip is valid Newick but not an overlapping input set.
    selected.write_text("((broken;")
    with pytest.raises(Exception):
        subject.run(args)
    assert hashes(args) == before
    selected.write_text("")
    with pytest.raises(ValueError, match="empty or invalid"):
        subject.run(args)
    assert hashes(args) == before


def test_failure_and_output_corruption_do_not_leave_partial_bundles(args, monkeypatch):
    subject.run(args)
    before = hashes(args)
    args.ranks = "order,family"
    original_plot = subject.plot
    def fail(*unused):
        raise RuntimeError("render failure")
    monkeypatch.setattr(subject, "plot", fail)
    with pytest.raises(RuntimeError, match="render failure"):
        subject.run(args)
    assert hashes(args) == before
    monkeypatch.setattr(subject, "plot", original_plot)
    subject.run(args)
    target = Path(args.output_dir) / "taxonomy_tree.svg"
    target.write_text("corruption")
    subject.run(args)
    assert "<svg" in target.read_text()


@pytest.mark.parametrize("alias", ["symlink", "hardlink", "directory"])
def test_output_aliases_preserve_inputs(args, alias):
    out = Path(args.output_dir)
    out.mkdir(parents=True)
    target = out / "species_taxonomy.tsv"
    original = Path(args.species_table).read_bytes()
    if alias == "symlink":
        target.symlink_to(args.species_table)
    elif alias == "hardlink":
        os.link(args.species_table, target)
    else:
        target.mkdir()
    with pytest.raises(ValueError):
        subject.run(args)
    assert Path(args.species_table).read_bytes() == original
    assert not (out / "provenance.json").exists()


def test_database_changes_invalidate_lineage_and_figures(args):
    subject.run(args)
    with sqlite3.connect(args.taxonomy_db) as conn:
        conn.execute("UPDATE species SET spname='Fixture family renamed' WHERE spname='Droseraceae'")
    subject.run(args)
    assert next(row for row in table(args) if row["species"] == "Dionaea_muscipula")["family"] == "Fixture family renamed"


def test_quoted_tip_aliases_keep_tree_names(args):
    args.species_tree = str(args.workspace / "quoted.nwk")
    Path(args.species_tree).write_text("('Dionaea muscipula':1,'Drosera capensis':1);")
    subject.run(args)
    rows = {row["species"]: row for row in table(args)}
    assert rows["Dionaea_muscipula"]["tree_tip"] == "Dionaea muscipula"
    assert rows["Dionaea_muscipula"]["tree_status"] == "mapped"
    assert (Path(args.output_dir) / "taxonomy_tree.nwk").read_text() == Path(args.species_tree).read_text()


def test_unclassified_descendants_do_not_assert_nonmonophyly(args):
    args.species_tree = str(args.workspace / "uncertain.nwk")
    Path(args.species_tree).write_text("((Dionaea_muscipula:1,Unknown_species:1):1,Drosera_capensis:1);")
    subject.run(args)
    group = next(row for row in table(args, "taxonomy_mapping.tsv") if row["name"] == "Droseraceae")
    assert group["status"] == "unresolved_membership"
    assert json.loads(group["unclassified_descendants"]) == ["Unknown_species"]


def test_late_publication_failure_restores_entire_previous_bundle(args, monkeypatch):
    import nwkit.output_transaction as transactions

    subject.run(args)
    before = hashes(args)
    args.ranks = "order,family"
    original = transactions.os.replace
    failed = False
    def replace(source, destination):
        nonlocal failed
        if str(destination) == str(Path(args.output_dir) / "taxonomy_tree.svg") and not failed:
            failed = True
            raise OSError("simulated publication failure")
        return original(source, destination)
    monkeypatch.setattr(transactions.os, "replace", replace)
    with pytest.raises(OSError, match="publication failure"):
        subject.run(args)
    assert failed
    assert hashes(args) == before


@pytest.mark.parametrize("newick", ["(No_match:1,Also_missing:1);", "(Dionaea_muscipula:1,Dionaea_muscipula:1);",
                                     "('Dionaea muscipula':1,Dionaea_muscipula:1);"])
def test_unusable_or_duplicate_tree_tips_fail_without_replacement(args, newick):
    subject.run(args)
    before = hashes(args)
    args.species_tree = str(args.workspace / "invalid.nwk")
    Path(args.species_tree).write_text(newick)
    with pytest.raises(ValueError):
        subject.run(args)
    assert hashes(args) == before


def test_core_input_taxonomy_skips_array_workers_and_dry_runs(args):
    text = (REPO / "workflow/core/gg_input_generation_core.sh").read_text()
    start = text.index("# Species taxonomy uses the current input set")
    snippet = text[start:text.index("\nfi", start) + 3]
    for mode, dry, download in [("array_prepare", "0", "0"), ("array_worker", "0", "0"),
                                ("single", "1", "0"), ("single", "0", "1")]:
        env = dict(os.environ, run_species_taxonomy="1", input_generation_mode=mode, dry_run=dry, download_only=download)
        # With -u, entering the body would fail because its variables are unset.
        result = subprocess.run(["bash", "-uc", snippet], env=env, text=True, capture_output=True)
        assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("core,mode", [("input_generation", "single"), ("input_generation", "array_finalize"),
                                      ("genome_evolution", "single"), ("gene_summary", "single")])
def test_core_stage_forwards_configuration_to_real_helper(args, core, mode):
    directory = args.workspace / "input/species_cds"
    directory.mkdir(parents=True)
    (directory / "Dionaea_muscipula_sample.fa").write_text(">gene\nATG\n")
    text = (REPO / f"workflow/core/gg_{core}_core.sh").read_text()
    start = text.index("# Species taxonomy uses the current input set")
    snippet = text[start:text.index("\nfi", start) + 3]
    setup = '''set -euo pipefail
ensure_ete_taxonomy_db() { export GG_TAXONOMY_DBFILE="$test_taxonomy_db"; }
effective_species_input_source_dir_path() { echo "$species_cds_dir"; }
'''
    env = dict(os.environ, gg_workspace_dir=str(args.workspace), gg_support_dir=str(SUPPORT),
               run_species_taxonomy="1", taxonomy_species_tree="auto", taxonomy_ranks="order,family",
               taxonomy_plot_clades="1", taxonomy_taxid_map="", input_generation_mode=mode, dry_run="0", download_only="0",
               species_cds_dir=str(directory), species_summary_output=args.species_table,
               test_taxonomy_db=args.taxonomy_db, MPLCONFIGDIR=str(args.workspace / "mpl"))
    result = subprocess.run(["bash", "-c", setup + snippet], env=env, text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr
    with (Path(args.output_dir) / "species_taxonomy.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1 and rows[0]["family"] == "Droseraceae"
    assert "genus" not in rows[0]
    assert json.loads((Path(args.output_dir) / "provenance.json").read_text())["fingerprint"]["plot_clades"] == 1


def test_provider_summary_duplicates_merge_or_report_conflicting_taxids(args):
    args.species_table = ""
    directory = args.workspace / "input/species_cds"
    directory.mkdir(parents=True)
    (directory / "Arabidopsis_thaliana_sample.fa").write_text(">gene\nATG\n")
    summary = args.workspace / "summary.tsv"
    args.species_summary = str(summary)
    summary.write_text("species_key\ttaxid\tprovider\nArabidopsis_thaliana\t3702\tfirst\n"
                       "Arabidopsis_thaliana\t3702\tsecond\nArabidopsis_thaliana\t\tthird\n")
    subject.run(args)
    assert len(table(args)) == 1
    assert table(args)[0]["taxid"] == "3702"
    with summary.open("a") as handle:
        handle.write("Arabidopsis_thaliana\t4530\tfourth\n")
    subject.run(args)
    assert table(args)[0]["resolution_status"] == "conflicting_taxids"
    assert table(args)[0]["taxid"] == ""
    overrides = args.workspace / "overrides.tsv"
    overrides.write_text("species\ttaxid\nArabidopsis_thaliana\t3702\n")
    args.taxid_map = str(overrides)
    subject.run(args)
    assert table(args)[0]["resolution_source"] == "override"
    assert table(args)[0]["taxid"] == "3702"


def test_clade_plot_toggle_preserves_taxonomy_data_and_invalidates_cache(args):
    subject.run(args)
    output = Path(args.output_dir)
    data_files = [name for name in subject.FILES if name.endswith((".tsv", ".nwk", ".nhx"))]
    original = {name: (output / name).read_bytes() for name in data_files}
    assert "clade (" not in (output / "taxonomy_tree.svg").read_text()
    assert "clade_3193" in (output / "species_taxonomy.tsv").read_text()
    for enabled in (1, 0):
        args.plot_clades = enabled
        subject.run(args)
        assert ("clade (" in (output / "taxonomy_tree.svg").read_text()) == bool(enabled)
        assert {name: (output / name).read_bytes() for name in data_files} == original
        assert json.loads((output / "provenance.json").read_text())["fingerprint"]["plot_clades"] == enabled


@pytest.mark.parametrize("payload", ["null", "[]", '{"fingerprint": null}', '{"outputs": []}'])
def test_malformed_provenance_rebuilds_outputs(args, payload):
    subject.run(args)
    path = Path(args.output_dir) / "provenance.json"
    path.write_text(payload)
    subject.run(args)
    assert json.loads(path.read_text())["input_count"] == 8


def test_auto_tree_rejects_broken_symlink(args):
    directory = args.workspace / "output/species_tree"
    directory.mkdir(parents=True)
    (directory / "dated_species_tree.nwk").symlink_to(directory / "missing.nwk")
    with pytest.raises(ValueError, match="empty or invalid"):
        subject.run(args)


@pytest.mark.parametrize("core", ["input_generation", "genome_evolution", "gene_summary"])
def test_taxonomy_paths_resolve_before_working_directory_changes(args, core):
    source = (REPO / f"workflow/core/gg_{core}_core.sh").read_text()
    start = source.index('run_species_taxonomy="')
    end = source.index("\nesac", start) + len("\nesac")
    setup = source[start:end]
    # Execute all top-level directory changes preceding configuration; a late
    # normalization would resolve against this scratch directory instead.
    prior_changes = [line for line in source[:start].splitlines() if line.startswith("cd ")]
    scratch = args.workspace / "scratch"
    scratch.mkdir()
    env = dict(os.environ, taxonomy_species_tree="tree.nwk", taxonomy_taxid_map="map.tsv", dir_tmp=str(scratch))
    command = "set -eu\n" + "\n".join(prior_changes) + "\n" + setup
    command += '\nprintf "%s\\n" "$taxonomy_species_tree" "$taxonomy_taxid_map"'
    result = subprocess.run(["bash", "-c", command], cwd=args.workspace, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == [str(args.workspace / "tree.nwk"), str(args.workspace / "map.tsv")]

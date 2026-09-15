import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from workflow.support import scaffold_taxonomy as scaffold


class Ncbi:
    def get_lineage(self, taxid):
        values = {1: [1], 2: [1, 2], 3: [1, 2, 3], 4: [1, 4], 5: [1, 4, 5]}
        if taxid not in values:
            raise ValueError("Unknown taxid")
        return values[taxid]

    def get_rank(self, ids):
        ranks = {1: "superkingdom", 2: "phylum", 3: "species", 4: "phylum", 5: "species"}
        return {i: ranks[i] for i in ids}


def build(species="Host_species"):
    gff = pd.DataFrame([
        dict(gene_id="x1", chromosome="s1", gff_transcript_id="t1"),
        dict(gene_id="x2", chromosome="s1", gff_transcript_id="t2"),
        dict(gene_id="h", chromosome="s1"),
        dict(gene_id="u", chromosome="s1"),
        dict(gene_id="missing_tax", chromosome="s1"),
        dict(gene_id="only", chromosome="s2"),
        dict(gene_id="trans", chromosome="s1", splice_mode="trans-splicing"),
        dict(gene_id="no_coordinates", chromosome=""),
    ])
    tax = pd.DataFrame({"gene_id": ["x1", "x2", "h", "u", "only"], "lca_taxid": [5, 5, 3, 1, 5]})
    return scaffold.build_tables(gff, tax, species, 3, scaffold.RankResolver(Ncbi()), {"t1": "locus_x", "t2": "locus_x"})


def test_rank_specific_composition_and_isoform_counting():
    genes, summaries = build()
    summary = summaries.loc[summaries.scaffold.eq("s1") & summaries["rank"].eq("phylum")].iloc[0]
    assert summary.total_count == 4  # x locus, h, u, missing
    assert summary.compatible_count == 1
    assert summary.incompatible_count == 1
    assert summary.unresolved_count == 2
    assert summary.cds_id_count == 3
    assert summary.classified_fraction == 0.5
    assert summary.compatible_fraction == 0.5
    assert summary.compatible_all_fraction == 0.25
    assert "trans" not in set(genes.gene_id)
    assert "no_coordinates" not in set(genes.gene_id)
    assert set(genes.loc[genes.gene_id.eq("u") & genes["rank"].eq("phylum"), "label"]) == {"unresolved"}
    assert set(genes.loc[genes["rank"].eq("class"), "label"]) == {"unresolved"}


def test_conflicting_isoforms_are_unresolved():
    gff = pd.DataFrame({"gene_id": ["a", "b"], "chromosome": ["s", "s"], "gff_transcript_id": ["t1", "t2"]})
    tax = pd.DataFrame({"gene_id": ["a", "b"], "lca_taxid": [3, 5]})
    genes, summaries = scaffold.build_tables(gff, tax, "Host_species", 3, scaffold.RankResolver(Ncbi()), {"t1": "L", "t2": "L"})
    assert set(genes.loc[genes["rank"].eq("phylum"), "label"]) == {"unresolved"}
    assert summaries.total_count.eq(1).all()
    with pytest.raises(ValueError, match="Duplicate"):
        scaffold.build_tables(pd.concat([gff, gff.assign(chromosome="other")]), tax, "Host_species", 3, scaffold.RankResolver(Ncbi()))


def test_explicit_gff3_and_gtf_loci(tmp_path):
    path = tmp_path / "input.gff"
    path.write_text("s\tx\tgene\t1\t9\t.\t+\t.\tID=g\n"
                    "s\tx\tmRNA\t1\t9\t.\t+\t.\tID=t;Parent=g\n"
                    "s\tx\tCDS\t1\t9\t.\t+\t0\tID=c;Parent=t\n"
                    's\tx\tCDS\t20\t29\t.\t+\t0\tgene_id "g2"; transcript_id "t2";\n')
    assert scaffold.gff_loci(path) == {"g": "g", "t": "g", "c": "g", "g2": "g2", "t2": "g2"}


def fixture_context(tmp_path):
    data, _ = build()
    other, _ = build("Donor_species")  # deliberately reuse the same gene/scaffold IDs
    pd.concat([data, other]).to_csv(tmp_path / "all_gene_taxonomy.tsv", sep="\t", index=False)
    tree = tmp_path / "species.nwk"
    tree.write_text("((Host_species:1)42:1,Donor_species:1)root;")
    genes = pd.DataFrame([
        dict(orthogroup="OG1", gene_id="x1", gene_taxon="Host species"),
        dict(orthogroup="OG1", gene_id="u", gene_taxon="Host species"),
        dict(orthogroup="OG1", gene_id="h", gene_taxon="Donor species"),
        dict(orthogroup="OG2", gene_id="only", gene_taxon="Host species"),
        dict(orthogroup="OG3", gene_id="absent", gene_taxon="Host species"),
    ])
    branches = pd.DataFrame([
        dict(orthogroup="OG1", generax_transfer="Y@Donor_species@42", candidate_genes="x1; u; h"),
        dict(orthogroup="OG2", generax_transfer="Y@Donor_species@Host_species", candidate_genes="only"),
        dict(orthogroup="OG3", generax_transfer="Y", candidate_genes="absent"),
    ])
    return branches, genes, tree


def test_recipient_only_unique_scaffolds_and_background(tmp_path):
    branches, genes, tree = fixture_context(tmp_path)
    branches, genes = scaffold.attach_context(branches, genes, tmp_path, tree)
    b = branches.iloc[0]
    assert b.host_scaffold_status == "measured"
    assert b.host_scaffold_recipient_gene_count == 2
    assert b.host_scaffold_count == 1
    assert b.host_scaffold_phylum_total_count == 4
    assert b.host_scaffold_background_phylum_total_count == 2  # h + missing; x2 excluded with x1
    assert b.host_scaffold_background_phylum_compatible_fraction == 1
    assert b.host_scaffold_background_phylum_classified_fraction == 0.5
    assert branches.iloc[1].host_scaffold_background_phylum_total_count == 0
    assert pd.isna(branches.iloc[1].host_scaffold_background_phylum_compatible_fraction)
    assert branches.iloc[2].host_scaffold_status == "recipient_unresolved"
    assert genes.iloc[-1].host_scaffold_status == "gene_not_mapped"
    assert genes.iloc[0].host_scaffold_count_unit == "gff_locus"


def test_cross_orthogroup_exclusion_and_partial_mapping(tmp_path):
    branches, genes, tree = fixture_context(tmp_path)
    # h is a candidate in a different OG: remove it from OG1's background too.
    genes.loc[len(genes)] = ["OG4", "h", "Host species"]
    genes.loc[len(genes)] = ["OG1", "absent", "Host species"]
    branches.loc[0, "candidate_genes"] += "; absent"
    branches, genes = scaffold.attach_context(branches, genes, tmp_path, tree)
    b = branches.iloc[0]
    assert b.host_scaffold_status == "partial"
    assert b.host_scaffold_background_phylum_total_count == 1
    assert b.host_scaffold_background_phylum_unresolved_count == 1
    assert pd.isna(b.host_scaffold_background_phylum_compatible_fraction)


def test_missing_input_is_not_zero_support(tmp_path):
    branches, genes, tree = fixture_context(tmp_path)
    branches, genes = scaffold.attach_context(branches, genes, "", tree)
    assert branches.host_scaffold_status.eq("missing_scaffold_taxonomy").all()
    assert genes.host_scaffold_phylum_compatible_fraction.isna().all()


def test_per_species_files_and_duplicate_rejection(tmp_path):
    branches, genes, tree = fixture_context(tmp_path)
    combined = tmp_path / "all_gene_taxonomy.tsv"
    data = pd.read_csv(combined, sep="\t")
    combined.unlink()
    for species, group in data.groupby("species"):
        group.to_csv(tmp_path / f"{species}_gene_taxonomy.tsv", sep="\t", index=False)
    measured, _ = scaffold.attach_context(branches, genes, tmp_path, tree)
    assert measured.iloc[0].host_scaffold_phylum_total_count == 4
    data.to_csv(combined, sep="\t", index=False)
    with pytest.raises(ValueError, match="Duplicate species"):
        scaffold.attach_context(branches, genes, tmp_path, tree)


def test_cli_raw_taxonomy_to_hgt_scorer(tmp_path):
    # A small real ETE-compatible SQLite taxonomy DB; no network/stub resolver.
    import sqlite3

    db = tmp_path / "taxonomy.sqlite"
    with sqlite3.connect(db) as connection:
        connection.execute("CREATE TABLE stats (version INT)")
        from ete4.ncbi_taxonomy.ncbiquery import DB_VERSION
        connection.execute("INSERT INTO stats VALUES (?)", (DB_VERSION,))
        connection.execute("CREATE TABLE species (taxid INT PRIMARY KEY, parent INT, spname TEXT COLLATE NOCASE, common TEXT, rank TEXT, track TEXT)")
        connection.execute("CREATE TABLE synonym (taxid INT, spname TEXT COLLATE NOCASE)")
        connection.execute("CREATE TABLE merged (taxid_old INT, taxid_new INT)")
        connection.executemany("INSERT INTO species VALUES (?, ?, ?, '', ?, ?)", [
            (1, 1, "root", "no rank", "1"),
            (2, 1, "Eukaryota", "superkingdom", "2,1"),
            (3, 2, "Arthropoda", "phylum", "3,2,1"),
            (4, 3, "Host species", "species", "4,3,2,1"),
        ])
    gff = tmp_path / "gff.tsv"
    pd.DataFrame({"gene_id": ["x", "h"], "chromosome": ["s", "s"]}).to_csv(gff, sep="\t", index=False)
    tax = tmp_path / "tax.tsv"
    tax.write_text("x\t2\tsuperkingdom\tEukaryota\t1\t1\t1\t1\t1;2\n"
                   "h\t4\tspecies\tHost species\t1\t1\t1\t1\t1;2;3;4\n")
    subprocess.run([sys.executable, str(Path(scaffold.__file__)), "--gff-info", str(gff), "--taxonomy", str(tax),
                    "--species", "Host_species", "--taxonomy-dbfile", str(db), "--gene-out", str(tmp_path / "host_gene_taxonomy.tsv"),
                    "--scaffold-out", str(tmp_path / "host_scaffold_taxonomy.tsv")], check=True)
    out = pd.read_csv(tmp_path / "host_scaffold_taxonomy.tsv", sep="\t")
    assert out.loc[out["rank"].eq("phylum"), "classified_fraction"].iloc[0] == 0.5
    branch_db = tmp_path / "branches.sqlite"
    with sqlite3.connect(branch_db) as connection:
        pd.DataFrame([
            dict(orthogroup="OG1", branch_id=2, node_name="n", gene_labels="x", num_leaf=1, so_event="S", taxon="", generax_event="H", generax_transfer="Y@Donor_species@Host_species"),
            dict(orthogroup="OG1", branch_id=1, node_name="x", gene_labels="x", num_leaf=1, so_event="L", taxon="Host species", generax_event="L", generax_transfer=""),
        ]).to_sql("branch", connection, index=False)
    tree = tmp_path / "species.nwk"
    tree.write_text("(Host_species,Donor_species)root;")
    outputs = {name: tmp_path / f"{name}.tsv" for name in ("branch", "gene", "orthogroup")}
    cmd = [sys.executable, str(Path(scaffold.__file__).with_name("score_hgt_candidates.py")), "--dbpath", str(branch_db),
           "--dir_scaffold_taxonomy", str(tmp_path), "--species_tree", str(tree), "--taxonomy_dbfile", ""]
    for name, path in outputs.items():
        cmd += [f"--{name}_out", str(path)]
    subprocess.run(cmd, check=True)
    result = pd.read_csv(outputs["branch"], sep="\t").iloc[0]
    assert result.host_scaffold_background_phylum_compatible_fraction == 1
    assert result.host_scaffold_background_phylum_total_count == 1

    # Execute the actual genome-annotation stage with pre-existing raw inputs,
    # then verify provenance reuse. No taxonomy search/download is permitted.
    from workflow.tests.test_hgt_end_to_end import _install_fake_conda
    workspace = tmp_path / "workspace"
    cds_dir = workspace / "input/species_cds"
    cds_dir.mkdir(parents=True)
    (cds_dir / "Host_species.fa").write_text(
        ">Host_species_x\nATGAAATAA\n>Host_species_h\nATGCCCTAA\n")
    raw_gff_dir = workspace / "input/species_gff"
    raw_gff_dir.mkdir()
    (raw_gff_dir / "Host_species.gff").write_text(
        "s\ttest\tgene\t1\t9\t.\t+\t.\tID=Host_species_x\n"
        "s\ttest\tgene\t20\t28\t.\t+\t.\tID=Host_species_h\n")
    info_dir = workspace / "output/species_gff_info"
    info_dir.mkdir(parents=True)
    pd.DataFrame({"gene_id": ["Host_species_x", "Host_species_h"], "chromosome": ["s", "s"]}).to_csv(
        info_dir / "Host_species_gff_info.tsv", sep="\t", index=False)
    tax_dir = workspace / "output/species_cds_mmseqs2taxonomy"
    tax_dir.mkdir()
    (tax_dir / "Host_species_mmseqs2taxonomy.tsv").write_text(
        tax.read_text().replace("x\t", "Host_species_x\t").replace("h\t", "Host_species_h\t"))
    db_dir = workspace / "downloads/ete_taxonomy"
    db_dir.mkdir(parents=True)
    shutil.copy2(db, db_dir / "taxa.sqlite")
    root = Path(__file__).resolve().parents[2]
    core = root / "workflow/core/gg_genome_annotation_core.sh"
    env = {key: value for key, value in os.environ.items() if not key.startswith(("CONDA", "MAMBA"))}
    env.update({flag: "0" for flag in set(re.findall(r"\brun_[a-z0-9_]+\b", core.read_text()))})
    env.update(gg_workspace_dir=str(workspace), GG_ARRAY_TASK_ID="1", run_scaffold_taxonomy="1", delete_tmp_dir="0",
               gg_support_dir=str(root / "workflow/support"))
    env["PATH"] = f"{_install_fake_conda(tmp_path)}{os.pathsep}{env['PATH']}"
    scaffold_output = workspace / "output/species_scaffold_taxonomy/Host_species_scaffold_taxonomy.tsv"
    mtimes = []
    for _ in range(2):
        completed = subprocess.run(["bash", str(core)], env=env, capture_output=True, text=True)
        assert completed.returncode == 0, completed.stdout + completed.stderr
        mtimes.append(scaffold_output.stat().st_mtime_ns)
    assert mtimes[0] == mtimes[1]
    assert pd.read_csv(scaffold_output, sep="\t").total_count.eq(2).all()
    assert (workspace / "output/artifact_provenance/genome_annotation/Host_species.scaffold_taxonomy.json").is_file()
    raw_taxonomy = tax_dir / "Host_species_mmseqs2taxonomy.tsv"
    raw_taxonomy.write_text(raw_taxonomy.read_text().replace("Host_species_h\t4\t", "Host_species_h\t2\t"))
    completed = subprocess.run(["bash", str(core)], env=env, capture_output=True, text=True)
    assert completed.returncode == 3  # Existing stop-on-stale policy applies.
    assert "Stale artifact detected" in completed.stdout + completed.stderr
    assert scaffold_output.stat().st_mtime_ns == mtimes[0]


def test_workflow_wiring_and_provenance():
    root = Path(__file__).resolve().parents[2]
    core = (root / "workflow/core/gg_genome_annotation_core.sh").read_text()
    assert core.index('task="Host scaffold taxonomy composition"') < core.index('task="Contaminated sequence removal from the CDS sequences"')
    assert '--input "cds_taxonomy=${file_sp_cds_mmseqs2taxonomy}"' in core
    assert '--input "aggregator=${gg_support_dir}/scaffold_taxonomy.py"' in core
    hgt = (root / "workflow/core/gg_hgt_core.sh").read_text()
    assert '--dir_scaffold_taxonomy "${gg_workspace_output_dir}/species_scaffold_taxonomy"' in hgt
    assert 'hgt_eval_provenance_args "species_tree" "${hgt_species_tree_path}"' in hgt

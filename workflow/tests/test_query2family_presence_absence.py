from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace

import pandas

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"


def load_module(name: str):
    path = SUPPORT_DIR / name
    spec = spec_from_file_location(path.stem, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_stat_branch(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    pandas.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_query2family_presence_absence_counts_leaf_species(tmp_path: Path):
    mod = load_module("gene_family_presence_absence.py")
    query_gene = tmp_path / "input" / "query_gene"
    query2family = tmp_path / "output" / "query2family"
    query_gene.mkdir(parents=True)
    (query_gene / "A").write_text("geneA\n", encoding="utf-8")
    (query_gene / "B").write_text("geneB\n", encoding="utf-8")
    (query_gene / ".hidden").write_text("ignored\n", encoding="utf-8")
    tree = tmp_path / "species.nwk"
    tree.write_text("((Species_one:1,Species_two:1):1,Species_three:1);", encoding="utf-8")

    write_stat_branch(
        query2family / "stat_branch" / "A_stat.branch.tsv",
        [
            {"branch_id": 0, "node_name": "Species_one_g1", "spnode_coverage": "Species_one", "so_event": "L"},
            {"branch_id": 1, "node_name": "Species_one_g2", "spnode_coverage": "Species_one", "so_event": "L"},
            {"branch_id": 2, "node_name": "Species_three_g1", "spnode_coverage": "Species_three", "so_event": "L"},
            {"branch_id": 3, "node_name": "internal", "spnode_coverage": "Species_one", "so_event": "D"},
        ],
    )
    write_stat_branch(
        query2family / "stat_branch" / "B_stat.branch.tsv",
        [
            {"branch_id": 0, "node_name": "Species_two_g1", "spnode_coverage": "Species_two", "so_event": "L"},
        ],
    )

    out_presence = tmp_path / "presence.tsv"
    out_copy_number = tmp_path / "copy.tsv"
    out_long = tmp_path / "long.tsv"
    mod.run(
        SimpleNamespace(
            dir_query2family=str(query2family),
            mode="query2family",
            dir_gene_family=str(query2family),
            dir_query_gene=str(query_gene),
            orthogroup_genecount="",
            species_tree=str(tree),
            out_presence=str(out_presence),
            out_copy_number=str(out_copy_number),
            out_long=str(out_long),
            out_plot_presence=str(tmp_path / "plot.tsv"),
            out_plot_copy_number=str(tmp_path / "copy.plot.tsv"),
            out_plot_long=str(tmp_path / "plot.long.tsv"),
            out_selection=str(tmp_path / "selection.tsv"),
            include_incomplete=0,
            max_families="auto",
            family_ids="",
            family_file="",
        )
    )

    presence = pandas.read_csv(out_presence, sep="\t", index_col=0)
    copy_number = pandas.read_csv(out_copy_number, sep="\t", index_col=0)
    assert presence.index.tolist() == ["Species_one", "Species_two", "Species_three"]
    assert presence.loc["Species_one", "A"] == 1
    assert presence.loc["Species_two", "A"] == 0
    assert presence.loc["Species_two", "B"] == 1
    assert copy_number.loc["Species_one", "A"] == 2
    assert copy_number.loc["Species_three", "A"] == 1

    long_df = pandas.read_csv(out_long, sep="\t")
    assert set(long_df["query"]) == {"A", "B"}
    assert set(long_df["status"]) == {"complete"}


def test_query2family_presence_absence_can_include_incomplete_queries(tmp_path: Path):
    mod = load_module("gene_family_presence_absence.py")
    query_gene = tmp_path / "input" / "query_gene"
    query2family = tmp_path / "output" / "query2family"
    query_gene.mkdir(parents=True)
    (query_gene / "A").write_text("geneA\n", encoding="utf-8")
    (query_gene / "C").write_text("geneC\n", encoding="utf-8")
    tree = tmp_path / "species.nwk"
    tree.write_text("(Species_one:1,Species_two:1);", encoding="utf-8")
    write_stat_branch(
        query2family / "stat_branch" / "A_stat.branch.tsv",
        [
            {"branch_id": 0, "node_name": "Species_one_g1", "spnode_coverage": "Species_one", "num_leaf": "1"},
        ],
    )

    out_presence = tmp_path / "presence.tsv"
    out_copy_number = tmp_path / "copy.tsv"
    out_long = tmp_path / "long.tsv"
    mod.run(
        SimpleNamespace(
            dir_query2family=str(query2family),
            mode="query2family",
            dir_gene_family=str(query2family),
            dir_query_gene=str(query_gene),
            orthogroup_genecount="",
            species_tree=str(tree),
            out_presence=str(out_presence),
            out_copy_number=str(out_copy_number),
            out_long=str(out_long),
            out_plot_presence=str(tmp_path / "plot.tsv"),
            out_plot_copy_number=str(tmp_path / "copy.plot.tsv"),
            out_plot_long=str(tmp_path / "plot.long.tsv"),
            out_selection=str(tmp_path / "selection.tsv"),
            include_incomplete=1,
            max_families="auto",
            family_ids="",
            family_file="",
        )
    )

    presence = pandas.read_csv(out_presence, sep="\t", index_col=0)
    assert presence.columns.tolist() == ["A", "C"]
    assert pandas.isna(presence.loc["Species_one", "C"])
    long_df = pandas.read_csv(out_long, sep="\t")
    assert "missing_stat_branch" in set(long_df["status"])


def test_gene_family_presence_absence_query_subset_keeps_requested_order(tmp_path: Path):
    mod = load_module("gene_family_presence_absence.py")
    query_gene = tmp_path / "input" / "query_gene"
    query2family = tmp_path / "output" / "query2family"
    query_gene.mkdir(parents=True)
    for query_id in ["A", "B", "C"]:
        (query_gene / query_id).write_text(f"{query_id}\n", encoding="utf-8")
        write_stat_branch(
            query2family / "stat_branch" / f"{query_id}_stat.branch.tsv",
            [
                {"branch_id": 0, "node_name": f"{query_id}_gene", "spnode_coverage": "Species_one", "so_event": "L"},
            ],
        )
    tree = tmp_path / "species.nwk"
    tree.write_text("(Species_one:1,Species_two:1);", encoding="utf-8")

    out_dir = tmp_path / "out"
    mod.run(
        SimpleNamespace(
            mode="query2family",
            dir_gene_family=str(query2family),
            dir_query_gene=str(query_gene),
            orthogroup_genecount="",
            species_tree=str(tree),
            out_presence=str(out_dir / "presence.tsv"),
            out_copy_number=str(out_dir / "copy.tsv"),
            out_long=str(out_dir / "long.tsv"),
            out_plot_presence=str(out_dir / "plot.tsv"),
            out_plot_copy_number=str(out_dir / "copy.plot.tsv"),
            out_plot_long=str(out_dir / "plot.long.tsv"),
            out_selection=str(out_dir / "selection.tsv"),
            include_incomplete=0,
            max_families="auto",
            family_ids="C,A",
            family_file="",
        )
    )

    selection = pandas.read_csv(out_dir / "selection.tsv", sep="\t")
    plot_presence = pandas.read_csv(out_dir / "plot.tsv", sep="\t", index_col=0)
    assert selection["family_id"].tolist() == ["C", "A"]
    assert plot_presence.columns.tolist() == ["C", "A"]


def test_gene_family_presence_absence_orthogroup_auto_caps_plot_to_100(tmp_path: Path):
    mod = load_module("gene_family_presence_absence.py")
    orthogroup = tmp_path / "output" / "orthogroup"
    genecount = tmp_path / "Orthogroups.GeneCount.selected.tsv"
    tree = tmp_path / "species.nwk"
    tree.write_text("(Sp1:1,Sp2:1);", encoding="utf-8")
    rows = []
    for i in range(1, 106):
        og = f"OG{i:07d}"
        rows.append({"Orthogroup": og, "Sp1": 1, "Sp2": 1, "Total": 2})
        write_stat_branch(
            orthogroup / "stat_branch" / f"{og}_stat.branch.tsv",
            [
                {"branch_id": 0, "node_name": f"{og}_sp1", "spnode_coverage": "Sp1", "so_event": "L"},
                {"branch_id": 1, "node_name": f"{og}_sp2", "spnode_coverage": "Sp2", "so_event": "L"},
            ],
        )
    pandas.DataFrame(rows).to_csv(genecount, sep="\t", index=False)

    out_dir = tmp_path / "out"
    mod.run(
        SimpleNamespace(
            mode="orthogroup",
            dir_gene_family=str(orthogroup),
            dir_query_gene="",
            orthogroup_genecount=str(genecount),
            species_tree=str(tree),
            out_presence=str(out_dir / "presence.tsv"),
            out_copy_number=str(out_dir / "copy.tsv"),
            out_long=str(out_dir / "long.tsv"),
            out_plot_presence=str(out_dir / "plot.tsv"),
            out_plot_copy_number=str(out_dir / "copy.plot.tsv"),
            out_plot_long=str(out_dir / "plot.long.tsv"),
            out_selection=str(out_dir / "selection.tsv"),
            include_incomplete=0,
            max_families="auto",
            family_ids="",
            family_file="",
        )
    )

    selection = pandas.read_csv(out_dir / "selection.tsv", sep="\t")
    presence = pandas.read_csv(out_dir / "presence.tsv", sep="\t", index_col=0)
    plot_presence = pandas.read_csv(out_dir / "plot.tsv", sep="\t", index_col=0)
    assert presence.shape[1] == 105
    assert plot_presence.shape[1] == 100
    assert selection["family_id"].iloc[0] == "OG0000001"
    assert selection["family_id"].iloc[-1] == "OG0000100"


def test_gene_family_database_id_strips_genegalleon_suffixes():
    mod = load_module("generate_orthogroup_database.py")

    assert mod.gene_family_id_from_path("/tmp/HOG000001_stat.branch.tsv") == "HOG000001"
    assert mod.gene_family_id_from_path("/tmp/AC_stat.tree.tsv") == "AC"
    assert mod.gene_family_id_from_path("/tmp/A.B.stat.branch.tsv") == "A.B"
    assert mod.gene_family_id_from_path("/tmp/A.B.stat.tree.tsv") == "A.B"
    assert mod.gene_family_id_from_path("/tmp/HOG000001.csubst_cb_3.tsv") == "HOG000001"
    assert mod.gene_family_id_from_path("/tmp/HOG000001_csubst_cb_2.tsv") == "HOG000001"

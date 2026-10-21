import subprocess
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace

import numpy
import pandas
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "orthogroup_selection.py"


def load_module():
    spec = spec_from_file_location("orthogroup_selection", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_get_df_gc_original_preserves_species_and_gene_order_for_ties():
    mod = load_module()
    args = SimpleNamespace(gene_size_quantiles="0.5")

    df_og_original = pandas.DataFrame(
        [
            {"Orthogroup": "OG1", "spA": "g2,g1", "spB": "g3"},
            {"Orthogroup": "OG2", "spA": numpy.nan, "spB": "g4,gx"},
        ]
    )
    df_gc_original = pandas.DataFrame(
        [
            {"Orthogroup": "OG1", "Total": 3},
            {"Orthogroup": "OG2", "Total": 2},
        ]
    )
    fx2tab = pandas.DataFrame(
        [
            {"#id": "g1", "length": 100},
            {"#id": "g2", "length": 100},
            {"#id": "g3", "length": 300},
            {"#id": "g4", "length": 400},
        ]
    )

    out = mod.get_df_gc_original(df_og_original, df_gc_original, fx2tab, args)
    assert out.loc[out["Orthogroup"] == "OG1", "geneid_0.5"].iloc[0] == "g2"
    assert out.loc[out["Orthogroup"] == "OG2", "geneid_0.5"].iloc[0] == "g4"


def test_get_df_gc_original_assigns_multiple_quantiles():
    mod = load_module()
    args = SimpleNamespace(gene_size_quantiles="0.25,0.75")

    df_og_original = pandas.DataFrame(
        [
            {"Orthogroup": "OG1", "spA": "g1,g2", "spB": "g3,g4"},
        ]
    )
    df_gc_original = pandas.DataFrame(
        [
            {"Orthogroup": "OG1", "Total": 4},
        ]
    )
    fx2tab = pandas.DataFrame(
        [
            {"#id": "g1", "length": 1},
            {"#id": "g2", "length": 2},
            {"#id": "g3", "length": 3},
            {"#id": "g4", "length": 4},
        ]
    )

    out = mod.get_df_gc_original(df_og_original, df_gc_original, fx2tab, args)
    assert out["geneid_0.25"].iloc[0] == "g2"
    assert out["geneid_0.75"].iloc[0] == "g3"


def test_get_species_protein_files_filters_and_sorts_fasta_paths(tmp_path):
    mod = load_module()
    # Deliberately create in reverse order so filesystem order is not the oracle.
    names = ["sp5.fna.gz", "sp4.fas", "sp3.fa.gz", "sp2.fasta.gz", "sp1.fa"]
    for name in names:
        (tmp_path / name).touch()
    (tmp_path / "README.txt").touch()
    (tmp_path / ".hidden.fa").touch()
    (tmp_path / "directory.fasta.gz").mkdir()

    assert mod.get_species_protein_files(str(tmp_path)) == sorted(names)


def test_format_command_for_log_escapes_control_characters():
    mod = load_module()
    command_log = mod._format_command_for_log(["mmseqs", "search", "clean", "bad\npath"])

    assert "\\n" in command_log
    assert "bad\npath" not in command_log


def test_normalize_mmseqs_split_memory_limit_rejects_control_characters():
    mod = load_module()

    with pytest.raises(ValueError, match="contains control characters"):
        mod.normalize_mmseqs_split_memory_limit("8G\n--bad")


def test_get_concatenated_fx2tab_runs_single_seqkit_call(monkeypatch, tmp_path):
    mod = load_module()
    file_a = tmp_path / "spA.fa"
    file_b = tmp_path / "spB.fa.gz"
    file_a.write_text(">a\nAT\n", encoding="utf-8")
    file_b.write_text("", encoding="utf-8")
    monkeypatch.setattr(mod, "get_species_protein_files", lambda _: [file_a.name, file_b.name])

    calls = []

    def fake_run_command(command, stdout_file=None, append=False, tool_name="tool", capture_output=False):
        calls.append(command)
        return b"#id\tlength\nA1\t100\nB1\t200\n"

    monkeypatch.setattr(mod, "run_command", fake_run_command)
    args = SimpleNamespace(dir_species_protein=str(tmp_path), ncpu=3)
    out = mod.get_concatenated_fx2tab(args)

    assert len(calls) == 1
    command = calls[0]
    assert command[:8] == ["seqkit", "fx2tab", "--threads", "3", "--length", "--name", "--header-line", "--only-id"]
    assert command[8:] == [str(file_a), str(file_b)]
    assert out.to_dict("records") == [{"#id": "A1", "length": 100}, {"#id": "B1", "length": 200}]


def test_get_concatenated_fx2tab_returns_empty_when_no_input(monkeypatch, tmp_path):
    mod = load_module()
    monkeypatch.setattr(mod, "get_species_protein_files", lambda _: [])
    args = SimpleNamespace(dir_species_protein=str(tmp_path), ncpu=1)
    out = mod.get_concatenated_fx2tab(args)
    assert out.empty
    assert list(out.columns) == ["#id", "length"]


@pytest.mark.parametrize("gene_ids", [["001", "1"], ["NA", "NULL", "nan", "001", "1"]])
def test_fx2tab_retains_literal_protein_ids_through_quantile_selection(monkeypatch, tmp_path, gene_ids, capsys):
    mod = load_module()
    (tmp_path / "Genus_species.fa").touch()
    raw = "#id\tlength\n" + "".join(f"{gene}\t{100 + i}\n" for i, gene in enumerate(gene_ids))
    monkeypatch.setattr(mod, "run_command", lambda *args, **kwargs: raw.encode())
    args = SimpleNamespace(dir_species_protein=str(tmp_path), ncpu=1, gene_size_quantiles="0.5")
    lengths = mod.get_concatenated_fx2tab(args)
    assert lengths["#id"].tolist() == gene_ids
    assert lengths["length"].dtype == numpy.dtype("int64")
    original = lengths.copy(deep=True)
    families = [f"OG{i}" for i in range(len(gene_ids))]
    orthogroups = pandas.DataFrame({"Orthogroup": families, "Genus_species": gene_ids})
    counts = pandas.DataFrame({"Orthogroup": families, "Genus_species": 1, "Total": 1})
    selected = mod.get_df_gc_original(orthogroups, counts, lengths, args)
    assert selected["geneid_0.5"].tolist() == gene_ids
    assert capsys.readouterr().err == ""
    pandas.testing.assert_frame_equal(lengths, original)


@pytest.mark.parametrize("ids", [["001", "1"], ["NA", "NULL", "nan", "001", "1"]])
def test_literal_ids_survive_annotation_and_orthogroup_selection(monkeypatch, tmp_path, ids, capsys):
    mod = load_module()
    monkeypatch.chdir(tmp_path)
    proteins = tmp_path / "proteins"
    proteins.mkdir()
    (proteins / "Genus_species.fa").write_text("".join(f">{gene}\nMA\n" for gene in ids))
    og_dir = tmp_path / "Orthogroups"
    og_dir.mkdir()
    membership = og_dir / "Orthogroups.tsv"
    count_path = og_dir / "Orthogroups.GeneCount.tsv"
    pandas.DataFrame({"Orthogroup": ids, "Genus_species": ids}).to_csv(membership, sep="\t", index=False)
    pandas.DataFrame({"Orthogroup": ids, "Genus_species": 1, "Total": 1}).to_csv(count_path, sep="\t", index=False)
    original = [path.read_bytes() for path in (membership, count_path)]
    db_prefix = tmp_path / "uniprot"
    for suffix in (".pin", ".phr", ".psq"):
        db_prefix.with_suffix(suffix).touch()

    def command(command, **kwargs):
        if command[:2] == ["seqkit", "fx2tab"]:
            return ("#id\tlength\n" + "".join(f"{gene}\t2\n" for gene in ids)).encode()
        if command[0] == "blastp":
            Path(command[command.index("-out") + 1]).write_text(
                "".join(f"{gene}\thit for {gene}\t1e-30\t100\n" for gene in ids))
        return b""

    monkeypatch.setattr(mod, "run_command", command)
    args = SimpleNamespace(
        dir_species_protein=str(proteins), dir_orthofinder_og=str(og_dir), gene_size_quantiles="0.5",
        ncpu=1, annotation_search_method="blastp", path_search_db=str(db_prefix), evalue="1e-2",
        min_gene_num=1, max_gene_num=10, min_species_num=1, min_percent_species_coverage=0,
        remove_unannotated=True,
    )
    annotated_path = mod.prepare_annotation(args)
    _, selected, species = mod.select_orthogroups(args, annotated_path)
    assert selected["Orthogroup"].tolist() == ids
    assert selected["geneid_0.5"].tolist() == ids
    assert selected["besthit_0.5"].tolist() == [f"hit for {gene}" for gene in ids]
    assert selected["Total"].dtype == numpy.dtype("int64")
    assert species == ["Genus_species"]
    result = pandas.read_csv(og_dir / "Orthogroups.selected.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert result["Orthogroup"].tolist() == ids
    assert result["Genus_species"].tolist() == ids
    assert [path.read_bytes() for path in (membership, count_path)] == original
    assert capsys.readouterr().err == ""


@pytest.mark.parametrize("suffix", [".tsv", ".tsv.gz"])
def test_identifier_readers_keep_empty_fields_and_numeric_missingness(tmp_path, monkeypatch, suffix):
    mod = load_module()
    path = tmp_path / ("counts" + suffix)
    frame = pandas.DataFrame({
        "Orthogroup": ["001", "NA", ""], "Genus_species": [1.0, numpy.nan, 2.0],
        "Total": [1.0, numpy.nan, 2.0], "geneid_0.5": ["NULL", "001", ""],
        "besthit_0.5": ["NA", "ordinary title", ""],
    })
    frame.to_csv(path, sep="\t", index=False)
    original = path.read_bytes()
    out = mod._load_orthogroup_table(path)
    assert out["Orthogroup"].iloc[:2].tolist() == ["001", "NA"]
    assert out["geneid_0.5"].iloc[:2].tolist() == ["NULL", "001"]
    assert pandas.isna(out["Orthogroup"].iloc[2])
    assert pandas.isna(out["geneid_0.5"].iloc[2])
    assert out["Genus_species"].dtype == numpy.dtype("float64")
    assert out["Total"].dtype == numpy.dtype("float64")
    assert out["Genus_species"].iloc[[0, 2]].tolist() == [1.0, 2.0]
    assert pandas.isna(out["Genus_species"].iloc[1])
    assert out["besthit_0.5"].isna().tolist() == [True, False, True]
    assert path.read_bytes() == original

    membership = tmp_path / ("genes" + suffix)
    pandas.DataFrame({"Orthogroup": ["001", "NA", ""], "sp": ["NULL", "001", ""]}).to_csv(
        membership, sep="\t", index=False)
    genes = mod._load_orthogroup_table(membership, gene_membership=True)
    assert genes["sp"].iloc[:2].tolist() == ["NULL", "001"]
    assert pandas.isna(genes["sp"].iloc[2])

    monkeypatch.setattr(mod, "get_species_protein_files", lambda _: ["Genus_species.fa"])
    monkeypatch.setattr(mod, "run_command", lambda *args, **kwargs: b"#id\tlength\nNA\t2\n\tNA\n")
    lengths = mod.get_concatenated_fx2tab(SimpleNamespace(dir_species_protein=str(tmp_path), ncpu=1))
    assert lengths["#id"].iloc[0] == "NA"
    assert pandas.isna(lengths["#id"].iloc[1])
    assert lengths["length"].dtype == numpy.dtype("float64")
    assert lengths["length"].iloc[0] == 2.0
    assert pandas.isna(lengths["length"].iloc[1])


def test_besthit_reader_preserves_ids_but_retains_missing_titles(tmp_path):
    mod = load_module()
    path = tmp_path / "hits.tsv"
    path.write_text("001\ttitle\t1e-10\t200\nNA\tNA\t1e-10\t200\n\ttitle\t1e-10\t200\n")
    out = mod._load_besthit_table(path)
    assert out["qseqid"].iloc[:2].tolist() == ["001", "NA"]
    assert pandas.isna(out["qseqid"].iloc[2])
    assert pandas.isna(out["stitle"].iloc[1])


def test_attach_besthits_maps_gene_columns_without_repeated_merges():
    mod = load_module()
    df_gc = pandas.DataFrame(
        [
            {"Orthogroup": "OG1", "geneid_0.25": "g1", "geneid_0.75": "g2"},
            {"Orthogroup": "OG2", "geneid_0.25": "g3", "geneid_0.75": "missing"},
        ]
    )
    df_diamond = pandas.DataFrame(
        [
            {"qseqid": "g1", "stitle": "hit1"},
            {"qseqid": "g2", "stitle": "hit2"},
            {"qseqid": "g3", "stitle": "hit3_first"},
            {"qseqid": "g3", "stitle": "hit3_second"},
        ]
    )
    out = mod._attach_besthits(df_gc, df_diamond, ["geneid_0.25", "geneid_0.75"])
    assert out.loc[out["Orthogroup"] == "OG1", "besthit_0.25"].iloc[0] == "hit1"
    assert out.loc[out["Orthogroup"] == "OG1", "besthit_0.75"].iloc[0] == "hit2"
    assert out.loc[out["Orthogroup"] == "OG2", "besthit_0.25"].iloc[0] == "hit3_first"
    assert pandas.isna(out.loc[out["Orthogroup"] == "OG2", "besthit_0.75"].iloc[0])


def test_annotation_search_does_not_reuse_existing_output(tmp_path, monkeypatch):
    mod = load_module()
    monkeypatch.chdir(tmp_path)
    (tmp_path / "tmp.blastp.out.tsv").write_text(
        "g1\tstale hit\t1e-50\t200\n",
        encoding="utf-8",
    )
    db_prefix = tmp_path / "uniprot"
    for extension in (".pin", ".phr", ".psq"):
        (tmp_path / ("uniprot" + extension)).write_text("db\n", encoding="utf-8")
    calls = []

    def fake_run_command(command, stdout_file=None, append=False, tool_name="tool", capture_output=False):
        calls.append(command)
        if command[0] == "blastp":
            Path(command[command.index("-out") + 1]).write_text(
                "g1\tfresh hit\t1e-20\t100\n",
                encoding="utf-8",
            )
        return b""

    monkeypatch.setattr(mod, "run_command", fake_run_command)
    monkeypatch.setattr(mod, "get_species_protein_files", lambda _path: [])
    df_gc = pandas.DataFrame([{"Orthogroup": "OG1", "geneid_0.5": "g1"}])
    args = SimpleNamespace(
        annotation_search_method="blastp",
        gene_size_quantiles="0.5",
        dir_species_protein=str(tmp_path),
        ncpu=1,
        path_search_db=str(db_prefix),
        evalue="1e-2",
    )
    out = mod.annotate_representative_genes(df_gc, args)
    assert "besthit_0.5" in out.columns
    assert out["besthit_0.5"].iloc[0] == "fresh hit"
    assert any(command[0] == "blastp" for command in calls)


def test_prepare_annotation_does_not_reuse_stale_quartile_table(tmp_path, monkeypatch):
    mod = load_module()
    monkeypatch.chdir(tmp_path)
    orthogroup_dir = tmp_path / "Orthogroups_filtered"
    orthogroup_dir.mkdir()
    pandas.DataFrame(
        [{"Orthogroup": "OG_FRESH", "spA": "geneA"}]
    ).to_csv(orthogroup_dir / "Orthogroups.tsv", sep="\t", index=False)
    pandas.DataFrame(
        [{"Orthogroup": "OG_FRESH", "spA": 1, "Total": 1}]
    ).to_csv(
        orthogroup_dir / "Orthogroups.GeneCount.tsv",
        sep="\t",
        index=False,
    )
    stale = tmp_path / "tmp.Orthogroups.GeneCount.quartile_genes.tsv"
    pandas.DataFrame(
        [{"Orthogroup": "OG_STALE", "spA": 1, "Total": 1}]
    ).to_csv(stale, sep="\t", index=False)
    monkeypatch.setattr(
        mod,
        "get_concatenated_fx2tab",
        lambda _args: pandas.DataFrame([{"#id": "geneA", "length": 10}]),
    )

    def add_fresh_quantile(df_og, df_gc, _fx2tab, _args):
        assert df_og["Orthogroup"].tolist() == ["OG_FRESH"]
        result = df_gc.copy()
        result["geneid_0.5"] = "geneA"
        return result

    monkeypatch.setattr(mod, "get_df_gc_original", add_fresh_quantile)
    monkeypatch.setattr(
        mod,
        "annotate_representative_genes",
        lambda frame, _args: frame.assign(**{"besthit_0.5": "fresh hit"}),
    )
    args = SimpleNamespace(
        dir_orthofinder_og=str(orthogroup_dir),
        gene_size_quantiles="0.5",
    )

    annotated_path = mod.prepare_annotation(args)

    annotated = pandas.read_csv(annotated_path, sep="\t")
    assert annotated["Orthogroup"].tolist() == ["OG_FRESH"]
    refreshed_quartile = pandas.read_csv(stale, sep="\t")
    assert refreshed_quartile["Orthogroup"].tolist() == ["OG_FRESH"]


def test_remove_tmp_files_unlinks_dangling_symlink(tmp_path, monkeypatch):
    mod = load_module()
    monkeypatch.chdir(tmp_path)
    dangling = tmp_path / "tmp.query.fasta"
    dangling.symlink_to(tmp_path / "missing-target")

    mod.remove_tmp_files()

    assert not dangling.exists()
    assert not dangling.is_symlink()


def test_load_besthit_table_reads_only_first_two_columns(tmp_path):
    mod = load_module()
    infile = tmp_path / "tmp.mmseqs2.out.tsv"
    infile.write_text(
        "geneA\thitA\t1e-50\t200\n"
        "geneB\thitB\t1e-20\t150\n",
        encoding="utf-8",
    )

    out = mod._load_besthit_table(str(infile))

    assert out.columns.tolist() == ["qseqid", "stitle"]
    assert out.to_dict("records") == [
        {"qseqid": "geneA", "stitle": "hitA"},
        {"qseqid": "geneB", "stitle": "hitB"},
    ]


def test_mmseqs_annotation_search_passes_split_memory_limit(tmp_path, monkeypatch):
    mod = load_module()
    monkeypatch.chdir(tmp_path)

    protein_dir = tmp_path / "proteins"
    protein_dir.mkdir()
    (protein_dir / "spA.fa").write_text(">geneA\nMA\n", encoding="utf-8")
    mmseqs_prefix = tmp_path / "uniprot.mmseqs"
    mmseqs_prefix.write_text("db\n", encoding="utf-8")
    (tmp_path / "uniprot.mmseqs.dbtype").write_text("0\n", encoding="utf-8")

    calls = []

    def fake_run_command(command, stdout_file=None, append=False, tool_name="tool", capture_output=False):
        calls.append(command)
        if stdout_file is not None:
            mode = "a" if append else "w"
            with open(stdout_file, mode, encoding="utf-8") as handle:
                handle.write(">geneA\nMA\n")
        if command[:2] == ["mmseqs", "convertalis"]:
            Path(command[5]).write_text("geneA\thitA\t1e-20\t100\n", encoding="utf-8")
        return b""

    monkeypatch.setattr(mod, "run_command", fake_run_command)
    args = SimpleNamespace(
        annotation_search_method="mmseqs2",
        gene_size_quantiles="0.5",
        dir_species_protein=str(protein_dir),
        ncpu=2,
        path_search_db=str(tmp_path / "uniprot"),
        evalue="1e-2",
        mmseqs_split_memory_limit="12G",
    )

    df_gc = pandas.DataFrame([{"Orthogroup": "OG1", "geneid_0.5": "geneA"}])
    out = mod.annotate_representative_genes(df_gc, args)

    search_calls = [call for call in calls if call[:2] == ["mmseqs", "search"]]
    assert len(search_calls) == 1
    assert "--split-memory-limit" in search_calls[0]
    limit_index = search_calls[0].index("--split-memory-limit")
    assert search_calls[0][limit_index + 1] == "12G"
    assert out.loc[0, "besthit_0.5"] == "hitA"


def test_remove_tmp_files_only_removes_known_selection_files(tmp_path, monkeypatch):
    mod = load_module()
    monkeypatch.chdir(tmp_path)
    known = tmp_path / "tmp.query.fasta"
    unrelated = tmp_path / "tmp.user.tsv"
    known.write_text("", encoding="utf-8")
    unrelated.write_text("keep\n", encoding="utf-8")

    mod.remove_tmp_files()

    assert not known.exists()
    assert unrelated.exists()


def test_resolve_mmseqs_db_prefix_rejects_multiline_path():
    mod = load_module()

    with pytest.raises(ValueError, match="contains control characters"):
        mod._resolve_mmseqs_db_prefix("createdb /tmp/uniprot.pep\n/tmp/uniprot.mmseqs")


def test_resolve_mmseqs_db_prefix_requires_existing_db_files(tmp_path):
    mod = load_module()
    prefix = tmp_path / "uniprot_sprot"

    with pytest.raises(FileNotFoundError, match="MMseqs2 search DB files were not found"):
        mod._resolve_mmseqs_db_prefix(str(prefix))

    mmseqs_prefix = tmp_path / "uniprot_sprot.mmseqs"
    mmseqs_prefix.write_text("db\n", encoding="utf-8")
    (tmp_path / "uniprot_sprot.mmseqs.dbtype").write_text("0\n", encoding="utf-8")

    assert mod._resolve_mmseqs_db_prefix(str(prefix)) == str(mmseqs_prefix)


def test_resolve_blast_db_prefix_allows_spaces_but_rejects_control_characters(tmp_path):
    mod = load_module()
    prefix = tmp_path / "uniprot sprot"

    with pytest.raises(ValueError, match="contains control characters"):
        mod._resolve_blast_db_prefix(str(prefix) + "\nlog")
    with pytest.raises(FileNotFoundError, match="BLASTP search DB files were not found"):
        mod._resolve_blast_db_prefix(str(prefix))

    for ext in (".pin", ".phr", ".psq"):
        (tmp_path / ("uniprot sprot" + ext)).write_text("db\n", encoding="utf-8")

    assert mod._resolve_blast_db_prefix(str(prefix)) == str(prefix)


@pytest.mark.parametrize("remove_unannotated", [False, True])
def test_all_computable_selection_includes_unannotated_and_large_families(tmp_path, remove_unannotated):
    mod = load_module()
    og = tmp_path / "orthogroups"
    og.mkdir()
    rows = [("large", 1001, 2), ("small", 4, 2), ("too_short", 3, 2), ("one_species", 10, 1)]
    counts, members = [], []
    for family, total, species in rows:
        second = 1 if species == 2 else 0
        counts.append({"Orthogroup": family, "A": total - second, "B": second,
                       "Total": total, "geneid_0.5": family + "_a0", "besthit_0.5": ""})
        members.append({"Orthogroup": family,
                        "A": ", ".join(f"{family}_a{i}" for i in range(total - second)),
                        "B": family + "_b0" if second else ""})
    annotated = og / "Orthogroups.GeneCount.annotated.tsv"
    pandas.DataFrame(counts).to_csv(annotated, sep="\t", index=False)
    pandas.DataFrame(members).to_csv(og / "Orthogroups.tsv", sep="\t", index=False)
    args = SimpleNamespace(dir_orthofinder_og=str(og), min_gene_num=4, max_gene_num=9930,
                           min_species_num=2, min_percent_species_coverage=0,
                           remove_unannotated=remove_unannotated, annotation_search_method="mmseqs2")
    before = annotated.read_bytes()
    _, selected, species = mod.select_orthogroups(args, str(annotated))
    expected = [] if remove_unannotated else ["large", "small"]
    assert selected["Orthogroup"].tolist() == expected
    membership = pandas.read_csv(og / "Orthogroups.selected.tsv", sep="\t", keep_default_na=False)
    assert membership["Orthogroup"].tolist() == expected
    assert species == ["A", "B"]
    assert annotated.read_bytes() == before


@pytest.mark.parametrize("failure", ["", "excluded", "changed_counts", "duplicates"])
def test_selection_prefix_preserves_indices_and_refuses_changed_inputs(tmp_path, failure):
    mod = load_module()
    og = tmp_path / "og"
    og.mkdir()
    rows = [{"Orthogroup": f"OG{i}", "A": i + 2, "B": 2, "Total": i + 4} for i in range(3)]
    annotated = og / "annotated.tsv"
    pandas.DataFrame(rows).to_csv(annotated, sep="\t", index=False)
    pandas.DataFrame([{"Orthogroup": r["Orthogroup"], "A": "a", "B": "b"} for r in rows]).to_csv(
        og / "Orthogroups.tsv", sep="\t", index=False)
    prefix_rows = [rows[2].copy(), rows[0].copy()]
    if failure == "excluded":
        prefix_rows[0]["Orthogroup"] = "absent"
    if failure == "changed_counts":
        prefix_rows[0]["A"] = 999
    if failure == "duplicates":
        prefix_rows.append(prefix_rows[0].copy())
    prefix = tmp_path / "prefix.tsv"
    pandas.DataFrame(prefix_rows).to_csv(prefix, sep="\t", index=False)
    args = SimpleNamespace(dir_orthofinder_og=str(og), min_gene_num=4, max_gene_num=9930,
                           min_species_num=2, min_percent_species_coverage=0,
                           remove_unannotated=False, annotation_search_method="mmseqs2",
                           selection_prefix=str(prefix))
    if failure:
        with pytest.raises(ValueError):
            mod.select_orthogroups(args, str(annotated))
        assert not (og / "Orthogroups.GeneCount.selected.tsv").exists()
        assert not (og / "Orthogroups.selected.tsv").exists()
    else:
        _, selected, _ = mod.select_orthogroups(args, str(annotated))
        assert selected["Orthogroup"].tolist() == ["OG2", "OG0", "OG1"]
        assert pandas.read_csv(og / "Orthogroups.selected.tsv", sep="\t")["Orthogroup"].tolist() == ["OG2", "OG0", "OG1"]


@pytest.mark.parametrize("failure", ["", "counts", "representative", "membership", "duplicate", "quantiles", "stale", "parameters", "uncovered", "native_source"])
def test_annotation_cache_is_bound_to_complete_native_membership(tmp_path, failure, monkeypatch):
    monkeypatch.syspath_prepend(str(SCRIPT_PATH.parent))
    mod = load_module()
    root = tmp_path / "native"
    root.mkdir()
    counts = pandas.DataFrame([{"Orthogroup": "OG1", "A": 2, "B": 2, "Total": 4}])
    counts.to_csv(root / "Orthogroups.GeneCount.tsv", sep="\t", index=False)
    pandas.DataFrame([{"Orthogroup": "OG1", "A": "a1, a2", "B": "b1, b2"}]).to_csv(
        root / "Orthogroups.tsv", sep="\t", index=False)
    cached = counts.assign(geneid_0_5="a1", besthit_0_5="known protein").rename(
        columns={"geneid_0_5": "geneid_0.5", "besthit_0_5": "besthit_0.5"})
    if failure == "counts":
        cached.loc[0, "Total"] = 99
    if failure == "representative":
        cached.loc[0, "geneid_0.5"] = "a_different_family_gene"
    if failure == "membership":
        counts.loc[0, "A"] = 3
        counts.loc[0, "Total"] = 5
        counts.to_csv(root / "Orthogroups.GeneCount.tsv", sep="\t", index=False)
        cached.loc[0, ["A", "Total"]] = [3, 5]
    if failure == "duplicate":
        cached = pandas.concat([cached, cached], ignore_index=True)
    if failure == "quantiles":
        cached = cached.drop(columns="besthit_0.5")
    source = tmp_path / "preserved-annotation.tsv"
    cached.to_csv(source, sep="\t", index=False)
    before = source.read_bytes()
    args = SimpleNamespace(dir_orthofinder_og=str(root), gene_size_quantiles="0.5", annotation_cache=str(source))
    _annotation_cache_proof(args, source, tmp_path)
    if failure == "stale":
        source.write_text(source.read_text() + "\n")
        before = source.read_bytes()
    if failure == "parameters":
        args.evalue = "1e-3"
    if failure == "uncovered":
        other = tmp_path / "foreign.tsv"
        other.write_bytes(source.read_bytes())
        args.annotation_cache = str(other)
    if failure == "native_source":
        (root / "Orthogroups.tsv").write_text((root / "Orthogroups.tsv").read_text().replace("a2", "a3"))
    target = root / "Orthogroups.GeneCount.annotated.tsv"
    if failure:
        with pytest.raises(ValueError):
            mod.prepare_annotation(args)
        assert not target.exists()
    else:
        assert mod.prepare_annotation(args) == str(target)
        assert target.read_bytes() == before
    assert source.read_bytes() == before


def test_additional_catalog_leaves_existing_family_provenance_current(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(SCRIPT_PATH.parent))
    mod = load_module()
    support = Path(__file__).resolve().parents[1] / "support"
    workspace = tmp_path / "workspace"
    old = workspace / "output/orthofinder/Orthogroups_filtered"
    new = workspace / "output/orthofinder/additional"
    old.mkdir(parents=True)
    new.mkdir()
    counts = pandas.DataFrame([{"Orthogroup": name, "A": 2, "B": 2, "Total": 4} for name in ["old", "new"]])
    members = pandas.DataFrame([{"Orthogroup": name, "A": f"{name}_a1, {name}_a2",
                                "B": f"{name}_b1, {name}_b2"} for name in ["old", "new"]])
    counts.iloc[:1].to_csv(old / "Orthogroups.GeneCount.selected.tsv", sep="\t", index=False)
    old_members = old / "Orthogroups.selected.tsv"
    members.iloc[:1].to_csv(old_members, sep="\t", index=False)
    counts.to_csv(new / "Orthogroups.GeneCount.tsv", sep="\t", index=False)
    members.to_csv(new / "Orthogroups.tsv", sep="\t", index=False)
    cache = old / "Orthogroups.GeneCount.annotated.tsv"
    counts.assign(**{"geneid_0.5": ["old_a1", "new_a1"], "besthit_0.5": ["hit", ""]}).to_csv(cache, sep="\t", index=False)
    store = workspace / "output/orthogroup"
    (store / "artifact_provenance").mkdir(parents=True)
    fasta = store / "cds_fasta/old.fa"
    fasta.parent.mkdir()
    fasta.write_text(">old_a1\nATG\n")
    contract = ["--manifest", str(store / "artifact_provenance/old.primary_fasta.json"),
                "--step", "primary_fasta", "--family-id", "old", "--logical-root", str(store),
                "--workspace-root", str(workspace), "--input", "orthogroups=" + str(old_members),
                "--output", "primary_fasta=" + str(fasta)]
    def provenance(command):
        return subprocess.run([sys.executable, str(support / "artifact_provenance.py"), command, *contract],
                              capture_output=True, text=True)
    recorded = provenance("record")
    assert recorded.returncode == 0, recorded.stdout + recorded.stderr
    assert provenance("needs-run").returncode == 1
    before = {p.name: p.read_bytes() for p in old.iterdir()}
    args = SimpleNamespace(dir_orthofinder_og=str(new), annotation_cache=str(cache), gene_size_quantiles="0.5",
                           selection_prefix=str(old / "Orthogroups.GeneCount.selected.tsv"), min_gene_num=4,
                           max_gene_num=9930, min_species_num=2, min_percent_species_coverage=0,
                           remove_unannotated=False, annotation_search_method="mmseqs2")
    _annotation_cache_proof(args, cache, workspace)
    annotated = mod.prepare_annotation(args)
    _, selected, _ = mod.select_orthogroups(args, annotated)
    assert selected["Orthogroup"].tolist() == ["old", "new"]
    assert {p.name: p.read_bytes() for p in old.iterdir()} == before
    current = provenance("needs-run")
    assert current.returncode == 1, current.stdout + current.stderr


def _annotation_cache_proof(args, cache, root):
    manifest = root / "annotation-source.json"
    args.annotation_search_method = "mmseqs2"
    args.evalue = "1e-2"
    args.annotation_cache_manifest = str(manifest)
    args.annotation_cache_workspace = str(root)
    args.annotation_cache_logical_root = str(root / "output/genome_evolution")
    Path(args.annotation_cache_logical_root).mkdir(parents=True, exist_ok=True)
    directory = Path(args.dir_orthofinder_og)
    command = [sys.executable, str(SCRIPT_PATH.parent / "artifact_provenance.py"), "record",
        "--manifest", str(manifest), "--step", "orthogroup_selection", "--family-id", "all_species",
        "--logical-root", args.annotation_cache_logical_root, "--workspace-root", str(root),
        "--input", "orthogroups=" + str(directory / "Orthogroups.tsv"),
        "--input", "gene_counts=" + str(directory / "Orthogroups.GeneCount.tsv"),
        "--output", "annotated=" + str(cache), "--parameter", "annotation_search_method=mmseqs2",
        "--parameter", "evalue=1e-2", "--parameter", "gene_size_quantiles=" + args.gene_size_quantiles]
    result = subprocess.run(command, text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr

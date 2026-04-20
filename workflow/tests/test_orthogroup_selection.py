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


def test_get_species_protein_files_includes_gzipped_fasta(tmp_path):
    mod = load_module()
    (tmp_path / "sp1.fa").write_text(">a\nAT\n", encoding="utf-8")
    (tmp_path / "sp2.fasta.gz").write_text("", encoding="utf-8")
    (tmp_path / "sp3.fa.gz").write_text("", encoding="utf-8")
    (tmp_path / "sp4.fas").write_text(">a\nAT\n", encoding="utf-8")
    (tmp_path / "sp5.fna.gz").write_text("", encoding="utf-8")
    (tmp_path / "README.txt").write_text("ignore", encoding="utf-8")

    files = mod.get_species_protein_files(str(tmp_path))
    assert set(files) == {"sp1.fa", "sp2.fasta.gz", "sp3.fa.gz", "sp4.fas", "sp5.fna.gz"}


def test_get_species_protein_files_ignores_hidden_entries_and_directories(tmp_path):
    mod = load_module()
    (tmp_path / "sp1.fa").write_text(">a\nAT\n", encoding="utf-8")
    (tmp_path / ".hidden.fa").write_text(">a\nAT\n", encoding="utf-8")
    (tmp_path / "sp2.fasta.gz").mkdir()

    files = mod.get_species_protein_files(str(tmp_path))
    assert files == ["sp1.fa"]


def test_get_species_protein_files_returns_sorted_file_names(tmp_path):
    mod = load_module()
    (tmp_path / "spB.fa").write_text(">a\nAT\n", encoding="utf-8")
    (tmp_path / "spA.fa").write_text(">a\nAT\n", encoding="utf-8")
    files = mod.get_species_protein_files(str(tmp_path))
    assert files == ["spA.fa", "spB.fa"]


def test_format_command_for_log_escapes_control_characters():
    mod = load_module()
    command_log = mod._format_command_for_log(["mmseqs", "search", "clean", "bad\npath"])

    assert "\\n" in command_log
    assert "bad\npath" not in command_log


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


def test_annotation_search_handles_empty_existing_output(tmp_path, monkeypatch):
    mod = load_module()
    monkeypatch.chdir(tmp_path)
    (tmp_path / "tmp.blastp.out.tsv").write_text("", encoding="utf-8")
    df_gc = pandas.DataFrame([{"Orthogroup": "OG1", "geneid_0.5": "g1"}])
    args = SimpleNamespace(
        annotation_search_method="blastp",
        gene_size_quantiles="0.5",
        dir_species_protein=str(tmp_path),
        ncpu=1,
        path_search_db="dummy",
        evalue="1e-2",
    )
    out = mod.annotate_representative_genes(df_gc, args)
    assert "besthit_0.5" in out.columns
    assert out["besthit_0.5"].iloc[0] == ""


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

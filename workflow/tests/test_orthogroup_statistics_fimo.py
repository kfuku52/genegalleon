import os
import subprocess
import sys
import types
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "orthogroup_statistics.py"


def load_module():
    if "kftools" not in sys.modules:
        pkg = types.ModuleType("kftools")
        pkg.__path__ = []
        sys.modules["kftools"] = pkg
    if "kftools.kfog" not in sys.modules:
        sys.modules["kftools.kfog"] = types.ModuleType("kftools.kfog")
    if "kftools.kfphylo" not in sys.modules:
        sys.modules["kftools.kfphylo"] = types.ModuleType("kftools.kfphylo")
    spec = spec_from_file_location("orthogroup_statistics", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_kftools_stub(stub_root):
    package_dir = stub_root / "kftools"
    package_dir.mkdir(parents=True)
    (package_dir / "__init__.py").write_text("", encoding="utf-8")
    (package_dir / "kfog.py").write_text("", encoding="utf-8")


@pytest.mark.parametrize("identifiers", [("001", "1"), ("NA", "NULL", "nan")])
def test_character_gff_literal_tip_names_join(tmp_path, identifiers):
    tree = tmp_path / "rooted.nwk"
    tree.write_text("(" + ",".join(f"{name}:1" for name in identifiers) + ")root;\n")
    traits = tmp_path / "traits.tsv"
    traits.write_text("gene_id\tnum_intron\n" + "".join(
        f"{name}\t{i}" + "\n" for i, name in enumerate(identifiers)))
    stub_root = tmp_path / "stub_packages"
    write_kftools_stub(stub_root)
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(filter(None, (str(stub_root), env.get("PYTHONPATH", ""))))
    proc = subprocess.run([
        sys.executable, str(SCRIPT_PATH), "--rooted_tree", str(tree), "--character_gff", str(traits),
    ], cwd=tmp_path, capture_output=True, text=True, env=env)
    assert proc.returncode == 0, proc.stderr
    branches = pandas.read_csv(tmp_path / "orthogroup.branch.tsv", sep="\t", converters={"node_name": str})
    tips = branches.set_index("node_name")
    assert [tips.loc[name, "num_intron"] for name in identifiers] == list(range(len(identifiers)))


@pytest.mark.parametrize("identifiers", [("001", "1"), ("NA", "NULL", "nan")])
@pytest.mark.parametrize("source", ["uniprot", "rpsblast", "expression"])
def test_branch_annotations_join_literal_tip_names(tmp_path, identifiers, source):
    tree = tmp_path / "rooted.nwk"
    tree.write_text("(" + ",".join(f"{name}:1" for name in identifiers) + ")root;\n")
    table = tmp_path / "annotation.tsv"
    if source == "uniprot":
        table.write_text("gene_id\tannotation_score\ttitle\n" + "".join(
            f"{name}\t{i + 1}\tNA\n" for i, name in enumerate(identifiers)))
        column = "annotation_score"
        expected = list(range(1, len(identifiers) + 1))
    elif source == "rpsblast":
        table.write_text("qacc\tstitle\tevalue\n" + "".join(
            f"{name}\tdomain{i + 1}\t1e-10\n" for i, name in enumerate(identifiers)))
        column = "pfam_domain"
        expected = [f"domain{i + 1}" for i in range(len(identifiers))]
    else:
        table.write_text("gene_id\ttissue1\ttissue2\n" + "".join(
            f"{name}\t{i + 1}\t{i + 3}\n" for i, name in enumerate(identifiers)))
        column = "expression_tissue1"
        expected = list(range(1, len(identifiers) + 1))
    stub_root = tmp_path / "stub_packages"
    write_kftools_stub(stub_root)
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(filter(None, (str(stub_root), env.get("PYTHONPATH", ""))))
    proc = subprocess.run([
        sys.executable, str(SCRIPT_PATH), "--rooted_tree", str(tree), "--" + source, str(table),
    ], cwd=tmp_path, capture_output=True, text=True, env=env)
    assert proc.returncode == 0, proc.stderr
    branches = pandas.read_csv(tmp_path / "orthogroup.branch.tsv", sep="\t", converters={"node_name": str})
    tips = branches.set_index("node_name")
    assert [tips.loc[name, column] for name in identifiers] == expected
    if source == "uniprot":
        assert branches["title"].isna().all()


def test_expression_preserves_implicit_gene_index(tmp_path):
    identifiers = ("Genus_species_gene1", "Other_species_gene2")
    tree = tmp_path / "rooted.nwk"
    tree.write_text("(" + ",".join(f"{name}:1" for name in identifiers) + ")root;\n")
    table = tmp_path / "expression.tsv"
    # Sample-only headers already make pandas infer the leading gene index.
    table.write_text("tissue1\ttissue2\nGenus_species_gene1\t1\t3\nOther_species_gene2\t2\t4\n")
    stub_root = tmp_path / "stub_packages"
    write_kftools_stub(stub_root)
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(filter(None, (str(stub_root), env.get("PYTHONPATH", ""))))
    proc = subprocess.run([
        sys.executable, str(SCRIPT_PATH), "--rooted_tree", str(tree), "--expression", str(table),
    ], cwd=tmp_path, capture_output=True, text=True, env=env)
    assert proc.returncode == 0, proc.stderr
    branches = pandas.read_csv(tmp_path / "orthogroup.branch.tsv", sep="\t").set_index("node_name")
    assert [branches.loc[name, "expression_tissue1"] for name in identifiers] == [1, 2]
    assert [branches.loc[name, "expression_tissue2"] for name in identifiers] == [3, 4]


def test_no_significant_pfam_hits_preserve_branch_table_schema(tmp_path):
    rooted_tree = tmp_path / "rooted.nwk"
    rooted_tree.write_text("(Species_A_gene1:1,Species_B_gene2:1)n0;\n", encoding="utf-8")
    rpsblast = tmp_path / "rpsblast.tsv"
    pandas.DataFrame(
        [
            {"qacc": "Species_A_gene1", "evalue": float("nan"), "stitle": None},
            {"qacc": "Species_B_gene2", "evalue": float("nan"), "stitle": None},
        ]
    ).to_csv(rpsblast, sep="\t", index=False)

    stub_root = tmp_path / "stub_packages"
    write_kftools_stub(stub_root)
    env = os.environ.copy()
    existing_pythonpath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = (
        f"{stub_root}{os.pathsep}{existing_pythonpath}"
        if existing_pythonpath
        else str(stub_root)
    )

    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--rooted_tree",
            str(rooted_tree),
            "--rpsblast",
            str(rpsblast),
        ],
        cwd=str(tmp_path),
        capture_output=True,
        text=True,
        env=env,
    )

    assert proc.returncode == 0, proc.stderr
    branch_table = pandas.read_csv(tmp_path / "orthogroup.branch.tsv", sep="\t")
    assert "pfam_domain" in branch_table.columns
    assert branch_table["pfam_domain"].isna().all()


def test_expression_gene_ids_join_to_tree_tips(tmp_path):
    rooted_tree = tmp_path / "rooted.nwk"
    rooted_tree.write_text("(Species_A_gene1:1,Species_B_gene2:1)n0;\n", encoding="utf-8")
    expression = tmp_path / "expression.tsv"
    expression.write_text(
        "gene_id\ttissue1_1\ttissue1_2\n"
        "Species_A_gene1\t2\t4\n"
        "Species_B_gene2\t6\t8\n",
        encoding="utf-8",
    )

    stub_root = tmp_path / "stub_packages"
    write_kftools_stub(stub_root)
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(filter(None, (str(stub_root), env.get("PYTHONPATH", ""))))

    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--rooted_tree",
            str(rooted_tree),
            "--expression",
            str(expression),
        ],
        cwd=str(tmp_path),
        capture_output=True,
        text=True,
        env=env,
    )

    assert proc.returncode == 0, proc.stderr
    branches = pandas.read_csv(tmp_path / "orthogroup.branch.tsv", sep="\t")
    tips = branches.set_index("node_name")
    assert tips.loc["Species_A_gene1", "expression_tissue1_1"] == 2
    assert tips.loc["Species_A_gene1", "expression_tissue1_2"] == 4
    assert tips.loc["Species_B_gene2", "expression_tissue1_1"] == 6
    assert tips.loc["Species_B_gene2", "expression_tissue1_2"] == 8


def test_load_fimo_hits_parses_modern_fimo_tsv(tmp_path):
    mod = load_module()
    infile = tmp_path / "fimo.tsv"
    infile.write_text(
        (
            "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n"
            "M1\tALT1\tgeneA\t10\t15\t+\t11.5\t1e-4\t1e-3\tACGTAC\n"
            "M2\t.\tgeneB\t20\t24\t-\t9.3\t1e-3\t2e-2\tTTGGCC\n"
        ),
        encoding="utf-8",
    )

    df = mod.load_fimo_hits(str(infile))
    assert df.shape[0] == 2
    assert df.columns.tolist() == [
        "motif_id",
        "motif_alt_id",
        "sequence_name",
        "start",
        "stop",
        "strand",
        "q-value",
    ]
    assert df.loc[0, "motif_alt_id"] == "ALT1"
    # "." should fall back to motif_id for consistent downstream grouping.
    assert df.loc[1, "motif_alt_id"] == "M2"


def test_fimo_annotations_preserve_interleaved_hit_order_and_literal_singletons(tmp_path):
    tree = tmp_path / "rooted.nwk"
    tree.write_text("(NA:1,001:1,nan:1,NULL:1)root;\n")
    fimo = tmp_path / "fimo.tsv"
    fimo.write_text(
        "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tq-value\n"
        "M1\t.\tNA\t11\t16\t+\t0.01\n"
        "M2\tALT2\t001\t21\t26\t-\t0.02\n"
        "M3\t.\tNA\t31\t36\t-\t0.03\n"
        "M4\tALT4\tnan\t41\t46\t+\t0.04\n"
        "M5\t.\t001\t51\t56\t+\t0.05\n"
        "M6\tALT6\tNULL\t61\t66\t-\t0.06\n"
        "invalid\t.\tNA\tbad\t99\t.\t0.001\n"
    )
    stub_root = tmp_path / "stub_packages"
    write_kftools_stub(stub_root)
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(filter(None, (str(stub_root), env.get("PYTHONPATH", ""))))
    proc = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "--rooted_tree", str(tree), "--fimo", str(fimo)],
        cwd=tmp_path, capture_output=True, text=True, env=env,
    )
    assert proc.returncode == 0, proc.stderr
    branches = pandas.read_csv(
        tmp_path / "orthogroup.branch.tsv", sep="\t", dtype=str, keep_default_na=False,
    ).set_index("node_name")
    columns = ["fimo_id", "fimo_alt_id", "fimo_start", "fimo_end", "fimo_strand", "fimo_qvalue"]
    expected = pandas.DataFrame.from_dict({
        "NA": ["M1; M3", "M1; M3", "11; 31", "16; 36", "+; -", "0.01; 0.03"],
        "001": ["M2; M5", "ALT2; M5", "21; 51", "26; 56", "-; +", "0.02; 0.05"],
        "nan": ["M4", "ALT4", "41", "46", "+", "0.04"],
        "NULL": ["M6", "ALT6", "61", "66", "-", "0.06"],
    }, orient="index", columns=columns)
    expected.index.name = "node_name"
    pandas.testing.assert_frame_equal(branches.loc[expected.index, columns], expected, check_exact=True)


def test_orthogroup_statistics_preserves_qualified_species_prefixes():
    mod = load_module()

    assert mod.extract_species_label("Dictyostelium_cf_discoideum_geneA") == "Dictyostelium_cf_discoideum"
    assert mod.extract_species_label("Bacillus_subtilis_subsp_subtilis_geneB") == "Bacillus_subtilis_subsp_subtilis"
    assert mod.scientific_name_from_label("Dictyostelium_cf_discoideum_geneA") == "Dictyostelium cf. discoideum"


def test_sequence_stats_validate_untrimmed_and_trimmed_alignments(tmp_path, monkeypatch):
    mod = load_module()
    untrimmed = tmp_path / "untrimmed.fasta"
    untrimmed.write_text(">geneA\nAAAAA\n>geneB\nA-AA-\n", encoding="utf-8")
    trimmed = tmp_path / "cleaned.fasta"
    trimmed.write_text(">geneA\nAAA\n>geneB\nA-A\n", encoding="utf-8")

    alignment_calls = []

    def strict_alignment_stats(path):
        alignment_calls.append(path)
        if path == str(untrimmed):
            return {"num_site": 5, "num_seq": 2, "len_max": 5, "len_min": 3}
        return {"num_site": 3, "num_seq": 2, "len_max": 3, "len_min": 2}

    monkeypatch.setattr(mod.kfog, "get_aln_stats", strict_alignment_stats, raising=False)

    out = mod.collect_sequence_stats(str(untrimmed), str(trimmed))

    assert alignment_calls == [str(untrimmed), str(trimmed)]
    assert out == {
        "original_num_site": 5,
        "original_num_seq": 2,
        "original_len_max": 5,
        "original_len_min": 3,
        "cleaned_num_site": 3,
        "cleaned_num_seq": 2,
        "cleaned_len_max": 3,
        "cleaned_len_min": 2,
    }


def test_sequence_stats_parser_names_alignment_inputs_explicitly():
    mod = load_module()

    args = mod.build_arg_parser().parse_args(["--untrimmed_aln", "before.fa", "--trimmed_aln", "after.fa"])

    assert args.untrimmed_aln == "before.fa"
    assert args.trimmed_aln == "after.fa"
    assert not hasattr(args, "unaligned_aln")
    assert not hasattr(args, "trimal_aln")


def test_new_unrooted_tree_preserves_numeric_iqtree_support():
    mod = load_module()

    tree = mod.new_unrooted_tree("((A:1,B:1)95:2,C:3);")
    internal = [node for node in tree.traverse() if sorted(leaf.name for leaf in mod.iter_leaves(node)) == ["A", "B"]][0]

    assert internal.support == 95.0
    assert internal.dist == 2.0


def test_new_unrooted_tree_accepts_internal_node_names_without_support():
    mod = load_module()

    tree = mod.new_unrooted_tree("((A:1,B:1)n3:2,C:3);")
    internal = [node for node in tree.traverse() if sorted(leaf.name for leaf in mod.iter_leaves(node)) == ["A", "B"]][0]

    assert internal.name == "n3"
    assert getattr(internal, "support", None) is None
    assert internal.dist == 2.0


def test_set_outgroup_compat_clears_root_distance_before_reroot():
    mod = load_module()

    tree = mod.new_tree("((A:1,B:1)n1:2,C:3)n0:1;", format=1)

    mod.set_outgroup_compat(tree, "A")

    assert sorted(leaf.name for leaf in mod.iter_leaves(tree)) == ["A", "B", "C"]


def test_set_outgroup_compat_handles_multiple_ete4_root_assertions():
    mod = load_module()

    tree = mod.new_tree("((A:1,B:1)95:2,C:3)100:1;", format=0)

    mod.set_outgroup_compat(tree, "A")

    assert sorted(leaf.name for leaf in mod.iter_leaves(tree)) == ["A", "B", "C"]


def test_load_fimo_hits_parses_legacy_fimo_txt_header(tmp_path):
    mod = load_module()
    infile = tmp_path / "fimo.txt"
    infile.write_text(
        (
            "#pattern name\tsequence name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched sequence\n"
            "MA0123.1\tgeneC\t30\t39\t+\t7.0\t1e-3\t4e-2\tACGTTGCAAA\n"
            "MA0567.1\tgeneD\t40\t47\t-\t6.1\t2e-3\t5e-2\tTTAACCGG\n"
        ),
        encoding="utf-8",
    )

    df = mod.load_fimo_hits(str(infile))
    assert df.shape[0] == 2
    assert set(df["sequence_name"].tolist()) == {"geneC", "geneD"}
    # legacy format has no alt id; should be mirrored from motif_id.
    assert df.loc[df["sequence_name"] == "geneC", "motif_alt_id"].iloc[0] == "MA0123.1"


def write_asr_fixture(path):
    path.write_text(
        "branch_id\tparent\tnode_class\tname\tnum_intron\tis_imputed\tp_intron_present\tp_intron_absent\n"
        "6\t2\tleaf\tD\tNA\tTrue\t0.1\t0.9\n"
        "0\t-1\troot\tn99\tNA\tFalse\t0.9\t0.1\n"
        "2\t0\tintnode\tn12\tNA\tFalse\t0.2\t0.8\n"
        "1\t0\tintnode\tn11\tNA\tFalse\t0.7\t0.3\n"
        "3\t1\tleaf\tA\t3\tFalse\t1\t0\n"
        "4\t1\tleaf\tB\t0\tFalse\t0\t1\n"
        "5\t2\tleaf\tC\t4\tFalse\t1\t0\n",
        encoding="utf-8",
    )


def test_load_asr_intron_branch_table_translates_native_ids_and_preserves_counts(tmp_path):
    mod = load_module()
    dated_tree = tmp_path / "dated.nwk"
    asr_intron = tmp_path / "asr.tsv"
    dated_tree.write_text("((A:1,B:1)n11:1,(C:1,D:1)n12:1)n99;\n", encoding="utf-8")
    write_asr_fixture(asr_intron)
    df_asr = mod.load_asr_intron_branch_table(str(asr_intron), str(dated_tree))
    # Different child ordering and internal labels must still join by clade ID.
    rooted = mod._ensure_branch_ids(mod.new_tree("((D:1,C:1)n2:1,(B:1,A:1)n1:1)n0;", format=1))
    df_branch = mod.pandas.DataFrame({
        "branch_id": [mod._get_node_label(node) for node in rooted.traverse()],
        "node_name": [node.name for node in rooted.traverse()],
    })
    merged = mod.pandas.merge(df_branch, df_asr.drop(columns=["node_name"]), on="branch_id", validate="one_to_one")
    assert len(merged) == 7
    by_name = merged.set_index("node_name")
    assert by_name.loc["n0", "intron_present"] == 0.9
    assert by_name.loc["n1", "intron_present"] == 0.7
    assert by_name.loc["n2", "intron_present"] == 0.2
    assert by_name.loc[["n0", "n1", "n2", "D"], "num_intron"].isna().all()
    assert by_name.loc["A", "num_intron"] == 3
    assert by_name.loc["B", "num_intron"] == 0
    assert by_name.loc["C", "num_intron"] == 4
    assert by_name.loc["D", "intron_is_imputed"]


@pytest.mark.parametrize("column,value,message", [
    ("branch_id", "99", "branch IDs"),
    ("name", "wrong", "does not match"),
    ("parent", "-1", "does not match"),
    ("p_intron_present", "1.2", "probabilities"),
    ("is_imputed", "unknown", "boolean"),
    ("num_intron", "-1", "non-negative integer"),
    ("num_intron", "1.5", "non-negative integer"),
    ("num_intron", "inf", "non-negative integer"),
    ("num_intron", "2", "imputation flags"),
    ("is_imputed", "False", "imputation flags"),
])
def test_load_asr_intron_branch_table_rejects_mismatched_or_invalid_results(tmp_path, column, value, message):
    mod = load_module()
    tree = tmp_path / "dated.nwk"
    tree.write_text("((A:1,B:1)n11:1,(C:1,D:1)n12:1)n99;", encoding="utf-8")
    table = tmp_path / "asr.tsv"
    write_asr_fixture(table)
    df = pandas.read_csv(table, sep="\t", dtype=str, keep_default_na=False)
    df.loc[0, column] = value
    df.to_csv(table, sep="\t", index=False)
    with pytest.raises(ValueError, match=message):
        mod.load_asr_intron_branch_table(str(table), str(tree))


def test_asr_node_validation_retains_input_row_error_precedence(tmp_path):
    mod = load_module()
    tree = tmp_path / 'dated.nwk'
    tree.write_text('((A:1,B:1)n11:1,(C:1,D:1)n12:1)n99;')
    table = tmp_path / 'asr.tsv'
    write_asr_fixture(table)
    frame = pandas.read_csv(table, sep='\t', dtype=str, keep_default_na=False)
    frame.loc[0, 'parent'] = '-1'
    frame.loc[1, 'name'] = 'incorrect'
    frame.loc[0, 'p_intron_present'] = '2'
    frame['additional_text'] = 'NA'
    frame.to_csv(table, sep='\t', index=False)
    with pytest.raises(ValueError) as error:
        mod.load_asr_intron_branch_table(str(table), str(tree))
    assert str(error.value) == 'Intron ASR node 6 does not match the dated tree.'


def test_asr_loading_retains_inferred_csv_index_and_parent_alignment(tmp_path):
    mod = load_module()
    tree = tmp_path / 'dated.nwk'
    tree.write_text('((A:1,B:1)n11:1,(C:1,D:1)n12:1)n99;')
    table = tmp_path / 'asr.tsv'
    write_asr_fixture(table)
    expected = mod.load_asr_intron_branch_table(str(table), str(tree))
    lines = table.read_text().splitlines()
    indices = [101 + 3 * i for i in range(len(lines) - 1)]
    table.write_text(lines[0] + '\n' + '\n'.join(
        f'{index}\t{line}' for index, line in zip(indices, lines[1:], strict=True)) + '\n')
    result = mod.load_asr_intron_branch_table(str(table), str(tree))
    assert result.index.tolist() == indices
    pandas.testing.assert_frame_equal(result.reset_index(drop=True), expected)
    # Duplicate inferred index labels must retain the rejection from parent
    # lookup, rather than silently admitting the ambiguous table.
    table.write_text(lines[0] + '\n' + '\n'.join('same\t' + line for line in lines[1:]) + '\n')
    with pytest.raises(ValueError, match='truth value of a Series is ambiguous'):
        mod.load_asr_intron_branch_table(str(table), str(tree))


def test_flatten_trait_variable_stats_builds_tree_info_keys():
    mod = load_module()
    df = mod.pandas.DataFrame(
        [
            {"trait": "salt", "variable": "omega", "coef": 1.5, "pvalue": 0.01},
            {"trait": "cold", "variable": "omega", "coef": -0.5, "pvalue": 0.20},
        ]
    )

    out = mod.flatten_trait_variable_stats(df, "pgls_geneTree_")

    assert out == {
        "pgls_geneTree_coef_salt_omega": 1.5,
        "pgls_geneTree_coef_cold_omega": -0.5,
        "pgls_geneTree_pvalue_salt_omega": 0.01,
        "pgls_geneTree_pvalue_cold_omega": 0.20,
    }


def test_load_synteny_summary_aggregates_per_node_metrics(tmp_path):
    mod = load_module()
    infile = tmp_path / "synteny.tsv"
    infile.write_text(
        (
            "node_name\tspecies\tdirection\toffset\tneighbor_gene\tgroup_id\tgroup_size\n"
            "geneA\tsp1\tupstream\t-1\tn1\tG1\t5\n"
            "geneA\tsp1\tdownstream\t2\tn2\tG2\t2\n"
            "geneB\tsp2\tupstream\t-3\tn3\tG1\t5\n"
            "geneB\tsp2\tdownstream\t1\tn4\tG3\t1\n"
            # duplicate row should not inflate counts
            "geneA\tsp1\tupstream\t-1\tn1\tG1\t5\n"
        ),
        encoding="utf-8",
    )

    out = mod.load_synteny_summary(str(infile))
    assert set(out["node_name"].tolist()) == {"geneA", "geneB"}
    row_a = out.loc[out["node_name"] == "geneA"].iloc[0]
    row_b = out.loc[out["node_name"] == "geneB"].iloc[0]

    assert row_a["synteny_edge_count"] == 2
    assert row_a["synteny_group_count"] == 2
    assert row_a["synteny_shared_group_count"] == 1
    assert row_a["synteny_shared_edge_count"] == 1
    assert row_a["synteny_shared_upstream_edge_count"] == 1
    assert row_a["synteny_shared_downstream_edge_count"] == 0
    assert row_a["synteny_min_abs_offset_shared"] == 1
    assert row_a["synteny_max_group_size"] == 5
    assert row_a["synteny_max_shared_group_tip_count"] == 2
    assert row_a["synteny_support_score"] == 0.5

    assert row_b["synteny_edge_count"] == 2
    assert row_b["synteny_group_count"] == 2
    assert row_b["synteny_shared_group_count"] == 1
    assert row_b["synteny_shared_edge_count"] == 1
    assert row_b["synteny_shared_upstream_edge_count"] == 1
    assert row_b["synteny_shared_downstream_edge_count"] == 0
    assert row_b["synteny_min_abs_offset_shared"] == 3
    assert row_b["synteny_max_group_size"] == 5
    assert row_b["synteny_max_shared_group_tip_count"] == 2
    assert row_b["synteny_support_score"] == 0.5


def test_load_synteny_summary_returns_empty_for_missing_required_columns(tmp_path):
    mod = load_module()
    infile = tmp_path / "synteny_missing.tsv"
    infile.write_text(
        (
            "node_name\toffset\tneighbor_gene\n"
            "geneA\t-1\tn1\n"
        ),
        encoding="utf-8",
    )
    out = mod.load_synteny_summary(str(infile))
    assert out.empty


@pytest.mark.parametrize("size_mode", ["absent", "all_missing", "partial_missing"])
def test_load_synteny_summary_preserves_order_missing_sizes_and_duplicate_first(tmp_path, size_mode):
    mod = load_module()
    infile = tmp_path / "synteny_mixed.tsv"
    columns = ["node_name", "group_id", "offset", "direction"]
    # The duplicate rows have large sizes which must not replace the first row.
    rows = [
        [" Z ", "P", "3", "?", "bad"],
        ["A", "S", "-2", " DOWNSTREAM ", "5"],
        ["Z", "S", "-1", "?", "bad"],
        ["A", "S", "1", "UPSTREAM", "7"],
        ["Z", "S", "-1", "upstream", "99"],
        ["B", "Q", "-2", "", "bad"],
        ["Z", "P", "3", "downstream", "100"],
    ]
    if size_mode == "absent":
        rows = [row[:4] for row in rows]
    else:
        columns.append("group_size")
        if size_mode == "all_missing":
            rows = [row[:4] + ["NA"] for row in rows]
    infile.write_text("\t".join(columns) + "\n" + "\n".join("\t".join(row) for row in rows) + "\n")

    out = mod.load_synteny_summary(str(infile))

    max_sizes = [0.0, 7.0, 0.0] if size_mode == "partial_missing" else [2.0, 2.0, 1.0]
    expected = pandas.DataFrame({
        "node_name": ["Z", "A", "B"],
        "synteny_edge_count": [2, 2, 1],
        "synteny_group_count": [2, 1, 1],
        "synteny_shared_group_count": [1, 1, 0],
        "synteny_shared_edge_count": [1, 2, 0],
        "synteny_shared_upstream_edge_count": [1, 1, 0],
        "synteny_shared_downstream_edge_count": [0, 1, 0],
        "synteny_min_abs_offset_shared": [1.0, 1.0, float("nan")],
        "synteny_max_group_size": max_sizes,
        "synteny_max_shared_group_tip_count": [2.0, 2.0, float("nan")],
        "synteny_support_score": [0.5, 1.0, 0.0],
    })
    pandas.testing.assert_frame_equal(out, expected, check_exact=True)


def test_load_synteny_summary_preserves_empty_shared_aggregate_types(tmp_path):
    mod = load_module()
    infile = tmp_path / "synteny_private.tsv"
    infile.write_text("node_name\tgroup_id\toffset\nZ\tP\t-2\nA\tQ\t1\n")

    out = mod.load_synteny_summary(str(infile))

    expected = pandas.DataFrame({
        "node_name": ["Z", "A"],
        "synteny_edge_count": [1, 1],
        "synteny_group_count": [1, 1],
        "synteny_shared_group_count": [0, 0],
        "synteny_shared_edge_count": [0, 0],
        "synteny_shared_upstream_edge_count": [0, 0],
        "synteny_shared_downstream_edge_count": [0, 0],
        "synteny_min_abs_offset_shared": [float("nan"), float("nan")],
        "synteny_max_group_size": [1.0, 1.0],
        "synteny_max_shared_group_tip_count": [float("nan"), float("nan")],
        "synteny_support_score": [0.0, 0.0],
    })
    pandas.testing.assert_frame_equal(out, expected, check_exact=True)


def test_load_synteny_summary_preserves_empty_schema_for_only_missing_node_names(tmp_path):
    mod = load_module()
    infile = tmp_path / "synteny_missing_names.tsv"
    infile.write_text("node_name\tgroup_id\toffset\nNA\tG\t-1\nNULL\tG\t2\n")

    out = mod.load_synteny_summary(str(infile))

    expected = pandas.DataFrame(columns=[
        "node_name",
        "synteny_edge_count",
        "synteny_group_count",
        "synteny_shared_group_count",
        "synteny_shared_edge_count",
        "synteny_shared_upstream_edge_count",
        "synteny_shared_downstream_edge_count",
        "synteny_min_abs_offset_shared",
        "synteny_max_group_size",
        "synteny_max_shared_group_tip_count",
        "synteny_support_score",
    ])
    pandas.testing.assert_frame_equal(out, expected, check_exact=True)


@pytest.mark.parametrize("row,column,value,message", [
    (1, "num_intron", "1", "imputation flags"),
    (1, "is_imputed", "True", "imputation flags"),
    (4, "p_intron_present", "0", "observed counts"),
])
def test_load_asr_intron_rejects_internally_inconsistent_observations(tmp_path, row, column, value, message):
    mod = load_module()
    tree = tmp_path / "dated.nwk"
    tree.write_text("((A:1,B:1)n11:1,(C:1,D:1)n12:1)n99;")
    table = tmp_path / "asr.tsv"
    write_asr_fixture(table)
    df = pandas.read_csv(table, sep="\t", dtype=str, keep_default_na=False)
    df.loc[row, column] = value
    if column == "p_intron_present":
        df.loc[row, "p_intron_absent"] = "1"
    df.to_csv(table, sep="\t", index=False)
    with pytest.raises(ValueError, match=message):
        mod.load_asr_intron_branch_table(str(table), str(tree))


@pytest.mark.parametrize("compressed", [False, True])
def test_promoter_coordinates_join_to_branch_annotations(tmp_path, compressed):
    import gzip

    rooted_tree = tmp_path / "rooted.nwk"
    rooted_tree.write_text("(terminal:1,isolated:1,short:1,noNs:1,empty:1)n0;\n")
    promoter = tmp_path / ("promoters.fa.gz" if compressed else "promoters.fa")
    text = ">terminal\nACNNnN\n>isolated\nNAnCNN\n>short\nN\n>noNs\nACG\n>empty\n"
    if compressed:
        with gzip.open(promoter, "wt") as handle:
            handle.write(text)
    else:
        promoter.write_text(text)
    stub_root = tmp_path / "stub_packages"
    write_kftools_stub(stub_root)
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(filter(None, (str(stub_root), env.get("PYTHONPATH", ""))))
    proc = subprocess.run([
        sys.executable, str(SCRIPT_PATH), "--rooted_tree", str(rooted_tree),
        "--promoter_fasta", str(promoter),
    ], cwd=tmp_path, capture_output=True, text=True, env=env)
    assert proc.returncode == 0, proc.stderr
    branches = pandas.read_csv(tmp_path / "orthogroup.branch.tsv", sep="\t", keep_default_na=False).set_index("node_name")
    expected = {"terminal": "3:6", "isolated": "1,3,5:6", "short": "1,1:6", "noNs": "3:6", "empty": "0:6"}
    assert branches.loc[list(expected), "promoter_N"].to_dict() == expected
    assert branches.loc[list(expected), "promoter_available"].tolist() == ["Y"] * len(expected)


@pytest.mark.parametrize('internal', ['', 'spAB'])
@pytest.mark.parametrize('root', ['', 'spRoot'])
@pytest.mark.parametrize('dated', [False, True])
def test_species_mapping_accepts_unnamed_internal_nodes_without_renaming_tree(tmp_path, internal, root, dated):
    species = tmp_path / 'species.nwk'
    original = f'((Species_a:1,Species_b:1){internal}:1,Species_c:2){root};\n'
    species.write_text(original)
    gene = tmp_path / 'gene.nwk'
    gene.write_text('((Species_a_g1:1,Species_b_g1:1)geneAB:1,Species_c_g1:2)geneRoot;\n')
    command = [sys.executable, str(SCRIPT_PATH), '--rooted_tree', str(gene), '--species_tree', str(species)]
    if dated:
        command += ['--dated_tree', str(gene)]
    stub_root = tmp_path / 'stub_packages'
    write_kftools_stub(stub_root)
    env = os.environ.copy()
    env['PYTHONPATH'] = os.pathsep.join(filter(None, (str(stub_root), env.get('PYTHONPATH', ''))))
    proc = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True, env=env)
    assert proc.returncode == 0, proc.stderr
    table = pandas.read_csv(tmp_path / 'orthogroup.branch.tsv', sep='\t', keep_default_na=False).set_index('node_name')
    assert table.loc['geneAB', 'spnode_coverage'] == internal
    assert table.loc['geneRoot', 'spnode_coverage'] == root
    assert table.loc['Species_a_g1', 'spnode_coverage'] == 'Species_a'
    if dated:
        assert table.loc['geneAB', 'spnode_age'] == internal
        assert table.loc['geneRoot', 'spnode_age'] == root
    assert species.read_text() == original

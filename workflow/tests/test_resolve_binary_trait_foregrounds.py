import subprocess
import sys
from importlib.util import module_from_spec, spec_from_file_location
from io import StringIO
from pathlib import Path

import pytest
from Bio import Phylo

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
SCRIPT_PATH = SUPPORT_DIR / "resolve_binary_trait_foregrounds.py"


def load_module():
    sys.path.insert(0, str(SUPPORT_DIR))
    try:
        spec = spec_from_file_location("resolve_binary_trait_foregrounds", SCRIPT_PATH)
        module = module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path.pop(0)


def parse_tree(newick):
    return Phylo.read(StringIO(newick), "newick")


def trait_rows():
    return [
        ["Alpha_one", "1", "2"],
        ["Beta_two", "1", "2"],
        ["Gamma_three", "0", "0"],
        ["Delta_four", "1", "7"],
        ["Epsilon_five", "0", "0"],
    ]


@pytest.mark.parametrize(
    "newick",
    [
        "(((Alpha_one,Beta_two),Gamma_three),(Delta_four,Epsilon_five));",
        "((Epsilon_five,Delta_four),(Gamma_three,(Beta_two,Alpha_one)));",
    ],
)
def test_resolver_assigns_stable_ids_to_maximal_foreground_clades(newick):
    mod = load_module()
    rows = trait_rows()

    resolved = mod.resolve_binary_foregrounds(
        ["species", "binary_trait", "manual_trait"],
        rows,
        parse_tree(newick),
    )

    assert resolved == [("binary_trait", 2)]
    assert [row[1] for row in rows] == ["1", "1", "0", "2", "0"]
    assert [row[2] for row in rows] == ["2", "2", "0", "7", "0"]


def test_resolver_leaves_nonbinary_and_nonnumeric_traits_unchanged():
    mod = load_module()
    rows = [
        ["Alpha_one", "1", "foreground"],
        ["Beta_two", "2", "background"],
    ]

    resolved = mod.resolve_binary_foregrounds(
        ["species", "lineage_id", "label"],
        rows,
        parse_tree("(Alpha_one,Beta_two);"),
    )

    assert resolved == []
    assert rows == [
        ["Alpha_one", "1", "foreground"],
        ["Beta_two", "2", "background"],
    ]


def test_resolver_rejects_foreground_species_absent_from_tree():
    mod = load_module()

    with pytest.raises(ValueError, match="foreground species absent from the species tree: Beta_two"):
        mod.resolve_binary_foregrounds(
            ["species", "binary_trait"],
            [["Alpha_one", "0"], ["Beta_two", "1"]],
            parse_tree("(Alpha_one,Gamma_three);"),
        )


def test_resolver_treats_unlisted_tree_tips_as_background():
    mod = load_module()
    rows = [["Alpha_one", "1"], ["Beta_two", "1"], ["Delta_four", "1"]]

    resolved = mod.resolve_binary_foregrounds(
        ["species", "binary_trait"],
        rows,
        parse_tree("(((Alpha_one,Beta_two),Gamma_three),Delta_four);"),
    )

    assert resolved == [("binary_trait", 2)]
    assert [row[1] for row in rows] == ["1", "1", "2"]


def test_resolver_prefers_an_exact_tree_tip_over_a_shortened_species_label():
    mod = load_module()
    rows = [["Arabidopsis_thaliana_Col_0", "1"], ["Arabidopsis_lyrata", "0"]]

    resolved = mod.resolve_binary_foregrounds(
        ["species", "binary_trait"],
        rows,
        parse_tree("(Arabidopsis_thaliana_Col_0,Arabidopsis_lyrata);"),
    )

    assert resolved == [("binary_trait", 1)]
    assert [row[1] for row in rows] == ["1", "0"]


def test_cli_writes_resolved_trait_table(tmp_path):
    trait_path = tmp_path / "species_trait.tsv"
    tree_path = tmp_path / "species_tree.nwk"
    output_path = tmp_path / "resolved.tsv"
    trait_path.write_text(
        "species\tbinary_trait\nAlpha_one\t1\nBeta_two\t1\nGamma_three\t0\nDelta_four\t1\n",
        encoding="utf-8",
    )
    tree_path.write_text("(((Alpha_one,Beta_two),Gamma_three),Delta_four);\n", encoding="utf-8")

    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--input",
            str(trait_path),
            "--species-tree",
            str(tree_path),
            "--output",
            str(output_path),
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert "Resolved binary foreground trait binary_trait: 2 lineage(s)." in completed.stdout
    assert output_path.read_text(encoding="utf-8").splitlines() == [
        "species\tbinary_trait",
        "Alpha_one\t1",
        "Beta_two\t1",
        "Gamma_three\t0",
        "Delta_four\t2",
    ]

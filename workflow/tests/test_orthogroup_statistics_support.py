import sys
import types
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

MODULE_PATH = Path(__file__).resolve().parents[1] / "support" / "orthogroup_statistics.py"
TARGET_MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "support"
    / "prepare_generax_ufboot_target.py"
)


def load_module():
    if "kftools" not in sys.modules:
        package = types.ModuleType("kftools")
        package.__path__ = []
        sys.modules["kftools"] = package
    if "kftools.kfog" not in sys.modules:
        sys.modules["kftools.kfog"] = types.ModuleType("kftools.kfog")
    spec = spec_from_file_location("orthogroup_statistics_support_test", MODULE_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_target_module():
    spec = spec_from_file_location("prepare_generax_ufboot_target_test", TARGET_MODULE_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def add_branch_ids(tree):
    for branch_id, node in enumerate(tree.traverse()):
        node.add_prop("branch_id", branch_id)
    return tree


def test_gff_join_rejects_duplicate_gene_rows_before_multiplying_branches(tmp_path):
    module = load_module()
    path = tmp_path / "gff.tsv"
    path.write_text("gene_id\tfeature_size\tnum_intron\nGene1\t6\t0\nGene1\t9\t1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="unique before branch join"):
        module.load_gff_gene_traits(path)
    path.write_text("gene_id\tfeature_size\tnum_intron\nGene1\t6\t0\n", encoding="utf-8")
    traits = module.load_gff_gene_traits(path)
    assert traits.columns.tolist() == ["node_name", "intron_feature_size"]
    assert traits.iloc[0]["node_name"] == "Gene1"


def test_prepared_generax_target_is_unrooted_and_preserves_tips_and_tree_length(
    tmp_path: Path,
):
    mod = load_module()
    target_mod = load_target_module()
    source = tmp_path / "generax.nwk"
    output = tmp_path / "generax.unrooted.nwk"
    source.write_text(
        "((a:1,b:2)n1:3,(c:4,(d:5,e:6)n2:7)n3:8)n0;\n",
        encoding="utf-8",
    )

    target_mod.prepare_target(source, output)
    tree = mod.new_tree(str(output), format=1)

    assert len(tree.get_children()) == 3
    assert {leaf.name for leaf in tree.leaves()} == {"a", "b", "c", "d", "e"}
    assert sum(float(node.dist) for node in tree.traverse() if not node.is_root) == pytest.approx(36)
    assert all(not str(node.name or "") for node in tree.traverse() if not node.is_leaf)


def test_maps_support_by_unrooted_split_and_preserves_legitimate_one_percent():
    mod = load_module()
    rooted = add_branch_ids(
        mod.new_tree("(((a:1,b:1):1,(c:1,d:1):1):1,(e:1,f:1):1);", format=1)
    )
    # A different root representation of the same three internal splits.
    support = mod.new_unrooted_tree(
        "((a:1,b:1)1:1,(c:1,d:1)77:1,(e:1,f:1)95:1);"
    )

    mapped, diagnostics = mod.map_internal_support_by_split(
        rooted, support, support_max=100
    )

    assert diagnostics == {
        "internal_split_count": 3,
        "supported_split_count": 3,
        "support_min": 1.0,
        "support_max": 95.0,
        "all_support_100": False,
    }
    assert 1.0 in mapped.values()
    assert 77.0 in mapped.values()
    assert 95.0 in mapped.values()
    assert rooted.props["branch_id"] not in mapped


def test_rejects_incompatible_support_topology_instead_of_silently_dropping_values():
    mod = load_module()
    rooted = add_branch_ids(
        mod.new_tree("(((a,b),(c,d)),(e,f));", format=1)
    )
    incompatible = mod.new_unrooted_tree(
        "((a,c)80,(b,d)90,(e,f)95);"
    )

    with pytest.raises(ValueError, match="incompatible unrooted topologies"):
        mod.map_internal_support_by_split(rooted, incompatible, support_max=100)


def test_rejects_partially_labelled_support_tree():
    mod = load_module()
    rooted = add_branch_ids(
        mod.new_tree("(((a,b),(c,d)),(e,f));", format=1)
    )
    partially_labelled = mod.new_unrooted_tree(
        "((a,b)80,(c,d),(e,f)95);"
    )

    with pytest.raises(ValueError, match="only a subset"):
        mod.map_internal_support_by_split(rooted, partially_labelled, support_max=100)


def test_explicit_support_mapping_rejects_fully_unlabelled_internal_splits():
    mod = load_module()
    rooted = add_branch_ids(
        mod.new_tree("(((a,b),(c,d)),(e,f));", format=1)
    )
    unlabelled = mod.new_unrooted_tree(
        "((a,b),(c,d),(e,f));"
    )

    with pytest.raises(ValueError, match="no explicit support labels"):
        mod.map_internal_support_by_split(
            rooted,
            unlabelled,
            support_max=100,
            require_support=True,
        )


def test_reports_all_one_hundred_without_rejecting_a_valid_distribution():
    mod = load_module()
    rooted = add_branch_ids(
        mod.new_tree("(((a,b),(c,d)),(e,f));", format=1)
    )
    support = mod.new_unrooted_tree(
        "((a,b)100,(c,d)100,(e,f)100);"
    )

    _mapped, diagnostics = mod.map_internal_support_by_split(
        rooted, support, support_max=100
    )

    assert diagnostics["all_support_100"] is True

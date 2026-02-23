from workflow.support.ete_compat import (
    backend_name,
    get_leaf_names,
    iter_ancestors,
    iter_descendants,
    iter_leaves,
    new_tree,
    node_is_leaf,
    node_is_root,
    write_tree,
)


def test_new_tree_supports_ete3_style_format_and_newick_kwargs():
    tree = new_tree(newick='("a a":1,b:1)root;', format=1, quoted_node_names=True)

    assert backend_name() == "ete4"
    assert node_is_root(tree)
    assert sorted(get_leaf_names(tree)) == ["a a", "b"]


def test_node_helpers_and_iterators_behave_consistently():
    tree = new_tree("((A:1,B:1)X:1,C:1)R;", format=1)
    node_x = [n for n in tree.traverse() if n.name == "X"][0]
    node_a = [n for n in tree.traverse() if n.name == "A"][0]

    assert not node_is_leaf(node_x)
    assert node_is_leaf(node_a)
    assert [n.name for n in iter_ancestors(node_a)] == ["X", "R"]
    assert [n.name for n in iter_descendants(node_x)] == ["A", "B"]
    assert sorted([n.name for n in iter_leaves(tree)]) == ["A", "B", "C"]


def test_write_tree_accepts_format_alias(tmp_path):
    tree = new_tree("(A:1,B:1)R;", format=1)
    outfile = tmp_path / "out.nwk"

    write_tree(tree, outfile=str(outfile), format=5)

    assert outfile.exists()
    tree_roundtrip = new_tree(str(outfile), format=1)
    assert sorted(get_leaf_names(tree_roundtrip)) == ["A", "B"]

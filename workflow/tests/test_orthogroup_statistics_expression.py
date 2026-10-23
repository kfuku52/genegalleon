"""Exact numeric and missing-tip contracts for nested expression annotation."""

import ast
import warnings
from pathlib import Path

import ete4
import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "support/orthogroup_statistics.py"


@pytest.fixture(scope="module")
def implementation():
    parsed = ast.parse(SCRIPT.read_text(), filename=str(SCRIPT))
    main = next(node for node in parsed.body if isinstance(node, ast.FunctionDef) and node.name == "main")
    functions = [node for node in parsed.body
                 if isinstance(node, ast.FunctionDef) and node.name == "node_is_leaf"]
    functions.extend(node for node in main.body if isinstance(node, ast.FunctionDef) and node.name in {
        "get_node_prop", "get_node_label", "compute_min_expression_corr_by_node",
    })
    namespace = {"numpy": np}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(SCRIPT), "exec"), namespace)
    return namespace["compute_min_expression_corr_by_node"], namespace["get_node_label"]


def make_tree(newick="((a,b)x,(c,d)y)root;", labels=None):
    tree = ete4.PhyloTree(newick, parser=1)
    for index, node in enumerate(tree.traverse()):
        node.add_prop("branch_id", index if labels is None else labels[index])
    return tree


def reference(tree, locations, matrix, label):
    # The prior algorithm deliberately indexes descendants by branch label.
    # Keeping it here also covers existing duplicate-label overwrite semantics.
    result, descendant_locations = {}, {}
    for node in tree.traverse(strategy="postorder"):
        key = label(node)
        if node.is_leaf:
            index = locations.get(node.name)
            descendant_locations[key] = None if index is None else [index]
            continue
        indices = []
        for child in node.get_children():
            selected = descendant_locations.get(label(child))
            if selected is None:
                descendant_locations[key] = None
                break
            indices.extend(selected)
        else:
            descendant_locations[key] = indices or None
            if indices:
                result[key] = matrix[np.ix_(indices, indices)].min()
    return result


def assert_exact(implementation, tree, locations, matrix):
    function, label = implementation

    def observe(call):
        result, failure = None, None
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always")
            try:
                result = call()
            except Exception as error:
                failure = error
        return result, failure, [(item.category, str(item.message)) for item in captured]

    expected, expected_error, expected_warnings = observe(lambda: reference(tree, locations, matrix, label))
    actual, actual_error, actual_warnings = observe(lambda: function(tree, locations, matrix))
    assert actual_warnings == expected_warnings
    if expected_error is not None:
        assert type(actual_error) is type(expected_error)
        assert str(actual_error) == str(expected_error)
        return
    assert actual_error is None
    assert [(type(key), key) for key in actual] == [(type(key), key) for key in expected]
    for key, value in expected.items():
        observed = actual[key]
        assert type(observed) is type(value)
        if value is np.ma.masked:
            assert observed is np.ma.masked
        elif isinstance(value, np.generic):
            assert observed.dtype == value.dtype
            assert observed.tobytes() == value.tobytes()
        else:
            assert observed == value
    return actual


@pytest.mark.parametrize("newick", [
    "((a,b)x,(c,d)y)root;",
    "(a,(b,(c,d)x)y)root;",
    "(a,b,c,d)root;",
    "(((a,b)x)y,(c,d)z)root;",
    "((((a,b)x)y)z)root;",
    "a;",
])
@pytest.mark.parametrize("missing", [(), ("a",), ("c",), ("a", "b", "c", "d")])
def test_every_clade_matches_original_bits_and_missing_tip_semantics(implementation, newick, missing):
    tree = make_tree(newick)
    locations = {name: index for index, name in enumerate(("a", "b", "c", "d")) if name not in missing}
    matrix = np.corrcoef(np.random.default_rng(1942).normal(size=(4, 13)))
    assert_exact(implementation, tree, locations, matrix)


def test_directed_cross_child_minima_and_permuted_duplicate_indices(implementation):
    tree = make_tree()
    # Opposite directions deliberately differ; a symmetric shortcut loses -0.9.
    matrix = np.array([[1, .5, .3, -.8], [-.9, 1, .2, .4],
                       [.7, .2, 1, .6], [.1, .2, .3, 1]], dtype=np.float64)
    for locations in ({"a": 2, "b": 0, "c": 3, "d": 1},
                      {"a": 2, "b": 0, "c": 2, "d": 1}):
        assert_exact(implementation, tree, locations, matrix)


@pytest.mark.parametrize("special", [
    np.array([0x7ff8000000000042], dtype=np.uint64).view(np.float64)[0],
    np.array([0xfff8000000000081], dtype=np.uint64).view(np.float64)[0],
    0.0, -0.0, np.inf, -np.inf,
])
@pytest.mark.parametrize("position", [(0, 0), (0, 1), (1, 0), (0, 2), (3, 1)])
def test_nan_payloads_signed_zero_and_infinities_keep_original_reduction(implementation, special, position):
    matrix = np.full((4, 4), 0.8)
    matrix[position] = special
    if special == 0:
        matrix[2, 3] = -special
    assert_exact(implementation, make_tree(), dict(zip("abcd", range(4), strict=True)), matrix)


@pytest.mark.parametrize("kind", ["float32", "int64", "complex128", "object", "fortran", "strided", "subclass", "masked"])
def test_other_matrix_types_and_layouts_keep_existing_contract(implementation, kind):
    matrix = np.arange(1, 17, dtype=np.float64).reshape(4, 4)
    if kind in {"float32", "int64", "complex128", "object"}:
        matrix = matrix.astype(kind)
    elif kind == "fortran":
        matrix = np.asfortranarray(matrix)
    elif kind == "strided":
        matrix = np.arange(1, 65, dtype=np.float64).reshape(8, 8)[::2, ::2]
    elif kind == "subclass":
        class CorrelationArray(np.ndarray):
            pass
        matrix = matrix.view(CorrelationArray)
    else:
        matrix = np.ma.array(matrix, mask=np.eye(4, dtype=bool))
    assert_exact(implementation, make_tree(), dict(zip("abcd", range(4), strict=True)), matrix)


@pytest.mark.parametrize("indices", [
    (-1, 0, 1, 2), (np.int64(0), 1, 2, 3), (True, False, True, False),
    (4, 1, 2, 3), (0.5, 1, 2, 3), (np.array([0]), 1, 2, 3),
])
def test_unusual_indices_preserve_values_or_exact_errors(implementation, indices):
    assert_exact(implementation, make_tree(), dict(zip("abcd", indices, strict=True)), np.full((4, 4), .4))


@pytest.mark.parametrize("matrix", [np.array(.4), np.array([.4, .8]), np.full((4, 3), .4)])
def test_unusual_shapes_preserve_original_error_behavior(implementation, matrix):
    assert_exact(implementation, make_tree(), dict(zip("abcd", range(4), strict=True)), matrix)


@pytest.mark.parametrize("labels", [
    ["NA", "001", "1", "NULL", "nan", "あ", "root"],
    [0, 1, 1, 0, 0, 1, 2],
    [0, 1, 2, 3, 4, 5, None],
])
def test_literal_and_duplicate_branch_labels_keep_existing_dictionary_behavior(implementation, labels):
    assert_exact(implementation, make_tree(labels=labels), dict(zip("abcd", range(4), strict=True)), np.full((4, 4), .4))


def test_literal_leaf_names_and_unused_expression_rows(implementation):
    tree = make_tree("((001,NA)x,(NULL,あ)y)root;")
    locations = {"unused": 0, "001": 4, "NA": 2, "NULL": 3, "あ": 1}
    matrix = np.corrcoef(np.random.default_rng(14).normal(size=(5, 9)))
    assert_exact(implementation, tree, locations, matrix)


def test_empty_catalog_and_none_index_skip_clade_output(implementation):
    tree = make_tree()
    assert assert_exact(implementation, tree, {}, np.empty((0, 0))) == {}
    assert_exact(implementation, tree, {"a": None, "b": 0, "c": 1, "d": 2}, np.full((3, 3), .4))


def large_tree(count=128):
    def subtree(start, stop):
        if stop - start == 1:
            return f"g{start}"
        middle = (start + stop) // 2
        return f"({subtree(start, middle)},{subtree(middle, stop)})n{start}_{stop}"
    return make_tree(subtree(0, count) + ";")


@pytest.mark.parametrize("bits", [
    0x7ff0000000000042, 0xfff0000000000081,
    0x7ff8000000000042, 0xfff8000000000081,
    0x0000000000000000, 0x8000000000000000,
    0x7ff0000000000000, 0xfff0000000000000,
])
@pytest.mark.parametrize("position", [(0, 0), (0, 1), (0, 64), (64, 0)])
@pytest.mark.parametrize("invalid", ["warn", "raise"])
def test_large_clades_preserve_special_bits_warnings_and_errors(implementation, bits, position, invalid):
    matrix = np.full((128, 128), .8)
    matrix[position] = np.array([bits], dtype=np.uint64).view(np.float64)[0]
    with np.errstate(invalid=invalid):
        assert_exact(implementation, large_tree(), {f"g{i}": i for i in range(128)}, matrix)


def test_large_clades_preserve_lower_finite_minimum_despite_other_special_values(implementation):
    matrix = np.full((128, 128), .8)
    matrix[0, 64], matrix[64, 0] = -.9, -.8
    matrix[2, 3], matrix[3, 2], matrix[4, 4] = 0., -0., np.inf
    matrix[5, 6] = np.nextafter(0., -1.)
    assert_exact(implementation, large_tree(), {f"g{i}": i for i in range(128)}, matrix)


@pytest.mark.parametrize("kind", ["symmetric", "asymmetric", "missing", "duplicate_labels", "duplicate_indices", "unary"])
def test_large_clades_match_all_node_bits_on_nonconstant_matrices(implementation, kind):
    rng = np.random.default_rng(2026)
    matrix = np.corrcoef(rng.normal(size=(128, 19)))
    tree = large_tree()
    locations = {f"g{i}": i for i in range(128)}
    if kind == "asymmetric":
        matrix = rng.uniform(-1., 1., size=(128, 128))
    elif kind == "missing":
        locations.pop("g17")
    elif kind == "duplicate_labels":
        for index, node in enumerate(tree.traverse()):
            node.add_prop("branch_id", index % 19)
    elif kind == "duplicate_indices":
        locations["g3"] = locations["g1"]
    elif kind == "unary":
        parent = ete4.PhyloTree()
        parent.add_prop("branch_id", 10000)
        parent.add_child(tree)
        tree = parent
    assert_exact(implementation, tree, locations, matrix)

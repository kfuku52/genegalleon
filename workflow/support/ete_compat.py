#!/usr/bin/env python3
# coding: utf-8

"""Compatibility helpers for ETE toolkit (ETE4).

These helpers provide call-level compatibility for common ETE3-style APIs.
"""

from functools import lru_cache


@lru_cache(maxsize=1)
def _load_backend():
    try:
        from ete4 import NCBITaxa, Tree

        return Tree, NCBITaxa, "ete4"
    except Exception as exc:
        raise ImportError("ete4 is required but unavailable. Install ete4 in conda base.") from exc


def new_tree(*args, **kwargs):
    tree_cls, _, _ = _load_backend()
    if ("newick" in kwargs) and (len(args) == 0):
        args = (kwargs.pop("newick"),)
    elif "newick" in kwargs:
        kwargs.pop("newick")
    if ("format" in kwargs) and ("parser" not in kwargs):
        kwargs["parser"] = kwargs.pop("format")
    kwargs.pop("quoted_node_names", None)
    return tree_cls(*args, **kwargs)


def ncbi_taxa(*args, **kwargs):
    _, ncbi_cls, _ = _load_backend()
    return ncbi_cls(*args, **kwargs)


def write_tree(tree, *args, **kwargs):
    if ("format" in kwargs) and ("parser" not in kwargs):
        kwargs["parser"] = kwargs.pop("format")
    return tree.write(*args, **kwargs)


def node_is_leaf(node):
    is_leaf_attr = getattr(node, "is_leaf")
    if callable(is_leaf_attr):
        return bool(is_leaf_attr())
    return bool(is_leaf_attr)


def node_is_root(node):
    is_root_attr = getattr(node, "is_root")
    if callable(is_root_attr):
        return bool(is_root_attr())
    return bool(is_root_attr)


def iter_leaves(node):
    if hasattr(node, "iter_leaves"):
        return node.iter_leaves()
    if hasattr(node, "leaves"):
        return node.leaves()
    return []


def iter_descendants(node):
    if hasattr(node, "iter_descendants"):
        return node.iter_descendants()
    if hasattr(node, "descendants"):
        return node.descendants()
    return []


def iter_ancestors(node):
    if hasattr(node, "iter_ancestors"):
        return node.iter_ancestors()
    if hasattr(node, "ancestors"):
        return node.ancestors()
    return []


def get_leaf_names(tree):
    if hasattr(tree, "get_leaf_names"):
        return tree.get_leaf_names()
    leaf_names_attr = getattr(tree, "leaf_names")
    return list(leaf_names_attr() if callable(leaf_names_attr) else leaf_names_attr)


def get_common_ancestor(tree, leaves):
    if hasattr(tree, "get_common_ancestor"):
        return tree.get_common_ancestor(leaves)
    if hasattr(tree, "common_ancestor"):
        return tree.common_ancestor(leaves)
    raise AttributeError("Tree object does not provide a common ancestor API.")


def backend_name():
    try:
        _, _, backend = _load_backend()
    except ImportError:
        return "unavailable"
    return backend

#!/usr/bin/env python3
# coding: utf-8

import argparse
import copy
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import json
import os
import re
import sys

import numpy
import pandas
import ete4
from kftools.kfog import *
from kftools.kfphylo import *


def new_tree(newick_or_path, format=1, quoted_node_names=False):
    _ = quoted_node_names
    if isinstance(newick_or_path, str) and os.path.exists(newick_or_path):
        with open(newick_or_path, 'r', encoding='utf-8') as handle:
            newick_or_path = handle.read().strip()
    return ete4.PhyloTree(newick_or_path, parser=format)


def iter_ancestors(node):
    return node.ancestors()


def iter_descendants(node):
    return node.descendants()


def iter_leaves(node):
    return (n for n in node.traverse() if n.is_leaf)


def node_is_leaf(node):
    return bool(node.is_leaf)


def node_is_root(node):
    return bool(node.is_root)


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--species_tree', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--unaligned_aln', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--trimal_aln', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--unrooted_tree', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--iqtree_model', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--rooted_tree', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--rooting_log', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--generax_nhx', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--notung_root_log', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--notung_reconcil_stats', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--dated_tree', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--dated_log', metavar='PATH', default='', type=str, help='')
    
    parser.add_argument('--targetp', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--rpsblast', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--uniprot', metavar='PATH', default='', type=str, help='')
    
    parser.add_argument('--l1ou_tree', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--l1ou_regime', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--l1ou_leaf', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--phylogeneticem_tree', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--phylogeneticem_regime', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--phylogeneticem_leaf', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--expression', metavar='PATH', default='', type=str, help='')
    
    parser.add_argument('--mapdnds_tree_dn', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--mapdnds_tree_ds', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--codeml_tsv', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--hyphy_dnds_json', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--hyphy_relax_json', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--hyphy_relax_reversed_json', metavar='PATH', default='', type=str, help='')
    
    parser.add_argument('--character_gff', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--scm_intron', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--fimo', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--promoter_fasta', metavar='PATH', default='', type=str, help='')
    
    parser.add_argument('--csubst_b', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--csubst_cb_stats', metavar='PATH', default='', type=str, help='')
    
    parser.add_argument('--gene_pgls_stats', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--species_pgls_stats', metavar='PATH', default='', type=str, help='')
    
    parser.add_argument('--clade_ortholog_prefix', metavar='STR', default='', type=str, help='')
    parser.add_argument('--ncpu', metavar='INT', default=1, type=int, help='Number of worker threads.')
    return parser


def load_fimo_hits(path):
    base_cols = ['motif_id', 'motif_alt_id', 'sequence_name', 'start', 'stop', 'strand', 'q-value']
    if (not path) or (not os.path.exists(path)) or (os.path.getsize(path) == 0):
        return pandas.DataFrame(columns=base_cols)
    try:
        df = pandas.read_csv(path, sep='\t', header=0, index_col=None, dtype=str)
    except Exception as exc:
        print('Failed to parse FIMO table {}: {}'.format(path, exc))
        return pandas.DataFrame(columns=base_cols)
    if df.empty:
        return pandas.DataFrame(columns=base_cols)

    colmap = {}
    for c in df.columns:
        c_norm = str(c).strip()
        if c_norm.startswith('#'):
            c_norm = c_norm.lstrip('#').strip()
        c_norm = c_norm.lower().replace(' ', '_').replace('-', '_')
        colmap[c] = c_norm
    df = df.rename(columns=colmap)

    def pick_col(candidates):
        for c in candidates:
            if c in df.columns:
                return c
        return None

    motif_col = pick_col(['motif_id', 'pattern_name'])
    motif_alt_col = pick_col(['motif_alt_id', 'motif_altid', 'alt_id'])
    seq_col = pick_col(['sequence_name'])
    start_col = pick_col(['start'])
    stop_col = pick_col(['stop', 'end'])
    strand_col = pick_col(['strand'])
    qvalue_col = pick_col(['q_value', 'qvalue'])
    if qvalue_col is None:
        # Legacy FIMO output can miss q-values when disabled; fall back to p-values.
        qvalue_col = pick_col(['p_value', 'pvalue'])

    required_missing = []
    if motif_col is None:
        required_missing.append('motif_id/pattern_name')
    if seq_col is None:
        required_missing.append('sequence_name')
    if start_col is None:
        required_missing.append('start')
    if stop_col is None:
        required_missing.append('stop/end')
    if qvalue_col is None:
        required_missing.append('q-value/p-value')
    if required_missing:
        print(
            'FIMO table is missing required columns ({}): {}'.format(
                ', '.join(required_missing), path
            )
        )
        return pandas.DataFrame(columns=base_cols)

    df_out = pandas.DataFrame(index=df.index.copy())
    df_out['motif_id'] = df[motif_col].astype(str).str.strip()
    if motif_alt_col is None:
        df_out['motif_alt_id'] = df_out['motif_id']
    else:
        motif_alt = df[motif_alt_col].astype(str).str.strip()
        is_missing_alt = motif_alt.str.lower().isin(['', 'na', 'nan', '.'])
        motif_alt = motif_alt.where(~is_missing_alt, df_out['motif_id'])
        df_out['motif_alt_id'] = motif_alt
    df_out['sequence_name'] = df[seq_col].astype(str).str.strip()
    df_out['start'] = pandas.to_numeric(df[start_col], errors='coerce')
    df_out['stop'] = pandas.to_numeric(df[stop_col], errors='coerce')
    if strand_col is None:
        df_out['strand'] = '.'
    else:
        df_out['strand'] = df[strand_col].astype(str).str.strip()
    df_out['q-value'] = pandas.to_numeric(df[qvalue_col], errors='coerce')

    keep = (
        (df_out['sequence_name'] != '') &
        df_out['start'].notna() &
        df_out['stop'].notna() &
        df_out['q-value'].notna() &
        (df_out['motif_id'] != '')
    )
    df_out = df_out.loc[keep, base_cols].copy()
    if df_out.empty:
        return df_out
    df_out['start'] = df_out['start'].astype(int)
    df_out['stop'] = df_out['stop'].astype(int)
    return df_out


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    args.ncpu = max(1, int(args.ncpu))
    params = vars(args).copy()
    
    print('Command: {}'.format(' '.join(sys.argv)), flush=True)
    print('\n')

    def open_text(path):
        if path.endswith('.gz'):
            return gzip.open(path, 'rt')
        return open(path, 'r')

    def add_dict_key_prefix(d, prefix):
        return {f'{prefix}_{k}': v for k, v in d.items()}

    def get_node_prop(node, key, default=None):
        if hasattr(node, key):
            return getattr(node, key)
        if hasattr(node, "props"):
            return node.props.get(key, default)
        return default

    def set_node_prop(node, key, value):
        if hasattr(node, "add_prop"):
            node.add_prop(key, value)
            return
        try:
            setattr(node, key, value)
            return
        except (AttributeError, TypeError):
            pass
        if hasattr(node, "props"):
            node.props[key] = value

    def get_node_label(node):
        label = get_node_prop(node, "branch_id", None)
        if label is None:
            raise AttributeError("Node branch_id is unavailable.")
        try:
            return int(label)
        except (TypeError, ValueError):
            return label

    def ensure_branch_ids(tree):
        nodes = list(tree.traverse())
        if len(nodes) == 0:
            return tree
        all_leaf_names = sorted([leaf.name for leaf in iter_leaves(tree)])
        if len(all_leaf_names) == 0:
            for i, node in enumerate(nodes):
                set_node_prop(node, "branch_id", i)
            return tree
        leaf_branch_ids = {leaf_name: (1 << i) for i, leaf_name in enumerate(all_leaf_names)}
        node_label_sum = {}
        for node in tree.traverse(strategy="postorder"):
            if node_is_leaf(node):
                node_label_sum[node] = leaf_branch_ids[node.name]
            else:
                node_mask = 0
                for child in node.get_children():
                    node_mask |= node_label_sum[child]
                node_label_sum[node] = node_mask
        branch_ids = [node_label_sum[node] for node in nodes]
        argsort_labels = numpy.argsort(branch_ids)
        label_ranks = numpy.empty_like(argsort_labels)
        label_ranks[argsort_labels] = numpy.arange(len(argsort_labels))
        for i, node in enumerate(nodes):
            set_node_prop(node, "branch_id", int(label_ranks[i]))
        return tree

    def get_leaf_names_compat(node):
        if hasattr(node, "get_leaf_names"):
            return list(node.get_leaf_names())
        if hasattr(node, "leaf_names"):
            leaf_names_attr = getattr(node, "leaf_names")
            return list(leaf_names_attr() if callable(leaf_names_attr) else leaf_names_attr)
        return [leaf.name for leaf in iter_leaves(node)]

    def get_common_ancestor_compat(tree, leaves):
        if hasattr(tree, "get_common_ancestor"):
            return tree.get_common_ancestor(leaves)
        if hasattr(tree, "common_ancestor"):
            return tree.common_ancestor(leaves)
        raise AttributeError("Tree object does not provide a common ancestor API.")

    def transfer_root(tree_to, tree_from, verbose=False):
        def clear_root_branch_property_compat(tree, prop_name):
            # ETE4 rejects reroot when the root branch itself has branch properties.
            # Keep behavior local to the root only; do not modify descendant branches.
            if hasattr(tree, "props"):
                tree.props.pop(prop_name, None)
            try:
                setattr(tree, prop_name, None)
            except Exception:
                pass

        def set_outgroup_compat(tree, outgroup):
            try:
                tree.set_outgroup(outgroup)
            except AssertionError as exc:
                msg = str(exc)
                # IQ-TREE can keep different support labels on both root children in unrooted trees.
                # ETE rejects rerooting in that case, so normalize only root-adjacent supports and retry.
                if 'inconsistent support at the root' in msg:
                    for child in tree.get_children():
                        child.support = 1.0
                    tree.set_outgroup(outgroup)
                    return
                # Newer ETE4 may also reject rerooting when the current root branch has support.
                if 'root has branch property: support' in msg:
                    clear_root_branch_property_compat(tree, 'support')
                    tree.set_outgroup(outgroup)
                    return
                raise

        tip_set_diff = set(get_leaf_names_compat(tree_to)) - set(get_leaf_names_compat(tree_from))
        if tip_set_diff:
            raise Exception('tree_to has more tips than tree_from. tip_set_diff = ' + str(tip_set_diff))
        subroot_leaves = [get_leaf_names_compat(node) for node in tree_from.get_children()]
        is_n0_bigger_than_n1 = len(subroot_leaves[0]) > len(subroot_leaves[1])
        ingroups = subroot_leaves[0] if is_n0_bigger_than_n1 else subroot_leaves[1]
        outgroups = subroot_leaves[0] if not is_n0_bigger_than_n1 else subroot_leaves[1]
        if verbose:
            print('outgroups:', outgroups)
        set_outgroup_compat(tree_to, ingroups[0])
        if len(outgroups) == 1:
            outgroup_ancestor = next(node for node in iter_leaves(tree_to) if node.name == outgroups[0])
        else:
            outgroup_ancestor = get_common_ancestor_compat(tree_to, outgroups)
        set_outgroup_compat(tree_to, outgroup_ancestor)
        subroot_to = tree_to.get_children()
        subroot_from = tree_from.get_children()
        total_subroot_length_to = sum(node.dist for node in subroot_to)
        total_subroot_length_from = sum(node.dist for node in subroot_from)
        dist_by_leafset = {
            frozenset(get_leaf_names_compat(node)): node.dist
            for node in subroot_from
        }
        for n_to in subroot_to:
            n_to_leafset = frozenset(get_leaf_names_compat(n_to))
            n_from_dist = dist_by_leafset.get(n_to_leafset)
            if n_from_dist is not None:
                n_to.dist = (n_from_dist / total_subroot_length_from) * total_subroot_length_to
        for n_to in tree_to.traverse():
            if n_to.name == '':
                n_to.name = tree_to.name
                tree_to.name = 'Root'
                break
        return tree_to

    def nwk2table(tree, attr='', age=False, parent=False, sister=False):
        tree_format = 0 if attr == 'support' else 1
        if isinstance(tree, str):
            tree = new_tree(tree, format=tree_format)
        tree = ensure_branch_ids(tree)
        nodes = list(tree.traverse())
        age_by_label = {}
        if ((attr == 'dist') and age):
            for node in tree.traverse(strategy='postorder'):
                nlabel = get_node_label(node)
                if node_is_leaf(node):
                    age_by_label[nlabel] = 0.0
                else:
                    first_child = node.get_children()[0]
                    age_by_label[nlabel] = age_by_label[get_node_label(first_child)] + first_child.dist
        rows = []
        for node in nodes:
            nlabel = get_node_label(node)
            if attr == '':
                attr_value = numpy.nan
            else:
                attr_value = get_node_prop(node, attr, numpy.nan)
            row = {
                'branch_id': nlabel,
                attr: attr_value,
            }
            if ((attr == 'dist') and age):
                row['age'] = age_by_label.get(nlabel, numpy.nan)
            if parent:
                if node_is_root(node):
                    row['parent'] = -1
                else:
                    row['parent'] = get_node_label(node.up)
            if sister:
                if node_is_root(node):
                    row['sister'] = -1
                else:
                    sisters = node.get_sisters()
                    row['sister'] = get_node_label(sisters[0]) if len(sisters) > 0 else -1
            rows.append(row)
        df = pandas.DataFrame.from_records(rows)
        return df.sort_values(by='branch_id').reset_index(drop=True)

    def node_gene2species(gene_tree, species_tree, is_ultrametric=False):
        gene_tree2 = copy.deepcopy(gene_tree)
        gene_tree2 = ensure_branch_ids(gene_tree2)
        for leaf in iter_leaves(gene_tree2):
            leaf_name_split = leaf.name.split("_", 2)
            if len(leaf_name_split) >= 2:
                leaf.name = leaf_name_split[0] + "_" + leaf_name_split[1]
        tip_set_diff = set(get_leaf_names_compat(gene_tree2)) - set(get_leaf_names_compat(species_tree))
        if tip_set_diff:
            sys.stderr.write(f"Warning. A total of {len(tip_set_diff)} species are missing in the species tree: {str(tip_set_diff)}\n")
        cn = ["branch_id", "spnode_coverage"]
        if is_ultrametric:
            cn.append("spnode_age")
        species_nodes = list(species_tree.traverse(strategy="postorder"))
        species_names = {sn: sn.name.replace('\'', '') for sn in species_nodes}
        species_leaf_node = {leaf.name: leaf for leaf in iter_leaves(species_tree)}
        species_depth = {}
        for sn in species_tree.traverse(strategy="preorder"):
            species_depth[sn] = 0 if node_is_root(sn) else (species_depth[sn.up] + 1)
        if is_ultrametric:
            species_age = {}
            for sn in species_nodes:
                if node_is_leaf(sn):
                    species_age[sn] = 0.0
                else:
                    first_child = sn.get_children()[0]
                    species_age[sn] = species_age[first_child] + first_child.dist
            species_up_age = {}
            for sn in species_nodes:
                if node_is_root(sn):
                    species_up_age[sn] = numpy.inf
                else:
                    species_up_age[sn] = species_age[sn.up]
        lca_cache = {}

        def pair_lca(node_a, node_b):
            if node_a is node_b:
                return node_a
            key = (node_a, node_b) if id(node_a) <= id(node_b) else (node_b, node_a)
            cached = lca_cache.get(key)
            if cached is not None:
                return cached
            a = node_a
            b = node_b
            depth_a = species_depth[a]
            depth_b = species_depth[b]
            while depth_a > depth_b:
                a = a.up
                depth_a -= 1
            while depth_b > depth_a:
                b = b.up
                depth_b -= 1
            while a is not b:
                a = a.up
                b = b.up
            lca_cache[key] = a
            return a

        gene_nodes = list(gene_tree2.traverse(strategy="postorder"))
        gene_coverage = {}
        gene_has_missing_species = {}
        if is_ultrametric:
            gene_age = {}
        for gn in gene_nodes:
            if node_is_leaf(gn):
                covered_species_node = species_leaf_node.get(gn.name)
                gene_coverage[gn] = covered_species_node
                gene_has_missing_species[gn] = covered_species_node is None
                if is_ultrametric:
                    gene_age[gn] = 0.0
                continue
            children = gn.get_children()
            if is_ultrametric:
                first_child = children[0]
                gene_age[gn] = gene_age[first_child] + first_child.dist
            has_missing_species = any(gene_has_missing_species[child] for child in children)
            gene_has_missing_species[gn] = has_missing_species
            if has_missing_species:
                gene_coverage[gn] = None
                continue
            covered_species_node = gene_coverage[children[0]]
            for child in children[1:]:
                covered_species_node = pair_lca(covered_species_node, gene_coverage[child])
            gene_coverage[gn] = covered_species_node

        rows = []
        for gn in gene_nodes:
            coverage_node = gene_coverage[gn]
            row = {
                "branch_id": get_node_label(gn),
                "spnode_coverage": "" if coverage_node is None else species_names[coverage_node],
            }
            if is_ultrametric:
                row["spnode_age"] = ""
                if coverage_node is not None:
                    gn_age = gene_age[gn]
                    current_species_node = coverage_node
                    while current_species_node is not None:
                        if (gn_age >= species_age[current_species_node]) and (gn_age < species_up_age[current_species_node]):
                            row["spnode_age"] = species_names[current_species_node]
                            break
                        current_species_node = current_species_node.up
            rows.append(row)
        return pandas.DataFrame(rows, columns=cn)

    def get_misc_node_statistics(tree_file, tax_annot=False):
        tree = new_tree(tree_file, format=1)
        tree = ensure_branch_ids(tree)
        cn1 = ["branch_id", "taxon", "taxid", "num_sp", "num_leaf", "so_event", "dup_conf_score"]
        cn2 = ["parent", "sister", "child1", "child2", "so_event_parent"]
        cn = cn1 + cn2
        nodes = list(tree.traverse())
        if tax_annot:
            raise NotImplementedError("tax_annot=True is unsupported in the ETE4 compatibility path.")
        sci_name_by_node = {}
        taxid_by_node = {}
        for node in nodes:
            taxid_by_node[node] = -999
            if node_is_leaf(node):
                name_split = node.name.split('_', 2)
                if len(name_split) >= 2:
                    sci_name_by_node[node] = name_split[0] + ' ' + name_split[1]
                else:
                    sci_name_by_node[node] = node.name
            else:
                sci_name_by_node[node] = ''
        species_index = {}
        species_mask_by_node = {}
        num_leaf_by_node = {}
        dup_conf_score_by_node = {}
        for node in tree.traverse(strategy="postorder"):
            if node_is_leaf(node):
                species_id = species_index.setdefault(sci_name_by_node[node], len(species_index))
                species_mask_by_node[node] = (1 << species_id)
                num_leaf_by_node[node] = 1
                dup_conf_score_by_node[node] = 0.0
                continue
            children = node.get_children()
            species_mask = 0
            num_leaf = 0
            for child in children:
                species_mask |= species_mask_by_node[child]
                num_leaf += num_leaf_by_node[child]
            species_mask_by_node[node] = species_mask
            num_leaf_by_node[node] = num_leaf
            if len(children) >= 2:
                child1_mask = species_mask_by_node[children[0]]
                child2_mask = species_mask_by_node[children[1]]
                union_mask = child1_mask | child2_mask
                union_count = union_mask.bit_count()
                if union_count == 0:
                    dup_conf_score_by_node[node] = 0.0
                else:
                    dup_conf_score_by_node[node] = (child1_mask & child2_mask).bit_count() / union_count
            else:
                dup_conf_score_by_node[node] = 0.0
        rows = []
        for node in nodes:
            label = get_node_label(node)
            children = node.get_children()
            if not node_is_root(node):
                parent = get_node_label(node.up)
                sisters = node.get_sisters()
                sister = get_node_label(sisters[0]) if len(sisters) > 0 else -999
            else:
                parent = -999
                sister = -999
            child1 = get_node_label(children[0]) if len(children) >= 1 else -999
            child2 = get_node_label(children[1]) if len(children) >= 2 else -999
            node_dup_conf_score = dup_conf_score_by_node.get(node, 0.0)
            if node_is_leaf(node):
                so_event = "L"
            else:
                so_event = "D" if node_dup_conf_score > 0 else "S"
            so_event_parent = "S"
            if (node.up is not None) and (dup_conf_score_by_node.get(node.up, 0.0) > 0):
                so_event_parent = "D"
            rows.append(
                {
                    "branch_id": label,
                    "taxon": str(sci_name_by_node[node]),
                    "taxid": taxid_by_node[node],
                    "num_sp": species_mask_by_node[node].bit_count(),
                    "num_leaf": num_leaf_by_node[node],
                    "so_event": so_event,
                    "dup_conf_score": node_dup_conf_score,
                    "parent": parent,
                    "sister": sister,
                    "child1": child1,
                    "child2": child2,
                    "so_event_parent": so_event_parent,
                }
            )
        return pandas.DataFrame(rows, columns=cn)
    
    def ou2table(regime_file, leaf_file, input_tree_file):
        df_regime = pandas.read_csv(regime_file, sep="\t")
        df_leaf = pandas.read_csv(leaf_file, sep="\t")
        regime_regimes = numpy.array([0, ] + df_regime['regime'].dropna().unique().tolist()).astype(int)
        leaf_regimes = df_leaf['regime'].dropna().unique().astype(int)
        is_same_regime = (set(regime_regimes) == set(leaf_regimes))
        assert is_same_regime, 'Regime numbers did not match between the regime and leaf files.'
        tree = new_tree(input_tree_file, format=1)
        tree = ensure_branch_ids(tree)
        nodes = list(tree.traverse())
        tissues = df_leaf.columns[3:].values
        if ('expectations' in df_leaf['param'].values):
            df_leaf.loc[(df_leaf['param'] == 'expectations'), 'param'] = 'mu'
        cn1 = ["branch_id", "regime", "is_shift", "num_child_shift"]
        cn2 = ["tau", "delta_tau", "delta_maxmu", "mu_complementarity"]
        cn3 = ["mu_" + tissue for tissue in tissues]
        cn = cn1 + cn2 + cn3

        name_to_regime = (
            df_regime.loc[df_regime['node_name'].notna(), ['node_name', 'regime']]
            .drop_duplicates(subset=['node_name'], keep='first')
            .set_index('node_name')['regime']
            .to_dict()
        )
        regime_by_label = {}
        for node in tree.traverse(strategy='preorder'):
            nlabel = get_node_label(node)
            inherited_regime = 0 if node_is_root(node) else regime_by_label[get_node_label(node.up)]
            regime_no = name_to_regime.get(node.name)
            regime_by_label[nlabel] = inherited_regime if regime_no is None else int(regime_no)

        is_mu = (df_leaf['param'] == 'mu')
        regime_cols = [c for c in df_leaf.columns if c not in ['node_name', 'param']]
        df_leaf_unique = df_leaf.loc[is_mu, regime_cols].copy()
        for col in df_leaf_unique.columns:
            if col != 'regime':
                df_leaf_unique[col] = pandas.to_numeric(df_leaf_unique[col], errors='coerce')
        df_leaf_unique = df_leaf_unique.drop_duplicates()
        if df_leaf_unique.shape[0] == 0:
            df_leaf_unique = pandas.DataFrame({'regime': [0]})
        else:
            df_leaf_unique = df_leaf_unique.groupby(by='regime', as_index=False).mean(numeric_only=True)
        if 'regime' not in df_leaf_unique.columns:
            df_leaf_unique['regime'] = 0
        for tissue in tissues:
            if tissue not in df_leaf_unique.columns:
                df_leaf_unique[tissue] = 0.0
        regime_labels = df_leaf_unique['regime'].to_numpy(dtype=int, copy=False)
        regime_mu_values = df_leaf_unique.loc[:, tissues].to_numpy(dtype=float, copy=False)
        regime_to_mu = {int(regime_no): mu_values for regime_no, mu_values in zip(regime_labels, regime_mu_values)}

        records = []
        for node in nodes:
            nlabel = get_node_label(node)
            node_regime = regime_by_label[nlabel]
            parent_regime = regime_by_label[get_node_label(node.up)] if not node_is_root(node) else node_regime
            row = {
                "branch_id": nlabel,
                "regime": node_regime,
                "is_shift": int((not node_is_root(node)) and (node_regime != parent_regime)),
                "num_child_shift": 0,
            }
            mu_values = regime_to_mu.get(int(node_regime))
            if mu_values is None:
                mu_values = numpy.zeros(len(tissues), dtype=float)
            for tissue, mu_value in zip(tissues, mu_values):
                row["mu_" + tissue] = float(mu_value)
            records.append(row)
        df = pandas.DataFrame.from_records(records, columns=cn)
        df[cn2] = df[cn2].astype(float)
        df[cn3] = df[cn3].astype(float)

        df["tau"] = calc_tau(df, cn3, unlog2=True, unPlus1=True)
        df = df.set_index("branch_id", drop=False)
        for node in nodes:
            nlabel = get_node_label(node)
            if not node_is_root(node):
                tau_up = float(df.at[get_node_label(node.up), "tau"])
                tau_my = float(df.at[nlabel, "tau"])
                df.at[nlabel, "delta_tau"] = tau_my - tau_up
                if int(df.at[nlabel, "is_shift"]) == 1:
                    sis_label = get_node_label(node.get_sisters()[0])
                    my_mu = df.loc[nlabel, cn3].to_numpy(dtype=float)
                    sis_mu = df.loc[sis_label, cn3].to_numpy(dtype=float)
                    my_maxmu = my_mu.max()
                    sis_maxmu = sis_mu.max()
                    delta_maxmu = my_maxmu - sis_maxmu
                    df.at[nlabel, "delta_maxmu"] = delta_maxmu
                    my_mu_unlog = numpy.clip(numpy.exp2(my_mu) - 1, a_min=0, a_max=None)
                    sis_mu_unlog = numpy.clip(numpy.exp2(sis_mu) - 1, a_min=0, a_max=None)
                    df.at[nlabel, "mu_complementarity"] = calc_complementarity(my_mu_unlog, sis_mu_unlog)
            if not node_is_leaf(node):
                node_regime = regime_by_label[nlabel]
                child_labels = [get_node_label(child) for child in node.get_children()]
                is_child1_shift = (node_regime != regime_by_label[child_labels[0]])
                is_child2_shift = (node_regime != regime_by_label[child_labels[1]])
                num_child_shift = sum([is_child1_shift, is_child2_shift])
                df.at[nlabel, "num_child_shift"] = num_child_shift
        return df.reset_index(drop=True)
    
    def get_N_coordinate(file):
        with open_text(file) as fh:
            seqs = fh.read().split('>')
        if seqs[0]=='':
            seqs = seqs[1::] # The first item is almost always empty ('')
        if len(seqs) == 0 or all((s.strip() == '') for s in seqs):
            return pandas.DataFrame({'node_name': [], 'promoter_N': [], 'promoter_available': []})
        seqnames = list()
        for i in range(len(seqs)):
            seqnames.append(re.sub('\n.*', '', seqs[i].replace('>','')))
            seqs[i] = re.sub('^[^\n]*\n', '', seqs[i])
            seqs[i] = re.sub('\n', '', seqs[i])
        max_seqlen = max([ len(seq) for seq in seqs ])
        out = list()
        for i in range(len(seqs)):
            N_positions = [nuc.span()[1] for nuc in re.finditer('[Nn]', seqs[i])]
            if len(N_positions)==0:
                N_slices = ['',]
            elif len(N_positions)==1:
                N_slices = [str(N_positions[0]),]
            else:
                N_slices = [str(N_positions[0]),]
                N_extension_flag = False
                for j in range(1, len(N_positions)):
                    if N_positions[j-1]==(N_positions[j]-1):
                        N_extension_flag = True
                        N_extension_until = str(N_positions[j])
                    else:
                        if N_extension_flag:
                            new_slice = N_slices[len(N_slices)-1]+':'+N_extension_until
                            N_slices[len(N_slices)-1] = new_slice
                        N_slices.append(str(N_positions[j]))
                        N_extension_flag = False
                # Flush the final contiguous N-run when it reaches the sequence end.
                if N_extension_flag:
                    new_slice = N_slices[len(N_slices)-1]+':'+N_extension_until
                    N_slices[len(N_slices)-1] = new_slice
            if (len(seqs[i]) < max_seqlen):
                slice_outside_scaffold = str(len(seqs[i]))+':'+str(max_seqlen)
                N_slices.append(slice_outside_scaffold)
            slice_str = re.sub('^,', '', ','.join(N_slices))
            out.append(slice_str)
        df_out = pandas.DataFrame({'node_name':seqnames, 'promoter_N':out})
        df_out.loc[:,'promoter_available'] = 'Y'
        return df_out
    
    def get_ancestor_and_descendant(tree, df_in):
        df_out = df_in.copy()
        tree = ensure_branch_ids(tree)
        ancestors_map = {}
        descendants_map = {}
        for node in tree.traverse():
            anocestor_nodes = [get_node_label(n) for n in iter_ancestors(node)]
            descendant_nodes = [get_node_label(n) for n in iter_descendants(node)]
            nlabel = get_node_label(node)
            ancestors_map[nlabel] = ','.join([str(n) for n in anocestor_nodes])
            descendants_map[nlabel] = ','.join([str(n) for n in descendant_nodes])
        df_out['ancestors'] = df_out['branch_id'].map(ancestors_map).fillna('')
        df_out['descendants'] = df_out['branch_id'].map(descendants_map).fillna('')
        return df_out
    
    def flatten_dict(data, col_name_prefix=''):
        flat_data = []
        for key, value in data.items():
            flat_record = {'id': key}  # or any other name for the key field
            for sub_key, sub_value in value.items():
                if isinstance(sub_value, dict):
                    for sub_sub_key, sub_sub_value in sub_value.items():
                        flat_record[f'{col_name_prefix}{sub_key}_{sub_sub_key}'] = sub_sub_value
                else:
                    flat_record[f'{col_name_prefix}{sub_key}'] = sub_value
            flat_data.append(flat_record)
        return flat_data

    def compute_min_expression_corr_by_node(gene_tree, exp_index_locs, corr_matrix):
        min_cor_by_label = {}
        exp_leaf_locs = {}
        for node in gene_tree.traverse(strategy='postorder'):
            nlabel = get_node_label(node)
            if node_is_leaf(node):
                idx = exp_index_locs.get(node.name)
                exp_leaf_locs[nlabel] = None if idx is None else [idx]
                continue
            child_locs = []
            has_missing = False
            for child in node.get_children():
                locs = exp_leaf_locs.get(get_node_label(child))
                if locs is None:
                    has_missing = True
                    break
                child_locs.extend(locs)
            if has_missing or (len(child_locs) == 0):
                exp_leaf_locs[nlabel] = None
                continue
            exp_leaf_locs[nlabel] = child_locs
            min_cor_by_label[nlabel] = corr_matrix[numpy.ix_(child_locs, child_locs)].min()
        return min_cor_by_label
    
    
    def build_hyphy_relax_branch_table(task):
        trait, json_data = task
        if len(json_data) == 0:  # Empty if no foreground lineage is included in the gene tree.
            return trait, None
        flat_dict = flatten_dict(data=json_data['branch attributes']['0'])
        df_tmp1 = pandas.DataFrame(flat_dict)
        df_tmp1 = df_tmp1.drop(labels=['original name'], axis=1)
        df_tmp1.columns = df_tmp1.columns.str.replace(' ', '_')
        df_tmp1 = df_tmp1.add_prefix('hyphy_relax_')
        df_tmp1 = df_tmp1.add_suffix('_' + trait)
        tree_txt = json_data['input']['trees']['0'] + ';'
        hyphy_tree = new_tree(tree_txt, format=1)
        df_tmp = nwk2table(tree=hyphy_tree, attr='name')
        df_tmp = pandas.merge(df_tmp, df_tmp1, left_on='name', right_on=f'hyphy_relax_id_{trait}', how='left')
        df_tmp = df_tmp.drop(labels=['name', f'hyphy_relax_id_{trait}'], axis=1)
        return trait, df_tmp
    
    
    def get_hyphy_relax_branch_tables(json_data_all_traits, ncpu):
        tasks = [(trait, json_data) for trait, json_data in sorted(json_data_all_traits.items())]
        if len(tasks) == 0:
            return []
        if ncpu > 1 and len(tasks) > 1:
            max_workers = min(ncpu, len(tasks))
            results = {}
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = {executor.submit(build_hyphy_relax_branch_table, task): task[0] for task in tasks}
                for future in as_completed(futures):
                    trait, df_tmp = future.result()
                    results[trait] = df_tmp
            return [(trait, results[trait]) for trait, _ in tasks]
        return [build_hyphy_relax_branch_table(task) for task in tasks]
    
    
    def extract_hyphy_relax_tree_metrics(json_data_all_traits):
        out = {}
        for trait, json_data in json_data_all_traits.items():
            if len(json_data) == 0:  # Empty if no foreground lineage is included in the gene tree.
                continue
            for key in json_data['test results'].keys():
                key_new = key.replace(' ', '_').replace('p-value', 'pvalue')
                out[f'hyphy_relax_{key_new}_{trait}'] = json_data['test results'][key]
        return out

    def merge_hyphy_relax_tables(relax_tables):
        tables = [df_tmp.set_index('branch_id') for _, df_tmp in relax_tables if df_tmp is not None]
        if len(tables) == 0:
            return None
        return pandas.concat(tables, axis=1, join='outer').reset_index()

    def concat_tables_on_key(tables, key):
        prepared = []
        for table in tables:
            if table is None or table.shape[0] == 0:
                continue
            if key not in table.columns:
                continue
            table_tmp = table.drop_duplicates(subset=[key], keep='first').set_index(key)
            prepared.append(table_tmp)
        if len(prepared) == 0:
            return None
        return pandas.concat(prepared, axis=1, join='outer').reset_index()

    def annotate_clade_signatures(tree, collect_leaf_labels=False):
        tree = ensure_branch_ids(tree)
        leaf_labels = {} if collect_leaf_labels else None
        for node in tree.traverse(strategy='postorder'):
            if node_is_leaf(node):
                node.add_prop('clade_sig', frozenset([node.name]))
                if collect_leaf_labels:
                    node.add_prop('_leaf_labels', [node.name])
            else:
                clade_sig = set()
                if collect_leaf_labels:
                    leaf_labels_this_node = []
                for child in node.get_children():
                    clade_sig.update(child.props['clade_sig'])
                    if collect_leaf_labels:
                        leaf_labels_this_node.extend(child.props['_leaf_labels'])
                node.add_prop('clade_sig', frozenset(clade_sig))
                if collect_leaf_labels:
                    node.add_prop('_leaf_labels', leaf_labels_this_node)
            if collect_leaf_labels:
                leaf_labels[get_node_label(node)] = '; '.join(node.props['_leaf_labels'])
        if collect_leaf_labels:
            for node in tree.traverse():
                node.props.pop('_leaf_labels', None)
            return tree, leaf_labels
        return tree
    
    
    hyphy_relax_json_data = None
    hyphy_relax_reversed_json_data = None
    numlabel_merge_tables = []
    node_left_merge_tables = []
    rooted_has_spnode_coverage = False
    needs_generax_parent = False
    
    print('Generating branch table.')
    df_branch = pandas.DataFrame({'branch_id':[0,]})
    if os.path.exists(params["species_tree"]):
        species_tree = new_tree(params['species_tree'], format=1)
    else:
        species_tree = None
    if os.path.exists(params["rooted_tree"]):
        rooted_tree = new_tree(params['rooted_tree'], format=1)
        rooted_tree = ensure_branch_ids(rooted_tree)
        rooted_tree, gene_labels = annotate_clade_signatures(rooted_tree, collect_leaf_labels=True)
        df_tmp_name = nwk2table(tree=rooted_tree, attr='name')
        df_tmp_name.columns = ['branch_id', 'node_name']
        df_tmp_dist = nwk2table(tree=rooted_tree, attr='dist')
        df_tmp_dist.columns = ['branch_id', 'bl_rooted']
        df_tmp_dist.loc[:, 'gene_labels'] = df_tmp_dist['branch_id'].map(gene_labels).fillna('')

        rooted_tables = [df_tmp_name, df_tmp_dist]
        if species_tree is not None:
            df_tmp_species = node_gene2species(rooted_tree, species_tree, is_ultrametric=False)
            rooted_has_spnode_coverage = ('spnode_coverage' in df_tmp_species.columns)
            rooted_tables.append(df_tmp_species)
        df_tmp = get_misc_node_statistics(tree_file=params["rooted_tree"], tax_annot=False)
        df_tmp = get_ancestor_and_descendant(tree=rooted_tree, df_in=df_tmp)
        rooted_tables.append(df_tmp)
        rooted_merged = concat_tables_on_key(rooted_tables, 'branch_id')
        if rooted_merged is not None:
            df_branch = pandas.merge(df_branch, rooted_merged, on='branch_id', how='outer')
    if os.path.exists(params["unrooted_tree"]) and os.path.exists(params["rooted_tree"]):
        # If a rooted tree was reconciled, unrooted tree may have inconsistent branch_id,
        # so branch matching cannot rely on it.
        nlabels = numpy.arange(len(list(rooted_tree.traverse())))
        df_tmp = pandas.DataFrame({'branch_id':nlabels})
        df_tmp.loc[:,'support_unrooted'] = numpy.nan
        df_tmp.loc[:,'bl_unrooted'] = numpy.nan
        unrooted = new_tree(params['unrooted_tree'], format=0)
        unrooted = transfer_root(tree_to=unrooted, tree_from=rooted_tree)
        unrooted = annotate_clade_signatures(unrooted)
        unrooted_by_clade = {unode.props.get('clade_sig'): unode for unode in unrooted.traverse()}
        for rnode in rooted_tree.traverse():
            unode = unrooted_by_clade.get(rnode.props.get('clade_sig'))
            if unode is not None:
                df_tmp.at[get_node_label(rnode), 'bl_unrooted'] = unode.dist
                if unode.support != 1.0:
                    df_tmp.at[get_node_label(rnode), 'support_unrooted'] = unode.support
        subroot_nodes = rooted_tree.get_children()
        subroot_nlabels = [get_node_label(n) for n in subroot_nodes]
        is_subroot_support = ~df_tmp.loc[subroot_nlabels,'support_unrooted'].isna()
        numlabel_merge_tables.append(df_tmp)
    if os.path.exists(params["dated_tree"]):
        df_tmp = nwk2table(tree=params["dated_tree"], attr='dist', age=True)
        df_tmp.columns = ['branch_id', 'bl_dated', 'age']
        numlabel_merge_tables.append(df_tmp)
    if os.path.exists(params["dated_tree"]) and os.path.exists(params["species_tree"]):
        gene_tree = new_tree(params['dated_tree'], format=1)
        df_tmp = node_gene2species(gene_tree, species_tree, is_ultrametric=True)
        if rooted_has_spnode_coverage and ('spnode_coverage' in df_tmp.columns):
            df_tmp = df_tmp.drop(labels='spnode_coverage', axis=1)
        numlabel_merge_tables.append(df_tmp)
    if (os.path.exists(params["generax_nhx"])):
        generax_tree = new_tree(params['generax_nhx'], format=1)
        generax_tree = ensure_branch_ids(generax_tree)
        df_tmp = nwk2table(tree=generax_tree, attr='S')
        for attr in ['H','D']:
            df_tmp1 = nwk2table(tree=generax_tree, attr=attr)
            df_tmp.loc[:,attr] = df_tmp1.loc[:,attr].values
        df_tmp.columns = ['branch_id','spnode_generax','generax_transfer','generax_event']
        is_N = (df_tmp.loc[:,'generax_event']=='N')
        is_Y = (df_tmp.loc[:,'generax_event']=='Y')
        is_H = df_tmp.loc[:,'generax_transfer'].str.startswith('Y')
        df_tmp.loc[is_N,'generax_event'] = 'S'
        df_tmp.loc[is_Y,'generax_event'] = 'D'
        df_tmp.loc[is_H,'generax_event'] = 'H'
        leaf_nlabels = {get_node_label(leaf) for leaf in iter_leaves(generax_tree)}
        df_tmp.loc[df_tmp['branch_id'].isin(leaf_nlabels), 'generax_event'] = 'L'
        numlabel_merge_tables.append(df_tmp)
        needs_generax_parent = True
    if (all([ os.path.exists(params[key]) for key in ['mapdnds_tree_ds','mapdnds_tree_dn'] ])):
        df_tmp_ds = nwk2table(tree=params["mapdnds_tree_ds"], attr='dist')
        df_tmp_ds.columns = ['branch_id', 'mapdnds_ds']
        df_tmp_dn = nwk2table(tree=params["mapdnds_tree_dn"], attr='dist')
        df_tmp_dn.columns = ['branch_id', 'mapdnds_dn']
        df_tmp = pandas.merge(df_tmp_ds, df_tmp_dn, on='branch_id', how='outer')
        df_tmp['mapdnds_omega'] = df_tmp['mapdnds_dn'] / df_tmp['mapdnds_ds']
        numlabel_merge_tables.append(df_tmp)
    if (all([ os.path.exists(params[key]) for key in ['hyphy_dnds_json'] ])):
        with open(params['hyphy_dnds_json'], 'r') as json_file:
            json_data = json.load(json_file)
        flat_dict = flatten_dict(data=json_data['branch attributes']['0'])
        df_tmp1 = pandas.DataFrame(flat_dict)
        df_tmp1 = df_tmp1.drop(labels=['original name',], axis=1)
        df_tmp1.columns = df_tmp1.columns.str.replace(' ', '_')
        df_tmp1.columns = df_tmp1.columns.str.replace('Confidence_Intervals_LB', 'omega_ci_lower_bound')
        df_tmp1.columns = df_tmp1.columns.str.replace('Confidence_Intervals_UB', 'omega_ci_upper_bound')
        df_tmp1.columns = df_tmp1.columns.str.replace('Confidence_Intervals_MLE', 'omega')
        df_tmp1.columns = df_tmp1.columns.str.replace('p-value', 'pvalue').str.replace('P-value', 'pvalue')
        df_tmp1.columns = 'hyphy_dnds_'+df_tmp1.columns
        tree_txt = json_data['input']['trees']['0']+';'
        hyphy_tree = new_tree(tree_txt, format=1)
        df_tmp = nwk2table(tree=hyphy_tree, attr='name')
        df_tmp = pandas.merge(df_tmp, df_tmp1, left_on='name', right_on='hyphy_dnds_id', how='left')
        df_tmp = df_tmp.drop(labels=['name', 'hyphy_dnds_id'], axis=1)
        numlabel_merge_tables.append(df_tmp)
    if (all([ os.path.exists(params[key]) for key in ['hyphy_relax_json'] ])):
        with open(params['hyphy_relax_json'], 'r') as json_file:
            hyphy_relax_json_data = json.load(json_file)
        relax_tables = get_hyphy_relax_branch_tables(hyphy_relax_json_data, params['ncpu'])
        relax_merged = merge_hyphy_relax_tables(relax_tables)
        if relax_merged is not None:
            numlabel_merge_tables.append(relax_merged)
    if (all([ os.path.exists(params[key]) for key in ['hyphy_relax_reversed_json'] ])):
        with open(params['hyphy_relax_reversed_json'], 'r') as json_file:
            hyphy_relax_reversed_json_data = json.load(json_file)
        relax_tables = get_hyphy_relax_branch_tables(hyphy_relax_reversed_json_data, params['ncpu'])
        relax_merged = merge_hyphy_relax_tables(relax_tables)
        if relax_merged is not None:
            numlabel_merge_tables.append(relax_merged)
    if (all([ os.path.exists(params[key]) for key in ['character_gff'] ])):
        df_tmp = pandas.read_csv(params['character_gff'], sep='\t',  header=0, index_col=None)
        df_tmp = df_tmp.drop(labels='num_intron', axis=1)
        df_tmp = df_tmp.rename(columns={'gene_id':'node_name'})
        df_tmp = df_tmp.rename(columns={'feature_size':'intron_feature_size'})
        df_branch = pandas.merge(df_branch, df_tmp, on='node_name', how='left')
    if (all([ os.path.exists(params[key]) for key in ['scm_intron'] ])):
        df_tmp = pandas.read_csv(params['scm_intron'], sep='\t',  header=0, index_col=None)
        df_tmp.columns = [ c.replace('leaf','node_name') for c in df_tmp.columns ]
        df_branch = pandas.merge(df_branch, df_tmp, on='node_name', how='outer')
        df_branch = df_branch.drop('intron_absent', axis=1)
        df_branch = compute_delta(df_branch, 'intron_present')
    if (os.path.exists(params["targetp"])):
        df_tmp = pandas.read_csv(params['targetp'], sep='\t',  header=None, index_col=None, comment='#')
        df_tmp.columns = ['node_name', 'targetp_prediction', 'targetp_noTP', 'targetp_SP', 'targetp_mTP', 'targetp_cTP', 'targetp_luTP', 'targetp_CS_position']
        node_left_merge_tables.append(df_tmp)
    if (os.path.exists(params["uniprot"])):
        df_tmp = pandas.read_csv(params['uniprot'], sep='\t',  header=0, index_col=None)
        df_tmp.index = df_tmp['gene_id']
        df_tmp = df_tmp.drop('gene_id', axis=1)
        df_tmp = df_tmp.reset_index().rename(columns={'index': 'node_name'})
        node_left_merge_tables.append(df_tmp)
    if (os.path.exists(params['rpsblast'])):
        df_tmp = pandas.read_csv(params['rpsblast'], sep='\t',  header=0, index_col=None, dtype = {'stitle':'object'})
        df_tmp = df_tmp.loc[(df_tmp['evalue']<=0.05),:]
        df_tmp = df_tmp.loc[:,['qacc','stitle']]
        df_tmp.loc[:,'stitle'] = df_tmp.loc[:,'stitle'].str.replace(',','|', 1, regex=False)
        df_tmp.loc[:,'stitle'] = df_tmp.loc[:,'stitle'].str.replace(',.*','', 1, regex=True)
        df_tmp.loc[:,'stitle'] = df_tmp.loc[:,'stitle'].str.replace('|',',', 1, regex=False)
        df_tmp.loc[:, 'stitle'] = df_tmp['stitle'].astype(str)
        df_tmp = df_tmp.groupby('qacc', sort=False)['stitle'].agg('; '.join).rename('pfam_domain')
        node_left_merge_tables.append(df_tmp.reset_index().rename(columns={'qacc': 'node_name'}))
    if (os.path.exists(params['fimo'])):
        df_tmp = load_fimo_hits(params['fimo'])
        if not df_tmp.empty:
            join_cols = ['motif_id', 'motif_alt_id', 'start', 'stop', 'strand', 'q-value']
            df_tmp_join = df_tmp.loc[:, ['sequence_name'] + join_cols].copy()
            df_tmp_join = df_tmp_join.astype({col_name: str for col_name in join_cols})
            df_tmp2 = (
                df_tmp_join.groupby('sequence_name', as_index=False, sort=False)[join_cols]
                .agg('; '.join)
                .rename(columns={
                    'motif_id': 'fimo_id',
                    'motif_alt_id': 'fimo_alt_id',
                    'start': 'fimo_start',
                    'stop': 'fimo_end',
                    'strand': 'fimo_strand',
                    'q-value': 'fimo_qvalue',
                })
            )
            df_tmp2['node_name'] = df_tmp2['sequence_name']
            node_left_merge_tables.append(df_tmp2)
    if (os.path.exists(params['promoter_fasta'])):
        df_tmp = get_N_coordinate(file=params['promoter_fasta'])
        node_left_merge_tables.append(df_tmp)
    if os.path.exists(params["l1ou_regime"]) and os.path.exists(params["l1ou_leaf"]):
        print('processing l1ou')
        df_tmp = ou2table(params['l1ou_regime'], params['l1ou_leaf'], params["dated_tree"])
        df_tmp.columns = [ 'l1ou_'+c if c!='branch_id' else c for c in df_tmp.columns ]
        numlabel_merge_tables.append(df_tmp)
    if os.path.exists(params["phylogeneticem_regime"]) and os.path.exists(params["phylogeneticem_leaf"]):
        print('processing PhylogeneticEM')
        df_tmp = ou2table(params['phylogeneticem_regime'], params['phylogeneticem_leaf'], params["dated_tree"])
        df_tmp.columns = [ 'phylogeneticem_'+c if c!='branch_id' else c for c in df_tmp.columns ]
        numlabel_merge_tables.append(df_tmp)
    if os.path.exists(params["expression"]) and os.path.exists(params["rooted_tree"]):
        col = 'clade_min_expression_pearsoncor'
        df_branch.loc[:,col] = numpy.nan
        try:
            df_exp = pandas.read_csv(params['expression'], sep='\t', index_col=0)
        except (FileNotFoundError, OSError, UnicodeDecodeError, ValueError, pandas.errors.EmptyDataError) as exc:
            print('Failed to read {}'.format(params['expression']))
            print('Reason: {}'.format(exc))
            print('Gene tree may not include genes from species with expression data.')
            df_exp = pandas.DataFrame()
        if not df_exp.empty:
            df_cor = numpy.corrcoef(df_exp)
            exp_index_locs = {name: i for i, name in enumerate(df_exp.index)}
            gene_tree = new_tree(params['rooted_tree'], format=1)
            gene_tree = ensure_branch_ids(gene_tree)
            min_cor_by_label = compute_min_expression_corr_by_node(gene_tree, exp_index_locs, df_cor)
            if len(min_cor_by_label) > 0:
                df_branch.loc[:, col] = df_branch['branch_id'].map(min_cor_by_label)
        df_exp.columns = 'expression_'+df_exp.columns
        node_left_merge_tables.append(df_exp.reset_index().rename(columns={'index': 'node_name'}))
    if (os.path.exists(params['csubst_b'])):
        df_tmp = pandas.read_csv(params['csubst_b'], sep='\t',  header=0, index_col=None)
        df_tmp = df_tmp.drop('branch_name', axis=1)
        df_tmp.columns = 'csubst_'+df_tmp.columns
        df_tmp.columns = df_tmp.columns.str.replace('csubst_branch_id','branch_id', regex=False)
        numlabel_merge_tables.append(df_tmp)
    merged_numlabel = concat_tables_on_key(numlabel_merge_tables, 'branch_id')
    if merged_numlabel is not None:
        df_branch = pandas.merge(df_branch, merged_numlabel, on='branch_id', how='outer')
    if needs_generax_parent and ('parent' in df_branch.columns) and ('generax_event' in df_branch.columns):
        parent_event_map = (
            df_branch.loc[:, ['branch_id', 'generax_event']]
            .drop_duplicates(subset=['branch_id'])
            .set_index('branch_id')['generax_event']
        )
        df_branch.loc[:, 'generax_event_parent'] = df_branch['parent'].map(parent_event_map)
    merged_node_left = concat_tables_on_key(node_left_merge_tables, 'node_name')
    if merged_node_left is not None:
        df_branch = pandas.merge(df_branch, merged_node_left, on='node_name', how='left')
    if (params['clade_ortholog_prefix'] != '') and os.path.exists(params["rooted_tree"]):
        newcol = params['clade_ortholog_prefix']+'closest_gene'
        ortholog_ids_by_label = {}
        for node in rooted_tree.traverse(strategy='postorder'):
            nlabel = get_node_label(node)
            if node_is_leaf(node):
                if node.name.startswith(params['clade_ortholog_prefix']):
                    ortholog_ids_by_label[nlabel] = [node.name.replace(params['clade_ortholog_prefix'], ''),]
                else:
                    ortholog_ids_by_label[nlabel] = []
            else:
                ids = []
                for child in node.get_children():
                    ids.extend(ortholog_ids_by_label.get(get_node_label(child), []))
                ortholog_ids_by_label[nlabel] = sorted(set(ids)) if len(ids) > 0 else []

        closest_by_label = {}
        closest_ids_by_label = {}
        for node in rooted_tree.traverse(strategy='preorder'):
            nlabel = get_node_label(node)
            own_ids = ortholog_ids_by_label.get(nlabel, [])
            if len(own_ids) > 0:
                closest_ids = own_ids
            elif node_is_root(node):
                closest_ids = []
            else:
                closest_ids = closest_ids_by_label.get(get_node_label(node.up), [])
            closest_ids_by_label[nlabel] = closest_ids
            if len(closest_ids) > 0:
                closest_by_label[nlabel] = ', '.join(closest_ids)
        df_branch.loc[:, newcol] = df_branch['branch_id'].map(closest_by_label).fillna('')
    if df_branch['node_name'].isnull().any():
        txt = ('Some inputs may have inconsistent tree sizes. '
               'Please double-check all files were generated under the same parameter set throughout the pipeline.')
        raise Exception(txt)
    df_branch.to_csv('orthogroup.branch.tsv', sep='\t', index=False)
    
    def branch2tree(df):
        out = dict()
        out['num_branch'] = df.shape[0]
        out['num_spe'] = (df['so_event'] == 'S').sum()
        out['num_dup'] = (df['so_event'] == 'D').sum()
        out['num_sp'] = df['num_sp'].max()
        if 'age' in df.columns:
            out['tree_height'] = df['age'].max()
        is_bl_col = df.columns.str.startswith('bl_')
        for col in df.columns[is_bl_col]:
            out[col] = df[col].sum()
        return out
    
    print('Generating tree table.')
    tree_info = dict()
    tree_info.update(branch2tree(df_branch))
    if os.path.exists(params['dated_log']):
        tree_info['dating_method'] = get_dating_method(params['dated_log'])
    if os.path.exists(params["unaligned_aln"]):
        tree_tmp = get_aln_stats(params['unaligned_aln'])
        tree_tmp = add_dict_key_prefix(tree_tmp, 'original')
        tree_info.update(tree_tmp)
    if os.path.exists(params["trimal_aln"]):
        tree_tmp = get_aln_stats(params['trimal_aln'])
        tree_tmp = add_dict_key_prefix(tree_tmp, 'cleaned')
        tree_info.update(tree_tmp)
    if os.path.exists(params["rooting_log"]):
        tree_tmp = get_root_stats(params['rooting_log'])
        tree_info.update(tree_tmp)
    if os.path.exists(params["notung_root_log"]):
        tree_tmp = get_notung_root_stats(params['notung_root_log'])
        tree_info.update(tree_tmp)
    if os.path.exists(params["notung_reconcil_stats"]):
        tree_tmp = get_notung_reconcil_stats(params['notung_reconcil_stats'])
        tree_info.update(tree_tmp)
    if os.path.exists(params["iqtree_model"]):
        tree_tmp = get_iqtree_model_stats(params['iqtree_model'])
        tree_info.update(tree_tmp)
    if os.path.exists(params["codeml_tsv"]):
        tmp = pandas.read_csv(params['codeml_tsv'], sep='\t', header=0, index_col=None)
        for col in tmp.columns:
            tree_info['codeml_'+col] = tmp.loc[:,col].values[0]
    for method in ['l1ou','phylogeneticem']:
        if os.path.exists(params[method+"_tree"]):
            tmp = pandas.read_csv(params[method+'_tree'], sep='\t')
            num_shift = tmp['num_shift'].values[0]
            tree_tmp = {method+'_num_shift':num_shift,}
            tree_info.update(tree_tmp)
        if os.path.exists(params[method+"_regime"]):
            try:
                tree_tmp = regime2tree(params[method+'_regime'])
            except ValueError as exc:
                print(
                    "Skipping {} regime summary due to invalid regime parameters: {}".format(
                        method, exc
                    ),
                    flush=True,
                )
            else:
                tree_tmp = add_dict_key_prefix(tree_tmp, method)
                tree_info.update(tree_tmp)
    if os.path.exists(params['csubst_cb_stats']):
        df_tmp = pandas.read_csv(params['csubst_cb_stats'], sep='\t',  header=0, index_col=None)
        df_tmp = df_tmp.loc[(df_tmp['arity']==2),:]
        for col in df_tmp.columns:
            tree_info['csubst_'+col] = df_tmp.loc[:,col].values[0]
    if os.path.exists(params['gene_pgls_stats']):
        df_tmp = pandas.read_csv(params['gene_pgls_stats'], sep='\t',  header=0, index_col=None)
        is_all_na = df_tmp.isnull().all().all()
        if is_all_na:
            print('Results of gene tree PGLS are all NA.', flush=True)
        else:
            df_tmp = df_tmp.loc[~df_tmp.isna().all(axis=1),:]
            stat_cols = df_tmp.columns.values.tolist()
            stat_cols.remove('trait')
            stat_cols.remove('variable')
            df_tmp['trait_var'] = '_' + df_tmp['trait'] + '_' + df_tmp['variable']
            for i in range(len(df_tmp)):
                for stat_col in stat_cols:
                    key = 'pgls_geneTree_'+stat_col+df_tmp['trait_var'][i]
                    tree_info[key] = df_tmp[stat_col][i]
    if os.path.exists(params['species_pgls_stats']):
        df_tmp = pandas.read_csv(params['species_pgls_stats'], sep='\t',  header=0, index_col=None)
        df_tmp = df_tmp.loc[~df_tmp.isna().all(axis=1),:]
        stat_cols = df_tmp.columns.values.tolist()
        stat_cols.remove('trait')
        stat_cols.remove('variable')
        df_tmp['trait_var'] = '_' + df_tmp['trait'] + '_' + df_tmp['variable']
        for i in range(len(df_tmp)):
            for stat_col in stat_cols:
                key = 'pgls_speciesTree_'+stat_col+df_tmp['trait_var'][i]
                tree_info[key] = df_tmp[stat_col][i]
    if (all([ os.path.exists(params[key]) for key in ['hyphy_relax_json'] ])):
        if hyphy_relax_json_data is None:
            with open(params['hyphy_relax_json'], 'r') as json_file:
                hyphy_relax_json_data = json.load(json_file)
        tree_info.update(extract_hyphy_relax_tree_metrics(hyphy_relax_json_data))
    if (all([ os.path.exists(params[key]) for key in ['hyphy_relax_reversed_json'] ])):
        if hyphy_relax_reversed_json_data is None:
            with open(params['hyphy_relax_reversed_json'], 'r') as json_file:
                hyphy_relax_reversed_json_data = json.load(json_file)
        tree_info.update(extract_hyphy_relax_tree_metrics(hyphy_relax_reversed_json_data))
    
    df_tree = pandas.DataFrame(tree_info, index=[0, ])
    df_tree.to_csv('orthogroup.tree.tsv', sep='\t', index=False)
    
    print('orthogroup_statistics: done!')


if __name__ == '__main__':
    main()

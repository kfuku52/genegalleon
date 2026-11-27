"""Separate species-branch identity from Newick support labels for ortholog plots."""

import argparse
from collections import Counter
from functools import reduce
from pathlib import Path

from Bio import Phylo
from reconciliation_database import branch_id


class SpeciesBranchMapping:
    def __init__(self, tree_path):
        self.tree = Phylo.read(tree_path, "newick")
        ordered = list(self.tree.find_clades(order="level"))
        tips = [node.name for node in ordered if node.is_terminal()]
        if not all(tips) or len(tips) != len(set(tips)):
            raise ValueError("Species mapping requires unique, nonempty tip names")
        counts = Counter(node.name for node in ordered if node.name)
        self.parents = {child: parent for parent in ordered for child in parent.clades}
        self.depths = {self.tree.root: 0}
        self.tips = {node.name: node for node in ordered if node.is_terminal()}
        self._lcas = {}
        for node in ordered[1:]:
            self.depths[node] = self.depths[self.parents[node]] + 1

        self.requires_coverage = False
        self.explicit_names = set(self.tips)
        for node in ordered:
            if node.is_terminal():
                continue
            if node.name and counts[node.name] == 1:
                self.explicit_names.add(node.name)
            else:
                node.name = branch_id([tip.name for tip in node.get_terminals()])
                self.requires_coverage = True
            # This tree carries mapping keys. The separate plot tree carries support.
            node.confidence = None
        if len({node.name for node in ordered}) != len(ordered):
            raise ValueError("Species mapping tree contains repeated clades or colliding names")

    def write(self, path):
        Phylo.write(self.tree, path, "newick", format_branch_length="%.17g")

    def common_ancestor(self, left, right):
        if left is None or right is None:
            return None
        if left is right:
            return left
        key = frozenset((left, right))
        if key not in self._lcas:
            a, b = left, right
            while a is not b:
                if self.depths[a] >= self.depths[b]:
                    a = self.parents[a]
                else:
                    b = self.parents[b]
            self._lcas[key] = a
        return self._lcas[key]

    def coverage_by_branch(self, rows, children, root):
        """Project complete descendant species onto their LCA in one bottom-up pass."""
        order = [root]
        for node in order:
            order.extend(children.get(node, ()))
        observed = {}
        for node in reversed(order):
            descendants = children.get(node)
            if descendants:
                observed[node] = reduce(self.common_ancestor, (observed[child] for child in descendants))
            else:
                observed[node] = self.tips.get(str(rows[node].get("spnode_coverage") or "").strip())
        return observed

    def use_coverage(self, tree_rows, source_rows, coverage):
        if not self.requires_coverage:
            return
        for row in tree_rows:
            node = int(row["node_id"])
            # Preserve explicit reconciled placement for duplication events.
            # Legacy coverage values, including numeric support and blanks, are
            # re-derived from the complete gene clade rather than displayed tips.
            reconciled = str(source_rows[node].get("spnode_generax") or "").strip()
            if row["event"] == "D" and reconciled and reconciled.lower() != "nan":
                if reconciled not in self.explicit_names:
                    raise ValueError(f"GeneRax species-node identity is missing or ambiguous: {reconciled}")
                row["mapped_species_node"] = reconciled
            elif row["event"] == "D":
                covered = coverage[node]
                row["mapped_species_node"] = "" if covered is None else covered.name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.source.resolve() == args.output.resolve():
        parser.error("The species mapping output must be separate from the scoring tree")
    from nwkit.output_transaction import output_transaction

    mapping = SpeciesBranchMapping(args.source)
    with output_transaction([args.output]) as staged:
        mapping.write(staged[args.output])


if __name__ == "__main__":
    main()

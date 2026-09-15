"""Read species-tree branch identifiers without interpreting them as support."""

from Bio import Phylo


def read_species_tree(path):
    # Bio.Phylo normally converts numeric internal names into float confidence
    # values, losing leading zeros and precision. This documented parser mode
    # takes confidence from comments instead and preserves ALL node names.
    # Confidence is not used in either HGT context or transfer plotting.
    tree = Phylo.read(path, "newick", comments_are_confidence=True)
    aliases = set()
    for node in tree.find_clades():
        if node.name is None:
            if node.is_terminal():
                raise ValueError("Species tree contains an unnamed terminal")
            continue
        key = node.name.strip().replace(" ", "_")
        if not key or key in aliases:
            raise ValueError(f"Ambiguous species-tree label: {node.name!r}")
        aliases.add(key)
    return tree

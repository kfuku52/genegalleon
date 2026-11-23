"""Translate in-frame ORFs for protein searches."""

import re

import numpy as np
from Bio.Data.CodonTable import ambiguous_generic_by_id
from Bio.Seq import translate
from cdskit.codonutil import get_codon_table_components

CANONICAL_DNA = re.compile("[ACGT]+")


def make_orf_peptide_translator(code):
    """Use a vectorized lookup for canonical codons and Biopython for other cases."""
    table = ambiguous_generic_by_id[code]
    components = get_codon_table_components(code)
    forward, stops = components["forward_table"], components["stop_codons"]

    def translate_general(sequence):
        return translate(sequence, table=table, gap="-").removesuffix("*")

    # Keep Biopython's context-dependent stop handling and warnings.
    if not stops.isdisjoint(forward):
        return translate_general

    bases = "ACGT"
    encode = np.zeros(256, dtype=np.uint8)
    for index, base in enumerate(bases):
        encode[ord(base)] = index
    amino_acids = np.array([
        ord(forward.get(a + b + c, "*"))
        for a in bases for b in bases for c in bases
    ], dtype=np.uint8)

    def translate_orf(sequence):
        upper = sequence.upper()
        if len(upper) % 3 or not CANONICAL_DNA.fullmatch(upper):
            return translate_general(sequence)
        codons = encode[np.frombuffer(upper.encode("ascii"), dtype=np.uint8)].reshape(-1, 3)
        indices = 16 * codons[:, 0] + 4 * codons[:, 1] + codons[:, 2]
        return amino_acids[indices].tobytes().decode("ascii").removesuffix("*")

    return translate_orf

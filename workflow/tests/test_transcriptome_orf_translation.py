import importlib.util
import warnings
from itertools import product
from pathlib import Path

import pytest
from Bio.Data.CodonTable import TranslationError, ambiguous_generic_by_id
from Bio.Seq import translate

SUPPORT = Path(__file__).resolve().parents[1] / "support"
spec = importlib.util.spec_from_file_location(
    "transcriptome_orf_translation", SUPPORT / "transcriptome_orf_translation.py"
)
translator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(translator)


@pytest.mark.parametrize("code", sorted(ambiguous_generic_by_id))
def test_orf_translation_matches_biopython_across_genetic_codes(code):
    canonical = "".join("".join(codon) for codon in product("ACGT", repeat=3))
    fast_translate = translator.make_orf_peptide_translator(code)
    table = ambiguous_generic_by_id[code]

    def reference(sequence):
        return translate(sequence, table=table, gap="-").removesuffix("*")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for sequence in (canonical, canonical.lower(), "", "ATGGCNTANTAA", "ATG---TAA", "AUGGCUUAA", "atga"):
            assert fast_translate(sequence) == reference(sequence)


def test_orf_translation_preserves_internal_stops_and_rejects_invalid_codons():
    fast_translate = translator.make_orf_peptide_translator(1)
    assert fast_translate("ATGTAATGGTAA") == "M*W"
    with pytest.raises(TranslationError):
        fast_translate("ATGZZZTAA")
